#!/usr/bin/env python3
"""Add a schema header to a pre-schema checkpoint, after verifying it earns one.

`check_checkpoint` refuses a checkpoint with no schema header, because resume
merges old records with new ones and a silent parameter change would produce a
certificate over columns that were never all tested the same way. That guard is
right, and it is also what strands a long run written by older code: the
records are fine, but nothing in the file says which parameters produced them.

This script closes that gap the honest way. It does NOT simply stamp a header
on. It re-derives the run's parameters from N and K, then checks that every
record in the file is consistent with them -- column ranges inside the phase
they are tagged with, pass indices inside the caps, primes in the right band --
and refuses if anything is off. A header is a claim about the records, so it is
only written once the records have been checked against it.

    python scripts/migrate_checkpoint.py --i 9 --chk PATH            # inspect
    python scripts/migrate_checkpoint.py --i 9 --chk PATH --write    # migrate

Safety:
  * never edits in place -- writes a new file and swaps it in
  * refuses a file that is still being appended to (a live run), because
    prepending under a writer loses records
  * keeps the original as PATH.pre-schema.bak
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bandii_kernel import (  # noqa: E402
    CHECKPOINT_SCHEMA,
    checkpoint_identity,
    make_fam,
)


def looks_live(path: Path, seconds: float = 6.0) -> bool:
    """Is something still appending to this file?"""
    a = path.stat().st_size
    time.sleep(seconds)
    return path.stat().st_size != a


def audit(path: Path, ident: dict) -> dict:
    """Check every record against the parameters the header would claim."""
    from family_sweep import z_cap

    fam_K, kmax, k_lo_z = ident["K"], ident["k_max"], ident["k_lo_z"]
    # The Z-jump's cap is the sweep's loop bound, not cap_z alone: every
    # small-k member runs to cap_z_small_k (15). Auditing against cap_z (12)
    # refused the legitimate z13/z14 records of every finished long run.
    cap_bii, cap_z = ident["cap_bii"], z_cap(ident)

    n = 0
    problems: list[str] = []
    tags: dict[str, int] = {}
    events: list[str] = []
    z_hi = 0

    with path.open(encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception as exc:
                problems.append(f"line {lineno}: not JSON ({exc})")
                continue

            if rec.get("event") == "schema":
                problems.append(f"line {lineno}: already has a schema header")
                continue
            if rec.get("event"):
                events.append(f"{rec['event']}:{rec.get('phase','')}")
                continue

            n += 1
            tag = rec.get("tag")
            if tag is None or "p" not in rec or "k_lo" not in rec:
                problems.append(f"line {lineno}: record is not resume-keyable")
                continue
            tags[tag] = tags.get(tag, 0) + 1
            k_lo, k_hi, p = int(rec["k_lo"]), int(rec["k_hi"]), int(rec["p"])

            if tag.startswith("bii"):
                lo, hi, cap = fam_K + 2, kmax, cap_bii
            elif tag.startswith("z"):
                lo, hi, cap = k_lo_z, fam_K - 1, cap_z
                z_hi = max(z_hi, k_hi)
            else:
                problems.append(f"line {lineno}: unknown tag {tag!r}")
                continue

            if not (lo <= k_lo <= k_hi <= hi):
                problems.append(
                    f"line {lineno}: {tag} columns [{k_lo},{k_hi}] outside "
                    f"phase range [{lo},{hi}]"
                )
            try:
                idx = int(tag[3:] if tag.startswith("bii") else tag[1:])
            except ValueError:
                problems.append(f"line {lineno}: tag {tag!r} has no pass index")
                continue
            if not (1 <= idx <= cap):
                problems.append(f"line {lineno}: {tag} pass {idx} outside cap {cap}")
            if p <= k_lo:
                problems.append(f"line {lineno}: p={p} not above k_lo={k_lo}")

    return {
        "records": n,
        "tags": tags,
        "events": events,
        "z_frontier": z_hi,
        "problems": problems,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--i", type=int, required=True)
    ap.add_argument("--chk", type=Path, required=True)
    ap.add_argument("--write", action="store_true",
                    help="actually migrate (default is inspect only)")
    ap.add_argument("--allow-live", action="store_true",
                    help="skip the still-being-written check (do not use)")
    args = ap.parse_args()

    if not args.chk.exists() or args.chk.stat().st_size == 0:
        print(f"  {args.chk}: missing or empty")
        return 1

    # THE identity family_sweep writes -- the same function, not a restated
    # dict -- so the header stamped here is the one the sweep checks against
    # on resume, n_chunks included. A restated copy drifted once (D3 added
    # n_chunks to the sweep; this file did not have it).
    from family_sweep import run_identity

    fam = make_fam(args.i)
    ident = run_identity(args.i)

    print(f"  checkpoint : {args.chk}")
    print(f"  size       : {args.chk.stat().st_size:,} bytes")
    print(f"  header would claim: {json.dumps(ident)}")
    print()

    rep = audit(args.chk, ident)
    print(f"  records         : {rep['records']:,}")
    print(f"  tags            : {rep['tags']}")
    print(f"  events          : {rep['events']}")
    if rep["z_frontier"]:
        span = fam.K - 1 - ident["k_lo_z"]
        done = rep["z_frontier"] - ident["k_lo_z"]
        print(f"  z frontier      : k={rep['z_frontier']:,} "
              f"({100 * done / span:.1f}% of round-1 span)")
    print()

    if rep["problems"]:
        print(f"  REFUSING: {len(rep['problems'])} inconsistencies")
        for p in rep["problems"][:15]:
            print(f"    {p}")
        if len(rep["problems"]) > 15:
            print(f"    ... and {len(rep['problems']) - 15} more")
        return 1

    print("  every record is consistent with those parameters.")

    if not args.write:
        print("\n  inspect only. re-run with --write to migrate.")
        return 0

    if not args.allow_live and looks_live(args.chk):
        print("\n  REFUSING: the file is still growing -- a run is writing to it.")
        print("  Stop the run first. Prepending under a live writer loses records.")
        return 1

    tmp = args.chk.with_suffix(args.chk.suffix + ".migrating")
    bak = args.chk.with_suffix(args.chk.suffix + ".pre-schema.bak")

    # Binary throughout. These files are written on Windows and carry CRLF; a
    # text-mode copy silently rewrites every one of them to LF, which turns
    # "prepend one line" into "rewrite 100 MB of a certificate-bearing
    # artifact". The body must come through byte for byte, so the header
    # matches the line ending already in use rather than imposing one.
    with args.chk.open("rb") as fh:
        head = fh.read(65536)
    nl = b"\r\n" if b"\r\n" in head else b"\n"

    with tmp.open("wb") as out:
        out.write(json.dumps(checkpoint_identity(**ident)).encode("utf-8") + nl)
        with args.chk.open("rb") as fh:
            while chunk := fh.read(1 << 20):
                out.write(chunk)
        out.flush()
        os.fsync(out.fileno())

    orig_size = args.chk.stat().st_size
    grew = tmp.stat().st_size - orig_size
    if grew <= 0:
        tmp.unlink()
        print(f"\n  REFUSING: migrated file did not grow ({grew:+,} bytes)")
        return 1

    os.replace(args.chk, bak)
    os.replace(tmp, args.chk)
    print(f"\n  migrated to schema {CHECKPOINT_SCHEMA}")
    print(f"  original kept at {bak.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
