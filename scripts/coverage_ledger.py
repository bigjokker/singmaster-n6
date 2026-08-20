#!/usr/bin/env python3
"""Q21: check that the coverage claims actually union to the full range.

The claim behind N(m) = 6 is that EVERY column k in [2, k_max], except the two
family columns {K, K+1}, is impossible. Until now no artifact asserted that.
Each one self-reported its own slice -- the witness table one range, the exact
scan another, the modular run a third -- with three different kinds of
evidence, and the union was assembled by a human reading three files.
`witness.coverage()` checks completeness of the range the sweep CLAIMED, which
is a different statement.

That is the same shape as the false-`clean` bugs this project has already had:
a local check passing while the global claim goes unchecked.

This script states the global claim and checks it:

    for each member i:
        witnessed(i)  ==  [2, k_max(i)]  \\  {K, K+1}

exactly -- no missing column, no extra column, and k_max recomputed from N and
K rather than read back from the file that is being checked.

    python scripts/coverage_ledger.py
    python scripts/coverage_ledger.py --json_out results/coverage_ledger.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import witness as W  # noqa: E402
from bandii_kernel import kmax_of, make_fam  # noqa: E402


def audit_member(i: int, path: Path) -> dict:
    ks, ps, meta = W.load(path)

    # Recompute the target range from N and K. Reading k_max out of the file
    # would make the check circular -- the file is the thing under audit.
    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    if (meta["N"], meta["K"]) != (fam.N, fam.K):
        return {"i": i, "ok": False,
                "error": f"file says N,K = {meta['N']},{meta['K']}; "
                         f"the family says {fam.N},{fam.K}"}

    expected = set(range(2, kmax + 1)) - {fam.K, fam.K + 1}
    got = {int(v) for v in ks}
    missing = sorted(expected - got)
    extra = sorted(got - expected)
    return {
        "i": i,
        "N": fam.N,
        "K": fam.K,
        "k_max": kmax,
        "k_max_source": "recomputed from N,K",
        "n_expected": len(expected),
        "n_witnessed": len(got),
        "n_missing": len(missing),
        "missing_sample": missing[:20],
        "n_extra": len(extra),
        "extra_sample": extra[:20],
        "family_columns_correctly_absent": not ({fam.K, fam.K + 1} & got),
        "sha256": meta["sha256"],
        "ok": not missing and not extra and not ({fam.K, fam.K + 1} & got),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    rows = []
    for f in sorted(args.results.glob("i*_witness.npz")):
        try:
            i = int(f.name.split("_")[0][1:])
        except ValueError:
            continue
        try:
            rows.append(audit_member(i, f))
        except Exception as exc:
            rows.append({"i": i, "ok": False, "error": repr(exc)})
    if not rows:
        print("  no witness tables found", flush=True)
        return 1
    rows.sort(key=lambda r: r["i"])

    print("  claim: witnessed(i) == [2, k_max] \\ {K, K+1}, exactly\n")
    print(f"  {'i':>3} {'k_max':>10} {'expected':>11} {'witnessed':>11} "
          f"{'missing':>8} {'extra':>7}  verdict")
    for r in rows:
        if "error" in r:
            print(f"  {r['i']:>3} {'-':>10} {'-':>11} {'-':>11} {'-':>8} {'-':>7}  "
                  f"ERROR {r['error'][:40]}")
            continue
        print(f"  {r['i']:>3} {r['k_max']:>10,} {r['n_expected']:>11,} "
              f"{r['n_witnessed']:>11,} {r['n_missing']:>8} {r['n_extra']:>7}  "
              f"{'COMPLETE' if r['ok'] else 'INCOMPLETE'}")
        if r["missing_sample"]:
            print(f"        missing: {r['missing_sample']}")
        if r["extra_sample"]:
            print(f"        extra:   {r['extra_sample']}")

    bad = [r for r in rows if not r.get("ok")]
    total = sum(r.get("n_witnessed", 0) for r in rows)
    print(f"\n  members audited {len(rows)}   incomplete {len(bad)}   "
          f"columns witnessed {total:,}")
    print("  RESULT", "ALL COMPLETE" if not bad else "GAPS PRESENT")
    print("\n  Note: completeness is about COVERAGE. That every witness is also")
    print("  VALID is a separate check -- scripts/witness.py verify.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(
            {"search": "coverage_ledger", "n_members": len(rows),
             "n_incomplete": len(bad), "total_columns": total, "members": rows},
            indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
