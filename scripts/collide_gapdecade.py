#!/usr/bin/env python3
"""Run the skipped collide decade (D5 hole). Does not start itself — the bats do.

Every recorded results/collide_*.json starts at min_m_for_digits(l, cap+1)
(62-digit values for l<10, 102-digit for l>=10). The classifier already calls
the PREVIOUS decade past 2017 (61 / 101 digits). Gap per l:

    [collide_frontier_m(l), min_m_for_digits(l, cap+1) - 1]

Computed from the engine, not copied from a table. Settled (k,l) skipped.
(3,5) and (4,5) are above COLLIDE_HARD_SKIP_M and are not in any pack.

Packs (cheapest first inside each):
  tonight  l = 20,19,18,17,16,9,15,8,14,13   ~1.9 h serial, ~15-25 min at 8 workers
  day      l = 12,7,11,6                     ~56 h serial leftover, ~7 h at 8 workers
  l10      l = 10                            ~163 h serial, ~20 h at 8 workers

Skip any pair whose json already exists (rerun is a resume). Do NOT start a
second copy of the same pack.

    python scripts/collide_gapdecade.py --pack tonight --dry-run
    python scripts/collide_gapdecade.py --pack tonight --workers 8
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import singmaster_intersect as si  # noqa: E402

# Cheapest first. Same order as docs/collide-decade-queue.md.
TONIGHT_L = (20, 19, 18, 17, 16, 9, 15, 8, 14, 13)
DAY_L = (12, 7, 11, 6)
L10_L = (10,)
PACKS = {"tonight": TONIGHT_L, "day": DAY_L, "l10": L10_L}

RESULTS = ROOT / "results"


def gap_range(l: int) -> tuple[int, int]:
    lo = si.collide_frontier_m(l)
    old = si.min_m_for_digits(l, si.frontier_digits_for_l(l) + 1)
    hi = old - 1
    if hi < lo:
        raise RuntimeError(f"l={l}: empty gap ({lo}..{hi})")
    if lo > si.COLLIDE_HARD_SKIP_M:
        raise RuntimeError(f"l={l}: frontier {lo} exceeds COLLIDE_HARD_SKIP_M")
    return lo, hi


def pairs_for(l: int) -> list[int]:
    ks = [k for k in range(2, l) if (k, l) not in si.SETTLED_KL]
    if not ks:
        raise RuntimeError(f"l={l}: no unsettled k")
    return ks


def json_path(k: int, l: int) -> Path:
    return RESULTS / f"collide_gapdecade_k{k}_l{l}.json"


def plan(pack: str) -> list[tuple[int, int, int, int]]:
    out = []
    for l in PACKS[pack]:
        lo, hi = gap_range(l)
        for k in pairs_for(l):
            out.append((k, l, lo, hi))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--pack", required=True, choices=sorted(PACKS))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    jobs = plan(args.pack)
    n_skip = sum(1 for k, l, *_ in jobs if json_path(k, l).exists())
    print(
        f"=== collide gap-decade  pack={args.pack}  pairs={len(jobs)}  "
        f"already={n_skip}  workers={args.workers} ===",
        flush=True,
    )
    t0 = time.time()
    n_run = n_fail = 0
    for k, l, lo, hi in jobs:
        dest = json_path(k, l)
        span = hi - lo + 1
        cmd = [
            sys.executable,
            str(ROOT / "singmaster_intersect.py"),
            "collide",
            "--k", str(k),
            "--l", str(l),
            "--min-m", str(lo),
            "--max-m", str(hi),
            "--workers", str(args.workers),
            "--json_out", str(dest),
        ]
        if dest.exists():
            print(f"  skip existing {dest.name}", flush=True)
            continue
        print(
            f"  C(n,{k})=C(m,{l})  m={lo}..{hi}  ({span:,} rows)  -> {dest.name}",
            flush=True,
        )
        if args.dry_run:
            print("    " + " ".join(cmd), flush=True)
            continue
        rc = subprocess.call(cmd, cwd=str(ROOT))
        if rc != 0:
            print(f"  FAILED rc={rc}  {dest.name}", flush=True)
            n_fail += 1
            return rc
        if not dest.exists():
            print(f"  FAILED no json written  {dest.name}", flush=True)
            return 1
        n_run += 1
    dt = time.time() - t0
    print(
        f"=== pack={args.pack} done  ran={n_run} skipped={n_skip} "
        f"failed={n_fail}  {dt:.1f}s ===",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
