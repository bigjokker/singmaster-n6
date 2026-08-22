#!/usr/bin/env python3
"""Where does a Z-jump job's time actually go, on the CURRENT kernel?

Grew out of the probe that first timed the production call rather than a proxy.
It pins one fact that decides what is worth optimising next:

    for a typical i=9 Z-jump job, r(p) is essentially the ENTIRE cost.

Measured 2026-08-22 over 120 real recorded jobs: r_two_digit_delta accounted
for 99.6% of scan_columns_general at 1 worker and 102% at 8 (i.e. all of it,
within noise). A median job scans 12 columns with g=12 -- about 72 half-scan
elements, microseconds of numpy -- while r(p) costs 0.05-0.43 s because the
Z-jump's remainders are generic: p is near k, not near N, so
min(k0, n0-k0, p-n0-1) stays O(p) and the Lucas product is a long loop.

WHAT THIS TEST MUST NOT DO is pin the 0.955 s median job time recorded in
results/i9_sweep.jsonl. That number came from the LIVE kernel, whose
scan_columns_general has no USE_WINDOWED_SCAN and builds a p-sized fact_table
on every job -- about a second at p ~ 1e7. The playground path replays the same
job in 0.087 s. Recorded seconds are not a budget for this kernel; the 11x gap
between them is Q2 not being wired on the machine that ran i=9.

So this asserts a SHARE, not a duration -- shares are machine-independent and
survive a faster box. It exists to fail loudly if someone makes r(p) cheap
(good news, update the threshold) or makes the scan dominant (a regression).

Jobs are rebuilt from the live-prime ladder rather than read from the 219 MB
checkpoint, so the test has no large-file dependency.

    python scripts/test_rp_cost.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import bandii_kernel as bk  # noqa: E402
from bandii_kernel import cells, first_live_after, live_intervals, make_fam  # noqa: E402

ok: list[str] = []
errors: list[str] = []

I_TEST = 9
N_JOBS = 16
# Measured 99.6%. A wide margin: this is a share of a real timing, and the
# point is to catch a REGIME change, not to police a few percent.
MIN_RP_SHARE = 0.70


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


def production_jobs(fam, n: int):
    """Real round-1 Z-jump buckets: contiguous columns sharing a live prime.

    bucket(p) = [previous live prime, p-1], the lower endpoint INCLUSIVE
    because first_live_after returns the first prime STRICTLY greater than x.
    """
    ivs = live_intervals(cells(fam), fam)
    out = []
    for anchor in range(2_000_000, 26_000_000, 1_500_000):
        prev = first_live_after(anchor, ivs, fam.D)
        p = first_live_after(prev, ivs, fam.D) if prev else None
        if not p:
            continue
        r = bk.r_of(p, N=fam.N, K=fam.K)
        if not r:                      # r(p)=0 kills nothing; the scan refuses it
            continue
        ks = list(range(prev, p))
        if 1 <= len(ks) <= 4096:
            out.append((p, ks))
        if len(out) >= n:
            break
    return out


def main() -> int:
    fam = make_fam(I_TEST)
    jobs = production_jobs(fam, N_JOBS)
    expect(len(jobs) >= 4, f"rebuilt {len(jobs)} production-shaped Z-jump jobs")
    if len(jobs) < 4:
        print("  RESULT FAIL (no jobs)")
        return 1

    t_full = t_r = 0.0
    cols = gs = 0
    for p, ks in jobs:
        t = time.perf_counter()
        bk.scan_columns_general(p, ks, N=fam.N, K=fam.K)
        t_full += time.perf_counter() - t
        t = time.perf_counter()
        bk.r_two_digit_delta(p, N=fam.N, K=fam.K)
        t_r += time.perf_counter() - t
        cols += len(ks)
        gs += p - ks[0]

    share = t_r / t_full if t_full else 0.0
    print(f"  {len(jobs)} jobs, {cols} columns, mean g {gs / len(jobs):,.0f}")
    print(f"  scan_columns_general {t_full * 1000:8.1f} ms")
    print(f"  r_two_digit_delta    {t_r * 1000:8.1f} ms   ({100 * share:.1f}%)")

    expect(share >= MIN_RP_SHARE,
           f"r(p) is the dominant cost of a Z-jump job "
           f"({100 * share:.1f}% >= {100 * MIN_RP_SHARE:.0f}% required)")
    expect(t_full > 0, "the production call was actually timed")

    # The scan itself must be the cheap half. If this ever inverts, the regime
    # changed and the optimisation target moved.
    expect(t_full - t_r < t_full * 0.5,
           f"the image scan is the cheap half "
           f"({100 * (t_full - t_r) / t_full:.1f}% of the job)")

    print("\n=== r(p) COST TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
