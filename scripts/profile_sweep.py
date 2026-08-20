#!/usr/bin/env python3
"""Where does a sweep's time actually go?

This exists because the question kept being answered by reasoning, and the
reasoning kept being wrong. Four in a row on this project:

  * r(p) was called the Band I bottleneck, then retracted, then un-retracted
    -- and the true figure was 13-15x, not the 26x the ratio suggested;
  * the "table-free" scan was ranked as the Z-jump fix and turned out to be
    2-26x SLOWER there than the table it was meant to replace;
  * the half-scan measured 2.2x in the kernel but 1.14x end-to-end, because
    the Z-jump is table-dominated and Band II is scan-dominated;
  * a folded engine scan measured 4-7x on a path that is never hot.

So: measure a representative slice of each phase, report the split, and
project the phase total. Costs seconds, not hours, because it samples rather
than sweeping.

    python scripts/profile_sweep.py --i 7
    python scripts/profile_sweep.py --i 9 --sample 200

Reports, per phase, the share going to the factorial table, the image scan,
and r(p) -- which is the number that decides which optimisation is worth
writing. A phase where the table is 80% will not be helped by a faster scan.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2  # noqa: E402

from bandii_kernel import (  # noqa: E402
    cells,
    fact_table,
    first_live_after,
    first_primes_above,
    kmax_of,
    live_intervals,
    make_fam,
    r_closed,
    r_two_digit_delta,
    scan_ks_full,
    scan_ks_half,
)


def _t(fn, *a, reps: int = 1):
    fn(*a)
    t = perf_counter()
    for _ in range(reps):
        out = fn(*a)
    return (perf_counter() - t) / reps, out


def _live_prime_count(ivs, k_lo: int, k_hi: int) -> float:
    """Live primes in [k_lo, k_hi], by integrating pi(x) over live intervals."""
    import math

    def pi(x):
        if x < 3:
            return 0.0
        L = math.log(x)
        return x / L * (1 + 1 / L + 2 / L**2)

    return sum(pi(min(hi, k_hi)) - pi(max(lo, k_lo))
               for lo, hi in ivs if max(lo, k_lo) < min(hi, k_hi))


def profile_bandii(fam, kmax: int, sample: int) -> dict:
    """Band II: one prime, every column. Table built once per chunk."""
    p = first_primes_above(fam.N2, fam.D, kmax, n=1)[0]
    t_tab, F = _t(fact_table, p)
    # r_closed/r_two_digit_delta default N,K to the i=8 constants; passing the
    # member's own values is not optional.
    rp = r_closed(p, N=fam.N, K=fam.K)
    t_r, _ = _t(lambda q: r_closed(q, N=fam.N, K=fam.K), p)

    lo, hi = fam.K + 2, kmax
    step = max(1, (hi - lo) // sample)
    ks = list(range(lo, hi + 1, step))
    t_half, _ = _t(scan_ks_half, F, p, rp, ks)
    t_full, _ = _t(scan_ks_full, F, p, rp, ks)

    n_cols = hi - lo + 1
    scan_total = t_half / len(ks) * n_cols
    return {
        "phase": "bandii",
        "p": p,
        "n_columns": n_cols,
        "n_primes": 1,
        "table_s": t_tab,
        "scan_s": scan_total,
        "r_s": t_r,
        "half_vs_full": t_full / t_half if t_half else None,
        "sampled": len(ks),
    }


def profile_zjump(fam, sample: int, k_lo: int) -> dict:
    """Z-jump: a prime per column, so the table is rebuilt constantly."""
    ivs = live_intervals(cells(fam), fam)
    hi = fam.K - 1
    if k_lo >= hi:
        return {"phase": "zjump", "n_columns": 0, "note": "empty"}
    step = max(1, (hi - k_lo) // sample)
    ks = list(range(k_lo, hi, step))

    t = perf_counter()
    ps = [first_live_after(k, ivs, fam.D) for k in ks]
    t_assign_per = (perf_counter() - t) / len(ks)
    ps = [q for q in ps if q]

    # Distinct primes = the LIVE primes in the range, not one per column:
    # consecutive columns share a prime. Counting them from a contiguous
    # window near k_lo is badly biased -- at small k the live primes are far
    # apart -- so integrate prime density over the live intervals instead.
    # That form recovers i=8's recorded 124,830 primes closely; the window
    # method gave 239,188, nearly 2x high.
    n_cols = hi - k_lo
    est_primes = _live_prime_count(ivs, k_lo, fam.K)

    mid = ps[len(ps) // 2]
    t_tab, F = _t(fact_table, mid)
    t_r, _ = _t(lambda q: r_two_digit_delta(q, N=fam.N, K=fam.K), mid)
    r = r_two_digit_delta(mid, N=fam.N, K=fam.K) or 1
    probe = [k for k in ks if k < mid][-8:] or [ks[0]]
    t_scan, _ = _t(scan_ks_half, F, mid, r, probe)

    return {
        "phase": "zjump",
        "n_columns": n_cols,
        "n_primes_est": est_primes,
        "cols_per_prime": n_cols / est_primes if est_primes else None,
        "table_s": t_tab * est_primes,
        "scan_s": t_scan / len(probe) * n_cols,
        "r_s": t_r * est_primes,
        "assign_s": t_assign_per * n_cols,
        "sample_p": mid,
    }


def report(rows: list[dict]) -> None:
    for r in rows:
        if not r.get("n_columns"):
            print(f"  {r['phase']}: empty")
            continue
        parts = {k[:-2]: r[k] for k in ("table_s", "scan_s", "r_s", "assign_s")
                 if r.get(k)}
        tot = sum(parts.values())
        print(f"\n  === {r['phase']} ===")
        print(f"    columns {r['n_columns']:,}"
              + (f"   distinct primes ~{r['n_primes_est']:,.0f}"
                 f"   ({r['cols_per_prime']:.1f} cols/prime)"
                 if r.get("n_primes_est") else "   1 prime per pass"))
        for name, v in sorted(parts.items(), key=lambda x: -x[1]):
            bar = "#" * int(40 * v / tot)
            print(f"    {name:8s} {v/3600:9.2f} core-h  {100*v/tot:5.1f}%  {bar}")
        print(f"    {'TOTAL':8s} {tot/3600:9.2f} core-h  "
              f"({tot/3600/8:.2f} h on 8 workers)")
        if r.get("half_vs_full"):
            print(f"    half-scan is {r['half_vs_full']:.2f}x the full scan here")
        dom = max(parts, key=parts.get)
        print(f"    -> {dom}-dominated: optimise that, not the others")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--i", type=int, required=True)
    ap.add_argument("--sample", type=int, default=400,
                    help="columns sampled per phase (cost, not accuracy, scales)")
    ap.add_argument("--k-lo", type=int, default=None,
                    help="Z-jump start; defaults to the member's K_EXACT+1")
    args = ap.parse_args()

    from family_sweep import K_EXACT

    fam = make_fam(args.i)
    kmax, logm = kmax_of(fam)
    k_lo = args.k_lo if args.k_lo is not None else K_EXACT.get(args.i, 2) + 1
    print(f"  i={args.i}  N={fam.N:,}  K={fam.K:,}  k_max={kmax:,}  "
          f"log10 m={logm:,.0f}")
    print(f"  sampling {args.sample} columns per phase", flush=True)
    rows = [profile_bandii(fam, kmax, args.sample),
            profile_zjump(fam, args.sample, k_lo)]
    report(rows)
    print("\n  (projections are sample x column-count; they are good to about a factor of 2;"
          "\n   the SPLIT is the reliable part, and the split is what to act on)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
