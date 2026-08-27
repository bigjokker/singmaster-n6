#!/usr/bin/env python3
"""Where does a sweep's time actually go?

This exists because the question kept being answered by reasoning, and the
reasoning kept being wrong. Then the tool itself was wrong, in both probes,
and was believed for two days.

REWRITTEN 2026-08-22 after `work_census.py` produced exact ground truth:

  profile_zjump was 8.5x-156x high, and its answer MOVED with --sample:
    --sample 20/100/400 -> 1183.75 / 270.25 / 64.50 core-h at i=9, against an
    exact 7.61. It timed 8 columns that were spread ~707,000 apart against
    `mid`, the first-live-prime of a DIFFERENT column, so the probe columns
    carried g in the millions where production's median g is 29.

  profile_bandii was a stable 1.7x high, which is why it looked trustworthy.
    `blk` started at lo = K+2 -- the smallest k, hence the LARGEST g -- so the
    marginal per-column cost was measured at g = 8,740,407 and applied to a
    population whose mean g is 5,121,476. Predicted 1.707x, measured 1.72x.

The fix is not a better sample. It is to stop sampling the thing that can be
computed exactly, and to separate the two quantities the old probes conflated:

  OPS   pure geometry: how many half-scan elements the phase costs.
        EXACT, closed form, no sampling and no timing. At a fixed prime the
        columns [a, b] contribute sum_{g} ceil(g/2) over a contiguous g range,
        and sum_{g=1}^{G} ceil(g/2) = floor((G+1)^2/4). Band II is one such
        range; the Z-jump is one per live-prime bucket.

  RATE  pure machine: how many elements per second this box does.
        Measured on production-shaped CONTIGUOUS blocks at several g scales,
        never at one end of a band, with the spread reported rather than hidden.

Consequently `--sample` no longer controls the answer. It controls the
precision of the modelled multi-pass correction (a ~1.07-1.13 factor) and
nothing else. `--check` gates the whole thing against `work_census.py`.

Scan work and wall clock are reported as TWO numbers, and the gap between
them is not constant: a faster scan is not a faster job, and no wall clock is
extrapolated to an unrun member.

RESOLVED 2026-08-25: the "6.11x at i=7 and 10.27x at i=8" this file used to
quote was a ratio against work_census.SCAN_RATE while that constant was still
the pre-Granlund-Montgomery 1.865e8. SCAN_RATE has since been re-measured on
the GM kernel (8.80e8 census-ops/s), so the denominator is current and the
old pair is simply superseded, not merely suspect. Note i=7's recorded wall
still predates GM (results/i7_sweep.json was never re-run), so its wall/scan
ratio compares a GM scan estimate against a pre-GM wall and reads high; i=8's
wall was regenerated under GM (7,508 s -> 174 s) and is the comparable one.
The rate this profiler MEASURES at run time is current; the fixed reference
it is compared against is not.

    python scripts/profile_sweep.py --i 9
    python scripts/profile_sweep.py --i 8 --check
    python scripts/profile_sweep.py --i 10
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2  # noqa: E402,F401

import sizelaw as S  # noqa: E402
from bandii_kernel import (  # noqa: E402
    cells,
    first_live_after,
    first_primes_above,
    kmax_of,
    live_intervals,
    make_fam,
    r_of,
    scan_ks_windowed,
)
from family_sweep import CAP_BII, CAP_Z, K_EXACT  # noqa: E402

# Elements the rate probe aims to touch per timed block: big enough to swamp
# call overhead, small enough that the probe stays in the "seconds" budget.
RATE_BLOCK_ELEMS = 4_000_000
# An anchor smaller than this cannot be timed meaningfully -- call overhead
# would swamp the scan and the 'rate' would be a measurement of Python.
MIN_RATE_ELEMS = 200_000
# The element rate genuinely varies with g -- short scans pay call overhead,
# long ones are memory-bandwidth bound -- so a spread across regimes is
# EXPECTED and the ops->time conversion is published as a band, never a point.
# This ceiling catches a machine that is not behaving like one machine at all
# (thermal throttling, another job competing), which is a different problem.
RATE_SPREAD_MAX = 6.0


def check(cond: bool, msg: str) -> None:
    """Guard. Not assert: `python -O` strips assert."""
    if not cond:
        raise RuntimeError(f"profile pre-flight failed: {msg}")


def _t(fn, *a, reps: int = 1):
    fn(*a)
    t = perf_counter()
    for _ in range(reps):
        out = fn(*a)
    return (perf_counter() - t) / reps, out


def _S(G):
    """sum_{g=1}^{G} ceil(g/2) = floor((G+1)^2/4). Exact, vectorised."""
    G = np.asarray(G, dtype=np.int64)
    return np.where(G <= 0, np.int64(0), ((G + 1) ** 2) // 4)


def _ops_range(p: int, k_lo: int, k_hi: int) -> int:
    """Exact half-scan elements for columns [k_lo, k_hi] scanned against p.

    g = p - k runs from p-k_hi up to p-k_lo, so the total is the difference of
    two triangular-ish sums. No loop over columns, no sampling.
    """
    if k_hi < k_lo:
        return 0
    return int(_S(p - k_lo)) - int(_S(p - k_hi - 1))


def zjump_buckets(fam, k_lo: int, k_hi: int):
    """Round-1 buckets as (prime, k_lo, k_hi), ascending.

    bucket(p) = [previous live prime, p-1]. The lower endpoint is INCLUSIVE of
    the previous live prime, because `first_live_after` returns the first prime
    STRICTLY greater than x. Taking prev+1 silently drops ~0.5% of columns.
    """
    ivs = live_intervals(cells(fam), fam)
    out = []
    lo = k_lo
    p = first_live_after(k_lo, ivs, fam.D)
    while p is not None and lo <= k_hi:
        hi = min(k_hi, p - 1)
        if hi >= lo:
            out.append((int(p), int(lo), int(hi)))
        lo = int(p)
        p = first_live_after(p, ivs, fam.D)
    return out, ivs


def measure_rate(anchors, reps: int = 2) -> dict:
    """Per-element rate for the COMPARISON, with the window build removed.

    scan_ks_windowed builds two O(g) factorial windows per call, in PURE
    PYTHON. Timing one block and dividing by its element count therefore
    charges the window as if it were scan work -- with a 2-column anchor at
    i=9 that dragged the apparent rate from ~2e8 down to 1.0e7, a 20x error
    in the wrong direction.

    So take a MARGINAL measurement: time n columns and 2n columns against the
    same prime. The window is built once either way, so it cancels in the
    difference and what remains is the per-element comparison cost -- which is
    exactly what the ops count models. This is the one idea worth keeping from
    the old profile_bandii; its defect was measuring at k=K+2 where g is
    largest, not the marginal trick itself.
    """
    rates = []
    for p, r, k0, k1 in anchors:
        g = p - k0
        if g <= 0 or k1 < k0:
            continue
        n1 = int(max(2, min(2048, RATE_BLOCK_ELEMS // max(1, g // 2))))
        n2 = 2 * n1
        if k0 + n2 - 1 > k1:                      # not enough columns share p
            continue
        e1 = _ops_range(p, k0, k0 + n1 - 1)
        e2 = _ops_range(p, k0, k0 + n2 - 1)
        if e2 - e1 < MIN_RATE_ELEMS:
            continue
        t1, _ = _t(scan_ks_windowed, p, r, list(range(k0, k0 + n1)), reps=reps)
        t2, _ = _t(scan_ks_windowed, p, r, list(range(k0, k0 + n2)), reps=reps)
        dt = t2 - t1
        if dt <= 0:
            continue
        rates.append({"p": p, "k_lo": k0, "g": g, "cols": f"{n1}->{n2}",
                      "elems": e2 - e1, "sec": dt, "rate": (e2 - e1) / dt})
    check(bool(rates), "no rate anchor had enough columns sharing one prime")
    vals = sorted(x["rate"] for x in rates)
    spread = vals[-1] / vals[0]
    return {"rate": float(np.median(vals)), "rate_lo": float(vals[0]),
            "rate_hi": float(vals[-1]), "spread": float(spread),
            "anchors": rates,
            "trustworthy": bool(spread <= RATE_SPREAD_MAX)}


def _bii_primes(fam, kmax: int) -> list[int]:
    return [int(x) for x in first_primes_above(fam.N2, fam.D, kmax, n=CAP_BII)]


def multipass_bandii(fam, kmax: int, primes, sample: int) -> dict:
    """Modelled (all passes)/(pass 1) ops ratio for Band II, from the size law.

    Band II g = p1 - k is LINEAR in k, so a uniform k-sample is unbiased here;
    there is no heavy tail to importance-weight against. Terminates when the
    expected survivor weight dies, NOT at CAP_BII -- production recorded nine
    passes at i=9, not fourteen.
    """
    lo, hi = fam.K + 2, kmax
    n = min(max(sample, 8), hi - lo + 1)
    ks = np.linspace(lo, hi, n).astype(np.int64)
    w1 = np.ceil((primes[0] - ks) / 2.0)
    alive = np.ones(n, dtype=np.float64)
    num = np.zeros(n, dtype=np.float64)
    passes = 0
    for p in primes:
        w = np.ceil((p - ks) / 2.0)
        num += alive * w
        passes += 1
        alive = alive * S.survival_vec(p - ks, p, ks)
        if float(alive.sum()) < 1e-9:
            break
    return {"factor": float(num.sum() / w1.sum()), "passes_modelled": passes,
            "sampled": int(n)}


def multipass_zjump(fam, buckets, ivs, sample: int, cap: int) -> dict:
    """Modelled (all rounds)/(round 1) ops ratio for the Z-jump.

    Z-jump g is heavy-tailed (i=9: median 29, mean 307,398), so buckets are
    sampled with probability proportional to their round-1 ops. A uniform
    sample would be dominated by the thousands of tiny buckets that carry
    almost none of the work.
    """
    ops = np.array([_ops_range(p, a, b) for p, a, b in buckets], dtype=np.float64)
    tot = ops.sum()
    check(tot > 0, "Z-jump round-1 ops is zero")
    n = min(max(sample, 8), len(buckets))
    rng = np.random.default_rng(12345)
    idx = rng.choice(len(buckets), size=n, replace=False, p=ops / tot)
    num_tot = den_tot = 0.0
    for j in idx:
        p1, a, b = buckets[j]
        # Columns WITHIN a bucket are not interchangeable: ops per column is
        # ceil((p-k)/2), so the bucket's work and its survival are both
        # concentrated at the small-k end. A single midpoint representative
        # measured 1.0037 against a true 1.0722 -- it missed the tail entirely.
        m = min(32, b - a + 1)
        ks = np.unique(np.linspace(a, b, m).astype(np.int64))
        den_tot += float(np.ceil((p1 - ks) / 2.0).sum()) * (b - a + 1) / len(ks)
        alive = np.ones(len(ks), dtype=np.float64)
        num = np.zeros(len(ks), dtype=np.float64)
        q = p1
        for _ in range(cap):
            if q is None:
                break
            num += alive * np.ceil((q - ks) / 2.0)
            alive = alive * S.survival_vec(q - ks, q, ks)
            if float(alive.max()) < 1e-9:
                break
            q = first_live_after(q, ivs, fam.D)
        num_tot += float(num.sum()) * (b - a + 1) / len(ks)
    check(den_tot > 0, "no usable Z-jump multipass sample")
    return {"factor": float(num_tot / den_tot), "sampled": int(n),
            "ops_weighted": True}


def profile_bandii(fam, kmax: int, sample: int) -> dict:
    primes = _bii_primes(fam, kmax)
    check(bool(primes), "no live Band II prime in (N/2, d]")
    p1 = primes[0]
    ops1 = _ops_range(p1, fam.K + 2, kmax)
    mp = multipass_bandii(fam, kmax, primes, sample)
    r1 = r_of(p1, N=fam.N, K=fam.K)
    check(bool(r1), f"r(p)=0 at Band II p1={p1}; that prime kills nothing")
    # anchors spanning the real g range, not one end of it
    lo, hi = fam.K + 2, kmax
    anchors = [(p1, r1, lo, hi), (p1, r1, (lo + hi) // 2, hi),
               (p1, r1, max(lo, hi - 8192), hi)]
    return {"phase": "bandii", "n_columns": kmax - fam.K - 1, "p": p1,
            "n_primes_available": len(primes),
            "ops_round1": ops1, "multipass": mp,
            "ops": int(round(ops1 * mp["factor"])), "anchors": anchors}


def profile_zjump(fam, sample: int, k_lo: int) -> dict:
    k_hi = fam.K - 1
    if k_lo > k_hi:
        return {"phase": "zjump", "n_columns": 0, "ops": 0, "note": "empty"}
    buckets, ivs = zjump_buckets(fam, k_lo, k_hi)
    check(bool(buckets), "no Z-jump bucket")
    covered = sum(b - a + 1 for _, a, b in buckets)
    check(covered == k_hi - k_lo + 1,
          f"buckets cover {covered} columns, band has {k_hi - k_lo + 1}")
    ops1 = sum(_ops_range(p, a, b) for p, a, b in buckets)
    cap = CAP_Z
    mp = multipass_zjump(fam, buckets, ivs, sample, cap)
    # rate anchors: biggest-ops bucket, median, and a small one
    # Rate anchors come from the buckets that carry the work. A bucket at the
    # 10th percentile of ops has g in single digits: timing it measures call
    # overhead, and mixing it in produced a 75,534x spread.
    order = sorted(range(len(buckets)), key=lambda j: _ops_range(*buckets[j]))
    picks = [order[-1], order[-max(2, len(order) // 100)],
             order[-max(3, len(order) // 20)]]
    anchors = []
    for j in dict.fromkeys(picks):
        p, a, b = buckets[j]
        r = r_of(p, N=fam.N, K=fam.K)
        if r:
            anchors.append((p, r, a, b))
    return {"phase": "zjump", "n_columns": k_hi - k_lo + 1,
            "n_buckets": len(buckets), "cols_per_bucket": covered / len(buckets),
            "ops_round1": ops1, "multipass": mp,
            "ops": int(round(ops1 * mp["factor"])), "anchors": anchors}


def profile(i: int, sample: int, k_lo: int | None = None) -> dict:
    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    k_lo = k_lo if k_lo is not None else K_EXACT.get(i, 2) + 1
    rows = [profile_bandii(fam, kmax, sample), profile_zjump(fam, sample, k_lo)]
    anchors = [a for r in rows for a in r.get("anchors", [])]
    rate = measure_rate(anchors)
    total = sum(r["ops"] for r in rows)
    out = {"i": i, "N": fam.N, "K": fam.K, "k_max": kmax, "sample": sample,
           "phases": rows, "rate": rate, "ops_total": total,
           "scan_core_h": total / rate["rate"] / 3600.0,
           "scan_core_h_lo": total / rate["rate_hi"] / 3600.0,
           "scan_core_h_hi": total / rate["rate_lo"] / 3600.0}
    rec = ROOT / "results" / f"i{i}_sweep.json"
    if rec.exists():
        d = json.loads(rec.read_text())
        partial = not (d.get("phases") or {}).get("bandii") and d.get("n_bii")
        if d.get("seconds") and d.get("workers") and not partial:
            out["recorded_wall_core_h"] = d["seconds"] * d["workers"] / 3600.0
        elif partial:
            out["record_partial"] = True
    return out


def report(r: dict) -> None:
    print(f"\n  i={r['i']}  N={r['N']:,}  K={r['K']:,}  k_max={r['k_max']:,}"
          f"   (--sample {r['sample']})")
    for ph in r["phases"]:
        if not ph.get("n_columns"):
            print(f"    {ph['phase']}: empty")
            continue
        mp = ph["multipass"]
        extra = (f"  {ph['n_buckets']:,} buckets, {ph['cols_per_bucket']:.1f} cols each"
                 if "n_buckets" in ph else f"  1 prime per pass")
        print(f"    === {ph['phase']} ===  {ph['n_columns']:,} columns{extra}")
        print(f"      ops round 1  {ph['ops_round1']:>20,d}   EXACT (closed form)")
        print(f"      x multipass  {mp['factor']:>20.6f}   modelled, "
              f"{mp['sampled']:,} sampled")
        print(f"      ops total    {ph['ops']:>20,d}")
    rate = r["rate"]
    flag = "" if rate["trustworthy"] else "   *** MACHINE NOT STABLE, DO NOT TRUST ***"
    print(f"    rate {rate['rate']:.4g} elem/s over {len(rate['anchors'])} anchors, "
          f"{rate['spread']:.2f}x across g regimes{flag}")
    print(f"    OPS        {r['ops_total']:,} elements"
          f"   (round 1 EXACT; x multipass MODELLED)")
    print(f"    SCAN WORK  {r['scan_core_h_lo']:.4f} - {r['scan_core_h_hi']:.4f} core-h"
          f"   (band, not a point: the rate depends on g)")
    if "recorded_wall_core_h" in r:
        w = r["recorded_wall_core_h"]
        print(f"    WALL CLOCK {w:.4f} core-h recorded  -> wall/scan "
              f"{w / r['scan_core_h_hi']:.2f}x - {w / r['scan_core_h_lo']:.2f}x")
        print("    These are two different numbers, and their ratio is not")
        print("    constant -- do not extrapolate it to an unrun member. (The")
        print("    once-quoted 6.11x/10.27x used a pre-Granlund-Montgomery")
        print("    SCAN_RATE; that constant was re-measured 2026-08-25.)")
    elif r.get("record_partial"):
        print("    WALL CLOCK not comparable: this record is partial (its `seconds`")
        print("    covers only a resumed leg). No ratio is printed.")
    else:
        print("    WALL CLOCK not recorded for this member; none is projected.")


def do_check(r: dict) -> bool:
    """Gate the profile against work_census.py, the exact ground truth."""
    import work_census as W
    try:
        c = W.census(r["i"])
    except Exception as exc:
        print(f"    CHECK skipped: census unavailable ({exc})")
        return True
    ok = True
    print(f"\n    --- check against work_census (i={r['i']}) ---")
    per = {p["phase"]: p for p in r["phases"]}
    for name in ("bandii", "zjump"):
        got, want = per[name]["ops"], c["phases"][name]["ops_exact"]
        err = abs(got - want) / want if want else 0.0
        good = err <= 0.02
        ok &= good
        print(f"    {'ok  ' if good else 'FAIL'} {name:7s} ops {got:>20,d} "
              f"vs census {want:>20,d}   {100 * err:6.2f}%")
    got, want = r["ops_total"], c["ops_exact"]
    err = abs(got - want) / want
    good = err <= 0.02
    ok &= good
    print(f"    {'ok  ' if good else 'FAIL'} {'total':7s} ops {got:>20,d} "
          f"vs census {want:>20,d}   {100 * err:6.2f}%")
    print(f"    tolerance 2% -- the PROFILER may carry a band; the census may not")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--i", type=int, required=True)
    ap.add_argument("--sample", type=int, default=400,
                    help="columns/buckets sampled for the MULTIPASS FACTOR only; "
                         "the ops count is exact and does not move with this")
    ap.add_argument("--k-lo", type=int, default=None,
                    help="Z-jump start; defaults to the member's K_EXACT+1")
    ap.add_argument("--check", action="store_true",
                    help="gate the result against work_census.py")
    ap.add_argument("--json_out", type=str, default=None)
    args = ap.parse_args()

    r = profile(args.i, args.sample, args.k_lo)
    report(r)
    ok = do_check(r) if args.check else True
    if args.json_out:
        p = Path(args.json_out)
        out = p if p.is_absolute() else ROOT / p
        out.write_text(json.dumps(r, indent=1))
        print(f"\n  wrote {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
