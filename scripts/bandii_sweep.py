#!/usr/bin/env python3
"""Band II p>N/2 image sweep for i=8. Spec: docs/bandii-spec.md.

Not a next-prime sweep from k. Live window (N/2, d]. Cap 14 primes.
Pre-flight §8 is mandatory and aborts the sweep on any failure.
Refuse if results/bandii_sweep.json exists. Resume from jsonl.

Does not touch Band I, Stage 4, E, or giant m.
"""

from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bandii_kernel import (  # noqa: E402
    CAP,
    D,
    K,
    KMAX,
    KMIN,
    LOG10_M,
    N,
    N2,
    NCOLS,
    P1,
    P2,
    PRIMES,
    R_P1,
    STRAGGLERS,
    equal_g_chunks,
    fact_table,
    log10_central,
    r_checked,
    r_closed,
    r_falling,
    r_from_F,
    scan_columns,
    scan_ks,
)


def load_done() -> list[dict]:
    from bandii_kernel import read_jsonl

    return read_jsonl(CHK)


def write_jsonl(rec: dict) -> None:
    from bandii_kernel import append_jsonl

    append_jsonl(CHK, rec)


def summarize(surv: list[dict]) -> dict:
    from bandii_kernel import summarize_survivors

    return summarize_survivors(surv)

OUT = ROOT / "results" / "bandii_sweep.json"
CHK = ROOT / "results" / "bandii_sweep.jsonl"
N_CHUNKS = 32
DEFAULT_WORKERS = 8
MAX_INLINE = 20_000

PREREGISTER = {
    1: {"alive": 1.026e5, "even": 0.658, "mean_k": 4_536_120, "band": [102_000, 103_200]},
    2: {"alive": 1.26e4, "even": 0.786, "mean_k": 4_450_918, "band": None},
    3: {"alive": 1816, "even": 0.875, "mean_k": 4_392_087, "band": None},
    4: {"alive": 289.9, "even": 0.930, "mean_k": 4_350_633, "band": [257, 324]},
    5: {"alive": 49.4, "even": 0.962, "mean_k": 4_320_257, "band": None},
    6: {"alive": 8.78, "even": 0.979, "mean_k": 4_297_132, "band": [3, 15]},
    7: {"alive": 1.61, "even": 0.989, "mean_k": 4_278_952, "band": None},
    8: {"alive": 0.300, "even": 0.994, "mean_k": 4_264_287, "band": None},
    10: {"alive": 0.0108, "even": 0.998, "mean_k": 4_242_079, "band": None},
    12: {"alive": 4.06e-4, "even": 1.000, "mean_k": 4_226_056, "band": None},
    14: {"alive": 1.56e-5, "even": 1.000, "mean_k": 4_213_949, "band": None},
}

CERTIFICATE = (
    "Every k in [4126649, 5182637] carries a modular impossibility "
    "certificate r(p) notin I_{p,k} for some prime p in {p1..pr}. "
    "Band II contains no extra representation of C(F_18 F_19, F_16 F_19). "
    "Unconditional."
)


def _job(payload: tuple) -> dict:
    p, r_expected, k_lo, k_hi, ks = payload
    t0 = time.time()
    r, rows = scan_columns(p, ks, r_expected=r_expected)
    return {
        "p": p,
        "r": r,
        "k_lo": k_lo,
        "k_hi": k_hi,
        "n_cols": len(ks),
        "n_survivors": len(rows),
        "survivors": rows,
        "seconds": round(time.time() - t0, 3),
    }


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: `python -O` strips assert, and every
    certificate downstream of this function is claimed unconditional."""
    if not cond:
        raise RuntimeError(f"pre-flight failed: {msg}")


def preflight() -> None:
    print("=== Band II pre-flight §8 ===", flush=True)
    t0 = time.time()

    from singmaster_intersect import binom_mod_lucas, fib

    check(N == 10_803_704 and K == 4_126_647, "N,K constants")
    check(int(fib(16) * fib(19)) == K, "K = F_16 F_19")
    check(int(fib(18) * fib(17)) == K + 1, "K+1 = F_18 F_17")
    check(K != int(fib(18) * fib(17)), "K is not F_18 F_17")
    check(D == 6_677_057 and KMAX == 5_182_637 and N2 == 5_401_852, "D,KMAX,N2")
    check(log10_central(KMAX) <= LOG10_M < log10_central(KMAX + 1), "KMAX brackets log10 m")
    check(KMAX - (K + 2) + 1 == 1_055_989 == NCOLS, "NCOLS")
    print("  8.1 constants ok", flush=True)

    import gmpy2

    for p in PRIMES:
        check(bool(gmpy2.is_prime(p)), f"{p} is prime")
        check(2 * p > N and p <= D and p > KMAX, f"{p} in Band II live window")
    print("  8.1 primes ok", flush=True)

    import numpy as np

    rng = np.random.default_rng(0)
    bad = 0
    ncase = 0
    for p in (11, 29, 101, 211, 1009):
        F = fact_table(p)
        check(int(F[p - 1]) == p - 1, f"Wilson at p={p}")
        for _ in range(200):
            k = int(rng.integers(1, p))
            rr = int(rng.integers(0, p))
            s = rr * int(F[k]) % p
            ker = bool(np.any(F[k:p] == (np.int64(s) * F[: p - k]) % np.int64(p)))
            brute = any(math.comb(n0, k) % p == rr for n0 in range(k, p))
            ncase += 1
            if ker != brute:
                bad += 1
    if bad:
        raise RuntimeError(f"8.2 kernel mismatches {bad}/{ncase}")
    print(f"  8.2 kernel {ncase} cases, 0 mismatches", flush=True)

    F = fact_table(P1)
    check(int(F[P1 - 1]) == P1 - 1, f"Wilson at p1={P1}")
    ra = r_from_F(F, P1)
    rb = r_falling(P1)
    rc = r_closed(P1)
    rl = int(binom_mod_lucas(N, K, P1))
    if not (ra == rb == rc == rl == R_P1):
        raise RuntimeError(f"8.4 r(p1) {ra} {rb} {rc} {rl} != {R_P1}")
    print(f"  8.4 r(p1)={ra} four ways ok", flush=True)

    hits = scan_ks(F, P1, ra, range(4_126_622, 4_126_647))
    got = {h["k"]: h["b"] for h in hits}
    if got != STRAGGLERS:
        raise RuntimeError(f"8.3 stragglers {got} != {STRAGGLERS}")
    r2 = r_closed(P2)
    F2 = fact_table(P2)
    r2t = r_checked(F2, P2)
    if r2t != r2:
        raise RuntimeError("8.3 r(p2) mismatch")
    still = scan_ks(F2, P2, r2t, STRAGGLERS.keys())
    if still:
        raise RuntimeError(f"8.3 still alive at p2: {still}")
    print("  8.3 stragglers exact four k and four b; p2 kills all", flush=True)

    spot = scan_ks(F, P1, ra, range(4_126_649, 4_126_689))
    if len(spot) != 9:
        raise RuntimeError(f"8.5 spot survivors {len(spot)} != 9")
    print("  8.5 hardest 40 columns: 9 survivors", flush=True)

    chunks = equal_g_chunks(KMIN, KMAX, P1, N_CHUNKS)
    covered = sum(hi - lo + 1 for lo, hi in chunks)
    if covered != NCOLS:
        raise RuntimeError(f"chunks cover {covered} != {NCOLS}")
    print(
        f"  chunks {N_CHUNKS}  first {chunks[0]}  last {chunks[-1]}  "
        f"{time.time()-t0:.1f}s",
        flush=True,
    )
    print("=== pre-flight passed ===", flush=True)


def chunk_key(prime_index: int, k_lo: int, k_hi: int) -> tuple[int, int, int]:
    return (prime_index, k_lo, k_hi)


def pass_n(prime_index: int, p: int, chunks: list[tuple[int, int]], prev: list[dict] | None, workers: int, done: set[tuple[int, int, int]]) -> list[dict]:
    r_expected = r_closed(p)
    jobs = []
    if prev is None:
        by_chunk_ks = None
    else:
        by_chunk_ks = { (lo, hi): [] for lo, hi in chunks }
        for s in prev:
            k = s["k"]
            for lo, hi in chunks:
                if lo <= k <= hi:
                    by_chunk_ks[(lo, hi)].append(k)
                    break
    pending = []
    survivors: list[dict] = []
    for lo, hi in chunks:
        key = chunk_key(prime_index, lo, hi)
        if key in done:
            continue
        if prev is None:
            ks = list(range(lo, hi + 1))
        else:
            ks = by_chunk_ks[(lo, hi)]
        pending.append((p, r_expected, lo, hi, ks))

    already = [rec for rec in load_done() if rec.get("prime_index") == prime_index]
    for rec in already:
        for s in rec.get("survivors") or []:
            if "g" not in s:
                s["g"] = p - s["k"]
                s["g_even"] = s["g"] % 2 == 0
            survivors.append(s)

    if not pending:
        print(f"  pass {prime_index} p={p} all chunks already done  alive={len(survivors)}", flush=True)
        return survivors

    print(
        f"  pass {prime_index} p={p} r={r_expected} d={2*p-N}  "
        f"jobs={len(pending)} workers={workers} prev={0 if prev is None else len(prev)}",
        flush=True,
    )
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        for rec in pool.imap_unordered(_job, pending):
            rec["prime_index"] = prime_index
            rec["delta"] = 2 * p - N
            rec["r"] = r_expected
            rows = rec["survivors"]
            if len(rows) > MAX_INLINE:
                rec["survivors"] = [{"k": x["k"], "b": x["b"]} for x in rows]
            write_jsonl(rec)
            survivors.extend(rows)
            print(
                f"    chunk [{rec['k_lo']},{rec['k_hi']}] cols={rec['n_cols']} "
                f"surv={rec['n_survivors']} {rec['seconds']}s  "
                f"alive_so_far={len(survivors)}",
                flush=True,
            )
    return survivors


def main() -> int:
    if "--preflight" in sys.argv:
        try:
            import numpy as np  # noqa: F401
        except ImportError:
            print("numpy is required. Not installed.", flush=True)
            return 1
        preflight()
        return 0

    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        print("numpy is required. Install numpy >= 1.24 before this bat.", flush=True)
        return 1

    workers = int(os.environ.get("BANDII_WORKERS", DEFAULT_WORKERS))
    workers = max(1, min(workers, N_CHUNKS))
    OUT.parent.mkdir(exist_ok=True)

    t0 = time.time()
    preflight()
    chunks = equal_g_chunks(KMIN, KMAX, P1, N_CHUNKS)
    done_recs = load_done()
    done_keys = {
        chunk_key(r["prime_index"], r["k_lo"], r["k_hi"])
        for r in done_recs
        if "prime_index" in r and "k_lo" in r
    }
    print(
        f"=== Band II sweep  cols={NCOLS}  chunks={N_CHUNKS}  "
        f"workers={workers}  cap={CAP}  resume={len(done_keys)} ===",
        flush=True,
    )
    print("preregister pass1 ~1.026e5  pass4 ~290  pass6 ~8.8  pass8 ~0.30", flush=True)

    prev: list[dict] | None = None
    passes = []
    last_complete = 0
    by_prime: dict[int, list[dict]] = {}
    for rec in done_recs:
        if "prime_index" in rec and "survivors" in rec:
            by_prime.setdefault(rec["prime_index"], []).append(rec)
    for pi in sorted(by_prime):
        ch_done = {(r["k_lo"], r["k_hi"]) for r in by_prime[pi]}
        if ch_done == set(chunks):
            last_complete = pi
            prev = []
            for r in by_prime[pi]:
                prev.extend(r["survivors"])

    start = last_complete + 1
    alive: list[dict] = prev or []
    for i, p in enumerate(PRIMES[:CAP], start=1):
        if i < start:
            continue
        if i > 1 and not alive and last_complete >= i - 1:
            print(f"  all dead before pass {i}. stopping.", flush=True)
            break
        src = None if i == 1 else alive
        alive = pass_n(i, p, chunks, src, workers, done_keys)
        summ = summarize(alive)
        pred = PREREGISTER.get(i)
        rec = {
            "prime_index": i,
            "p": p,
            "n_alive": summ["n"],
            "even_g_frac": summ["even"],
            "mean_k": summ["mean_k"],
            "pred": pred,
        }
        passes.append(rec)
        print(
            f"  PASS {i} alive={summ['n']} even={summ['even']} mean_k={summ['mean_k']}"
            + (
                f"  pred {pred['alive']} even {pred['even']} mean_k {pred['mean_k']}"
                if pred
                else ""
            ),
            flush=True,
        )
        done_keys = {
            chunk_key(r["prime_index"], r["k_lo"], r["k_hi"])
            for r in load_done()
            if "prime_index" in r and "k_lo" in r
        }
        if not alive:
            break

    clean = len(alive) == 0
    payload = {
        "search": "bandii_sweep",
        "N": N,
        "K": K,
        "d": D,
        "k_min": KMIN,
        "k_max": KMAX,
        "n_columns": NCOLS,
        "primes": PRIMES[:CAP],
        "cap": CAP,
        "n_chunks": N_CHUNKS,
        "workers": workers,
        "preregister": PREREGISTER,
        "passes": passes,
        "n_final_alive": len(alive),
        "final_survivors": alive,
        "clean": clean,
        "certificate": CERTIFICATE if clean else None,
        "seconds": round(time.time() - t0, 3),
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(flush=True)
    print(
        f"wrote {OUT}  passes={len(passes)}  final_alive={len(alive)}  "
        f"clean={clean}  {payload['seconds']}s",
        flush=True,
    )
    if clean:
        print("Band II closed. Unconditional certificates for all 1055989 columns.", flush=True)
    else:
        print("ALIVE AT CAP. Log and hand back. Do not extend. Do not touch Band I.", flush=True)
        for s in alive:
            print(f"  k={s['k']} g={s['g']} even={s['g_even']} b={s['b']}", flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
