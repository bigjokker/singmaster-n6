#!/usr/bin/env python3
"""Stage-1 next-prime obstruction sweep.

For a Fibonacci member i and each extra column k in [kmin, kmax]
(skipping the two family columns), test the next `nprimes` primes p > k.
Record the least killing prime q(k), the gap q(k)-k, and how many
consecutive primes after k survive (r(k)).

This is a finite census, not a theorem. Family columns must stay possible.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import gmpy2

from singmaster_intersect import (
    binom_mod_lucas,
    column_possible,
    fib_member,
)


# How often to drop dead r(p) entries. The walk window is a few hundred
# primes wide, so pruning this rarely costs nothing and bounds the table.
PRUNE_EVERY = 1000


def next_primes_after(n: int, count: int) -> list[int]:
    p = int(n)
    out: list[int] = []
    while len(out) < count:
        p = int(gmpy2.next_prime(p))
        out.append(p)
    return out


class RCache:
    """Memoised r(p) = C(N,K) mod p for one fixed (N,K).

    r(p) does not depend on the column k, but the walk is per-column, so
    without a cache every (k,p) pair recomputes the same Lucas product.
    Consecutive columns share their next primes: Stage 3 (i=8,
    k=10^5..10^6) ran ~900k column-walks against ~34.7k live primes, i.e.
    it evaluated each r(p) about 26 times.

    The sweep visits columns in increasing k and only ever asks about
    p > k, so an entry with p <= k can never be asked for again; `prune`
    drops those and keeps the table at the width of the walk window
    rather than the width of the whole k-range.
    """

    def __init__(self, N: int, K: int) -> None:
        self.N = N
        self.K = K
        self._r: dict[int, int] = {}
        self.lookups = 0
        self.misses = 0

    def r(self, p: int) -> int:
        self.lookups += 1
        val = self._r.get(p)
        if val is None:  # r(p)=0 is a legitimate cached value, not a miss
            val = binom_mod_lucas(self.N, self.K, p)
            self._r[p] = val
            self.misses += 1
        return val

    def prune(self, k: int) -> None:
        """Forget primes no later column can ask about (p <= k)."""
        if self._r:
            self._r = {q: v for q, v in self._r.items() if q > k}

    def stats(self) -> dict:
        return {
            "lookups": self.lookups,
            "evaluations": self.misses,
            "reuse": round(self.lookups / self.misses, 3) if self.misses else None,
            "resident": len(self._r),
        }


def walk_until_kill(
    N: int,
    K: int,
    k: int,
    *,
    nprimes: int,
    max_gap: int,
    rcache: RCache | None = None,
) -> dict:
    """Least p>k that obstructs, or None if both budgets expire.

    nprimes<=0 means no prime-count cap (still stop if p-k exceeds max_gap).
    `rcache` only memoises r(p); with or without it the rows are identical.
    """
    p = k
    r = 0
    while True:
        p = int(gmpy2.next_prime(p))
        if max_gap > 0 and p - k > max_gap:
            return {
                "k": k,
                "q": None,
                "gap": None,
                "r": r,
                "survived_all": True,
                "stop": "max_gap",
            }
        if nprimes > 0 and r >= nprimes:
            return {
                "k": k,
                "q": None,
                "gap": None,
                "r": r,
                "survived_all": True,
                "stop": "nprimes",
            }
        m_mod = rcache.r(p) if rcache is not None else binom_mod_lucas(N, K, p)
        if not column_possible(m_mod, k, p):
            return {
                "k": k,
                "q": p,
                "gap": p - k,
                "r": r,
                "survived_all": False,
                "stop": "killed",
            }
        r += 1


def _load_checkpoint(path: Path) -> dict[int, dict]:
    done: dict[int, dict] = {}
    if not path.exists():
        return done
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            done[int(rec["k"])] = rec
    return done


def sweep(
    i: int,
    kmin: int,
    kmax: int,
    nprimes: int,
    checkpoint: Path | None = None,
    max_gap: int = 0,
) -> dict:
    mem = fib_member(i, compute_m=False)
    N, K = mem.n, mem.k
    family = {K, K + 1}
    done = _load_checkpoint(checkpoint) if checkpoint is not None else {}
    if done:
        print(f"  resume: {len(done)} columns already in {checkpoint}", flush=True)
    chk = checkpoint.open("a", encoding="utf-8") if checkpoint is not None else None
    rcache = RCache(N, K)
    rows = []
    unkilled = []
    t0 = time.time()
    total = kmax - kmin + 1
    try:
        for k in range(kmin, kmax + 1):
            if k in family:
                continue
            if k in done:
                rec = done[k]
            else:
                rec = walk_until_kill(
                    N, K, k, nprimes=nprimes, max_gap=max_gap, rcache=rcache
                )
                if chk is not None:
                    chk.write(json.dumps(rec) + "\n")
                    if (k - kmin + 1) % 200 == 0:
                        chk.flush()
            rows.append(rec)
            if rec.get("survived_all"):
                unkilled.append(k)
            done_n = k - kmin + 1
            if done_n % PRUNE_EVERY == 0:
                rcache.prune(k)
            if done_n % 5000 == 0 or k == kmax:
                dt = time.time() - t0
                frac = done_n / total
                eta = (dt / frac - dt) if frac > 0 else 0
                print(
                    f"  k={k}  {100*frac:.1f}%  elapsed={dt:.0f}s  eta={eta:.0f}s",
                    flush=True,
                )
    finally:
        if chk is not None:
            chk.flush()
            chk.close()

    killed = [row for row in rows if row["q"] is not None]
    gaps = [row["gap"] for row in killed]
    rs = [row["r"] for row in rows]
    return {
        "search": "nextprime_sweep",
        "i": i,
        "N": N,
        "K": K,
        "k_range": [kmin, kmax],
        "nprimes": nprimes,
        "max_gap_cap": max_gap,
        "n_columns": len(rows),
        "n_unkilled": len(unkilled),
        "unkilled": unkilled,
        "max_gap": max(gaps) if gaps else None,
        "max_r": max(rs) if rs else None,
        "gap_hist": dict(sorted(Counter(gaps).items())),
        "r_hist": dict(sorted(Counter(rs).items())),
        "seconds": round(time.time() - t0, 3),
        "r_cache": rcache.stats(),
        "rows": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--i", type=int, default=8)
    ap.add_argument("--kmin", type=int, default=401)
    ap.add_argument("--kmax", type=int, default=10000)
    ap.add_argument(
        "--nprimes",
        type=int,
        default=20,
        help="Max consecutive surviving primes. 0 = no count cap (use --max-gap).",
    )
    ap.add_argument(
        "--max-gap",
        type=int,
        default=0,
        help="Stop if p-k exceeds this. 0 = no gap cap. Need at least one cap.",
    )
    ap.add_argument("--json_out", type=Path, default=None)
    ap.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help="JSONL of finished rows; resume if the file already exists",
    )
    args = ap.parse_args()
    if args.nprimes <= 0 and args.max_gap <= 0:
        raise SystemExit("need --nprimes > 0 or --max-gap > 0 so a never-kill cannot hang")
    print(
        f"=== NEXTPRIME i={args.i} k={args.kmin}..{args.kmax} "
        f"nprimes={args.nprimes} max_gap={args.max_gap} ===",
        flush=True,
    )
    chk = args.checkpoint
    if chk is None and args.json_out is not None:
        chk = args.json_out.with_suffix(".jsonl")
    rep = sweep(
        args.i,
        args.kmin,
        args.kmax,
        args.nprimes,
        checkpoint=chk,
        max_gap=args.max_gap,
    )
    print(
        f"columns={rep['n_columns']}  unkilled={rep['n_unkilled']}  "
        f"max_gap={rep['max_gap']}  max_r={rep['max_r']}  "
        f"{rep['seconds']}s",
        flush=True,
    )
    rc = rep["r_cache"]
    print(
        f"r(p): {rc['lookups']} lookups, {rc['evaluations']} evaluations, "
        f"reuse={rc['reuse']}x",
        flush=True,
    )
    print(f"r_hist={rep['r_hist']}", flush=True)
    print(f"gap_hist (first 20)={dict(list(rep['gap_hist'].items())[:20])}", flush=True)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        slim = dict(rep)
        # keep rows; 10k records is small
        args.json_out.write_text(json.dumps(slim), encoding="utf-8")
        print(f"wrote {args.json_out}", flush=True)
    return 0 if rep["n_unkilled"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
