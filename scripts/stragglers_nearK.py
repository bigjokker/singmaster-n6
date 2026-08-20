#!/usr/bin/env python3
"""Kill the 25 Band I k with no prime in (k, K].

Those k are 4126622..4126646. Next prime after them is 4126651 > K,
which divides m. First possible kill is p > N/2.

Uses two-digit/three-digit Lucas, not giant m.
Does not touch nextprime_sweep.py.

p in (N/2, d]:  C(N,K) ≡ C(N-p, K) (mod p), can be nonzero.
p in (d, N]:    C(N,K) ≡ 0, cannot kill.
p > N:          single-digit C(N,K) mod p, can kill again.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2

from bandii_kernel import image_j, inv_table  # noqa: E402

from singmaster_intersect import binom_mod_lucas

N = 10_803_704
K = 4_126_647
D = N - K
KMIN, KMAX = 4_126_622, 4_126_646  # inclusive; not family columns
MAX_PRIMES = 40
OUT = ROOT / "results" / "stragglers_nearK.json"


def main() -> int:
    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    ks = list(range(KMIN, KMAX + 1))
    remaining = set(ks)
    rows = {
        k: {"k": k, "q": None, "g": None, "r": 0, "kill_j": None, "survived": []}
        for k in ks
    }
    t0 = time.time()
    p = int(gmpy2.next_prime(N // 2))
    nprimes = 0
    print(
        f"=== stragglers k={KMIN}..{KMAX}  first p>{N//2} = {p} ===",
        flush=True,
    )
    while remaining and nprimes < MAX_PRIMES:
        nprimes += 1
        zone = (
            "(N/2,d] can-kill"
            if p <= D
            else ("(d,N] forced-zero" if p <= N else "p>N can-kill")
        )
        t1 = time.time()
        mmod = binom_mod_lucas(N, K, p)
        print(
            f"  p={p}  zone={zone}  alpha={N//p} beta={K//p}  "
            f"mmod={mmod}  left={len(remaining)}  lucas={time.time()-t1:.2f}s",
            flush=True,
        )
        if mmod == 0:
            for k in remaining:
                rows[k]["r"] += 1
                rows[k]["survived"].append({"p": p, "g": p - k, "j": None, "zero": True})
            p = int(gmpy2.next_prime(p))
            continue
        gmax = max(p - k for k in remaining)
        inv = inv_table(p, gmax)
        newly = []
        for k in sorted(remaining):
            j = image_j(mmod, k, p, inv)
            if j is None:
                rows[k]["q"] = p
                rows[k]["g"] = p - k
                rows[k]["kill_j"] = None
                newly.append(k)
                print(f"    KILL k={k}  p={p}  g={p-k}  mmod={mmod}", flush=True)
            else:
                rows[k]["r"] += 1
                rows[k]["survived"].append(
                    {"p": p, "g": p - k, "j": j, "zero": False}
                )
                print(f"    live k={k}  p={p}  g={p-k}  j={j}", flush=True)
        remaining -= set(newly)
        p = int(gmpy2.next_prime(p))

    unkilled = sorted(remaining)
    payload = {
        "search": "stragglers_nearK",
        "N": N,
        "K": K,
        "d": D,
        "k_range": [KMIN, KMAX],
        "first_p": int(gmpy2.next_prime(N // 2)),
        "n_columns": len(ks),
        "n_killed": len(ks) - len(unkilled),
        "n_unkilled": len(unkilled),
        "unkilled": unkilled,
        "nprimes_used": nprimes,
        "seconds": round(time.time() - t0, 3),
        "rows": [rows[k] for k in ks],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        f"wrote {OUT}  killed={payload['n_killed']}  "
        f"unkilled={payload['n_unkilled']}  {payload['seconds']}s",
        flush=True,
    )
    return 0 if not unkilled else 3


if __name__ == "__main__":
    raise SystemExit(main())
