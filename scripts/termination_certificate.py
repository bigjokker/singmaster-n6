#!/usr/bin/env python3
"""Q14: prove that the census MUST terminate for a given column.

The worry Q14 names: if (x)_k - k!m were intersective -- a root modulo every
prime but no integer root -- column k would survive every prime without
representing m, and the census would never finish. That is the only logical
gap between the method and a guarantee.

The general question is open for k >= 6 (see docs/q14-intersective.md). But it
does not have to be settled to close the gap FOR A GIVEN COLUMN, because:

    (x)_k - k!m irreducible over Q
      => Gal acts transitively on the k roots
      => Jordan: a transitive group of degree >= 2 has a derangement
      => Chebotarev: the primes whose Frobenius is that derangement have
         positive density, and for each of them the polynomial has NO root
      => a killing prime exists, with density >= 1/k (Cameron-Cohen)

and irreducibility over Q has a cheap certificate: if the polynomial is
irreducible modulo a single prime q, it is irreducible over Q.

So the tool is: find q with (x)_k - k!m irreducible mod q. Since m mod q is
just Lucas, this never builds m -- which matters, since m has 3.1 million
digits at i=8 and 21 million at i=9.

The density of such q is roughly 1/k (the proportion of k-cycles), so a few k
primes suffice. Cost is O(k^2 log q) per prime, so this is practical to about
k = few hundred -- which is exactly the regime where columns are hardest to
kill, since small columns sit at g/p up to 0.94 and survive many live primes.
Large columns die almost at once and never need a certificate.

    python scripts/termination_certificate.py --i 8 --kmax 60
    python scripts/termination_certificate.py --i 9 --k 11
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2  # noqa: E402

from singmaster_intersect import (  # noqa: E402
    _pgcd,
    _pmul,
    binom_mod_lucas,
    falling_poly,
    fib_member,
)


def _pow_q(g, f, q):
    """g^q mod f over F_q."""
    e, base, acc = q, list(g), [1]
    while e:
        if e & 1:
            acc = _pmul(acc, base, f, q)
        base = _pmul(base, base, f, q)
        e >>= 1
    return acc


def irreducible_mod_q(k: int, r: int, q: int) -> bool:
    """Is (x)_k - k! r irreducible over F_q?

    Distinct-degree test: a reducible degree-k polynomial has an irreducible
    factor of degree <= k/2, and gcd(x^(q^d) - x, f) collects exactly the
    factors whose degree divides d. So if that gcd is trivial for every
    d <= k/2, f is irreducible. Requires q > k, else k! vanishes and f
    degenerates to (x)_k, which splits completely.
    """
    if q <= k:
        raise ValueError(f"need q > k (got q={q}, k={k})")
    f = falling_poly(k, q)
    kf = 1
    for i in range(1, k + 1):
        kf = kf * i % q
    f[0] = (f[0] - kf * r) % q
    cur = _pow_q([0, 1], f, q)
    for d in range(1, k // 2 + 1):
        h = list(cur)
        while len(h) < 2:
            h.append(0)
        h[1] = (h[1] - 1) % q
        while len(h) > 1 and h[-1] == 0:
            h.pop()
        if h == [0]:
            # f divides x^(q^d) - x with d < k: every factor has degree <= d
            return False
        if len(_pgcd(f, h, q)) > 1:
            return False
        cur = _pow_q(cur, f, q)
    return True


def certify(N: int, K: int, k: int, tries: int = 600) -> dict:
    """Find q proving (x)_k - k!m irreducible over Q, hence termination."""
    t0 = perf_counter()
    q = int(gmpy2.next_prime(max(k, 2)))
    live = dead = 0
    for _ in range(tries):
        r = int(binom_mod_lucas(N, K, q))
        if r:
            live += 1
            if irreducible_mod_q(k, r, q):
                return {
                    "k": k,
                    "certified": True,
                    "q": q,
                    "r": r,
                    "live_primes_tried": live,
                    "dead_primes_skipped": dead,
                    "seconds": round(perf_counter() - t0, 3),
                    "conclusion": (
                        f"(x)_{k} - {k}!m is irreducible mod {q}, hence over Q; "
                        f"by Jordan its Galois group has a derangement, so a "
                        f"prime with no root exists (density >= 1/{k})."
                    ),
                }
        else:
            dead += 1  # r=0: f = (x)_k splits completely, tells us nothing
        q = int(gmpy2.next_prime(q))
    return {
        "k": k,
        "certified": False,
        "live_primes_tried": live,
        "dead_primes_skipped": dead,
        "seconds": round(perf_counter() - t0, 3),
        "conclusion": (
            "no irreducible reduction found. That is NOT evidence of "
            "intersectivity -- it may simply be reducible, in which case the "
            "factor degrees are the thing to look at."
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--i", type=int, required=True)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--kmax", type=int, default=40)
    ap.add_argument("--tries", type=int, default=600)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    mem = fib_member(args.i, compute_m=False)
    N, K = mem.n, mem.k
    ks = [args.k] if args.k else [
        k for k in range(2, args.kmax + 1) if k not in (K, K + 1)
    ]
    print(f"  i={args.i}  N={N:,}  K={K:,}   (m never built)")
    print(f"  {'k':>5} {'certifying q':>14} {'live q':>8} {'dead q':>8} {'time':>8}")
    rows = []
    for k in ks:
        res = certify(N, K, k, args.tries)
        rows.append(res)
        q = f"{res['q']:,}" if res["certified"] else "NONE"
        print(f"  {k:>5} {q:>14} {res['live_primes_tried']:>8} "
              f"{res['dead_primes_skipped']:>8} {res['seconds']:>7.2f}s")
    n_ok = sum(1 for r in rows if r["certified"])
    print(f"\n  certified {n_ok}/{len(rows)} columns: for each, a killing prime")
    print(f"  provably exists, so the census had to terminate -- independently")
    print(f"  of whether it did.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(
            {"search": "termination_certificate", "i": args.i, "N": N, "K": K,
             "n_certified": n_ok, "n_columns": len(rows), "rows": rows},
            indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0 if n_ok == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
