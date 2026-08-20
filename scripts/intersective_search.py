#!/usr/bin/env python3
"""Q14: can (x)_k - c be intersective?

A polynomial is INTERSECTIVE if it has a root modulo every prime but no
rational root. For this project the question is not decorative:

    column k survives prime p   <=>  (x)_k - k!m  has a root mod p
    column k represents m       <=>  (x)_k - k!m  has an integer root

so an intersective (x)_k - c is exactly a column that survives every prime
without representing m -- a census that never terminates. It is the only
logical gap between the method and a guarantee that it always finishes.

WHAT IS PROVED ELSEWHERE (see docs/q14-intersective.md):
  k <= 4  never intersective, elementary.
  k = 5   never intersective. A 2+3 split forces 5a^4-10a^2+9 to be square,
          which is the Pell equation v^2-5u^2=4 with u=a^2-1, so a^2 = F_2j+1;
          F_n+1 is a square only for n=0,4,6, giving c in {0,+-210,+-2160},
          and each of those is killed by an explicit small prime.

WHAT THIS SCRIPT DOES: searches k >= 6, where no such clean reduction is
available, using two necessary conditions in increasing cost --

  (1) rad(k!) | c.  For p <= k the product of k consecutive residues always
      contains a multiple of p, so (x)_k vanishes identically as a FUNCTION
      on F_p and f_c = -c there. This alone cuts the search by 30x (k=5,6)
      to 2310x (k=12,13).

  (2) c mod p must lie in I_p = {(x)_k mod p} for every p > k. Just above k
      the image is tiny -- |I_p| = 2 at p = k+1 when that is prime -- so the
      first few primes above k are brutally selective.

Survivors are then factored over Q. A survivor is interesting only if it has
no rational root and every irreducible factor has degree >= 2; anything else
is either a genuine falling-factorial value or ruled out by Jordan.

    python scripts/intersective_search.py --kmax 14 --bound 1000000000
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def radical_factorial(k: int) -> int:
    from sympy import isprime

    r = 1
    for q in range(2, k + 1):
        if isprime(q):
            r *= q
    return r


def falling_image(k: int, p: int) -> set[int]:
    """{(x)_k mod p : x}. For p <= k this is {0}."""
    img = set()
    for v in range(p):
        prod = 1
        for i in range(k):
            prod = prod * (v - i) % p
        img.add(prod)
    return img


def sieve_primes_for(k: int, n: int) -> list[int]:
    from sympy import isprime

    out, p = [], k + 1
    while len(out) < n:
        if isprime(p):
            out.append(p)
        p += 1
    return out


def search_k(k: int, bound: int, n_primes: int, chunk: int) -> dict:
    import numpy as np

    rad = radical_factorial(k)
    primes = sieve_primes_for(k, n_primes)
    allowed = {}
    for p in primes:
        img = falling_image(k, p)
        a = np.zeros(p, dtype=bool)
        for v in img:
            a[v] = True
        allowed[p] = a

    jmax = bound // rad
    survivors: list[int] = []
    t0 = perf_counter()
    for lo in range(-jmax, jmax + 1, chunk):
        hi = min(lo + chunk, jmax + 1)
        j = np.arange(lo, hi, dtype=np.int64)
        c = j * rad
        mask = np.ones(c.size, dtype=bool)
        for p in primes:
            mask &= allowed[p][c % p]
            if not mask.any():
                break
        survivors.extend(int(v) for v in c[mask])
    return {
        "k": k,
        "rad_k_factorial": rad,
        "bound": bound,
        "n_candidates": 2 * jmax + 1,
        "sieve_primes": primes,
        "n_survivors": len(survivors),
        "survivors": survivors,
        "sieve_seconds": round(perf_counter() - t0, 2),
    }


def classify(k: int, c: int) -> dict:
    """Factor (x)_k - c over Q and say what kind of object it is."""
    from sympy import Poly, factor_list, symbols

    x = symbols("x")
    f = Poly(1, x)
    for i in range(k):
        f = f * Poly(x - i, x)
    f = f - Poly(c, x)
    facs = [(Poly(g, x).degree(), Poly(g, x).as_expr()) for g, _ in
            factor_list(f.as_expr(), x)[1]]
    degs = sorted(d for d, _ in facs)
    has_rational_root = 1 in degs
    return {
        "c": c,
        "degrees": degs,
        "rational_root": has_rational_root,
        "irreducible": degs == [k],
        # only these could possibly be intersective
        "candidate": (not has_rational_root) and degs != [k] and min(degs) >= 2,
        "factors": [str(e) for _, e in facs] if degs != [k] else None,
    }


def first_killing_prime(k: int, c: int, limit: int = 100000) -> int | None:
    """Smallest p with no root of (x)_k - c mod p. None means none found."""
    from sympy import isprime

    p = 2
    while p < limit:
        if isprime(p):
            hit = False
            for v in range(p):
                prod = 1
                for i in range(k):
                    prod = prod * (v - i) % p
                if prod == c % p:
                    hit = True
                    break
            if not hit:
                return p
        p += 1
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kmin", type=int, default=5)
    ap.add_argument("--kmax", type=int, default=14)
    ap.add_argument("--bound", type=int, default=10**8)
    ap.add_argument("--primes", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    report = {"search": "intersective", "bound": args.bound,
              "n_sieve_primes": args.primes, "k": {}}
    print(f"  |c| <= {args.bound:,}, {args.primes} sieve primes above k\n")
    print(f"  {'k':>3} {'rad(k!)':>8} {'candidates':>14} {'survive':>8} "
          f"{'reducible':>10} {'CANDIDATE':>10} {'sieve s':>8}")
    for k in range(args.kmin, args.kmax + 1):
        res = search_k(k, args.bound, args.primes, args.chunk)
        cls = [classify(k, c) for c in res["survivors"]]
        red = [d for d in cls if not d["irreducible"]]
        cand = [d for d in cls if d["candidate"]]
        res["classified"] = cls
        res["n_reducible"] = len(red)
        res["candidates"] = cand
        for d in cand:
            d["first_killing_prime"] = first_killing_prime(k, d["c"])
        report["k"][k] = res
        print(f"  {k:>3} {res['rad_k_factorial']:>8} {res['n_candidates']:>14,} "
              f"{res['n_survivors']:>8} {len(red):>10} {len(cand):>10} "
              f"{res['sieve_seconds']:>8.1f}")
        for d in cand:
            kp = d["first_killing_prime"]
            print(f"        c={d['c']}  degrees {d['degrees']}  "
                  + (f"killed by p={kp}" if kp else "NO KILLING PRIME FOUND"))

    n_cand = sum(len(r["candidates"]) for r in report["k"].values())
    n_open = sum(1 for r in report["k"].values() for d in r["candidates"]
                 if d.get("first_killing_prime") is None)
    print(f"\n  factorisation candidates (no rational root, all factors >= 2): {n_cand}")
    print(f"  of those, still unkilled: {n_open}")
    print("  RESULT", "no intersective polynomial found" if n_open == 0
          else "OPEN CANDIDATE -- investigate")
    report["n_candidates"] = n_cand
    report["n_unkilled"] = n_open
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
