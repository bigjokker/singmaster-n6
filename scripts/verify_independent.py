#!/usr/bin/env python3
"""Q10: re-check certificates through a different implementation.

witness.py shares no code path with the sweep, but both it and the sweep rest
on the same author's arithmetic. This checks a sample of certificates using
sympy instead -- a separate codebase, and more importantly a DIFFERENT PIECE
OF MATHEMATICS for the central claim.

Four routes, chosen by cost, in decreasing independence:

  brute    math.comb(x, k) % p for every x in [0, p). No cleverness at all.
           Exact, but O(p) big-integer binomials, so only for small p.

  poly     sympy factors f(x) = (x)_k - k! r over GF(p) and asks for linear
           factors. Cantor-Zassenhaus -- no relation to walking an image.
           "Column k represents r mod p" is exactly "f has a root in F_p".

  factor   F[k+b] == r F[k] F[b] with F[x] = x! mod p, written fresh here.
           A different identity from the upper-index negation witness.py
           walks, so it cross-checks witness.py -- but it is the identity the
           sweep uses, so it does not independently certify the sweep.

  lucas    C(N,K) mod p rebuilt with sympy.binomial on the base-p digits,
           against witness.lucas_mod_pure. Applied to every certificate.

What this does and does not buy: sympy is a different library with different
authors and algorithms, so it catches implementation error in witness.py. It
is not a different LANGUAGE or a different CPU, so it does not catch a shared
misunderstanding of the mathematics -- for that the `brute` route is the real
control, since it assumes nothing beyond the definition of C(n,k).

    python scripts/verify_independent.py --file results/i7_witness.npz --sample 40
    python scripts/verify_independent.py --file results/i3_witness.npz --all --brute
"""

from __future__ import annotations

import argparse
import sys
from math import comb
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import witness as W  # noqa: E402


def sympy_is_prime(p: int) -> bool:
    from sympy import isprime

    return bool(isprime(p))


def sympy_lucas(N: int, K: int, p: int) -> int:
    """C(N,K) mod p via Lucas, with sympy's binomial for each digit."""
    from sympy import binomial

    r, n, k = 1, int(N), int(K)
    while n or k:
        ni, ki = n % p, k % p
        if ki > ni:
            return 0
        r = r * int(binomial(ni, ki)) % p
        n //= p
        k //= p
    return r


def sympy_has_root(k: int, r: int, p: int) -> bool:
    """Does (x)_k - k! r have a root in F_p? sympy factors over GF(p).

    Nothing here walks an image or touches a factorial table: the question
    "is r in column k's image" is being answered as a polynomial
    factorisation, which is a different theorem and a different algorithm.
    """
    from sympy import GF, Poly, factor_list, symbols

    x = symbols("x")
    f = Poly(1, x, domain=GF(p))
    for i in range(k):
        f = f * Poly(x - i, x, domain=GF(p))
    kf = 1
    for i in range(1, k + 1):
        kf = kf * i % p
    f = f - Poly(kf * r % p, x, domain=GF(p))
    _, factors = factor_list(f.as_expr(), x, modulus=p)
    return any(Poly(fac, x, domain=GF(p)).degree() == 1 for fac, _ in factors)


def brute_has_root(k: int, r: int, p: int) -> bool:
    """math.comb over the whole domain. Assumes only the definition."""
    return any(comb(x, k) % p == r % p for x in range(p))


def factorial_has_root(k: int, r: int, p: int) -> bool:
    """C(k+b,k) = r for some b, via factorials. Written fresh, uses no repo code.

    This is a DIFFERENT identity from the one witness.py walks. witness.py
    uses C(k+j,k) = (-1)^j C(p-k-1,j), i.e. upper-index negation, and carries
    two running products with no inverses. Here the test is the direct
    definition C(k+b,k) = (k+b)! / (k! b!), i.e.

        F[k+b] == r * F[k] * F[b]   (mod p),  F[x] = x! mod p,

    which is how the sweep kernel does it. So agreement between this and
    witness.py cross-checks witness.py against a second derivation; it does
    not certify the sweep, which uses this same identity (vectorised).

    O(p) to build the table, O(g) to test -- affordable where brute force is
    not, and the only route available at Band II sizes.
    """
    F = [1] * p
    acc = 1
    for x in range(1, p):
        acc = acc * x % p
        F[x] = acc
    s = r % p * F[k] % p
    return any(F[k + b] == s * F[b] % p for b in range(p - k))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--file", type=Path, required=True)
    ap.add_argument("--sample", type=int, default=25)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--brute-max-p", type=int, default=40000,
                    help="use the math.comb control when p is at most this")
    ap.add_argument("--poly-max-k", type=int, default=60,
                    help="use sympy GF(p) factorisation when k is at most this")
    args = ap.parse_args()

    import numpy as np

    ks, ps, meta = W.load(args.file)
    N, K = int(meta["N"]), int(meta["K"])
    idx = np.arange(ks.size)
    if not args.all and args.sample < ks.size:
        idx = np.random.default_rng(0).choice(ks.size, args.sample, replace=False)
        idx.sort()

    print(f"  {args.file.name}: i={meta.get('i')} N={N:,} K={K:,}")
    print(f"  checking {len(idx)} of {ks.size:,} certificates through sympy",
          flush=True)

    bad, n_brute, used = [], 0, {}
    t0 = perf_counter()
    for j, i in enumerate(idx):
        k, p = int(ks[i]), int(ps[i])
        ours = W.check_witness(N, K, k, p)

        if not sympy_is_prime(p):
            bad.append((k, p, "sympy says p is not prime"))
            continue
        r_sym = sympy_lucas(N, K, p)
        if r_sym != ours.get("r"):
            bad.append((k, p, f"r disagrees: ours {ours.get('r')} sympy {r_sym}"))
            continue
        if r_sym == 0:
            bad.append((k, p, "r=0 certifies nothing"))
            continue
        # Pick the membership route by cost. Brute force is the strongest
        # control (it assumes only the definition of C(n,k)) and is cheap when
        # p is small. sympy's GF(p) factorisation is cheap only for small k --
        # Cantor-Zassenhaus on a degree-300 polynomial is far slower than the
        # O(p) walk it is meant to check, which is the same crossover Q7 found.
        route = None
        if p <= args.brute_max_p:
            route = "brute"
            n_brute += 1
            hit = brute_has_root(k, r_sym, p)
        elif k <= args.poly_max_k:
            route = "sympy-poly"
            hit = sympy_has_root(k, r_sym, p)
        else:
            route = "factorial"
            hit = factorial_has_root(k, r_sym, p)
        used[route] = used.get(route, 0) + 1
        if hit:
            bad.append((k, p, f"{route} finds a root: the column DOES represent m"))
            continue
        if not ours["ok"]:
            bad.append((k, p, f"witness.py rejected what the check accepts: {ours['reason']}"))
            continue
        if (j + 1) % 10 == 0:
            print(f"    {j+1}/{len(idx)}  bad={len(bad)}  "
                  f"{perf_counter()-t0:.0f}s", flush=True)

    print(f"\n  checked      {len(idx)}")
    print(f"  membership routes used: {used}")
    if used.get("factorial"):
        print(f"    NOTE {used['factorial']} used the factorial identity, which is")
        print("         the one the sweep kernel uses. That cross-checks")
        print("         witness.py against a second derivation, but does not")
        print("         independently certify the sweep itself.")
    print(f"  disagreements {len(bad)}")
    for row in bad[:10]:
        print(f"    k={row[0]} p={row[1]}: {row[2]}")
    print(f"  {perf_counter()-t0:.1f}s")
    print("  RESULT", "AGREE" if not bad else "DISAGREE")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
