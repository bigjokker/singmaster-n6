#!/usr/bin/env python3
"""Q25: is (x)_6 - c ever intersective?  No. This verifies every step.

Q14 settled k <= 5 and left k >= 6 open, predicting k = 6 "will take three
curves" -- one per split of the sextic (2+4, 3+3, 2+2+2). It takes less,
because for even k the problem is not really about the sextic at all.

With t = 2x - 5 the roots 0..5 become +-1, +-3, +-5, so

    2^6 (x)_6 = (t^2-1)(t^2-9)(t^2-25) = R(t^2),   R(u) = u^3-35u^2+259u-225

and therefore  2^6 f_c = g(t^2)  with

    g(u) = u^3 - 35u^2 + 259u - (225 + 64c).

THE KEY EQUIVALENCE. Since t ranges over all of F_p and t -> t^2 has image
exactly the squares,

    f_c has a root mod p  <=>  g has a root u mod p that is a SQUARE mod p.

(u = 0 is a square, so it is included.) So intersectivity of the sextic is a
question about a CUBIC plus a quadratic-residue condition, and the case
analysis runs over the factorisation of g over Q -- a trichotomy -- rather
than over the five splits of f_c.

  CASE 1  g irreducible over Q.
          Gal acts transitively on 3 roots, Jordan gives a derangement,
          Chebotarev gives a positive density of p at which g has NO root
          mod p at all. A fortiori no square root. Dead.

  CASE 2  g = (u - beta) * q(u),  q irreducible quadratic, disc Delta.
          beta is an integer (monic, rational root). beta is a non-square,
          else t^2 = beta gives f_c a rational root. Delta is a non-square
          because q is irreducible. Both non-squares => Q(sqrt beta,
          sqrt Delta) supplies a positive density of p with BOTH symbols
          -1. At such p, q has no root mod p, so beta is the only root of
          g mod p -- and it is a non-residue. Dead.
          This is an INFINITE family and needs no finiteness argument.

  CASE 3  g splits completely over Q.
          Then e1 = 35 and e2 = 259 force -3s^2 + 140s - 1036 to be a
          perfect square, a BOUNDED conic. The only solution is
          {beta} = {1,9,25}, i.e. c = 0, which is (x)_6 itself and has
          rational roots. Dead.

So (x)_6 - c is never intersective, matching k <= 5.

Case 3 is where an intersective example would have to live if one existed:
three rational beta_i, none a square, but with beta_1 beta_2 beta_3 a square
so that the all-non-residue pattern is forbidden. That is exactly the "odd
square-relation" of q14-intersective.md section 4, and here beta_1 beta_2
beta_3 = 225 + 64c. The conic simply has no room for it.

    python scripts/k6_intersective.py
    python scripts/k6_intersective.py --brute 300000 --json_out results/k6.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

x, t, u, c = sp.symbols("x t u c")
BASE = sp.prod([x - i for i in range(6)])
RAD6 = 30  # rad(6!) = rad(720) = 2*3*5


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: python -O strips assert."""
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def issq(v) -> bool:
    v = int(v)
    return v >= 0 and sp.integer_nthroot(v, 2)[1]


def g_poly(cc):
    return sp.Poly(u**3 - 35 * u**2 + 259 * u - (225 + 64 * cc), u)


def factors_of(cc) -> list:
    return [p for p, _ in sp.Poly(BASE - cc, x).factor_list()[1]]


def degs_of(cc) -> tuple:
    return tuple(sorted(p.degree() for p in factors_of(cc)))


def kill_prime(cc, hi=20000):
    """First prime at which NO irreducible factor of f_c has a root."""
    facs = factors_of(cc)
    for p in sp.primerange(2, hi):
        if all(not any(int(f.eval(v)) % p == 0 for v in range(p)) for f in facs):
            return p
    return None


def step0() -> None:
    print("  [0] reduction, and the key equivalence")
    gg = sp.expand(sp.expand((BASE - c).subs(x, (t + 5) / 2)) * 64)
    R = sp.expand(sp.prod([u - (2 * j - 1) ** 2 for j in (1, 2, 3)]))
    check(sp.simplify(gg - (R.subs(u, t**2) - 64 * c)) == 0, "2^6 f_c == R(t^2) - 64c")
    check(sp.expand(R) == u**3 - 35 * u**2 + 259 * u - 225, "R(u)")
    # verify the equivalence directly, over many (c, p)
    bad = 0
    for cc in range(-400, 401, 7):
        for p in (7, 11, 13, 17, 19, 23, 29, 31, 37):
            lhs = any(int(sp.prod([v - i for i in range(6)]) - cc) % p == 0 for v in range(p))
            sq = {(v * v) % p for v in range(p)}
            gv = [v for v in range(p) if (v**3 - 35 * v**2 + 259 * v - (225 + 64 * cc)) % p == 0]
            rhs = any(v in sq for v in gv)
            if lhs != rhs:
                bad += 1
    check(bad == 0, f"equivalence failed in {bad} cases")
    print("      2^6 f_c = g(t^2);  f_c has a root mod p  <=>  g has a SQUARE root mod p")
    print("      equivalence checked on 115 x 9 = 1035 (c,p) pairs: 0 failures")


def case3() -> list:
    print("  [CASE 3] g splits completely over Q")
    sols = []
    for s in range(0, 60):
        D = -3 * s * s + 140 * s - 1036
        if D < 0:
            continue
        r, exact = sp.integer_nthroot(D, 2)
        if not exact or (s + r) % 2:
            continue
        b1, b2, b3 = (s + r) // 2, (s - r) // 2, 35 - s
        if b1 * b2 + b3 * (b1 + b2) != 259:
            continue
        e3 = b1 * b2 * b3
        if (e3 - 225) % 64 == 0:
            sols.append(((b1, b2, b3), (e3 - 225) // 64))
    cs = {cc for _, cc in sols}
    check(cs == {0}, f"case 3 must give only c=0, got {cs}")
    check(degs_of(0) == (1,) * 6, "c=0 splits into linear factors")
    print("      e1=35, e2=259  =>  disc = -3s^2+140s-1036 must be a square")
    print("      disc >= 0 only for s in [9.22, 37.44] -- a BOUNDED conic, 28 integers")
    print(f"      solutions: {sorted({tuple(sorted(b)) for b,_ in sols})} -> c = 0 only")
    print("      c=0 is (x)_6 itself: rational roots, so not a candidate. Dead.")
    return [{"c": 0, "case": 3}]


def case2(lo: int, hi: int) -> dict:
    print(f"  [CASE 2] g = (u-beta) * irreducible quadratic   (beta in [{lo},{hi}])")
    n = 0
    worst = 0
    for beta in range(lo, hi + 1):
        num = (beta - 1) * (beta - 9) * (beta - 25)
        if num % 64:
            continue
        cc = num // 64
        if cc == 0:
            continue
        check(g_poly(cc).eval(beta) == 0, f"beta={beta} is a root of g")
        # the complementary quadratic and its discriminant
        Delta = -3 * beta * beta + 70 * beta + 189
        if issq(Delta):
            continue  # g splits completely -> case 3 territory
        if issq(beta):
            continue  # t^2=beta rational -> f_c has a rational root
        if cc % RAD6:
            continue  # rad(6!) | c is necessary for intersectivity
        n += 1
        hit = None
        for p in sp.primerange(5, 20000):
            if beta % p == 0 or Delta % p == 0:
                continue
            if (
                sp.legendre_symbol(beta % p, p) == -1
                and sp.legendre_symbol(Delta % p, p) == -1
            ):
                hit = p
                break
        check(hit is not None, f"beta={beta}: no prime with both symbols -1")
        facs = factors_of(cc)
        check(
            not any(int(f.eval(v)) % hit == 0 for f in facs for v in range(hit)),
            f"beta={beta}: predicted prime {hit} did not kill c={cc}",
        )
        worst = max(worst, hit)
    print(f"      candidates (beta and Delta both non-squares, rad(6!)|c): {n}")
    print(f"      for each, the prime predicted by (beta|p)=(Delta|p)=-1 kills it")
    print(f"      largest predicted prime needed: p = {worst}")
    print("      no finiteness used -- this closes an INFINITE family")
    return {"n": n, "max_prime": int(worst)}


def case1(lim: int) -> dict:
    print("  [CASE 1] g irreducible over Q  -> Jordan + Chebotarev")
    n = 0
    worst = 0
    for cc in range(-lim, lim + 1, RAD6):
        if not g_poly(cc).is_irreducible:
            continue
        n += 1
        if n > 200:
            break
        kp = kill_prime(cc)
        check(kp is not None, f"c={cc}: g irreducible but no killing prime found")
        worst = max(worst, kp)
    print(f"      sampled {n-1 if n>200 else n} values of c with g irreducible")
    print(f"      every one has a killing prime (largest {worst}), as Jordan predicts")
    return {"sampled": n, "max_prime": int(worst)}


def brute(lim: int) -> dict:
    print(f"  [CROSS-CHECK] brute force, c = 30j, |c| <= {lim}")
    seen = {}
    inter = []
    for cc in range(-lim, lim + 1, RAD6):
        dg = degs_of(cc)
        seen[dg] = seen.get(dg, 0) + 1
        if len(dg) > 1 and min(dg) >= 2:
            if kill_prime(cc) is None:
                inter.append(cc)
    check(not inter, f"possible intersective c: {inter}")
    print(f"      patterns: {seen}")
    print(f"      candidates with no killing prime: 0")
    return {"limit": lim, "patterns": {str(k): v for k, v in seen.items()}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--beta_lo", type=int, default=-3000)
    ap.add_argument("--beta_hi", type=int, default=3000)
    ap.add_argument("--brute", type=int, default=60000)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q25: (x)_6 - c is never intersective ===")
    step0()
    r3 = case3()
    r2 = case2(args.beta_lo, args.beta_hi)
    r1 = case1(args.brute)
    rb = brute(args.brute)
    print()
    print("  RESULT: the trichotomy on g is exhaustive and every branch is closed.")
    print("    g irreducible -> Jordan/Chebotarev kills it")
    print("    g = linear x quadratic -> (beta|p)=(Delta|p)=-1 kills it (infinite family)")
    print("    g splits -> bounded conic -> c=0 only, which has rational roots")
    print("  Hence (x)_6 - c is NEVER intersective, extending Q14 from k<=5 to k<=6.")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "search": "k6_intersective",
                    "claim": "(x)_6 - c is never intersective",
                    "method": "trichotomy on the factorisation of g(u)=R(u)-64c",
                    "case1_g_irreducible": r1,
                    "case2_linear_times_quadratic": r2,
                    "case3_g_splits": r3,
                    "brute_force": rb,
                    "intersective_found": 0,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
