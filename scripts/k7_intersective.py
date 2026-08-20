#!/usr/bin/env python3
"""Q26: is (x)_7 - c ever intersective?  No -- but the proof has one gap.

Q25 settled k=6 by a trick that is unavailable here. For EVEN k the roots are
symmetric about (k-1)/2, so f_c is a polynomial in t^2 and the degree halves
(6 -> 3). For ODD k the centred polynomial is ODD rather than even:

    (x)_7 = y (y^2-1)(y^2-4)(y^2-9) = y^7 - 14y^5 + 49y^3 - 36y,   y = x-3.

The analogous identity does exist -- f_c(y) f_c(-y) = T(y^2) with
T(v) = c^2 - v S(v)^2, S(v) = (v-1)(v-4)(v-9) -- and the same equivalence
holds (f_c has a root mod p iff T has a SQUARE root mod p). But deg T = 7,
not 3. Odd k gets NO degree reduction, which is exactly why it is the hard
side. So k=7 goes the k=5 way: coefficient matching, and curves.

Oddness does give one thing: f_c(-y) = -f_{-c}(y), so the whole problem is
symmetric under c -> -c, and it suffices to search c >= 0.

Intersective needs f_c reducible over Q with every irreducible factor of
degree >= 2 (else Jordan + Chebotarev supply a prime with no root). The
partitions of 7 into parts >= 2 are 2+5, 3+4 and 2+2+3. But 2+2+3 is a 2+5
with a reducible quintic, so there are only TWO branches to analyse:

  BRANCH A   f_c has a quadratic factor y^2+ay+b.
             Eliminating the quintic's coefficients leaves ONE condition,
             free of c, in (a,b) -- and it involves a only through a^2:

               Phi(A,b) = -A^3 + 5A^2 b + 14A^2 - 6A b^2 - 42A b - 49A
                          + b^3 + 14b^2 + 49b + 36 = 0,     A = a^2.

             Its leading form -A^3+5A^2b-6Ab^2+b^3 is an irreducible cubic
             form, so the curve has 3 distinct points at infinity and
             Siegel's theorem gives FINITELY many integral points.
             Found: c = 0 (trivial) and c = +-17472, +-459648.

  BRANCH B   f_c has a cubic factor and no quadratic factor (3+4).
             Two conditions survive; the first is linear in the cubic's
             constant term d, giving a plane curve C(A,b) of degree 5 whose
             leading form is squarefree with 5 distinct roots -- Siegel
             again. The linear solve has a DEGENERATE locus 3a^2-2b-14 = 0
             which must be handled separately; it is not empty, and it is
             where the only nontrivial 3+4 solution lives.
             Found: c = 0 (trivial) and c = +-896.

Every nontrivial candidate then dies for the same reason: rad(7!) = 210 does
not divide c. In fact 5 divides none of 896, 17472, 459648. Since (x)_7
vanishes identically on F_5 (any 7 consecutive residues contain a multiple of
5), f_c == -c mod 5, which is nonzero. No root mod 5.

THE GAP. Siegel proves both curves have finitely many integral points, and
Baker's method makes that effective, but this script does not carry out the
effective computation -- it searches |a| up to a bound and cross-checks
against brute-force factorisation of (x)_7 - c over a range of c. The two
methods agree exactly. Closing the gap properly means running a Thue solver
(PARI/GP's `thue`) on the two curves. Until then this is "no counterexample
exists below the search bound, and only finitely many can exist anywhere",
not a finished theorem.

    python scripts/k7_intersective.py
    python scripts/k7_intersective.py --a_max 3000 --brute 400000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

x, y, A, a, b, d, c = sp.symbols("x y A a b d c")
BASE = sp.prod([x - i for i in range(7)])
RAD7 = 210  # rad(7!) = 2*3*5*7

PHI = (
    -(A**3) + 5 * A**2 * b + 14 * A**2 - 6 * A * b**2 - 42 * A * b - 49 * A
    + b**3 + 14 * b**2 + 49 * b + 36
)
CURVE34 = (
    2 * A**5 - 9 * A**4 * b - 70 * A**4 + 12 * A**3 * b**2 + 210 * A**3 * b
    + 833 * A**3 - A**2 * b**3 - 98 * A**2 * b**2 - 1225 * A**2 * b - 3792 * A**2
    - 9 * A * b**4 - 168 * A * b**3 - 882 * A * b**2 - 432 * A * b + 4179 * A
    + 4 * b**5 + 112 * b**4 + 1176 * b**3 + 5632 * b**2 + 11620 * b + 7056
)
E2 = (
    a**5 - 4 * a**3 * b - 14 * a**3 + 3 * a**2 * d + 3 * a * b**2 + 28 * a * b
    + 49 * a - 2 * b * d - 14 * d
)
E1 = (
    a**4 * b - a**3 * d - 3 * a**2 * b**2 - 14 * a**2 * b + 4 * a * b * d
    + 14 * a * d + b**3 + 14 * b**2 + 49 * b - d**2 + 36
)
E0 = (
    a**4 * d - 3 * a**2 * b * d - 14 * a**2 * d + 2 * a * d**2 + b**2 * d
    + 14 * b * d + c + 49 * d
)


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: python -O strips assert."""
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def degs_of(cc) -> tuple:
    return tuple(sorted(p.degree() for p, _ in sp.Poly(BASE - cc, x).factor_list()[1]))


def step0() -> None:
    print("  [0] setup: odd k has no degree reduction")
    P = sp.expand(BASE.subs(x, y + 3))
    check(sp.expand(P - (y**7 - 14 * y**5 + 49 * y**3 - 36 * y)) == 0, "centred form")
    check(sp.simplify(P.subs(y, -y) + P) == 0, "centred (x)_7 is ODD in y")
    S = (A - 1) * (A - 4) * (A - 9)
    prod = sp.expand((P - c) * (P.subs(y, -y) - c))
    T = sp.expand(c**2 - A * S**2)
    check(sp.expand(prod - T.subs(A, y**2)) == 0, "f_c(y)f_c(-y) = T(y^2)")
    check(sp.Poly(T, A).degree() == 7, "deg T = 7, so NO degree reduction")
    print("      f_c(y)f_c(-y) = T(y^2), T(v) = c^2 - v S(v)^2, deg T = 7 (not 3)")
    print("      f_c(-y) = -f_{-c}(y), so the problem is symmetric in c -> -c")


def branch_a(a_max: int) -> list:
    print(f"  [A] f_c has a QUADRATIC factor   (covers 2+5 and 2+2+3), |a| <= {a_max}")
    F3 = -(A**3) + 5 * A**2 * b - 6 * A * b**2 + b**3
    fl = sp.factor_list(sp.expand(F3.subs(A, 1)))[1]
    check(len(fl) == 1 and fl[0][1] == 1, "leading cubic form squarefree/irreducible")
    check(sp.Poly(fl[0][0], b).degree() == 3, "3 distinct points at infinity")
    out = {}
    for av in range(0, a_max + 1):
        Av = av * av
        for root in sp.polys.polyroots.roots(sp.Poly(PHI.subs(A, Av), b)):
            if not root.is_Integer:
                continue
            bv = int(root)
            cv = (
                av**5 * bv - 4 * av**3 * bv**2 - 14 * av**3 * bv
                + 3 * av * bv**3 + 28 * av * bv**2 + 49 * av * bv
            )
            if cv:
                out[abs(cv)] = (av, bv)
    check(set(out) == {17472, 459648}, f"branch A nontrivial c, got {sorted(out)}")
    print("      leading form irreducible cubic -> 3 points at infinity -> Siegel: FINITE")
    for cv, (av, bv) in sorted(out.items()):
        print(f"      c = +-{cv:<8} (a={av}, b={bv})  split {degs_of(cv)}")
    return [{"c": int(v), "branch": "A", "split": list(degs_of(v))} for v in sorted(out)]


def branch_b(a_max: int) -> list:
    print(f"  [B] f_c has a CUBIC factor (3+4), |a| <= {a_max}")
    F5 = 2 * A**5 - 9 * A**4 * b + 12 * A**3 * b**2 - A**2 * b**3 - 9 * A * b**4 + 4 * b**5
    fl = sp.factor_list(sp.expand(F5.subs(A, 1)))[1]
    check(all(m == 1 for _, m in fl), "leading quintic form squarefree")
    check(sum(sp.Poly(p, b).degree() for p, _ in fl) == 5, "5 distinct points at infinity")
    out = {}
    # generic locus: 3a^2-2b-14 != 0
    for av in range(0, a_max + 1):
        Av = av * av
        for root in sp.polys.polyroots.roots(sp.Poly(CURVE34.subs(A, Av), b)):
            if not root.is_Integer:
                continue
            bv = int(root)
            den = 3 * Av - 2 * bv - 14
            if den == 0:
                continue
            num = -av * (Av - 3 * bv - 7) * (Av - bv - 7)
            if num % den:
                continue
            dv = num // den
            cv = int(sp.solve(E0.subs({a: av, b: bv, d: dv}), c)[0])
            if cv:
                out[abs(cv)] = ("generic", av, bv, dv)
    # DEGENERATE locus: 3a^2-2b-14 = 0, which forces num = 0 as well
    ndeg = 0
    for av in range(-a_max, a_max + 1):
        if (3 * av * av - 14) % 2:
            continue
        bv = (3 * av * av - 14) // 2
        if -av * (av * av - 3 * bv - 7) * (av * av - bv - 7) != 0:
            continue
        ndeg += 1
        for dv in sp.solve(E1.subs({a: av, b: bv}), d):
            if not dv.is_Integer:
                continue
            dv = int(dv)
            cv = int(sp.solve(E0.subs({a: av, b: bv, d: dv}), c)[0])
            if cv:
                out[abs(cv)] = ("degenerate", av, bv, dv)
    check(set(out) == {896}, f"branch B nontrivial c, got {sorted(out)}")
    print("      leading form squarefree, 5 points at infinity -> Siegel: FINITE")
    print(f"      degenerate locus 3a^2-2b-14=0 handled separately ({ndeg} points);")
    print("      it is NOT empty and is where the only nontrivial 3+4 solution lives")
    for cv, info in sorted(out.items()):
        print(f"      c = +-{cv:<8} {info}  split {degs_of(cv)}")
    return [{"c": int(v), "branch": "B", "split": list(degs_of(v))} for v in sorted(out)]


def step_kill(cands: list) -> None:
    print("  [C] every nontrivial candidate fails rad(7!) = 210")
    for cv in sorted({r["c"] for r in cands}):
        check(cv % RAD7 != 0, f"c={cv} unexpectedly divisible by 210")
        check(cv % 5 != 0, f"c={cv} unexpectedly divisible by 5")
        vals = {int((sp.prod([v - i for i in range(7)]) - cv) % 5) for v in range(5)}
        check(0 not in vals, f"c={cv} has a root mod 5")
        print(f"      c=+-{cv:<8} {sp.factorint(cv)}   5 does not divide it -> no root mod 5")
    print("      ((x)_7 vanishes identically on F_5, so f_c == -c != 0 there)")


def step_brute(lim: int, cands: list) -> dict:
    print(f"  [D] brute-force cross-check: ALL integer c in [0,{lim}]")
    found = []
    for cc in range(0, lim + 1):
        dg = degs_of(cc)
        if len(dg) > 1 and min(dg) >= 2:
            found.append(cc)
    predicted = sorted({r["c"] for r in cands if r["c"] <= lim})
    check(found == predicted, f"brute {found} != curves {predicted}")
    print(f"      reducible with all parts >= 2: {found}")
    print(f"      curve analysis predicted     : {predicted}   AGREE")
    return {"limit": lim, "found": found}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a_max", type=int, default=400)
    ap.add_argument("--brute", type=int, default=20000)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q26: (x)_7 - c is never intersective (one gap, stated) ===")
    step0()
    ca = branch_a(args.a_max)
    cb = branch_b(min(args.a_max, 400))
    cands = ca + cb
    step_kill(cands)
    rb = step_brute(args.brute, cands)
    print()
    print("  RESULT: the only c with f_c reducible and all factors of degree >= 2")
    print("          are c = +-896, +-17472, +-459648, and 5 divides none of them,")
    print("          so none has a root mod 5. Not intersective.")
    print("  GAP:    Siegel gives finiteness on both curves; this script does not")
    print("          carry out the effective (Baker/Thue) computation, so the lists")
    print("          are complete only up to the search bound. PARI/GP `thue` on")
    print("          the two curves would close it.")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "search": "k7_intersective",
                    "claim": "(x)_7 - c is never intersective, modulo effective Thue",
                    "candidates": cands,
                    "brute_force": rb,
                    "a_max": args.a_max,
                    "intersective_found": 0,
                    "gap": "Siegel gives finiteness; effective Thue not carried out",
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
