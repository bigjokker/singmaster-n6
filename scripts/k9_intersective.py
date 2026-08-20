#!/usr/bin/env python3
"""Q29: is (x)_9 - c ever intersective?  No -- modulo the usual effective step.

k=9 is ODD, so Q25's degree-halving does not apply (see Q26): with y = x-4,

    (x)_9 = y(y^2-1)(y^2-4)(y^2-9)(y^2-16) = y^9 - 30y^7 + 273y^5 - 820y^3 + 576y

is ODD in y, which gives only the symmetry f_c(-y) = -f_{-c}(y) -- enough to
restrict to c >= 0, not enough to reduce the degree.

BRANCH REDUCTION. The partitions of 9 into parts >= 2 are 9, 2+7, 3+6, 4+5,
2+2+5, 2+3+4, 3+3+3, 2+2+2+3 -- eight shapes. But a shape containing a 2 is
just "f_c has a quadratic factor", and one containing a 3 (but no 2) is "f_c
has a cubic factor". So there are only THREE branches:

  A  a quadratic factor      covers 2+7, 2+2+5, 2+3+4, 2+2+2+3
  B  a cubic, no quadratic   covers 3+6, 3+3+3
  C  a quartic, no quad/cubic covers 4+5

(and irreducible, which Jordan + Chebotarev kill with no computation).

BRANCH A. Eliminating the septic's coefficients leaves ONE condition, free of
c, involving a only through a^2 -- the chord curve for P(y1)=P(y2). Its
constant part factors as -(b+1)(b+4)(b+9)(b+16), the same shape as k=7's
-(w+1)(w+4)(w+9). The leading form is squarefree with 4 distinct roots, so
there are >= 3 points at infinity and Siegel gives finiteness. Nontrivial
integer points:

    c = +-176774400   (2,7)   fails rad(9!) = 210
    c = +-2630880     (2,7)   PASSES rad(9!) -- the first candidate anywhere
                              in this ladder to do so -- and dies at p = 13.

BRANCH B. Here the surviving y^2 condition is QUADRATIC in the cubic's
constant term, not linear, so the elimination needs a resultant rather than a
substitution. It has degree 18 in a and 9 in b. Every integer point found has
c = 0, so 3+6 and 3+3+3 contribute nothing.

BRANCH C. Three conditions in four unknowns. p0 is linear in the first, and
eliminating p1 by resultant gives a curve with TWO components:

    (2p2 - 3p3^2 + 30)^4    -- the degenerate locus, where the p0-solve is
                               invalid and which must be handled separately
    a main component of degree 14 in p2, 24 in p3

The degenerate locus is empty and the main component's integer points all have
c = 0, so 4+5 contributes nothing either.

A STRUCTURAL NOTE ON 4+5, worth recording even though the branch is empty.
Killing a 4+5 needs Frobenius to be a derangement of BOTH factors at once. If
the two splitting fields were independent that has density
(9/24)(44/120) ~ 0.14 > 0, and checking the sign character across all 9
transitive subgroups of S4 and all 20 of S5 finds NO blocked pair. But the
sign character is not the only shared quotient: C4 and F20 can fuse over C4,
and there

    C4's derangements avoid the identity coset,
    F20's derangements lie entirely inside it,

so the two requirements are incompatible and the kill is blocked. That is a
SECOND obstruction mechanism, distinct from the n=3 one that produced Q27's
genus-3 curve. It is moot here only because no 4+5 solution exists at all.

THE GAP is the same as Q26/Q27/Q28: Siegel gives finiteness on the branch
curves, but the effective computation was not carried out. Searches cover
|a| <= 300 (A), |a| <= 200 (B), |p3| <= 120 (C), cross-checked by factoring
(x)_9 - c for every multiple of 210 up to 4,200,000, which returns exactly the
one case the curves predict.

    python scripts/k9_intersective.py
    python scripts/k9_intersective.py --a_max 300 --p3_max 120 --brute 20000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

x, y, a, b, d, c = sp.symbols("x y a b d c")
p0, p1, p2, p3 = sp.symbols("p0 p1 p2 p3")
BASE = sp.prod([x - i for i in range(9)])
P = y**9 - 30 * y**7 + 273 * y**5 - 820 * y**3 + 576 * y
RAD9 = 210  # rad(9!) = 2*3*5*7


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: python -O strips assert."""
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def degs_of(cv) -> tuple:
    return tuple(sorted(p.degree() for p, _ in sp.Poly(BASE - cv, x).factor_list()[1]))


def kill_prime(cv, hi=20000):
    facs = [f for f, _ in sp.Poly(BASE - cv, x).factor_list()[1]]
    for p in sp.primerange(2, hi):
        if all(not any(int(f.eval(v)) % p == 0 for v in range(p)) for f in facs):
            return p
    return None


def step0() -> None:
    print("  [0] odd k: no degree reduction")
    check(sp.expand(BASE.subs(x, y + 4) - P) == 0, "centred form")
    check(sp.simplify(P.subs(y, -y) + P) == 0, "centred (x)_9 is ODD")
    print("      f_c(-y) = -f_{-c}(y): restricts to c >= 0, does NOT halve the degree")
    ps = [q for q in sp.utilities.iterables.partitions(9) if all(k >= 2 for k in q)]
    check(len(ps) == 8, f"expected 8 partitions, got {len(ps)}")
    print(f"      {len(ps)} shapes collapse to THREE branches (a 2 subsumes, a 3 subsumes)")


def branch_a(a_max: int) -> list:
    print(f"  [A] a quadratic factor  (2+7, 2+2+5, 2+3+4, 2+2+2+3), |a| <= {a_max}")
    Phi = (
        -(a**8) + 7 * a**6 * b + 30 * a**6 - 15 * a**4 * b**2 - 150 * a**4 * b
        - 273 * a**4 + 10 * a**2 * b**3 + 180 * a**2 * b**2 + 819 * a**2 * b
        + 820 * a**2 - b**4 - 30 * b**3 - 273 * b**2 - 820 * b - 576
    )
    check(
        sp.expand(sp.factor(-(b**4 + 30 * b**3 + 273 * b**2 + 820 * b + 576))
                  + (b + 1) * (b + 4) * (b + 9) * (b + 16)) == 0,
        "constant part is -(b+1)(b+4)(b+9)(b+16)",
    )
    A = sp.symbols("A")
    L = -(A**4) + 7 * A**3 * b - 15 * A**2 * b**2 + 10 * A * b**3 - b**4
    fl = sp.factor_list(sp.expand(-L.subs(A, 1)))[1]
    check(all(m == 1 for _, m in fl), "leading form squarefree")
    check(sum(sp.Poly(f, b).degree() for f, _ in fl) == 4, "4 distinct points at infinity")
    out = {}
    for av in range(0, a_max + 1):
        poly = sp.Poly(Phi.subs(a, av), b)
        if poly.is_zero or poly.degree() < 1:
            continue
        for r in poly.ground_roots():
            if not r.is_Integer:
                continue
            bv = int(r)
            cv = (
                av**7 * bv - 6 * av**5 * bv**2 - 30 * av**5 * bv + 10 * av**3 * bv**3
                + 120 * av**3 * bv**2 + 273 * av**3 * bv - 4 * av * bv**4
                - 90 * av * bv**3 - 546 * av * bv**2 - 820 * av * bv
            )
            if cv:
                out[abs(cv)] = (av, bv)
    check(set(out) == {2630880, 176774400}, f"branch A nontrivial c: {sorted(out)}")
    print("      leading form squarefree, 4 points at infinity -> Siegel: FINITE")
    res = []
    for cv in sorted(out):
        kp = kill_prime(cv)
        rad = cv % RAD9 == 0
        check(degs_of(cv) == (2, 7), f"c={cv} should be (2,7)")
        if rad:
            check(kp is not None, f"c={cv} passes rad but has no killing prime")
        res.append({"c": int(cv), "rad_ok": rad, "kill_prime": kp})
        note = "PASSES rad(9!)" if rad else "fails rad(9!)"
        print(f"      c = +-{cv:<10} (2,7)  {note}, dies at p={kp}")
    return res


def branch_b(a_max: int) -> dict:
    print(f"  [B] a cubic, no quadratic  (3+6, 3+3+3), |a| <= {a_max}")
    R2 = (
        a**7 - 6 * a**5 * b - 30 * a**5 + 5 * a**4 * d + 10 * a**3 * b**2
        + 120 * a**3 * b + 273 * a**3 - 12 * a**2 * b * d - 90 * a**2 * d
        - 4 * a * b**3 - 90 * a * b**2 - 546 * a * b + 3 * a * d**2 - 820 * a
        + 3 * b**2 * d + 60 * b * d + 273 * d
    )
    R1 = (
        a**6 * b - a**5 * d - 5 * a**4 * b**2 - 30 * a**4 * b + 8 * a**3 * b * d
        + 30 * a**3 * d + 6 * a**2 * b**3 + 90 * a**2 * b**2 + 273 * a**2 * b
        - 3 * a**2 * d**2 - 9 * a * b**2 * d - 120 * a * b * d - 273 * a * d
        - b**4 - 30 * b**3 - 273 * b**2 + 3 * b * d**2 - 820 * b + 30 * d**2 - 576
    )
    R0 = -(
        a**6 * d - 5 * a**4 * b * d - 30 * a**4 * d + 4 * a**3 * d**2
        + 6 * a**2 * b**2 * d + 90 * a**2 * b * d + 273 * a**2 * d - 6 * a * b * d**2
        - 60 * a * d**2 - b**3 * d - 30 * b**2 * d - 273 * b * d + d**3 - 820 * d
    )
    check(sp.Poly(R2, d).degree() == 2, "y^2 is QUADRATIC in d -- needs a resultant")
    Res = sp.expand(sp.resultant(sp.Poly(R2, d), sp.Poly(R1, d)) / -3)
    check(sp.Poly(Res, a).degree() == 18, "resultant degree 18 in a")
    n = nz = 0
    for av in range(0, a_max + 1):
        Pb = sp.Poly(Res.subs(a, av), b)
        if Pb.is_zero or Pb.degree() < 1:
            continue
        for r in Pb.ground_roots():
            if not r.is_Integer:
                continue
            bv = int(r)
            for dr in sp.Poly(R2.subs({a: av, b: bv}), d).ground_roots():
                if not dr.is_Integer:
                    continue
                dv = int(dr)
                if sp.expand(R1.subs({a: av, b: bv, d: dv})) != 0:
                    continue
                n += 1
                if int(R0.subs({a: av, b: bv, d: dv})) != 0:
                    nz += 1
    check(nz == 0, f"branch B has {nz} nontrivial points")
    print(f"      resultant is degree 18 x 9; {n} integer points, ALL with c = 0")
    return {"points": n, "nontrivial": nz}


def branch_c(p3_max: int) -> dict:
    print(f"  [C] a quartic, no quad/cubic  (4+5), |p3| <= {p3_max}")
    q = sp.symbols("q0:5")
    F = P - c
    Q4 = y**4 + p3 * y**3 + p2 * y**2 + p1 * y + p0
    Q5 = y**5 + sum(q[i] * y**i for i in range(5))
    D = sp.Poly(sp.expand(Q4 * Q5 - F), y)
    co = [D.coeff_monomial(y**i) for i in range(9)]
    sol = sp.solve([co[i] for i in range(8, 3, -1)], list(q), dict=True)[0]
    E3 = sp.expand(co[3].subs(sol))
    E2 = sp.expand(co[2].subs(sol))
    E1 = sp.expand(co[1].subs(sol))
    E0 = sp.expand(co[0].subs(sol))
    den = E3.coeff(p0, 1)
    num = sp.expand(-(E3 - p0 * den))
    E2b = sp.expand(sp.numer(sp.together(E2.subs(p0, sp.cancel(num / den)))))
    E1b = sp.expand(sp.numer(sp.together(E1.subs(p0, sp.cancel(num / den)))))
    Res = sp.expand(sp.resultant(sp.Poly(E2b, p1), sp.Poly(E1b, p1)))
    comps = sp.factor_list(Res)[1]
    degen = [f for f, _ in comps if sp.Poly(f, p2).degree() == 1]
    big = [f for f, _ in comps if sp.Poly(f, p2).degree() > 1][0]
    check(len(degen) == 1, "expected one degenerate component")
    print(f"      degenerate locus {degen[0]} = 0  (handled separately)")
    print(f"      main component: degree {sp.Poly(big,p2).degree()} in p2, "
          f"{sp.Poly(big,p3).degree()} in p3")
    rest = sp.expand(E3 - p0 * den)
    ndeg = 0
    for v3 in range(-200, 201):
        if (3 * v3 * v3 - 30) % 2:
            continue
        v2 = (3 * v3 * v3 - 30) // 2
        if sp.expand(rest.subs({p2: v2, p3: v3})) != 0:
            continue
        ndeg += 1
    check(ndeg == 0, f"degenerate locus not empty: {ndeg}")
    n = nz = 0
    for v3 in range(0, p3_max + 1):
        Pp = sp.Poly(big.subs(p3, v3), p2)
        if Pp.is_zero or Pp.degree() < 1:
            continue
        for r in Pp.ground_roots():
            if not r.is_Integer:
                continue
            v2 = int(r)
            dv = int(den.subs({p2: v2, p3: v3}))
            if dv == 0:
                continue
            A_ = sp.Poly(E2b.subs({p2: v2, p3: v3}), p1)
            B_ = sp.Poly(E1b.subs({p2: v2, p3: v3}), p1)
            if A_.is_zero or B_.is_zero:
                continue
            for rr in sp.gcd(A_, B_).ground_roots():
                if not rr.is_Integer:
                    continue
                v1 = int(rr)
                nv = int(num.subs({p1: v1, p2: v2, p3: v3}))
                if nv % dv:
                    continue
                v0 = nv // dv
                if sp.expand(E1.subs({p0: v0, p1: v1, p2: v2, p3: v3})) != 0:
                    continue
                n += 1
                if int(sp.solve(E0.subs({p0: v0, p1: v1, p2: v2, p3: v3}), c)[0]) != 0:
                    nz += 1
    check(nz == 0, f"branch C has {nz} nontrivial points")
    print(f"      degenerate locus EMPTY; main component {n} points, ALL with c = 0")
    return {"points": n, "nontrivial": nz}


def obstruction_note() -> dict:
    """The (C4, F20) blocking configuration for 4+5 -- moot here, but real."""
    from itertools import permutations

    print("  [C'] a second obstruction mechanism, distinct from n=3")

    def close(gens, n):
        G = {tuple(range(n))}
        fr = [tuple(range(n))]
        while fr:
            xx = fr.pop()
            for g in gens:
                z = tuple(xx[g[i]] for i in range(n))
                if z not in G:
                    G.add(z)
                    fr.append(z)
        return G

    C4 = close([(1, 2, 3, 0)], 4)
    F20 = close([(1, 2, 3, 4, 0), (0, 2, 4, 1, 3)], 5)
    C5 = close([(1, 2, 3, 4, 0)], 5)
    D4 = {g for g in C4 if all(g[i] != i for i in range(4))}
    D5 = {g for g in F20 if all(g[i] != i for i in range(5))}
    check(len(F20) == 20 and len(C4) == 4, "group orders")
    check(tuple(range(4)) not in D4, "C4 derangements avoid the identity")
    check(D5 <= C5, "F20 derangements lie in C5, the kernel of F20 -> C4")
    print("      C4's derangements avoid the identity coset of C4;")
    print("      F20's derangements lie entirely INSIDE it (they generate C5).")
    print("      Fused over C4 the two demands are incompatible -> kill blocked.")
    print("      Moot for k=9 (branch C is empty) but a real second mechanism.")
    return {"blocking_pair": "C4 x F20 fused over C4"}


def brute(lim: int) -> dict:
    print(f"  [D] cross-check: factor (x)_9 - c for every c = 210j, 0 <= j <= {lim}")
    found = []
    for j in range(lim + 1):
        cv = 210 * j
        dg = degs_of(cv)
        if len(dg) > 1 and min(dg) >= 2:
            found.append(cv)
    check(found == [2630880], f"brute force found {found}")
    print(f"      reducible with all parts >= 2: {found}   AGREES with the curves")
    return {"limit": 210 * lim, "found": found}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a_max", type=int, default=300)
    ap.add_argument("--b_max", type=int, default=200)
    ap.add_argument("--p3_max", type=int, default=120)
    ap.add_argument("--brute", type=int, default=20000)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q29: (x)_9 - c is never intersective (one gap, stated) ===")
    step0()
    ra = branch_a(args.a_max)
    rb = branch_b(args.b_max)
    rc = branch_c(args.p3_max)
    ro = obstruction_note()
    rd = brute(args.brute)
    print()
    print("  RESULT: the only c with f_c reducible and all factors of degree >= 2")
    print("          are c = +-2630880 and +-176774400, both (2,7). The second")
    print("          fails rad(9!); the first passes it and dies at p = 13.")
    print("  GAP:    Siegel gives finiteness on all three branch curves; the")
    print("          effective computation was not carried out.")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "search": "k9_intersective",
                    "claim": "(x)_9 - c is never intersective, modulo effective Siegel",
                    "branch_A": ra,
                    "branch_B": rb,
                    "branch_C": rc,
                    "second_obstruction": ro,
                    "brute_force": rd,
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
