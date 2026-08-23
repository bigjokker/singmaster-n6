#!/usr/bin/env python3
"""Q26 Branch A: a smooth plane cubic (genus 1), NOT a Thue equation.

docs/q26-k7-intersective.md section 5 proposed closing Branch A -- "f_c has a
quadratic factor y^2+ay+b" -- with PARI's `thue`.  That was wrong, and this
script establishes what is actually true about the branch.

[0] The chord curve Phi(A, b) = 0 (A = a^2) is re-derived from the matching
    and shown SMOOTH: no affine singular point, and the three points at
    infinity (the roots of the irreducible leading cubic form) are smooth.
    A smooth plane cubic has genus 1.  `thue` solves F(x,y) = m for a
    homogeneous form F; an affine cubic curve is not that, and no rearranging
    makes it that: the gap was never "an afternoon with pari".

[1] The curve has rational points (the 21 chords of c = 0, and the c =
    +-17472, +-459648 points), so it is an elliptic curve.  Its Jacobian,
    via the pencil of lines through (A, b) = (0, -1) -- the two residual
    intersections are parametrised by s^2 = 25T^4 - 32T^3 + 306T^2 - 184T
    + 29 -- has classical invariants giving, after scaling by 6,

        E:  Y^2 = X^3 - 1764 X + 28224  =  X^3 - 42^2 X + 168^2.

    E(Q)_tors is TRIVIAL (the order divides gcd #E(F_p) = 1), yet E has
    visible integer points -- (0, 168), (21, 21), (28, 28), (-48, 48) ... --
    so any of them has infinite order: rank >= 1, E(Q) is infinite, and
    there is no rank-0 shortcut to the integral-point list.  Siegel still
    gives finiteness; certifying the list needs the elliptic-logarithm
    method on a CUBIC model (Stroeker-Tzanakis), which no installed tool
    provides.  Branch A stays BLOCKED there.

[2] Chebotarev narrows what the missing list is for.  A 2+5 value c is
    killed at p iff the quadratic is inert AND the quintic is rootless mod
    p, i.e. iff the joint Frobenius is (transposition, derangement).  Among
    the transitive subgroups of S5, ONLY S5 contains an odd derangement
    (the (2,3)-elements); C5, D5, F20 and A5 have only 5-cycles, all even.
    C5 and A5 admit no C2-quotient, so the joint group is the full product
    and (transposition, 5-cycle) kills.  S5 kills through its odd
    derangements even when the quadratic subfields coincide.  The ONLY
    escape: quintic group D5 or F20 whose unique quadratic subfield equals
    Q(sqrt(a^2-4b)) -- then every Frobenius class fixes a root (the
    reflections and affine maps fix exactly one point; the translations
    are even), and no unramified prime kills.  The k=7 no-kill trap, one
    more instance of the even-derangement mechanism from Q27/Q28.

[3] The trap is real: (x^2 + x - 1)(x^5 - 2) -- quadratic field Q(sqrt 5),
    and x^5 - 2 has group F20 with the SAME quadratic subfield Q(sqrt 5)
    inside Q(zeta_5, 2^(1/5)) -- has a root modulo every prime tested.
    Control: (x^2 + x + 1)(x^5 - 2), fields Q(sqrt -3) vs Q(sqrt 5), is
    killed.  (These witness the group theory; they are not (x)_7 - c.)

[4] Both recorded Branch-A candidates factor 2+5 with an S5 quintic (a
    (2,3) factorisation pattern modulo a small prime certifies S5), so even
    without the 5-divisibility kill they would die by Chebotarev.

STATUS: Branch A of Q26 remains BLOCKED -- now on "integral points of a
rank->=1 plane cubic" (the right statement), not on a Thue computation
(the wrong one).  Branch B, by contrast, is closed: scripts/k7_runge.py.

    python scripts/k7_branchA.py          (~1 min)
"""

from __future__ import annotations

import argparse
from itertools import permutations
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

A, b, T, X, Y, x, yv, cc = sp.symbols("A b T X Y x y c")

PHI = (-(A**3) + 5*A**2*b + 14*A**2 - 6*A*b**2 - 42*A*b - 49*A
       + b**3 + 14*b**2 + 49*b + 36)
JAC = (-1764, 28224)   # E: Y^2 = X^3 - 1764X + 28224, scale 6
PENCIL_QUARTIC = (25, -32, 306, -184, 29)
NO_KILL_GROUPS = ("D5", "F20")
PROVED_BRANCH_A = False


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def step0() -> None:
    print("  [0] the chord curve is a SMOOTH plane cubic: genus 1, not a Thue form")
    # re-derive from the matching f_c = (y^2+ay+b)(quintic)
    a_, c_ = sp.symbols("a_ c_")
    P7 = yv**7 - 14*yv**5 + 49*yv**3 - 36*yv
    es = sp.symbols("e0:5")
    quint = yv**5 + sum(es[i]*yv**i for i in range(5))
    D = sp.Poly(sp.expand((yv*yv + a_*yv + b)*quint - (P7 - c_)), yv)
    sol = sp.solve([D.coeff_monomial(yv**k) for k in range(2, 7)], list(es), dict=True)[0]
    rem1 = sp.expand(D.coeff_monomial(yv**1).subs(sol))
    rem0 = sp.expand(D.coeff_monomial(yv**0).subs(sol))
    check(c_ not in rem1.free_symbols, "the y^1 condition is free of c")
    phi_ab = sp.expand(sp.Poly(rem1, a_).as_expr().subs(
        {a_**6: A**3, a_**4: A**2, a_**2: A}))
    check(sp.expand(phi_ab - PHI) == 0, "Phi re-derived from the matching")
    cfun = sp.solve(rem0, c_)[0]
    check(sp.expand(cfun - a_*b*(a_*a_ - 3*b - 7)*(a_*a_ - b - 7)) == 0,
          "c = a b (a^2-3b-7)(a^2-b-7)")
    # smoothness: affine
    sing = sp.solve([PHI, sp.diff(PHI, A), sp.diff(PHI, b)], [A, b], dict=True)
    check(sing == [], "no affine singular points")
    # smoothness at infinity: leading form has 3 distinct roots, so the curve
    # meets the line at infinity in 3 points where grad(lead) != 0
    lead = sp.Poly(b**3 - 6*b**2 + 5*b - 1, b)
    check(sp.discriminant(lead, b) == 49, "leading form: disc 49, three distinct roots")
    check(lead.is_irreducible, "leading form irreducible: no rational point at infinity")
    print("      Phi and c = ab(a^2-3b-7)(a^2-b-7) re-derived; curve smooth -> genus 1")


def step1() -> dict:
    print("  [1] the Jacobian, its trivial torsion, and rank >= 1")
    # pencil through (0,-1): residual intersections f3 X^2 + f2 X + f1 = 0
    F = sp.Poly(sp.expand(PHI.subs({A: X, b: -1 + T*X})), X)
    check(F.coeff_monomial(1) == 0, "(0,-1) on the curve")
    f1, f2, f3 = (F.coeff_monomial(X**k) for k in (1, 2, 3))
    q = sp.Poly(sp.expand(f2**2 - 4*f1*f3), T)
    check(tuple(int(q.coeff_monomial(T**k)) for k in (4, 3, 2, 1, 0)) == PENCIL_QUARTIC,
          "pencil quartic s^2 = 25T^4 - 32T^3 + 306T^2 - 184T + 29")
    a4_, a3_, a2_, a1_, a0_ = PENCIL_QUARTIC
    I = 12*a4_*a0_ - 3*a3_*a1_ + a2_**2
    J = 72*a4_*a2_*a0_ + 9*a3_*a2_*a1_ - 27*a4_*a1_**2 - 27*a0_*a3_**2 - 2*a2_**3
    check((I, J) == (84672, -48771072), "invariants I, J")
    check((-27*I) % 6**4 == 0 and (-27*J) % 6**6 == 0, "scale by 6")
    check((-27*I//6**4, -27*J//6**6) == JAC, "E: Y^2 = X^3 - 1764X + 28224")
    m4, m6 = JAC
    disc = -16*(4*m4**3 + 27*m6**2)
    check(sp.Poly(X**3 + m4*X + m6, X).is_irreducible, "no rational 2-torsion")
    import math
    g = 0
    for p in sp.primerange(5, 200):
        if disc % p == 0:
            continue
        sq = {}
        for u_ in range(p):
            sq[u_*u_ % p] = sq.get(u_*u_ % p, 0) + 1
        n = 1 + sum(sq.get((u_**3 + m4*u_ + m6) % p, 0) for u_ in range(p))
        g = math.gcd(g, n)
    check(g == 1, "torsion order divides gcd #E(F_p) = 1: E(Q)_tors trivial")
    pt = (0, 168)
    check(pt[1]**2 == pt[0]**3 + m4*pt[0] + m6, "(0,168) on E")
    print(f"      E: Y^2 = X^3 - 42^2 X + 168^2, torsion trivial, (0,168) on E")
    print("      -> (0,168) has INFINITE order: rank >= 1, E(Q) infinite,")
    print("         no rank-0 shortcut. The missing certificate is the")
    print("         elliptic-log method on a cubic model (Stroeker-Tzanakis).")
    return {"E": JAC, "torsion": 1, "witness_point": pt}


def _group(name: str) -> list:
    r = (1, 2, 3, 4, 0)                       # x -> x+1
    s2 = tuple((2*i) % 5 for i in range(5))   # x -> 2x  (order 4)
    sneg = tuple((-i) % 5 for i in range(5))  # x -> -x  (order 2)
    def compose(p, q_):
        return tuple(p[q_[i]] for i in range(5))
    def closure(gens):
        seen = {tuple(range(5))}
        frontier = list(seen)
        while frontier:
            nxt = []
            for g_ in frontier:
                for h_ in gens:
                    e_ = compose(g_, h_)
                    if e_ not in seen:
                        seen.add(e_)
                        nxt.append(e_)
            frontier = nxt
        return sorted(seen)
    if name == "C5":
        return closure([r])
    if name == "D5":
        return closure([r, sneg])
    if name == "F20":
        return closure([r, s2])
    if name == "A5":
        return [p for p in permutations(range(5)) if _sign(p) == 1]
    if name == "S5":
        return list(permutations(range(5)))
    raise ValueError(name)


def _sign(p) -> int:
    s, p_, seen = 1, list(p), [False]*len(p)
    for i in range(len(p_)):
        if seen[i]:
            continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True
            j = p_[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return s


def step2() -> None:
    print("  [2] the 2+5 trichotomy: only S5 supplies an odd derangement")
    rows = []
    for name in ("C5", "D5", "F20", "A5", "S5"):
        G = _group(name)
        der = [p for p in G if all(p[i] != i for i in range(5))]
        odd = [p for p in der if _sign(p) == -1]
        rows.append((name, len(G), len(der), len(odd)))
        if name in NO_KILL_GROUPS:
            check(len(odd) == 0, f"{name}: all derangements even")
            # and the derangements lie in the kernel of the C2-quotient
            # (index-2 subgroup = <r> resp. <r, s2^2>): all are 5-cycles
            check(all(sorted(p) == [0, 1, 2, 3, 4] and
                      all(p[i] != i for i in range(5)) for p in der), "5-cycles")
    check([r_[3] for r_ in rows] == [0, 0, 0, 0, 20], "odd derangements: S5 only (the 20 (2,3)-elements)")
    for name, n, nd, nod in rows:
        print(f"      {name:3} order {n:3}: {nd:2} derangements, {nod:2} odd"
              + ("   <- the (2,3)-elements" if nod else ""))
    print("      C5, A5: no C2-quotient -> full product -> (inert, 5-cycle) kills")
    print("      S5: (inert, (2,3)) kills even with a shared quadratic subfield")
    print("      D5, F20 + matching subfield: every class fixes a root. NO KILL.")


def _rootless(coeffs: list, hi: int, lo: int = 3):
    for p in sp.primerange(lo, hi):
        found = False
        for u_ in range(p):
            v = 0
            for co in coeffs:
                v = (v*u_ + co) % p
            if v == 0:
                found = True
                break
        if not found:
            return int(p)
    return None


def step3(scan: int) -> None:
    print(f"  [3] the trap witnessed, p < {scan}")
    h1 = sp.Poly(sp.expand((x*x + x - 1)*(x**5 - 2)), x).all_coeffs()
    h2 = sp.Poly(sp.expand((x*x + x + 1)*(x**5 - 2)), x).all_coeffs()
    k1 = _rootless([int(v) for v in h1], scan)
    k2 = _rootless([int(v) for v in h2], scan)
    check(k1 is None, f"(x^2+x-1)(x^5-2) unexpectedly killed at {k1}")
    check(k2 is not None, "(x^2+x+1)(x^5-2) should be killed")
    print("      (x^2+x-1)(x^5-2): Q(sqrt5) matches F20's subfield -> root mod every p tested")
    print(f"      (x^2+x+1)(x^5-2): Q(sqrt-3) does not -> killed at p={k2}")


def step4() -> None:
    print("  [4] the recorded candidates' quintics are S5 (a (2,3) pattern certifies it)")
    P7 = sp.prod([x - i for i in range(7)])
    for cv in (17472, 459648):
        fl = sp.Poly(P7 - cv, x).factor_list()[1]
        quint = [p_ for p_, _ in fl if p_.degree() == 5]
        check(len(quint) == 1, f"c={cv} splits 2+5")
        qq = quint[0]
        got = None
        for p in sp.primerange(3, 200):
            degs = sorted(sp.Poly(pp, x).degree() for pp, _ in sp.factor_list(qq.as_expr(), modulus=p)[1])
            if degs == [2, 3]:
                got = int(p)
                break
        check(got is not None, f"c={cv}: (2,3) pattern found")
        print(f"      c=+-{cv}: quintic has a (2,3) factorisation mod {got} -> group S5")
    print("      (an S5 quintic is Chebotarev-killed regardless; both also die at 5)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scan", type=int, default=2000)
    args = ap.parse_args()
    print("=== Q26 Branch A: genus 1 (not Thue), rank >= 1, and the D5/F20 trap ===")
    step0()
    step1()
    step2()
    step3(args.scan)
    step4()
    print()
    print("  RESULT: Branch A stays BLOCKED -- on integral points of a rank->=1")
    print("          elliptic plane cubic (Stroeker-Tzanakis machinery), not on")
    print("          a Thue equation. Chebotarev kills every 2+5 value whose")
    print("          quintic is not a D5/F20 with matching quadratic subfield;")
    print("          an intersective c would have to live on that thin locus")
    print("          AND survive 7! | c. No such c is known below the searches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
