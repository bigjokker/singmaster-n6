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
    there is no rank-0 shortcut to the integral-point list.  Magma Rank(E)
    on 2026-08-24 printed "3 true" (torsion order 1), so rank = 3
    unconditionally.  Siegel still gives finiteness; certifying the Phi-list
    needs the elliptic-logarithm method on a CUBIC model (Stroeker-Tzanakis),
    which no installed tool provides.  Branch A stays BLOCKED there.

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

[5] (2026-08-24) THEOREM: the REAL half of Branch A closes elementarily.
    If a^2-4b > 0 the two roots are real, and a real chord of P7 forces both
    endpoints inside [-Y0, Y0], Y0 = 3.1042... -- beyond that P7 exceeds its
    largest local maximum (95.8419... < 96) and is monotone.  So |a| <= 6 and
    |b| <= 9, and enumerating that box gives exactly the 21 chord points, ALL
    with c = 0.  Consequence: every surviving Branch-A point has a^2-4b < 0,
    and the entire c = 0 trap is now a theorem rather than an obstacle.

[6] THEOREM: Phi is neither a Thue form nor a Desboves cubic.  The old
    argument here ("smooth, hence not Thue") was INVALID -- N(x,y) = m z^3 is
    itself a smooth plane cubic.  The correct test: a cubic is affinely a
    binary-form equation only if the Hessian of its homogenisation vanishes
    identically on the line at infinity.  For Phi it is 392(b-2A)(2b-A)^2,
    not identically 0 and not even divisible by the leading form, so no point
    at infinity is an inflection and no affine change of variables (which
    preserves the line at infinity) can produce a Thue equation.  A projective
    change would be a different affine chart and would not preserve Phi(Z).
    And psi_3 of E is irreducible, so E has no rational 3-isogeny and Phi is
    not Q-equivalent to a Desboves cubic either.

[7] The explicit birational map sigma: Phi -> E, X = NX/A^2, Y = NY/A^3, with
    sigma^-1(O) = (0,-9) an AFFINE rational point.  This is what rules out the
    obvious shortcut: sigma sends the genuine integral point (A,b) = (4,13) --
    carrying the real value c = +-17472 -- to (1345/4, 48959/8), which is on E
    but NOT integral.  IntegralPoints(E) therefore returns a demonstrably wrong
    list: of the 21 known integral points of Phi (all with |A| <= 95), 8 have
    integral image, 10 have finite non-integral image, and 3 are poles of
    sigma (A = 0).
    The reason is Galois: the three points at infinity form a single C_3-orbit,
    so every rational effective divisor supported there has degree divisible by
    3, and no Q-rational function of degree 1 or 2 has poles only at infinity.
    Over K = Q(theta) the orbit splits and there IS an integral quartic model;
    Magma has no number-field IntegralQuarticPoints.  The block is arithmetic.

[8] THEOREM: 7! | c removes F20 from the no-kill door.  |c| >= 5040 exceeds
    every critical value of P7, so f_c has one real root; by [5] the quadratic
    is complex, so the quintic carries it, giving two conjugate pairs, an even
    complex conjugation, and disc(quintic) > 0 -- a REAL quadratic subfield,
    while the door demands the imaginary Q(sqrt(a^2-4b)).  Only D5 survives,
    and D5 < A5 makes "disc(quintic) is a perfect square" a new necessary
    condition.

STATUS: Branch A of Q26 remains BLOCKED -- on integral points of a rank-3
elliptic plane cubic (Magma Rank 3 true, 2026-08-24), with Thue, Desboves,
Runge, and every Q-rational Weierstrass/quartic model now excluded BY PROOF
rather than by inspection.  Branch B is closed: scripts/k7_runge.py.

    python scripts/k7_branchA.py          (~1 s, measured)
    python scripts/k7_branchA.py --magma  (the block to paste, if wanted)
"""

from __future__ import annotations

import argparse
from itertools import permutations
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

A, b, T, X, Y, x, yv, cc = sp.symbols("A b T X Y x y c")
z = sp.Symbol("z")          # the homogenising variable for step [6]

PHI = (-(A**3) + 5*A**2*b + 14*A**2 - 6*A*b**2 - 42*A*b - 49*A
       + b**3 + 14*b**2 + 49*b + 36)
JAC = (-1764, 28224)   # E: Y^2 = X^3 - 1764X + 28224, scale 6
PENCIL_QUARTIC = (25, -32, 306, -184, 29)
NO_KILL_GROUPS = ("D5", "F20")
PROVED_BRANCH_A = False

# ---- [5] the real half (2026-08-24): a^2-4b > 0 forces c = 0 ----
CRIT_BOUND = 96          # every |critical value| of P7 is < 96 (max 95.8419...)
REAL_A_MAX, REAL_B_MAX = 6, 9     # from |y_i| <= Y0 = 3.1042..., P7(Y0) = 96
REAL_HALF_CLOSED = True           # the only real-chord points are the 21 c = 0 ones

# ---- [6] Phi is neither a Thue form nor a Desboves cubic (both PROVED) ----
# Hessian of the homogenised Phi, restricted to z = 0.  A cubic is affinely
# equivalent over Q to a binary-form equation N(x,y) = m ONLY IF this vanishes
# IDENTICALLY; here it does not, and it is not even divisible by the leading
# form, so no point at infinity is an inflection.
HESSIAN_AT_INFINITY = "392*(b - 2A)*(2b - A)^2"
NOT_THUE = True
PSI3 = (3, 0, -10584, 338688, -3111696)   # 3-division polynomial of E
NOT_DESBOVES = True      # PSI3 irreducible over Q: E has no rational 3-isogeny

# ---- [7] the explicit birational map sigma: Phi -> E ----
# X = NX(A,b)/A^2, Y = NY(A,b)/A^3; sigma^{-1}(O) = (A,b) = (0,-9).
SIGMA_POLE = (0, -9)
SIGMA_DENOMS = ("A^2", "A^3")
NX = (5*A**3 - 25*A**2*b - 22*A**2 + 30*A*b**2 + 127*A*b + 97*A
      - 5*b**3 - 30*b**2 - 45*b - 20)
NY = (-8*A**4 + 65*A**3*b - 14*A**3 - 173*A**2*b**2 - 196*A**2*b - 23*A**2
      + 158*A*b**3 + 833*A*b**2 + 1192*A*b + 517*A
      - 25*b**4 - 175*b**3 - 375*b**2 - 325*b - 100)
# a genuine integral point of Phi whose image is NOT integral on E:
SIGMA_WITNESS = ((4, 13), (sp.Rational(1345, 4), sp.Rational(48959, 8)), 17472)
# Recorded, not pinned: the O_K-quartic identity is not re-derived here.
# kappa = b - theta*A over K = Q(theta), theta^3-6theta^2+5theta-1 = 0.
K_QUARTIC = ("(-3t^2+12t+16) k^4 + (-56t^2+252t+168) k^3 + (-294t^2+1470t+588) k^2"
             " + (-432t^2+3100t+652) k + (385t^2+1246t+385)")
K_QUARTIC_LEADING_IS_SQUARE = "(-6t^2+33t-14)^2"

# ---- [8] F20 falls out of the no-kill door once 7! | c ----
SURVIVING_NO_KILL_GROUP = "D5"

# ---- [9] known integral points of Phi, and Magma rank ----
# Complete list with |A| <= 95, found by factoring Phi(A, b) in b over Z.
# A larger sweep is not in this repo; do not cite |a| <= 3162 as coverage.
PHI_INTEGRAL_POINTS = (
    (-95, -66), (-44, -15), (-39, -30), (-15, -6), (-11, -12), (-2, -6),
    (0, -1), (0, -4), (0, -9),
    (1, 0), (1, -2), (1, -6),
    (4, 13), (4, 0), (4, -3),
    (9, 38), (9, 2), (9, 0),
    (13, 4), (16, 3), (25, 6),
)
SEARCH_A_MAX = 95
SEARCH_POINTS = 21
# sigma images of those 21: integral on E, finite non-integral, poles (A = 0)
SIGMA_IMAGE_COUNTS = (8, 10, 3)
SEARCH_C_VALUES = (0, 17472, 459648)   # |c| over the A = a^2 subset
# Magma 2026-08-24: TorsionSubgroup order 1; Rank(E) printed "3 true".
MAGMA_COMMAND = (
    "E := EllipticCurve([-1764, 28224]); TorsionSubgroup(E); Rank(E);"
)
MAGMA_RANK = (3, True)
MAGMA_TORSION_ORDER = 1
RANK_LOWER = 3
RANK_GENERATORS = ((-48, 48), (60, 372), (490, 10808))  # on E; Magma certifies rank 3

BLOCKED_ON = (
    "the elliptic-logarithm method on a CUBIC model (Stroeker-Tzanakis, Math. "
    "Comp. 72 (2003) 1917-1933).  No Magma routine does it: IntegralPoints and "
    "SIntegralPoints take Weierstrass models, IntegralQuarticPoints takes "
    "y^2 = quartic over Z, SIntegralDesbovesPoints takes Desboves cubics -- and "
    "Phi is none of these (steps [6], [7])."
)

MAGMA_BLOCK = """// Q26 Branch A -- Jacobian data.  Rank already run 2026-08-24:
//   E := EllipticCurve([-1764, 28224]);
//   TorsionSubgroup(E);   // Abelian Group of order 1
//   Rank(E);              // 3 true
// This does NOT close Branch A (no Magma routine lists Phi(Z)).
E := EllipticCurve([-1764, 28224]);
E;
TorsionSubgroup(E);
Rank(E);
// Demonstration only -- do not treat IntegralPoints as the Phi-list.
Q := E ! [1345/4, 48959/8];
Q;
Denominator(Q[1]);            // 4: sigma(4,13) is not integral on E
"""


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


P7Y = yv**7 - 14*yv**5 + 49*yv**3 - 36*yv     # (x)_7 shifted by 3: y(y^2-1)(y^2-4)(y^2-9)


def step5() -> None:
    """The real half of Branch A closes, elementarily."""
    print("  [5] THEOREM: a^2 - 4b > 0 forces c = 0 (the real half is finished)")
    crit = sp.real_roots(sp.Poly(sp.diff(P7Y, yv), yv))
    check(len(crit) == 6, "P7' has six real critical points")
    check(all(abs(sp.N(P7Y.subs(yv, r), 40)) < CRIT_BOUND for r in crit),
          f"every |critical value| < {CRIT_BOUND}")
    # |c| > every critical value => P7(y) = c has ONE real root, so no real chord
    check(len(sp.real_roots(sp.Poly(P7Y - CRIT_BOUND, yv))) == 1,
          "P7(y) = 96 has a single real root")
    Y0 = max(sp.N(r, 30) for r in sp.real_roots(sp.Poly(P7Y - CRIT_BOUND, yv)))
    check(int(sp.floor(2*Y0)) == REAL_A_MAX and int(sp.floor(Y0**2)) == REAL_B_MAX,
          f"|a| <= {REAL_A_MAX}, |b| <= {REAL_B_MAX}")
    hits = []
    for av in range(-REAL_A_MAX, REAL_A_MAX + 1):
        for bv in range(-REAL_B_MAX, REAL_B_MAX + 1):
            if av*av - 4*bv <= 0:
                continue
            if sp.expand(PHI.subs({A: av*av, b: bv})) == 0:
                hits.append((av, bv, av*bv*(av*av - 3*bv - 7)*(av*av - bv - 7)))
    chords = {(-(i + j), i*j) for i in range(-3, 4) for j in range(i + 1, 4)}
    check({(h[0], h[1]) for h in hits} == chords and len(hits) == 21,
          "exactly the 21 chord points")
    check({h[2] for h in hits} == {0}, "every real-chord value is c = 0")
    print(f"      a real chord of P7 needs both roots in [-Y0, Y0], Y0 = {sp.N(Y0, 8)}")
    print(f"      (beyond it P7 exceeds its largest local max {CRIT_BOUND}, and P7 is")
    print(f"      monotone there), so |a| <= {REAL_A_MAX} and |b| <= {REAL_B_MAX}.")
    print("      Enumerating that box: exactly the 21 chords, ALL with c = 0.")
    print("      => every surviving Branch-A point has a^2 - 4b < 0: the quadratic")
    print("         factor is complex, and the whole c = 0 trap is now a theorem.")


def step6() -> None:
    """Phi is not a Thue form and not a Desboves cubic -- both proved."""
    print("  [6] Phi is NEITHER a Thue form NOR a Desboves cubic (proved, not assumed)")
    F = sp.expand(PHI.subs({A: A/z, b: b/z}) * z**3)
    H = sp.Matrix(3, 3, lambda i, j: sp.diff(F, [A, b, z][i], [A, b, z][j]))
    hess0 = sp.Poly(sp.expand(H.det().subs(z, 0)), A, b)
    check(sp.expand(hess0.as_expr() - 392*(b - 2*A)*(2*b - A)**2) == 0,
          f"Hessian|_(z=0) = {HESSIAN_AT_INFINITY}")
    check(hess0.as_expr() != 0, "Hessian|_(z=0) is NOT identically zero")
    lead = sp.Poly(-(A**3) + 5*A**2*b - 6*A*b**2 + b**3, A, b)
    _, r = sp.div(hess0, lead)
    check(sp.expand(r.as_expr()) != 0, "not even divisible by the leading form")
    # a Thue curve N(x,y) = m z^3 has Hessian|_(z=0) == 0 identically; control:
    xx, yy = sp.symbols("xx yy")
    G = xx**3 + 2*yy**3 - 7*z**3
    HG = sp.Matrix(3, 3, lambda i, j: sp.diff(G, [xx, yy, z][i], [xx, yy, z][j]))
    check(sp.expand(HG.det().subs(z, 0)) == 0,
          "control: a genuine Thue cubic has Hessian|_(z=0) == 0")
    m4, m6 = JAC
    psi3 = sp.Poly(3*X**4 + 6*m4*X**2 + 12*m6*X - m4**2, X)
    check(tuple(int(v) for v in psi3.all_coeffs()) == PSI3, "3-division polynomial")
    check(psi3.is_irreducible and not psi3.ground_roots(),
          "psi_3 irreducible: E has NO rational 3-isogeny")
    print(f"      Hessian at infinity = {HESSIAN_AT_INFINITY}, not identically 0 and")
    print("      not divisible by the leading form: no point at infinity is an")
    print("      inflection, so NO affine change of variables turns Phi into")
    print("      N(x,y) = m.  PARI's `thue` can never apply -- that is now a")
    print("      theorem, not the (invalid) inference 'smooth, hence not Thue'.")
    print("      psi_3 irreducible => no rational 3-isogeny => Phi is not")
    print("      Q-equivalent to a Desboves cubic either, so Magma's")
    print("      SIntegralDesbovesPoints is out as well.")


def sigma_of(Av: int, bv: int):
    """Image of (A,b) under sigma, or None at a pole (A = 0)."""
    if Av == 0:
        return None
    return (sp.Rational(int(NX.subs({A: Av, b: bv})), Av**2),
            sp.Rational(int(NY.subs({A: Av, b: bv})), Av**3))


def classify_sigma_images(pts=PHI_INTEGRAL_POINTS) -> tuple[int, int, int]:
    n_int = n_non = n_pole = 0
    m4, m6 = JAC
    for Av, bv in pts:
        img = sigma_of(Av, bv)
        if img is None:
            n_pole += 1
            continue
        Xv, Yv = img
        check(Yv**2 == Xv**3 + m4*Xv + m6, f"sigma({Av},{bv}) on E")
        if Xv.q == 1 and Yv.q == 1:
            n_int += 1
        else:
            n_non += 1
    return (n_int, n_non, n_pole)


def step7() -> None:
    """The birational map, and why no Q-model hands Phi's problem to Magma."""
    print("  [7] the explicit map sigma: Phi -> E, and the model obstruction")
    m4, m6 = JAC
    W = sp.expand(NY**2 - NX**3 - m4*NX*A**4 - m6*A**6)
    _, r = sp.div(sp.Poly(W, b), sp.Poly(PHI, b))
    check(sp.expand(r.as_expr()) == 0, "sigma lands on E identically modulo Phi")
    (Av, bv), (Xv, Yv), cw = SIGMA_WITNESS
    check(sp.expand(PHI.subs({A: Av, b: bv})) == 0, "(4,13) is on Phi")
    check(sigma_of(Av, bv) == (Xv, Yv), "sigma(4,13) = (1345/4, 48959/8)")
    check(Yv**2 == Xv**3 + m4*Xv + m6, "that image really lies on E")
    check(Xv.q != 1 and Yv.q != 1, "but it is NOT an integral point of E")
    a_w = sp.integer_nthroot(Av, 2)[0]
    check(a_w*bv*(a_w*a_w - 3*bv - 7)*(a_w*a_w - bv - 7) == cw,
          f"and it carries the real Branch-A value c = {cw}")
    counts = classify_sigma_images()
    check(counts == SIGMA_IMAGE_COUNTS,
          f"sigma images of the 21 points: {counts} != {SIGMA_IMAGE_COUNTS}")
    print(f"      X = NX/{SIGMA_DENOMS[0]},  Y = NY/{SIGMA_DENOMS[1]};  "
          f"sigma^-1(O) = {SIGMA_POLE}, an AFFINE rational point.")
    print(f"      Witness: the genuine Phi-point (A,b) = {(Av, bv)} (c = +-{cw})")
    print(f"      maps to ({Xv}, {Yv}) -- on E, NOT integral.  So IntegralPoints(E)")
    print("      provably returns a WRONG list.  Of the 21 known integral points")
    print(f"      of Phi (all |A| <= {SEARCH_A_MAX}): {counts[0]} integral image,")
    print(f"      {counts[1]} finite non-integral, {counts[2]} poles (A = 0).")
    print("      This is Q28's Sage trap, at k=7.")
    print("      The reason is Galois: the three points at infinity of Phi form one")
    print("      C_3-orbit (leading form irreducible, disc 49), so every rational")
    print("      effective divisor supported there has degree divisible by 3.  No")
    print("      Q-rational function of degree 1 or 2 has poles only at infinity,")
    print("      hence NO Q-rational Weierstrass or quartic model has the SAME")
    print("      integral points as Phi.  A model whose integral points merely")
    print("      CONTAIN them exists -- rescale (X,Y) -> (u^2 X, u^3 Y) -- but")
    print("      choosing u needs the points' denominators, i.e. the answer.")
    print("      Over K = Q(theta) the orbit splits and the obstruction vanishes:")
    print(f"      kappa = b - theta A gives the INTEGRAL quartic s^2 = {K_QUARTIC},")
    print(f"      leading coefficient {K_QUARTIC_LEADING_IS_SQUARE}.  That is exactly")
    print("      IntegralQuarticPoints' input shape -- but over O_K, and Magma has")
    print("      no number-field version.  The block is arithmetic, not geometric.")


def step8() -> None:
    """7! | c removes F20 from the no-kill door."""
    print("  [8] THEOREM: 7! | c leaves only D5 in the no-kill door (F20 is out)")
    check(5040 > CRIT_BOUND, "7! exceeds every critical value of P7")
    for cv in (5040, -5040):
        check(len(sp.real_roots(sp.Poly(P7Y - cv, yv))) == 1,
              f"f_c has exactly one real root at c = {cv}")
    # a quintic with one real root has two conjugate pairs -> disc > 0
    for poly in (x**5 - 2, x**5 + x + 1, x**5 - x - 1):
        p_ = sp.Poly(poly, x)
        check(len(sp.real_roots(p_)) == 1 and sp.discriminant(p_, x) > 0,
              "control: one real root => disc > 0")
    P7x = sp.prod([x - i for i in range(7)])
    for cv in (17472, 459648):
        q5 = [p_ for p_, _ in sp.Poly(P7x - cv, x).factor_list()[1]
              if p_.degree() == 5][0]
        d5 = sp.discriminant(q5, x)
        check(len(sp.real_roots(q5)) == 1 and d5 > 0,
              f"c={cv}: quintic has one real root, disc > 0")
        check(not sp.integer_nthroot(int(d5), 2)[1],
              f"c={cv}: disc not a square -> not D5 (it is S5, step [4])")
    print("      An intersective c needs 7! | c, so |c| >= 5040 > every critical")
    print("      value of P7: f_c has exactly ONE real root.  By [5] the surviving")
    print("      region has a^2-4b < 0, so the quadratic contributes none and the")
    print("      QUINTIC has that single real root -- two conjugate pairs, so")
    print("      complex conjugation is a (2,2)-element, hence EVEN, hence")
    print("      disc(quintic) > 0 and Q(sqrt disc) is REAL.  F20's unique")
    print("      quadratic subfield is Q(sqrt disc) (F20 cap A5 = D5), but the door")
    print("      demands it equal Q(sqrt(a^2-4b)), which is IMAGINARY.  So F20 is")
    print("      impossible and only D5 survives -- and D5 < A5 forces")
    print("      disc(quintic) to be a perfect SQUARE, a new necessary condition.")
    print("      Both recorded candidates fail it (they are S5, as step [4] found).")


def step9() -> None:
    """What is known, what is searched, and exactly what is missing."""
    print("  [9] the known list, Magma rank 3 true, and the named wall")
    m4, m6 = JAC
    for (Xv, Yv) in RANK_GENERATORS:
        check(Yv*Yv == Xv**3 + m4*Xv + m6, f"({Xv},{Yv}) on E")
    check(len(PHI_INTEGRAL_POINTS) == SEARCH_POINTS, "21 recorded points")
    check(max(abs(Av) for Av, _ in PHI_INTEGRAL_POINTS) == SEARCH_A_MAX,
          f"every recorded point has |A| <= {SEARCH_A_MAX}")
    for Av, bv in PHI_INTEGRAL_POINTS:
        check(sp.expand(PHI.subs({A: Av, b: bv})) == 0, f"({Av},{bv}) on Phi")
    check(MAGMA_RANK == (3, True) and MAGMA_TORSION_ORDER == 1,
          "Magma Rank(E) was 3 true, torsion order 1")
    print(f"      {SEARCH_POINTS} integral points of Phi, all |A| <= {SEARCH_A_MAX}")
    print("      (listed in PHI_INTEGRAL_POINTS; a larger sweep is not in this repo).")
    print(f"      A = a^2 subset gives |c| in {set(SEARCH_C_VALUES)},")
    print("      i.e. c = 0 and the two recorded values, both dying at 5.")
    print(f"      Magma 2026-08-24: {MAGMA_COMMAND}")
    print(f"      TorsionSubgroup order {MAGMA_TORSION_ORDER}; Rank(E) printed "
          f"{MAGMA_RANK[0]} {'true' if MAGMA_RANK[1] else 'false'}.")
    print("      MISSING: " + BLOCKED_ON.replace("  ", " "))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scan", type=int, default=2000)
    ap.add_argument("--magma", action="store_true", help="print the Magma block")
    args = ap.parse_args()
    if args.magma:
        print(MAGMA_BLOCK)
        return 0
    print("=== Q26 Branch A: not Thue, not Desboves, real half closed, still BLOCKED ===")
    step0()
    step1()
    step2()
    step3(args.scan)
    step4()
    step5()
    step6()
    step7()
    step8()
    step9()
    print()
    print("  RESULT: Branch A stays BLOCKED -- on integral points of a rank-3")
    print("          elliptic plane cubic (Magma Rank 3 true; Stroeker-Tzanakis")
    print("          machinery), not on a Thue equation and not on a Desboves")
    print("          cubic: both are now excluded by proof, and no Q-rational")
    print("          Weierstrass or quartic model preserves integrality (one")
    print("          Galois orbit at infinity).  NEW: the real half a^2-4b > 0")
    print("          is CLOSED -- only the 21 chords, all c = 0 -- so every")
    print("          survivor has a complex quadratic factor; and 7! | c removes")
    print("          F20, leaving D5 alone in the no-kill door.  21 known")
    print("          integral points of Phi, all |A| <= 95. PROVED_BRANCH_A = False.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
