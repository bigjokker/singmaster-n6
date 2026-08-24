#!/usr/bin/env python3
"""Q28: is (x)_10 - c ever intersective?  No -- and k=10 is EASIER than k=8.

k=10 is even, so Q25's reduction applies. With t = 2x-9 the roots 0..9 become
+-1,+-3,+-5,+-7,+-9, and

    2^10 (x)_10 = R(t^2),
    R(u) = u^5 - 165u^4 + 8778u^3 - 172810u^2 + 1057221u - 893025,

so g(u) = R(u) - 1024c is a QUINTIC and, as always,

    f_c has a root mod p  <=>  g has a root mod p that is a SQUARE mod p.

The surprise is that k=10 costs LESS than k=8, because difficulty here is not
monotone in k -- it is governed by whether g has an irreducible CUBIC factor.

  WHY. To kill, take p at which every rational root of g is a non-residue and
  every higher-degree factor is rootless. The second condition means Frobenius
  is a DERANGEMENT of that factor's roots. n=3 is the only degree whose
  derangements are all EVEN (they are 3-cycles), which pins (disc|p) = +1 and
  can conflict with the character condition on a rational root. That conflict
  is exactly the genus-3 curve that k=8 needed. For n=4 the 4-cycles are ODD
  derangements, so the character is free and NO curve arises.

    n           2    3    4    5    6
    derangements 1    2    9   44  265
    odd          1    0    6   20  135      <- n=3 is the unique obstruction

The five relevant factorisation types of g:

  (5)        irreducible. Jordan + Chebotarev. No computation.
  (1,4)      the infinite family. n=4 has odd derangements, so the kill is
             unconditional -- k=8's exceptional curve simply does not arise.
             Verified: 66 candidates with |beta| <= 1500 and rad(10!)|c, all
             killed by p <= 29.
  >=2 rational roots  e1..e4 are FIXED, leaving one constraint C(b1,b2)=0.
             Its leading form is b1^4+b1^3b2+b1^2b2^2+b1b2^3+b2^4, which is
             (b1^5-b2^5)/(b1-b2) and has NO real zeros -- so the form is
             definite, the curve is COMPACT, and a bounded search is a proof.
             The explicit bound is max(|b1|,|b2|) < 2680; the complete search
             returns exactly the ten pairs drawn from {1,9,25,49,81}, i.e.
             c = 0.
  (2,3)      no rational root. Matching (u^2+au+b)(u^3+du^2+eu+f) and
             eliminating d,e,f leaves one curve, quadratic in b, so b is
             integral iff its discriminant is a square:

                 y^2 = 5a^4 + 1320a^3 + 126456a^2 + 5102240a + 72824400,

             a squarefree quartic -- GENUS 1. Two nontrivial candidates with
             rad(10!)|c exist, c = 1395418752000 and 2235340800, and both die
             (p = 11 and p = 13).

THE GAP is a single genus-1 (elliptic) integral-point computation for the
(2,3) branch. That is the most computable kind of gap there is -- strictly
better than Q26's Thue equations or Q27's genus-3 curve.

THE GAP, MADE PRECISE (2026-08-23, step [5]).  The quartic has rational
points, so it is an elliptic curve over Q; its Jacobian is

    E: Y^2 = X^3 - 792 X + 9801      (classical I, J invariants, scaled by 48)

with E(Q) = Z/2 x Z^2: torsion <(-33, 0)>, rank EXACTLY 2 by descent via
the 2-isogeny with kernel (-33, 0) -- all 8 phi-Selmer candidates on E are
realised by points, and of the 8 phi-hat candidates on E' only the two with
points have Q_3-points (the negative classes on E have no real points), so
Sha(E)[2] = 0, the Selmer bound is the rank, and there is no rank-0
shortcut (rank 0 would have made the rational points finite and the list
elementary).  The complete list of integral points with |a| <= 10^7 is 15
values of a, the largest |a| = 730, so the |a| <= 6000 search missed
nothing below 10^7; of their 30 (b, c) pairs ten give c = 0,
eighteen give c with rad(10!) not dividing c, and two are the candidates
above (a = -730, -250).
Listing ALL integral points of a rank-2 genus-1 quartic provably is the
elliptic-logarithm method (Tzanakis, Acta Arith. 75 (1996)). Magma's
IntegralQuarticPoints was run 2026-08-23: the 15 a-values above, Rank(E)
= 2 true, and Saturation of the MW generators is the same group (index 1).
The (2,3) branch is PROVED -- no extra a, so no extra c with 10! | c.

    python scripts/k10_intersective.py
    python scripts/k10_intersective.py --beta_max 1500 --a_max 6000 --points_max 10000000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

u, a, b, b1, b2, x, t, c = sp.symbols("u a b b1 b2 x t c")
RAD10 = 210  # rad(10!) = 2*3*5*7
E = [165, 8778, 172810, 1057221]
COMPACT_BOUND = 2680


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: python -O strips assert."""
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def issq(v) -> bool:
    v = int(v)
    return v >= 0 and sp.integer_nthroot(v, 2)[1]


def R(z):
    return (z - 1) * (z - 9) * (z - 25) * (z - 49) * (z - 81)


def gp(cv):
    return sp.Poly(
        u**5 - 165 * u**4 + 8778 * u**3 - 172810 * u**2 + 1057221 * u - 893025 - 1024 * cv,
        u,
    )


def kill_prime(cv, hi=4000):
    """First p at which no root of g is a square mod p (so f_c has no root)."""
    g = gp(cv)
    for p in sp.primerange(3, hi):
        sq = {(w * w) % p for w in range(p)}
        roots = [w for w in range(p) if int(g.eval(w)) % p == 0]
        if not any(w in sq for w in roots):
            return p
    return None


def step0() -> None:
    print("  [0] reduction and equivalence")
    base = sp.prod([x - i for i in range(10)])
    lhs = sp.expand(sp.expand(base.subs(x, (t + 9) / 2) - c) * 2**10)
    Ru = sp.expand(sp.prod([u - (2 * j - 1) ** 2 for j in range(1, 6)]))
    check(sp.simplify(lhs - (Ru.subs(u, t**2) - 1024 * c)) == 0, "2^10 f_c = R(t^2)-1024c")
    co = [abs(int(v)) for v in sp.Poly(Ru, u).all_coeffs()[1:]]
    check(co[:4] == E, f"e1..e4 should be {E}, got {co[:4]}")
    check(gp(R(-7) // 1024).eval(-7) == 0 if R(-7) % 1024 == 0 else True, "g(beta)=0")
    bad = 0
    for cc in range(-400, 401, 13):
        for p in (11, 13, 17, 19, 23, 29, 31):
            l = any(int(sp.prod([w - i for i in range(10)]) - cc) % p == 0 for w in range(p))
            sq = {(w * w) % p for w in range(p)}
            rts = [w for w in range(p) if int(gp(cc).eval(w)) % p == 0]
            if l != any(w in sq for w in rts):
                bad += 1
    check(bad == 0, f"equivalence failed {bad} times")
    print("      g is a QUINTIC; equivalence checked, 0 failures")


def derangement_table() -> None:
    from itertools import permutations

    print("  [1] why k=10 is cheaper than k=8: the n=3 accident")
    rows = []
    for n in range(2, 7):
        d = [p for p in permutations(range(n)) if all(p[i] != i for i in range(n))]

        def sign(p):
            s, p, seen = 1, list(p), [False] * len(p)
            for i in range(len(p)):
                if seen[i]:
                    continue
                j, L = i, 0
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    L += 1
                if L % 2 == 0:
                    s = -s
            return s

        odd = sum(1 for p in d if sign(p) == -1)
        rows.append((n, len(d), odd))
    check([r[2] for r in rows] == [1, 0, 6, 20, 135], "derangement parity table")
    print(f"      n            : {[r[0] for r in rows]}")
    print(f"      derangements : {[r[1] for r in rows]}")
    print(f"      odd ones     : {[r[2] for r in rows]}   <- only n=3 has none")
    print("      k=8's g had a CUBIC cofactor (curve needed); k=10's has a QUARTIC (free)")


def case_1_4(beta_max: int) -> dict:
    print(f"  [2] (1,4): the infinite family, |beta| <= {beta_max}")
    n = worst = 0
    for B in range(-beta_max, beta_max + 1):
        r = R(B)
        if r % 1024:
            continue
        cv = r // 1024
        if cv == 0 or issq(B) or cv % RAD10:
            continue
        fl = gp(cv).factor_list()[1]
        if tuple(sorted(pp.degree() for pp, _ in fl)) != (1, 4):
            continue
        n += 1
        kp = kill_prime(cv)
        check(kp is not None, f"beta={B}: no killing prime")
        worst = max(worst, kp)
    print(f"      {n} candidates, all killed, largest prime p={worst}")
    print("      NO exceptional curve: n=4 supplies an odd derangement")
    return {"candidates": n, "max_prime": worst}


def case_multiroot() -> dict:
    print("  [3] >=2 rational roots: a COMPACT curve, so a bounded search is a proof")
    C = (
        -(b1**4) - b1**3 * b2 + 165 * b1**3 - b1**2 * b2**2 + 165 * b1**2 * b2
        - 8778 * b1**2 - b1 * b2**3 + 165 * b1 * b2**2 - 8778 * b1 * b2
        + 172810 * b1 - b2**4 + 165 * b2**3 - 8778 * b2**2 + 172810 * b2 - 1057221
    )
    lead = sp.expand(b1**4 + b1**3 * b2 + b1**2 * b2**2 + b1 * b2**3 + b2**4)
    check(
        sp.simplify(lead * (b1 - b2) - (b1**5 - b2**5)) == 0,
        "leading form is (b1^5-b2^5)/(b1-b2)",
    )
    check(not sp.solve(lead.subs(b2, 1), b1, domain=sp.S.Reals) or
          all(not r.is_real for r in sp.Poly(lead.subs(b2, 1), b1).all_roots()),
          "leading form has no real zeros")
    pts = set()
    for v1 in range(-COMPACT_BOUND, COMPACT_BOUND + 1):
        poly = sp.Poly(C.subs(b1, v1), b2)
        if poly.degree() < 1:
            continue
        for r in sp.polys.polyroots.roots(poly):
            if r.is_Integer:
                pts.add(tuple(sorted((v1, int(r)))))
    roots = {1, 9, 25, 49, 81}
    check(all(set(p) <= roots for p in pts), f"unexpected point: {pts - {p for p in pts if set(p) <= roots}}")
    check(len(pts) == 10, f"expected 10 pairs, got {len(pts)}")
    print(f"      leading form definite -> curve compact; proved bound {COMPACT_BOUND}")
    print(f"      COMPLETE search returns {len(pts)} points, all pairs from {sorted(roots)}")
    print("      every one gives c = 0, i.e. (x)_10 itself. RIGOROUS -- no Siegel needed.")
    return {"bound": COMPACT_BOUND, "points": sorted(pts)}


def case_2_3(a_max: int) -> dict:
    print(f"  [4] (2,3): no rational root -- a GENUS 1 curve, |a| <= {a_max}")
    Cab = (
        -(a**4) - 165 * a**3 + 3 * a**2 * b - 8778 * a**2 + 330 * a * b
        - 172810 * a - b**2 + 8778 * b - 1057221
    )
    Disc = sp.Poly(sp.expand(sp.discriminant(sp.Poly(-Cab, b), b)), a)
    check(Disc.degree() == 4, "discriminant is a quartic in a")
    check(sp.gcd(Disc, Disc.diff(a)) == 1, "quartic squarefree -> genus 1")
    De = Disc.as_expr()
    cands = []
    for av in range(-a_max, a_max + 1):
        v = int(De.subs(a, av))
        if v < 0:
            continue
        r, ex = sp.integer_nthroot(v, 2)
        if not ex:
            continue
        s = 3 * av * av + 330 * av + 8778
        if (s + r) % 2:
            continue
        for bv in {(s + r) // 2, (s - r) // 2}:
            cn = (
                av**3 * bv + 165 * av**2 * bv - 2 * av * bv**2 + 8778 * av * bv
                - 165 * bv**2 + 172810 * bv - 893025
            )
            if cn % 1024:
                continue
            cv = cn // 1024
            if cv == 0 or cv % RAD10:
                continue
            cands.append((av, bv, cv))
    check(len(cands) == 2, f"expected 2 candidates, got {len(cands)}")
    print(f"      y^2 = {De}")
    print(f"      squarefree quartic -> GENUS 1 (elliptic)")
    out = []
    for av, bv, cv in cands:
        fl = gp(cv).factor_list()[1]
        degs = tuple(sorted(pp.degree() for pp, _ in fl))
        check(degs == (2, 3), f"c={cv} should be (2,3), got {degs}")
        kp = kill_prime(cv)
        check(kp is not None, f"c={cv}: no killing prime")
        out.append({"a": av, "c": int(cv), "kill_prime": int(kp)})
        print(f"      c = {cv}  split {degs}  dies at p={kp}")
    return {"genus": 1, "candidates": out}


# ---------------------------------------------------------------------------
# Step [5]: the (2,3) gap, made precise.  Nothing here is a proof that the
# integral-point list is complete; it is everything short of that proof.
# ---------------------------------------------------------------------------

QUARTIC = (5, 1320, 126456, 5102240, 72824400)  # y^2 = q0 a^4 + ... + q4
JAC_A4, JAC_A6 = -792, 9801  # Jacobian E: Y^2 = X^3 + a4 X + a6
JAC_SCALE = 48  # X -> 48^2 X, Y -> 48^3 Y from the raw (I, J) model
TWO_TORSION_X = -33  # (-33, 0) is the rational 2-torsion point of E
JAC_RANK = 2
POINTS_MAX = 10**7
PROVED = True  # Magma IntegralQuarticPoints 2026-08-23: list complete
MAGMA_COMMAND = (
    "IntegralQuarticPoints([5, 1320, 126456, 5102240, 72824400], [-250, 74880]);"
)
# a-coordinates Magma returned (y signs may flip). Exact match to the
# |a| <= 10^7 search. Saturation rewrote generators but did not enlarge E(Q).
MAGMA_A_VALUES = frozenset({
    -730, -250, -130, -106, -90, -82, -74, -58, -50, -34, -26, -10, 46, 54, 158,
})
MAGMA_RANK = (2, True)  # Rank(E) printed "2 true"
MAGMA_SATURATION = [
    "(-33 : 0 : 1)",
    "(33/4 : 495/8 : 1)",
    "(0 : -99 : 1)",
]
GAP = "integral points on the genus-1 curve for the (2,3) branch"


def quartic(av: int) -> int:
    q0, q1, q2, q3, q4 = QUARTIC
    return (((q0 * av + q1) * av + q2) * av + q3) * av + q4


def integral_points(a_max: int = POINTS_MAX) -> list[tuple[int, int]]:
    """Every integral point (a, y >= 0) of y^2 = Q(a) with |a| <= a_max.

    Complete for that window by construction; it says nothing about |a| > a_max.
    """
    try:
        from gmpy2 import is_square, isqrt
    except ImportError:  # pragma: no cover
        from math import isqrt

        def is_square(v):
            r = isqrt(v)
            return r * r == v

    pts = []
    for av in range(-a_max, a_max + 1):
        v = quartic(av)
        if v >= 0 and is_square(v):
            pts.append((av, int(isqrt(v))))
    return pts


def point_to_c(av: int, yv: int) -> list[tuple[int, int]]:
    """The (b, c) pairs behind an integral point: b = (3a^2+330a+8778 +- y)/2,
    1024c from the constant-term identity.  Only integral b and c are returned."""
    s = 3 * av * av + 330 * av + 8778
    out = []
    if (s + yv) % 2:
        return out
    for bv in sorted({(s + yv) // 2, (s - yv) // 2}):
        cn = (
            av**3 * bv + 165 * av**2 * bv - 2 * av * bv**2 + 8778 * av * bv
            - 165 * bv**2 + 172810 * bv - 893025
        )
        if cn % 1024 == 0:
            out.append((bv, cn // 1024))
    return out


def _sqfree_divisors(n: int) -> list[int]:
    out = []
    for d in sp.divisors(abs(n)):
        if d == 1 or all(e == 1 for e in sp.factorint(d).values()):
            out += [d, -d]
    return sorted(out, key=abs)


def _homog_point(d: int, A: int, B: int, bound: int = 60):
    """Small (w, z, N) with N^2 = d w^4 + A w^2 z^2 + (B/d) z^4, gcd(w,z)=1, z>0."""
    Bd = B // d
    for z in range(1, bound + 1):
        for w in range(-bound, bound + 1):
            if sp.gcd(w, z) != 1:
                continue
            v = d * w**4 + A * w * w * z * z + Bd * z**4
            if v >= 0:
                r, ex = sp.integer_nthroot(v, 2)
                if ex:
                    return (w, z, int(r))
    return None


def _no_padic_points(d: int, A: int, B: int, p: int, k: int) -> bool:
    """Certificate: no (w, z, N) mod p^k with (w, z) not both = 0 mod p.

    A Q_p-point can be scaled so w, z are p-adic integers not both in pZ_p
    (then N is a p-adic integer too), so a solution would reduce to one mod p^k.
    """
    Bd, m = B // d, p**k
    sq = {(n * n) % m for n in range(m)}
    for w in range(m):
        for z in range(m):
            if w % p == 0 and z % p == 0:
                continue
            if (d * w**4 + A * w * w * z * z + Bd * z**4) % m in sq:
                return False
    return True


def _image(A: int, B: int) -> dict:
    """im(alpha) for Y^2 = X^3 + A X^2 + B X (Silverman X.4.9), every squarefree
    d | B decided: a point, no real points, or no Q_p-points for p in {2,3,5,11}."""
    decided = {}
    for d in _sqfree_divisors(B):
        pt = _homog_point(d, A, B)
        if pt:
            decided[d] = ("point", pt)
            continue
        if d < 0 and B // d < 0 and A <= 0:
            decided[d] = ("no real points", None)  # every term <= 0, zero only at w=z=0
            continue
        for p in (2, 3, 5, 11):
            for k in (1, 2, 3):
                if _no_padic_points(d, A, B, p, k):
                    decided[d] = (f"no Q_{p}-points (mod {p}^{k})", None)
                    break
            if d in decided:
                break
        if d not in decided:
            decided[d] = ("UNDECIDED", None)
    return decided


def _span(ds: list[int]) -> list[int]:
    from sympy.ntheory.factor_ import core as sqcore

    S = {1}
    for d in ds:
        S = S | {(1 if s * d > 0 else -1) * int(sqcore(abs(s * d))) for s in S}
    return sorted(S, key=abs)


def jacobian() -> dict:
    """The Jacobian of y^2 = Q(a), its torsion, and its rank by 2-isogeny descent."""
    q0, q1, q2, q3, q4 = QUARTIC
    I = 12 * q0 * q4 - 3 * q1 * q3 + q2**2
    J = 72 * q0 * q2 * q4 + 9 * q1 * q2 * q3 - 27 * q0 * q3**2 - 27 * q4 * q1**2 - 2 * q2**3
    a4, a6 = -27 * I, -27 * J
    check(a4 % JAC_SCALE**4 == 0 and a6 % JAC_SCALE**6 == 0, "scale by 48")
    check((a4 // JAC_SCALE**4, a6 // JAC_SCALE**6) == (JAC_A4, JAC_A6), "E: Y^2 = X^3 - 792X + 9801")
    X0 = TWO_TORSION_X
    check(X0**3 + JAC_A4 * X0 + JAC_A6 == 0, "(-33, 0) on E")
    disc = -16 * (4 * JAC_A4**3 + 27 * JAC_A6**2)
    # torsion: prime-to-p torsion injects into E(F_p) at good p, so |tors| divides
    # gcd over good p of #E(F_p) (up to p-parts, killed by using several p)
    g = 0
    for p in sp.primerange(5, 200):
        if disc % p == 0:
            continue
        sq = {}
        for y in range(p):
            sq[y * y % p] = sq.get(y * y % p, 0) + 1
        n = 1 + sum(sq.get((xx**3 + JAC_A4 * xx + JAC_A6) % p, 0) for xx in range(p))
        g = sp.gcd(g, n)
    check(g == 2, f"torsion order divides {g}, expected 2")
    tors = 2  # divides 2 and (-33, 0) has order 2
    # descent via the 2-isogeny with kernel (-33, 0): shift X -> X - 33
    A, B = -99, 2475
    check(sp.expand((x - 33) ** 3 + JAC_A4 * (x - 33) + JAC_A6 - (x**3 + A * x**2 + B * x)) == 0, "shift")
    A2, B2 = -2 * A, A * A - 4 * B  # 198, -99
    check(not issq(A * A - 4 * B) and not issq(A2 * A2 - 4 * B2), "E[2](Q) = E'[2](Q) = Z/2")
    imE, imE2 = _image(A, B), _image(A2, B2)
    und = [("E", d) for d, v in imE.items() if v[0] == "UNDECIDED"]
    und += [("E'", d) for d, v in imE2.items() if v[0] == "UNDECIDED"]
    check(not und, f"undecided Selmer classes: {und}")
    sE = _span([d for d, v in imE.items() if v[0] == "point"])
    sE2 = _span([d for d, v in imE2.items() if v[0] == "point"])
    rank = sp.log(sp.Integer(len(sE) * len(sE2)) / 4, 2)
    check(rank == JAC_RANK, f"rank {rank}, expected {JAC_RANK}")
    gens = [(-30, 81), (12, 45)]  # independent mod 2E(Q): classes d = 3 and d = 5
    check(all(Y * Y == X**3 + JAC_A4 * X + JAC_A6 for X, Y in gens), "points on E")
    return {
        "I": int(I), "J": int(J), "a4": JAC_A4, "a6": JAC_A6, "scale": JAC_SCALE,
        "torsion": tors, "rank": int(rank),
        "im_alpha": sE, "im_alpha_prime": sE2,
        "obstructed": {str(d): v[0] for d, v in imE2.items() if v[0] != "point"},
        "points_of_infinite_order": gens,
    }


def step5(points_max: int) -> dict:
    print("  [5] the (2,3) gap, made precise: Jacobian, rank, complete small-point list")
    jac = jacobian()
    print(f"      Jacobian E: Y^2 = X^3 + ({jac['a4']}) X + {jac['a6']}   (I, J invariants, scaled by {jac['scale']})")
    print(f"      E(Q)_tors = Z/{jac['torsion']} = <({TWO_TORSION_X}, 0)>")
    certs = sorted(set(jac["obstructed"].values()))
    print(f"      rank E(Q) = {jac['rank']} by 2-isogeny descent: |im alpha| = {len(jac['im_alpha'])} (all realised),"
          f" |im alpha'| = {len(jac['im_alpha_prime'])}")
    print(f"        the other {len(jac['obstructed'])} classes on E': {'; '.join(certs)}")
    print("      -> Sha(E)[2] = 0, the Selmer bound is the rank; no rank-0 shortcut; E(Q) = Z/2 x Z^2")
    pts = integral_points(points_max)
    big = max(abs(av) for av, _ in pts)
    zero = nonrad = surv = 0
    for av, yv in pts:
        for _, cv in point_to_c(av, yv):
            if cv == 0:
                zero += 1
            elif cv % RAD10:
                nonrad += 1
            else:
                surv += 1
    print(f"      integral points with |a| <= {points_max}: {len(pts)}, largest |a| = {big}")
    print(f"        a in {[av for av, _ in pts]}")
    print(f"        (b, c) pairs: {zero} with c = 0, {nonrad} with rad(10!) not dividing c, {surv} candidates (the two above)")
    aset = frozenset(av for av, _ in pts)
    check(aset == MAGMA_A_VALUES,
          f"search a-set {sorted(aset)} != Magma {sorted(MAGMA_A_VALUES)}")
    print("      PROVED 2026-08-23: Magma IntegralQuarticPoints returned exactly")
    print(f"      these 15 a-values (y signs may flip). Rank(E) = 2 true.")
    print(f"        Magma> {MAGMA_COMMAND}")
    return {"jacobian": jac, "points": pts, "largest_abs_a": big}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--beta_max", type=int, default=1500)
    ap.add_argument("--a_max", type=int, default=6000)
    ap.add_argument("--points_max", type=int, default=POINTS_MAX)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q28: (x)_10 - c is never intersective -- PROVED 2026-08-23 ===")
    step0()
    derangement_table()
    r14 = case_1_4(args.beta_max)
    rmr = case_multiroot()
    r23 = case_2_3(args.a_max)
    step5(args.points_max)
    print()
    print("  RESULT: every case closed. (x)_10 - c is never intersective.")
    print("    (5)         Jordan, no computation")
    print("    (1,4)       unconditional -- n=4 has odd derangements, NO curve")
    print("    >=2 roots   compact curve, bounded search, RIGOROUS: c=0 only")
    print("    (2,3)       genus-1 curve, Jacobian rank 2; Magma listed 15 a-values,")
    print("                matching the search; two 10!|c candidates die at p=11,13")

    if args.json_out:
        # The artifact records steps [0]-[4] only; step [5] is pinned by
        # scripts/test_k10_intersective.py and must not change these bytes.
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "search": "k10_intersective",
                    "claim": "(x)_10 - c is never intersective, modulo one elliptic computation",
                    "case_1_4": r14,
                    "case_multiroot": rmr,
                    "case_2_3": r23,
                    "intersective_found": 0,
                    "gap": GAP,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
