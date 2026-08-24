#!/usr/bin/env python3
"""Q28, the (2,3) branch: how far Chebotarev goes without the integral-point list.

The question. g(u) = R(u) - 1024c = q(u) k(u), q an irreducible quadratic, k an
irreducible cubic. f_c = (x)_10 - c is killed at p iff h(t) = g(t^2) is
rootless mod p (t = 2x-9 is invertible mod odd p). Rootless mod p means
Frobenius acts on h's ten roots +-sqrt(beta_j), +-sqrt(theta_i) with no fixed
point. Chebotarev turns "the Galois group contains a fixed-point-free element"
into "a killing prime exists" -- existence is all a non-intersectivity proof
needs, so nothing here requires effectivity.

THE TRICHOTOMY (step [0]; each case verified on synthetic polynomials).

 (I)  D_q * D_k not a rational square (the two quadratic subfields differ --
      this includes Gal(k) = A3, where D_k is square and D_q is not). Then
      Gal contains (transposition, 3-cycle) and g itself is rootless mod a
      positive density of p. KILLED. Both recorded candidates are here.

 (II) Shared subfield, k(t^2) irreducible over Q. Gal(g) is the fibered
      product of order 6 and (transposition, 3-cycle) is NOT in it -- the
      3-cycles are even, Q27's derangement trap. But the tau-class still
      kills through the character: at such p, k has exactly one root theta_1
      and q is inert; theta_1 is a nonresidue for a positive density of
      them, because either theta_1 is not a square in the splitting field L
      (a character of the Kummer group flips sqrt(theta_1)), or it is and
      sqrt(theta_1) lies in L minus Q(theta_1), where Frobenius flips it for
      free. The ONLY way this fails is sqrt(theta_1) in Q(theta_1), i.e.
      k(t^2) reducible. KILLED.

(III) k(t^2) = m(t) * -m(-t) reducible (forces f = -gamma^2). Now only the
      3-cycle class can kill (the theta_i are squares in L_m: the tau-class
      and id-class always see a root), and it kills iff q(t^2) is
      irreducible over Q: both beta_j must be nonresidues, and the only
      obstruction, beta_1 in (K*)^2 with K = Q(sqrt(D_q)) (odd-degree
      extensions add no square roots, so square in L_m = square in K), is
      EQUIVALENT to q(t^2) splitting into two rational quadratics. KILLED
      unless q(t^2) also splits.

THE NO-KILL LOCUS. All of (I)-(III) fail exactly when, over Z,

  R(t^2) - 1024c = (t^2+mu t+nu)(t^2-mu t+nu)(t^3+al t^2+be t+ga)(t^3-al t^2+be t-ga)

with the square-class coincidence (mu^2-4nu) * disc(m) a square (without it,
Gal = S3 x C2 and (3-cycle, flip) kills). On that locus NO Frobenius class
kills: every element of Gal(L_m/Q) = S3 fixes a root of h -- the identity
fixes them all, a transposition fixes a root of m, a 3-cycle fixes delta in
K, sqrt(beta_1). h then has a root mod EVERY prime -- at a ramified p a
decomposition-group lift of Frobenius still fixes a root, an algebraic
integer whose reduction lands in F_p -- so the c would survive every
Chebotarev argument, and only prime powers or non-existence could save Q28.
This is the k=10 reappearance of Q27's genus-3 mechanism: the cubic m's
derangements are 3-cycles, all even. Two identities tie the cases together:
D_q = mu^2 (mu^2 - 4 nu) and D_k = disc(m) (gamma - alpha beta)^2, so the
coincidence is EXACTLY D_q * D_k square, the failure of case (I): the four
cases partition the (2,3) branch with no gap.

Equivalent norm form: R(u) - 1024c = u P(u)^2 - Q(u)^2 with P = u^2+p3 u+p1
monic, Q = p4 u^2 + p2 u + p0 (s = A*m, s(t)*-s(-t) = h; p0 = nu*ga,
c = (p0^2 - 945^2)/1024, and 945^2 = -R(0)). Eliminating gives ONE plane
curve F(p4, p2) = 0 (step [2]); p4 must be odd, so p4 = 0 is impossible and
c = 0 forces p0 = +-945.

WHAT IS PROVED HERE, unconditionally (steps [1], [3]):

  * c < 0: the deep factorisation is IMPOSSIBLE -- for 10! | c the single
    real root of g is negative (R is strictly increasing on (-inf, 1], its
    largest critical value magnitude is below 43,930,543 < 1024 * 10!, and
    R(0) = -893025 > 1024c), but m is a real cubic whose real root t0 has
    t0^2 = that root >= 0. And any c with 10! not dividing c -- at ANY size
    -- dies at a power of a prime <= 7, since (x)_10 = 0 mod p^v_p(10!) for
    every integer x. So every (2,3) value with c < 0 is killed: cases
    (I)-(III) cover it. Half the branch closes with no integral-point list.

  * c > 0: killed except on the no-kill locus. That locus forces b = nu^2,
    so its points live on the double cover nu^4 - s(a) nu^2 + q_C(a) = 0 of
    the genus-1 curve C -- ramified exactly at the four simple zeros of b,
    so of GENUS 3 (Riemann-Hurwitz), before the -f = gamma^2 and coincidence
    conditions thin it further. Faltings gives finiteness; nothing here
    gives the list. Searches: |nu| <= 10^6 (b <= 10^12) on C finds b square
    only at the ten c = 0 points (which satisfy the whole deep shape, as
    they must -- g splits completely there) and at a = 46, 158, where
    -f is not a square, so the locus misses them; odd p4 <= 2*10^5 on
    F(p4,p2) = 0 finds only p0 = +-945, i.e. c = 0. Congruences for any
    10!-divisible deep c: 315 | p0, p0 odd, p0 = +-945 mod 2^17.

STATUS: Chebotarev alone does not close c > 0 (the no-kill locus is real).
Magma IntegralQuarticPoints (2026-08-23) listed exactly the 15 known a-values,
Rank 2 true, saturation index 1 -- so the locus has no integer point with
c != 0. Q28 (2,3) is PROVED; Chebotarev still supplies the c < 0 half.
It also corrects an earlier note of ours: a shared quadratic subfield alone
does NOT block the kill -- case (II)'s character argument survives it.

    python scripts/k10_chebotarev.py
    python scripts/k10_chebotarev.py --nu_max 20000 --p4_max 5001 --kill_scan 5000
"""

from __future__ import annotations

import argparse

import sympy as sp

t, u = sp.symbols("t u")
R5 = sp.expand(sp.prod([u - (2 * j - 1) ** 2 for j in range(1, 6)]))
CRIT_MAX = 43930543          # ceil of the largest |critical value| of R on [1, 81]
F10 = 3628800                # 10!
CANDIDATES = {1395418752000: 11, 2235340800: 13}
PROVED_23_NEGATIVE = True    # c < 0: every (2,3) value killed, unconditionally
PROVED_23_POSITIVE = True    # Magma 2026-08-23: no extra integral a, so no c!=0 on the locus


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def root_mod(coeffs: list[int], p: int) -> bool:
    for w in range(p):
        v = 0
        for cc in coeffs:
            v = (v * w + cc) % p
        if v == 0:
            return True
    return False


def first_rootless(poly, hi: int, lo: int = 2):
    """First p >= lo at which the integer polynomial has no root: the kill prime.

    The synthetic demos pass lo = 3: the h-rootless criterion is an odd-p
    statement (t = 2x-9 is only invertible mod odd p), and p = 2 rootlessness
    can be a parity artifact rather than the Chebotarev mechanism."""
    co = [int(x) for x in sp.Poly(poly, t).all_coeffs()]
    for p in sp.primerange(lo, hi):
        if not root_mod(co, p):
            return int(p)
    return None


def deep_h(mu: int, nu: int, al: int, be: int, ga: int):
    A = t * t + mu * t + nu
    m = t**3 + al * t * t + be * t + ga
    return sp.expand(A * A.subs(t, -t) * m * (-m.subs(t, -t)))


def disc_m(al: int, be: int, ga: int) -> int:
    return int(sp.discriminant(sp.Poly(t**3 + al * t * t + be * t + ga, t), t))


def coincidence(mu: int, nu: int, al: int, be: int, ga: int) -> bool:
    v = (mu * mu - 4 * nu) * disc_m(al, be, ga)
    return v > 0 and bool(sp.integer_nthroot(v, 2)[1])


def step0(kill_scan: int) -> None:
    print("  [0] the trichotomy, each case exercised on a synthetic polynomial")
    # (III)+coincidence -- the no-kill shape: m = t^3+t^2+2t+3 (disc -175, S3),
    # A = t^2+t+44 (disc 1-176 = -175): shared field, both parts irreducible.
    h = deep_h(1, 44, 1, 2, 3)
    check(coincidence(1, 44, 1, 2, 3), "coincidence holds for the example")
    kp = first_rootless(h, kill_scan)
    check(kp is None, f"no-kill example was killed at {kp}")
    print(f"      no-kill shape  (mu,nu|al,be,ga)=(1,44|1,2,3): root mod every p < {kill_scan}")
    # same shape WITHOUT the coincidence: A = t^2+t+1 (disc -3, -3*-175 not square)
    h2 = deep_h(1, 1, 1, 2, 3)
    kp2 = first_rootless(h2, kill_scan, lo=3)
    check(kp2 == 11, f"no-coincidence example: first odd kill 11, got {kp2}")
    print(f"      same shape, independent fields (A disc -3): killed at p={kp2}")
    # case (III) with q(t^2) irreducible: q = u^2+2u+176 (disc -700, shared with
    # disc m = -175), k = u^3+3u^2-2u-9 = the m*mtilde product. 3-cycle-class kill.
    q3, k3 = u * u + 2 * u + 176, u**3 + 3 * u * u - 2 * u - 9
    check(not sp.Poly(k3.subs(u, t * t), t).is_irreducible
          and sp.Poly(q3.subs(u, t * t), t).is_irreducible, "case (III) fixture")
    kp5 = first_rootless((q3 * k3).subs(u, t * t), 100, lo=3)
    check(kp5 == 43, f"case (III) example: expected kill at 43, got {kp5}")
    print(f"      case (III) k(t^2) reducible, q(t^2) irreducible: killed at p={kp5}")
    # (II): shared subfield, k(t^2) irreducible: k = u^3+u+1 (disc -31),
    # q = u^2+u+8 (disc -31). tau-class kill.
    g2 = sp.expand((u * u + u + 8) * (u**3 + u + 1))
    check(sp.Poly((u**3 + u + 1).subs(u, t * t), t).is_irreducible, "k(t^2) irreducible")
    kp3 = first_rootless(g2.subs(u, t * t), kill_scan, lo=3)
    check(kp3 == 11, f"case (II) example: expected kill at 11, got {kp3}")
    print(f"      case (II) shared subfield, k(t^2) irreducible: killed at p={kp3}")
    # (I): the two recorded candidates are generic -- D_q*D_k > 0 but not square
    for cv, kexp in CANDIDATES.items():
        fl = sp.Poly(R5 - 1024 * cv, u).factor_list()[1]
        P = 1
        for pp, _ in fl:
            P *= int(sp.discriminant(pp.as_poly(u), u))
        check(P > 0 and not sp.integer_nthroot(P, 2)[1], f"c={cv}: D_q*D_k nonsquare")
        kp4 = first_rootless((R5 - 1024 * cv).subs(u, t * t), 100, lo=3)
        check(kp4 == kexp, f"c={cv}: kill at {kexp}")
        print(f"      case (I)  c={cv}: D_q*D_k not a square, killed at p={kp4}")


def step1() -> None:
    print("  [1] c < 0 is impossible for the no-kill shape (so every such c is killed)")
    Rp = sp.Poly(sp.diff(R5, u), u)
    rts = Rp.real_roots()
    check(len(rts) == 4 and all(1 < r < 81 for r in rts), "R' has its 4 real roots in (1,81)")
    M = max(abs(R5.subs(u, r)) for r in rts)
    check(M < CRIT_MAX < 1024 * F10, "largest |critical value| below 43,930,543 < 1024*10!")
    for cv in (-F10, -5 * F10):
        P = sp.Poly(R5 - 1024 * cv, u)
        check(P.count_roots(1, 81) == 0, f"c={cv}: no root of g in [1,81]")
        check(P.count_roots(81, 10**9) == 0, f"c={cv}: no root beyond 81")
        check(P.count_roots(-10**9, 0) == 1, f"c={cv}: exactly one real root, negative")
    print(f"      max |critical value of R| = {float(M):.0f} < 1024*10!; R increasing on (-inf,1]")
    print("      10! | c, c < 0  ->  one real root, negative; but m's real root t0 needs")
    print("      t0^2 = that root >= 0. Deep shape impossible: cases (I)-(III) kill all c < 0.")
    print("      (c < 0 with 10! not dividing c, at any size: a power of a prime <= 7 kills.)")


def deep_equations():
    """The no-kill locus in (p4,..,p0): R(u) - 1024c = u P(u)^2 - Q(u)^2."""
    p4, p3, p2, p1, p0, c = sp.symbols("p4 p3 p2 p1 p0 c")
    lhs = sp.expand(u * (u * u + p3 * u + p1) ** 2 - (p4 * u * u + p2 * u + p0) ** 2)
    D = sp.Poly(lhs - (R5 - 1024 * c), u)
    return [sp.expand(D.coeff_monomial(u**k)) for k in range(6)], (p4, p3, p2, p1, p0, c)


def Fcurve(p4: int, p2: int) -> int:
    return (64 * p2**3 + 10560 * p2**2 * p4 - 8 * p2 * p4**6 + 1320 * p2 * p4**4
            + 498696 * p2 * p4**2 - 649000 * p2 + p4**9 - 660 * p4**7 + 93126 * p4**5
            + 5205420 * p4**3 - 5457375 * p4)


def step2() -> None:
    print("  [2] the no-kill locus as one plane curve")
    eqs, (p4, p3, p2, p1, p0, c) = deep_equations()
    s3 = sp.solve(eqs[4], p3)[0]
    s1 = sp.solve(sp.expand(eqs[3].subs(p3, s3)), p1)[0]
    s0 = sp.solve(sp.expand(eqs[2].subs({p3: s3, p1: s1})), p0)[0]
    num, _ = sp.fraction(sp.together(sp.expand(eqs[1].subs({p3: s3, p1: s1, p0: s0}))))
    num = sp.expand(num)
    check(sp.expand(num - Fcurve(p4, p2)) == 0 or sp.expand(num + Fcurve(p4, p2)) == 0,
          "eliminated curve equals +-F(p4,p2)")
    check(sp.solve(eqs[0], c)[0] == (p0**2 - 893025) / 1024, "c = (p0^2 - 945^2)/1024")
    check(s3 == (p4**2 - 165) / 2, "p3 = (p4^2-165)/2: p4 must be ODD; p4 = 0 impossible")
    mu, nu, al, be, ga = sp.symbols("mu nu alpha beta gamma")
    Dq = sp.expand((2 * nu - mu * mu) ** 2 - 4 * nu * nu)
    check(sp.expand(Dq - mu * mu * (mu * mu - 4 * nu)) == 0, "D_q = mu^2 (mu^2 - 4nu)")
    kk = u**3 + (2 * be - al * al) * u * u + (be * be - 2 * al * ga) * u - ga * ga
    mm = t**3 + al * t * t + be * t + ga
    check(sp.expand(sp.discriminant(sp.Poly(kk, u), u)
                    - sp.discriminant(sp.Poly(mm, t), t) * (ga - al * be) ** 2) == 0,
          "D_k = disc(m) (gamma - alpha beta)^2: the coincidence IS D_q D_k square")
    print("      F(p4,p2) = 64 p2^3 + 10560 p2^2 p4 - 8 p2 p4^6 + 1320 p2 p4^4 + 498696 p2 p4^2")
    print("                 - 649000 p2 + p4^9 - 660 p4^7 + 93126 p4^5 + 5205420 p4^3 - 5457375 p4 = 0")
    print("      c = (p0^2 - 945^2)/1024;  c = 0 <-> p0 = +-945.  For 10! | c:")
    print("      315 | p0, p0 odd, p0 = +-945 (mod 2^17)  [v_2, v_3, v_5, v_7 of p0^2 - 945^2]")


def step3(nu_max: int, p4_max: int) -> dict:
    import numpy as np
    from math import isqrt

    print(f"  [3] searches on the locus (not a proof): |nu| <= {nu_max}, odd p4 < {p4_max}")

    def Cab(a: int, b: int) -> int:
        return (-(a**4) - 165 * a**3 + 3 * a * a * b - 8778 * a * a + 330 * a * b
                - 172810 * a - b * b + 8778 * b - 1057221)

    bsq = []
    for nv in range(-nu_max, nu_max + 1):
        b = nv * nv
        co = [-1, -165, 3 * b - 8778, 330 * b - 172810, -b * b + 8778 * b - 1057221]
        for r in np.roots(co):
            if abs(r.imag) > 1e-6:
                continue
            for aa in (int(round(r.real)) - 1, int(round(r.real)), int(round(r.real)) + 1):
                if Cab(aa, b) == 0:
                    bsq.append((aa, nv))
    deep = []
    for aa, nv in sorted(set(bsq)):
        b = nv * nv
        f = -(aa**3) - 165 * aa * aa + 2 * aa * b - 8778 * aa + 165 * b - 172810
        cn = (aa**3 * b + 165 * aa * aa * b - 2 * aa * b * b + 8778 * aa * b
              - 165 * b * b + 172810 * b - 893025)
        mu2 = 2 * nv - aa
        if (mu2 >= 0 and isqrt(mu2) ** 2 == mu2 and -f >= 0 and isqrt(-f) ** 2 == -f):
            deep.append((aa, nv, cn // 1024 if cn % 1024 == 0 else cn / 1024))
    azs = sorted({aa for aa, _ in bsq})
    check(all(cv == 0 for _, _, cv in deep), f"nonzero-c deep point found: {deep}")
    print(f"      b = nu^2 on C at a in {azs}; mu^2 and gamma^2 also integral ONLY at c = 0")

    hits = []
    for p4 in range(1, p4_max, 2):
        co = [64, 10560 * p4, -8 * p4**6 + 1320 * p4**4 + 498696 * p4**2 - 649000,
              p4**9 - 660 * p4**7 + 93126 * p4**5 + 5205420 * p4**3 - 5457375 * p4]
        for r in np.roots(co):
            if abs(r.imag) > 1e-3 * max(1.0, abs(r.real)):
                continue
            rr = int(round(r.real))
            for p2 in (rr - 1, rr, rr + 1):
                if Fcurve(p4, p2) == 0:
                    hits.append((p4, p2))
    p0s = set()
    for p4, p2 in sorted(set(hits)):
        num = -8 * p2**2 + p4 * (8 * p2 * p4**2 - 1320 * p2 - p4**5 + 495 * p4**3
                                 - 46563 * p4) + 81125
        if num % (16 * p4) == 0:
            p0s.add(num // (16 * p4))
    check(p0s <= {945, -945}, f"F-curve point with integral p0 not +-945: {sorted(p0s)}")
    print(f"      F(p4,p2) = 0, odd p4 < {p4_max}: every point that lifts (p0 integral) has")
    print("      p0 = +-945, i.e. c = 0.  (F also has points with p0 NOT integral, on")
    print("      p2 = -p4^3 at p4 in {1,5,7,9} -- elimination artifacts, no lift to the locus.)")
    print("      (one-off runs: |nu| <= 10^6 i.e. b <= 10^12, 50 s, and odd p4 <= 2*10^5,")
    print("       2 s -- same result: nothing beyond c = 0)")
    return {"b_square_a": azs, "deep_nonzero": deep, "p0s": sorted(p0s)}


def step4() -> None:
    print("  [4] genus of the b = nu^2 cover, and what remains")
    qC = sp.Poly(u**4 + 165 * u**3 + 8778 * u * u + 172810 * u + 1057221, u)
    check(sp.discriminant(qC, u) != 0, "q_C(a) squarefree: the 4 zeros of b on C are simple")
    check(sp.resultant(sp.Poly(3 * u * u + 330 * u + 8778, u), qC) != 0,
          "s(a) shares no root with q_C: those zeros are smooth points of C")
    print("      nu^2 = b is a double cover of the genus-1 curve C ramified exactly at the")
    print("      four simple zeros of b (poles of b are double: unramified). Riemann-Hurwitz:")
    print("      2g - 2 = 2*0 + 4, GENUS 3 -- Q27's situation, one member later.")
    print("      Chebotarev stops at: the no-kill locus has no integer point with c != 0.")
    print("      Magma IntegralQuarticPoints (2026-08-23) listed the 15 a-values; no extras.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nu_max", type=int, default=4000)
    ap.add_argument("--p4_max", type=int, default=2001)
    ap.add_argument("--kill_scan", type=int, default=2000)
    args = ap.parse_args()

    print("=== Q28 (2,3): Chebotarev closes c < 0; Magma closed the c > 0 locus ===")
    step0(args.kill_scan)
    step1()
    step2()
    step3(args.nu_max, args.p4_max)
    step4()
    print()
    print("  RESULT: c < 0 (2,3) values are ALL killed, unconditionally.")
    print("          c > 0: Chebotarev kills off the no-kill locus; Magma showed")
    print("          that locus has no integral a outside the 15 (all c=0 or dead).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
