#!/usr/bin/env python3
"""Q29 Branch A CLOSED by Runge: the chord quartic's integer points, provably.

Branch A of k=9 -- f_c has a quadratic factor y^2+ay+b -- reduces (matching,
re-derived in step [0]) to the chord curve, a plane QUARTIC in A = a^2:

  Phi9(A,b) = -A^4 + 7A^3 b + 30A^3 - 15A^2 b^2 - 150A^2 b - 273A^2
              + 10A b^3 + 180A b^2 + 819A b + 820A
              - b^4 - 30b^3 - 273b^2 - 820b - 576 = 0,

smooth, genus 3 (NOT k=7's genus-1 cubic: the outline's warning is real).
Its leading form SPLITS over Q:

  L = -(A - b) * (A^3 - 6A^2 b + 9A b^2 - b^3),

a linear factor times the cyclic cubic form of discriminant 81 -- the field
Q(2cos(2pi/9)), k = 9 announcing itself the way Q(2cos(2pi/7)) did in Q26
Branch B.  Two coprime factors mean RUNGE'S METHOD applies, and it closes
the branch unconditionally:

  * the four slopes at infinity are s = 1 (rational) and the three conjugate
    roots beta of s^3 - 9s^2 + 6s - 1 (all real);
  * along the rational channel, b = A - 10 + 9/A + 90/A^2 + ... , so the
    INTEGER m = b - A satisfies m + 10 in (0, 1) for A >= A_FAR: no integer
    value exists.  (Even the trap value m = -10 is empty for ALL A:
    Phi9(A, A-10) = -27(A^2 - 10A - 12), discriminant 148, not a square.)
  * along the cyclic channels, the degree-4 integer polynomial

      W1 = A^3 b - 6A^2 b^2 - 20A^2 b - 9A^2 + 9A b^3 + 70A b^2 + 118A b
           + 90A - b^4 - 20b^3 - 82b^2 - 180b

    has NO positive powers of A in its branch expansion (verified
    symbolically on every run); its limits are omega = -270beta^2 + 2241beta
    - 594, numerically 18.827, 311.712, -573.538 (trace exactly -243 =
    -3^5), with distances 0.173, 0.288, 0.462 to the nearest integers --
    while W1 takes INTEGER values at integer points and is certified within
    err << 0.17 of a limit once A >= A_FAR.  No integer point escapes.
  * a cell-wise domination argument ([F]) shows every real point with
    A >= A_FAR lies in one of the four channel tubes, so the two traps
    cover everything; a CRT modular sweep ([S]) settles 0 <= A < A_FAR.

RESULT (the theorem): the complete list of integer points of Phi9 with
A >= 0 is the sweep list; those with A = a^2 a perfect square give exactly

    c = 0,  c = +-2630880,  c = +-176774400,

i.e. the two known Branch-A candidates and nothing else.  Both die
(p = 13 and p = 7, re-verified here).  PROVED_BRANCH_A_K9 = True.

Everything is exact: rational interval arithmetic end to end, algebraic
numbers enter only through certified rational enclosures.  Modelled on
scripts/k7_runge.py (Q26 Branch B), which pioneered the pattern.

    python scripts/k9_branchA_runge.py           full certificate + sweep
    python scripts/k9_branchA_runge.py --quick   certificate only, no sweep
"""

from __future__ import annotations

import argparse
import itertools
import json
import time
from fractions import Fraction as Q
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

A, b, beta, t, X, y, a_, c_ = sp.symbols("A b beta t X y a c")

PHI9 = (-(A**4) + 7*A**3*b + 30*A**3 - 15*A**2*b**2 - 150*A**2*b
        - 273*A**2 + 10*A*b**3 + 180*A*b**2 + 819*A*b
        + 820*A - b**4 - 30*b**3 - 273*b**2 - 820*b - 576)
K3F = A**3 - 6*A**2*b + 9*A*b**2 - b**3       # cyclic cubic form, disc 81
MP = beta**3 - 9*beta**2 + 6*beta - 1          # slopes of the K3F channels
W1 = (A**3*b - 6*A**2*b**2 - 20*A**2*b - 9*A**2 + 9*A*b**3 + 70*A*b**2
      + 118*A*b + 90*A - b**4 - 20*b**3 - 82*b**2 - 180*b)
OMEGA = -270*beta**2 + 2241*beta - 594         # W1's limit along each channel

SLOPE_BOUND = 16      # [L]: no real point with |b| > max(16A, 1000)
B_FLOOR = 1000
DELTA = [Q(7, 100), Q(7, 100), Q(1, 4), Q(1)]  # tube radii in s = b/A units
CELL = Q(1, 32)       # cell width for the [F] domination sweep
SWEEP_MODS = (128, 243, 125, 343, 121)         # CRT moduli for [S]

PROVED_BRANCH_A_K9 = True
BRANCH_A_C_VALUES = (0, 2630880, 176774400)
KNOWN_KILLS = {2630880: 13, 176774400: 7}


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def red(e):
    return sp.rem(sp.Poly(sp.expand(e), beta), sp.Poly(MP, beta)).as_expr()


def kinv(al):
    inv = sp.invert(sp.Poly(sp.expand(al), beta, domain="QQ"),
                    sp.Poly(MP, beta, domain="QQ"))
    return red(inv.as_expr())


# --------------------- exact interval arithmetic over Q ----------------------

def iv(x):
    x = Q(x)
    return (x, x)


def iadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def imul(u, v):
    ps = (u[0]*v[0], u[0]*v[1], u[1]*v[0], u[1]*v[1])
    return (min(ps), max(ps))


def ipow(u, n):
    r = iv(1)
    for _ in range(n):
        r = imul(r, u)
    return r


def iabsmax(u):
    return max(abs(u[0]), abs(u[1]))


def iabsmin(u):
    if u[0] <= 0 <= u[1]:
        return Q(0)
    return min(abs(u[0]), abs(u[1]))


def enclose_roots(bits: int = 100) -> list:
    """certified rational enclosures of the three real roots of MP."""
    f = sp.Poly(MP, beta)
    out = []
    for lo, hi in ((Q(1, 4), Q(1, 3)), (Q(2, 5), Q(1, 2)), (Q(8), Q(9))):
        flo = sp.sign(f.as_expr().subs(beta, sp.Rational(lo)))
        fhi = sp.sign(f.as_expr().subs(beta, sp.Rational(hi)))
        check(flo * fhi < 0, "root bracketed")
        for _ in range(bits):
            mid = (lo + hi) / 2
            fm = sp.sign(f.as_expr().subs(beta, sp.Rational(mid)))
            if fm == 0:
                lo, hi = mid - Q(1, 10**30), mid + Q(1, 10**30)
                break
            if fm == flo:
                lo = mid
            else:
                hi = mid
        out.append((lo, hi))
    check(out[0][1] < out[1][0] < out[2][0], "roots ordered and separated")
    return out


def eval_K(expr, om):
    """interval value of a Q(beta)-element on the enclosure om."""
    p = sp.Poly(red(expr), beta)
    r = iv(0)
    for (k,), co in p.terms():
        r = iadd(r, imul(iv(Q(int(sp.numer(co)), int(sp.denom(co)))), ipow(om, k)))
    return r


def tbound(expr, om, T):
    """upper bound for |expr(t)| over |t| <= T, K-coefficients on om."""
    p = sp.Poly(sp.expand(expr), t)
    tot = Q(0)
    for (k,), co in p.terms():
        tot += iabsmax(eval_K(co, om)) * Q(T) ** k
    return tot


# ------------------------------ the certificate ------------------------------

def certificate(verbose: bool = True) -> dict:
    log = print if verbose else (lambda *aa, **kk: None)

    # ---- [0] Phi9 re-derived from the matching; leading form split; W1 check
    P9 = y**9 - 30*y**7 + 273*y**5 - 820*y**3 + 576*y
    check(sp.expand(sp.prod([x_ - i for i in range(9)]).subs(x_, y + 4) - P9) == 0
          if (x_ := sp.Symbol("x_")) else False, "(x)_9 centred")
    es = sp.symbols("e0:7")
    sept = y**7 + sum(es[i]*y**i for i in range(7))
    D = sp.Poly(sp.expand((y*y + a_*y + b)*sept - (P9 - c_)), y)
    sol = sp.solve([D.coeff_monomial(y**k) for k in range(2, 9)], list(es), dict=True)[0]
    rem1 = sp.expand(D.coeff_monomial(y).subs(sol))
    check(c_ not in rem1.free_symbols, "y^1 condition free of c")
    phi_re = sp.expand(sp.Poly(rem1, a_).as_expr().subs(
        {a_**8: A**4, a_**6: A**3, a_**4: A**2, a_**2: A}))
    check(sp.expand(phi_re - PHI9) == 0, "Phi9 re-derived from the matching")
    cfun = sp.expand(-sp.solve(D.coeff_monomial(1).subs(sol), c_)[0])
    log("  [0]  Phi9 re-derived; c = -e0*b determined by (a, b)")
    check(sp.expand(PHI9 + (A - b)*K3F
                    - (30*A**3 - 150*A**2*b - 273*A**2 + 180*A*b**2 + 819*A*b
                       + 820*A - 30*b**3 - 273*b**2 - 820*b - 576)) == 0,
          "leading form is exactly -(A-b)*K3F")
    lead1 = sp.Poly(K3F.subs({A: 1, b: beta}), beta)
    check(sp.expand(lead1.as_expr() + MP) == 0, "K3F(1, beta) = -MP")
    check(sp.Poly(MP, beta).is_irreducible and sp.discriminant(sp.Poly(MP, beta), beta) == 81,
          "cyclic cubic, disc 81: the field Q(2cos(2pi/9))")
    sing = sp.solve([PHI9, sp.diff(PHI9, A), sp.diff(PHI9, b)], [A, b], dict=True)
    check(sing == [], "Phi9 smooth: genus 3")
    G = sp.expand(PHI9 + (A - b)*K3F)

    # ---- branch expansions (exact, in Q(beta)); rational channel in Q
    ts = sp.symbols("s0:6")
    bser = beta/X + sum(ts[i]*X**i for i in range(6))
    Pser = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        PHI9.subs({A: 1/X, b: bser}) * X**4)))), X)
    vals, order_of = {}, {}
    for k in range(0, Pser.degree() + 1):
        co = red(Pser.coeff_monomial(X**k).subs(vals))
        if co == 0:
            continue
        free = [s for s in ts if s in co.free_symbols]
        check(bool(free), f"series inconsistency at X^{k}")
        s = free[0]
        lin = sp.Poly(co, s)
        check(lin.degree() == 1, "series step is linear")
        al, be = lin.all_coeffs()
        vals[s] = red(-red(be) * kinv(red(al)))
        order_of[s] = k
        if len(vals) == 6:
            break
    check([order_of[ts[i]] for i in range(6)] == [1, 2, 3, 4, 5, 6], "t_i at X^{i+1}")
    log("  [ser] cyclic tail t0..t5 solved exactly in Q(2cos(2pi/9))")

    rs = sp.symbols("r0:6")
    rser = 1/X + sum(rs[i]*X**i for i in range(6))
    Rser = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        PHI9.subs({A: 1/X, b: rser}) * X**4)))), X)
    rvals = {}
    for k in range(0, Rser.degree() + 1):
        co = sp.expand(Rser.coeff_monomial(X**k).subs(rvals))
        if co == 0:
            continue
        free = [s for s in rs if s in co.free_symbols]
        check(bool(free), f"rational series inconsistency at X^{k}")
        s = free[0]
        lin = sp.Poly(co, s)
        check(lin.degree() == 1, "rational series step is linear")
        al, be = lin.all_coeffs()
        rvals[s] = sp.Rational(sp.expand(-be/al))
        if len(rvals) == 6:
            break
    check(rvals[rs[0]] == -10 and rvals[rs[1]] == 9 and rvals[rs[2]] == 90,
          "rational channel: b = A - 10 + 9/A + 90/A^2 + ...")
    log("  [ser] rational tail: b = A - 10 + 9/A + 90/A^2 + 630/A^3 + ...")
    lin_line = sp.factor(sp.expand(PHI9.subs(b, A - 10)))
    check(sp.expand(lin_line + 27*(A**2 - 10*A - 12)) == 0
          and not sp.integer_nthroot(148, 2)[1],
          "Phi9(A, A-10) = -27(A^2-10A-12): no integer roots (disc 148 nonsquare)")

    # ---- W1: no positive powers along the cyclic channels; limits omega
    PWX = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        W1.subs({A: 1/X, b: bser.subs(vals)}) * X**4)))), X)
    for k in range(0, 4):
        check(red(PWX.coeff_monomial(X**k)) == 0, f"W1 bounded: X^{k} coefficient")
    w0 = red(PWX.coeff_monomial(X**4))
    check(sp.expand(w0 - OMEGA) == 0, "W1 limit is omega = -270b^2+2241b-594")
    log("  [W1] no positive powers of A along any cyclic branch; limit omega")

    oms = enclose_roots()
    slopes = [oms[0], oms[1], (Q(1), Q(1)), oms[2]]
    delta = [DELTA[0], DELTA[1], DELTA[2], DELTA[3]]

    # trace pin: sum of omega over the three roots is -243 = -3^5
    e1, e2 = Q(9), Q(6)                       # elementary symmetric of the roots
    p2 = e1*e1 - 2*e2                          # power sum beta^2
    check(-270*p2 + 2241*e1 - 3*594 == -243, "trace(omega) = -243 = -3^5")

    # ---- [L] no real point with |b| >= max(16A, 1000)  (A >= 0)
    P = sp.Poly(PHI9, A, b)
    bud = Q(0)
    for (p, q), co in P.terms():
        if (p, q) == (0, 4):
            continue
        d = p + q
        bud += abs(co) * Q(1, SLOPE_BOUND**p) * Q(1, B_FLOOR**(4 - d))
    check(bud < 1, f"[L] budget {float(bud):.4f} < 1")
    log(f"  [L]  no real point with |b| >= max({SLOPE_BOUND}A, {B_FLOOR})"
        f"   (budget {float(bud):.3f} < 1)")

    # ---- [F] cell-wise domination: A_FAR
    # outside the four tubes and with |s| <= 16, every cell certifies
    # prodmin * A^4 > Gam_cell(A) for A >= A_far_cell; A_FAR = max over cells.
    Gterms = sp.Poly(G, A, b).terms()

    def gam_cell(smax):
        g = {}
        for (p, q), co in Gterms:
            d = p + q
            g[d] = g.get(d, Q(0)) + abs(co) * Q(smax) ** q
        return g

    def prod_lo(s):
        tot = Q(1)
        for (lo, hi) in slopes:
            if lo <= s <= hi:
                return Q(0)
            tot *= min(abs(s - lo), abs(s - hi))
        return tot

    tube_iv = [(lo - dd, hi + dd) for (lo, hi), dd in zip(slopes, delta)]
    cells = []
    s_pts = [Q(-SLOPE_BOUND)]
    s_cur = Q(-SLOPE_BOUND)
    while s_cur < SLOPE_BOUND:
        s_next = min(s_cur + CELL, Q(SLOPE_BOUND))
        cells.append((s_cur, s_next))
        s_cur = s_next
    A_FAR = 0
    for (sl, sr) in cells:
        if any(sl < thi and sr > tlo for (tlo, thi) in tube_iv):
            # clip the cell to the part outside every tube
            segs = [(sl, sr)]
            for (tlo, thi) in tube_iv:
                segs = [seg for pair in
                        ([(max(x0, thi), x1)] if x0 >= tlo else
                         [(x0, min(x1, tlo)), (max(x0, thi), x1)]
                         for (x0, x1) in segs)
                        for seg in pair if seg[0] < seg[1]]
        else:
            segs = [(sl, sr)]
        for (x0, x1) in segs:
            pl = Q(1)
            for (lo, hi) in slopes:
                cands = []
                for xx in (x0, x1):
                    if lo <= xx <= hi:
                        cands.append(Q(0))
                    else:
                        cands.append(min(abs(xx - lo), abs(xx - hi)))
                if x0 <= lo and hi <= x1:
                    cands.append(Q(0))
                pl *= min(cands)
            check(pl > 0, f"cell [{float(x0)},{float(x1)}] outside tubes has prod > 0")
            g = gam_cell(max(abs(x0), abs(x1)))
            hi_ = 2
            while pl * Q(hi_)**4 <= sum(g.get(d_, Q(0)) * Q(hi_)**d_ for d_ in range(4)):
                hi_ *= 2
            lo_ = hi_ // 2
            while lo_ + 1 < hi_:
                mid = (lo_ + hi_) // 2
                if pl * Q(mid)**4 > sum(g.get(d_, Q(0)) * Q(mid)**d_ for d_ in range(4)):
                    hi_ = mid
                else:
                    lo_ = mid
            A_FAR = max(A_FAR, hi_)
    check(A_FAR < 100000, f"A_FAR = {A_FAR} manageable")
    for i in range(4):
        check(delta[i] * A_FAR > 1, "tube wide enough at A_FAR to matter")
    log(f"  [F]  A_FAR = {A_FAR}: for A >= A_FAR every real point lies in a tube")
    A_COVER = A_FAR   # coverage threshold; traps below may need a larger A_FAR


    def far_analysis(A_FAR):
        # ---- [T] per-channel tube radii r_i (|b - s_i A| <= r_i for A >= A_FAR)
    # Load-bearing subtlety (verified in adversarial review): the narrowing
    # iteration starts from r = delta*A_FAR, but for A > A_FAR a tube point
    # can have |t| up to delta*A > delta*A_FAR.  The Gb bound is nevertheless
    # valid for ALL A >= A_FAR because every monomial t^m A^k of G (and of
    # dW1/db) along a channel line has m + k <= 3, so each term scales by
    # (A/A_FAR)^(m+k-3) <= 1.  This holds automatically here since deg G = 3;
    # any adaptation of this pattern to a curve whose sub-leading part has
    # higher degree along the line must re-check it or become unsound.
        Gs = sp.Poly(sp.expand(G.subs(b, beta*A + t)), A)
        g_of = {k: red(Gs.coeff_monomial(A**k)) for k in range(4)}
        tubes = []
        for i in (0, 1, 3):
            om = oms[i if i < 2 else 2]
            r = Q(delta[i]) * A_FAR
            for _ in range(6):
                Gb = sum(tbound(g_of[k], om, r) * Q(1, A_FAR) ** (3 - k) for k in range(4))
                prod_lo_ = Q(1)
                for j in range(4):
                    if j == i:
                        continue
                    (lo, hi) = slopes[j]
                    d_ = iabsmin((om[0] - hi, om[1] - lo))
                    d_ = d_ - Q(r, A_FAR)
                    check(d_ > 0, "tube inside separation")
                    prod_lo_ *= d_
                r_new = Gb / prod_lo_
                if r_new >= r:
                    break
                r = r_new
            r = r + Q(1, 2)
            check(r < delta[i] * A_FAR, "tube radius below the [F] coverage radius")
            tubes.append((i, r))
        # rational channel tube
        Gr = sp.Poly(sp.expand(G.subs(b, A + t)), A)
        gr_of = {k: sp.expand(Gr.coeff_monomial(A**k)) for k in range(4)}
        rrat = Q(delta[2]) * A_FAR
        om1 = (Q(1), Q(1))
        for _ in range(6):
            Gb = sum(tbound(gr_of[k], om1, rrat) * Q(1, A_FAR) ** (3 - k) for k in range(4))
            prod_lo_ = Q(1)
            for j in (0, 1, 3):
                (lo, hi) = slopes[j]
                d_ = min(abs(1 - lo), abs(1 - hi)) - Q(rrat, A_FAR)
                check(d_ > 0, "rational tube inside separation")
                prod_lo_ *= d_
            r_new = Gb / prod_lo_
            if r_new >= rrat:
                break
            rrat = r_new
        rrat = rrat + Q(1, 2)
        check(rrat < delta[2] * A_FAR, "rational tube consistent")
        log(f"  [T]  tube radii: cyclic {[float(r) for _, r in tubes]},"
            f" rational {float(rrat):.2f}")

        # ---- [N]+[W] cyclic channels: Newton residual and the W1 trap
        hsc = sp.Poly(sp.expand(PHI9.subs(b, beta*A + t)), A)
        hs = {k: red(hsc.coeff_monomial(A**k)) for k in range(5)}
        check(hs[4] == 0 or red(hs[4]) == 0, "leading form vanishes along cyclic line")
        dh = {k: sp.diff(hs[k], t) for k in range(4)}
        tail6 = [vals[ts[i]] for i in range(6)]
        results = []
        for idx, (i, r) in enumerate(tubes):
            om = oms[idx if idx < 2 else 2]
            # d/db Phi9 has A^3 coefficient with slope bounded below on the tube
            F_lo = iabsmin(eval_K(sp.Poly(dh[3], t).coeff_monomial(1), om))
            check(F_lo > 0, "dPhi9/db A^3 term bounded away from 0")
            dlo = F_lo - (tbound(dh[2], om, r)/A_FAR + tbound(dh[1], om, r)/Q(A_FAR)**2
                          + tbound(dh[0], om, r)/Q(A_FAR)**3
                          + tbound(sp.diff(sp.expand(hsc.coeff_monomial(A**4)
                                                     if hsc.degree() >= 4 else 0), t), om, r))
            check(dlo > 0, "dPhi9/db dominated by its A^3 term on the tube")
            tail_expr = sum(tail6[k] * X**k for k in range(6))
            Pres = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
                PHI9.subs({A: 1/X, b: beta/X + tail_expr}) * X**4)))), X)
            rb = Q(0)
            for k in range(Pres.degree() + 1):
                co = red(Pres.coeff_monomial(X**k))
                if co == 0:
                    continue
                check(k >= 7, f"residual order (X^{k})")
                rb += iabsmax(eval_K(co, om)) * Q(1, A_FAR) ** (k - 7)
            rho = rb / dlo          # |b* - (betaA + tail(1/A))| <= rho / A^6
            tail_mag = sum(iabsmax(eval_K(tail6[k], om)) * Q(1, A_FAR)**k for k in range(6))
            check(tail_mag + rho/Q(A_FAR)**6 <= r, "tail point and root inside the tube")
            # W1 trap: |W1 - omega| <= wtail/A + |dW1/db| * rho/A^6
            wtail = Q(0)
            for k in range(5, PWX.degree() + 1):
                co = red(PWX.coeff_monomial(X**k))
                if co != 0:
                    wtail += iabsmax(eval_K(co, om)) * Q(1, A_FAR) ** (k - 5)
            dW = sp.Poly(sp.expand(sp.diff(W1, b).subs(b, beta*A + t)), A)
            DWrho = sum(tbound(red(dW.coeff_monomial(A**k)), om, r) * rho / Q(A_FAR) ** (6 - k)
                        for k in range(4))
            err = wtail / A_FAR + DWrho
            w0i = eval_K(OMEGA, om)
            import math as _m
            fl = _m.floor(w0i[0])
            check(_m.floor(w0i[1]) == fl, "omega enclosure contains no integer")
            dist = min(w0i[0] - fl, Q(fl + 1) - w0i[1])
            check(err < dist, f"cyclic channel {idx}: trap err {float(err):.2e} < "
                              f"dist(omega, Z) {float(dist):.3f}")
            results.append((float(w0i[0]), float(err), float(dist)))
            log(f"  [N/W] cyclic {idx}: omega ~ {float(w0i[0]):+.3f}, err {float(err):.1e}"
                f" < dist {float(dist):.3f} -> NO integer W1 value, NO points")

        # ---- [M] rational channel: m = b - A has m + 10 in (0, 1)
        rtail6 = [rvals[rs[i]] for i in range(6)]
        dr_lo = None
        # Newton residual for the rational branch
        hr = sp.Poly(sp.expand(PHI9.subs(b, A + t)), A)
        hrs = {k: sp.expand(hr.coeff_monomial(A**k)) for k in range(5)}
        check(sp.expand(hrs[4]) == 0, "leading form vanishes along s = 1")
        dhr = {k: sp.diff(hrs[k], t) for k in range(4)}
        F_lo = abs(Q(sp.Rational(sp.Poly(dhr[3], t).coeff_monomial(1))))
        dlo = F_lo - sum(tbound(dhr[k], om1, rrat) / Q(A_FAR) ** (3 - k) for k in range(3))
        check(dlo > 0, "dPhi9/db dominated on the rational tube")
        rtail_expr = sum(rtail6[k] * X**k for k in range(6))
        Rres = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
            PHI9.subs({A: 1/X, b: 1/X + rtail_expr}) * X**4)))), X)
        rb = Q(0)
        for k in range(Rres.degree() + 1):
            co = sp.expand(Rres.coeff_monomial(X**k))
            if co == 0:
                continue
            check(k >= 7, f"rational residual order (X^{k})")
            rb += abs(Q(sp.Rational(co))) * Q(1, A_FAR) ** (k - 7)
        rho_r = rb / dlo
        # m + 10 = 9/A + 90/A^2 + ... +- rho_r/A^6: certify strictly in (0, 1)
        tail_hi = sum(abs(Q(rvals[rs[i]])) * Q(1, A_FAR)**i for i in range(1, 6))
        m_hi = tail_hi + rho_r / Q(A_FAR)**5
        m_lo = Q(9, A_FAR) - sum(abs(Q(rvals[rs[i]])) * Q(1, A_FAR)**i
                                 for i in range(2, 6)) - rho_r / Q(A_FAR)**5
        check(m_lo > 0, f"rational channel: m + 10 >= {float(m_lo):.2e} > 0")
        check(m_hi < 1, f"rational channel: m + 10 <= {float(m_hi):.2e} < 1")
        log(f"  [M]  rational channel: m + 10 in ({float(m_lo):.1e}, {float(m_hi):.1e})"
            f" (0,1) -> m = b - A is never an integer for A >= A_FAR")

        return A_FAR, tubes, rrat, results

    A_FAR = A_COVER
    for _ in range(10):
        try:
            A_FAR, tubes, rrat, results = far_analysis(A_FAR)
            break
        except RuntimeError as ex:
            log(f'       (A_FAR = {A_FAR}: {ex} -- doubling)')
            A_FAR *= 2
    else:
        check(False, 'trap escalation did not converge')
    log(f'  [far] all traps certified at A_FAR = {A_FAR}')
    # tie the hard-coded c_of to the matching-derived cfun (c is odd in a:
    # cfun(a, b) = -c_of(a, b) is the same +-orbit)
    av, bv = sp.symbols("av bv")
    check(sp.expand(cfun.subs({a_: av, b: bv})
                    + c_of(av, bv)) == 0
          or sp.expand(cfun.subs({a_: av, b: bv}) - c_of(av, bv)) == 0,
          "c_of matches the matching-derived c (up to the a -> -a sign)")
    return {"A_FAR": A_FAR, "tubes": [float(r) for _, r in tubes],
            "rational_tube": float(rrat), "traps": results}


# --------------------------------- [S] sweep ---------------------------------

def sweep(A_max: int, verbose: bool = True) -> list:
    """every integer point of Phi9 with 0 <= A < A_max, |b| <= max(16A, 1000).
    CRT modular sieve, complete by construction: prod(SWEEP_MODS) exceeds the
    b-range, so each allowed residue combination is a unique candidate."""
    t0 = time.time()
    Pc = sp.Poly(PHI9, b)
    cpolys = [sp.Poly(cq, A) for cq in Pc.all_coeffs()]

    def coeffs_at(Av):
        return [int(cp.eval(Av)) for cp in cpolys]

    from math import prod
    M = prod(SWEEP_MODS)
    check(M > 2 * max(SLOPE_BOUND * A_max, B_FLOOR) + 1, "CRT range covers |b|")
    tables = {}
    for m_ in SWEEP_MODS:
        tab = []
        for rA in range(m_):
            cs = [cq % m_ for cq in coeffs_at(rA)]
            ok = tuple(r_ for r_ in range(m_)
                       if ((((cs[0]*r_ + cs[1]) * r_ + cs[2]) * r_ + cs[3]) * r_
                           + cs[4]) % m_ == 0)
            tab.append(ok)
        tables[m_] = tab
    basis = []
    for m_ in SWEEP_MODS:
        Mi = M // m_
        basis.append(Mi * pow(Mi, -1, m_))
    pts = []
    for Av in range(0, A_max):
        allowed = []
        empty = False
        for m_ in SWEEP_MODS:
            ok = tables[m_][Av % m_]
            if not ok:
                empty = True
                break
            allowed.append(ok)
        if empty:
            continue
        cs = None
        Bmax = max(SLOPE_BOUND * Av, B_FLOOR)
        for combo in itertools.product(*allowed):
            r_ = sum(cq * e_ for cq, e_ in zip(combo, basis)) % M
            cand = r_ if r_ <= Bmax else r_ - M
            if cand < -Bmax or cand > Bmax:
                continue
            if cs is None:
                cs = coeffs_at(Av)
            v = ((((cs[0]*cand + cs[1])*cand + cs[2])*cand + cs[3])*cand + cs[4])
            if v == 0:
                pts.append((Av, cand))
    pts = sorted(set(pts))
    if verbose:
        print(f"  [S]  sweep 0 <= A < {A_max}: {len(pts)} integer points"
              f"   ({time.time()-t0:.0f}s)")
    return pts


# ------------------------------- [C] the c-list -------------------------------

def c_of(av, bv):
    return (av**7*bv - 6*av**5*bv**2 - 30*av**5*bv + 10*av**3*bv**3
            + 120*av**3*bv**2 + 273*av**3*bv - 4*av*bv**4 - 90*av*bv**3
            - 546*av*bv**2 - 820*av*bv)


def c_list(pts, verbose: bool = True) -> dict:
    import math
    sq = []
    for (Av, bv) in pts:
        r = math.isqrt(Av)
        if r*r == Av:
            cv = c_of(r, bv)
            sq.append((Av, bv, r, cv))
    cs = sorted({abs(cv) for _, _, _, cv in sq})
    if verbose:
        print(f"  [C]  points with A = a^2: {len(sq)}; |c| values: {cs}")
    check(tuple(cs) == BRANCH_A_C_VALUES, f"the Branch-A c-list: {cs}")
    BASE = sp.prod([sp.Symbol('x') - i for i in range(9)])
    xs = sp.Symbol('x')
    for cv, kp in KNOWN_KILLS.items():
        facs = [f for f, _ in sp.Poly(BASE - cv, xs).factor_list()[1]]
        check(tuple(sorted(f.degree() for f in facs)) == (2, 7), f"c={cv} is (2,7)")
        for p_ in sp.primerange(2, kp):
            check(any(any(int(f.eval(v)) % p_ == 0 for v in range(p_)) for f in facs),
                  f"c={cv}: root mod {p_}")
        check(all(not any(int(f.eval(v)) % kp == 0 for v in range(kp)) for f in facs),
              f"c={cv}: killed at {kp}")
        if verbose:
            print(f"       c = +-{cv}: (2,7), dies at p = {kp}")
    return {"square_points": len(sq), "c_values": cs}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true", help="certificate only")
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q29 Branch A: CLOSED by Runge (leading form splits off the")
    print("    disc-81 cyclic cubic; four channels, two traps, one sweep) ===")
    cert = certificate()
    out = {"A_FAR": cert["A_FAR"], "tubes": cert["tubes"],
           "rational_tube": cert["rational_tube"], "traps": cert["traps"]}
    if not args.quick:
        pts = sweep(cert["A_FAR"])
        cl = c_list(pts)
        out["sweep_points"] = len(pts)
        out["square_points"] = cl["square_points"]
        out["c_values"] = [int(v) for v in cl["c_values"]]
        print()
        print("  RESULT: Branch A of Q29 is a THEOREM.  The complete list of")
        print("          integer points of Phi9 with A >= 0 is the sweep list;")
        print("          A = a^2 points give exactly c = 0, +-2630880 (dies at 13),")
        print("          +-176774400 (dies at 7).  No other 2+7 / 2+2+5 / 2+3+4 /")
        print("          2+2+2+3 value of (x)_9 - c exists.  PROVED, no Siegel.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "claim": "Q29 Branch A (quadratic factor) closed by Runge",
            "curve": "chord quartic Phi9(A,b), genus 3, leading form -(A-b)*K3(disc 81)",
            "W1_limits_trace": -243,
            **out,
        }, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
