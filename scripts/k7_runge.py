#!/usr/bin/env python3
"""Q26 Branch B closed: the 3+4 curve's integral points, complete, by Runge.

docs/q26-k7-intersective.md reduced "f_c = (x)_7 - c has a cubic factor and no
quadratic factor" to integer points (A, b) = (a^2, b) of the plane quintic
C34(A, b) = 0 (generic locus) plus a degenerate locus 3a^2 - 2b - 14 = 0 that
is FINITE by pure algebra (step [D]: a in {0, +-2}; it carries c = -+896).
The doc left the generic locus "complete only to |a| <= 400" -- Siegel
finiteness with no certificate -- and proposed PARI's `thue` to close it.

`thue` was never the right tool, but something better is true: RUNGE'S METHOD
applies, because the quintic's leading form factors over Q as

    F5 = Q2 * K3,  Q2 = 2A^2 - 5Ab + 4b^2   (positive definite, disc -7),
                   K3 = A^3 - 2A^2 b - Ab^2 + b^3,

and K3(1, w) = w^3 - w^2 - 2w + 1 is the cyclic cubic of discriminant 49 --
the field Q(2cos(2pi/7)); k = 7 announcing itself.  Q2 definite means the two
Q2-branches at infinity are complex, so every REAL point runs along one of
K3's three real branches b = om*A + t0 + t1/A + ...  The polynomial

    W2 = A^4 - 10A^3 - 5A^2 b^2 - 36A^2 b - 49A^2 - Ab^3 + 10Ab^2 + 119Ab
         + 238A + 2b^4 + 32b^3 + 175b^2 + 404b

is BOUNDED along all three branches (its expansion has no positive power of
A -- re-derived and verified symbolically below, never assumed).  Hence on
integer points with A large the integer m = W2(A, b) is trapped within 1/2
of one of three conjugate limits w0(om); each possible m cuts the quintic in
a finite Bezout set; below the explicit threshold a modular sieve enumerates
every integer point directly.  All inequalities are exact rational
arithmetic; omega enters only through certified rational enclosures.

The chain (A >= 0 suffices, since A = a^2):
  [L] real roots have |b| <= max(3A, 1000)         (coefficient budget < 4)
  [F] for A >= A_FAR every real root is within sep/2 of one branch line
  [T] per-branch tube radius sharpens to O(10) via the exact expansion of
      G = C34 - Q2*K3 along the line (its A^4-coefficient is t-free)
  [N] the Newton residual of the 5-term branch tail pins the root within rho/A^5
  [W] |W2 - w0(om)| < 1/2  =>  m = W2 lies in an explicit 3-element-ish set
  [E] Res_b(C34, W2 - m): no integer root A >= A_FAR for any candidate m
  [S] modular-sieve sweep of 0 <= A < A_FAR: every remaining point
  [D] the degenerate locus is exactly a in {0, +-2}, by algebra

RESULT: the complete list of integer points of C34 with A >= 0, hence the
complete Branch-B list: c = 0 and c = -+896.  Branch B is CLOSED,
unconditionally -- no Thue solver, no external CAS, no unproved bound.  (Branch A,
the quadratic-factor branch, is a smooth plane cubic -- genus 1, NOT a Thue
equation; see scripts/k7_branchA.py.)

    python scripts/k7_runge.py            (~3-4 min, single core)
    python scripts/k7_runge.py --quick    (certificate only, no sweep)
"""

from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import sympy as sp
from sympy import Rational as Q

ROOT = Path(__file__).resolve().parents[1]

A, b, w, t, X, aa, dd, cc = sp.symbols("A b w t X a d c")

C34 = (2*A**5 - 9*A**4*b - 70*A**4 + 12*A**3*b**2 + 210*A**3*b + 833*A**3 - A**2*b**3
       - 98*A**2*b**2 - 1225*A**2*b - 3792*A**2 - 9*A*b**4 - 168*A*b**3 - 882*A*b**2
       - 432*A*b + 4179*A + 4*b**5 + 112*b**4 + 1176*b**3 + 5632*b**2 + 11620*b + 7056)
Q2F = 2*A**2 - 5*A*b + 4*b**2
K3F = A**3 - 2*A**2*b - A*b**2 + b**3
W2 = (A**4 - 10*A**3 - 5*A**2*b**2 - 36*A**2*b - 49*A**2 - A*b**3 + 10*A*b**2
      + 119*A*b + 238*A + 2*b**4 + 32*b**3 + 175*b**2 + 404*b)
MINPOL = w**3 - w**2 - 2*w + 1
E2 = (aa**5 - 4*aa**3*b - 14*aa**3 + 3*aa**2*dd + 3*aa*b**2 + 28*aa*b + 49*aa
      - 2*b*dd - 14*dd)
E1 = (aa**4*b - aa**3*dd - 3*aa**2*b**2 - 14*aa**2*b + 4*aa*b*dd + 14*aa*dd
      + b**3 + 14*b**2 + 49*b - dd**2 + 36)
E0 = (aa**4*dd - 3*aa**2*b*dd - 14*aa**2*dd + 2*aa*dd**2 + b**2*dd + 14*b*dd + cc
      + 49*dd)
SWEEP_MODS = [32, 27, 25, 7, 11, 13]  # pairwise coprime, lcm 21,621,600


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def red(e):
    e = sp.expand(e)
    return sp.rem(e, MINPOL, w) if w in getattr(e, "free_symbols", set()) else e


def kinv(al):
    u_, v_, g_ = sp.gcdex(sp.Poly(al, w), sp.Poly(MINPOL, w))
    check(g_.degree() == 0, "field inverse exists")
    return red(u_.as_expr() / g_.as_expr())


# --------------------- exact interval arithmetic over Q ----------------------

def iv(x):
    x = Q(x)
    return (x, x)


def iadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def imul(u, v):
    p = (u[0]*v[0], u[0]*v[1], u[1]*v[0], u[1]*v[1])
    return (min(p), max(p))


def ipow(u, n):
    r = iv(1)
    for _ in range(n):
        r = imul(r, u)
    return r


def iabsmax(u):
    return max(abs(u[0]), abs(u[1]))


def iabsmin(u):
    return Q(0) if u[0] <= 0 <= u[1] else min(abs(u[0]), abs(u[1]))


def enclose_roots(bits: int = 100) -> list:
    """Certified rational enclosures of the three real roots of MINPOL."""
    f = sp.Poly(MINPOL, w)
    out = []
    for lo0, hi0 in ((-2, -1), (0, 1), (1, 2)):
        lo, hi = Q(lo0), Q(hi0)
        check(f.eval(lo) * f.eval(hi) < 0, "root bracketed")
        for _ in range(bits):
            mid = (lo + hi) / 2
            v = f.eval(mid)
            if v == 0:
                lo = hi = mid
                break
            if f.eval(lo) * v < 0:
                hi = mid
            else:
                lo = mid
        out.append((lo, hi))
    return out


def eval_K(expr, om):
    p = sp.Poly(red(expr), w)
    r = iv(0)
    for (k,), co in p.terms():
        r = iadd(r, imul(iv(co), ipow(om, k)))
    return r


def tbound(poly_t, om, T):
    """Upper bound for |poly(t)| over |t| <= T, K-coefficients on interval om."""
    r = Q(0)
    for (k,), co in sp.Poly(poly_t, t).terms():
        r += iabsmax(eval_K(co, om)) * Q(T) ** k
    return r


# ------------------------------ the certificate ------------------------------

def certificate(verbose: bool = True) -> dict:
    log = print if verbose else (lambda *a, **k: None)

    check(sp.expand(sp.Poly(C34, A, b).homogeneous_order() or 0) == 0 or True, "n/a")
    check(sp.expand((C34 - Q2F*K3F) - (C34 - sp.expand(Q2F*K3F))) == 0, "G defined")
    Gpoly = sp.expand(C34 - Q2F*K3F)
    check(sp.Poly(Gpoly, A, b).total_degree() == 4, "leading form is exactly Q2*K3")

    # expansions along b = w*A + t
    PE = sp.Poly(sp.expand(C34.subs(b, w*A + t)), A)
    hs = {k: red(PE.coeff_monomial(A**k)) for k in range(6)}
    check(hs[5] == 0, "leading form vanishes along the line")
    PW = sp.Poly(sp.expand(W2.subs(b, w*A + t)), A)
    check(red(PW.coeff_monomial(A**4)) == 0, "W2 leading form vanishes at (1, omega)")

    # branch tail t0..t4, re-derived from scratch
    ts = sp.symbols("s0:5")
    bser = w/X + sum(ts[i]*X**i for i in range(5))
    Pser = sp.Poly(sp.expand(C34.subs({A: 1/X, b: bser}) * X**5), X)
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
        if len(vals) == 5:
            break
    check(sp.simplify(red(vals[ts[0]] + 2*(w + 2))) == 0, "t0 = -2(w+2)")
    check([order_of[ts[i]] for i in range(5)] == [1, 2, 3, 4, 5], "t_i fixed at X^{i+1}")

    # W2 bounded along the branch, and its limit w0 (exact, in K)
    PWX = sp.Poly(sp.expand(W2.subs({A: 1/X, b: bser.subs(vals)}) * X**4), X)
    for k in range(0, 4):
        check(red(PWX.coeff_monomial(X**k)) == 0, f"W2 bounded: X^{k} coefficient")
    w0 = red(PWX.coeff_monomial(X**4))
    log("  [W2] no positive powers of A along any branch; limit w0 =", w0)

    oms = enclose_roots()
    sep = {}
    for i in range(3):
        for j in range(3):
            if i != j:
                sep[i, j] = iabsmin((oms[i][0]-oms[j][1], oms[i][1]-oms[j][0]))
    sep_min = min(sep.values())
    half = sep_min / 2

    # ---- [L]
    P = sp.Poly(C34, A, b)
    bud5 = sum(abs(co)*Q(1, 3**p) for (p, q), co in P.terms()
               if p + q == 5 and (p, q) != (0, 5))
    bud4 = sum(abs(co)*Q(1, 3**p) for (p, q), co in P.terms() if p + q == 4)
    budlow = sum(abs(co)*Q(1, 3**p) for (p, q), co in P.terms() if p + q <= 3)
    budget = bud5 + Q(bud4, 1000) + Q(budlow, 1000**2)
    check(budget < 4, f"[L] budget {float(budget)} < 4")  # ~3.87
    log(f"  [L]  no real root with |b| >= max(3A, 1000)   (budget {float(budget):.3f} < 4)")

    # ---- [F]
    gam = {}
    for (p, q), co in sp.Poly(Gpoly, A, b).terms():
        gam[p + q] = gam.get(p + q, Q(0)) + abs(co) * Q(3) ** q

    def Gam(av):
        return sum(gam.get(d_, 0) * Q(av) ** d_ for d_ in range(5))

    hi_ = 2
    while Q(7, 16) * half**3 * Q(hi_)**5 <= Gam(hi_):
        hi_ *= 2
    lo_ = hi_ // 2
    while lo_ + 1 < hi_:
        mid = (lo_ + hi_) // 2
        if Q(7, 16) * half**3 * Q(mid)**5 > Gam(mid):
            hi_ = mid
        else:
            lo_ = mid
    A_FAR = hi_
    log(f"  [F]  A_FAR = {A_FAR}")

    # ---- [T] per-branch tubes (all bounds monotone: valid for every A >= A_FAR)
    Gs = sp.Poly(sp.expand(Gpoly.subs(b, w*A + t)), A)
    g_of = {k: red(Gs.coeff_monomial(A**k)) for k in range(5)}
    check(sp.Poly(g_of[4], t).degree() <= 0, "G's A^4 coefficient is t-free")
    Q2s = sp.Poly(sp.expand(Q2F.subs(b, w*A + t)), A)
    q_of = {k: red(Q2s.coeff_monomial(A**k)) for k in range(3)}
    tubes = []
    for i in range(3):
        om = oms[i]
        r = Q(16, 7) * Gam(A_FAR) / (Q(A_FAR)**2 * ((sep_min - half) * A_FAR)**2)
        for _ in range(4):
            Gb = sum(tbound(g_of[k], om, r) * Q(1, A_FAR) ** (4 - k) for k in range(5))
            q2lo = (iabsmin(eval_K(q_of[2], om)) - tbound(q_of[1], om, r)/A_FAR
                    - tbound(q_of[0], om, r)/Q(A_FAR)**2)
            # Q2 >= (7/16)A^2 globally (complete the square in b), so:
            q2lo = max(q2lo, Q(7, 16))
            prod_lo = Q(1)
            for j in range(3):
                if j != i:
                    pl = sep[i, j] - Q(r, A_FAR)
                    check(pl > 0, "tube inside the separation")
                    prod_lo *= pl
            r_new = Gb / (q2lo * prod_lo)
            if r_new >= r:
                break
            r = r_new
        r = r + Q(1, 2)  # slack: any bound above a valid one is valid
        check(r < half * A_FAR, "tube consistent with [F]")
        tubes.append(r)
        log(f"  [T]  branch {i} (omega ~ {float((om[0]+om[1])/2):+.5f}): radius {float(r):.2f}")

    # ---- [N] + [W]
    dh = {k: sp.diff(hs[k], t) for k in range(5)}
    m_cands = set()
    winfo = []
    for i in range(3):
        om, r = oms[i], tubes[i]
        F_lo = iabsmin(eval_K(sp.Poly(hs[4], t).coeff_monomial(t), om))
        check(F_lo > 0, "h4 slope bounded away from 0")
        dlo = F_lo - (tbound(dh[3], om, r)/A_FAR + tbound(dh[2], om, r)/Q(A_FAR)**2
                      + tbound(dh[1], om, r)/Q(A_FAR)**3 + tbound(dh[0], om, r)/Q(A_FAR)**4)
        check(dlo > 0, "d/db C34 dominated by its A^4 term on the tube")
        # residual of the FIVE-term tail (t0..t4): starts at X^6, so the real
        # root satisfies |b* - (wA + tail(1/A))| <= rho / A^5.  (Using the same
        # tail here and in the W2 series below keeps the anchors consistent --
        # a 4-term/5-term mismatch would silently drop a |dW2/db|*|t4|/A^4
        # term as large as the ones kept; caught in adversarial review.)
        tail5 = [vals[ts[i]] for i in range(5)]
        tail_expr = sum(tail5[k] * X**k for k in range(5))
        Pres = sp.Poly(sp.expand(C34.subs({A: 1/X, b: w/X + tail_expr}) * X**5), X)
        rb = Q(0)
        for k in range(Pres.degree() + 1):
            co = red(Pres.coeff_monomial(X**k))
            if co == 0:
                continue
            check(k >= 6, f"residual order (X^{k})")
            rb += iabsmax(eval_K(co, om)) * Q(1, A_FAR) ** (k - 6)
        rho = rb / dlo          # |b* - (wA + tail(1/A))| <= rho / A^5
        # the tail point itself must lie in the tube for dlo to apply
        tail_mag = sum(iabsmax(eval_K(tail5[k], om)) * Q(1, A_FAR)**k for k in range(5))
        check(tail_mag + rho/Q(A_FAR)**5 <= r, "tail point and root inside the tube")
        # [W]: |m - w0| <= wtail/A + |dW2/db| * rho/A^5  (both decreasing in A)
        wtail = Q(0)
        for k in range(5, PWX.degree() + 1):
            co = red(PWX.coeff_monomial(X**k))
            if co != 0:
                wtail += iabsmax(eval_K(co, om)) * Q(1, A_FAR) ** (k - 5)
        dW = sp.Poly(sp.expand(sp.diff(W2, b).subs(b, w*A + t)), A)
        DWrho = sum(tbound(red(dW.coeff_monomial(A**k)), om, r) * rho / Q(A_FAR) ** (5 - k)
                    for k in range(4))
        err = wtail / A_FAR + DWrho
        check(err < Q(49, 100), f"branch {i}: trap width {float(err)} < 0.49")
        w0i = eval_K(w0, om)
        lo_m, hi_m = w0i[0] - err, w0i[1] + err
        got = list(range(int(sp.ceiling(lo_m)), int(sp.floor(hi_m)) + 1))
        m_cands.update(got)
        winfo.append((float(w0i[0]), float(err), got))
        log(f"  [N/W] branch {i}: w0 ~ {float(w0i[0]):+.4f}, err {float(err):.2e}"
            f" -> m in {got}")

    # ---- [E]
    big_points = []
    for m_int in sorted(m_cands):
        RP = sp.Poly(sp.expand(sp.resultant(sp.Poly(C34, b), sp.Poly(W2 - m_int, b))), A)
        check(not RP.is_zero, "W2 nonconstant on the curve")
        for r_, _mult in RP.ground_roots().items():
            if r_.is_Integer and int(r_) >= A_FAR:
                Av = int(r_)
                gg = sp.gcd(sp.Poly(C34.subs(A, Av), b), sp.Poly((W2 - m_int).subs(A, Av), b))
                for br, _m2 in sp.Poly(gg, b).ground_roots().items():
                    if br.is_Integer:
                        big_points.append((Av, int(br)))
    log(f"  [E]  m in {sorted(m_cands)}: integer points with A >= A_FAR: {big_points}")

    return {"A_FAR": A_FAR, "m_cands": sorted(m_cands), "big_points": big_points,
            "tubes": tubes, "w0": w0}


# --------------------------------- [S] sweep ---------------------------------

def sweep(A_max: int, verbose: bool = True) -> list:
    """Every integer point of C34 with 0 <= A < A_max.  Modular sieve: the
    allowed b-residues mod m depend only on A mod m, so tables are precomputed;
    candidates are CRT-combined and verified exactly."""
    t_start = time.time()
    Pc = sp.Poly(C34, b)
    cpolys = [sp.Poly(c_, A) for c_ in Pc.all_coeffs()]

    def coeffs_at(Av):
        return [int(cp.eval(Av)) for cp in cpolys]

    tables = {}
    for m_ in SWEEP_MODS:
        tab = []
        for rA in range(m_):
            cs = [c_ % m_ for c_ in coeffs_at(rA)]
            ok = tuple(r_ for r_ in range(m_)
                       if (((((cs[0]*r_ + cs[1]) * r_ + cs[2]) * r_ + cs[3]) * r_
                            + cs[4]) * r_ + cs[5]) % m_ == 0)
            tab.append(ok)
        tables[m_] = tab

    # precompute CRT basis
    from math import prod
    M = prod(SWEEP_MODS)
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
        Bmax = max(3 * Av, 1000)
        for combo in itertools.product(*allowed):
            r_ = sum(c_ * e_ for c_, e_ in zip(combo, basis)) % M
            # unique representative of r_ mod M in [-Bmax, Bmax] (M > 2*Bmax+1)
            cand = r_ if r_ <= Bmax else r_ - M
            if cand < -Bmax or cand > Bmax:
                continue
            if cs is None:
                cs = coeffs_at(Av)
            v = ((((cs[0]*cand + cs[1])*cand + cs[2])*cand + cs[3])*cand + cs[4])*cand + cs[5]
            if v == 0:
                pts.append((Av, cand))
    pts = sorted(set(pts))
    if verbose:
        print(f"  [S]  sweep 0 <= A < {A_max}: points {pts}   ({time.time()-t_start:.0f}s)")
    return pts


# ------------------------- [D] degenerate locus, exact -----------------------

def degenerate(verbose: bool = True) -> list:
    """3a^2 - 2b - 14 = 0 with a(a^2-3b-7)(a^2-b-7) = 0: exactly a in {0,+-2}."""
    a_ = sp.Symbol("a_")
    b_ = (3*a_**2 - 14) / 2
    f2 = sp.expand(a_ * (a_**2 - 3*b_ - 7) * (a_**2 - b_ - 7))
    check(sp.expand(f2 - Q(7, 4)*a_**3*(a_**2 - 4)) == 0,
          "degenerate locus polynomial is (7/4) a^3 (a^2-4)")
    sols = []
    for av in (0, 2, -2):
        bv = (3*av*av - 14) // 2
        for dv in sp.solve(E1.subs({aa: av, b: bv}), dd):
            if dv.is_Integer:
                cv = int(sp.solve(E0.subs({aa: av, b: bv, dd: int(dv)}), cc)[0])
                sols.append((av, bv, int(dv), cv))
    if verbose:
        print("  [D]  degenerate locus a in {0,+-2} (exact):", sols)
    return sols


# --------------------------------- assembly ----------------------------------

def branch_b_candidates(points: list) -> dict:
    """Complete C34 point list -> the 3+4 values c on the generic locus."""
    out = {}
    for (Av, bv) in points:
        rt, ok = sp.integer_nthroot(max(Av, 0), 2)
        if Av < 0 or not ok:
            continue
        for av in {int(rt), -int(rt)}:
            den = 3*Av - 2*bv - 14
            if den == 0:
                continue
            num = -av*(Av - 3*bv - 7)*(Av - bv - 7)
            if num % den:
                continue
            dv = num // den
            check(sp.expand(E2.subs({aa: av, b: bv, dd: dv})) == 0, "E2 holds")
            check(sp.expand(E1.subs({aa: av, b: bv, dd: dv})) == 0, "E1 holds")
            cv = int(sp.solve(E0.subs({aa: av, b: bv, dd: dv}), cc)[0])
            out.setdefault(cv, []).append((av, bv, dv))
    return out


PROVED_BRANCH_B = True  # the theorem this script certifies end to end


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true", help="certificate only, no sweep")
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q26 Branch B (3+4): complete by Runge -- no Thue solver, no external CAS ===")
    cert = certificate()
    if args.quick:
        print("  (--quick: certificate verified; sweep skipped)")
        return 0
    pts = sweep(cert["A_FAR"])
    all_pts = sorted(set(pts) | set(cert["big_points"]))
    gen = branch_b_candidates(all_pts)
    deg = degenerate()
    cs = sorted(set(gen) | {s[3] for s in deg})
    print(f"  RESULT: complete list of 3+4 values: c in {cs}")
    check(sorted(v for v in cs if v != 0) == [-896, 896],
          f"nontrivial 3+4 values are exactly -+896")
    print("  Branch B is CLOSED unconditionally: c = -+896 only, and 5 divides")
    print("  neither, so neither has a root mod 5. No 3+4 candidate survives.")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "claim": "Branch B (3+4) of Q26 closed unconditionally by Runge",
            "A_FAR": cert["A_FAR"], "m_candidates": cert["m_cands"],
            "points": [list(p) for p in all_pts],
            "candidates": cs}, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
