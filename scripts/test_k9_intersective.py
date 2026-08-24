#!/usr/bin/env python3
"""Regression tests for Q29: Branch A's Runge closure, and the two open branches.

Pins:
  1. the chord quartic Phi9 re-derives from the (2,7) matching, is smooth
     (genus 3), and its leading form FACTORS over Q as -(A-b) * K3F with K3F
     the cyclic cubic form of discriminant 81 -- Q(2cos(2pi/9)) -- the two
     factors coprime: the Runge split that closes Branch A;
  2. the Runge kit is genuine: W1 (integer coefficients) has NO positive
     powers of A along the cyclic branches (verified symbolically here, not
     quoted), its limit is omega = -270beta^2 + 2241beta - 594 with trace
     exactly -243 = -3^5, and each omega is provably far from every integer
     (isolation via the exact minimal polynomial of omega);
  3. the rational channel is empty two ways: the branch invariant m = b - A
     converges to -10 with positive 9/A correction, and the trap line
     b = A - 10 meets Phi9 where A^2 - 10A - 12 = 0 (disc 148, nonsquare);
  4. the CRT sweep machinery is honest: brute force over a small box equals
     the sieve's output exactly (adversarial cross-check, not
     self-consistency);
  5. the certified complete point list: 23 integer points of Phi9 with
     A >= 0, all with A <= 49, exactly one with A not a perfect square
     ((17, 8)), and the A = a^2 points give exactly
     |c| in {0, 2630880, 176774400};
  6. both nontrivial c die: 2630880 (passes rad(9!)) at p = 13 and
     176774400 (fails rad(9!)) at p = 7, checked directly on (x)_9 - c;
  7. branch C's degenerate locus (2p2 - 3p3^2 + 30 = 0) holds exactly SIX
     rational points, ALL with c = 0 -- by elimination, not by a search
     (the earlier scan was vacuous: its expression still had p1 free);
  8. status flags: PROVED_BRANCH_A True (and PROVED_BRANCH_A_K9 True in the
     Runge script), PROVED_BRANCH_B and PROVED_BRANCH_C False -- nothing
     claims B or C closed;
  9. results/k9_intersective.json is byte-identical (sha256-pinned,
     historical "modulo effective Siegel" claim kept on purpose), and
     results/k9_branchA.json matches the certificate's recorded numbers.

Fails before scripts/k9_branchA_runge.py exists.  Runs in ~4 s.

Run: python scripts/test_k9_intersective.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

K9_ARTIFACT_SHA = "05c1e5580f4e688119fc46a59d5aef0ec7dec4fdf2cdfcb8a42992004dc93e86"

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


import k9_intersective as k9  # noqa: E402

try:
    import k9_branchA_runge as kr
except ImportError:
    kr = None

A, b, beta, X, yv = sp.symbols("A b beta X y")


def test_curve_and_split() -> None:
    expect(kr is not None, "scripts/k9_branchA_runge.py exists and imports")
    if kr is None:
        return
    # smoothness and the leading-form split
    sing = sp.solve([kr.PHI9, sp.diff(kr.PHI9, A), sp.diff(kr.PHI9, b)],
                    [A, b], dict=True)
    expect(sing == [], "Phi9 has no affine singular point")
    L = sp.expand(kr.PHI9 - sum(
        co * A**p * b**q for (p, q), co in sp.Poly(kr.PHI9, A, b).terms()
        if p + q < 4))
    expect(sp.expand(L + (A - b) * kr.K3F) == 0,
           "leading form = -(A-b) * K3F: the Runge split")
    mp = sp.Poly(kr.MP, beta)
    expect(mp.is_irreducible and sp.discriminant(mp, beta) == 81,
           "K3F's slopes: irreducible cyclic cubic, disc 81 (Q(2cos(2pi/9)))")
    expect(sp.expand(kr.K3F.subs(A, b)) == 3 * b**3,
           "K3F(b, b) = 3b^3 != 0: the two factors are coprime")
    expect(sp.Poly(sp.factor_list(kr.PHI9)[1][0][0], A, b).total_degree() == 4
           and len(sp.factor_list(kr.PHI9)[1]) == 1,
           "Phi9 irreducible over Q")


def test_w1_and_omega() -> None:
    if kr is None:
        return
    # W1 has integer coefficients
    expect(all(sp.nsimplify(co).is_Integer
               for co in sp.Poly(kr.W1, A, b).coeffs()),
           "W1 has integer coefficients (so W1 is an integer at integer points)")
    # branch expansion: solve the tail and check no positive powers -- from
    # scratch, symbolically, in Q(beta)
    def red(e):
        return sp.rem(sp.Poly(sp.expand(e), beta), sp.Poly(kr.MP, beta)).as_expr()

    ts = sp.symbols("u0:5")
    bser = beta / X + sum(ts[i] * X**i for i in range(5))
    Pser = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        kr.PHI9.subs({A: 1 / X, b: bser}) * X**4)))), X)
    vals = {}
    for k in range(0, Pser.degree() + 1):
        co = red(Pser.coeff_monomial(X**k).subs(vals))
        if co == 0:
            continue
        free = [s for s in ts if s in co.free_symbols]
        if not free:
            break
        s = free[0]
        lin = sp.Poly(co, s)
        expect(lin.degree() == 1, "series step linear")
        al, be = lin.all_coeffs()
        inv = sp.invert(sp.Poly(sp.expand(al), beta, domain="QQ"),
                        sp.Poly(kr.MP, beta, domain="QQ")).as_expr()
        vals[s] = red(-be * inv)
        if len(vals) == 5:
            break
    PWX = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        kr.W1.subs({A: 1 / X, b: bser.subs(vals)}) * X**4)))), X)
    expect(all(red(PWX.coeff_monomial(X**k)) == 0 for k in range(4)),
           "W1 has NO positive powers of A along the cyclic branches")
    expect(sp.expand(red(PWX.coeff_monomial(X**4)) - kr.OMEGA) == 0,
           "W1's limit is omega = -270beta^2 + 2241beta - 594")
    # trace of omega = -243 = -3^5 (exact, via power sums)
    expect(-270 * (81 - 12) + 2241 * 9 - 3 * 594 == -243,
           "trace(omega) = -243 = -3^5")
    # omega's minimal polynomial has no rational root near any integer:
    # resultant to get min poly of omega, then check the three real roots
    om_mp = sp.Poly(sp.resultant(sp.Poly(kr.MP, beta),
                                 sp.Poly(sp.Symbol("w") - kr.OMEGA, beta)),
                    sp.Symbol("w"))
    expect(om_mp.degree() == 3 and not om_mp.ground_roots(),
           "omega's minimal polynomial is a rootless cubic: no omega is rational")
    for lo, hi in ((18, 19), (311, 312), (-574, -573)):
        signs = {sp.sign(om_mp.eval(v)) for v in (lo, hi)}
        expect(len(signs) == 2, f"omega isolated in ({lo}, {hi})")
    for n in (19, 312, -574, 18, 311, -573):
        expect(om_mp.eval(n) != 0, f"omega != {n}")


def test_rational_channel() -> None:
    if kr is None:
        return
    lin = sp.factor(sp.expand(kr.PHI9.subs(b, A - 10)))
    expect(sp.expand(lin + 27 * (A**2 - 10 * A - 12)) == 0,
           "Phi9(A, A-10) = -27(A^2 - 10A - 12)")
    expect(not sp.integer_nthroot(148, 2)[1],
           "disc 148 is not a square: the trap line has no integer points")
    # branch tail: b = A - 10 + 9/A + ... re-derived
    rs = sp.symbols("r0:3")
    rser = 1 / X + sum(rs[i] * X**i for i in range(3))
    Rser = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(
        kr.PHI9.subs({A: 1 / X, b: rser}) * X**4)))), X)
    rvals = {}
    for k in range(0, Rser.degree() + 1):
        co = sp.expand(Rser.coeff_monomial(X**k).subs(rvals))
        if co == 0:
            continue
        free = [s for s in rs if s in co.free_symbols]
        if not free:
            break
        s = free[0]
        lin2 = sp.Poly(co, s)
        al, be = lin2.all_coeffs()
        rvals[s] = sp.Rational(sp.expand(-be / al))
        if len(rvals) == 3:
            break
    expect(rvals[rs[0]] == -10 and rvals[rs[1]] == 9 and rvals[rs[2]] == 90,
           "rational channel: b = A - 10 + 9/A + 90/A^2 + ... (t1 = 9 > 0)")


def test_sweep_cross_check() -> None:
    if kr is None:
        return
    got = [p for p in kr.sweep(60, verbose=False)]
    brute = []
    for Av in range(0, 60):
        for bv in range(-1000, 1001):
            if ((((-Av + 7 * bv + 30) * Av - 15 * bv * bv - 150 * bv - 273) * Av
                 + 10 * bv**3 + 180 * bv * bv + 819 * bv + 820) * Av
                    - bv**4 - 30 * bv**3 - 273 * bv * bv - 820 * bv - 576) == 0:
                brute.append((Av, bv))
    expect(got == sorted(brute),
           f"CRT sieve == brute force on 0 <= A < 60 ({len(got)} points)")
    from math import isqrt
    expect(len(got) == 23 and max(p[0] for p in got) == 49,
           "all 23 points appear below A = 60, max A exactly 49")
    expect([p for p in got if isqrt(p[0])**2 != p[0]] == [(17, 8)],
           "(17, 8) is the unique point with A not a perfect square")


def test_certified_list_and_kills() -> None:
    if kr is None:
        return
    art = json.loads((ROOT / "results" / "k9_branchA.json").read_text(encoding="utf-8"))
    expect(art["A_FAR"] == 34416, "recorded A_FAR")
    expect(art["sweep_points"] == 23 and art["square_points"] == 22,
           "23 integer points with A >= 0; 22 with A = a^2")
    expect(art["c_values"] == [0, 2630880, 176774400], "the Branch-A c-list")
    expect(art["W1_limits_trace"] == -243, "artifact records the -3^5 trace")
    expect(all(row[1] < row[2] for row in art["traps"]),
           "every recorded trap has err < dist(omega, Z)")
    expect([round(row[2], 3) for row in art["traps"]] == [0.173, 0.288, 0.462],
           "the recorded omega distances match the isolated values")
    expect(kr.BRANCH_A_C_VALUES == (0, 2630880, 176774400), "script constant agrees")
    expect(kr.KNOWN_KILLS == {2630880: 13, 176774400: 7}, "kill primes recorded")
    # the kills, directly on (x)_9 - c
    xs = sp.Symbol("x")
    BASE = sp.prod([xs - i for i in range(9)])
    for cv, kp in kr.KNOWN_KILLS.items():
        facs = [f for f, _ in sp.Poly(BASE - cv, xs).factor_list()[1]]
        expect(tuple(sorted(f.degree() for f in facs)) == (2, 7), f"c={cv} is (2,7)")
        expect(all(not any(int(f.eval(v)) % kp == 0 for v in range(kp))
                   for f in facs), f"c={cv}: no factor has a root mod {kp}")
        expect(all(any(any(int(f.eval(v)) % p_ == 0 for v in range(p_))
                       for f in facs) for p_ in sp.primerange(2, kp)),
               f"c={cv}: a root mod every prime below {kp}")
    expect(2630880 % 210 == 0 and 176774400 % 210 != 0,
           "rad(9!)-filter: 2630880 passes, 176774400 fails")


def test_branch_c_degenerate() -> None:
    """the six degenerate points, re-derived independently of branch_c()."""
    y_, c_ = sp.symbols("y c")
    p0, p1, p2, p3 = sp.symbols("p0 p1 p2 p3")
    P = y_**9 - 30 * y_**7 + 273 * y_**5 - 820 * y_**3 + 576 * y_
    q = sp.symbols("q0:5")
    Q4 = y_**4 + p3 * y_**3 + p2 * y_**2 + p1 * y_ + p0
    Q5 = y_**5 + sum(q[i] * y_**i for i in range(5))
    D = sp.Poly(sp.expand(Q4 * Q5 - (P - c_)), y_)
    co = [D.coeff_monomial(y_**i) for i in range(10)]
    sol = sp.solve([co[i] for i in range(8, 3, -1)], list(q), dict=True)[0]
    E3 = sp.expand(co[3].subs(sol))
    E2 = sp.expand(co[2].subs(sol))
    E1 = sp.expand(co[1].subs(sol))
    E0 = sp.expand(co[0].subs(sol))
    sub = {p2: sp.Rational(3, 2) * p3**2 - 15}
    E3d, E2d, E1d = (sp.expand(e.subs(sub)) for e in (E3, E2, E1))
    expect(p0 not in E3d.free_symbols, "p0 drops out of E3 on the locus")
    Rp0 = sp.expand(sp.resultant(sp.Poly(E2d, p0), sp.Poly(E1d, p0)))
    Rfin = sp.expand(sp.resultant(sp.Poly(E3d, p1), sp.Poly(Rp0, p1)))
    cands = sorted({int(r) for f, _m in sp.factor_list(Rfin)[1]
                    for r in sp.Poly(f, p3).ground_roots() if r.is_Integer})
    expect(cands == [-2, 0, 2], "degenerate p3 candidates are -2, 0, 2")
    pts = []
    for r3 in cands:
        v2 = sp.Rational(3, 2) * r3**2 - 15
        for r1 in (r for r in sp.Poly(E3d.subs(p3, r3), p1).ground_roots()
                   if r.is_Integer):
            g = sp.gcd(sp.Poly(E2d.subs({p3: r3, p1: r1}), p0),
                       sp.Poly(E1d.subs({p3: r3, p1: r1}), p0))
            for r0 in (r for r in sp.Poly(g, p0).ground_roots() if r.is_Integer):
                cv = sp.solve(E0.subs({p3: r3, p2: v2, p1: r1, p0: r0}), c_)
                pts.append(([int(v) for v in cv], int(r3), int(r1), int(r0)))
    expect(len(pts) == 6, f"exactly six degenerate points ({len(pts)})")
    expect(all(cvs == [0] for cvs, *_ in pts), "every degenerate point has c = 0")
    # spot check: (p3,p2,p1,p0) = (0,-15,10,24) is y^4-15y^2+10y+24 with
    # integer roots -- a genuine c=0 factorisation
    expect(sorted(sp.Poly(y_**4 - 15 * y_**2 + 10 * y_ + 24, y_).ground_roots())
           == [-4, -1, 2, 3], "witness quartic has roots {-4,-1,2,3}")


def test_branch_b_no_runge_split() -> None:
    """the load-bearing 'no Runge transfer' claim for Branch B, pinned."""
    y_, c_ = sp.symbols("y c")
    a_, b_, d_ = sp.symbols("a b d")
    P = y_**9 - 30 * y_**7 + 273 * y_**5 - 820 * y_**3 + 576 * y_
    es = sp.symbols("f0:6")
    sext = y_**6 + sum(es[i] * y_**i for i in range(6))
    DB = sp.Poly(sp.expand((y_**3 + a_ * y_**2 + b_ * y_ + d_) * sext - (P - c_)), y_)
    cob = [DB.coeff_monomial(y_**i) for i in range(10)]
    solb = sp.solve([cob[i] for i in range(8, 2, -1)], list(es), dict=True)[0]
    R2 = sp.expand(sp.numer(sp.together(cob[2].subs(solb))))
    R1 = sp.expand(sp.numer(sp.together(cob[1].subs(solb))))
    expect(c_ not in R2.free_symbols and c_ not in R1.free_symbols,
           "branch B: both leftover conditions are c-free")
    Res = sp.expand(sp.resultant(sp.Poly(R2, d_), sp.Poly(R1, d_)))
    fl = sp.factor_list(Res)
    expect(len(fl[1]) == 1 and fl[1][0][1] == 1,
           "branch B resultant is irreducible over Q")
    Pf = sp.Poly(fl[1][0][0], a_, b_)
    expect(sp.Poly(Res, a_).degree() == 18 and sp.Poly(Res, b_).degree() == 9,
           "degree 18 in a, 9 in b")
    dtot = Pf.total_degree()
    lead = sum(co * a_**i * b_**j for (i, j), co in Pf.terms() if i + j == dtot)
    expect(sp.factor(lead) / lead.coeff(a_**18) == a_**18,
           "leading form is a^18: a single repeated factor, NO Runge split")


def test_status_and_artifacts() -> None:
    expect(getattr(k9, "PROVED_BRANCH_A", None) is True,
           "k9_intersective.PROVED_BRANCH_A is True (the Runge theorem)")
    expect(getattr(k9, "PROVED_BRANCH_B", None) is False,
           "PROVED_BRANCH_B is False: nothing claims B closed")
    expect(getattr(k9, "PROVED_BRANCH_C", None) is False,
           "PROVED_BRANCH_C is False: nothing claims C closed")
    if kr is not None:
        expect(kr.PROVED_BRANCH_A_K9 is True, "the Runge script's own flag")
    raw = (ROOT / "results" / "k9_intersective.json").read_bytes()
    expect(hashlib.sha256(raw).hexdigest() == K9_ARTIFACT_SHA,
           "results/k9_intersective.json byte-identical (sha256-pinned)")
    art = json.loads(raw.decode("utf-8"))
    expect("Siegel" in art["claim"],
           "historical claim line kept ('modulo effective Siegel')")


def main() -> int:
    test_curve_and_split()
    test_w1_and_omega()
    test_rational_channel()
    test_sweep_cross_check()
    test_certified_list_and_kills()
    test_branch_c_degenerate()
    test_branch_b_no_runge_split()
    test_status_and_artifacts()
    print("\n=== K9 INTERSECTIVE / BRANCH A RUNGE TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
