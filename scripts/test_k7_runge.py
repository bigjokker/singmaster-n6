#!/usr/bin/env python3
"""Regression tests for Q26's two branches: Runge closure (B) and the corrected gap (A).

Pins the load-bearing facts of scripts/k7_runge.py and scripts/k7_branchA.py:

  1. the Runge split is exact: C34's leading form is Q2 * K3 with Q2 positive
     definite (disc -7) and K3(1,w) = w^3 - w^2 - 2w + 1, disc 49 (cyclic);
  2. the certificate holds end to end: A_FAR = 111184, all three tube radii
     below 10, W2's three branch limits (trace exactly -54) are bounded away
     from every integer, so NO integer point has A >= A_FAR;
  3. a small sweep reproduces the complete point list -- every one of the 19
     points of C34 (A >= 0) has A <= 36, so the certified list is closed;
  4. the degenerate locus is exactly a in {0, +-2} and carries c = -+896;
     the generic locus carries only c = 0: Branch B's values are {0, -+896};
  5. Branch A is a SMOOTH plane cubic (genus 1) -- not a Thue equation --
     with Jacobian Y^2 = X^3 - 1764X + 28224, trivial torsion, and (0,168)
     of infinite order. Magma Rank(E) 2026-08-24 printed "3 true";
  6. the 2+5 trichotomy: only S5 has odd derangements (20 of them); the
     D5/F20 trap is witnessed by (x^2+x-1)(x^5-2) having a root mod every
     prime tested while the mismatched control dies at 11;
  7. nothing overclaims: PROVED_BRANCH_B True, PROVED_BRANCH_A False,
     results/k7_intersective.json is byte-identical (sha256-pinned; its gap
     line is historical), and results/k7_runge.json matches the script.

Added 2026-08-24, when Branch A was re-attacked and stayed blocked:

  8. "not a Thue equation" is now a THEOREM, not an inference from
     smoothness (which is invalid -- N(x,y) = m z^3 is a smooth cubic too):
     the Hessian of the homogenised Phi restricted to z = 0 is
     392(b-2A)(2b-A)^2, neither identically zero nor divisible by the
     leading form, so no point at infinity is an inflection.  A genuine
     Thue cubic is pinned as the control;
  9. "not a Desboves cubic" too: psi_3 of E is irreducible, so E has no
     rational 3-isogeny and Magma's SIntegralDesbovesPoints cannot apply;
 10. the REAL half of Branch A is closed: a^2-4b > 0 forces |a| <= 6 and
     |b| <= 9 (a real chord cannot escape [-Y0,Y0], P7's largest local
     maximum being under 96), and that box holds exactly the 21 chord
     points, every one with c = 0;
 11. the explicit map sigma: Phi -> E is an identity modulo Phi, and it
     sends the genuine integral point (4,13) -- carrying c = +-17472 -- to
     (1345/4, 48959/8), on E but NOT integral.  Of the 21 known integral
     points of Phi (all |A| <= 95): 8 integral image, 10 finite
     non-integral, 3 poles.  This is the pin that stops anyone "closing"
     Branch A with IntegralPoints(E);
 12. 7! | c removes F20 from the no-kill door: |c| >= 5040 leaves f_c one
     real root, which the complex quadratic cannot hold, so disc(quintic)
     > 0 and its quadratic subfield is real while the door needs an
     imaginary one.  Only D5 survives.  Magma Rank is (3, True);
     PROVED_BRANCH_A stays False.

Fails before scripts/k7_runge.py and scripts/k7_branchA.py exist.
Runtime ~11 s, measured (the full certificate, but only a 41-value sweep).

Run: python scripts/test_k7_runge.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


try:
    import k7_runge as kr
except ImportError:
    kr = None
try:
    import k7_branchA as ka
except ImportError:
    ka = None

A, b, w, x = sp.symbols("A b w x")

EXPECTED_POINTS = [(0, -9), (0, -7), (0, -4), (0, -1), (1, -9), (1, -6), (1, -4),
                   (1, -2), (4, -9), (4, -5), (4, -3), (4, -1), (9, -4), (9, -1),
                   (9, 2), (16, 1), (16, 3), (25, 6), (36, 11)]


def test_runge_split() -> None:
    expect(kr is not None, "scripts/k7_runge.py exists and imports")
    if kr is None:
        return
    lead = sum(term for term in sp.Poly(kr.C34, A, b).as_expr().as_ordered_terms()
               if sp.Poly(term, A, b).total_degree() == 5)
    expect(sp.expand(lead - kr.Q2F * kr.K3F) == 0, "leading form = Q2 * K3, exactly")
    expect(sp.discriminant(sp.Poly(kr.Q2F.subs(A, 1), b), b) == -7,
           "Q2 positive definite (disc -7): the Q2 branches are complex")
    expect(sp.expand(kr.K3F.subs({A: 1, b: w}) - kr.MINPOL) == 0,
           "K3(1, w) IS the branch minimal polynomial w^3 - w^2 - 2w + 1")
    expect(sp.discriminant(sp.Poly(kr.K3F.subs(A, 1), b), b) == 49,
           "disc(K3) = 49: the cyclic cubic of Q(2cos(2pi/7))")


def test_certificate() -> None:
    if kr is None:
        return
    cert = kr.certificate(verbose=False)
    expect(cert["A_FAR"] == 111184, "A_FAR = 111184")
    expect(cert["m_cands"] == [], "no integer is trapped: W2's limits avoid Z")
    expect(cert["big_points"] == [], "no integer point with A >= A_FAR")
    expect(all(r < 10 for r in cert["tubes"]), "all tube radii below 10")
    w0 = cert["w0"]
    expect(sp.expand(w0 - (-294*w**2 + 69*w + 449)) == 0, "w0 = -294w^2 + 69w + 449")
    # trace of w0 over the cyclic cubic: tr(w)=1, tr(w^2)=5
    expect(-294*5 + 69*1 + 3*449 == -54, "trace of the three limits is -54")


def test_sweep_and_lists() -> None:
    if kr is None:
        return
    pts = kr.sweep(41, verbose=False)
    expect(pts == EXPECTED_POINTS, f"complete point list (A <= 40): {len(pts)} points")
    expect(max(p[0] for p in pts) == 36, "every point has A <= 36")
    gen = kr.branch_b_candidates(pts)
    expect(set(gen) == {0}, "generic locus carries only c = 0")
    deg = kr.degenerate(verbose=False)
    expect(sorted({s[3] for s in deg}) == [-896, 0, 896],
           "degenerate locus carries exactly c in {0, -+896}")
    expect((2, -1, 14, -896) in deg and (-2, -1, -14, 896) in deg,
           "the c = -+896 factorisations at (a,b,d) = (+-2,-1,+-14)")
    for cv in (896, -896):
        expect(cv % 5 != 0, f"5 does not divide {cv}: killed at 5")


def test_branch_a_geometry() -> None:
    expect(ka is not None, "scripts/k7_branchA.py exists and imports")
    if ka is None:
        return
    sing = sp.solve([ka.PHI, sp.diff(ka.PHI, A), sp.diff(ka.PHI, b)], [A, b], dict=True)
    expect(sing == [], "Phi = 0 has no affine singular point")
    expect(sp.discriminant(sp.Poly(b**3 - 6*b**2 + 5*b - 1, b), b) == 49,
           "leading form separable: 3 smooth points at infinity -> genus 1")
    expect(ka.JAC == (-1764, 28224), "Jacobian E: Y^2 = X^3 - 42^2 X + 168^2")
    m4, m6 = ka.JAC
    expect(168**2 == 0**3 + m4*0 + m6, "(0,168) on E")
    expect(sp.Poly(x**3 + m4*x + m6, x).is_irreducible, "no rational 2-torsion")
    expect(ka.PENCIL_QUARTIC == (25, -32, 306, -184, 29), "pencil quartic")
    # torsion trivial: gcd over two good primes with coprime counts suffices
    import math
    g = 0
    disc = -16*(4*m4**3 + 27*m6**2)
    for p in [11, 13, 17, 19, 23, 29]:
        if disc % p == 0:
            continue
        sq = {}
        for u_ in range(p):
            sq[u_*u_ % p] = sq.get(u_*u_ % p, 0) + 1
        n = 1 + sum(sq.get((u_**3 + m4*u_ + m6) % p, 0) for u_ in range(p))
        g = math.gcd(g, n)
    expect(g == 1, "gcd #E(F_p) = 1: torsion trivial, so (0,168) has infinite order")


def test_branch_a_not_thue_not_desboves() -> None:
    """The two exclusions that used to be assertions and are now theorems."""
    if ka is None:
        return
    z = sp.Symbol("z")
    F = sp.expand(ka.PHI.subs({A: A/z, b: b/z}) * z**3)
    H = sp.Matrix(3, 3, lambda i, j: sp.diff(F, [A, b, z][i], [A, b, z][j]))
    hess0 = sp.expand(H.det().subs(z, 0))
    expect(sp.expand(hess0 - 392*(b - 2*A)*(2*b - A)**2) == 0,
           "Hessian of the homogenised Phi on z=0 is 392(b-2A)(2b-A)^2")
    expect(hess0 != 0, "it is NOT identically zero -- so Phi is not affinely Thue")
    lead = sp.Poly(-(A**3) + 5*A**2*b - 6*A*b**2 + b**3, A, b)
    _, r = sp.div(sp.Poly(hess0, A, b), lead)
    expect(sp.expand(r.as_expr()) != 0,
           "and not divisible by the leading form: no infinite point is an inflection")
    # control: a genuine Thue cubic DOES have Hessian identically 0 at infinity
    xx, yy = sp.symbols("xx yy")
    G = xx**3 + 2*yy**3 - 7*z**3
    HG = sp.Matrix(3, 3, lambda i, j: sp.diff(G, [xx, yy, z][i], [xx, yy, z][j]))
    expect(sp.expand(HG.det().subs(z, 0)) == 0,
           "control: x^3+2y^3=7 has Hessian|_(z=0) identically 0")
    expect(ka.NOT_THUE is True, "k7_branchA records NOT_THUE")
    m4, m6 = ka.JAC
    psi3 = sp.Poly(3*x**4 + 6*m4*x**2 + 12*m6*x - m4**2, x)
    expect(tuple(int(v) for v in psi3.all_coeffs()) == ka.PSI3, "psi_3 of E")
    expect(psi3.is_irreducible and not psi3.ground_roots(),
           "psi_3 irreducible: no rational 3-isogeny, so Phi is not Desboves")
    expect(ka.NOT_DESBOVES is True, "k7_branchA records NOT_DESBOVES")


def test_branch_a_real_half_closed() -> None:
    """NEW theorem: a^2-4b > 0 admits only the 21 chords, all with c = 0."""
    if ka is None:
        return
    yv = sp.Symbol("yv")
    P7 = yv**7 - 14*yv**5 + 49*yv**3 - 36*yv
    crit = sp.real_roots(sp.Poly(sp.diff(P7, yv), yv))
    expect(len(crit) == 6, "P7 has six real critical points")
    expect(all(abs(sp.N(P7.subs(yv, r), 40)) < ka.CRIT_BOUND for r in crit),
           f"every |critical value| < {ka.CRIT_BOUND}")
    expect(len(sp.real_roots(sp.Poly(P7 - ka.CRIT_BOUND, yv))) == 1,
           "P7(y) = 96 has one real root, so |c| > 96 admits no real chord")
    Y0 = max(sp.N(r, 30) for r in sp.real_roots(sp.Poly(P7 - ka.CRIT_BOUND, yv)))
    expect(int(sp.floor(2*Y0)) == ka.REAL_A_MAX
           and int(sp.floor(Y0**2)) == ka.REAL_B_MAX,
           f"|a| <= {ka.REAL_A_MAX} and |b| <= {ka.REAL_B_MAX} on the real half")
    hits = []
    for av in range(-ka.REAL_A_MAX, ka.REAL_A_MAX + 1):
        for bv in range(-ka.REAL_B_MAX, ka.REAL_B_MAX + 1):
            if av*av - 4*bv > 0 and sp.expand(ka.PHI.subs({A: av*av, b: bv})) == 0:
                hits.append((av, bv, av*bv*(av*av - 3*bv - 7)*(av*av - bv - 7)))
    chords = {(-(i + j), i*j) for i in range(-3, 4) for j in range(i + 1, 4)}
    expect(len(hits) == 21 and {(h[0], h[1]) for h in hits} == chords,
           "the real half is exactly the 21 chord points")
    expect({h[2] for h in hits} == {0}, "and every one of them has c = 0")
    expect(ka.REAL_HALF_CLOSED is True, "k7_branchA records REAL_HALF_CLOSED")


def test_branch_a_sigma_and_the_wrong_list() -> None:
    """sigma lands on E, and IntegralPoints(E) provably misses a Phi-point."""
    if ka is None:
        return
    m4, m6 = ka.JAC
    W = sp.expand(ka.NY**2 - ka.NX**3 - m4*ka.NX*A**4 - m6*A**6)
    _, r = sp.div(sp.Poly(W, b), sp.Poly(ka.PHI, b))
    expect(sp.expand(r.as_expr()) == 0, "sigma: Phi -> E is an identity modulo Phi")
    (Av, bv), (Xv, Yv), cw = ka.SIGMA_WITNESS
    expect(sp.expand(ka.PHI.subs({A: Av, b: bv})) == 0, "(4,13) lies on Phi")
    expect(ka.sigma_of(Av, bv) == (Xv, Yv), "sigma(4,13) = (1345/4, 48959/8)")
    expect(Yv**2 == Xv**3 + m4*Xv + m6, "the image lies on E")
    expect(Xv.q != 1 and Yv.q != 1,
           "but is NOT an integral point of E: IntegralPoints(E) misses it")
    aw = sp.integer_nthroot(Av, 2)[0]
    expect(aw*bv*(aw*aw - 3*bv - 7)*(aw*aw - bv - 7) == cw,
           f"and that point carries the real Branch-A value c = {cw}")
    pa, pb = ka.SIGMA_POLE
    expect(sp.expand(ka.PHI.subs({A: pa, b: pb})) == 0,
           "sigma^-1(O) = (0,-9) is an affine rational point of Phi")
    expect(len(ka.PHI_INTEGRAL_POINTS) == 21, "21 recorded integral points of Phi")
    expect(max(abs(p[0]) for p in ka.PHI_INTEGRAL_POINTS) == 95,
           "every recorded point has |A| <= 95")
    expect(ka.classify_sigma_images() == (8, 10, 3),
           "sigma images: 8 integral, 10 finite non-integral, 3 poles")
    expect(ka.SIGMA_IMAGE_COUNTS == (8, 10, 3), "SIGMA_IMAGE_COUNTS matches")


def test_branch_a_f20_excluded() -> None:
    """NEW theorem: 7! | c leaves D5 alone in the no-kill door."""
    if ka is None:
        return
    yv = sp.Symbol("yv")
    P7 = yv**7 - 14*yv**5 + 49*yv**3 - 36*yv
    expect(5040 > ka.CRIT_BOUND, "7! exceeds every critical value of P7")
    for cv in (5040, -5040):
        expect(len(sp.real_roots(sp.Poly(P7 - cv, yv))) == 1,
               f"f_c has exactly one real root at c = {cv}")
    for poly in (x**5 - 2, x**5 + x + 1):
        p_ = sp.Poly(poly, x)
        expect(len(sp.real_roots(p_)) == 1 and sp.discriminant(p_, x) > 0,
               f"control: {poly} has one real root and positive discriminant")
    P7x = sp.prod([x - i for i in range(7)])
    for cv in (17472, 459648):
        q5 = [p_ for p_, _ in sp.Poly(P7x - cv, x).factor_list()[1]
              if p_.degree() == 5][0]
        d5 = sp.discriminant(q5, x)
        expect(len(sp.real_roots(q5)) == 1 and d5 > 0,
               f"c={cv}: quintic has one real root, disc > 0 (real quadratic subfield)")
        expect(not sp.integer_nthroot(int(d5), 2)[1],
               f"c={cv}: disc is not a square, so the quintic is not D5")
    expect(ka.SURVIVING_NO_KILL_GROUP == "D5", "only D5 survives the door")
    expect(ka.NO_KILL_GROUPS == ("D5", "F20"),
           "the group-theoretic trichotomy itself is unchanged (both still no-kill)")


def test_trichotomy() -> None:
    if ka is None:
        return
    rows = []
    for name in ("C5", "D5", "F20", "A5", "S5"):
        G = ka._group(name)
        der = [p for p in G if all(p[i] != i for i in range(5))]
        odd = [p for p in der if ka._sign(p) == -1]
        rows.append((name, len(G), len(der), len(odd)))
    expect([(r[1], r[2], r[3]) for r in rows] ==
           [(5, 4, 0), (10, 4, 0), (20, 4, 0), (60, 24, 0), (120, 44, 20)],
           "derangement table: odd ones exist in S5 only (20 of them)")
    h1 = sp.Poly(sp.expand((x*x + x - 1)*(x**5 - 2)), x).all_coeffs()
    expect(ka._rootless([int(v) for v in h1], 500) is None,
           "(x^2+x-1)(x^5-2): root mod every p < 500 -- the D5/F20 trap is real")
    h2 = sp.Poly(sp.expand((x*x + x + 1)*(x**5 - 2)), x).all_coeffs()
    expect(ka._rootless([int(v) for v in h2], 100) == 11,
           "(x^2+x+1)(x^5-2): mismatched subfield, killed at 11")


def test_status_and_artifacts() -> None:
    if kr is not None:
        expect(kr.PROVED_BRANCH_B is True, "Branch B claimed closed (it is)")
    if ka is not None:
        expect(ka.PROVED_BRANCH_A is False, "Branch A NOT claimed closed")
        expect(ka.MAGMA_RANK == (3, True), "Magma Rank(E) was 3 true")
        expect(ka.MAGMA_TORSION_ORDER == 1, "Magma TorsionSubgroup had order 1")
        expect("EllipticCurve([-1764, 28224])" in ka.MAGMA_COMMAND,
               "MAGMA_COMMAND names this Weierstrass model")
    import hashlib
    raw = (ROOT / "results" / "k7_intersective.json").read_bytes()
    expect(hashlib.sha256(raw).hexdigest() == "6f4ba80b505e7559d9f7ca1e861e17634e940ade91c1ceec60d9dd0aebe8081b",
           "historical artifact byte-identical (sha256 pinned)")
    old = json.loads(raw.decode("utf-8"))
    expect(old["gap"] == "Siegel gives finiteness; effective Thue not carried out",
           "its gap line predates the closure, deliberately")
    expect([r["c"] for r in old["candidates"]] == [17472, 459648, 896],
           "historical candidate list unchanged")
    newp = ROOT / "results" / "k7_runge.json"
    expect(newp.exists(), "results/k7_runge.json exists")
    if newp.exists():
        art = json.loads(newp.read_text(encoding="utf-8"))
        expect(art["A_FAR"] == 111184 and art["m_candidates"] == [],
               "artifact records the empty trap")
        expect([tuple(p) for p in art["points"]] == EXPECTED_POINTS,
               "artifact records the complete 19-point list")
        expect(art["candidates"] == [-896, 0, 896], "artifact candidate list")


def main() -> int:
    test_runge_split()
    test_certificate()
    test_sweep_and_lists()
    test_branch_a_geometry()
    test_branch_a_not_thue_not_desboves()
    test_branch_a_real_half_closed()
    test_branch_a_sigma_and_the_wrong_list()
    test_branch_a_f20_excluded()
    test_trichotomy()
    test_status_and_artifacts()
    print("\n=== K7 RUNGE / BRANCH A TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
