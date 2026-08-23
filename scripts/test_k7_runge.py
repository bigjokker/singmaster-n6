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
     of infinite order: rank >= 1, no rank-0 shortcut;
  6. the 2+5 trichotomy: only S5 has odd derangements (20 of them); the
     D5/F20 trap is witnessed by (x^2+x-1)(x^5-2) having a root mod every
     prime tested while the mismatched control dies at 11;
  7. nothing overclaims: PROVED_BRANCH_B True, PROVED_BRANCH_A False,
     results/k7_intersective.json is byte-identical (sha256-pinned; its gap
     line is historical), and results/k7_runge.json matches the script.

Fails before scripts/k7_runge.py and scripts/k7_branchA.py exist.
Runtime ~40 s (the full certificate, but only a 41-value sweep).

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
