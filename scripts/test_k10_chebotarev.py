#!/usr/bin/env python3
"""Regression tests for Q28's Chebotarev pass on the (2,3) branch.

Pins the trichotomy and the two unconditional results of
scripts/k10_chebotarev.py, each recomputed here independently where cheap:

  1. the no-kill shape is real: a synthetic (2,2,3,3) product with the
     square-class coincidence has a root mod every prime tested, while the
     same product without the coincidence, the shared-subfield case with
     k(t^2) irreducible, and both recorded candidates all die;
  2. c < 0 closes: R is strictly increasing on (-inf, 1], its largest
     critical value magnitude is under 43,930,543 < 1024*10!, and for
     negative multiples of 10! the quintic has exactly one real root,
     negative -- incompatible with the real root of the cubic factor m;
  3. the no-kill locus is the plane curve F(p4, p2) = 0 with p4 odd,
     c = (p0^2 - 945^2)/1024, and (for 10! | c) 315 | p0 and
     p0 = +-945 mod 2^17; small sweeps find only c = 0 on it;
  4. Chebotarev still does not kill the locus by itself; Magma did.
     PROVED and PROVED_23_POSITIVE are True; the committed json is the
     historical pre-Magma artifact.

Fails before scripts/k10_chebotarev.py exists. Runs in a few seconds.

Run: python scripts/test_k10_chebotarev.py
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
    import k10_chebotarev as kc
except ImportError:
    kc = None

import k10_intersective as k10  # noqa: E402


def test_trichotomy_examples() -> None:
    expect(kc is not None, "scripts/k10_chebotarev.py exists and imports")
    if kc is None:
        return
    # the no-kill example: coincidence holds and no prime below 500 kills it
    h = kc.deep_h(1, 44, 1, 2, 3)
    expect(kc.coincidence(1, 44, 1, 2, 3), "example (1,44|1,2,3): (mu^2-4nu)*disc(m) is a square")
    expect((1 - 4 * 44) == kc.disc_m(1, 2, 3) == -175, "both discriminants are -175")
    expect(kc.first_rootless(h, 500) is None, "no-kill shape: root mod every p < 500")
    t = sp.Symbol("t")
    q_u = sp.Poly(sp.expand((t * t + t + 44) * (t * t - t + 44)).subs(t**2, sp.Symbol("u")), sp.Symbol("u"))
    expect(q_u.is_irreducible, "its q is an irreducible quadratic in u")
    # without the coincidence the same shape dies at an odd prime (11)
    expect(kc.first_rootless(kc.deep_h(1, 1, 1, 2, 3), 500, lo=3) == 11,
           "same shape, independent fields: first odd kill at 11")
    # case (III) with q(t^2) irreducible: 3-cycle-class kill at 43
    u0 = sp.Symbol("u")
    t0 = sp.Symbol("t")
    h5 = ((u0 * u0 + 2 * u0 + 176) * (u0**3 + 3 * u0 * u0 - 2 * u0 - 9)).subs(u0, t0 * t0)
    expect(kc.first_rootless(sp.expand(h5), 100, lo=3) == 43,
           "case (III), q(t^2) irreducible: killed at 43")
    # the identities behind the coincidence
    mu, nu, al, be, ga = sp.symbols("mu nu al be ga")
    expect(sp.expand((2 * nu - mu * mu) ** 2 - 4 * nu * nu
                     - mu * mu * (mu * mu - 4 * nu)) == 0, "D_q = mu^2(mu^2-4nu)")
    kk = sp.Poly(u0**3 + (2 * be - al * al) * u0 * u0 + (be * be - 2 * al * ga) * u0 - ga * ga, u0)
    mm = sp.Poly(t0**3 + al * t0 * t0 + be * t0 + ga, t0)
    expect(sp.expand(sp.discriminant(kk, u0) - sp.discriminant(mm, t0) * (ga - al * be) ** 2) == 0,
           "D_k = disc(m)(gamma - alpha*beta)^2: coincidence = D_q*D_k square")
    # shared subfield with k(t^2) irreducible dies (tau-class)
    u = sp.Symbol("u")
    g2 = sp.expand((u * u + u + 8) * (u**3 + u + 1))
    expect(int(sp.discriminant(sp.Poly(u * u + u + 8, u), u)) ==
           int(sp.discriminant(sp.Poly(u**3 + u + 1, u), u)) == -31, "shared subfield fixture")
    expect(kc.first_rootless(g2.subs(u, t * t), 100, lo=3) == 11, "case (II): killed at 11")
    # the recorded candidates are generic and still die at 11 and 13
    for cv, kp in kc.CANDIDATES.items():
        P = 1
        for pp, _ in sp.Poly(kc.R5 - 1024 * cv, u).factor_list()[1]:
            P *= int(sp.discriminant(pp.as_poly(u), u))
        expect(P > 0 and not sp.integer_nthroot(P, 2)[1], f"c={cv}: D_q*D_k > 0, not a square")
        expect(kc.first_rootless((kc.R5 - 1024 * cv).subs(u, t * t), 40, lo=3) == kp,
               f"c={cv}: killed at {kp}")


def test_c_negative_closes() -> None:
    if kc is None:
        return
    u = sp.Symbol("u")
    R = kc.R5
    rts = sp.Poly(sp.diff(R, u), u).real_roots()
    expect(len(rts) == 4 and all(1 < r < 81 for r in rts), "R' has 4 real roots, all in (1,81)")
    M = max(abs(R.subs(u, r)) for r in rts)
    expect(M < 43930543 and not M < 43930542, "largest |critical value| in (43930542, 43930543)")
    expect(1024 * 3628800 > 43930543, "1024*10! clears the critical range")
    P = sp.Poly(R - 1024 * (-3628800), u)
    expect(P.count_roots(1, 81) == 0 and P.count_roots(81, 10**9) == 0
           and P.count_roots(-10**9, 0) == 1,
           "c = -10!: exactly one real root and it is negative")
    expect(kc.PROVED_23_NEGATIVE is True, "the c < 0 half of the branch is closed")


def test_no_kill_locus_curve() -> None:
    if kc is None:
        return
    # a known c = 0 point of F: s with roots {1,-3,5,-7,9} has (p4,p2,p0)=(-5,-230,-945);
    # the sweep's normalised sign gives (5,-230,945) -- both must lie on F.
    expect(kc.Fcurve(5, -230) == 0 and kc.Fcurve(-5, 230) == 0, "(5,-230) is on F (a c=0 point)")
    expect(kc.Fcurve(1, 14) == 0 and kc.Fcurve(3, -82) == 0, "more c=0 points on F")
    expect(kc.Fcurve(5, -231) != 0, "F is not identically zero nearby")
    # p4 even is impossible: p3 = (p4^2-165)/2
    expect((2**2 - 165) % 2 != 0 and (1 - 165) % 2 == 0, "p3 integral forces p4 odd")
    # congruences for 10! | c: 315 | p0, p0 odd, p0 = +-945 mod 2^17
    good = 945 + 2**17 * 3  # = 945 mod 2^17 and odd
    expect(good % 2 == 1 and (good * good - 945 * 945) % 2**18 == 0,
           "p0 = 945 mod 2^17 gives 2^18 | p0^2 - 945^2")
    bad = 945 + 2**16
    expect((bad * bad - 945 * 945) % 2**18 != 0, "p0 = 945 + 2^16 does not")
    p0x = 945 + 2**17  # odd, = 945 mod 2^17, but 315 does not divide it
    cx = (p0x * p0x - 945 * 945) // 1024
    expect((p0x * p0x - 945 * 945) % 1024 == 0 and p0x % 315 != 0 and cx % 81 != 0,
           "p0 = 945 + 2^17: c is an integer but 3^4 does not divide it -- 315 | p0 is real")
    # F has integral points whose p0 is NOT integral: elimination artifacts on p2 = -p4^3
    expect(kc.Fcurve(1, -1) == 0 and kc.Fcurve(9, -729) == 0,
           "spurious F points at p2 = -p4^3 exist")
    num1 = -8 * 1 + 1 * (8 * (-1) * 1 - 1320 * (-1) - 1 + 495 - 46563) + 81125
    expect(num1 % 16 != 0, "their p0 = num/(16 p4) is not an integer at (1,-1)")
    # small sweeps: only c = 0 (step3 raises on any nonzero-c deep point)
    r = kc.step3(1200, 301)
    expect(r["deep_nonzero"] == [] or all(cv == 0 for *_, cv in r["deep_nonzero"]),
           "nu-sweep: no deep point with c != 0")
    expect(set(r["p0s"]) <= {945, -945}, "F-sweep: every point has p0 = +-945 (c = 0)")
    expect(46 in r["b_square_a"] and 158 in r["b_square_a"],
           "the two b-square, non-deep integral points at a = 46, 158 are seen")


def test_genus_three_cover() -> None:
    if kc is None:
        return
    u = sp.Symbol("u")
    qC = sp.Poly(u**4 + 165 * u**3 + 8778 * u * u + 172810 * u + 1057221, u)
    expect(sp.discriminant(qC, u) != 0, "b's four zeros on C are simple (q_C squarefree)")
    expect(sp.resultant(sp.Poly(3 * u * u + 330 * u + 8778, u), qC) == -5350630125159,
           "s(a) shares no root with q_C: the branch points are smooth on C")
    expect(int(2 * 1 - 2) * 2 + 4 == 2 * 3 - 2, "Riemann-Hurwitz: 2g-2 = 2*0 + 4 -> genus 3")


def test_status_proved_by_magma_not_chebotarev() -> None:
    if kc is not None:
        expect(kc.PROVED_23_POSITIVE is True, "c > 0 closed by Magma's 15 a-values")
    expect(k10.PROVED is True, "k10_intersective.PROVED True after Magma")
    art = json.loads((ROOT / "results" / "k10_intersective.json").read_text(encoding="utf-8"))
    expect("modulo" in art["claim"] and art["gap"] == k10.GAP,
           "committed json is still the pre-Magma artifact")


def main() -> int:
    test_trichotomy_examples()
    test_c_negative_closes()
    test_no_kill_locus_curve()
    test_genus_three_cover()
    test_status_proved_by_magma_not_chebotarev()
    print("\n=== K10 CHEBOTAREV TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
