#!/usr/bin/env python3
"""Regression tests for Q28's (2,3) branch: what is established, and what is not.

docs/q28-k10-intersective.md reduces the (2,3) factorisation of
g(u) = R(u) - 1024c to integral points on the genus-1 quartic

    y^2 = 5a^4 + 1320a^3 + 126456a^2 + 5102240a + 72824400.

These pins fix the facts the write-up now rests on, each recomputed here
independently of scripts/k10_intersective.py where that is cheap:

  1. the curve is what the matching says it is (re-derived from the
     coefficient equations, not copied from the script's Cab);
  2. the complete list of integral points with |a| <= 10^7 is 15 values of
     a, none beyond |a| = 730 -- so the |a| <= 6000 search missed nothing;
  3. the Jacobian is E: Y^2 = X^3 - 792X + 9801 with E(Q) = Z/2 x Z^2,
     rank EXACTLY 2 by 2-isogeny descent (no undecided Selmer class), so
     there is no rank-0 shortcut to a finite rational-point list;
  4. the two candidates still die, checked directly on (x)_10 - c mod p:
     c = 1395418752000 at p = 11 and c = 2235340800 at p = 13, each with
     a root mod every smaller prime;
  5. the script does NOT claim the list complete: PROVED is False, the
     artifact's gap string is unchanged, and the named Magma call carries a
     point that really lies on the curve.

Fails on the pre-step-[5] script (no jacobian / integral_points / PROVED).
Runs in a few seconds; touches nothing under results/.

Run: python scripts/test_k10_intersective.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import k10_intersective as k10  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


Q = (5, 1320, 126456, 5102240, 72824400)
EXPECTED_A = [-730, -250, -130, -106, -90, -82, -74, -58, -50, -34, -26, -10, 46, 54, 158]
CANDIDATES = {1395418752000: 11, 2235340800: 13}


def test_curve_rederives_from_the_matching() -> None:
    u, a, b, d, e, f, c = sp.symbols("u a b d e f c")
    R = sp.expand(sp.prod([u - (2 * j - 1) ** 2 for j in range(1, 6)]))
    prod = sp.expand((u**2 + a * u + b) * (u**3 + d * u**2 + e * u + f))
    P = sp.Poly(prod - (R - 1024 * c), u)
    eq = {n: P.coeff_monomial(u**n) for n in range(5)}
    sol = {d: sp.solve(eq[4], d)[0]}
    sol[e] = sp.solve(eq[3].subs(sol), e)[0]
    sol[f] = sp.solve(eq[2].subs(sol), f)[0]
    C = sp.Poly(sp.expand(eq[1].subs(sol)), b)
    expect(C.degree() == 2, "u^1 constraint is quadratic in b")
    D = sp.Poly(sp.discriminant(C, b), a)
    expect([int(v) for v in D.all_coeffs()] == list(Q), f"discriminant is the quartic {Q}")
    expect(sp.gcd(D, D.diff(a)) == 1, "quartic squarefree (genus 1)")
    expect(sp.factor_list(D.as_expr())[1][0][0].as_poly(a).degree() == 4, "irreducible over Q")
    expect(tuple(getattr(k10, "QUARTIC", ())) == Q, "script's QUARTIC matches")


def test_integral_points_complete_to_1e7() -> None:
    fn = getattr(k10, "integral_points", None)
    expect(callable(fn), "k10_intersective exposes integral_points(a_max)")
    if not callable(fn):
        return
    pts = fn(10**7)
    q = k10.quartic
    expect([av for av, _ in pts] == EXPECTED_A, f"15 integral points, a = {EXPECTED_A}")
    expect(all(yv >= 0 and yv * yv == q(av) for av, yv in pts), "each is on the curve")
    expect(max(abs(av) for av, _ in pts) == 730, "largest |a| is 730: the 6000 search missed nothing below 1e7")
    pairs = [(av, bv, cv) for av, yv in pts for bv, cv in k10.point_to_c(av, yv)]
    expect(len(pairs) == 30, "every point gives two integral (b, c) pairs")
    surv = sorted(cv for _, _, cv in pairs if cv and cv % 3628800 == 0)
    expect(surv == sorted(CANDIDATES), f"exactly the two candidates have 10! | c: {surv}")
    expect(sum(1 for _, _, cv in pairs if cv == 0) == 10, "ten pairs give c = 0")
    expect(all(cv % k10.RAD10 for _, _, cv in pairs if cv and cv not in CANDIDATES),
           "every other nonzero c fails rad(10!) | c")
    # independent recount with a different square test
    from math import isqrt
    n = sum(1 for av in range(-1500, 1501) if q(av) >= 0 and isqrt(q(av)) ** 2 == q(av))
    expect(n == 15, "math.isqrt recount over |a| <= 1500 agrees: 15")


def test_jacobian_rank_two() -> None:
    fn = getattr(k10, "jacobian", None)
    expect(callable(fn), "k10_intersective exposes jacobian()")
    if not callable(fn):
        return
    j = fn()
    expect((j["a4"], j["a6"]) == (-792, 9801), "Jacobian E: Y^2 = X^3 - 792X + 9801")
    expect(j["I"] == 155713536 and j["J"] == -4439704338432, "I, J invariants")
    expect(j["torsion"] == 2, "E(Q)_tors = Z/2")
    expect(j["rank"] == 2 and k10.JAC_RANK == 2, "rank 2")
    expect(j["im_alpha"] == [1, 3, 5, 11, 15, 33, 55, 165], "im alpha: all 8 positive classes")
    expect(j["im_alpha_prime"] == [1, -11], "im alpha': {1, -11}")
    expect(sorted(int(d) for d in j["obstructed"]) == [-33, -3, -1, 3, 11, 33]
           and all("Q_3" in v for v in j["obstructed"].values()),
           "the six other classes of E' have no 3-adic points")
    # the descent's ingredients, independently: shift, isogenous curve, 2-torsion
    X = sp.symbols("X")
    expect(sp.expand((X - 33) ** 3 - 792 * (X - 33) + 9801 - (X**3 - 99 * X**2 + 2475 * X)) == 0,
           "shift by 33 gives Y^2 = X^3 - 99X^2 + 2475X")
    expect(99 * 99 - 4 * 2475 == -99 and not k10.issq(198 * 198 + 4 * 99),
           "E and E' each have exactly one rational 2-torsion point (discs -99, 39600)")
    # 3-adic obstruction, by hand, for d = 3 on E': 3w^4 + 198w^2z^2 - 33z^4 = N^2
    found = any((3 * w**4 + 198 * w * w * z * z - 33 * z**4) % 9 in {n * n % 9 for n in range(9)}
                for w in range(9) for z in range(9) if w % 3 or z % 3)
    expect(not found, "d = 3 on E': no solution mod 9 with (w, z) not both divisible by 3")
    for Xp, Yp in j["points_of_infinite_order"]:
        expect(Yp * Yp == Xp**3 - 792 * Xp + 9801 and (Xp, Yp) != (-33, 0),
               f"({Xp}, {Yp}) is on E and is not the torsion point")


def test_two_kills_regenerate() -> None:
    for cv, kp in CANDIDATES.items():
        expect(k10.kill_prime(cv) == kp, f"script: kill_prime({cv}) == {kp}")
        degs = tuple(sorted(pp.degree() for pp, _ in k10.gp(cv).factor_list()[1]))
        expect(degs == (2, 3), f"c={cv} splits (2,3)")
        expect(cv % 3628800 == 0, f"c={cv} passes the real necessary condition 10! | c")

        def has_root(p):
            return any(sp.prod([w - i for i in range(10)]) % p == cv % p for w in range(p))

        expect(not has_root(kp), f"(x)_10 - {cv} has no root mod {kp} (direct, not via g)")
        expect(all(has_root(p) for p in sp.primerange(2, kp)),
               f"(x)_10 - {cv} has a root mod every prime below {kp}")


def test_status_is_not_proved() -> None:
    expect(getattr(k10, "PROVED", None) is False, "PROVED is False: the list is not certified complete")
    art = json.loads((ROOT / "results" / "k10_intersective.json").read_text(encoding="utf-8"))
    expect(art["gap"] == getattr(k10, "GAP", None), "artifact gap string unchanged")
    expect("modulo" in art["claim"], "artifact claim still says 'modulo'")
    expect(art["case_2_3"]["candidates"] == [{"a": -730, "c": 1395418752000, "kill_prime": 11},
                                             {"a": -250, "c": 2235340800, "kill_prime": 13}],
           "artifact's two candidates and kills")
    cmd = getattr(k10, "MAGMA_COMMAND", "")
    expect("IntegralQuarticPoints" in cmd and str(list(Q)) in cmd, "Magma call names this quartic")
    expect("[-250, 74880]" in cmd and 74880**2 == sum(qi * (-250) ** (4 - i) for i, qi in enumerate(Q)),
           "its base point lies on the curve")


def main() -> int:
    test_curve_rederives_from_the_matching()
    test_integral_points_complete_to_1e7()
    test_jacobian_rank_two()
    test_two_kills_regenerate()
    test_status_is_not_proved()
    print("\n=== K10 INTERSECTIVE TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
