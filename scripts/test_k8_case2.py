#!/usr/bin/env python3
"""Regression tests for the Q27 closure: the genus-3 curve is solved, Q27 proved.

Pins the load-bearing facts of scripts/k8_case2.py, each recomputed here
independently where cheap:

  1. the curve: disc((R(u)-R(B))/(u-B)) = -16*P6(B), y^2 = -B*P6(B) is
     squarefree of degree 7 (genus 3), P6(0) = -3^7 * 5^2 * 173;
  2. the descent: -B = delta*m^2 with delta a signed squarefree divisor of
     3*5*173 (16 classes), n = y/m, n^2 = delta*P6(-delta*m^2);
  3. beta >= 0 is compact (P6 <= 0 only on a bounded window) and the
     complete list is (0,0),(1,3840),(9,20736),(25,19200),(49,5376);
  4. all eight positive classes are empty by congruence: mod 7 for
     {1,15,519,865}, mod 5^4 for {5,2595}, mod 2^14 for 3, mod 2^9 for 173
     (2-adic depths 13 and 6 that a 2^7 battery misses); delta = 1 also
     dies by the Runge squeeze, kept as an independent cross-check;
  5. hence the curve's integral points are EXACTLY the five degenerate
     ones: no dangerous beta exists, Case 2 closes, and Q27 is PROVED;
  6. the second-road mod-8 family filter: min v2((x)_8) = 7, the exact
     v2(R(beta)) residue table (4 / 8 / 4 for beta == 3, 5, 7 mod 8), and
     the two delta == 7 mod 8 classes are {15, 519};
  7. the historical artifact results/k8_intersective.json is byte-identical
     (sha256-pinned) and the known candidates' record is unchanged;
     results/k8_case2.json's key fields are spot-checked.

Fails before scripts/k8_case2.py exists.  Runs in ~15 s.

Run: python scripts/test_k8_case2.py
"""

from __future__ import annotations

import hashlib
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
    import k8_case2 as kc
except ImportError:
    kc = None

u, B, m = sp.symbols("u B m")


def test_curve_and_descent() -> None:
    expect(kc is not None, "scripts/k8_case2.py exists and imports")
    if kc is None:
        return
    Rq = (u - 1) * (u - 9) * (u - 25) * (u - 49)
    qq = sp.Poly(sp.expand(sp.cancel((Rq - Rq.subs(u, B)) / (u - B))), u)
    expect(sp.expand(sp.discriminant(qq) + 16 * kc.P6) == 0, "disc(q) = -16*P6(B)")
    F = sp.Poly(sp.expand(-B * kc.P6), B)
    expect(F.degree() == 7 and sp.gcd(F, F.diff(B)) == 1 and (7 - 1) // 2 == 3,
           "y^2 = -B*P6(B): squarefree degree 7, genus 3")
    expect(-kc.p6_int(0) == 3**7 * 5**2 * 173, "P6(0) = -3^7*5^2*173")
    expect(sorted(kc.DELTAS_POS) == sorted(sp.divisors(3 * 5 * 173)),
           "the 8 positive squarefree divisor classes")
    # the descent's engine: for p outside {3,5,173}, p | B => p does not divide P6(B)
    for p_ in (2, 7, 11, 13, 17):
        expect(kc.p6_int(0) % p_ != 0, f"p={p_}: P6(0) is a unit, so v_p(-B) is even")


def test_compact_side() -> None:
    if kc is None:
        return
    rr = sp.Poly(kc.P6, B).real_roots()
    expect(len(rr) == 4 and all(-2 < r < 50 for r in rr), "P6: 4 real roots in (-2,50)")
    expect(kc.p6_int(-2) > 0 and kc.p6_int(50) > 0, "P6 > 0 at -2 and 50")
    pts = []
    for Bv in range(0, 51):
        val = -Bv * kc.p6_int(Bv)
        if val >= 0 and sp.integer_nthroot(val, 2)[1]:
            pts.append((Bv, int(sp.integer_nthroot(val, 2)[0])))
    expect(pts == [(0, 0), (1, 3840), (9, 20736), (25, 19200), (49, 5376)],
           "complete B >= 0 list: the five degenerate points")
    expect((0 - 1) * (0 - 9) * (0 - 25) * (0 - 49) % 256 != 0,
           "beta = 0 gives non-integral c")


def test_class_certificates() -> None:
    if kc is None:
        return
    expect(kc.DEAD_MOD == {1: 7, 15: 7, 519: 7, 865: 7, 5: 625, 2595: 625,
                           3: 16384, 173: 512},
           "the eight congruence certificates (delta = 1 included: P6(-m^2) == 3 mod 7)")
    for d, N in sorted(kc.DEAD_MOD.items()):
        sq = {(n_ * n_) % N for n_ in range(N)}
        found = any((d * kc.p6_mod(-d * mm * mm, N)) % N in sq for mm in range(N))
        expect(not found, f"delta={d}: empty mod {N} (re-scanned)")
    # the 2-adic depth facts behind the two deep certificates
    vals3 = {kc.v2(3 * kc.p6_int(-3 * mm * mm)) for mm in range(1, 60, 2)}
    expect(vals3 == {13}, "delta=3, odd m: v2 = 13 exactly (odd -> never a square)")
    u173 = {(173 * kc.p6_int(-173 * mm * mm) >> 6) % 8 for mm in range(1, 60, 2)}
    expect(all(kc.v2(173 * kc.p6_int(-173 * mm * mm)) == 6 for mm in range(1, 30, 2))
           and u173 == {5}, "delta=173, odd m: value = 2^6 * (5 mod 8), not a square")
    # Runge for delta = 1
    S = m**6 + 63 * m**4 + 651 * m**2 + 269
    Rm = sp.expand(sp.expand(kc.P6.subs(B, -m**2)) - S**2)
    expect(sp.expand(Rm - (336**2 * m**4 + 5429760 * m**2 - 9531136)) == 0,
           "Runge residual = 336^2 m^4 + 5429760 m^2 - 9531136")
    expect(2 * 56644 - 112770 == 518 and 518 * 56644 - 5428458 > 0,
           "squeeze dominates for |m| >= 238")
    expect(336**2 * 238**4 + 5429760 * 238**2 - 9531136 > 0,
           "lower half of the squeeze at the threshold")
    expect(all(not (kc.p6_int(-mm * mm) >= 0 and sp.integer_nthroot(kc.p6_int(-mm * mm), 2)[1])
               for mm in range(0, 238)), "delta=1: nothing below the squeeze threshold")


def test_mod8_filter() -> None:
    if kc is None:
        return
    worst = min(kc.v2(sp.prod([x_ - i for i in range(8)]))
                for x_ in range(-300, 301) if all(x_ - i != 0 for i in range(8)))
    expect(worst == 7, "min v2((x)_8) = 7 = v2(8!): intersective needs 2^7 | c")
    R8 = kc.R8
    expect({kc.v2(R8(b_)) for b_ in range(3 - 400, 3 + 401, 8)} == {4}
           and {kc.v2(R8(b_)) for b_ in range(7 - 400, 7 + 401, 8)} == {4},
           "beta == 3, 7 mod 8: v2(R) = 4 -> c not an integer")
    expect({kc.v2(R8(b_)) for b_ in range(5 - 400, 5 + 401, 8)} == {8},
           "beta == 5 mod 8: v2(R) = 8 -> c an ODD integer -> no root mod 2")
    expect(min(kc.v2(R8(b_)) for b_ in range(1 - 400, 1 + 402, 8) if R8(b_) != 0) >= 8,
           "beta == 1 mod 8: v2(R) >= 8 (the only viable residue)")
    expect(all(R8(b_) % 2 == 1 for b_ in range(-40, 41, 2)), "even beta: R odd")
    expect([d for d in kc.DELTAS_POS if (-d) % 8 == 1] == [15, 519],
           "negative beta == 1 mod 8 needs delta in {15, 519}: both empty mod 7")


def test_candidates_and_flags() -> None:
    if kc is not None:
        expect(kc.PROVED_CASE2 is True and kc.Q27_PROVED is True,
               "the flags claim exactly the theorem proved here")
    raw = (ROOT / "results" / "k8_intersective.json").read_bytes()
    expect(hashlib.sha256(raw).hexdigest()
           == "b400f9986704514c6a5c1c318352a5f220ec144e263945289f042081b9ed4a0b",
           "historical artifact byte-identical (sha256-pinned)")
    old = json.loads(raw.decode("utf-8"))
    expect(old["case_2_direct"] == {"candidates": 52, "dangerous": 0, "max_prime": 29},
           "the recorded candidates: 52 checked, 0 dangerous, all killed by p <= 29")
    expect(old["gap"] == "effective integral points on the genus-3 curve not computed",
           "its gap line predates the closure, deliberately")
    newp = ROOT / "results" / "k8_case2.json"
    expect(newp.exists(), "results/k8_case2.json exists")
    if newp.exists():
        art = json.loads(newp.read_text(encoding="utf-8"))
        expect(art["q27_proved"] is True, "artifact records the theorem")
        expect([tuple(p) for p in art["curve_integral_points"]]
               == [(0, 0), (1, 3840), (9, 20736), (25, 19200), (49, 5376)],
               "artifact records the complete curve point list")
        expect(art["class_certificates"]["3"] == "empty mod 16384"
               and art["class_certificates"]["173"] == "empty mod 512",
               "the two deep 2-adic certificates recorded")
    # the compact side's four nonzero points are the squares 1, 9, 25, 49:
    # excluded because t^2 = beta would give f_c a rational root
    expect(all(sp.integer_nthroot(b_, 2)[1] for b_ in (1, 9, 25, 49)),
           "the only dangerous-curve betas are perfect squares: degenerate")


def main() -> int:
    test_curve_and_descent()
    test_compact_side()
    test_class_certificates()
    test_mod8_filter()
    test_candidates_and_flags()
    print("\n=== K8 CASE 2 TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
