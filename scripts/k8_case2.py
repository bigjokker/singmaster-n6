#!/usr/bin/env python3
"""Q27 Case 2 CLOSED -- and with it, (x)_8 - c is NEVER intersective. Theorem.

docs/q27-k8-intersective.md left one gap: the dangerous locus of Case 2
(g = (u-beta) * irreducible cubic, c = R(beta)/256) is the genus-3 curve

    y^2 = -B * P6(B),
    P6(B) = B^6 - 126B^5 + 5271B^4 - 82564B^3 + 570591B^2 - 5779998B - 9458775,

with Siegel finiteness but no effective computation, and integral points
claimed (search to |B| <= 20000) to be exactly {0, 1, 9, 25, 49}.  It turns
out no effective machinery is needed: a divisor-class descent closes the
CURVE outright, and the family's own 2-adic structure gives a second road
that needs only a fraction of the certificates.  No external CAS, no
Chabauty, no search-as-proof.

THE DESCENT.  P6(B) = P6(0) mod B and P6(0) = -9458775 = -3^7 * 5^2 * 173,
so no prime outside {3, 5, 173} divides both B and P6(B) (2 never does:
P6(0) is odd).  At an integral point, v_p(-B) is therefore even for every
p outside {3, 5, 173}: writing -B = delta * m^2 with delta squarefree,
delta | 3*5*173 up to sign (16 classes), m | y, and n = y/m satisfies
n^2 = delta * P6(-delta * m^2).  Class by class:

  [C] delta <= 0  (i.e. beta >= 0, 8 classes): COMPACT.  P6 <= 0 only on
      (-1.41, 26.44) u (48.67, 49.68), so B = |delta| m^2 <= 49: the
      complete list of curve points with B >= 0 is
      (B,|y|) = (0,0), (1,3840), (9,20736), (25,19200), (49,5376).
  [M] delta in {1, 15, 519, 865}: EMPTY mod 7 -- for delta == 1 mod 7 this
      is the single computation P6(-m^2) == m^12 - m^6 + 3 == 3 (mod 7), a
      non-residue; delta = 1 ALSO dies by an independent Runge squeeze
      (S = m^6+63m^4+651m^2+269, S^2 < P6(-m^2) < (S+1)^2 for |m| >= 238,
      nothing below), kept as a cross-check.  delta in {5, 2595}: EMPTY
      mod 5^4.  delta = 3: EMPTY mod 2^14 -- for odd m the value
      3*P6(-3m^2) has 2-adic valuation exactly 13, odd, never a square;
      for even m it is 3 mod 8.  delta = 173: EMPTY mod 2^9 -- the value
      is 2^6 * (5 mod 8) for odd m, 5 mod 8 for even m.  (The 2-adic
      obstructions sit at depth 13 and 6: a congruence battery capped at
      2^7 wrongly calls both classes alive.  Certificates are full (m, n)
      scans mod 2^14 and 2^9, re-run on every execution.)

  THEOREM (the curve).  The integral points of y^2 = -B*P6(B) are exactly
  (0,0), (1,+-3840), (9,+-20736), (25,+-19200), (49,+-5376).  All are
  degenerate: 1, 9, 25, 49 are perfect squares (f_c would have a rational
  root) and B = 0 gives c = 11025/256, not an integer.  So NO dangerous
  beta exists, and the doc's Chebotarev kill applies to every genuine
  Case-2 value: CASE 2 IS CLOSED.

CROSS-CHECK: the mod-8 family filter closes Case 2 without the curve
theorem.  Intersective needs roots mod every power of 2, and (x)_8 == 0
mod 2^7 identically, so 2^7 | c.  For odd beta the exact 2-adic valuation
of R(beta) = (beta-1)(beta-9)(beta-25)(beta-49) is pinned by beta mod 8:
v2 = 4 for beta == 3, 7 (c not an integer), v2 = 8 for beta == 5 (c an ODD
integer: f_c is odd for every x -- no root even mod 2), and only
beta == 1 mod 8 survives; even beta gives R odd (c not an integer).  A
negative dangerous beta == 1 mod 8 would need -beta = delta m^2 with
delta == 7 mod 8, i.e. delta in {15, 519} -- both empty mod 7 -- and the
positive ones are the compact list.  Same conclusion.

CONSEQUENCE.  With Cases 1 and 4 (Jordan / two non-residue discriminants),
Case 5 ({1,9,25,49} only, c = 0) and Case 3 (empty) closed in the doc:

    (x)_8 - c IS NEVER INTERSECTIVE.  Unconditionally.  Q27 is PROVED,
    joining k <= 6.  The "modulo effective Siegel" label is retired.

    python scripts/k8_case2.py          (~1 s)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

u, B, m = sp.symbols("u B m")
P6 = (B**6 - 126*B**5 + 5271*B**4 - 82564*B**3 + 570591*B**2
      - 5779998*B - 9458775)
P6_COEFFS = [1, -126, 5271, -82564, 570591, -5779998, -9458775]
DELTAS_POS = (1, 3, 5, 15, 173, 519, 865, 2595)   # squarefree divisors of 3*5*173
DEAD_MOD = {1: 7, 15: 7, 519: 7, 865: 7, 5: 625, 2595: 625, 3: 16384, 173: 512}
# delta = 3 dies mod 2^14 and delta = 173 mod 2^9: the obstruction sits at
# 2-adic depth 13 resp. 6, which is why a battery capped at 2^7 called them
# alive -- valuation-deep obstructions need moduli past the valuation.
KNOWN_POINTS = [(0, 0), (1, 3840), (9, 20736), (25, 19200), (49, 5376)]
RUNGE_S = (1, 63, 651, 269)
RUNGE_M0 = 238
PROVED_CASE2 = True    # the theorem this script certifies
Q27_PROVED = True      # all five cases closed; see docs/q27-k8-intersective.md


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def p6_int(x: int) -> int:
    v = 0
    for cc in P6_COEFFS:
        v = v * x + cc
    return v


def p6_mod(x: int, N: int) -> int:
    v = 0
    for cc in P6_COEFFS:
        v = (v * x + cc) % N
    return v


def R8(z: int) -> int:
    return (z - 1) * (z - 9) * (z - 25) * (z - 49)


def v2(n: int) -> int:
    check(n != 0, "v2 of 0")
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def step0() -> None:
    print("  [0] the curve, re-derived, and the descent")
    Rq = (u - 1) * (u - 9) * (u - 25) * (u - 49)
    qq = sp.Poly(sp.expand(sp.cancel((Rq - Rq.subs(u, B)) / (u - B))), u)
    check(sp.expand(sp.discriminant(qq) + 16 * P6) == 0, "disc(q) = -16 * P6(B)")
    F = sp.Poly(sp.expand(-B * P6), B)
    check(F.degree() == 7 and sp.gcd(F, F.diff(B)) == 1, "degree 7, squarefree")
    check((F.degree() - 1) // 2 == 3, "genus 3")
    check(sp.Poly(P6, B).is_irreducible, "P6 irreducible")
    check(-p6_int(0) == 3**7 * 5**2 * 173, "P6(0) = -3^7 * 5^2 * 173")
    check([int(v) for v in sp.Poly(P6, B).all_coeffs()] == P6_COEFFS, "coefficients")
    # descent sanity: divisors list is complete
    divs = sorted({d for d in sp.divisors(3 * 5 * 173)})
    check(tuple(sorted(DELTAS_POS)) == tuple(divs), "the 8 squarefree divisors")
    print("      y^2 = -B*P6(B): squarefree, genus 3.  P6(0) = -3^7*5^2*173, so the")
    print("      squarefree part of -B divides 3*5*173: -B = delta*m^2, 16 classes,")
    print("      and n = y/m satisfies n^2 = delta*P6(-delta*m^2).")


def step_compact() -> list:
    print("  [C] beta >= 0 is COMPACT: the complete list")
    rr = sp.Poly(P6, B).real_roots()
    check(len(rr) == 4 and all(-2 < r < 50 for r in rr), "4 real roots, all in (-2,50)")
    check(p6_int(-2) > 0 and p6_int(50) > 0, "P6 > 0 at the window ends")
    pts = []
    for Bv in range(0, 51):
        val = -Bv * p6_int(Bv)
        if val >= 0 and sp.integer_nthroot(val, 2)[1]:
            pts.append((Bv, int(sp.integer_nthroot(val, 2)[0])))
    check(pts == KNOWN_POINTS, f"points with B >= 0: {pts}")
    check(all(sp.integer_nthroot(b_, 2)[1] for b_, _ in pts if b_ > 0),
          "1, 9, 25, 49 are perfect squares (rational root of f_c)")
    check(R8(0) % 256 != 0, "beta = 0: c = 11025/256 not an integer")
    print(f"      P6 <= 0 only on (-1.41,26.44) u (48.67,49.68) -> B <= 49; complete:")
    print(f"      (B,|y|) = {pts}; all degenerate")
    return pts


def step_mod8() -> None:
    print("  [8] the mod-8 verdict: a viable dangerous beta must be == 1 mod 8")
    # (x)_8 == 0 mod 2^7 identically (exhaustive check of the 2-adic valuation
    # of 8 consecutive integers, x mod 256 suffices: v2 depends on x mod 2^k
    # through the multiples of 2, 4, ... below 2^8 > 8)
    worst = min(v2(sp.prod([x_ - i for i in range(8)])) if all(x_ - i != 0 for i in range(8)) else 99
                for x_ in range(-260, 261))
    check(worst == 7, f"min v2((x)_8) over a full period is {worst}, = v2(8!) = 7")
    # exact valuations of R(beta) by beta mod 8 (odd beta): each factor's v2 is
    # exact because beta - r == 2 or 4 mod 8 pins it
    for res, expected in ((3, 4), (5, 8), (7, 4)):
        vals = {v2(R8(bb)) for bb in range(res - 800, res + 801, 8) if R8(bb) != 0}
        check(vals == {expected}, f"beta == {res} mod 8: v2(R) = {vals}, expected {expected}")
    seen1 = {v2(R8(bb)) for bb in range(1 - 800, 802, 8) if R8(bb) != 0}
    check(min(seen1) >= 8, f"beta == 1 mod 8: v2(R) >= 8 (saw {sorted(seen1)[:4]}...)")
    # even beta: R odd
    check(all(R8(bb) % 2 == 1 for bb in range(-50, 51, 2)), "even beta: R(beta) odd")
    print("      (x)_8 == 0 mod 2^7 for every x, so intersective needs 2^7 | c, and")
    print("      c = R(beta)/256:  beta even or == 3,7 mod 8  -> c not an integer;")
    print("      beta == 5 mod 8 -> v2(R) = 8 exactly -> c ODD -> f_c odd: no root mod 2.")


def step_neg1mod8() -> None:
    print("  [7] negative beta == 1 mod 8 does not exist on the curve")
    # -beta = delta m^2, m odd (even m makes beta even), so beta == -delta mod 8
    need = [d for d in DELTAS_POS if (-d) % 8 == 1]
    check(need == [15, 519], f"delta == 7 mod 8: {need}")
    for d in need:
        N = 7
        sq = {(n_ * n_) % N for n_ in range(N)}
        ok = any((d * p6_mod(-d * mm * mm, N)) % N in sq for mm in range(N))
        check(not ok, f"delta={d} should be empty mod 7")
        print(f"      delta = {d:3}: n^2 = {d}*P6(-{d} m^2) has NO solution mod 7 -- class empty")
    print("      so every dangerous beta is in {0, 1, 9, 25, 49}: all degenerate.")


def step_allclasses() -> None:
    print("  [M] every negative-beta class dies: the curve is COMPLETELY SOLVED")
    # Runge on delta = 1
    s6, s4, s2, s0 = RUNGE_S
    S = s6 * m**6 + s4 * m**4 + s2 * m**2 + s0
    Rm = sp.expand(sp.expand(P6.subs(B, -m**2)) - S**2)
    check(sp.expand(Rm - (336**2 * m**4 + 5429760 * m**2 - 9531136)) == 0, "Runge residual")
    check(sp.expand(2 * S + 1 - Rm
                    - (2*m**6 - 112770*m**4 - 5428458*m**2 + 9531675)) == 0, "squeeze")
    check(2 * 56644 - 112770 == 518 and 518 * 56644 - 5428458 > 0, "domination")
    check(336**2 * 238**4 + 5429760 * 238**2 - 9531136 > 0,
          "lower half of the squeeze at the threshold (increasing in |m| beyond)")
    found = [mm for mm in range(0, RUNGE_M0)
             if p6_int(-mm * mm) >= 0 and sp.integer_nthroot(p6_int(-mm * mm), 2)[1]]
    check(found == [], "delta = 1 empty below the squeeze threshold")
    print(f"      delta = 1: doubly dead -- mod 7 (below) AND the Runge squeeze")
    print(f"      (S^2 < P6(-m^2) < (S+1)^2 for |m| >= {RUNGE_M0}; nothing below)")
    # the seven congruence certificates, re-run in full every time
    for d, N in sorted(DEAD_MOD.items()):
        sq = {(n_ * n_) % N for n_ in range(N)}
        ok = any((d * p6_mod(-d * mm * mm, N)) % N in sq for mm in range(N))
        check(not ok, f"delta={d} should be empty mod {N}")
        print(f"      delta = {d:4}: no (m, n) mod {N:>5} -- class EMPTY")
    print("      (delta 3 and 173 die at 2-adic depth 13 and 6: for odd m,")
    print("       v2(3*P6(-3m^2)) = 13 -- odd, never a square -- and")
    print("       173*P6(-173m^2) = 2^6*(5 mod 8); a 2^7-capped battery misses both.)")


def step_theorem() -> None:
    print("  [T] the curve theorem, and Q27")
    print("      INTEGRAL POINTS OF y^2 = -B*P6(B):  exactly")
    print("      (0,0), (1,+-3840), (9,+-20736), (25,+-19200), (49,+-5376).")
    print("      All degenerate -> NO dangerous beta exists -> the Chebotarev kill")
    print("      applies to every genuine Case-2 value. Case 2 is CLOSED.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q27 Case 2: closed. (x)_8 - c is never intersective -- THEOREM ===")
    step0()
    pts = step_compact()
    step_mod8()
    step_neg1mod8()
    step_allclasses()
    step_theorem()
    print()
    print("  RESULT: the genus-3 curve is solved outright -- its integral points are")
    print("          the five degenerate ones -- so no dangerous beta exists and")
    print("          Case 2 closes.  Independently, the mod-8 family filter alone")
    print("          closes Case 2 (steps [8] + [7] + [C]).  With Cases 1, 3, 4, 5:")
    print("          (x)_8 - c IS NEVER INTERSECTIVE.  Unconditionally.  Q27 PROVED.")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "claim": "(x)_8 - c is never intersective -- PROVED; the genus-3 curve is solved",
            "curve_integral_points": [list(p) for p in pts],
            "class_certificates": {"1": "Runge squeeze",
                                   **{str(k): f"empty mod {v}" for k, v in sorted(DEAD_MOD.items())},
                                   "negative delta": "compact: P6 <= 0 only on a bounded window"},
            "mod8_family_filter": {"1": "only viable residue", "3": "c not integral",
                                   "5": "c odd: no root mod 2", "7": "c not integral"},
            "q27_proved": True,
        }, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
