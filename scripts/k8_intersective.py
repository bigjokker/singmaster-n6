#!/usr/bin/env python3
"""Q27: is (x)_8 - c ever intersective?  No -- PROVED (gap closed 2026-08-23).

k=8 is EVEN, so Q25's reduction applies: with t = 2x-7 the roots 0..7 become
+-1,+-3,+-5,+-7, and

    2^8 (x)_8 = (t^2-1)(t^2-9)(t^2-25)(t^2-49) = R(t^2),
    R(u) = u^4 - 84u^3 + 1974u^2 - 12916u + 11025,

so 2^8 f_c = g(t^2) with g(u) = R(u) - 256c, and the same equivalence holds:

    f_c has a root mod p  <=>  g has a root mod p that is a SQUARE mod p.

But g is a QUARTIC now, not a cubic, so the trichotomy of Q25 becomes five
cases on the factorisation of g over Q. Two need no computation:

  CASE 1  g irreducible. Transitive => Jordan => a derangement => Chebotarev
          gives a positive density of p with no root of g at all. Dead.

  CASE 4  g = two irreducible quadratics. Both discriminants are non-squares,
          so a positive density of p has BOTH non-residues (1/4, or 1/2 if
          their product is a square), and there g has no root at all. Dead.
          No condition needed -- this case can never be intersective.

Two are finite, and both collapse to c = 0:

  CASE 5  g splits over Q. e1=84, e2=1974, e3=12916 are FIXED, so with
          s = beta1+beta2 the product is forced,
          p(s) = (-s^3+84s^2-1974s+12916)/(84-2s). What bounds s is
          INTEGRALITY: s even with (s-42) | 2048 (22 values; the old
          "disc >= 0 bounds s to [10,74]" reason was wrong -- disc >= 0
          holds on [9,35] u [43,81] -- caught in the 2026-08-23 review;
          the enumeration was complete). Only solution: {1,9,25,49}, c = 0.

  CASE 3  exactly two rational roots. Same p(s) relation, same bound. All six
          solutions have c = 0 and a reducible remainder, i.e. they are
          Case 5 in disguise. Case 3 proper is EMPTY.

That leaves one infinite family, exactly as k=6 did:

  CASE 2  g = (u - beta) q(u), q an irreducible cubic.
          c = R(beta)/256, and q is the divided difference
          q(u) = (R(u)-R(beta))/(u-beta), so its coefficients -- and hence
          disc(q) -- are POLYNOMIALS in beta:

              disc(q) = -16 * P6(beta),
              P6(B) = B^6-126B^5+5271B^4-82564B^3+570591B^2-5779998B-9458775.

          The kill: take p with (beta|p) = -1 and Frobenius a 3-cycle on q,
          so q has no root mod p and beta is a non-residue -- no root of g is
          a square, so f_c has no root. Such p exist with positive density
          UNLESS Q(sqrt beta) is the quadratic subfield of q's splitting
          field, i.e. unless beta*disc(q) is a square. Since -16 is a square
          times -1, that says

              y^2 = -beta * P6(beta),

          a degree-7 squarefree hyperelliptic curve of GENUS 3. By Siegel it
          has finitely many INTEGRAL points, and its integral points are
          beta in {0, 1, 9, 25, 49} -- the roots of R, plus 0. Every one is
          either a perfect square (excluded: t^2 = beta would give f_c a
          rational root) or gives c = 0. So NO dangerous beta exists and the
          Chebotarev kill always applies.

THE GAP IS CLOSED -- 2026-08-23, scripts/k8_case2.py: the genus-3 curve is
SOLVED by an elementary descent (P6(0) = -3^7*5^2*173 splits the integral
points into 16 classes; compactness, a Runge squeeze and congruences -- two
of them at 2-adic depth 13 and 6 -- kill every class), so the integral
points are exactly the degenerate ones and Q27 is a THEOREM. This script's
search is retained as the historical record and cross-check.

    python scripts/k8_intersective.py
    python scripts/k8_intersective.py --beta_max 4000 --curve_max 20000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

x, u, t, c, B = sp.symbols("x u t c B")
BASE = sp.prod([x - i for i in range(8)])
RAD8 = 210  # rad(8!) = 2*3*5*7
P6 = (
    B**6 - 126 * B**5 + 5271 * B**4 - 82564 * B**3
    + 570591 * B**2 - 5779998 * B - 9458775
)


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: python -O strips assert."""
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def issq(v) -> bool:
    v = int(v)
    return v >= 0 and sp.integer_nthroot(v, 2)[1]


def R(z):
    return (z - 1) * (z - 9) * (z - 25) * (z - 49)


def step0() -> None:
    print("  [0] reduction and equivalence")
    g = sp.expand(sp.expand(BASE.subs(x, (t + 7) / 2) - c) * 2**8)
    Ru = sp.expand(sp.prod([u - (2 * j - 1) ** 2 for j in (1, 2, 3, 4)]))
    check(sp.simplify(g - (Ru.subs(u, t**2) - 256 * c)) == 0, "2^8 f_c = R(t^2)-256c")
    check(
        sp.expand(Ru) == u**4 - 84 * u**3 + 1974 * u**2 - 12916 * u + 11025, "R(u)"
    )
    bad = 0
    for cc in range(-300, 301, 11):
        for p in (11, 13, 17, 19, 23, 29, 31, 37, 41):
            lhs = any(int(sp.prod([v - i for i in range(8)]) - cc) % p == 0 for v in range(p))
            sq = {(v * v) % p for v in range(p)}
            gv = [
                v for v in range(p)
                if (v**4 - 84 * v**3 + 1974 * v**2 - 12916 * v + (11025 - 256 * cc)) % p == 0
            ]
            if lhs != any(v in sq for v in gv):
                bad += 1
    check(bad == 0, f"equivalence failed {bad} times")
    print("      g(u) = R(u)-256c is a QUARTIC; equivalence checked, 0 failures")


def cases_1_and_4() -> None:
    print("  [1,4] g irreducible / two irreducible quadratics -- no computation")
    print("      irreducible: transitive -> Jordan -> derangement -> Chebotarev")
    print("      two quadratics: both discs non-square, so a positive density of p")
    print("        has both non-residues and g has no root at all. Never intersective.")


def case_5_and_3() -> list:
    print("  [5,3] g splits / exactly two rational roots -- both BOUNDED")
    out5, out3, lo, hi = [], [], None, None
    for sv in range(-3000, 3001):
        if sv == 42:
            continue
        num, den = -(sv**3) + 84 * sv**2 - 1974 * sv + 12916, 84 - 2 * sv
        if num % den:
            continue
        pv = num // den
        D = sv * sv - 4 * pv
        if D < 0:
            continue
        lo = sv if lo is None else lo
        hi = sv
        r, ex = sp.integer_nthroot(D, 2)
        if not ex or (sv + r) % 2:
            continue
        b1, b2 = (sv + r) // 2, (sv - r) // 2
        rest = 84 - sv
        b34p = 1974 - pv - sv * rest
        if (11025 - pv * b34p) % 256:
            continue
        cv = (11025 - pv * b34p) // 256
        Dq = rest * rest - 4 * b34p
        (out5 if (Dq >= 0 and issq(Dq)) else out3).append((sorted([b1, b2]), cv))
    check(lo is not None and (lo, hi) == (10, 74), f"s range should be [10,74], got {lo},{hi}")
    check({cv for _, cv in out5} == {0}, f"case 5 must give c=0 only, got {out5}")
    check(not out3, f"case 3 proper must be empty, got {out3}")
    print(f"      feasible s range [{lo},{hi}]; the true bound is integrality:")
    print("      p(s) integral only for s even with (s-42) | 2048 (22 values)")
    print(f"      case 5: only betas {{1,9,25,49}}, c = 0 (which has rational roots)")
    print("      case 3: EMPTY -- every solution has a reducible remainder")
    return [{"c": 0, "case": 5}]


def case_2_curve(curve_max: int) -> dict:
    print("  [2] g = (u-beta) x irreducible cubic -- the infinite family")
    qq = sp.Poly(sp.expand(sp.cancel((R(u) - R(B)) / (u - B))), u)
    D = sp.factor(sp.discriminant(qq))
    check(sp.simplify(D - (-16 * P6)) == 0, "disc(q) = -16*P6(beta)")
    check(sp.Poly(P6, B).is_irreducible, "P6 irreducible")
    check(P6.subs(B, 0) != 0, "P6 coprime to B")
    F = sp.Poly(sp.expand(-B * P6), B)
    check(F.degree() == 7, "curve has degree 7")
    check(sp.gcd(F, F.diff(B)) == 1, "curve squarefree")
    genus = (F.degree() - 1) // 2
    check(genus == 3, f"genus should be 3, got {genus}")
    pts = []
    Fe = F.as_expr()
    for bv in range(-curve_max, curve_max + 1):
        val = int(Fe.subs(B, bv))
        if val < 0:
            continue
        r, ex = sp.integer_nthroot(val, 2)
        if ex:
            pts.append(bv)
    check(set(pts) == {0, 1, 9, 25, 49}, f"unexpected dangerous beta: {pts}")
    print(f"      disc(q) = -16*P6(beta) -- a POLYNOMIAL in beta (q is a divided difference)")
    print(f"      dangerous <=> y^2 = -beta*P6(beta): squarefree, degree 7, GENUS {genus}")
    print(f"      Siegel => finitely many integral points; found {pts}")
    print("      all are roots of R (c=0) or perfect squares (excluded). NO dangerous beta.")
    return {"genus": genus, "integral_points": pts, "curve_max": curve_max}


def case_2_verify(beta_max: int) -> dict:
    print(f"  [2b] direct check of the family, |beta| <= {beta_max}")
    n = dang = worst = 0
    for bv in range(-beta_max, beta_max + 1):
        r = R(bv)
        if r % 256:
            continue
        cv = r // 256
        if cv == 0 or issq(bv) or cv % RAD8:
            continue
        g = sp.Poly(u**4 - 84 * u**3 + 1974 * u**2 - 12916 * u + (11025 - 256 * cv), u)
        fl = g.factor_list()[1]
        if tuple(sorted(pp.degree() for pp, _ in fl)) != (1, 3):
            continue
        n += 1
        q = [pp for pp, _ in fl if pp.degree() == 3][0]
        Dq = int(sp.discriminant(q))
        if issq(bv * Dq):
            dang += 1
        kp = None
        for p in sp.primerange(5, 5000):
            if all(
                not any(int(f.eval(v)) % p == 0 for v in range(p))
                for f, _ in sp.Poly(BASE - cv, x).factor_list()[1]
            ):
                kp = p
                break
        check(kp is not None, f"beta={bv}: no killing prime")
        worst = max(worst, kp)
    check(dang == 0, f"{dang} dangerous beta found in range")
    print(f"      candidates: {n}   dangerous: {dang}   largest killing prime: {worst}")
    return {"candidates": n, "dangerous": dang, "max_prime": worst}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--beta_max", type=int, default=1200)
    ap.add_argument("--curve_max", type=int, default=20000)
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("=== Q27: (x)_8 - c is never intersective (historical run; gap closed 2026-08-23) ===")
    step0()
    cases_1_and_4()
    r53 = case_5_and_3()
    rc = case_2_curve(args.curve_max)
    r2 = case_2_verify(args.beta_max)
    print()
    print("  RESULT: all five cases closed.")
    print("    1,4  no computation needed -- Jordan / both-discs-non-residue")
    print("    3,5  bounded (s in [10,74]); case 3 empty, case 5 only c=0")
    print("    2    dangerous locus is a genus-3 curve whose integral points are")
    print("         exactly the roots of R, all excluded. Chebotarev always applies.")
    print("  GAP:  closed 2026-08-23 (scripts/k8_case2.py): the curve is solved by")
    print("        descent; the points above are ALL of them. Q27 is a theorem.")

    if args.json_out:
        # FROZEN TEXT: results/k8_intersective.json is sha256-pinned by
        # test_k8_case2.py and must regenerate byte-identically; the closure
        # lives in scripts/k8_case2.py. Do not "update" these strings.
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "search": "k8_intersective",
                    "claim": "(x)_8 - c is never intersective, modulo effective Siegel",
                    "cases_1_4": "closed with no computation",
                    "cases_3_5": r53,
                    "case_2_curve": rc,
                    "case_2_direct": r2,
                    "intersective_found": 0,
                    "gap": "effective integral points on the genus-3 curve not computed",
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
