#!/usr/bin/env python3
"""Regression tests for the Q28 wall audit (the Q27 battery-depth lesson applied).

Pins:
  1. the calibration: Q27's delta = 3 class is ALIVE mod 2^7 and DEAD mod
     2^14 through this harness (v2 = 13 exactly for odd m) -- proving the
     2^7 cap was the bug at k=8 and that the harness can catch such kills;
  2. the battery is genuinely deep: the p=2 branch witnesses run at
     precision 2^64-headroom (>= 2^40 resolved), far past any 2^7 cap;
  3. the fast lift-tree (linear Hensel lift) produces IDENTICAL survivor
     sets to digit brute force at four (p, depth) points -- the deep table
     is not an artifact of the lifting method;
  4. mid-depth tree counts (raw tuples + resolved-valuation counts) are
     pinned to live computation, and DEEP_RECORD carries the same shape
     one-off deep runs, re-derivable with `python scripts/k10_deep.py --deep`;
  5. the die/survive table: the k=10 locus SURVIVES at p = 2, 3, 5, 7 and
     the spot-primes 11, 13, 31 -- exact branch witnesses obey the linear
     laws v2(cnum) = t+14, v3 = v5 = v7 = t+1, v11 = v13 = v31 = t, with
     both side square-conditions CERTIFIED (sqclass_decisive: resolved even
     valuation and resolved square unit -- not just consistent);
  6. the prune trap, both halves: step4's pruned tree really reports the
     false death at depth 2^19 AND the exact branch carries a certified
     point (v2 = 26) at the same offset a = -10 + 2^12;
  7. the 160 degenerate c = 0 tuples (ten a-values) satisfy the system
     exactly, and both witness bases are among them;
  8. the recorded candidates c = 1395418752000 and 2235340800 are OFF the
     locus (D_q*D_k nonsquare) and still die at 11 and 13 -- no false
     2-adic kill is claimed for them;
  9. WALL_STANDS is True (p-adics do not close the locus). Magma, not this
     battery, set PROVED_23_POSITIVE. results/k10_intersective.json is
     byte-identical (sha256-pinned, historical "modulo" claim), and
     results/k10_deep.json's verdict is pinned in FULL -- a prefix-only pin
     once let the artifact go stale against its own generator.

Fails before scripts/k10_deep.py exists.  Runs in ~5 s.

Run: python scripts/test_k10_deep.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

# sha256 of results/k10_intersective.json as committed before this lab:
K10_ARTIFACT_SHA = "8bc906a427faf963d3441164e007d0d142bd4045b7ec0733a86ff3681d4270cf"

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


try:
    import k10_deep as kd
except ImportError:
    kd = None
import k10_chebotarev as kc  # noqa: E402


def test_calibration() -> None:
    expect(kd is not None, "scripts/k10_deep.py exists and imports")
    if kd is None:
        return
    for k, expect_alive in ((7, True), (14, False)):
        N = 2**k
        sq = {(n_*n_) % N for n_ in range(N)}
        alive = any(kd.q27_delta3_value(mm) % N in sq for mm in range(N))
        expect(alive == expect_alive,
               f"Q27 delta=3 through the harness: alive mod 2^{k} is {expect_alive}")
    expect({kd.v_p(kd.q27_delta3_value(mm), 2) for mm in range(1, 31, 2)} == {13},
           "the Q27 kill was the exact odd valuation 13 -- the cap was the bug")


def test_degenerate_tuples() -> None:
    if kd is None:
        return
    deg = kd.degenerate_tuples()
    expect(len(deg) == 160, "160 degenerate c = 0 tuples")
    expect(len({t[0] for t in deg}) == 10, "over ten a-values (3b's ten c = 0 points)")
    expect(all(all(v == 0 for v in kd._eq_vals(*t)) for t in deg),
           "every degenerate tuple satisfies the system exactly over Z")
    expect(kd.WITNESS_BASE in deg, "the witness base is a degenerate tuple")
    expect((-10, -3, 315, -7, -53) in deg, "so is the prune-trap base")


def test_sqclass_decisive() -> None:
    if kd is None:
        return
    expect(kd.sqclass_decisive(2**4 * 17**2, 2, 30) is True,
           "decisive: 2^4*17^2 certified a Z_2-square")
    expect(kd.sqclass_decisive(0, 2, 30) is False,
           "decisive: zero-to-precision is NOT certified (sqclass_ok would pass it)")
    expect(kd.sqclass_ok(0, 2, 30) is True,
           "necessary-condition check keeps its permissive semantics")
    expect(kd.sqclass_decisive(2**3, 2, 30) is False, "decisive: odd valuation fails")
    expect(kd.sqclass_decisive(2**28, 2, 30) is False,
           "decisive: unresolved unit head fails even with even valuation")
    expect(kd.sqclass_decisive(5**2 * 6, 5, 10) is True,
           "decisive: 5^2*6 certified a Z_5-square (6 = 1 mod 5, a residue)")
    expect(kd.sqclass_decisive(2, 5, 10) is False, "decisive: non-residue unit fails")


def test_tree_fast_equals_slow() -> None:
    if kd is None:
        return
    for p, k in ((2, 6), (3, 4), (5, 3), (7, 2)):
        fast = kd._tree_fast(p, k, collect=True)[2]
        slow = kd._tree_slow(p, k)
        expect(fast == slow and len(fast) > 0,
               f"p={p}, depth {k}: linear Hensel lift == digit brute force "
               f"({len(fast)} tuples)")


def test_tree_alive() -> None:
    if kd is None:
        return
    for p, k, want in ((2, 8, (245760, 0)), (3, 5, (76324, 4)),
                       (5, 4, (132900, 119620)), (7, 3, (12320, 6160))):
        leaves, frozen, _ = kd._tree_fast(p, k)
        expect((leaves, frozen) == want,
               f"p={p}, depth {k}: {want[0]} raw tuples, {want[1]} with resolved "
               f"v_{p}(cnum) -- live computation matches the pin")
    expect(kd.DEEP_RECORD == {2: (10, 6881280, 0), 3: (9, 6181920, 5724000),
                              5: (6, 3322500, 3309220), 7: (5, 603680, 597520)},
           "the deep table (raw tuples; re-derive with k10_deep.py --deep)")
    expect(kd.DEEP_RECORD[2][2] == 0,
           "p=2 resolved-count 0 is structural: the filter forces cnum = 0 mod 2^k "
           "for k <= 18, so v2 resolution needs the witnesses, not the tree")


def test_witness_laws() -> None:
    if kd is None:
        return
    for p, K, ts in ((2, 64, (14, 20)), (3, 40, (6, 9)), (5, 34, (4, 7)),
                     (7, 30, (3, 6)), (11, 20, (2, 4)), (13, 18, (2, 4)),
                     (31, 16, (1, 3))):
        for t in ts:
            r = kd.branch_witness(kd.WITNESS_BASE, 1, t, p, K)
            expect(isinstance(r, dict), f"p={p}, t={t}: branch witness exists")
            if isinstance(r, dict):
                expect(r["v"] == t + kd.WITNESS_LAW[p],
                       f"p={p}, t={t}: v_{p}(cnum) = t + {kd.WITNESS_LAW[p]} exactly")
                expect(r["T1"] and r["T2"],
                       f"p={p}, t={t}: side squares certified (decisive)")
                expect(r["v"] >= kd.REQ_V.get(p, 0), f"p={p}, t={t}: inside the 10!-window")
    r = kd.branch_witness(kd.WITNESS_BASE, 1, 20, 2, 64)
    expect(isinstance(r, dict) and r["prec"] >= 40,
           "the 2-adic battery resolves past 2^40 -- not capped at 2^7")


def test_prune_trap() -> None:
    if kd is None:
        return
    died_at = kd.step4()
    expect(died_at == 19,
           "step4 runs live: the pruned tree really reports the false death at 2^19")
    r = kd.branch_witness((-10, -3, 315, -7, -53), 1, 12, 2, 64)
    expect(isinstance(r, dict) and r["v"] == 26,
           "the offset the pruned tree called dead carries a certified point (v2 = 26)")


def test_candidates_and_status() -> None:
    if kd is not None:
        expect(kd.WALL_STANDS is True, "the verdict: the wall stands")
    expect(kc.PROVED_23_POSITIVE is True,
           "PROVED_23_POSITIVE True because Magma listed the 15 a-values, not this battery")
    raw = (ROOT / "results" / "k10_intersective.json").read_bytes()
    expect(hashlib.sha256(raw).hexdigest() == K10_ARTIFACT_SHA,
           "results/k10_intersective.json byte-identical (sha256-pinned)")
    art = json.loads(raw.decode("utf-8"))
    expect("modulo" in art["claim"], "historical k10 artifact claim unchanged")
    u, t_ = sp.symbols("u t")
    R5 = sp.expand(sp.prod([u - (2*j - 1)**2 for j in range(1, 6)]))
    for cv, kp in ((1395418752000, 11), (2235340800, 13)):
        P_ = 1
        for pp, _ in sp.Poly(R5 - 1024*cv, u).factor_list()[1]:
            P_ *= int(sp.discriminant(pp.as_poly(u), u))
        expect(P_ > 0 and not sp.integer_nthroot(P_, 2)[1],
               f"c={cv}: off the locus (D_q*D_k nonsquare)")
        co = [int(x) for x in sp.Poly(sp.expand((R5 - 1024*cv).subs(u, t_*t_)), t_).all_coeffs()]
        killed = None
        for p_ in sp.primerange(3, 40):
            if not any(sum(c2 * pow(w, i, p_) for i, c2 in enumerate(reversed(co))) % p_ == 0
                       for w in range(p_)):
                killed = int(p_)
                break
        expect(killed == kp, f"c={cv}: still killed at {kp} -- no false 2-adic kill")
    newp = ROOT / "results" / "k10_deep.json"
    expect(newp.exists(), "results/k10_deep.json exists")
    if newp.exists() and kd is not None:
        art2 = json.loads(newp.read_text(encoding="utf-8"))
        expect(art2["verdict"] == ("SURVIVE at 2, 3, 5, 7 (spot: 11, 13, 31); "
                                   "p-adic wall stands; Magma closed the list 2026-08-23"),
               "artifact verdict is the FULL current string (a prefix pin let a "
               "stale 'PROVED_23_POSITIVE unchanged (False)' tail survive once)")
        expect("False" not in art2["verdict"],
               "artifact verdict does not contradict the repo's proof state")
        expect(art2["witness_laws"]["2"] == "v_2(cnum) = t + 14", "the 2-adic law recorded")
        expect(art2["witness_laws"]["11"] == "v_11(cnum) = t + 0", "the spot-law recorded")
        expect(art2["tree_alive"] == {str(p): {"depth": d, "raw_tuples": n, "v_resolved": f}
                                      for p, (d, n, f) in kd.DEEP_RECORD.items()},
               "artifact deep table matches DEEP_RECORD")


def main() -> int:
    test_calibration()
    test_degenerate_tuples()
    test_sqclass_decisive()
    test_tree_fast_equals_slow()
    test_tree_alive()
    test_witness_laws()
    test_prune_trap()
    test_candidates_and_status()
    print("\n=== K10 DEEP-BATTERY TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
