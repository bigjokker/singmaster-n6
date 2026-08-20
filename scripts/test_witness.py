#!/usr/bin/env python3
"""Regression tests for the certificate layer.

The verifier is the only thing standing between "the program printed True" and
"the theorem is checkable", so most of these are adversarial: they feed it
certificates that must be REJECTED. A checker that always says VALID would
pass a test suite made only of good inputs.

Ground truth is exact integer arithmetic (math.comb), not any part of the
sweep. Run: python scripts/test_witness.py
"""

from __future__ import annotations

import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2  # noqa: E402

import witness as W  # noqa: E402
from bandii_kernel import make_fam  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


I_TEST = 3  # small enough to brute-force every image exactly


def _fam():
    fam = make_fam(I_TEST)
    return fam, int(gmpy2.bincoef(fam.N, fam.K))


def test_against_exact_arithmetic() -> None:
    """Every stored witness, checked against math.comb over the whole domain."""
    path = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not path.exists():
        errors.append(f"{path.name} missing; run family_sweep --i {I_TEST} first")
        return
    fam, m = _fam()
    ks, ps, meta = W.load(path)
    disagree = []
    for k, p in zip(ks.tolist(), ps.tolist()):
        in_image = any(math.comb(x, k) % p == m % p for x in range(p))
        claimed_kill = W.check_witness(fam.N, fam.K, k, p)["ok"]
        if claimed_kill == in_image:  # claimed_kill means "not in image"
            disagree.append((k, p))
    expect(
        not disagree,
        f"all {len(ks)} witnesses agree with brute-force exact arithmetic",
    )
    if disagree:
        errors.append(f"  first disagreements: {disagree[:5]}")


def test_rejects_bad_certificates() -> None:
    """The verifier must refuse anything that is not a genuine obstruction."""
    fam, m = _fam()
    N, K = fam.N, fam.K

    # a prime that does NOT obstruct column K (K is representable by construction)
    p_live = int(gmpy2.next_prime(K))
    while int(m % p_live) == 0:
        p_live = int(gmpy2.next_prime(p_live))
    expect(
        W.check_witness(N, K, K, p_live)["ok"] is False,
        "rejects a prime that does not obstruct the column",
    )

    expect(
        W.check_witness(N, K, 400, 359)["ok"] is False,
        "rejects p <= k (the image characterisation needs p > k)",
    )
    expect(W.check_witness(N, K, 1, 359)["ok"] is False, "rejects k < 2")

    # r(p) = 0 is always in the image, so a dead prime certifies nothing
    p_dead = 2
    while True:
        p_dead = int(gmpy2.next_prime(p_dead))
        if p_dead > K and int(m % p_dead) == 0:
            break
    res = W.check_witness(N, K, 5, p_dead)
    expect(
        res["ok"] is False and "certifies nothing" in res["reason"],
        f"rejects a dead prime r(p)=0 (p={p_dead})",
    )


def test_rejects_composite_p() -> None:
    """Job 2 verifier spec, item (i): p must be prime.

    Both legs of the argument assume it -- Lucas for r(p), and the claim that
    I_{p,k} is the complete image of column k. A composite modulus proves
    nothing. Most composites happen to crash or be rejected anyway, which is
    luck, not a check: the verifier has to establish primality itself.
    """
    fam, _ = _fam()
    bad = []
    composites = [12, 25, 91, 1001, 5401853 * 3, 10**6 + 3]
    for c in composites:
        if gmpy2.is_prime(c):
            continue
        try:
            res = W.check_witness(fam.N, fam.K, 5, c)
        except Exception as exc:
            bad.append((c, f"raised {type(exc).__name__} instead of refusing"))
            continue
        if res["ok"] or "not prime" not in res["reason"]:
            bad.append((c, res.get("reason")))
    expect(not bad, f"composite p is refused as not prime ({len(composites)} tried)")
    if bad:
        errors.append(f"  {bad[:3]}")


def test_primality_routine() -> None:
    """The verifier carries its own Miller-Rabin; it must be exact."""
    wrong = [n for n in range(2, 20000) if W.is_prime_pure(n) != bool(gmpy2.is_prime(n))]
    for x in (10**6, 10**7, 5401852, 37024873, 2**31):
        q = int(gmpy2.next_prime(x))
        if not W.is_prime_pure(q):
            wrong.append(q)
    expect(not wrong, "is_prime_pure agrees with gmpy2 on 0..20000 and large primes")


def test_r_routes_agree() -> None:
    """Item (ii): r must be recomputed by every route that applies, all agreeing."""
    from bandii_kernel import D, K, N

    m = int(gmpy2.bincoef(N, K))
    counts = {"two_digit": 0, "delta": 0}
    mis = []
    p = int(gmpy2.next_prime(N // 2))
    while p <= D and counts["delta"] < 25:
        exact = int(m % p)
        if W.lucas_mod_pure(N, K, p) != exact:
            mis.append(("lucas", p))
        for name, fn in (("two_digit", W.r_two_digit_pure), ("delta", W.r_delta_pure)):
            v = fn(N, K, p)
            if v is not None:
                counts[name] += 1
                if v != exact:
                    mis.append((name, p))
        p = int(gmpy2.next_prime(p))
    expect(
        not mis and counts["delta"] > 0 and counts["two_digit"] > 0,
        f"all r(p) routes equal m mod p on the live window {counts}",
    )
    expect(
        W.r_two_digit_pure(N, K, 3187) is None,
        "two-digit route declines below sqrt(N) instead of returning a false r",
    )


def test_certificate_is_bound_to_m() -> None:
    """The verifier must actually read N and K, not just (k, p).

    Note what this can NOT be tested by: at these g/p the kill rate is close
    to 1, so a certificate for m usually obstructs a different m' as well,
    purely by chance. "Still valid for m+1" is therefore not evidence of a
    bug. The binding that matters is that the r being tested really is
    m mod p, which is exact and checkable against gmpy2.
    """
    path = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not path.exists():
        return
    fam, m = _fam()
    ks, ps, _ = W.load(path)
    sample = ps[:: max(1, len(ps) // 50)].tolist()
    wrong = [p for p in sample if W.lucas_mod_pure(fam.N, fam.K, p) != int(m % p)]
    expect(
        not wrong,
        f"the r(p) the verifier tests is exactly m mod p ({len(sample)} primes, "
        f"checked against exact gmpy2)",
    )

    # and the verdict must genuinely depend on m for at least some columns
    pairs = list(zip(ks.tolist(), ps.tolist()))
    flipped = sum(
        1
        for k, p in pairs
        if W.check_witness(fam.N, fam.K, k, p)["ok"]
        != W.check_witness(fam.N + 1, fam.K, k, p)["ok"]
    )
    expect(
        flipped > 0,
        f"the verdict depends on (N,K): {flipped}/{len(pairs)} certificates "
        f"change verdict for a different m",
    )


def test_coverage_detects_a_hole() -> None:
    """Dropping one column must make the table incomplete."""
    path = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not path.exists():
        return
    ks, ps, meta = W.load(path)
    full = W.coverage(ks, meta)
    expect(full["complete"], f"stored table is complete ({full['n_expected']} columns)")

    holed = ks[1:]
    hole = W.coverage(holed, meta)
    expect(
        not hole["complete"] and hole["n_missing"] == 1,
        f"a single dropped witness is detected (missing={hole['n_missing']})",
    )

    import numpy as np

    extra = np.append(ks, np.int64(10**9))
    exc = W.coverage(extra, meta)
    expect(
        not exc["complete"] and exc["n_extra"] == 1,
        "a witness outside the claimed range is detected",
    )


def test_tamper_detection() -> None:
    """The stored sha256 must catch an edited witness table."""
    src = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not src.exists():
        return
    import numpy as np

    with tempfile.TemporaryDirectory() as td:
        ks, ps, meta = W.load(src)
        bad = Path(td) / "tampered.npz"
        ps2 = ps.copy()
        ps2[0] = int(ps2[0]) + 2  # swap in a different prime, keep the old hash
        np.savez_compressed(bad, k=ks, p=ps2, meta=json.dumps(meta))
        try:
            W.load(bad)
            expect(False, "tampered witness table was accepted")
        except RuntimeError:
            expect(True, "tampered witness table is rejected by its own sha256")


def test_checkpoint_formats_agree() -> None:
    """family_sweep, bandii_sweep and zjump label passes differently.

    The three writers tag chunks as "bii3"/"z2", prime_index, and round. The
    builder must reconstruct the same witnesses from all three shapes.
    """
    chk = ROOT / "results" / f"i{I_TEST}_sweep.jsonl"
    if not chk.exists():
        errors.append(f"{chk.name} missing; cannot test checkpoint adapters")
        return
    rows = W.read_jsonl(chk)

    bii = W._pass_tag("bii")
    z = W._pass_tag("z")
    tagged = [r for r in rows if str(r.get("tag", "")).startswith("bii")]
    expect(bool(tagged), "found tagged Band II chunk records to relabel")

    relabelled = []
    for r in tagged:
        r2 = dict(r)
        r2["prime_index"] = int(r["tag"][3:])
        r2.pop("tag")
        relabelled.append(r2)
    same = all(
        bii(a) == bii(b) for a, b in zip(tagged, relabelled)
    )
    expect(same, "bandii_sweep 'prime_index' records resolve to the same pass")

    zt = [r for r in rows if str(r.get("tag", "")).startswith("z")]
    zrel = []
    for r in zt:
        r2 = dict(r)
        r2["round"] = int(r["tag"][1:])
        r2.pop("tag")
        zrel.append(r2)
    expect(
        all(z(a) == z(b) for a, b in zip(zt, zrel)),
        "zjump 'round' records resolve to the same round",
    )
    expect(
        bii({"tag": "z3", "p": 5}) is None and z({"tag": "bii3", "p": 5}) is None,
        "the two phase adapters do not read each other's records",
    )


def test_rebuild_is_deterministic() -> None:
    """Rebuilding from the same checkpoint must give byte-identical witnesses."""
    chk = ROOT / "results" / f"i{I_TEST}_sweep.jsonl"
    src = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not (chk.exists() and src.exists()):
        return
    with tempfile.TemporaryDirectory() as td:
        again = Path(td) / "again.npz"
        meta = W._build_family(I_TEST, chk, again)
        _, _, orig = W.load(src)
        expect(
            meta["sha256"] == orig["sha256"],
            "rebuilding from the checkpoint reproduces the same sha256",
        )


def main() -> int:
    test_against_exact_arithmetic()
    test_rejects_bad_certificates()
    test_rejects_composite_p()
    test_primality_routine()
    test_r_routes_agree()
    test_certificate_is_bound_to_m()
    test_coverage_detects_a_hole()
    test_tamper_detection()
    test_checkpoint_formats_agree()
    test_rebuild_is_deterministic()
    print("\n=== WITNESS TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
