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

import hashlib
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


def test_ghost_census_claim_holds() -> None:
    """The census claim, checked exactly on real certificates.

    "c = k! C(N,K) is outside (x)_k(F_p)" must be literally true for the
    recorded pairs, computed from the definition with no repo arithmetic.
    """
    path = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not path.exists():
        return
    fam, m = _fam()
    ks, ps, _ = W.load(path)
    bad = []
    for k, p in list(zip(ks.tolist(), ps.tolist()))[:12]:
        c = math.factorial(k) * m
        # (x)_k == c mod p for some x?  falling factorial, straight from
        # the definition -- no image walk, no factorial table
        hit = any(
            (math.prod(range(x - k + 1, x + 1)) - c) % p == 0 for x in range(p)
        )
        if hit:
            bad.append((k, p))
    expect(not bad, "ghost-census claim holds exactly on sampled certificates")


def test_independent_verifier_agrees() -> None:
    """The sympy/brute-force cross-check must agree on real certificates."""
    path = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not path.exists():
        return
    try:
        import verify_independent as V
    except ImportError as exc:
        errors.append(f"verify_independent not importable: {exc}")
        return
    fam, m = _fam()
    ks, ps, _ = W.load(path)
    bad = []
    for k, p in list(zip(ks.tolist(), ps.tolist()))[:10]:
        r = V.sympy_lucas(fam.N, fam.K, p)
        if r != int(m % p):
            bad.append((k, p, "sympy Lucas disagrees with gmpy2"))
        elif V.brute_has_root(k, r, p) or V.factorial_has_root(k, r, p):
            bad.append((k, p, "an independent route finds a root"))
    expect(not bad, f"independent routes agree on sampled certificates {bad[:2]}")


def test_fill_leaves_no_stale_label() -> None:
    """The label must describe the file it is attached to.

    `fill` carried the ORIGINAL build's n_unresolved/unresolved forward through
    `{**meta, ...}` without rewriting them, so a table that had since been
    completed kept advertising the very columns it now held. i=9 shipped to two
    trees -- including the published one -- claiming n_unresolved=4 for
    k=87/399/553/1281 while holding all four at p=191/421/557/1321.

    A referee reading that caption would conclude the theorem had a hole.

    Both paths are pinned: the early return for an already-complete table, and
    the normal path, which must write what is ACTUALLY still missing rather
    than inherit a stale value. The arrays must never move -- only the caption.
    """
    import json

    import numpy as np

    i = 3
    tmp = Path(tempfile.mkdtemp(prefix="label_test_"))
    try:
        path = tmp / f"i{i}_witness.npz"
        src = ROOT / "results" / f"i{i}_witness.npz"
        if not src.exists():
            ok.append("i3_witness.npz absent; stale-label test skipped")
            return
        ks, ps, meta = W.load(src)
        before = hashlib.sha256(ks.tobytes() + ps.tobytes()).hexdigest()

        # a COMPLETE table wearing a stale caption -- exactly i=9's shape
        liar = {**meta, "n_unresolved": 4, "unresolved": [87, 399, 553, 1281]}
        W.save(path, liar, {int(a): int(b) for a, b in zip(ks, ps)})
        z = np.load(path, allow_pickle=False)
        expect(json.loads(str(z["meta"]))["n_unresolved"] == 4,
               "stale-label fixture really does carry the false caption")

        res = W.fill_small_gaps(path, i)
        ks2, ps2, meta2 = W.load(path)
        after = hashlib.sha256(ks2.tobytes() + ps2.tobytes()).hexdigest()

        expect(meta2.get("n_unresolved") == 0 and meta2.get("unresolved") == [],
               f"fill rewrites a stale label on a complete table "
               f"(n_unresolved={meta2.get('n_unresolved')}, "
               f"unresolved={meta2.get('unresolved')})")
        expect(after == before,
               "fill rewrote the LABEL and not the TABLE (array sha256 unchanged)")
        expect(np.array_equal(ks, ks2) and np.array_equal(ps, ps2),
               "every witness row survives a relabel byte-for-byte")
        expect(meta2.get("sha256") == meta.get("sha256"),
               "the array digest in meta is untouched, so the file's own "
               "self-check still passes")
        expect(bool(res.get("relabelled")),
               f"the relabel is REPORTED, not silent (got {res.get('relabelled')!r})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_repair_replaces_a_false_certificate() -> None:
    """A table can carry a certificate that is simply FALSE, and repair must
    fix exactly that row -- or refuse.

    This is not hypothetical. i=8 shipped with k=1021 recorded as killed by
    p=3517, a prime it merely SURVIVED: the Z-jump's extended small-k rounds
    drop survivors with k >= SMALL_K instead of killing them
    (family_sweep.py:365-368), and the builder reads "absent from this round's
    survivors" as "killed by this round's prime". Sampled verification
    (5,000 of 5,182,634) never looked at it.

    Two properties are pinned, and the second is the one that matters: repair
    must RAISE rather than write out a table it could not actually fix.
    """
    import json

    import numpy as np

    i = 3
    fam = make_fam(i)
    tmp = Path(tempfile.mkdtemp(prefix="repair_test_"))
    try:
        src = ROOT / "results" / f"i{i}_witness.npz"
        if not src.exists():
            ok.append(f"i{i}_witness.npz absent; repair test skipped")
            return
        ks, ps, meta = W.load(src)
        path = tmp / f"i{i}_witness.npz"

        # corrupt ONE row to a prime that does not obstruct its column
        have = {int(a): int(b) for a, b in zip(ks, ps)}
        victim = int(ks[len(ks) // 2])
        good_p = have[victim]
        bad_p = next(q for q in range(victim + 1, victim + 4000)
                     if W.is_prime_pure(q)
                     and not W.check_witness(fam.N, fam.K, victim, q)["ok"])
        have[victim] = bad_p
        W.save(path, meta, have)
        expect(not W.check_witness(fam.N, fam.K, victim, bad_p)["ok"],
               f"fixture really does carry a FALSE certificate (k={victim}, p={bad_p})")

        res = W.repair_invalid(path, i, [victim])
        ks2, ps2, meta2 = W.load(path)
        have2 = {int(a): int(b) for a, b in zip(ks2, ps2)}

        expect(len(res["repaired"]) == 1 and res["repaired"][0]["k"] == victim,
               f"repair fixes exactly the corrupted row ({res['repaired']})")
        expect(W.check_witness(fam.N, fam.K, victim, have2[victim])["ok"],
               f"the replacement certificate actually holds "
               f"(p={have2[victim]})")
        untouched = sum(1 for kk in have if kk != victim and have2.get(kk) == have[kk])
        expect(untouched == len(have) - 1,
               f"every other row is untouched ({untouched}/{len(have)-1})")
        expect(meta2.get("n_repaired") == 1 and meta2["repaired"][0]["old_p"] == bad_p,
               "the repair is recorded in meta, so provenance is not laundered")
        expect(np.array_equal(ks, ks2), "the column set itself is unchanged")

        # ADVERSARIAL: if no obstructing prime exists, refuse -- do not write a
        # table that still carries a known-false certificate.
        have[victim] = bad_p
        W.save(path, meta, have)
        before = path.read_bytes()
        import singmaster_intersect as SI
        real_op = SI.obstructing_prime
        SI.obstructing_prime = lambda *a, **kw: None
        try:
            W.repair_invalid(path, i, [victim])
            expect(False, "repair RAISES when a false certificate cannot be replaced")
        except RuntimeError as exc:
            expect("FALSE certificate" in str(exc),
                   f"repair raises rather than shipping an unrepaired table "
                   f"({str(exc)[:70]})")
        finally:
            SI.obstructing_prime = real_op
        expect(path.read_bytes() == before,
               "and it leaves the file byte-identical rather than half-written")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_ledger_catches_a_certificate_naming_another_table() -> None:
    """A certificate names its witness table BY DIGEST. Check the digest.

    Nothing did. Any post-sweep change to a table -- fill adding engine
    witnesses, repair replacing a false row -- left the certificate naming a
    table that no longer existed, and it was never noticed. Measured
    2026-08-22: five of seven certified members were stale, and i=8's named a
    digest matching no table that had ever been at that path.

    A referee following such a certificate verifies the wrong file, or one that
    is absent. Worse, if a stale table is still lying around it verifies
    cleanly and the reader concludes the claim is backed.
    """
    import json

    import coverage_ledger as CL

    i = 3
    src = ROOT / "results" / f"i{i}_witness.npz"
    sweep_src = ROOT / "results" / f"i{i}_sweep.json"
    if not (src.exists() and sweep_src.exists()):
        ok.append(f"i{i} artifacts absent; certificate-binding test skipped")
        return
    tmp = Path(tempfile.mkdtemp(prefix="certbind_test_"))
    try:
        shutil.copy(src, tmp / src.name)
        rep = json.loads(sweep_src.read_text(encoding="utf-8"))
        cert = rep.get("certificate") or ""
        if "sha256 " not in cert:
            ok.append(f"i{i} has no certificate digest; binding test skipped")
            return

        # honest copy first: the check must PASS on a correctly bound pair
        (tmp / sweep_src.name).write_text(json.dumps(rep), encoding="utf-8")
        good = CL.audit_member(i, tmp / src.name)
        expect(good["certificate_names_this_table"] is True and good["ok"],
               "ledger accepts a certificate that names the table on disk")

        # now point the certificate at a different table
        import re
        bad_cert = re.sub(r"sha256 [0-9a-f]+", "sha256 " + "0" * 16, cert)
        (tmp / sweep_src.name).write_text(
            json.dumps({**rep, "certificate": bad_cert}), encoding="utf-8")
        bad = CL.audit_member(i, tmp / src.name)
        expect(bad["certificate_names_this_table"] is False,
               f"ledger DETECTS a certificate naming another table "
               f"({bad['certificate_sha256']})")
        expect(bad["ok"] is False,
               "and a member whose certificate names another table is not COMPLETE")
        expect(bad["n_missing"] == 0 and bad["n_extra"] == 0,
               "coverage itself is still reported as intact -- the two failures "
               "are kept distinct rather than conflated")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_deferred_column_is_not_credited_a_kill() -> None:
    """A column the sweep never scanned must NOT get a witness.

    This is the k=1021 defect, reduced to its smallest form. The builders used
    to infer "absent from this round's survivors" => "killed by this round's
    prime". A DEFERRED column is absent too: past CAP_Z the sweep drops big-k
    columns from `current` entirely, so they stop appearing in records without
    ever having been tested.

    Nothing downstream caught it. The coverage ledger checks that every column
    HAS a witness, not that the witness is VALID, so it reported
    "missing 0, ok=True" over a certificate that check_witness rejects. i=8
    shipped k=1021 -> p=3517 that way, and re-minted it on every rebuild.

    Coverage decides killed-vs-untested where the records can prove it, and
    the count identity catches the residue a sparse span can hide. Both are
    exercised, because each covers the other's blind spot:
    the declared-alive list is exact but the sweep json truncates it above
    100, and the count identity needs no list but only fires when the records
    carry n_cols. And two non-regression cases, since a guard that also
    refuses legitimate builds is not a fix.
    """
    from bandii_kernel import cells, first_live_after, live_intervals

    fam = make_fam(5)
    ivs = live_intervals(cells(fam), fam)
    prev = first_live_after(300, ivs, fam.D)
    p = first_live_after(prev, ivs, fam.D)
    while p - prev < 2:
        prev, p = p, first_live_after(p, ivs, fam.D)
    a, b = prev, prev + 1          # same bucket => same next live prime
    tag = W._pass_tag("z")

    def rec(k_hi, n_cols, survivors):
        return {"tag": "z1", "p": p, "k_lo": a, "k_hi": k_hi,
                "n_cols": n_cols, "survivors": survivors,
                "seconds": 0.0, "n_survivors": len(survivors)}

    # b deferred and OUTSIDE every record's span: the round provably never
    # touched it, so it is handed to the engine rather than crashing the build
    w, alive = W.build_zjump([rec(a, 1, [])], a, b, ivs, fam.D, tag)
    expect(w == {a: p} and alive == {b},
           "a column no record covers is left unresolved, not credited "
           "and not a hard failure")

    # b deferred but INSIDE a sparse span -- coverage cannot separate them, so
    # the count identity is what catches it. This is the residue case, and the
    # reason both guards stay.
    sparse = rec(b, 1, [])          # span [a,b] but only 1 column scanned
    try:
        w, _ = W.build_zjump([sparse], a, b, ivs, fam.D, tag)
        expect(False, f"sparse-span deferral credited a kill: {w}")
    except RuntimeError as exc:
        expect("never scanned" in str(exc),
               "count identity catches a deferral hidden in a sparse span")

    # b deferred, and the sweep json lists it as still alive
    w, alive = W.build_zjump([rec(a, 1, [])], a, b, ivs, fam.D, tag,
                             declared_alive=frozenset({b}))
    expect(b not in w and b in alive and w.get(a) == p,
           "a column the run declares alive is left unresolved, not credited")

    # non-regression: both genuinely scanned and killed
    w, alive = W.build_zjump([rec(b, 2, [])], a, b, ivs, fam.D, tag)
    expect(w == {a: p, b: p} and not alive,
           "a genuine two-column kill still builds")

    # non-regression: a legacy record with no n_cols must not false-alarm
    legacy = rec(b, 2, [])
    del legacy["n_cols"]
    try:
        w, _ = W.build_zjump([legacy], a, b, ivs, fam.D, tag)
        expect(w == {a: p, b: p},
               "a legacy record without n_cols builds as before")
    except RuntimeError as exc:
        expect(False, f"count identity false-alarmed on a legacy record: {exc}")

    # and Band II carries the same guard
    brec = {"tag": "bii1", "p": p, "prime_index": 1, "k_lo": a, "k_hi": a,
            "n_cols": 1, "survivors": [], "seconds": 0.0, "n_survivors": 0}
    w, alive = W.build_bandii([brec], a, b, [p], W._pass_tag("bii"))
    expect(w == {a: p} and b in alive,
           "Band II leaves an uncovered column unresolved too")


def test_fill_and_repair_bind_to_the_table_directory() -> None:
    """fill/repair of a COPY must never touch the tracked sweep record.

    retarget_certificate used to default its target to ROOT/results/
    i{i}_sweep.json -- resolved from the family index, not from the table
    being edited. So fill or repair on a temp copy (which is exactly what the
    suites do) rewrote the TRACKED proof record whenever the copy's digest
    diverged from the shipped table, appending a `certificate_retargeted`
    entry describing a fill or repair that never happened to the shipped
    table, and pointing the certificate at a digest that existed only in a
    temp dir. The suites were a silent no-op only because their fixtures
    happened to reproduce the shipped bytes.

    The sweep record a table change may retarget is the one NEXT TO that
    table, or nothing. Run against a stand-in ROOT holding a copy of the
    tracked record, so a regression shows on the copy and can never reach
    results/. The real tracked files are digest-guarded in main() as well.
    """
    import json

    i = 3
    fam = make_fam(i)
    src = ROOT / "results" / f"i{i}_witness.npz"
    tracked = ROOT / "results" / f"i{i}_sweep.json"
    if not src.exists() or not tracked.exists():
        ok.append("i3 artifacts absent; directory-binding test skipped")
        return
    stand = Path(tempfile.mkdtemp(prefix="standin_root_"))
    tmp = Path(tempfile.mkdtemp(prefix="bind_test_"))
    real_root = W.ROOT
    try:
        (stand / "results").mkdir()
        decoy = stand / "results" / f"i{i}_sweep.json"
        shutil.copyfile(tracked, decoy)
        decoy_before = decoy.read_bytes()
        W.ROOT = stand

        ks, ps, meta = W.load(src)
        have = {int(a): int(b) for a, b in zip(ks, ps)}
        # Make the copy DIVERGE from the shipped table: one column gets a
        # different (still valid) witness, and one column is dropped so fill
        # has something to add. The rebuilt digest then differs from the one
        # the certificate names, which is the case that used to leak.
        victim = int(ks[len(ks) // 3])
        alt = next(q for q in range(have[victim] + 1, have[victim] + 5000)
                   if W.is_prime_pure(q) and W.check_witness(fam.N, fam.K, victim, q)["ok"])
        have[victim] = alt
        dropped = int(ks[len(ks) // 2])
        del have[dropped]
        path = tmp / f"i{i}_witness.npz"
        W.save(path, meta, have)

        res = W.fill_small_gaps(path, i)
        expect(res["added"] == 1 and not res["unresolved"],
               f"fixture: fill re-adds the dropped column (added {res['added']})")
        expect(res["sha256"] != meta["sha256"],
               "fixture: the filled copy really has a different digest from the shipped table")
        expect(decoy.read_bytes() == decoy_before,
               "fill of a temp copy leaves ROOT/results/i3_sweep.json byte-identical")
        expect(res.get("retargeted") is None,
               f"with no sweep record next to the table, nothing is retargeted "
               f"(got {res.get('retargeted')!r})")

        # Now put a sweep record NEXT TO the table: that one, and only that
        # one, is what a table change may retarget.
        sibling = tmp / f"i{i}_sweep.json"
        shutil.copyfile(tracked, sibling)
        ks2, ps2, meta2 = W.load(path)
        have2 = {int(a): int(b) for a, b in zip(ks2, ps2)}
        victim2 = int(ks2[len(ks2) // 4])
        bad_p = next(q for q in range(victim2 + 1, victim2 + 4000)
                     if W.is_prime_pure(q)
                     and not W.check_witness(fam.N, fam.K, victim2, q)["ok"])
        have2[victim2] = bad_p
        W.save(path, meta2, have2)
        rep = W.repair_invalid(path, i, [victim2])
        expect(len(rep["repaired"]) == 1, "fixture: repair replaced the corrupted row")
        sib = json.loads(sibling.read_text(encoding="utf-8"))
        expect(rep.get("retargeted") is not None
               and sib.get("certificate_retargeted")
               and rep["sha256"].startswith(sib["certificate_retargeted"][-1]["to"]),
               "repair retargets the sweep record NEXT TO the table it edited")
        expect(decoy.read_bytes() == decoy_before,
               "repair of a temp copy leaves ROOT/results/i3_sweep.json byte-identical")
    finally:
        W.ROOT = real_root
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.rmtree(stand, ignore_errors=True)


def test_build_family_clips_engine_band_at_kmax() -> None:
    """The engine band below the Z-jump start must stop at k_max.

    i=2 is the one member whose k_max (49) lies below K_EXACT (200). The
    builder searched obstructing primes for k=50..200 -- columns no claim
    contains -- wrote witnesses for them, claimed [2,200], and reported k=63
    and 65 (above k_max) as unresolved, so a rebuild exited 2 and a fresh
    sweep of i=2 would have withheld a certificate the member earns.

    The shipped results/i2_witness.npz (46 rows, digest ab08157a...) is the
    reference; a rebuild from the checkpoint must reproduce it exactly.
    """
    import coverage_ledger as CL
    from bandii_kernel import kmax_of

    i = 2
    chk = ROOT / "results" / f"i{i}_sweep.jsonl"
    src = ROOT / "results" / f"i{i}_witness.npz"
    if not chk.exists() or not src.exists():
        ok.append("i2 artifacts absent; k_max clip test skipped")
        return
    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    _, _, shipped = W.load(src)
    tmp = Path(tempfile.mkdtemp(prefix="clip_test_"))
    try:
        path = tmp / f"i{i}_witness.npz"
        meta = W._build_family(i, chk, path)
        ks, _ps, meta = W.load(path)
        expect(meta["n_unresolved"] == 0 and not meta["unresolved"],
               f"i=2 rebuild leaves nothing unresolved "
               f"(n_unresolved={meta['n_unresolved']}, {meta['unresolved']})")
        expect(int(ks.max()) <= kmax and ks.size == shipped["n_witnesses"],
               f"i=2 rebuild emits no column above k_max={kmax} "
               f"({ks.size} rows, max k {int(ks.max())})")
        expect(meta["sha256"] == shipped["sha256"],
               f"i=2 rebuild reproduces the shipped table "
               f"({meta['sha256'][:16]} vs {shipped['sha256'][:16]})")
        cov = W.coverage(ks, meta)
        expect(cov["complete"],
               f"rebuilt i=2 claim is complete and exact "
               f"(missing {cov['n_missing']}, extra {cov['n_extra']})")
        aud = CL.audit_member(i, path)
        expect(aud.get("ok") and aud["n_extra"] == 0 and aud["n_missing"] == 0,
               f"ledger accepts the rebuilt i=2 table (extra {aud.get('n_extra')})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_verifier_rejects_the_j0_image_entry() -> None:
    """r(p)=1 is C(k,k): the image entry at j=0, present for EVERY column.

    No adversarial case fed the verifier a prime with r(p)=1, so an image
    walk that started at j=1 would have accepted a false certificate for
    every column at that prime and passed every suite. Pin j=0 explicitly.
    """
    fam, _m = _fam()
    p = 787
    r = W.lucas_mod_pure(fam.N, fam.K, p)
    expect(r == 1, f"fixture: r({p}) = 1 for i={I_TEST} (got {r})")
    flips = []
    for k in range(2, 60):
        res = W.check_witness(fam.N, fam.K, k, p)
        if res["ok"] or res.get("j") != 0:
            flips.append((k, res))
    expect(not flips,
           f"check_witness rejects p={p} for every k<p with the hit at j=0 "
           f"({len(flips)} wrong verdicts {flips[:2]})")


def test_bandii_count_identity_catches_a_sparse_span() -> None:
    """Band II: a column absent from the survivors inside a SPARSE span is
    not a kill unless the record's n_cols says it was scanned.

    The Z-jump half of this identity is pinned in
    test_deferred_column_is_not_credited_a_kill; the Band II half was not,
    and removing it passed every suite. This is the k=1021 mechanism in the
    Band II phase.
    """
    p = 359
    sparse = {"tag": "bii1", "p": p, "prime_index": 1, "k_lo": 300, "k_hi": 301,
              "n_cols": 1, "survivors": [], "seconds": 0.0, "n_survivors": 0}
    try:
        w, alive = W.build_bandii([sparse], 300, 301, [p], W._pass_tag("bii"))
        expect(False, f"Band II sparse-span record credited {w}")
    except RuntimeError as exc:
        expect("never scanned" in str(exc),
               "Band II count identity refuses to credit a column the record "
               "did not scan")
    dense = {**sparse, "n_cols": 2}
    w, alive = W.build_bandii([dense], 300, 301, [p], W._pass_tag("bii"))
    expect(w == {300: p, 301: p} and not alive,
           "and a dense Band II record still credits both columns")


def test_verify_reports_invalid_rows() -> None:
    """witness.verify is the command the README hands a referee, and nothing
    called it. Pin: a table with one FALSE row is NOT VALID under a full
    check, the bad row is named, and a sampled check labels itself sampled.
    """
    fam, _m = _fam()
    src = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not src.exists():
        ok.append("i3_witness.npz absent; verify test skipped")
        return
    tmp = Path(tempfile.mkdtemp(prefix="verify_test_"))
    try:
        ks, ps, meta = W.load(src)
        have = {int(a): int(b) for a, b in zip(ks, ps)}
        victim = int(ks[len(ks) // 5])
        bad_p = next(q for q in range(victim + 1, victim + 4000)
                     if W.is_prime_pure(q)
                     and not W.check_witness(fam.N, fam.K, victim, q)["ok"])
        have[victim] = bad_p
        path = tmp / f"i{I_TEST}_witness.npz"
        W.save(path, meta, have)

        full = W.verify(path, sample=None, workers=1)
        expect(full["valid"] is False and full["n_invalid"] == 1
               and full["invalid"][0]["k"] == victim,
               f"full verify reports the one false certificate "
               f"(valid={full['valid']}, n_invalid={full['n_invalid']})")
        expect(full["coverage"]["complete"] and not full["sampled"],
               "full verify: coverage complete, not sampled")
        part = W.verify(path, sample=5, workers=1, seed=0)
        expect(part["sampled"] is True and part["n_checked"] == 5,
               "a sampled verify says so and checks exactly the sample")

        clean = W.verify(src, sample=20, workers=1, seed=1)
        expect(clean["valid"] is True and clean["n_invalid"] == 0,
               "the shipped table passes a sampled verify")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_ledger_recomputes_kmax_from_the_family() -> None:
    """coverage_ledger's headline property -- k_max from N,K, never from the
    file under audit -- had no test. A table that drops its top column AND
    shrinks its own claimed range to match must still be INCOMPLETE.
    """
    import coverage_ledger as CL
    from bandii_kernel import kmax_of

    fam, _m = _fam()
    src = ROOT / "results" / f"i{I_TEST}_witness.npz"
    if not src.exists():
        ok.append("i3_witness.npz absent; ledger k_max test skipped")
        return
    kmax, _ = kmax_of(fam)
    tmp = Path(tempfile.mkdtemp(prefix="ledger_kmax_"))
    try:
        ks, ps, meta = W.load(src)
        have = {int(a): int(b) for a, b in zip(ks, ps) if int(a) != kmax}
        liar = {**meta, "claimed_ranges": [[2, kmax - 1]], "k_max": kmax - 1,
                "excluded": [fam.K, fam.K + 1]}
        path = tmp / f"i{I_TEST}_witness.npz"
        W.save(path, liar, have)
        ks2, _ps2, meta2 = W.load(path)
        expect(W.coverage(ks2, meta2)["complete"],
               "fixture: the table's OWN claim is self-consistent")
        aud = CL.audit_member(I_TEST, path)
        expect(aud["k_max"] == kmax and aud["ok"] is False
               and aud["n_missing"] == 1 and aud["missing_sample"] == [kmax],
               f"ledger recomputes k_max={kmax} from N,K and reports the top "
               f"column missing (ok={aud['ok']}, k_max={aud['k_max']}, "
               f"missing={aud['missing_sample']})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# Tracked artifacts that no test in this file may modify. The fill/repair
# tests operate on temp copies; if a table change ever reaches one of these
# again, fail loudly rather than let git status be the only witness.
GUARDED = ("i3_sweep.json", "i5_sweep.json", "i3_witness.npz", "i5_witness.npz",
           "i2_witness.npz", "i2_sweep.json")


def _digest_guarded() -> dict:
    out = {}
    for name in GUARDED:
        p = ROOT / "results" / name
        out[name] = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    return out


def main() -> int:
    before = _digest_guarded()
    test_fill_and_repair_bind_to_the_table_directory()
    test_build_family_clips_engine_band_at_kmax()
    test_verifier_rejects_the_j0_image_entry()
    test_bandii_count_identity_catches_a_sparse_span()
    test_verify_reports_invalid_rows()
    test_ledger_recomputes_kmax_from_the_family()
    test_fill_leaves_no_stale_label()
    test_ledger_catches_a_certificate_naming_another_table()
    test_deferred_column_is_not_credited_a_kill()
    test_repair_replaces_a_false_certificate()
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
    test_ghost_census_claim_holds()
    test_independent_verifier_agrees()
    after = _digest_guarded()
    changed = [n for n in GUARDED if before[n] != after[n]]
    expect(not changed,
           "the suite left every tracked results artifact byte-identical"
           + (f" -- MUTATED: {changed}" if changed else ""))
    print("\n=== WITNESS TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
