#!/usr/bin/env python3
"""Regression tests for the sweep drivers.

singmaster_intersect.py's `sanity` covers the engine. Nothing covered the
sweep drivers, which is where the two "false clean certificate" bugs lived.
Run: python scripts/test_sweeps.py
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import bandii_kernel as bk  # noqa: E402
import family_sweep as fsw  # noqa: E402
import nextprime_sweep as nps  # noqa: E402

from singmaster_intersect import binom_mod_lucas, fib_member  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


# Artifacts the sweep tests used to overwrite in the REAL tree.
GUARDED = ("i3_sweep.json", "i3_sweep.jsonl", "i3_witness.npz")


def _digest_guarded() -> dict:
    """sha256 of every tracked artifact the sweep tests could touch."""
    out = {}
    for name in GUARDED:
        p = ROOT / "results" / name
        out[name] = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    return out


@contextlib.contextmanager
def sandboxed_results():
    """Point family_sweep's output at a throwaway results/ directory.

    Three of these tests call `fsw.main()`, which WRITES results/i{i}_sweep.json,
    its jsonl checkpoint, and -- on a clean run -- results/i{i}_witness.npz.
    Aimed at the real tree they deleted and regenerated tracked artifacts on
    every run. Worse, two of them unlink i3_sweep.json at the END, leaving it
    MISSING until the third regenerates it: a crash in between lost the file
    outright, and the proof object was being rewritten too.

    It went unnoticed because the sweep is deterministic, so the only field
    that ever changed was `seconds` -- a one-line diff. It stopped being
    cosmetic once work_census.py began gating against results/i{i}_sweep.json:
    the suite was mutating the ground truth it exists to protect.

    family_sweep resolves ROOT at call time in both `paths()` and the witness
    writer, so redirecting the module global covers both. Tests that read the
    REAL recorded runs use this module's own ROOT and are unaffected.
    """
    tmp = Path(tempfile.mkdtemp(prefix="sweep_test_"))
    (tmp / "results").mkdir()
    real = fsw.ROOT
    fsw.ROOT = tmp
    try:
        yield tmp
    finally:
        fsw.ROOT = real
        shutil.rmtree(tmp, ignore_errors=True)


def test_phase_count() -> None:
    """A phase_complete record is not proof the phase ended empty."""
    done = [
        {"event": "phase_complete", "phase": "bandii", "n_alive": 7},
        {"event": "phase_complete", "phase": "zjump", "n_alive": 0, "n_nolive": 3},
    ]
    expect(fsw._phase_count(done, "bandii", "n_alive") == 7,
           "_phase_count reads a non-zero survivor count back off the checkpoint")
    expect(fsw._phase_count(done, "zjump", "n_alive") == 0,
           "_phase_count returns 0 for a genuinely empty phase")
    expect(fsw._phase_count(done, "zjump", "n_nolive") == 3,
           "_phase_count recovers n_nolive across a resume")
    expect(fsw._phase_count([], "bandii", "n_alive") == 0,
           "_phase_count returns 0 when the phase never ran")


def test_resume_of_unfinished_phase_is_not_clean() -> None:
    """Resuming a phase that ended with survivors must not certify.

    Write a checkpoint whose Band II phase_complete records 5 survivors,
    then resume. Before the fix the resume branch hardcoded an empty list
    and the run emitted an unconditional certificate anyway.
    """
    i = 3
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    # the checkpoint must carry a matching schema header or resume refuses
    fam = bk.make_fam(i)
    kmax, _ = bk.kmax_of(fam)
    fsw.write_jsonl(chk, bk.checkpoint_identity(
        i=i, N=fam.N, K=fam.K, k_max=kmax,
        k_lo_z=fsw.K_EXACT.get(i, 2) + 1,
        cap_bii=fsw.CAP_BII, cap_z=fsw.CAP_Z,
    ))
    for rec in (
        {"event": "phase_complete", "phase": "bandii", "n_alive": 5},
        {"event": "phase_complete", "phase": "zjump", "n_alive": 0, "n_nolive": 2},
    ):
        fsw.write_jsonl(chk, rec)

    sys.argv = ["family_sweep.py", "--i", str(i)]
    fsw.main()
    rep = json.loads(out.read_text(encoding="utf-8"))
    expect(rep["n_bii_alive"] == 5,
           f"resumed Band II carries its recorded survivors (n_bii_alive={rep['n_bii_alive']})")
    expect(rep["n_z_nolive"] == 2,
           f"resumed Z-jump carries its recorded n_nolive ({rep['n_z_nolive']})")
    expect(rep["clean"] is False, "a resumed unfinished phase does not certify")
    expect(rep["certificate"] is None, "no certificate after resuming an unfinished phase")
    expect(rep["z_survivors"] == [] and rep["bii_survivors"] == [],
           "resume emits no fabricated survivor records")
    for f in (out, chk):
        f.unlink(missing_ok=True)


def test_nolive_is_not_clean() -> None:
    """A column with no live prime was never tested; it is not killed.

    Patch first_live_after so one Z-jump column can never be assigned a
    prime. Before the fix that column silently vanished from the survivor
    list and the run still emitted an unconditional certificate.
    """
    i = 3
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)

    real = bk.first_live_after
    doomed = 250

    def fake(x, ivs, d=bk.D):
        # k=250 starts its walk at x=250 and can never find a prime
        if x == doomed:
            return None
        return real(x, ivs, d)

    fsw.first_live_after = fake
    try:
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
    finally:
        fsw.first_live_after = real

    rep = json.loads(out.read_text(encoding="utf-8"))
    expect(rep["n_z_nolive"] == 1,
           f"untestable column is counted (n_z_nolive={rep['n_z_nolive']})")
    expect(doomed in rep["z_nolive"],
           f"the untestable column is named in z_nolive ({rep['z_nolive']})")
    expect(rep["clean"] is False,
           "clean is False when a column could not be tested")
    expect(rep["certificate"] is None,
           "no certificate is emitted when a column could not be tested")
    # the pre-fix predicate would have passed this run
    expect(rep["n_bii_alive"] == 0 and rep["n_z_alive"] == 0,
           "and the OLD predicate (survivors only) would have called this clean "
           "-- which is exactly the bug")

    for f in (out, chk):
        f.unlink(missing_ok=True)


def test_clean_run_still_clean() -> None:
    """The same member, unpatched, must still certify."""
    i = 3
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    sys.argv = ["family_sweep.py", "--i", str(i)]
    fsw.main()
    rep = json.loads(out.read_text(encoding="utf-8"))
    expect(rep["clean"] is True and rep["n_z_nolive"] == 0,
           "an unobstructed i=3 sweep still certifies clean")
    expect(rep["certificate"] is not None, "i=3 certificate still emitted")


def test_rcache_matches_uncached() -> None:
    """The r(p) cache is an optimisation: rows must be byte-identical.

    r(p)=C(N,K) mod p does not depend on k, so caching it per prime is
    sound -- but only if nothing else in the walk depends on k through
    it. This pins that: same range, cache on and off, same rows.
    """
    mem = fib_member(8, compute_m=False)
    N, K = mem.n, mem.k
    ks = [k for k in range(2000, 2120) if k not in (K, K + 1)]
    plain = [nps.walk_until_kill(N, K, k, nprimes=20, max_gap=0) for k in ks]
    rc = nps.RCache(N, K)
    cached = [
        nps.walk_until_kill(N, K, k, nprimes=20, max_gap=0, rcache=rc)
        for k in ks
    ]
    expect(plain == cached, "cached and uncached walks produce identical rows")
    expect(rc.lookups > rc.misses,
           f"the cache actually hits (reuse={rc.stats()['reuse']}x over {len(ks)} columns)")
    expect(all(rc._r[q] == binom_mod_lucas(N, K, q) for q in rc._r),
           "every cached r(p) equals a fresh Lucas evaluation")


def test_rcache_matches_the_recorded_runs() -> None:
    """Anchor the cached walk to ground truth, not to the uncached path.

    test_rcache_matches_uncached compares current code against current code:
    a regression touching both legs passes it. These two results files are
    from the original Stage-2 runs and are independent of anything the cache
    could break.

    Two annotation fields have to be excluded, and they are worth naming
    because they look like failures: `stop` was added to walk_until_kill
    after both files were written, and `needed_more_than_20` was written by
    an older vintage into the nprimes=20 run only. Every substantive field
    -- k, q, gap, r, survived_all -- must match exactly.
    """
    mem = fib_member(8, compute_m=False)
    N, K = mem.n, mem.k
    ANNOTATIONS = ("stop", "needed_more_than_20")
    checked = 0
    for name, lo, hi in (
        ("nextprime_i8_k401-10000.json", 401, 1200),
        ("nextprime_i8_k10001-100000.json", 10001, 10800),
    ):
        path = ROOT / "results" / name
        if not path.exists():
            ok.append(f"{name} absent; skipped ground-truth check")
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        rows = {r["k"]: r for r in payload["rows"]}
        nprimes = payload["nprimes"]
        rc = nps.RCache(N, K)
        bad = []
        for k in range(lo, hi):
            if k not in rows:
                continue
            got = nps.walk_until_kill(N, K, k, nprimes=nprimes, max_gap=0, rcache=rc)
            want = dict(rows[k])
            for f in ANNOTATIONS:
                got.pop(f, None)
                want.pop(f, None)
            checked += 1
            if got != want:
                bad.append((k, got, want))
        expect(
            not bad,
            f"cached walk reproduces {name} exactly "
            f"(nprimes={nprimes}, reuse={rc.stats()['reuse']}x)",
        )
        if bad:
            errors.append(f"  first: k={bad[0][0]} got {bad[0][1]} want {bad[0][2]}")
    expect(checked > 1000, f"ground-truth check covered {checked} recorded columns")


def test_rcache_zero_is_cached() -> None:
    """r(p)=0 is a value, not a miss.

    Dead primes are the cheap case, so a `if not val` test here would be
    invisible in a timing but would recompute them forever.
    """
    mem = fib_member(8, compute_m=False)
    N, K = mem.n, mem.k
    zero_p = next(
        q for q in nps.next_primes_after(2000, 400) if binom_mod_lucas(N, K, q) == 0
    )
    rc = nps.RCache(N, K)
    expect(rc.r(zero_p) == 0, f"r({zero_p})=0 for i=8")
    rc.r(zero_p)
    rc.r(zero_p)
    expect(rc.misses == 1,
           f"a zero r(p) is evaluated once, not per lookup (misses={rc.misses})")


def test_rcache_prune_keeps_answers() -> None:
    """Pruning drops only primes no later column can ask about."""
    mem = fib_member(8, compute_m=False)
    N, K = mem.n, mem.k
    rc = nps.RCache(N, K)
    ps = nps.next_primes_after(5000, 12)
    for q in ps:
        rc.r(q)
    cut = ps[5]
    rc.prune(cut)
    expect(all(q > cut for q in rc._r),
           "prune drops every entry with p <= k")
    expect(set(rc._r) == {q for q in ps if q > cut},
           "prune keeps every entry a later column can still reach")
    before = rc.misses
    expect(rc.r(ps[-1]) == binom_mod_lucas(N, K, ps[-1]) and rc.misses == before,
           "a surviving entry is still served from the cache after a prune")


def test_checkpoint_schema_guard() -> None:
    """Resume must refuse a checkpoint from different code or parameters.

    Resume merges old records with new ones, so a silent schema or parameter
    change yields a certificate over columns that were never all tested the
    same way. Refusing is cheaper than trusting it.
    """
    import tempfile

    fam = bk.make_fam(3)
    ident = dict(i=3, N=fam.N, K=fam.K, k_max=342, k_lo_z=201,
                 cap_bii=14, cap_z=12)
    with tempfile.TemporaryDirectory() as td:
        good = Path(td) / "good.jsonl"
        bk.append_jsonl(good, bk.checkpoint_identity(**ident))
        try:
            bk.check_checkpoint(good, **ident)
            expect(True, "matching schema header is accepted")
        except SystemExit as exc:
            expect(False, f"matching header rejected: {exc}")

        bad = Path(td) / "bad.jsonl"
        bk.append_jsonl(bad, bk.checkpoint_identity(**{**ident, "k_max": 999}))
        try:
            bk.check_checkpoint(bad, **ident)
            expect(False, "mismatched k_max was accepted")
        except SystemExit:
            expect(True, "a changed run parameter is refused")

        legacy = Path(td) / "legacy.jsonl"
        bk.append_jsonl(legacy, {"tag": "bii1", "p": 359, "survivors": []})
        try:
            bk.check_checkpoint(legacy, **ident)
            expect(False, "headerless legacy checkpoint was accepted")
        except SystemExit:
            expect(True, "a checkpoint with no schema header is refused")

        expect(bk.check_checkpoint(Path(td) / "absent.jsonl", **ident) is None,
               "a missing checkpoint is fine (fresh run)")


def main() -> int:
    before = _digest_guarded()
    with sandboxed_results():
        test_rcache_matches_uncached()
        test_rcache_matches_the_recorded_runs()
        test_rcache_zero_is_cached()
        test_rcache_prune_keeps_answers()
        test_checkpoint_schema_guard()
        test_phase_count()
        test_resume_of_unfinished_phase_is_not_clean()
        test_nolive_is_not_clean()
        test_clean_run_still_clean()
    # Guard the sandbox itself. If the redirect is ever removed, this fails
    # loudly instead of quietly regenerating tracked artifacts again.
    after = _digest_guarded()
    changed = [n for n in GUARDED if before[n] != after[n]]
    expect(not changed,
           "the suite left every tracked results artifact byte-identical"
           + (f" -- MUTATED: {changed}" if changed else ""))
    print("\n=== SWEEP TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
