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
    import witness as _w

    tmp = Path(tempfile.mkdtemp(prefix="sweep_test_"))
    (tmp / "results").mkdir()
    real = fsw.ROOT
    real_w = _w.ROOT
    fsw.ROOT = tmp
    # witness resolves its own ROOT for CLI defaults; redirect it too, so a
    # fill/repair/retarget reached through the sweep can only ever see the
    # sandbox. (The binding itself is by table directory -- see
    # witness.retarget_certificate -- this is the belt to that brace.)
    _w.ROOT = tmp
    try:
        yield tmp
    finally:
        fsw.ROOT = real
        _w.ROOT = real_w
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


def test_certificate_basis_follows_k_exact() -> None:
    """The certificate sentence must state what actually closed the columns
    below the Z-jump's start, per member -- not a hardcoded 'exact k<=200'.

    K_EXACT is {2..7: 200 (exact intersect), 9: 80 (modular)} and i=8 has no
    exact band at all (k_lo_z=3, k=2 from the engine). The emitted sentence
    said 'Together with exact k<=200' for every member, which is false for
    i=8 and i=9. Latent only because both of their shipped records carry
    certificate=null -- but a future clean run would have shipped it.
    """
    basis = getattr(fsw, "certificate_basis", None)
    expect(callable(basis), "family_sweep exposes certificate_basis(i, k_lo_z)")
    if not callable(basis):
        return
    b3 = basis(3, fsw.K_EXACT[3] + 1)
    b9 = basis(9, fsw.K_EXACT[9] + 1)
    b8 = basis(8, fsw.K_EXACT.get(8, 2) + 1)
    expect("exact" in b3 and "200" in b3,
           f"i=3 basis names the exact k<=200 run ({b3!r})")
    expect("modular" in b9 and "80" in b9 and "200" not in b9,
           f"i=9 basis names the modular k<=80 run, not exact k<=200 ({b9!r})")
    expect("200" not in b8 and "exact" not in b8,
           f"i=8 basis claims no exact band ({b8!r})")
    out, _chk = fsw.paths(3)
    if out.exists():
        rep = json.loads(out.read_text(encoding="utf-8"))
        cert = rep.get("certificate") or ""
        expect(b3 in cert,
               "the i=3 sweep record's certificate carries the per-member basis clause")


def test_sparse_records_carry_their_column_list() -> None:
    """Q30. A chunk record whose columns are not contiguous must say WHICH
    columns it scanned; a contiguous one must not (99.8% of i=8's records
    are contiguous and their [k_lo, k_hi] is already exact).

    Without the list, a sparse record's span over-covers, and a column the
    sweep deferred that falls inside it can only be told from a killed one
    by the count identity -- which stops the build rather than repairing it.
    """
    fam = bk.make_fam(5)
    ivs = bk.live_intervals(bk.cells(fam), fam)
    a = 300
    p = bk.first_live_after(a, ivs, fam.D)
    r = int(binom_mod_lucas(fam.N, fam.K, p))
    expect(r != 0, f"fixture: p={p} is live for i=5")
    # same bucket, a gap inside: [a, a+1, a+3] spans 4 columns but scans 3
    rec = fsw._job(("z", p, [a, a + 1, a + 3], fam.N, fam.K, r))
    expect(rec["n_cols"] == 3 and rec["k_lo"] == a and rec["k_hi"] == a + 3,
           "fixture: the chunk is sparse (3 columns over a span of 4)")
    expect(rec.get("ks") == [a, a + 1, a + 3],
           f"a sparse chunk record carries its exact column list (got {rec.get('ks')})")
    dense = fsw._job(("z", p, [a, a + 1, a + 2], fam.N, fam.K, r))
    expect("ks" not in dense,
           "a contiguous chunk record carries no column list (its span is exact)")
    bii = fsw._job(("bii", bk.first_primes_above(fam.N2, fam.D, bk.kmax_of(fam)[0])[0],
                    [fam.K + 2, fam.K + 5], fam.N, fam.K, None))
    expect(bii.get("ks") == [fam.K + 2, fam.K + 5],
           "Band II sparse chunks carry the list as well")


def test_resume_refuses_a_changed_chunk_partition() -> None:
    """D3. The chunk partition is part of a run's identity; resuming under a
    different one must be REFUSED, not merged.

    done_keys are (tag, p, k_lo, k_hi) -- chunk boundaries. With N_CHUNKS
    absent from the checkpoint header, a resume under N_CHUNKS=16 of a pass
    written at 32 passed check_checkpoint, every new chunk looked pending and
    was re-scanned, AND run_jobs recovered the old records' survivors of the
    same tag: every pass-1 survivor twice (327 -> 654 at i=5), a false
    ESCALATE, doubled n_alive in the record, clean=True. The only thing that
    stopped a certificate was the witness builder's count identity -- a guard
    built for a different bug. check_checkpoint's docstring already promised
    to refuse exactly this merge.

    Three things pinned: a changed partition against a header that carries
    the key is refused; the same partition resumes and reproduces the
    uninterrupted run; a legacy header without the key still resumes at the
    default partition (the 90f105c compatibility rule -- no schema bump).
    """
    i = 5
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    sys.argv = ["family_sweep.py", "--i", str(i)]
    fsw.main()
    base = json.loads(out.read_text(encoding="utf-8"))
    base_alive1 = base["phases"]["bandii"][0]["n"]
    expect(base["clean"] is True and base["certificate"], "fixture: baseline i=5 run is clean")

    lines = chk.read_text(encoding="utf-8").splitlines()
    header = lines[0]
    bii1 = [ln for ln in lines[1:] if '"tag": "bii1"' in ln]
    expect(json.loads(header).get("event") == "schema" and len(bii1) == fsw.N_CHUNKS,
           f"fixture: pass 1 was written as {len(bii1)} chunks at N_CHUNKS={fsw.N_CHUNKS}")

    def cut(hdr: str) -> None:
        out.unlink(missing_ok=True)
        chk.write_text("\n".join([hdr] + bii1) + "\n", encoding="utf-8")

    saved = fsw.N_CHUNKS
    try:
        # 1. changed partition against a header that records it: refuse
        cut(header)
        fsw.N_CHUNKS = 16
        refused = False
        try:
            fsw.main()
        except SystemExit as exc:
            refused = "n_chunks" in str(exc)
        if refused:
            expect(True, "a resume under a different N_CHUNKS is refused by the header")
        else:
            rep = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {}
            got = rep.get("phases", {}).get("bandii", [{}])[0].get("n")
            expect(False, f"a resume under a different N_CHUNKS was ACCEPTED "
                          f"(pass-1 n_alive {got} vs {base_alive1}, clean={rep.get('clean')})")

        # 2. same partition: resumes and reproduces the uninterrupted run
        fsw.N_CHUNKS = saved
        cut(header)
        fsw.main()
        rep = json.loads(out.read_text(encoding="utf-8"))
        expect(rep["phases"]["bandii"][0]["n"] == base_alive1
               and rep["clean"] is True and rep["certificate"] == base["certificate"],
               f"the same partition resumes and reproduces the run "
               f"(pass-1 n_alive {rep['phases']['bandii'][0]['n']} vs {base_alive1})")

        # 3. legacy header without the key: still resumes at the default
        legacy = json.loads(header)
        legacy.pop("n_chunks", None)
        cut(json.dumps(legacy))
        fsw.main()
        rep = json.loads(out.read_text(encoding="utf-8"))
        expect(rep["phases"]["bandii"][0]["n"] == base_alive1 and rep["clean"] is True,
               "a legacy header without n_chunks still resumes at the default partition")
    finally:
        fsw.N_CHUNKS = saved
        for f in (out, chk):
            f.unlink(missing_ok=True)


def test_cap_survivors_are_judged_at_phase_end() -> None:
    """D4 wiring. When a phase ends with columns alive at its cap, the sweep
    record's escalation block must carry one judged row per survivor --
    Lambda over the primes it actually faced, peers x Lambda against the
    same threshold -- for BOTH phases. Caps of 1 at i=4 leave both phases
    with survivors; every one faced exactly one prime, so all are ordinary,
    and the point here is the wiring, not the verdict (test_sizelaw pins
    the verdicts).
    """
    i = 4
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    saved = (fsw.CAP_BII, fsw.CAP_Z, fsw.CAP_Z_SMALL_K)
    try:
        fsw.CAP_BII, fsw.CAP_Z, fsw.CAP_Z_SMALL_K = 1, 1, 1
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
    finally:
        fsw.CAP_BII, fsw.CAP_Z, fsw.CAP_Z_SMALL_K = saved
    rep = json.loads(out.read_text(encoding="utf-8"))
    expect(rep["clean"] is False and rep["n_bii_alive"] > 0 and rep["n_z_alive"] > 0,
           f"fixture: caps of 1 leave survivors in both phases "
           f"(bii {rep['n_bii_alive']}, z {rep['n_z_alive']})")
    for phase, n_alive, peers in (("bandii", rep["n_bii_alive"], rep["n_bii"]),
                                  ("zjump", rep["n_z_alive"], rep["n_z"])):
        v = rep["escalation"][phase]
        rows = v.get("cap_survivors") or []
        expect(v.get("n_cap_survivors") == n_alive and v.get("cap_peers") == peers,
               f"{phase}: every cap survivor is judged ({v.get('n_cap_survivors')} of "
               f"{n_alive}) against the {peers} columns that entered the phase")
        expect(rows and all(r["run"] == 1 and "lambda" in r and "expected" in r
                            and r["expected"] == peers * r["lambda"] for r in rows),
               f"{phase}: each row carries run, Lambda and expected = peers x Lambda")
        expect(all(r["escalate"] is False for r in rows) and not v.get("escalating_survivors")
               and v["escalate"] is False,
               f"{phase}: one-prime survivors are ordinary, phase verdict does not fire")
    for f in (out, chk):
        f.unlink(missing_ok=True)


def test_import_failure_does_not_lose_the_sweep_json() -> None:
    """A finished sweep must write its record even if `import witness` fails.

    The except handler around the witness build referenced `wpath`, which was
    assigned AFTER the import inside the same try -- so an ImportError turned
    into UnboundLocalError in the handler, main() crashed, and a finished
    multi-hour sweep had no json at all, despite the comment 'never lose a
    finished sweep to this'.
    """
    i = 3
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    saved = sys.modules.get("witness")
    sys.modules["witness"] = None            # makes `import witness` raise
    try:
        sys.argv = ["family_sweep.py", "--i", str(i)]
        try:
            rc = fsw.main()
        except Exception as exc:            # noqa: BLE001
            rc = f"raised {type(exc).__name__}: {exc}"
    finally:
        if saved is None:
            sys.modules.pop("witness", None)
        else:
            sys.modules["witness"] = saved
    expect(rc == 0 and out.exists(),
           f"a failed `import witness` still writes the sweep json (rc={rc!r}, "
           f"json exists={out.exists()})")
    if out.exists():
        rep = json.loads(out.read_text(encoding="utf-8"))
        expect(rep.get("witness_ok") is False and rep.get("certificate") is None
               and rep.get("clean") is True,
               "and the record says witness_ok=False, certificate withheld, "
               "sweep verdict preserved")


def test_m_route_matches_lucas_route_exactly() -> None:
    """r(p) from m must produce the SAME sweep as r(p) from Lucas.

    The Z-jump's r(p) now comes from reducing m = C(N,K) against the whole
    round's prime set with a product/remainder tree in the parent, instead of
    two-digit Lucas per job. That is an optimisation, so the only acceptable
    evidence is identical output -- not a plausible argument, and not a faster
    number.

    Measured at i=9 on the real 990,683-prime ladder: 2.07 us/prime against
    3.14 ms for direct m %% p and ~86 ms for Lucas, i.e. 23.82 core-h of r(p)
    becomes ~2.9 s. None of that matters if a single survivor moves.
    """
    i = 5
    out, chk = fsw.paths(i)
    saved = fsw.USE_M_FOR_RP
    reports = {}
    try:
        for flag in (False, True):
            fsw.USE_M_FOR_RP = flag
            fsw._M_KEY = None                   # force a rebuild, do not reuse
            for f in (out, chk):
                f.unlink(missing_ok=True)
            sys.argv = ["family_sweep.py", "--i", str(i)]
            fsw.main()
            reports[flag] = json.loads(out.read_text(encoding="utf-8"))
    finally:
        fsw.USE_M_FOR_RP = saved
        fsw._M_KEY = None
        for f in (out, chk):
            f.unlink(missing_ok=True)

    lucas, viam = reports[False], reports[True]
    diff = _diff_ignoring_timing(lucas, viam)
    expect(not diff,
           "the m route reproduces the Lucas route exactly"
           + (f" -- {len(diff)} differences, first: {diff[0]}" if diff else ""))
    expect(lucas.get("certificate") == viam.get("certificate")
           and lucas.get("certificate") is not None,
           "and emits a byte-identical certificate")
    expect(lucas.get("witness", {}).get("sha256")
           == viam.get("witness", {}).get("sha256"),
           "and the witness table it produces is byte-identical "
           f"({str(lucas.get('witness', {}).get('sha256'))[:16]} vs "
           f"{str(viam.get('witness', {}).get('sha256'))[:16]})")
    # m must not leak into any artifact
    blob = json.dumps(viam)
    expect(len(blob) < 200_000 and "bincoef" not in blob,
           "m does not leak into the sweep record")


def test_incremental_done_keys_equals_full_rederivation() -> None:
    """The set maintained in place must equal the one read off the whole file.

    done_keys used to be re-derived from the entire checkpoint after every
    round, and run_jobs re-parsed it again on entry -- about 42 full passes at
    i=9 over a file reaching 219 MB, at ~5.7x file size in RAM. That is roughly
    1.25 GB of Python objects in the parent, 24 times, on a machine already
    crashed once by memory pressure.

    Maintaining it incrementally is only safe if it is EQUAL, so VERIFY_DONE_KEYS
    re-derives and asserts equality every round. That is the expensive path by
    design: it is switched on here and off in production.

    Exercised on a fresh run AND a resume, because the resume path is where the
    two false-clean certificates lived and where a missing key would silently
    re-run or silently skip a chunk.
    """
    i = 5
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)

    saved = fsw.VERIFY_DONE_KEYS
    fsw.VERIFY_DONE_KEYS = True
    try:
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
        full = json.loads(out.read_text(encoding="utf-8"))
        expect(full.get("clean") is True,
               "fresh run with VERIFY_DONE_KEYS on is still clean "
               "(the equality assertion did not fire)")

        lines = chk.read_text(encoding="utf-8").splitlines()
        chk.write_text("\n".join(lines[: max(2, int(len(lines) * 0.6))]) + "\n",
                       encoding="utf-8")
        out.unlink(missing_ok=True)
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
        got = json.loads(out.read_text(encoding="utf-8"))
        expect(got.get("clean") is True,
               "resumed run with VERIFY_DONE_KEYS on is still clean")
        expect(got.get("certificate") == full.get("certificate"),
               "and the resumed certificate is byte-identical to the "
               "uninterrupted one")
    except RuntimeError as exc:
        expect(False, f"incremental done_keys diverged: {str(exc)[:120]}")
    finally:
        fsw.VERIFY_DONE_KEYS = saved
        for f in (out, chk):
            f.unlink(missing_ok=True)


def test_extended_rounds_defer_rather_than_drop() -> None:
    """A survivor the extended rounds do not serve must still be REPORTED.

    Past CAP_Z the Z-jump serves only the small-k tail. It used to do that with
        current = [c for c in current if int(c["k"]) < SMALL_K]
    which does not merely stop testing the large-k survivors -- it deletes them
    from the accounting. They reach neither z_none nor z_left, so n_z_alive can
    fall to 0 and the run certifies over columns it never killed.

    That is not hypothetical. i=8's round 12 held k=145 and k=1021; round 13
    recorded n=1, mean_k=145; k=1021 was dropped; clean went True; and the
    witness builder recorded p=3517 -- a prime k=1021 had merely SURVIVED -- as
    its killer. Sampled verification (5,000 of 5,182,634) never looked.

    Forced here by shrinking CAP_Z so ordinary survivors fall into the extended
    rounds, which is the only way to reach this branch on a small member.
    """
    i = 5
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)

    saved = (fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K)
    fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K = 1, 3, 1000
    try:
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
    finally:
        fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K = saved

    rep = json.loads(out.read_text(encoding="utf-8"))
    n_alive = rep.get("n_z_alive", 0)
    big = [int(c["k"]) for c in (rep.get("z_survivors") or [])
           if int(c["k"]) >= 1000]
    expect(n_alive > 0,
           f"deferred large-k survivors are counted, not dropped "
           f"(n_z_alive={n_alive})")
    expect(rep.get("clean") is False,
           f"a deferred survivor blocks clean exactly as an ordinary one does "
           f"(clean={rep.get('clean')!r})")
    expect(rep.get("certificate") is None,
           "and no certificate is emitted over a column that was never killed")
    # it must be SURVIVORS blocking clean, not untestable columns -- the
    # nolive path is a different mechanism with a different meaning
    expect(rep.get("n_z_nolive", 0) == 0,
           f"clean is blocked by deferred SURVIVORS, not by nolive columns "
           f"(n_z_nolive={rep.get('n_z_nolive')})")

    for f in (out, chk):
        f.unlink(missing_ok=True)


def test_failed_witness_build_withholds_the_certificate() -> None:
    """A certificate names a witness table; if the build failed, do not emit it.

    Before the fix the `except` around the witness build printed the error,
    left clean=True, and still emitted a certificate reading
    "results/i{i}_witness.npz (sha256 n/a)". The dangerous case is a STALE npz
    from an earlier run: a referee who opens the named path verifies a table
    this sweep never wrote, and it passes.

    `clean` is deliberately NOT flipped. It is the SWEEP's verdict -- every
    column testable, nothing surviving -- and a write failure does not make a
    clean sweep dirty. clean=True with certificate=None is the honest state.
    """
    import witness as _w

    i = 3
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)

    # a stale table from a notional earlier run, at the exact path the
    # certificate would name
    wpath = fsw.ROOT / "results" / f"i{i}_witness.npz"
    wpath.parent.mkdir(parents=True, exist_ok=True)
    wpath.write_bytes(b"stale-not-this-runs-output")
    stale_before = wpath.read_bytes()

    real_build = _w._build_family

    def boom(*a, **kw):
        raise RuntimeError("simulated witness build failure")

    _w._build_family = boom
    try:
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
    finally:
        _w._build_family = real_build

    rep = json.loads(out.read_text(encoding="utf-8"))
    expect(rep.get("certificate") is None,
           f"a failed witness build emits NO certificate "
           f"(got {str(rep.get('certificate'))[:40]!r})")
    expect(rep.get("witness_ok") is False,
           f"witness_ok records the failure (got {rep.get('witness_ok')!r})")
    expect(rep.get("clean") is True,
           "the SWEEP's own verdict is preserved -- a write failure does not "
           f"make a clean sweep dirty (clean={rep.get('clean')!r})")
    expect("stale" in (rep.get("certificate_withheld") or "").lower()
           or "NOT written" in (rep.get("certificate_withheld") or ""),
           f"the withholding reason names the stale-file hazard "
           f"({str(rep.get('certificate_withheld'))[:90]!r})")
    expect(wpath.read_bytes() == stale_before,
           "the stale table is left untouched, not silently overwritten")

    for f in (out, chk):
        f.unlink(missing_ok=True)
    wpath.unlink(missing_ok=True)


def _diff_ignoring_timing(a, b, path: str = "") -> list[str]:
    """Recursive compare, ignoring wall-clock fields only.

    `seconds` legitimately differs between an uninterrupted run and a resumed
    one. Everything else -- every survivor, every round count, the certificate
    string itself -- must be identical, or the resume changed the claim.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k == "seconds":
                continue
            if k not in a or k not in b:
                out.append(f"{path}.{k}: present in only one")
                continue
            out += _diff_ignoring_timing(a[k], b[k], f"{path}.{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        out = []
        for n, (x, y) in enumerate(zip(a, b)):
            out += _diff_ignoring_timing(x, y, f"{path}[{n}]")
        return out
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def test_resume_reproduces_an_uninterrupted_run() -> None:
    """A run killed mid-round and resumed must emit the SAME certificate.

    The resume path is where both "false clean certificate" bugs lived. The
    two tests below cover a resumed unfinished phase and an untestable column,
    but neither checks the property that matters most: that resuming AGREES
    WITH not having been interrupted.

    Truncating the checkpoint is exactly what a crash leaves on disk -- some
    job records written, no output json -- so it exercises the real resume
    branch rather than a mock of it. Two cut points, one mid-phase and one
    later, because the Band II and Z-jump resume paths are separate code.
    """
    i = 5
    out, chk = fsw.paths(i)
    for f in (out, chk):
        f.unlink(missing_ok=True)
    sys.argv = ["family_sweep.py", "--i", str(i)]
    fsw.main()
    full = json.loads(out.read_text(encoding="utf-8"))
    lines = chk.read_text(encoding="utf-8").splitlines()
    expect(full.get("clean") is True and full.get("certificate"),
           f"baseline i={i} run is clean and certified (for the resume comparison)")

    # `escalation` and `phases` are rebuilt only for phases this run actually
    # executed. A resume that skips an already-complete phase cannot
    # reconstruct them from the checkpoint, so it emits evaluated=False with a
    # note saying why. That is DECLARED, not silent -- so assert the
    # declaration rather than equality, and pin everything else exactly.
    #
    # This is not cosmetic: i=9 resumed, and its recorded json carries NO
    # escalation block for either phase. The certificates are unaffected (each
    # kill has an independently checkable witness), but the size-law anomaly
    # trigger did not run on that member and a human checked the Band II curve
    # by hand instead.
    LOST_ON_RESUME = {"escalation", "phases"}

    for frac in (0.35, 0.75):
        keep = max(2, int(len(lines) * frac))
        chk.write_text("\n".join(lines[:keep]) + "\n", encoding="utf-8")
        out.unlink(missing_ok=True)
        sys.argv = ["family_sweep.py", "--i", str(i)]
        fsw.main()
        got = json.loads(out.read_text(encoding="utf-8"))
        pct = int(frac * 100)

        claim_full = {k: v for k, v in full.items() if k not in LOST_ON_RESUME}
        claim_got = {k: v for k, v in got.items() if k not in LOST_ON_RESUME}
        diff = _diff_ignoring_timing(claim_full, claim_got)
        expect(not diff,
               f"resume from a {pct}% checkpoint reproduces the uninterrupted CLAIM"
               + (f" -- {len(diff)} differences, first: {diff[0]}" if diff else ""))
        expect(got.get("certificate") == full.get("certificate"),
               f"resume from {pct}% emits a byte-identical certificate")

        # whatever escalation it does NOT evaluate, it must say so
        undeclared = [
            ph for ph, blk in (got.get("escalation") or {}).items()
            if isinstance(blk, dict) and not blk.get("evaluated")
            and "resumed" not in (blk.get("note") or "")
        ]
        expect(not undeclared,
               f"resume from {pct} declares every escalation it could not evaluate"
               + (f" -- silent for {undeclared}" if undeclared else ""))


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
        test_certificate_basis_follows_k_exact()
        test_sparse_records_carry_their_column_list()
        test_cap_survivors_are_judged_at_phase_end()
        test_import_failure_does_not_lose_the_sweep_json()
        test_resume_reproduces_an_uninterrupted_run()
        test_resume_refuses_a_changed_chunk_partition()
        test_m_route_matches_lucas_route_exactly()
        test_incremental_done_keys_equals_full_rederivation()
        test_extended_rounds_defer_rather_than_drop()
        test_failed_witness_build_withholds_the_certificate()
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
