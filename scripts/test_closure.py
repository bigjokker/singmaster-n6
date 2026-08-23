#!/usr/bin/env python3
"""End-to-end regression for the path that shipped a false certificate.

The pieces are tested elsewhere -- the builder in test_witness, the scan in
test_kernel, the sweep in test_sweeps. What no test drove is the SEQUENCE that
produced i=8's k=1021 -> p=3517:

    sweep ends clean=False with columns still alive
      -> build must leave them UNRESOLVED, not invent a killer
      -> fill must close them with the engine
      -> ledger must read COMPLETE
      -> and every filled witness must actually verify

Past CAP_Z the sweep defers big-k survivors: it stops testing them, so they
stop appearing in the records. The builder used to read that absence as
"killed by this round's prime" and fabricate a witness. Nothing downstream
disagreed -- the coverage ledger checks that every column HAS a witness, not
that the witness is VALID -- so it reported "missing 0, ok=True" over a
certificate check_witness rejects.

i=8 takes ~3 minutes, too slow for a suite, so the same shape is forced on a
small family by capping the Z rounds at 1 and setting SMALL_K a quarter into
the Z range, so most round-1 survivors defer while the rest keep running and
generate the later rounds. That split matters: defer EVERYTHING and the round
loop simply ends, leaving no later round to mis-credit, and the test passes on
the broken builder. It deliberately produces MORE than 100 survivors,
because the sweep json only lists them when there are <= 100. Above that it
writes [], `declared_alive` comes back empty, and coverage is the only guard
left that can fire. That is the blind spot this test exists to cover.

WHAT THIS DOES AND DOES NOT PIN. Against the pre-fix builder (8e64945) this
exits 1, so it is a real regression test -- but it fails by RAISING from the
old `p in seen` guard, not by minting a false witness. At i=5 the deferred
column's rebuilt prime is absent from the round's prime set, so the old guard
happens to catch it; at i=8 that prime WAS present, the guard passed, and the
certificate went out silently. Forcing that coincidence is data-dependent, so
the SILENT fabrication is pinned synthetically in test_witness.py
(test_deferred_column_is_not_credited_a_kill), where the records are built so
the prime is in `seen`. This file pins the other half: that the whole closure
path completes, and completes CORRECTLY, when columns are deferred.

Run: python scripts/test_closure.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import coverage_ledger as CL  # noqa: E402
import family_sweep as fsw  # noqa: E402
import test_sweeps as T  # noqa: E402
import witness as W  # noqa: E402

I_TEST = 5
ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


def test_deferred_columns_close_through_the_engine() -> None:
    from bandii_kernel import make_fam

    # SMALL_K must SPLIT the survivors, not defer all of them. If everything
    # defers, `current` empties, the round loop breaks, and there is no later
    # round for the old builder to credit -- the bug's precondition is gone and
    # the test passes on broken code. (It did, on the first version of this
    # file.) A quarter into the Z range defers most of them while leaving
    # enough columns running to generate rounds 2+, which is the i=8 shape:
    # k=1021 deferred while small-k columns carried on into rounds 13-14.
    fam = make_fam(I_TEST)
    k_lo_z = fsw.K_EXACT.get(I_TEST, 2) + 1
    split = k_lo_z + (fam.K - k_lo_z) // 4

    caps = (fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K)
    with T.sandboxed_results() as tmp:
        try:
            # The loop is bounded by cap_z but the deferral tests `rnd >
            # CAP_Z`, so deferral only ever happens when CAP_Z_SMALL_K > CAP_Z
            # -- the default 15 > 12 is why i=8 deferred at rounds 13-14. Set
            # them equal and deferral never fires at all.
            fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K = 1, 3, split
            out, chk = fsw.paths(I_TEST)
            sys.argv = ["family_sweep.py", "--i", str(I_TEST)]
            t = time.perf_counter()
            fsw.main()
            sweep_s = time.perf_counter() - t
        finally:
            fsw.CAP_Z, fsw.CAP_Z_SMALL_K, fsw.SMALL_K = caps

        rep = json.loads(out.read_text(encoding="utf-8"))
        n_alive = int(rep.get("n_z_alive") or 0)
        expect(rep.get("clean") is False and n_alive > 0,
               f"forced deferral leaves the sweep unclean "
               f"({n_alive} columns still alive, {sweep_s:.1f}s)")
        if not n_alive:
            return

        # Above 100 the json writes [] -- so declared_alive is empty and the
        # builder cannot lean on the run's own survivor list.
        listed = rep.get("z_survivors") or []
        expect(not listed and n_alive > 100,
               f"the survivor list is truncated ({n_alive} alive, "
               f"{len(listed)} listed), so coverage is the only guard left")

        wpath = tmp / "results" / f"i{I_TEST}_witness.npz"
        meta = W._build_family(I_TEST, chk, wpath)
        expect(meta["n_unresolved"] == n_alive,
               f"builder leaves every deferred column unresolved "
               f"({meta['n_unresolved']} of {n_alive}) rather than "
               f"crediting a prime that never tested it")

        # and none of them slipped into the table with a fabricated witness
        ks, _ps, _ = W.load(wpath)
        got = set(int(x) for x in ks)
        unresolved = set(meta.get("unresolved") or [])
        expect(unresolved and not (unresolved & got),
               f"no deferred column carries a witness in the built table "
               f"({len(unresolved)} unresolved, 0 present)")

        fl = W.fill_small_gaps(wpath, I_TEST)
        expect(fl["added"] == n_alive and not fl["unresolved"],
               f"the engine closes every deferred column "
               f"(added {fl['added']}, {len(fl['unresolved'])} left)")

        aud = CL.audit_member(I_TEST, wpath)
        expect(aud["ok"] and aud["n_missing"] == 0 and aud["n_extra"] == 0,
               f"ledger COMPLETE after fill ({aud['n_witnessed']:,} witnessed, "
               f"missing {aud['n_missing']}, extra {aud['n_extra']})")

        # Coverage complete is not validity. Check the filled columns really do
        # kill -- this is the assertion the original defect would have failed.
        ks2, ps2, _ = W.load(wpath)
        pos = {int(k): int(p) for k, p in zip(ks2, ps2)}
        sample = sorted(unresolved)[:40]
        bad = [(k, pos.get(k)) for k in sample
               if k not in pos
               or not W.check_witness(rep["N"], rep["K"], k, pos[k])["ok"]]
        expect(not bad and sample,
               f"every sampled engine-filled witness verifies "
               f"({len(sample)} checked, {len(bad)} invalid {bad[:2]})")


# Tracked artifacts this suite must leave untouched. fill runs on the sandbox
# table only; if a table change ever reaches results/ again, fail loudly.
GUARDED = (f"i{I_TEST}_sweep.json", f"i{I_TEST}_sweep.jsonl", f"i{I_TEST}_witness.npz")


def _digest_guarded() -> dict:
    import hashlib

    out = {}
    for name in GUARDED:
        p = ROOT / "results" / name
        out[name] = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
    return out


def main() -> int:
    t = time.perf_counter()
    before = _digest_guarded()
    test_deferred_columns_close_through_the_engine()
    after = _digest_guarded()
    changed = [n for n in GUARDED if before[n] != after[n]]
    expect(not changed,
           "the suite left every tracked results artifact byte-identical"
           + (f" -- MUTATED: {changed}" if changed else ""))
    print(f"\n=== CLOSURE TESTS (i={I_TEST}, {time.perf_counter()-t:.1f}s) ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
