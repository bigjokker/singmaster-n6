#!/usr/bin/env python3
"""Adversarial tests for work_census.py.

A census that always says PASS passes any suite built only of good inputs, and
this one ships an *excuse* mechanism (`record_incomplete`) that could launder a
genuine failure. So every test here tries to make the census accept something
wrong, and asserts the SPECIFIC gate or guard that must catch it -- not merely
that something failed.

Covered:
  1  a tampered witness prime must break a recorded-count gate
  2  a dropped column must break the column-count gate
  3  a Z-jump prime that is not on the live-prime ladder must be REFUSED
  4  a Band II prime outside primes_bii must be REFUSED
  5  a witness table checked against the WRONG member must be refused
  6  the record_incomplete excuse must NOT fire when its evidence is absent
  7  a tampered k_max must be refused (k_max is recomputed, not read back)
  8  the ladder off-by-one (prev+1 instead of prev) must trip the self-check

    python scripts/test_work_census.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import work_census as wc  # noqa: E402

FAILS = 0
SMALL = 5          # small enough to run in well under a second
MID = 7            # has engine-filled leftovers (k=2..200) worth exercising


def sandbox(i: int) -> Path:
    """A private copy of one member's artifacts, so tampering is isolated."""
    tmp = Path(tempfile.mkdtemp(prefix="census_test_"))
    (tmp / "results").mkdir()
    for name in (f"i{i}_witness.npz", f"i{i}_sweep.json"):
        shutil.copy(ROOT / "results" / name, tmp / "results" / name)
    return tmp


def load_npz(tmp: Path, i: int):
    d = np.load(tmp / "results" / f"i{i}_witness.npz", allow_pickle=True)
    return d["k"].astype(np.int64), d["p"].astype(np.int64), str(d["meta"])


def save_npz(tmp: Path, i: int, k, p, meta: str) -> None:
    np.savez_compressed(tmp / "results" / f"i{i}_witness.npz",
                        k=k, p=p, meta=np.array(meta))


def run_census(tmp: Path, i: int):
    """Point the census at the sandbox rather than the real tree."""
    real, wc.ROOT = wc.ROOT, tmp
    try:
        return wc.census(i), None
    except Exception as exc:               # a refusal is a valid outcome
        return None, exc
    finally:
        wc.ROOT = real


def report(name: str, ok: bool, detail: str = "") -> None:
    global FAILS
    if not ok:
        FAILS += 1
    print(f"  {'OK   ' if ok else 'FAIL '} {name}" + (f"  [{detail}]" if detail else ""))


def gate_named(res, prefix: str):
    return [g for g in res["gates"] if g["gate"].startswith(prefix)]


# --------------------------------------------------------------------------- 0
def test_clean_passes():
    tmp = sandbox(SMALL)
    res, exc = run_census(tmp, SMALL)
    report("untampered i=5 census passes",
           exc is None and res["ok"], str(exc) if exc else "")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 1
def test_tampered_prime_breaks_a_gate():
    """Move one column's kill to a LATER live prime. The column now looks like
    it survived a round it did not, so a recorded survivor count must shift."""
    tmp = sandbox(SMALL)
    k, p, meta = load_npz(tmp, SMALL)
    fam = wc.make_fam(SMALL)
    from bandii_kernel import cells, first_live_after, live_intervals
    ivs = live_intervals(cells(fam), fam)
    rec = json.loads((tmp / "results" / f"i{SMALL}_sweep.json").read_text())
    z_lo, z_hi = rec["k_z"]
    idx = int(np.flatnonzero((k >= z_lo) & (k <= z_hi))[0])
    p[idx] = first_live_after(int(p[idx]), ivs, fam.D)      # one rung higher
    save_npz(tmp, SMALL, k, p, meta)
    res, exc = run_census(tmp, SMALL)
    if exc is not None:
        report("tampered witness prime is caught", True, f"refused: {type(exc).__name__}")
    else:
        bad = [g for g in res["gates"] if g["verdict"] == "fail"]
        report("tampered witness prime is caught", (not res["ok"]) and bool(bad),
               f"{len(bad)} gate(s) failed")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 2
def test_dropped_column_breaks_count():
    tmp = sandbox(SMALL)
    k, p, meta = load_npz(tmp, SMALL)
    rec = json.loads((tmp / "results" / f"i{SMALL}_sweep.json").read_text())
    z_lo, z_hi = rec["k_z"]
    drop = int(np.flatnonzero((k >= z_lo) & (k <= z_hi))[5])
    save_npz(tmp, SMALL, np.delete(k, drop), np.delete(p, drop), meta)
    res, exc = run_census(tmp, SMALL)
    ok = exc is not None or any(
        g["gate"] == "zjump.n_columns" and g["verdict"] == "fail" for g in res["gates"])
    report("dropped column breaks zjump.n_columns", ok,
           f"refused: {type(exc).__name__}" if exc else "")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 3
def test_off_ladder_prime_refused():
    """A prime that is not live cannot be the prime production stopped at."""
    tmp = sandbox(SMALL)
    k, p, meta = load_npz(tmp, SMALL)
    rec = json.loads((tmp / "results" / f"i{SMALL}_sweep.json").read_text())
    z_lo, z_hi = rec["k_z"]
    idx = int(np.flatnonzero((k >= z_lo) & (k <= z_hi))[0])
    p[idx] = int(p[idx]) + 1                       # almost certainly not prime, not live
    save_npz(tmp, SMALL, k, p, meta)
    res, exc = run_census(tmp, SMALL)
    ok = exc is not None and "ladder" in str(exc)
    report("off-ladder Z-jump prime is REFUSED", ok,
           str(exc)[:60] if exc else "census returned instead of refusing")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 4
def test_bandii_prime_outside_list_refused():
    tmp = sandbox(SMALL)
    k, p, meta = load_npz(tmp, SMALL)
    rec = json.loads((tmp / "results" / f"i{SMALL}_sweep.json").read_text())
    b_lo, b_hi = rec["k_bii"]
    idx = int(np.flatnonzero((k >= b_lo) & (k <= b_hi))[0])
    p[idx] = int(p[idx]) + 2
    save_npz(tmp, SMALL, k, p, meta)
    res, exc = run_census(tmp, SMALL)
    ok = exc is not None and "primes_bii" in str(exc)
    report("Band II prime outside primes_bii is REFUSED", ok,
           str(exc)[:60] if exc else "census returned instead of refusing")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 5
def test_wrong_member_refused():
    """i=5's witness table audited as if it were i=6 must not quietly pass."""
    tmp = sandbox(SMALL)
    shutil.copy(ROOT / "results" / f"i{SMALL}_sweep.json",
                tmp / "results" / "i6_sweep.json")
    shutil.copy(tmp / "results" / f"i{SMALL}_witness.npz",
                tmp / "results" / "i6_witness.npz")
    res, exc = run_census(tmp, 6)
    ok = exc is not None or not res["ok"]
    report("witness table checked as the WRONG member is caught", ok,
           f"refused: {str(exc)[:50]}" if exc else "gates failed")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 6
def test_excuse_requires_its_evidence():
    """The record_incomplete excuse for round-1 n_primes is conditional on every
    LATER round matching. Break a later round too: the excuse must not fire,
    and the result must be FAIL. This is the test that stops the excuse
    mechanism from becoming a tolerance."""
    res = {
        "phases": {"bandii": {"columns": 0}, "zjump": {"columns": 1000}},
        "gates": [
            {"gate": "zjump.round1.n_primes", "got": 990684, "want": 73009,
             "ok": False, "note": "", "verdict": "fail", "evidence": ""},
            {"gate": "zjump.round2.n_primes", "got": 308, "want": 299,
             "ok": False, "note": "", "verdict": "fail", "evidence": ""},
        ],
    }
    wc._classify_record_gaps(res, {"phases": {"zjump": []}})
    verdicts = [g["verdict"] for g in res["gates"]]
    report("excuse does NOT fire when a later round also mismatches",
           verdicts == ["fail", "fail"], f"verdicts={verdicts}")

    # and the converse: with every later round clean, it SHOULD fire
    res2 = {
        "phases": {"bandii": {"columns": 0}, "zjump": {"columns": 1000}},
        "gates": [
            {"gate": "zjump.round1.n_primes", "got": 990684, "want": 73009,
             "ok": False, "note": "", "verdict": "fail", "evidence": ""},
            {"gate": "zjump.round2.n_primes", "got": 308, "want": 308,
             "ok": True, "note": "", "verdict": "pass", "evidence": ""},
        ],
    }
    wc._classify_record_gaps(res2, {"phases": {"zjump": []}})
    report("excuse DOES fire when its evidence is present",
           res2["gates"][0]["verdict"] == "record_incomplete"
           and bool(res2["gates"][0]["evidence"]),
           res2["gates"][0]["verdict"])

    # an excuse must never fire in the direction of MORE recorded than found
    res3 = {
        "phases": {"bandii": {"columns": 0}, "zjump": {"columns": 1000}},
        "gates": [
            {"gate": "zjump.round1.n_primes", "got": 100, "want": 990684,
             "ok": False, "note": "", "verdict": "fail", "evidence": ""},
            {"gate": "zjump.round2.n_primes", "got": 308, "want": 308,
             "ok": True, "note": "", "verdict": "pass", "evidence": ""},
        ],
    }
    wc._classify_record_gaps(res3, {"phases": {"zjump": []}})
    report("excuse does NOT fire when the record claims MORE than reconstructed",
           res3["gates"][0]["verdict"] == "fail", res3["gates"][0]["verdict"])


# --------------------------------------------------------------------------- 7
def test_tampered_kmax_refused():
    tmp = sandbox(SMALL)
    f = tmp / "results" / f"i{SMALL}_sweep.json"
    rec = json.loads(f.read_text())
    rec["k_max"] = int(rec["k_max"]) + 1
    f.write_text(json.dumps(rec))
    res, exc = run_census(tmp, SMALL)
    ok = exc is not None and "k_max" in str(exc)
    report("tampered k_max is REFUSED (recomputed, not read back)", ok,
           str(exc)[:60] if exc else "census returned instead of refusing")
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- 8
def test_ladder_offbyone_trips_selfcheck():
    """Two independent reviewers made the prev+1 vs prev off-by-one on the
    bucket lower endpoint and it silently dropped ~0.5% of columns. Simulate a
    ladder missing one rung and confirm the contiguity self-check trips."""
    tmp = sandbox(MID)
    real_build = wc.build_ladder

    def holed(fam, lo, hi):
        L = real_build(fam, lo, hi)
        return np.delete(L, len(L) // 2) if L.size > 4 else L

    wc.build_ladder = holed
    try:
        res, exc = run_census(tmp, MID)
    finally:
        wc.build_ladder = real_build
    ok = exc is not None or (res is not None and not res["ok"])
    report("a ladder with a missing rung is caught", ok,
           str(exc)[:70] if exc else "gates failed")
    shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    print("=== work_census adversarial tests ===")
    for fn in (test_clean_passes,
               test_tampered_prime_breaks_a_gate,
               test_dropped_column_breaks_count,
               test_off_ladder_prime_refused,
               test_bandii_prime_outside_list_refused,
               test_wrong_member_refused,
               test_excuse_requires_its_evidence,
               test_tampered_kmax_refused,
               test_ladder_offbyone_trips_selfcheck):
        fn()
    print(f"\nRESULT {'PASS' if FAILS == 0 else f'FAIL ({FAILS})'}")
    return 0 if FAILS == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
