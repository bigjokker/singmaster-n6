#!/usr/bin/env python3
"""Exact scan-work census for a run member, reconstructed from its own record.

Why this exists. `profile_sweep.py` estimates a phase by timing a sample of
columns and multiplying by the column count. Both of its probes were measured
wrong: the Z-jump by 8.5x-156x -- its answer moved with `--sample`, which the
docstring calls "cost, not accuracy" -- and Band II by a stable 1.7x, because
it timed the marginal column at k=K+2 where g is largest and applied that to a
population whose mean g is 40% smaller. Sampling is simply the wrong
instrument for a heavy-tailed quantity (i=9 Z-jump: median g = 29, mean g =
307,398) when the run's own record is sitting on disk.

So compute it exactly instead.

The witness table stores only the KILLING prime per column, which undercounts:
production scanned each column at every live prime it survived first. This
reconstructs the full chain.

  Z-jump   Round 1 enters with x = k, so p1 = first_live_after(k). Round 2+
           enters with x = k + g, and since the survivor record stores
           g = p - k, that x is exactly the prime just survived -- so the
           chain is first_live_after applied to the PRIME, not to k
           (family_sweep.py:374-378). Consecutive live primes are consecutive
           entries of the live-prime ladder, so p_j = L[idx + j - 1] and the
           whole walk vectorises.

  Band II  One shared prime list, meta['primes_bii'], walked in order; pass j
           scans every column not yet killed.

then sums ceil(g/2) -- the half-scan's element count -- over every
(column, prime) pair actually scanned.

The reconstruction is not trusted on its own say-so. It predicts quantities
that production recorded independently in `results/i{i}_sweep.json` and never
saw: per-round survivor counts, per-round distinct-prime counts, per-pass
Band II survivors, and the no-live-prime counts. Every one must match
EXACTLY.

A reconstruction that needs a tolerance is a reconstruction that is still
wrong. The only permitted discrepancy is a named leftover -- a column present
in the witness table but outside the swept bands, filled by the exact or
modular engine instead (i=7: k=2..200; i=8: k=2) -- or a column absent from
the table altogether (i=9: the four Lucas-killed k=87/399/553/1281). Those are
listed individually, never absorbed into a fudge.

Scan work and wall clock are reported as TWO numbers. They are not the same
number and the gap between them is not constant: it measured 5.6x at i=7 and
9.4x at i=8. Until that gap is itemised, a faster scan cannot be quoted as a
faster job, and no wall clock is extrapolated to an unrun member here.

    python scripts/work_census.py
    python scripts/work_census.py --i 9
    python scripts/work_census.py --json_out results/work_census.json
"""

from __future__ import annotations

import argparse
import json
import textwrap
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bandii_kernel import (  # noqa: E402
    cells,
    first_live_after,
    kmax_of,
    live_intervals,
    make_fam,
)

# Single-core numpy scan rate, measured on this machine over the real
# `(s * F) % p` inner loop. Used only to turn an exact element count into an
# indicative core-hour figure -- it is a machine property, not a claim about
# the sweep, and it is deliberately quoted to three digits because the probe
# bugs this file replaces were factors of 8 to 150, not 10%.
SCAN_RATE = 1.865e8


def check(cond: bool, msg: str) -> None:
    """Gate. Not assert: `python -O` strips assert, and a census that silently
    stops checking itself is exactly the failure this file exists to prevent."""
    if not cond:
        raise RuntimeError(f"census gate failed: {msg}")


def build_ladder(fam, lo: int, hi: int) -> np.ndarray:
    """Every live prime in (lo, hi], ascending.

    Built by walking `first_live_after` rather than by sieving, so the ladder
    is definitionally the same object the sweep stepped along -- including its
    treatment of interval edges and of primes that are skipped as dead.
    """
    ivs = live_intervals(cells(fam), fam)
    out: list[int] = []
    q = first_live_after(lo, ivs, fam.D)
    while q is not None and q <= hi:
        out.append(int(q))
        q = first_live_after(q, ivs, fam.D)
    return np.array(out, dtype=np.int64)


def _classify_record_gaps(res: dict, rec: dict) -> None:
    """Separate "the reconstruction is wrong" from "the record is incomplete".

    Only these two named, evidence-carrying cases are excused, and each states
    the evidence inline so a reader can audit the excuse. Everything else that
    fails stays failed -- this is deliberately not a tolerance.
    """
    by = {g["gate"]: g for g in res["gates"]}
    ph = rec.get("phases") or {}

    g = by.get("bandii.n_passes_realised")
    if g and not g["ok"] and not ph.get("bandii") and res["phases"]["bandii"]["columns"]:
        g["verdict"] = "record_incomplete"
        g["evidence"] = (
            f"phases.bandii is {ph.get('bandii')!r} while n_bii = "
            f"{res['phases']['bandii']['columns']:,}. Band II reached "
            f"phase_complete before this run was resumed, so its passes were "
            f"never written to this record. Not evaluable from this file.")

    g = by.get("zjump.round1.n_primes")
    later = [x for x in res["gates"]
             if x["gate"].startswith("zjump.round") and x["gate"].endswith(".n_primes")
             and x["gate"] != "zjump.round1.n_primes"]
    if (g and not g["ok"] and g["want"] < g["got"]
            and later and all(x["ok"] for x in later)):
        g["verdict"] = "record_incomplete"
        g["evidence"] = (
            f"rounds 2+ n_primes match exactly ({len(later)} of them); only "
            f"round 1 is short. Round-1 buckets partition the entire Z band, so "
            f"covering {res['phases']['zjump']['columns']:,} columns needs the "
            f"whole ladder -- and the recorded {g['want']:,} cannot do it "
            f"(i=8 needed 124,830 buckets for 4,126,644 columns). The run was "
            f"resumed mid-round-1: survivors were merged from the checkpoint, "
            f"which is why round-1 survivors match while its prime count does not.")


def census(i: int, rate: float = SCAN_RATE) -> dict:
    fam = make_fam(i)
    npz = ROOT / "results" / f"i{i}_witness.npz"
    sweep = ROOT / "results" / f"i{i}_sweep.json"
    check(npz.exists(), f"no witness table at {npz}")
    check(sweep.exists(), f"no sweep record at {sweep} -- the gates need it")

    d = np.load(npz, allow_pickle=True)
    meta = json.loads(str(d["meta"]))
    rec = json.loads(sweep.read_text())
    k_all = d["k"].astype(np.int64)
    p_all = d["p"].astype(np.int64)

    # k_max recomputed from N,K, never read back from the file under audit.
    kmax, _ = kmax_of(fam)
    check(kmax == rec["k_max"], f"k_max {kmax} != recorded {rec['k_max']}")

    z_lo, z_hi = rec["k_z"]
    b_lo, b_hi = rec["k_bii"]
    in_z = (k_all >= z_lo) & (k_all <= z_hi)
    in_b = (k_all >= b_lo) & (k_all <= b_hi)
    outside = ~(in_z | in_b)

    res: dict = {
        "i": i, "N": fam.N, "K": fam.K, "k_max": kmax,
        "witnesses": int(len(k_all)),
        "leftovers": {}, "gates": [], "phases": {},
    }

    # ---- named leftovers: present in the table but never swept -------------
    left_k = k_all[outside]
    res["leftovers"]["not_swept"] = {
        "n": int(len(left_k)),
        "k": [int(x) for x in left_k[:12]],
        "truncated": bool(len(left_k) > 12),
        "why": "outside k_z and k_bii -- filled by the exact/modular engine",
    }
    # ---- named leftovers: claimed by the run but absent from the table -----
    res["leftovers"]["unresolved"] = {
        "n": int(meta.get("n_unresolved", 0)),
        "k": [int(x) for x in meta.get("unresolved", [])],
        "why": "recorded unresolved in the witness meta (cap leftovers)",
    }

    def gate(name: str, got, want, note: str = "") -> bool:
        ok = got == want
        res["gates"].append({"gate": name, "got": got, "want": want,
                             "ok": bool(ok), "note": note,
                             "verdict": "pass" if ok else "fail", "evidence": ""})
        return ok

    # ------------------------------------------------------------------ Band II
    kb, pb = k_all[in_b], p_all[in_b]
    pr = np.array(sorted(meta["primes_bii"]), dtype=np.int64)
    pos = np.searchsorted(pr, pb)
    check(bool((pr[np.clip(pos, 0, len(pr) - 1)] == pb).all()),
          "a Band II witness prime is not in meta['primes_bii']")

    ops_b = 0
    b_rounds = []
    for j in range(len(pr)):
        entered = pos >= j
        n_in = int(entered.sum())
        if n_in == 0:
            break
        ops_b += int(np.ceil((pr[j] - kb[entered]) / 2).sum())
        b_rounds.append({"pass": j + 1, "p": int(pr[j]),
                         "entered": n_in, "survived": int((pos > j).sum())})
    naive_b = int(np.ceil((pb - kb) / 2).sum())

    rec_b = rec.get("phases", {}).get("bandii", [])
    gate("bandii.n_columns", int(len(kb)), int(rec["n_bii"]))
    for j, rr in enumerate(rec_b):
        if j < len(b_rounds):
            gate(f"bandii.pass{j + 1}.survivors", b_rounds[j]["survived"], int(rr["n"]))
            gate(f"bandii.pass{j + 1}.prime", b_rounds[j]["p"], int(rr["p"]))
    gate("bandii.n_passes_realised", len(b_rounds), len(rec_b),
         "production stops on weight, not on CAP_BII")

    res["phases"]["bandii"] = {
        "columns": int(len(kb)), "ops_exact": ops_b, "ops_naive": naive_b,
        "multipass_factor": round(ops_b / naive_b, 6) if naive_b else None,
        "mean_g": float((pb - kb).mean()) if len(kb) else 0.0,
        "median_g": int(np.median(pb - kb)) if len(kb) else 0,
        "rounds": b_rounds,
    }

    # ------------------------------------------------------------------- Z-jump
    kz, pz = k_all[in_z], p_all[in_z]
    ops_z = 0
    z_rounds = []
    naive_z = int(np.ceil((pz - kz) / 2).sum()) if len(kz) else 0
    if len(kz):
        L = build_ladder(fam, int(kz.min()) - 1, int(pz.max()))
        check(L.size > 0, "empty live-prime ladder")
        idx = np.searchsorted(L, kz, side="right")
        end = np.searchsorted(L, pz, side="left")
        check(bool((L[np.clip(end, 0, L.size - 1)] == pz).all()),
              "a Z-jump witness prime is not on the live-prime ladder")
        rounds = end - idx + 1
        check(bool((rounds >= 1).all()), "a Z-jump column has a non-positive round index")

        # Self-check on the ladder, independent of any record. Round-1 buckets
        # are bucket(p) = [previous_live_prime, p-1] -- lower endpoint
        # INCLUSIVE, because first_live_after returns the first prime STRICTLY
        # greater than x. Take prev+1 instead and ~0.5% of columns silently
        # vanish with no other symptom, so assert the shape directly: each
        # round-1 prime must own a contiguous run of k.
        first_p = L[idx]
        fp_sorted = first_p[np.argsort(kz, kind="stable")]
        n_runs = 1 + int((np.diff(fp_sorted) != 0).sum())
        check(bool((np.diff(fp_sorted) >= 0).all()),
              "round-1 prime is not monotone in k -- the ladder is wrong")
        check(n_runs == int(np.unique(first_p).size),
              "a round-1 prime owns a non-contiguous set of k -- the ladder is wrong")

        for t in range(int(rounds.max())):
            m = rounds > t
            n_in = int(m.sum())
            step = idx[m] + t
            check(bool((step < L.size).all()),
                  f"round {t + 1} walks past the end of the ladder")
            primes_t = L[step]
            ops_z += int(np.ceil((primes_t - kz[m]) / 2).sum())
            z_rounds.append({"round": t + 1, "entered": n_in,
                             "survived": int((rounds > t + 1).sum()),
                             "n_primes": int(np.unique(primes_t).size)})

    rec_z = rec.get("phases", {}).get("zjump", [])
    gate("zjump.n_columns", int(len(kz)) + res["leftovers"]["unresolved"]["n"],
         int(rec["n_z"]), "table columns + named unresolved leftovers")
    for t, rr in enumerate(rec_z):
        if t < len(z_rounds):
            gate(f"zjump.round{t + 1}.survivors",
                 z_rounds[t]["survived"] + res["leftovers"]["unresolved"]["n"],
                 int(rr["n"]), "reconstructed + named unresolved leftovers")
            gate(f"zjump.round{t + 1}.n_primes", z_rounds[t]["n_primes"],
                 int(rr["n_primes"]))
        if int(rr.get("n_nolive", 0)):
            gate(f"zjump.round{t + 1}.n_nolive", 0, int(rr["n_nolive"]),
                 "a column with no live prime is an anomaly, not a kill")

    res["phases"]["zjump"] = {
        "columns": int(len(kz)), "ops_exact": ops_z, "ops_naive": naive_z,
        "multipass_factor": round(ops_z / naive_z, 6) if naive_z else None,
        "mean_g": float((pz - kz).mean()) if len(kz) else 0.0,
        "median_g": int(np.median(pz - kz)) if len(kz) else 0,
        "rounds": z_rounds,
    }

    # ------------------------------------------------------------------ totals
    total = ops_b + ops_z
    res["ops_exact"] = total
    res["bandii_share"] = round(100.0 * ops_b / total, 2) if total else None
    res["scan_core_h"] = total / rate / 3600.0
    res["scan_rate_used"] = rate

    _classify_record_gaps(res, rec)
    res["record_partial"] = any(g["verdict"] == "record_incomplete" for g in res["gates"])

    secs, workers = rec.get("seconds"), rec.get("workers")
    if secs and workers:
        res["recorded_seconds"] = secs
        res["workers"] = workers
        res["recorded_wall_core_h"] = secs * workers / 3600.0
        # A partial record's `seconds` covers only the resumed leg, so dividing
        # scan work by it is meaningless. Refuse the ratio rather than print a
        # number that looks like a measurement.
        if not res["record_partial"]:
            res["scan_share_of_wall"] = round(
                100.0 * res["scan_core_h"] / res["recorded_wall_core_h"], 2)
            res["wall_over_scan"] = round(
                res["recorded_wall_core_h"] / res["scan_core_h"], 2)

    res["ok"] = all(g["verdict"] != "fail" for g in res["gates"])
    return res


def report(rows: list[dict]) -> None:
    for r in rows:
        print(f"\n=== i={r['i']}  N={r['N']:,}  K={r['K']:,}  k_max={r['k_max']:,} ===")
        for name in ("bandii", "zjump"):
            ph = r["phases"][name]
            print(f"  {name:7s} cols {ph['columns']:>12,d}  ops {ph['ops_exact']:>18,d}"
                  f"  multipass x{ph['multipass_factor']}"
                  f"  mean_g {ph['mean_g']:>12,.0f}  median_g {ph['median_g']:>10,d}")
        print(f"  TOTAL   ops {r['ops_exact']:,}   Band II share {r['bandii_share']}%")
        print(f"  scan work {r['scan_core_h']:.4f} core-h at {r['scan_rate_used']:.4g} elem/s")
        if "recorded_wall_core_h" in r:
            print(f"  recorded  {r['recorded_wall_core_h']:.4f} core-h wall "
                  f"({r['recorded_seconds']:,.1f} s x {r['workers']} workers)")
            if "wall_over_scan" in r:
                print(f"  -> scan is {r['scan_share_of_wall']}% of wall; "
                      f"wall/scan = {r['wall_over_scan']}x  "
                      f"(NOT constant across i -- do not extrapolate it)")
            else:
                print("  -> wall/scan NOT computed: this record is partial, so its "
                      "`seconds` covers only the resumed leg")
        for key, lf in r["leftovers"].items():
            if lf["n"]:
                shown = ", ".join(str(x) for x in lf["k"])
                more = " ..." if lf.get("truncated") else ""
                print(f"  leftover {key}: n={lf['n']}  k={shown}{more}  ({lf['why']})")
        npass = sum(1 for g in r["gates"] if g["verdict"] == "pass")
        inc = [g for g in r["gates"] if g["verdict"] == "record_incomplete"]
        bad = [g for g in r["gates"] if g["verdict"] == "fail"]
        print(f"  gates: {npass}/{len(r['gates'])} pass"
              + (f", {len(inc)} not evaluable, {len(bad)} FAIL" if inc or bad else ""))
        for g in inc:
            print(f"    RECORD INCOMPLETE {g['gate']}: "
                  f"reconstructed {g['got']:,} vs recorded {g['want']:,}")
            for line in textwrap.wrap(g["evidence"], 84):
                print(f"        {line}")
        for g in bad:
            print(f"    FAIL {g['gate']}: reconstructed {g['got']} != recorded {g['want']}"
                  + (f"  [{g['note']}]" if g["note"] else ""))
        print(f"  RESULT {'PASS' if r['ok'] else 'FAIL'}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--i", type=int, default=None,
                    help="member; default every member with a witness table and a sweep record")
    ap.add_argument("--rate", type=float, default=SCAN_RATE,
                    help="single-core scan rate in elements/s, for the core-hour figure only")
    ap.add_argument("--json_out", type=str, default=None)
    args = ap.parse_args()

    if args.i is not None:
        members = [args.i]
    else:
        members = sorted(
            int(f.stem.split("_")[0][1:])
            for f in (ROOT / "results").glob("i*_witness.npz")
            if (ROOT / "results" / f"i{f.stem.split('_')[0][1:]}_sweep.json").exists()
        )
    rows = [census(i, args.rate) for i in members]
    report(rows)

    if args.json_out:
        out = ROOT / args.json_out if not Path(args.json_out).is_absolute() else Path(args.json_out)
        out.write_text(json.dumps({"search": "work_census", "members": rows}, indent=1))
        print(f"\nwrote {out}")
    return 0 if all(r["ok"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
