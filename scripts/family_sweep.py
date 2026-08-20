#!/usr/bin/env python3
"""Replay Band II + Z-jump on a Fibonacci family member. No giant m.

Default i=7. Exact i=7 already scanned k<=200; this job does
k=201..k_max except {K, K+1}.

Not a next-prime sweep from k through (k, N/2]. Cap 14 Band II,
cap 12 Z-jump live primes. Refuse if results/i{i}_sweep.json exists.
"""

from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2

from sizelaw import RoundLedger  # noqa: E402

from bandii_kernel import (  # noqa: E402
    Fam,
    append_jsonl as write_jsonl,
    cells,
    chunk_ks,
    equal_g_chunks,
    check_checkpoint,
    checkpoint_identity,
    first_live_after,
    first_primes_above,
    live_intervals,
    read_jsonl as load_done,
    summarize_survivors as summarize,
    fact_table,
    kmax_of,
    make_fam,
    r_closed,
    r_from_F,
    r_two_digit,
    scan_columns,
    scan_columns_general,
)

CAP_BII = 14
CAP_Z = 12
# Regime-aware extension of the Z-jump cap. Small columns are the HARDEST in
# this whole programme, not the easiest: at k < 10^3 most nearby primes are
# dead (Kummer carry), so the first LIVE prime sits far above k, g/p runs to
# 0.94 and the per-prime survival rate to 0.61. Runs of 8-10 are ordinary
# there. A uniform 12 is tight at the bottom and generous everywhere above,
# so give the bottom a few more rounds rather than reporting a fat-tail
# survivor as an anomaly. See sizelaw.py and docs/zjump-spec.md.
CAP_Z_SMALL_K = 15
SMALL_K = 10**3
N_CHUNKS = 32
DEFAULT_WORKERS = 8
K_EXACT = {2: 200, 3: 200, 4: 200, 5: 200, 6: 200, 7: 200, 9: 80}
# i=2..7: exact k_extra=200 in fibonacci_i1-7.json
# i=9: modular k<=80 all impossible. i=1 is 3003 (N=8), skip.
# i=8 closed by the dedicated pipeline.



def _job(payload: tuple) -> dict:
    kind, p, ks, N, K, r_expected = payload
    t0 = time.time()
    if kind == "bii":
        r, rows = scan_columns(p, ks, r_expected=r_expected, N=N, K=K)
    else:
        r, rows = scan_columns_general(p, ks, N=N, K=K)
    return {
        "p": p,
        "r": r,
        "k_lo": int(ks[0]),
        "k_hi": int(ks[-1]),
        "n_cols": len(ks),
        "n_survivors": len(rows),
        "survivors": rows,
        "seconds": round(time.time() - t0, 3),
    }


def _phase_count(done: list[dict], phase: str, field: str) -> int:
    """A counter off the last phase_complete record for `phase`.

    phase_complete means the phase ran to its cap, not that it ended empty.
    Resuming must read these counts back; assuming zero turns an unfinished
    phase into a clean certificate.
    """
    n = 0
    for rec in done:
        if rec.get("event") == "phase_complete" and rec.get("phase") == phase:
            n = int(rec.get(field) or 0)
    return n


def paths(i: int) -> tuple[Path, Path]:
    return ROOT / "results" / f"i{i}_sweep.json", ROOT / "results" / f"i{i}_sweep.jsonl"


def run_jobs(jobs: list[tuple], workers: int, chk: Path, tag: str, done_keys: set) -> list[dict]:
    pending = []
    for job in jobs:
        _kind, p, ks, _N, _K, _r = job
        key = (tag, p, int(ks[0]), int(ks[-1]))
        if key not in done_keys:
            pending.append(job)
    surv = []
    for rec in load_done(chk):
        if rec.get("tag") == tag:
            surv.extend(rec.get("survivors") or [])
    if not pending:
        print(f"  {tag} all chunks done  alive={len(surv)}", flush=True)
        return surv
    print(f"  {tag} jobs={len(pending)} workers={workers}", flush=True)
    nprint = 0
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        for rec in pool.imap_unordered(_job, pending):
            rec["tag"] = tag
            write_jsonl(chk, rec)
            surv.extend(rec["survivors"])
            nprint += 1
            fat = rec["n_cols"] >= 500 or rec["seconds"] >= 1.0
            if fat or nprint % 50 == 0:
                print(
                    f"    p={rec['p']} cols={rec['n_cols']} surv={rec['n_survivors']} "
                    f"{rec['seconds']}s  alive={len(surv)}",
                    flush=True,
                )
    return surv


def check(cond: bool, msg: str) -> None:
    """Pre-flight guard. Not assert: `python -O` strips assert, and every
    certificate downstream of this function is claimed unconditional."""
    if not cond:
        raise RuntimeError(f"pre-flight failed: {msg}")


def preflight(fam: Fam, kmax: int, primes: list[int]) -> list[tuple[int, int]]:
    print(f"=== i={fam.i} pre-flight ===", flush=True)
    t0 = time.time()
    check(fam.K != int(gmpy2.fib(2 * fam.i + 2) * gmpy2.fib(2 * fam.i + 1)), "K is F_{2i} F_{2i+3}")
    check(fam.N == int(gmpy2.fib(2 * fam.i + 2) * gmpy2.fib(2 * fam.i + 3)), "N = F_{2i+2} F_{2i+3}")
    check(fam.K == int(gmpy2.fib(2 * fam.i) * gmpy2.fib(2 * fam.i + 3)), "K = F_{2i} F_{2i+3}")
    check(fam.K1 == fam.K + 1, "K1 = K+1")
    for p in primes:
        check(bool(gmpy2.is_prime(p)), f"{p} is prime")
        check(2 * p > fam.N and p <= fam.D and p > kmax, f"{p} in Band II live window")
    print(
        f"  N={fam.N} K={fam.K} d={fam.D} N/2={fam.N2} kmax={kmax} "
        f"p1={primes[0]} dlt={2*primes[0]-fam.N}",
        flush=True,
    )

    F = fact_table(primes[0])
    ra = r_from_F(F, primes[0], N=fam.N, K=fam.K)
    rc = r_closed(primes[0], N=fam.N, K=fam.K)
    rt = r_two_digit(F, primes[0], N=fam.N, K=fam.K)
    if not (ra == rc == rt):
        raise RuntimeError(f"r(p1) {ra} {rc} {rt}")
    if int(F[-1]) != primes[0] - 1:
        raise RuntimeError("Wilson")
    print(f"  r(p1)={ra} table=closed=two-digit  Wilson ok", flush=True)

    from singmaster_intersect import binom_mod_lucas

    if int(binom_mod_lucas(fam.N, fam.K, primes[0])) != ra:
        raise RuntimeError("lucas(p1) mismatch")

    import numpy as np

    rng = np.random.default_rng(1)
    bad = 0
    for p in (11, 29, 101, 211, 1009):
        Fp = fact_table(p)
        check(int(Fp[p - 1]) == p - 1, f"Wilson at p={p}")
        for _ in range(100):
            k = int(rng.integers(1, p))
            rr = int(rng.integers(0, p))
            s = rr * int(Fp[k]) % p
            ker = bool(np.any(Fp[k:p] == (np.int64(s) * Fp[: p - k]) % np.int64(p)))
            brute = any(math.comb(n0, k) % p == rr for n0 in range(k, p))
            if ker != brute:
                bad += 1
    if bad:
        raise RuntimeError(f"kernel mismatches {bad}")
    print(f"  kernel 500 cases, 0 mismatches  {time.time()-t0:.1f}s", flush=True)

    windows = cells(fam)
    ivs = live_intervals(windows, fam)
    k0 = K_EXACT.get(fam.i, 2) + 1
    if k0 < fam.K:
        p0 = first_live_after(k0, ivs, fam.D)
        if p0 is None or p0 <= k0:
            raise RuntimeError(f"no live prime after {k0}")
        print(f"  live after k={k0}: {p0}  intervals={len(ivs)}", flush=True)
    else:
        print(f"  Band I extra already in exact k<= {k0-1}; Z-jump empty", flush=True)
    print("=== pre-flight passed ===", flush=True)
    return ivs


def main() -> int:
    i = 7
    if "--i" in sys.argv:
        i = int(sys.argv[sys.argv.index("--i") + 1])
    # i=8 was originally closed by a stitched pipeline (Stages 1-3 by
    # next-prime walk, hang-guards, stragglers, Band II, Z-jump remnant) whose
    # checkpoints are gone, so its 4.27M certificates exist nowhere. Running
    # it here re-derives them by ONE uniform method over k=3..k_max, which is
    # a better proof object than the original stitching. Outputs are
    # i8_sweep.json / i8_witness.npz and do not collide with the historical
    # bandii_sweep.json or zjump.json.

    out, chk = paths(i)
    if "--preflight" in sys.argv:
        fam = make_fam(i)
        kmax, _ = kmax_of(fam)
        primes = first_primes_above(fam.N2, fam.D, kmax)
        preflight(fam, kmax, primes)
        return 0

    if out.exists():
        print(f"{out} already exists. Not rerunning.", flush=True)
        return 2
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        print("numpy is required.", flush=True)
        return 1

    workers = int(os.environ.get("I7_WORKERS", os.environ.get("FAMILY_WORKERS", DEFAULT_WORKERS)))
    workers = max(1, min(workers, 16))
    out.parent.mkdir(exist_ok=True)

    t0 = time.time()
    fam = make_fam(i)
    kmax, logm = kmax_of(fam)
    primes = first_primes_above(fam.N2, fam.D, kmax, n=max(CAP_BII, 16))
    if not primes:
        raise RuntimeError("no live Band II primes in (N/2, d]")
    ivs = preflight(fam, kmax, primes)
    k_lo_z = K_EXACT.get(i, 2) + 1
    n_z = max(0, fam.K - k_lo_z)  # k_lo_z .. K-1, or 0 if exact already covers
    n_bii = kmax - (fam.K + 2) + 1
    print(
        f"=== i={i} sweep  Z {k_lo_z}..{fam.K-1} ({n_z})  "
        f"BII {fam.K+2}..{kmax} ({n_bii})  workers={workers} ===",
        flush=True,
    )
    print(f"    log10 m={logm:.4f}  p1={primes[0]}", flush=True)

    # cap_z alone does not describe the loop: the Z-jump runs to
    # CAP_Z_SMALL_K when k_lo_z < SMALL_K (which is every member, since the
    # exact band never reaches 1000), so a header claiming cap_z=12 for a run
    # that tested up to 15 rounds under-describes it. Record both. Adding a
    # key is backward compatible -- check_checkpoint only compares keys that
    # are present in the recorded header, so pre-existing headers still pass.
    ident = dict(i=i, N=fam.N, K=fam.K, k_max=kmax, k_lo_z=k_lo_z,
                 cap_bii=CAP_BII, cap_z=CAP_Z,
                 cap_z_small_k=CAP_Z_SMALL_K, small_k=SMALL_K)
    check_checkpoint(chk, **ident)
    if not chk.exists() or chk.stat().st_size == 0:
        write_jsonl(chk, checkpoint_identity(**ident))
    done = load_done(chk)
    done_keys = {
        (r["tag"], r["p"], r["k_lo"], r["k_hi"])
        for r in done
        if "tag" in r and "p" in r and "k_lo" in r
    }
    complete = {r["phase"] for r in done if r.get("event") == "phase_complete"}

    phases = {}

    # --- Band II ---
    ledger_bii = RoundLedger("bandii")
    ledger_z = RoundLedger("zjump")
    if "bandii" not in complete:
        alive = None
        for pi, p in enumerate(primes[:CAP_BII], start=1):
            tag = f"bii{pi}"
            if alive is None:
                if n_bii < 2000:
                    jobs = [
                        (
                            "bii",
                            p,
                            list(range(fam.K + 2, kmax + 1)),
                            fam.N,
                            fam.K,
                            r_closed(p, N=fam.N, K=fam.K),
                        )
                    ]
                else:
                    chunks = equal_g_chunks(fam.K + 2, kmax, p, N_CHUNKS)
                    jobs = [
                        ("bii", p, list(range(lo, hi + 1)), fam.N, fam.K, r_closed(p, N=fam.N, K=fam.K))
                        for lo, hi in chunks
                    ]
            else:
                if not alive:
                    break
                ks_alive = [row["k"] for row in alive]
                r_p = r_closed(p, N=fam.N, K=fam.K)
                jobs = [
                    ("bii", p, ch, fam.N, fam.K, r_p)
                    for ch in chunk_ks(ks_alive, p, N_CHUNKS if len(ks_alive) >= 2000 else 1)
                ]
            entered = (
                list(range(fam.K + 2, kmax + 1)) if alive is None
                else [row["k"] for row in alive]
            )
            exp = ledger_bii.expect((k_, p) for k_ in entered)
            surv = run_jobs(jobs, workers, chk, tag, done_keys)
            sm = summarize(surv)
            judged = ledger_bii.record(pi, len(entered), exp, sm["n"])
            phases.setdefault("bandii", []).append(
                {"prime_index": pi, "p": p, **sm, "expected": round(exp, 4),
                 "escalate": judged["escalate"]}
            )
            if judged["escalate"]:
                print(f"  ** ESCALATE ** BII pass {pi}: {judged['reason']} "
                      f"(expected {exp:.3g}, observed {sm['n']})", flush=True)
            print(f"  BII pass {pi} p={p} alive={sm['n']} (size law expected "
                  f"{exp:.4g}) even={sm['even']} mean_k={sm['mean_k']}", flush=True)
            alive = surv
            write_jsonl(chk, {"event": "round_complete", "phase": "bandii", "pass": pi, "n_alive": sm["n"]})
            done_keys = {
                (r["tag"], r["p"], r["k_lo"], r["k_hi"])
                for r in load_done(chk)
                if "tag" in r and "p" in r and "k_lo" in r
            }
            if not alive:
                break
        write_jsonl(chk, {"event": "phase_complete", "phase": "bandii", "n_alive": 0 if not alive else len(alive)})
        bii_left = [] if not alive else alive
        n_bii_left = len(bii_left)
    else:
        # Band II can exhaust CAP_BII with survivors and still write
        # phase_complete, so read the recorded count rather than assuming 0.
        bii_left = []
        n_bii_left = _phase_count(done, "bandii", "n_alive")
        print(f"  bandii phase already complete (recorded alive={n_bii_left})", flush=True)

    # --- Z-jump Band I remnant ---
    cap_z = CAP_Z          # bound even when the phase is skipped or resumed
    if "zjump" not in complete and n_z > 0:
        current: list = [{"k": k} for k in range(k_lo_z, fam.K)]
        zrounds = []
        z_none: list[int] = []
        cap_z = CAP_Z_SMALL_K if k_lo_z < SMALL_K else CAP_Z
        for rnd in range(1, cap_z + 1):
            if not current:
                break
            if rnd > CAP_Z:
                # extra rounds serve the small-k tail only
                current = [c for c in current if int(c["k"]) < SMALL_K]
                if not current:
                    break
            buckets = defaultdict(list)
            none = []
            for it in current:
                if isinstance(it, dict) and "g" in it:
                    k, x = int(it["k"]), int(it["k"] + it["g"])
                else:
                    k = int(it["k"] if isinstance(it, dict) else it)
                    x = k
                p = first_live_after(x, ivs, fam.D)
                if p is None:
                    none.append(k)
                else:
                    buckets[p].append(k)
            if none:
                # A column with no remaining live prime was never tested. It is
                # NOT killed -- carry it to the anomaly list so `clean` is false.
                print(f"  Z round {rnd} no live prime: {len(none)} (anomaly)", flush=True)
                z_none.extend(none)
            jobs = []
            for p, ks in sorted(buckets.items()):
                for ch in chunk_ks(ks, p, N_CHUNKS if len(ks) >= 2000 else 1):
                    jobs.append(("z", p, ch, fam.N, fam.K, None))
            exp = ledger_z.expect(
                (k_, p_) for p_, ks_ in buckets.items() for k_ in ks_
            )
            n_entered = sum(len(v) for v in buckets.values())
            surv = run_jobs(jobs, workers, chk, f"z{rnd}", done_keys)
            sm = summarize(surv)
            judged = ledger_z.record(rnd, n_entered, exp, sm["n"])
            zrounds.append({"round": rnd, "n_primes": len(buckets), **sm,
                            "n_nolive": len(none), "expected": round(exp, 4),
                            "escalate": judged["escalate"]})
            print(
                f"  Z round {rnd} alive={sm['n']} (size law expected {exp:.4g}) "
                f"even={sm['even']} mean_k={sm['mean_k']}",
                flush=True,
            )
            if judged["escalate"]:
                print(f"  ** ESCALATE ** Z round {rnd}: {judged['reason']} "
                      f"(expected {exp:.3g}, observed {sm['n']})", flush=True)
            current = surv
            write_jsonl(chk, {"event": "round_complete", "phase": "zjump", "round": rnd, "n_alive": sm["n"]})
            done_keys = {
                (r["tag"], r["p"], r["k_lo"], r["k_hi"])
                for r in load_done(chk)
                if "tag" in r and "p" in r and "k_lo" in r
            }
        write_jsonl(
            chk,
            {
                "event": "phase_complete",
                "phase": "zjump",
                "n_alive": len(current),
                "n_nolive": len(z_none),
            },
        )
        phases["zjump"] = zrounds
        z_left = current
        n_z_left = len(z_left)
        n_z_nolive = len(z_none)
    elif "zjump" in complete:
        z_left, z_none = [], []
        n_z_left = _phase_count(done, "zjump", "n_alive")
        n_z_nolive = _phase_count(done, "zjump", "n_nolive")
        print(
            f"  zjump phase already complete "
            f"(recorded alive={n_z_left} nolive={n_z_nolive})",
            flush=True,
        )
    else:
        print("  zjump skipped (exact already covers Band I extra)", flush=True)
        z_left, z_none = [], []
        n_z_left = n_z_nolive = 0
        write_jsonl(chk, {"event": "phase_complete", "phase": "zjump", "skipped": True})

    # A sweep is clean only if nothing survived AND every column was testable.
    # A column with no live prime left was never tested; it is not a kill.
    clean = n_bii_left == 0 and n_z_left == 0 and n_z_nolive == 0

    # The Z-jump starts at K_EXACT[i]+1, so the small-k columns below it are
    # closed by `modular` and their run statistic is never logged. That is the
    # regime with the LONGEST runs (i=9: k=11 run 8, k=29 run 7, k=45 run 5),
    # so a readout that omits it reports the Band II maximum as if it were the
    # maximum. Census it here, cheaply, and record it alongside.
    smallk = None
    if k_lo_z > 2:
        try:
            from sizelaw import assess, live_run, run_lambda

            rows = []
            for kk in range(2, min(k_lo_z, 400)):
                if kk in (fam.K, fam.K + 1):
                    continue
                surv, kill, dead = live_run(fam.N, fam.K, kk, cap=20)
                if surv:
                    lam = run_lambda(kk, surv)
                    rows.append(assess(kk, surv,
                                       expected=(k_lo_z - 2) * lam, observed=1))
            rows.sort(key=lambda r: r["expected"])
            longest = max(rows, key=lambda r: r["run"]) if rows else None
            smallk = {
                "k_range": [2, min(k_lo_z, 400) - 1],
                "n_with_run": len(rows),
                "longest_run": longest and {"k": longest["k"], "run": longest["run"],
                                            "expected": longest["expected"]},
                "most_surprising": rows[0] if rows else None,
                "escalate": any(r["escalate"] for r in rows),
                "note": "closed by `modular`, not by the Z-jump; run length here "
                        "is NOT comparable with Band II -- see sizelaw.py",
            }
            if longest:
                print(f"  small-k census k<{min(k_lo_z,400)}: longest run "
                      f"k={longest['k']} run={longest['run']} "
                      f"(expected {longest['expected']:.3g}); most surprising "
                      f"k={rows[0]['k']} run={rows[0]['run']}", flush=True)
            if smallk["escalate"]:
                print("  ** ESCALATE ** small-k census", flush=True)
        except Exception as exc:
            print(f"  small-k census failed: {exc!r}", flush=True)
            smallk = {"error": repr(exc)}

    # Emit the per-column witness table while the checkpoint still exists.
    # Without this the certificate below asserts a fact about every column and
    # records the witness for none of them, and the only way to re-check it is
    # to run the whole sweep again.
    witness_meta = None
    if clean:
        try:
            import witness as _witness

            wpath = ROOT / "results" / f"i{i}_witness.npz"
            witness_meta = _witness._build_family(i, chk, wpath)
            print(
                f"  witnesses {witness_meta['n_witnesses']} -> {wpath.name}  "
                f"sha256 {witness_meta['sha256'][:16]}...",
                flush=True,
            )
            if witness_meta["n_unresolved"]:
                print(
                    f"  WARNING {witness_meta['n_unresolved']} columns have no "
                    f"witness; certificate withheld",
                    flush=True,
                )
                clean = False
        except Exception as exc:  # never lose a finished sweep to this
            print(f"  witness build failed: {exc!r}", flush=True)
            witness_meta = {"error": repr(exc)}
    payload = {
        "search": f"i{i}_sweep",
        "i": i,
        "N": fam.N,
        "K": fam.K,
        "d": fam.D,
        "k_max": kmax,
        "log10_m": round(logm, 6),
        "k_z": [k_lo_z, fam.K - 1],
        "k_bii": [fam.K + 2, kmax],
        "n_z": n_z,
        "cap_z": cap_z,
        "cap_z_small_k": CAP_Z_SMALL_K,
        "n_bii": n_bii,
        "primes_bii": primes[:CAP_BII],
        "workers": workers,
        "phases": phases,
        "n_bii_alive": n_bii_left,
        "n_z_alive": n_z_left,
        "n_z_nolive": n_z_nolive,
        "z_nolive": z_none[:100],
        "witness": witness_meta,
        "small_k_census": smallk,
        "escalation": {
            "bandii": ledger_bii.verdict(),
            "zjump": ledger_z.verdict(),
        },
        "bii_survivors": bii_left if n_bii_left and n_bii_left <= 100 else [],
        "z_survivors": z_left if n_z_left and n_z_left <= 100 else [],
        "clean": clean,
        "certificate": (
            f"Every extra k in [2, k_max] except {{K, K+1}} for i={i} "
            f"has r(p) notin I_{{p,k}}. Together with exact k<=200 and "
            f"two family columns, N(C(N,K))=6. Not Singmaster. "
            f"Witness prime for every column: results/i{i}_witness.npz "
            f"(sha256 {witness_meta['sha256'][:16] if witness_meta and 'sha256' in witness_meta else 'n/a'}); "
            f"re-check with: python scripts/witness.py verify "
            f"--file results/i{i}_witness.npz"
        )
        if clean
        else None,
        "seconds": round(time.time() - t0, 3),
    }
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(flush=True)
    print(
        f"wrote {out}  clean={clean}  bii_left={n_bii_left}  "
        f"z_left={n_z_left}  {payload['seconds']}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
