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
    iter_jsonl as stream_done,
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
    scan_ks_windowed,
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
# Re-derive done_keys from the whole checkpoint each round and assert it
# equals the incrementally maintained set. Costs a full parse per round,
# so it is for tests, not production.
VERIFY_DONE_KEYS = False
# r(p) for a WHOLE ROUND at once, by reducing m = C(N,K) against that round's
# primes with a product/remainder tree, in the PARENT.
#
# Band II primes have a tiny third Lucas index (2p-N-1), so r there is a
# handful of multiplies and is already computed parent-side. Z-jump primes sit
# near k, their remainders are generic, and min(k0, n0-k0, p-n0-1) stays O(p):
# measured over all 990,683 i=9 round-1 primes the median min-index is
# 1,032,763 and the Lucas bill is 23.82 core-h -- 99.7% of a job that scans a
# median of 12 columns.
#
# MEASURED at i=9 on the real ladder: build m 0.86 s (8.9 MB), tree 2.05 s for
# 990,683 primes = 2.07 us/prime. Direct m %% p is 3.14 ms/prime, so the tree is
# 1,517x that and ~30,000x the Lucas route: 23.82 core-h becomes ~2.9 s.
# Amdahl caps the member at ~1.6x (62 -> ~38 core-h) since the Band II scan is
# untouched; the case for it is i=10, where the Lucas r(p) bill is order 1e3
# core-hours.
#
# m stays in the PARENT. It is never pickled, never written to the npz or the
# json, and check_witness never sees it -- the proof rule that a certificate is
# re-derivable from (N, K, k, p) alone is unchanged.
USE_M_FOR_RP = True
# Every Nth prime of a round is ALSO reduced by Lucas and compared. m is one
# integer, so if it is right it is right for every p -- but a silent GMP fault
# would poison an entire run, and check_witness only samples.
RP_CROSSCHECK = 1000

_M = None
_M_KEY = None


def _m(N: int, K: int):
    """m = C(N,K), built once in the parent. gmpy2.bincoef, not math.comb --
    measured 282x apart, which is why the old "never build m" cost premise held
    for the naive binomial and does not for GMP.

    KEYED on (N, K). Caching on "is it None" was wrong: one process sweeps a
    single member so it never bit in production, but the test suite runs i=3
    before i=5 and the i=5 rounds silently reduced against m_3. The 1-in-1000
    Lucas cross-check caught it on its first run, which is the entire argument
    for keeping that check in the live sweep.
    """
    global _M, _M_KEY
    if _M_KEY != (N, K):
        _M = gmpy2.bincoef(N, K)
        _M_KEY = (N, K)
    return _M


# Primes per tree. The whole ladder in one tree is fastest but the product
# tree holds every level at once, and at i=10 scale that MEASURED a 2,286 MB
# peak working set -- in the parent, while eight workers run, on a machine with
# a memory-crash history. Blocking to 8 pieces measured 766 MB for 27.0 s
# against 17.7 s, i.e. 3x less memory for 1.5x time, with byte-identical
# output. Nine seconds against i=10's order-1e3 core-hour Lucas bill is free.
# 1e6 leaves i=9 (990,683 primes) in a single tree and splits i=10 into seven.
TREE_BLOCK = 1_000_000


def remainder_tree(m, primes: list[int], block: int = TREE_BLOCK) -> list[int]:
    """m mod p for every p. Product tree up, remainder tree down.

    Blocked at `block` primes to bound peak memory; the blocks are independent
    and concatenate, so the result is identical to one big tree.
    """
    if block and len(primes) > block:
        out: list[int] = []
        for s in range(0, len(primes), block):
            out.extend(remainder_tree(m, primes[s:s + block], block=0))
        return out
    level = [gmpy2.mpz(int(x)) for x in primes]
    levels = [level]
    while len(level) > 1:
        level = [level[i] * level[i + 1] if i + 1 < len(level) else level[i]
                 for i in range(0, len(level), 2)]
        levels.append(level)
    cur = [m % levels[-1][0]]
    for d in range(len(levels) - 2, -1, -1):
        below, out = levels[d], []
        for i, r in enumerate(cur):
            a = 2 * i
            out.append(r % below[a])
            if a + 1 < len(below):
                out.append(r % below[a + 1])
        cur = out
    return [int(x) for x in cur]


def round_rp(fam, primes: list[int]) -> dict:
    """{p: r(p)} for one round's primes. MUST be called per round -- rounds 2+
    scan a different prime set, derived from the previous round's survivors."""
    if not USE_M_FOR_RP or not primes:
        return {}
    from singmaster_intersect import binom_mod_lucas

    t0 = time.time()
    rmap = dict(zip(primes, remainder_tree(_m(fam.N, fam.K), primes)))
    checked = 0
    for j in range(0, len(primes), RP_CROSSCHECK or len(primes)):
        q = primes[j]
        check(rmap[q] == int(binom_mod_lucas(fam.N, fam.K, q)),
              f"r(p) from m disagrees with Lucas at p={q}")
        checked += 1
    print(f"  r(p) for {len(primes):,} primes in {time.time()-t0:.1f}s "
          f"({checked} Lucas cross-checks)", flush=True)
    return rmap
SMALL_K = 10**3
N_CHUNKS = 32
# MEASURED at i=8, three arms, byte-identical survivors in every one:
#     numpy@8   wall 772.6 s   Band II cpu 5071.5 s
#     GM@8      wall 168.8 s   Band II cpu  826.6 s     kernel alone   4.58x
#     GM@16     wall 163.6 s   Band II cpu 1347.5 s     +workers       1.03x
# The kernel is the whole win; 8 -> 16 workers buys 3% for 2.9x the memory
# (2,345 MB vs ~800 MB) and 1.63x MORE cpu for identical work. A 40-column
# microbenchmark had predicted 1.61x from the worker count; it did not survive
# real chunks -- the third component ratio in this project to fail that way.
#
# The cause of the cpu inflation is NOT established. "E-cores are slower here"
# was an inference, and pinning the kernel to individual cores does not support
# it (0.177-0.267 Gelem/s with no bimodal split, core 7 slowest and core 12
# fastest). The competing explanation is L3/memory contention -- Band II
# streams ~5 MB factorial windows per worker against a 36 MB L3 -- which more
# cores of any kind would not fix. Untested either way.
DEFAULT_WORKERS = 8
K_EXACT = {2: 200, 3: 200, 4: 200, 5: 200, 6: 200, 7: 200, 9: 80}
# i=2..7: exact k_extra=200 in fibonacci_i1-7.json
# i=9: modular k<=80 all impossible. i=1 is 3003 (N=8), skip.
# i=8 closed by the dedicated pipeline.
# What closed each member's pre-Z band, for the certificate sentence. The
# sentence used to hardcode "exact k<=200" for every i, which is false for i=9
# (modular k<=80) and i=8 (no pre-Z band at all: the Z-jump starts at k=3).
K_EXACT_KIND = {2: "exact intersect", 3: "exact intersect", 4: "exact intersect",
                5: "exact intersect", 6: "exact intersect", 7: "exact intersect",
                9: "modular"}


def certificate_basis(i: int, k_lo_z: int) -> str:
    """One clause stating what closed the columns below the Z-jump's start.

    The witness table carries an engine witness for every column k < k_lo_z
    regardless (witness._build_family), so the certificate rests on the table
    alone; this clause names the INDEPENDENT run that also closed that band,
    per member, rather than asserting an exact k<=200 result that only
    i=2..7 have.
    """
    if k_lo_z <= 2:
        return "No column lies below the Z-jump's start."
    kind = K_EXACT_KIND.get(i)
    if kind is None:
        return (f"Columns k<{k_lo_z} carry engine witnesses in the table; "
                f"no separate pre-Z run exists for i={i}.")
    return (f"Columns k<={K_EXACT[i]} were also closed independently by the "
            f"{kind} run, and carry engine witnesses in the table.")



def _job(payload: tuple) -> dict:
    kind, p, ks, N, K, r_expected = payload
    t0 = time.time()
    if kind == "bii":
        r, rows = scan_columns(p, ks, r_expected=r_expected, N=N, K=K)
    elif r_expected is not None:
        # The parent already reduced m against this round's primes. Scan only.
        # scan_columns_general is deliberately NOT given a "skip Lucas" branch:
        # it stays the no-r path, and this is a different route to the same r
        # rather than a shortcut through the same one. scan_ks_windowed still
        # refuses r=0, and a wrong r cannot mint a certificate -- check_witness
        # re-derives r from N and K by Lucas and rejects it.
        r = int(r_expected)
        rows = scan_ks_windowed(p, r, ks)
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


def run_jobs(jobs: list[tuple], workers: int, chk: Path, tag: str,
             done_keys: set, tags_on_disk: set | None = None) -> list[dict]:
    """Run this round's chunks, recovering anything already on disk.

    `done_keys` is updated IN PLACE as each chunk lands, so the caller never has
    to re-derive it from the whole checkpoint. That re-derivation, plus the
    survivor recovery below, used to re-parse the jsonl twice per round: at i=9
    about 42 full passes over a file that reaches 219 MB. Measured 16.6 MB/s
    and ~5.7x file size in RAM, so roughly 1.25 GB of Python objects allocated
    in the parent, 24 times, on a machine that has already been crashed once by
    memory pressure. The time was never the problem (~185 s against 62 core-h);
    the allocation was.
    """
    pending = []
    for job in jobs:
        _kind, p, ks, _N, _K, _r = job
        key = (tag, p, int(ks[0]), int(ks[-1]))
        if key not in done_keys:
            pending.append(job)
    surv = []
    # Only look for recoverable survivors when a chunk of THIS tag is already
    # on disk. On a fresh round there is nothing to recover and the parse is
    # pure cost -- it read the whole file to build an empty list.
    if tag in (tags_on_disk if tags_on_disk is not None
               else {k[0] for k in done_keys}):
        # Stream, filtering to THIS tag. Materialising the whole checkpoint to
        # recover one round's survivors is what cost ~1.25 GB per call at i=9;
        # peak here scales with the round's survivor list instead of the file.
        for rec in stream_done(chk):
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
            # Keep done_keys current in place. The record is on disk first, so
            # a crash between the two lines loses only the in-memory copy --
            # the next run's startup parse recovers it from the file.
            done_keys.add((tag, int(rec["p"]), int(rec["k_lo"]), int(rec["k_hi"])))
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
    # Stream the checkpoint rather than materialising it. Only two small
    # things are needed here -- the set of finished chunk keys, and which
    # phases are complete -- but the records carry whole survivor lists, so
    # holding them all costs ~5.7x the file size in Python objects (about
    # 1.25 GB at i=9's 219 MB) in the parent while the workers run.
    done_keys = set()
    tags_on_disk: set = set()
    complete = set()
    phase_counts: list[dict] = []
    for r in stream_done(chk):
        if "tag" in r:
            # Tracked SEPARATELY from done_keys. Inferring "is this tag on
            # disk" from done_keys would silently depend on every tagged
            # record also carrying p/k_lo -- true today, but a record type
            # with a tag and no p would then be invisible, the recovery guard
            # would read False, and that tag's survivors would vanish.
            tags_on_disk.add(r["tag"])
        if "tag" in r and "p" in r and "k_lo" in r:
            done_keys.add((r["tag"], r["p"], r["k_lo"], r["k_hi"]))
        elif r.get("event") == "phase_complete":
            complete.add(r["phase"])
            phase_counts.append(r)
    done = phase_counts

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
            surv = run_jobs(jobs, workers, chk, tag, done_keys, tags_on_disk)
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
            # done_keys is maintained in place by run_jobs. Re-deriving it
            # here meant a second full parse of the growing checkpoint every
            # round; VERIFY_DONE_KEYS re-derives and compares instead, so the
            # equivalence is testable without paying for it in production.
            if VERIFY_DONE_KEYS:
                _recs = load_done(chk)
                _full = {
                    (r["tag"], r["p"], r["k_lo"], r["k_hi"])
                    for r in _recs
                    if "tag" in r and "p" in r and "k_lo" in r
                }
                check(_full == done_keys,
                      f"incremental done_keys diverged from the checkpoint: "
                      f"missing {sorted(_full - done_keys)[:5]}, "
                      f"extra {sorted(done_keys - _full)[:5]}")
                # Keys are the easy half. The survivors are what feed `alive`,
                # the size-law ledger and the n_alive written to the record --
                # and the ledger CANNOT catch an under-recovery, because
                # escalate() is a one-sided upper-tail test, so observed far
                # BELOW expected reads as "ordinary".
                _on_disk = {int(c["k"]) for r in _recs if r.get("tag") == tag
                            for c in (r.get("survivors") or [])}
                _in_hand = {int(c["k"]) for c in surv}
                check(_on_disk == _in_hand,
                      f"recovered survivors for {tag} diverge from the "
                      f"checkpoint: missing {sorted(_on_disk - _in_hand)[:5]}, "
                      f"extra {sorted(_in_hand - _on_disk)[:5]}")
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
        # Survivors the extended rounds do not serve. They are SET ASIDE, never
        # dropped -- see the loop below.
        z_deferred: list = []
        cap_z = CAP_Z_SMALL_K if k_lo_z < SMALL_K else CAP_Z
        for rnd in range(1, cap_z + 1):
            if not current:
                break
            if rnd > CAP_Z:
                # Extra rounds serve the small-k tail only. Columns with
                # k >= SMALL_K are NOT killed here, so they must still be
                # REPORTED as survivors rather than filtered out of existence.
                #
                # Dropping them is what let i=8 certify k=1021, a column it
                # never killed: k=145 and k=1021 both survived round 12, round
                # 13 filtered k=1021 away, n_z_alive fell to 0, clean went True,
                # and the witness builder recorded the last prime k=1021 had
                # merely SURVIVED as the prime that killed it.
                keep = [c for c in current if int(c["k"]) < SMALL_K]
                z_deferred.extend(c for c in current if int(c["k"]) >= SMALL_K)
                current = keep
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
            rmap = round_rp(fam, sorted(buckets))
            for p, ks in sorted(buckets.items()):
                for ch in chunk_ks(ks, p, N_CHUNKS if len(ks) >= 2000 else 1):
                    jobs.append(("z", p, ch, fam.N, fam.K, rmap.get(p)))
            exp = ledger_z.expect(
                (k_, p_) for p_, ks_ in buckets.items() for k_ in ks_
            )
            n_entered = sum(len(v) for v in buckets.values())
            surv = run_jobs(jobs, workers, chk, f"z{rnd}", done_keys, tags_on_disk)
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
            # done_keys is maintained in place by run_jobs. Re-deriving it
            # here meant a second full parse of the growing checkpoint every
            # round; VERIFY_DONE_KEYS re-derives and compares instead, so the
            # equivalence is testable without paying for it in production.
            if VERIFY_DONE_KEYS:
                _recs = load_done(chk)
                _full = {
                    (r["tag"], r["p"], r["k_lo"], r["k_hi"])
                    for r in _recs
                    if "tag" in r and "p" in r and "k_lo" in r
                }
                check(_full == done_keys,
                      f"incremental done_keys diverged from the checkpoint: "
                      f"missing {sorted(_full - done_keys)[:5]}, "
                      f"extra {sorted(done_keys - _full)[:5]}")
                # Keys are the easy half. The survivors are what feed `alive`,
                # the size-law ledger and the n_alive written to the record --
                # and the ledger CANNOT catch an under-recovery, because
                # escalate() is a one-sided upper-tail test, so observed far
                # BELOW expected reads as "ordinary".
                _tag = f"z{rnd}"
                _on_disk = {int(c["k"]) for r in _recs if r.get("tag") == _tag
                            for c in (r.get("survivors") or [])}
                _in_hand = {int(c["k"]) for c in surv}
                check(_on_disk == _in_hand,
                      f"recovered survivors for {_tag} diverge from the "
                      f"checkpoint: missing {sorted(_on_disk - _in_hand)[:5]}, "
                      f"extra {sorted(_in_hand - _on_disk)[:5]}")
        write_jsonl(
            chk,
            {
                "event": "phase_complete",
                "phase": "zjump",
                "n_alive": len(current) + len(z_deferred),
                "n_nolive": len(z_none),
                "n_deferred": len(z_deferred),
            },
        )
        phases["zjump"] = zrounds
        # Deferred columns are survivors. They block `clean` exactly as an
        # ordinary survivor does; the extended rounds simply never tested them.
        z_left = current + z_deferred
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
    # SEPARATE from `clean`. `clean` is the SWEEP's verdict -- every column was
    # testable and nothing survived. Whether the proof object could be written
    # is a different question, and conflating them would let a clean sweep be
    # reported as dirty (or, before this, a failed build be reported as
    # certified). A run may legitimately be clean=True with certificate=None.
    witness_ok = True
    withheld = None
    # Bound BEFORE the try: the except handler reads it. With it assigned
    # after `import witness` inside the try, an ImportError became an
    # UnboundLocalError in the handler, main() crashed, and a finished sweep
    # had no json at all -- the exact loss the handler exists to prevent.
    wpath = ROOT / "results" / f"i{i}_witness.npz"
    if clean:
        try:
            import witness as _witness

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
            # ...but never certify one either. Before this, a failed build left
            # clean=True and still emitted a certificate reading
            # "results/i{i}_witness.npz (sha256 n/a)". The dangerous case is a
            # STALE npz from an earlier run: a referee who opens the named path
            # verifies a table this sweep did not write, and it passes.
            witness_ok = False
            witness_meta = {"error": repr(exc)}
            stale = wpath.exists()
            withheld = (
                f"witness build failed: {exc!r}"
                + ("; results/{}_witness.npz EXISTS but was NOT written by this "
                   "run -- do not verify it as this run's proof object"
                   .format(f"i{i}") if stale else "")
            )
            print(f"  witness build failed: {exc!r}", flush=True)
            print("  certificate WITHHELD: the sweep may be clean, but its "
                  "witness table was not produced", flush=True)
            if stale:
                print(f"  WARNING {wpath.name} exists and is STALE -- it is not "
                      f"this run's output", flush=True)
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
        "witness_ok": witness_ok,
        "certificate_withheld": withheld,
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
            f"has r(p) notin I_{{p,k}}. {certificate_basis(i, k_lo_z)} "
            f"Together with the two family columns, N(C(N,K))=6. Not Singmaster. "
            f"Witness prime for every column: results/i{i}_witness.npz "
            f"(sha256 {witness_meta['sha256'][:16] if witness_meta and 'sha256' in witness_meta else 'n/a'}); "
            f"re-check with: python scripts/witness.py verify "
            f"--file results/i{i}_witness.npz"
        )
        # A certificate names a witness table. Emit it only when that table was
        # actually written by THIS run -- `clean` alone is not enough, because a
        # failed build leaves the named path absent or stale.
        if clean and witness_ok
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
