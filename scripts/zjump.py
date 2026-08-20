#!/usr/bin/env python3
"""Band I Z-jump remnant for i=8.

Worklist: Stage-3 hang-guards (89195) plus k=1000001..4126621
(stragglers 4126622..4126646 already certified). Skip family K, K+1.

For each k, jump to the first LIVE prime p>k (NONE / PART-lower /
p>N/2). Digit-0 and (K, N/2] are skipped analytically. Test at most
12 live primes with the factorial image kernel. Cap 12: anything
still up is an anomaly, not a representation.

Not a next-prime walk from k. Not until-kill. Not Band II.
Refuse if results/zjump.json exists. Resume from jsonl.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


from bandii_kernel import (  # noqa: E402
    D,
    K,
    N,
    P1,
    R_P1,
    cells,
    chunk_ks,
    fact_table,
    first_live_after,
    live_intervals,
    live_primes,
    r_closed,
    r_two_digit,
    scan_columns_general,
    scan_ks,
)


def load_done() -> list[dict]:
    from bandii_kernel import read_jsonl

    return read_jsonl(CHK)


def write_jsonl(rec: dict) -> None:
    from bandii_kernel import append_jsonl

    append_jsonl(CHK, rec)


def summarize(rows: list[dict]) -> dict:
    from bandii_kernel import summarize_survivors

    return summarize_survivors(rows, with_g=True)

OUT = ROOT / "results" / "zjump.json"
CHK = ROOT / "results" / "zjump.jsonl"
CAP = 12
N_CHUNKS = 32
DEFAULT_WORKERS = 8
MAX_INLINE = 20_000

# Stage-3 hang-guard runs (inclusive). Count must be 89195.
HANG_RUNS = [
    (515813, 520186),
    (607003, 615518),
    (615520, 615520),
    (687773, 700252),
    (741883, 751696),
    (834629, 880328),
    (953861, 962170),
]
HANG_N = 89195
K_LO_NEW = 1_000_001
K_HI_NEW = 4_126_621  # stragglers start at 4126622


def hang_ks() -> list[int]:
    out = []
    for a, b in HANG_RUNS:
        out.extend(range(a, b + 1))
    return out


def worklist() -> list[int]:
    ks = hang_ks()
    ks.extend(range(K_LO_NEW, K_HI_NEW + 1))
    return ks


def _job(payload: tuple) -> dict:
    p, ks = payload
    t0 = time.time()
    r, rows = scan_columns_general(p, ks)
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


def preflight() -> dict:
    print("=== Z-jump pre-flight ===", flush=True)
    t0 = time.time()
    from singmaster_intersect import binom_mod_lucas

    hg = hang_ks()
    if len(hg) != HANG_N:
        raise RuntimeError(f"hang-guard count {len(hg)} != {HANG_N}")
    n_new = K_HI_NEW - K_LO_NEW + 1
    print(f"  hang-guards {len(hg)}  1e6-K {n_new}  total {len(hg)+n_new}", flush=True)

    windows = cells()
    ivs = live_intervals(windows)
    w20 = next(w for w in windows if w["a"] == 20 and w["b"] == 7)
    if w20["kind"] != "FULL" or w20["z_last"] != 540185:
        raise RuntimeError(f"(20,7) {w20}")
    p_hang = first_live_after(515813, ivs)
    if p_hang is None or p_hang <= 540185:
        raise RuntimeError(f"hang k=515813 first live {p_hang}")
    if p_hang - 515813 <= 20000:
        raise RuntimeError(f"hang jump too short {p_hang - 515813}")
    print(f"  (20,7) FULL z_last=540185  k=515813 -> {p_hang}  g={p_hang-515813}", flush=True)

    p513 = first_live_after(513593, ivs)
    if p513 != 514499:
        raise RuntimeError(f"k=513593 first live {p513} != 514499")
    F = fact_table(p513)
    r = r_two_digit(F, p513)
    hits = scan_ks(F, p513, r, [513593])
    if hits:
        raise RuntimeError(f"k=513593 survived {p513}: {hits}")
    print("  k=513593 jumps to 514499 and dies (image miss)", flush=True)

    ps = live_primes(268733, 3, ivs)
    if ps[:3] != [270097, 270121, 270131]:
        raise RuntimeError(f"k=268733 live {ps[:3]}")
    F1 = fact_table(270097)
    r1 = r_two_digit(F1, 270097)
    h1 = scan_ks(F1, 270097, r1, [268733])
    F2 = fact_table(270121)
    r2 = r_two_digit(F2, 270121)
    h2 = scan_ks(F2, 270121, r2, [268733])
    F3 = fact_table(270131)
    r3 = r_two_digit(F3, 270131)
    h3 = scan_ks(F3, 270131, r3, [268733])
    if not h1 or not h2 or h3:
        raise RuntimeError(f"k=268733 hits {h1} {h2} {h3}")
    if h1[0]["b"] != 589 or h2[0]["b"] != 196:
        raise RuntimeError(f"k=268733 b {h1[0]['b']} {h2[0]['b']}")
    print("  k=268733 270097/589, 270121/196, kill 270131", flush=True)

    Fp = fact_table(P1)
    if r_two_digit(Fp, P1) != R_P1 or r_closed(P1) != R_P1:
        raise RuntimeError("r_two_digit(p1) != 1275205")
    if int(binom_mod_lucas(N, K, P1)) != R_P1:
        raise RuntimeError("lucas(p1) mismatch")
    print("  r_two_digit(p1)=1275205 matches closed form and lucas", flush=True)

    for p in (270097, 514499, P1):
        F = fact_table(p)
        if r_two_digit(F, p) != int(binom_mod_lucas(N, K, p)):
            raise RuntimeError(f"r vs lucas at {p}")
        if r_two_digit(F, p) == 0:
            raise RuntimeError(f"r=0 at supposed live {p}")
    print(f"  r_two_digit == lucas on sample live primes  {time.time()-t0:.1f}s", flush=True)
    print("=== pre-flight passed ===", flush=True)
    return {"windows": windows, "ivs": ivs}


def run_round(
    rnd: int,
    buckets: dict[int, list[int]],
    workers: int,
    done_keys: set[tuple],
) -> list[dict]:
    jobs = []
    for p, ks in sorted(buckets.items()):
        chunks = chunk_ks(ks, p, N_CHUNKS if len(ks) >= 2000 else 1)
        for ch in chunks:
            key = (rnd, p, int(ch[0]), int(ch[-1]))
            if key in done_keys:
                continue
            jobs.append((p, ch))

    survivors: list[dict] = []
    for rec in load_done():
        if rec.get("round") == rnd:
            survivors.extend(rec.get("survivors") or [])

    if not jobs:
        print(f"  round {rnd} all chunks already done  alive={len(survivors)}", flush=True)
        return survivors

    print(
        f"  round {rnd}  primes={len(buckets)}  jobs={len(jobs)}  "
        f"workers={workers}  cols={sum(len(v) for v in buckets.values())}",
        flush=True,
    )
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        for rec in pool.imap_unordered(_job, jobs):
            rec["round"] = rnd
            rows = rec["survivors"]
            if rnd == 1 and len(rows) > MAX_INLINE:
                rec["survivors"] = [{"k": x["k"], "g": x["g"], "g_even": x["g_even"], "b": x["b"]} for x in rows]
            write_jsonl(rec)
            survivors.extend(rows)
            print(
                f"    p={rec['p']} [{rec['k_lo']},{rec['k_hi']}] "
                f"cols={rec['n_cols']} surv={rec['n_survivors']} "
                f"{rec['seconds']}s  alive_so_far={len(survivors)}",
                flush=True,
            )
    return survivors


def assign_buckets(items: list, ivs: list[tuple[int, int]], rnd: int) -> tuple[dict[int, list[int]], list[int]]:
    """items is list of k (round 1) or survivor dicts (later)."""
    buckets: dict[int, list[int]] = defaultdict(list)
    no_prime: list[int] = []
    for it in items:
        if isinstance(it, dict):
            k = int(it["k"])
            last_p = int(it.get("last_p") or (it["k"] + it["g"] if "g" in it else it["k"]))
            # after a hit at p=k+g, next live after that p
            x = last_p
        else:
            k = int(it)
            x = k
        p = first_live_after(x, ivs)
        if p is None:
            no_prime.append(k)
            continue
        buckets[p].append(k)
    return buckets, no_prime


def main() -> int:
    if "--preflight" in sys.argv:
        try:
            import numpy as np  # noqa: F401
        except ImportError:
            print("numpy is required. Not installed.", flush=True)
            return 1
        preflight()
        return 0

    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        print("numpy is required. Install numpy >= 1.24 before this bat.", flush=True)
        return 1

    workers = int(os.environ.get("ZJUMP_WORKERS", DEFAULT_WORKERS))
    workers = max(1, min(workers, 16))
    OUT.parent.mkdir(exist_ok=True)

    t0 = time.time()
    geo = preflight()
    ivs = geo["ivs"]
    ks = worklist()
    n_hang = HANG_N
    n_new = K_HI_NEW - K_LO_NEW + 1
    print(
        f"=== Z-jump  hang={n_hang}  1e6-K={n_new}  total={len(ks)}  "
        f"cap={CAP}  workers={workers} ===",
        flush=True,
    )

    done_recs = load_done()
    done_keys = {
        (r["round"], r["p"], r["k_lo"], r["k_hi"])
        for r in done_recs
        if "round" in r and "p" in r and "k_lo" in r
    }
    complete_rounds = [r["round"] for r in done_recs if r.get("event") == "round_complete"]
    last_complete = max(complete_rounds) if complete_rounds else 0
    if last_complete:
        alive_rows = []
        for r in done_recs:
            if r.get("round") == last_complete and "survivors" in r:
                for s in r.get("survivors") or []:
                    s = dict(s)
                    if "g" in s:
                        s["last_p"] = s["k"] + s["g"]
                    alive_rows.append(s)
        start_round = last_complete + 1
        print(
            f"  resume after round {last_complete}  alive={len(alive_rows)}",
            flush=True,
        )
    else:
        alive_rows = [{"k": k} for k in ks]
        start_round = 1
        if done_keys:
            print(f"  resume incomplete round 1  done_chunks={len(done_keys)}", flush=True)

    rounds = []
    anomalies = []
    current = alive_rows
    for rnd in range(start_round, CAP + 1):
        if not current:
            print(f"  all dead before round {rnd}. stopping.", flush=True)
            break
        if rnd == 1 and current and "g" not in current[0]:
            buckets, none = assign_buckets([c["k"] for c in current], ivs, rnd)
        else:
            buckets, none = assign_buckets(current, ivs, rnd)
        if none:
            print(f"  round {rnd} no further live prime: {len(none)} k (anomaly)", flush=True)
            anomalies.extend(none)
        surv = run_round(rnd, buckets, workers, done_keys)
        for s in surv:
            s["last_p"] = s["k"] + s["g"]
        sm = summarize(surv)
        rec = {"round": rnd, "n_primes": len(buckets), **sm}
        rounds.append(rec)
        print(
            f"  ROUND {rnd} alive={sm['n']} even={sm['even']} "
            f"mean_k={sm['mean_k']} mean_g={sm['mean_g']}",
            flush=True,
        )
        current = surv
        write_jsonl({"event": "round_complete", "round": rnd, "n_alive": sm["n"]})
        done_keys = {
            (r["round"], r["p"], r["k_lo"], r["k_hi"])
            for r in load_done()
            if "round" in r and "p" in r and "k_lo" in r
        }

    if current:
        anomalies.extend(int(s["k"]) if isinstance(s, dict) else int(s) for s in current)

    clean = len(anomalies) == 0
    payload = {
        "search": "zjump",
        "N": N,
        "K": K,
        "cap": CAP,
        "n_hang": n_hang,
        "n_1e6_K": n_new,
        "n_columns": len(ks),
        "k_ranges": {
            "hang_runs": HANG_RUNS,
            "new": [K_LO_NEW, K_HI_NEW],
        },
        "workers": workers,
        "rounds": rounds,
        "n_anomalies": len(anomalies),
        "anomalies": anomalies[:1000],
        "clean": clean,
        "certificate": (
            "Every leftover Band I k (Stage-3 hang-guards and "
            "1000001..4126621) carries r(p) notin I_{p,k} for some live "
            "prime after a Z-jump. Together with Stage 1-3 kills, "
            "stragglers, modular small-k, and Band II, "
            "N(C(F_18 F_19, F_16 F_19))=6 exactly. Not Singmaster."
        )
        if clean
        else None,
        "seconds": round(time.time() - t0, 3),
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(flush=True)
    print(
        f"wrote {OUT}  rounds={len(rounds)}  anomalies={len(anomalies)}  "
        f"clean={clean}  {payload['seconds']}s",
        flush=True,
    )
    if clean:
        print(
            "Band I remnant closed. Hang-guards + 1e6-K all killed.",
            flush=True,
        )
    else:
        print("ANOMALIES at cap or no live prime. Log and hand back. Do not extend.", flush=True)
        print("  first 20:", anomalies[:20], flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
