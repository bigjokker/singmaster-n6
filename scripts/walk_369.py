#!/usr/bin/env python3
"""Walk the 369 fat-image triples past prime 3.

Pre-registered test of Claude's size law (2026-08-19). Not a fishing trip.
Image clause only. Stops at kill, cell end (nhi), digit-0, or MAX_RUN.
Digit-0 / leaving (alpha,beta) is NOT image-run length.

Predictions (from 369 triples, independence-after-size):
  survive prime 4 ~44  (2sigma ~32-58; ~40 from (3,1))
  survive prime 5 ~6
  survive prime 6 ~0.8
  survive prime 7 ~0.1
  max_run 5 or 6
  7 mild; 8 real correlation (Q1 back); 100+ quadruples = law incomplete

Requires results/fat_image_hunt.json. Checkpoint jsonl after each k.
Resume if jsonl exists. Refuse if the final json exists.
Does not touch nextprime_sweep. No giant m. No 1e6-K until-kill.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gmpy2

from bandii_kernel import (  # noqa: E402
    image_j,
    inv_table,
    lucas_digits as lucas,
)

N = 10_803_704
K = 4_126_647
SRC = ROOT / "results" / "fat_image_hunt.json"
OUT = ROOT / "results" / "walk_369.json"
CHK = ROOT / "results" / "walk_369.jsonl"
MAX_RUN = 12

PREREGISTER = {
    "survive_4": 44,
    "survive_4_2sigma": [32, 58],
    "survive_4_from_3_1": 40,
    "survive_5": 6,
    "survive_6": 0.8,
    "survive_7": 0.1,
    "max_run": "5 or 6",
    "run_7": "mild",
    "run_8": "real correlation; Q1 back",
    "quadruples_100plus": "size law incomplete",
    "even_g_n3": 0.878,
    "measured_even_g_n3": 0.900,
}


def row_id(key: str, k: int) -> str:
    return f"{key}|{k}"


def load_done() -> dict[str, dict]:
    done = {}
    if not CHK.exists():
        return done
    with CHK.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            done[row_id(rec["key"], rec["k"])] = rec
    return done


def walk_one(trip: dict, region: dict) -> dict:
    k = trip["k"]
    hits = list(trip["hits"])
    nlo, nhi = region["nlo"], region["nhi"]
    a0, b0 = region["a"], region["b"]
    last_p = hits[-1]["p"]
    extra = []
    n_run = len(hits)
    p = int(gmpy2.next_prime(last_p))
    stop = None
    kill_p = None
    left_p = None
    while p <= nhi and n_run < MAX_RUN:
        a, b, prod = lucas(p)
        if a != a0 or b != b0:
            stop = "left_cell"
            left_p = p
            break
        if prod == 0:
            stop = "digit0"
            left_p = p
            break
        g = p - k
        inv = inv_table(p, g)
        j = image_j(prod, k, p, inv)
        if j is None:
            stop = "killed"
            kill_p = p
            break
        extra.append({"p": p, "g": g, "j": j, "prod": prod})
        n_run += 1
        p = int(gmpy2.next_prime(p))
    if stop is None:
        if n_run >= MAX_RUN:
            stop = "max_run"
        elif p > nhi:
            stop = "cell_end"
        else:
            stop = "killed"
    g0 = hits[0]["g"]
    return {
        "key": trip["key"],
        "k": k,
        "g0": g0,
        "even_g": g0 % 2 == 0,
        "n_run": n_run,
        "stop": stop,
        "kill_p": kill_p,
        "left_p": left_p,
        "nlo": nlo,
        "nhi": nhi,
        "extra": extra,
    }


def summarize(rows: list[dict]) -> dict:
    stops = Counter(r["stop"] for r in rows)
    max_run = max((r["n_run"] for r in rows), default=0)
    survive = {}
    for n in range(4, MAX_RUN + 1):
        survive[str(n)] = sum(1 for r in rows if r["n_run"] >= n)
    by_cell = {}
    for r in rows:
        cell = by_cell.setdefault(
            r["key"],
            {"n": 0, "max_run": 0, "survive_4": 0, "survive_5": 0, "survive_6": 0, "survive_7": 0},
        )
        cell["n"] += 1
        cell["max_run"] = max(cell["max_run"], r["n_run"])
        if r["n_run"] >= 4:
            cell["survive_4"] += 1
        if r["n_run"] >= 5:
            cell["survive_5"] += 1
        if r["n_run"] >= 6:
            cell["survive_6"] += 1
        if r["n_run"] >= 7:
            cell["survive_7"] += 1
    s4 = [r for r in rows if r["n_run"] >= 4]
    even4 = sum(1 for r in s4 if r["even_g"])
    return {
        "n_rows": len(rows),
        "max_run": max_run,
        "survive": survive,
        "n_quadruples": survive.get("4", 0),
        "stops": dict(stops),
        "even_g_survive_4": None if not s4 else round(even4 / len(s4), 4),
        "n_even_survive_4": even4,
        "by_cell": by_cell,
    }


def main() -> int:
    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    if not SRC.exists():
        print(f"missing {SRC}", flush=True)
        return 1

    src = json.loads(SRC.read_text(encoding="utf-8"))
    regions = {r["key"]: r for r in src["regions"]}
    triples = src["triples"]
    if len(triples) != 369:
        print(f"warning: expected 369 triples, got {len(triples)}", flush=True)

    done = load_done()
    t0 = time.time()
    print(
        f"=== walk 369  resume={len(done)}/{len(triples)}  "
        f"MAX_RUN={MAX_RUN} ===",
        flush=True,
    )
    print(
        "preregister  s4~44  s5~6  s6~0.8  s7~0.1  max_run 5-6",
        flush=True,
    )

    chk = CHK.open("a", encoding="utf-8")
    rows = []
    try:
        for i, trip in enumerate(triples, 1):
            rid = row_id(trip["key"], trip["k"])
            if rid in done:
                rec = done[rid]
                rows.append(rec)
                continue
            region = regions[trip["key"]]
            rec = walk_one(trip, region)
            chk.write(json.dumps(rec) + "\n")
            chk.flush()
            rows.append(rec)
            done[rid] = rec
            if i % 10 == 0 or rec["n_run"] >= 4 or rec["stop"] != "killed":
                alive4 = sum(1 for r in rows if r["n_run"] >= 4)
                alive5 = sum(1 for r in rows if r["n_run"] >= 5)
                print(
                    f"  {i}/{len(triples)} {trip['key']} k={trip['k']} "
                    f"n_run={rec['n_run']} {rec['stop']}  "
                    f"s4={alive4} s5={alive5}",
                    flush=True,
                )
    finally:
        chk.close()

    meas = summarize(rows)
    payload = {
        "search": "walk_369",
        "N": N,
        "K": K,
        "source": str(SRC.name),
        "n_triples": len(triples),
        "max_run_cap": MAX_RUN,
        "preregister": PREREGISTER,
        "measured": meas,
        "seconds": round(time.time() - t0, 3),
        "rows": rows,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(flush=True)
    print(
        f"wrote {OUT}  max_run={meas['max_run']}  "
        f"s4={meas['survive'].get('4', 0)}  "
        f"s5={meas['survive'].get('5', 0)}  "
        f"s6={meas['survive'].get('6', 0)}  "
        f"s7={meas['survive'].get('7', 0)}  "
        f"{payload['seconds']}s",
        flush=True,
    )
    print(f"stops {meas['stops']}", flush=True)
    print("preregister s4~44 s5~6 s6~0.8 s7~0.1 max_run 5-6", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
