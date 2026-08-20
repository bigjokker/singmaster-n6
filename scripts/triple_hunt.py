#!/usr/bin/env python3
"""NONE-window image-run hunt: first three primes, look for a triple.

Does NOT rerun the 1e5-1e6 census (already: max run 2, 0 triples).
New territory: two-digit NONE cells with P_hi > 1e6.

Exhaustive under-Z k only when g_max <= G_EXHAUST (else three
canonical k's only). Fat slabs like (4,1) Zw~475k are O(g^2) if
done fully — that is not this job.

Does not touch nextprime_sweep.py. No giant m.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


from bandii_kernel import (  # noqa: E402
    cells,
    first_primes,
    image_j,
    inv_table,
    preceding_z,
    r_of,
)

N = 10_803_704
K = 4_126_647
D = N - K
P_NEW = 1_000_000
G_EXHAUST = 25_000
OUT = ROOT / "results" / "triple_hunt_p1e6-K.json"


def walk3(k: int, primes: list[int], prods: list[int]) -> dict:
    hits = []
    invs: dict[int, list[int]] = {}
    for p, prod in zip(primes, prods):
        if p <= k:
            continue
        g = p - k
        if p not in invs:
            invs[p] = inv_table(p, g)
        elif len(invs[p]) <= g:
            invs[p] = inv_table(p, g)
        j = image_j(prod, k, p, invs[p])
        if j is None:
            return {"k": k, "n": len(hits), "hits": hits, "kill_p": p}
        hits.append({"p": p, "g": g, "j": j, "prod": prod})
    return {"k": k, "n": len(hits), "hits": hits, "kill_p": None}


def main() -> int:
    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    t0 = time.time()
    windows = cells()
    none_new = [
        (i, w)
        for i, w in enumerate(windows)
        if w["kind"] == "NONE" and w["phi"] > P_NEW
    ]
    print(
        f"=== triple hunt  NONE cells with P_hi>{P_NEW}  "
        f"n={len(none_new)}  exhaust if g<={G_EXHAUST} ===",
        flush=True,
    )
    reports = []
    triples = []
    max_run = 0
    best = None
    for i, w in none_new:
        primes = first_primes(w["plo"], w["phi"], 3)
        if len(primes) < 1:
            continue
        p0 = primes[0]
        z = preceding_z(windows, i)
        k_lo = (z["z_first"] - 1) if z else max(2, w["plo"] - 1)
        g_max = p0 - k_lo
        ks = {p0 - 1}
        if z:
            ks.add(z["z_last"])
            ks.add(z["z_first"] - 1)
        exhaustive = g_max <= G_EXHAUST
        if exhaustive:
            ks.update(range(k_lo, p0))
        ks = {k for k in ks if 2 <= k < p0}
        prods = [r_of(p) for p in primes]
        local_max = 0
        local_best = None
        n_ge2 = n_ge3 = 0
        for k in sorted(ks):
            rec = walk3(k, primes, prods)
            local_max = max(local_max, rec["n"])
            if rec["n"] >= 2:
                n_ge2 += 1
            if rec["n"] >= 3:
                n_ge3 += 1
                triples.append(
                    {
                        "window": (w["a"], w["b"]),
                        "k": k,
                        "hits": rec["hits"],
                        "kill_p": rec["kill_p"],
                    }
                )
            if rec["n"] > local_max - 1 and rec["n"] >= local_max:
                local_best = rec
        if local_max > max_run:
            max_run = local_max
            best = {"window": (w["a"], w["b"]), **(local_best or {})}
        print(
            f"  ({w['a']},{w['b']}) plo={w['plo']} phi={w['phi']} "
            f"p0={p0} g_max={g_max} ks={len(ks)} exh={exhaustive} "
            f"max_run={local_max} ge2={n_ge2} ge3={n_ge3}",
            flush=True,
        )
        reports.append(
            {
                "a": w["a"],
                "b": w["b"],
                "plo": w["plo"],
                "phi": w["phi"],
                "primes": primes,
                "g_max": g_max,
                "n_k": len(ks),
                "exhaustive": exhaustive,
                "max_run": local_max,
                "n_ge2": n_ge2,
                "n_ge3": n_ge3,
                "best": local_best,
            }
        )
    payload = {
        "search": "triple_hunt_p1e6-K",
        "N": N,
        "K": K,
        "p_new": P_NEW,
        "g_exhaust": G_EXHAUST,
        "n_windows": len(reports),
        "max_run": max_run,
        "n_triples": len(triples),
        "triples": triples,
        "best": best,
        "seconds": round(time.time() - t0, 3),
        "windows": reports,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        f"wrote {OUT}  windows={len(reports)}  max_run={max_run}  "
        f"triples={len(triples)}  {payload['seconds']}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
