#!/usr/bin/env python3
"""Overnight: exhaustive first-three-prime image runs on fat non-zero cells.

This is the remaining Band I census that the bounded triple hunt skipped.
It covers NONE cells AND PART-lower (digit-0 silent) with g_max > 25000
and P_hi > 1e6, largest g first (worst rho is (3,1) lower, ~0.176).

A triple here is usable: it kills two-to-three as a hoped lemma.
A clean max-run=2 on these cells is the first exhaustive look at the
high-rho end of Band I.

O(g^2) per first prime. (3,1) g~4.75e5 is the long piece.
Checkpoint jsonl after each region. Resume if that file exists.
Does not touch nextprime_sweep.py. No giant m. No 1e6-K until-kill.
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
P_CUT = 1_000_000
G_MIN = 25_000
OUT = ROOT / "results" / "fat_image_hunt.json"
CHK = ROOT / "results" / "fat_image_hunt.jsonl"


def surviving_g(prod: int, p: int, gmin: int, gmax: int) -> list[tuple[int, int]]:
    """All (g, j) in [gmin, gmax] with prod in I_g."""
    inv = inv_table(p, max(gmax, 1))
    hits = []
    tmark = time.time()
    for g in range(gmin, gmax + 1):
        if g % 25000 == 0:
            print(
                f"      g={g}/{gmax}  hits={len(hits)}  "
                f"{time.time()-tmark:.0f}s",
                flush=True,
            )
            tmark = time.time()
        n = g - 1
        c = 1
        sign = 1
        found = None
        for j in range(g):
            if (sign * c) % p == prod:
                found = j
                break
            if j + 1 >= g:
                break
            c = c * (n - j) % p
            c = c * inv[j + 1] % p
            sign = -sign
        if found is not None:
            hits.append((g, found))
    return hits


def worklist(windows: list[dict]) -> list[dict]:
    jobs = []
    for i, w in enumerate(windows):
        if w["kind"] == "NONE":
            nlo, nhi = w["plo"], w["phi"]
            tag = "NONE"
        elif w["kind"] == "PART":
            nlo, nhi = w["plo"], w["zlo"] - 1
            tag = "PART-lower"
            if nlo > nhi:
                continue
        else:
            continue
        if nhi <= P_CUT:
            continue
        nlo = max(nlo, 2)
        primes = first_primes(nlo, nhi, 3)
        if not primes:
            continue
        z = preceding_z(windows, i)
        k_lo = (z["z_first"] - 1) if z else max(2, nlo - 1)
        p0 = primes[0]
        g_max = p0 - k_lo
        if g_max <= G_MIN:
            continue
        rho = g_max / p0
        jobs.append(
            {
                "key": f"{w['a']},{w['b']}:{tag}",
                "a": w["a"],
                "b": w["b"],
                "tag": tag,
                "nlo": nlo,
                "nhi": nhi,
                "primes": primes,
                "k_lo": k_lo,
                "g_max": g_max,
                "rho": round(rho, 5),
                "prev_z": None
                if z is None
                else {
                    "a": z["a"],
                    "b": z["b"],
                    "z_first": z["z_first"],
                    "z_last": z["z_last"],
                },
            }
        )
    jobs.sort(key=lambda j: -j["g_max"])
    return jobs


def walk_rest(k: int, primes: list[int], prods: list[int], j0: int) -> dict:
    hits = [{"p": primes[0], "g": primes[0] - k, "j": j0, "prod": prods[0]}]
    for p, prod in zip(primes[1:], prods[1:]):
        g = p - k
        inv = inv_table(p, g)
        j = image_j(prod, k, p, inv)
        if j is None:
            return {"k": k, "n": len(hits), "hits": hits, "kill_p": p}
        hits.append({"p": p, "g": g, "j": j, "prod": prod})
    return {"k": k, "n": len(hits), "hits": hits, "kill_p": None}


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
            done[rec["key"]] = rec
    return done


def main() -> int:
    if OUT.exists():
        print(f"{OUT} already exists. Not rerunning.", flush=True)
        return 2
    t0 = time.time()
    windows = cells()
    jobs = worklist(windows)
    done = load_done()
    print(f"=== fat image hunt  {len(jobs)} regions  resume={len(done)} ===", flush=True)
    for j in jobs:
        mark = "done" if j["key"] in done else "todo"
        print(
            f"  {mark} {j['key']} g_max={j['g_max']} rho={j['rho']} "
            f"p0={j['primes'][0]} k_lo={j['k_lo']}",
            flush=True,
        )

    reports = []
    triples = []
    max_run = 0
    best = None
    chk = CHK.open("a", encoding="utf-8")
    try:
        for job in jobs:
            if job["key"] in done:
                rec = done[job["key"]]
                reports.append(rec)
                max_run = max(max_run, rec.get("max_run") or 0)
                triples.extend(rec.get("triples") or [])
                if rec.get("best") and rec["max_run"] >= max_run:
                    best = rec["best"]
                continue
            print(f"\n-- {job['key']} g_max={job['g_max']} --", flush=True)
            primes = job["primes"]
            t1 = time.time()
            prods = [r_of(p) for p in primes]
            print(f"   lucas {[ (p, pr) for p, pr in zip(primes, prods) ]}", flush=True)
            p0 = primes[0]
            prod0 = prods[0]
            gmin = 1
            gmax = job["g_max"]
            sg = surviving_g(prod0, p0, gmin, gmax)
            print(f"   first-survivors={len(sg)}  {time.time()-t1:.1f}s", flush=True)
            local_max = 0
            local_best = None
            n_ge2 = n_ge3 = 0
            loc_triples = []
            examples = []
            for g, j0 in sg:
                k = p0 - g
                rec = walk_rest(k, primes, prods, j0)
                nrun = rec["n"]
                if nrun >= 2:
                    n_ge2 += 1
                if nrun >= 3:
                    n_ge3 += 1
                    loc_triples.append(
                        {
                            "key": job["key"],
                            "k": k,
                            "hits": rec["hits"],
                            "kill_p": rec["kill_p"],
                        }
                    )
                if nrun > local_max:
                    local_max = nrun
                    local_best = rec
                if nrun >= 2 and len(examples) < 6:
                    examples.append(rec)
            row = {
                "key": job["key"],
                "a": job["a"],
                "b": job["b"],
                "tag": job["tag"],
                "nlo": job["nlo"],
                "nhi": job["nhi"],
                "primes": primes,
                "prods": prods,
                "k_lo": job["k_lo"],
                "g_max": job["g_max"],
                "rho": job["rho"],
                "n_first_survivors": len(sg),
                "max_run": local_max,
                "n_ge2": n_ge2,
                "n_ge3": n_ge3,
                "triples": loc_triples,
                "best": local_best,
                "examples_run_ge2": examples,
                "seconds": round(time.time() - t1, 3),
            }
            chk.write(json.dumps(row) + "\n")
            chk.flush()
            reports.append(row)
            triples.extend(loc_triples)
            if local_max > max_run:
                max_run = local_max
                best = {"key": job["key"], **(local_best or {})}
            print(
                f"   max_run={local_max} ge2={n_ge2} ge3={n_ge3} "
                f"{row['seconds']}s",
                flush=True,
            )
    finally:
        chk.close()

    payload = {
        "search": "fat_image_hunt",
        "N": N,
        "K": K,
        "p_cut": P_CUT,
        "g_min": G_MIN,
        "n_regions": len(reports),
        "max_run": max_run,
        "n_triples": len(triples),
        "triples": triples,
        "best": best,
        "seconds": round(time.time() - t0, 3),
        "regions": reports,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        f"\nwrote {OUT}  regions={len(reports)}  max_run={max_run}  "
        f"triples={len(triples)}  {payload['seconds']}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
