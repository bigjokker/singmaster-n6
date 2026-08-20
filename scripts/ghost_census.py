#!/usr/bin/env python3
"""Q9: the census read as a test of the falling-factorial local-global question.

WHAT A KILLED COLUMN CERTIFIES. Column k of m = C(N,K) is killed by p when
m is outside the image of C(x,k) mod p. Multiplying through by k!, that says

    c_k := k! * m   is not in   (x)_k (F_p),

hence not in (x)_k(Z_p), hence c_k is NOT in the intersection over all primes.
So every certificate in a witness table is a certified counter-instance to
"c is locally a falling factorial everywhere" -- a certified NON-ghost.

A ghost is a c that lies in (x)_k(Z_p) for every p but is not (a)_k for any
integer a. Whether any exists is open; it is the only logical gap between
this project's method and a guarantee that the method always terminates
(see Q14). Every column this project kills is one more c that is not one.

WHAT IT IS WORTH, HONESTLY. The count is large and the c are structured --
k! times one fixed enormous binomial coefficient, with k sweeping six orders
of magnitude -- and structured families are where a ghost would have to live,
since a random c is overwhelmingly not a ghost for cheap reasons. But this
census was not designed to hunt ghosts: it tests the c that the Singmaster
search happens to produce, and it stops at the first prime that kills, so it
establishes "not a ghost" rather than probing how close any c came. Read it
as a broad negative result, not as a targeted search.

    python scripts/ghost_census.py
    python scripts/ghost_census.py --json_out results/ghost_census.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import witness as W  # noqa: E402


def census_one(path: Path) -> dict:
    import numpy as np

    ks, ps, meta = W.load(path)
    g = ps - ks
    return {
        "file": path.name,
        "i": meta.get("i"),
        "N": meta["N"],
        "K": meta["K"],
        "n_values": int(ks.size),
        "k_min": int(ks.min()),
        "k_max": int(ks.max()),
        "k_odd_frac": round(float((ks % 2 == 1).mean()), 4),
        "p_min": int(ps.min()),
        "p_max": int(ps.max()),
        "g_median": int(np.median(g)),
        "sum_g": int(g.sum()),
        "sha256": meta["sha256"],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    rows = []
    for f in sorted(args.results.glob("i*_witness.npz")):
        try:
            rows.append(census_one(f))
        except Exception as exc:
            print(f"  {f.name}: unreadable ({exc})", flush=True)
    if not rows:
        print("  no witness tables found; run family_sweep first", flush=True)
        return 1
    rows.sort(key=lambda r: (r["i"] is None, r["i"]))

    total = sum(r["n_values"] for r in rows)
    print("  Certified non-ghosts: values c = k! * C(N,K) shown to lie outside")
    print("  (x)_k(F_p) for an explicit prime p.\n")
    print(f"  {'member':>7} {'columns':>12} {'k range':>22} {'p range':>25}")
    for r in rows:
        print(f"  {'i=' + str(r['i']):>7} {r['n_values']:>12,} "
              f"{str(r['k_min']) + '..' + str(r['k_max']):>22} "
              f"{str(r['p_min']) + '..' + str(r['p_max']):>25}")
    print(f"  {'TOTAL':>7} {total:>12,}")

    kmin = min(r["k_min"] for r in rows)
    kmax = max(r["k_max"] for r in rows)
    print(f"\n  distinct c values certified : {total:,}")
    print(f"  degrees k covered           : {kmin} .. {kmax:,}")
    print(f"  ghosts found                : 0")
    print(f"  (0 is the only possible answer here: a surviving column would be")
    print(f"   an unresolved anomaly, not a ghost -- it would have to survive")
    print(f"   EVERY prime to be one, and the sweeps stop at a cap.)")

    payload = {
        "search": "ghost_census",
        "claim": (
            "For each (k, p) recorded, c = k! * C(N,K) is not in (x)_k(F_p), "
            "hence not in the intersection over all primes, hence not a ghost."
        ),
        "caveats": [
            "Not a targeted ghost hunt: these are the c the Singmaster search "
            "produces, not candidates chosen for being ghost-like.",
            "Each column stops at its first killing prime, so the census "
            "records that c fails locally somewhere, not how nearly it passed.",
            "A column surviving its cap would be an unresolved anomaly, not a "
            "ghost; ghosthood needs failure at every prime.",
        ],
        "n_values": total,
        "k_range": [kmin, kmax],
        "ghosts_found": 0,
        "members": rows,
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\n  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
