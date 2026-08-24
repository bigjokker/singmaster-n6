#!/usr/bin/env python3
"""Regression tests for results/ghost_census.json (Q9): it cannot go stale silently.

The census json went stale once -- committed with seven members and the
PRE-repair i=8 digest (b4c02030...), so its count was not a statement about
the tables on disk.  These pins compare the json LIVE against W.load of
every witness table, so any table change (repair, extension, new member)
fails the suite until the census is regenerated:

  1. the json lists exactly the members i = 2..9, one per witness table on
     disk -- no table is missing from the census and no census member lacks
     a table;
  2. for EVERY member, n_values and sha256 equal what W.load reads from
     results/i{i}_witness.npz right now (nothing is hardcoded: if a table
     is repaired again, this fails until the json is regenerated);
  3. n_values totals the members, k_range matches the tables' true min/max,
     and ghosts_found == 0;
  4. the three honesty caveats survive in the payload (not a targeted hunt;
     first-killing-prime only; a surviving column is an unresolved anomaly,
     not a ghost);
  5. the i=8 member is NOT the pre-repair table (digest b4c02030...).

The claim itself (c = k!*C(N,K) outside (x)_k(F_p) for the recorded (k,p))
is checked on samples by scripts/test_witness.py::
test_ghost_census_claim_holds; this suite deliberately does not duplicate
that.

Runs in ~3 s.

Run: python scripts/test_ghost_census.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import witness as W  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


def main() -> int:
    path = ROOT / "results" / "ghost_census.json"
    expect(path.exists(), "results/ghost_census.json exists")
    art = json.loads(path.read_text(encoding="utf-8"))

    tables = sorted(ROOT.glob("results/i*_witness.npz"))
    members = {m["i"]: m for m in art["members"]}
    expect(sorted(members) == list(range(2, 10)),
           f"members are exactly i = 2..9 (got {sorted(members)})")
    expect(len(art["members"]) == 8, "n_members == 8")
    expect(len(tables) == len(members),
           "one census member per witness table on disk, none missing")

    total = 0
    kmin, kmax = None, None
    for t in tables:
        ks, ps, meta = W.load(t)
        i = meta.get("i")
        expect(i in members, f"{t.name}: table has a census member")
        if i not in members:
            continue
        m = members[i]
        expect(m["n_values"] == int(ks.size),
               f"i={i}: census n_values {m['n_values']:,} == table rows {ks.size:,}")
        expect(m["sha256"] == meta["sha256"],
               f"i={i}: census sha256 matches the table on disk "
               f"({meta['sha256'][:12]}...)")
        expect(m["k_min"] == int(ks.min()) and m["k_max"] == int(ks.max()),
               f"i={i}: census k range matches the table")
        total += int(ks.size)
        kmin = int(ks.min()) if kmin is None else min(kmin, int(ks.min()))
        kmax = int(ks.max()) if kmax is None else max(kmax, int(ks.max()))

    expect(art["n_values"] == total,
           f"census total {art['n_values']:,} == sum over tables {total:,}")
    expect(art["k_range"] == [kmin, kmax],
           f"census k_range == tables' true range [{kmin}, {kmax:,}]")
    expect(art["ghosts_found"] == 0,
           "ghosts_found == 0 (the only possible answer on killed columns)")

    cav = " ".join(art.get("caveats", []))
    expect(len(art.get("caveats", [])) == 3, "three caveats present")
    expect("Not a targeted ghost hunt" in cav,
           "caveat: not a targeted hunt survives")
    expect("first killing prime" in cav,
           "caveat: first-killing-prime-only survives")
    expect("unresolved anomaly" in cav and "every prime" in cav,
           "caveat: surviving column is an anomaly, not a ghost")

    expect(not members[8]["sha256"].startswith("b4c02030"),
           "i=8 member is NOT the pre-repair table (digest b4c02030)")

    print("\n=== GHOST CENSUS TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
