#!/usr/bin/env python3
"""Regression tests for the family-wide termination certificates (Q14 §6).

results/termination_i{2..9}.json now cover every member: all 232 columns
k <= 30 certified -- for each, (x)_k - k!m is irreducible mod one prime q,
hence over Q, hence Jordan + Chebotarev give a killing prime of density
>= 1/k.  That upgrades ~230 columns from "a killing prime was found" to "a
killing prime HAD to exist", family-wide (2026-08-25; i=8 committed
2026-08-20 and reproduced content-identically before the rest were added).

Pins, live where it matters:
  1. all eight artifacts exist, i = 2..9, each with 29/29 certified, 232
     rows total, and N/K matching make_fam;
  2. EVERY row re-verified independently of termination_certificate.py:
     q prime (sympy), r == C(N,K) mod q recomputed through witness's
     lucas_mod_pure (a different code path from the generator), and
     (x)_k - k!·r irreducible mod q re-checked through sympy's GF(q)
     factorisation -- not by trusting the recorded 'certified' flag;
  3. the spot q-values quoted in the 2026-08-25 lab note: i=8 k=2 -> 227,
     k=3 -> 59; i=9 k=11 -> 631, k=29 -> 433;
  4. no row claims a dead prime: r != 0 mod q for every row (r = 0 would
     certify nothing -- 0 lies in every column image);
  5. the i=9 rows include the four Lucas-fill long-run columns' small
     ones (k = 11, 29): their killers now provably had to exist --
     existence, not location (the sub-sqrt(N) location question is Q31,
     and effective location bounds are Q17, blocked).

Runs in ~15 s (232 sympy irreducibility checks, degree <= 30, small q).

Run: python scripts/test_termination.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from bandii_kernel import make_fam  # noqa: E402
from witness import lucas_mod_pure  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


SPOT_Q = {(8, 2): 227, (8, 3): 59, (9, 11): 631, (9, 29): 433}

x = sp.Symbol("x")


def main() -> int:
    total = 0
    for i in range(2, 10):
        path = ROOT / "results" / f"termination_i{i}.json"
        expect(path.exists(), f"results/termination_i{i}.json exists")
        if not path.exists():
            continue
        art = json.loads(path.read_text(encoding="utf-8"))
        fam = make_fam(i)
        expect(art["i"] == i and art["N"] == fam.N and art["K"] == fam.K,
               f"i={i}: artifact N, K match make_fam")
        rows = art["rows"]
        expect(art["n_columns"] == 29 and art["n_certified"] == 29
               and len(rows) == 29,
               f"i={i}: 29/29 columns certified")
        expect(sorted(r["k"] for r in rows) == list(range(2, 31)),
               f"i={i}: rows cover exactly k = 2..30")
        bad = []
        for r in rows:
            k, q, rr = r["k"], r["q"], r["r"]
            if not r["certified"]:
                bad.append((k, "not certified"))
                continue
            if not sp.isprime(q):
                bad.append((k, f"q={q} not prime"))
                continue
            if rr % q == 0:
                bad.append((k, "r = 0: dead prime certifies nothing"))
                continue
            if lucas_mod_pure(fam.N, fam.K, q) != rr:
                bad.append((k, f"r != C(N,K) mod {q} by lucas_mod_pure"))
                continue
            kfact_r = (sp.factorial(k) * rr) % q
            poly = sp.Poly(sp.prod([x - j for j in range(k)]) - kfact_r,
                           x, modulus=q)
            fl = poly.factor_list()[1]
            if not (len(fl) == 1 and fl[0][1] == 1
                    and fl[0][0].degree() == k):
                bad.append((k, f"NOT irreducible mod {q}"))
        expect(not bad,
               f"i={i}: every certificate re-verified independently "
               f"(prime q, Lucas r, sympy GF(q) irreducibility)"
               + (f" -- failures: {bad[:3]}" if bad else ""))
        total += len(rows)
    expect(total == 232, f"232 rows family-wide (got {total})")
    for (i, k), q in SPOT_Q.items():
        art = json.loads((ROOT / "results" / f"termination_i{i}.json")
                         .read_text(encoding="utf-8"))
        row = {r["k"]: r for r in art["rows"]}[k]
        expect(row["q"] == q, f"spot pin: i={i} k={k} certified at q={q}")

    print("\n=== TERMINATION CERTIFICATE TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
