#!/usr/bin/env python3
"""Regression tests for the image-scan kernel.

The half-scan is an exact identity, not an approximation, so the bar is
byte-identical output against the reference full scan -- survivor set and
reported b alike. Anything less would silently change every recorded
certificate and break comparability with past runs.

Run: python scripts/test_kernel.py
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import bandii_kernel as K  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


def test_half_equals_full() -> None:
    """Same survivors, same b, over every column at several primes."""
    rng = random.Random(1)
    bad = 0
    n = 0
    for p in (11, 29, 101, 211, 1009, 5003, 10007):
        F = K.fact_table(p)
        ks = list(range(2, p))
        for _ in range(40):
            r = rng.randrange(1, p)
            n += 1
            if K.scan_ks_full(F, p, r, ks) != K.scan_ks_half(F, p, r, ks):
                bad += 1
    expect(not bad, f"half-scan == full scan on {n} (p,r) cases, all columns")


def test_against_exact_arithmetic() -> None:
    """Both scans against math.comb, including the reported b."""
    p = 211
    F = K.fact_table(p)
    wrong = 0
    for r in range(1, p):
        got = {x["k"]: x["b"] for x in K.scan_ks_half(F, p, r, range(2, p))}
        for k in range(2, p):
            hits = [b for b in range(p - k) if math.comb(k + b, k) % p == r]
            if bool(hits) != (k in got):
                wrong += 1
            elif hits and got[k] != hits[0]:
                wrong += 1
    expect(not wrong, f"half-scan matches math.comb over every (r,k) at p={p}")


def test_reported_b_is_the_least_hit() -> None:
    """b must be the FIRST full-scan hit, not merely some hit.

    For odd k the -s branch finds half-index b whose full-scan index is
    g-1-b, which decreases as b grows -- so the last -s hit gives the
    smallest full index. Getting that backwards passes a survivor-set test
    and still records the wrong witness.
    """
    p = 1009
    F = K.fact_table(p)
    checked = 0
    wrong = []
    for r in (3, 17, 200, 511, 1008):
        for row in K.scan_ks_half(F, p, r, range(2, 200)):
            k, b, g = row["k"], row["b"], row["g"]
            hits = [x for x in range(g) if math.comb(k + x, k) % p == r]
            checked += 1
            if not hits or b != hits[0]:
                wrong.append((k, r, b, hits[:3]))
    expect(
        not wrong and checked > 50,
        f"reported b is the least full-scan hit ({checked} survivors checked)",
    )
    if wrong:
        errors.append(f"  wrong b: {wrong[:3]}")


def test_odd_and_even_k_both_covered() -> None:
    """Both parities must actually be exercised, not just one."""
    p = 1009
    F = K.fact_table(p)
    par = {0: 0, 1: 0}
    for r in range(1, 60):
        for row in K.scan_ks_half(F, p, r, range(2, 300)):
            par[row["k"] % 2] += 1
    expect(par[0] > 20 and par[1] > 20,
           f"both k parities exercised (even {par[0]}, odd {par[1]})")


def test_r_zero_is_refused() -> None:
    """r=0 lies in every image, so it certifies nothing.

    The mask scans cannot see that -- a factorial table has no zero entry --
    so they would report "no survivors", i.e. every column killed. That is a
    false certificate, so both scans must refuse instead.
    """
    F = K.fact_table(101)
    for name, fn in (("full", K.scan_ks_full), ("half", K.scan_ks_half)):
        try:
            fn(F, 101, 0, [3, 4, 5])
            errors.append(f"{name} scan accepted r=0 instead of refusing")
        except ValueError:
            ok.append(f"{name} scan refuses r(p)=0 rather than reporting all-killed")


def test_straggler_benchmark() -> None:
    """The four near-K stragglers, with their recorded witness indices."""
    F = K.fact_table(K.P1)
    got = {h["k"]: h["b"] for h in K.scan_ks_half(F, K.P1, K.R_P1, K.STRAGGLERS.keys())}
    expect(got == K.STRAGGLERS,
           f"half-scan reproduces the straggler witnesses at p={K.P1}")
    if got != K.STRAGGLERS:
        errors.append(f"  got {got} want {K.STRAGGLERS}")


def test_toggle_round_trips() -> None:
    """USE_HALF_SCAN must select the implementation, for A/B checking."""
    F = K.fact_table(1009)
    ks = list(range(2, 400))
    old = K.USE_HALF_SCAN
    try:
        K.USE_HALF_SCAN = True
        a = K.scan_ks(F, 1009, 77, ks)
        K.USE_HALF_SCAN = False
        b = K.scan_ks(F, 1009, 77, ks)
    finally:
        K.USE_HALF_SCAN = old
    expect(a == b, "USE_HALF_SCAN toggles implementation, output unchanged")
    expect(K.USE_HALF_SCAN is True, "half-scan is the default")


def main() -> int:
    test_half_equals_full()
    test_against_exact_arithmetic()
    test_reported_b_is_the_least_hit()
    test_odd_and_even_k_both_covered()
    test_r_zero_is_refused()
    test_straggler_benchmark()
    test_toggle_round_trips()
    print("\n=== KERNEL TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
