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


def test_gm_scan_equals_numpy_scan() -> None:
    """The division-free scan must be byte-identical to the numpy one.

    numpy does the test as (s*F) %% p, and int64 %% compiles to scalar idiv --
    there is no SIMD integer division on x86, which is why production sat at
    ~0.167 Gelem/s. Granlund-Montgomery replaces it: with pinv = p^-1 mod 2^64,
    p | d iff d*pinv <= floor((2^64-1)/p). One multiply, one compare.

    MEASURED end to end at i=8, three arms, byte-identical survivors in every
    one: numpy@8 wall 772.6 s, GM@8 168.8 s (4.58x), GM@16 163.6 s. But a
    faster kernel that moved a single survivor would be worthless, so what is
    pinned here is EQUALITY, not speed -- every record, including the witness
    index b, and both k parities (the odd branch scans the -s side downward,
    since a hit at half-index j is a full-scan hit at g-1-j).
    """
    saved = K.USE_GM_SCAN
    try:
        for i in (5, 6, 7):
            fam = K.make_fam(i)
            kmax, _ = K.kmax_of(fam)
            p = int(K.first_primes_above(fam.N2, fam.D, kmax, n=1)[0])
            r = K.r_closed(p, N=fam.N, K=fam.K)
            lo = fam.K + 2
            ks = [k for k in range(lo, min(lo + 90, kmax + 1))]
            if not ks:
                continue
            K.USE_GM_SCAN = False
            want = K.scan_ks_windowed(p, r, ks)
            K.USE_GM_SCAN = True
            got = K.scan_ks_windowed(p, r, ks)
            odd = sum(1 for d in want if d["k_odd"])
            expect(got == want,
                   f"GM scan == numpy scan at i={i} "
                   f"({len(want)} survivors of {len(ks)} columns, {odd} odd-k)")
            if got != want:
                bad = [(a, b) for a, b in zip(got, want) if a != b][:2]
                errors.append(f"  first divergence: {bad}")
            expect(any(d["k_odd"] for d in want) or not want,
                   f"i={i} sample exercises the odd-k negated branch "
                   f"({odd} odd-k survivors)")
    finally:
        K.USE_GM_SCAN = saved

    # and it must refuse r(p)=0 exactly as the numpy path does
    K.USE_GM_SCAN = True
    try:
        fam = K.make_fam(5)
        kmax, _ = K.kmax_of(fam)
        p = int(K.first_primes_above(fam.N2, fam.D, kmax, n=1)[0])
        try:
            K.scan_ks_windowed(p, 0, [fam.K + 2])
            expect(False, "GM path refuses r(p)=0")
        except ValueError as exc:
            expect("r(p)=0" in str(exc) and "kills nothing" in str(exc),
                   f"GM path refuses r(p)=0 like the numpy path ({str(exc)[:46]})")
    finally:
        K.USE_GM_SCAN = saved


def test_iv_starts_survives_id_reuse() -> None:
    """first_live_after must not answer from a cache keyed by a dead id.

    _iv_starts cached interval starts under id(ivs) with only a len() guard
    and held no reference to the list, so once the list was freed CPython
    could hand the same address to a different list of the same length --
    and the bisect then ran on the OLD starts. The wrong answer is a wrong
    live prime (or None) for every column that asks; it is masked today
    only because no two families have the same interval count.
    """
    def ask(ivs):
        return K.first_live_after(0, ivs, 10**6)

    hits = 0
    wrong = []
    for trial in range(200):
        a = [(100 + 10 * j, 105 + 10 * j) for j in range(50)]
        ask(a)
        id_a = id(a)
        del a
        b = [(5000 + 10 * j, 5005 + 10 * j) for j in range(50)]
        if id(b) == id_a:
            hits += 1
        got = ask(b)
        if got != 5003:
            wrong.append((trial, got))
        expect_starts = [lo for lo, _hi in b]
        if K._iv_starts(b) != expect_starts:
            wrong.append((trial, "stale starts"))
        del b
    expect(not wrong,
           f"first_live_after is correct across {200} freed-and-reallocated "
           f"interval lists ({hits} id reuses observed, {len(wrong)} wrong "
           f"{wrong[:2]})")


def test_scanners_refuse_primes_past_the_word_size() -> None:
    """The int64 kernel needs s*F[b] < p^2 < 2^63; above it numpy wraps
    silently and the scan can report a FALSE KILL. The bound was stated in
    a docstring and enforced nowhere.
    """
    import gmpy2

    ok_p = 37_024_873                       # i=9 Band II p1
    bad_p = int(gmpy2.next_prime(3_037_000_499))
    expect(bad_p * bad_p >= 2**63 and ok_p * ok_p < 2**63,
           f"fixture: {bad_p} is past the int64 bound, {ok_p} inside it")
    expect(K._check_r(5, ok_p) == 5, "a prime inside the bound is accepted")
    checks = [
        ("_check_r", lambda: K._check_r(5, bad_p)),
        ("fact_window", lambda: K.fact_window(bad_p, 0, 3)),
        ("scan_ks_windowed", lambda: K.scan_ks_windowed(bad_p, 5, [bad_p - 50])),
    ]
    refused = 0
    for name, fn in checks:
        try:
            fn()
            expect(False, f"{name} accepted p={bad_p} past the int64 bound")
        except ValueError as exc:
            refused += 1
            expect("2^63" in str(exc) or "int64" in str(exc),
                   f"{name} refuses p={bad_p} past the int64 bound")
    # fact_table(bad_p) without a guard would allocate 24 GB and loop 3e9
    # times; only probe it once the cheap entry points are known to refuse.
    if refused == len(checks):
        try:
            K.fact_table(bad_p)
            expect(False, f"fact_table accepted p={bad_p} past the int64 bound")
        except ValueError:
            expect(True, f"fact_table refuses p={bad_p} before allocating")
    else:
        expect(False, "fact_table bound not probed: the cheap entry points do not refuse")


def main() -> int:
    test_iv_starts_survives_id_reuse()
    test_scanners_refuse_primes_past_the_word_size()
    test_gm_scan_equals_numpy_scan()
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
