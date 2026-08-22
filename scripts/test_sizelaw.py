#!/usr/bin/env python3
"""Regression tests for the size law and the escalation trigger.

Two things have to hold. The model must reproduce measurements it was not
fitted to -- there are no free parameters, so any drift is a real error. And
the trigger must FIRE: a trigger that never escalates is as useless as one
that escalates on run length, which is the bug this replaced.

Run: python scripts/test_sizelaw.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import sizelaw as S  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


def test_image_size_against_measurement() -> None:
    """Four image sizes measured at real fat-cell parameters (p odd, so g even
    <=> k odd). Zero fitted parameters: agreement is the whole claim."""
    rows = [
        ("(3,1) g_max", 475282, 2700967, 435720),
        ("(3,1) g_max-1", 475281, 2700967, 227550),
        ("(2,1) g_max-1", 262708, 3601237, 253496),
        ("(5,2) g_max", 131355, 1800619, 64397),
    ]
    worst = 0.0
    for _name, g, p, meas in rows:
        k = 3 if g % 2 == 0 else 4  # any odd / even k reproduces the parity
        worst = max(worst, abs(S.image_size(g, p, k) / meas - 1.0))
    expect(worst < 0.002, f"modelled |I| matches all four measurements to {worst:.2%}")


def test_proved_bound_is_never_violated() -> None:
    """image_bound is a theorem. Brute-force small cases; a breach is a bug."""
    breaches = []
    for p in (11, 29, 101, 211, 1009):
        for k in (2, 3, 4, 5, 8, 17):
            if k >= p:
                continue
            img = {math.comb(x, k) % p for x in range(p)}
            nonzero = len(img - {0})
            if nonzero > S.image_bound(p - k, k):
                breaches.append((p, k, nonzero, S.image_bound(p - k, k)))
    expect(not breaches, "proved image bound holds on every brute-forced (p,k)")
    if breaches:
        errors.append(f"  breaches: {breaches[:5]}")


def test_parity_fold_is_real() -> None:
    """Even k folds 2-to-1, odd k does not.

    Checked against the model, not a hand-picked ratio: at g ~ p the two
    parities tend to 1-1/e = 0.632 and 1-e^(-1/2) = 0.394 of p, a ratio of
    1.606, and finite p sits just under that. The doc's Theorem 3 table
    measures 0.6277/0.3939 at p = 10007.
    """
    p = 1009
    worst = 0.0
    frac = {}
    for k in (5, 6, 7, 8, 11, 12):
        img = {math.comb(x, k) % p for x in range(p)}   # includes 0 = C(x,k), x<k
        frac[k] = len(img) / p
        worst = max(worst, abs(S.image_size(p - k, p, k) / len(img) - 1.0))
    expect(worst < 0.02, f"exact image sizes match the model to {worst:.2%} at p=1009")
    expect(
        all(frac[o] > 1.5 * frac[e] for o, e in ((5, 6), (7, 8), (11, 12))),
        f"odd k has ~1.6x the image of even k: "
        f"{ {k: round(v, 4) for k, v in frac.items()} }",
    )
    endpoints = (1 - 1 / math.e, 1 - math.exp(-0.5))
    expect(
        abs(frac[11] - endpoints[0]) < 0.02 and abs(frac[12] - endpoints[1]) < 0.02,
        f"the two parity endpoints are 1-1/e and 1-e^-1/2 "
        f"({frac[11]:.4f}, {frac[12]:.4f} vs {endpoints[0]:.4f}, {endpoints[1]:.4f})",
    )


def test_preregistration_regenerates() -> None:
    """The recorded Band II pre-registration must fall out of the model."""
    from bandii_kernel import CAP, KMAX, KMIN, PRIMES
    from bandii_sweep import PREREGISTER

    rows = S.expected_alive(KMIN, KMAX, PRIMES[:CAP])
    worst, worst_at = 0.0, None
    for row in rows:
        pre = PREREGISTER.get(row["prime_index"])
        if not pre:
            continue
        rel = abs(row["alive"] - pre["alive"]) / pre["alive"]
        if rel > worst:
            worst, worst_at = rel, row["prime_index"]
        if abs(row["even_g_frac"] - pre["even"]) > 0.005:
            errors.append(
                f"even-g fraction pass {row['prime_index']}: model "
                f"{row['even_g_frac']:.3f} vs recorded {pre['even']:.3f}"
            )
    expect(
        worst < 0.01,
        f"Band II pre-registration regenerates (worst pass {worst_at}, {worst:.2%})",
    )


def test_poisson_tail() -> None:
    expect(abs(S.poisson_tail(1.0, 1) - (1 - math.exp(-1))) < 1e-12,
           "poisson_tail(1,1) = 1 - e^-1")
    expect(S.poisson_tail(0.5, 0) == 1.0, "poisson_tail(.,0) = 1")
    expect(S.poisson_tail(2.0, 20) < 1e-10, "poisson_tail is tiny far out")


def test_trigger_fires() -> None:
    """A trigger that never escalates is worthless. Make it escalate."""
    quiet = S.escalate(expected=0.5, observed=1)
    expect(not quiet["escalate"], "does not fire on an ordinary survivor (E=0.5, n=1)")

    rare = S.escalate(expected=0.005, observed=1)
    expect(
        rare["escalate"] and "below threshold" in rare["reason"],
        "FIRES when a column survives a round that expected almost none",
    )

    flood = S.escalate(expected=2.0, observed=25)
    expect(
        flood["escalate"] and "above the size law" in flood["reason"],
        "FIRES when far more columns survive than the law allows",
    )

    led = S.RoundLedger("test")
    led.record(1, 1000, 40.0, 38)
    led.record(2, 38, 0.004, 1)
    v = led.verdict()
    expect(
        v["escalate"] and v["escalating_rounds"] == [2],
        f"ledger escalates on the right round ({v['escalating_rounds']})",
    )


def test_run_length_is_not_the_criterion() -> None:
    """The bug this replaced: judging by run length across regimes.

    i=9 k=11 has a LONGER run than the i=8 fat-cell record column, but is
    three orders of magnitude more likely. Any criterion monotone in run
    length gets this pair backwards.
    """
    from bandii_kernel import D, cells, live_intervals, make_fam

    fam9 = make_fam(9)
    surv9, kill9, _ = S.live_run(fam9.N, fam9.K, 11, cap=20)
    lam9 = S.run_lambda(11, surv9)

    ivs = live_intervals(cells())
    chain8 = S.zjump_chain(2227205, ivs, D, 6)
    lam8 = S.run_lambda(2227205, chain8)

    expect(len(surv9) == 8 and kill9 == 449,
           f"i=9 k=11 reproduces: run {len(surv9)}, killed at {kill9}")
    expect(
        len(surv9) > len(chain8) and lam9 > 100 * lam8,
        f"the LONGER run (i=9 k=11, {len(surv9)}) is {lam9/lam8:,.0f}x more "
        f"likely than the shorter one (i=8 k=2227205, {len(chain8)})",
    )


def test_fat_cell_chain_matches_the_record() -> None:
    """docs/band-I.md: record k=2227205, kill 2701099. Rebuild it."""
    from bandii_kernel import D, cells, live_intervals

    ivs = live_intervals(cells())
    chain = S.zjump_chain(2227205, ivs, D, 9)
    expect(
        len(chain) > 6 and chain[6] == 2701099,
        f"reconstructed chain kills k=2227205 at 2701099 "
        f"after 6 survived primes (got {chain[:7]})",
    )


def test_round_model_tracks_a_real_sweep() -> None:
    """Round by round against i=7's actual Z-jump, from its checkpoint."""
    chk = ROOT / "results" / "i7_sweep.jsonl"
    if not chk.exists():
        ok.append("i7_sweep.jsonl absent; skipped round-model check")
        return
    from bandii_kernel import cells, first_live_after, live_intervals, make_fam
    from family_sweep import K_EXACT
    from witness import read_jsonl

    fam = make_fam(7)
    ivs = live_intervals(cells(fam), fam)
    rows = read_jsonl(chk)

    def surv_of(rnd):
        out = {}
        for r in rows:
            if str(r.get("tag", "")) == f"z{rnd}":
                for s in r.get("survivors") or []:
                    out[int(s["k"])] = int(s["k"]) + int(s["g"])
        return out

    alive = {k: k for k in range(K_EXACT[7] + 1, fam.K)}
    fired = []
    for rnd in range(1, 5):  # rounds with enough counts to be meaningful
        obs = surv_of(rnd)
        exp = 0.0
        for k, lp in alive.items():
            p = first_live_after(lp, ivs, fam.D)
            if p is not None:
                exp += S.survival(p - k, p, k)
        verdict = S.escalate(exp, len(obs))
        if verdict["escalate"]:
            fired.append((rnd, exp, len(obs)))
        if rnd == 1:
            expect(
                abs(len(obs) / exp - 1.0) < 0.05,
                f"i=7 Z round 1: model {exp:.0f} vs observed {len(obs)} "
                f"({len(obs)/exp:.2f}x)",
            )
        alive = obs
    expect(
        not fired,
        f"no false escalation on the clean i=7 sweep (fired: {fired})",
    )


def test_accuracy_bound_by_regime() -> None:
    """Q23: the model's error must stay inside the bound the docs claim.

    Two regimes with opposite characters, and the trigger is safe because they
    are anti-correlated: the model is loosest where the headroom is enormous,
    and near-exact where the headroom is thin. A change that made the model
    worse at Band II scale would be dangerous even though the same error is
    harmless at small-k scale, so both are pinned separately.

    Direction matters too: OVERestimating |I| inflates Lambda and can mask an
    anomaly. The proved involution bound caps that -- p(1-(1-1/p)^M) < M -- so
    this also asserts the model never exceeds the theorem.
    """
    def exact(k, p):
        vals = {0, 1}; v = 1
        for b in range(p - k - 1):
            v = v * (k + b + 1) % p * pow(b + 1, -1, p) % p
            vals.add(v)
        return len(vals)

    worst_small = 0.0
    for k, p in ((11, 191), (29, 283), (45, 239), (247, 347)):
        e = exact(k, p); m = S.image_size(p - k, p, k)
        worst_small = max(worst_small, m / e, e / m)
        expect(m < S.preimage_count(p - k, k) + 1e-9,
               f"model stays under the proved bound at k={k}, p={p}")
    expect(worst_small < 1.35,
           f"small-k census regime within its documented 1.263x bound "
           f"(measured {worst_small:.3f}x)")

    worst_big = 0.0
    for k, p in ((4126649, 5401853), (5182636, 5401853)):
        e = exact(k, p); m = S.image_size(p - k, p, k)
        worst_big = max(worst_big, m / e, e / m)
    expect(worst_big < 1.002,
           f"Band II regime within its documented 1.00014x bound "
           f"(measured {worst_big:.5f}x) -- this is the thin-headroom regime")


def test_survival_vec_matches_scalar() -> None:
    """`survival_vec` must be the SAME law as `survival`, elementwise.

    The profiler sweeps millions of columns and needs the vector form. If the
    two ever drift, the cost model and the escalation trigger would silently
    disagree about the same column -- so pin them together here rather than
    trusting that two copies of a formula stay equal.
    """
    import numpy as np
    rng = np.random.default_rng(4)
    g = rng.integers(0, 5_000_000, size=4000, dtype=np.int64)
    k = rng.integers(2, 40_000_000, size=4000, dtype=np.int64)
    p = (k + g).astype(np.int64)
    got = S.survival_vec(g, p, k)
    want = np.array([S.survival(int(a), int(b), int(c))
                     for a, b, c in zip(g, p, k)])
    worst = float(np.max(np.abs(got - want)))
    if worst > 1e-12:
        errors.append(f"survival_vec differs from survival by {worst:.3g}")
    else:
        ok.append(f"survival_vec == survival on 4000 random (g,p,k), max |diff| {worst:.3g}")

    mg = S.preimage_count_vec(g, k)
    mw = np.array([S.preimage_count(int(a), int(c)) for a, c in zip(g, k)])
    if not np.array_equal(mg, mw):
        errors.append("preimage_count_vec differs from preimage_count")
    else:
        ok.append("preimage_count_vec == preimage_count on the same 4000")

    # g <= 0 is the degenerate branch the scalar form special-cases
    z = S.survival_vec(np.array([0, -3]), np.array([101, 101]), np.array([5, 6]))
    if not np.allclose(z, 1.0):
        errors.append(f"survival_vec mishandles g<=0: {z}")
    else:
        ok.append("survival_vec returns 1.0 for g <= 0, matching the scalar branch")


def test_poisson_tail_survives_band_ii_scale() -> None:
    """The trigger must be able to fire where the counts are large.

    poisson_tail used to start from term = math.exp(-expected), which
    underflows to 0.0 above ~745. The cdf then accumulated nothing and the tail
    returned 1.0 -- "ordinary" -- for every round with expected above that.
    Measured before the fix:

        escalate(expected=102,563, observed=500,000) -> escalate: False

    A five-fold excess of survivors read as ordinary, and that covered Band II
    passes 1-4 and Z-jump rounds 1-2 at every member. The anomaly detector was
    inoperative precisely where almost every column lives.
    """
    import math

    def brute(lam, k):
        """Direct summation. Ground truth ONLY below the underflow cliff."""
        t = math.exp(-lam)
        cdf = t
        for i in range(1, k):
            t *= lam / i
            cdf += t
        return max(0.0, min(1.0, 1.0 - cdf))

    worst = 0.0
    for lam in (0.001, 0.5, 1, 5, 20, 100, 300, 700):
        for k in (1, 2, 5, 10, 50, 120, 400, 900):
            got, want = S.poisson_tail(lam, k), brute(lam, k)
            if want > 1e-12:
                worst = max(worst, abs(got - want) / want)
    if worst > 1e-6:
        errors.append(f"poisson_tail disagrees with direct summation by {worst:.3g}")
    else:
        ok.append(f"poisson_tail matches direct summation below the cliff "
                  f"(worst rel. diff {worst:.2g})")

    # the regression itself
    fired = S.escalate(102563.32, 500000)
    if not fired["escalate"]:
        errors.append("a 5x survivor excess at Band II scale still reads as ordinary")
    else:
        ok.append("a 5x survivor excess at Band II scale now escalates "
                  f"(tail={fired['poisson_tail']:.3g})")

    quiet = S.escalate(102563.32, 102754)
    if quiet["escalate"]:
        errors.append(f"an ordinary Band II pass now escalates spuriously: {quiet}")
    else:
        ok.append(f"an ordinary Band II pass stays ordinary "
                  f"(tail={quiet['poisson_tail']:.3g})")

    # monotone in observed, for a fixed expected
    tails = [S.poisson_tail(5000.0, n) for n in (4000, 4900, 5000, 5100, 6000)]
    if any(b > a + 1e-12 for a, b in zip(tails, tails[1:])):
        errors.append(f"poisson_tail is not monotone decreasing in observed: {tails}")
    else:
        ok.append("poisson_tail is monotone decreasing in observed")


def test_recorded_bandii_passes_stay_ordinary() -> None:
    """The fix must make the trigger WORK without making it trigger-happy.

    Every recorded Band II pass of i=8 was, by the project's own reading,
    ordinary and on the pre-registered curve. If the repaired tail escalates
    any of them, the model and the trigger now disagree about data that has
    already been accepted -- which would be a worse failure than the underflow.
    """
    rows = [(102600.0, 102754), (12600.0, 12478), (1816.0, 1782),
            (290.0, 293), (49.4, 48), (8.78, 7), (1.61, 1), (0.30, 0)]
    fired = [(e, o) for e, o in rows if S.escalate(e, o)["escalate"]]
    if fired:
        errors.append(f"repaired tail escalates recorded-ordinary i=8 Band II "
                      f"passes: {fired}")
    else:
        ok.append(f"all {len(rows)} recorded i=8 Band II passes still read "
                  f"ordinary under the repaired tail")


def main() -> int:
    test_image_size_against_measurement()
    test_survival_vec_matches_scalar()
    test_poisson_tail_survives_band_ii_scale()
    test_recorded_bandii_passes_stay_ordinary()
    test_accuracy_bound_by_regime()
    test_proved_bound_is_never_violated()
    test_parity_fold_is_real()
    test_preregistration_regenerates()
    test_poisson_tail()
    test_trigger_fires()
    test_run_length_is_not_the_criterion()
    test_fat_cell_chain_matches_the_record()
    test_round_model_tracks_a_real_sweep()
    print("\n=== SIZE LAW TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
