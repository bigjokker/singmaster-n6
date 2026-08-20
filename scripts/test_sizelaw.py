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


def main() -> int:
    test_image_size_against_measurement()
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
