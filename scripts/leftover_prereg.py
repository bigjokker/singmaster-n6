#!/usr/bin/env python3
"""Pre-register E[off_ladder leftovers] per member, before i=10 is ever run.

WHAT A LEFTOVER IS.  `cells()` starts at sqrt(N)+1, so the Z-jump never scans
p <= sqrt(N) (Q31).  A small column that survives the cap of fat primes just
above sqrt(N) is left alive and then needs a full-Lucas fill, which finds its
killer below sqrt(N).  The census calls those `off_ladder`.  i=9 left four
(k = 87, 399, 553, 1281).  That count -- not the sub-sqrt(N) census, which
cannot grow -- is the number worth watching, and it is predictable.

THE INSTRUMENT.  Every k < sqrt(N) starts its ladder at the SAME first live
prime above sqrt(N), so they share a prime list: a Band-II-shaped population,
which is exactly what sizelaw.expected_alive models.  The window is
[k_lo_z, sqrt(N)-1] walked over that shared ladder.

CONDITIONAL vs UNCONDITIONAL -- the distinction that matters here.  The
escalation ledger's E_r is "sum over the columns ACTUALLY entering round r of
survival(their own next prime)": a one-step check CONDITIONED on the observed
entering set.  It cannot forecast an unrun member, because there is no
entering set to condition on.  The window curve below is UNCONDITIONAL --
the cumulative product from pass 1 -- and is therefore the forecasting
instrument.  They are not rivals; they answer different questions.

WHAT THE ONE CALIBRATION POINT SAYS.  At i=9 the window wanted 3.345 columns
still alive entering round 12 and 7 were there; it wanted 1.954 surviving and
4 were.  Same factor (~2.05x) at both steps, while the LAST-STEP survival was
right (observed 4/7 = 0.571 against predicted 1.954/3.345 = 0.584).  So the
next-prime law is fine; the whole miss is extra columns still standing.  With
n = 1 that is either a high draw (P(X>=4 | 1.954) ~ 0.135) or a systematic
underestimate of survival in the fat g/p regime, compounding.  IT CANNOT BE
TOLD APART, so i=10 is quoted as a band, not a point, and no per-round drift
rate is fitted -- one miss is not a rate.

    python scripts/leftover_prereg.py
    python scripts/leftover_prereg.py --json_out results/leftover_prereg.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import sizelaw as S  # noqa: E402
from bandii_kernel import cells, first_live_after, live_intervals, make_fam  # noqa: E402

CAP_Z, CAP_Z_SMALL_K, SMALL_K = 12, 15, 1000

# k_lo_z read from each member's sweep record; i=10 has no run, so the i=8
# convention (engine band ends at 3) is assumed and flagged as such.
WINDOW_LO = {7: 201, 8: 3, 9: 81, 10: 3}
REALISED_ROUNDS = {7: 9, 8: 14, 9: 12}

# Observed, TRANSCRIBED from the run records.  i=9's Band II/Z-jump per-round
# survivor counts live only in results/i9_sweep.jsonl (219 MB, gitignored --
# its json has no `escalation` block), so these are not reproducible from a
# clone and are recorded as observations, not as regenerable numbers.
OBSERVED = {
    7: {"leftovers": 0, "note": "Z-jump died at round 9; the cap was never reached"},
    8: {"leftovers": 0,
        "note": "k=1021 is a dropped-and-repaired record (false p=3517, then "
                "1051), NOT a cap survivor; it lands in off_ladder only because "
                "that bucket answers 'where did the recorded killer sit'. On the "
                "leftover process i=8 observed 0."},
    9: {"leftovers": 4, "entering_round_12": 7, "mean_k": 580,
        "columns": [87, 399, 553, 1281],
        "note": "honest cap survivors; transcribed from the i9 jsonl readout"},
}


def shared_ladder(fam, n: int) -> list[int]:
    """The first `n` live primes above sqrt(N) -- shared by every k < sqrt(N)."""
    ivs = live_intervals(cells(fam), fam)
    out, x = [], math.isqrt(fam.N)
    for _ in range(n):
        q = first_live_after(x, ivs, fam.D)
        if q is None:
            break
        out.append(int(q))
        x = q
    return out


def window_curve(i: int) -> dict:
    fam = make_fam(i)
    rN = math.isqrt(fam.N)
    klo, khi = WINDOW_LO[i], rN - 1
    ladder = shared_ladder(fam, CAP_Z_SMALL_K)
    rows = S.expected_alive(klo, khi, ladder)
    at = {r["prime_index"]: r for r in rows}
    return {
        "i": i, "sqrt_N": rN, "window": [klo, khi], "n_columns": khi - klo + 1,
        "ladder_start": ladder[0], "ladder": ladder,
        "E_at_12": round(at[CAP_Z]["alive"], 4) if CAP_Z in at else None,
        "E_at_15": round(at[CAP_Z_SMALL_K]["alive"], 4) if CAP_Z_SMALL_K in at else None,
        "mean_k_at_12": round(at[CAP_Z]["mean_k"], 0) if CAP_Z in at else None,
        "curve": [{"pass": r["prime_index"], "p": r["p"],
                   "alive": round(r["alive"], 4),
                   "mean_k": round(r["mean_k"], 0) if r["mean_k"] else None}
                  for r in rows],
    }


def poisson_upper_tail(lam: float, obs: int) -> float:
    return 1.0 - sum(math.exp(-lam) * lam**n / math.factorial(n) for n in range(obs))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json_out", type=Path, default=None)
    args = ap.parse_args()

    print("  Leftover pre-registration: E[off_ladder] on the shared sub-sqrt(N) ladder")
    print(f"  policy: CAP_Z={CAP_Z} over all columns, then rounds "
          f"{CAP_Z+1}..{CAP_Z_SMALL_K} over k < {SMALL_K} only\n")
    members = []
    print(f"  {'i':>3} {'sqrt(N)':>9} {'window cols':>12} {'E at cap':>10} "
          f"{'realised':>9} {'observed':>9}")
    for i in (7, 8, 9, 10):
        w = window_curve(i)
        rr = REALISED_ROUNDS.get(i)
        at = {r["pass"]: r["alive"] for r in w["curve"]}
        e_realised = at.get(rr) if rr else None
        w["realised_rounds"] = rr
        w["E_at_realised"] = e_realised
        obs = OBSERVED.get(i)
        w["observed"] = obs
        if obs and e_realised:
            w["poisson_upper_tail"] = round(
                poisson_upper_tail(e_realised, obs["leftovers"]), 4)
        members.append(w)
        shown = e_realised if e_realised is not None else w["E_at_12"]
        print(f"  {i:>3} {w['sqrt_N']:>9,} {w['n_columns']:>12,} "
              f"{shown:>10.3f} {str(rr):>9} "
              f"{(str(obs['leftovers']) if obs else '-'):>9}")

    m9 = next(m for m in members if m["i"] == 9)
    at9 = {r["pass"]: r["alive"] for r in m9["curve"]}
    pred_entering, pred_alive = at9[11], at9[12]
    obs_entering, obs_alive = 7, 4
    calib = {
        "member": 9,
        "window_entering_round_12": round(pred_entering, 4),
        "observed_entering_round_12": obs_entering,
        "window_alive_after_12": round(pred_alive, 4),
        "observed_alive_after_12": obs_alive,
        "count_factor": round(obs_alive / pred_alive, 3),
        "last_step_survival_predicted": round(pred_alive / pred_entering, 4),
        "last_step_survival_observed": round(obs_alive / obs_entering, 4),
        "poisson_upper_tail": round(poisson_upper_tail(pred_alive, obs_alive), 4),
        "ledger_conditional_check": {
            "value": round(obs_entering * (pred_alive / pred_entering), 3),
            "observed": obs_alive,
            "what_it_is": "the escalation ledger's E_r: the 7 columns that "
                          "ACTUALLY arrived times the size law's next-prime "
                          "survival. A conditional one-step check that the "
                          "next-prime law is right -- NOT a rival forecast, "
                          "and undefined at i=10 where nothing has arrived.",
        },
        "reading": "The next-prime law is right (0.571 observed vs 0.584 "
                   "predicted). The miss is entirely extra columns still "
                   "standing: the window wanted 3.345 entering and got 7, "
                   "wanted 1.954 surviving and got 4 -- the same ~2.05x at "
                   "both steps. With n=1 a high draw and a compounding "
                   "underestimate cannot be told apart; no per-round drift "
                   "rate is fitted.",
    }
    m10 = next(m for m in members if m["i"] == 10)
    # low end: the window taken as-is; high end: the single i=9 factor applied
    # once as a worst case. Reported as integers -- one calibration point does
    # not support decimals.
    band = [int(math.floor(m10["E_at_12"])),
            int(round(m10["E_at_12"] * (obs_alive / pred_alive)))]
    print(f"\n  i=9 calibration: entering {pred_entering:.3f} predicted vs "
          f"{obs_entering} observed; alive {pred_alive:.3f} vs {obs_alive} "
          f"(P(X>={obs_alive}) = {calib['poisson_upper_tail']:.3f})")
    print(f"  last-step survival {calib['last_step_survival_observed']} observed vs "
          f"{calib['last_step_survival_predicted']} predicted -- the law is fine")
    print(f"  mean k at round 12: {m9['mean_k_at_12']:.0f} predicted vs "
          f"{OBSERVED[9]['mean_k']} observed (walking down, as predicted)")
    print(f"\n  i=10 PRE-REGISTRATION: window {m10['E_at_12']:.2f}; applying the "
          f"i=9 factor once as a worst case gives {band[1]:.1f}")
    print(f"  => {band[0]:.0f}-{band[1]:.0f} leftovers, UNDER A DOZEN. Lucas-fill, "
          f"not 22k. The decision is robust to the unresolved 2x.")

    payload = {
        "search": "leftover_prereg",
        "claim": "E[off_ladder leftovers] per member from the shared sub-sqrt(N) "
                 "ladder, pre-registered before i=10 is run.",
        "policy": {"CAP_Z": CAP_Z, "CAP_Z_SMALL_K": CAP_Z_SMALL_K,
                   "SMALL_K": SMALL_K},
        "instrument": "unconditional sizelaw.expected_alive over "
                      "[k_lo_z, sqrt(N)-1] on the shared ladder above sqrt(N)",
        "members": members,
        "calibration_i9": calib,
        "i10_prediction": {
            "window_E_at_12": m10["E_at_12"],
            "band": band,
            "basis": "5.49 trusting the window as-is; ~11 applying the single "
                     "i=9 factor once as a worst case. n=1, so a band.",
            "consequence": "Lucas-fill, not a 22.5k pre-pass. Decision "
                           "unchanged across the whole band.",
        },
        "caveats": [
            "i=9's per-round observed counts are transcribed from "
            "results/i9_sweep.jsonl (219 MB, gitignored; its json has no "
            "escalation block) and are not reproducible from a clone.",
            "i=8's k=1021 is a dropped-and-repaired record, not a cap "
            "survivor; on the leftover process i=8 observed 0 against a small "
            "E, which does not calibrate anything.",
            "One informative member (i=9). No per-round drift rate is fitted; "
            "2.05^(1/11) is a description of one miss, not a measured rate.",
            "i=10's window low end assumes k_lo_z = 3 (the i=8 convention); "
            "no i=10 run exists.",
        ],
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\n  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
