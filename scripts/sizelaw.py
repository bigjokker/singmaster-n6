#!/usr/bin/env python3
"""The size law, and the escalation trigger built on it.

A column that survives its prime cap is an anomaly only if surviving was
unlikely. Run *length* does not measure that, because the survival rate per
live prime varies by two orders of magnitude across the regimes this project
scans:

    regime              g/p        survival per live prime
    small k             up to 0.94   up to 0.61
    fat cell / Band II  <= 0.236     <= 0.16

So a run of 8 at i=9 k=11 (expected count 65 over the small-k columns alone)
is four orders of magnitude *less* surprising than a run of 6 in an i=8 fat
cell (expected count ~1). The old pre-registered trigger -- "escalate on a
run of 8" -- fires on the first and not the second, i.e. exactly backwards.

The corrected trigger is on expected count, not length. Per column,

    Lambda = prod_j |I_{p_j,k}| / p_j          (probability of that run)

but Lambda alone is not the verdict: it has to be weighed against how many
columns had a chance to do it. Multiplying by a whole phase's column count is
wrong by orders of magnitude -- most Band I columns sit at g ~ 9 and die at
once, while the fat cells sit at g ~ 5e5 -- so the multiplier is computed
round by round instead:

    E_r = sum over the columns ACTUALLY entering round r of
          survival(their own next prime)
    escalate iff (observed_r >= 1 and E_r < THRESHOLD)
              or (Poisson upper tail P(X >= observed_r | E_r) < TAIL_ALPHA)

The sweep has those columns and their primes in hand when it builds each
round's buckets, so the expectation costs nothing and needs no hand-tuned
regime size.

THE MODEL. I_{p,k} is the value set of f_k(x) = (x)_k / k! over F_p. The
involution (k-1-x)_k = (-1)^k (x)_k gives an exact pre-image count:

    k odd   no fold          M = g + 1          (|I| <= g+1)
    k even  exact 2-to-1     M = ceil(g/2) + 1  (|I| <= ceil(g/2)+1)

Those bounds are theorems, uniform in p and k. The only further loss is
birthday collisions in F_p, so |I| ~ p(1 - (1-1/p)^M). Zero fitted parameters.

The parity that matters is k's, not g's. They are equivalent for odd p
(g = p-k even <=> k odd), but k-parity is the invariant: it is a property of
the polynomial f_k, fixed once and for all, while g changes with every prime.

ACCURACY, BY REGIME (Q23; test_sizelaw.py pins these). The model's error and
the trigger's headroom are anti-correlated, in the safe direction:

    regime                 p        model/exact   smallest expected   headroom
    Band II / Z-jump   ~5.4e6         1.00014x         0.206            21x
    small-k census     ~2e2           1.263x          15.3            1527x

So the model is loosest exactly where the trigger has three orders of
magnitude of room, and is accurate to 0.014% exactly where the room is thin.
Measured: 0 verdict flips over all 102 columns of the i=9 small-k census when
Lambda is recomputed from EXACT image sizes; a flip there would need Lambda
overestimated by 1527x against a measured worst case of 1.263x.

The often-quoted "0.2%" is a large-p figure and should not be quoted globally.
Note also which direction is dangerous: OVERestimating |I| inflates Lambda and
can mask a real anomaly, while underestimating only causes a spurious
escalation. The proved involution bound caps the model from above --
p(1-(1-1/p)^M) < M -- so the overestimate can never exceed the theorem.

VALIDATION (scripts/test_sizelaw.py re-checks all of it):
  * four measured image sizes at the real fat-cell parameters, to 0.15%
  * the Band II pre-registration, 102600/12600/1816/289.9/49.4/8.78/1.61/0.30
    and its even-g fractions, both confirmed by the real run
  * i=9 k=11: run 8, Lambda 2.2e-2 -- ordinary, not Q1
  * the i=7 Z-jump round by round: predicted 8724/643/79/20/6.9/2.0 against
    observed 8783/613/87/23/5/1, every round inside Poisson noise
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

# Escalate when a survivor's expected count falls below this. 0.01 = "we would
# not have expected to see even one such column in this whole regime".
THRESHOLD = 0.01

# Poisson upper-tail alpha for "far more survived than the law allows".
TAIL_ALPHA = 1e-3


# ---------------------------------------------------------------------------
# image size
# ---------------------------------------------------------------------------

def preimage_count(g: int, k: int) -> int:
    """M: pre-images of f_k on its non-vanishing locus, before collisions.

    Exact, from the involution. k odd is fold-free; k even folds 2-to-1.
    """
    return g + 1 if k % 2 else (g + 1) // 2 + 1


def image_bound(g: int, k: int) -> int:
    """PROVED upper bound on |I_{p,k}|, uniform in p and k (Theorem 2).

    Not a model. A measured image larger than this is a bug, not a surprise.
    """
    return preimage_count(g, k)


def image_size(g: int, p: int, k: int) -> float:
    """MODELLED |I_{p,k}|: the proved pre-image count, minus birthday losses."""
    if g <= 0:
        return 1.0
    return p * -math.expm1(preimage_count(g, k) * math.log1p(-1.0 / p))


def survival(g: int, p: int, k: int) -> float:
    """Probability column k survives prime p under the size law."""
    if g <= 0:
        return 1.0
    return -math.expm1(preimage_count(g, k) * math.log1p(-1.0 / p))


# ---------------------------------------------------------------------------
# the trigger
# ---------------------------------------------------------------------------

def run_lambda(k: int, primes) -> float:
    """Lambda: probability of surviving every prime in `primes`."""
    lam = 1.0
    for p in primes:
        lam *= survival(int(p) - int(k), int(p), int(k))
    return lam


def poisson_tail(expected: float, observed: int) -> float:
    """P(X >= observed) for X ~ Poisson(expected). Upper tail, no scipy."""
    if observed <= 0:
        return 1.0
    if expected <= 0:
        return 0.0
    term = math.exp(-expected)
    cdf = term
    for i in range(1, observed):
        term *= expected / i
        cdf += term
    return max(0.0, min(1.0, 1.0 - cdf))


def escalate(expected: float, observed: int, threshold: float = THRESHOLD,
             tail: float = TAIL_ALPHA) -> dict:
    """Judge one round's survivors against what the size law expected.

    Two independent grounds to escalate:
      * we did not expect even one such column (expected < threshold) yet saw one
      * far more survived than expected (Poisson upper tail below `tail`)

    `expected` must be summed over the columns that ACTUALLY entered the round,
    each with the prime it was ACTUALLY tested against. That is what makes the
    verdict regime-correct without a hand-supplied multiplier -- a single rate
    applied to a whole phase is wrong by orders of magnitude, because most Band
    I columns sit at g ~ 9 while the fat cells sit at g ~ 5e5.
    """
    P = poisson_tail(expected, observed)
    unexpected = observed >= 1 and expected < threshold
    excess = P < tail
    return {
        "expected": expected,
        "observed": int(observed),
        "poisson_tail": P,
        "threshold": threshold,
        "escalate": bool(unexpected or excess),
        "reason": (
            "expected count below threshold, yet a column survived"
            if unexpected
            else ("survivor count far above the size law" if excess else "ordinary")
        ),
    }


def assess(k: int, primes, expected: float, observed: int = 1,
           threshold: float = THRESHOLD) -> dict:
    """Report for one column that survived `primes`, in its round's context."""
    primes = [int(x) for x in primes]
    return {
        "k": int(k),
        "k_parity": "odd" if k % 2 else "even",
        "run": len(primes),
        "primes": primes[:16],
        "g_over_p": [round((x - k) / x, 4) for x in primes[:8]],
        "lambda": run_lambda(k, primes),
        **escalate(expected, observed, threshold),
    }


class RoundLedger:
    """Expected vs observed survivors, round by round, during a sweep.

    Cheap: the expectation for a round is a sum over exactly the columns the
    round is about to test, with exactly their primes -- data the sweep has
    in hand when it builds its buckets.
    """

    def __init__(self, phase: str, threshold: float = THRESHOLD):
        self.phase = phase
        self.threshold = threshold
        self.rounds: list[dict] = []

    def expect(self, pairs) -> float:
        """Sum survival over (k, p) pairs about to be tested.

        Sorted first: the caller's order comes from imap_unordered, so an
        unsorted float sum differs in the last digit between runs and makes
        results non-reproducible for no reason.
        """
        terms = sorted((int(k), int(p)) for k, p in pairs)
        return float(math.fsum(survival(p - k, p, k) for k, p in terms))

    def record(self, rnd: int, entered: int, expected: float, observed: int) -> dict:
        row = {
            "round": rnd,
            "entered": int(entered),
            **escalate(expected, observed, self.threshold),
        }
        self.rounds.append(row)
        return row

    def verdict(self) -> dict:
        firing = [r for r in self.rounds if r["escalate"]]
        return {
            "phase": self.phase,
            "threshold": self.threshold,
            "rounds": self.rounds,
            # No rounds means the phase was resumed from a checkpoint and never
            # re-run here, so nothing was judged. Absence of evidence is not
            # "no escalation"; say which it is.
            "evaluated": bool(self.rounds),
            "escalate": bool(firing),
            "escalating_rounds": [r["round"] for r in firing],
            "note": (
                "Escalation is on expected count, not run length. A long run in "
                "a high-g/p regime is ordinary; see scripts/sizelaw.py."
                if self.rounds
                else "not evaluated in this run: the phase was empty, skipped, "
                "or resumed from a checkpoint"
            ),
        }


# ---------------------------------------------------------------------------
# prime chains, per regime
# ---------------------------------------------------------------------------

def live_run(N: int, K: int, k: int, cap: int = 40):
    """Small-k regime: walk primes above k, skipping dead ones, until a kill.

    Dead primes (r(p)=0, a Kummer carry) test nothing and must not count
    toward run length -- that conflation is what made Stage 1's max_r=89
    look alarming.
    """
    import gmpy2

    from witness import image_hit_tablefree, lucas_mod_pure

    p = int(k)
    survived: list[int] = []
    dead = 0
    for _ in range(cap * 40):
        p = int(gmpy2.next_prime(p))
        r = lucas_mod_pure(N, K, p)
        if r == 0:
            dead += 1
            continue
        if image_hit_tablefree(p, r, k) is None:
            return survived, p, dead
        survived.append(p)
        if len(survived) >= cap:
            return survived, None, dead
    return survived, None, dead


def zjump_chain(k: int, ivs, d: int, n: int):
    """Z-jump regime: the live primes column k would be tested against."""
    from bandii_kernel import first_live_after

    out: list[int] = []
    x = k
    for _ in range(n):
        p = first_live_after(x, ivs, d)
        if p is None:
            break
        out.append(p)
        x = p
    return out


# ---------------------------------------------------------------------------
# pre-registration
# ---------------------------------------------------------------------------

def expected_alive(kmin: int, kmax: int, primes) -> list[dict]:
    """Per-pass expected survivors over a contiguous column block.

    This is how the Band II pre-registration was produced. Regenerating it
    is an integrity check on the model, not a replacement for the recorded
    prediction: a pre-registration is only worth something if it was fixed
    before the run.
    """
    import numpy as np

    k = np.arange(kmin, kmax + 1, dtype=np.int64)
    k_odd = (k % 2 == 1)
    alive = np.ones(k.size)
    rows = []
    for j, p in enumerate(primes, start=1):
        g = p - k
        M = np.where(k_odd, g + 1, (g + 1) // 2 + 1)
        alive = alive * -np.expm1(M * np.log1p(-1.0 / int(p)))
        tot = float(alive.sum())
        rows.append(
            {
                "prime_index": j,
                "p": int(p),
                "alive": tot,
                "even_g_frac": (float(alive[k_odd].sum()) / tot) if tot > 0 else None,
                "mean_k": (float((alive * k).sum()) / tot) if tot > 0 else None,
            }
        )
        if tot < 1e-6:
            break
    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="assess one column's live run")
    r.add_argument("--i", type=int, required=True)
    r.add_argument("--k", type=int, required=True)
    r.add_argument("--cap", type=int, default=40)
    r.add_argument("--peers", type=int, default=3000,
                   help="columns in the same regime that could have done the "
                        "same; the honest multiplier for a one-off assessment")

    s = sub.add_parser("scan", help="small-k census ranked by surprise")
    s.add_argument("--i", type=int, required=True)
    s.add_argument("--kmax", type=int, default=3000)
    s.add_argument("--kmin", type=int, default=2)
    s.add_argument("--cap", type=int, default=20)
    s.add_argument("--json_out", type=Path, default=None)

    p = sub.add_parser("predict", help="regenerate a Band II pre-registration")
    p.add_argument("--i", type=int, default=8)
    p.add_argument("--check", action="store_true",
                   help="compare against the recorded pre-registration")

    args = ap.parse_args()
    from bandii_kernel import make_fam

    if args.cmd == "run":
        fam = make_fam(args.i)
        survived, kill, dead = live_run(fam.N, fam.K, args.k, args.cap)
        lam = run_lambda(args.k, survived)
        rep = assess(args.k, survived, expected=args.peers * lam, observed=1)
        rep.update({"i": args.i, "kill_prime": kill, "dead_primes_skipped": dead,
                    "peers": args.peers})
        print(json.dumps(rep, indent=2))
        return 0

    if args.cmd == "scan":
        fam = make_fam(args.i)
        rows = []
        for k in range(args.kmin, args.kmax + 1):
            if k in (fam.K, fam.K + 1):
                continue
            survived, kill, dead = live_run(fam.N, fam.K, k, args.cap)
            if not survived:
                continue
            lam = run_lambda(k, survived)
            rows.append(assess(k, survived,
                               expected=(args.kmax - args.kmin + 1) * lam,
                               observed=1))
        rows.sort(key=lambda x: x["expected"])
        longest = max(rows, key=lambda x: x["run"]) if rows else None
        print(f"  i={args.i}  k={args.kmin}..{args.kmax}  columns with run>=1: "
              f"{len(rows)}", flush=True)
        if longest:
            print(f"  longest run      k={longest['k']} run={longest['run']} "
                  f"expected={longest['expected']:.3g}", flush=True)
            print(f"  most surprising  k={rows[0]['k']} run={rows[0]['run']} "
                  f"expected={rows[0]['expected']:.3g}", flush=True)
            print(f"  -> {'ESCALATE' if rows[0]['escalate'] else 'nothing to escalate'}",
                  flush=True)
        if args.json_out:
            args.json_out.write_text(
                json.dumps({"search": "sizelaw_scan", "i": args.i,
                            "k_range": [args.kmin, args.kmax],
                            "n_with_run": len(rows), "rows": rows[:500]}, indent=2),
                encoding="utf-8",
            )
            print(f"  wrote {args.json_out}", flush=True)
        return 0

    from bandii_kernel import CAP, KMAX, KMIN, PRIMES

    rows = expected_alive(KMIN, KMAX, PRIMES[:CAP])
    print(f"  {'pass':>4} {'p':>9} {'alive':>10} {'even_g':>7} {'mean_k':>10}")
    for row in rows:
        print(f"  {row['prime_index']:>4} {row['p']:>9} {row['alive']:>10.4g} "
              f"{row['even_g_frac']:>7.3f} {row['mean_k']:>10.0f}")
    if args.check:
        from bandii_sweep import PREREGISTER

        bad = 0
        for row in rows:
            pre = PREREGISTER.get(row["prime_index"])
            if not pre:
                continue
            if abs(row["alive"] - pre["alive"]) > 0.02 * max(pre["alive"], 1e-9):
                print(f"    MISMATCH pass {row['prime_index']}: model "
                      f"{row['alive']:.4g} vs recorded {pre['alive']:.4g}", flush=True)
                bad += 1
        print(f"  pre-registration check: {'OK' if not bad else str(bad)+' mismatches'}")
        return 0 if not bad else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
