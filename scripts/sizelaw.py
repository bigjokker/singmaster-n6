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


def preimage_count_vec(g, k):
    """Vectorised `preimage_count`. Same formula, numpy arrays in and out."""
    import numpy as np
    g = np.asarray(g, dtype=np.int64)
    k = np.asarray(k, dtype=np.int64)
    return np.where(k % 2 == 1, g + 1, (g + 1) // 2 + 1)


def survival_vec(g, p, k):
    """Vectorised `survival`, for callers sweeping millions of columns.

    Exists so the profiler does not carry a second copy of the size law --
    a re-implementation that drifted would make the model and the escalation
    trigger silently disagree. `test_sizelaw.py` pins this elementwise
    against the scalar form.
    """
    import numpy as np
    g = np.asarray(g, dtype=np.int64)
    M = preimage_count_vec(g, k)
    out = -np.expm1(M * np.log1p(-1.0 / np.asarray(p, dtype=np.float64)))
    return np.where(g <= 0, 1.0, out)


# ---------------------------------------------------------------------------
# the trigger
# ---------------------------------------------------------------------------

def run_lambda(k: int, primes) -> float:
    """Lambda: probability of surviving every prime in `primes`."""
    lam = 1.0
    for p in primes:
        lam *= survival(int(p) - int(k), int(p), int(k))
    return lam


def _gammainc_lower_reg(a: float, x: float) -> float:
    """Regularized lower incomplete gamma P(a,x). Series / continued fraction.

    The whole point is the prefactor exp(-x + a ln x - lgamma(a)), formed in
    LOG space. Those three terms very nearly cancel for a ~ x, which is exactly
    the regime the escalation trigger lives in -- computing exp(-x) on its own
    first, as the old summation did, throws the answer away before the
    cancellation can happen.
    """
    if x <= 0.0 or a <= 0.0:
        return 0.0
    lg = -x + a * math.log(x) - math.lgamma(a)
    if lg < -745.0:                      # prefactor is zero to double precision
        return 0.0 if x < a else 1.0
    pref = math.exp(lg)
    if x < a + 1.0:
        ap, term, total = a, 1.0 / a, 1.0 / a
        for _ in range(100_000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        return min(1.0, max(0.0, total * pref))
    tiny = 1e-300                        # modified Lentz for Q(a,x)
    b, c, d = x + 1.0 - a, 1.0 / tiny, 1.0 / max(x + 1.0 - a, tiny)
    h = d
    for i in range(1, 100_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return min(1.0, max(0.0, 1.0 - pref * h))


def poisson_tail(expected: float, observed: int) -> float:
    """P(X >= observed) for X ~ Poisson(expected). Upper tail, no scipy.

    P(X >= k) = P(k, lambda), the regularized LOWER incomplete gamma.

    The previous form summed the pmf from i=0 and started with
    term = math.exp(-expected). That underflows to 0.0 for expected > ~745, so
    the cdf accumulated nothing and the tail returned 1.0 -- "ordinary" -- for
    EVERY round with expected above that. Measured before the fix:

        escalate(expected=102,563, observed=500,000) -> escalate: False

    A five-fold excess of survivors read as ordinary. That covered Band II
    passes 1-4 and Z-jump rounds 1-2 at every member, i.e. essentially every
    column the project has ever swept: the anomaly detector was silently
    inoperative exactly where the counts are large. It was also O(observed),
    which is 100k+ iterations at Band II scale.
    """
    if observed <= 0:
        return 1.0
    if expected <= 0:
        return 0.0
    return _gammainc_lower_reg(float(observed), float(expected))


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
        self.survivors: list[dict] = []
        self.peers: int | None = None

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

    def cap_survivors(self, entries, peers: int) -> list[dict]:
        """Judge every column still alive at the phase's cap. (D4)

        The per-round rule is memoryless: with ONE column entering a round,
        E_r is that column's own next-prime survival -- 0.04..0.21 in Band
        II, ~0.63 at small k -- never below THRESHOLD, and P(X >= 1 | E_r) is
        never below TAIL_ALPHA. So a lone column surviving every pass to the
        cap reads 'ordinary' in every round, although its run has Lambda
        ~ 3e-10 (i=8 Band II k=4126649 over 14 passes) and the whole phase
        expected ~3e-4 such columns. That is the motivating anomaly, and the
        per-round ledger cannot see it.

        This judges each cap survivor with the SAME machinery the small-k
        census already uses: Lambda over the primes the column ACTUALLY
        faced, weighed against how many columns had a chance to do it --
        expected = peers x Lambda, observed = 1, the same THRESHOLD and
        TAIL_ALPHA. No new rule. Run length never enters: i=9's four Z-jump
        columns alive at round 12/12 have Lambda ~ 4e-3 over 28.3M peers,
        expected ~1e5, and stay ordinary, exactly as their round did.

        `entries` are (k, primes_faced) pairs; `peers` is the number of
        columns that entered the phase. The per-round records are untouched;
        verdict() fires on either.
        """
        self.peers = int(peers)
        rows = []
        for k, primes in entries:
            primes = [int(p) for p in primes]
            lam = run_lambda(int(k), primes)
            rows.append(assess(int(k), primes, expected=self.peers * lam,
                               observed=1, threshold=self.threshold))
        rows.sort(key=lambda r: r["k"])
        self.survivors.extend(rows)
        return rows

    def cap_unassessable(self, k: int, reason: str) -> dict:
        """A cap survivor whose chain cannot be reconstructed (its recorded
        last prime is not on the live ladder). The record and the ladder
        disagree, so no Lambda exists to judge; that is not 'ordinary', it is
        a flag -- record it as an escalation with the reason, and keep going
        so the sweep's json is still written."""
        row = {"k": int(k), "run": None, "primes": [], "lambda": None,
               "expected": None, "observed": 1, "escalate": True,
               "reason": f"unassessable: {reason}"}
        self.survivors.append(row)
        return row

    def verdict(self) -> dict:
        firing = [r for r in self.rounds if r["escalate"]]
        firing_caps = [r for r in self.survivors if r["escalate"]]
        return {
            "phase": self.phase,
            "threshold": self.threshold,
            "rounds": self.rounds,
            # No rounds means the phase was resumed from a checkpoint and never
            # re-run here, so nothing was judged. Absence of evidence is not
            # "no escalation"; say which it is.
            "evaluated": bool(self.rounds),
            "escalate": bool(firing or firing_caps),
            "escalating_rounds": [r["round"] for r in firing],
            # Columns alive at the cap, each judged by peers x Lambda over the
            # primes it actually faced. Empty when nothing reached the cap.
            "cap_peers": self.peers,
            "n_cap_survivors": len(self.survivors),
            "cap_survivors": self.survivors[:100],
            "escalating_survivors": [r["k"] for r in firing_caps],
            "note": (
                "Escalation is on expected count, not run length: per round "
                "over the columns entering it, and per column alive at the cap "
                "(peers x Lambda over the primes it faced). A long run in a "
                "high-g/p regime is ordinary; see scripts/sizelaw.py."
                if self.rounds
                else "not evaluated in this run: the phase was empty, skipped, "
                "or resumed from a checkpoint"
            ),
        }


# ---------------------------------------------------------------------------
# prime chains, per regime
# ---------------------------------------------------------------------------

class _Exhausted:
    """Sentinel: the walk's prime budget ran out before `cap` live primes were
    seen and before a kill. Distinct from None (= the column reached the cap)
    and from a prime (= killed). Falsy, so `if kill:` still means "killed"."""

    def __bool__(self) -> bool:
        return False

    def __repr__(self) -> str:
        return "EXHAUSTED"


EXHAUSTED = _Exhausted()


def live_run(N: int, K: int, k: int, cap: int = 40, budget: int | None = None):
    """Small-k regime: walk primes above k, skipping dead ones, until a kill.

    Dead primes (r(p)=0, a Kummer carry) test nothing and must not count
    toward run length -- that conflation is what made Stage 1's max_r=89
    look alarming.

    Returns (survived, kill, dead) with THREE distinct second slots:
      a prime   -- the column was killed there, after `survived`;
      None      -- the column survived `cap` live primes (reached the cap);
      EXHAUSTED -- the walk's budget (`budget` primes, default cap*40) ran
                   out first. `survived` is then a PARTIAL run, possibly
                   empty: the column was not tested to the cap and must not
                   be assessed as if it had been. Under a forced-zero slab
                   the first live prime can be tens of thousands of primes
                   away (i=8 k=2227205: 2700967), so a small budget sees
                   only dead primes; use live_ladder for Band I columns.
    """
    import gmpy2

    from witness import image_hit_tablefree, lucas_mod_pure

    p = int(k)
    survived: list[int] = []
    dead = 0
    n_budget = int(budget) if budget is not None else cap * 40
    for _ in range(n_budget):
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
    return survived, EXHAUSTED, dead


def live_ladder(fam, k: int, cap: int = 40):
    """Band I regime: walk the member's LIVE-PRIME LADDER above k, testing
    each rung, until a kill or `cap` rungs. The ladder skips forced-zero
    slabs analytically (cells / live_intervals), so this reaches the first
    live prime of a fat-cell column in one step where a prime-by-prime walk
    sees only dead primes. Same return shape as live_run; EXHAUSTED here
    means the ladder ran out (no live prime left below d), which is the
    sweep's "no live prime" anomaly, not a budget.
    """
    from bandii_kernel import cells, first_live_after, live_intervals
    from witness import image_hit_tablefree, lucas_mod_pure

    ivs = live_intervals(cells(fam), fam)
    survived: list[int] = []
    x = int(k)
    for _ in range(cap):
        p = first_live_after(x, ivs, fam.D)
        if p is None:
            return survived, EXHAUSTED, 0
        r = lucas_mod_pure(fam.N, fam.K, p)
        if r == 0:
            # the ladder promised a live prime; r=0 means ladder and family
            # disagree -- refuse rather than count it either way
            raise RuntimeError(f"ladder prime {p} has r(p)=0 for k={k}")
        if image_hit_tablefree(p, r, k) is None:
            return survived, p, 0
        survived.append(p)
        x = p
    return survived, None, 0


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


def chain_to(k: int, last_p: int, ivs, d: int, limit: int = 256) -> list[int]:
    """The live primes column k faced, from the first after k through last_p.

    A survivor record stores k and g = p_last - k, so the chain a cap
    survivor actually walked is the live ladder in (k, k+g]. Reconstruct it
    from the ladder rather than guessing a count from the cap. A last_p that
    is not ON the ladder means the record and the ladder disagree; refuse
    rather than assess a chain the column never walked.
    """
    from bandii_kernel import first_live_after

    out: list[int] = []
    x = int(k)
    last_p = int(last_p)
    while len(out) < limit:
        p = first_live_after(x, ivs, d)
        if p is None or p > last_p:
            break
        out.append(p)
        if p == last_p:
            return out
        x = p
    raise ValueError(
        f"k={k}: last prime {last_p} is not on the live ladder after k "
        f"(walked {len(out)} primes, reached {out[-1] if out else None})"
    )


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


def predict_member(i: int, cap: int | None = None):
    """The Band II pre-registration curve for member i, from the model alone.

    Window and primes are derived from (N, K) exactly as family_sweep derives
    them: columns [K+2, k_max], the first CAP_BII live primes in (N/2, d]
    above k_max. Nothing is read from a run record. Returns (rows, meta).

    This is what `predict --i N` prints. It used to import i=8's module
    constants regardless of N, so `predict --i 9` printed i=8's curve and
    `--check` scored it against i=8's table.
    """
    from bandii_kernel import first_primes_above, kmax_of, make_fam
    from family_sweep import CAP_BII

    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    n = int(cap) if cap else CAP_BII
    primes = first_primes_above(fam.N2, fam.D, kmax, n=max(n, 16))[:n]
    kmin = fam.K + 2
    rows = expected_alive(kmin, kmax, primes)
    return rows, {"i": int(i), "N": fam.N, "K": fam.K, "kmin": kmin,
                  "kmax": kmax, "primes": primes}


def check_member(i: int, rows, meta, tol: float = 0.02) -> dict:
    """Is the regenerated curve the one this member is on record with?

    i=8: against the pre-registration that was fixed BEFORE its Band II run
    (scripts/bandii_sweep.py PREREGISTER, also results/bandii_sweep.json).
    Any other member: against its run record results/i{i}_sweep.json --
    the window k_bii and the prime list primes_bii must be the ones derived
    here, and the model over THAT window must equal `rows`. Where the record
    carries per-pass observed survivors they are reported beside the model,
    as information; the recorded run is an outcome, not a prediction, and is
    not scored. A member is never scored against another member's table.
    """
    out = {"i": int(i), "against": None, "ok": False, "mismatches": 0,
           "notes": [], "preregister_compared": False}
    if int(i) == 8:
        from bandii_sweep import PREREGISTER

        out["against"] = "preregister"
        out["preregister_compared"] = True
        for row in rows:
            pre = PREREGISTER.get(row["prime_index"])
            if not pre:
                continue
            if abs(row["alive"] - pre["alive"]) > tol * max(pre["alive"], 1e-9):
                out["mismatches"] += 1
                out["notes"].append(
                    f"pass {row['prime_index']}: model {row['alive']:.4g} vs "
                    f"recorded {pre['alive']:.4g}")
        out["ok"] = out["mismatches"] == 0
        return out

    rec_path = ROOT / "results" / f"i{i}_sweep.json"
    if not rec_path.exists():
        out["against"] = None
        out["notes"].append(f"no run record {rec_path.name}; nothing to check "
                            f"against -- the model is printed, not scored")
        out["ok"] = False
        return out
    rec = json.loads(rec_path.read_text(encoding="utf-8"))
    out["against"] = "sweep_record"
    k_bii = [int(x) for x in rec.get("k_bii", [])]
    primes_rec = [int(x) for x in rec.get("primes_bii", [])]
    if k_bii != [meta["kmin"], meta["kmax"]]:
        out["mismatches"] += 1
        out["notes"].append(f"window: derived [{meta['kmin']}, {meta['kmax']}] vs "
                            f"recorded {k_bii}")
    if primes_rec[:len(meta["primes"])] != meta["primes"]:
        out["mismatches"] += 1
        out["notes"].append("primes_bii in the record differ from the derived ladder")
    if not out["mismatches"]:
        ref = expected_alive(k_bii[0], k_bii[1], primes_rec)
        for a, b in zip(rows, ref):
            if abs(a["alive"] - b["alive"]) > 1e-9 * max(1.0, b["alive"]):
                out["mismatches"] += 1
                out["notes"].append(f"pass {a['prime_index']}: model over the derived "
                                    f"window {a['alive']:.6g} != over the recorded "
                                    f"window {b['alive']:.6g}")
    observed = {int(r["prime_index"]): int(r["n"])
                for r in (rec.get("phases") or {}).get("bandii") or []}
    out["observed"] = observed
    if not observed:
        out["notes"].append("the record carries no per-pass Band II survivors "
                            "(resumed leg); model printed without an observed column")
    out["ok"] = out["mismatches"] == 0
    return out


def run_report(i: int, k: int, cap: int, peers: int, ladder: bool = False,
               budget: int | None = None) -> dict:
    """What `run` prints. Assesses ONLY a column that was actually tested to
    a kill or to the cap; an exhausted walk is reported as untested, never as
    'ordinary'."""
    from bandii_kernel import make_fam

    fam = make_fam(i)
    if ladder:
        survived, kill, dead = live_ladder(fam, k, cap)
    else:
        survived, kill, dead = live_run(fam.N, fam.K, k, cap, budget=budget)
    base = {"i": int(i), "k": int(k), "cap": int(cap), "peers": int(peers),
            "walk": "ladder" if ladder else "next-prime",
            "dead_primes_skipped": dead}
    if kill is EXHAUSTED:
        n_budget = int(budget) if budget is not None else cap * 40
        why = ("no live prime left on the ladder below d"
               if ladder else
               f"walked {n_budget} primes above k and met {len(survived)} live "
               f"prime(s) and {dead} dead; the first live prime of a Band I "
               f"column can sit tens of thousands of primes up. Re-run with "
               f"--ladder (walks the member's live ladder) or a larger --budget.")
        return {**base, "tested": False, "run": len(survived), "primes": survived,
                "lambda": None, "expected": None, "escalate": None,
                "reason": f"NOT TESTED: {why}", "kill_prime": None}
    lam = run_lambda(k, survived)
    rep = assess(k, survived, expected=peers * lam, observed=1)
    rep.update(base)
    rep.update({"tested": True, "kill_prime": kill,
                "reached_cap": kill is None})
    return rep


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
    r.add_argument("--ladder", action="store_true",
                   help="walk the member's live-prime ladder (skips forced-zero "
                        "slabs analytically) instead of every prime above k; "
                        "what a Band I / fat-cell column needs")
    r.add_argument("--budget", type=int, default=None,
                   help="primes to walk before giving up (default cap*40)")

    s = sub.add_parser("scan", help="small-k census ranked by surprise")
    s.add_argument("--i", type=int, required=True)
    s.add_argument("--kmax", type=int, default=3000)
    s.add_argument("--kmin", type=int, default=2)
    s.add_argument("--cap", type=int, default=20)
    s.add_argument("--json_out", type=Path, default=None)

    p = sub.add_parser("predict", help="regenerate a member's Band II curve "
                                        "from the model alone")
    p.add_argument("--i", type=int, default=8)
    p.add_argument("--check", action="store_true",
                   help="i=8: compare against its pre-registration; other i: "
                        "against the window and primes in results/i{i}_sweep.json")

    args = ap.parse_args()
    from bandii_kernel import make_fam

    if args.cmd == "run":
        rep = run_report(args.i, args.k, args.cap, args.peers,
                         ladder=args.ladder, budget=args.budget)
        print(json.dumps(rep, indent=2))
        return 0 if rep["tested"] else 2

    if args.cmd == "scan":
        fam = make_fam(args.i)
        rows = []
        untested: list[int] = []
        for k in range(args.kmin, args.kmax + 1):
            if k in (fam.K, fam.K + 1):
                continue
            survived, kill, dead = live_run(fam.N, fam.K, k, args.cap)
            if kill is EXHAUSTED:
                # not tested to a kill or to the cap: do not rank it, say so
                untested.append(k)
                continue
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
        if untested:
            print(f"  NOT TESTED (walk budget exhausted before a kill or the cap): "
                  f"{len(untested)} column(s), e.g. {untested[:8]} -- excluded "
                  f"from the ranking, not called ordinary", flush=True)
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
                            "n_with_run": len(rows), "n_untested": len(untested),
                            "untested": untested[:500], "rows": rows[:500]}, indent=2),
                encoding="utf-8",
            )
            print(f"  wrote {args.json_out}", flush=True)
        return 0

    rows, meta = predict_member(args.i)
    print(f"  i={args.i}  Band II columns [{meta['kmin']}, {meta['kmax']}]  "
          f"{len(meta['primes'])} live primes from p1={meta['primes'][0]}")
    print(f"  {'pass':>4} {'p':>9} {'alive':>10} {'even_g':>7} {'mean_k':>10}")
    for row in rows:
        print(f"  {row['prime_index']:>4} {row['p']:>9} {row['alive']:>10.4g} "
              f"{row['even_g_frac']:>7.3f} {row['mean_k']:>10.0f}")
    if args.check:
        chk = check_member(args.i, rows, meta)
        for note in chk["notes"]:
            print(f"    {note}", flush=True)
        obs = chk.get("observed") or {}
        if obs:
            print(f"  {'pass':>4} {'model':>10} {'observed':>9}   (observed is the "
                  f"recorded run: an outcome, not scored)")
            for row in rows:
                if row["prime_index"] in obs:
                    print(f"  {row['prime_index']:>4} {row['alive']:>10.4g} "
                          f"{obs[row['prime_index']]:>9}")
        what = {"preregister": "pre-registration (fixed before the run)",
                "sweep_record": f"results/i{args.i}_sweep.json window and primes",
                None: "nothing (no record for this member)"}[chk["against"]]
        print(f"  check against {what}: "
              f"{'OK' if chk['ok'] else str(chk['mismatches']) + ' mismatches'}")
        return 0 if chk["ok"] else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
