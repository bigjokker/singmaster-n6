# Q11 — What is the claim, and is more compute the way to get it?

Answered 2026-08-20, with costs from `scripts/profile_sweep.py` and rates from
the recorded `nearby` / `collide` runs. Short version: **after i=9, stop the
family ladder.** It is the most expensive thing this project can do and the
least informative per hour, and the alternatives are 10x cheaper and aimed at
where a new result could actually be.

---

## 1. What is established

**Proved, effectively, by this machine:**
`N(C(F_{2i+2}F_{2i+3}, F_{2i}F_{2i+3})) = 6` for i = 2..8, with i=9 running.
Every extra column carries a modular impossibility certificate, and since Q4
those certificates are recorded and independently checkable.

**Not proved, and not approached:** Singmaster's conjecture; that 3003 is the
only N=8; that the Lind/Singmaster/Tovey family is the only infinite one. The
README is already correct about this and should stay that way.

**What the 2022 interior theorem does to the picture.** MRSTT prove at most two
left-half solutions in `exp(log^{2/3+eps} n) <= k <= n/2` for sufficiently
large t. The family attains that bound, so morally the interior is settled --
but the threshold is ineffective, and the one explicit exponent puts it past
`log N > 10^2274` against `log N_8 = 16.2` (see `interior-2022.md`). **So the
census's unique contribution is not the conclusion; it is that the conclusion
is effective.** That observation drives everything below: effectiveness is
bought by covering more ground, not by climbing further up one family.

---

## 2. What one more family member buys

Decreasing, and now quantifiably so. Each member is structurally identical to
the last, and since the size law was implemented (`sizelaw.py`) the outcome is
*predicted before the run with zero fitted parameters* -- it reproduced the
Band II pre-registration (102600 / 12600 / 1816 / 289.9 / ...) and tracked the
i=7 Z-jump round by round inside Poisson noise.

A member therefore tests the method, not the mathematics. The method is now
well tested: seven members, ~884k recorded certificates before i=8's
re-derivation, and an escalation trigger that has never fired.

Concretely, i=10 would answer: *"does the size law continue to hold for a
family member whose behaviour it already predicts?"* That is worth something.
It is not worth 65 days.

---

## 3. The cost curve

Columns grow as `phi^(4i)` and per-column work as `phi^(4i)` again, so cost
grows as **phi^8 = 47x per family step**. Measured base: i=9 is 268 core-hours
(Band II 31.7, Z-jump 236.7).

| i | N | columns | core-hours | on 8 workers |
|---|---:|---:|---:|---:|
| 8 | 10,803,704 | 5,182,635 | 6 | 0.7 h |
| 9 | 74,049,690 | 35,522,327 | 268 | **33.5 h** |
| 10 | 507,544,127 | 243,473,668 | 12,609 | **65.7 days** |
| 11 | 3,478,759,200 | 1,668,793,361 | 592,359 | 8.5 years |

Optimisation does not rescue this. The i=9 Z-jump is 90.6% factorial-table
construction, so even solving Q2 outright -- say a 10x cheaper table -- takes
i=10 from 65.7 days to roughly 36. The half-scan (2x, exact) applies mostly to
Band II, which is only ~12% of the cost. **The curve is the problem, not the
constants.**

---

## 4. What the same machine time buys elsewhere

Recorded rates, all consistent at ~1.6e5 (k,pair)/s:

| run | coverage | cost |
|---|---|---:|
| `nearby_k333k-2M_de6` | k=333k..2M, d,e<=6 | 287 s |
| `nearby_k2M-8M_de8` | k=2M..8M, d,e<=8 | 2,227 s |
| `nearby_k8M-1B` | k=8M..1e9, d,e<=6 | 209,016 s (2.4 days) |
| `collide_l20` | 162 pairs, 80k m-steps each | 39 s |

Extrapolating at the measured rate:

- **`nearby` to k=1e9 with d,e<=10** (97 unsettled pairs instead of 33):
  about **7 days**. One tenth of i=10.
- **`collide` to l<=30 with 1e7 m-steps per pair** (426 pairs): about
  **3.6 hours**.

And unlike the family, these regions are where a new result could live. BBW
2017 settled `(d,e)` in {(1,1),(1,2),(2,1)} and the nine elliptic `(k,l)`
pairs; everything else is open. A hit there is a new collision. A hit on the
family would be a new N=8 -- but the family is precisely the case the 2022
theorem morally covers and where the structure makes a surprise least likely.

**The honest caveat on `nearby`:** its nulls are *sampled*, not exhaustive
(Q12). Extending its range without fixing that buys more sampled coverage, not
proof. That makes Q12 the better companion to this recommendation than raw
range extension.

---

## 5. Recommendation

1. **Finish i=9.** It is running, most of the cost is sunk, and it is the
   first member where Band II overtakes small-k as the dominant regime -- so
   it is a genuine test of the size law's two-regime prediction.
2. **Do not start i=10** unless something below changes.
3. **Spend the next block on breadth and on the proof object**, in this order:
   Q12 (make `nearby` exhaustive, or characterise exactly what it misses),
   then extended `nearby` / `collide` range, then Q10/Q9 hardening. Total cost
   of all of it is a small fraction of one family step.
4. **Re-open i=10 only if** (a) i=9's escalation trigger fires, meaning the
   size law failed and the ladder is suddenly informative again; (b) Q2 is
   solved so decisively that the cost drops by more than an order of
   magnitude; or (c) a structural reason emerges to expect member 10 to differ
   -- which, per Q14, would most likely come from the intersectivity question,
   not from more computation.

---

## 6. What would change this answer

This is a judgement about marginal value, so it is worth naming its failure
modes. The recommendation is wrong if:

- the size law is fitting rather than predicting, and i=9 breaks it -- watch
  the escalation ledger, not the run length;
- `nearby`'s sampling is missing solutions systematically rather than
  randomly, in which case extending its range is worthless and Q12 is not
  optional but prerequisite;
- the goal is a *publication about the family specifically*, in which case one
  more member may be worth 65 days for presentation reasons that have nothing
  to do with the mathematics. That is a legitimate reason, but it should be
  named as such rather than dressed up as evidence.
