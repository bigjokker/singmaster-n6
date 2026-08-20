# Open questions

Ranked by what they move, not by how interesting they are. Each carries a
difficulty and a recommended reasoning budget for a follow-up session:

- **HARD** — bounded, known shape, mostly execution. A normal session settles it.
- **EXTRA** — needs real derivation or a design decision with trade-offs, or the
  answer might be "no" and proving that is the work.
- **MAX** — reserved. Only where a positive answer changes what the project *is*.

Written 2026-08-20, after the witness layer, the size-law trigger and the
half-scan landed. Measurements quoted are from this machine.

---

## Correction to the earlier ranking

**The "table-free scan" is not a sweep optimisation.** It was #4 on the previous
list on the theory that killing the O(p) factorial table would speed up the
Z-jump. Measured at p=300007:

| columns per prime | table + numpy | table-free (pure Python) |
|---|---:|---:|
| 1 | 0.015 s | 0.036 s |
| 4 | 0.012 s | 0.166 s |
| 22 | 0.032 s | 0.528 s |
| 100 | 0.100 s | 2.547 s |

Table-free loses 2–26x, even at one column per prime, because a vectorised
numpy scan beats a Python loop by ~50x per element. It is the right tool for
the *verifier* — O(1) memory, checkable on a laptop from five integers, which
is what `witness.py` uses it for — and the wrong tool for the sweep.

The real question underneath it is Q2 below.

---

## Tier A — execution, high value, mostly bounded

### Q1. Where does the time actually go? *(difficulty: low · **HARD**)* — DONE 2026-08-20
`scripts/profile_sweep.py` samples each phase and reports the split between factorial table, image scan, r(p) and prime assignment. The split is decisive and was not guessable: **i=9 Band II is 100% scan; the i=9 Z-jump is 90.6% table** (214 of 237 core-hours). Its distinct-prime count integrates prime density over the live intervals -- a window sample near the range start is ~2x high -- and recovers i=8's recorded 124,830 primes to 0.2%. Wall-clock projection is good to about a factor of 2; the split is the reliable part.

*Original text:*
There is still no profiling harness. `Claude-Answer.txt` was wrong twice about
the bottleneck (claimed r(p) dominated Stage 3, retracted, then un-retracted),
and this session found the Z-jump is table-dominated rather than scan-dominated
only by measuring. Everything below is guesswork until a profile exists.
Deliverable: per-phase, per-function wall-clock for a real i=7 run, plus a
`--profile` flag.

### Q2. Can `fact_table` be made cheap, or avoided? *(difficulty: medium · **EXTRA**)* — MEASURED 2026-08-20
`scripts/profile_sweep.py --i 9` now quantifies it: the i=9 Z-jump is **90.6% factorial-table construction** (214 of 237 core-hours, ~989,000 distinct live primes), against Band II which is 100% scan. So this is the single biggest remaining cost in the pipeline, and the three routes in the original text are still the options. Not yet fixed.

*Original text:*
It is a pure-Python loop over p — **282 MB and several seconds per prime at
i=9**, rebuilt for every distinct prime. At i=7's Z-jump that is 27,622 primes
at ~22 columns each, making the table **29%** of the work; at i=9 it is worse.
Three routes, and the question is which wins:
- **Vectorise the modular prefix product.** No numpy modular `cumprod` exists,
  but a log-depth scan (pairwise multiply + reduce, ~25 passes at i=9) is all
  vectorised. Is it actually faster than the Python loop, and by how much?
- **Build only what is touched.** The half-scan reads `F[0:half]` and
  `F[k:k+half]` — that is `g` entries, and Band II has `g ≈ 0.24p`. The running
  product still costs O(k), but memory could drop **282 MB → 67 MB** per worker
  at i=9. Worth it on its own for an 8-worker box.
- **Amortise it** — see Q3.

### Q3. Is per-column prime selection the right schedule? *(difficulty: medium · **EXTRA**)*
The Z-jump gives every column its own next live prime, which maximises kill rate
per test (smallest g) but also maximises the number of distinct primes, and each
one costs a table. Band II does the opposite: one shared prime for every column,
one table. There is a schedule between them — batch columns onto shared primes,
trading a slightly worse per-test kill rate for far fewer tables. Given a cost
model (Q1) and the size law (already implemented, `sizelaw.survival`), what is
the optimal batching, and how much does it save at i=9/i=10 scale?

### Q4. Re-derive i=8's witnesses. *(difficulty: low · **HARD**)*
i=8 is the headline result and the **only member with no witness table** — its
`bandii_sweep.jsonl` and `zjump.jsonl` were deleted, so its 4.27M certificates
exist nowhere. The `build --i 8` path is implemented and its format adapters are
tested, but it has no data. With the half-scan, Band II should re-run in roughly
half the recorded 1235 s. Question: cheapest complete re-derivation covering
Band II, the Z-jump remnant, the stragglers, and the modular small-k band — and
does the rebuilt table verify?

### Q5. Wire the banked r(p) work. *(difficulty: low · **HARD**)* — DONE 2026-08-20
The delta identity is wired via `bandii_kernel.r_two_digit_delta`, used by `lucas_digits`/`r_of`. It reproduces all nine published Band I r(p) values and matches exact `m mod p` on 400 primes. Measured **7,554-14,123x** on cell-bottom primes and **1.47x aggregate** over the primes the i=8 Z-jump actually uses (the doc said 1.44x). The per-prime ratio is below the doc's 87,568x because the deferred inverse below already sped the baseline 4-5x, so the two overlap. Deferred inverse landed in both `binom_mod_prime` copies and in `r_closed`: **3.1-5.0x**, approaching the doc's 5.5x at large lower index.

*Original text:*
The live tree has `r_two_digit_delta` and `binom_mod_deferred`, written and
spot-checked, both unwired. The δ-identity
`C(n₀,k₀) ≡ (−1)^{k₀} C(k₀+δ−1, δ−1)` with `δ = (α+1)p − N` gives a third index,
so the true cost is `min(k₀, n₀−k₀, δ−1)` — up to **87,568×** on cell-bottom
primes (p=3601237: six multiplications instead of half a million), ~1.44×
aggregate. Deferring the modular inverse is a further 5.5×. Both need wiring and
a regression test; neither changes any result.

### Q19. Cache `r(p)` in `nextprime_sweep`. *(difficulty: trivial · **HARD**) — DONE 2026-08-20*
*(numbered Q19 because it is the next free id, not because it ranks last.)*
`binom_mod_lucas(N, K, p)` sat inside `walk_until_kill`, i.e. inside the
per-column loop, so every `(k, p)` pair recomputed the same value. `r(p)` does
not depend on `k` and consecutive columns share their next prime, so the
recomputation factor is roughly the local prime gap.

Landed as `RCache` — a dict keyed by `p`, pruned of `p <= k` every 1000 columns
so the table stays at the width of the walk window (93 entries resident after
2600 columns, not 371). Rows are byte-identical with the cache on and off;
`scripts/test_sweeps.py` pins that, plus that `r(p)=0` is cached as a value
rather than treated as a miss, plus that a prune drops only dead entries.

Measured on this machine, i=8:

| range | reuse factor | wall-clock speedup |
|---|---:|---:|
| `k=401..10000` (Stage 2) | 33.2× | **7.3×** |
| `k=100001..100200` | 13.8× | **13.7×** |
| `k=999001..1000000` | 15.2× | **15.1×** |

**The estimate was 26×; the measurement is 13–15× where it matters.** Two
things separate the reuse factor from the speedup. At small `k` the column test
is a real share of the work, so removing all the `r(p)` cost cannot buy more
than 7×, whatever the reuse. At Stage-3 `k` the reverse holds — nearly every
column dies on its first prime, so lookups ≈ columns and reuse ≈ the prime gap
(~14 at `10⁶`), and since `r(p)` is then almost the whole cost the speedup
tracks the reuse almost exactly. The document's 26 came from dividing 900k
column-walks by 34.7k *live* primes; the dead primes it left out are the ones
that return 0 before any product, so they inflate the ratio without carrying
any time.

So the retraction-of-the-retraction is right in substance — `r(p)` recomputation
really was essentially all of Stage 3 — and wrong in the constant.

### Q6. Fix the engine's own scan. *(difficulty: low · **HARD**) — DONE 2026-08-20*
**Landed, but the payoff is smaller and narrower than estimated.** The folded,
inverse-free scan is 4-7x on a large-\(g\) membership test (42x on one k=5
case where the fold hits early), verified on 57k exhaustive cases and against
the recorded Stage-2 sweep (400 columns, 0 mismatches in q(k) and r(k)). But
the earlier "20-40x" was the document's estimate, not a measurement, and
CPython's `pow(x,-1,p)` is C-level and cheaper than that assumed.

More importantly the path is **not hot**: `column_possible` uses the cached
image for p <= 4000, and above that the only caller with large g is a deep
single membership test. `nextprime_sweep` has g ~ 15 (it takes the *next*
prime above k), so it measured 0.70x — i.e. unchanged, with r(p) dominating
exactly as Q-beta concluded. The fix is worth keeping (strictly faster, O(1)
memory, makes large --pmax modular runs practical) but it moved no current
workload. For k=3 and k=4 specifically, Q7's closed form is ~1700x and is the
real answer; 4-7x is the best available for k >= 5, where no closed form exists.

*Original text:*
`_column_possible_scan` in `singmaster_intersect.py` still takes a modular
inverse per step and does no folding. Both fixes are now proven in this repo
(the two-modmul recurrence in `witness.image_hit_tablefree`, the involution in
`bandii_kernel.scan_ks_half`). This is the path `extra_reps`, `intersect` and
`collide` hammer. Expect 20–40×.

### Q7. Add the k=3 and k=4 closed-form membership tests. *(difficulty: low-medium · **HARD**)* — DONE 2026-08-20, by a better route
The doc's radical criteria for k=3,4 were verified independently (exact against the brute image at p=101,211,1009, densities 2/3 and 3/8 as predicted) but NOT used. A polynomial gcd is strictly better: `gcd(x^p - x, (x)_k - k! m)` is nontrivial iff the column represents m, needs no radicals and no case analysis, and therefore works for **every k**, not just those where the Galois group is solvable (k in {1,2,3,4,6,8}). Exact on 13,599 (p,k,r) triples for k=2..11. Measured **1239x at k=3**, 847x at k=4, 86x at k=16, 2.9x at k=100, and 0.2x at k=400 -- so it is dispatched only when `k^2 log p < g`, crossover near k=150.

*Original text:*
`Claude-Answer.txt` derives and verifies both (k=3 against 3595 brute-force
cases, k=4 against every c at 11 primes). `column_possible` short-circuits k=2
with the QR test and falls through to an O(g) scan for k=3,4 — roughly **1700×**
per call at p~10⁶. Twenty lines. Also worth recording *why* it stops at k=4:
closed-form criteria exist exactly for k ∈ {1,2,3,4,6,8}, and k=5 is where S_k
arrives.

---

## Tier B — design and scope

### Q8. Should the caps be adaptive? *(difficulty: medium · **EXTRA**)*
Caps are fixed at 12 (Z-jump) and 14 (Band II). Now that Λ is computed per
round, the natural rule is "keep testing until the expected survivor count drops
below the threshold" — which is tighter at the bottom (the document recommends
~15 live primes at k<10³) and looser above. Does an adaptive cap change any
outcome, and does it cost or save time?

### Q9. Make the ghost census a first-class output. *(difficulty: low-medium · **HARD**)* — DONE 2026-08-20
`scripts/ghost_census.py` reads the witness tables and reports the census with its claim stated exactly: for each recorded (k,p), c = k! C(N,K) is outside (x)_k(F_p), hence outside the intersection over all primes, hence a certified NON-ghost. 884,236 values over k = 41..756,136 before i=8's re-derivation. Verified exactly on a spot-check (k=201, p=211, c a 582-digit number).

The output carries its own caveats, because the honest reading is narrower than 'largest test in existence': this is not a targeted ghost hunt (these are the c the Singmaster search happens to produce, not ghost-like candidates), each column stops at its FIRST killing prime so the census records that c fails somewhere rather than how nearly it passed, and 'ghosts found: 0' is the only possible answer -- a surviving column would be an unresolved anomaly, since ghosthood needs failure at every prime.

*Original text:*
Every killed column is a certified non-ghost: a value `c = m·k!` with a prime
witnessing `c ∉ (x)_k(F_p)`. Across i=2..8 that is ~6.07M values of c, and i=9
takes it past 42M — which `Claude-Answer.txt` calls *"the largest test of this
conjecture in existence"* and *"a genuine contribution to a question nobody has
data on."* Right now it is an unrecorded by-product. Should the witness tables
be indexed and reported as a ghost census in their own right?

### Q10. An independent reimplementation of the verifier. *(difficulty: low-medium · **HARD**)* — DONE 2026-08-20
`scripts/verify_independent.py`. No PARI/GP, Sage or Magma on this machine, so sympy plus two from-scratch controls, with the route chosen by cost:

* **brute** -- `math.comb(x,k) % p` over the whole domain. Assumes only the definition; the real control. Used when p is small.
* **sympy-poly** -- factor (x)_k - k! r over GF(p) and look for a linear factor. Different theorem, different algorithm. Cheap only for small k: Cantor-Zassenhaus on degree ~300 is slower than the walk it checks, the same crossover Q7 found.
* **factorial** -- F[k+b] == r F[k] F[b], written fresh, for the large-k/large-p certificates the other two cannot reach.
* **lucas** -- r(p) rebuilt with `sympy.binomial`, on every certificate.

Result: **AGREE** on all 140 i=3 certificates and on samples from i=4..i=7, with full membership coverage (no certificate left at primality-and-r only). Stated limits: sympy is a different library, not a different language or machine, so it catches implementation error rather than a shared misunderstanding; and the `factorial` route uses the same identity as the sweep, so it cross-checks witness.py without independently certifying the sweep.

*Original text:*
`witness.py` shares no code path with the sweep, but both rest on the same
author's Lucas implementation. A second implementation in PARI/GP or Sage,
checking a sample of certificates, closes the last shared-assumption gap and is
the kind of thing a referee asks for. Cheap, and it either finds nothing or
finds something important.

### Q11. What is the claim, and is more compute the way to get it? *(difficulty: judgment · **HARD**)* — ANSWERED 2026-08-20
Full analysis in [`q11-what-is-the-claim.md`](q11-what-is-the-claim.md). **Recommendation: after i=9, stop the family ladder.** Cost grows as phi^8 = 47x per step, so i=10 is **65.7 days on 8 workers** and i=11 is 8.5 years; optimisation cannot rescue a curve like that (even a 10x cheaper factorial table only takes i=10 to ~36 days). Meanwhile the whole `nearby` sweep over k=2M..8M cost 37 minutes, and extending it to k=1e9 with d,e<=10 is about 7 days -- a tenth of one family step, aimed at (d,e) pairs BBW 2017 left open, where a hit would be a new result. The family's own outcome is now predicted before the run with zero fitted parameters, so another member tests the method rather than the mathematics. Named conditions for re-opening i=10 are in the doc.

*Original text:*
i=9 costs 1–3 days; i=10 costs ~47× that. Each member adds one more instance of
a pattern nobody doubts. Meanwhile `nearby` and `collide` are **sampled, not
exhaustive**, and have not been touched in this session's work. An honest
accounting of what each extra family member buys, against what the same machine
time buys elsewhere, before committing to i=10.

### Q12. Can the nearby / collide searches be made exhaustive? *(difficulty: medium · **EXTRA**)*
`nearby_solutions` brackets roots with a geometric sample plus targeted probes
near the attractor `c·k`; the README correctly calls the nulls sampled. What
would a certified-exhaustive version cost — a proven bracketing argument for the
degree-(d+e) residual, or an interval/Sturm method — and is it reachable for
k up to 10⁷?

### Q13. Version control, and the two trees. *(difficulty: trivial · **HARD**)*
Still not a git repository. `Desktop\Singmaster` and `Desktop\Claude-Singmaster`
have already diverged silently once (the sub-√N guard existed in one and not the
other for a day). This is the cheapest risk reduction available and it is
overdue.

---

## Tier C — research; would change what the tool can claim

### Q14. Can `(x)_k − c` be intersective? *(difficulty: open, well-posed · **MAX**)*
No rational root, but a root modulo every prime. This is **the only logical gap
between the census method and a guarantee that it always terminates**: Frobenius
+ Jordan give a killer prime with density ≥ 1/k *provided* F_k has a derangement
in its Galois group, which fails exactly for intersective F_k. For k=2,
`x²−x−2m` is irreducible unless 8m+1 is square, so never. For even k the
`R(z²)` decomposition blocks the standard construction, and a direct sieve over
|c| ≤ 10⁷ for k ≤ 12 found **zero** ghosts. Odd k ≥ 5 is open. A positive answer
would mean the method can fail; a negative one upgrades every certificate from
"worked" to "provably had to work."

### Q15. A certificate shorter than the search. *(difficulty: research · **EXTRA**)*
Full verification currently costs about what the search cost (≈ 2·Σg modmuls).
Is there a certificate for a *range* of columns rather than one each — a single
prime plus an argument covering an interval of k, or a batched aggregate — that
a referee could check in minutes instead of hours? This is what would turn 4.3M
certificates into a publishable object rather than a large file.

### Q16. Prove the birthday deficit above `g > p^{1/3}`. *(difficulty: research · **EXTRA**)*
Theorem B covers ~77% of generic Band I columns exactly, and **0%** of the fat
cells and Band II — precisely the columns the pre-registrations depend on. The
counting argument provably runs out at `p^{1/3}`, which is where the phenomenon
starts. Improves the write-up, not the kill guarantees; worth knowing whether
`p^{1/3}` is the true barrier or just this argument's.

### Q17. Effective Chebotarev without the discriminant. *(difficulty: research, likely blocked · **EXTRA**)*
The qualitative half of "a killer prime exists" is done. The effective half
needs `disc(F_k)`, which contains m (3.1M digits), so even under GRH the
least-prime bound dwarfs anything testable. The question worth asking is not
"can we make it effective" but **"can the dependence on disc be avoided for this
special shape"** — and if not, a clean statement of why, so it stops being
re-attempted.

### Q18. Is 3003 alone, and is the family the only infinite one? *(difficulty: open — do not budget)*
The actual prize, and neither MAX nor any other effort level will settle these
in a session. They are listed so they are not confused with the tractable
questions above. If they are to be attacked, it should be through a named
sub-question with a decidable answer — e.g. Q14 — not head-on.

---

## Notes on effort levels

Most of Tier A is **HARD**: the shape is known and the work is execution plus
regression tests. Q2, Q3, Q8, Q12, Q15, Q16, Q17 are **EXTRA** because each
involves either a real derivation or a design trade-off where the honest answer
may be "this does not help," and establishing *that* is the deliverable.

Only **Q14** is marked **MAX**, because it is the one question whose answer
changes the standing of every certificate the project has produced.

A pattern worth carrying into any of these: four times now — the r(p)
bottleneck, the table-free scan, the parity-fold test threshold, and Q19's 26×
that measured 13–15× — a plausible claim survived reasoning and died on
measurement. Q19 adds a sharper version of the rule: a *reuse* factor is not a
*speedup*. Measure the clock, not the counter.
