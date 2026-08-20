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

**Superseded the same day, by its own consequence.** That split is what
motivated Q2, and Q2 then removed the factorial table from the hot path — so
the 90.6% figure now describes a pipeline that no longer exists. Both phases
are scan-dominated today; see Q24 for the re-derived budget, which also
corrected two measurement errors in the profiler itself. The *method* stands
and the ~2x projection caveat stands; only the number is historical.

*Original text:*
There is still no profiling harness. `Claude-Answer.txt` was wrong twice about
the bottleneck (claimed r(p) dominated Stage 3, retracted, then un-retracted),
and this session found the Z-jump is table-dominated rather than scan-dominated
only by measuring. Everything below is guesswork until a profile exists.
Deliverable: per-phase, per-function wall-clock for a real i=7 run, plus a
`--profile` flag.

### Q2. Can `fact_table` be made cheap, or avoided? *(difficulty: medium · **EXTRA**)* — SOLVED 2026-08-20
**Solved by Wilson's theorem.** `k!(p-1-k)! = (-1)^(k+1)` mod p, so `k!` mod p costs `g-1` multiplications rather than `k`. Since the scan reads only O(g) entries anyway, the p-sized table goes away entirely. Measured on real i=8 workloads: **1476x fewer multiplications on the Z-jump** (median g/p ~ 0 there -- the old code built a multi-million entry table to run a fifty-step test) and 1.75x on Band II. End to end on real Z-jump buckets: **11.8x, identical output**; memory 282 MB/worker at i=9 becomes O(g). Implemented as `fact_at` / `fact_window` / `scan_ks_windowed` behind `USE_WINDOWED_SCAN`. Full write-up in [`extra-questions.md`](extra-questions.md).

*The paragraph that used to sit here* said the Z-jump was 90.6% table, that
this was the single biggest remaining cost, and that it was "not yet fixed".
That was written before the Wilson route landed and is exactly what the route
removed. Kept only as the motivation: the profile is what identified the
target, and the target is gone. Q24 re-derived the budget afterwards and
downgraded the Band II half of Q2's claim from 1.75x to a wash.

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

### Q3. Is per-column prime selection the right schedule? *(difficulty: medium · **EXTRA**)* — RESOLVED 2026-08-20
By Q2, not by measurement of the original trade-off. Batching existed to amortise an O(p) table; the table is now O(g_max + spread of k), which makes first-live-prime optimal on BOTH counts (it minimises g, hence both scan and table) and makes batching actively worse (spread-out k widen the window). No change.

*Original text:*
The Z-jump gives every column its own next live prime, which maximises kill rate
per test (smallest g) but also maximises the number of distinct primes, and each
one costs a table. Band II does the opposite: one shared prime for every column,
one table. There is a schedule between them — batch columns onto shared primes,
trading a slightly worse per-test kill rate for far fewer tables. Given a cost
model (Q1) and the size law (already implemented, `sizelaw.survival`), what is
the optimal batching, and how much does it save at i=9/i=10 scale?

### Q4. Re-derive i=8's witnesses. *(difficulty: low · **HARD**)* — DONE 2026-08-20
**Re-derived, and the table is bigger than the one that was lost.** One uniform
pass over \(k=2..k_{\max}=5{,}182{,}637\) — Band II, the Z-jump remnant, the
stragglers and the small-\(k\) band all by the same method rather than four —
returned `clean=True` in 7,508 s on 4 workers and wrote **5,182,634
certificates** to `results/i8_witness.npz`. That is every extra column except
\(\{K,K+1\}\), so i=8 is now covered by a single evidence type instead of a
union of three, and `coverage_ledger.py` reports COMPLETE against a \(k_{\max}\)
recomputed from \(N,K\).

The count is larger than the 4.27M quoted below because the lost runs had
excluded the bands that other artifacts covered; re-deriving in one pass
absorbed them. Verification samples clean, and the run's jsonl carries the
schema header, so this table is reproducible in a way its predecessor was not.

*Original text:*
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

### Q8. Should the caps be adaptive? *(difficulty: medium · **EXTRA**)* — DECLINED 2026-08-20
Solves a problem that has never occurred: no column has reached a cap (i=8 Band II died at pass 8/14, its Z-jump at 7/12, i=7's at 9/12). Caps are already regime-aware and the escalation ledger already distinguishes a surprising survivor from a long one. The one real gain -- certifying a fat-tail survivor rather than flagging it -- is better served by Q14's termination certificate, which proves a killing prime must exist without testing more primes.

*Original text:*
Caps are fixed at 12 (Z-jump) and 14 (Band II). Now that Λ is computed per
round, the natural rule is "keep testing until the expected survivor count drops
below the threshold" — which is tighter at the bottom (the document recommends
~15 live primes at k<10³) and looser above. Does an adaptive cap change any
outcome, and does it cost or save time?

### Q9. Make the ghost census a first-class output. *(difficulty: low-medium · **HARD**)* — DONE 2026-08-20
`scripts/ghost_census.py` reads the witness tables and reports the census with its claim stated exactly: for each recorded (k,p), c = k! C(N,K) is outside (x)_k(F_p), hence outside the intersection over all primes, hence a certified NON-ghost. **6,067,902 values over k = 2..5,182,637**, re-run after Q4 rebuilt i=8's table (it read 884,236 over k = 41..756,136 before that). Verified exactly on a spot-check (k=201, p=211, c a 582-digit number).

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

### Q12. Can the nearby / collide searches be made exhaustive? *(difficulty: medium · **EXTRA**)* — SOLVED 2026-08-20
Yes, and cheaply. The residual has EXACTLY ONE root on [k+d+e, inf) because d/dn log g < 0 there (every term of the d+e-term sum exceeds every term of the d-term sum), so bracket-and-bisect is complete. Verified: exactly one sign change over 216 (k,d,e); **0 differences from the sampled version across 7,200 pairs**, so the recorded nulls were not missing anything -- now provably. Also **2.0-2.6x faster**. `collide` was already exhaustive on its range; its limit is range, not method.

*Original text:*
`nearby_solutions` brackets roots with a geometric sample plus targeted probes
near the attractor `c·k`; the README correctly calls the nulls sampled. What
would a certified-exhaustive version cost — a proven bracketing argument for the
degree-(d+e) residual, or an interval/Sturm method — and is it reachable for
k up to 10⁷?

### Q13. Version control, and the two trees. *(difficulty: trivial · **HARD**)* — DONE 2026-08-20
**Local repository initialised**, deliberately with no remote — linking to
GitHub is a separate decision and is being deferred. `.gitattributes` pins
`*.py` and `*.md` to `eol=lf` before the first commit, because
`pathlib.write_text` on Windows had been silently converting LF to CRLF and
would otherwise have baked that into history. `.gitignore` excludes the live
checkpoints: `results/i8_sweep.jsonl` was committed once by mistake and had to
be untracked — a running job's jsonl changes under the repository and does not
belong in it.

This does not by itself fix the two-tree divergence that motivated the
question; it makes divergence *detectable*, which is the part that was missing.
The live tree `Desktop\Singmaster` remains the authority while its i=9 job runs.

*Original text:*
Still not a git repository. `Desktop\Singmaster` and `Desktop\Claude-Singmaster`
have already diverged silently once (the sub-√N guard existed in one and not the
other for a day). This is the cheapest risk reduction available and it is
overdue.

---

## Tier C — research; would change what the tool can claim

### Q14. Can `(x)_k − c` be intersective? *(difficulty: open, well-posed · **MAX**)* — ANSWERED 2026-08-20
Full answer in [`q14-intersective.md`](q14-intersective.md). **Never for k <= 5** -- and k=5, the first open case, is settled here: a 2+3 split forces 5a^4-10a^2+9 to be square, which is the Pell equation v^2-5u^2=4 with u=a^2-1, so a^2 = F_2j+1; F_n+1 is a square only for n=0,4,6, leaving c in {0,+-210,+-2160}, each killed by an explicit small prime. Open for k >= 6, but constrained by rad(k!)|c, reducibility with all factors >= 2, and (even k) an odd square-relation among the beta_i -- and unobserved across 180 million values of c (k=5..15, |c| <= 1e9).

**The gap is closed operationally regardless.** Irreducibility of (x)_k - k!m over Q implies, by Jordan and Chebotarev, that a killing prime exists with density >= 1/k -- and irreducibility has a one-prime certificate computed by Lucas without ever building m. `scripts/termination_certificate.py` certifies all 29 columns k<=30 of i=8, and all four of i=9's long-run columns k=11,29,40,45. Each such column is upgraded from 'a killing prime was found' to 'a killing prime had to exist'.

*Original text:*
No rational root, but a root modulo every prime. This is **the only logical gap
between the census method and a guarantee that it always terminates**: Frobenius
+ Jordan give a killer prime with density ≥ 1/k *provided* F_k has a derangement
in its Galois group, which fails exactly for intersective F_k. For k=2,
`x²−x−2m` is irreducible unless 8m+1 is square, so never. For even k the
`R(z²)` decomposition blocks the standard construction, and a direct sieve over
|c| ≤ 10⁷ for k ≤ 12 found **zero** ghosts. Odd k ≥ 5 is open. A positive answer
would mean the method can fail; a negative one upgrades every certificate from
"worked" to "provably had to work."

### Q15. A certificate shorter than the search. *(difficulty: research · **EXTRA**)* — ANSWERED (no) 2026-08-20
No, structurally. The claim is universal ('no j < g works'), and the natural short proof -- gcd(x^p - x, (x)_k - k!r) = 1 -- is the polynomial route Q7 measured at O(k^2 log p), which LOSES to the O(g) scan above k ~ 150 and by orders of magnitude at Band II's k ~ 0.76p. What is available instead, and already built: per-certificate independence, parallel verification, sampling with separately-checked coverage, and O(g) time in O(1) memory (~32 ms for the hardest i=7 column).

*Original text:*
Full verification currently costs about what the search cost (≈ 2·Σg modmuls).
Is there a certificate for a *range* of columns rather than one each — a single
prime plus an argument covering an interval of k, or a batched aggregate — that
a referee could check in minutes instead of hours? This is what would turn 4.3M
certificates into a publishable object rather than a large file.

### Q16. Prove the birthday deficit above `g > p^{1/3}`. *(difficulty: research · **EXTRA**)* — DECLINED 2026-08-20
Declined on value, not only difficulty. The payoff is a LOWER bound on |I|, which sharpens pre-registrations but cannot improve a kill guarantee -- a larger image makes killing less likely. The size law already matches measurement to 0.014% at Band II scale (Q23), so a proof would replace an accurate model with a weaker inequality.

*Original text:*
Theorem B covers ~77% of generic Band I columns exactly, and **0%** of the fat
cells and Band II — precisely the columns the pre-registrations depend on. The
counting argument provably runs out at `p^{1/3}`, which is where the phenomenon
starts. Improves the write-up, not the kill guarantees; worth knowing whether
`p^{1/3}` is the true barrier or just this argument's.

### Q17. Effective Chebotarev without the discriminant. *(difficulty: research, likely blocked · **EXTRA**)* — ANSWERED 2026-08-20
Still blocked -- disc(F_k) contains m, 3.1M digits, so GRH bounds are unreachable and the shape does not obviously help. But **Q14 removed the need**: an effective bound would say 'checking to X suffices', whereas the census only needs to know it will stop, and Q14's per-column irreducibility certificate gives exactly that (existence, plus a density >= 1/k rate).

*Original text:*
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

---

## Tier A2 — from the 2026-08-20 first-principles audit

Full findings in [`audit-2026-08-20.md`](audit-2026-08-20.md). The pipeline's
mathematics re-derived clean throughout (image characterisation over 146k
checks, cell geometry over 5,594 primes, Lemma A with no live prime skipped,
four membership paths in agreement). These are what the audit turned up.

### Q19. Should `first_live_after` bisect into the interval list? *(low - **HARD**)* — DONE 2026-08-20
Implemented. Intervals verified DISJOINT at i=8 and i=9 (0 overlaps, 0 containments), which is what makes bisect provably correct -- with overlapping intervals it could land past the container and wrongly return None. Measured **292.5 -> 2.7 us/call**, identical answers on 300 consecutive k plus 400 spread over [100, K).

*Original text:*
It rescans from index 0 on every call. Measured **142x** (292.5 -> 2.1 us/call)
with identical answers: 2.3 h of single-threaded parent time at i=9 becomes
1.2 minutes. That is the stretch during which a running i=9 job looks idle.

### Q20. Should `binom_mod_lucas` take the delta-identity branch? *(low - **HARD**)* — DONE (re-posed) 2026-08-20
**The question was mis-posed and the measurement re-posed it.** delta only wins where it is small: at nextprime_sweep's primes (k~10^6) delta-1 = 196,328 against min(k0,n0-k0) = 126,635, so delta LOSES. The right fix is general and lives one level down, in `binom_mod_prime`: n = -(p-n) mod p gives C(n,k) = (-1)^k C(p-n+k-1, p-n-1), a THIRD lower index p-n-1, and taking min(k, n-k, p-n-1) never loses. Verified on 515,362 exact checks and 203 primes of the full Lucas chain. Band II Lucas: **100+ms -> 0.001ms**; 3890x mean on Band II primes, 940x on fat-cell primes, 1.0x where it does not apply.

The EXTRA sub-question is answered NO: `witness.binom_mod_prime_pure` is deliberately left naive. The verifier's independence rests on full Lucas being a different computation from the delta route, and the negation identity would make them algebraically the same. The verifier stays slow on purpose.

*Original text:*
Measured **33,000-35,000x** at i=9 two-digit primes, same answer.
`r_two_digit_delta` already exists and the sweeps use it; `binom_mod_lucas`
does not, and it is the r(p) route for `modular`, `nextprime_sweep`, and
`witness.lucas_mod_pure`. One sub-question is EXTRA: the verifier uses generic
Lucas as its *independent* route, so making delta primary would cost
independence -- the honest answer may be "keep the verifier slow on purpose".

### Q21. Build a machine-checked coverage ledger. *(low-medium - **HARD**)* — DONE 2026-08-20
`scripts/coverage_ledger.py` states the global claim and checks it: witnessed(i) == [2, k_max] \ {K, K+1}, exactly, with k_max **recomputed from N and K** rather than read back from the file under audit. On first run it found the predicted gap in 6 of 7 members -- 199 missing columns for i=3,4,6,7, 37 for i=2, 1 for i=8. After Q22 all seven report **COMPLETE, 6,067,902 columns, 0 missing, 0 extra**.

*Original text:*
**The one correctness gap found.** `N(m)=6` claims every k in [2,kmax] except
{K,K+1}, but each artifact self-reports only its own slice and the union is
assembled by reading three files with three different evidence types.
`witness.coverage()` checks completeness of the range the sweep *claimed*,
which is not the same statement. Same shape as the false-`clean` bugs: a local
check passing while the global claim goes unchecked.

### Q22. Should the Z-jump start at k=2? *(trivial - **HARD**)* — DONE (re-posed) 2026-08-20
**No -- measurement killed the original plan.** k=2 in the sweep kernel costs a full g~p scan (410 ms at i=8); the O(1) QR shortcut lives in the engine's `column_possible`, not in the kernel. So extending the Z-jump down would be expensive, not free. Instead the witness builder now takes those columns from the engine's modular scan, and `witness.py fill` back-fills existing tables without re-running any sweep -- legitimate precisely because a witness is checkable from (N,K,k,p) however it was found. All 199-per-member added witnesses verify.

*Original text:*
i=8's re-derivation starts at k=3, leaving k=2 to a separate artifact. k=2 is
decided by the QR test in O(1), so including it should be free.

### Q23. Is the size law accurate enough at small primes? *(medium - **EXTRA**)* — ANSWERED 2026-08-20
**Answered: no code change. The bound is the deliverable, and it is now pinned by a test.**

The two error directions are not symmetric. OVERestimating |I| inflates Lambda and can mask a real anomaly; underestimating only causes a spurious escalation. So the question is whether the model can overestimate enough to matter -- and the proved involution bound already caps it from above, since p(1-(1-1/p)^M) < M.

Measured against EXACT image sizes:

| regime | p | model/exact | smallest expected | headroom to threshold |
|---|---:|---:|---:|---:|
| Band II / Z-jump | ~5.4e6 | **1.00014x** | 0.206 | 21x |
| small-k census | ~2e2 | **1.263x** | 15.3 | **1527x** |

**The error and the headroom are anti-correlated, in the safe direction.** The model is loosest exactly where the trigger has three orders of magnitude of room, and accurate to 0.014% exactly where the room is thin (late Band II rounds, where expected counts fall to ~0.2).

Recomputing Lambda from exact image sizes for **all 102 columns** of the i=9 small-k census gives **0 verdict flips**. A flip would need Lambda overestimated by 1527x against a measured worst case of 1.263x -- and the per-prime errors scatter both ways (+-9%) rather than compounding, which is why my ~50%-over-a-run guess in the audit was wrong.

Recorded in `sizelaw.py`'s docstring by regime, pinned by `test_accuracy_bound_by_regime`, and the global "0.2%" claim is corrected wherever it was used as a general figure.

*Original text:*
The 0.2% figure is a large-p measurement. Against exact images the error is
5.6% at p~10^3 and 1.9% at p~10^4, and compounds over a run. Irrelevant for
Band II and the Z-jump; not obviously irrelevant for the small-k census.
Likely answer: no verdict changes, and the deliverable is a documented accuracy
bound by regime rather than a code change.

### Q24. Re-derive the phase budget now that the table is gone. *(low - **HARD**)* — DONE 2026-08-20
`profile_sweep.py` now times the windowed path and reports the old table cost as historical only. Two corrections came out of it: (1) timing the windowed scan on a SPREAD sample charges the whole factorial window to `sample` columns instead of a real chunk, so the profiler now splits fixed per-chunk cost from marginal per-column cost; (2) **Q2's "1.75x on Band II" was a multiplication count and is a wash in wall-clock** (1.01x over chunks of 2k/10k/33k columns, identical output). Z-jump: the pre-Q2 table path would have added ~208 core-h at i=9 that is now gone. The extrapolated totals remain good only to a factor of ~2, so Q11's i=10 figure should be quoted as an order of magnitude, not a number.

*Original text:*
`profile_sweep.py` still times `fact_table`, which Q2 removed from the hot
path, so the "90.6% table" split and the 29.6 h i=9 Z-jump estimate describe
the *old* pipeline. This also feeds Q11, whose "i=10 is 65.7 days" rests on
pre-Q2 numbers.

**No MAX from this audit.** Every gap found is bounded engineering or
documentation accuracy; Q21 is the most structurally important and is still a
scripting job.

---

## Tier C2 — follow-on research

### Q25. Is `(x)_6 - c` ever intersective? *(research · **EXTRA**)* — ANSWERED (no) 2026-08-20
**No.** Q14 settled `k <= 5` and forecast that `k = 6` "will take three
curves", one per split of the sextic. That route works but is the wrong
decomposition. For even `k` the roots are symmetric about `(k-1)/2`, so with
`t = 2x-5` we get `2^6 f_c = g(t^2)` with `g(u) = u^3-35u^2+259u-(225+64c)`,
and since `t^2` ranges over exactly the squares:

> `f_c` has a root mod `p`  <=>  `g` has a root mod `p` that is a SQUARE mod `p`.

A sextic question becomes a **cubic** question, and the case analysis is a
trichotomy on `g` rather than five splits of `f_c`:

- **`g` irreducible** -- Jordan + Chebotarev give a prime with no root at all.
- **`g = (u-beta) x irreducible quadratic`** -- the only infinite case. Both
  `beta` and `Delta = -3beta^2+70beta+189` are forced non-squares (the first
  else `f_c` has a rational root, the second else the quadratic splits), so a
  positive density of primes has both Legendre symbols `-1`; there the
  quadratic has no root and `beta` is a non-residue, so `f_c` has none.
  **Closes an infinite family with no finiteness argument.**
- **`g` splits over Q** -- `e1=35, e2=259` force `-3s^2+140s-1036` to be a
  square, a BOUNDED conic (28 integers), whose only solution is `{1,9,25}`,
  i.e. `c=0`. This is exactly where an odd square-relation (Q14 section 4)
  would have to live, and there is no room for one.

Verified in `scripts/k6_intersective.py`: the equivalence on 1,035 `(c,p)`
pairs, 273 Case-2 candidates all killed by the predicted prime (largest
`p=113`), and a brute-force cross-check over `|c| <= 300000` finding no
unexplained reducible case and no candidate without a killing prime.

**Smallest open degree moves from 6 to 7.** The reduction is general for even
`k` (`k=8` gives a quartic `g` with `e1=84`), and it inverts the intuition:
**even `k` is the easy side**, because odd `k` has no `t^2` reduction at all.
That is why `k=5` needed a Pell equation and `k=6` needed a Legendre symbol.

Operationally this changes nothing -- our `c = k!m` has 3.1M digits and the
census still relies on Q14 section 6's per-column termination certificate.
Full write-up in [`q25-k6-intersective.md`](q25-k6-intersective.md).

### Q26. Is `(x)_7 - c` ever intersective? *(research · **EXTRA**)* — ANSWERED (no, modulo effective Thue) 2026-08-20
**No counterexample exists, and only finitely many can.** Weaker than Q25 and
labelled so: Q25 is a proof, this is a proof modulo an effective Thue
computation not carried out.

Odd `k` gets NO degree reduction. The odd analogue of Q25's trick exists --
`f_c(y)f_c(-y) = T(y^2)` with `T(v) = c^2 - v S(v)^2` -- but `deg T = 7`, not
3. Even `k` maps `2m -> m`; odd `k` maps `2m+1 -> 2m+1`. That is the whole
difference between `k=6` and `k=7`.

Only TWO branches, since `2+2+3` is a `2+5` with a reducible quintic:

- **quadratic factor** -- eliminating the quintic's coefficients leaves one
  condition free of `c`, involving `a` only via `a^2`: the chord curve
  `Phi(A,b)=0` for `P(y1)=P(y2)`. Leading form is an irreducible cubic form,
  so 3 points at infinity and Siegel gives finiteness. Nontrivial points:
  `c = +-17472, +-459648`.
- **cubic factor (3+4)** -- two conditions; the first is linear in `d`, giving
  a degree-5 plane curve with a squarefree leading form (5 points at infinity,
  Siegel again). **The degenerate locus `3a^2-2b-14=0` must be handled
  separately** -- my first pass skipped it and wrongly concluded 3+4 was
  empty. It is where the only nontrivial solution lives: `c = +-896`.

All of `896 = 2^7*7`, `17472 = 2^6*3*7*13`, `459648 = 2^7*3^3*7*19` fail
`rad(7!) = 210`, and specifically **5 divides none of them**. `(x)_7` vanishes
identically on `F_5`, so `f_c == -c != 0` there: no root mod 5.

Cross-validated: brute-force factorisation for EVERY integer `c` in
`[0,400000]` finds exactly the two candidates in that range that the curves
predict. **The gap** is that Siegel gives finiteness but the effective
(Baker/Thue) computation was not run, so the lists are complete only to
`|a| <= 400`. PARI/GP `thue` on the two curves would close it.

Next target is `k = 8`, not `k = 9`: even, so Q25's reduction applies, giving
a quartic `g` with `e1 = 84`, where Case 3 becomes a surface rather than a
bounded conic. Full write-up in [`q26-k7-intersective.md`](q26-k7-intersective.md).

### Q27. Is `(x)_8 - c` ever intersective? *(research · **EXTRA**)* — ANSWERED (no, modulo effective Siegel) 2026-08-20
**No.** Q26 predicted `k=8` would be hard because "Case 3 becomes a surface
rather than a bounded conic". **That prediction was wrong** -- Cases 3 and 5
both stayed bounded, and the difficulty moved into Case 2.

`k=8` is even so Q25's reduction applies, but `g = R(u)-256c` is a QUARTIC, so
the trichotomy becomes five cases. Two need no computation at all:
**g irreducible** (Jordan) and **g = two irreducible quadratics** (both discs
non-square, so a positive density of `p` has both non-residues and `g` has no
root at all -- unconditionally dead, no side condition).

Two are bounded. `e1=84, e2=1974, e3=12916` are FIXED, so two roots summing to
`s` have their product forced, and `disc = s^2-4p(s) ~ -s^2` confines `s` to
`[10,74]`. **Case 5** gives only `{1,9,25,49}`, i.e. `c=0`; **Case 3 proper is
EMPTY**.

That leaves **Case 2**, `g = (u-beta) x irreducible cubic` -- infinite. The
kill needs a prime with `(beta|p) = -1` and Frobenius a 3-cycle on `q`, which
exists unless `beta*disc(q)` is a square. The key: `q` is the divided
difference `(R(u)-R(beta))/(u-beta)`, so **disc(q) is a POLYNOMIAL in beta**,
`disc(q) = -16*P6(beta)`. The dangerous condition is therefore

    y^2 = -beta*P6(beta),   squarefree, degree 7, GENUS 3

whose integral points are exactly `beta in {0,1,9,25,49}` -- the roots of `R`
plus 0. All are perfect squares (excluded) or give `c=0`. **No dangerous beta
exists**, so Chebotarev always applies. Direct check: 197 candidates with
`|beta| <= 4000`, zero dangerous, all killed by `p <= 29`.

**Gap**: Siegel gives finiteness and Baker makes it effective, but the
effective computation was not run (search covers `|beta| <= 20000`). Better
than Q26's gap: the points found are exactly the geometrically forced ones,
not merely "whatever the search turned up".

Next target is **k=10, not k=9**: odd `k=9` has no reduction and seven
partitions to grind, while `k=10` is even with a quintic `g` where the same
divided-difference trick applies. Full write-up in
[`q27-k8-intersective.md`](q27-k8-intersective.md).

### Q28. Is `(x)_10 - c` ever intersective? *(research · **EXTRA**)* — ANSWERED (no, modulo one elliptic computation) 2026-08-20
**No -- and k=10 is CHEAPER than k=8.** Difficulty in this family is **not
monotone in k**; it is governed by whether `g` has an irreducible CUBIC factor.

The kill always needs Frobenius to be a DERANGEMENT of each higher-degree
factor. **n=3 is the only degree whose derangements are all EVEN** (3-cycles,
lying in A_3), which pins `(disc|p)=+1` and can conflict with the character
condition on a rational root -- that conflict is exactly Q27's genus-3 curve.
At k=10 the cofactor is a QUARTIC, 4-cycles are odd derangements, the
character is free, and **no exceptional curve arises at all**.

| n | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| derangements | 1 | 2 | 9 | 44 | 265 |
| odd ones | 1 | **0** | 6 | 20 | 135 |

Cases: **(5)** Jordan, no computation. **(1,4)** unconditional, 66 candidates
all killed by `p<=29`. **>=2 rational roots** -- the constraint curve has
leading form `(b1^5-b2^5)/(b1-b2)`, which has NO real zeros, so the curve is
COMPACT and a bounded search is a PROOF; the bound is 2680 and the complete
search returns exactly the ten pairs from `{1,9,25,49,81}`, i.e. `c=0`. (My
first pass searched to 400 and would have been reported as complete -- the
bound had to be computed.) **(2,3)** reduces to `y^2 = 5a^4+1320a^3+126456a^2
+5102240a+72824400`, squarefree, **GENUS 1**; exactly two candidates survive
`rad(10!)|c` and both die, at `p=11` and `p=13`.

**Gap**: one elliptic integral-point computation. The mildest in the ladder --
easier than Q26's Thue or Q27's genus 3.

**Corrected ordering by difficulty: 6 < 10 < 7 < 8.** Scaling: case count is
`p(k/2)`; where a curve IS needed its genus is `(m-1)(m-2)/2`, m=k/2 -- so
1,3,6,10,15 for k=6..14. Full write-up in [`q28-k10-intersective.md`](q28-k10-intersective.md).

### Q29. Is `(x)_9 - c` ever intersective? *(research · **EXTRA**)* — ANSWERED (no, modulo effective Siegel) 2026-08-20
**No.** Odd `k`, so no degree reduction -- the centred `(x)_9` is ODD in `y`,
which buys only `f_c(-y) = -f_{-c}(y)` (restrict to `c >= 0`).

The eight partitions of 9 into parts >=2 collapse to **three branches**: a
shape containing a 2 is "has a quadratic factor", one containing a 3 but no 2
is "has a cubic factor".

- **A (quadratic factor)** -- the chord curve for `P(y1)=P(y2)`, free of `c`,
  involving `a` only via `a^2`. Constant part factors as
  `-(b+1)(b+4)(b+9)(b+16)`, echoing k=7's `-(w+1)(w+4)(w+9)`. Leading form
  squarefree, 4 points at infinity, Siegel applies. Nontrivial:
  `c = +-176774400` (fails rad, dies p=7) and **`c = +-2630880 = 2^5*3^4*5*7*29`,
  which PASSES rad(9!)** -- the first candidate anywhere in this ladder to do
  so -- and dies at `p=13`.
- **B (cubic, no quadratic)** -- the y^2 condition is QUADRATIC in the cubic's
  constant term, so this needs a resultant (degree 18 x 9). All 42 integer
  points have `c=0`.
- **C (quartic, no quad/cubic)** -- resultant splits into a DEGENERATE locus
  `(2p2-3p3^2+30)^4` (empty; this is the trap that cost me `c=+-896` in Q26)
  and a main component of degree 14 x 24 whose 65 points all have `c=0`.

**A second obstruction mechanism, recorded though moot here.** Killing a 4+5
needs Frobenius to be a derangement of BOTH factors. The sign character blocks
no pair across all 9 transitive subgroups of S4 and all 20 of S5 -- I briefly
concluded 4+5 was unconditionally dead, which was wrong. `C4` and `F20` fuse
over `C4`, and there C4's derangements AVOID the identity coset while F20's lie
entirely INSIDE it (they generate the kernel C5). Incompatible -> kill blocked.
Distinct from the `n=3` mechanism of Q27: that one is a single factor whose
derangements all lie in `A_n`; this is two factors in incompatible cosets of a
shared quotient. **Lesson, now twice: check every transitive group, not the
generic one.**

Cross-checked by factoring `(x)_9 - c` for every multiple of 210 up to
4,200,000: exactly one reducible-all-parts>=2 case, as predicted.

**Gap**: Siegel finiteness on all three branch curves, effective computation
not run. Full write-up in [`q29-k9-intersective.md`](q29-k9-intersective.md).

**Ladder by gap severity: 6 < 10 < 7 < 9 < 8** -- still nothing to do with the
size of `k`.
