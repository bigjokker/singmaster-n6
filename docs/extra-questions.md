# The EXTRA questions — Q2, Q3, Q8, Q12, Q15, Q16, Q17, and Q18

Answered 2026-08-20. EXTRA was defined as "the honest answer may be *this does
not help*, and establishing that is the deliverable". Two of these earned an
implementation; four did not; one resolved itself as a consequence of another.
Reasoning first, code only where it survived the reasoning.

| | question | verdict |
|---|---|---|
| Q2 | make `fact_table` cheap or avoid it | **implemented — 11.8x measured** |
| Q3 | is per-column prime selection right | **resolved by Q2: yes, don't change it** |
| Q8 | adaptive caps | **declined** — solves a problem that has never occurred |
| Q12 | exhaustive nearby / collide | **implemented — provably exhaustive, and 2.4x faster** |
| Q15 | a certificate shorter than the search | **no** — with the argument |
| Q16 | prove the birthday deficit above p^(1/3) | **declined** — blocked, and low value even if solved |
| Q17 | effective Chebotarev without the discriminant | **blocked, and Q14 removed the need** |
| Q18 | is 3003 alone / is the family the only infinite one | not attackable head-on; see below |

---

## Q2 — `fact_table`, solved by Wilson's theorem

**The finding.** The scan reads only `F[0:half]` and `F[k:k+half]` with
`half = ceil(g/2)` — that is `O(g)` entries. `fact_table` built all `p` of
them. And the reason it seemed unavoidable is that reaching `F[k]` looks like
`k` sequential multiplications.

It is not. Wilson's theorem splits `(p-1)!` at `k`:

    k! * (p-1-k)! = (-1)^(k+1)  (mod p),   so   k! = (-1)^(k+1) / (g-1)!

with `g = p-k`. **So `k!` mod `p` costs `g-1` multiplications, not `k`** —
and in the Band II and Z-jump regimes `k` is close to `p`, so `g` is the small
parameter. Verified against direct computation at `p` up to `10^6`.

**Why it matters so much.** Measured on real i=8 workloads, the fraction of
the table actually touched:

| phase | median `g/p` | multiplications: table | windowed | ratio |
|---|---:|---:|---:|---:|
| Band II | 0.138 | 5,401,853 | 3,078,401 | 1.75x |
| Z-jump (sampled) | **0.000** | 4,780,976,329 | 3,239,050 | **1476x** |

The Z-jump number is the striking one: the first live prime usually sits just
above `k`, so `g` is tiny, and the old code built a multi-million-entry
factorial table to run a fifty-step test.

**Implemented** as `fact_at`, `fact_window`, `scan_ks_windowed` in
`bandii_kernel.py`, wired into both `scan_columns` and `scan_columns_general`
behind `USE_WINDOWED_SCAN` for A/B. Output is identical record for record
(175 random `(p,r,ks)` cases; the zjump preflight benchmark `k=268733`,
`p=270097`, `b=589` reproduces exactly; every witness table sha256 unchanged).

**Measured end to end: 11.8x on real i=8 Z-jump buckets**, identical output.

**Correction (Q24, same day).** The 1.75x for Band II in the table above
counts MULTIPLICATIONS, and it does not survive a wall-clock A/B: on
realistic Band II chunks of 2,000 / 10,000 / 33,000 columns the two paths
measure 1.02x, 0.98x and 1.01x -- a wash, with identical output. Band II
amortises a single table over an entire chunk and the numpy scan dominates
it, so there is nothing there for the windowed build to win back. **The
11.8x is entirely a Z-jump effect**, and that is the honest scope of Q2.

One safety note: `scan_columns` used `r_checked` to cross-check `r` against the
factorial table. With no table that check is gone, so it now cross-checks two
independent table-free routes instead (`r_closed` against
`r_two_digit_delta`), preserving the property.

---

## Q3 — per-column prime selection: resolved by Q2, don't change it

The question was whether to batch columns onto shared primes to amortise the
factorial table, trading a worse per-test kill rate for fewer tables.

**Q2 dissolves it.** The amortisation argument existed only because the table
cost `O(p)` regardless of how many columns used it. Now it costs
`O(g_max + spread of k)`. That flips the trade-off twice over:

- picking the *first* live prime above `k` minimises `g`, which now minimises
  the table cost as well as the scan cost — the two objectives agree;
- batching columns with spread-out `k` onto one prime *widens* the window
  `[min k, max k + half]`, so batching is now actively worse.

So the current schedule is optimal on both counts, more clearly than before.
No change. This is the cleanest kind of EXTRA outcome: the work went into
establishing that the intervention is unnecessary.

---

## Q8 — adaptive caps: declined

The proposal: instead of a fixed cap (14 Band II, 12/15 Z-jump), keep testing
until the size law says a survivor would be surprising.

**It solves a problem that has never occurred.** i=8 Band II died at pass 8 of
14; its Z-jump at round 7 of 12; i=7's Z-jump at round 9. No column has ever
reached a cap. The caps are already regime-aware (15 below `k=1000`, where
`g/p` runs to 0.94 and long runs are ordinary), and the escalation ledger
already reports whether a survivor is surprising rather than merely long.

The one real gain would be turning a hypothetical fat-tail survivor from
"unresolved" into "certified" instead of flagging it. But that gain is now
better served by Q14's termination certificate, which proves a killing prime
*must* exist for the column without testing any more primes at all. Adaptive
caps would be strictly weaker and add a control-flow branch to the phase loop.

**Revisit if** a column ever does reach a cap — at which point the certificate
is the first thing to run, not more rounds.

---

## Q12 — exhaustive `nearby`: implemented, and it was cheap

The nulls were "sampled, not exhaustive" because `nearby_solutions` bracketed
roots with a geometric walk plus probes near the attractor `c*k`.

**That caveat can be removed, because the residual has exactly one root.**
With `f(n) = (k+e)_e (n)_d - (n-k)_{d+e}` on `n >= k+d+e`, put
`g(n) = (k+e)_e (n)_d / (n-k)_{d+e}`. Then

    d/dn log g = sum_{i<d} 1/(n-i) - sum_{j<d+e} 1/(n-k-j),

and for `n > k+d+e` every term of the second sum exceeds every term of the
first (`n-k-j < n-i` whenever `k > i-j`, and `k >= d` here) while the second
sum has `d+e > d` terms. So `log g` is strictly decreasing; `g` runs from a
large value at `n = k+d+e` down to `0`, crossing `1` exactly once. Verified:
over 216 `(k,d,e)` combinations with `d,e <= 6` and `k` to `10^5`, `f` changes
sign **exactly once** on the region, never twice.

So a bracket-and-bisect is complete. `nearby_solutions_exhaustive` implements
it and `run_nearby` now uses it.

- Recovers all five family members tested (`i=1..5`, the `(1,1)` case).
- **0 differences from the sampled version across 7,200 `(k,d,e)` pairs** — so
  the recorded nulls were not missing anything, and that is now provable
  rather than hoped.
- **2.0–2.6x faster** than the sampled walk, since a bisection beats ~120
  probe evaluations.

This matters because Q11 recommends spending the next machine block on
`nearby` breadth rather than another family member, and that recommendation
was explicitly conditional on the sampling being fixed first.

`collide` is untouched: it enumerates `m` directly over an interval with an
exact incremental `C(m,l)`, so it is already exhaustive on the range it
covers. Its limitation is range, not method.

---

## Q15 — a certificate shorter than the search: no

Full verification currently costs about what the search cost (`~2 * sum g`).
The question was whether a certificate could cover a *range* of columns, or
otherwise be checked in materially less time.

**The answer is no, and the reason is structural.** The claim being certified
is universal — "no `j < g` satisfies the congruence" — and the natural
candidate for a short proof of such a statement is the polynomial route:

    r not in I_{p,k}   <=>   gcd(x^p - x, (x)_k - k! r) = 1  in F_p[x].

That is a genuine certificate, and it is *slower* here, not faster: Q7
measured the polynomial route at `O(k^2 log p)` against the scan's `O(g)`,
with the crossover near `k = 150` and a 0.2x ratio by `k = 400`. In the Band
II regime `k ~ 0.76p`, so the polynomial certificate loses by orders of
magnitude. Batching does not help either: one prime against many columns is
already what Band II does, and the per-column work is unchanged.

What *is* available, and already implemented: verification is embarrassingly
parallel and independent per certificate, sampling gives calibrated
confidence with complete coverage checked separately, and a single
certificate is `O(g)` in `O(1)` memory — about 32 ms for the hardest i=7
column. So the practical goal ("a referee can check any one line cheaply") is
met even though the theoretical goal ("check all of them cheaply") is not.

---

## Q16 — birthday deficit above `p^(1/3)`: declined

Theorem B covers `g <~ p^(1/3)`, which is ~77% of generic Band I columns and
**0%** of the fat cells and Band II — precisely the columns the
pre-registrations rest on. The counting argument provably runs out at
`p^(1/3)`, which is also where the phenomenon starts, so the barrier is real
rather than an artefact of one proof.

**Declined on value, not only on difficulty.** The payoff would be a *lower*
bound on `|I|`, which sharpens pre-registrations but cannot improve a kill
guarantee — a larger image makes killing *less* likely, not more. Since the
size law already agrees with measurement to 0.014% at Band II scale (Q23),
a proof would replace an accurate model with a weaker inequality.
Worth doing for a write-up that wants everything proved; not worth doing to
make the census better.

---

## Q17 — effective Chebotarev: blocked, and Q14 removed the need

The qualitative half ("a killing prime exists") is settled. The effective half
needs `disc(F_k)`, which contains `m` — 3.1 million digits at i=8 — so even
under GRH the least-prime bound is astronomically past anything testable. Asking
whether the dependence on `disc` can be avoided *for this special shape* is
reasonable, but the shape does not obviously help: the discriminant of
`(x)_k - c` genuinely involves `c`, and `c = k!m` is the large object.

**The more useful observation is that the project does not need it.** An
effective bound would let us say "checking up to `X` suffices", turning the
census into a finite verification. Q14 gives something different and, for this
purpose, sufficient: a per-column proof that a killing prime *exists*, from a
one-prime irreducibility certificate — with density `>= 1/k` by
Cameron–Cohen, which is a rate even if not a bound. The census does not need
to know in advance when it will stop; it needs to know that it will, and now
it does, column by column.

---

## Q18 — is 3003 alone, is the family the only infinite one?

Deliberately left unbudgeted, and that stands. These are the actual open
problems, and no reasoning budget settles them in a session; listing them at
MAX would only spend effort confirming that.

What can be said honestly:

- **What is known.** 3003 is the only known `N = 8`; the Lind/Singmaster/Tovey
  family is the only known infinite `N = 6` family; no `N = 5` or `N = 7` is
  known at all. Blokhuis–Brouwer–de Weger (2017) closed the small-row and
  small-value regions; MRSTT (2022) closed the interior for sufficiently large
  `t`, ineffectively (see `interior-2022.md`).
- **What this project adds.** Effectiveness in a range the 2022 theorem cannot
  reach, plus **6,067,902 recorded certificates** — every extra column of every
  member i=2..8, checked to union to exactly \([2,k_{\max}]\setminus\{K,K+1\}\)
  — and a ghost census of the same size as a by-product.
- **How to attack it.** Not head-on. Through named sub-questions with decidable
  answers, which is exactly what Q14 was: it did not settle Singmaster, but it
  did settle `k <= 5` completely and closed the termination gap operationally.
  The pattern generalises — pick a sub-question whose answer is a theorem, not
  a search.

The natural next sub-question in that style is `k = 6` for intersectivity, by
the same route that settled `k = 5`: enumerate the `2+4`, `3+3` and `2+2+2`
splits, reduce each to a Diophantine condition, and check the finitely many
survivors. `k = 5` took one Pell equation. The even-`k` obstruction in
`q14-intersective.md` §4 should dispose of `2+2+2` outright.
