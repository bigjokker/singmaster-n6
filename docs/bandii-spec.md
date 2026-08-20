# Band II sweep spec (i=8)

Claude, 2026-08-19. Implemented in scripts/bandii_kernel.py and
scripts/bandii_sweep.py. Bat: scripts/run-bandii-sweep.bat.

Pre-flight ran on this box before the bat was handed over (numpy 2.2.6,
168.5 M elem/s, stragglers exact). **This file does not start a job.**
The bat refuses if 
esults/bandii_sweep.json exists.

# BAND II SWEEP — SPEC

For Grok to implement. Written 2026-08-19. Claude.
Target file if committed: `docs/bandii-spec.md`.

**This is a spec, not a start order.** Nothing here starts a job. The bat is
written only after Sam says to run it, and it must refuse if
`results/bandii_sweep.json` already exists (same guard as fat-image).

Everything numeric below was computed and checked while writing this file.
The kernel was validated end-to-end against the live-machine straggler log
(§8.3) — it reproduces all four survivors **and their witness j values**.

---

## 0. What this is / is not

**Is:** an exhaustive, certificate-producing modular sweep of every Band II
column `k ∈ [K+2, k_max]` for `i=8`, using primes just above `N/2`.

**Is not:** a next-prime sweep from `k`. The refresher's Band II ban stands
and is untouched: for `k` in Band II and `k < p ≤ N/2`, the theorem gives
`p | m`, so `0 ∈ I_{p,k}` and no such prime can ever kill. This spec never
tests one.

The regime here is different in kind, not degree: `α = ⌊N/p⌋ = 1`,
`β = ⌊K/p⌋ = 0`. There are no Z-slabs, no digit-0 silence, no NONE/PART
classification, and `r(p) ≠ 0` for every prime in the window. That regime
does not exist anywhere in Band I.

---

## 1. The prize

For any prime `p > k`, Lucas on the single digit of `k` gives the complete
column image over **all** `n`, including `n ≥ p`:

```
{ C(n,k) mod p : n ≥ k }  =  {0} ∪ { C(k+j, j) mod p : 0 ≤ j < p−k }  =  I_{p,k}
```

so `r(p) ∉ I_{p,k}` is an **unconditional certificate** that column `k` never
represents `m`. Not a sample. Not a null.

Band II is one of the two regions still open for `i=8`. With Band II and the
Band I remnant closed, every `k ∈ [2, k_max] \ {K, K+1}` carries a
certificate, and that is a theorem:

> **N( C(F₁₈F₁₉, F₁₆F₁₉) ) = 6, exactly.**

That is the first theorem-shaped end of the i=8 line. This sweep is one of
the two jobs standing between here and there.

---

## 2. Fixed objects and the window

```
N      = F_18 F_19 = 10803704
K      = F_16 F_19 = 4126647          (NEVER F_18 F_17 — that is K+1)
d      = N − K     = 6677057
N/2    = 5401852
k_max  = 5182637                      largest k with C(2k,k) ≤ m
m      = C(N,K),  3120256 digits,  log10 m = 3120255.2212
```

`k_max` re-derived independently here by bisection on lgamma:
`log10 C(2k,k) = 3120254.781` at `k = 5182637`, `3120255.383` at `k+1`.
**Agrees with the refresher.** Assert this in pre-flight (§8.1).

**Columns swept:** `k = K+2 … k_max` = `4126649 … 5182637`, **1,055,989
columns**. `K` and `K+1` are the family and are never tested.

**Live prime window:** `p ∈ (N/2, d] = (5401852, 6677057]`, ≈81,500 primes.

Why exactly that window, for every Band II `k`:

| prime range | α, β | r(p) | usable? |
|---|---|---|---|
| `k < p ≤ N/2` | ≥2 | 0 (theorem) | no — 0 ∈ image, cannot kill |
| `N/2 < p ≤ d` | 1, 0 | `C(N−p, K) mod p`, **never 0** | **yes — all of them** |
| `d < p ≤ N` | 1, 0 | 0 (`K > N−p`, borrow) | no |

`r(p) ≠ 0` on the middle band because both `N−p < p` and `K < p`, so no
factor of `p` can enter `C(N−p,K)`.

Note the first prime is tight: `2·5401853 − N = 2`. Assert `2p > N` and
`p ≤ d` and `p > k_max` for every prime used.

---

## 3. The kernel

Do **not** build Pascal rows. Use factorials.

`C(n₀,k) ≡ r (mod p)` with `F(x) := x! mod p`:

```
F(n₀) ≡ r · F(k) · F(n₀−k)   (mod p)
```

Substituting `n₀ = k + b`:

> **Column `k` survives prime `p`  ⟺  ∃ b ∈ [0, p−k−1] with**
> **`F[k+b] ≡ s·F[b] (mod p)`, where `s = r(p)·F[k] mod p`.**

`b` here is exactly the `j` that fat-image logged. `n₀ = k+b` ranges over
`[k, p−1]`, which is precisely the image range. No signs anywhere.

**Half of that range is redundant.** The involution `(k−1−x)_k = (−1)^k (x)_k`
sends `b ↦ g−1−b` on `[0,g)` and multiplies the value by `(−1)^k`, so:

| | upper half | test on `b ∈ [0, ⌈g/2⌉)` |
|---|---|---|
| `k` even | repeats | `F[k+b] ≡ s·F[b]` |
| `k` odd | negates | that, or `F[k+b] ≡ −s·F[b]` |

The `−s` branch costs a **subtraction**, not a second multiply, since
`(−s)F[b] = p − sF[b]` for `r ≠ 0` — so both parities do half the modular
multiplications, on half the temporary memory. Exact, not approximate.

The reported `b` is unchanged: the least full-scan hit is
`min(first s-hit, g−1 − last (−s)-hit)`, because a `−s` hit at half-index `b`
is a full hit at `g−1−b`, which *decreases* as `b` grows. So the survivor set
**and** every recorded witness index are byte-identical to the full scan.
Measured ~2.2× at Band II parameters. `scan_ks_full` remains as the reference;
`USE_HALF_SCAN = False` restores it for A/B.

`r(p) = 0` is refused rather than scanned: `0 = C(x,k)` for every `x < k`, so
such a prime certifies nothing, and a factorial-table mask cannot represent
that — it would report every column killed.

Vectorized, this is **one multiply, one modulo, one compare per element**:

```python
mask = Fa[k:p] == (s * Fa[0:p-k]) % p
```

Overflow: `s < p` and `Fa < p`, so `s·Fa[b] < p² ≈ 2.92e13 < 2^63`. int64 is
safe. Do not use uint64, do not use float.

**Measured throughput: 162 M elem/s, single core, numpy 2.4, int64.**

Validated by brute force against `math.comb` on 1000 random `(p,k,r)` for
`p ∈ {11, 29, 101, 211, 1009}`: **0 mismatches.** Reproduce this in
pre-flight.

### numpy is a hard requirement

The stack is Python + gmpy2; gmpy2 has no vectorization and there is no
pure-Python path to 7.9e11 element-ops. `numpy >= 1.24` must be importable
before the bat is written. If numpy is not on the box, that is the blocker
to clear first — nothing else in this spec matters until it is.

---

## 4. r(p), three independent ways — compute all three, assert equal

Getting `r(p)` wrong once makes every "kill" at that prime worthless, so it
is worth three computations and an assert. All three verified at
`p = 5401853`, all give **r = 1275205**.

**(a) Factorial table.** `r = F[N−p] · F[K]⁻¹ · F[N−p−K]⁻¹ mod p`.

**(b) Direct falling factorial.** `min(K, N−p−K) = N−p−K` (1,275,204 at p₁),
so `r = (∏_{i=0}^{N−p−K−1}(N−p−i)) · F[N−p−K]⁻¹ mod p`. ~1.3 s in Python.

**(c) Closed form — cheap, and the real check.** Write `δ = 2p − N`. Then
`N − p = p − δ`, and

```
C(p−δ, j) ≡ (−1)^j C(j+δ−1, j)   (mod p)
```

hence

> **`r(p) = (−1)^K · C(K+δ−1, δ−1) mod p`,  δ = 2p − N.**

The lower index is `δ−1`, not `K`. For the first 16 primes `δ ≤ 422`, so
this costs a few hundred modmuls — no table needed. At `p₁`, `δ = 2`,
`K` odd, so `r = −(K+1) = 5401853 − 4126648 = 1275205`. ✓

Also cross-check (a) against the engine's existing `binom_mod_lucas(N,K,p)`
for at least the first prime. That reuses the code path with 4114 validated
trials behind it.

---

## 5. Algorithm

Prime-major, one factorial table live at a time. Simple, resumable,
memory-safe.

```
build F for p₁                         # ~2 s, 43 MB int64
r₁ ← §4 (three ways, assert equal)
pass 1: all 1,055,989 columns          # ≈7.89e11 elements
        → survivors S₁ (expect ~1.03e5)
free F, build F for p₂
pass 2: columns in S₁                  # ≈9e10 elements
        → S₂ (expect ~1.26e4)
pass 3: columns in S₂                  # ≈1e10
...
pass r: columns in S_{r−1}
stop at r = 14 (§10)
```

Total ≈ `1.13 × pass 1`. Pass 1 is 99% of the cost; everything after is
noise.

**Memory:** one `int64[p]` table = 43 MB per worker. With 8 workers, ~350 MB
plus per-column temporaries of ≤10 MB each. Given the earlier 119 GB
incident: do **not** cache tables across primes, do **not** hold more than
one table per worker, and do **not** keep the boolean masks — reduce each to
`(any, argmax)` immediately and discard.

**Parallelism.** Cost per column is `∝ (p−k)`, which is linear in `k`, so
equal-column-count chunks are badly unbalanced. Split `[K+2, k_max]` into
`4·W` contiguous chunks of **equal Σg** (closed form: chunk boundaries at
equal steps of `(p−k)²`). One chunk = one work unit = one checkpoint line.

**Windows spawn:** each worker rebuilds `F` itself (~2 s). Do not pickle the
43 MB array to workers.

---

## 6. Prime list

The first 16 primes above `N/2`. All satisfy `2p > N`, `p ≤ d`, `p > k_max`.

```
p1  = 5401853    δ=2      p9  = 5401993    δ=282
p2  = 5401861    δ=18     p10 = 5401999    δ=294
p3  = 5401867    δ=30     p11 = 5402003    δ=302
p4  = 5401897    δ=90     p12 = 5402011    δ=318
p5  = 5401901    δ=98     p13 = 5402021    δ=338
p6  = 5401951    δ=198    p14 = 5402051    δ=398
p7  = 5401969    δ=234    p15 = 5402057    δ=410
p8  = 5401973    δ=242    p16 = 5402063    δ=422
```

Hard-code these and re-verify primality at load. The cap is 14; p15/p16 are
listed only so an anomaly escalation has somewhere to go without a new spec.
An escalation means the size law failed, not that a run was long; the test is
in `scripts/sizelaw.py` and each pass records expected-vs-observed.

---

## 7. Pre-registered kill curve

Same size law that produced walk-369 (42/4/1/0, max_run 6), even-g 0.929,
and the held-out stragglers (3.97 predicted / 4 measured). Zero fitted
parameters. Computed over all 1,055,989 columns with the actual primes
above.

| after r primes | expected still alive | frac g even | mean surviving k |
|---:|---:|---:|---:|
| 1 | 1.026 × 10⁵ | 0.658 | 4 536 120 |
| 2 | 1.26 × 10⁴ | 0.786 | 4 450 918 |
| 3 | 1816 | 0.875 | 4 392 087 |
| 4 | 289.9 | 0.930 | 4 350 633 |
| 5 | 49.4 | 0.962 | 4 320 257 |
| 6 | 8.78 | 0.979 | 4 297 132 |
| 7 | 1.61 | 0.989 | 4 278 952 |
| 8 | 0.300 | 0.994 | 4 264 287 |
| 10 | 0.0108 | 0.998 | 4 242 079 |
| 12 | 4.06 × 10⁻⁴ | 1.000 | 4 226 056 |
| 14 | 1.56 × 10⁻⁵ | 1.000 | 4 213 949 |

Three things are pre-registered, not just the count:

1. **Counts.** Poisson bands: pass 1 in `[102 000, 103 200]`, pass 4 in
   `[257, 324]`, pass 6 in `[3, 15]`.
2. **Parity.** Survivors should skew even-`g` exactly as tabled — `g` even
   means `|I_g| ≈ g`, `g` odd means `≈ g/2` (Pascal reflection). Prime gaps
   are even so parity is invariant along a column's chain.
3. **Location.** Survivors should pile up near `K+2` (largest `g/p`), with
   the mean surviving `k` walking *down* toward `K` as `r` grows. This is
   the sharpest test: it is a prediction about *which* columns survive, not
   how many.

---

## 8. Pre-flight checks — mandatory, all cheap, run before any sweep

Abort on any failure. Total runtime under a minute.

### 8.1 Constants
```
assert N == 10803704 and K == 4126647 and K == fib(16)*fib(19)
assert K != fib(18)*fib(17)                      # that is K+1
assert d == 6677057 and k_max == 5182637
assert log10 C(2·k_max, k_max) <= log10 m < log10 C(2k_max+2, k_max+1)
assert kmax - (K+2) + 1 == 1055989
```

### 8.2 Kernel identity vs brute force
For `p ∈ {11, 29, 101, 211, 1009}`, 200 random `(k, r)` each: the kernel
`∃b: F[k+b] ≡ r·F[k]·F[b]` must equal
`any(comb(n0,k) % p == r for n0 in range(k,p))`.
**Expected: 1000 cases, 0 mismatches.** Also assert Wilson: `F[p−1] == p−1`.

### 8.3 Straggler regression — the one that matters

Run the kernel on `k = 4126622 … 4126646` at `p = 5401853`. These are Band I
columns, but every prime here exceeds them so the same kernel applies, and
the machine already logged the answer.

**Required output, exactly:**

| k | g | g even | witness b |
|---|---|---|---|
| 4126624 | 1275229 | no | 273671 |
| 4126638 | 1275215 | no | 268500 |
| 4126642 | 1275211 | no | 2006 |
| 4126643 | 1275210 | yes | 554555 |

and nothing else survives. The four `b` values are the `j = 273671, 268500,
2006, 554555` in the refresher's straggler line. Then at `p = 5401861` all
four must die.

I ran this: 2.2 s total including the table build, exact match on all four
`k` **and** all four witnesses. If the implementation does not reproduce
this table byte-for-byte, it is wrong; do not proceed.

### 8.4 r(p) triple check
§4 (a), (b), (c) must agree at `p₁`, and `r(p₁) == 1275205`. Also check (a)
vs `binom_mod_lucas(N,K,p₁)`.

### 8.5 Live spot check
`k = 4126649 … 4126688` (the 40 hardest columns) at `p₁`: expect ~8.6
survivors by the size law. I measured **9**.

---

## 9. Output

Follow the fat-image convention. `results/bandii_sweep.jsonl` is the
checkpoint, `results/bandii_sweep.json` is the record.

**jsonl, one line per (prime, chunk):**
```json
{"prime_index":1,"p":5401853,"delta":2,"r":1275205,
 "k_lo":4126649,"k_hi":4160000,"n_cols":33352,
 "n_survivors":7331,
 "survivors":[{"k":4126649,"g":1275204,"g_even":true,"b":41277}, ...],
 "seconds":38.4}
```

Log **`g` and parity of everyone still up**, at every stage, plus the
witness `b`. Cap the inline `survivors` array at pass 1 (7331 per chunk is
fine to write; if a chunk blows past 20k, write `k` + `b` only). From pass 2
on, log every survivor in full.

**json record:** constants, prime list, per-pass counts, per-pass observed
even-`g` fraction and mean surviving `k` next to the §7 predictions, total
seconds, and the full final survivor list (expected: empty).

`certificate` field, on a clean run:

> Every k in [4126649, 5182637] carries a modular impossibility certificate
> `r(p) ∉ I_{p,k}` for some prime p ∈ {p₁…p_r}. Band II contains no extra
> representation of C(F₁₈F₁₉, F₁₆F₁₉). Unconditional.

---

## 10. Stop rules

- **Cap at 14 primes.** Do not auto-extend. Same discipline as walk-369.
- **Clean finish** = zero columns alive after ≤14. Write the record and stop.
- **Anything alive at 14** → stop, log, hand back. Do not extend, do not
  start an exact check, do not touch Band I.
- **Counts far off §7** (e.g. ~300 alive at 12 instead of 4×10⁻⁴): that is
  **not** a broken size law. The likeliest cause is that `r(p) = C(N−p, K)`
  is not behaving like a generic residue in this window — which would be a
  new arithmetic fact about `m` and more interesting than the sweep. Log the
  survivors with `g`, parity, `b`, and `k`; the `k`-distribution and the
  parity fraction will say immediately which it is. Hand back, do not chase.
- **Parity fraction badly off** (e.g. 0.5 at pass 3 instead of 0.875) with
  counts on target: suspect an implementation bug in `g`, not a discovery.

A column surviving all 14 primes is **not** a representation. It means the
image filter failed there. Escalation, if it ever happens, is more primes
from the same window (≈81,500 available), then exact — but that is a new
decision, not part of this run.

---

## 11. Cost

| | |
|---|---|
| pass-1 elements | 7.89 × 10¹¹ |
| measured throughput | 162 M elem/s / core (numpy 2.4, int64) |
| pass 1, single core | 1.35 h |
| all 14 passes, single core | ~1.5 h |
| **8 workers** | **~12–20 min** |
| F table build | ~2 s per prime per worker |
| peak RSS | ~350 MB at 8 workers |

This is a short job, not an overnight. Checkpoint per chunk anyway.

---

## 12. Not in scope

This sweep does **not** touch:

- **Band I remnant** — the 89,195 Stage-3 hang-guard columns, and
  `10⁶ … K` (~3.13 M columns, never swept). That is the second job, by
  Z-jump: skip to the first NONE/PART-lower prime past the Q3 slab, test at
  most 12; a column still up is an anomaly only if its round's expected
  survivor count was below \(10^{-2}\) (see `zjump-spec.md`, Escalation —
  run length alone is the wrong trigger). The hang-guard `p−k > 20000` is
  obsolete — walk-369 established the image tail is size-driven and maxes at
  6.
- **The same kernel serves it.** For a Band I column, `s = r(p)·F[k]` with
  `r(p)` the two-digit `C(α,β)C(n₀,k₀)`; everything else is identical. Group
  by prime, one `F` table per prime, same `(any, argmax)` reduction. Worth
  writing the kernel once, in a module both jobs import.
- i=9, i=10, nearby ≥10⁹, Stage 4, E, Q1, more fat-image. Finish i=8.

---

## 13. Do not

- Do not next-prime-sweep Band II from `k`. The ban stands; this spec never
  does it.
- Do not use `K = F₁₈F₁₇`. That is `K+1`.
- Do not mark `K` or `K+1` impossible.
- Do not test primes `≤ N/2` or `> d` against Band II — provably `m ≡ 0`,
  guaranteed survivor, pure waste.
- Do not build `m`.
- Do not use `math.comb` in the hot path (Lucas digits use
  `_binom_mod_prime`).
- Do not cache factorial tables across primes.
- Do not run without §8.3 passing.
- Do not call a clean pass "Singmaster proved". It closes Band II for `i=8`.
  The theorem needs Band I too.