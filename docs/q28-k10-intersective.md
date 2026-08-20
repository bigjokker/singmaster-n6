# Q28 — `(x)_10 − c` is never intersective (one genus-1 gap)

Worked out 2026-08-20. Verified by `scripts/k10_intersective.py`; artifact
`results/k10_intersective.json`.

**The headline is not the result, it is the cost.** `k=10` is *cheaper* than
`k=8`. Difficulty in this family is **not monotone in `k`** — it is governed
entirely by whether `g` acquires an irreducible **cubic** factor.

---

## 1. Why k=8 was hard and k=10 is not

The kill is always the same: find `p` at which every rational root of `g` is a
non-residue and every higher-degree factor is rootless. "Rootless" means
Frobenius is a **derangement** of that factor's roots.

| n | derangements | odd ones |
|---|---|---|
| 2 | 1 | 1 |
| **3** | **2** | **0** |
| 4 | 9 | 6 |
| 5 | 44 | 20 |
| 6 | 265 | 135 |

**`n=3` is the only degree whose derangements are all even.** They are
3-cycles, which lie in `A₃`, so `(disc|p) = +1` is forced. If `Q(√β)` happens
to be the cubic's quadratic subfield, the character condition on `β` then
conflicts with rootlessness, and the kill fails on an exceptional locus —
which is precisely the genus-3 curve Q27 needed for `k=8`, whose cofactor was
a cubic.

At `k=10` the cofactor in the infinite case is a **quartic**, and 4-cycles are
odd derangements. The character is free, and **no exceptional curve arises at
all**. The genus-6 curve one might have expected (`disc(q)` has degree 12 in
`β`) is simply never invoked.

---

## 2. The four cases

With `t = 2x−9`, `2^10 (x)_10 = R(t²)` and

\[
R(u)=u^5-165u^4+8778u^3-172810u^2+1057221u-893025,\qquad g(u)=R(u)-1024c,
\]

so `g` is a quintic, and `f_c` has a root mod `p` iff `g` has a root that is a
square mod `p`.

**(5) irreducible.** Jordan + Chebotarev. No computation. This is 8000 of the
8001 values scanned with `rad(10!) | c`.

**(1,4) — the infinite family.** Unconditional, per §1. Verified on 66
candidates with `|β| ≤ 1500`; all die by `p ≤ 29`.

**≥2 rational roots — and this one is fully rigorous.** `e₁..e₄` are fixed, so
two rational roots satisfy a single constraint `C(b₁,b₂)=0`. Its leading form
is

\[
b_1^4+b_1^3b_2+b_1^2b_2^2+b_1b_2^3+b_2^4=\frac{b_1^5-b_2^5}{b_1-b_2},
\]

which has **no real zeros** (they would need `b₁⁵=b₂⁵` with `b₁≠b₂`). A
definite leading form makes the curve **compact**, so a bounded search is a
*proof*, not a sample. The explicit bound is `max(|b₁|,|b₂|) < 2680`, and the
complete search to it returns exactly the ten pairs from `{1,9,25,49,81}` —
all giving `c=0`.

> Worth flagging: my first pass searched only to 400 and would have been
> presented as complete. The bound had to be computed, not assumed.

**(2,3) — no rational root.** Matching `(u²+au+b)(u³+du²+eu+f)` and
eliminating `d,e,f` leaves one curve, *quadratic in `b`*, so `b` is integral
iff its discriminant is a square:

\[
y^2 = 5a^4+1320a^3+126456a^2+5102240a+72824400,
\]

a squarefree quartic — **genus 1**. Exactly two candidates with `rad(10!) | c`
survive, and both die:

| `c` | split | killed at |
|---|---|---|
| 1,395,418,752,000 | (2,3) | `p = 11` |
| 2,235,340,800 | (2,3) | `p = 13` |

---

## 3. The gap

**One elliptic integral-point computation**, for the (2,3) branch. Siegel gives
finiteness; the search covers `|a| ≤ 6000`.

That is the mildest gap in the whole ladder — strictly easier than Q26's Thue
equations or Q27's genus-3 curve, and squarely in the range where standard
software returns a provably complete answer.

---

## 4. Ladder, corrected

| `k` | status | gap |
|---|---|---|
| ≤6 | **proved** | — |
| 7 | modulo effective Thue | two curves |
| 8 | modulo effective Siegel | genus 3 |
| **10** | modulo one elliptic computation | **genus 1** |

The ordering by difficulty is `6 < 10 < 7 < 8`, not `6 < 7 < 8 < 10`. Anyone
extending this should target the even `k` whose partitions avoid a 3 in the
cofactor, and should expect the *odd* `k` and the *cubic-cofactor* cases to be
where the work actually is.

Scaling, for the record: the case count is `p(k/2)` (3, 5, 7, 11, 15, 22, … for
`k = 6, 8, 10, 12, 14, 16`), and where an exceptional curve *is* needed its
genus is `(m−1)(m−2)/2` with `m = k/2` — 1, 3, 6, 10, 15 for `k = 6..14`.
Quadratic growth in the genus is the real ceiling, and it only engages on
cubic-cofactor cases.

As with Q25–Q27, none of this touches the census: `c = k!m` has 3.1 million
digits, and termination is certified per column by Q14 §6.
