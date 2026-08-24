# Q29 — `(x)_9 − c` is never intersective (one gap, stated)

Worked out 2026-08-20. Verified by `scripts/k9_intersective.py`; artifact
`results/k9_intersective.json`.

Same shape as Q26: the case analysis is complete and every branch closes,
but Siegel's finiteness is not backed by the effective computation. (`k=8`
and `k=10` have both since been proved outright — 2026-08-23:
`scripts/k8_case2.py` by descent, and Q28 by Magma's
`IntegralQuarticPoints`.)

---

## 1. Odd k, again no reduction

With `y = x−4`,

\[
(x)_9=y(y^2-1)(y^2-4)(y^2-9)(y^2-16)=y^9-30y^7+273y^5-820y^3+576y,
\]

which is **odd** in `y`. That yields `f_c(−y) = −f_{−c}(y)`, enough to restrict
to `c ≥ 0` and nothing more. Q25's degree-halving needs `f_c` to be a
polynomial in `t²`, which only happens for even `k`.

## 2. Eight shapes, three branches

Partitions of 9 into parts ≥ 2: `9, 2+7, 3+6, 4+5, 2+2+5, 2+3+4, 3+3+3,
2+2+2+3`. But a shape containing a 2 is just "`f_c` has a quadratic factor",
and one containing a 3 but no 2 is "`f_c` has a cubic factor". So:

| branch | means | covers |
|---|---|---|
| A | a quadratic factor | 2+7, 2+2+5, 2+3+4, 2+2+2+3 |
| B | a cubic, no quadratic | 3+6, 3+3+3 |
| C | a quartic, no quadratic or cubic | 4+5 |

plus irreducible, killed by Jordan + Chebotarev with no computation.

## 3. Branch A — where everything lives

Eliminating the septic's coefficients leaves one condition, free of `c`,
involving `a` only through `a²` — the chord curve for `P(y₁)=P(y₂)`. Its
constant part factors as

\[
-(b+1)(b+4)(b+9)(b+16),
\]

the same shape as k=7's `−(w+1)(w+4)(w+9)`. The leading form is squarefree with
four distinct roots, so there are ≥3 points at infinity and Siegel applies.

| `c` | split | `rad(9!) = 210 ∣ c`? | dies at |
|---|---|---|---|
| ±176,774,400 | (2,7) | no | `p = 7` |
| **±2,630,880** | (2,7) | **yes** | `p = 13` |

`2,630,880 = 2⁵·3⁴·5·7·29` is **the first candidate anywhere in this ladder to
pass the `rad(k!)` filter.** Every earlier k had its survivors eliminated for
free by that congruence. Here the filter does not help and the polynomial has
to be killed on its merits:

```
x² − 10x + 45  ·  x⁷ − 26x⁶ + 241x⁵ − 956x⁴ + 2044x³ − 3824x² − 12096x − 58464
```

## 4. Branches B and C — empty

**B.** The surviving y² condition is **quadratic** in the cubic's constant
term, not linear, so this needs a resultant rather than a substitution — degree
18 in `a`, 9 in `b`. All 42 integer points have `c = 0`.

**C.** Three conditions in four unknowns; `p₀` is linear in the first, and
eliminating `p₁` by resultant splits the locus into

- `(2p₂ − 3p₃² + 30)⁴` — the **degenerate locus**, where the `p₀`-solve is
  invalid and which must be handled separately (this is exactly the trap that
  cost me the `c = ±896` case in Q26). It is empty here.
- a main component of degree 14 in `p₂`, 24 in `p₃`, whose 65 integer points
  all have `c = 0`.

## 5. A second obstruction mechanism

Worth recording even though branch C is empty.

Killing a 4+5 needs Frobenius to be a derangement of **both** factors at once.
If the splitting fields were independent that has density
`(9/24)(44/120) ≈ 0.14 > 0`. I checked the sign character across all 9
transitive subgroups of `S₄` and all 20 of `S₅` — **no blocked pair**, and I
briefly concluded 4+5 was unconditionally dead.

That was wrong, because the sign character is not the only shared quotient.
`C₄` and `F₂₀` can fuse over `C₄`, and there:

- `C₄`'s derangements are its three non-identity elements — they **avoid** the
  identity coset;
- `F₂₀`'s derangements are its four elements of order 5, which generate `C₅` =
  the kernel of `F₂₀ ↠ C₄` — they lie **entirely inside** it.

So one demands image ≠ 0 and the other forces image = 0. Incompatible: the kill
is blocked.

This is a **second obstruction mechanism, distinct from the `n=3` one** that
produced Q27's genus-3 curve. That one came from a *single* factor whose
derangements all lie in `Aₙ`; this one comes from *two* factors whose
derangement sets sit in incompatible cosets of a shared quotient. It is moot at
k=9 only because no 4+5 solution exists.

The general lesson, and it has now cost me twice: **check every transitive
group, not the generic one.** Both obstructions are invisible if you assume
`Gal = Sₙ`.

## 6. Cross-check and the gap

Factoring `(x)_9 − c` for **every** multiple of 210 up to 4,200,000 returns
exactly one reducible case with all parts ≥ 2 — `c = 2,630,880` — which is
precisely what the curve analysis predicts.

**Gap:** Siegel gives finiteness on all three branch curves; the effective
computation was not carried out. Searches cover `|a| ≤ 300` (A), `|a| ≤ 200`
(B), `|p₃| ≤ 120` (C).

## 7. Ladder

| `k` | status | gap |
|---|---|---|
| ≤6 | **proved** | — |
| 7 | 3+4 branch proved (Runge); 2+5 modulo one elliptic computation | two curves |
| 8 | **proved** (2026-08-23: descent) | — |
| **9** | modulo effective Siegel | three branch curves |
| 10 | **proved** (Magma 2026-08-23: 15 \(a\)-values) | — |

Difficulty ordering is `6 < 10 < 7 < 9 < 8` by gap severity, which continues to
have nothing to do with the size of `k`. What actually drives cost is (i) odd
vs even, and (ii) which obstruction mechanisms the factor degrees admit.

As with Q25–Q28, none of this touches the census: `c = k!m` has 3.1 million
digits, and termination is certified per column by Q14 §6.
