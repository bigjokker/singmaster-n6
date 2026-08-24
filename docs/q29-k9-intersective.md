# Q29 — `(x)_9 − c` is never intersective (Branch A proved; B, C stated)

Worked out 2026-08-20; **Branch A closed by Runge's method 2026-08-24**
(`scripts/k9_branchA_runge.py`, §3b below — the branch where every
nontrivial candidate lives, including the ladder's first `rad(k!)`-passer).
Verified by `scripts/k9_intersective.py`; artifacts
`results/k9_intersective.json` (historical) and `results/k9_branchA.json`;
pinned by `scripts/test_k9_intersective.py`.

Branches B and C are the remaining shape-of-Q26 gap: their case analysis is
complete and every point found has `c = 0`, but Siegel's finiteness is not
backed by an effective computation. (`k=8` and `k=10` have both since been
proved outright — 2026-08-23: `scripts/k8_case2.py` by descent, and Q28 by
Magma's `IntegralQuarticPoints`.)

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

### 3b. Branch A CLOSED (2026-08-24): Runge, because the leading form splits

In `A = a²` the chord condition is a plane **quartic** — not k=7's cubic —

```text
Φ₉(A,b) = −A⁴ + 7A³b + 30A³ − 15A²b² − 150A²b − 273A² + 10Ab³ + 180Ab²
          + 819Ab + 820A − b⁴ − 30b³ − 273b² − 820b − 576,
```

irreducible and smooth: **genus 3**. That would be Q27-territory — except its
leading form **factors over `Q`**:

```text
L = −(A − b) · (A³ − 6A²b + 9Ab² − b³),
```

a linear factor times the cyclic cubic form of discriminant **81** — the field
`Q(2cos(2π/9))`, `k = 9` announcing itself exactly as `Q(2cos(2π/7))` did in
Q26 Branch B. Two coprime factors mean **Runge's method** applies
(`scripts/k9_branchA_runge.py`, exact rational interval arithmetic end to
end, modelled on `scripts/k7_runge.py`):

- **Four real channels** at infinity: slope `1` (rational) and the three
  conjugate roots `β ≈ 0.2831, 0.4260, 8.2909` of `s³−9s²+6s−1`.
- **The rational channel traps an integer.** Along it
  `b = A − 10 + 9/A + 90/A² + 630/A³ + …`, so `m = b − A ∈ Z` satisfies
  `m + 10 ∈ (0, 1)` once `A ≥ A_FAR`: no integer value exists. Even the trap
  line itself is empty for every `A`: `Φ₉(A, A−10) = −27(A²−10A−12)`,
  discriminant `148`, not a square.
- **The cyclic channels trap `W₁`.** The degree-4 integer polynomial
  `W₁ = A³b − 6A²b² − 20A²b − 9A² + 9Ab³ + 70Ab² + 118Ab + 90A − b⁴ − 20b³
  − 82b² − 180b` has no positive powers of `A` along any cyclic branch
  (re-verified symbolically on every run); its limits are
  `ω = −270β² + 2241β − 594 ≈ 18.827, 311.712, −573.538` (trace exactly
  `−243 = −3⁵`), at distances `0.173, 0.288, 0.462` from the nearest
  integers, while the certified trap widths at `A_FAR` are
  `0.038, 0.175, 0.024`. An integer point would make `W₁` an integer within
  the trap — impossible.
- **Coverage**: a cell-wise domination argument shows every real point with
  `A ≥ A_FAR = 34416` lies in one of the four channel tubes; a CRT modular
  sweep (complete by construction — the moduli product exceeds the `b`-range)
  settles `0 ≤ A < A_FAR` in about a minute.

**The theorem**: the integer points of `Φ₉` with `A ≥ 0` are exactly
**23 points, all with `A ≤ 49`** — twenty chords (`c = 0`), the two
candidate points `(4, 21)` and `(9, 62)` carrying `c = ±2,630,880` and
`c = ±176,774,400`, and one point `(17, 8)` with `A` not a square. So the
complete list of Branch-A values is `c = 0, ±2630880, ±176774400`; the two
nontrivial ones die at `13` and `7`. **Every shape of 9 containing a 2 —
`2+7`, `2+2+5`, `2+3+4`, `2+2+2+3` — is now a theorem.** No Siegel, no
Magma. Artifact: `results/k9_branchA.json`.

## 4. Branches B and C — only `c = 0` found; still open

**B.** The surviving y² condition is **quadratic** in the cubic's constant
term, not linear, so this needs a resultant rather than a substitution — degree
18 in `a`, 9 in `b`. All 42 integer points found have `c = 0`. Classified
2026-08-24: the resultant curve is **irreducible with leading form `a¹⁸`** —
a single repeated factor, so the Runge split that closed Branch A provably
does not transfer to this projection.

**C.** Three conditions in four unknowns; `p₀` is linear in the first, and
eliminating `p₁` by resultant splits the locus into

- `(2p₂ − 3p₃² + 30)⁴` — the **degenerate locus**, where the `p₀`-solve is
  invalid and which must be handled separately (this is exactly the trap that
  cost me the `c = ±896` case in Q26). **Correction (2026-08-24): it is NOT
  empty**, and the original emptiness check was vacuous (it tested an
  expression that still had `p₁` free, so it tested nothing). Full
  elimination — `Res_{p₁}(E₃|_{deg}, Res_{p₀}(E₂|_{deg}, E₁|_{deg}))` has
  integer roots `p₃ ∈ {0, ±2}` only — shows the locus holds exactly **six
  rational points**: `(p₃,p₂,p₁,p₀) = (±2,−9,∓2,8), (±2,−9,∓18,0),
  (0,−15,±10,24)`. **Every one has `c = 0`** (e.g. `y⁴−15y²+10y+24` has
  roots `{−4,−1,2,3}` — a chord factorisation), which is the statement that
  matters. Proved by algebra, not by a search; pinned in
  `scripts/test_k9_intersective.py`.
- a main component of degree 14 in `p₂`, 24 in `p₃`, whose 65 integer points
  found all have `c = 0`.

## 5. A second obstruction mechanism

Worth recording even though branch C carries only `c = 0`.

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
k=9 only because no nontrivial 4+5 value was found.

The general lesson, and it has now cost me twice: **check every transitive
group, not the generic one.** Both obstructions are invisible if you assume
`Gal = Sₙ`.

## 6. Cross-check and the gap

Factoring `(x)_9 − c` for **every** multiple of 210 up to 4,200,000 returns
exactly one reducible case with all parts ≥ 2 — `c = 2,630,880` — which is
precisely what the curve analysis predicts.

**Gap, narrowed 2026-08-24:** Branch A — the source of that one case, and of
every nontrivial candidate — is **closed** (§3b, Runge). What remains is
Branches B and C: Siegel gives finiteness on their curves, but no effective
tool applies (B's projection has the non-split leading form `a¹⁸`; C's main
component has degree `(14, 24)`). Searches cover `|a| ≤ 200` (B),
`|p₃| ≤ 120` (C), and find only `c = 0`.

## 7. Ladder

| `k` | status | gap |
|---|---|---|
| ≤6 | **proved** | — |
| 7 | 3+4 branch proved (Runge); 2+5 modulo one elliptic computation | two curves |
| 8 | **proved** (2026-08-23: descent) | — |
| **9** | **Branch A proved** (2026-08-24: Runge, the disc-81 split); B and C modulo effective Siegel | two curves |
| 10 | **proved** (Magma 2026-08-23: 15 \(a\)-values) | — |

Difficulty ordering is `6 < 10 < 7 < 9 < 8` by the severity of the gap each
`k` originally posed, which continues to
have nothing to do with the size of `k`. What actually drives cost is (i) odd
vs even, and (ii) which obstruction mechanisms the factor degrees admit.
(With Branch A closed, k=9's remaining gap is smaller than k=7's: two
searched-empty branches versus one live rank-3 cubic.)

As with Q25–Q28, none of this touches the census: `c = k!m` has 3.1 million
digits, and termination is certified per column by Q14 §6.
