# Q26 — `(x)_7 − c` is never intersective (with one gap, stated)

Worked out 2026-08-20, immediately after Q25; Branch B closed and the gap
restated 2026-08-23 (§5). Verified by `scripts/k7_intersective.py`,
`scripts/k7_runge.py` and `scripts/k7_branchA.py`; artifacts
`results/k7_intersective.json` (historical) and `results/k7_runge.json`;
pinned by `scripts/test_k7_runge.py`.

**Status: weaker than Q25, and deliberately labelled so.** Q25 is a proof.
This is a proof *modulo one elliptic integral-point computation* — Branch A
below. The original text said "modulo an effective Thue computation"; §5
explains why that framing was wrong on both branches, in opposite
directions: Branch B needed no solver at all and is now **closed**
unconditionally by Runge's method, while Branch A was never a Thue equation
— it is a smooth genus-1 cubic, the wall Q28's quartic *was*. That one fell
to Magma's `IntegralQuarticPoints` (2026-08-23); no analogous routine exists
for a cubic model, so Branch A is now the harder of the two.

---

## 1. Why k=7 is harder, precisely

Q25 worked because for **even** `k` the roots are symmetric about `(k−1)/2`,
so `f_c` is a polynomial in `t²` and the degree **halves**: a sextic question
became a cubic question.

For **odd** `k` the centred polynomial is *odd*, not even:

\[
(x)_7=y(y^2-1)(y^2-4)(y^2-9)=y^7-14y^5+49y^3-36y,\qquad y=x-3.
\]

The analogous identity does still exist. Writing `S(v) = (v−1)(v−4)(v−9)`,

\[
f_c(y)\,f_c(-y)=T(y^2),\qquad T(v)=c^2-v\,S(v)^2,
\]

and the same equivalence holds — `f_c` has a root mod `p` iff `T` has a
**square** root mod `p`. But \(\deg T = 7\), not 3.

> Even `k`: \(2m \to m\). Odd `k`: \(2m+1 \to 2m+1\). **No reduction.**

That is the whole difference, and it is why `k=5` needed a Pell equation,
`k=6` needed a Legendre symbol, and `k=7` needs curves.

Oddness does buy one thing: `f_c(−y) = −f_{−c}(y)`, so everything is symmetric
under `c → −c` and it suffices to search `c ≥ 0`.

---

## 2. Only two branches, not three

Intersective requires `f_c` reducible with every irreducible factor of degree
`≥ 2`. The partitions of 7 into parts `≥ 2` are `2+5`, `3+4`, `2+2+3`.

But **`2+2+3` is a `2+5` with a reducible quintic.** So asking "does `f_c` have
a quadratic factor?" covers both, and only two branches need analysis.

### Branch A — `f_c` has a quadratic factor

Matching `f_c = (y²+ay+b)(y⁵+…)` and eliminating the five quintic coefficients
leaves exactly one condition, **free of `c`**, and it involves `a` only through
`a²`. With `A = a²`:

\[
\Phi(A,b)=-A^3+5A^2b+14A^2-6Ab^2-42Ab-49A+b^3+14b^2+49b+36=0.
\]

`c` is then determined by `(a,b)`. Geometrically `Φ = 0` is the classical
*chord* condition `P(y₁) = P(y₂)` with `σ = y₁+y₂ = −a`, `π = y₁y₂ = b`.

Its leading form `−A³+5A²b−6Ab²+b³` is an **irreducible cubic form**
(`b³−6b²+5b−1` has no rational root), so the curve has three distinct points
at infinity and **Siegel's theorem gives finitely many integral points**.

Integer points with `|a| ≤ 400`: `c = 0` (trivial) and

| `c` | `(a,b)` | split |
|---|---|---|
| `±17472` | `(2, 13)` | `2+5` |
| `±459648` | `(3, 38)` | `2+5` |

### Branch B — `f_c` has a cubic factor and no quadratic factor (`3+4`)

Matching `f_c = (y³+ay²+by+d)(y⁴+…)` leaves **two** conditions free of `c`.
The first is linear in `d`:

\[
d\,(3a^2-2b-14)=-a\,(a^2-3b-7)(a^2-b-7),
\]

and substituting into the second gives a plane curve `C(A,b) = 0` of degree 5,
whose leading form `2A⁵−9A⁴b+12A³b²−A²b³−9Ab⁴+4b⁵` is squarefree with five
distinct roots — Siegel again.

**The degenerate locus matters — and it is finite by pure algebra.** When
`3a²−2b−14 = 0` the linear solve is invalid and must be handled separately.
My first pass skipped it and reported "no nontrivial 3+4 exists", which was
wrong: the locus is not empty, and it is exactly where the only nontrivial
solution lives. (No Siegel needed here: substituting `b = (3a²−14)/2` into
the vanishing numerator gives `(7/4)·a³(a²−4) = 0`, so `a ∈ {0, ±2}` exactly
— three points, no search.) The solution:

\[
(a,b,d)=(2,-1,14)\ \Longrightarrow\ c=-896,
\]

giving `(x)_7 + 896 = (x³+…)(x⁴+…)`. The generic branch contributes nothing
but `c = 0`.

---

## 3. Every candidate dies for one reason

The complete nontrivial list is `c = ±896, ±17472, ±459648`, and

| `c` | factorisation | `5 ∣ c`? |
|---|---|---|
| 896 | `2⁷·7` | no |
| 17472 | `2⁶·3·7·13` | no |
| 459648 | `2⁷·3³·7·19` | no |

`rad(7!) = 210` divides none of them, and specifically **5 divides none of
them**. Since any 7 consecutive residues mod 5 contain a multiple of 5,
`(x)_7` vanishes *identically* on \(\mathbf{F}_5\), so `f_c ≡ −c ≢ 0`. No root
mod 5.

> **`(x)_7 − c` is never intersective.**

---

## 4. Cross-validation

Two independent methods agree exactly. Brute-force factorisation of `(x)_7 − c`
for **every** integer `c` in `[0, 400000]` finds exactly two reducible cases
with all parts `≥ 2` — `c = 896` and `c = 17472` — which is precisely what the
curve analysis predicts in that range (`459648` lies beyond it and was found on
the curve). That the chord curve and blind factorisation produce the same list
is the main reason to believe the case analysis is not missing a branch.

---

## 5. The gap — restated 2026-08-23: half closed, half corrected

The 2026-08-20 text under this heading proposed closing both branches with
PARI/GP's `thue`. That was wrong twice over, and in opposite directions.

**Branch B is CLOSED — no solver was ever needed.** The quintic's leading
form factors over `Q`:

```text
F5 = (2A² − 5Ab + 4b²) · (A³ − 2A²b − Ab² + b³) = Q2 · K3
```

`Q2` is positive definite (disc `−7`) and `K3(1, w) = w³ − w² − 2w + 1` is
the cyclic cubic of discriminant `49` — the field `Q(2cos(2π/7))`, `k = 7`
announcing itself. Two coprime factors mean **Runge's method** applies, and
the definiteness of `Q2` makes it sharp: the two `Q2`-branches at infinity
are complex, so every real point of the curve runs along one of `K3`'s three
real branches, and the degree-4 polynomial `W2` (constructed so its
expansion along each branch has no positive power of `A` — a property
re-verified symbolically on every run) is bounded there. On an integer point
with `A ≥ 111184`, the integer `m = W2(A,b)` would have to lie within that
branch's certified error of that branch's limit. The three conjugate limits
are `w0 = −294ω² + 69ω + 449 ≈ −94.199, +421.478, −381.278` (trace exactly
`−54`), their distances to the nearest integer are `0.199, 0.478, 0.278`,
and the certified per-branch errors are `0.0073, 0.121, 0.0006` — each
smaller than its own margin. So there are no integer points at all with
`A ≥ 111184`; a modular-sieve sweep over `0 ≤ A` below the threshold
(`A ≥ 0` suffices: `A = a²`) finds exactly nineteen points, all with
`A ≤ 36`. Complete list of 3+4 values: `c = 0` and `c = ±896`. Both die at
5. `scripts/k7_runge.py` is the certificate — exact rational interval
arithmetic end to end, ~3–4 minutes, no external CAS; an earlier draft of
its final inequality mixed a 4-term Newton anchor with a 5-term series tail
and was caught in adversarial review — the committed chain uses one anchor
throughout. Artifact: `results/k7_runge.json`.

**Branch A was never a Thue equation.** `Φ(A,b) = 0` is a *smooth* plane
cubic — no affine singular point, three smooth points at infinity — hence
**genus 1**: `thue` solves `F(x,y) = m` for homogeneous `F` and simply does
not apply. Via the pencil of lines through `(0,−1)` the Jacobian is

```text
E:  Y² = X³ − 1764·X + 28224  =  X³ − 42²·X + 168²
```

with **trivial torsion** (gcd of `#E(F_p)` is 1) and visible integer points
— `(0,168)`, `(21,21)`, `(28,28)` — each therefore of **infinite order**:
rank ≥ 1, `E(Q)` is infinite, and there is no rank-0 shortcut. Siegel still
gives finiteness of the *integral* points of `Φ`; certifying that list needs
the elliptic-logarithm method on a **cubic model** (Stroeker–Tzanakis), and
no installed tool provides it — Magma's `IntegralPoints` takes Weierstrass
models, and the birational map does not preserve integrality. **Branch A
stays blocked there** — the analogue of Q28's rank-2 quartic, which Magma
has since solved (2026-08-23). The quartic case had a routine; this one has
none, which is why Q28 closed and Branch A did not.

**Chebotarev narrows what the missing list is for** (`scripts/k7_branchA.py`
§2). A 2+5 value is killed at `p` iff the quadratic is inert and the quintic
rootless — Frobenius `(transposition, derangement)`. Among the transitive
subgroups of `S₅` only `S₅` has an *odd* derangement (its twenty
`(2,3)`-elements); `C₅, D₅, F₂₀, A₅` have only 5-cycles, all even. `C₅` and
`A₅` admit no `C₂`-quotient, so the full product appears and
`(inert, 5-cycle)` kills; `S₅` kills through `(2,3)` even with a shared
subfield. The one escape: **quintic group `D₅` or `F₂₀` whose unique
quadratic subfield is `Q(√(a²−4b))`** — then every class fixes a root and no
unramified prime kills. The trap is real:
`(x²+x−1)(x⁵−2)` — `Q(√5)` on both sides — has a root modulo every prime
tested, while the mismatched control `(x²+x+1)(x⁵−2)` dies at 11. This is
the k=7 instance of the even-derangement mechanism of Q27/Q28 (the pair
witnesses the group theory; neither factor product is of the form
`(x)_7 − c`). Both recorded candidates' quintics are `S₅` (a `(2,3)` pattern
mod 37 resp. 103 certifies it), so they were never in the trap. The
trichotomy assumes the quintic irreducible; a `2+2+3` counterexample would
need its own no-kill configuration among two quadratics and a cubic, not
analysed here — but every `2+2+3` value lies on the Branch-A curve too, so
the enumeration below covers it.

The honest statement now: **the 3+4 branch is a theorem; the 2+5 branch has
no counterexample below the search bound, only finitely many can exist
anywhere, and any counterexample with irreducible quintic would need a
`D₅/F₂₀` quintic with matched subfield *and* `7! | c`.** The remaining
certificate is the integral-point list of one rank-≥1 elliptic cubic —
the elliptic-log method for general genus-1 models: Stroeker–Tzanakis,
*Computing all integer solutions of a genus 1 equation*, Math. Comp. 72
(2003), 1917–1933. (Magma's only non-Weierstrass cubic routine,
`SIntegralDesbovesPoints`, takes special Desboves-form cubics only, which
`Φ` is not.)

---

## 6. Where this leaves the general question

- `k ≤ 4` — Q14, elementary.
- `k = 5` — Q14, one Pell equation.
- `k = 6` — Q25, proved outright via the `t²` reduction.
- `k = 7` — this note: the 3+4 branch is a theorem (Runge, §5); the 2+5
  branch modulo one elliptic integral-point computation.
- `k = 8` — Q27: **proved** (2026-08-23). The predicted "surface" never
  materialised; the difficulty was a genus-3 curve in Case 2, since solved
  by descent (`scripts/k8_case2.py`). (This bullet said "open" long after
  Q27 existed — a stale line, fixed with the closure.)

The pattern is now clear and worth stating: **even `k` is the easy side.** The
odd cases need Diophantine geometry (Pell at 5, Siegel at 7) while the even
cases need only character arguments on a halved degree. Anyone continuing this
should do `k = 8` before `k = 9`.

Operationally, as with Q25, this changes nothing for the census: our `c = k!m`
has 3.1 million digits, and termination is certified per column by Q14 §6.
