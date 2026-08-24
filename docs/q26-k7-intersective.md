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

**Branch A was never a Thue equation — and here is the actual proof.** The
2026-08-23 text argued "`Φ` is a *smooth* plane cubic, hence genus 1, hence
`thue` does not apply." **That inference is invalid**: a Thue curve
`N(x,y) = m z³` is itself a smooth plane cubic of genus 1, so smoothness
excludes nothing. The correct test is projective. A cubic is affinely
equivalent over `Q` to a binary-form equation `N(x,y) = m` **only if** the
Hessian of its homogenisation vanishes identically on the line at infinity
(equivalently: all three infinite points are inflections *and* their tangents
are concurrent). For `Φ`,

```text
Hessian(Φ)|_{z=0} = 392·(b − 2A)·(2b − A)²
```

which is not identically zero and is not even divisible by the leading form
`b³ − 6Ab² + 5A²b − A³` — the roots `b = 2A`, `b = A/2` are not roots of
`θ³−6θ²+5θ−1` (its values there are `−7` and `1/8`). So **no** point at
infinity is an inflection, and no change of variables over `Q` can turn `Φ`
into a Thue equation. (Affine: a projective change would move the line at
infinity and would not preserve `Φ(Z)`.) PARI's `thue` is excluded by
theorem, not by inspection. (A genuine Thue cubic is carried as the control in
`test_k7_runge.py`: for `x³+2y³=7` the Hessian *does* vanish identically at
infinity.) The same pass also excludes the other special-form escape:
`ψ₃` of the Jacobian is irreducible over `Q`, so `E` has no rational
3-isogeny, `Φ` is not `Q`-equivalent to a Desboves cubic, and Magma's
`SIntegralDesbovesPoints` cannot apply either.

Via the pencil of lines through `(0,−1)` the Jacobian is

```text
E:  Y² = X³ − 1764·X + 28224  =  X³ − 42²·X + 168²
```

with **trivial torsion** (gcd of `#E(F_p)` is 1; Magma `TorsionSubgroup`
2026-08-24: abelian group of order 1) and visible integer points
— `(0,168)`, `(21,21)`, `(28,28)` — each therefore of **infinite order**:
rank ≥ 1, `E(Q)` is infinite, and there is no rank-0 shortcut. Magma
`Rank(E)` the same day printed **`3 true`**, so rank = 3 unconditionally —
which matters, because the size of a Stroeker–Tzanakis reduction lattice is
`rank + 2`.

```text
> E := EllipticCurve([-1764, 28224]);
> E;
Elliptic Curve defined by y^2 = x^3 - 1764*x + 28224 over Rational Field
> TorsionSubgroup(E);
Abelian Group of order 1
> Rank(E);
3 true
```

**Why the obvious shortcut is wrong, concretely.** Write the birational map
out. Projecting from `(0,−1)` and reducing gives `σ : Φ → E`,

```text
X = N_X(A,b)/A²,   Y = N_Y(A,b)/A³,   σ⁻¹(O) = (A,b) = (0,−9)
```

(the identity `Y² = X³ − 1764X + 28224` holds modulo `Φ`; both are pinned).
`σ⁻¹(O)` is an **affine** rational point, not one of the points at infinity —
it has to be, since the infinite points are irrational. So `X∘σ` has a pole
on the affine curve and the map does not preserve integrality. The failure is
not theoretical: the genuine integral point `(A,b) = (4,13)`, which carries
the real Branch-A value `c = ±17472`, maps to

```text
σ(4,13) = (1345/4, 48959/8)  ∈ E(Q),  NOT an integral point of E
```

Of the 21 known integral points of `Φ` (all with `|A| ≤ 95`; listed in
`scripts/k7_branchA.py` as `PHI_INTEGRAL_POINTS`), **8** have integral image
on `E`, **10** have finite non-integral image, and **3** are poles of `σ`
(`A = 0`). So `IntegralPoints(E)` returns a demonstrably **wrong** list —
the k=7 instance of exactly the trap that `Sage`'s `E.integral_points()`
would have been on Q28.

**The obstruction, stated properly.** The three points at infinity of `Φ` are
a single Galois orbit: the leading form is irreducible with discriminant `49`,
so `C₃` acts simply transitively. Hence every Galois-stable effective divisor
supported there has degree divisible by 3, and **no `Q`-rational function of
degree 1 or 2 has poles only at infinity**. A Weierstrass model needs a
degree-2 function (`x`, pole divisor `2·O`) and a quartic model needs one too
(`T`); neither exists over `Q`. Therefore no `Q`-rational Weierstrass or
quartic model has the *same* integral points as `Φ`, and `IntegralPoints`,
`SIntegralPoints` and `IntegralQuarticPoints` cannot be handed this problem
on any such model.

Two honest caveats, both established while checking the above. First, a model
whose integral points merely *contain* `Φ`'s always exists — rescale
`(X,Y) → (u²X, u³Y)` to clear denominators — so the barrier is not
impossibility but **circularity**: choosing `u`, or choosing `S` for
`SIntegralPoints(E,S)`, requires the denominators of the very points one is
trying to find. Second, the obstruction is about the **ground field**, not the
curve. Over `K = Q(θ)`, `θ³−6θ²+5θ−1 = 0`, the orbit splits, `κ = b − θA`
becomes a degree-2 function, and `Φ` becomes a genuinely **integral quartic
over `O_K = Z[θ]`**:

```text
s² = (−3θ²+12θ+16)κ⁴ + (−56θ²+252θ+168)κ³ + (−294θ²+1470θ+588)κ²
     + (−432θ²+3100θ+652)κ + (385θ²+1246θ+385),   s = 2α(κ)A + β(κ)
```

with leading coefficient the square `(−6θ²+33θ−14)²`, and `κ, s ∈ Z[θ]` at
every known integral point of `Φ`. That is precisely
`IntegralQuarticPoints`' input shape — but over `O_K`, and Magma has no
number-field version. (`K` is the same cyclic cubic field of conductor 7 that
Branch B's `K3` produced: `θ = 1/(2+t)`, `t = 2cos(2π/7)`. Both branches of
`k=7` live over it.)

Siegel gives finiteness of the integral points of `Φ`; certifying the list
needs the elliptic-logarithm method on a **cubic model** (Stroeker–Tzanakis),
and no installed or online tool provides it. **Branch A stays blocked
there** — the analogue of Q28's rank-2 quartic, which Magma solved on
2026-08-23. The quartic case had a routine; this one has none, which is why
Q28 closed and Branch A did not.

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

**Two new theorems (2026-08-24) that shrink what is left.**

*The real half of Branch A is closed, elementarily.* If `a² − 4b > 0` the two
roots `y₁ ≠ y₂` are real, so `c` is a **real chord** of
`P(y) = y(y²−1)(y²−4)(y²−9)`. `P`'s six critical values are
`±95.8419…, ±23.1490…, ±12.3588…`, all under `96` in absolute value, and `P`
is strictly monotone beyond the outermost critical point. So a chord value
satisfies `|c| ≤ 95.85`, both endpoints lie in `[−Y₀, Y₀]` with
`Y₀ = 3.1042…` (the root of `P = 96`), and therefore

```text
|a| = |y₁+y₂| ≤ 2Y₀ < 7   and   |b| = |y₁y₂| ≤ Y₀² < 10.
```

Enumerating that finite box on `Φ` returns **exactly the 21 chord points,
every one with `c = 0`**. So the whole `c = 0` locus — the obstruction that
makes any congruence or descent argument vacuous, exactly as in Q28 — lives
in the real region, and the real region is now finished by hand. **Every
surviving Branch-A point has `a² − 4b < 0`**: the quadratic factor is complex.

*`7! | c` removes `F₂₀` from the no-kill door.* An intersective `c` needs
`7! | c`, so `|c| ≥ 5040`, which exceeds every critical value of `P`; hence
`f_c` has exactly **one** real root. By the previous theorem the quadratic
factor is complex, so the **quintic** carries that single real root — two
conjugate pairs, so complex conjugation acts as a `(2,2)`-element, which is
*even*, so `disc(quintic) > 0` and `Q(√disc)` is **real**. But `F₂₀`'s unique
quadratic subfield *is* `Q(√disc)` (`F₂₀ ∩ A₅ = D₅`, index 2), while the door
demands it equal `Q(√(a²−4b))`, which is now **imaginary**. Contradiction:
**`F₂₀` is impossible**, and only `D₅` survives. Since `D₅ ⊂ A₅`, that adds a
new necessary condition — `disc(quintic)` must be a perfect **square** — which
both recorded candidates fail (they are `S₅`, as above). The
group-theoretic trichotomy itself is unchanged; what changed is that the
arithmetic of `7! | c` closes one of its two doors.

The honest statement now: **the 3+4 branch is a theorem; on the 2+5 branch
the real half is a theorem too, and what remains is the complex-quadratic
half, where no counterexample exists among the 21 known integral points of
`Φ` (all `|A| ≤ 95`), only finitely many can exist anywhere, and any
counterexample would need a `D₅` quintic with matched imaginary subfield,
square discriminant, *and* `7! | c`.** The remaining certificate is the
integral-point list of one rank-3 elliptic cubic — the elliptic-log method
for general genus-1 models: Stroeker–Tzanakis,
*Computing all integer solutions of a genus 1 equation*, Math. Comp. 72
(2003), 1917–1933. Nothing available implements it: Magma's `IntegralPoints`
and `SIntegralPoints` take Weierstrass models, `IntegralQuarticPoints` takes
`y² = quartic` over `Z`, and `SIntegralDesbovesPoints` takes Desboves cubics —
and `Φ` is provably none of the three. Searching is not a substitute: 21
integral points of `Φ` are known, **all with `|A| ≤ 95`**, the 14 with `A` a
square giving `|c| ∈ {0, 17472, 459648}` — nothing new. A larger sweep is
not in this repository.

One more piece of geometry worth recording: the curve actually needed is not
`Φ` but its double cover `Ψ : Φ(a², b) = 0`, since `A` must be a perfect
square. `A` has three simple zeros (`b = −1, −4, −9`) and three simple poles
(the points at infinity) on `Φ`, so `Ψ → Φ` is ramified at all six and
Riemann–Hurwitz gives **genus 4**. Faltings therefore bounds the *rational*
Branch-A solutions, not merely the integral ones — strictly stronger than
Siegel, and equally ineffective. `Ψ` is bielliptic over a positive-rank
elliptic quotient, which is exactly the configuration in which elliptic-curve
Chabauty does **not** apply.

---

## 6. Where this leaves the general question

- `k ≤ 4` — Q14, elementary.
- `k = 5` — Q14, one Pell equation.
- `k = 6` — Q25, proved outright via the `t²` reduction.
- `k = 7` — this note: the 3+4 branch is a theorem (Runge, §5); the 2+5
  real half is a theorem too; the complex half is rank 3 (Magma `3 true`)
  and still needs Stroeker–Tzanakis on the cubic.
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
