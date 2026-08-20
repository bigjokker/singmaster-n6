# Q26 — `(x)_7 − c` is never intersective (with one gap, stated)

Worked out 2026-08-20, immediately after Q25. Verified by
`scripts/k7_intersective.py`; artifact `results/k7_intersective.json`.

**Status: weaker than Q25, and deliberately labelled so.** Q25 is a proof.
This is a proof *modulo an effective Thue computation I did not run*. The
structure is complete and both curves are provably finite by Siegel; what is
missing is the certificate that the integer-point lists are exhaustive rather
than merely unbeaten below the search bound. That gap is mechanical, not
conceptual — see §5.

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

**The degenerate locus matters.** When `3a²−2b−14 = 0` the linear solve is
invalid and must be handled separately. My first pass skipped it and reported
"no nontrivial 3+4 exists", which was wrong: the locus is not empty, and it is
exactly where the only nontrivial solution lives —

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

## 5. The gap, stated plainly

Siegel proves both curves have finitely many integral points, and Baker's
method makes the bound effective. **I did not carry out the effective
computation.** So the lists above are complete up to `|a| ≤ 400` and
cross-checked against brute force to `c ≤ 400000`, but not proved exhaustive.

Closing it is mechanical: run PARI/GP's `thue` on the two curves. Branch A is
a genuine cubic Thue equation and `thue` handles it directly; Branch B is
degree 5 and may need `thue` on the leading form plus a bounded search. Neither
is research — it is an afternoon with the right tool, which this machine does
not currently have installed.

Until then the honest statement is: **no counterexample exists below the search
bound, and only finitely many can exist anywhere.**

---

## 6. Where this leaves the general question

- `k ≤ 4` — Q14, elementary.
- `k = 5` — Q14, one Pell equation.
- `k = 6` — Q25, proved outright via the `t²` reduction.
- `k = 7` — this note, modulo effective Thue.
- `k = 8` — open, and now the natural target: even, so the reduction applies
  and it becomes a **quartic** `g` with `e₁ = 1+9+25+49 = 84`. The analogue of
  Q25's Case 3 is a surface rather than a bounded conic, which is the first
  place the even method needs a genuinely new idea.

The pattern is now clear and worth stating: **even `k` is the easy side.** The
odd cases need Diophantine geometry (Pell at 5, Siegel at 7) while the even
cases need only character arguments on a halved degree. Anyone continuing this
should do `k = 8` before `k = 9`.

Operationally, as with Q25, this changes nothing for the census: our `c = k!m`
has 3.1 million digits, and termination is certified per column by Q14 §6.
