# Q27 — `(x)_8 − c` is never intersective (one gap, stated)

Worked out 2026-08-20, after Q25 (`k=6`) and Q26 (`k=7`). Verified by
`scripts/k8_intersective.py`; artifact `results/k8_intersective.json`.

**Status.** Same shape as Q26: the case analysis is complete and every branch
closes, but one step rests on Siegel's theorem without the effective
computation being carried out. Stronger than Q26 in one respect — the
exceptional locus here turns out to be *explicitly* the degenerate points,
which is a much sharper landing than "no counterexample below the bound".

Q26 predicted this would be the hard case: "the analogue of Q25's Case 3 is a
surface rather than a bounded conic." **That prediction was wrong.** Case 3
and Case 5 both stayed bounded; the difficulty moved somewhere else entirely
— into a genus-3 curve in Case 2.

---

## 1. The reduction still works

`k=8` is even, so Q25 applies. With `t = 2x−7` the roots `0..7` become
`±1,±3,±5,±7`:

\[
2^8(x)_8=(t^2-1)(t^2-9)(t^2-25)(t^2-49)=R(t^2),
\]
\[
R(u)=u^4-84u^3+1974u^2-12916u+11025,\qquad g(u)=R(u)-256c,
\]

and `f_c` has a root mod `p` iff `g` has a root mod `p` that is a **square**
mod `p`. Checked on 495 `(c,p)` pairs, 0 failures.

The cost of `k=8` over `k=6` is that `g` is a **quartic**, so the trichotomy
becomes five cases.

---

## 2. Two cases need no computation

**Case 1 — `g` irreducible.** Transitive ⟹ Jordan ⟹ a derangement ⟹
Chebotarev gives a positive density of primes at which `g` has no root at all.

**Case 4 — `g` = two irreducible quadratics.** Both discriminants are
non-squares, so a positive density of primes has *both* as non-residues
(density 1/4, or 1/2 if their product is a square), and at those primes `g`
has no root at all. **This case can never be intersective, unconditionally** —
no side condition is required.

---

## 3. Two cases are bounded, and both collapse to `c = 0`

The symmetric functions `e₁ = 84`, `e₂ = 1974`, `e₃ = 12916` are **fixed** —
only `e₄ = 11025 − 256c` varies. So if two roots have sum `s`, their product is
forced:

\[
p(s)=\frac{-s^3+84s^2-1974s+12916}{84-2s},\qquad
\mathrm{disc}=s^2-4p(s)=-\frac{s^3-126s^2+3948s-25832}{s-42}.
\]

For large `|s|` that behaves like `−s²`, so `disc ≥ 0` **bounds `s` to
[10, 74]** — 65 integers.

- **Case 5** (`g` splits over **Q**): the only solution is `{1,9,25,49}`,
  i.e. `c = 0`, which is `(x)_8` itself and has rational roots.
- **Case 3** (exactly two rational roots): all six solutions have `c = 0` and
  a *reducible* remainder, so they are Case 5 in disguise. **Case 3 proper is
  empty.**

---

## 4. Case 2 is the whole difficulty — and it has a clean answer

`g = (u−β)·q(u)` with `q` an irreducible cubic. Here `c = R(β)/256`, an
infinite one-parameter family.

**The kill.** Choose `p` with `(β|p) = −1` and Frobenius a 3-cycle on `q`.
Then `q` has no root mod `p`, so `β` is the only root of `g` mod `p`, and it is
a non-residue — no root of `g` is a square, so `f_c` has no root mod `p`.

Such `p` exist with positive density **unless** `Q(√β)` is the quadratic
subfield of `q`'s splitting field, i.e. unless `β·disc(q)` is a square. (If
`Gal(q) = C₃` there is no quadratic subfield at all and the kill is automatic.)

**The key observation.** `q` is not arbitrary — it is the *divided difference*

\[
q(u)=\frac{R(u)-R(\beta)}{u-\beta},
\]

so its coefficients, and hence its discriminant, are **polynomials in β**:

\[
\operatorname{disc}(q)=-16\,P_6(\beta),
\]
\[
P_6(B)=B^6-126B^5+5271B^4-82564B^3+570591B^2-5779998B-9458775,
\]

with `P₆` irreducible and coprime to `B`. Since `−16` is `−1` times a square,
the dangerous condition `β·disc(q) = □` becomes

> \(y^2=-\beta\,P_6(\beta)\) — squarefree, degree 7, **genus 3**.

By Siegel, a curve of genus ≥ 1 has finitely many **integral** points. Its
integral points are

\[
\beta\in\{0,\;1,\;9,\;25,\;49\}
\]

— exactly the roots of `R`, plus `0`. Every one is excluded: `1, 9, 25, 49` are
perfect squares (so `t² = β` would give `f_c` a rational root), and they give
`c = 0`; `β = 0` gives `R(0)/256 = 11025/256`, not an integer.

**So no dangerous `β` exists, and the Chebotarev kill always applies.**

Direct check: 52 candidates with `|β| ≤ 1200` satisfying `rad(8!) | c`, zero
dangerous, every one killed by a prime ≤ 29. (Extended to `|β| ≤ 4000`: 197
candidates, still zero dangerous.)

---

## 5. The gap

Siegel gives finiteness of integral points on the genus-3 curve, and Baker's
method makes that bound effective — but **the effective computation was not
carried out**; the search covers `|β| ≤ 20000`.

This is a better position than Q26's, though. There the found points were
"whatever the search turned up". Here they are exactly `{0} ∪ {roots of R}` —
the points forced by the geometry — which is what one expects the complete
answer to be. Closing it properly means Chabauty–Coleman or an effective
Baker bound on `y² = −B·P₆(B)`.

---

## 6. Where the ladder stands

| `k` | status | what it took |
|---|---|---|
| ≤4 | proved | elementary |
| 5 | proved | one Pell equation |
| 6 | **proved** | `t²` reduction + a Legendre symbol |
| 7 | modulo effective Thue | two curves, Siegel |
| 8 | modulo effective Siegel | genus-3 curve, explicit exceptional set |

The even/odd asymmetry from Q25 holds up: even `k` reduces the degree by half
and the work is character arguments, odd `k` does not reduce and the work is
Diophantine geometry. But `k=8` shows the even side is not free — the
*discriminant of the cofactor* is where the difficulty hides, and it grows
with `k`.

**For `k=9`** (odd, no reduction, degree 9 with partitions `2+7`, `3+6`,
`4+5`, `2+2+5`, `2+3+4`, `3+3+3`, `2+2+2+3`) the branch count alone makes the
Q26 approach painful. The better next target is `k=10`: even, `g` is a
**quintic**, and the Case-2 analogue would be a divided difference again, so
the same "disc is a polynomial in β" trick applies. Whether its exceptional
curve stays as clean as genus 3 is the question worth asking.

As with Q25 and Q26, none of this changes the census: `c = k!m` has 3.1
million digits, and termination is certified per column by Q14 §6.
