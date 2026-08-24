# Q27 — `(x)_8 − c` is never intersective. **Proved.**

Worked out 2026-08-20, after Q25 (`k=6`) and Q26 (`k=7`); the last gap closed
2026-08-23 (§5). Verified by `scripts/k8_intersective.py` and
`scripts/k8_case2.py`; artifacts `results/k8_intersective.json` (historical)
and `results/k8_case2.json`; pinned by `scripts/test_k8_case2.py`.

**Status: THEOREM.** The 2026-08-20 text rested one step on Siegel's theorem
without the effective computation. That step is now closed outright — the
genus-3 curve of §4 is *solved*: its integral points are exactly the five
degenerate ones, by an elementary descent whose every class dies by
compactness or an explicit congruence (two hiding at 2-adic depth 13 and 6),
with a Runge squeeze as an independent cross-check on one class. No Siegel,
no Chabauty, no external CAS — only sympy polynomial arithmetic and integer
scans. `k = 8` joins `k ≤ 6` as fully proved.

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

What actually bounds `s` is **integrality**, not the sign of `disc` (the
2026-08-20 text claimed "`disc ≥ 0` bounds `s` to `[10, 74]`", which is
wrong — `disc ≥ 0` holds on `[9, 35] ∪ [43, 81]` and fails inside
`[36, 41]`): from `s³−84s²+1974s−12916 = (s−42)(s²−42s+210) − 4096`, `p(s)`
is an integer only if `s` is even and `(s−42) | 2048` — twenty-two values of
`s`, all within the scan; intersecting with `disc ≥ 0` leaves
`s ∈ {10, 26, 34, 44, 46, 50, 58, 74}`. (Caught in the 2026-08-23
adversarial review; the enumeration itself was complete and its conclusions
stand.)

- **Case 5** (`g` splits over **Q**): the only solution is `{1,9,25,49}`,
  i.e. `c = 0`, which is `(x)_8` itself and has rational roots.
- **Case 3** (exactly two rational roots): all six solutions have `c = 0` and
  a *reducible* remainder, so they are Case 5 in disguise. **Case 3 proper is
  empty.**

---

## 4. Case 2 is the whole difficulty — and it has a clean answer

`g = (u−β)·q(u)` with `q` an irreducible cubic. Here `c = R(β)/256`, an
infinite one-parameter family.

**The kill.** (Here `β` is always a *non-square*: an odd square `β = s²`
gives `f_c` the rational root `(s+7)/2` — excluded by the definition of
intersective — and an even `β` makes `R(β)` odd, `c` non-integral. So the
choice `(β|p) = −1` is always available.) Choose `p` with `(β|p) = −1` and
Frobenius a 3-cycle on `q`.
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

By Siegel, a curve of genus ≥ 1 has finitely many **integral** points — and
since 2026-08-23 no appeal to Siegel is needed: §5 *solves* the curve by an
elementary descent. Its integral points are exactly

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

## 5. The gap, closed (2026-08-23): the curve is solved by descent

Siegel was never needed. Since `P₆(B) ≡ P₆(0) (mod B)` and
`P₆(0) = −9458775 = −3⁷·5²·173`, no prime outside `{3, 5, 173}` can divide
both `B` and `P₆(B)` (and 2 never divides `P₆(B)` when `2 | B`: `P₆(0)` is
odd). At an integral point of `y² = −B·P₆(B)`, the valuation `v_p(−B)` is
therefore **even** for every `p ∉ {3, 5, 173}`: writing `−B = δm²` with `δ`
squarefree, `δ` divides `3·5·173` up to sign — **sixteen classes** — and
`n = y/m` is an integer with

```text
n² = δ · P₆(−δ m²).
```

Every class dies elementarily (`scripts/k8_case2.py`, certificates re-run on
every execution):

| class | fate |
|---|---|
| `δ < 0` (8 classes, i.e. `β ≥ 0`) | **compact**: `n² ≥ 0` forces `P₆(B) ≤ 0`, which holds only inside the outer bounds `(−1.41, 26.44) ∪ (48.67, 49.68)` (roots `≈ −1.4001, 26.4361, 48.6792, 49.6796`); so `B ≤ 49` and the complete list is `(B, |y|) = (0,0), (1,3840), (9,20736), (25,19200), (49,5376)` |
| `δ = 1, 15, 519, 865` | empty **mod 7** — for `δ ≡ 1 (mod 7)` this is the single computation `P₆(−m²) ≡ m¹² − m⁶ + 3 ≡ 3 (mod 7)`, a non-residue. `δ = 1` also dies by an independent **Runge squeeze** (`S = m⁶+63m⁴+651m²+269`, `P₆(−m²) − S² = 336²m⁴ + 5429760m² − 9531136`, and `S² < P₆(−m²) < (S+1)²` for `|m| ≥ 238`; nothing below) — kept as a cross-check |
| `δ = 5, 2595` | empty **mod 5⁴** |
| `δ = 3` | empty **mod 2¹⁴**: for odd `m` the value `3·P₆(−3m²)` has 2-adic valuation exactly **13** — odd, never a square; for even `m` it is `3 mod 8` |
| `δ = 173` | empty **mod 2⁹**: the value is `2⁶·(5 mod 8)` for odd `m`, `5 mod 8` for even `m` |

> **Theorem.** The integral points of `y² = −B·P₆(B)` are exactly
> `(0,0), (1,±3840), (9,±20736), (25,±19200), (49,±5376)`.

All five are degenerate (§4), so **no dangerous `β` exists and the
Chebotarev kill applies to every genuine Case-2 value. Case 2 is closed.**

Two remarks worth keeping. First, the two deep certificates were nearly
missed: a congruence battery capped at `2⁷` declares the classes `δ = 3` and
`δ = 173` *alive* — the obstructions sit at 2-adic depth 13 and 6, past the
cap. Valuation-deep obstructions need moduli beyond the valuation. Second,
there is a second road that needs only three of the certificates — the mod-7
kills of `δ = 15, 519` and the compact scan — and none of the Runge, mod-5⁴,
or deep 2-adic ones: since
`(x)_8 ≡ 0 (mod 2⁷)` identically, an intersective `c` needs `2⁷ | c`, and
the exact 2-adic valuation of `R(β) = (β−1)(β−9)(β−25)(β−49)` is pinned by
`β mod 8` — it is 4 for `β ≡ 3, 7` (so `c = R(β)/256` is not an integer)
and exactly 8 for `β ≡ 5` (so `c` is an **odd** integer and `f_c` is odd for
every `x`: no root even mod 2). Only `β ≡ 1 (mod 8)` survives, negative such
`β` would need `δ ∈ {15, 519}` (both empty mod 7), and positive ones are the
compact list. Either road ends the same way.

---

## 6. Where the ladder stands

| `k` | status | what it took |
|---|---|---|
| ≤4 | proved | elementary |
| 5 | proved | one Pell equation |
| 6 | **proved** | `t²` reduction + a Legendre symbol |
| 7 | 3+4 branch proved (Runge); 2+5 modulo one elliptic computation | two curves |
| 8 | **proved** | genus-3 curve *solved*: descent into 16 classes, all elementary |

The even/odd asymmetry from Q25 holds up: even `k` reduces the degree by half
and the work is character arguments, odd `k` does not reduce and the work is
Diophantine geometry. But `k=8` shows the even side is not free — the
*discriminant of the cofactor* is where the difficulty hides, and it grows
with `k`.

**For `k=9`** (odd, no reduction, degree 9) the branch count makes the Q26
approach painful. **`k=10` is now proved** (Q28, Magma 2026-08-23): even,
quintic `g`, no cubic-cofactor curve, (2,3) list complete. Remaining
intersective leftover is `k=9` (Siegel) and Q26 Branch A (genus-1 cubic).

As with Q25 and Q26, none of this changes the census: `c = k!m` has 3.1
million digits, and termination is certified per column by Q14 §6.
