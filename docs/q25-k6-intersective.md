# Q25 — `(x)_6 − c` is never intersective

Worked out 2026-08-20, as a follow-on to Q14. Verified by
`scripts/k6_intersective.py`; artifact `results/k6_intersective.json`.

**Status.** This is new, derived in a single session, and the argument below
should be read critically before it is quoted anywhere. Every step is checked
computationally, and the load-bearing step (Case 2) is a short Chebotarev
argument that a referee can check by hand in a few minutes. But it has not
been reviewed by anyone yet.

---

## What Q14 expected, and what actually happened

Q14 closed `k ≤ 5` and left `k ≥ 6` open, with this forecast:

> `k = 5` took one Pell equation; `k = 6` will take three curves, and the
> even-`k` obstruction of §4 should kill the `2+2+2` case outright.

The forecast was reasonable and it is what I started on — enumerate the `2+4`,
`3+3` and `2+2+2` splits of the sextic, reduce each to a Diophantine
condition. That route works: `3+3` collapses to `a | 189` (five values of
`c`), and `2+2+2` collapses to two finite conics (one value, `c = 2240`).

But it is the wrong decomposition. **The problem is not about the sextic.**

---

## 1. The reduction: a sextic question is a cubic question

For even `k` the roots `0..k−1` are symmetric about `(k−1)/2`, so with
`t = 2x − (k−1)` they become `±1, ±3, …, ±(k−1)`. At `k = 6`:

\[
2^6\,(x)_6=(t^2-1)(t^2-9)(t^2-25)=R(t^2),\qquad
R(u)=u^3-35u^2+259u-225,
\]

hence

\[
2^6 f_c = g(t^2),\qquad g(u)=u^3-35u^2+259u-(225+64c).
\]

**The key equivalence.** As `t` runs over \(\mathbf{F}_p\), `t²` runs over
exactly the squares. Therefore

> \(f_c\) has a root mod \(p\)  ⟺  \(g\) has a root mod \(p\) **that is a
> square mod \(p\)**.

(`u = 0` is a square, so it is included.) Checked on 1,035 \((c,p)\) pairs
with zero failures.

So the question is about a **cubic together with a quadratic-residue
condition**, and the natural case analysis runs over the factorisation of `g`
over **Q** — a trichotomy — not over the five splits of `f_c`. Degree halves.

---

## 2. The trichotomy

### Case 1 — `g` irreducible over **Q**

Its Galois group acts transitively on three roots, Jordan gives a
fixed-point-free element, and Chebotarev turns that into a positive density of
primes at which `g` has **no root at all** mod `p` — a fortiori no square
root. Density `1/3` for `S₃`, `2/3` for `C₃`.

This is the generic case: 3,986 of the 4,001 values `c = 30j` with
`|c| ≤ 60000`. Sampled 200 of them; all have a killing prime, the largest
being `p = 19`.

### Case 2 — `g = (u − β)·q(u)`, `q` an irreducible quadratic

This is the only infinite case, and the one that matters.

`β` is an integer (rational root of a monic integer cubic), and
`c = (β−1)(β−9)(β−25)/64`, so the family is parametrised by a single integer.
Two forced non-squares:

- **`β` is a non-square**, else `t² = β` gives `f_c` a rational root;
- **`Δ := −3β² + 70β + 189` is a non-square**, since it is the discriminant of
  the complementary quadratic and `q` is irreducible.

Now consider \(\mathbf{Q}(\sqrt\beta,\sqrt\Delta)\). Since both are
non-squares there is a positive density of primes with

\[
\left(\tfrac{\beta}{p}\right)=\left(\tfrac{\Delta}{p}\right)=-1
\]

— density `1/4` if `βΔ` is a non-square, `1/2` if `βΔ` is a square (the two
fields then coincide). At any such `p`:

- `q` has no root mod `p`, because `Δ` is a non-residue;
- so `β` is the **only** root of `g` mod `p`;
- and `β` is a non-residue.

No root of `g` mod `p` is a square, so by the equivalence **`f_c` has no root
mod `p`.** Not intersective. ∎

**This needs no finiteness.** It closes an infinite family in one stroke,
which is why the "three curves" forecast was pessimistic.

Verified on all 273 candidates with `|β| ≤ 3000` satisfying `rad(6!) | c`: the
prime predicted by the two Legendre symbols kills every one, largest needed
`p = 113`.

### Case 3 — `g` splits completely over **Q**

Then `β₁+β₂+β₃ = 35` and `Σβᵢβⱼ = 259`. Writing `s = β₁+β₂` gives
`β₁β₂ = s²−35s+259`, so `β₁,β₂` are integers only if

\[
-3s^2+140s-1036 \ \text{is a perfect square.}
\]

The leading coefficient is negative, so the left side is non-negative only for
`s ∈ [9.22, 37.44]` — **a bounded conic, 28 integers to test.** The only
solution is `{β} = {1, 9, 25}`, i.e. `c = 0`, which is `(x)_6` itself and has
rational roots.

**This is where a counterexample would have had to live.** For the covering to
work you would need three rational `βᵢ`, none a square, but with `β₁β₂β₃` a
square so that the all-non-residue pattern `(−1,−1,−1)` is impossible. That is
exactly the "odd square-relation" of `q14-intersective.md` §4, and here
`β₁β₂β₃ = 225 + 64c`. The conic leaves no room for it.

---

## 3. Result

> **`(x)_6 − c` is never intersective.** Q14's answer extends from `k ≤ 5` to
> `k ≤ 6`.

Cross-checked by brute force over `c = 30j`, `|c| ≤ 60000` (and `≤ 300000` in
a longer run): every reducible case with all factors of degree `≥ 2` is
accounted for by the trichotomy, and none lacks a killing prime.

The two splits Q14 named are subsumed rather than needed:

| Q14's split | where it lands | outcome |
|---|---|---|
| `3+3` | `g` irreducible (Case 1) | `c ∈ {0, 14850, 77350, 14731200, 11084273200}`; the two with `rad(6!) \| c` die at `p = 7` |
| `2+2+2` | `g` has a rational root (Case 2) | only `c = 2240`, and `3 ∤ 2240`, so no root mod 3 |
| `2+4` | Case 2 | the infinite family, killed by the two-symbol argument |

Both finite lists were computed independently and agree, which is why the
first version of the script is worth keeping in the history: it is a
consistency check on the cleaner argument.

---

## 4. What this suggests for general even `k`

The reduction is not special to 6. For any even `k`, with `m = k/2`,

\[
2^k f_c = g(t^2),\qquad g(u)=\prod_{j=1}^{m}\bigl(u-(2j-1)^2\bigr)-2^k c,
\]

and the same equivalence holds: **`f_c` has a root mod `p` iff `g` has a
square root mod `p`.** So every even `k` reduces to a degree-`m` problem —
`k = 8` to a quartic, `k = 10` to a quintic.

Two structural facts carry over unchanged:

1. **Any irreducible factor of `g` of degree ≥ 2 that must supply the root is
   vulnerable to Jordan/Chebotarev.**
2. **Rational roots of `g` must be non-squares** (else `f_c` gets a rational
   root), so each contributes a quadratic character that can be forced to
   `−1`.

Intersectivity therefore requires the characters to be *dependent* — an odd
square-relation, §4 again. For `k = 6` the relation would have to hold among
three integers pinned by `e₁ = 35, e₂ = 259`, and no such triple exists.

**The natural next question (`k = 8`) is now well-posed and materially
easier than `k = 6` looked before this reduction:** it is a quartic `g` with
`e₁ = 1+9+25+49 = 84`, and the analogue of Case 3 is a surface rather than a
conic. That is the first place where the method might genuinely need more than
elementary bounding.

Note the inversion this produces relative to intuition: **even `k` is the
easy side.** Odd `k` has no `t²` reduction at all — which is why `k = 5`
needed a Pell equation while `k = 6` needed a Legendre symbol.

---

## 5. What it means for the project

Nothing operationally — and that is worth stating plainly. Our `c = k!·m` has
3.1 million digits at `i=8`, so no result about small `c` reaches it, and
Q14 §6 already closed the termination gap *per column* via the one-prime
irreducibility certificate. That remains the mechanism the census actually
relies on.

What this adds is to the general question rather than to the tool: the
smallest genuinely open degree for `(x)_k − c` intersectivity moves from 6 to
7, and the even case now has a uniform method rather than a case-by-case
hunt.
