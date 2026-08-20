# Q14 — Can `(x)_k − c` be intersective?

Answered 2026-08-20. Every claim below was verified computationally; the
scripts are `intersective_search.py` and `termination_certificate.py`.

**Why it matters.** For this project the question is not decorative:

| | |
|---|---|
| column `k` survives prime `p` | `(x)_k − k!m` has a root mod `p` |
| column `k` represents `m` | `(x)_k − k!m` has an integer root |

A polynomial is **intersective** if it has a root modulo every prime but no
rational root. So an intersective `(x)_k − k!m` is exactly a column that
survives every prime without representing `m` — a census that never
terminates. It is the only logical gap between the method and a guarantee
that it always finishes.

**Answer in one line.** Never for `k ≤ 5` (proved, and `k = 5` — the first
open case — is settled here). Open for `k ≥ 6`, but constrained by three
necessary conditions and unobserved across 180 million values of `c`. And
**operationally the gap is closed anyway**: termination can be certified per
column, cheaply, without settling the general question.

---

## 1. The reduction

**Intersective ⟹ reducible over Q with every irreducible factor of degree ≥ 2.**
If `f` is irreducible of degree ≥ 2 its Galois group acts transitively on the
roots; Jordan's theorem (1872) gives a fixed-point-free element; Chebotarev
turns that into a positive density of primes with no root. A degree-1 factor
would be a rational root.

**`rad(k!) | c` is necessary.** For `p ≤ k`, any `k` consecutive residues
contain a multiple of `p`, so `(x)_k` vanishes *identically as a function* on
`F_p` and `f_c ≡ −c` there. A root then requires `p | c`. Verified for every
prime `p ≤ k` at `k = 5, 7, 12`. This alone cuts the search space by 30× at
`k = 5,6` and 2310× at `k = 12,13`.

---

## 2. `k ≤ 4`: never intersective

- `k ≤ 3`: no rational root forces irreducibility (a cubic with no rational
  root is irreducible), so Jordan applies.
- `k = 4`: the only split with all parts ≥ 2 is `2 + 2`, and **two irreducible
  quadratics can never cover every prime**. Covering needs `χ_p(d₁) = 1` or
  `χ_p(d₂) = 1` for all `p`. If `d₁, d₂, d₁d₂` are all non-squares then
  `Gal(Q(√d₁,√d₂)/Q) = (Z/2)²` and the pattern `(−1,−1)` has density 1/4; if
  `d₁d₂` *is* a square the two fields coincide and `(−1,−1)` has density 1/2.
  Either way a bad prime exists. (Checked on 300 random pairs: a witness prime
  every time.)

This also recovers the standard fact that intersective polynomials need degree
≥ 5, with `(x²+x+1)(x³−2)` attaining it.

---

## 3. `k = 5`: never intersective — the first open case, settled

Degree 5 with no rational root and all factors ≥ 2 leaves only a `2 + 3` split.
Centring on `y = x − 2` gives `(x)_5 = y⁵ − 5y³ + 4y`. Matching

    y⁵ − 5y³ + 4y − c = (y² + ay + b)(y³ − ay² + ey + f)

and eliminating `e, f` leaves one condition on `(a, b)`:

    b² + b(5 − 3a²) + (a² − 1)(a² − 4) = 0,   discriminant  5a⁴ − 10a² + 9.

So a `2+3` split exists exactly when `5a⁴ − 10a² + 9` is a perfect square. Now

    5a⁴ − 10a² + 9 = 5(a² − 1)² + 4,

so with `u = a² − 1` this is the **Pell equation `v² − 5u² = 4`**, whose
solutions are `(v, u) = (L_{2j}, F_{2j})` since `L_n² − 5F_n² = 4(−1)ⁿ`.
Therefore

> `(x)_5 − c` is reducible  ⟺  `a² = F_{2j} + 1` for some `j`
> — a Fibonacci number one less than a perfect square.

`F_n + 1` is a perfect square only for `n = 0, 4, 6` (Cohn; Robbins,
*Fibonacci and Lucas numbers of the forms w²−1, w³±1*). Verified here to
`n ≤ 5000`. Hence `a ∈ {0, ±1, ±2, ±3}` and `c ∈ {0, ±210, ±2160}`.

Brute-force factorisation over `|c| ≤ 6000` confirms the list exactly, and
also turns up `c = ±120, ±720, ±2520` with a `1 + 4` split — those are
`(x)_5` at `x = 5, 6, 7`, i.e. genuine falling-factorial values with a
rational root, so not candidates.

The four real candidates each split as (irreducible quadratic)(irreducible
cubic) and each dies at a small prime:

| `c` | factorisation | killed at |
|---|---|---|
| ±210 | `(y²∓2y+7)(y³±2y²−8y∓30)` | `p = 19` |
| ±2160 | `(y²∓3y+20)(y³±3y²−16y∓108)` | `p = 7` |

**So `(x)_5 − c` is never intersective.** The obstruction turns out to be a
Fibonacci condition — in a project about the Fibonacci family.

---

## 4. `k` even: an extra structural obstruction

The roots `0..k−1` centre to `t = 2x − k + 1` taking the odd values
`±1, ±3, …, ±(k−1)`, so

    (x)_k = 2^(−k) ∏_{j=1}^{k/2} (t² − (2j−1)²),

verified symbolically at `k = 2,4,6,8`. Writing `β_i` for the roots of
`R(u) − 2^k c`, we get `f_c = 2^(−k) ∏_i (t² − β_i)`.

The "flip every sign" element `ε` (trivial on `M = Q(β_i)`, sending
`√β_i ↦ −√β_i`) is a **derangement** on the `k` roots, provided no `β_i = 0`.
And `β_i = 0` would need `c = (−1)^{k/2}((k−1)!!)²/2^k`, which is never an
integer because `(k−1)!!` is odd. By Kummer theory `ε ∈ Gal(L/M)` unless some
**odd-size** subset `S` has `∏_{i∈S} β_i` a square in `M`. Therefore

> `k` even and `f_c` intersective ⟹ an odd square-relation among the `β_i`.

That is precisely what the standard `(x²−a)(x²−b)(x²−ab)` construction
arranges — the three discriminants multiply to a square, an odd relation — and
precisely what cannot be arranged freely here, because the `β_i` are the roots
of one fixed polynomial `R(u) − 2^k c`. This is the sense in which the even
case is harder to break than the odd case.

---

## 5. Search: `k = 5..15`, `|c| ≤ 10⁹`

`scripts/intersective_search.py`. Sieve by `rad(k!) | c`, then by
`c mod p ∈ I_p = {(x)_k mod p}` for the 40 smallest primes above `k` — the
first of which is brutally selective (`|I_p| = 2` at `p = k+1` when prime).

| k | rad(k!) | candidates | survivors | with rational root | **intersective candidates** |
|---|---:|---:|---:|---:|---:|
| 5 | 30 | 66,666,667 | 125 | 123 | **0** |
| 6 | 30 | 66,666,667 | 30 | 30 | **0** |
| 7–10 | 210 | 9,523,809 each | 33/11/13/4 | all | **0** |
| 11–12 | 2310 | 865,801 | 5/2 | all | **0** |
| 13–15 | 30030 | 66,601 | 1/1/1 | all | **0** |

**Every survivor is a genuine falling-factorial value.** (The two `k=5`
survivors that are irreducible rather than rational-rooted are not
intersective either — irreducible means Jordan applies.) Zero intersective
polynomials across ~180 million values of `c`.

---

## 6. What this means for the project — the part that actually closes the gap

Our `c = k!·m` with `m = C(N,K)` has 3.1 million digits at `i=8`. **No search
will ever reach it.** So the general question, even if settled, would not by
itself certify our columns. But it does not need to be:

> `(x)_k − k!m` irreducible over Q
> ⟹ Gal transitive ⟹ (Jordan) a derangement exists
> ⟹ (Chebotarev) a prime with no root exists, density ≥ 1/k
> ⟹ **the census provably terminates for that column.**

and irreducibility over Q has a one-prime certificate: irreducible mod a
single `q` implies irreducible over Q. Since `m mod q` is just Lucas, this
never builds `m`.

`scripts/termination_certificate.py` does it. Results:

- **all 29 columns `k ≤ 30` of `i=8`: certified**, in milliseconds each;
- **`i=9`'s four long-run columns — `k = 11, 29, 40, 45` — all certified**
  (`q = 631, 433, 563, 1549`). These are exactly the columns whose runs of
  8, 7, 4 and 5 would have fired the old run-length trigger.

Cost is `O(k² log q)` per prime with about `k` primes needed (the density of
irreducible reductions is ~`1/k`), so it is practical to a few hundred — which
is exactly the regime where columns are hardest to kill, since small columns
sit at `g/p` up to 0.94. Large columns die at once and never need it.

**So every column we can certify is upgraded from "a killing prime was found"
to "a killing prime had to exist".** That is the upgrade Q14 was posed to
deliver, obtained per column rather than in general.

---

## 7. What remains open, stated precisely

For `k ≥ 6` it is not proved that `(x)_k − c` is never intersective. Any
counterexample must satisfy, simultaneously:

1. `rad(k!) | c`;
2. `f_c` reducible over Q with every irreducible factor of degree ≥ 2;
3. the covering condition — the Galois group is the union of its point
   stabilisers on the roots;
4. for even `k`, an odd square-relation among the `β_i` (§4).

None was found below `10⁹`. The natural next step, if anyone wants the general
theorem, is `k = 6` by the same route as `k = 5`: enumerate the `2+4`, `3+3`
and `2+2+2` splits, reduce each to a Diophantine condition, and check the
finitely many survivors. `k = 5` took one Pell equation; `k = 6` will take
three curves, and the even-`k` obstruction of §4 should kill the `2+2+2` case
outright.
