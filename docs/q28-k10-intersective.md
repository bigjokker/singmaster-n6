# Q28 — `(x)_10 − c` is never intersective (one genus-1 gap)

Worked out 2026-08-20; the gap made precise 2026-08-23 (§3); the Chebotarev
door tried the same day (§3b: half the branch closes, the rest is a locus on
a genus-3 cover). Verified by `scripts/k10_intersective.py` and
`scripts/k10_chebotarev.py`; artifact `results/k10_intersective.json`; pinned
by `scripts/test_k10_intersective.py` and `scripts/test_k10_chebotarev.py`.

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

## 3. The gap — made precise (2026-08-23)

**One elliptic integral-point computation**, for the (2,3) branch, and it is
**still open**: not proved here. This section records everything short of that
computation, so that whoever has the tool can finish it in one call — and
says exactly which tool, because it is not the one a Sage user would reach for.

**The curve is an elliptic curve over Q.** The quartic has rational points
(`a = −250, −730`, and the ten `c = 0` points), so it is an elliptic curve,
not merely a genus-1 curve. From the classical invariants of the binary
quartic, `I = 155713536` and `J = −4439704338432`, scaled by `48`, its
Jacobian is

\[
E:\ Y^2 = X^3 - 792X + 9801 = (X+33)(X^2-33X+297).
\]

**Torsion.** `E(Q)_tors = Z/2 = ⟨(−33, 0)⟩`: the quadratic factor has
discriminant `−99`, and the order of the torsion divides
`gcd_p #E(F_p) = 2` over the good primes `p < 200`.

**Rank: exactly 2**, by descent via the 2-isogeny with kernel `(−33, 0)`.
Shifting `X → X − 33` gives `E: Y² = X³ − 99X² + 2475X` and the isogenous
`E′: Y² = X³ + 198X² − 99X`. Of the sixteen squarefree `d | 2475`, the eight
positive ones are all realised by explicit points (e.g. `(−30, 81)` and
`(12, 45)` on `E` are the classes `d = 3` and `d = 5`) and the eight negative
ones have no real points; of the eight squarefree `d | −99` on `E′`, only
`d = 1` and `d = −11` (the torsion class) have points, and the other six have
**no 3-adic points** — no solution modulo 9 (for `±3, ±33`) or modulo 27 (for
`−1, 11`) with `(w, z)` not both divisible by 3. So
`2^r = |im α| · |im α′| / 4 = 8 · 2 / 4` and `r = 2` — the `4` because
`E(Q)[2] = E′(Q)[2] = Z/2` (discriminants `−99` and `39600` are not squares),
and `im α`, `im α′` are the realised classes `E(Q)/φ̂(E′(Q))` and
`E′(Q)/φ(E(Q))`. Every `φ`-Selmer and `φ̂`-Selmer class is decided — a point,
or no points over `R` or `Q₃` — so `Ш(E)[φ] = Ш(E′)[φ̂] = 0`, hence
`Ш(E)[2] = 0` and the Selmer bound *is* the rank: **`E(Q) ≅ Z/2 × Z²`.**

**What rank 2 means for the gap.** Rank 0 would have made the rational points
finite and the integral-point list elementary; that door is closed. With rank
2 the complete list is an elliptic-logarithm computation — Mordell–Weil
generators, elliptic logarithms, a Baker-type bound, LLL reduction — which is
the quartic case of Stroeker–Tzanakis: N. Tzanakis, *Solving elliptic
diophantine equations by estimating linear forms in elliptic logarithms. The
case of quartic equations*, Acta Arith. 75 (1996), 165–190. Siegel alone gives
finiteness, not a bound a search can meet, and the leading coefficient `5` is
not a square, so Runge's method does not apply either (the two points at
infinity are conjugate over `Q(√5)`). A bigger search is not a proof.

**The complete small list, which is what the tool would have to confirm.**
Every integral point with `|a| ≤ 10⁷` — fifteen values of `a`:

\[
-730,\ -250,\ -130,\ -106,\ -90,\ -82,\ -74,\ -58,\ -50,\ -34,\ -26,\ -10,\ 46,\ 54,\ 158,
\]

the largest `|a|` being `730`. So the `|a| ≤ 6000` search missed nothing
below `10⁷` — and nothing here says anything about `|a| > 10⁷`. Each
point gives two integral `(b, c)` pairs; of the thirty, ten give `c = 0`,
eighteen give a `c` not divisible by `rad(10!)`, and two are the candidates of
§2, both killed. (The necessary condition is in fact `10! | c`, since
`(x)_10 ≡ 0 (mod 10!)` for every integer `x`; both candidates satisfy it,
`c = 384540 · 10!` and `c = 616 · 10!`, so the script's weaker `rad(10!)`
filter hides nothing.)

**The missing tool.** The routine that does this is Magma's
`IntegralQuarticPoints` (linear forms in elliptic logarithms, after
Stroeker–Tzanakis 1994 and Tzanakis 1996); since neither the leading nor the
constant coefficient is a square, it is the two-argument form that takes a
rational point:

```text
> IntegralQuarticPoints([5, 1320, 126456, 5102240, 72824400], [-250, 74880]);
```

(Use a current Magma: a V2.24-3 release note records a p-adic precision fix
in this routine for high-rank curves. The handbook documents the
one-argument form as returning points modulo negation; compare `a`-values
whichever sign convention the output uses. Neither was tested here.) Sage
has no equivalent:
`E.integral_points()` works only on a Weierstrass model, and the integral
points of `E` are not the integral points of the quartic — the birational
map has denominators. PARI/GP has no integral-point routine at all
(`hyperellratpoints` is a height-bounded search). None of Magma, Sage or
PARI is installed on this machine, there is no WSL, and no Windows-native
wheel provides the routine; nothing was installed in this pass. If that
Magma call returns exactly the fifteen `a`-values above **and** Magma
reports the Mordell–Weil basis it used as unconditional and saturated
(`MordellWeilShaInformation` / `Saturation` on the Jacobian), Q28 is
**proved**. If it returns anything else, that point is the result — write
it down.

The rank Magma finds must be 2 — that much is certified here, independently
of Magma. Saturation of its basis is Magma's certificate, not ours: the
elliptic-logarithm method lists the integral points in the subgroup
generated by the basis it is given, so an unsaturated basis can miss points.

---

## 3b. The Chebotarev door (tried 2026-08-23): half the branch closes, the rest is Q27's trap

Could a direct Chebotarev kill make the integral-point list unnecessary? For
`p` odd, `f_c` dies at `p` iff `h(t) = g(t²)` is rootless mod `p`, i.e. iff
Frobenius fixes none of `h`'s ten roots `±√β_j, ±√θ_i` — and Chebotarev turns
"the Galois group contains a fixed-point-free element" into "a killing prime
exists". Existence is all a non-intersectivity proof needs; no effectivity
anywhere. The answer is a trichotomy (each case exercised on a synthetic
polynomial in `scripts/k10_chebotarev.py`):

1. **`D_q·D_k` not a square** (the two factors' quadratic subfields differ;
   this includes `Gal(k) = A₃`). Then `Gal ∋ (transposition, 3-cycle)` and `g`
   itself is rootless mod a positive density of `p`. Killed. *Both recorded
   candidates are here* — `D_q·D_k > 0` but not square — consistent with
   their deaths at 11 and 13.
2. **Shared subfield, `k(t²)` irreducible.** `(transposition, 3-cycle)` is
   now excluded — the 3-cycles are even, exactly Q27's derangement trap —
   but the transposition class still kills through the quadratic character:
   `k` has one root `θ₁` mod such `p`, and either `θ₁` is not a square in
   the splitting field (a Kummer character flips `√θ₁`) or `√θ₁ ∈ L∖Q(θ₁)`
   and Frobenius flips it for free. The only escape is `√θ₁ ∈ Q(θ₁)`, which
   is precisely `k(t²)` reducible. Killed.
3. **`k(t²) = m(t)·(−m(−t))` reducible** (forces `f = −γ²`). Only the
   3-cycle class can still kill (the `θ_i` are squares in `L_m`, so the
   other classes always see a root), and it does iff `q(t²)` is irreducible
   over `Q`: the one obstruction, `β₁ ∈ (K*)²` with `K = Q(√D_q)` — an
   odd-degree extension adds no square roots, so square in `L_m` equals
   square in `K` — is *equivalent* to `q(t²)` splitting into two rational
   quadratics. Killed unless `q(t²)` also splits.

**The no-kill locus.** All three fail exactly when, over `Z`,

\[
R(t^2) - 1024c = (t^2{+}\mu t{+}\nu)(t^2{-}\mu t{+}\nu)(t^3{+}\alpha t^2{+}\beta t{+}\gamma)(t^3{-}\alpha t^2{+}\beta t{-}\gamma)
\]

with `(μ²−4ν)·disc(m)` a square (without that coincidence, `Gal = S₃ × C₂`
and (3-cycle, flip) kills). Two identities tie the cases together:
`D_q = μ²(μ²−4ν)` and `D_k = disc(m)·(γ−αβ)²`, so the coincidence is
*exactly* `D_q·D_k` square — the failure of case 1 — and the four cases
partition the branch with no gap. On the locus **no Frobenius class kills**:
every element of `Gal = S₃` fixes a root of `h`, so `h` has a root mod
**every** prime (at a ramified `p`, a decomposition-group lift of Frobenius
still fixes a root — an algebraic integer whose reduction lands in `F_p`) —
this is where an intersective `c` would have to live (it would still owe
roots mod prime powers), and it is the k=10 reappearance of Q27's mechanism
through the cubic `m`. Witness that
the obstruction is real: `(t²+t+44)(t²−t+44)(t³+t²+2t+3)(t³−t²+2t−3)` (both
discriminants `−175`) has a root modulo every prime tested to 2000.
Equivalently, the locus is `R(u) − 1024c = u·P(u)² − Q(u)²` with
`P = u²+p₃u+p₁` monic, `Q = p₄u²+p₂u+p₀`; eliminating gives one plane curve
`F(p₄,p₂) = 0` with `p₄` odd forced, and `c = (p₀² − 945²)/1024` where
`945² = −R(0)` — so `c = 0 ⟺ p₀ = ±945`.

**What this proves, unconditionally.**

- **Every `(2,3)` value with `c < 0` is killed.** For `10! | c` the quintic
  has exactly one real root and it is negative (`R` is strictly increasing on
  `(−∞,1]`, its largest critical value magnitude `43,930,542.03 < 1024·10!`,
  and `R(0) = −893025 > 1024c`); but `m` is a real cubic whose real root
  `t₀` needs `t₀² =` that root `≥ 0`. So the deep factorisation is
  impossible and cases 1–3 kill everything. Half the branch closes with no
  integral-point list at all. (For `c < 0` with `10! ∤ c` — at any size — a
  power of a prime `≤ 7` already kills, since `(x)_10 ≡ 0 mod p^{v_p(10!)}`
  for every integer `x`.)
- **For `c > 0`, killed except on the no-kill locus**, and that locus forces
  `b = ν²`: its points lie on the double cover `ν⁴ − s(a)ν² + q_C(a) = 0` of
  the genus-1 curve `C`, ramified exactly at the four simple zeros of `b` —
  **genus 3** by Riemann–Hurwitz, before `f = −γ²` and the coincidence thin
  it further. Faltings gives finiteness, nothing gives the list: the same
  wall as Q27, one member later.

**Searches on the locus** (evidence, not proof): every point of `C` with
`b = ν²`, `|ν| ≤ 10⁶` (i.e. `b ≤ 10¹²`), is either one of the ten `c = 0`
points — which satisfy the whole deep shape, as they must, since `g` splits
completely there — or `a = 46, 158`, where `−f` is not a square and the
locus misses them. On `F(p₄,p₂) = 0`, every integral point with odd
`p₄ ≤ 2·10⁵` that lifts to the locus (`p₀ ∈ Z`) has `p₀ = ±945`, i.e.
`c = 0`; four spurious points on `p₂ = −p₄³` (`p₄ ∈ {1,5,7,9}`) have
`p₀ ∉ Z` — elimination artifacts. Any `10!`-divisible deep `c`
must have `315 | p₀`, `p₀` odd, and `p₀ ≡ ±945 (mod 2¹⁷)`.

**Outcome: Q28 stays BLOCKED.** The missing statement narrows to *"the
no-kill locus has no integer point with `c ≠ 0`"* — implied by the §3 Magma
call, since deep points are integral points of the rank-2 quartic with `b` a
perfect square. One correction to the flag that prompted this section (the
previous commit's report): a shared quadratic subfield alone does **not**
block the kill — case 2's character argument survives it. The block needs
the full `(2,2,3,3)` factorisation, whose coincidence condition is, by the
identities above, that same shared-subfield condition resurfacing at the
bottom of the tree.

## 4. Ladder, corrected

| `k` | status | gap |
|---|---|---|
| ≤6 | **proved** | — |
| 7 | modulo effective Thue | two curves |
| 8 | modulo effective Siegel | genus 3 |
| **10** | modulo one elliptic computation (Jacobian of rank 2; needs Magma's `IntegralQuarticPoints`, not done). §3b: `(2,3)` with `c < 0` closed unconditionally; `c > 0` reduced to a no-kill locus on a genus-3 cover | **genus 1** |

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
