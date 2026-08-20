# The 2022 interior theorem, and why the census is still the tool

Matomäki–Radziwiłł–Shao–Tao–Teräväinen, *Singmaster's conjecture in the
interior of Pascal's triangle*, QJMath **73** (2022) 1137–1177
([arXiv:2106.03335](https://arxiv.org/abs/2106.03335)).

**Theorem 1.3.** For `0 < ε < 1` and `t` sufficiently large depending on `ε`,
there are at most two solutions of `C(n,m) = t` with
`exp(log^{2/3+ε} n) ≤ m ≤ n/2`.

Three conventions decide whether it touches this project:

- The largeness hypothesis is on **`t`**, not on `n`. For us `t = C(N,K)` is
  enormous — the favourable direction.
- The region condition is **per solution, against that solution's own `n`**. A
  hypothetical third representation `(n₃,k₃)` with small `k₃` has `n₃ ≫ N`, so
  it is tested against `n₃`. Reading the condition against `N` is the mistake
  that makes the theorem look applicable.
- The family's two solutions are `(N,K)` and `(N−1,K+1)` — the sibling sits in
  row `N−1`. Both must be tested separately. The engine's trivial pair sits at
  `m = 1` and the theorem never sees it.

## It does not apply here

The display inequality is satisfied comfortably: at `i=8`,
`exp((log N)^{2/3}) ≈ 603` against `K = 4126647`, and both `K` and `K+1` are
below `N/2`. The same holds at `i=9` for any fixed `ε ≲ 0.33`. **That is
necessary and cheap, and it is not the question.**

The question is `t₀(ε)`. Remark 1.7 states the bounds are effective but that no
attempt was made to optimise them and they *"will likely be too large to be of
use in numerical verification."* No number appears in the paper. From the one
explicit tolerance in §3 — `δ = log^{−1/1000} P` against a contradiction margin
of `1/72` — the argument cannot close before roughly

    log N > (72C)^{1000/(2/3+ε/2)}

which with an absolute-floor `C = 1` is `log N > 10^2274` at `ε = 0.3`, rising
to `10^2784` as `ε → 0`. Here `log N₈ = 16.2` and `log N₉ = 18.1`. No plausible
sharpening of a hidden constant closes that gap, and no effective `i` can be
named without redoing §§3–5 with every constant tracked — a paper-sized job the
authors deliberately did not do.

## What that means for the tool

**The census is the only effective instrument in this range.** The alternative
route — Remark 1.5's finiteness for fixed `2 ≤ m < m′` — rests on Siegel's
theorem and is completely ineffective, and it lives in exactly the small-`k`
regime the census occupies.

Worth registering the counterfactual: *if* the theorem did apply at `i=8`, every
third representation would have to satisfy `k < exp((log n)^{2/3+ε})` for its
own row, giving a cutoff `k*` solving `k*(log k*)^{1/(2/3+ε)} = log t`. At
`log t₈ = 7.18×10⁶` that is `k* ≈ 1.7×10⁵` (ε→0) to `5.0×10⁵` (ε=0.3) — already
inside Stages 1–3, which would have made Band II and the `10⁶..K` remnant
redundant. It does not apply, so those sweeps are what carry the result.

Nor does 2017 cover the edge: at `k* ≈ 1.7×10⁵` the relevant rows are
`n ~ 10^23`, far past BBW's `n ≤ 10⁶`, while `t ~ 10^{3120255}` is far past
their `t ≤ 10^60`.

## Consistency with the code

Lemma A of the Z-jump notes is what `bandii_kernel` implements, and the two
places it constrains the code are both enforced:

- `r(p) = 0` cannot kill — `scan_ks` refuses such a prime rather than reporting
  every column killed, and `witness.check_witness` rejects the certificate.
- The digit-0 criterion `r(p)=0 ⟺ p(α−β) > d` is valid only for `p > √N`;
  below that `scan_columns_general` falls back to full Lucas, and
  `witness.r_two_digit_pure` returns `None` rather than a false `r`.

Lemma B's warning — *"observed maxima are search caps, not maxima … a column
still alive at the cap is an anomaly to be reported, never a column to be
declared killed"* — is enforced by `sizelaw.RoundLedger` (escalation on expected
count, not run length) and by `family_sweep` refusing to certify when any column
was untestable.
