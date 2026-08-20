# Band I Z-jump spec (i=8)

2026-08-19. Implemented in `scripts/zjump.py`. Bat: `scripts/run-zjump.bat`.
Kernel: `scripts/bandii_kernel.py` (`scan_columns_general`, two-digit `r(p)`).

**This file does not start a job.** The bat refuses if `results/zjump.json`
exists. The Band II ban (do not next-prime from \(k\) through \((k,N/2]\))
is untouched. Hang-guard \(p-k>20000\) is obsolete.

## What this is

The leftover Band I columns for \(i=8\):

| set | columns | why leftover |
|---|---:|---|
| Stage-3 hang-guards | 89195 | six fat Z slabs, gap cap 20000 |
| \(10^6..K\) minus stragglers | 3126621 | never swept |
| **total** | **3215816** | |

Hang-guard runs (inclusive):
`[515813,520186]`, `[607003,615518]`, `{615520}`,
`[687773,700252]`, `[741883,751696]`, `[834629,880328]`,
`[953861,962170]`.

New range: \(k=1000001..4126621\). Stragglers \(4126622..4126646\)
already certified at \(p>N/2\). Family \(K,K+1\) never tested.

Already closed, **not** in this job: exact/modular \(k\le 400\),
Stage 1–2, Stage 3 killed rows, stragglers, Band II.

## Method

For two-digit \(p\), digit-0 / Z-slab primes have \(r(p)=0\) and cannot
kill. PART-lower and NONE are live (\(r(p)\ne 0\)). After \(K\), primes
in \((K,N/2]\) cannot kill (Band II zero theorem). Live again on
\((N/2,d]\).

For each leftover \(k\):

1. Jump to the first live prime \(p>k\) (cell formula, Q3 Z-width).
2. Image test: \(k\) survives iff \(r(p)\in I_{p,k}\), same factorial
   kernel as Band II, with
   \(r(p)=C(\alpha,\beta)C(n_0,k_0)\bmod p\).
3. At most **12** live primes. A column still up at the cap is an anomaly
   only if the size law did not expect it: see the escalation rule below.

## Escalation (supersedes "run length")

Run length is **not** comparable across regimes. Survival per live prime runs
from \(\le 0.16\) in a fat cell to \(0.61\) at small \(k\), so a run of 8
at \(i=9,k=11\) is ordinary (expected count \(\approx 66\)) while a run of 6
in an \(i=8\) fat cell is a thousand times rarer. The old pre-registered
"escalate on a run of 8" fires on the first and is silent on the second.

Escalate on **expected count**, computed per round over the columns that
actually entered it:

$$E_r=\sum_{k\ \text{entering round } r}\frac{|I_{p_k,k}|}{p_k},\qquad
\text{escalate iff } (n_r\ge 1 \text{ and } E_r<10^{-2})
\text{ or } \Pr[X\ge n_r\mid X\sim\mathrm{Poi}(E_r)]<10^{-3}.$$

`scripts/sizelaw.py` implements this and `family_sweep.py` records it per
round. A single phase-wide multiplier is wrong by orders of magnitude: most
Band I columns sit at \(g\approx 9\), the fat cells at \(g\approx 5\times10^5\).

Do not walk Z primes. Do not next-prime-sweep Band II from \(k\).

## Pre-registered checks (pre-flight)

- Hang-guard count \(=89195\).
- \((\alpha,\beta)=(20,7)\) is FULL, \(z_\mathrm{last}=540185\).
  \(k=515813\) first live \(p>540185\), gap \(>20000\).
- \(k=513593\) first live \(=514499\) and dies there.
- \(k=268733\): \(270097/j=589\), \(270121/j=196\), kill \(270131\).
- \(r(p_1)=1275205\) via two-digit table, closed form, and
  `binom_mod_lucas`.

## Stop rules

- Cap 12 live primes. Do not auto-extend.
- Clean = zero leftover \(k\) alive. That closes the Band I remnant.
  Together with earlier Band I kills, stragglers, and Band II, that is
  \(N(C(F_{18}F_{19},F_{16}F_{19}))=6\) exactly. **Not Singmaster.**
- Alive at 12: log \(k,g\), parity, \(b\), hand back. Do not start exact.
  Do not touch Band II again.

## Cost

Most of \(10^6..K\) die at the first live prime with tiny \(g\) (prime
gap). The expensive buckets are under-Z, especially (4,1)\(\to\)(3,1),
same geometry as fat-image. Expect hours, not overnight, at 8 workers.
Checkpoint per \((p,\mathrm{chunk})\).

## Do not

- Do not walk primes inside a Z-slab.
- Do not next-prime Band II from \(k\).
- Do not use \(K=F_{18}F_{17}\).
- Do not mark \(K\) or \(K+1\) impossible.
- Do not build \(m\).
- Do not cache \(F\) across primes.
- Do not call a clean run “Singmaster proved”.
