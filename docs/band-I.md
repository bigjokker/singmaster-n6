# i=8 extra columns (Band I and Band II)

This does **not** prove Singmaster’s conjecture. It proves there is no extra
left-half column for one value: \(m=C(F_{18}F_{19},F_{16}F_{19})\).
The readable statement is [`i8-N6.md`](i8-N6.md).

Fixed objects (do not use \(K=F_{18}F_{17}\); that is \(K+1\)):

\[
N=F_{18}F_{19}=10803704,\quad
K=F_{16}F_{19}=4126647,\quad
d=N-K=6677057,\quad
N/2=5401852,\quad
k_{\max}=5182637.
\]

Band I is \(2\le k<K\). Band II is \(K+2\le k\le k_{\max}\).
Family columns \(K,K+1\) are never marked impossible.
Two-digit primes: \(3287\le p\le K\).
Write \(\alpha=\lfloor N/p\rfloor\), \(\beta=\lfloor K/p\rfloor\),
\(n_0=N-\alpha p\), \(k_0=K-\beta p\), \(g=p-k\),
\(r(p)=C(N,K)\bmod p\).

## Theorem

**Image.** For any prime \(p>k\),

\[
\{C(n,k)\bmod p:n\ge k\}
=\{0\}\cup\{C(k+j,j)\bmod p:0\le j<p-k\}=I_{p,k}.
\]

So \(r(p)\notin I_{p,k}\) is an unconditional certificate that column \(k\)
never represents \(m\).

**No extra column.** Every
\(k\in[2,k_{\max}]\setminus\{K,K+1\}\) has such a prime. Therefore

\[
N\bigl(C(F_{18}F_{19},F_{16}F_{19})\bigr)=6
\]

exactly (the known family, no extra left-half). Not Singmaster.

Coverage:

| \(k\) | job | result |
|---|---|---|
| \(2..400\) | exact \(k\le 300\) + modular \(k\le 400\) | all impossible |
| \(401..10^5\) | Stage 1–2 next-prime | all killed |
| \(10^5..10^6\) | Stage 3 + Z-jump hang-guards | all killed |
| \(10^6..4126621\) | Z-jump | all killed |
| \(4126622..4126646\) | stragglers, \(p>N/2\) | 25/25 killed |
| \(K,K+1\) | family | not tested |
| \(4126649..k_{\max}\) | Band II, \(p>N/2\) | 1,055,989/1,055,989 killed |

## Lemmas

**Digit-0.** For two-digit \(p\), \(C(N,K)\equiv 0\pmod p\) iff \(p(\alpha-\beta)>d\).

**NONE nonvanishing.** If \(p(\alpha-\beta)\le d\) then
\(C(N,K)\equiv C(\alpha,\beta)C(n_0,k_0)\not\equiv 0\pmod p\).
PART-lower (\(p<z_\mathrm{lo}\)) is the same image clause.

**Kill test (no giant \(m\)).** \(I_g=\{(-1)^j\binom{g-1}{j}\bmod p:0\le j<g\}\)
(equivalently \(\binom{k+j}{j}\), **not** \(\binom{k}{j}\)).
A live prime kills iff \(r(p)\notin I_g\).

**Z-width.** Cell
\(P_\mathrm{lo}=\max(\lfloor N/(\alpha+1)\rfloor+1,\lfloor K/(\beta+1)\rfloor+1)\),
\(P_\mathrm{hi}=\min(\lfloor N/\alpha\rfloor,\lfloor K/\beta\rfloor)\),
\(z_\mathrm{lo}=\lfloor d/(\alpha-\beta)\rfloor+1\).
FULL / PART / NONE from where \(z_\mathrm{lo}\) sits.
Live primes are NONE, PART-lower, and \(N/2<p\le d\). Digit-0 and
\(K<p\le N/2\) cannot kill.

**Pascal size of \(I_g\).** Coefficient set of \((1-X)^{g-1}\) over
\(\mathbf{F}_p\). Only systematic collision is reflection \(j\leftrightarrow n-j\):
\(g\) odd \(\Rightarrow \lvert I_g\rvert\approx(g+1)/2\);
\(g\) even \(\Rightarrow \lvert I_g\rvert\approx g\);
then birthday \(p(1-(1-1/p)^M)\). Prime gaps even \(\Rightarrow\) \(g\)-parity
frozen along a chain.

**Band II zero block.** \(K<p\le N/2\Rightarrow p\mid m\). Those primes cannot
kill. The Band II sweep starts at \(p=5401853\), never at \(k\).

**Band II \(r(p)\).** On \(N/2<p\le d\), \(\alpha=1,\beta=0\),
\(r(p)=C(N-p,K)\equiv(-1)^K C(K+\delta-1,\delta-1)\pmod p\) with
\(\delta=2p-N\). At \(p_1=5401853\), \(\delta=2\), \(r=1275205\).

## Census (do not rerun)

**3(b)** on 38 NONE windows, \(10^5\le p\le 10^6\): max run 2, 0 triples,
0 whole-window. [`zeromap-p1e5-1e6.md`](zeromap-p1e5-1e6.md).

**Fat-image**, first 3 primes of fat NONE+PART-lower, \(P_\mathrm{hi}>10^6\),
\(g_{\max}>25000\): 369 triples, hunt cap 3.
[`../results/fat_image_hunt.json`](../results/fat_image_hunt.json).
Two-to-three is false as a lemma. Size law predicts the counts
(\(+1.4\sigma\) pooled).

**Walk-369:** those 369 past prime 3. Pre-registered 44/6/0.8/0.1 vs
measured 42/4/1/0, `max_run=6`. Record \(k=2227205\), kill 2701099.

*Run length depends on the prime-selection rule, so never compare two.*
\(i=7,k=487\) has a Z-jump run of 8 (first live prime 1259) and a
next-prime run of 0 (dies at \(p=491\)). Both are correct; they are
different statistics. Compare expected counts instead — see
[`zjump-spec.md`](zjump-spec.md), Escalation.
[`../results/walk_369.json`](../results/walk_369.json).

**Band II** \(p>N/2\), 1,055,989 columns, cap 14, factorial kernel, 1235 s:
died at prime 8. Pass 1/4/6 inside Poisson bands. Last live \(k=4155257\).
[`../results/bandii_sweep.json`](../results/bandii_sweep.json).
Spec: [`bandii-spec.md`](bandii-spec.md).

**Z-jump** remnant, 3,215,816 columns (89,195 hang-guards +
\(10^6..4126621\)), cap 12, 3102 s: 0 anomalies, died at live-prime 7.
Rounds 4–6 are the walk-369 tail (same 42/4/1, same \(k=2227205\)).
[`../results/zjump.json`](../results/zjump.json).
Spec: [`zjump-spec.md`](zjump-spec.md).

**Stragglers** \(k=4126622..4126646\): 25/25 killed at \(p>N/2\).
[`../results/stragglers_nearK.json`](../results/stragglers_nearK.json).

**Stage 3** hang-guard \(p-k>20000\) is obsolete (Z-jump). Summary:
[`../results/nextprime_i8_k100001-1000000_summary.json`](../results/nextprime_i8_k100001-1000000_summary.json).

## Heuristic (closed as a search)

Long \(q(k)-k\) is Z-slab plus a short image tail (max 6 at \(\rho=0.176\)).
\(G(k)=\lceil(\ln k)^2\rceil\) is false. Independence-after-size is the
right model; the Q1 character sum is the wrong shape.

Do not start Stage 4 until-kill, Band II next-prime from \(k\), exact \(i=10\),
or nearby \(10^9..2\cdot 10^9\).
