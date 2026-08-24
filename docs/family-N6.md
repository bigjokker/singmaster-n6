# \(N=6\) for the Fibonacci family, \(i=2\) through \(i=9\)

This note records a theorem about eight binomial values. It does **not**
prove Singmaster's conjecture.

[`i8-N6.md`](i8-N6.md) is the detailed \(i=8\) case and is not superseded;
this note states the family claim, the lemmas common to every member, and
the exhaustion and provenance for each. Where the two overlap, the lemmas
are the same lemmas.

## 1. Statement

Let \(F_n\) be the Fibonacci sequence with \(F_1=F_2=1\). For \(i\ge1\) put

\[
N_i=F_{2i+2}F_{2i+3},\qquad
K_i=F_{2i}F_{2i+3},\qquad
m_i=C(N_i,K_i).
\]

(The complementary family column is \(K_i+1=F_{2i+1}F_{2i+2}\). Never
confuse \(K\) with \(F_{2i+2}F_{2i+1}\); that is \(K+1\).)

**Family identity.** \(C(N,K)=C(N-1,K+1)\). Equivalently, and checkable
without constructing \(m\),

\[
N(K+1)=(N-K)(N-K-1),
\]

an identity of the family for every \(i\ge1\); it is verified here by exact
integer arithmetic for \(i=1,\dots,9\). So each \(m_i\) has two off-centre
left-half representations, \((N,K)\) and \((N-1,K+1)\).

Write \(k_{\max}\) for the largest \(k\) with \(C(2k,k)\le m\).

> **Theorem A (no extra column).** For each \(i=2,\dots,9\) there is no
> pair \((n,k)\) with \(2\le k\le n/2\), \(k\notin\{K_i,K_i+1\}\) and
> \(C(n,k)=m_i\).
>
> **Theorem B (multiplicity).** For each \(i=2,\dots,9\),
> \(N(m_i)=6\).

B follows from A together with two facts proved in §2: that \(C(n,k)\) is
strictly increasing in \(n\) for \(n\ge k\ge1\), so each of the columns
\(K\) and \(K+1\) contains **at most one** \(n\) with \(C(n,k)=m\) — and it
contains exactly one, by the family identity — and that there is no central
solution \(C(2k,k)=m\). The six appearances are then exactly: the two
trivial \(C(m,1)=C(m,m-1)=m\), and the mirror pair of each of the two
left-half representations.

The objects, recomputed from \(N\) and \(K\) rather than read back from any
file under audit (`bandii_kernel.make_fam`, `kmax_of`):

| \(i\) | \(N\) | \(K\) | \(\log_{10}m\) | \(k_{\max}\) | extra columns | certificates |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 104 | 39 | 28.7869 | 49 | 46 | 46 |
| 3 | 714 | 272 | 204.5487 | 342 | 339 | 339 |
| 4 | 4,895 | 1,869 | 1,411.6618 | 2,347 | 2,344 | 2,344 |
| 5 | 33,552 | 12,815 | 9,687.7814 | 16,094 | 16,091 | 16,091 |
| 6 | 229,970 | 87,840 | 66,415.5955 | 110,318 | 110,315 | 110,315 |
| 7 | 1,576,239 | 602,069 | 455,236.2642 | 756,136 | 756,133 | 756,133 |
| 8 | 10,803,704 | 4,126,647 | 3,120,255.2212 | 5,182,637 | 5,182,634 | 5,182,634 |
| 9 | 74,049,690 | 28,284,464 | 21,386,569.3411 | 35,522,329 | 35,522,326 | 35,522,326 |

"Extra columns" is \(\bigl|[2,k_{\max}]\setminus\{K,K+1\}\bigr|=k_{\max}-3\)
(both family columns lie in range for every member). "Certificates" is the
row count of `results/i{i}_witness.npz`. That the two agree, for every
member, is the coverage claim of §6; it is machine-checked by
`coverage_ledger.py`, not eyeballed. **41,590,228 columns in total.**

Coverage is not the same as having re-checked every certificate. §6 gives
the exact fraction of each table that has been independently verified; it
is 100% for \(i=2..6\) and a sample for the three largest. A reader who
wants the strongest available statement per member should read §6 before
quoting this table.

\(i=1\) is \(C(15,5)=3003\), which has \(N=8\) — the sporadic value, not a
counterexample to anything here, and excluded from the theorem. Its third
representation \(C(78,2)=3003\) is recorded in
`results/fibonacci_i1-7.json` (`status: known_3003`).

## 2. What a certificate is, and which \(k\) can host a solution

**Lemma 1 (monotonicity).** For fixed \(k\ge1\), \(C(n,k)\) is strictly
increasing in \(n\) on \(n\ge k\). *Proof.* \(C(n+1,k)-C(n,k)=C(n,k-1)>0\)
for \(n\ge k-1\). \(\square\) Consequently each column \(k\) contains at
most one \(n\) with \(C(n,k)=m\): a column is killed or hosts exactly one
solution, never two.

**Lemma 2 (complete image).** Let \(p\) be prime and \(k<p\). Then

\[
\{C(n,k)\bmod p:n\ge k\}
=\{0\}\cup\bigl\{C(k+j,j)\bmod p:0\le j<p-k\bigr\}.
\]

Write \(I_{p,k}\) for the right-hand side. With \(g=p-k\),

\[
C(k+j,j)\equiv(-1)^j C(g-1,j)\pmod p\qquad(0\le j<g).
\]

(The identity uses \(C(k+j,j)\), not \(C(k,j)\).)

*Proof sketch.* Lucas on the single digit of \(k\): since \(k<p\),
\(C(n,k)\equiv C(n_0,k)\pmod p\) with \(n_0=n\bmod p\). As \(n_0\) runs
through \(\{k,\dots,p-1\}\) one gets \(C(k+j,k)=C(k+j,j)\) for
\(0\le j<g\), together with \(0\) when \(n_0<k\). For \(n\ge p\) the same
residues repeat. \(\square\)

**Corollary (certificate).** Write \(r(p)=m\bmod p\). If \(p>k\) and
\(r(p)\notin I_{p,k}\) then no \(n\) satisfies \(C(n,k)=m\). This is a
statement about the whole column, not a sample, and it is unconditional.

Three hypotheses are load-bearing, and all three are checked per
certificate rather than assumed: \(p\) is **prime** (Lucas needs it, and so
does the claim that \(I_{p,k}\) is the *complete* image); \(p>k\); and
\(r(p)\neq0\), because \(0\in I_{p,k}\) always — \(0=C(x,k)\) for every
\(x<k\) — so a prime dividing \(m\) certifies nothing.

**Which columns need a certificate.** A left-half solution has \(k\le n/2\),
hence \(n\ge2k\) and, by Lemma 1, \(C(n,k)\ge C(2k,k)\). For
\(k>k_{\max}\), \(C(2k,k)>m\), so no left-half solution exists there. The
bracket

\[
C(2k_{\max},k_{\max})\ \le\ m\ <\ C(2k_{\max}+2,k_{\max}+1)
\]

holds exactly for every \(i=2,\dots,9\) (exact integer arithmetic; 2.6 s at
\(i=9\)). So extra left-half columns can occur only for

\[
2\le k\le k_{\max},\qquad k\notin\{K,K+1\},
\]

and \(k=1\) is the trivial pair, outside this range.

**No central solution.** A central representation is \(C(2k,k)=m\), the case
\(n=2k\). By cases, exhaustively: \(k=1\) gives \(C(2,1)=2\neq m\); for
\(k\in[2,k_{\max}]\setminus\{K,K+1\}\) the column certificate rules out
*every* \(n\), hence \(n=2k\) with the rest; for \(k>k_{\max}\),
\(C(2k,k)>m\) by definition of \(k_{\max}\); and the two family columns fall
to Lemma 1, since \(2K<N\) gives \(C(2K,K)<C(N,K)=m\) and \(2(K+1)<N-1\)
gives \(C(2K+2,K+1)<C(N-1,K+1)=m\). Both inequalities hold for every
\(i=2,\dots,9\) (\(K/N\to\varphi^{-2}\approx0.382\)). That is every case,
so no central solution exists — and with Lemma 1 fixing the family columns
at one solution each, Theorem B follows from Theorem A.

## 3. Where a prime can kill

Fix a member; write \(d=N-K\), \(\alpha=\lfloor N/p\rfloor\),
\(\beta=\lfloor K/p\rfloor\). Numbers in this section are \(i=8\)'s, as the
worked example.

**Two-digit range.** For \(\sqrt N<p\le K\) the base-\(p\) expansions of
\(N\) and \(K\) have two digits and Lucas gives
\(r(p)\equiv C(\alpha,\beta)C(n_0,k_0)\) with \(n_0=N-\alpha p\),
\(k_0=K-\beta p\). At \(i=8\) that range is \(3287\le p\le4126647\).

**Digit-0.** \(m\equiv0\pmod p\) iff \(p(\alpha-\beta)>d\). Such a prime
cannot kill, since \(0\in I_{p,k}\).

**Z-width.** For a cell with \(\beta\ge1\), the prime range of
\((\alpha,\beta)\) is
\(P_{\mathrm{lo}}=\max(\lfloor N/(\alpha+1)\rfloor+1,\lfloor K/(\beta+1)\rfloor+1)\),
\(P_{\mathrm{hi}}=\min(\lfloor N/\alpha\rfloor,\lfloor K/\beta\rfloor)\),
and \(z_{\mathrm{lo}}=\lfloor d/(\alpha-\beta)\rfloor+1\) splits it into a
forced-zero part and a live part (FULL / PART / NONE by where
\(z_{\mathrm{lo}}\) falls). The \(\beta=0\) regime — \(p>K\) — is Band II
below, handled separately.

**Band II zero block.** For \(K<p\le N/2\), every prime divides \(m\), so
none can kill. *Why:* \(p>K\) gives \(\beta=0\) and \(k_0=K\); \(p\le N/2\)
gives \(\alpha\ge2\), and \(3K>N\) for every member (\(N/K\to\varphi^2
\approx2.618<3\)) forces \(\alpha=2\) exactly, so \(n_0=N-2p<N-2K<K=k_0\)
and the Lucas digit \(C(n_0,k_0)\) vanishes. This is why Band II is never
next-prime-swept from \(k\): the whole block is provably useless.

**Band II live window.** On \(N/2<p\le d\) one has \(\alpha=1\),
\(\beta=0\), and \(r(p)=C(N-p,K)\not\equiv0\). With \(\delta=2p-N\),

\[
r(p)\equiv(-1)^K C(K+\delta-1,\delta-1)\pmod p ,
\]

whose lower index \(\delta-1\) is tiny at the bottom of the window. At
\(i=8\), \(p_1=5401853\), \(\delta=2\), \(K\) odd, so
\(r(p_1)=-(K+1)=1275205\). At \(i=9\) the first Band II prime is
\(37024873\).

**The two-digit formula is false below \(\sqrt N\).** If \(p\le\sqrt N\)
then \(N\) has three or more base-\(p\) digits and
\(C(\alpha,\beta)C(n_0,k_0)\) drops one, so it is not \(r(p)\). Full Lucas
is correct at any \(p\) and is what the verifier and the fills of §7 use.

**Three reasons the sweep does not test every prime above every column,**
and they are not equivalent. Only the first is a proof that nothing was
missed:

1. *provably cannot kill* — digit-0 primes and the Band II zero block are
   skipped analytically, because \(0\in I_{p,k}\) there;
2. *not tested, by construction* — the cell geometry starts at
   \(\lfloor\sqrt N\rfloor+1\) (at \(i=9\), never below 8606), so
   sub-\(\sqrt N\) primes are never scanned;
3. *not tested, by budget* — each column is tested against a capped number
   of live primes.

Cases 2 and 3 leave columns without a certificate rather than with a false
one. That is exactly what happened at \(i=8\) and \(i=9\), and §7 says what
closed those columns instead.

## 4. The size law (calibration, not evidence)

**The size law is not part of the proof.** Theorems A and B rest on §2 and
§6: every extra column carries a certificate, and the certificates are
checkable. A column that survived its cap would simply have no certificate,
and the coverage check would report it missing — which is what makes a
*bounded* search safe to run in the first place. The size law decides where
to spend compute and whether an outcome is surprising. It never supplies a
kill. Read this section as calibration.

**Image size.** \(I_{p,k}\) is the value set of \(f_k(x)=(x)_k/k!\) over
\(\mathbf F_p\). The involution \((k-1-x)_k=(-1)^k(x)_k\) gives an exact
pre-image count, hence a proved upper bound:

\[
|I_{p,k}|\le
\begin{cases}
g+1, & k\text{ odd (no fold)},\\[2pt]
\lceil g/2\rceil+1, & k\text{ even (exact 2-to-1 fold)}.
\end{cases}
\]

Both are theorems, uniform in \(p\) and \(k\). The model then subtracts
birthday collisions, \(|I|\approx p\bigl(1-(1-1/p)^M\bigr)\) with \(M=g+1\)
or \(\lceil g/2\rceil+1\) — that step is a **heuristic**, an assumption that
the coefficients behave like random residues, and it is the only part of
§4 that is not proved. It has no fitted parameters. The parity that matters
is \(k\)'s, not \(g\)'s: they agree for odd \(p\), but \(k\)-parity is a
property of the polynomial while \(g\) changes with every prime. Survival
per live prime is \(|I_{p,k}|/p\): at \(i=8\)'s first Band II prime that is
\(0.2103\) at \(k=K+2\) and \(0.0398\) at \(k=k_{\max}\).

**Two different expectations, and they must not be confused.** The
*pre-registration* is the forward model, computed before the run from
\(k\)-range and prime list alone (`results/bandii_sweep.json`,
`preregister`). The *per-round* \(E_r\) recorded in
`results/i{i}_sweep.json` is a one-step-ahead conditional: round \(r\)'s
expectation is computed over the columns that actually entered it, i.e.
from the previous round's observed survivors. Both are useful; only the
first is a prediction of the whole curve. For \(i=8\) Band II:

| pass | pre-registered | per-round \(E_r\) | observed |
|---:|---:|---:|---:|
| 1 | 102,600 | 102,563 | 102,754 |
| 2 | 12,600 | 12,624 | 12,478 |
| 3 | 1,816 | 1,805 | 1,782 |
| 4 | 289.9 | 284.3 | 293 |
| 5 | 49.4 | 49.52 | 48 |
| 6 | 8.78 | 8.90 | 7 |
| 7 | 1.61 | 1.394 | 1 |
| 8 | 0.300 | 0.2061 | 0 |

Observed even-\(g\) runs \(0.657\to1.000\) against a pre-registered
\(0.658\to0.994\), and the mean surviving \(k\) walks *down* toward \(K\) as
predicted — a statement about *which* columns survive, not how many.

\(i=9\) is the same shape one member up. Its Band II forward model
regenerates in seconds —
`sizelaw.expected_alive(*rep["k_bii"], rep["primes_bii"])` gives
702,981 / 86,359 / 12,444 / 1,987 / 338.6 / 60.2 / 11.0 / 2.05 / 0.388 for
passes 1..9 — against measured
703,000 / 86,212 / 12,411 / 1,946 / 318 / 47 / 6 / 1 / 0, death at pass 9,
`n_bii_alive = 0`. Two caveats on that measured row: it comes from the run's
jsonl events, tabulated in `Singmaster-REFRESHER.txt`, because
`results/i9_sweep.json` has no Band II phase record (§7); and that jsonl is
219 MB and gitignored, so the row is not reproducible from a clone. The
count `n_bii_alive = 0` *is* in the json.

**Escalation.** Run length is not the criterion: survival per live prime
spans two orders of magnitude across regimes, so a run of 8 at small \(k\)
and a run of 6 in a fat cell are not comparable. Two tests, both on
expected count, both at the same thresholds (\(10^{-2}\), and \(10^{-3}\)
on the Poisson upper tail):

* *per round* — \(E_r\) over the columns entering round \(r\) with the
  primes they were actually tested against;
* *per column alive at the cap* — \(\Lambda=\prod_j|I_{p_j,k}|/p_j\) over
  the primes that column actually faced, judged as
  \(\text{peers}\times\Lambda\) against the columns that entered the phase.

The second exists because the first is blind to its motivating case: once a
single column enters a round, \(E_r\) is just that column's own next-prime
survival (0.04–0.21 in Band II), never below the threshold, so a lone column
surviving every pass would read "ordinary" every time. Calibration on real
columns:

| case | run | \(\Lambda\) | peers × \(\Lambda\) | verdict |
|---|---:|---:|---:|---|
| \(i=8\) Band II \(k=4{,}126{,}649\), all 14 passes | 14 | \(3.31\times10^{-10}\) | \(3.49\times10^{-4}\) | **would escalate** |
| \(i=9\) \(k=11\) | 8 | \(2.19\times10^{-2}\) | 1.73 | ordinary |
| \(i=9\) \(k=87\) at the cap | 12 | \(3.8\times10^{-3}\) | \(1.07\times10^{5}\) | ordinary |
| \(i=9\) \(k=1281\) at the cap | 12 | \(1.28\times10^{-3}\) | \(3.61\times10^{4}\) | ordinary |

The first row is hypothetical: no Band II column ever reached the cap. The
last two are the real \(i=9\) leftovers of §7 — long runs, and unremarkable
ones. (Their \(\Lambda\) is recomputed from the live-prime ladder, which is
derivable from \(N,K\); the per-round record of which primes each faced
lives only in the gitignored jsonl.)

**What the ledgers actually say.** No Band II pass and no Z-jump round
escalated, in any member whose record carries an escalation ledger — that
is \(i=2\) through \(i=8\); \(i=9\)'s record has none (§7). Two flags did
fire, both in the *small-\(k\) census*, a separate diagnostic that walks
\(k<200\) outside the sweep's bands:

* \(i=2\), \(k=65\): a run of 20 live primes, \(\Lambda=6.0\times10^{-8}\).
  This column lies **above \(k_{\max}=49\)**, where \(C(2k,k)>m\) and no
  left-half solution can exist at all; the theorem claims nothing there.
* \(i=7\), \(k=149\): a run of 8 live primes, \(\Lambda=4.2\times10^{-5}\),
  expected \(8.4\times10^{-3}\), just under the \(10^{-2}\) threshold. This
  column *is* in range, and it is **killed**: it dies at \(p=347\), which
  is its row in the witness table, and that certificate verifies.

Neither is a surviving column, and neither touches the theorem. Both are
what the trigger is for — flagging a run longer than the law expected — in
the regime where the model is loosest (documented model/exact ratio 1.263x
at small \(p\), against 1.00014x at Band II scale).

## 5. Exhaustion, per member

Every member is closed by three bands and the witness table is their union.
`results/i{i}_witness.npz` is the proof object; `results/i{i}_sweep.json` is
the run record. They are different things, and for \(i=8\) and \(i=9\) they
say different things — see §7.

| \(i\) | small-\(k\) band (rows) | Z-jump band (columns) | Band II band (columns) | run outcome |
|---:|---|---|---|---|
| 2 | \(k\le38\) (37) | — none | 41–49 (9) | Band II died pass 2; clean, certificate |
| 3 | \(k\le200\) (199) | 201–271 (71) | 274–342 (69) | Z 2 rounds, Band II 3 passes; clean, certificate |
| 4 | \(k\le200\) (199) | 201–1,868 (1,668) | 1,871–2,347 (477) | Z 3, Band II 3; clean, certificate |
| 5 | \(k\le200\) (199) | 201–12,814 (12,614) | 12,817–16,094 (3,278) | Z 4, Band II 6; clean, certificate |
| 6 | \(k\le200\) (199) | 201–87,839 (87,639) | 87,842–110,318 (22,477) | Z 6, Band II 6; clean, certificate |
| 7 | \(k\le200\) (199) | 201–602,068 (601,868) | 602,071–756,136 (154,066) | Z 9, Band II 7; clean, certificate |
| 8 | \(k=2\) (1) | 3–4,126,646 (4,126,644) | 4,126,649–5,182,637 (1,055,989) | Band II died pass 8; Z ran 14 rounds and left \(k=1021\) alive; **not clean, no certificate** |
| 9 | \(k\le80\) (79) | 81–28,284,463 (28,284,383) | 28,284,466–35,522,329 (7,237,864) | Band II died pass 9; Z stopped at its cap of 12 with four columns alive; **not clean, no certificate** |

The three counts sum exactly to the certificate count in every row (e.g.
\(79+28{,}284{,}383+7{,}237{,}864=35{,}522{,}326\) at \(i=9\)). The band
columns are *ranges of \(k\)*, not claims about which engine closed each
one: within the \(i=9\) Z band, four columns were closed by full-Lucas fills
rather than by the Z-jump, and at \(i=8\) one Z-band column was closed by
the engine after repair (§7). \(i=2\) is the degenerate shape —
\(k_{\max}=49\) is below the exact scan's reach, so there is no Z-jump band
at all and the sweep contributed only the nine Band II columns.

Two notes on reading the \(i=8\) row against its json: round 14 of the
Z-jump records `entered: 1, observed: 0`, which looks like the last column
dying, but \(k=1021\) was *dropped* from the accounting rather than killed —
§7. And three cap numbers are in play for that member: `cap_z: 15` in the
record (the small-\(k\) extension), 14 rounds actually run, and the "cap 12"
quoted in [`i8-N6.md`](i8-N6.md) §5, which is the base Z-jump cap before the
extension. They are consistent; they are not the same number.

**The small-\(k\) band.** Those rows carry witnesses from the engine —
`singmaster_intersect.py`'s `obstructing_prime`, the modular route that
tests small primes directly and never builds \(m\). Whether a *second,
independent* artifact also closes that band differs by member, and the
distinction matters:

* \(i=2..7\) — genuinely independent, and of a different evidence type:
  `results/fibonacci_i1-7.json` is an **exact** inversion, `k_extra: 200`,
  `also_central: true`, `extra_reps: []` for each of \(i=2..7\). It settles
  \(k\le200\) by construction rather than by obstruction.
* \(i=9\) — `results/modular_i9_k80.json` (79 columns \(k=2..80\), all
  `impossible`) is **the same computation recorded twice**, not a check:
  all 79 witness primes in the table are exactly this file's, \(k=2\to191\)
  and so on. It documents the fills; it does not corroborate them.
* \(i=8\) — the small-\(k\) band is the single column \(k=2\), whose
  witness \(p=227\) is likewise `modular_i8_k400.json`'s. That file does
  provide independent cover for \(k=3..400\), where its primes differ from
  the table's in 398 of 399 rows.

## 6. The certificate layer

`check_witness(N, K, k, p)` establishes four things, in order, from four
integers and nothing else:

1. \(p\) is prime, by its own deterministic Miller–Rabin, and \(p>k\ge2\);
2. \(r\) really is \(m\bmod p\), recomputed from \(N\) and \(K\) by every
   route that applies to that \(p\) — full Lucas always, two-digit when
   \(p^2>N\), the \(\delta\) closed form on the \(\alpha=1\) window — all
   agreeing;
3. \(r\neq0\);
4. \(r\notin I_{p,k}\), by walking the whole image.

`check_witness` itself shares no code path with the sweep: no factorial
table, no numpy, no gmpy2, no image set, walking
\(C(k+j,k)=(-1)^jC(p-k-1,j)\) as two modular multiplications per step in
\(O(1)\) memory. (The surrounding `witness.py verify` harness does use numpy
to read the `.npz`; the arithmetic that decides a verdict does not.) So the
builder is untrusted: propose the wrong prime and verification fails.

Three properties, checked separately, and the distinction is the whole
reason this section exists:

* **validity** — every stored \((k,p)\) really is an obstruction;
* **coverage** — the witnessed set equals
  \([2,k_{\max}]\setminus\{K,K+1\}\) exactly, with \(k_{\max}\) recomputed
  from \(N,K\) rather than read back from the file under audit;
* **binding** — a sweep record's certificate names *this* table by digest.

**Coverage** holds for all eight members: 0 missing, 0 extra, 41,590,228
columns. **Binding** holds for \(i=2..7\) and fails for \(i=8\) and
\(i=9\), whose runs ended `clean=false` and so emit no certificate; the
ledger reports that as its own state, *coverage complete, UNBOUND*, and
exits 1. Unbound is not a coverage hole and not a bad row — it means no run
record vouches for the file, so the closure rests on the table and the
ledger.

**Validity** is the one to be careful about. It is established in full for
the five small members and sampled for the three large ones:

| \(i\) | table sha256 | rows | \(\sum g\) (verifier steps) | verified |
|---:|---|---:|---:|---|
| 2 | `ab08157ad3407871` | 46 | 329 | **all**, 0 invalid |
| 3 | `028be3a55d945913` | 339 | 6,526 | **all**, 0 invalid |
| 4 | `64d3642c2346bd7e` | 2,344 | 222,074 | **all**, 0 invalid |
| 5 | `cbdbd736ff61d969` | 16,091 | 9,586,090 | **all**, 0 invalid |
| 6 | `fa1d304c3f462f50` | 110,315 | 442,960,748 | **all**, 0 invalid |
| 7 | `1ede909367ad5be5` | 756,133 | 20,748,416,901 | 20,000 sampled (2.6%), 0 invalid |
| 8 | `596dbf47543fa60a` | 5,182,634 | 974,216,503,629 | 5,000 sampled (0.097%), 0 invalid |
| 9 | `e375de140d32cf77` | 35,522,326 | 45,763,123,617,761 | 2,000 sampled (0.0056%), 0 invalid |

So 129,135 certificates have been checked exhaustively and 27,000 more by
sample: about 0.375% of the 41.6 million. Nothing invalid has ever been
found by a sample. The verifier is deliberately slow — its independence is
the point — and the cost is linear in \(\sum g\): at the rate measured on
the full \(i=6\) run (about \(1.3\times10^{7}\) image-steps per second
across 8 workers), a full \(i=8\) verification is roughly 21 hours of
wall-clock on this machine (~170 core-hours) and a full \(i=9\) roughly
40 days (~8,000 core-hours). \(i=8\) is an overnight job that has not been
run; \(i=9\) is not a laptop job.

## 7. Provenance exceptions

Where the proof object and the run record diverge. Stated here, not in a
footnote, because a reader checking the claim will meet them.

**\(i=8\), \(k=1021\): a false certificate, found and repaired.** Past its
cap the Z-jump's extended small-\(k\) rounds *dropped* survivors with
\(k\ge10^3\) instead of killing them, so \(k=1021\) — which had merely
survived \(p=3517\) — was recorded as killed by it. That row was false:
\(r(3517)=2662\) lies in \(I_{3517,1021}\) at \(j=14\). It was repaired to
\(p=1051\), which verifies (\(r=541\), \(g=30\)), and the drop was fixed at
the source, after which the run's record regenerated as `clean=false` with
the column listed alive. The column is impossible either way — the repaired
certificate proves it, independently of how it was found. The npz metadata
carries the repair rather than laundering it.

**\(i=8\), \(k=2\).** One engine fill: the Z-jump starts at \(k=3\) because
the kernel has no \(O(1)\) route for \(k=2\); the engine does (the
\(8m+1\) quadratic-residue test) and gives \(p=227\).

**\(i=9\), 83 engine fills.** Columns \(k=2..80\) (79) sit below the
Z-jump's start and are closed by the modular engine. The other four are the
interesting ones: \(k=87,399,553,1281\) were still alive when the Z-jump
stopped at 12 live primes. They are not anomalies (§4) and not unkilled.
The sweep could not kill them because its two-digit \(r(p)\) is false below
\(\lfloor\sqrt N\rfloor=8605\) and its cell geometry never scans there; each
dies at the *first* live prime below that floor, under full Lucas:

| \(k\) | \(p\) | \(r(p)\) |
|---:|---:|---:|
| 87 | 191 | 161 |
| 399 | 421 | 346 |
| 553 | 557 | 145 |
| 1281 | 1,321 | 434 |

All four verify through the same `check_witness` as every other row. On the
"cap of 12": `results/i9_sweep.json` records 12 rounds but no `cap_z` field,
so the number itself comes from `family_sweep.CAP_Z = 12` and the run log.
It is corroborated by the ladder: exactly 12 live primes lie in
\((8605,8779]\), every one of the four survivors has last prime 8779, and
the next live prime would be 8803 — the run stopped at a budget, not at
exhaustion.

**\(i=9\)'s run record is a resumed leg.** `results/i9_sweep.json` carries
the Z-jump phase but no Band II *per-pass* rows and no escalation ledger:
the job resumed after Band II had already written `phase_complete`, and
that path records the totals without the per-pass detail. The totals are
there (`n_bii: 7237864`, `n_bii_alive: 0`, the prime list); what is missing
is the pass-by-pass curve and the ledger. Nothing in the theorem depends on
either — it is a gap in the record, and it is why §4 regenerates the
prediction rather than quoting the file.

**Termination certificates, for a subset.** 29 columns \(k\le30\) of
\(i=8\) (`results/termination_i8.json`, 29 rows, \(k=2..30\)) carry a
stronger statement than "a killing prime was found": if \((x)_k-k!\,m\) is
irreducible over \(\mathbf Q\) then Jordan gives a derangement and
Chebotarev a killing prime of density \(\ge1/k\), so a killing prime *had to
exist*. Irreducibility has a one-prime certificate computed by Lucas
without ever building \(m\). The same is recorded for four \(i=9\) columns
— \(k=11,29,40,45\), the long-run columns of its small-\(k\) census, **not**
the four cap columns above — in `docs/open-questions.md`; no
`termination_i9.json` is committed, so that half is a claim in a doc rather
than an artifact in `results/`. This upgrades those columns; it does not
weaken the others, which stand on their own certificates.

## 8. What this is not

* **Not Singmaster's conjecture**, which bounds \(N(t)\) for every \(t\).
  This is eight values.
* **Not a statement about \(i\ge10\)**, and not a reason to start one. Cost
  grows as \(\varphi^{8}\approx47\times\) per member and the outcome is
  predicted before the run by a law with no fitted parameters; see
  [`q11-what-is-the-claim.md`](q11-what-is-the-claim.md).
* **Not "3003 is the only \(N=8\)"**, and not "this is the only infinite
  family". Neither is touched.
* **Not the `nearby` or `collide` results.** Different searches over
  different objects; a null there is not part of this theorem. (Nor are the
  recorded `collide` nulls contiguous with the 2017 bound: they start at
  62-digit values, one decade above the first value the classifier calls
  past. The scan start was corrected on 2026-08-23; the artifacts were not
  re-run.)
* **Not `results/coverage_ledger.json`.** Still stale: the ledger json
  predates the coverage/binding split, and its \(i=8\) row still asserts a
  certificate digest that the current \(i=8\) record does not contain. Run
  `coverage_ledger.py`; do not quote that file. (`results/ghost_census.json`
  was in the same boat until 2026-08-24; it has since been regenerated from
  the current tables — eight members, digests matching this section's table,
  live-pinned by `scripts/test_ghost_census.py` — and IS now citable, with
  its own caveats.)
* **Not reproducible from a clone alone.** Every `results/i*_sweep.jsonl` is
  gitignored, and those checkpoints are the only record of which prime
  killed which column. The witness tables are committed and can be
  *verified* from a clone (§9); they cannot be *rebuilt* from one.

## 9. How to check this

All read-only. Timings are wall-clock on this machine with the default 8
workers, not core-time.

```text
python scripts/coverage_ledger.py
```

Expect the eight-row table of §1, `0 missing, 0 extra` on every row,
`columns witnessed 41,590,228`, i=2..7 COMPLETE, i=8 and i=9
`coverage complete, UNBOUND`, and

```text
  RESULT UNBOUND MEMBERS: i=8, i=9
```

with **exit status 1** — the expected outcome, not a failure: the script
exits non-zero while any member is unbound.

Full verification of the five small members — every certificate, no sample:

```text
python scripts/witness.py verify --file results/i3_witness.npz
```

`certificates checked 339   invalid=0` → `RESULT VALID`, about 0.12 s. The
same for \(i=2\) (46, instant), \(i=4\) (2,344, ~0.15 s), \(i=5\) (16,091,
~0.73 s) and \(i=6\) (110,315, ~34 s) — all VALID, 0 invalid.

The three large members are sampled, per §6:

```text
python scripts/witness.py verify --file results/i7_witness.npz --sample 20000
python scripts/witness.py verify --file results/i8_witness.npz --sample 5000
python scripts/witness.py verify --file results/i9_witness.npz --sample 2000
```

→ `RESULT VALID (coverage complete; certificates sampled)`, on the order of
1, 2 and 4 minutes respectively (they vary with machine load). Coverage is
complete in all three; the sampling is over *validity* only.

The repaired \(i=8\) row, both directions — the two lines a referee should
run if they run nothing else:

```text
python scripts/witness.py one --N 10803704 --K 4126647 --k 1021 --p 1051
```

→ `"ok": true, "r": 541, "g": 30`, exit 0. And the prime it replaced:

```text
python scripts/witness.py one --N 10803704 --K 4126647 --k 1021 --p 3517
```

→ `"ok": false, "reason": "r(p) IS in the image", "r": 2662, "j": 14`,
exit 1.

The four \(i=9\) cap columns, the same way — e.g.

```text
python scripts/witness.py one --N 74049690 --K 28284464 --k 87 --p 191
```

→ `"ok": true, "r": 161`. Likewise \(399/421\), \(553/557\),
\(1281/1321\).

Nothing above rebuilds a table, writes into `results/`, or starts a sweep.
