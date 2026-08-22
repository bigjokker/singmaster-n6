# Singmaster intersect

Search past the Blokhuis–Brouwer–de Weger (2017) envelope for extra binomial
representations of the Lind/Singmaster/Tovey \(N=6\) Fibonacci family, plus
unsettled nearby-row and column-pair collisions.

This does **not** prove Singmaster's conjecture. It does prove
\(N\bigl(C(F_{18}F_{19},F_{16}F_{19})\bigr)=6\) exactly: every extra left-half
column \(2\le k\le k_{\max}\), \(k\notin\{K,K+1\}\), carries a modular kill
certificate — and since 2026-08-20 every one of those certificates is
**recorded and independently checkable**, for i=2 through i=8.

## Engine

`singmaster_intersect.py` (needs `gmpy2`):

| Command | Question |
|---|---|
| `intersect` | Exact extra-rep on a Fibonacci member (builds \(m\)) |
| `modular` | Lucas / image obstruction: prove a column cannot represent \(m_i\) **without** building \(m\) |
| `nearby` | \(C(n,k)=C(n-d,k+e)\) for unsettled \((d,e)\) — exhaustive, not sampled |
| `collide` | Finite \(m\)-slice of \(C(n,k)=C(m,l)\) |
| `sanity` | Catalog, classifier, Lucas, image, 3003 tripwire |

```text
python singmaster_intersect.py sanity
```

`impossible` + a witness prime is a column-level proof. `possible` only means
the prime list did not kill it.

### Tests

| Suite | Covers |
|---|---|
| `singmaster_intersect.py sanity` | the engine |
| `scripts/test_kernel.py` | the image-scan kernel |
| `scripts/test_sweeps.py` | the sweep drivers' certificate logic |
| `scripts/test_witness.py` | the certificate verifier |
| `scripts/test_sizelaw.py` | the size law and the escalation trigger |
| `scripts/test_work_census.py` | the exact scan-work census and its gates |
| `scripts/test_rp_cost.py` | that r(p), not the scan, is a Z-jump job's cost |

```text
python scripts/test_sweeps.py
```

## Certificates

A sweep proves *every* extra column has a killing prime. `scripts/witness.py`
records **which** prime, per column, and re-checks it independently.

| Command | Question |
|---|---|
| `witness build --i N` | Rebuild the witness table from a run's checkpoint |
| `witness fill --i N` | Add engine witnesses for columns the table lacks |
| `witness verify --file F` | Re-check every certificate, and that none is missing |
| `witness one --N .. --K .. --k .. --p ..` | Check a single certificate |

```text
python scripts/witness.py verify --file results/i8_witness.npz --sample 5000
```

`check_witness` establishes four things, in order: **(i)** \(p\) is prime — by
its own Miller–Rabin, because both Lucas and the claim that \(I_{p,k}\) is the
*complete* image assume it; **(ii)** \(r\) really is \(m \bmod p\), recomputed
from \(N,K\) by every route that applies, all agreeing; **(iii)** \(r\neq 0\),
since 0 lies in every column image and a dead prime certifies nothing;
**(iv)** \(r\) is outside \(I_{p,k}\), by walking the whole image.

The verifier shares no code path with the sweep — no factorial table, no numpy,
no gmpy2 — walking \(C(k+j,k)=(-1)^jC(p-k-1,j)\) as two modular multiplications
per step in **O(1) memory**. So the builder is untrusted: propose the wrong
prime and verification fails. The hardest i=7 certificate checks in 32 ms.

`verify` reports two independent things and needs both: no certificate invalid,
**and** no claimed column missing one. A clean `family_sweep` writes
`results/i{N}_witness.npz` automatically and names it in the certificate.

### Coverage

Validity is not completeness. `scripts/coverage_ledger.py` states the global
claim and checks it — `witnessed(i) == [2, k_max] \ {K, K+1}`, exactly, with
\(k_{\max}\) **recomputed from \(N,K\)** rather than read back from the file
under audit:

```text
python scripts/coverage_ledger.py
```

All seven members report COMPLETE: **6,067,902 columns, 0 missing, 0 extra**.

### Independent re-checks

`scripts/verify_independent.py` re-checks certificates through sympy plus two
from-scratch controls, choosing the route by cost: `math.comb` brute force
where \(p\) is small, GF(\(p\)) factorisation where \(k\) is small, and a
freshly written factorial identity elsewhere. It agrees on every certificate
tried. Its limits are stated in the file — sympy is a different library, not a
different language, and the factorial route shares the sweep's identity.

`scripts/termination_certificate.py` proves a column *must* die: if
\((x)_k - k!m\) is irreducible over \(\mathbf{Q}\), Jordan gives a derangement
and Chebotarev a killing prime of density \(\ge 1/k\). Irreducibility has a
one-prime certificate, computed by Lucas without ever building \(m\). All 29
columns \(k\le 30\) of i=8 and all four of i=9's long-run columns are
certified — upgrading them from *a killing prime was found* to *a killing prime
had to exist*.

`scripts/ghost_census.py` reads the witness tables as a test of the
falling-factorial local–global question: each recorded \((k,p)\) certifies that
\(c = k!\,C(N,K)\) lies outside \((x)_k(\mathbf{F}_p)\), hence is not a
"ghost". The output carries its own caveats, worth reading before quoting the
count.

## Results worth keeping

Witness tables are the proof object; the sweep JSONs are the run records.

| File | Claim |
|---|---|
| `results/coverage_ledger.json` | Every member covers \([2,k_{\max}]\setminus\{K,K+1\}\) exactly. 6,067,902 columns |
| `results/i8_witness.npz` | i=8: 5,182,634 per-column certificates, \(k=2..5{,}182{,}637\) |
| `results/i7_witness.npz` | i=7: 756,133 certificates |
| `results/i2..i6_witness.npz` | 46 / 339 / 2,344 / 16,091 / 110,315 certificates |
| `results/i8_sweep.json` | i=8 by one uniform method, \(k=3..k_{\max}\) plus Band II. clean, 7508 s on 4 workers |
| `results/i7_sweep.json` | i=7 Band II + Z-jump, \(N=6\) (95 s) |
| `results/i2..i6_sweep.json` | same pipeline, \(N=6\) each (0.7–7.9 s) |
| `results/fibonacci_i8_k300.json` | Exact: i=8 has no extra left-half with \(2\le k\le 300\), no central. ~16.3 h |
| `results/modular_i8_k400.json` | Modular: i=8, \(k=2..400\) all impossible |
| `results/modular_i9_k80.json` | Modular: i=9, \(k=2..80\) all impossible |
| `results/modular_i10_k20.json` | Modular: i=10, \(k=2..20\) all impossible |
| `results/nearby_k8M-1B.json` | Nearby, \(k=8\text{M}..10^9\), \(d,e\le 6\), `new_hits=0` |
| `results/nearby_k2M-8M_de8.json` | Nearby, \(k=2\text{M}..8\text{M}\), \(d,e\le 8\), `new_hits=0` |
| `results/intersective_search.json` | \(k=5..15\), \(|c|\le 10^9\): every survivor a genuine falling-factorial value, 0 intersective |
| `results/k6_intersective.json` | \((x)_6-c\) is **never** intersective — proved, not searched |
| `results/k7_intersective.json` | \((x)_7-c\): only \(c=\pm896,\pm17472,\pm459648\) even factor, and 5 divides none. Modulo effective Thue |
| `results/k8_intersective.json` | \((x)_8-c\): five cases, all closed; the exceptional locus is a genus-3 curve with only degenerate points |
| `results/k10_intersective.json` | \((x)_{10}-c\): all cases closed; only a genus-1 gap remains. Difficulty is not monotone in \(k\) |
| `results/k9_intersective.json` | \((x)_9-c\): three branches, all closed; \(c=\pm2630880\) is the first to pass \(\mathrm{rad}(k!)\) and dies at \(p=13\) |
| `results/termination_i8.json` | 29 i=8 columns with a proof that a killing prime must exist |
| `results/ghost_census.json` | Certified non-ghosts, with its own caveats |
| `results/bandii_sweep.json` | Historical Band II \(p>N/2\): 1,055,989 columns, died at prime 8 |
| `results/zjump.json` | Historical Band I remnant: 3,215,816 columns, 0 anomalies |
| `results/walk_369.json` | 369 triples walked past prime 3: max run 6, counts match the size law |
| `results/stragglers_nearK.json` | 25 near-\(K\) Band I \(k\) all killed at \(p>N/2\) |

The historical `bandii_sweep.json` / `zjump.json` runs kept no per-column
witness, and their checkpoints are gone — which is why i=8 was re-derived in
one uniform pass. Do not build \(m_{10}\) (~147 million digits) to answer
small-\(k\) extra-rep; modular already killed those columns.

## i=8 extra columns

Every extra \(k\in[2,k_{\max}]\setminus\{K,K+1\}\) has an unconditional
modular certificate, all 5,182,634 of them recorded and verifiable.
\(N(C(F_{18}F_{19},F_{16}F_{19}))=6\). Not Singmaster.

- [`docs/i8-N6.md`](docs/i8-N6.md) — the \(N=6\) theorem, as a note
- [`docs/band-I.md`](docs/band-I.md) — lemmas and census behind it
- [`docs/bandii-spec.md`](docs/bandii-spec.md) — Band II \(p>N/2\) sweep
- [`docs/zjump-spec.md`](docs/zjump-spec.md) — Band I Z-jump remnant
- [`docs/zeromap-p1e5-1e6.md`](docs/zeromap-p1e5-1e6.md) — 136 digit-windows, 38 NONE

## Escalation

A survivor at the prime cap is an anomaly only if the size law did not expect
it. **Run length is not the criterion** — survival per live prime is 0.61 at
small \(k\) and under 0.16 in a fat cell, so a run of 8 at \(i=9,k=11\) is
ordinary while a run of 6 in an \(i=8\) fat cell is a thousand times rarer.

`scripts/sizelaw.py` computes the expected survivor count per round, summed
over the columns that actually entered it with the primes they were actually
tested against, and `family_sweep.py` records expected-vs-observed for every
pass. Escalate when a column survives a round that expected almost none
(\(E_r<10^{-2}\)), or when far more survive than the law allows (Poisson upper
tail \(<10^{-3}\)).

```text
python scripts/sizelaw.py run --i 9 --k 11
python scripts/sizelaw.py predict --check
```

The law has no fitted parameters: \(|I_{p,k}|\) is the proved involution bound
(\(g+1\) for odd \(k\), \(\lceil g/2\rceil+1\) for even) minus birthday
collisions. Its accuracy and the trigger's headroom are **anti-correlated in
the safe direction** — 1.00014x against exact images at Band II scale where
headroom is 21x, and 1.263x at small-\(k\) census scale where headroom is
1527x. Recomputing from exact image sizes flips **no** verdict.

## Profiling

Several performance claims on this project survived reasoning and died on
measurement, so the pipeline carries a profiler rather than an argument:

```text
python scripts/profile_sweep.py --i 9
```

It samples each phase and reports the split between the image scan, r(p) and
prime assignment. The split is the part to trust; wall-clock projection is
good to about a factor of two, so quote \(i=10\) as an order of magnitude, not
a number. Since the factorial table was removed (below), **both phases are now
scan-dominated**.

The distinct-prime count integrates prime density over the live intervals
rather than sampling a window near the range start, which is biased ~2x high;
the integrated form recovers i=8's recorded 124,830 primes to 0.2%.

## Checkpoints

Every run's jsonl opens with a schema header pinning the version and the run's
parameters (`i`, `N`, `K`, `k_max`, `k_lo_z`, caps). Resume verifies it and
refuses on mismatch, or on a checkpoint with no header at all. Resume merges
old records with new ones, so a silent format or parameter change would produce
a certificate over columns that were never all tested the same way.

A run's jsonl is the only record of which prime killed which column, so a clean
`family_sweep` builds the witness table from it before it stops mattering.

## Notes

### Documents

- [`docs/open-questions.md`](docs/open-questions.md) — the ranked question list, Q1–Q24, with outcomes
- [`docs/audit-2026-08-20.md`](docs/audit-2026-08-20.md) — first-principles audit of the pipeline
- [`docs/extra-questions.md`](docs/extra-questions.md) — Q2/Q3/Q8/Q12/Q15/Q16/Q17/Q18, with the declines and why
- [`docs/q14-intersective.md`](docs/q14-intersective.md) — can \((x)_k-c\) be intersective? Settled for \(k\le 5\)
- [`docs/q25-k6-intersective.md`](docs/q25-k6-intersective.md) — \(k=6\) settled too: for even \(k\) the sextic is a **cubic** in \(t^2\)
- [`docs/q26-k7-intersective.md`](docs/q26-k7-intersective.md) — \(k=7\), modulo one effective-Thue step: odd \(k\) gets **no** reduction
- [`docs/q27-k8-intersective.md`](docs/q27-k8-intersective.md) — \(k=8\): the difficulty hides in a **genus-3 curve**, whose points are exactly the degenerate ones
- [`docs/q28-k10-intersective.md`](docs/q28-k10-intersective.md) — \(k=10\) is **cheaper than \(k=8\)**: difficulty tracks cubic factors, not \(k\)
- [`docs/q29-k9-intersective.md`](docs/q29-k9-intersective.md) — \(k=9\): eight shapes collapse to three branches; a **second** obstruction mechanism
- [`docs/q11-what-is-the-claim.md`](docs/q11-what-is-the-claim.md) — what another family member buys, against what it costs
- [`docs/interior-2022.md`](docs/interior-2022.md) — the MRSTT interior theorem, and why it does not reach i=8/9
- [`docs/modular-spec.txt`](docs/modular-spec.txt) — Lucas/modular layer: what a certificate is, what not to rebuild
- [`docs/campaign-log.txt`](docs/campaign-log.txt) — original search campaign

### Invariants worth not breaking

- A sweep is `clean` only if nothing survived **and** every column was actually
  testable. A column with no live prime left is an anomaly, not a kill.
- A scan with \(r(p)=0\) is refused, not run: 0 lies in every column image, so
  such a prime certifies nothing, and a mask scan would silently report every
  column killed.
- Pre-flight guards use `check()`, not `assert` — `python -O` strips `assert`,
  and every certificate downstream of a pre-flight is claimed unconditional.
- The proved image bound \(|I|\le g+1\) (odd \(k\)) / \(\lceil g/2\rceil+1\)
  (even) is asserted where the image is built. It is a theorem, so a violation
  is a bug.
- Settled theorem pairs \((k,l)\) and nearby \((1,1),(1,2),(2,1)\) are not
  re-sieved. The classifier uses `is_fibonacci_pair`, not "rows differ by 1".
- `witness.binom_mod_prime_pure` is deliberately left naive. The verifier's
  independence rests on full Lucas being a *different computation* from the
  delta route; optimising it would make them algebraically the same.

### Kernel

- Shared Band I / cell-geometry helpers (`cells`, `live_intervals`,
  `first_live_after`, `image_j`, `inv_table`, `binom_mod_prime`, `chunk_ks`,
  jsonl I/O) live in [`scripts/bandii_kernel.py`](scripts/bandii_kernel.py).
  They were copy-pasted into six scripts and had begun to drift; do not
  re-inline them.
- The image scan uses the involution \((k-1-x)_k=(-1)^k(x)_k\) to test only
  \(b\in[0,\lceil g/2\rceil)\): even \(k\) repeats on the upper half, odd \(k\)
  negates, and the \(-s\) branch is a subtraction rather than a second multiply.
  Exact, ~2.2x measured, byte-identical output including the witness index
  \(b\). `USE_HALF_SCAN=False` switches back for A/B.
- The scan builds only the \(O(g)\) factorial entries it reads, not all \(p\).
  Wilson's theorem gives \(k! = (-1)^{k+1}/(g-1)!\), so reaching \(F[k]\) costs
  \(g-1\) multiplications rather than \(k\). **11.8x on real Z-jump buckets**;
  a wash on Band II, which amortises one table over a whole chunk.
  `USE_WINDOWED_SCAN=False` switches back.
- `binom_mod_prime` takes `min(k, n-k, p-n-1)` — three lower indices, not two.
  The third comes from \(n\equiv-(p-n)\), giving
  \(C(n,k)=(-1)^kC(p-n+k-1,\,p-n-1)\), which is tiny when \(n\) is near \(p\).
  Band II Lucas: 100+ ms to 0.001 ms.
- `first_live_after` bisects into the interval list instead of rescanning from
  index 0 — safe because `live_intervals` returns disjoint intervals. 292.5 to
  2.7 µs/call, which is 2.3 h of serial parent time at i=9 down to about a
  minute.
- `_column_possible_scan` walks the image as \(\{(-1)^jC(g-1,j)\}\) keeping
  \(A_j=(g-1)_j\) and \(B_j=j!\) — two multiplications per step, no modular
  inverse — and folds by the same involution. `_column_possible_scan_ref` is
  the superseded version, kept because sanity checks the two agree.
- `column_possible` dispatches by cost: QR test for \(k=2\), cached image for
  \(p\le 4000\), polynomial gcd when \(k^2\log p < g\), else the folded scan.
  The gcd route needs no radicals, so it works for every \(k\) — 1239x at
  \(k=3\), 0.2x at \(k=400\), crossover near \(k=150\).
- `nextprime_sweep` memoises \(r(p)\) per prime (`RCache`), not per \((k,p)\)
  pair. 7.3x on the Stage-2 range, 13–15x on Stage-3-style \(k\). Rows are
  identical with it on or off, and `test_sweeps.py` pins that against the
  original recorded runs.
