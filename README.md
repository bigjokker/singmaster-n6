# Singmaster intersect

Search past the Blokhuis–Brouwer–de Weger (2017) envelope for extra binomial
representations of the Lind/Singmaster/Tovey \(N=6\) Fibonacci family, plus
unsettled nearby-row and column-pair collisions.

This does **not** prove Singmaster's conjecture. It does prove
\(N\bigl(C(F_{18}F_{19},F_{16}F_{19})\bigr)=6\) exactly: every extra left-half
column \(2\le k\le k_{\max}\), \(k\notin\{K,K+1\}\), carries a modular kill
certificate — and every one of those certificates is **recorded and
independently checkable**: for i=2 through i=8 since 2026-08-20, and for i=9
since 2026-08-22, when its table was harvested (`c061eeb`).

Two of those members ship *unbound*: i=8 and i=9 have complete, verifiable
witness tables, but their sweep records carry `certificate: null`, because
neither run ended clean (i=8 dropped a column it had not killed; i=9 left
four alive at its cap). That is the honest state, and the coverage ledger
now says so rather than printing COMPLETE — see **Coverage** below.

## Engine

`singmaster_intersect.py` (needs `gmpy2`):

| Command | Question |
|---|---|
| `intersect` | Exact extra-rep on a Fibonacci member (builds \(m\)) |
| `modular` | Lucas / image obstruction: prove a column cannot represent \(m_i\) **without** building \(m\) |
| `nearby` | \(C(n,k)=C(n-d,k+e)\) for unsettled \((d,e)\) — exhaustive, not sampled |
| `collide` | Finite \(m\)-slice of \(C(n,k)=C(m,l)\). The past-2017 start now abuts the bound (first 61-digit value; 101 for \(l\ge10\)); the recorded `results/collide_*.json` were run from one decade higher (62 / 102 digits) and do **not** cover that first decade |
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
| `scripts/test_rp_cost.py` | the r(p)-vs-scan split of a Z-jump job on the *no-`r_expected`* fallback path |
| `scripts/test_collide.py` | the collide frontier abuts the 2017 bound (first 61 / 101-digit value, not one decade up) |
| `scripts/test_k10_intersective.py` | Q28's (2,3) branch: the curve, its 15 small integral points, the rank-2 Jacobian, the two kills, and that nothing claims the list complete |
| `scripts/test_k10_chebotarev.py` | Q28's Chebotarev pass: the kill trichotomy, the no-kill shape's reality, c<0 closed, the locus on its genus-3 cover, still BLOCKED |
| `scripts/test_k7_runge.py` | Q26: Branch B's Runge certificate (empty trap, 19 points, c = 0, +-896 only) and Branch A's corrected geometry (genus 1, rank >= 1, the D5/F20 trap) |
| `scripts/test_k8_case2.py` | Q27 proved: the genus-3 curve's 16 descent classes all die; the five degenerate points; the mod-8 filter; the 2-adic depth-13/6 certificates |

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

It reports two facts per member, separately, and requires both:
*coverage-complete* (the witnessed set is exactly the claimed set) and
*bound* (a sweep certificate names this table by digest). Eight members,
**41,590,228 columns, 0 missing, 0 extra**; i=2..7 COMPLETE AND BOUND; i=8
and i=9 *coverage complete, UNBOUND*, since their sweep records carry no
certificate. Unbound is neither a coverage hole nor a clean bill, so the
script prints it as its own state and **exits 1**.

`results/coverage_ledger.json` in the tree predates that split and still
shows the old single verdict; the script is the current statement, not that
file.

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

**`results/ghost_census.json` is stale and has not been regenerated.** It
lists seven members and 6,067,902 values, so it predates i=9; worse, its i=8
entry names digest `b4c02030`, the table as it stood *before* the k=1021
repair — the one that held the false 1021 → 3517 row. So its count is not a
statement about the current i=8 table and should not be quoted as one.

## Results worth keeping

Witness tables are the proof object; the sweep JSONs are the run records.

| File | Claim |
|---|---|
| `results/coverage_ledger.json` | Stale: written before coverage and binding were separated. Run `scripts/coverage_ledger.py` for the current statement (8 members, 41,590,228 columns, i=8/i=9 UNBOUND) |
| `results/i9_witness.npz` | i=9: 35,522,326 per-column certificates, \(k=2..35{,}522{,}329\); 83 engine-filled (k=2..80 and the four Lucas columns) |
| `results/i8_witness.npz` | i=8: 5,182,634 per-column certificates, \(k=2..5{,}182{,}637\); two rows are not from the sweep (engine fill k=2, repair 1021: 3517 → 1051) |
| `results/i7_witness.npz` | i=7: 756,133 certificates |
| `results/i2..i6_witness.npz` | 46 / 339 / 2,344 / 16,091 / 110,315 certificates |
| `results/i8_sweep.json` | i=8 by one uniform method, \(k=3..k_{\max}\) plus Band II. **`clean=false`**, `certificate=null`, `n_z_alive=1` (k=1021), 174 s on 8 workers — the honest record of the regenerated run |
| `results/i9_sweep.json` | i=9 Band II + Z-jump. **`clean=false`**, `certificate=null`, four columns alive at Z-jump cap 12 (k=87/399/553/1281), later killed by full Lucas below \(\sqrt N\) |
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
| `results/k7_intersective.json` | \((x)_7-c\): only \(c=\pm896,\pm17472,\pm459648\) even factor, and 5 divides none. Historical: see `k7_runge.json` |
| `results/k7_runge.json` | \((x)_7-c\), 3+4 branch CLOSED by Runge: all 19 quintic-curve points have \(A\le36\); complete list \(c\in\{0,\pm896\}\). The 2+5 branch is a rank-\(\ge1\) genus-1 cubic (never a Thue equation), still blocked |
| `results/k8_intersective.json` | \((x)_8-c\): five cases, all closed; the exceptional locus is a genus-3 curve with only degenerate points. Historical: see `k8_case2.json` |
| `results/k8_case2.json` | \((x)_8-c\) is never intersective -- **PROVED**: the genus-3 curve is solved by descent (16 classes: compactness, Runge, congruences at 2-adic depth 13/6) |
| `results/k10_intersective.json` | \((x)_{10}-c\): four cases closed, the (2,3) case reduced to one Magma `IntegralQuarticPoints` call (Jacobian of rank 2; not run). Difficulty is not monotone in \(k\) |
| `results/k9_intersective.json` | \((x)_9-c\): three branches, all closed; \(c=\pm2630880\) is the first to pass \(\mathrm{rad}(k!)\) and dies at \(p=13\) |
| `results/termination_i8.json` | 29 i=8 columns with a proof that a killing prime must exist |
| `results/ghost_census.json` | Certified non-ghosts — **stale**: 7 members, and its i=8 digest is the pre-repair table |
| `results/bandii_sweep.json` | Historical Band II \(p>N/2\): 1,055,989 columns, died at prime 8 |
| `results/zjump.json` | Historical Band I remnant: 3,215,816 columns, 0 anomalies |
| `results/walk_369.json` | 369 triples walked past prime 3: max run 6, counts match the size law |
| `results/stragglers_nearK.json` | 25 near-\(K\) Band I \(k\) all killed at \(p>N/2\) |

The historical `bandii_sweep.json` / `zjump.json` runs kept no per-column
witness, and their checkpoints are gone — which is why i=8 was re-derived in
one uniform pass. Do not build \(m_{10}\) (~147 million digits) to answer
small-\(k\) extra-rep; modular already killed those columns.

## Family i=2..9

\(N\bigl(C(F_{2i+2}F_{2i+3},\,F_{2i}F_{2i+3})\bigr)=6\) for every
\(i=2,\dots,9\): 41,590,228 extra columns, each with a recorded,
independently checkable modular certificate, and no central solution.
i=8 and i=9 are *coverage complete, UNBOUND* — their runs ended
`clean=false`, so no sweep certificate names their tables, which does not
touch the theorem.

- [`docs/family-N6.md`](docs/family-N6.md) — the family \(N=6\) note:
  statement, lemmas, per-member exhaustion, the provenance exceptions, and
  the commands to re-check it

## i=8 extra columns

Every extra \(k\in[2,k_{\max}]\setminus\{K,K+1\}\) has an unconditional
modular certificate, all 5,182,634 of them recorded and verifiable.
\(N(C(F_{18}F_{19},F_{16}F_{19}))=6\). Not Singmaster.

The claim rests on the **table**, not on the sweep record: two of those rows
did not come from the sweep (k=2 filled by the engine, k=1021 repaired
3517 → 1051 after the sweep credited a prime it had only survived), and
`results/i8_sweep.json` is correspondingly `clean=false` with no certificate.
Coverage is complete and every sampled certificate verifies; the ledger
therefore reports i=8 as *coverage complete, UNBOUND*.

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

That rule is per round, and it is blind to one case: when a single column
enters a round, \(E_r\) is just that column's own next-prime survival
(0.04–0.21 in Band II), never below \(10^{-2}\) — so a lone column
surviving *every* pass to the cap read "ordinary" in every round, though its
run is the anomaly the trigger exists for. So each column still alive at a
phase's cap is judged again, as a column: \(\Lambda\) over the primes it
actually faced, against the columns that entered the phase, with the same
threshold. i=8 Band II \(k=4{,}126{,}649\) through 14 passes would fire
(\(\Lambda=3.3\times10^{-10}\), peers \(\times\) \(\Lambda=3.5\times10^{-4}\));
i=9's \(k=11\) run of 8 and its four columns at the Z-jump cap do not. The
phase verdict is the per-round fire **or** the cap-survivor fire.

```text
python scripts/sizelaw.py run --i 9 --k 11
python scripts/sizelaw.py run --i 8 --k 2227205 --cap 8 --ladder
python scripts/sizelaw.py predict --i 9 --check
```

`predict --i N` regenerates member N's own Band II curve from \((N,K)\)
alone; `--check` scores i=8 against the pre-registration fixed before its
run and any other member against the window and primes its run record
carries (a member is never scored against another member's table). `run`
assesses a column only if it was tested to a kill or to the cap: a walk
whose prime budget runs out inside a forced-zero slab is reported as
**NOT TESTED** with exit 2, never as "ordinary" — use `--ladder` for a
Band I / fat-cell column, which walks the member's live-prime ladder.

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

It reports each phase's scan work; wall-clock projection is good to about a
factor of two, so quote \(i=10\) as an order of magnitude, not a number.
What the code makes checkable without a new measurement: a production
Z-jump worker computes **no \(r(p)\) at all**. The parent reduces \(m\)
against a whole round's primes at once (`USE_M_FOR_RP`, a remainder tree),
so `_job` takes its `r_expected` branch and calls only the scan; Band II's
\(r(p)\) is a handful of multiplies, also parent-side. Note what that means
for `scripts/test_rp_cost.py`: it measures the *other* branch, the
no-`r_expected` fallback, where r(p) is indeed almost the whole cost. No
artifact in the tree currently measures a per-phase scan-vs-r(p) split on
the GM kernel, so none is quoted here.

**Stale, and knowingly so:** `work_census.SCAN_RATE = 1.865e8` was measured on
the old numpy `(s*F) % p` loop, and production has run the
Granlund–Montgomery kernel since `8e64945`. Every figure derived from that
constant is stale with it — the census's core-hour column, its "5.6x / 9.4x",
and the profiler's "6.11x at i=7 and 10.27x at i=8". The visible symptom is
i=8 printing "scan is 209.66% of wall", which is the constant, not the run.
Both files say so in place. Re-measuring is a job, not a doc fix, so no
replacement number is quoted here.

## Checkpoints

A run's jsonl opens with a schema header pinning the version and the run's
parameters (`i`, `N`, `K`, `k_max`, `k_lo_z`, the caps, and since 2026-08-23
`n_chunks`, because the chunk partition *is* the resume key). Resume refuses
on mismatch, or on a checkpoint with no header at all: it merges old records
with new ones, so a silent format or parameter change would produce a
certificate over columns that were never all tested the same way.

Two limits worth knowing. The comparison covers only keys present in **both**
the header and the current parameters, so a checkpoint written before a key
existed still resumes — deliberate, and why `n_chunks` could be added without
a schema bump. And not every jsonl in the tree has a header:
`results/i9_sweep.jsonl`, the file i=9's table was built from, predates the
schema and starts with a data record; `scripts/migrate_checkpoint.py` exists
to stamp one on such a file after checking its records against the
parameters, and `witness.py build` does not require one.

A record's coverage is exact from `[k_lo, k_hi]` only when the chunk is
contiguous. Since 2026-08-23 a *sparse* chunk also carries `ks`, its exact
column list, so the builder never infers coverage from a span that
over-covers; checkpoints written before that keep the over-covering span and
are read exactly as before.

A run's jsonl is the only record of which prime killed which column, so a clean
`family_sweep` builds the witness table from it before it stops mattering.

## Notes

### Documents

- [`docs/open-questions.md`](docs/open-questions.md) — the ranked question list, Q1–Q30, with outcomes
- [`docs/audit-2026-08-20.md`](docs/audit-2026-08-20.md) — first-principles audit of the pipeline
- [`docs/extra-questions.md`](docs/extra-questions.md) — Q2/Q3/Q8/Q12/Q15/Q16/Q17/Q18, with the declines and why
- [`docs/q14-intersective.md`](docs/q14-intersective.md) — can \((x)_k-c\) be intersective? Settled for \(k\le 5\)
- [`docs/q25-k6-intersective.md`](docs/q25-k6-intersective.md) — \(k=6\) settled too: for even \(k\) the sextic is a **cubic** in \(t^2\)
- [`docs/q26-k7-intersective.md`](docs/q26-k7-intersective.md) — \(k=7\), 3+4 branch closed by Runge, 2+5 modulo one elliptic integral-point step: odd \(k\) gets **no** reduction
- [`docs/q27-k8-intersective.md`](docs/q27-k8-intersective.md) — \(k=8\), **proved** (2026-08-23): the genus-3 curve is solved by descent — its points are exactly the degenerate ones
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
