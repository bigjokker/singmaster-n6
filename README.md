# Singmaster intersect

Search past the Blokhuis–Brouwer–de Weger (2017) envelope for extra binomial representations of the Lind/Singmaster/Tovey \(N=6\) Fibonacci family, plus unsettled nearby-row and column-pair collisions.

This does **not** prove Singmaster’s conjecture. It does prove
\(N\bigl(C(F_{18}F_{19},F_{16}F_{19})\bigr)=6\) exactly: every extra
left-half column \(2\le k\le k_{\max}\), \(k\notin\{K,K+1\}\), has a
modular kill certificate.

## Engine

`singmaster_intersect.py` (needs `gmpy2`):

| Command | Question |
|---|---|
| `intersect` | Exact extra-rep on a Fibonacci member (builds \(m\)) |
| `modular` | Lucas / image obstruction: prove a column cannot represent \(m_i\) **without** building \(m\) |
| `nearby` | Sampled \(C(n,k)=C(n-d,k+e)\) for unsettled \((d,e)\) |
| `collide` | Finite \(m\)-slice of \(C(n,k)=C(m,l)\) |
| `sanity` | Catalog, classifier, Lucas, image, 3003 tripwire |

```text
python singmaster_intersect.py sanity
python singmaster_intersect.py modular --imin 9 --imax 9 --kextra 80
```

Tests: `singmaster_intersect.py sanity` covers the engine,
`scripts/test_sweeps.py` covers the sweep drivers' certificate logic,
`scripts/test_witness.py` covers the certificate verifier,
`scripts/test_sizelaw.py` covers the size law and the escalation trigger, and
`scripts/test_kernel.py` covers the image-scan kernel.

```text
python scripts/test_sweeps.py
```

## Certificates

A sweep proves *every* extra column has a killing prime. `scripts/witness.py`
records which prime, per column, and re-checks it independently.

| Command | Question |
|---|---|
| `witness build --i N` | Rebuild the witness table from a run's checkpoint |
| `witness verify --file F` | Re-check every certificate, and that none is missing |
| `witness one --N .. --K .. --k .. --p ..` | Check a single certificate |

```text
python scripts/witness.py verify --file results/i7_witness.npz --sample 4000
```

The verifier shares no code path with the sweep: it re-derives
\(r(p)=C(N,K) mod p\) by Lucas from \(N,K\) alone and walks the column
image using \(C(k+j,k)=(-1)^jC(p-k-1,j)\) — two modular multiplications per
step, no factorial table, no inverses, **O(1) memory**. So the builder is
untrusted: if it proposes the wrong prime, verification fails. One Band II
certificate checks in about a second from five integers.

`verify` reports two independent things, and needs both: no certificate is
invalid, **and** no claimed column is missing one. A clean `family_sweep` writes
`results/i{N}_witness.npz` automatically and names it in the certificate.

`impossible` + a witness prime is a column-level proof. `possible` only means the prime list did not kill it. Nearby nulls are **sampled**, not exhaustive.

## Results worth keeping

| File | Claim |
|---|---|
| `results/fibonacci_i8_k300.json` | Exact: i=8 has no extra left-half with \(2\le k\le 300\), no central. \(N=6\). ~16.3 h |
| `results/modular_i8_k400.json` | Modular: i=8, \(k=2..400\) all impossible |
| `results/modular_i9_k80.json` | Modular: i=9, \(k=2..80\) all impossible |
| `results/modular_i10_k20.json` | Modular: i=10, \(k=2..20\) all impossible |
| `results/nearby_k2M-8M_de8.json` | Sampled nearby, \(k=2\text{M}..8\text{M}\), \(d,e\le 8\), `new_hits=0` |
| `results/nextprime_i8_k100001-1000000_summary.json` | Stage 3 census summary (i=8, \(k=10^5..10^6\)). Not a theorem |
| `results/stragglers_nearK.json` | 25 near-\(K\) Band I \(k\) all killed at \(p>N/2\) |
| `results/triple_hunt_p1e6-K.json` | NONE-window image runs, \(P_\mathrm{hi}>10^6\): max run 2, 0 triples (fat cells not exhaustive) |
| `results/fat_image_hunt.json` | Fat NONE+PART-lower, first 3 primes: 369 triples, hunt cap max run 3 |
| `results/walk_369.json` | Those 369 walked past prime 3: max run 6, counts match size law |
| `results/bandii_sweep.json` | Band II \(p>N/2\): 1,055,989 columns, all killed by prime 8 |
| `results/zjump.json` | Band I remnant Z-jump: 3,215,816 columns, 0 anomalies, tail = walk-369 |
| `results/i7_sweep.json` | i=7 Band II + Z-jump: all extra \(k\) killed, \(N=6\) (110.7 s) |
| `results/i6_sweep.json` | i=6 same pipeline, \(N=6\) (4.7 s) |
| `results/i5_sweep.json` | i=5 same pipeline, \(N=6\) (1.6 s) |
| `results/i4_sweep.json` | i=4 same pipeline, \(N=6\) (0.9 s) |
| `results/i3_sweep.json` | i=3 same pipeline, \(N=6\) (0.8 s) |
| `results/i2_sweep.json` | i=2 Band II only (exact already covered Band I), \(N=6\) |

Do not build \(m_{10}\) (~147 million digits) to answer small-\(k\) extra-rep. Modular already killed those columns.

## i=8 extra columns

Every extra \(k\in[2,k_{\max}]\setminus\{K,K+1\}\) has an unconditional
modular certificate. \(N(C(F_{18}F_{19},F_{16}F_{19}))=6\). Not Singmaster.

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
(\(E_r<10^{-2}\)), or when far more survive than the law allows (Poisson
upper tail \(<10^{-3}\)).

```text
python scripts/sizelaw.py run --i 9 --k 11
python scripts/sizelaw.py predict --check
```

The law has no fitted parameters: \(|I_{p,k}|\) is the proved involution
bound (\(g+1\) for odd \(k\), \(\lceil g/2
ceil+1\) for even) minus
birthday collisions. It reproduces the recorded Band II pre-registration and
the i=7 Z-jump round by round.


## Profiling

Four performance claims on this project survived reasoning and died on
measurement, so the pipeline now carries a profiler rather than an argument:

```text
python scripts/profile_sweep.py --i 9
```

It samples each phase and reports the split between factorial-table
construction, the image scan, r(p), and prime assignment. The split is the
part to trust (wall-clock projection is good to about a factor of two), and
it is decisive: at i=9 Band II is 100% scan, while the Z-jump is **90.6%
table**. A faster scan does nothing for the second.

The distinct-prime count integrates prime density over the live intervals
rather than sampling a window near the range start, which is biased low at
small k; the integrated form recovers i=8's recorded 124,830 primes to 0.2%.

## Checkpoints

Every run's jsonl opens with a schema header pinning the version and the
run's parameters (`i`, `N`, `K`, `k_max`, `k_lo_z`, caps). Resume verifies it
and refuses on mismatch, or on a checkpoint with no header at all. Resume
merges old records with new ones, so a silent format or parameter change
would produce a certificate over columns that were never all tested the same
way.

## Notes

- [`docs/modular-spec.txt`](docs/modular-spec.txt) — Lucas/modular layer: what a certificate is, how to run scans, what not to rebuild
- [`docs/interior-2022.md`](docs/interior-2022.md) — the MRSTT interior theorem, and why it does not reach i=8/9
- [`docs/open-questions.md`](docs/open-questions.md) — ranked open questions, with difficulty and recommended effort
- [`docs/campaign-log.txt`](docs/campaign-log.txt) — original search campaign (exact i=8 \(k\le 300\), nearby, collide)
- [`scripts/run-i8-k300.bat`](scripts/run-i8-k300.bat) — historical wrapper for the exact i=8 \(k\le 300\) job (already finished)
- Settled theorem pairs \((k,l)\) and nearby \((1,1),(1,2),(2,1)\) are not re-sieved
- Classifier uses `is_fibonacci_pair`, not “rows differ by 1”
- The image scan uses the involution \((k-1-x)_k=(-1)^k(x)_k\) to test only
  \(b\in[0,\lceil g/2
ceil)\): even \(k\) repeats on the upper half, odd
  \(k\) negates, and the \(-s\) branch is a subtraction rather than a second
  multiply. Exact, ~2.2x measured, byte-identical output including the
  witness index \(b\). `scan_ks_full` is kept as the reference;
  `USE_HALF_SCAN=False` switches back for A/B
- A scan with \(r(p)=0\) is refused, not run: 0 lies in every column image,
  so such a prime certifies nothing, and a mask scan would silently report
  every column killed
- `_column_possible_scan` walks the image as \(\{(-1)^jC(g-1,j)\}\) keeping
  \(A_j=(g-1)_j\) and \(B_j=j!\) — two multiplications per step, no modular
  inverse — and folds by the same involution as the sweep kernel, so only
  \(j<\lceil g/2
ceil\) is tested. 4-7x on a large-\(g\) membership test.
  `_column_possible_scan_ref` is the superseded version, kept because sanity
  checks the two agree
- Lucas digits use `binom_mod_prime`, not `math.comb`. Large-\(p\) column tests scan and do not cache residue images
- `nextprime_sweep` memoises \(r(p)=C(N,K) mod p\) per prime (`RCache`), not per
  \((k,p)\) pair: it does not depend on \(k\), and consecutive columns share their
  next prime. 7.3x on the Stage-2 range, 13-15x on Stage-3-style \(k\). The table
  is pruned of \(p\le k\) every 1000 columns, so it stays the width of the walk
  window. Rows are identical with it on or off, and `test_sweeps.py` pins that
- Shared Band I / cell-geometry helpers (`cells`, `live_intervals`, `first_live_after`,
  `image_j`, `inv_table`, `binom_mod_prime`, `chunk_ks`, jsonl I/O) live in
  [`scripts/bandii_kernel.py`](scripts/bandii_kernel.py). They were copy-pasted into six
  scripts and had begun to drift; do not re-inline them
- Pre-flight guards use `check()`, not `assert` — `python -O` strips `assert`, and every
  certificate downstream of a pre-flight is claimed unconditional
- A sweep is `clean` only if nothing survived **and** every column was actually testable.
  A column with no live prime left is an anomaly, not a kill

