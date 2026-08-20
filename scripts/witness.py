#!/usr/bin/env python3
"""Per-column kill certificates: build them, then check them independently.

The sweeps prove "every extra column k has a prime p with r(p) not in I_{p,k}"
but record only aggregate pass counts. The witness primes exist transiently in
the run's jsonl -- a column absent from pass j+1's survivors died at pass j+1's
prime -- and that file is gitignored. So until now the headline certificate was
re-checkable only by re-running the whole sweep.

This module closes that. Two halves, deliberately separated:

  build   reconstructs (k, p) pairs from a run's jsonl. UNTRUSTED. It knows the
          sweep's prime-selection rule, so a bug here proposes a wrong p.
  verify  re-derives r(p) = C(N,K) mod p by Lucas from N and K alone, then
          walks the image of column k. TRUSTED. It shares no code path with
          the sweep: no factorial table, no numpy, no gmpy2, no image set.

A wrong builder therefore cannot manufacture a certificate -- the verifier
rejects it. Coverage is checked separately: the witness set must equal the
claimed column set exactly, or the theorem has a hole rather than a bad row.

check_witness establishes four things, in order:
  (i)   p is prime (its own Miller-Rabin) and p > k >= 2 -- both Lucas and the
        claim that I_{p,k} is the COMPLETE image assume primality, so a
        certificate naming a composite modulus proves nothing;
  (ii)  r really is m mod p, recomputed from N and K by every route that
        applies to this p -- full Lucas always, two-digit when p^2 > N, the
        delta closed form on the alpha=1 live window -- all agreeing;
  (iii) r != 0, since 0 lies in every column image and a dead prime certifies
        nothing;
  (iv)  r is outside I_{p,k}, by walking the whole image.

The walk uses C(k+j,k) = (-1)^j C(p-k-1,j) (mod p), maintained as A_j = (g-1)_j
and B_j = j!: two modular multiplications per step, no inverses, no table,
O(1) memory. One Band II certificate checks in well under a second on a laptop
from (N, K, k, p) alone.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

SCHEMA = 1


# ---------------------------------------------------------------------------
# the trusted half: pure Python, sharing no code path with the sweep
# ---------------------------------------------------------------------------

def binom_mod_prime_pure(n: int, k: int, p: int) -> int:
    """C(n,k) mod p for prime p, 0 <= k <= n < p. One inverse, not k of them."""
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    if k == 0:
        return 1
    num = den = 1
    for i in range(k):
        num = num * (n - i) % p
        den = den * (i + 1) % p
    return num * pow(den, -1, p) % p


def lucas_mod_pure(n: int, k: int, p: int) -> int:
    """C(n,k) mod p by Lucas. Independent of every r(p) form the sweep uses."""
    if k < 0 or n < 0 or k > n:
        return 0
    result = 1
    nn, kk = int(n), int(k)
    while nn > 0 or kk > 0:
        ni, ki = nn % p, kk % p
        if ki > ni:
            return 0
        result = result * binom_mod_prime_pure(ni, ki, p) % p
        nn //= p
        kk //= p
    return result


def image_hit_tablefree(p: int, r: int, k: int) -> int | None:
    """Least j with r = (-1)^j C(g-1,j) mod p, else None. O(g) time, O(1) memory.

    I_{p,k} = {C(x,k) mod p : x} = {0} u {(-1)^j C(g-1,j) : 0 <= j < g}, g = p-k.
    Returns None exactly when r lies outside the image, i.e. exactly when p
    certifies column k impossible.
    """
    g = p - int(k)
    if g <= 0:
        return None
    rp = int(r) % p
    if rp == 0:
        return None
    A = B = 1
    sign = 1
    for j in range(g):
        if A == (sign * rp * B) % p:
            return j
        A = A * ((g - 1 - j) % p) % p
        B = B * ((j + 1) % p) % p
        sign = -sign
    return None


_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime_pure(n: int) -> bool:
    """Deterministic Miller-Rabin. Exact for n < 3.3e24, well past any p here.

    Not optional. Both legs of the argument -- Lucas for r(p), and the claim
    that I_{p,k} is the complete image of column k -- assume p prime. A
    certificate naming a composite modulus proves nothing, so the verifier
    must establish primality rather than take it on trust from the builder.
    """
    n = int(n)
    if n < 2:
        return False
    for q in _MR_BASES:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def r_two_digit_pure(N: int, K: int, p: int) -> int | None:
    """r(p) by two-digit Lucas. Valid only when p*p > N; else None.

    Below sqrt(N) the product C(alpha,beta) C(n0,k0) is simply false -- N has
    three or more base-p digits and the formula drops one.
    """
    if p * p <= N:
        return None
    a, b = N // p, K // p
    n0, k0 = N - a * p, K - b * p
    if not (0 <= b <= a and 0 <= k0 <= n0):
        return 0
    return binom_mod_prime_pure(a, b, p) * binom_mod_prime_pure(n0, k0, p) % p


def r_delta_pure(N: int, K: int, p: int) -> int | None:
    """r(p) by the delta closed form, for the alpha=1 live window; else None.

    For N/2 < p <= N one has alpha=1, beta=0, n0 = N-p = p-delta with
    delta = 2p-N, and C(n0,K) = (-1)^K C(K+delta-1, delta-1) mod p. Cheap when
    delta is small, and a genuinely different arithmetic route to the same r.
    """
    if not (N // 2 < p <= N) or K >= p:
        return None
    dlt = 2 * p - N
    if dlt < 1:
        return None
    c = binom_mod_prime_pure(K + dlt - 1, dlt - 1, p)
    return (-c) % p if K % 2 else c


def check_witness(N: int, K: int, k: int, p: int) -> dict:
    """Is p a valid impossibility certificate for column k of m = C(N,K)?

    Self-contained: five integers in, verdict out. This is the whole
    referee-facing surface of the project. It establishes, in order:
      (i)   p is prime and p > k >= 2;
      (ii)  r really is m mod p, recomputed from N and K by every route that
            applies to this p, all agreeing;
      (iii) r != 0, since a dead prime certifies nothing;
      (iv)  r is outside I_{p,k}, by walking the whole image.
    """
    k, p = int(k), int(p)
    if k < 2:
        return {"ok": False, "reason": "k < 2"}
    if p <= k:
        return {"ok": False, "reason": f"p={p} <= k={k}; the image test needs p > k"}
    if not is_prime_pure(p):
        return {"ok": False, "reason": f"p={p} is not prime; Lucas and I_p,k both need it"}
    r = lucas_mod_pure(N, K, p)
    routes = {"lucas": r}
    for name, fn in (("two_digit", r_two_digit_pure), ("delta", r_delta_pure)):
        alt = fn(N, K, p)
        if alt is not None:
            routes[name] = alt
    if len(set(routes.values())) != 1:
        return {"ok": False, "reason": f"r(p) routes disagree: {routes}", "routes": routes}
    if r == 0:
        # 0 = C(x,k) for every x < k, so 0 is always in the image.
        return {"ok": False, "reason": "r(p)=0 certifies nothing", "r": 0}
    j = image_hit_tablefree(p, r, k)
    if j is not None:
        return {"ok": False, "reason": "r(p) IS in the image", "r": r, "j": j}
    return {"ok": True, "r": r, "g": p - k, "routes": sorted(routes)}


# ---------------------------------------------------------------------------
# the untrusted half: reconstruct (k, p) from a run's checkpoint
# ---------------------------------------------------------------------------

def read_jsonl(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def _survivors_of(recs: list[dict]) -> dict[int, int]:
    """k -> last prime tested, across a set of chunk records."""
    out: dict[int, int] = {}
    for rec in recs:
        p = int(rec["p"])
        for s in rec.get("survivors") or []:
            k = int(s["k"])
            out[k] = int(s["k"]) + int(s["g"]) if "g" in s else p
    return out


def build_bandii(rows, k_lo, k_hi, primes, pass_of):
    """Band II: pass j tests every still-alive column against primes[j-1]."""
    alive = set(range(k_lo, k_hi + 1))
    witness: dict[int, int] = {}
    for j, p in enumerate(primes, start=1):
        recs = [r for r in rows if pass_of(r) == j]
        if not recs:
            break
        survivors = set(_survivors_of(recs))
        for k in alive - survivors:
            witness[k] = p
        alive &= survivors
        if not alive:
            break
    return witness, alive


def build_zjump(rows, k_lo, k_hi, ivs, d, round_of):
    """Z-jump: each column jumps to the first LIVE prime above its last one."""
    from bandii_kernel import first_live_after

    last_p = {k: k for k in range(k_lo, k_hi + 1)}
    witness: dict[int, int] = {}
    untestable: set[int] = set()
    rnd = 1
    while last_p:
        recs = [r for r in rows if round_of(r) == rnd]
        if not recs:
            break
        survivors = _survivors_of(recs)
        seen = {int(r["p"]) for r in recs}
        for k, lp in last_p.items():
            if k in survivors:
                continue
            p = first_live_after(lp, ivs, d)
            if p is None:
                untestable.add(k)
                continue
            if p not in seen:
                raise RuntimeError(
                    f"rebuilt prime {p} for k={k} in round {rnd} is absent from "
                    f"the checkpoint's prime set: builder and run disagree"
                )
            witness[k] = p
        last_p = dict(survivors)
        rnd += 1
    return witness, set(last_p) | untestable


# ---------------------------------------------------------------------------
# store
# ---------------------------------------------------------------------------

def save(path: Path, meta: dict, witness: dict[int, int]) -> dict:
    import numpy as np

    ks = np.array(sorted(witness), dtype=np.int64)
    ps = np.array([witness[int(k)] for k in ks], dtype=np.int64)
    digest = hashlib.sha256(ks.tobytes() + ps.tobytes()).hexdigest()
    meta = {**meta, "schema": SCHEMA, "n_witnesses": int(ks.size), "sha256": digest}
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, k=ks, p=ps, meta=json.dumps(meta))
    return meta


def load(path: Path):
    import numpy as np

    z = np.load(path, allow_pickle=False)
    ks, ps = z["k"], z["p"]
    meta = json.loads(str(z["meta"]))
    digest = hashlib.sha256(ks.tobytes() + ps.tobytes()).hexdigest()
    if meta.get("sha256") != digest:
        raise RuntimeError(f"{path}: witness table does not match its own sha256")
    return ks, ps, meta


# ---------------------------------------------------------------------------
# coverage: a witness per row is not enough, every claimed row needs one
# ---------------------------------------------------------------------------

def coverage(ks, meta: dict) -> dict:
    expected: set[int] = set()
    for lo, hi in meta["claimed_ranges"]:
        expected |= set(range(int(lo), int(hi) + 1))
    expected -= {int(x) for x in meta.get("excluded", [])}
    got = {int(x) for x in ks}
    missing = sorted(expected - got)
    extra = sorted(got - expected)
    return {
        "n_expected": len(expected),
        "n_witnessed": len(got),
        "n_missing": len(missing),
        "missing_sample": missing[:20],
        "n_extra": len(extra),
        "extra_sample": extra[:20],
        "complete": not missing and not extra,
    }


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------

def _check_chunk(args):
    N, K, pairs = args
    bad = []
    for k, p in pairs:
        res = check_witness(N, K, k, p)
        if not res["ok"]:
            bad.append({"k": int(k), "p": int(p), **res})
    return len(pairs), bad


def verify(path: Path, sample: int | None, workers: int, seed: int = 0) -> dict:
    import numpy as np

    ks, ps, meta = load(path)
    N, K = int(meta["N"]), int(meta["K"])
    cov = coverage(ks, meta)

    idx = np.arange(ks.size)
    if sample is not None and sample < ks.size:
        idx = np.random.default_rng(seed).choice(ks.size, size=sample, replace=False)
        idx.sort()
    pairs = [(int(ks[i]), int(ps[i])) for i in idx]
    total_g = sum(p - k for k, p in pairs)

    print(
        f"  verifying {len(pairs)} of {ks.size} certificates  "
        f"(sum g = {total_g:,} steps)",
        flush=True,
    )
    t0 = time.time()
    bad: list[dict] = []
    done = 0
    if workers <= 1 or len(pairs) < 64:
        n, b = _check_chunk((N, K, pairs))
        done, bad = n, b
    else:
        import multiprocessing as mp

        size = max(1, len(pairs) // (workers * 8))
        tasks = [
            (N, K, pairs[i : i + size]) for i in range(0, len(pairs), size)
        ]
        ctx = mp.get_context("spawn")
        with ctx.Pool(workers) as pool:
            for n, b in pool.imap_unordered(_check_chunk, tasks):
                done += n
                bad.extend(b)
                if done % max(1, len(pairs) // 20) < size:
                    print(
                        f"    {done}/{len(pairs)}  bad={len(bad)}  "
                        f"{time.time()-t0:.0f}s",
                        flush=True,
                    )
    dt = time.time() - t0
    return {
        "search": "witness_verify",
        "file": str(path),
        "N": N,
        "K": K,
        "source": meta.get("source"),
        "n_total": int(ks.size),
        "n_checked": len(pairs),
        "sampled": sample is not None and sample < ks.size,
        "sum_g": total_g,
        "coverage": cov,
        "n_invalid": len(bad),
        "invalid": bad[:50],
        "valid": not bad and cov["complete"],
        "seconds": round(dt, 3),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def fill_small_gaps(path: Path, i: int) -> dict:
    """Add engine-derived witnesses for any column the table is missing.

    The sweep only ever claimed [k_lo_z, K-1] u [K+2, k_max]; columns below the
    Z-jump's start were closed by the engine's modular scan and lived in a
    separate file. This puts them in the table, so one artifact carries the
    whole claim. It re-runs no sweep -- a witness is checkable from (N,K,k,p)
    regardless of how it was found, which is exactly the property that makes
    this legitimate rather than a shortcut.
    """
    from singmaster_intersect import obstructing_prime, primes_upto

    from bandii_kernel import kmax_of, make_fam

    ks, ps, meta = load(path)
    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    have = {int(a): int(b) for a, b in zip(ks, ps)}
    want = set(range(2, kmax + 1)) - {fam.K, fam.K + 1}
    missing = sorted(want - set(have))
    if not missing:
        return {"i": i, "added": 0, "unresolved": [], "sha256": meta["sha256"]}

    small_primes = primes_upto(50_000)
    added, unresolved = 0, []
    for k in missing:
        q = obstructing_prime(fam.N, fam.K, k, small_primes)
        if q is None:
            unresolved.append(k)
        else:
            have[k] = int(q)
            added += 1
    meta = {**meta, "claimed_ranges": [[2, kmax]], "excluded": [fam.K, fam.K + 1],
            "n_filled_from_engine": added}
    out = save(path, meta, have)
    return {"i": i, "added": added, "unresolved": unresolved,
            "sha256": out["sha256"], "n_witnesses": out["n_witnesses"]}


def _pass_tag(prefix: str):
    """Pass/round index out of a checkpoint record, across the three writers.

    family_sweep tags chunks "bii3"/"z2"; bandii_sweep uses prime_index;
    zjump uses round. Same reconstruction either way.
    """

    def of(rec: dict):
        if prefix == "bii" and "prime_index" in rec:
            return int(rec["prime_index"])
        if prefix == "z" and "round" in rec and "p" in rec:
            return int(rec["round"])
        tag = str(rec.get("tag", ""))
        if tag.startswith(prefix):
            rest = tag[len(prefix):]
            if rest.isdigit():
                return int(rest)
        return None

    return of


def build_i8(chk_bandii: Path, chk_zjump: Path, out: Path) -> dict:
    """i=8's dedicated pipeline: bandii_sweep.jsonl + zjump.jsonl."""
    from bandii_kernel import CAP, D, K, KMAX, KMIN, N, PRIMES, cells, live_intervals
    from zjump import HANG_RUNS, K_HI_NEW, K_LO_NEW

    witness: dict[int, int] = {}
    unresolved: set[int] = set()
    claimed: list[list[int]] = []

    if chk_bandii.exists():
        w, alive = build_bandii(
            read_jsonl(chk_bandii), KMIN, KMAX, PRIMES[:CAP], _pass_tag("bii")
        )
        witness.update(w)
        unresolved |= alive
        claimed.append([KMIN, KMAX])

    if chk_zjump.exists():
        ivs = live_intervals(cells())
        rows = read_jsonl(chk_zjump)
        for lo, hi in HANG_RUNS + [(K_LO_NEW, K_HI_NEW)]:
            w, alive = build_zjump(rows, lo, hi, ivs, D, _pass_tag("z"))
            witness.update(w)
            unresolved |= alive
            claimed.append([lo, hi])

    if not claimed:
        raise RuntimeError(
            f"no i=8 checkpoints found ({chk_bandii.name}, {chk_zjump.name}). "
            "These jsonl files are the only record of which prime killed which "
            "column; without them the sweeps must be repeated."
        )
    meta = {
        "source": f"{chk_bandii.name}+{chk_zjump.name}",
        "i": 8,
        "N": N,
        "K": K,
        "k_max": KMAX,
        "claimed_ranges": claimed,
        "excluded": [],
        "n_unresolved": len(unresolved),
        "unresolved": sorted(unresolved)[:100],
        "primes_bii": PRIMES[:CAP],
    }
    return save(out, meta, witness)


def _build_family(i: int, chk: Path, out: Path) -> dict:
    """Build from a family_sweep / bandii_sweep / zjump checkpoint."""
    from bandii_kernel import (
        cells,
        first_primes_above,
        kmax_of,
        live_intervals,
        make_fam,
    )
    from family_sweep import CAP_BII, K_EXACT

    fam = make_fam(i)
    kmax, _ = kmax_of(fam)
    rows = read_jsonl(chk)
    primes = first_primes_above(fam.N2, fam.D, kmax, n=max(CAP_BII, 16))[:CAP_BII]

    w_bii, alive_bii = build_bandii(
        rows,
        fam.K + 2,
        kmax,
        primes,
        _pass_tag("bii"),
    )
    k_lo_z = K_EXACT.get(i, 2) + 1
    if k_lo_z < fam.K:
        ivs = live_intervals(cells(fam), fam)
        w_z, alive_z = build_zjump(
            rows,
            k_lo_z,
            fam.K - 1,
            ivs,
            fam.D,
            _pass_tag("z"),
        )
    else:
        w_z, alive_z = {}, set()

    # Columns below the Z-jump's start are closed by the engine's modular
    # scan, not by the sweep, so the witness table used to have a hole there
    # ([2,200] for i=3..7, [2,80] for i=9, k=2 for i=8) and the full claim was
    # only assembled by reading a second file. Extending the Z-jump down is NOT
    # the fix: the sweep kernel has no O(1) route for k=2, so k=2 there costs a
    # full g~p scan (410 ms at i=8). The engine does have one, so take the
    # witnesses from there and make the table self-contained.
    w_small: dict[int, int] = {}
    unresolved_small: set[int] = set()
    if k_lo_z > 2:
        from singmaster_intersect import obstructing_prime, primes_upto

        small_primes = primes_upto(50_000)
        for kk in range(2, k_lo_z):
            if kk in (fam.K, fam.K + 1):
                continue
            q = obstructing_prime(fam.N, fam.K, kk, small_primes)
            if q is None:
                unresolved_small.add(kk)
            else:
                w_small[kk] = q

    witness = {**w_bii, **w_z, **w_small}
    claimed = [[fam.K + 2, kmax]]
    if k_lo_z < fam.K:
        claimed.append([k_lo_z, fam.K - 1])
    if w_small or unresolved_small:
        claimed.append([2, k_lo_z - 1])
    meta = {
        "source": str(chk.name),
        "i": i,
        "N": fam.N,
        "K": fam.K,
        "k_max": kmax,
        "claimed_ranges": claimed,
        "excluded": [],
        "n_unresolved": len(alive_bii) + len(alive_z) + len(unresolved_small),
        "unresolved": sorted(alive_bii | alive_z | unresolved_small)[:100],
        "n_small_k_from_engine": len(w_small),
        "primes_bii": primes,
    }
    return save(out, meta, witness)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="reconstruct witnesses from a run's jsonl")
    b.add_argument("--i", type=int, required=True)
    b.add_argument("--checkpoint", type=Path, default=None)
    b.add_argument("--out", type=Path, default=None)

    v = sub.add_parser("verify", help="independently re-check a witness table")
    v.add_argument("--file", type=Path, required=True)
    v.add_argument("--sample", type=int, default=None,
                   help="check a random subset; omit to check every certificate")
    v.add_argument("--workers", type=int, default=8)
    v.add_argument("--seed", type=int, default=0)
    v.add_argument("--json_out", type=Path, default=None)

    fl = sub.add_parser("fill", help="add engine witnesses for columns the table lacks")
    fl.add_argument("--i", type=int, required=True)
    fl.add_argument("--file", type=Path, default=None)

    o = sub.add_parser("one", help="check a single certificate from five integers")
    o.add_argument("--N", type=int, required=True)
    o.add_argument("--K", type=int, required=True)
    o.add_argument("--k", type=int, required=True)
    o.add_argument("--p", type=int, required=True)

    args = ap.parse_args()

    if args.cmd == "one":
        t0 = time.time()
        res = check_witness(args.N, args.K, args.k, args.p)
        print(json.dumps({**res, "seconds": round(time.time() - t0, 3)}, indent=2))
        return 0 if res["ok"] else 1

    if args.cmd == "fill":
        f = args.file or ROOT / "results" / f"i{args.i}_witness.npz"
        res = fill_small_gaps(f, args.i)
        print(f"  {f.name}: added {res['added']} engine witnesses"
              + (f", UNRESOLVED {res['unresolved'][:10]}" if res["unresolved"] else ""))
        print(f"  sha256 {res['sha256'][:32]}...")
        return 0 if not res["unresolved"] else 2

    if args.cmd == "build":
        out = args.out or ROOT / "results" / f"i{args.i}_witness.npz"
        if args.i == 8 and args.checkpoint is None:
            meta = build_i8(
                ROOT / "results" / "bandii_sweep.jsonl",
                ROOT / "results" / "zjump.jsonl",
                out,
            )
            print(f"wrote {out}", flush=True)
            print(f"  witnesses  {meta['n_witnesses']}", flush=True)
            print(f"  unresolved {meta['n_unresolved']}", flush=True)
            print(f"  sha256     {meta['sha256'][:32]}...", flush=True)
            return 0 if meta["n_unresolved"] == 0 else 2
        chk = args.checkpoint or ROOT / "results" / f"i{args.i}_sweep.jsonl"
        if not chk.exists():
            print(f"missing checkpoint {chk}", flush=True)
            print("  the jsonl is the ONLY record of which prime killed which", flush=True)
            print("  column. Without it the run must be repeated.", flush=True)
            return 1
        meta = _build_family(args.i, chk, out)
        print(f"wrote {out}", flush=True)
        print(f"  witnesses  {meta['n_witnesses']}", flush=True)
        print(f"  unresolved {meta['n_unresolved']}", flush=True)
        print(f"  sha256     {meta['sha256'][:32]}...", flush=True)
        return 0 if meta["n_unresolved"] == 0 else 2

    rep = verify(args.file, args.sample, args.workers, args.seed)
    cov = rep["coverage"]
    print(flush=True)
    print(f"  columns claimed   {cov['n_expected']}", flush=True)
    print(f"  witnesses present {cov['n_witnessed']}   missing={cov['n_missing']} "
          f"extra={cov['n_extra']}", flush=True)
    print(f"  certificates checked {rep['n_checked']}   invalid={rep['n_invalid']}"
          f"   {rep['seconds']}s", flush=True)
    if rep["invalid"]:
        for badrow in rep["invalid"][:10]:
            print(f"    INVALID k={badrow['k']} p={badrow['p']}: {badrow['reason']}",
                  flush=True)
    if cov["missing_sample"]:
        print(f"    MISSING {cov['missing_sample']}", flush=True)
    verdict = "VALID" if rep["valid"] else "NOT VALID"
    if rep["sampled"] and rep["valid"]:
        verdict = "VALID (coverage complete; certificates sampled)"
    print(f"  RESULT {verdict}", flush=True)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}", flush=True)
    return 0 if rep["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
