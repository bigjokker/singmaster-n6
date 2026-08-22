#!/usr/bin/env python3
"""Band II image kernel: factorial table, r(p), survive-test.

Scans use the involution half-range (see scan_ks_half); scan_ks_full is the
reference implementation and USE_HALF_SCAN switches between them.

Column k survives prime p iff exists b in [0, p-k-1] with
    F[k+b] ≡ r(p)·F[k]·F[b]  (mod p)
which is C(k+b, k) ≡ r(p) (mod p). b is the fat-image j.

Do not build Pascal rows. Do not cache F across primes.
int64 only: s·F[b] < p² < 2^63.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass

import gmpy2

try:
    import numpy as np
except ImportError:  # pragma: no cover
    # Only the factorial-table scanners (fact_table / scan_ks) need numpy.
    # The pure-integer helpers below must stay importable without it.
    np = None

N = 10_803_704
K = 4_126_647
D = 6_677_057
N2 = 5_401_852
KMAX = 5_182_637
KMIN = K + 2  # 4126649
NCOLS = KMAX - KMIN + 1  # 1055989
LOG10_M = 3120255.2212


@dataclass(frozen=True)
class Fam:
    i: int
    N: int
    K: int

    @property
    def D(self) -> int:
        return self.N - self.K

    @property
    def N2(self) -> int:
        return self.N // 2

    @property
    def p_two(self) -> int:
        return int(math.isqrt(self.N)) + 1

    @property
    def K1(self) -> int:
        return self.K + 1


def make_fam(i: int) -> Fam:
    import gmpy2

    n = int(gmpy2.fib(2 * i + 2) * gmpy2.fib(2 * i + 3))
    k = int(gmpy2.fib(2 * i) * gmpy2.fib(2 * i + 3))
    return Fam(i=i, N=n, K=k)


def kmax_of(fam: Fam) -> tuple[int, float]:
    """Largest k with C(2k,k) <= C(N,K), via lgamma. Also log10 m."""
    logm = (
        math.lgamma(fam.N + 1)
        - math.lgamma(fam.K + 1)
        - math.lgamma(fam.N - fam.K + 1)
    ) / math.log(10)

    def log_central(k: int) -> float:
        return (math.lgamma(2 * k + 1) - 2 * math.lgamma(k + 1)) / math.log(10)

    lo, hi = 2, fam.N
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if log_central(mid) <= logm:
            lo = mid
        else:
            hi = mid - 1
    return lo, logm

PRIMES = [
    5_401_853,
    5_401_861,
    5_401_867,
    5_401_897,
    5_401_901,
    5_401_951,
    5_401_969,
    5_401_973,
    5_401_993,
    5_401_999,
    5_402_003,
    5_402_011,
    5_402_021,
    5_402_051,
    5_402_057,
    5_402_063,
]
CAP = 14
P1 = PRIMES[0]
P2 = PRIMES[1]
R_P1 = 1_275_205

STRAGGLERS = {
    4_126_624: 273_671,
    4_126_638: 268_500,
    4_126_642: 2_006,
    4_126_643: 554_555,
}


def delta(p: int) -> int:
    return 2 * p - N


def fact_table(p: int) -> np.ndarray:
    """F[i] = i! mod p, i = 0..p-1. Python loop, ~0.3 s at p1."""
    F = np.empty(p, dtype=np.int64)
    F[0] = 1
    acc = 1
    for i in range(1, p):
        acc = acc * i % p
        F[i] = acc
    return F


def r_from_F(F: np.ndarray, p: int, *, N: int = N, K: int = K) -> int:
    n0 = N - p
    return (
        int(F[n0])
        * pow(int(F[K]), -1, p)
        * pow(int(F[n0 - K]), -1, p)
        % p
    )


def r_falling(p: int, *, N: int = N, K: int = K) -> int:
    n0 = N - p
    kk = n0 - K
    num = 1
    for i in range(kk):
        num = num * (n0 - i) % p
    den = 1
    for i in range(1, kk + 1):
        den = den * i % p
    return num * pow(den, -1, p) % p


def r_closed(p: int, *, N: int = N, K: int = K) -> int:
    """(-1)^K C(K+δ-1, δ-1) mod p, δ = 2p-N. Lower index is δ-1."""
    dlt = 2 * p - N
    n = K + dlt - 1
    kk = dlt - 1
    num = den = 1
    for i in range(kk):
        num = num * (n - i) % p
        den = den * (i + 1) % p
    c = num * pow(den, -1, p) % p if kk else 1
    if K % 2:
        c = (-c) % p
    return c


def r_checked(F: np.ndarray, p: int, falling: bool = False, *, N: int = N, K: int = K) -> int:
    ra = r_from_F(F, p, N=N, K=K)
    rc = r_closed(p, N=N, K=K)
    if ra != rc:
        raise RuntimeError(f"r(p) table {ra} != closed {rc} at p={p}")
    if falling:
        rb = r_falling(p, N=N, K=K)
        if ra != rb:
            raise RuntimeError(f"r(p) table {ra} != falling {rb} at p={p}")
    return ra


# Half-scan is exact, not an approximation, so it is on by default. Flip to
# False only to A/B against the reference implementation.
USE_HALF_SCAN = True

# Build only the O(g) factorial entries the scan reads, instead of all p of
# them. Exact; flip to False to A/B against the table-based path.
USE_WINDOWED_SCAN = True


def _check_r(r: int, p: int) -> int:
    """r=0 would make the scan report every column killed, which is false.

    0 = C(x,k) for every x < k, so 0 lies in every column image and a prime
    with r(p)=0 certifies nothing. The mask-based scans cannot see that (the
    factorial table has no zero entry), so they would silently return "no
    survivors" -- a false certificate. Refuse instead.
    """
    rp = int(r) % p
    if rp == 0:
        raise ValueError(
            f"scan called with r(p)=0 at p={p}: 0 is in every column image, "
            f"so this prime kills nothing. Do not scan it."
        )
    return rp


def scan_ks_full(F: np.ndarray, p: int, r: int, ks) -> list[dict]:
    """Reference scan: test every b in [0, g). Kept to check the half-scan."""
    out: list[dict] = []
    p64 = np.int64(p)
    rp = _check_r(r, p)
    for k in ks:
        k = int(k)
        s = rp * int(F[k]) % p
        n = p - k
        left = F[k:p]
        right = (np.int64(s) * F[:n]) % p64
        eq = left == right
        hit = bool(eq.any())
        if hit:
            b = int(eq.argmax())
            g = n
            # g_even is kept for the recorded jsonl format and the
            # pre-registration tables; k_odd is the same bit, but it is the
            # INVARIANT one -- a property of the polynomial f_k, fixed once,
            # where g changes with every prime. (p odd => g even <=> k odd.)
            out.append({"k": k, "g": g, "g_even": g % 2 == 0,
                        "k_odd": k % 2 == 1, "b": b})
        del eq
    return out


def scan_ks_half(F: np.ndarray, p: int, r: int, ks) -> list[dict]:
    """Exact 2x scan via the involution (k-1-x)_k = (-1)^k (x)_k.

    On the domain x = k+b, b in [0,g), the involution sends b -> g-1-b and
    multiplies the value by (-1)^k, so the upper half of the column carries
    no new information:

        k even   values repeat    -> scan b in [0, ceil(g/2)) against s
        k odd    values negate    -> scan the same range against s and -s

    The -s branch costs a subtraction, not a second multiply, because
    (-s)F[b] = p - (s F[b]) for r != 0. So both parities do half the modular
    multiplications of the full scan, on half the temporary memory.

    Output is byte-identical to scan_ks_full, b included. The least full-scan
    hit is min(first s-hit, g-1 - last (-s)-hit): -s hits at half-index b
    correspond to full-scan hits at g-1-b, which decreases as b grows, so the
    LAST -s hit gives the smallest full index.
    """
    out: list[dict] = []
    p64 = np.int64(p)
    rp = _check_r(r, p)
    for k in ks:
        k = int(k)
        g = p - k
        if g <= 0:
            continue
        half = (g + 1) // 2
        s = rp * int(F[k]) % p
        left = F[k : k + half]
        right = (np.int64(s) * F[:half]) % p64
        eq = left == right
        b = int(eq.argmax()) if bool(eq.any()) else None
        if k % 2 == 1:
            # Negate in place: `right` is dead once eq is built, and reusing
            # the buffer is worth ~1.2x on odd columns at Band II sizes.
            np.subtract(p64, right, out=right)
            eq2 = left == right
            if bool(eq2.any()):
                last = eq2.size - 1 - int(eq2[::-1].argmax())
                cand = g - 1 - last
                b = cand if b is None else min(b, cand)
            del eq2
        if b is not None:
            # g_even is kept for the recorded jsonl format and the
            # pre-registration tables; k_odd is the same bit, but it is the
            # INVARIANT one -- a property of the polynomial f_k, fixed once,
            # where g changes with every prime. (p odd => g even <=> k odd.)
            out.append({"k": k, "g": g, "g_even": g % 2 == 0,
                        "k_odd": k % 2 == 1, "b": b})
        del eq
    return out


def fact_at(p: int, k: int) -> int:
    """k! mod p in min(k, p-k) multiplications.

    Wilson's theorem splits (p-1)! at k:  k! (p-1-k)! = (-1)^(k+1) (mod p),
    so with g = p-k,  k! = (-1)^(k+1) / (g-1)!.  When k is close to p -- which
    is the whole Band II and Z-jump regime -- that turns O(k) into O(g).
    """
    g = p - k
    if k <= g:
        acc = 1
        for i in range(1, k + 1):
            acc = acc * i % p
        return acc
    acc = 1
    for i in range(1, g):
        acc = acc * i % p                      # (g-1)!
    v = pow(acc, -1, p)
    return v if (k + 1) % 2 == 0 else (-v) % p


def fact_window(p: int, lo: int, n: int):
    """[lo!, (lo+1)!, ..., (lo+n-1)!] mod p, without ever building all of p."""
    import numpy as np

    out = np.empty(n, dtype=np.int64)
    acc = fact_at(p, lo)
    for j in range(n):
        out[j] = acc
        acc = acc * (lo + j + 1) % p
    return out


def scan_ks_windowed(p: int, r: int, ks) -> list[dict]:
    """scan_ks without the p-sized factorial table.

    The scan reads only F[0:half] and F[k:k+half] with half = ceil(g/2), which
    is O(g) entries -- but fact_table built all p of them. Getting F at the
    high offset needs F[lo] for lo = min(ks), and Wilson delivers that in
    p-lo-1 multiplications instead of lo.

    Cost drops from p to about 1.5*g_max + (spread of ks). Measured on real
    i=8 workloads: 1.75x on Band II, 1476x on the Z-jump, where the first live
    prime usually sits just above k so g/p is near zero and the old code built
    a multi-million entry table to run a fifty-step test.

    Output is identical to scan_ks, record for record.
    """
    import numpy as np

    ks = [int(v) for v in ks]
    if not ks:
        return []
    rp = _check_r(r, p)
    lo = min(ks)
    halves = {k: (p - k + 1) // 2 for k in ks}
    hi = max(k + halves[k] for k in ks)
    H = max(halves.values())
    F_low = fact_window(p, 0, H)
    F_hi = fact_window(p, lo, hi - lo + 1)
    p64 = np.int64(p)
    out: list[dict] = []
    for k in ks:
        g = p - k
        if g <= 0:
            continue
        half = halves[k]
        s = rp * int(F_hi[k - lo]) % p
        left = F_hi[k - lo : k - lo + half]
        right = (np.int64(s) * F_low[:half]) % p64
        eq = left == right
        b = int(eq.argmax()) if bool(eq.any()) else None
        if k % 2 == 1:
            np.subtract(p64, right, out=right)
            eq2 = left == right
            if bool(eq2.any()):
                last = eq2.size - 1 - int(eq2[::-1].argmax())
                cand = g - 1 - last
                b = cand if b is None else min(b, cand)
            del eq2
        if b is not None:
            out.append({"k": k, "g": g, "g_even": g % 2 == 0,
                        "k_odd": k % 2 == 1, "b": b})
        del eq
    return out


def scan_ks(F: np.ndarray, p: int, r: int, ks) -> list[dict]:
    """Survivors in ks. Half-scan by default; identical output either way."""
    if USE_HALF_SCAN:
        return scan_ks_half(F, p, r, ks)
    return scan_ks_full(F, p, r, ks)


def scan_columns(
    p: int,
    ks,
    r_expected: int | None = None,
    *,
    N: int = N,
    K: int = K,
) -> tuple[int, list[dict]]:
    """Check Band II r, scan. Returns (r, survivors). Worker entry."""
    if USE_WINDOWED_SCAN:
        # Two independent table-free routes to r, so dropping the factorial
        # table does not drop the cross-check r_checked used to provide.
        r = r_closed(p, N=N, K=K)
        alt = r_two_digit_delta(p, N=N, K=K)
        if r != alt:
            raise RuntimeError(f"r(p) closed {r} != delta {alt} at p={p}")
        if r_expected is not None and r != r_expected:
            raise RuntimeError(f"r(p) {r} != expected {r_expected} at p={p}")
        return r, scan_ks_windowed(p, r, ks)
    F = fact_table(p)
    r = r_checked(F, p, falling=False, N=N, K=K)
    if r_expected is not None and r != r_expected:
        raise RuntimeError(f"r(p) {r} != expected {r_expected} at p={p}")
    return r, scan_ks(F, p, r, ks)


def r_two_digit(F: np.ndarray, p: int, *, N: int = N, K: int = K) -> int:
    """C(N,K) mod p via two-digit Lucas from the factorial table.

    Works for every p > k_max of Band I leftover and for p > N/2.
    C(α,β) C(n0,k0); α = N//p < p on the two-digit range.
    """
    a, b = N // p, K // p
    n0, k0 = N - a * p, K - b * p
    if not (0 <= b <= a and 0 <= k0 <= n0):
        return 0
    cab = int(F[a]) * pow(int(F[b]), -1, p) * pow(int(F[a - b]), -1, p) % p
    c0 = int(F[n0]) * pow(int(F[k0]), -1, p) * pow(int(F[n0 - k0]), -1, p) % p
    return cab * c0 % p


def scan_columns_general(p: int, ks, *, N: int = N, K: int = K) -> tuple[int, list[dict]]:
    """Image scan. Two-digit r(p) only when p^2 > N; else full Lucas.

    cells() already skips p <= sqrt(N), so this branch is a guard. Do not
    use C(alpha,beta) C(n0,k0) below sqrt(N) -- that formula is false: N
    then has three or more base-p digits and the product drops one.
    """
    if p * p <= N:
        from singmaster_intersect import binom_mod_lucas

        r = int(binom_mod_lucas(N, K, p))
    elif USE_WINDOWED_SCAN:
        r = r_two_digit_delta(p, N=N, K=K)          # table-free
    else:
        r = r_two_digit(fact_table(p), p, N=N, K=K)
    if r == 0:
        raise RuntimeError(f"live prime {p} has r=0 (Z / digit-0); do not scan")
    if USE_WINDOWED_SCAN:
        return r, scan_ks_windowed(p, r, ks)
    return r, scan_ks(fact_table(p), p, r, ks)


def cumulative_g(k_end: int, kmin: int, p: int) -> int:
    if k_end < kmin:
        return 0
    m = k_end - kmin + 1
    return m * (2 * p - kmin - k_end) // 2


def equal_g_chunks(kmin: int, kmax: int, p: int, n_chunks: int) -> list[tuple[int, int]]:
    """Contiguous [lo, hi] with equal Σ(p-k). Covers kmin..kmax exactly."""
    total = cumulative_g(kmax, kmin, p)
    edges = [kmin]
    for i in range(1, n_chunks):
        target = total * i // n_chunks
        lo, hi = kmin, kmax
        while lo < hi:
            mid = (lo + hi) // 2
            if cumulative_g(mid, kmin, p) < target:
                lo = mid + 1
            else:
                hi = mid
        if lo <= edges[-1]:
            lo = edges[-1] + 1
        if lo > kmax - (n_chunks - i) + 1:
            lo = kmax - (n_chunks - i) + 1
        edges.append(lo)
    edges.append(kmax + 1)
    chunks = []
    for i in range(n_chunks):
        lo = edges[i]
        hi = edges[i + 1] - 1
        if hi < lo:
            raise RuntimeError(f"empty chunk {i}: {lo}>{hi}")
        chunks.append((lo, hi))
    if chunks[0][0] != kmin or chunks[-1][1] != kmax:
        raise RuntimeError("chunk coverage failed")
    for i in range(1, n_chunks):
        if chunks[i][0] != chunks[i - 1][1] + 1:
            raise RuntimeError("chunk gap")
    return chunks


def log10_central(k: int) -> float:
    return (math.lgamma(2 * k + 1) - 2 * math.lgamma(k + 1)) / math.log(10)


# ---------------------------------------------------------------------------
# shared Band I / cell-geometry helpers
#
# These used to be copy-pasted into fat_image_hunt / triple_hunt / walk_369 /
# stragglers_nearK / zjump / family_sweep. Identical maths in every copy, but
# they had already begun to drift (triple_hunt's cells() had dropped the "zlo"
# key). One home, so a fix lands once.
# ---------------------------------------------------------------------------

FAM8 = Fam(i=8, N=N, K=K)


def binom_mod_prime(n: int, k: int, p: int) -> int:
    """C(n,k) mod p for prime p and 0 <= k <= n < p. Modular product, not exact comb."""
    if k < 0 or k > n:
        return 0
    # THREE lower indices, not two. Beyond k and n-k there is the negation
    # identity: n = -(p-n) mod p, so
    #     C(n,k) = (-1)^k C(p-n+k-1, p-n-1),
    # whose lower index is p-n-1. That is tiny exactly when n is close to p,
    # which is the whole Band II and Z-jump regime -- measured 3890x there and
    # 940x on fat-cell primes, while never losing elsewhere (the min is taken).
    d1 = p - n - 1
    kk = min(k, n - k)
    if 0 <= d1 < kk:
        num = den = 1
        for i in range(d1):
            num = num * (k + d1 - i) % p
            den = den * (i + 1) % p
        c = num * pow(den, -1, p) % p if d1 else 1
        return (-c) % p if k % 2 else c
    k = kk
    if k == 0:
        return 1
    # One inverse at the end, not one per step: pow(i+1,-1,p) inside the loop
    # was ~5x the cost of the two multiplications it sat between.
    num = den = 1
    for i in range(k):
        num = num * (n - i) % p
        den = den * (i + 1) % p
    return num * pow(den, -1, p) % p


def r_two_digit_delta(p: int, *, N: int = N, K: int = K) -> int:
    """Table-free r(p) using the cheapest of THREE lower indices.

    Standard two-digit Lucas costs min(k0, n0-k0). But n0 = p - delta with
    delta = (alpha+1)p - N, and negating the upper index gives

        C(n0, k0) = C(p-delta, k0) = (-1)^k0 C(k0+delta-1, delta-1)  (mod p)

    whose lower index is delta-1. delta is small at the BOTTOM of a cell,
    which is exactly where the Z-jump lands and where every fat-image prime
    sits, so the true cost is min(k0, n0-k0, delta-1). r_closed is the
    alpha=1 case of this; this is the general form, valid in every cell.

    Requires p*p > N (two digits). Returns 0 on a Kummer carry.
    """
    a, b = N // p, K // p
    n0, k0 = N - a * p, K - b * p
    if not (0 <= b <= a and 0 <= k0 <= n0):
        return 0
    cab = binom_mod_prime(a, b, p)
    dlt = (a + 1) * p - N          # = p - n0
    if 0 <= dlt - 1 <= min(k0, n0 - k0):
        c0 = binom_mod_prime(k0 + dlt - 1, dlt - 1, p)
        if k0 % 2:
            c0 = (-c0) % p
    else:
        c0 = binom_mod_prime(n0, k0, p)
    return cab * c0 % p


def delta_of(p: int, *, N: int = N) -> int:
    """delta = (alpha+1)p - N, the third lower index. Small at a cell bottom."""
    return (N // p + 1) * p - N


def lucas_digits(p: int, *, N: int = N, K: int = K) -> tuple[int, int, int]:
    """Two-digit Lucas: (alpha, beta, C(N,K) mod p). Requires p*p > N.

    The product goes through r_two_digit_delta, which picks the cheapest of
    the three lower indices; the digits themselves are unchanged.
    """
    a, b = N // p, K // p
    return a, b, r_two_digit_delta(p, N=N, K=K)


def r_of(p: int, *, N: int = N, K: int = K) -> int:
    """C(N,K) mod p by two-digit Lucas, table-free. Requires p*p > N."""
    return lucas_digits(p, N=N, K=K)[2]


def inv_table(p: int, upto: int) -> list[int]:
    """inv[i] = i^-1 mod p for i = 1..upto, in O(upto) with no pow() calls."""
    inv = [0] * (upto + 1)
    if upto >= 1:
        inv[1] = 1
    for i in range(2, upto + 1):
        inv[i] = (p - (p // i) * inv[p % i] % p) % p
    return inv


def image_j(prod: int, k: int, p: int, inv: list[int]) -> int | None:
    """j with prod = (-1)^j C(g-1,j) (mod p), or None = kill.

    Uses C(k+j,k) = (-1)^j C(p-k-1,j) (mod p), so the whole image of column
    k mod p is walked in g = p-k steps with two modmuls each.
    """
    g = p - k
    if g <= 0:
        return None
    n = g - 1
    c = 1
    sign = 1
    for j in range(g):
        if (sign * c) % p == prod:
            return j
        if j + 1 >= g:
            break
        c = c * (n - j) % p
        c = c * inv[j + 1] % p
        sign = -sign
    return None


def cells(fam: Fam = FAM8) -> list[dict]:
    """(alpha,beta) prime windows with their Z (forced-zero) sub-range."""
    out = []
    p2 = fam.p_two
    amax = fam.N // p2 + 1
    for a in range(2, amax + 1):
        for b in range(1, a):
            plo = max(fam.N // (a + 1) + 1, fam.K // (b + 1) + 1)
            phi = min(fam.N // a, fam.K // b)
            if plo > phi or phi < p2 or plo > fam.K:
                continue
            s = a - b
            zlo = fam.D // s + 1
            if zlo <= plo:
                kind = "FULL"
            elif zlo <= phi:
                kind = "PART"
            else:
                kind = "NONE"
            z_first = max(plo, zlo) if kind != "NONE" else None
            z_last = phi if kind != "NONE" and z_first is not None else None
            if kind != "NONE" and z_first is not None and z_first > phi:
                kind = "NONE"
                z_first = z_last = None
            out.append(
                {
                    "a": a,
                    "b": b,
                    "plo": plo,
                    "phi": phi,
                    "kind": kind,
                    "zlo": zlo,
                    "z_first": z_first,
                    "z_last": z_last,
                }
            )
    out.sort(key=lambda w: w["plo"])
    return out


def live_intervals(windows: list[dict], fam: Fam = FAM8) -> list[tuple[int, int]]:
    """Prime ranges where r(p) can be non-zero: NONE cells, PART-lower, p>N/2."""
    ivs = []
    for w in windows:
        if w["kind"] == "NONE":
            ivs.append((w["plo"], w["phi"]))
        elif w["kind"] == "PART" and w["plo"] <= w["zlo"] - 1:
            ivs.append((w["plo"], w["zlo"] - 1))
    ivs.append((fam.N2 + 1, fam.D))
    ivs.sort()
    return ivs


_IV_CACHE: dict[int, list[int]] = {}


def _iv_starts(ivs) -> list[int]:
    """Cached list of interval left endpoints, keyed by identity of `ivs`."""
    key = id(ivs)
    got = _IV_CACHE.get(key)
    if got is None or len(got) != len(ivs):
        got = [lo for lo, _hi in ivs]
        _IV_CACHE[key] = got
    return got


def first_live_after(x: int, ivs: list[tuple[int, int]], d: int = D) -> int | None:
    """First live prime after x. Bisects to the right interval rather than
    rescanning from index 0.

    The old linear scan cost O(#intervals) per call -- 292 us at i=9's 7161
    intervals, which is 2.3 h of single-threaded parent time over 28.3M
    columns, all of it before any worker starts. Bisecting first is 142x.

    Safe because live_intervals returns DISJOINT intervals sorted by lo
    (verified: 0 overlaps and 0 containments at i=8 and i=9), so the last
    interval with lo <= p is the only one that can contain p.
    """
    p = int(gmpy2.next_prime(x))
    n = len(ivs)
    i = bisect.bisect_right(_iv_starts(ivs), p) - 1
    if i < 0:
        i = 0
    while p <= d and i < n:
        lo, hi = ivs[i]
        if p > hi:
            i += 1
            continue
        if p < lo:
            p = int(gmpy2.next_prime(lo - 1))
            i = max(i, bisect.bisect_right(_iv_starts(ivs), p) - 1)
            continue
        return int(p)
    return None


def live_primes(k: int, n: int, ivs: list[tuple[int, int]], d: int = D) -> list[int]:
    out = []
    x = k
    for _ in range(n):
        p = first_live_after(x, ivs, d)
        if p is None:
            break
        out.append(p)
        x = p
    return out


def first_primes(lo: int, hi: int, n: int) -> list[int]:
    ps = []
    p = int(gmpy2.next_prime(lo - 1))
    while p <= hi and len(ps) < n:
        ps.append(int(p))
        p = int(gmpy2.next_prime(p))
    return ps


def first_primes_above(n2: int, d: int, kmax: int, n: int = 16) -> list[int]:
    out = []
    p = int(gmpy2.next_prime(n2))
    while len(out) < n and p <= d:
        if p > kmax:
            out.append(int(p))
        p = int(gmpy2.next_prime(p))
    return out


def preceding_z(windows: list[dict], idx: int) -> dict | None:
    for j in range(idx - 1, -1, -1):
        w = windows[j]
        if w["kind"] in {"FULL", "PART"} and w["z_first"] is not None:
            return w
    return None


def chunk_ks(ks: list[int], p: int, n_chunks: int) -> list[list[int]]:
    """Split ks into n_chunks contiguous pieces of roughly equal sum(p-k)."""
    ks = sorted(ks)
    if not ks:
        return []
    if len(ks) < 2000 or n_chunks <= 1:
        return [ks]
    weights = [p - k for k in ks]
    total = sum(weights)
    if total <= 0:
        return [ks]
    chunks, acc, start, t = [], 0, 0, 1
    for i, w in enumerate(weights):
        acc += w
        if t < n_chunks and acc >= total * t / n_chunks:
            chunks.append(ks[start : i + 1])
            start = i + 1
            t += 1
    if start < len(ks):
        chunks.append(ks[start:])
    return [c for c in chunks if c]


CHECKPOINT_SCHEMA = 2
# 1: original chunk records {tag|prime_index|round, p, k_lo, k_hi, survivors}
# 2: survivor records additionally carry "k_odd" (2026-08-20)


def checkpoint_identity(**params) -> dict:
    """The header record that pins a checkpoint to one run's parameters."""
    return {"event": "schema", "version": CHECKPOINT_SCHEMA, **params}


def check_checkpoint(path, **params) -> None:
    """Refuse to resume a checkpoint written by different code or parameters.

    Resume merges old records with new ones. If the schema or the run's
    parameters changed in between, that merge is silent and the result is a
    certificate over a set of columns that were never all tested the same
    way. Cheaper to refuse and re-run than to trust it.
    """
    import json

    if not path.exists() or path.stat().st_size == 0:
        return
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("event") != "schema":
                break
            bad = {k: (v, params.get(k)) for k, v in rec.items()
                   if k not in ("event",) and k in params and params[k] != v}
            if rec.get("version") != CHECKPOINT_SCHEMA:
                bad["version"] = (rec.get("version"), CHECKPOINT_SCHEMA)
            if bad:
                raise SystemExit(
                    f"{path.name}: checkpoint does not match this run "
                    f"(recorded vs now: {bad}). Delete it and start over, "
                    f"or resume with the code/parameters that wrote it."
                )
            return
    raise SystemExit(
        f"{path.name}: checkpoint has no schema header, so it predates "
        f"schema {CHECKPOINT_SCHEMA} and cannot be safely resumed. "
        f"Delete it and start over."
    )


def iter_jsonl(path):
    """Yield each record in a checkpoint jsonl without materialising the file.

    read_jsonl builds a list of every record, and a chunk record carries its
    whole survivor list. At i=9 the checkpoint reaches 219 MB, which measures
    at about 5.7x that in Python objects -- roughly 1.25 GB held in the parent
    while eight workers are running, on a machine that has already been crashed
    once by memory pressure.

    Callers that only need a few fields per record (done_keys, phase_complete
    events, the survivors of ONE tag) should stream instead, so peak scales
    with what they keep rather than with the file.

    Deliberately NOT guarded against a malformed line: json.loads raises, so a
    record truncated by a kill or a power loss stops the resume instead of
    being silently skipped. read_jsonl has the same property and it is load
    bearing -- a skipped record would drop a chunk from the worklist and the
    run would certify over columns it never scanned.
    """
    if not path.exists():
        return
    import json

    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                yield json.loads(line)


def read_jsonl(path) -> list[dict]:
    """Every record in a checkpoint jsonl, or [] if it does not exist."""
    rows: list[dict] = []
    if not path.exists():
        return rows
    import json

    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def append_jsonl(path, rec: dict) -> None:
    import json

    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
        fh.flush()


def summarize_survivors(rows: list[dict], with_g: bool = False) -> dict:
    """Survivor census. with_g adds mean_g (zjump's shape)."""
    n = len(rows)
    if not n:
        out = {"n": 0, "even": None, "mean_k": None}
        if with_g:
            out["mean_g"] = None
        return out
    out = {
        "n": n,
        # identical numbers; "even" is the recorded name, "odd_k" the meaningful one
        "even": round(sum(1 for s in rows if s.get("g_even")) / n, 4),
        "odd_k": round(sum(1 for s in rows if s.get("k_odd", s.get("g_even"))) / n, 4),
        "mean_k": round(sum(s["k"] for s in rows) / n),
    }
    if with_g:
        out["mean_g"] = round(sum(s["g"] for s in rows) / n)
    return out
