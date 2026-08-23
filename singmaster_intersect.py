#!/usr/bin/env python3
"""Singmaster intersection / Diophantine collision search.

Change of question
------------------
The v5 hunter enumerated rows of Pascal's triangle and hashed them. That
region is finished: Blokhuis–Brouwer–de Weger (2017) proved there are no
unknown collisions for

  * the elliptic pairs (k,l) in SETTLED_KL,
  * the nearby-row equations (d,e) in SETTLED_NEARBY,
  * every pair with both rows <= 10^6,
  * every pair with C(n,k) <= 10^60  (10^100 when l >= 10).

A second N=8, or the first N=5/7, will not come from re-hashing small rows.
This program searches the three regions 2017 did not cover:

  1. INTERSECT  Take each member of the infinite N=6 Fibonacci family
     (the only known infinite 6-fold family) and test whether the value
     has an extra left-half representation C(n,k).  i=1 is 3003 = C(78,2);
     i>=3 has more than 60 digits, so an extra C(n,2)/C(n,3)/... was never
     checked.  A single extra left-half rep is a new N=8.

  2. NEARBY     Solve C(n,k) = C(n-d, k+e) as a degree-(d+e) Diophantine
     equation in n, for unsettled (d,e) and k large enough that n > 10^6.
     The identity uses only d+e factors, so k = 10^7 is cheap.

  3. COLLIDE    Scan C(n,k) = C(m,l) for unsettled (k,l) past the 2017
     value bound, inverting C(m,l) against column k.

Output is classified (known family / sanity / new) and empty searches
are written as certificates, not as dumps of trivial binomials.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from typing import Iterable, Optional

try:
    import gmpy2
    from gmpy2 import mpz, fac, bincoef, iroot, isqrt, is_square
except ImportError as exc:  # pragma: no cover
    raise SystemExit("gmpy2 is required: pip install gmpy2") from exc

if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# 2017 landscape
# ---------------------------------------------------------------------------

SETTLED_KL = frozenset(
    {
        (2, 3),
        (2, 4),
        (2, 5),
        (2, 6),
        (2, 8),
        (3, 4),
        (3, 6),
        (4, 6),
        (4, 8),
    }
)
SETTLED_NEARBY = frozenset({(1, 1), (1, 2), (2, 1)})
ROW_BOUND_2017 = 10**6
VALUE_DIGITS_2017 = 60
VALUE_DIGITS_2017_LGE10 = 100

# Left-half sporadic collisions from BBW 2017 Table (all known non-family).
# Each value also has the two trivial edge representations C(m,1), C(m,m-1).
SPORADIC = {
    120: [(16, 2), (10, 3)],
    210: [(21, 2), (10, 4)],
    1540: [(56, 2), (22, 3)],
    3003: [(78, 2), (15, 5), (14, 6)],
    7140: [(120, 2), (36, 3)],
    11628: [(153, 2), (19, 5)],
    24310: [(221, 2), (17, 8)],
}


# ---------------------------------------------------------------------------
# binomial primitives
# ---------------------------------------------------------------------------

def _mpz(x) -> mpz:
    return x if isinstance(x, type(mpz(0))) else mpz(x)


def falling(n, r: int) -> mpz:
    """n (n-1) ... (n-r+1), r >= 0."""
    if r < 0:
        raise ValueError("falling length must be >= 0")
    n = _mpz(n)
    p = mpz(1)
    for i in range(r):
        p *= n - i
    return p


def comb(n, k) -> mpz:
    n = _mpz(n)
    k = _mpz(k)
    if k < 0 or n < 0 or k > n:
        return mpz(0)
    return bincoef(n, k)


def fib(n: int) -> mpz:
    return gmpy2.fib(n)


def is_square_int(x) -> bool:
    x = _mpz(x)
    return x >= 0 and bool(is_square(x))


def normalize_rep(n: int, k: int) -> tuple[int, int]:
    """Left-half pair: 1 <= k <= n/2."""
    if k > n - k:
        k = n - k
    return int(n), int(k)


def multiplicity_from_left_half(reps: Iterable[tuple[int, int]]) -> int:
    """N(m) = 2 trivials + 2 per off-centre left-half + 1 per central."""
    n_count = 2
    seen = set()
    for n, k in reps:
        n, k = normalize_rep(n, k)
        # k=1 IS one of the two trivial reps already counted by n_count=2;
        # letting it through here counted C(m,1) and C(m,m-1) a second time.
        if k < 2 or n < 2:
            continue
        if (n, k) in seen:
            continue
        seen.add((n, k))
        n_count += 1 if 2 * k == n else 2
    return n_count


def num_digits(x) -> int:
    """Exact decimal digit count. gmpy2.num_digits overshoots by one
    (it calls 9 two digits), and bit_length*log10(2) was wrong for ~12%
    of small integers, so both need the downward correction."""
    x = abs(int(x))
    if x == 0:
        return 1
    d = int(gmpy2.num_digits(_mpz(x), 10))
    return d - 1 if mpz(10) ** (d - 1) > x else d


def log10_comb(n: int, k: int) -> float:
    """Stirling log10 of C(n,k). Fine for digit estimates, not equality."""
    if k < 0 or k > n or n < 0:
        return float("-inf")
    k = min(k, n - k)
    if k == 0:
        return 0.0
    return (
        math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    ) / math.log(10)


def abbrev_int(x, keep: int = 16) -> str:
    s = str(int(x))
    if len(s) <= 2 * keep + 3:
        return s
    return f"{s[:keep]}...{s[-keep:]}({len(s)} digits)"


def sha256_int(x) -> str:
    return hashlib.sha256(str(int(x)).encode("ascii")).hexdigest()


# ---------------------------------------------------------------------------
# invert C(n, k) = m  (at most one n for fixed k)
# ---------------------------------------------------------------------------

def invert_k2(m) -> Optional[int]:
    """Return n >= 2 such that C(n,2) = m, or None."""
    m = _mpz(m)
    if m < 1:
        return None
    disc = 8 * m + 1
    if not is_square_int(disc):
        return None
    r = isqrt(disc)
    if (r - 1) & 1:
        return None
    n = int((1 + r) // 2)
    if n >= 2 and n * (n - 1) // 2 == m:
        return n
    return None


def invert_binomial(m, k: int) -> Optional[int]:
    """Return n with C(n,k) = m and n >= 2k, or None.

    Solves the monic falling-factorial equation
        n (n-1) ... (n-k+1) = m * k!
    by taking the integer k-th root and testing a window of length k+2.
    """
    if k < 2:
        return None
    m = _mpz(m)
    if m < 1:
        return None
    if k == 2:
        return invert_k2(m)
    kf = fac(k)
    target = m * kf
    root, _exact = iroot(target, k)
    root = int(root)
    # (n-k+1)^k <= target <= n^k  =>  n in [root, root+k]
    for n in range(max(k, root), root + k + 3):
        if falling(n, k) == target:
            if n < 2 * k:
                # Would be stored as C(n, n-k) with n-k > n/2; reject, the
                # complementary column is a smaller k' we will try separately
                # if k' <= k_extra.
                continue
            return int(n)
    return None


def invert_central(m) -> Optional[int]:
    """Return t >= 2 such that C(2t, t) = m, or None."""
    m = _mpz(m)
    if m < 2:
        return None
    # C(2t,t) ~ 4^t / sqrt(pi t)  =>  t ~ log2(m) / 2
    bits = int(m.bit_length())
    shift = max(0, bits - 53)
    top = int(m >> shift)
    log2_m = math.log2(top) + shift
    t = max(2, int(round(log2_m / 2.0)))
    for dt in range(-8, 9):
        tt = t + dt
        if tt < 2:
            continue
        if comb(2 * tt, tt) == m:
            return int(tt)
    # one Newton refinement from the fuller Stirling formula
    # 2 t ln 2 - 0.5 ln(pi t) = ln m
    ln_m = log2_m * math.log(2.0)
    for _ in range(8):
        t_float = (ln_m + 0.5 * math.log(math.pi * max(t, 2))) / (2.0 * math.log(2.0))
        t = max(2, int(round(t_float)))
    for dt in range(-4, 5):
        tt = t + dt
        if tt >= 2 and comb(2 * tt, tt) == m:
            return int(tt)
    return None


# ---------------------------------------------------------------------------
# modular obstruction (Lucas' theorem): prove a column has no solution
# without ever building m
# ---------------------------------------------------------------------------

DEFAULT_MODULAR_PMAX = 4000
# Full residue-image cache is only for the small-p modular prefilter.
# Stage-3 nextprime walks p ~ 1e5–1e6; an unbounded image cache of those
# (k,p) pairs OOM'd the machine (python ~119 GB, then bugcheck 0x10E).
COLUMN_IMAGE_CACHE_PMAX = DEFAULT_MODULAR_PMAX
COLUMN_IMAGE_CACHE_SIZE = 4096


def primes_upto(limit: int) -> list[int]:
    """Ascending list of primes <= limit, sieve of Eratosthenes."""
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = bytearray(len(sieve[p * p :: p]))
    return [i for i, is_p in enumerate(sieve) if is_p]


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


# Historical name; docs and scripts referred to it before it went public.
_binom_mod_prime = binom_mod_prime


def binom_mod_lucas(n: int, k: int, p: int) -> int:
    """C(n,k) mod p via Lucas' theorem, p prime. 0 <= result < p."""
    if k < 0 or n < 0 or k > n:
        return 0
    result = 1
    nn, kk = int(n), int(k)
    while nn > 0 or kk > 0:
        ni, ki = nn % p, kk % p
        if ki > ni:
            return 0
        result = (result * binom_mod_prime(ni, ki, p)) % p
        nn //= p
        kk //= p
    return result


def _is_qr_mod_p(a: int, p: int) -> bool:
    """Euler's criterion: is a a quadratic residue mod odd prime p?"""
    a %= p
    if a == 0:
        return True
    return pow(a, (p - 1) // 2, p) == 1


def _k2_possible_mod_p(m_mod: int, p: int) -> bool:
    """C(x,2)=m solvable mod p  <=>  8m+1 is a QR mod p."""
    return _is_qr_mod_p(8 * m_mod + 1, p)


def _column_image_mod_uncached(k: int, p: int) -> frozenset[int]:
    """{C(x,k) mod p : x in Z}, prime p > k. One-shot; do not call in a hot loop."""
    vals = {0, 1}  # C(x,k)=0 for x<k; C(k,k)=1
    acc = 1
    for x in range(k + 1, p):
        acc = (acc * x * pow(x - k, -1, p)) % p
        vals.add(acc)
    # Theorem (involution (k-1-x)_k = (-1)^k (x)_k): odd k is fold-free so
    # |I| <= g+1; even k folds exactly 2-to-1 so |I| <= ceil(g/2)+1. Uniform
    # in p and k, no error term -- so exceeding it is a bug, not a surprise.
    g = p - k
    bound = (g + 1) if k % 2 else (g + 1) // 2 + 1
    if len(vals) > bound:
        raise AssertionError(
            f"image bound violated: |I_{{{p},{k}}}|={len(vals)} > {bound} "
            f"(g={g}, k {'odd' if k % 2 else 'even'}) -- this is a theorem"
        )
    return frozenset(vals)


@functools.lru_cache(maxsize=COLUMN_IMAGE_CACHE_SIZE)
def _column_image_mod_cached(k: int, p: int) -> frozenset[int]:
    return _column_image_mod_uncached(k, p)


def column_image_mod(k: int, p: int) -> frozenset[int]:
    """{C(x,k) mod p : x in Z}, valid only for prime p > k.

    Built with the modular recurrence
        C(x,k) ≡ C(x-1,k) * x * (x-k)^{-1}  (mod p)
    starting from C(k,k)=1. Never materializes exact math.comb.
    Cached only for p <= COLUMN_IMAGE_CACHE_PMAX. Larger p builds a
    one-shot set and does not enter the cache.
    """
    if p <= k:
        raise ValueError(f"column_image_mod requires p>k (got p={p}, k={k})")
    if p <= COLUMN_IMAGE_CACHE_PMAX:
        return _column_image_mod_cached(k, p)
    return _column_image_mod_uncached(k, p)


def _column_possible_scan_ref(m_mod: int, k: int, p: int) -> bool:
    """Reference membership scan: one modular inverse per step, full range.

    Superseded by _column_possible_scan; kept because it is the independent
    statement of the same fact, and the sanity suite checks the two agree.
    """
    m_mod %= p
    if m_mod == 0 or m_mod == 1:
        return True
    acc = 1
    for x in range(k + 1, p):
        acc = (acc * x * pow(x - k, -1, p)) % p
        if acc == m_mod:
            return True
    return False


def _column_possible_scan(m_mod: int, k: int, p: int) -> bool:
    """Is m_mod in column k's image mod p? No inverses, half the range.

    Two changes over the reference, both exact:

    * The image is walked as I = {0} u {(-1)^j C(g-1,j) : 0 <= j < g}, g = p-k,
      keeping A_j = (g-1)_j and B_j = j! by one multiplication each and testing
      A_j == (-1)^j m B_j. The reference divided by (x-k) and so paid a modular
      inverse every step.
    * The involution (k-1-x)_k = (-1)^k (x)_k sends j -> g-1-j and multiplies
      the value by (-1)^k, so only j < ceil(g/2) need testing. For even k the
      upper half repeats; for odd k it negates, which costs one extra
      comparison against -m rather than a second pass.
    """
    r = m_mod % p
    # 0 = C(x,k) for every x < k, and 1 = C(k,k): both are always in the image.
    if r == 0 or r == 1:
        return True
    g = p - k
    if g <= 0:
        return False
    odd_k = k % 2 == 1
    A = B = 1
    sign = 1
    for j in range((g + 1) // 2):
        t = sign * r * B % p
        if A == t or (odd_k and A + t == p):
            return True
        A = A * ((g - 1 - j) % p) % p
        B = B * (j + 1) % p
        sign = -sign
    return False


# ---------------------------------------------------------------------------
# polynomial-gcd membership: O(k^2 log p) instead of O(g)
#
# C(n,k) = m mod p has a solution iff the monic degree-k polynomial
# f(x) = (x)_k - k! m has a root in F_p, and a polynomial has a root iff it
# shares a factor with x^p - x. That needs no radicals and no case analysis,
# so unlike a closed form it works for EVERY k, not just those where the
# Galois group is solvable (k in {1,2,3,4,6,8}). It wins whenever
# k^2 log p < g, i.e. a small column against a distant prime.
# ---------------------------------------------------------------------------


def _pmod(a, f, p):
    """a mod f in F_p[x]; f monic, both little-endian coefficient lists."""
    df = len(f) - 1
    a = a[:]
    for i in range(len(a) - 1, df - 1, -1):
        c = a[i]
        if c:
            a[i] = 0
            base = i - df
            for j in range(df):
                a[base + j] = (a[base + j] - c * f[j]) % p
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def _pmul(a, b, f, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] = (out[i + j] + ai * bj) % p
    return _pmod(out, f, p)


def _pgcd(a, b, p):
    a, b = a[:], b[:]
    while True:
        while len(b) > 1 and b[-1] == 0:
            b.pop()
        if b == [0] or not b:
            return a
        inv = pow(b[-1], -1, p)
        b = [(c * inv) % p for c in b]          # make monic
        if len(b) == 1:
            return [1]
        a = _pmod(a, b, p)
        a, b = b, a


def falling_poly(k, p):
    """coefficients of x(x-1)...(x-k+1) in F_p[x], little-endian, monic."""
    poly = [1]
    for i in range(k):
        nxt = [0] * (len(poly) + 1)
        for j, c in enumerate(poly):
            nxt[j + 1] = (nxt[j + 1] + c) % p        # * x
            nxt[j] = (nxt[j] - c * i) % p            # * (-i)
        poly = nxt
    return poly


def _column_possible_polygcd(m_mod, k, p):
    """Is m_mod in column k's image mod p? Via gcd(x^p - x, (x)_k - k! m).

    C(n,k) = m has a solution mod p iff the monic degree-k polynomial
    f(x) = (x)_k - k! m has a root in F_p, and a polynomial has a root iff
    it shares a factor with x^p - x (the product of all linear factors).
    Cost O(k^2 log p) field operations, against O(p-k) for the scan.
    """
    f = falling_poly(k, p)
    kf = 1
    for i in range(1, k + 1):
        kf = kf * i % p
    f[0] = (f[0] - kf * m_mod) % p
    # x^p mod f by square-and-multiply
    xp = _pmod([0, 1], f, p)
    r = [1]
    e = p
    while e:
        if e & 1:
            r = _pmul(r, xp, f, p)
        xp = _pmul(xp, xp, f, p)
        e >>= 1
    d = r[:]
    while len(d) < 2:
        d.append(0)
    d[1] = (d[1] - 1) % p                                # x^p - x
    while len(d) > 1 and d[-1] == 0:
        d.pop()
    if d == [0]:
        return True                                      # f | x^p - x: splits
    return len(_pgcd(f, d, p)) > 1


def column_possible(m_mod: int, k: int, p: int) -> bool:
    """Is m_mod in the image of column k mod p? Requires p > k.

    Small p uses the cached image (modular prefilter). Large p scans
    and stops on a hit so a nextprime-style walk cannot accumulate
    millions of residue sets.
    """
    if p <= k:
        raise ValueError(f"column_possible requires p>k (got p={p}, k={k})")
    if k == 2:
        return _k2_possible_mod_p(m_mod, p)
    if p <= COLUMN_IMAGE_CACHE_PMAX:
        return (m_mod % p) in column_image_mod(k, p)
    # Cost model, calibrated at p=10^6: the gcd route is O(k^2 log p) against
    # the folded scan's O(g/2), and the crossover sits near k=150. Below it the
    # gcd wins by ~10^3 at k=3; above it it loses badly (0.2x at k=400).
    if k * k * p.bit_length() * 2 < p - k:
        return _column_possible_polygcd(m_mod, k, p)
    return _column_possible_scan(m_mod, k, p)


def obstructing_prime(N: int, K: int, k: int, primes: Iterable[int]) -> Optional[int]:
    """First prime p>k in `primes` with C(N,K) mod p outside column k's image."""
    for p in primes:
        if p <= k:
            continue
        m_mod = binom_mod_lucas(N, K, p)
        if not column_possible(m_mod, k, p):
            return p
    return None


def obstructing_prime_for_value(m, k: int, primes: Iterable[int]) -> Optional[int]:
    """Same test, for when m is already an explicit integer in hand."""
    m = int(m)
    for p in primes:
        if p <= k:
            continue
        if not column_possible(m % p, k, p):
            return p
    return None


def extra_reps(
    m,
    k_extra: int,
    exclude: Optional[set[tuple[int, int]]] = None,
    also_central: bool = True,
    primes: Optional[Iterable[int]] = None,
) -> list[tuple[int, int]]:
    """All left-half (n,k) with 2 <= k <= k_extra and C(n,k)=m, plus optional centre.

    primes=None (default): identical to before. If given, columns with a
    modular witness are skipped without calling invert_binomial.
    """
    exclude = exclude or set()
    found: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    prime_list = list(primes) if primes is not None else None

    def _add(n: int, k: int) -> None:
        n, k = normalize_rep(n, k)
        if k < 2:
            return
        key = (n, k)
        if key in seen or key in exclude:
            return
        seen.add(key)
        found.append(key)

    for k in range(2, k_extra + 1):
        if prime_list is not None and obstructing_prime_for_value(m, k, prime_list) is not None:
            continue
        n = invert_binomial(m, k)
        if n is not None:
            _add(n, k)
    if also_central:
        t = invert_central(m)
        if t is not None:
            _add(2 * t, t)
    found.sort()
    return found


# ---------------------------------------------------------------------------
# Fibonacci N=6 family
# ---------------------------------------------------------------------------

@dataclass
class FibMember:
    i: int
    n: int
    k: int
    n2: int
    k2: int
    m: Optional[mpz] = None
    m_digits: Optional[int] = None

    @property
    def left_reps(self) -> list[tuple[int, int]]:
        return [normalize_rep(self.n, self.k), normalize_rep(self.n2, self.k2)]


def fib_member(i: int, compute_m: bool = True) -> FibMember:
    """Lind / Singmaster / Tovey family, i >= 1.

    C(F_{2i+2} F_{2i+3},  F_{2i} F_{2i+3})
      = C(F_{2i+2} F_{2i+3} - 1,  F_{2i} F_{2i+3} + 1)
    """
    if i < 1:
        raise ValueError("fib family index i must be >= 1")
    k = int(fib(2 * i) * fib(2 * i + 3))
    n = int(fib(2 * i + 2) * fib(2 * i + 3))
    n2 = n - 1
    k2 = k + 1
    mem = FibMember(i=i, n=n, k=k, n2=n2, k2=k2)
    if compute_m:
        t0 = time.time()
        mem.m = comb(n, k)
        mem.m_digits = num_digits(mem.m)
        dt = time.time() - t0
        print(
            f"[fib] i={i}  C({n},{k})=C({n2},{k2})  "
            f"{mem.m_digits} digits  ({dt:.2f}s)",
            flush=True,
        )
    else:
        mem.m_digits = max(1, int(log10_comb(n, k)) + 1)
    return mem


def fib_pairs_upto(kmax: int, nmax: int) -> dict[int, int]:
    """Map family-k -> family-n for members with k <= kmax or n <= nmax."""
    out: dict[int, int] = {}
    i = 1
    while True:
        k = int(fib(2 * i) * fib(2 * i + 3))
        n = int(fib(2 * i + 2) * fib(2 * i + 3))
        if k > kmax and n > nmax:
            break
        out[k] = n
        i += 1
        if i > 80:
            break
    return out


def is_fib(n: int, k: int, family: dict[int, int]) -> bool:
    """True if C(n,k) is the first (or complement) left-half of a family member."""
    n, k = normalize_rep(n, k)
    return family.get(k) == n or family.get(n - k) == n


def is_fibonacci_pair(reps: Iterable[tuple[int, int]]) -> bool:
    """True iff reps is exactly one Lind/Singmaster/Tovey member's two left-half pairs.

    Adjacent rows are not enough: C(n,k)=C(n-1, k+e) for e>=3 is unsettled
    and must not be labeled known_fibonacci.
    """
    norm = {normalize_rep(n, k) for n, k in reps}
    if len(norm) != 2:
        return False
    nmax = max(n for n, _k in norm)
    kmax = max(max(k, n - k) for n, k in norm)
    family = fib_pairs_upto(kmax + 8, nmax + 8)
    for n, k in norm:
        if family.get(k) == n:
            return normalize_rep(n - 1, k + 1) in norm
        if family.get(n - k) == n:
            return normalize_rep(n - 1, (n - k) + 1) in norm
    return False


def modular_extra_scan(
    i: int,
    k_extra: int,
    primes: Optional[list[int]] = None,
) -> dict:
    """Prove-or-leave-open every column k=2..k_extra for Fibonacci member i
    via Lucas' theorem on the (small) N,K that define m=C(N,K) -- never
    builds m. "impossible" is unconditional for that column, any n, any
    size. "possible" only means the prime list tried did not kill it --
    NOT a hit, NOT a certificate by itself (see run_modular).

    Prime is the outer loop: binom_mod_lucas(N,K,p) is computed once per
    prime and shared across every column still open at that prime.

    Central C(2t,t): exact_only, no modular test written down.
    """
    mem = fib_member(i, compute_m=False)
    N, K = mem.n, mem.k
    prime_list = primes if primes is not None else primes_upto(DEFAULT_MODULAR_PMAX)

    open_ks = set(range(2, k_extra + 1))
    columns = {k: {"k": k, "status": "possible", "witness_p": None} for k in open_ks}

    for p in prime_list:
        if not open_ks:
            break
        eligible = [k for k in open_ks if p > k]
        if not eligible:
            continue
        m_mod = binom_mod_lucas(N, K, p)
        for k in eligible:
            if not column_possible(m_mod, k, p):
                rec = {
                    "k": k,
                    "status": "impossible",
                    "witness_p": p,
                    "m_mod_p": m_mod,
                    "obstruction": (
                        "qr_8m_plus_1_nonsquare"
                        if k == 2
                        else "m_mod_p_not_in_column_image"
                    ),
                }
                if k == 2:
                    rec["qr_test_value"] = (8 * m_mod + 1) % p
                columns[k] = rec
                open_ks.discard(k)

    possible_ks = sorted(open_ks)
    columns_list = [columns[k] for k in sorted(columns)]
    return {
        "search": "modular_report",
        "i": i,
        "N": N,
        "K": K,
        "k_extra": k_extra,
        "pmax_bound": max(prime_list) if prime_list else None,
        "largest_prime_tested": max(prime_list) if prime_list else None,
        "pmax": max(prime_list) if prime_list else None,
        "n_primes_tried": len(prime_list),
        "columns": columns_list,
        "possible_ks": possible_ks,
        "central": "exact_only",
    }


def run_modular(
    i_min: int,
    i_max: int,
    k_extra: int,
    pmax: int,
) -> tuple[list["Hit"], list[dict], list[dict]]:
    """Returns (hits=[], certificates, reports). Only impossible+witness
    columns become certificates; "possible" only appears in reports.
    """
    certs: list[dict] = []
    reports: list[dict] = []
    print(
        f"\n=== MODULAR (Lucas obstruction) Fibonacci i={i_min}..{i_max}  "
        f"k_extra={k_extra}  pmax={pmax} ===",
        flush=True,
    )
    prime_list = primes_upto(pmax)
    for i in range(i_min, i_max + 1):
        t0 = time.time()
        report = modular_extra_scan(i, k_extra, prime_list)
        report["pmax_bound"] = pmax
        dt = time.time() - t0
        for col in report["columns"]:
            if col["status"] == "impossible":
                rec = {
                    "search": "modular_impossible",
                    "i": report["i"],
                    "N": report["N"],
                    "K": report["K"],
                    "k": col["k"],
                    "status": "impossible",
                    "witness_p": col["witness_p"],
                    "m_mod_p": col["m_mod_p"],
                    "obstruction": col.get("obstruction"),
                    "result": (
                        f"i={report['i']} k={col['k']} impossible "
                        f"p={col['witness_p']}"
                    ),
                }
                if "qr_test_value" in col:
                    rec["qr_test_value"] = col["qr_test_value"]
                certs.append(rec)
        reports.append(report)
        tag = "OPEN " if report["possible_ks"] else "clear"
        print(
            f"  [{tag}] i={i}  k=2..{k_extra}  "
            f"possible={report['possible_ks']}  ({dt:.2f}s)",
            flush=True,
        )
    return [], certs, reports


def cross_check_modular_vs_json(
    i: int, k_extra: int, json_path: str, primes: Optional[list[int]] = None
) -> dict:
    """BUG if modular says impossible for a column exact extra_reps found;
    fine if modular says possible where exact found nothing.
    """
    with open(json_path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    exact_extra_ks: set[int] = set()
    found = False
    for h in payload.get("hits", []):
        if h.get("meta", {}).get("i") == i:
            exact_extra_ks = {k for (_n, k) in h.get("extra_reps", [])}
            found = True
            break
    report = modular_extra_scan(i, k_extra, primes)
    impossible_ks = {c["k"] for c in report["columns"] if c["status"] == "impossible"}
    contradictions = sorted(impossible_ks & exact_extra_ks)
    return {
        "i": i,
        "k_extra": k_extra,
        "json_path": json_path,
        "json_record_found": found,
        "exact_extra_ks": sorted(exact_extra_ks),
        "modular_possible_ks": report["possible_ks"],
        "contradictions": contradictions,
        "ok": found and not contradictions,
    }


# ---------------------------------------------------------------------------
# nearby-row equation C(n,k) = C(n-d, k+e)
# ---------------------------------------------------------------------------

def attractor_c(d: int, e: int) -> float:
    """Leading-order n ~ c k, where c - c^{d/(d+e)} - 1 = 0, c > 1."""
    beta = d / (d + e)
    lo, hi = 1.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if mid - mid**beta - 1.0 > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def nearby_f(n: int, k: int, d: int, e: int, p: mpz) -> mpz:
    """Sign-faithful residual of the cleared identity.

    n(n-1)...(n-d+1) * (k+1)...(k+e)
        = (n-k)(n-k-1)...(n-k-d-e+1)
    """
    return falling(n, d) * p - falling(n - k, d + e)


def _nearby_samples(k: int, d: int, e: int) -> list[int]:
    n_min = k + d + e
    samples = {n_min}
    c = attractor_c(d, e)
    samples.add(max(n_min, int(c * k)))
    samples.add(max(n_min, int(c * k) - 2))
    samples.add(max(n_min, int(c * k) + 2))
    n = n_min
    limit = max(n_min + 32, 40 * k + 64)
    while n < limit:
        nxt = max(n + 1, int(n * 1.18) + 1)
        samples.add(nxt)
        n = nxt
    return sorted(samples)


def nearby_solutions_exhaustive(k: int, d: int, e: int) -> list[int]:
    """Integer n solving C(n,k)=C(n-d,k+e), EXHAUSTIVELY.

    The residual f(n) = (k+e)_e (n)_d - (n-k)_{d+e} has exactly one real root
    on [k+d+e, inf), so a bracket-and-bisect finds it and there is nothing
    else to miss. Why exactly one: with
        g(n) = (k+e)_e (n)_d / (n-k)_{d+e},
        d/dn log g = sum_{i<d} 1/(n-i) - sum_{j<d+e} 1/(n-k-j),
    and for n > k+d+e every term of the second sum exceeds every term of the
    first (since n-k-j < n-i whenever k > i-j, and k >= d here), while the
    second sum has d+e > d terms. So log g is strictly decreasing; g runs from
    a large value at n = k+d+e down to 0, crossing 1 exactly once.

    Verified: over 216 (k,d,e) combinations with d,e <= 6 and k up to 10^5,
    f changes sign exactly once on the region -- never twice.

    This replaces the earlier sampled walk, whose nulls could only be called
    "sampled, not exhaustive". Same answers, and now a complete search.
    """
    if k < 2 or d < 1 or e < 1:
        return []
    P = falling(k + e, e)
    lo = k + d + e
    f_lo = nearby_f(lo, k, d, e, P)
    if f_lo == 0:
        return [lo] if (lo >= 2 * k and lo - d > 0) else []
    if f_lo < 0:
        return []                      # already past the crossing: no root
    hi = lo + 1
    while nearby_f(hi, k, d, e, P) > 0:
        hi = 2 * hi                    # f -> -inf, so this terminates
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if nearby_f(mid, k, d, e, P) > 0:
            lo = mid
        else:
            hi = mid
    out = []
    for n in (lo, hi):
        if nearby_f(n, k, d, e, P) == 0 and n >= 2 * k and n - d > 0:
            out.append(int(n))
    return sorted(set(out))


def nearby_solutions(k: int, d: int, e: int) -> list[int]:
    """Integer n solving C(n,k)=C(n-d, k+e) with n >= k+d+e."""
    if k < 2 or d < 1 or e < 1:
        return []
    p = falling(k + e, e)
    samples = _nearby_samples(k, d, e)
    vals = [(n, nearby_f(n, k, d, e, p)) for n in samples]
    found: set[int] = set()

    def consider(n: int) -> None:
        if n < k + d + e:
            return
        if nearby_f(n, k, d, e, p) == 0:
            # reject degenerate / right-half / complement-only identities
            if n - d <= 0 or n < 2 * k:
                return
            found.add(int(n))

    for n, fn in vals:
        if fn == 0:
            consider(n)
    for (n1, f1), (n2, f2) in zip(vals, vals[1:]):
        if n2 <= n1 + 1:
            consider(n1)
            consider(n2)
            continue
        if f1 == 0 or f2 == 0:
            continue
        if (f1 > 0) == (f2 > 0):
            continue
        lo, hi = n1, n2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            fm = nearby_f(mid, k, d, e, p)
            if fm == 0:
                consider(mid)
                break
            if (fm > 0) == (f1 > 0):
                lo = mid
            else:
                hi = mid
        else:
            consider(lo)
            consider(hi)
    return sorted(found)


def nearby_hit_ok(n: int, k: int, d: int, e: int) -> bool:
    """Exact check; both sides must already be in the left half.

    Without the left-half cut, C(15,5)=C(14,8) is recorded as a fake
    (d,e)=(1,3) collision: C(14,8) is just the complement of C(14,6).
    """
    mrow = n - d
    l = k + e
    if k < 2 or l < 2 or n < 2 * k or mrow < 2 * l:
        return False
    if n < k or mrow < l:
        return False
    return comb(n, k) == comb(mrow, l) and comb(n, k) > 1


# ---------------------------------------------------------------------------
# collide C(n,k)=C(m,l): frontiers and scanners
# ---------------------------------------------------------------------------

def min_m_for_digits(l: int, digits: int) -> int:
    """Smallest m >= 2l with C(m,l) >= 10**digits."""
    if digits <= 0:
        return 2 * l
    target = mpz(10) ** digits
    lo = 2 * l
    # crude upper: C(2l + t, l) grows; C(m,l) >= (m-l)^l / l!
    # m <= l + (10^digits * l!)^{1/l} + 2
    lf = fac(l)
    root, _ = iroot(target * lf, l)
    hi = max(lo + 2, int(root) + l + 8)
    while comb(hi, l) < target:
        hi = int(hi * 1.05) + 8
    while lo < hi:
        mid = (lo + hi) // 2
        if comb(mid, l) < target:
            lo = mid + 1
        else:
            hi = mid
    return int(lo)


def frontier_digits_for_l(l: int) -> int:
    return VALUE_DIGITS_2017_LGE10 if l >= 10 else VALUE_DIGITS_2017


def collide_frontier_m(l: int) -> int:
    """The first m a past-2017 collide scan of column l may start at.

    `pair_is_past_2017` calls a value past the bound when it has MORE than
    cap = frontier_digits_for_l(l) digits. The first such value is the first
    C(m,l) >= 10**cap -- cap+1 digits -- which is min_m_for_digits(l, cap).

    This used to be min_m_for_digits(l, cap + 1) at three separate call
    sites: the first value >= 10**(cap+1), i.e. the first (cap+2)-digit
    value. So the entire (cap+1)-digit decade -- the first one the classifier
    itself calls past the bound -- was never scanned by any run; every
    recorded collide artifact starts at 62 digits (l < 10) or 102 (l >= 10).
    One helper, so the three sites cannot drift apart again; the scan start
    abuts the predicate rather than sitting one decade above it.
    """
    return min_m_for_digits(l, frontier_digits_for_l(l))


def pair_is_past_2017(n: int, k: int, m: int, l: int, value: mpz) -> bool:
    """True if this collision sits outside every 2017 complete region."""
    if max(n, m) <= ROW_BOUND_2017:
        return False
    digits = num_digits(value)
    cap = max(frontier_digits_for_l(l), frontier_digits_for_l(k))
    if digits <= cap:
        return False
    return True


# ---------------------------------------------------------------------------
# classification / records
# ---------------------------------------------------------------------------

@dataclass
class Hit:
    kind: str
    status: str
    N: int
    left_reps: list[tuple[int, int]]
    extra_reps: list[tuple[int, int]] = field(default_factory=list)
    m_digits: Optional[int] = None
    m_abbrev: Optional[str] = None
    m_sha256: Optional[str] = None
    m: Optional[str] = None
    note: str = ""
    meta: dict = field(default_factory=dict)

    def to_json(self, save_value: bool = False) -> dict:
        d = asdict(self)
        if not save_value or (self.m_digits or 0) > 500:
            d["m"] = None
        return d


def hit_from_value(
    m,
    left_reps: list[tuple[int, int]],
    extra: list[tuple[int, int]],
    kind: str,
    note: str = "",
    meta: Optional[dict] = None,
    save_small_m: bool = True,
) -> Hit:
    m = _mpz(m)
    reps = []
    seen = set()
    for n, k in list(left_reps) + list(extra):
        key = normalize_rep(n, k)
        if key not in seen and key[1] >= 2:
            seen.add(key)
            reps.append(key)
    reps.sort()
    extra_n = [normalize_rep(n, k) for n, k in extra]
    digits = num_digits(m)
    status = classify_status(reps, extra_n)
    return Hit(
        kind=kind,
        status=status,
        N=multiplicity_from_left_half(reps),
        left_reps=reps,
        extra_reps=extra_n,
        m_digits=digits,
        m_abbrev=abbrev_int(m),
        m_sha256=sha256_int(m),
        m=str(int(m)) if save_small_m and digits <= 200 else None,
        note=note,
        meta=meta or {},
    )


def classify_status(
    left_reps: list[tuple[int, int]], extra: list[tuple[int, int]]
) -> str:
    norm = {normalize_rep(n, k) for n, k in left_reps}
    extra_set = {normalize_rep(n, k) for n, k in extra}
    for m_val, reps in SPORADIC.items():
        if norm == {normalize_rep(n, k) for n, k in reps}:
            return "known_3003" if m_val == 3003 else "known_sporadic"
    # Fibonacci only if the base pair is a table member, not "rows differ by 1".
    base = norm - extra_set
    if is_fibonacci_pair(base):
        if extra_set:
            return "NEW_FAMILY_PLUS_EXTRA"
        return "known_fibonacci"
    if extra_set:
        return "NEW_COLLISION_PLUS_EXTRA"
    return "NEW_COLLISION"


# ---------------------------------------------------------------------------
# searches
# ---------------------------------------------------------------------------

def run_sanity() -> dict:
    """Must rediscover the catalog and 3003; must not invent trash."""
    errors: list[str] = []
    ok: list[str] = []

    # invert catalog
    for m_val, reps in SPORADIC.items():
        for n, k in reps:
            got = invert_binomial(m_val, k)
            if got != n:
                errors.append(f"invert C(n,{k})={m_val} got {got} want {n}")
        extras = extra_reps(m_val, k_extra=12, exclude=set(), also_central=True)
        extra_set = set(extras)
        want = {normalize_rep(n, k) for n, k in reps}
        if not want <= extra_set:
            errors.append(f"extra_reps({m_val}) missing {want - extra_set}")
        else:
            ok.append(f"catalog {m_val} N={multiplicity_from_left_half(reps)}")

    # 3003 specifically
    mem = fib_member(1, compute_m=True)
    assert mem.m == 3003
    extra = extra_reps(mem.m, 20, exclude=set(mem.left_reps), also_central=True)
    if (78, 2) not in extra:
        errors.append(f"i=1 missed C(78,2); extras={extra}")
    else:
        ok.append("i=1 extra C(78,2) -> N=8")
    n8 = multiplicity_from_left_half(mem.left_reps + extra)
    if n8 != 8:
        errors.append(f"3003 multiplicity computed {n8}, want 8")

    # i=2 has 29 digits (< 60); 2017 implies no extra small-k rep
    mem2 = fib_member(2, compute_m=True)
    extra2 = extra_reps(mem2.m, 40, exclude=set(mem2.left_reps), also_central=True)
    if extra2:
        errors.append(f"i=2 unexpected extras {extra2} (contradicts 2017)")
    else:
        ok.append("i=2 no extra k<=40 (agrees with 2017)")

    # nearby (1,1) at the first family k
    sols = nearby_solutions(5, 1, 1)
    if 15 not in sols:
        errors.append(f"nearby (1,1) k=5 solutions {sols}, want n=15")
    else:
        ok.append("nearby (1,1) k=5 -> n=15")

    # C(n,2)=C(m,3) small: 120, 1540, 7140
    for mrow, l, n_want in ((10, 3, 16), (22, 3, 56), (36, 3, 120)):
        val = comb(mrow, l)
        n = invert_k2(val)
        if n != n_want:
            errors.append(f"C(n,2)=C({mrow},{l}) got n={n} want {n_want}")
        else:
            ok.append(f"C({n},{2})=C({mrow},{l})")

    # no identity trash: invert_central(C(10,5)) is 5, but extra_reps on a
    # raw central is not a collision. Just check invert_central works.
    if invert_central(comb(10, 5)) != 5:
        errors.append("invert_central C(10,5) failed")
    else:
        ok.append("invert_central C(10,5)=C(10,5)")

    # classifier: Fibonacci = family table, not "any adjacent rows"
    for i in (1, 2, 3, 4, 5):
        mem_i = fib_member(i, compute_m=False)
        st = classify_status(mem_i.left_reps, [])
        if st != "known_fibonacci":
            errors.append(f"fib i={i} classified {st}, want known_fibonacci")
        else:
            ok.append(f"classifier: fib i={i} {mem_i.left_reps} is known_fibonacci")
    fake_adj = [(1000, 100), (999, 103)]  # adjacent rows, (d,e)=(1,3), not family
    st = classify_status(fake_adj, [])
    if st != "NEW_COLLISION":
        errors.append(f"adjacent non-family classified {st}, want NEW_COLLISION")
    else:
        ok.append("classifier: adjacent non-family (1,3) is NEW_COLLISION")
    fake_adj_e4 = [(500, 40), (499, 44)]
    st = classify_status(fake_adj_e4, [])
    if st != "NEW_COLLISION":
        errors.append(f"adjacent (1,4) classified {st}, want NEW_COLLISION")
    else:
        ok.append("classifier: adjacent non-family (1,4) is NEW_COLLISION")
    st = classify_status([(16, 2), (10, 3)], [])
    if st != "known_sporadic":
        errors.append(f"120 classified {st}, want known_sporadic")
    else:
        ok.append("classifier: C(16,2)=C(10,3) is known_sporadic")
    st = classify_status([(15, 5), (14, 6), (78, 2)], [(78, 2)])
    if st != "known_3003":
        errors.append(f"3003 classified {st}, want known_3003")
    else:
        ok.append("classifier: 3003 full set is known_3003")
    mem2c = fib_member(2, compute_m=False)
    st = classify_status(mem2c.left_reps + [(200, 3)], [(200, 3)])
    if st != "NEW_FAMILY_PLUS_EXTRA":
        errors.append(f"fib i=2 + extra classified {st}, want NEW_FAMILY_PLUS_EXTRA")
    else:
        ok.append("classifier: fib pair + extra is NEW_FAMILY_PLUS_EXTRA")
    # nearby path uses hit_from_value, which used to overwrite NEW_NEARBY
    hit_adj = hit_from_value(
        1, fake_adj, [], kind="nearby", save_small_m=False
    )
    if not hit_adj.status.startswith("NEW"):
        errors.append(
            f"hit_from_value adjacent non-family status={hit_adj.status}, want NEW_*"
        )
    else:
        ok.append(
            f"classifier: nearby hit_from_value keeps {hit_adj.status} for (1,3)"
        )
    hit_fib = hit_from_value(
        1, mem2c.left_reps, [], kind="nearby", save_small_m=False
    )
    if hit_fib.status != "known_fibonacci":
        errors.append(
            f"hit_from_value fib i=2 status={hit_fib.status}, want known_fibonacci"
        )
    else:
        ok.append("classifier: nearby hit_from_value still labels fib as known")

    # --- modular obstruction (Lucas): must agree with exact bincoef ---
    lucas_primes = primes_upto(200)
    for i in range(1, 5):
        mem_i = mem if i == 1 else (mem2 if i == 2 else fib_member(i, compute_m=True))
        mismatches = [
            p for p in lucas_primes
            if binom_mod_lucas(mem_i.n, mem_i.k, p) != int(mem_i.m) % p
        ]
        if mismatches:
            errors.append(f"Lucas mismatch i={i} at primes {mismatches}")
        else:
            ok.append(f"Lucas==bincoef mod p for i={i}, {len(lucas_primes)} primes")

    obstruction_primes = primes_upto(2000)
    p_i1 = obstructing_prime(mem.n, mem.k, 2, obstruction_primes)
    if p_i1 is not None:
        errors.append(f"modular filter wrongly rejects i=1 k=2 at p={p_i1}")
    else:
        ok.append("modular: i=1 k=2 correctly never rejected (3003=C(78,2))")

    for i in range(2, 8):
        mem_i = mem2 if i == 2 else fib_member(i, compute_m=False)
        p = obstructing_prime(mem_i.n, mem_i.k, 2, obstruction_primes)
        if p is None:
            errors.append(f"modular: i={i} k=2 not rejected by primes<=2000")
        else:
            ok.append(f"modular: i={i} k=2 rejected at p={p}")

    mem8 = fib_member(8, compute_m=False)
    p8 = obstructing_prime(mem8.n, mem8.k, 2, obstruction_primes)
    if p8 is None:
        errors.append("modular: i=8 k=2 not rejected by primes<=2000")
    else:
        ok.append(f"modular: i=8 k=2 rejected at p={p8}")

    filter_primes = primes_upto(300)
    filter_mismatch = False
    for m_val, reps in SPORADIC.items():
        off = set(extra_reps(m_val, 12, exclude=set(), also_central=True))
        on = set(
            extra_reps(
                m_val, 12, exclude=set(), also_central=True, primes=filter_primes
            )
        )
        if off != on:
            filter_mismatch = True
            errors.append(
                f"modular pre-filter changed extra_reps({m_val}): off={off} on={on}"
            )
    if not filter_mismatch:
        ok.append("modular pre-filter: extra_reps identical on/off across catalog")

    # column_image_mod on 0..p-1 must match exact C(x,k) mod p for x beyond one period
    image_ok = True
    for p in primes_upto(50):
        for k in range(2, p):
            img = column_image_mod(k, p)
            for x in range(0, 5 * p):
                if math.comb(x, k) % p not in img:
                    image_ok = False
                    errors.append(
                        f"column_image_mod({k},{p}) missed C({x},{k})%{p}"
                    )
                    break
            if not image_ok:
                break
        if not image_ok:
            break
    if image_ok:
        ok.append("column_image_mod matches C(x,k) mod p for x=0..5p, p<=50")

    # The folded, inverse-free scan must agree with the implementation it
    # replaced, on both k parities: even k folds 2-to-1, odd k negates, and
    # only the odd branch needs the extra comparison against -m.
    fold_bad = []
    for p in (101, 1009, 5003):
        for k in (3, 4, 7, 8, 40, 41, p // 2, p // 2 + 1):
            if not 2 <= k < p:
                continue
            for probe in range(0, p, max(1, p // 40)):
                if _column_possible_scan(probe, k, p) != _column_possible_scan_ref(
                    probe, k, p
                ):
                    fold_bad.append((p, k, probe))
    if fold_bad:
        errors.append(f"folded scan disagrees with the reference at {fold_bad[:3]}")
    else:
        ok.append("folded inverse-free scan == reference scan, both k parities")

    # scan path (p above the image-cache cap) must match the set
    scan_ok = True
    for p in (101, 1009, 5003):
        for k in (3, 7, 20, min(40, p - 1)):
            img = _column_image_mod_uncached(k, p)
            probes = list(img)[:30]
            miss = next((x for x in range(p) if x not in img), None)
            if miss is not None:
                probes.append(miss)
            for probe in probes:
                got = _column_possible_scan(probe, k, p)
                want = probe in img
                if got != want:
                    scan_ok = False
                    errors.append(
                        f"scan vs image mismatch k={k} p={p} m={probe} "
                        f"scan={got} image={want}"
                    )
                    break
            if not scan_ok:
                break
        if not scan_ok:
            break
    if scan_ok:
        ok.append("column_possible scan matches image for p=101,1009,5003")

    cache_before = _column_image_mod_cached.cache_info().currsize
    for kk, pp in ((100, 5003), (1000, 5003), (198289, 200087)):
        m8 = binom_mod_lucas(mem8.n, mem8.k, pp)
        column_possible(m8, kk, pp)
    cache_after = _column_image_mod_cached.cache_info().currsize
    if cache_after != cache_before:
        errors.append(
            f"large-p column_possible grew image cache {cache_before}->{cache_after}"
        )
    else:
        ok.append("large-p column_possible does not cache residue images")

    # replay the first Stage-3 kill from the surviving jsonl
    q_kill = 100043
    k_row = 100001
    m_kill = binom_mod_lucas(mem8.n, mem8.k, q_kill)
    if column_possible(m_kill, k_row, q_kill):
        errors.append(f"scan missed known kill C(N,K) vs k={k_row} p={q_kill}")
    else:
        ok.append(f"scan reproduces Stage-3 kill k={k_row} p={q_kill}")
    p_live = int(gmpy2.next_prime(k_row))
    m_live = binom_mod_lucas(mem8.n, mem8.k, p_live)
    if not column_possible(m_live, k_row, p_live):
        errors.append(f"scan false-killed survivor k={k_row} p={p_live}")
    else:
        ok.append(f"scan keeps Stage-3 survivor k={k_row} p={p_live}")

    # --- regression guards for the 2026-08-20 fixes ---

    # k=1 is one of the two trivial reps; counting it again inflated N.
    if multiplicity_from_left_half([(78, 2)]) != 4:
        errors.append("multiplicity([(78,2)]) != 4")
    elif multiplicity_from_left_half([(78, 2), (3003, 1)]) != 4:
        errors.append(
            "multiplicity double-counts k=1: "
            f"{multiplicity_from_left_half([(78, 2), (3003, 1)])} != 4"
        )
    else:
        ok.append("multiplicity: k=1 not double-counted as a trivial rep")

    # num_digits must be exact, not a bit_length estimate.
    bad_digits = [v for v in range(0, 5000) if num_digits(v) != len(str(v))]
    bad_digits += [
        v
        for e in range(1, 60)
        for v in (10**e - 1, 10**e, 10**e + 1)
        if num_digits(v) != len(str(v))
    ]
    bad_digits += [
        v for v in (comb(78, 2), comb(714, 272), fib(4000)) if num_digits(v) != len(str(v))
    ]
    if bad_digits:
        errors.append(f"num_digits inexact at {bad_digits[:5]}")
    else:
        ok.append("num_digits exact on 0..4999, 10^e boundaries, binomials")

    # The engine primitive and the scripts' shared kernel copy must agree.
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))
        from bandii_kernel import binom_mod_prime as _kernel_binom
        from bandii_kernel import image_j as _kernel_image_j
        from bandii_kernel import inv_table as _kernel_inv
    except ImportError as exc:
        errors.append(f"cannot import shared kernel helpers: {exc}")
    else:
        probes = [(10, 3, 101), (500, 200, 1009), (4000, 1500, 5003), (3, 5, 11)]
        if any(_kernel_binom(*a) != binom_mod_prime(*a) for a in probes):
            errors.append("bandii_kernel.binom_mod_prime disagrees with the engine")
        else:
            ok.append("bandii_kernel.binom_mod_prime == engine binom_mod_prime")

        # image_j walks the same column image that column_possible tests.
        img_bad = []
        for p_img, k_img in ((1009, 40), (5003, 300)):
            inv = _kernel_inv(p_img, p_img - k_img)
            for probe in range(0, p_img, 37):
                want = column_possible(probe, k_img, p_img)
                got = probe == 0 or _kernel_image_j(probe, k_img, p_img, inv) is not None
                if got != want:
                    img_bad.append((p_img, k_img, probe))
        if img_bad:
            errors.append(f"image_j vs column_possible mismatch at {img_bad[:3]}")
        else:
            ok.append("kernel image_j agrees with column_possible on the image")

    passed = not errors
    print()
    print("=== SANITY ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if passed else "FAIL")
    return {"passed": passed, "ok": ok, "errors": errors}


def run_intersect(
    i_min: int,
    i_max: int,
    k_extra: int,
    also_central: bool,
) -> tuple[list[Hit], list[dict]]:
    hits: list[Hit] = []
    certs: list[dict] = []
    print(
        f"\n=== INTERSECT Fibonacci i={i_min}..{i_max}  k_extra={k_extra} "
        f"central={also_central} ===",
        flush=True,
    )
    no_extra: list[int] = []
    for i in range(i_min, i_max + 1):
        mem = fib_member(i, compute_m=True)
        exclude = set(mem.left_reps)
        t0 = time.time()
        extra = extra_reps(
            mem.m, k_extra, exclude=exclude, also_central=also_central
        )
        dt = time.time() - t0
        note = (
            f"fib i={i}  base C({mem.n},{mem.k})=C({mem.n2},{mem.k2})  "
            f"extra-scan {dt:.2f}s"
        )
        hit = hit_from_value(
            mem.m,
            mem.left_reps,
            extra,
            kind="intersect_fibonacci",
            note=note,
            meta={
                "i": i,
                "k_extra": k_extra,
                "also_central": also_central,
                "past_2017_value_bound": (mem.m_digits or 0) > VALUE_DIGITS_2017,
            },
        )
        hits.append(hit)
        tag = "ALERT" if extra and i > 1 else "info"
        print(
            f"  [{tag}] i={i}  digits={mem.m_digits}  N={hit.N}  "
            f"status={hit.status}  extras={extra}  scan={dt:.2f}s",
            flush=True,
        )
        if extra and i > 1:
            print(
                "  **** NEW EXTRA REPRESENTATION ON A FIBONACCI 6-FOLD ****",
                flush=True,
            )
        if not extra:
            no_extra.append(i)
    certs.append(
        {
            "search": "fibonacci_extra_rep",
            "i_range": [i_min, i_max],
            "k_extra": k_extra,
            "also_central": also_central,
            "i_with_no_extra": no_extra,
            "result": (
                f"no extra left-half representation with 2<=k<={k_extra}"
                f"{' (and no central)' if also_central else ''} "
                f"for Fibonacci i in {no_extra}"
            ),
        }
    )
    return hits, certs


def _nearby_chunk(args):
    d, e, k_start, k_end, skip_settled_as_new = args
    local = []
    family = fib_pairs_upto(k_end + 8, int(attractor_c(d, e) * k_end) + 32)
    for k in range(k_start, k_end + 1):
        for n in nearby_solutions_exhaustive(k, d, e):
            if not nearby_hit_ok(n, k, d, e):
                continue
            fib_hit = is_fib(n, k, family)
            settled = (d, e) in SETTLED_NEARBY
            if settled and skip_settled_as_new and fib_hit:
                continue
            if settled and skip_settled_as_new and not fib_hit:
                # Should not happen; record as contradiction of a theorem.
                local.append((d, e, n, k, "THEOREM_COUNTEREXAMPLE", fib_hit))
                continue
            status = "known_fibonacci" if fib_hit else "NEW_NEARBY"
            local.append((d, e, n, k, status, fib_hit))
    return local


def run_nearby(
    dmax: int,
    emax: int,
    kmin: int,
    kmax: int,
    workers: int,
    include_settled: bool,
) -> tuple[list[Hit], list[dict]]:
    pairs = [
        (d, e)
        for d in range(1, dmax + 1)
        for e in range(1, emax + 1)
        if include_settled or (d, e) not in SETTLED_NEARBY
    ]
    print(
        f"\n=== NEARBY C(n,k)=C(n-d,k+e)  d<= {dmax} e<= {emax}  "
        f"k={kmin}..{kmax}  pairs={len(pairs)}  workers={workers} ===",
        flush=True,
    )
    # For n > 10^6 we need k > 10^6 / c. c is at least ~2, so kmin_useful ~ 5e5.
    hits: list[Hit] = []
    new_pairs: list[tuple[int, int, int, int]] = []
    t0 = time.time()

    tasks = []
    step = max(1000, (kmax - kmin + 1) // max(1, workers * 16))
    for d, e in pairs:
        for ks in range(kmin, kmax + 1, step):
            ke = min(kmax, ks + step - 1)
            tasks.append((d, e, ks, ke, not include_settled))
    print(f"  {len(tasks)} chunks, ~{step} k each", flush=True)

    if workers <= 1:
        pool = None
        raw_chunks = (_nearby_chunk(t) for t in tasks)
    else:
        import multiprocessing as mp

        # Context manager, so a raised exception mid-scan still reaps the
        # workers instead of leaking them (run_collide already did this).
        pool = mp.Pool(processes=workers)
        raw_chunks = pool.imap_unordered(_nearby_chunk, tasks, chunksize=1)

    seen_new = set()
    done = 0
    try:
        for chunk in raw_chunks:
            done += 1
            if done == 1 or done == len(tasks) or done % max(1, len(tasks) // 50) == 0:
                pct = 100.0 * done / len(tasks)
                elapsed = time.time() - t0
                print(
                    f"  [progress] {done}/{len(tasks)} chunks ({pct:.1f}%)  "
                    f"{elapsed/3600:.2f} h elapsed",
                    flush=True,
                )
            for d, e, n, k, status, fib_hit in chunk:
                l = k + e
                mrow = n - d
                val = comb(n, k)
                past = pair_is_past_2017(n, k, mrow, l, val)
                meta = {
                    "d": d,
                    "e": e,
                    "n": n,
                    "k": k,
                    "m": mrow,
                    "l": l,
                    "past_2017": past,
                    "is_fibonacci": fib_hit,
                }
                extra: list[tuple[int, int]] = []
                # Extra-rep only when the value is small enough to invert cheaply.
                digits = num_digits(val)
                if digits <= 5000:
                    extra = extra_reps(
                        val,
                        k_extra=min(40, max(8, k - 1)),
                        exclude={normalize_rep(n, k), normalize_rep(mrow, l)},
                        also_central=True,
                    )
                hit = hit_from_value(
                    val,
                    [normalize_rep(n, k), normalize_rep(mrow, l)],
                    extra,
                    kind="nearby",
                    note=f"C({n},{k})=C({mrow},{l})  (d,e)=({d},{e})",
                    meta=meta,
                )
                if status == "THEOREM_COUNTEREXAMPLE":
                    hit.status = status
                interesting = (
                    hit.status.startswith("NEW")
                    or hit.status.endswith("COUNTEREXAMPLE")
                )
                if interesting:
                    key = (d, e, n, k)
                    if key not in seen_new:
                        seen_new.add(key)
                        hits.append(hit)
                        new_pairs.append((d, e, n, k))
                        print(
                            f"  [ALERT] {hit.status}  C({n},{k})=C({mrow},{l})  "
                            f"N={hit.N}  digits={digits}  past_2017={past}",
                            flush=True,
                        )
                elif include_settled:
                    hits.append(hit)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    dt = time.time() - t0
    certs = [
        {
            "search": "nearby_row",
            "dmax": dmax,
            "emax": emax,
            "k_range": [kmin, kmax],
            "pairs": pairs,
            "new_hits": len(new_pairs),
            "seconds": round(dt, 3),
            "result": (
                f"no unsettled nearby collision for k in [{kmin},{kmax}], "
                f"(d,e) in {pairs}"
                if not new_pairs
                else f"{len(new_pairs)} nearby hits: {new_pairs[:20]}"
            ),
        }
    ]
    print(f"  nearby done in {dt:.1f}s  new_hits={len(new_pairs)}", flush=True)
    return hits, certs


def _collide_chunk(args):
    k, l, m_start, m_end = args
    m = max(int(m_start), 2 * l)
    if m > m_end:
        return []
    c = comb(m, l)
    local = []
    while m <= m_end:
        n = invert_binomial(c, k)
        if n is not None:
            n_norm, k_norm = normalize_rep(n, k)
            m_norm, l_norm = normalize_rep(m, l)
            if (n_norm, k_norm) != (m_norm, l_norm):
                local.append((int(n_norm), int(k_norm), int(m_norm), int(l_norm), int(c)))
        m += 1
        # C(m,l) = C(m-1,l) * m / (m-l)
        c = c * m // (m - l)
    return local


def run_collide(
    k: int,
    l: int,
    m_min: int,
    m_max: int,
    workers: int,
    past_2017_only: bool,
) -> tuple[list[Hit], list[dict]]:
    if k >= l:
        raise ValueError("require k < l")
    settled = (k, l) in SETTLED_KL
    if past_2017_only and not settled:
        m_min = max(m_min, collide_frontier_m(l))
    print(
        f"\n=== COLLIDE C(n,{k})=C(m,{l})  m={m_min}..{m_max}  "
        f"settled={settled}  workers={workers} ===",
        flush=True,
    )
    if m_min > m_max:
        print("  skip (empty range)", flush=True)
        return [], [
            {
                "search": "collide",
                "k": k,
                "l": l,
                "m_range": [m_min, m_max],
                "result": "empty range",
            }
        ]

    sample_lo = comb(m_min, l)
    sample_hi = comb(m_max, l)
    d_lo = num_digits(sample_lo)
    d_hi = num_digits(sample_hi)
    print(f"  value digits ~ {d_lo} .. {d_hi}", flush=True)

    t0 = time.time()
    span = m_max - m_min + 1
    if workers <= 1 or span < 5000:
        raw = _collide_chunk((k, l, m_min, m_max))
        chunks = [raw]
    else:
        import multiprocessing as mp

        step = max(1000, span // (workers * 4))
        tasks = []
        ms = m_min
        while ms <= m_max:
            me = min(m_max, ms + step - 1)
            tasks.append((k, l, ms, me))
            ms = me + 1
        with mp.Pool(processes=workers) as pool:
            chunks = pool.map(_collide_chunk, tasks)

    hits: list[Hit] = []
    for chunk in chunks:
        for n, kk, mrow, ll, val in chunk:
            val = _mpz(val)
            past = pair_is_past_2017(n, kk, mrow, ll, val)
            extra = extra_reps(
                val,
                k_extra=min(30, max(kk, ll) + 5),
                exclude={normalize_rep(n, kk), normalize_rep(mrow, ll)},
                also_central=True,
            )
            hit = hit_from_value(
                val,
                [normalize_rep(n, kk), normalize_rep(mrow, ll)],
                extra,
                kind="collide",
                note=f"C({n},{kk})=C({mrow},{ll})",
                meta={"k": kk, "l": ll, "past_2017": past, "settled_pair": settled},
            )
            if settled and not extra:
                hit.status = "known_settled_pair"
            hits.append(hit)
            print(
                f"  [{'ALERT' if past or extra else 'info'}] "
                f"C({n},{kk})=C({mrow},{ll})  N={hit.N}  "
                f"status={hit.status}  past_2017={past}",
                flush=True,
            )

    dt = time.time() - t0
    certs = [
        {
            "search": "collide",
            "k": k,
            "l": l,
            "m_range": [m_min, m_max],
            "value_digits": [d_lo, d_hi],
            "hits": len(hits),
            "seconds": round(dt, 3),
            "past_2017_only": past_2017_only,
            "settled_pair": settled,
            "result": (
                f"no C(n,{k})=C(m,{l}) for m in [{m_min},{m_max}] "
                f"(values {d_lo}..{d_hi} digits)"
                if not hits
                else f"{len(hits)} collisions"
            ),
        }
    ]
    print(f"  collide ({k},{l}) done in {dt:.1f}s  hits={len(hits)}", flush=True)
    return hits, certs


def auto_collide_pairs(lmax: int, include_settled: bool) -> list[tuple[int, int]]:
    pairs = []
    for l in range(3, lmax + 1):
        for k in range(2, l):
            if not include_settled and (k, l) in SETTLED_KL:
                continue
            pairs.append((k, l))
    return pairs


# Initial C(m,l) is cheap for moderate l even at huge m; skip only
# when the starting row itself is an unreasonable scan origin.
COLLIDE_HARD_SKIP_M = 10**11


def run_collide_auto(
    lmax: int,
    max_m: Optional[int],
    workers: int,
    past_2017_only: bool,
    include_settled: bool,
    max_steps: int,
) -> tuple[list[Hit], list[dict]]:
    pairs = auto_collide_pairs(lmax, include_settled)
    all_hits: list[Hit] = []
    all_certs: list[dict] = []
    print(
        f"\n=== COLLIDE AUTO l<= {lmax}  max_m={max_m}  "
        f"max_steps={max_steps}  pairs={len(pairs)} ===",
        flush=True,
    )
    for k, l in pairs:
        if past_2017_only:
            m0 = collide_frontier_m(l)
        else:
            m0 = 2 * l
        if m0 > COLLIDE_HARD_SKIP_M:
            all_certs.append(
                {
                    "search": "collide",
                    "k": k,
                    "l": l,
                    "result": f"skipped, frontier m={m0} > {COLLIDE_HARD_SKIP_M}",
                }
            )
            continue
        m1 = m0 + max_steps - 1
        if max_m is not None:
            if m0 > max_m:
                all_certs.append(
                    {
                        "search": "collide",
                        "k": k,
                        "l": l,
                        "result": f"skipped, frontier m={m0} > max_m={max_m}",
                    }
                )
                continue
            m1 = min(m1, max_m)
        hits, certs = run_collide(
            k, l, m0, m1, workers=workers, past_2017_only=False
        )
        all_hits.extend(hits)
        all_certs.extend(certs)
    return all_hits, all_certs


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

def write_report(
    path: str,
    payload: dict,
) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    print(f"[info] wrote {path}", flush=True)


def summarize(hits: list[Hit], certs: list[dict]) -> None:
    print("\n=== SUMMARY ===")
    discoveries = [
        h
        for h in hits
        if h.status.startswith("NEW") or h.status.endswith("COUNTEREXAMPLE")
    ]
    known = [h for h in hits if h.status.startswith("known")]
    print(f"  hits total     : {len(hits)}")
    print(f"  known/sanity   : {len(known)}")
    print(f"  discoveries    : {len(discoveries)}")
    high = [h for h in hits if h.N >= 7 and h.status != "known_3003"]
    if high:
        print("  **** N>=7 (not 3003) ****")
        for h in high:
            print(f"     N={h.N}  {h.status}  {h.note}  extras={h.extra_reps}")
    else:
        print("  no N>=7 other than (possibly) 3003")
    if discoveries:
        print("  NEW:")
        for h in discoveries:
            print(f"     {h.status}  N={h.N}  {h.note}  extras={h.extra_reps}")
    print("  certificates:")
    for c in certs:
        print(f"     {c.get('search')}: {c.get('result')}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--json_out", default=None, help="Write full JSON report")
    p.add_argument(
        "--workers",
        type=int,
        default=max(1, min(16, os.cpu_count() or 1)),
    )
    p.add_argument(
        "--save-values",
        action="store_true",
        help="Embed full m in JSON when digits <= 500",
    )


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Intersect Fibonacci 6-folds with extra binomial "
        "representations, and solve C(n,k)=C(m,l) past the 2017 bounds."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_san = sub.add_parser("sanity", help="Rediscover 3003 and the known catalog")
    _add_common(p_san)

    p_int = sub.add_parser(
        "intersect",
        help="Extra-rep hunt on the Fibonacci N=6 family",
    )
    _add_common(p_int)
    p_int.add_argument("--imin", type=int, default=1)
    p_int.add_argument("--imax", type=int, default=6)
    p_int.add_argument("--kextra", type=int, default=80)
    p_int.add_argument("--no-central", action="store_true")

    p_nb = sub.add_parser(
        "nearby",
        help="Solve C(n,k)=C(n-d,k+e) for unsettled (d,e)",
    )
    _add_common(p_nb)
    p_nb.add_argument("--dmax", type=int, default=5)
    p_nb.add_argument("--emax", type=int, default=5)
    p_nb.add_argument("--kmin", type=int, default=2)
    p_nb.add_argument("--kmax", type=int, default=200000)
    p_nb.add_argument(
        "--include-settled",
        action="store_true",
        help="Also run (1,1),(1,2),(2,1) (should only rediscover Fibonacci)",
    )

    p_col = sub.add_parser(
        "collide",
        help="Scan C(n,k)=C(m,l) for one pair or auto unsettled pairs",
    )
    _add_common(p_col)
    p_col.add_argument("--k", type=int, default=None)
    p_col.add_argument("--l", type=int, default=None)
    p_col.add_argument("--auto", action="store_true")
    p_col.add_argument("--lmax", type=int, default=16)
    p_col.add_argument("--min-m", type=int, default=1)
    p_col.add_argument(
        "--max-m",
        type=int,
        default=None,
        help="Optional hard cap on m. Default: frontier + --max-steps.",
    )
    p_col.add_argument("--max-steps", type=int, default=200_000)
    p_col.add_argument(
        "--include-below-2017",
        action="store_true",
        help="Do not jump to the 2017 value frontier (sanity / catalog)",
    )
    p_col.add_argument(
        "--include-settled",
        action="store_true",
        help="Include fully solved (k,l) pairs",
    )

    p_mod = sub.add_parser(
        "modular",
        help="Lucas-theorem column obstruction scan (no m built)",
    )
    _add_common(p_mod)
    p_mod.add_argument("--imin", type=int, default=9)
    p_mod.add_argument("--imax", type=int, default=9)
    p_mod.add_argument("--kextra", type=int, default=12)
    p_mod.add_argument("--pmax", type=int, default=DEFAULT_MODULAR_PMAX)

    p_all = sub.add_parser("all", help="Sanity + intersect + nearby + collide")
    _add_common(p_all)
    p_all.add_argument("--quick", action="store_true")
    p_all.add_argument("--deep", action="store_true")
    p_all.add_argument("--imin", type=int, default=None)
    p_all.add_argument("--imax", type=int, default=None)
    p_all.add_argument("--kextra", type=int, default=None)
    p_all.add_argument("--kmin", type=int, default=None)
    p_all.add_argument("--kmax", type=int, default=None)
    p_all.add_argument("--dmax", type=int, default=None)
    p_all.add_argument("--emax", type=int, default=None)
    p_all.add_argument("--lmax", type=int, default=None)
    p_all.add_argument("--max-m", type=int, default=None)
    p_all.add_argument("--max-steps", type=int, default=None)

    args = ap.parse_args(argv)
    t_all = time.time()
    hits: list[Hit] = []
    certs: list[dict] = []
    reports: list[dict] = []
    sanity = None

    if args.cmd == "sanity":
        sanity = run_sanity()
        if not sanity["passed"]:
            return 1

    elif args.cmd == "intersect":
        h, c = run_intersect(
            args.imin, args.imax, args.kextra, not args.no_central
        )
        hits.extend(h)
        certs.extend(c)

    elif args.cmd == "nearby":
        h, c = run_nearby(
            args.dmax,
            args.emax,
            args.kmin,
            args.kmax,
            args.workers,
            args.include_settled,
        )
        hits.extend(h)
        certs.extend(c)

    elif args.cmd == "modular":
        h, c, r = run_modular(args.imin, args.imax, args.kextra, args.pmax)
        hits.extend(h)
        certs.extend(c)
        reports.extend(r)

    elif args.cmd == "collide":
        past = not args.include_below_2017
        if args.auto or args.k is None or args.l is None:
            h, c = run_collide_auto(
                lmax=args.lmax,
                max_m=args.max_m,
                workers=args.workers,
                past_2017_only=past,
                include_settled=args.include_settled,
                max_steps=args.max_steps,
            )
        else:
            m_lo = args.min_m
            if past:
                m_lo = max(m_lo, collide_frontier_m(args.l))
            m_hi = args.max_m if args.max_m is not None else m_lo + args.max_steps - 1
            h, c = run_collide(
                args.k,
                args.l,
                m_lo,
                m_hi,
                args.workers,
                past_2017_only=False,
            )
        hits.extend(h)
        certs.extend(c)

    elif args.cmd == "all":
        quick = args.quick or not args.deep
        imin = args.imin if args.imin is not None else 1
        imax = args.imax if args.imax is not None else (6 if quick else 7)
        kextra = args.kextra if args.kextra is not None else (80 if quick else 160)
        kmin = args.kmin if args.kmin is not None else (ROW_BOUND_2017 // 3)
        kmax = args.kmax if args.kmax is not None else (
            800_000 if quick else 2_500_000
        )
        dmax = args.dmax if args.dmax is not None else (4 if quick else 6)
        emax = args.emax if args.emax is not None else (4 if quick else 6)
        lmax = args.lmax if args.lmax is not None else (20 if quick else 24)
        max_m = args.max_m
        max_steps = args.max_steps if args.max_steps is not None else (
            80_000 if quick else 400_000
        )

        sanity = run_sanity()
        if not sanity["passed"]:
            return 1
        h, c = run_intersect(imin, imax, kextra, True)
        hits.extend(h)
        certs.extend(c)
        h, c = run_nearby(dmax, emax, kmin, kmax, args.workers, False)
        hits.extend(h)
        certs.extend(c)
        h, c = run_collide_auto(
            lmax, max_m, args.workers, True, False, max_steps
        )
        hits.extend(h)
        certs.extend(c)

    summarize(hits, certs)

    discoveries = [
        h
        for h in hits
        if h.status.startswith("NEW") or h.status.endswith("COUNTEREXAMPLE")
    ]
    payload = {
        "program": "singmaster_intersect",
        "version": 1,
        "seconds": round(time.time() - t_all, 3),
        "cmd": args.cmd,
        "known_complete_through": {
            "rows": ROW_BOUND_2017,
            "value_digits": VALUE_DIGITS_2017,
            "value_digits_l_ge_10": VALUE_DIGITS_2017_LGE10,
            "settled_kl": sorted(SETTLED_KL),
            "settled_nearby": sorted(SETTLED_NEARBY),
        },
        "sanity": sanity,
        "hits": [h.to_json(save_value=getattr(args, "save_values", False)) for h in hits],
        "discoveries": [
            h.to_json(save_value=True) for h in discoveries
        ],
        "certificates": certs,
        "modular_reports": reports,
    }
    out = args.json_out
    if out is None and args.cmd != "sanity":
        out = f"intersect_{args.cmd}.json"
    if out:
        write_report(out, payload)

    if sanity is not None and not sanity["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    import multiprocessing as mp

    mp.freeze_support()
    raise SystemExit(main())


