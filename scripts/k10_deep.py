#!/usr/bin/env python3
"""Q28's wall, audited with Q27's battery-depth lesson: it is REAL, not shallow.

Q27's genus-3 curve fell because two descent classes that looked alive mod 2^7
were dead at 2-adic depth 13 and 6 (scripts/k8_case2.py).  This lab applies
the same lesson to Q28's remaining wall -- the no-kill locus of the (2,3)
branch (docs/q28-k10-intersective.md section 3b): over Z,

    R(t^2) - 1024c = (t^2+mu t+nu)(t^2-mu t+nu)(t^3+al t^2+be t+ga)(t^3-al t^2+be t-ga)

with (mu^2-4nu)*disc(m) a square, whose unproved missing statement is "no
integer point with c != 0".  Verdict, in one line: THE WALL STANDS -- the
locus is p-adically alive at EVERY depth at p = 2, 3, 5, 7 (spot-certified
at 11, 13, 31 as well), with certified analytic branch witnesses, so no
congruence supported on these primes can close it, and the section-3
elliptic machinery (Magma's IntegralQuarticPoints) remains the missing
tool.  PROVED_23_POSITIVE stays False.

WHY THE LESSON DOES NOT TRANSFER (the structural finding).  Q27's classes
were 0-dimensional families n^2 = delta*P6(-delta*m^2) over a FREE integer
parameter m, and the value's 2-adic valuation was m-uniform (exactly 13),
so one deep modulus killed a whole class.  Q28's no-kill locus is a
1-dimensional p-adic CURVE passing through 160 genuine integer points --
the c = 0 factorisations of R(t^2) itself, 160 tuples (a,nu,ga,al,be)
over ten a-values (the ten c = 0 points of section 3b) -- and c is a
nonconstant analytic function vanishing there.  Along each branch,
offsetting a = a0 + p^t gives certified Z_p-solutions of the full system
with

    v_p(cnum) = t + const_p        (cnum = p0^2 - 945^2 = 1024c),

e.g. v2 = t + 14, v3 = v5 = v7 = t + 1 (and v11 = v13 = v31 = t) on the
branch through (a, nu, ga, al, be) = (-130, -63, -15, -9, 23).  The
valuation SWEEPS every value as t grows: at each certified prime p and
every modulus p^k, the locus carries Z_p-points with v_p(cnum) finite,
10!-compatible, and taking every value >= const_p -- so no congruence at
these primes can force c = 0.  A square polynomial value on a
free-parameter family can hide a fixed deep obstruction; a curve through
honest points cannot.  (The certified statement is for the listed primes;
nothing suggests other primes behave differently -- every prime beyond the
10!-primes where the construction was pointed gave the same finite law
v_p = t -- but it was run only there.)

The lab also records a methodological trap that nearly produced a fake kill
here: a lift-tree with an aggressive branch-neighborhood prune reported the
2-adic fiber at a = -10 + 2^12 "dead at depth 2^19" -- but the exact branch
computation shows that very offset carries a certified point with
v2(cnum) = 26 (the prune window was tighter than the branch's derivative
allowed).  Deep batteries can produce false kills as well as miss real
ones; only exact-valuation laws or unpruned certificates settle a class.

Structure of this script:
  [0] setup: the coincidence simplifies -- mu^2 = 2nu - a exactly, so
      (mu^2-4nu) = -(a+2nu): the coincidence is -(a+2nu)*disc(m) = square,
      and the whole system needs only (a, nu, ga, al, be).  A Q27-style
      delta-split does NOT exist: on the b = nu^2 cover, b*(s-b) = q_C(a)
      makes b the SQUARE ITSELF, not a cofactor of a y^2-product, so no
      valuation-evenness is forced on anything (the finite kernel
      Res(s, q_C) = 3^10 * 11^2 * 748871 is computed for the record).
      The 160 degenerate c = 0 tuples are generated and checked exactly.
  [1] calibration: the battery harness reproduces Q27's delta = 3 class --
      ALIVE mod 2^7, DEAD mod 2^14 -- proving a 2^7-capped battery was the
      bug there, and that this harness would catch such a kill here.
  [2] the filtered lift-tree battery (10!-filters v2 >= 18, v3 >= 4,
      v5 >= 2, v7 >= 1), a linear Hensel lift with the same semantics as
      digit brute force (cross-checked in the test): alive at every prime
      at every depth reached.  DEEP_RECORD holds the deep table, raw
      5-tuple counts, reproducible with --deep.  The equation set is
      invariant under the mirror (ga, al) -> (-ga, -al) (it swaps the two
      cubic factors), so survivor sets are mirror-closed; counts here are
      raw tuples, NOT mirror orbits.
  [3] THE MAIN RESULT: exact p-adic branch witnesses at p = 2, 3, 5, 7
      plus spot-checks at 11, 13, 31 (integer Hensel square roots + a
      Newton solve for alpha, precision up to 2^56), with the linear
      valuation laws above and the two square-class side conditions
      certified decisively (nonzero, even valuation, resolved unit that is
      a p-adic square).  SURVIVE at every prime.
  [4] the prune-artifact demonstration (pinned).
  [5] the recorded candidates c = 1395418752000 (kill 11) and 2235340800
      (kill 13) are OFF the locus (D_q*D_k nonsquare) and still die.

    python scripts/k10_deep.py           (fast steps, seconds)
    python scripts/k10_deep.py --deep    (re-derives DEEP_RECORD, ~75 s)
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]

REQ_V = {2: 18, 3: 4, 5: 2, 7: 1}       # v_p(cnum) needed for 10! | c
TREE_DEPTH = {2: 8, 3: 4, 5: 3, 7: 2}   # fast default depths
DEEP_DEPTHS = {2: 10, 3: 9, 5: 6, 7: 5}  # depths of the recorded deep runs
# raw 5-tuple survivor counts (and how many have v_p(cnum) resolved finite
# below the depth), reproducible with --deep; sets are mirror-closed, these
# are NOT orbit counts:
DEEP_RECORD = {2: (10, 6881280, 0), 3: (9, 6181920, 5724000),
               5: (6, 3322500, 3309220), 7: (5, 603680, 597520)}
WITNESS_LAW = {2: 14, 3: 1, 5: 1, 7: 1,  # v_p(cnum) = t + const on the branch
               11: 0, 13: 0, 31: 0}      # spot-checks outside the 10!-window
WITNESS_BASE = (-130, -63, -15, -9, 23)
WALL_STANDS = True                       # PROVED_23_POSITIVE stays False
P6_Q27 = [1, -126, 5271, -82564, 570591, -5779998, -9458775]


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise RuntimeError(f"FAILED: {msg}")


def v_p(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


# ---------------------------- the locus equations ----------------------------

def eqs_ok(a, nu, ga, al, be, M) -> bool:
    b = nu * nu % M
    if (-a**4 - 165*a**3 + 3*a*a*b - 8778*a*a + 330*a*b - 172810*a
            - b*b + 8778*b - 1057221) % M:
        return False
    if (ga*ga - a**3 - 165*a*a + 2*a*b - 8778*a + 165*b - 172810) % M:
        return False
    if (2*be + 165 + a - al*al) % M:
        return False
    if (be*be - 2*al*ga - (a*a + 165*a - b + 8778)) % M:
        return False
    return True


def sqclass_ok(x, p, k) -> bool:
    """could x (known mod p^k) be a p-adic square?  (necessary condition)"""
    x %= p**k
    if x == 0:
        return True
    v = v_p(x, p)
    if v % 2:
        return False
    u = x // p**v
    head = k - v
    if p == 2:
        if head >= 3:
            return u % 8 == 1
        return u % 4 == 1 if head == 2 else True
    return pow(u % p, (p-1)//2, p) == 1 if head >= 1 else True


def sqclass_decisive(x, p, k) -> bool:
    """is x (known mod p^k) CERTAINLY a p-adic square?  True only when the
    valuation is resolved (nonzero mod p^k), even, and the unit head is
    resolved far enough to certify a Z_p-square by Hensel: u = 1 mod 8 with
    >= 3 unit bits at p = 2, or a quadratic residue unit at odd p."""
    x %= p**k
    if x == 0:
        return False
    v = v_p(x, p)
    if v % 2:
        return False
    u = x // p**v
    head = k - v
    if p == 2:
        return head >= 3 and u % 8 == 1
    return head >= 1 and pow(u % p, (p-1)//2, p) == 1


def degenerate_tuples() -> list:
    out = set()
    base = [1, 3, 5, 7, 9]
    for signs in product([1, -1], repeat=5):
        roots = [s*r for s, r in zip(signs, base)]
        for pi in combinations(range(5), 2):
            pr = [roots[i] for i in pi]
            tr = [roots[i] for i in range(5) if i not in pi]
            nu = pr[0]*pr[1]
            mu = -(pr[0] + pr[1])
            al = -(tr[0] + tr[1] + tr[2])
            be = tr[0]*tr[1] + tr[0]*tr[2] + tr[1]*tr[2]
            ga = -tr[0]*tr[1]*tr[2]
            out.add((2*nu - mu*mu, nu, ga, al, be))
    return sorted(out)


# --------------------------------- steps -------------------------------------

def step0() -> None:
    print("  [0] setup: the coincidence needs no mu, and no Q27-split exists")
    a, nu, mu = sp.symbols("a nu mu")
    check(sp.expand((2*nu - a) - 4*nu - (-(a + 2*nu))) == 0,
          "mu^2 - 4nu = -(a + 2nu) given mu^2 = 2nu - a")
    s = 3*a**2 + 330*a + 8778
    qC = a**4 + 165*a**3 + 8778*a**2 + 172810*a + 1057221
    r = int(sp.resultant(sp.Poly(s, a), sp.Poly(qC, a)))
    check(abs(r) == 3**10 * 11**2 * 748871, "Res(s, q_C) = 3^10 * 11^2 * 748871")
    check(sp.Poly(qC, a).is_irreducible, "q_C irreducible: no rational split of the cover")
    deg = degenerate_tuples()
    check(len(deg) == 160, "160 degenerate c = 0 tuples")
    check(len({t[0] for t in deg}) == 10, "over ten a-values (the ten c = 0 points)")
    check(all(all(v == 0 for v in _eq_vals(*t)) for t in deg),
          "every degenerate tuple satisfies the system exactly over Z")
    check(WITNESS_BASE in deg, "the witness base is one of them")
    check((-10, -3, 315, -7, -53) in deg, "so is the prune-trap base")
    print("      mu^2 = 2nu - a exactly, so the coincidence is -(a+2nu)*disc(m) = square;")
    print("      on the cover, b*(s-b) = q_C(a) makes b the square ITSELF, not a")
    print("      cofactor: no evenness is forced anywhere -- no delta-classes arise.")
    print("      160 degenerate c = 0 tuples (ten a-values) satisfy the system exactly.")


def q27_delta3_value(m: int) -> int:
    x = -3 * m * m
    v = 0
    for cc in P6_Q27:
        v = v * x + cc
    return 3 * v


def step1() -> None:
    print("  [1] calibration: Q27's delta = 3 class through this harness")
    for k, expect_alive in ((7, True), (14, False)):
        N = 2**k
        sq = {(n_*n_) % N for n_ in range(N)}
        alive = any((q27_delta3_value(mm)) % N in sq for mm in range(N))
        check(alive == expect_alive, f"delta=3 alive mod 2^{k} should be {expect_alive}")
    vals = {v_p(q27_delta3_value(mm), 2) for mm in range(1, 41, 2)}
    check(vals == {13}, "the kill is the exact odd valuation 13 for odd m")
    print("      ALIVE mod 2^7, DEAD mod 2^14 (v2 = 13, odd): a 2^7-capped battery")
    print("      was the bug at k=8.  The same harness now probes k=10's locus.")


def _eq_vals(a, nu, ga, al, be):
    b = nu * nu
    return (
        -a**4 - 165*a**3 + 3*a*a*b - 8778*a*a + 330*a*b - 172810*a
        - b*b + 8778*b - 1057221,
        ga*ga - a**3 - 165*a*a + 2*a*b - 8778*a + 165*b - 172810,
        2*be + 165 + a - al*al,
        be*be - 2*al*ga - (a*a + 165*a - b + 8778),
    )


_GRADS = None


def _grads():
    """Jacobian entry functions of the four equations, derived by sympy
    (not by hand) so the linear lift below cannot drift from eqs_ok."""
    global _GRADS
    if _GRADS is None:
        xs = sp.symbols("a nu ga al be")
        a, nu, ga, al, be = xs
        b = nu * nu
        eqs = [
            -a**4 - 165*a**3 + 3*a*a*b - 8778*a*a + 330*a*b - 172810*a
            - b*b + 8778*b - 1057221,
            ga*ga - a**3 - 165*a*a + 2*a*b - 8778*a + 165*b - 172810,
            2*be + 165 + a - al*al,
            be*be - 2*al*ga - (a*a + 165*a - b + 8778),
        ]
        _GRADS = [[sp.lambdify(xs, sp.diff(e, x), "math") for x in xs]
                  for e in eqs]
    return _GRADS


def _solve_template(J, p):
    """RREF of the 4x5 matrix J over F_p with the row-op matrix T recorded:
    returns (E, T, pivots, free_cols, zero_rows) so that J d = r has the
    solutions read off from T r for every right-hand side r."""
    E = [row[:] for row in J]
    T = [[int(i == j) for j in range(4)] for i in range(4)]
    pivots = []
    r_ = 0
    for c in range(5):
        pr = next((i for i in range(r_, 4) if E[i][c] % p), None)
        if pr is None:
            continue
        E[r_], E[pr] = E[pr], E[r_]
        T[r_], T[pr] = T[pr], T[r_]
        inv = pow(E[r_][c], -1, p)
        E[r_] = [v * inv % p for v in E[r_]]
        T[r_] = [v * inv % p for v in T[r_]]
        for i in range(4):
            if i != r_ and E[i][c]:
                f = E[i][c]
                E[i] = [(E[i][j] - f*E[r_][j]) % p for j in range(5)]
                T[i] = [(T[i][j] - f*T[r_][j]) % p for j in range(4)]
        pivots.append((r_, c))
        r_ += 1
        if r_ == 4:
            break
    piv_cols = {c for _, c in pivots}
    free_cols = [c for c in range(5) if c not in piv_cols]
    return E, T, pivots, free_cols, list(range(r_, 4))


def _child_filters_ok(y, p, k) -> bool:
    if ((y[1]*y[2])**2 - 893025) % p**min(k, REQ_V[p]):
        return False
    if not sqclass_ok(2*y[1] - y[0], p, k):
        return False
    dm = (18*y[3]*y[4]*y[2] - 4*y[3]**3*y[2] + y[3]*y[3]*y[4]*y[4]
          - 4*y[4]**3 - 27*y[2]*y[2])
    return sqclass_ok((-(y[0] + 2*y[1])) * dm, p, k)


def _tree_fast(p, kmax, collect=False):
    """DFS over the filtered lift-tree via exact linear Hensel lifting:
    for k >= 2, E(x + d p^(k-1)) = E(x) + p^(k-1) J(x) d (mod p^k), so the
    lifts of a node are the F_p-solutions of J d = -E(x)/p^(k-1) -- the
    same node set as digit brute force (_tree_slow pins this), without the
    p^5 trial factor.  Returns (leaves, frozen, sorted_leaf_list|None)
    where frozen counts leaves with v_p(cnum) resolved below kmax."""
    grads = _grads()
    tmpl = {}
    roots = [t for t in product(range(p), repeat=5) if eqs_ok(*t, p)]
    leaves = frozen = 0
    out = [] if collect else None
    Mk = p**kmax
    stack = [(t, 1) for t in roots]
    while stack:
        x, k = stack.pop()
        if k == kmax:
            leaves += 1
            if ((x[1]*x[2])**2 - 893025) % Mk:
                frozen += 1
            if collect:
                out.append(x)
            continue
        cls = tuple(v % p for v in x)
        t_ = tmpl.get(cls)
        if t_ is None:
            J = [[int(g(*cls)) % p for g in row] for row in grads]
            t_ = tmpl[cls] = _solve_template(J, p)
        E, T, pivots, free_cols, zero_rows = t_
        Mp = p**k
        ev = _eq_vals(*x)
        r = [(-(e // Mp)) % p for e in ev]
        rt = [sum(T[i][j]*r[j] for j in range(4)) % p for i in range(4)]
        if any(rt[i] for i in zero_rows):
            continue
        kk = k + 1
        for assign in product(range(p), repeat=len(free_cols)):
            d = [0]*5
            for c, val in zip(free_cols, assign):
                d[c] = val
            for ri, ci in reversed(pivots):
                s_ = rt[ri]
                for j in free_cols:
                    if j > ci and E[ri][j]:
                        s_ -= E[ri][j]*d[j]
                d[ci] = s_ % p
            y = tuple(x[i] + d[i]*Mp for i in range(5))
            if _child_filters_ok(y, p, kk):
                stack.append((y, kk))
    return leaves, frozen, (sorted(out) if collect else None)


def _tree_slow(p, kmax):
    """digit brute force, the reference semantics for _tree_fast (test-pinned)."""
    nodes = {t for t in product(range(p), repeat=5) if eqs_ok(*t, p)}
    for k in range(2, kmax + 1):
        M, Mp = p**k, p**(k-1)
        new = set()
        for (a, nu, ga, al, be) in nodes:
            for da, dn, dg, dl, db in product(range(p), repeat=5):
                t = (a + da*Mp, nu + dn*Mp, ga + dg*Mp, al + dl*Mp, be + db*Mp)
                if eqs_ok(*t, M) and _child_filters_ok(t, p, k):
                    new.add(t)
        nodes = new
    return sorted(nodes)


def step2(depths: dict) -> dict:
    print("  [2] the 10!-filtered lift-tree battery on the k=10 locus")
    table = {}
    for p in (2, 3, 5, 7):
        kmax = depths[p]
        leaves, frozen, _ = _tree_fast(p, kmax)
        check(leaves > 0, f"p={p}: tree unexpectedly dead at depth {kmax}")
        table[p] = (kmax, leaves, frozen)
        dk, dn_, df_ = DEEP_RECORD[p]
        print(f"      p={p}: ALIVE at p^{kmax} ({leaves} raw tuples, {frozen} with")
        print(f"           v_{p}(cnum) resolved finite; deep record {p}^{dk}: "
              f"{dn_} tuples, {df_} resolved -- reproduce with --deep)")
    return table


# ------------------------- exact p-adic branch witnesses ----------------------

def _sqrt_padic(n, target, p, K):
    M = p**K
    n %= M
    if n == 0:
        return None, "zero to precision", 0
    v = v_p(n, p)
    if v % 2:
        return None, f"v odd ({v})", 0
    u = n // p**v
    if p == 2:
        if u % 8 != 1:
            return None, f"unit {u % 8} mod 8", 0
        x = 1
        for j in range(3, K - v):
            if (x*x - u) % (1 << (j+1)):
                x += 1 << (j-1)
            if (x*x - u) % (1 << (j+1)):
                return None, "stall", 0
    else:
        if pow(u % p, (p-1)//2, p) != 1:
            return None, "unit non-residue", 0
        x = next(c for c in range(1, p) if (c*c - u) % p == 0)
        for j in range(1, K - v):
            r = (x*x - u) % p**(j+1)
            if r:
                x = (x - (r // p**j) * pow(2*x % p, -1, p) * p**j) % p**(j+1)
    root = x * p**(v//2)
    Mr = p**(K - v//2)
    m = p**(K - v//2 - 2)
    if (root - target) % m and ((-root) - target) % m == 0:
        root = (-root) % Mr
    elif (root - target) % m and ((-root) - target) % m:
        if v_p(((-root) - target) % Mr, p) > v_p((root - target) % Mr, p):
            root = (-root) % Mr
    return root % Mr, v // 2, K - v // 2


def branch_witness(base, u_, t, p, K):
    """Follow the locus branch through a degenerate point at offset a = a0 + u*p^t.
    Returns dict with v_p(cnum) etc., or a string naming the failing step.
    T1/T2 use sqclass_decisive, so True means the mu-square 2nu - a and the
    coincidence -(a+2nu)*disc(m) are CERTIFIED Z_p-squares (resolved even
    valuation, resolved square unit), not merely square-class-consistent:
    a dict with T1 and T2 True is a certified Z_p-point of the full system."""
    a0, nu0, ga0, al0, be0 = base
    M = p**K
    a = a0 + u_ * p**t
    s = 3*a*a + 330*a + 8778
    D = 5*a**4 + 1320*a**3 + 126456*a*a + 5102240*a + 72824400
    y0 = 2*(nu0*nu0) - (3*a0*a0 + 330*a0 + 8778)
    y, l1, prec1 = _sqrt_padic(D % M, y0 % M, p, K)
    if y is None:
        return f"y = sqrt(D): {l1}"
    b = s + y
    if p == 2:
        if b % 2:
            return "s + y odd"
        b //= 2
    else:
        b = b * pow(2, -1, p**prec1) % p**prec1
    nu, l2, prec2 = _sqrt_padic(b % p**prec1, nu0 % M, p, prec1)
    if nu is None:
        return f"nu = sqrt(b): {l2}"
    f = -a**3 - 165*a*a + 2*a*b - 8778*a + 165*b - 172810
    ga, l3, prec3 = _sqrt_padic((-f) % p**prec1, ga0 % M, p, prec1)
    if ga is None:
        return f"gamma = sqrt(-f): {l3}"
    P = min(prec2, prec3, K - 8)
    Mp = p**P
    d_, e_ = -165 - a, a*a + 165*a - b + 8778
    al = al0 % Mp

    def G(al):
        return ((d_ + al*al)**2 - 8*al*ga - 4*e_) % Mp

    def Gp(al):
        return (4*al*(d_ + al*al) - 8*ga) % Mp

    hv = v_p(Gp(al), p) if Gp(al) else 99
    for _ in range(P + 2):
        g = G(al)
        if g == 0:
            break
        vg = v_p(g, p)
        if vg <= 2*hv:
            return f"alpha Newton: v(G) = {vg} <= 2 v(G') = {2*hv}"
        al = (al - (g // p**hv) * pow(Gp(al) // p**hv, -1, Mp // p**hv)) % Mp
    if p == 2:
        if (d_ + al*al) % 2:
            return "beta not integral"
        be = (d_ + al*al) // 2
    else:
        be = ((d_ + al*al) * pow(2, -1, Mp)) % Mp
    p0 = (nu * ga) % Mp
    cn = (p0*p0 - 893025) % Mp
    vC = v_p(cn, p)
    if vC >= P - 4:
        return "cnum unresolved at precision"
    t1 = sqclass_decisive(2*nu - a, p, min(prec2, 30))
    dm = 18*al*be*ga - 4*al**3*ga + al*al*be*be - 4*be**3 - 27*ga*ga
    t2 = sqclass_decisive((-(a + 2*nu)) * dm, p, min(P, 30))
    return {"v": vC, "prec": P, "T1": t1, "T2": t2}


def step3() -> dict:
    print("  [3] exact p-adic branch witnesses: the locus is alive at EVERY depth")
    out = {}
    for p, K, ts in ((2, 64, (14, 18, 22)), (3, 40, (6, 8, 10)),
                     (5, 34, (4, 6, 8)), (7, 30, (3, 5, 7)),
                     (11, 20, (2, 4)), (13, 18, (2, 4)), (31, 16, (1, 3))):
        rows = []
        for t in ts:
            r = branch_witness(WITNESS_BASE, 1, t, p, K)
            check(isinstance(r, dict), f"p={p}, t={t}: witness failed: {r}")
            check(r["T1"] and r["T2"], f"p={p}, t={t}: side squares fail")
            check(r["v"] == t + WITNESS_LAW[p], f"p={p}: law v = t + {WITNESS_LAW[p]}")
            check(r["v"] >= REQ_V.get(p, 0), f"p={p}: witness inside the 10!-window")
            rows.append((t, r["v"], r["prec"]))
        out[p] = rows
        law = WITNESS_LAW[p]
        note = "" if p in REQ_V else " (10!-window imposes nothing at this p)"
        print(f"      p={p}: v_{p}(cnum) = t + {law} exactly at t = {ts}, side squares")
        print(f"           certified, precision to {p}^{rows[0][2]}: c != 0 at every "
              f"depth{note}.")
    print("      No congruence at these primes can force c = 0 on the locus.")
    return out


def step4() -> int:
    print("  [4] the prune trap: an aggressive tree prune fakes a kill here")
    # the pruned tree at a = -10 + 2^12 dies at depth 2^19 (pinned below) ...
    base = (-10, -3, 315, -7, -53)
    t = 12
    a = base[0] + 2**t
    tgt = base[1:]
    sols = [tuple(v % 2 for v in tgt)]
    died_at = None
    for k in range(1, 30):
        M, Mp = 2**(k+1), 2**k
        j = min(k, max(t - 3, 3))
        new = []
        for x in sols:
            for d in product(range(2), repeat=4):
                y_ = tuple((xi + di*Mp) for xi, di in zip(x, d))
                if all(v % M == 0 for v in
                       (lambda a_, x_: [(-a_**4 - 165*a_**3 + 3*a_*a_*(x_[0]*x_[0]) - 8778*a_*a_
                                         + 330*a_*(x_[0]*x_[0]) - 172810*a_ - (x_[0]*x_[0])**2
                                         + 8778*(x_[0]*x_[0]) - 1057221) % M,
                                        (x_[1]*x_[1] - a_**3 - 165*a_*a_ + 2*a_*(x_[0]*x_[0])
                                         - 8778*a_ + 165*(x_[0]*x_[0]) - 172810) % M,
                                        (2*x_[3] + 165 + a_ - x_[2]*x_[2]) % M,
                                        (x_[3]*x_[3] - 2*x_[2]*x_[1]
                                         - (a_*a_ + 165*a_ - (x_[0]*x_[0]) + 8778)) % M])(a, y_)):
                    if all((y_[i] - tgt[i]) % 2**min(k+1, j) == 0 for i in range(4)):
                        new.append(y_)
        sols = new[:600]
        if not sols:
            died_at = k + 1
            break
    check(died_at == 19, "the pruned tree reports a death at depth 2^19")
    r = branch_witness(base, 1, t, 2, 64)
    check(isinstance(r, dict) and r["v"] == t + 14,
          "but the exact branch has a certified point at the same offset")
    print(f"      pruned tree: 'dead at depth 2^{died_at}' -- yet the exact branch")
    print(f"      carries a certified Z_2-point with v2(cnum) = {t + 14} at that very")
    print("      offset.  The prune window was tighter than the branch derivative")
    print("      allows.  Deep batteries can produce FALSE kills; only exact")
    print("      valuation laws or unpruned certificates settle a class.")
    return died_at


def step5() -> None:
    print("  [5] the recorded candidates are off the locus and still die")
    u = sp.Symbol("u")
    t_ = sp.Symbol("t")
    R5 = sp.expand(sp.prod([u - (2*j - 1)**2 for j in range(1, 6)]))
    for cv, kp in ((1395418752000, 11), (2235340800, 13)):
        fl = sp.Poly(R5 - 1024*cv, u).factor_list()[1]
        P_ = 1
        for pp, _ in fl:
            P_ *= int(sp.discriminant(pp.as_poly(u), u))
        check(P_ > 0 and not sp.integer_nthroot(P_, 2)[1], f"c={cv}: D_q*D_k nonsquare")
        co = [int(x) for x in sp.Poly(sp.expand((R5 - 1024*cv).subs(u, t_*t_)), t_).all_coeffs()]
        killed = None
        for p_ in sp.primerange(3, 40):
            if not any(sum(c2 * pow(w, i, p_) for i, c2 in enumerate(reversed(co))) % p_ == 0
                       for w in range(p_)):
                killed = int(p_)
                break
        check(killed == kp, f"c={cv}: kill prime {kp}")
        print(f"      c = {cv}: D_q*D_k nonsquare (case I, not on the locus), killed at {killed}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split(chr(10))[0])
    ap.add_argument("--json_out", type=Path, default=None)
    ap.add_argument("--deep", action="store_true",
                    help="re-derive DEEP_RECORD (raw counts, minutes)")
    args = ap.parse_args()

    print("=== Q28 wall audit: the Q27 battery-depth lesson, applied. Verdict: REAL ===")
    step0()
    step1()
    table = step2(TREE_DEPTH)
    if args.deep:
        deep = step2(DEEP_DEPTHS)
        for p, row in deep.items():
            check(row == DEEP_RECORD[p], f"p={p}: deep run reproduces DEEP_RECORD, got {row}")
        print("      --deep: DEEP_RECORD reproduced exactly.")
    wit = step3()
    step4()
    step5()
    print()
    print("  RESULT: SURVIVE -- certified analytic branches with c != 0 and")
    print("          10!-compatible valuations at every depth, at p = 2, 3, 5, 7")
    print("          and spot-checks 11, 13, 31.  Unlike Q27's free-parameter")
    print("          classes (uniform valuations, secretly empty), this locus is")
    print("          a p-adic curve THROUGH the 160 c = 0 points, and v_p(c)")
    print("          sweeps upward along it: no congruence at these primes can")
    print("          close it.  The wall stands; the section-3 Magma call remains")
    print("          the missing tool.  PROVED_23_POSITIVE stays False.")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps({
            "claim": "Q28's no-kill locus survives deep p-adic batteries: the wall is real",
            "calibration": "Q27 delta=3: alive mod 2^7, dead mod 2^14 (v2 = 13)",
            "tree_alive": {str(p): {"depth": d, "raw_tuples": n, "v_resolved": f}
                           for p, (d, n, f) in DEEP_RECORD.items()},
            "tree_note": "raw 5-tuple counts (mirror-closed sets, not orbits); "
                         "reproduce with: python scripts/k10_deep.py --deep",
            "witness_base": list(WITNESS_BASE),
            "witness_laws": {str(p): f"v_{p}(cnum) = t + {c}" for p, c in WITNESS_LAW.items()},
            "witnesses": {str(p): [[t, v, pr] for (t, v, pr) in rows]
                          for p, rows in wit.items()},
            "prune_trap": "pruned tree fakes a death at a = -10 + 2^12 (depth 2^19); "
                          "exact branch is alive there with v2(cnum) = 26",
            "verdict": "SURVIVE at 2, 3, 5, 7 (spot: 11, 13, 31); "
                       "PROVED_23_POSITIVE unchanged (False)",
        }, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
