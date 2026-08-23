#!/usr/bin/env python3
"""Regression tests for the collide frontier (D5).

`collide` scans C(n,k)=C(m,l) "past the 2017 bound": BBW 2017 settles every
value with at most VALUE_DIGITS_2017 = 60 digits (100 for l >= 10), so
`pair_is_past_2017` calls a value past the bound when num_digits(C) > cap.
The scan therefore has to START at the first C with cap+1 digits, i.e. at
C >= 10**cap, i.e. at min_m_for_digits(l, cap).

It started one decade later. Three call sites each used
min_m_for_digits(l, frontier_digits_for_l(l) + 1), which is the first value
>= 10**(cap+1) -- the first (cap+2)-digit value -- so the whole (cap+1)-digit
decade, which the classifier itself calls past the bound, was never scanned
by any run. Every recorded collide artifact shows value_digits starting at
62 (l < 10) or 102 (l >= 10).

The pins here fail on that code and pass once the three sites share one
helper, `collide_frontier_m(l)`, that abuts the classifier. Nothing here
runs a long collide: the run-level pin scans 21 values of C(m,l).

Run: python scripts/test_collide.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import singmaster_intersect as si  # noqa: E402

ok: list[str] = []
errors: list[str] = []


def expect(cond: bool, msg: str) -> None:
    (ok if cond else errors).append(msg)


ROW = si.ROW_BOUND_2017 + 1          # any n past the 2017 row bound


def test_frontier_start_abuts_the_classifier() -> None:
    """For each regime, the frontier m is the first m whose C(m,l) the
    classifier calls past 2017 -- and its predecessor is not past."""
    fn = getattr(si, "collide_frontier_m", None)
    expect(callable(fn), "singmaster_intersect exposes collide_frontier_m(l)")
    if not callable(fn):
        return
    for l, cap in ((3, si.VALUE_DIGITS_2017), (10, si.VALUE_DIGITS_2017_LGE10)):
        expect(si.frontier_digits_for_l(l) == cap,
               f"fixture: frontier_digits_for_l({l}) == {cap}")
        m = fn(l)
        m_ref = si.min_m_for_digits(l, cap)
        m_old = si.min_m_for_digits(l, cap + 1)
        d_m = si.num_digits(si.comb(m, l))
        d_prev = si.num_digits(si.comb(m - 1, l))
        expect(m == m_ref,
               f"l={l}: frontier m == min_m_for_digits({l}, {cap}) "
               f"(got {m}, ref {m_ref})")
        expect(d_m == cap + 1 and d_prev == cap,
               f"l={l}: C(m,{l}) has {cap+1} digits and C(m-1,{l}) has {cap} "
               f"(got {d_m} / {d_prev})")
        expect(m < m_old,
               f"l={l}: the frontier start is strictly below the old "
               f"min_m_for_digits({l}, {cap+1}) = {m_old} -- the skipped decade")
        # the predicate the scan must abut, both sides of the line
        expect(si.pair_is_past_2017(ROW, 2, m, l, si.comb(m, l)) is True,
               f"l={l}: pair_is_past_2017 is True at the frontier value")
        expect(si.pair_is_past_2017(ROW, 2, m - 1, l, si.comb(m - 1, l)) is False,
               f"l={l}: pair_is_past_2017 is False on the {cap}-digit predecessor")


def test_run_collide_starts_in_the_first_past_decade() -> None:
    """A tiny past-2017 run on an unsettled pair must report value_digits
    whose low end is cap+1 (61 / 101), not cap+2 (62 / 102)."""
    fn = getattr(si, "collide_frontier_m", None)
    for k, l in ((3, 5), (2, 10)):
        expect((k, l) not in si.SETTLED_KL, f"fixture: ({k},{l}) is an unsettled pair")
        cap = si.frontier_digits_for_l(l)
        m_front = fn(l) if callable(fn) else si.min_m_for_digits(l, cap)
        hits, certs = si.run_collide(k, l, 2 * l, m_front + 20,
                                     workers=1, past_2017_only=True)
        c = certs[0]
        expect(c.get("m_range", [None])[0] == m_front,
               f"({k},{l}): past-2017 run starts at the frontier m={m_front} "
               f"(got m_range {c.get('m_range')})")
        expect(c.get("value_digits", [None])[0] == cap + 1,
               f"({k},{l}): value_digits low end is {cap+1}, not {cap+2} "
               f"(got {c.get('value_digits')})")
        expect(c.get("past_2017_only") is True and c.get("settled_pair") is False,
               f"({k},{l}): certificate records past_2017_only=True, unsettled")
        expect(not hits, f"({k},{l}): no collision in the 21 values scanned")


def test_sanity_unaffected() -> None:
    """The engine's own suite must still pass with the frontier moved."""
    import io
    import contextlib

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rep = si.run_sanity()
    expect(rep.get("passed") is True and not rep.get("errors"),
           f"singmaster_intersect sanity still passes "
           f"({len(rep.get('ok', []))} checks, {len(rep.get('errors', []))} errors)")


def main() -> int:
    test_frontier_start_abuts_the_classifier()
    test_run_collide_starts_in_the_first_past_decade()
    test_sanity_unaffected()
    print("\n=== COLLIDE TESTS ===")
    for line in ok:
        print("  OK   ", line)
    for line in errors:
        print("  FAIL ", line)
    print("  RESULT", "PASS" if not errors else "FAIL")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
