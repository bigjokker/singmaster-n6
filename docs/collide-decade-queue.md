# The skipped collide decade: the command queue (prepared 2026-08-24)

**Status 2026-08-25.** Tonight pack (`l=20..13`, 127 pairs, ~17 min) and
day pack (`l=12,7,11,6`, 25 pairs, ~6.0 h) both finished: **152 json,
0 hits, 0 discoveries**, every `m_range` matching the engine gap
`[collide_frontier_m(l), old-start−1]`, digits 61 (`l<10`) / 101
(`l>=10`). Artifacts: `results/collide_gapdecade_k*_l*.json`.

Remaining: **l=10** (8 pairs, ~20 h at 8 workers). No bat; the finished
tonight/day wrappers were removed.

```text
python scripts\collide_gapdecade.py --pack l10 --workers 8
```

Existing json are skipped. Do not start a second copy. Do not paste the
PowerShell loop below — ranges are computed live in
`scripts/collide_gapdecade.py`.

**The hole (D5):** every recorded `results/collide_*.json` starts at the
*second* past-2017 value decade — 62-digit values for `l < 10`, 102-digit for
`l >= 10` — because three call sites used `min_m_for_digits(l, cap + 1)`.
The classifier calls the *first* decade (61 / 101 digits) past the bound, and
`collide_frontier_m(l)` now abuts it (pinned by `scripts/test_collide.py`).
That first decade has never been scanned as a job, for any pair.

**Coverage check (2026-08-24):** every recorded `m_range` low end equals the
old frontier `min_m_for_digits(l, cap+1)` exactly, for all 160 recorded
(k, l) pairs — so the gap is exactly `[collide_frontier_m(l), recorded_start - 1]`
per pair, nothing more, and the commands below cannot overlap a recorded
range (`--max-m` is `recorded_start - 1`).

**Cost model:** q11 §4 rate, ~1.6e5 (k,pair)-checks/s serial. `--workers W`
divides wall-clock by roughly W (multiprocessing pool). Times below are
SERIAL; divide by your worker count.

## The queue, cheapest tier first (stop after any tier; each is complete)

| tier | l | pairs | m-range (per pair) | checks | serial time |
|---|---|---|---|---|---|
| 1 | 20 | k=2..19 (18) | 830,446..931,774 | 1.8M | 11 s |
| 2 | 19 | k=2..18 (17) | 1,453,370..1,640,617 | 3.2M | 20 s |
| 3 | 18 | k=2..17 (16) | 2,714,483..3,084,909 | 5.9M | 37 s |
| 4 | 17 | k=2..16 (15) | 5,473,925..6,267,902 | 11.9M | 74 s |
| 5 | 16 | k=2..15 (14) | 12,093,138..13,964,936 | 26.2M | 2.7 min |
| 6 | 9 | k=2..8 (7) | 19,249,445..24,861,612 | 39.3M | 4.1 min |
| 7 | 15 | k=2..14 (13) | 29,814,904..34,761,624 | 64.3M | 6.7 min |
| 8 | 8 | k=3,5,6,7 (4) | 119,039,222..158,741,352 | 158.8M | 16.5 min |
| 9 | 14 | k=2..13 (12) | 84,006,892..99,024,687 | 180.2M | 18.8 min |
| 10 | 13 | k=2..12 (11) | 279,070,468..333,147,804 | 594.9M | 62 min |
| 11 | 12 | k=2..11 (10) | 1,139,448,627..1,380,473,525 | 2.41B | 4.2 h |
| 12 | 7 | k=2..6 (5) | 1,259,932,332..1,750,670,296 | 2.45B | 4.3 h |
| 13 | 11 | k=2..10 (9) | 6,052,339,028..7,461,606,434 | 12.68B | 22.0 h |
| 14 | 6 | k=5 (1) | 29,937,951,658..43,942,903,516 | 14.0B | 24.3 h |
| 15 | 10 | k=2..9 (8) | 45,287,286,886..57,013,316,290 | 93.8B | 162.9 h |

Tiers 1–10 (everything through l=13): **~1.9 h serial**. Through tier 14
(everything but l=10): **~56.6 h serial, ~7 h at --workers 8**. Tier 15 —
the 101-digit decade at l=10 — is the monster: **~163 h serial, ~20 h at
--workers 8**; it is 74% of the whole hole. Grand total 126.4B checks,
219.5 h serial.

Settled pairs are excluded (l=8 skips k=2,4; l=6 is k=5 only; l>=9 has no
settled pairs). k < l always. Each (k,l) writes its own json under a name
no recorded artifact uses (`collide_gapdecade_k{k}_l{l}.json`).

## The command (paste into PowerShell; runs the tiers in the order above)

```powershell
$decade = @{
  20 = @(830446, 931774);        19 = @(1453370, 1640617)
  18 = @(2714483, 3084909);      17 = @(5473925, 6267902)
  16 = @(12093138, 13964936);     9 = @(19249445, 24861612)
  15 = @(29814904, 34761624);     8 = @(119039222, 158741352)
  14 = @(84006892, 99024687);    13 = @(279070468, 333147804)
  12 = @(1139448627, 1380473525); 7 = @(1259932332, 1750670296)
  11 = @(6052339028, 7461606434); 6 = @(29937951658, 43942903516)
  10 = @(45287286886, 57013316290)
}
$pairs = @{
  20 = 2..19; 19 = 2..18; 18 = 2..17; 17 = 2..16; 16 = 2..15
   9 = 2..8;  15 = 2..14;  8 = @(3,5,6,7); 14 = 2..13; 13 = 2..12
  12 = 2..11;  7 = 2..6;  11 = 2..10;  6 = @(5); 10 = 2..9
}
$order = @(20,19,18,17,16,9,15,8,14,13,12,7,11,6,10)
foreach ($l in $order) {
  $lo = $decade[$l][0]; $hi = $decade[$l][1]
  foreach ($k in $pairs[$l]) {
    python singmaster_intersect.py collide --k $k --l $l --min-m $lo --max-m $hi --workers 8 --json_out "results/collide_gapdecade_k${k}_l${l}.json"
  }
}
```

To run only one tier, set `$order = @(<that l>)`. To stop before the l=10
monster, drop `10` from `$order`.

**Still unscanned after this queue, deliberately not proposed:** the (3,5)
and (4,5) pairs — their frontier m = 4,128,917,917,336 exceeds
`COLLIDE_HARD_SKIP_M` = 1e11, and their first decade alone is ~2.4e12
checks (~4,100 h serial); every recorded run skipped them entirely. And all
l > 20. Those are jobs of a different order, not this hole.

A hit from any tier is a new past-2017 collision — the headline event. A
clean sweep extends the recorded nulls down to the exact 2017 boundary,
closing D5's "abuts but unscanned" caveat pair by pair.
