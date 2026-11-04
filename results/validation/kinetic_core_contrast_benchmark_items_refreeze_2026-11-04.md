# Contrast benchmark items re-frozen (2026-11-04): v6 items with the shipped engine's current refusals

- Why: the scorer re-derives the item sets and refuses to score if they differ from the committed record. Two things
  changed since the items were frozen (79138999, 2026-11-03):
  1. PR #47 (CF-106) made the acrylamide lane refuse programmes peaking inside 60 - 120 C. Atlas item 56 (the
     acrylamide temperature leg: Asn 200 + Glc 100 mM, 0.1 M phosphate pH 7, 100 vs 130 C, 21 min) moves from
     answerable to refused ("ACRYLAMIDE BELOW ITS FIT WINDOW (CF-106)"). Atlas coverage: answerable 4 -> 3,
     refused_other 1 -> 2. The keyed family, the S items (31) and the detection curve are unchanged.
  2. M1 (EX-061) added two Config keys, off by default (`fit.dimer_encounter_ceiling`, `fit.mft_decay_prior`), so
     the full Config digest changes (40b74647... -> 09d122f9...); the digest over the pre-existing keys is unchanged.
- No rule, item definition, prereg or hand table changed. The frozen v1 score (`kinetic_core_contrast_benchmark_score_v1`)
  stays the "before" record; it scored item 56 as answerable, which an acrylamide-only item cannot affect for M1's
  sulfur-lane side effects (sec. 8 of the M1 prereg).
