# Wave B9 pre-registration — the fit / validate split

*Written 2026-09-03 before any B9 number existed. Generator:
`scripts/generators/generate_kinetic_core_b9_fit.py`.*

## 1. The rule this wave applies

Owner's principle (2026-09-03): **primary evidence fits the model; causally secondary evidence
validates it.** Rate constants, activation energies, fed-intermediate yields, conversions and
within-study ratios are causes and belong in the objective. The concentration of a target odorant in a
full precursor system is the consequence the model exists to predict, and belongs in validation.

B8's objective (62 rows) contains eight rows that are consequences and are also scored panel bundles:
the Hofmann & Schieberle 1998 Table 1 levels of FFT and MFT at pH 5 for ribose, xylose, glucose and
fructose + cysteine (`hofmann_{ribose,xylose,glucose,fructose}_{FFT,MFT}`). Step B3 found this; step 2
of the post-retirement plan declared it (`kinetic_core_b8_fit_targets.json`) and moved the xylose bundle
off the hold-out panel. B9 removes the cause.

## 2. What changes, and what does not

| | B8 | B9 |
| --- | --- | --- |
| objective rows | 62 | **54** (the eight level rows removed) |
| network, T-structure, species | B8 | unchanged |
| free set | 23 coordinates | unchanged |
| declared bands (`full_bounds`) | B8 | **unchanged**, and every bound the optimum lands on is listed in the report as `active_bounds` |
| pH weighting | W-HALF, E = 0.70 | unchanged |
| optimiser / budget / starts | `least_squares` trf, 600 evaluations, 2 starts | unchanged; start 0 = B8's frozen optimum, start 1 = the perturbation protocol |
| rows kept that are also levels | — | the in-situ NF and furan-2-aldehyde levels in the ribose pot (Hofmann T5): they are intermediates no bundle scores |

Nothing else is edited. The B8 generator is frozen (`WAVES.md`); B9 imports it and removes rows.

## 3. What B9 gives back

The four Hofmann pH-5 bundles become **out-of-sample** for the kinetic core. The xylose bundle returns
to the hold-out panel (`data/benchmarks/external_validation/maillard_path/`); the ribose, glucose and
fructose bundles stay on the trust loop but are no longer flagged `in_core_fit`. The scorecard's
out-of-sample row count rises from 40 to 49, i.e. every scored row.

## 4. Expectations and the falsifier, stated in advance

- **Expected:** the eight returned rows will score worse than under B8, because B8 fitted them. B8's
  in-fit rows were 4 of 9 within 3x; a fair expectation for B9 on the same eight rows is 2 to 4 within 3x.
- **Expected:** the hold-out panel (16 bundles) should not get worse than B8's 3 of 28 rows within 3x
  by more than one row; the removed rows carried information mostly about the same conditions the
  hold-out probes (145 °C, 20 min, pH 5).
- **Expected:** the absolute scale stays identified from the fed-intermediate rows (yields per mmol fed);
  the Laplace covariance at the B9 optimum should identify at least as many coordinates as B8's 18 of 23.
- **The falsifier:** if the four returned bundles score **0 of 8** within 3x AND the hold-out panel loses
  two or more rows, the fed-intermediate rows do not carry the absolute scale and the split has cost more
  than it bought. B9 then does NOT ship as the engine's sulfur report; the report is kept, the finding is
  recorded, and the next question is which level measurement would be a legitimate fit anchor (a
  time course, not an end point).
- **Bands:** B9 does not widen any band. If the optimum lands on the same bounds as B8 (profile: nine
  active), the report says so and the owner decides per bound whether it is a measurement (Gigl 2021's
  thiol-sink range) or a search convenience. That decision is a wave of its own.

## 5. What is read and what is not

The fit reads its row table (literals in the frozen B8 / B2.3 generators) and nothing under
`data/benchmarks/`. The hold-out guard's fit-isolation check (check 2) runs on the B9 generator like
any other. The scorecard, envelope and directional artifacts are regenerated AFTER the fit is frozen,
never during it.
