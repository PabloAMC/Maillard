# Pre-registration: calibrate on the user's own data (`maillard calibrate`), 2026-09-08

## 1. What it is

A verb that takes a user's scored measurements (the document `maillard score` reads) and writes a
**calibration**: a per-laboratory overlay on the shipped model that the other verbs can apply with
`--calibration`. The shipped parameters, the panel scorecard and every tracked artifact are
untouched; a calibration lives under `results/user/<lab>/` and is never read by a generator.

## 2. The rule it keeps

The repository fits primary evidence and validates on levels. A user's pots are levels. The verb
therefore splits what a level may teach from what only a contrast may teach:

- **Levels fit the laboratory's response factor**: one multiplicative factor per compound, the
  offset between what this laboratory measures and what the model predicts at the shipped
  parameters. It is a property of the measurement chain (extraction, calibration, vessel), not of
  the chemistry, so a level may set it.
- **Contrasts fit the kinetics**: log-ratios between the user's own pots (different time,
  temperature, pH, recipe), from which the response factor cancels. Only these move rate
  constants, and only the coordinates the contrasts can identify (a Jacobian rank test), each
  pulled toward its shipped value by its shipped uncertainty (a maximum-a-posteriori fit, not a
  free one). Everything else stays exactly at the shipped value and the card says so.

## 3. What runs

1. Records are read from the measurement document. Each record is tagged `role: fit` or
   `role: validate`; untagged records with four or more in the document are split by holding out
   every second record in sorted order of (temperature, time, pH, recipe), chosen before any fit
   and written to the card. Fewer than four: no hold-out, and the card says the calibration is
   unvalidated.
2. Stage one: response factors from the fit records' levels at the shipped parameters.
3. Stage two: contrasts between consecutive fit records that measure the same compound; candidate
   coordinates are the lane's fitted coordinates with a shipped uncertainty (the envelope's priors
   table); the Jacobian of the contrast residuals selects the identifiable ones; a bounded
   least-squares fit with the prior as a residual moves them; the posterior uncertainty comes from
   the same Jacobian.
4. Stage one again, at the fitted kinetics.
5. The validate records are scored at the shipped parameters and with the calibration; the card
   reports both, per compound and as a median fold error.

The calibration records the laboratory, the base wave, the matrix of the records, the factors with
their sigma and row counts, the overrides with prior and posterior, the fit and validate record
names, and provenance. Applying it to a spec whose matrix differs from the records' prints a
warning on the answer.

## 4. What counts as success, declared before the code runs

- **T1, synthetic recovery (unit tier, a stand-in model).** With a known threefold response factor
  and a known shift in one coordinate, the two stages recover both within their reported sigma,
  and touch no other coordinate.
- **T2, synthetic recovery (scientific tier, the real engine).** Measurements generated from the
  engine at 145 °C with one sulfur coordinate shifted and a response factor of three on MFT:
  the factor is recovered within 0.1 dex and the shifted coordinate within twice its posterior
  sigma; the hold-out median fold error after calibration is below the one before.
- **T3, the hold-out is never fitted.** A runtime guard records every record a residual read;
  a validate record in that set is a test failure.
- **T4, the shipped numbers do not move.** No file under `scripts/generators` or `src/kinetic_core`
  other than the calibration and user-scoring modules names `results/user`; the panel scorecard
  computed with the calibration modules imported equals the tracked one.
- **T5, a real laboratory as a test.** Yiltirak 2026's four-temperature ladder as a pretend
  laboratory: the 100 and 120 °C pots fit, the 110 and 130 °C pots validate. The hold-out median
  fold error with the calibration must be below the one without. The sign of the temperature trend
  is not expected to flip: a response factor cannot fix a structural miss, and the card must
  say which of the two it did.

Failure of T1 to T4 means the verb does not ship. Failure of T5 ships with the number stated,
since it is a test of the model, not of the verb.

## 5. What it will not do

It will not move the shipped parameters, write under `data/` or `results/validation/`, fit a
coordinate from levels, fit a coordinate the contrasts cannot identify, or widen a prior. The
envelope's interval does not yet carry the factor's sigma; the card prints it.

## 6. Outcome (2026-09-08, first run)

Shipped: `src/kinetic_core/calibration.py` (the overlay), `src/kinetic_core/user_fit.py` (the two
stages), the `calibrate` verb and `--calibration` on `compare`, `predict` and `score`;
`tests/unit/test_calibration.py` (a stand-in model) and `tests/scientific/test_calibration_engine.py`
(the engine).

- **T1 passed** on the stand-in: a threefold factor and a 0.4 dex shift recovered, the other
  coordinates untouched, the hold-out never read.
- **T2 passed with a lesson.** The first synthetic set shifted the norfuraneol-to-thiol rate at
  145 °C, and the fit left it to the response factor: at that temperature the rate scales MFT almost
  uniformly, so no contrast can see it, and a uniform offset *is* a response factor by the rule in
  section 2. The set was changed to a time series with the FFT decay rate shifted, which changes
  the shape. With the coordinate named by the caller (`--coordinate`), the fit recovers the factor
  (3.41 for a true 3.0) and the shift (−0.14 ± 0.06 for a true −0.20, the prior pulling back as
  designed); the hold-out median fold error goes from 2.4 to 1.03, six of six rows within 3x. Left
  to choose among all sixteen calibratable sulfur coordinates, the fit reaches the same hold-out
  (1.02) through two other coordinates: one time series at one temperature cannot tell a sink rate
  from the osone and cysteine rates. The card lists the chosen and the unchosen, and the verb gained
  `--coordinate` so a laboratory that knows where its pot differs can say so.
- **T3 passed**: the runtime log of records read during the fit never contains a validate record.
- **T4 passed**: no generator or engine module names `results/user`; the plain engine and the
  engine under a null calibration return the same concentrations.
- **T5 passed, and says what a factor can and cannot do.** Yiltirak 2026 as a pretend laboratory
  (100 and 120 °C fit, 110 and 130 °C validate): the factors are 0.04 on MFT and 0.002 on FFT, the
  model being 25 to 500 times too high for that laboratory; the contrasts moved the furfural decay
  and the furfural-to-FFT rates; the hold-out median fold error falls from 115 to 1.2 and three of
  four rows land within 3x. The one that does not is MFT at 130 °C: measured 1.71, calibrated 6.5,
  because the measured ladder falls with temperature and the calibrated model still rises. A
  response factor moves the level; it cannot fix the sign of a trend, and the card shows that.

Cost: 75 to 130 engine evaluations per calibration, under a minute. Not yet done: the envelope does
not carry the factor's sigma (the card prints it); the overlay type is not yet shared with a matrix
declaration (roadmap, section 6).
