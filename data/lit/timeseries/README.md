# `data/lit/timeseries/` — harvested kinetic time series

Concentration-vs-time data for Maillard model systems, harvested from primary sources
so that the reaction-network model can eventually be calibrated and tested against
*measured trajectories* rather than against endpoints or against numbers of unclear origin.

---

## Read this first

**This changed on 2026-08-27 (audit Wave S3). Two of these files are now a FIT CORPUS.**

The previous version of this section said "nothing in this directory is wired into the
model, the calibration, the registry, the priors, or any test", and required this README to
change if that ever stopped being true. It has stopped being true, and here is the exact
scope of it:

| File | Status |
|---|---|
| `martins2005_glucose_glycine_80_100_120C_pH68.yml` | **IN THE FIT CORPUS** (glucose, DFG, 3-DG, 1-DG series) |
| `martins2003_DFG_amadori_degradation.yml` | **IN THE FIT CORPUS**, pH 6.8 DFG series only; the pH 5.5 series is held out |
| `martins2005_glucose_glycine_100C_pH68.yml` | not fitted — it is the same experiment as the 100 °C series of the file above, plotted in a second figure. Used only to measure the read-off reproducibility floor |
| `brands_sugar_casein_120C_pH68.yml` | not fitted — endpoints, not trajectories |

What consumes them: `scripts/generators/generate_trunk_rate_calibration.py`, which fits the
eight-step trunk in `src/trunk_kinetics.py` and writes
`results/validation/trunk_rate_calibration_refit.{json,md}`. The fit is machine-declared in
that artifact's `fit_target_files` block and is read by `scripts/ci/fit_target_gate.py`
(which was extended in the same wave to strip `.yml`/`.yaml`, not only `.json`, from a
declared fit target — before that a YAML corpus declaration was accepted but inert).

**What is still true, and matters more:** *no shipped prediction changed.* The fitted
constants live in a dedicated calibration lane that no module under `src/` imports —
`tests/scientific/test_wave_s3_trunk_kinetics.py` asserts exactly that, so the day someone
wires it into a prediction path, a test fails and says so. The individual file headers still
say "NOTHING IN THIS FILE IS WIRED INTO THE MODEL OR THE CALIBRATION"; for the two fitted
files that sentence is now **superseded by this table**, and the fit is the reason.

**This repository has been burned by invented literature values**, which is why the
paragraph below still stands and why the fit was built to be auditable rather than
convenient.

**This repository has been burned by invented literature values.** See `AUDIT.md` at the
repository root: a 2026-08 audit found roughly 30–45% citation contamination across the
reference set, anchors tuned to benchmarks that did not exist, and a validation headline
that was circular. That history is why every file here carries its own provenance block,
why "unverifiable" is an acceptable recorded status, and why the verification pass
described below re-derived the numbers *independently* instead of re-reading the harvest.

**Do not promote anything here into calibration without reading its `residual_uncertainty`
block.** Three of the four files are figure read-offs. A figure read-off recovers where the
authors' plotting program placed a marker — not the number in the authors' spreadsheet.

---

## Verification status

Verified 2026-08-27 (audit Wave Q). Every file was checked against the **retrieved source
document**, not against the harvest notes.

| File | System | Points | Access | Status | Evidence |
|---|---|---:|---|---|---|
| `brands_sugar_casein_120C_pH68.yml` | sugar/casein, 120 °C, pH 6.8 | 60 | `full_text` | **verified_table** | Brands (2002) WUR thesis, Tables 6.1 / 6.2 / 3.1, re-read from 4× image renders |
| `martins2003_DFG_amadori_degradation.yml` | DFG degradation, 100/120 °C, pH 5.5/6.8 | 106 | `mixed` | **verified_table + verified_figure_shape** | Martins (2003) WUR thesis, Table 4.1.1 exact; Figs 4.1.1/4.1.2 re-extracted |
| `martins2005_glucose_glycine_100C_pH68.yml` | glucose/glycine, 100 °C, pH 6.8 | 183 | `figure_readoff` | **verified_figure_shape** | Martins (2003) thesis Fig 5.9; all 183 points reproduced exactly |
| `martins2005_glucose_glycine_80_100_120C_pH68.yml` | glucose/glycine, 80/100/120 °C, pH 6.8 | 418 | `figure_readoff` | **corrected** | Martins (2003) thesis Fig 5.10; 415/418 reproduced, **3 phantom points removed** |

Point counts include replicate markers listed separately, `t0_markers`, `ambiguous_*`
blocks, and the 32-value glycine-yield table.

### Citations

All six DOIs referenced across these files resolve on CrossRef with exact
title/journal/volume/page/author agreement:

| DOI | Work |
|---|---|
| `10.1016/j.foodchem.2004.04.006` | Martins & van Boekel 2005, Food Chem 90(1-2):257-269 |
| `10.1016/s0008-6215(03)00173-3` | Martins, Marcelis & van Boekel 2003, Carbohydr Res 338(16):1651-1663 (Part I) |
| `10.1016/s0008-6215(03)00174-5` | Martins & van Boekel 2003, Carbohydr Res 338(16):1665-1678 (Part II) |
| `10.1021/jf9907586` | Brands, Alink, van Boekel & Jongen 2000, JAFC 48:2271-2275 |
| `10.1021/jf010789c` | Brands, Wedzicha & van Boekel 2002, JAFC 50:1178-1183 |
| `10.1021/jf001430b` | Brands & van Boekel 2001, JAFC 49:4667-4675 (referenced as a *target*, not a source) |

**On the missing DOI in the Brands file.** The Brands file has no top-level DOI because its
source is a PhD thesis, which has none; it is openly downloadable at
<https://edepot.wur.nl/199005>. That is not a gap. The two thesis chapters supplying every
number in it were published as journal articles and both of those DOIs are recorded and
verified. Note specifically that these data are **not** from Brands & van Boekel (2001)
`10.1021/jf001430b` — that is a *different* chapter of the same thesis, and the natural next
retrieval target, not the origin of anything currently in the file.

---

## Corrections made in Wave Q

**Data corrections (one file):**

- `martins2005_glucose_glycine_80_100_120C_pH68.yml`, glycine 120 °C: removed
  `[10, 196.1]`, `[15, 179.5]`, `[20, 185.1]`. An exhaustive marker census of Figure 5.10
  panel B found 66 markers at 48 distinct positions; the series claimed 15 values at 5
  sampling times where only 12 exist. Each removed value duplicated an adjacent time point,
  i.e. they were padding to make every time show three replicates. **This was the only data
  defect found in 767 values.**

**Data additions:**

- `brands_sugar_casein_120C_pH68.yml`: added the six values of the `a/b'` microanalysis
  column of Table 3.1, which the original harvest had skipped.

**Metadata corrections:**

- Brands thesis ISBN `90-5808-579-4` → `90-5808-591-0`; unsupported "137 pp" removed
  (the scan is 133 PDF pages ending at printed p. 127, and no page count is stated).
- Brands Table 6.1 page 101 → 99; Table 3.1 page 36 → 37.
- Brands scan described as "200 dpi" → measured 1-bit, ≈1310×1889 px/page (≈160 dpi).
- Martins thesis ISBN `90-5808-923-4` → `90-5808-823-5` (both Martins files).
- DFG file: the excluded non-data geometry was described as "a filled-diamond series";
  no diamond markers exist on that page. The fits are stroked polylines. Exclusion was
  correct; the description was not.
- Fig 5.10 file: the dotted fit curve whose dots mimic data is in **panel C** (3-DG at
  80 °C), not panel B.
- Fig 5.10 file: read-off agreement with Fig 5.9 restated from "under 1% of full scale"
  to the measured worst cases (up to 2.2%, for DFG at 60–120 min).

**Verified and left alone:** every other value. Notably, the DFG file's flag that the thesis
Materials and Methods contains an internal inconsistency (`0.237 g, 10 mmol` — 0.237 g of
DFG is 1.0 mmol, not 10) was confirmed verbatim: the error is the thesis's, correctly
flagged rather than silently fixed.

---

## How the figure files were verified

The three Martins files are read-offs from vector figures in the thesis PDF. Verification
did **not** re-read the harvest; it re-derived the numbers with a separate content-stream
parser sharing no code or calibration with the original pass:

1. Marker classification by **path geometry**, not bounding-box heuristics — filled closed
   paths give diamond/triangle/circle/square; `+`, `×` and `*` are drawn as 4 or 6 short
   segments radiating from a common start point, and that shared start point *is* the datum,
   recovered with zero centring error.
2. Axis calibration re-derived from the drawn plot frames and their tick marks, including
   the **dual axes** in Figure 5.10 (fructose and DFG are on right-hand scales — the single
   easiest way to misread that figure).
3. PDF clipping paths (operator `n`) and model-fit polylines identified and excluded.
4. Marker **multiplicity** recorded per position. This is what made the three phantom
   glycine values detectable: a position with zero drawn markers cannot be data.
5. Cross-check: Figure 5.9 and the 100 °C series of Figure 5.10 are the same experiment
   plotted twice. They agree at every shared time point for all ten species.

Independent corroboration of qualitative claims came from the thesis body — e.g. the
formic-vs-acetic assignment in Fig 5.9 panel C, which the caption cannot settle, is fixed by
the text: *"acetic acid was always formed in higher concentrations than formic acid"*.

---

## What a human should pull next

Highest value first:

1. **Brands & van Boekel (2001), `10.1021/jf001430b`.** The multiresponse concentration data
   for glucose/fructose-casein at 120 °C (glucose, fructose, formic acid, acetic acid, lysine
   loss, melanoidins vs time) exist in the thesis only as figures in a bitonal scan, which is
   unreadable. Someone with ACS access should check whether the publisher PDF has vector
   figures; if so the same extraction used on the Martins figures recovers the whole set.
   This is the single biggest gap: the Brands file currently has browning and mutagenicity
   endpoints, not trajectories.
2. **Martins thesis Chapter 4 figures 4.1.3–4.1.7** (PDF pp. 78–82) and the Chapter 4.2
   model-fit figures (PDF pp. 93–105): 1-DG, 3-DG, methylglyoxal, glucose, mannose, formic
   and acetic acid and melanoidins during DFG degradation, at four conditions. These are
   vector figures in the same PDF and extract with the identical method. They were not done
   for time, **not** because they are inaccessible.
3. **The 90 °C and 110 °C glucose/glycine series.** The thesis ran five temperatures and
   fitted rate constants to all five, but plots only three. Recovering the other two turns
   this into a five-temperature Arrhenius data set. Requires the authors' raw data from
   Wageningen Food Quality and Design.
4. **The underlying concentration tables for Figures 5.9 and 5.10.** Never published, in
   either the thesis or the Food Chemistry paper. Until someone obtains them, the vector
   extraction here is the best available reconstruction — and it is a reconstruction.

---

## Conventions

- `access:` — `full_text` (source retrieved and read), `figure_readoff`, `mixed`,
  or `unverifiable`.
- `data_quality:` per series — `table` (transcribed digits) or `figure_readoff`.
- `verification:` per series — `status` is one of `verified_table`,
  `verified_figure_shape`, `corrected`, `unverifiable`; plus `date` and a `note` saying
  what was checked against what.
- `file_verification:` — the same, for the file as a whole, plus `residual_uncertainty`.
- Replicate markers are listed **separately and never averaged**. Where replicates were
  drawn at identical coordinates they collapse to one value, so the number of values at a
  time point is a *lower bound* on the number of replicates.
- Corrections cite where the correct value came from (which table or figure of which
  document, and which PDF page).
