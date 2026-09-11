# Pre-registration: wave B37, sixteen sources arrive (written 2026-09-11, BEFORE any probe was run)

## 1. Why

Wave B36 read seven papers and ended by naming, in the PR record and in `docs/guides/EXPERIMENTS.md`,
the measurements this model still lacks. The owner then downloaded **sixteen** papers against that
list in one go. This wave reads all sixteen, writes a dossier for each, and states — before running
anything — what each one can and cannot do to the model.

| file | paper | what it was fetched for |
|---|---|---|
| `mittelmaier2010.pdf` | Mittelmaier et al. 2011, Anal. Bioanal. Chem. 399:1689 | **the fed 3-deoxyglucosone experiment**: pure 3-DG in water at 120 °C, pH 5, with 3,4-DGE followed against time |
| `frankel1993.pdf` | Frankel 1993, JAOCS 70:767 | **hexanal-specific activation energies** for thermal decomposition of oxidized oils |
| `tazi2009.pdf` | Tazi et al. 2009, Food Chem. 115:958 | **a measured Q10 table** for lipoxidation in a food matrix, 60–130 °C, three water activities |
| `yu2012.pdf` | Yu, Meng, Ramaswamy & Boye, J. Food Process. Preserv. | **measured residence time** of an SPI feed in a twin-screw extruder |
| `chen2010.pdf` | Chen, Wei, Zhang & Ojokoh 2010, J. Food Eng. 96:208 | the same, for pure SPI |
| `ebert2021.pdf` | Ebert et al. 2022, J. Sci. Food Agric. | pea isolate **and** its texturates: the unheated column |
| `Zhou2026.pdf` | Zhou et al. 2026, Front. Nutr. | soybean flour → meal → isolate → extrudate by isotope dilution |
| `Xu2024.pdf` | Xu et al. 2024, Food Chem. 437:137924 | soy isolate heated against an unheated control |
| `Kong2024.pdf` | Kong et al. 2024, Food Chem. 445:138795 | the same, with an untreated control column |
| `Liu2025.pdf` | Liu et al. 2025, Food Biophys. 20:190 | a soy isolate fraction from unheated to 180 °C |
| `cai2021.pdf` | Cai et al. 2021, Food Chem. 340:127880 | roasted soybeans with an unroasted control |
| `zhang2012.pdf` | Zhang et al. 2012, JAFC | soymilk raw against UHT |
| `moisio2015.pdf` | Moisio et al. 2015, Eur. Food Res. Technol. | rye bran extruded over a temperature series |
| `Wang2026b.pdf` | Wang et al. 2026, Food Chem. 521:149913 | xylose–cysteine thiols and disulfides against time |
| `Zhai2023c.pdf` | Zhai et al. 2023, Food Chem. 404:134420 | the same, fed with deoxyosone fragments |
| `hernandez2023.pdf` | Hernandez et al. 2023, Molecules 28:3151 | the source of the PBMA identity row, previously off disk |

**This wave fits nothing and moves no constant.** It reads, records, corrects provenance, and runs
declared probes whose predictions are written below before they are run.

## 2. What was seen while reading, declared here because reading precedes this file

1. **Mittelmaier runs the exact experiment `EXPERIMENTS.md` asks for.** ~200 µM of pure 3-DG,
   0.5 mL, sealed gas-tight vial, glucose-free PD model at pH 5, 120 °C, sampled at 0/10/20/30/60/90/
   120 min. 3,4-DGE peaks at **26.7 µM at 30 min** — 13.4 % of the charge. The time courses
   themselves are Fig. 5 and are figure-only; six maxima are printed as numbers.
2. **And it shows the step is reversible.** Feeding 3,4-DGE regenerates 3-DG (26.9 µM at 30 min).
   It also shows a sink the model does not have: 3-deoxygalactosone, the C4 epimer, which reaches
   **26 % of the fed 3-DG by 60 min**. The trunk's `r_tdg_ddg` is one-way and has no epimer.
3. **The engine cannot be charged with 3-deoxyglucosone.** `PRECURSOR_ALIASES` carries methylglyoxal,
   glyoxal, glucosone, diacetyl, norfuraneol and the Amadori compounds — wave B13's "dicarbonyl
   trio" — but not TDG, DDG or ODG. The fed experiment cannot be run against the model as the code
   stands.
4. **Frankel prints a hexanal-specific activation energy**: n-6 seed oils at 130–160 °C, 27.2 / 29.2 /
   27.3 kcal/mol (113.8 / 122.2 / 114.2 kJ/mol). **And it prints the reason a single barrier will not
   do**: the apparent Ea rises with the decomposition window in every oil, which the author reads as
   a mechanism change.
5. **Tazi prints a Q10 table** resolved by temperature and water activity: 3.3 → 2.4 at aw 0.38 and
   **2.0 → 1.6 at aw 0.72** across 60–70 → 120–130 °C, with Ea 114 → 65 kJ/mol (chemiluminescence)
   and 100 → 62 (TBARS). The markers are bulk lipoxidation extent, not hexanal.
6. **Chen 2010's residence times are figure-only**; Yu's are a full 18-run table, but on a 20 %-SPI
   corn-flour feed rather than the 90 %-SPI feed the model's extrusion row describes.

## 3. Predictions, written before the probes are run

- **P1 — the fed 3-DG pot, once it can be charged.** Adding `"3-deoxyglucosone": "TDG"` to
  `PRECURSOR_ALIASES` is inert (no bundle, benchmark or claim charges it). Charging 200 µM 3-DG at
  120 °C, pH 5, water, and integrating to 120 min, I predict:
  **(a)** the model's peak 3,4-DGE is **below** the printed 26.7 µM, by a factor between **5 and 50**
  — the bracket set by B34's 32× at 121 °C and B36's 7–10× at 90–110 °C;
  **(b)** the model's 3,4-DGE rises **monotonically** across 0–120 min rather than peaking at 30 min,
  because the trunk has no reverse reaction, no epimer sink, and a slower exit than entry.
  If (b) is wrong — if the model does peak — the mechanism is closer than claimed and I will say so.
- **P2 — the lipid lane's Q10 against a measurement.** `Q10_ASSUMPTION` is a constant with a declared
  default and band, referenced at 25 °C. Evaluated at the midpoints of Tazi's four intervals, I
  predict the model's **default sits above the measured Q10 at aw 0.57 and 0.72 at every interval**,
  and that the gap **widens with temperature**, because a constant Q10 cannot reproduce the fall that
  an Arrhenius barrier produces.
- **P3 — Frankel's barrier expressed as a Q10.** Converting 113.8–122.2 kJ/mol at 145 °C, I predict
  the implied Q10 lands **between 2.0 and 2.5**, i.e. inside the model's 2–3 band at cooking
  temperature while the same barrier implies a Q10 above 4 at the model's 25 °C reference. The point
  of the probe is the discrepancy between those two readings of one band, not the band itself.
- **P4 — nothing moves.** No constant, no fit row, no target, no benchmark value. Predicted headline
  after the run: **unchanged** — panel 10/45, out-of-sample 9/44, refused 32, hold-out 5/31, envelope
  16/44 with 1 not evaluable. If any count moves, this wave stops and reports it.

## 4. What is built

1. **Sixteen extraction dossiers**, one per PDF.
2. **The provenance correction for `resconi_2023_pbma_beef_identity_benchmark`**, whose vessel note
   says its source is not on disk. It is. Corrected through the generator, prior claim retained and
   labelled, as B34, B35 and B36 did.
3. **One inert code change**: 3-deoxyglucosone becomes chargeable, so that a fed-dicarbonyl pot can
   be expressed at all. Guarded by a test that no existing pot charges it and that every frozen
   artifact is unchanged.
4. **The three probes above, recorded with their predictions and outcomes.**
5. **A list of what each paper unlocks and what wave would have to do it** — none of it done here,
   because each is a fit or a ship decision that needs its own pre-registration and, in the lipid
   lane's case, a frozen before/after pair. Naming them is the deliverable; doing them is not.

## 5. Outcome (written 2026-09-11, after the probes)

**Three predictions held, one was refuted, and the refuted one is the more interesting.**

### P1 — the fed 3-deoxyglucosone pot

Charged 200 µM 3-DG at 120 °C, pH 5, water, and integrated. The model's 3,4-DGE:

| t (min) | model 3,4-DGE (µM) | % of charge | model 3-DG left (µM) |
|---:|---:|---:|---:|
| 5 | 4.24 | 2.12 | 138.8 |
| **10** | **5.28** | **2.64** | 96.3 |
| 20 | 4.15 | 2.08 | 46.4 |
| 30 | 2.49 | 1.24 | 22.3 |
| 60 | 0.35 | 0.17 | 2.0 |
| 120 | 0.005 | 0.002 | 0.03 |

- **P1(a) HELD, at the optimistic edge.** The model's peak is **5.28 µM against the printed 26.7 µM
  — 5.05× low**. The prediction was "between 5 and 50", and the answer landed on the boundary. It is
  the third independent measurement of the same deficit, and the *smallest* of the three: 32× on
  Leitzen 2021 at 121 °C from glucose, 7–10× on Zhang 2021 at 90–110 °C from glucose, 5× here from
  pure 3-DG. **The deficit shrinks as the pot gets closer to the step.** That ordering is itself
  evidence: part of what B34 read as a slow `k_tdg_ddg` is upstream of 3-DG, not in the step.
- **P1(b) REFUTED.** I predicted the model's 3,4-DGE would rise monotonically because the trunk has
  no reverse reaction and no epimer sink. It does not. **It peaks at about 10 minutes and then
  collapses** — to 0.002 % of the charge by 120 minutes, three and a half orders below its own peak.
  The paper's pot still holds 3,4-DGE at 120 min and its 3-DG is only partly gone. **The model's
  error is not that it makes too little 3,4-DGE; it is that it destroys the whole 3-deoxy pool far
  too fast.** At 60 minutes the model has 1 % of the fed 3-DG left, where the paper still has enough
  for a quarter of it to have epimerised to 3-deoxygalactosone. A fit that moved `k_tdg_ddg` alone
  to chase the peak height would be fitting the wrong constant.

### P2 — the lipid lane's Q10 against a measured one

`Q10_ASSUMPTION` is a constant **2.449** (the geometric mean of its declared 2–3 band), referenced at
25 °C. Against Tazi 2009's measured table:

| interval midpoint (°C) | model | measured, aw 0.38 (CL/TBARS) | aw 0.57 | aw 0.72 | model above the moist rows? |
|---:|---:|---|---|---|---|
| 65 | 2.45 | 3.3 / 2.8 | 2.2 / 1.9 | 2.0 / 1.9 | yes |
| 85 | 2.45 | 2.9 / 2.5 | 2.1 / 1.8 | 1.8 / 1.8 | yes |
| 105 | 2.45 | 2.6 / 2.3 | 1.9 / 1.7 | 1.7 / 1.7 | yes |
| 125 | 2.45 | 2.4 / 2.1 | 1.8 / 1.6 | 1.6 / 1.6 | yes |

**HELD on both counts.** The model's default is above every measured value at aw 0.57 and 0.72, and
the gap widens monotonically with temperature (0.45 → 0.65 → 0.75 → 0.85 against the wettest column).
The model's band bottom of 2.0 is still above the measured 1.6 at cooking temperature in a moist
matrix.

### P3 — Frankel's barrier read as a Q10

| oil, 130–160 °C | Ea | Q10 at 145 °C | Q10 at 25 °C |
|---|---:|---:|---:|
| soybean | 113.8 kJ/mol | **2.15** | 4.44 |
| safflower | 122.2 kJ/mol | **2.27** | 4.95 |
| canola | 114.2 kJ/mol | **2.15** | 4.46 |

**HELD** — the prediction was 2.0 to 2.5. And the probe's real point, which is the reciprocal: the
model's constant Q10 of 2.449 corresponds to **Ea = 68.4 kJ/mol if read at its 25 °C reference and
118.1 kJ/mol if read at 120 °C**; its band bottom of 2.0 spans 52.9 to 91.3 kJ/mol and its top of 3.0
spans 83.9 to 144.8. **A constant Q10 is not a constant barrier.** Frankel's measured hexanal barrier
(113.8–122.2 kJ/mol, bulk oil) sits inside the 120 °C reading and far above the 25 °C one; Tazi's
moist-matrix barrier (61–65 kJ/mol) sits inside the 25 °C reading and far below the 120 °C one. The
two measurements do not agree with each other, and the model's formulation is precisely what hides
the disagreement.

### P4 — nothing moved

**HELD.** No constant, no fit row, no target, no measured value. The 3-deoxy aliases are inert:
no bundle, benchmark, claim or fit row charges any of them, and a test holds it.

## 6. What the sixteen unlocked, and which wave would have to do it

Named here, done nowhere, because each is a fit or a ship decision that needs its own
pre-registration:

1. **The 3-DG limb, refit against a fed pot** (Mittelmaier). The finding above says the target is the
   3-deoxy pool's total lifetime, not the single step. Needs a reversible step, or a declared reason
   not to have one, before any constant moves.
2. **The lipid lane's temperature term** (Frankel and Tazi). Replacing a constant Q10 with a barrier
   would move every lipid row, so it needs a frozen before/after pair and a ship rule, like ENV-B34.
3. **The extrusion residence time** (Yu). The fit row's 25 s is unsourced and 40 s is measured on a
   different feed; changing a FIT row's conditions invalidates a frozen report.
4. **A dry-roast temperature ladder** (Cai): an unroasted zero and four temperatures, with pyrazines
   and furfural quantified — the only external Maillard series with a genuine blank.
5. **A soymilk raw-versus-UHT row** (Zhang 2012): the best-quantified unheated/heated pair of the
   sixteen, and two compounds that rise rather than fall.
6. **The thiol sink** (Zhai): a printed seven-point time course at 100 °C with both thiols measured
   free, on a TTCA charge the engine already accepts.

And two requests, because the numbers exist and did not arrive: **Table S1 of Wang et al. 2026**
(the three canonical dimers as concentrations) and **Supplementary Table 1 of Liu et al. 2025** (a
four-temperature moist-heat ladder with a blank).
