# Pre-registration: wave B36, seven sources arrive (written 2026-09-11, BEFORE any bundle was edited or any scorer run)

## 1. Why

The 2026-09-11 reading-list row in `scripts/generators/WAVES.md` named twelve bundles whose cited
source was not on disk and gave the owner nine DOIs. The owner downloaded seven PDFs:

| file | paper | bundles that cite it |
|---|---|---|
| `data/articles/Ma2024.pdf` | Ma, Fu, Cheng & Liu 2024, Int. J. Mol. Sci. 25:8668 | `acrylamide_spi_extrusion_130C_ACSRef3` (a FIT row) |
| `data/articles/li2026.pdf` | Li, Dai, Mao, An, Bai & Kaur 2026, Foods 15:912 | `external_validation_li_2026_spi_wg_hme_control` |
| `data/articles/lin2021.pdf` | Lin, Chan, Kao & Sung 2021, Polymers 13:1901 — the paper the three "Chang 2021" bundles score | `mp_holdout_glucose_asparagine_180C_{10min,30min,30min_water}_Chang2021` |
| `data/articles/lin2022.pdf` | Lin, Ting, Ndraha, Hsiao & Sung 2022, Polymers 14:1565 | `mp_holdout_fructose_asparagine_180C_Lin2022` |
| `data/articles/ye2024.pdf` | Ye et al. 2024, Foods 13:2836 | `mp_holdout_glucose_asparagine_180C_Ye2024` |
| `data/articles/singh2021.pdf` | Singh, Shi, Magreault, Kitts, Jarzębski, Siejak & Pratap-Singh 2021, Molecules 26:4104 | `pea_isolate_40C_PratapSingh2021`, `soy_isolate_40C_PratapSingh2021` |
| `data/articles/koelsch1991.pdf` | Koelsch, Downes & Labuza 1991, J. Food Sci. 56:816 | none (a lipid-rate candidate named by the search) |

Three named sources did **not** arrive and their bundles keep their "not on disk" notes, which are
true: Zhang et al. 2021 (Food Sci. Nutr. 9:290, the aqueous 3-DG/3,4-DGE time courses B34 asked for) **[WRONG -- see the correction at the end of §5: it has been on disk since 2026-09-07 as `zhang2020.pdf`]**,
Schibilsky 2019 (TU Berlin dissertation, two hold-out bundles) and Hernandez et al. 2023 (Molecules
28:3151, the PBMA identity row).

This wave does one thing: **reads the seven papers and checks every number and every condition the
nine bundles transcribed second-hand against the print.** No constant moves. Nothing is fitted.

## 2. What was already seen while reading, before this file was written

The papers were read before this pre-registration could be written, so the following were seen and
are declared here rather than presented as predictions.

1. **Li 2026, Table 2, nonanal.** The bundle scores 72.66 µg/kg, and its 2026-08-27 correction note
   says "the Nonanal row reads 72.66 ± 1.46". The printed table (page 12, rendered and read as an
   image) has **74.37 ± 0.11** in the HMPE-0 min control column and 72.66 ± 1.46 in the HMPE-20 min
   column. The August correction fixed the row (decanal → nonanal) and took the wrong column. The
   other three scored values are in the control column as they should be: hexanal 605.64 ± 6.50,
   1-hexanol 20.04 ± 0.66, 2-pentylfuran 5625.80 ± 63.75.
2. **Ma 2024 prints no residence time.** The bundle's vessel note says "the bundle's 25 s is the
   residence time". The Materials and Methods give ten heating zones, screw speed, feed rates and
   moisture, and no residence time anywhere in the paper. The 25 s is an assumption and will be
   labelled one.
3. **Ma 2024, Figure 2D (read as an image).** The 130 °C bar reads ≈150 µg/kg, the 150 °C bar ≈120,
   the 170 °C bar ≈82 and the **unextruded control ≈38 µg/kg**. The bundle's 150 is confirmed. The
   control bar is a starting level the bundle does not declare (see §4, item 4).
4. **Ye 2024 never states the reactant molarity.** Already recorded in the bundle's
   `precursor_concentration_provenance` (0.2 M taken from the paper's cited method, Knol 2005); the
   print confirms the absence.
5. **Every other transcribed number matches the print**: Singh 2021 Table 1 hexanal (pea 1138.00 ±
   297.30, soy 1621.71 ± 159.69) and the prose 2-pentylfuran hexanal-equivalents (pea 638 ± 49, soy
   2492 ± 199); Lin 2021 prose (28, 912, 1459 ppb; 832 ppb; HMF 7 ppm); Lin 2022 prose (1859 ppb; HMF
   12.28 ppm); Ye 2024 prose (140.58 ± 13.92 µmol/mol Asn). Every verbatim quotation the bundles
   carry from the Europe PMC XML is in the PDF word for word.

## 3. Predictions, written before any generator runs

- **P1 (provenance).** After the wave, no bundle whose source is on disk says "NOT ON DISK" except in
  a retained, labelled prior note. The vessel provenance class of the Ma, Li, three Chang, Lin 2022
  and Ye bundles moves from `repo_verbatim_methods_quote` to `primary_source_pdf`; the Singh bundles
  stay `not_applicable` (no cook); the Li buffer stays `buffer_unknown` because the print gives no
  medium pH for the extrusion blend.
- **P2 (the one value change).** Nonanal on `external_validation_li_2026_spi_wg_hme_control` moves
  72.66 → 74.37 µg/kg, a factor 1.0235 on that row's ratio. No 3× verdict can flip unless the row's
  current ratio lies within 2.35 % of 3 or of 1/3. Predicted headline counts after the run:
  **unchanged** — panel 10/45, out-of-sample 9/44, refused 32, hold-out 5/31, envelope 16/44 with 1
  not evaluable, families headspace 9 / extraction 36. If any count moves, the wave stops and says so.
- **P3 (no fit input changes).** No constant moves, no fit row's value or conditions change, the
  frozen fit reports are untouched, `fit_target_gate` and `holdout_guard` stay green.
- **P4 (the acrylamide fit row).** Its 150 µg/kg stays; the 25 s residence time stays as a labelled
  assumption; the zone profile from the print is recorded in the vessel note and nowhere executable,
  because a benchmark bundle is an isothermal hold (`panel.py`) and the row already scores 4 247× low.
- **P5 (nothing added).** No new target: the candidates are named in §4 with the reason each is not
  a validation.

## 4. What is built

1. **Seven extraction dossiers**, one per PDF, so that every number below has a dossier:
   `ma2024_extraction.md`, `li2026_extraction.md`, `lin2021_extraction.md`, `lin2022_extraction.md`,
   `ye2024_extraction.md`, `singh2021_extraction.md`, `koelsch1991_extraction.md`.
2. **Provenance notes corrected through their generators** (`complete_benchmark_vessel_fields.py`,
   `complete_benchmark_buffer_fields.py`), the false claims retained and labelled superseded, as B34
   and B35 did. The Ma note gains the print's zone profile (first five zones 80/80/85/90/100 °C; the
   130 °C arm sets the last five to 110/120/130/130/130 °C; 150 rpm; 6.0 kg/h raw + 2.57 kg/h water;
   30 % moisture; SPI:corn starch 9:1).
3. **The nonanal correction** on the Li 2026 bundle with a dated correction note, and the bundle's
   frozen hash in `tests/unit/test_kinetic_core_b2_3.py` re-pinned with the reason.
4. **Named and not acted on.** Ma 2024's unextruded control already carries ≈38 µg/kg acrylamide
   (Figure 2D). Declaring it as a carried level is the same cure Trikusuma took, but this is a FIT row
   whose B3 report is frozen, the value is a figure read, and on a row that scores 4 247× low a 25 %
   shift of the formed amount changes no verdict. Recorded in the dossier and in
   `docs/guides/EXPERIMENTS.md` as the row's second named debt.
5. **Candidates the trunk carries and does not get**, with the reason: Li 2026's control column also
   prints 2,5-dimethylpyrazine 134.30 ± 4.36 µg/kg, but the bundle charges no free amino acid and the
   engine's reachability rule would refuse the target by name, so it would be a refusal row, not a
   validation; furfural is "–" (not detected) in the control with no detection limit printed; Singh
   2021's hexanol is n.d. for pea and soy. The lipid products in both tables that the trunk does not
   carry (heptanal, benzaldehyde, 1-octen-3-ol, 2-nonanone) are listed in the dossiers.
6. **One test amended.** `test_kinetic_core_b2_3.py`'s rule for `buffer_unknown` demanded the words
   "NOT ON DISK" in the note. That conflated "the paper is absent" with "the paper is silent". It now
   accepts either, so the Li 2026 buffer note can say the true thing: the paper is on disk and does
   not state the medium.

## 5. Outcome (written 2026-09-11, after the run)

**Every prediction held.**

- **P1 held.** Nine bundles' notes now open with the on-disk reading; the old claims sit after them,
  labelled superseded. Six vessel blocks are `primary_source_pdf`; the Singh pair stay
  `not_applicable`; the Li buffer stays `buffer_unknown` with the note saying the print is silent.
- **P2 held.** Nonanal on the Li 2026 row moved 72.66 → 74.37 µg/kg. Its fold error moved 3.04× →
  3.11× against a prediction of 23.9 µg/kg: outside the threefold band before and after, so the row's
  verdict did not change. Headline counts after the run: panel **10/45**, out-of-sample **9/44**,
  refused **32**, hold-out 5/31, envelope **16/44** with 1 not evaluable — all as predicted. The one
  visible README change is the lipid-lane median on the protein-matrix line, 2.79× → 2.82×.
- **P3 held.** No constant moved; `fit_target_gate` and `holdout_guard` green; the frozen fit
  reports untouched.
- **P4 held.** The acrylamide row keeps 150 µg/kg and 130 °C / 25 s; the zone profile and the
  "no residence time" finding live in its vessel note and dossier.
- **P5 held.** No target added.

**One thing the pre-registration did not foresee.** Running the buffer-note generator to apply the
corrections regressed **seven** bundles whose buffer blocks B34, B35 and the Yiltirak reading had
edited in place without updating the generator (Steinhagen, liu_2023, four Yiltirak rows, li_2026).
Caught by `git status` before anything was committed. The generator now carries those blocks
verbatim (`_EDITED_IN_PLACE`), its `--check` reports drift instead of always passing, and a unit test
runs it — the guard the vessel generator has had since R1 and the buffer generator never had. The
rule for the future is one line: a buffer or vessel note is edited in its generator, never in the
bundle.

**What the seven papers did not give.** No rate, no barrier, no time course usable by the trunk:
Koelsch 1991 is one temperature; the Lin and Ye papers keep their time courses in figures; Ma 2024's
acrylamide is bars. The reading list's top item, Zhang 2021 (doi 10.1002/fsn3.1995), did not arrive
and remains the one download that could move a constant. **[WRONG, corrected below.]**

**Correction, 2026-09-11, same day.** CORRECTION 2026-09-11 (same day, on the owner's word): THIS PAPER HAS BEEN ON DISK SINCE 2026-09-07 as `data/articles/zhang2020.pdf`, dossier `zhang2020_extraction.md` (the file is named by its received date, the venue year is 2021), and the audit matched it by year and missed it -- the same error B34 found on the Steinhagen row. The dossier already records that its absolute levels do not mass-balance and that only orderings and time shapes are used (directional claim DIC-01). So the claim that it "did not arrive" was false twice over: it had arrived four days earlier and had already been read. What it can give the aqueous ask is not a rate but a within-study ratio, and it was probed the same day (the engine on 0.3 M glucose in water, pH 6.5, 6 h, the DIC-01 pot; nothing fitted, nothing changed):

| T (°C) | model 3-DG (µg/L) | model 3,4-DDG (µg/L) | model ratio 3,4-DDG / 3-DG | Zhang 2021 Table 1, 6 h, same ratio | model ÷ paper |
|---:|---:|---:|---:|---:|---:|
| 90 | 5264 | 106 | 0.020 | 0.212 | 0.095 |
| 95 | 7062 | 168 | 0.024 | 0.231 | 0.103 |
| 100 | 9298 | 260 | 0.028 | 0.241 | 0.116 |
| 105 | 11928 | 391 | 0.033 | 0.259 | 0.126 |
| 110 | 14659 | 560 | 0.038 | 0.277 | 0.138 |

The paper's ratio is unit-free: both quinoxalines come from one derivatisation and one LC-MS/MS run, so the calibration error that breaks its absolute levels cancels if it is common to the two, which is the assumption this ratio rests on. Read that way, **a second aqueous laboratory says the same thing B34 found on Leitzen 2021: the model's 3-deoxyglucosone → 3,4-dideoxyglucosone step is too slow — 7–10× at 90–110 °C over 6 h here, about 30× at 121 °C over 18 min there** (where 3-DG itself is right to 11 %). The temperature trend disagrees too: the paper's ratio rises 1.3× from 90 to 110 °C, the model's 1.9×. Under the owner's rule a within-study ratio is FIT evidence, so the next wave is a refit of `k_tdg_ddg` (and its barrier) against these two ratios and the Leitzen hold-out, pre-registered before any constant moves. Not done here.
