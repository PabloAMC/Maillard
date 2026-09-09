# Cheng, Rouseff, Li & Wu 2020 — EXTRACTION (methanethiol in heated mandarin juice; methionine → methional → MeSH / DMDS / DMTS in a model juice at 100 °C, with and without ascorbic acid)
### A juice paper. The one time-series experiment (0–60 min at 100 °C) is reported as three chromatograms; what survives as numbers is the juice levels, the model-juice recipe, and the direction of the ascorbate effect.

**Source on disk:** `data/articles/cheng2020.pdf` (owner's download, 2026-09-08). Read from the `pdftotext`
text layer in the scratchpad; Tables 1–3 are clean except four Table 2 rows where the cultivar code fused
with the first digit of the next column (e.g. "Dafen D3 4 ± 31 8 ± 0.8 47" = Dafen, D, 34 ± 3, 18 ± 0.8,
47 %) — resolved by the printed % loss/gain column in every case. Figure 3 (the methionine-degradation
chromatograms) is FIGURE-ONLY. The Supporting Information (PCA plot, calibration slopes and R²) is not
on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Methanethiol, an Off-Flavor Produced from the Thermal Treatment of Mandarin Juices: A Study of Citrus Sulfur Volatiles" |
| Authors | Yujiao Cheng, Russell Rouseff (Southwest Univ. Chongqing / Univ. Florida), Guijie Li, Houjiu Wu |
| Venue | J. Agric. Food Chem. 2020, 68, 1030–1037 |
| DOI | 10.1021/acs.jafc.9b06647 |
| System | 27 mandarin juices (10 cultivars) fresh vs heated 60–90 s at 100 °C; model juices with 1 mM methionine; 1 mM aqueous precursors |
| Analytes | H2S, COS, MeSH, CS2, DMS quantified in juice; MeSH, DMDS, DMTS, methional followed in the model juices (Fig. 3 only) |

## 1. Why it matters

Programme 6 (roadmap §5c) needs methional's chain to methanethiol and on to dimethyl disulfide and
dimethyl trisulfide. This paper is the only one on disk that runs that chain from methionine at 100 °C
in a food-like aqueous matrix and varies the one ingredient that turns out to dominate it — ascorbic
acid. Its qualitative result is sharp: methionine in water alone gives a little MeSH, DMDS and DMTS and
no methional; in a sugar/citrate model juice methional is the main product with less MeSH/DMDS/DMTS;
with 1.42 mM ascorbic acid all four appear at their highest levels in a quarter of the heating time.
The authors read ascorbic acid / dehydroascorbic acid as a ready-made α-dicarbonyl Strecker agent.
For B17's oxidant finding the same ingredient is also a pro-oxidant (Chin & Lindsay 1994), and this
design cannot separate the two roles. In the real juice, 60–90 s at 100 °C takes MeSH from
not-detected to 279–917 ng/L while H2S falls 3–36 %.

## 2. Methods as they matter to a model

- **Juice heating:** 5 mL juice in a 20 mL glass vial with 1 cm stir bar, Teflon-faced silicone septum
  cap (sealed; ~15 mL air headspace), water bath at **100 °C**, 350 rpm, **60 s** (initial samples) or
  **90 s** ("to mimic thermal processing"), then ice bath. Which samples got 60 s vs 90 s is not stated
  per cultivar.
- **Model juice MJ1** (per 100 mL, Milli-Q water): fructose 2.5 g (139 mmol/L), glucose 2.5 g (139 mmol/L),
  sucrose 5.0 g (146 mmol/L), citric acid 1.0 g (52 mmol/L), tripotassium citrate 0.5 g (16 mmol/L, MW
  306.4), **methionine 14.9 mg = 0.999 mmol/L** (MW 149.21). **pH not stated** (a 52:16 citric
  acid:citrate mix is roughly pH 3 — my estimate, not printed).
- **MJ2** = MJ1 + **25 mg ascorbic acid / 100 mL = 1.42 mmol/L** (MW 176.12).
- **Control:** 1 mM methionine in water.
- **Model heating:** ~5 mL in 20 mL vial, 100 °C water bath, 350 rpm; MJ1 and control at **0, 1, 5, 10,
  15, 30, 60 min**; MJ2 at **0, 1, 5, 10, 15 min**; ice bath. ⚠ These time series were run but **no
  numbers from them are printed**; the paper shows one chromatogram per matrix (Fig. 3: control 60 min,
  MJ1 60 min, MJ2 15 min).
- **Precursor screen:** 1 mM aqueous solutions of methionine, S-methylmethionine sulfonium (MMS),
  cysteine (10 mM in the text), thiamine, glutathione; 15–60 min at 100 °C.
- **Quantification:** static headspace SPME. Juice: 5 mL in 20 mL vial, 40 °C / 20 min equilibration,
  2 cm DVB/CAR/PDMS fibre 30 min. Models: 1 cm fibre, 15 min. Internal standards EMS (8.42 mg/L stock)
  and IPDS (0.943 mg/L stock), 1.5 µL into 5 mL. GC Agilent 7890B, Rtx-Wax 60 m × 0.25 mm × 0.25 µm,
  35 °C 6 min → 203 °C at 7 °C/min, hold 10 min; effluent split to a **pulsed flame photometric
  detector** (sulfur mode, square-root output) and an MSD in SIM (two ions per compound, Table 1).
  **External calibration:** eight-level standards of MeSH, DMS, DMTS and methional added to the
  sugar/citrate model juice (no methionine); COS and H2S estimated against the EMS response (no
  standards). Triplicate. ⚠ DMDS is listed as identified (Table 1) and named in the abstract as
  quantitated, but it is **not in the calibration list and has no column in Table 2**.
- **Identification:** LRI on wax and DB-5 within 1 % of standards on both PFPD and SIM, plus SIM ion
  ratio.
- **Sensory:** triangle test, 18 assessors, Satsuma juice spiked with MeSH to 680.3 ng/L (OAV 34):
  10/18 correct (p < 0.05); 9 of 10 preferred the control.
- **Conversions used below:** MeSH MW 48.11 (1 ng/L = 0.0208 nmol/L); H2S 34.08; DMS 62.13; CS2 76.14;
  COS 60.07.

## 3. Tables re-typed

### Table 1. "Identification and Characteristics of Sulfur Volatiles in Mandarin Juices"

Columns: boiling point (°C); LRI DB-5 sample / standard; LRI WAX sample / standard; SIM ions (m/z);
identification; odor descriptor; odor threshold (µg/L; a = in air, w = in water); reference.

| compound | bp | DB-5 sample / std | WAX sample / std | ions | ID | descriptor | threshold |
|---|---:|---|---|---|---|---|---|
| COS | −50 | — / — | 479 / — | 44, 60 | PFPD, SIM | repulsive | 55 a |
| H2S | −60 | 371 / — | 536 / — | 34, 33 | PFPD, SIM | rotten eggs | 10 w |
| MeSH | 6 | 429 / 425 | 676 / 671 | 48, 47 | PFPD, SIM, Std | rotten cabbage | **0.02 w** |
| CS2 | 46 | 551 / 553 | 721 / 721 | 76, 44 | PFPD, SIM, Std | sweet, chemical | 210 a |
| DMS | 37.3 | 523 / 514 | 736 / 734 | 62, 47 | PFPD, SIM, Std | cabbage, sulfurous | 0.33 w |
| DMDS | 110 | 724 / 725 | 1081 / 1084 | 94, 79 | PFPD, SIM, Std | onion, garlic | 12 w |
| DMTS | 170 | 909 / 910 | 1396 / 1410 | 126, 79 | PFPD, SIM, Std | onion, cabbage | 0.01 w |
| methional | 165 | 963 / 963 | 1457 / 1462 | 104, 48 | PFPD, SIM, Std | cooked potato | 0.2 w |

(The DB-5 LRI cells for COS and H2S are blank in print; the paper says DB-5 did not resolve COS from H2S.)

### Table 2. "Average Concentrations (±SD) of Sulfur Volatiles in Different Mandarin Juice Cultivars before and after Heating" — ng/L; ND = not detected; bold in print = Satsuma group (D, N, Y, I)

| cultivar | code | COS before / after / % loss | H2S before / after / % loss | MeSH before / after / % gain | CS2 before / after / % gain | DMS before / after / % gain |
|---|---|---|---|---|---|---|
| Clementine | C | ND / ND / 0.0 | 1187 ± 101 / 1003 ± 69 / 16 | ND / **477 ± 14** / 100 | ND / ND / 0.0 | ND / 464 ± 30 / 100 |
| Gonggan | G | ND / ND / 0.0 | 513 ± 33 / 397 ± 15 / 23 | ND / **279 ± 16** / 100 | 142 ± 3 / 385 ± 46 / 171 | 768 ± 38 / 1026 ± 50 / 34 |
| Nanfengmiju | F | ND / ND / 0.0 | 380 ± 12 / 271 ± 18 / 29 | ND / **307 ± 3** / 100 | 177 ± 10 / 305 ± 2 / 72 | 454 ± 47 / 706 ± 41 / 56 |
| Zhoupigan | Z | 25 ± 0.6 / 24 ± 1 / 4.9 | 597 ± 14 / 514 ± 31 / 14 | ND / **453 ± 20** / 100 | ND / 103 ± 2 / 100 | ND / 392 ± 43 / 100 |
| Xinshengxi ponkan | P | ND / ND / 0.0 | 866 ± 79 / 763 ± 39 / 12 | ND / **385 ± 16** / 100 | ND / 66 ± 0.6 / 100 | 166 ± 27 / 471 ± 16 / 184 |
| Xingyidahongpao | X | ND / ND / 0.0 | 751 ± 98 / 596 ± 50 / 21 | ND / **400 ± 19** / 100 | 48 ± 6 / 175 ± 6 / 264 | 438 ± 22 / 1097 ± 117 / 150 |
| Dafen (Satsuma) | D | 34 ± 3 / 18 ± 0.8 / 47 | 1178 ± 29 / 841 ± 22 / 29 | ND / **825 ± 24** / 100 | 140 ± 3 / 412 ± 41 / 194 | 1719 ± 42 / 6488 ± 430 / 277 |
| Nichinan (Satsuma) | N | 39 ± 0.3 / 24 ± 2 / 38 | 1453 ± 124 / 928 ± 87 / 36 | ND / **917 ± 7** / 100 | 116 ± 1 / 133 ± 4 / 15 | 4355 ± 208 / 7667 ± 945 / 76 |
| Yanjiangmigan (Satsuma) | Y | 9 ± 1 / ND / 100 | 1703 ± 77 / 1446 ± 81 / 15 | ND / **686 ± 42** / 100 | 55 ± 7 / 59 ± 3 / 7 | 5397 ± 368 / 8287 ± 117 / 54 |
| Iwasaki wase (Satsuma) | I | ND / ND / 0.0 | 930 ± 43 / 901 ± 38 / 3 | ND / **379 ± 12** / 100 | 41 ± 9 / 129 ± 2 / 216 | 4557 ± 108 / 9130 ± 591 / 100 |

Text figures on the same data: MeSH average in heated juice **511 ng/L** (= 10.6 nmol/L); heated ranges
H2S 271–1446, DMS 392–9130, MeSH 279–917 ng/L; DMS average 1785 → 3573 ng/L (+87 %); Satsuma DMS
average 7893 vs 693 ng/L for the other six; "total sulfur volatiles almost doubled". DMDS and DMTS
"were only observed in storage samples or samples heated for longer than 90 s or in commercial samples".

### Table 3. "Average OAV of Sulfur Volatiles ... before and after Heating" (concentration / threshold from Table 1)

| cultivar | code | COS b / a | H2S b / a | MeSH b / a | CS2 b / a | DMS b / a |
|---|---|---|---|---|---|---|
| Clementine | C | ND / ND | <1 / <1 | ND / 23.9 | ND / ND | ND / 1.4 |
| Gonggan | G | ND / ND | <1 / <1 | ND / 14 | <1 / <1 | 2.3 / 3.1 |
| Nanfengmiju | F | ND / ND | <1 / <1 | ND / 15.4 | <1 / <1 | 1.4 / 2.1 |
| Zhoupigan | Z | <1 / <1 | <1 / <1 | ND / 22.6 | ND / <1 | ND / 1.1 |
| Xinshengxi ponkan | P | ND / ND | <1 / <1 | ND / 19.2 | ND / <1 | <1 / 1.4 |
| Xingyidahongpao | X | ND / ND | <1 / <1 | ND / 20 | <1 / <1 | 1.3 / 3.3 |
| Dafen | D | <1 / <1 | <1 / <1 | ND / 41.2 | <1 / <1 | 5.2 / 19.7 |
| Nichinan | N | <1 / <1 | <1 / <1 | ND / 45.8 | <1 / <1 | 13.2 / 23.2 |
| Yanjiangmigan | Y | <1 / ND | <1 / <1 | ND / 34.3 | <1 / <1 | 16.3 / 25.1 |
| Iwasaki wase | I | ND / ND | <1 / <1 | ND / 19 | <1 / <1 | 13.9 / 27.7 |

Averages in text: MeSH OAV 25.5 (Satsuma 35.1, others 19.2); DMS OAV 10.8 (Satsuma 23.9, others 2.1).

### Figure 3 — FIGURE-ONLY (three PFPD chromatograms; no axis values read)

Printed captions and the text's reading of them, verbatim in substance:

| panel | system | heating | what the text says |
|---|---|---|---|
| 3A | 1.0 mM methionine in water | 100 °C, 60 min | "only limited amounts of methanethiol, DMDS, and DMTS are formed. No methional is formed." |
| 3B | 1.0 mM methionine in MJ1 (sugars + citric acid) | 100 °C, 60 min | "methional is preferentially formed with reduced amounts of methanethiol, DMDS, and DMTS" |
| 3C | 1.0 mM methionine in MJ1 + 1.42 mM ascorbic acid | 100 °C, **15 min** | "produces the largest amounts of MeSH, DMDS, DMTS, and methional in only 1/4 the amount of heating time" |

Other precursor results (text only): thiamine (1 mM, 30 min) → "small but measurable" H2S and
2-methylfuran-3-thiol; cysteine 10 mM, 30 min → trace CS2 and H2S only; methionine → no measurable DMS;
MMS → DMS "almost exclusively", "the most thermally unstable" precursor.

## 4. Numbers the repository can use

Registry keys from `data/keys/compounds.yml`: `methanethiol`, `methional`, `dimethyl_disulfide`,
`dimethyl_trisulfide`, `hydrogen_sulfide` exist. Methionine, dimethyl sulfide, COS, CS2, S-methylmethionine:
not in registry.

| quantity | value | unit | conditions | source | evidence class | registry key |
|---|---|---|---|---|---|---|
| MeSH in fresh juice | ND, all 10 cultivars | ng/L | HS-SPME GC-PFPD, external curve | Table 2 | level_only | methanethiol |
| MeSH after heating | 279–917; mean 511 (10.6 nmol/L) | ng/L | 100 °C water bath, 60–90 s, sealed 20 mL vial, 5 mL juice | Table 2, text | level_only | methanethiol |
| H2S loss on heating | 3–36 % (mean of the ten rows 19.8 %, my arithmetic) | % | same | Table 2 | within_study_ratio | hydrogen_sulfide |
| DMS gain on heating | +87 % on the all-juice average (1785 → 3573) | ng/L | same | text | within_study_ratio | not in registry |
| DMDS, DMTS in 60–90 s juice | not observed (only after > 90 s, storage, commercial) | — | same | text | level_only | dimethyl_disulfide, dimethyl_trisulfide |
| MeSH odor threshold used | 0.02 | µg/L water | Guadagni 1963 | Table 1 | — | methanethiol |
| DMTS odor threshold used | 0.01 | µg/L water | Buttery 1976 | Table 1 | — | dimethyl_trisulfide |
| methional threshold used | 0.2 | µg/L water | Buttery 1971 | Table 1 | — | methional |
| Met in water → products | MeSH, DMDS, DMTS "limited"; methional none | qualitative | 1 mM Met, 100 °C, 60 min | Fig. 3A | figure_only | — |
| Met in sugar/citrate juice → products | methional preferred; MeSH/DMDS/DMTS reduced vs 3A | qualitative | 1 mM Met, MJ1, 100 °C, 60 min, pH unstated | Fig. 3B | figure_only | — |
| Met + 1.42 mM ascorbic acid → products | all four highest, at 15 min | qualitative | MJ2, 100 °C, 15 min | Fig. 3C | figure_only | — |
| MeSH sensory detection in juice | 10/18 triangle at 680 ng/L over DMS OAV 13.3 | — | 40 µL cups, orbital 10 min | text | level_only | methanethiol |

Nothing here is a rate. The methional → MeSH step is asserted (β-elimination, ref 27) but not measured;
the ascorbate acceleration is on the whole chain from methionine.

## 5. Flags

1. **The kinetic experiment exists but its numbers were never printed.** MJ1 and the control were
   sampled at 0, 1, 5, 10, 15, 30, 60 min and MJ2 at 0–15 min, with an eight-point calibration for
   MeSH, DMTS and methional in the same matrix; only three chromatograms appear. The SI has calibration
   slopes only. If the authors' data can be obtained, this becomes a methional/MeSH/DMDS/DMTS time series
   at 100 °C with and without 1.42 mM ascorbate — worth a request.
2. **Fig. 3 panels are not comparable by eye:** different heating times (60, 60, 15 min), PFPD in
   square-root mode, 1 cm fibre / 15 min extraction. Do not read intensities off them.
3. **Model-juice pH is not stated.** Citric acid 52 mM + tripotassium citrate 16 mM is well below the
   pH 6 at which Schutte 1972 found the methional β-elimination to need a catalyst; the "MeSH in water
   alone" result (3A) is at the unbuffered pH of 1 mM methionine.
4. **DMDS was identified but not quantified**, in juice or in models; DMTS was calibrated but is
   reported only in Fig. 3.
5. **Ascorbate's two roles are confounded.** The paper attributes the effect to Strecker activity
   (ascorbic acid / dehydroascorbic acid as α-dicarbonyl). Chin & Lindsay 1994 show ascorbate + Fe(III) is
   a pro-oxidant for MeSH → DMDS/DMTS at 30 °C; Schutte 1972 (Table II) found dehydroascorbic acid does
   NOT catalyse the methional → MeSH elimination at pH 6. The 15-min MJ2 chromatogram cannot say which
   step ascorbate accelerates.
6. **Sealed 20 mL vial with 15 mL air headspace** at 100 °C: the oxidant available is that headspace
   oxygen plus dissolved O2; not measured or varied.
7. **Static HS-SPME with external calibration in a sugar matrix**: levels depend on the juice's
   partitioning; the "% gain 100" entries mean "from ND", not a measured doubling.
8. Juice heating was 60–90 s: relevant to pasteurisation, not to a cook; the 8-h Sawamura values
   (H2S 300 µg/L, DMS 560 µg/L) quoted for comparison are 100–1000× higher.
