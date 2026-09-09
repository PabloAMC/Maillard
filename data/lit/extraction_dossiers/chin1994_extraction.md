# Chin & Lindsay 1994 — EXTRACTION (methanethiol oxidation to dimethyl disulfide and dimethyl trisulfide at 30 °C, pH 6.3: Cu(II), Fe(III), ascorbate, H2S, H2O2, aerobic vs anaerobic)
### The one paper on disk that varies the oxidant on methanethiol directly. One temperature, one pH, one table, four figures.

**Source on disk:** `data/articles/chin1994.pdf` (owner's download, 2026-09-08). Read from the `pdftotext`
text layer in the scratchpad; the layer is OCR-grade (scanned Food Chemistry 1994): "30 rain" = 30 min,
"10-6 µg CuCl2" = 10.6 µg (checked: 1 ppm Cu in 5 mL = 5 µg Cu = 10.6 µg CuCl2), "({ 15%)" = "(< 15 %)".
Table 1 re-typed. Figures 1–4 (all the time courses) are FIGURE-ONLY; only the text's statements about
them are carried. Fig. 5 (mechanism) transcribed as reactions.

## 0. Identity

| field | value |
|---|---|
| Title | "Ascorbate and transition-metal mediation of methanethiol oxidation to dimethyl disulfide and dimethyl trisulfide" |
| Authors | Hsi-Wen Chin, Robert C. Lindsay (Food Science, Univ. Wisconsin–Madison) |
| Venue | Food Chemistry 49 (1994) 387–392; received 19 March 1993, accepted 19 August 1993 |
| DOI | none printed (Elsevier: 10.1016/0308-8146(94)90009-4 — not verified against the PDF) |
| System | MeSH 1–2 ppm in 50 mM sodium phosphate pH 6.3, 30 °C, sealed serum vials; Cu(II), Fe(III), ascorbate, H2S, H2O2, phytate, benzoate |
| Analytes | MeSH, DMDS, DMTS by headspace GC-FPD |
| Companions | Chin & Lindsay 1993a (J. Food Sci. 58, 835), 1993b (cabbage sulfoxide lyase mechanisms) |

## 1. Why it matters

B17 §6 found the sulfur lane's disulfide channel is oxidant-limited: with the ambient oxidant pool
spent, the model holds 10–200× too little disulfide. Programme 6 will add methional → MeSH → DMDS/DMTS
onto that same pool. This paper is the direct measurement of what oxidises methanethiol in water:
Cu(II) at 1 ppm removes 70 % of 2 ppm MeSH in 30 min with air and 30 % without air; Fe(III) at 4 ppm
removes < 15 % in 300 min; the metal-free control is slow. DMDS carries only 43–77 % of the MeSH lost;
the rest goes to non-volatile species (sulfenic/sulfinic/sulfonic acids, Cu(I) mercaptide). DMTS is
NOT formed from MeSH plus metal alone — it needs hydrogen sulfide — and ascorbate + Fe(III) (via
H2O2 / hydroxyl radical, the authors argue) turns MeSH + H2S into DMTS about seven times faster than
Fe(III) alone. So: order of magnitude of the metal-catalysed rate at 30 °C; the split of MeSH loss
between disulfide and non-volatile sinks; the H2S requirement for the trisulfide; and the pro-oxidant
role of ascorbate. No temperature dependence anywhere.

## 2. Methods as they matter to a model

- **Vessels:** 120 mL serum vials with Mininert valves (gas-tight, re-sealable for headspace sampling).
- **Experiment A — metals (Fig. 1, Table 1).** 5 mL sodium phosphate **50 mM, pH 6.3**, demetallised over
  Chelex 100; purged 15 min at 30 mL/min with purified air (aerobic) or nitrogen (anaerobic), sealed
  under gas flow; **30 °C** water bath. Metal from fresh stock: **FeCl3 58 µg = 4 ppm Fe(III)** (20 µg Fe
  in 5 mL = **71.6 µmol/L**) or **CuCl2 10.6 µg = 1 ppm Cu(II)** (5 µg Cu = **15.7 µmol/L**), added just
  before the thiol. Start: **MeSH 10 µg = 2 ppm** (0.208 µmol; **41.6 µmol/L**) from a chilled
  N2-saturated methanolic stock via 10 µL gas-tight syringe. Headspace ≈ 115 mL: an air-purged vial holds
  roughly 1 mmol O2 against 0.2 µmol MeSH, so the aerobic arms are oxygen-saturated by a factor of
  ~5000 (my arithmetic). Residual O2 in the N2-purged arms is not measured.
- **Sampling A:** 0.3 mL headspace periodically → Varian 3700 GC, flame photometric detector, 6 ft × 2 mm
  glass column, 40/60 Carbopack B HT100, 40 °C isothermal for MeSH. DMDS and DMTS: 2 mL headspace
  **after 300 min**, oven 40 °C 1 min → 180 °C at 20 °C/min, hold 5 min. He 30 mL/min; injector and
  detector 200 °C.
- **Experiment B — ascorbate (Figs 2–4).** **50 mL** phosphate 50 mM pH 6.3 in the 120 mL vial, purged
  with purified air 20 min, sealed; **ascorbate 450 ppm** (2.56 mmol/L, MW 176.12), **Fe(III) 4 ppm**
  (71.6 µmol/L), **H2S 1 ppm** (29.3 µmol/L; from acidified Na2S), equilibrated at 30 °C; start with
  **MeSH 1 ppm** (20.8 µmol/L). 4 mL headspace sampled at times for DMDS and DMTS. Variants: ascorbate
  alone (450 ppm); ascorbate + phytate 1 mM; **H2O2 100 ppm** (2.94 mmol/L) with MeSH + H2S, no
  ascorbate/Fe; ascorbate + Fe(III) + benzoate 10 mM; Cu(II) in place of Fe(III) (data not shown);
  Table 1 row 4: Cu(II) 1 ppm + **ascorbate 10 ppm** (56.8 µmol/L), anaerobic, MeSH 2 ppm.
- **Quantification:** standard curves of authentic MeSH, DMDS, DMTS in the same buffer; **square root
  of peak area vs total mass (µg) in the vial** is linear (FPD sulfur response). Duplicate analyses, CV
  < 5 %. So "ppm" is the total mass in the vial divided by the liquid volume, with the headspace/liquid
  partition folded into the calibration; the reported quantity is the whole-vial inventory, not the
  dissolved concentration.
- **Conversions used below:** MeSH MW 48.11; DMDS 94.20; DMTS 126.26. 1 ppm MeSH = 20.8 µmol/L; 1 ppm
  DMDS = 10.6 µmol/L = 21.2 µmol/L MeSH-equivalents. Because MW(DMDS) ≈ 2 × MW(MeSH) (94.2 vs 96.2),
  the paper's mass-% conversion equals the sulfur-% conversion within 2 %.

## 3. Tables re-typed

### Table 1. "Fe(III) and Cu(II)-catalyzed conversion of methanethiol (2 ppm) to dimethyl disulfide after 300 min at 30 °C"

| sample | loss of CH3SH (ppm) | production of CH3SSCH3 (ppm) | % conversion (production / loss × 100) |
|---|---:|---:|---:|
| Fe(III), aerobic | 0.22 | 0.17 | 77 |
| Cu(II), aerobic | 1.64 | 0.84 | 51 |
| Cu(II), anaerobic | 1.06 | 0.46 | 43 |
| Cu(II) + ascorbate (10 ppm), anaerobic | 0.87 | 0.20 | 23 |

Molar re-statement (mine): Fe aerobic lost 4.6 µM, made 1.8 µM DMDS (3.6 µM MeSH-eq, 79 %); Cu aerobic
lost 34.1 µM, made 8.9 µM DMDS (17.8 µM MeSH-eq, 52 %); Cu anaerobic lost 22.0 µM, made 4.9 µM DMDS
(9.8 µM MeSH-eq, 44 %); Cu + ascorbate anaerobic lost 18.1 µM, made 2.1 µM DMDS (4.2 µM MeSH-eq, 23 %).
In the anaerobic Cu arm the unaccounted 12 µM MeSH is close to the 15.7 µM Cu charged — the authors'
Cu(I) mercaptide reading (CH3S–Cu) is numerically consistent. No control-row entry in Table 1.

### Statements about the figures (FIGURE-ONLY curves; text values only)

| figure | system | statement |
|---|---|---|
| Fig. 1 | MeSH 2 ppm, 30 °C, pH 6.3 | Cu(II) 1 ppm: **−70 % in 30 min aerobic, −30 % in 30 min anaerobic**. Fe(III) 4 ppm: **< 15 % removed over 300 min**, aerobic or anaerobic. Controls: "rates of disappearance ... were low", air > nitrogen. |
| Fig. 1 / text | metal experiment, 300 min | **only DMDS** detected as a volatile product; **DMTS not detected** at any time (no H2S present). |
| text | Fe(III) aerobic ± Fe | equal DMDS "after complete disappearance" of MeSH (data not shown) — Fe changes the rate, not the DMDS yield. |
| Fig. 2 | MeSH 1 ppm + H2S 1 ppm, aerobic, 30 °C, 330 min | **DMTS with ascorbate 450 ppm + Fe(III) 4 ppm ≈ 7 × DMTS with Fe(III) alone**; DMTS the predominant product; DMDS "very low and similar" across systems. Cu(II) gave similar DMTS enhancement (data not shown). |
| Fig. 3 | same, 330 min | ascorbate alone: **DMTS −45 %** vs ascorbate + Fe(III); ascorbate + phytate 1 mM: **DMTS −75 %** vs ascorbate alone (catalytically active metal required; Chelex left traces). |
| Fig. 4 | MeSH 1 ppm + H2S 1 ppm + H2O2 100 ppm, aerobic | DMTS formed; **DMDS > DMTS** with H2O2, whereas DMTS > DMDS with ascorbate + Fe(III). |
| text | ascorbate + Fe(III) + benzoate 10 mM, 330 min | **DMDS −42 %, DMTS −50 %** (hydroxyl-radical scavenger; data not shown). |

Fig. 2–4 y-axes are labelled 0–80 (units not legible in the text layer; the captions say "rates of
formation" in ppb-type units) — not read.

### Fig. 5 — proposed mechanism (AH2 = ascorbate, A = dehydroascorbate)

(a) CH3SH + Fe(III) → CH3S• + Fe(II) + H+; (b) AH2 + Fe(III) → AH• + Fe(II) + H+; (c) AH2 + O2 →(Fe(III))
A + H2O2; (d) Fe(II) + H2O2 → Fe(III) + •OH + OH− (Fenton); (e) CH3S• + •OH → CH3SOH (methanesulfenic
acid); (f) 2 CH3SOH + H2S → CH3SSSCH3 + 2 H2O; (g) 2 CH3SOH → CH3S(O)SCH3 + H2O (methyl
methanethiosulfinate); (h) 2 CH3S(O)SCH3 + H2S → CH3SSSCH3 + 2 CH3SOH.

## 4. Numbers the repository can use

Registry keys: `methanethiol`, `dimethyl_disulfide`, `dimethyl_trisulfide`, `hydrogen_sulfide` exist.
Cu(II), Fe(III), ascorbate: not compound-registry items (process levers).

| quantity | value | unit | conditions | source | evidence class | registry key |
|---|---|---|---|---|---|---|
| MeSH loss, Cu(II) 1 ppm, aerobic | 70 % at 30 min; 82 % (1.64/2.00) at 300 min | % of 41.6 µM | 50 mM phosphate pH 6.3, 30 °C, air-saturated sealed vial | Fig. 1 text, Table 1 | within_study_ratio (two time points) | methanethiol |
| apparent first-order k, Cu aerobic, first 30 min | ≈ 0.040 (t½ ≈ 17 min) | min⁻¹ | same; **my two-point derivation**; the 300-min point (82 %) shows the loss is far from first-order with constant k — a fast Cu-limited phase then slow | derived | within_study_ratio | methanethiol |
| MeSH loss, Cu(II) 1 ppm, anaerobic (N2) | 30 % at 30 min; 53 % at 300 min | % | same, N2-purged | Fig. 1 text, Table 1 | within_study_ratio | methanethiol |
| MeSH loss, Fe(III) 4 ppm, aerobic or anaerobic | < 15 % at 300 min (Table 1: 11 %) | % | same | text, Table 1 | within_study_ratio | methanethiol |
| apparent first-order k, Fe(III) aerobic | ≈ 4 × 10⁻⁴ (t½ ≈ 30 h) | min⁻¹ | same; my derivation from 11 % / 300 min | derived | within_study_ratio | methanethiol |
| MeSH loss, no metal (Chelexed buffer) | "low"; air > N2 | qualitative | same | Fig. 1 text | figure_only | methanethiol |
| DMDS yield on MeSH lost | 77 / 51 / 43 / 23 % (Fe aer / Cu aer / Cu anaer / Cu + asc anaer) | % (mass ≈ sulfur) | 300 min, 30 °C | Table 1 | within_study_ratio | dimethyl_disulfide |
| DMDS at 300 min | 0.17 / 0.84 / 0.46 / 0.20 | ppm (1.8 / 8.9 / 4.9 / 2.1 µM) | same rows | Table 1 | level_only | dimethyl_disulfide |
| DMTS without H2S | not detected through 300 min | — | MeSH + Cu or Fe, air or N2 | text | level_only | dimethyl_trisulfide |
| DMTS enhancement by ascorbate 450 ppm + Fe(III) 4 ppm | ≈ 7 × Fe(III) alone | ratio | MeSH 1 ppm + H2S 1 ppm, air, 30 °C, 330 min | Fig. 2 text | within_study_ratio | dimethyl_trisulfide |
| DMTS, ascorbate alone vs ascorbate + Fe | 0.55 × | ratio | same, 330 min | Fig. 3 text | within_study_ratio | dimethyl_trisulfide |
| DMTS, ascorbate + phytate 1 mM vs ascorbate alone | 0.25 × | ratio | same | Fig. 3 text | within_study_ratio | dimethyl_trisulfide |
| benzoate 10 mM on ascorbate + Fe system | DMDS × 0.58, DMTS × 0.50 | ratio | same, 330 min | text (data not shown) | within_study_ratio | both |
| H2O2 100 ppm (2.9 mM) as oxidant | DMTS formed; DMDS > DMTS | qualitative | MeSH 1 ppm + H2S 1 ppm, air, 30 °C | Fig. 4 | figure_only | both |
| MeSH odor threshold quoted | 0.02 ppb water | — | Hansen 1992 | intro | — | methanethiol |

No measured barrier (single temperature), no explicit order in thiol or oxygen (single MeSH level per
experiment; O2 only as air vs N2).

## 5. Flags

1. **30 °C only, pH 6.3 only, 50 mM phosphate.** Nothing on temperature; any use at cook temperature is
   an extrapolation the paper does not support.
2. **Every time course is a figure**; the numbers above are the text's two-point statements and one
   300-min table. Fig. 1 would give the full MeSH decay curves for six arms if re-read from a clean scan
   (not attempted here — the scan's text layer is OCR-grade).
3. **The Cu(II) arm is not a clean catalytic rate.** 1 ppm Cu (15.7 µM) against 41.6 µM MeSH: aerobic
   loss is 2.2 MeSH per Cu (turnover with O2), anaerobic 1.4 per Cu with ~12 µM unaccounted ≈ the Cu
   charge (mercaptide). Fitting a single first-order constant to it would be wrong in kind.
4. **DMDS is only half the MeSH sink even with air and Cu.** 23–57 % of lost MeSH is not recovered as
   any headspace volatile; the authors name CH3SOH / CH3SO2H / CH3SO3H and CH3S–Cu(I). A model that
   closes MeSH loss onto DMDS alone will overpredict DMDS about twofold under Cu, and by 1.3× under Fe.
5. **DMTS needs a second sulfur.** With MeSH as the only sulfur source there was no DMTS in 300 min at
   all; DMTS appeared only when H2S (1 ppm) was present, and then mostly under ascorbate + metal or
   H2O2. For the repo: the trisulfide should draw on the H2S pool, not on MeSH alone (Yu 1995 and Zhang
   2023 point the same way).
6. **Ascorbate is a pro-oxidant here, via H2O2/•OH with a metal**, not an antioxidant; its 450 ppm
   (2.56 mM) is a cabbage level. Cheng 2020 saw ascorbate accelerate the whole methionine chain at
   100 °C; this paper says part of that can be oxidation, not just Strecker.
7. **Anaerobic ≠ oxygen-free.** N2 purge 15 min at 30 mL/min into a 120 mL vial, then Chelexed buffer;
   residual O2 not measured. The 30 % anaerobic loss with Cu is consistent with Cu(II) itself as the
   electron acceptor, but a trace-O2 contribution cannot be excluded.
8. **"ppm" is a whole-vial inventory** (5 mL liquid, 115 mL headspace; or 50 mL / 70 mL), calibrated
   with standards in the same geometry. Dissolved concentrations are lower than the ppm values by the
   headspace partition, which for MeSH at 30 °C is large.
9. Duplicates with CV < 5 % on the analytical side; no statement of independent replicate vials.
