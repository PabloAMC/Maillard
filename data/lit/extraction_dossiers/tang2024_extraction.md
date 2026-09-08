# Tang, Teodorowicz, Boeren, Wichers & Hettinga 2024 — EXTRACTION (soy and pea protein + glucose: lysine, furosine, CML and CEL after wet heating 85 °C / 30 min and dry heating 60 °C / aw 0.6 / 48 h)
### A bioactivity paper whose §3.1 carries the glycation-marker end points Programme 7 wants. The absolute furosine, CML and CEL values live in a bar chart; only the lysine levels and the changes quoted in the text are numbers.

**Source on disk:** `data/articles/Tang2024.pdf` (owner's download, 2026-09-08). Read from the `pdftotext`
text layer in the scratchpad; clean. Fig. 2 (the four glycation markers, 2 proteins × 3 treatments) is a
raster image; page 4 was re-extracted with `pypdf` in layout mode to confirm that no bar values are
embedded — FIGURE-ONLY. The LC-MS/MS method for the markers is given only by reference to Tang, Wichers
& Hettinga 2023 (Food Hydrocolloids 136, 108251), not on disk. Figures 3–6 and Table 1 (peptidomics,
sRAGE, antimicrobial) are outside the repository's question and are summarised in one line each.

## 0. Identity

| field | value |
|---|---|
| Title | "sRAGE-binding and antimicrobial bioactivities of soy and pea protein after heating and in vitro infant digestion" |
| Authors | Jiaying Tang, Malgorzata Teodorowicz, Sjef Boeren, Harry J. Wichers, Kasper A. Hettinga (Wageningen) |
| Venue | Food Research International 183 (2024) 114224; received 4 Oct 2023, accepted 11 Mar 2024; CC BY |
| DOI | 10.1016/j.foodres.2024.114224 |
| Systems | soy protein (SP, 77.9 % protein) and pea protein (PP, 78.8 %) from fresh seed, + glucose 4-fold (w/w); NT / wet-heated (W) / dry-heated (D); then in vitro infant digestion |
| Markers | lysine, furosine, CML, CEL by LC-MS/MS; independent duplicates |

## 1. Why it matters

Programme 7 (roadmap §5d) makes the isolate's protein-bound lysine a slow Maillard reactant with CML
and CEL as the markers. This paper gives the starting lysine per gram of protein for a soy and a pea
protein (84.2 and 101.7 mg/g protein), and the end-point change in lysine, furosine, CML and CEL after
two heat treatments with glucose in large excess. Dry heating at 60 °C / aw 0.6 / 48 h takes 47–49 % of
the lysine and adds ~46 mg furosine per g protein; wet heating at 85 °C / 30 min at 1.2 % protein takes
8–22 % of the lysine and adds 1.4–1.9 mg/g furosine. These are single end points, not series, and the
CML/CEL absolutes are figure-only, but the lysine numbers are exactly what the matrix table needs to
charge lysine sites, and the wet/dry contrast is a within-study ratio on the early-stage marker.

## 2. Methods as they matter to a model

- **Raw material:** fresh soybean (Glycine max, 39.3 % protein w/w) and fresh pea (Pisum sativum, 23.3 %
  protein w/w; retail, brand HAK), protein by Dumas with **N × 5.7**. Protein preparations per Tang 2022/
  2023: **SP 77.9 % protein, PP 78.8 % protein (w/w)**. Not commercial isolates.
- **Sugar:** "proteins were mixed with 4-fold glucose (w/w) to mimic plant-based IF" — read as glucose :
  protein = 4 : 1 by mass. ⚠ Whether "protein" here is the protein mass or the preparation mass is not
  stated; the factor 0.78 between them is the ambiguity.
- **Wet heating (W):** protein 1.2 % (w/v) (glucose therefore 4.8 % w/v ≈ 266 mmol/L if 4× the protein
  mass), water bath **85 °C, 30 min**. Solvent: 10 mM PBS is named only for the dry-heated
  reconstitution; the wet and NT solvents are presumably the same PBS (not stated). Lysine sites at
  1.2 % protein (12 g/L): SP 6.9 mmol/L, PP 8.4 mmol/L (from §3 lysine, my conversion) → glucose :
  lysine ≈ 32–39 : 1.
- **Dry heating (D):** the protein + glucose powder in a desiccator at **60 °C, water activity 0.6, 48 h**,
  then reconstituted in 10 mM PBS to 1.2 % protein (w/v).
- **Non-treated (NT):** 1.2 % protein (w/v), unheated.
- **Replicates:** "independent duplicates" for all treatments; Fig. 2 error bars = SD of duplicates;
  ANOVA + Duncan.
- **Markers:** lysine, furosine (Nε-2-furoylmethyllysine), CML, CEL by **LC-MS/MS "according to our
  previously established method (Tang et al., 2023)"**; hydrolysis, internal standards and calibration are
  not restated. Units: **mg/g protein for lysine and furosine; mg/100 g protein for CML and CEL.**
- **Time structure:** one end point per treatment (NT, W, D). **No time series; nothing here is a
  rate.**
- **Downstream (not extracted):** infant in vitro digestion (gastric pH 5.3, 60 min, pepsin 268 U/mL;
  intestinal pH 6.6, 60 min, trypsin 16 U/mL in pancreatin), peptidomics with +162 (hexose), +58 (CML),
  +72 (CEL) at K and +54 (MG-H), +40 (G-H) at R; sRAGE inhibition ELISA; growth of E. cloacae and
  S. epidermidis; CAMPR4 AMP prediction.
- **Conversions used below:** lysine MW 146.19 → 1 mg/g = 6.84 µmol/g; furosine C12H18N2O4 MW 254.28
  → 1 mg/g = 3.93 µmol/g; CML 204.22; CEL 218.25. Per gram of preparation: multiply per-protein values by
  0.779 (SP) or 0.788 (PP).

## 3. Tables re-typed

There is no printed table of marker values. Everything quoted in §3.1 is collected here; blank cells
are FIGURE-ONLY (Fig. 2).

### §3.1 glycation markers as stated in text

| marker (unit) | protein | NT | W (85 °C / 30 min, 1.2 % w/v) | D (60 °C / aw 0.6 / 48 h) |
|---|---|---:|---:|---:|
| lysine (mg/g protein) | SP | **84.2** | "decreased by ~22 %" → ≈ 65.7 | "decreased by ~47 %" → ≈ 44.6 |
| lysine (mg/g protein) | PP | **101.7** | "~8 %" → ≈ 93.6 | "~49 %" → ≈ 51.9 |
| furosine (mg/g protein) | SP | figure-only | NT + 1.4 | NT + ~46 |
| furosine (mg/g protein) | PP | figure-only | NT + 1.9 | NT + ~46 |
| CML (mg/100 g protein) | SP | figure-only | "no increase" | ~12 × NT |
| CML (mg/100 g protein) | PP | figure-only | ~2 × NT | ~48 × NT |
| CEL (mg/100 g protein) | SP | figure-only | ~2 × NT ("around 2-fold") | ~9 × NT |
| CEL (mg/100 g protein) | PP | figure-only | ~2 × NT | ~16 × NT |

Text also: "despite PP having a higher initial lysine content than SP, after dry heating, their lysine
levels became comparable" (≈ 44.6 vs 51.9 by the percentages); "dry heating contributed more to glycation
than wet heating, and glycation of PP was more extensive than SP".

### Derived (my arithmetic)

| quantity | SP | PP |
|---|---:|---:|
| lysine, NT, µmol/g protein | 576 | 696 |
| lysine, NT, mg/g preparation (× 0.779 / 0.788) | 65.6 | 80.1 |
| lysine lost, W, µmol/g protein | ≈ 127 | ≈ 56 |
| lysine lost, D, µmol/g protein | ≈ 271 | ≈ 341 |
| furosine gained, W, µmol/g protein | 5.5 | 7.5 |
| furosine gained, D, µmol/g protein | ≈ 181 | ≈ 181 |
| furosine gained / lysine lost, W | 0.043 | 0.13 |
| furosine gained / lysine lost, D | 0.67 | 0.53 |

The literature convention furosine × ~3.1 ≈ Amadori-lysine (acid hydrolysis converts ~32 % of
fructoselysine to furosine) is NOT applied by the paper and is not applied here; note that applying it
to the D rows would give 0.56 mmol/g Amadori-lysine, more than the lysine lost, which says either the
"lysine" measured is not simply unmodified lysine after hydrolysis (Amadori partly reverts to lysine in
acid) or the factor does not hold for this method — see Flag 4.

### Other results (one line each; not for the model)

Glycated peptides after digestion: relative intensity ~28 % (SP-D) and 29 % (PP-D) of total at I60;
relative number ~10 % and 20 %; > 80 % of glycated peptides K-modified, hexose (+162) dominant (~80 %
SP-D, ~95 % PP-D); MG-H and G-H at R in SP, MG-H only in PP. sRAGE inhibition before digestion > 50 %
for all, SP-D ~50 %, PP-D ~53 %; after digestion < 18 % (SP −3 to 4 %, PP 6–18 %). E. cloacae growth
rate up 28–59 % on D digests; S. epidermidis suppressed by PP digests regardless of heating. 519 (SP) and
133 (PP) predicted AMPs, 9 and 6 glycated (Table 1).

## 4. Numbers the repository can use

Registry keys: `reactive_lysine` (the closest key for protein lysine), `furosine`, `cml`, `cel`; marker sets
`furosine_cml`, `cml_cel`. Free lysine as a molecule: not in registry.

| quantity | value | unit | conditions | source | evidence class | registry key |
|---|---|---|---|---|---|---|
| lysine, soy protein, unheated | 84.2 (576 µmol/g protein; 65.6 mg/g preparation) | mg/g protein | SP 77.9 % protein from fresh soybean; LC-MS/MS | §3.1 | level_only | reactive_lysine |
| lysine, pea protein, unheated | 101.7 (696 µmol/g protein; 80.1 mg/g preparation) | mg/g protein | PP 78.8 % protein from fresh pea | §3.1 | level_only | reactive_lysine |
| lysine loss, wet | SP ~22 %, PP ~8 % | % of NT | + glucose 4:1 w/w, 1.2 % protein, 85 °C, 30 min | §3.1 | within_study_ratio (single end point) | reactive_lysine |
| lysine loss, dry | SP ~47 %, PP ~49 % | % of NT | powder, 60 °C, aw 0.6, 48 h | §3.1 | within_study_ratio (single end point) | reactive_lysine |
| furosine gain, wet | SP +1.4, PP +1.9 | mg/g protein | 85 °C, 30 min | §3.1 | level_only (difference; NT absolute figure-only) | furosine |
| furosine gain, dry | ~ +46 (both) | mg/g protein | 60 °C, aw 0.6, 48 h | §3.1 | level_only (difference) | furosine |
| CML, wet / NT | SP 1.0 ×; PP ~2 × | ratio | 85 °C, 30 min | §3.1 | within_study_ratio | cml |
| CML, dry / NT | SP ~12 ×; PP ~48 × | ratio | 60 °C, aw 0.6, 48 h | §3.1 | within_study_ratio | cml |
| CEL, wet / NT | ~2 × (both) | ratio | 85 °C, 30 min | §3.1 | within_study_ratio | cel |
| CEL, dry / NT | SP ~9 ×; PP ~16 × | ratio | 60 °C, aw 0.6, 48 h | §3.1 | within_study_ratio | cel |
| CML, CEL, furosine absolutes (NT, W, D) | bar chart | mg/100 g protein; mg/g protein | — | Fig. 2 | figure_only | cml, cel, furosine |
| average furosine formation rate, dry | ≈ 0.96 mg/g protein/h (≈ 3.8 µmol/g/h) over 48 h | mg/g/h | 60 °C, aw 0.6; **two-point average, my derivation, not a measured rate** | derived | within_study_ratio | furosine |
| average lysine loss rate, wet | SP ≈ 0.62, PP ≈ 0.27 mg/g protein/min over 30 min | mg/g/min | 85 °C, 1.2 % w/v, glucose 4:1; two-point average | derived | within_study_ratio | reactive_lysine |

None of these is a measured rate; each treatment is one end point in duplicate.

## 5. Flags

1. **FIGURE-ONLY absolutes.** The NT furosine, CML and CEL levels — the baselines the fold changes
   multiply — are only in Fig. 2. Without them "~48 × NT" for pea CML is not a level. If the authors'
   Tang 2023 paper (Food Hydrocolloids) tabulates the same materials, it is the source to retrieve.
2. **Single end points, duplicates.** Two heating regimes at different temperatures, water contents and
   times (85 °C / 30 min / dilute solution vs 60 °C / aw 0.6 / 48 h / powder): they cannot be combined
   into a rate law or a barrier; the per-hour numbers in §4 are averages over the whole treatment.
3. **The glucose ratio is large** (≈ 32–39 glucose per lysine site in the wet system, 266 mM glucose):
   the lysine loss is pseudo-first-order in lysine at near-saturating sugar, an infant-formula ratio
   rather than a savoury-recipe one. A pea-isolate recipe with a few percent sugar sits below this.
4. **What "lysine" is analytically is not restated.** LC-MS/MS after hydrolysis (method by reference).
   If acid hydrolysis was used, part of the Amadori product reverts to lysine and part becomes furosine,
   so "lysine" overstates unmodified lysine and the furosine/lysine-loss ratios in §3 are method-bound.
5. **The dry regime is a storage-type glycation, not a cook.** 60 °C for 48 h at aw 0.6 is the classic
   dry-glycation protocol; the CML/CEL multipliers it produces (9–48 ×) say nothing about a 100–180 °C
   cook of an isolate at high water. For Programme 7 they bound the marker response of these two proteins,
   not its rate at cooking temperature.
6. **PP is more reactive than SP** at every marker after dry heating (lysine −49 vs −47 %, CML 48 vs 12 ×,
   CEL 16 vs 9 ×) and less after wet heating (lysine −8 vs −22 %) — the ordering flips with the regime,
   so a single "pea vs soy reactivity" factor would be wrong.
7. **"~46 mg/g protein" is one rounded figure given for both proteins**; the paper does not say they
   were equal to that precision.
8. Wet-heating solvent and the exact basis of "4-fold glucose" are unstated (§2).
9. The materials are laboratory preparations from fresh seed (77.9 / 78.8 % protein), not commercial
   isolates; commercial soy/pea isolates run 80–90 % protein and have been through their own heat
   history, so their NT furosine and CML would differ.
