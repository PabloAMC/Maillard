# Zhang 2024b — EXTRACTION (methionine + thiamine + xylose 15 + 15 + 15 mg/mL, phosphate pH 4.9, 115 C, up to 60 min, ± cysteine / glutathione / cystine 15 mg/mL; MFT, MMFT, MFT-MFT, MeSH, DMDS, DMTS; dose series 0.5-75 mg/mL Cys or GSH; text-layer re-extraction under the no-figures rule)
### Identity correction first: this is the SAME paper as the on-disk `Zhang2024_extraction.md` (wave Z3, 2026-08-28), not a companion to it. That dossier's numbers are 300-dpi figure read-offs; this one keeps only what is printed.

**Source on disk:** `data/articles/Zhang2024.pdf` (10 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Zhang2024.txt`); Table 1 (sensory difference matrix) came through
clean. Figures 1, 2, 3, 4, 5 were not read and no value is taken from them. The Supporting
Information (Fig. S1 dynamic profiles, S2 single-factor sweeps, S3 response surfaces, S4
metabolomics; Tables S1-S3 design, S4 SIM ions, S5 calibration curves, **S6 zero-order fits, S7
linear regressions**) is NOT on disk.

**On the brief's premise.** The brief describes `Zhang2024_extraction.md` as "the same group's JAFC
2024 paper". It is not: its header reads "Zhang, Kuang, Wang & Cao 2024 (10.1016/j.foodres.2024.114149)
— Wave Z3 extraction 2026-08-28", i.e. this Food Research International paper. No JAFC 2024 paper by
Zhang, Wang and Cao is on disk or in the repository's registries; the only "Zhang JAFC 2024" there
(10.1021/acs.jafc.4c05736) is Zhang, Cui, Xia et al. (Jiangnan; DiXyl-Lys-ARP degradation), a
different group and topic. The BTBU group's JAFC paper is 2023 (71, 13043; `zhang2023_extraction.md`,
the storage study). "What is new against it" below therefore means: what this re-extraction changes
relative to the wave-Z3 dossier of the same paper.

## 0. Identity

| field | value |
|---|---|
| Title | "Metabolomics reveals factors affecting the radical reaction of sulfides during thermal processing for meaty aroma" |
| Authors | Zeyu Zhang, Huiyu Kuang, Bei Wang*, Yanping Cao* (Beijing Technology and Business University) |
| Venue | Food Research International 182 (2024) 114149. Received 30 November 2023, revised 15 February 2024, accepted 17 February 2024, online 20 February 2024 |
| DOI | 10.1016/j.foodres.2024.114149 |
| Naming | MFT = 2-methyl-3-furanthiol; MMFT = 2-methyl-3-(methyldisulfanyl)furan (the MFT-MeSH mixed disulfide); MFT-MFT = bis(2-methyl-3-furyl) disulfide; MFFT = methyl furfuryl disulfide (named, not measured); MeSH, DMDS, DMTS; VB1 = thiamine; Xyl = xylose; Met; Cys, GSH, GCys = cystine; IPDS = isopropyl disulfide (internal standard); "second-order thermal procedure" = a two-stage temperature programme, not a kinetic order; "simple model" = VB1 + Xyl; "multi-component model" = Met + VB1 + Xyl; R1 = DMTS concentration, R2 = MMFT / MFT-MFT ratio (the two RSM responses) |
| Companions | `Zhang2024_extraction.md` (same paper, figure-based); Zhang, Wang & Cao 2023 JAFC 71, 13043 (`zhang2023_extraction.md`; the isotope evidence that MeSH + MFT -> MMFT competes with DMTS, cited here); Zhang 2023 FRI 172, 113200 (storage; cited by Pan 2025); Chin & Lindsay 1994 (`chin1994_extraction.md`, the oxidant study the B19 log leans on) |

## 1. Why it matters

Two of the repository's open questions run through this paper. (i) B17's finding
(`results/validation/kinetic_core_b17_prereg.md` section 6, T3): the model's MFT-dimer share is ten to
two hundred times below what "Zhang 2024's cysteine arm (measured 8.7 %)" shows, and the diagnosis is
that the thiol-disulfide channel is oxidant-limited. That 8.7 % came from the wave-Z3 dossier's
read-off of Fig. 1a (0.115 / 1.34 ng/mL). Under the current rule it is figure-only, and this dossier
says what the printed text supports instead. (ii) The methionine chain: the paper's system has Met
as the MeSH source and thiamine + xylose as the MFT source, so MeSH is partitioned between DMDS /
DMTS (self-oxidation) and MMFT (capture by MFT) — the one measured competition between the two sinks
of methanethiol in the corpus. The printed kinetic content is thin (five numbers without units for
MMFT, and their ratios to a control), but the ratios are unit-free and usable, and the paper's
statements about the oxidant are what B17 asked for.

## 2. Methods as they matter to a model

- **Simple model (2.2).** "The binary mixture of VB1 at 15.0 mg/mL (w/v) and Xyl in a ratio 1:1 was
  solubilized in pH 4.9 ± 0.1 phosphate buffered solution. The mixture was made to react for 60 min,
  then immediately cooled in ice water." Additives Cys, GSH, GCys "of equal weight as VB1 and Xyl",
  i.e. 15 mg/mL each. **Temperature not stated in 2.2** (115 C by analogy with 2.3). Results: Fig. 1
  (figure-only).
- **Multi-component model (2.3).** "Met (15.0 mg/mL, w/v) as the source of MeSH mixed with VB1 and
  Xyl as a ratio of 1:1:1, which dissolved in phosphate buffered solution (pH 4.9 ± 0.1). For
  treatment groups, Cys (15.0 mg/mL), GSH (15.0 mg/mL), and GCys (15.0 mg/mL) were mixed into pH 4.9
  phosphate buffered solution, respectively. The mixture solutions were incubated for up to 60 min
  at 115 C." Buffer concentration, vessel, fill, headspace and atmosphere are not stated. The
  "dynamic profiles" (Fig. S1) and the regression "before and after 90 min" imply runs longer than
  60 min whose duration is nowhere printed (Flags 4).
  **Molar conversions (mine):** Met 15 mg/mL = **100.5 mmol/L** (M 149.21); Xyl 15 mg/mL = **99.9
  mmol/L** (150.13); VB1 15 mg/mL = **44.5 mmol/L** if thiamine hydrochloride (337.27) or 45.8 if the
  mononitrate (327.36) — the salt is not stated; Cys 15 mg/mL = **123.8 mmol/L** (121.16); GSH 15
  mg/mL = **48.8 mmol/L** (307.32); GCys 15 mg/mL = **62.4 mmol/L** (240.30; 124.8 mmol/L of cysteine
  equivalents). The three sulfur additives are therefore not equimolar: on a thiol basis Cys 124,
  GSH 49, cystine 0 (oxidised) mmol/L.
- **Dose series (2.6) and RSM (2.7-2.8).** "The ternary mixture of Met, VB1 and Xyl (1:3:3, w/w) mixed
  with a series of Cys (or GSH) solutions ... (0.5, 2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 62.5 and
  75.0 mg/mL) ... incubated at 115 C for 60 min"; then Cys:GSH mixtures 10:0 ... 0:10. Section 2.5
  fixes the ternary's methionine: the Cys series "from 0.5 mg/mL to 100.0 mg/mL in models, as a
  ratio of 0.1, 0.5, 1, 2, ..., and 20 with Met concentration" -> **Met = 5 mg/mL = 33.5 mmol/L** in
  the 1:3:3 ternary, with VB1 and Xyl at 15 mg/mL; consistent with 3.3, where the Cys/Met ratio 6
  "corresponds to the optimal Cys used level of 30.0 mg/mL". So the Fig. 2 dose series has one third
  of the methionine of the 1:1:1 model. RSM: Box-Behnken, four factors (Cys:GSH ratio, initial
  reaction temperature, second reaction time, second-stage temperature), 27 runs, three centre
  replicates; responses R1 = DMTS, R2 = MMFT/MFT-MFT.
- **Thermal process flavourings (2.4).** Pork longissimus + ingredients, 115 C 60 min (control);
  treatment with Cys 2.5 mg/mL + GSH 6.0 mg/mL, 115 C 50 min then 95 C 10 min; storage 112 d.
  Real-food arm, semi-quantitative (see below).
- **Quantification (2.10).** 10 mL sample + 1 µL IPDS "0.0943 ng/mL" (as printed; 9.4e-5 ng per vial
  is not a usable internal-standard amount — mg/mL is the likely intent, Flags 5) in a 40-mL vial,
  45 C 30 min, DVB/CAR/PDMS 30 min at 45 C, splitless, DB-WAX 60 m x 0.25 mm x 0.25 µm, 40 C (2 min)
  -> 190 C at 3 C/min -> 230 C at 10 C/min, **SIM** (ions Table S4). Concentrations in ng/mL from
  **external standard curves of eight levels** (Table S5) as mass ratio vs peak-area ratio to the IS.
  Verification test (flavourings): semi-quantitative, area ratio to IS. Triplicates, randomised.
- **Statistics.** ANOVA, Duncan p < 0.05; RSM in Design-Expert 12; metabolomics OPLS-DA, VIP > 1,
  p < 0.05.

## 3. Tables re-typed and printed numbers

### Table 1. "The ratio of Cys to Met corresponding to differences between individual and group after Cys added into the model and thermal process flavorings"

Forced-choice ascending method, twelve panellists, "+" = difference detected. Cys/Met ratios 0.1,
0.5, 1, 2, 4, 6, 8, 10, 12.5, 15, 17.5, 20. Model (1:3:3 ternary, 115 C 60 min): first "+" at ratio 6
for panellists 2, 3, 5, 6, 7, 9, 11, 12; at 4 for 8 and 10; at 8 for 1 and 4. Thermal process
flavourings: first "+" at ratio 8 for 2, 3, 4, 6, 7, 9, 10, 12; at 6 for 5 and 8; at 10 for 1 and 11.
Authors' reading: the additive becomes sensorially detectable at Cys/Met = 6 in the model (30.0
mg/mL) and 8 in the flavouring; keep Cys/Met < 6. (Sensory; no chemistry.)

### The kinetic sentences (3.1, p. 4) — the only rate numbers in the paper

"The dynamic analysis results (Table S6) showed that the zero-order kinetic equations of MMFT levels
in the Cys and GSH treatment groups had better-fitted data. The two groups' respective reaction rate
constants were 0.0028 and 0.0031 ... linear regression (shown in Table S7) ... the two treatments and
control group showed interesting two-stage trends before and after 90 min of thermal processing ...
during the first 90 min of thermal processing, the slope of MMFT level in Cys and GSH groups reached
0.0043 and 0.0033, which are about 3 times and 2.5 times the slope in the control group,
respectively. After that, the slope in the Cys group dropped to a comparable level to that in the
control group. The slope of MMFT in the GSH treatment group is 0.0032, which is 4 times that of the
control group."

| quantity | Cys arm | GSH arm | control |
|---|---:|---:|---:|
| zero-order k of MMFT, whole run | 0.0028 | 0.0031 | not printed |
| slope of MMFT, first 90 min | 0.0043 (≈ 3 x control) | 0.0033 (≈ 2.5 x control) | ≈ 0.0014 / 0.0013 (mine, from the two ratios) |
| slope of MMFT, after 90 min | "comparable" to control | 0.0032 (4 x control) | ≈ 0.0008 (mine, 0.0032 / 4) |

**Units are not printed**; MMFT is reported in ng/mL and the paper's times are in minutes, so ng
mL-1 min-1 is the natural reading (1 ng/mL/min of MMFT, M 160.25, = 6.2e-3 µmol L-1 min-1), but
Table S6/S7 are not on disk to confirm, and the run length behind them is not stated. The ratios
(3, 2.5, 4, "comparable") are unit-free.

### Other printed numbers

- Dose series (3.2, describing Fig. 2, values figure-only): MeSH rises then falls with additive
  amount, falling above 40 mg/mL; DMDS "decreased gradually with the increase of Cys and GSH", no
  difference between the two; DMTS rises then falls, lower in GSH than Cys; MFT rises with
  additive; MFT-MFT "similar to the changes of MFT exhibiting an increasing trend with the addition
  of additives"; MMFT no difference between arms; the ratio MMFT / MFT-MFT "decreased" with
  additive and is "significantly larger" in the Cys arm.
- Single-factor sweeps (3.4.1, Fig. S2): DMTS lower and MMFT/MFT-MFT higher at lower initial
  temperature, at second-stage reaction time < 20 min and second-stage temperature < 105 C.
- RSM: R2 0.905 (DMTS) and 0.873 (ratio); p 0.0004 and 0.0019; lack of fit not significant.
  **Optimum: Cys:GSH 3:7, initial 115 C, second reaction time 10 min, second-stage 95 C.**
  Verification in the model: **MMFT / MFT-MFT measured 1.422 vs predicted 1.294**; DMTS measured
  within **7.0 %** of predicted (neither value printed).
- Flavourings: DMTS lower and MMFT/MFT-MFT higher in the treatment group (Fig. 4a, figure-only);
  DMTS lower over 112 d storage (Fig. 4c, semi-quantitative, figure-only); aroma scores Fig. 4b
  (figure-only).
- Metabolomics: 232 differential metabolites, 86 KEGG pathways; thiamine phosphate down-regulated,
  Cys up-regulated in the treatment; interpretation only.

## 4. Kinetic numbers the repository can use

Registry mapping: MFT -> `2_methyl_3_furanthiol`; MFT-MFT -> `bis_2_methyl_3_furyl_disulfide`
(alias "mft disulfide"); methanethiol -> `methanethiol`; DMDS -> `dimethyl_disulfide`; DMTS ->
`dimethyl_trisulfide`; MMFT (2-methyl-3-(methyldisulfanyl)furan), MFFT, methionine, thiamine (the
registry has the lever `thiamine_availability`, not the molecule), xylose, cysteine, glutathione,
cystine, isopropyl disulfide -> not in registry.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| MMFT formation slope, first 90 min, Cys / GSH arm relative to control | 3 / 2.5 | — | Met + VB1 + Xyl 100 + 45 + 100 mmol/L, phosphate pH 4.9, 115 C, + Cys 124 or GSH 49 mmol/L; run length not printed | zero order (linear in t) | text 3.1, Table S7 | within_study_ratio |
| MMFT formation slope after 90 min, GSH / control; Cys / control | 4; ≈ 1 | — | same | zero order | text 3.1 | within_study_ratio |
| MMFT zero-order k, Cys / GSH | 0.0028 / 0.0031 | not printed (ng mL-1 min-1 presumed) | same | zero order | text 3.1, Table S6 | measured_rate — DO NOT INGEST until the unit and run length are confirmed |
| MMFT slopes, first 90 min, Cys / GSH; after 90 min, GSH | 0.0043 / 0.0033; 0.0032 | not printed | same | linear | text 3.1, Table S7 | measured_rate (same caveat) |
| control MMFT slopes (mine) | ≈ 0.0014 (first 90 min), ≈ 0.0008 (after) | not printed | same | — | derived from the printed ratios | derived_assumption |
| MMFT / MFT-MFT at the RSM optimum (model) | 1.422 measured, 1.294 predicted | — | 1:3:3 ternary (Met 33.5 mmol/L), Cys:GSH 3:7, 115 C then 95 C / 10 min | — | text 3.4.2 | level_only (concentration ratio; response factors differ) |
| DMTS at the optimum vs prediction | within 7.0 % | — | same | — | text 3.4.2 | level_only (values not printed) |
| ordering of the MFT-dimer share by additive | GCys > Cys ≈ GSH ("GCys group exhibited the lowest MFT percentage"; "an unclear distinction ... Cys plays an equal role ... with GSH") | — | simple and multi-component models, 115 C, 60 min | — | text 3.1 | qualitative (the percentages themselves are Fig. 1b: figure_only) |
| direction of each volatile with additive dose | MeSH up then down (> 40 mg/mL); DMDS monotone down; DMTS up then down (GSH < Cys); MFT up; MFT-MFT up with MFT; MMFT flat; MMFT/MFT-MFT down, Cys > GSH | — | 1:3:3 ternary, Cys or GSH 0.5-75 mg/mL, 115 C, 60 min | — | text 3.2 | qualitative |
| sensory detection threshold of added Cys | Cys/Met = 6 (model, = 30 mg/mL), 8 (flavouring) | — | — | — | Table 1 | threshold (sensory) |
| MFT, MMFT, MFT-MFT levels and percentages (simple model); dose-response concentrations of six volatiles; flavouring levels; time courses (Fig. S1); single-factor sweeps (Fig. S2) | — | ng/mL | — | — | Figs. 1, 2, 4, S1, S2 | figure_only |

## What is new against `Zhang2024_extraction.md`

1. **Identity.** Same paper; the brief's "JAFC 2024" premise is wrong (see the note above). The two
   dossiers should be read together, with this one governing which numbers may enter a fit.
2. **Admissibility.** The wave-Z3 dossier's sections 2 (Fig. 1a/b), 3 (Fig. 2a-g, including the
   flux arithmetic) and 5 (Fig. 4) are 300-dpi read-offs. Under the repository's rule ("nothing from
   figures is transcribed") they are figure_only. In particular the **8.7 % dimer share of the Cys
   arm and 54.2 % of the GCys arm** that B17 section 6 (T3) uses as measured targets are read-offs;
   the printed text supports only the ordering GCys > Cys ≈ GSH and "MFT-MFT increasing with
   additive". The B17 target should be re-labelled figure-derived, or the values requested from the
   authors (raw data "upon reasonable request").
3. **Molar conversions** for every pot, the non-equimolar thiol loading (Cys 124 vs GSH 49 mmol/L),
   and the inference that the dose-response / RSM ternary carries **Met 5 mg/mL = 33.5 mmol/L**
   (from 2.5 and 3.3), one third of the 1:1:1 model — the old dossier treats the ternary's Met as
   unstated.
4. **The oxidant, as printed, for B17 section 6.** The paper measures no oxidant, adds none, and does
   not describe the reaction vessel's atmosphere. What it states: MFT-MFT is "generated by the
   oxidation of MFT"; the cystine arm (no free thiol) has the largest dimer share, attributed to
   "MFT being readily oxidized by free radicals"; Cys and GSH act "as an antioxidant against the
   attack from hydroxyl and peroxyl radicals", GSH the stronger; DMDS, the "MeSH oxidation
   product", falls monotonically with added thiol; and yet MFT-MFT *rises* with added thiol in step
   with MFT (3.2). So at pH 4.9 / 115 C the dimer tracks the thiol across a 150-fold additive range
   rather than saturating: by this paper the oxidant is not exhausted in the additive arms —
   consistent with an oxidant supplied continuously (sugar-derived radicals, or air in an
   undescribed headspace) rather than a fixed pool. The paper cannot say which. This is the opposite
   sign from B17's model, where the pool runs out; it supports B17's own suggestion of "a larger
   reservoir ... or a second oxidant" and says nothing about its size.
5. **The MeSH partition.** The one printed quantitative statement about the competition between
   MeSH's two sinks is the MMFT slope ratios (3 x and 2.5 x control in the first 90 min): adding a
   thiol that raises MFT raises the capture of MeSH into MMFT threefold while DMDS falls. As
   within-study ratios these can enter a fit of the MeSH -> MMFT vs MeSH -> DMDS branching once the
   engine has an MMFT species (not in the registry today).
6. **Two typographic hazards** not in the old dossier: the IPDS internal standard "0.0943 ng/mL"
   (nine orders below a usable amount) and the unstated temperature of the simple model (2.2).
7. **Evidence classes and registry keys** in the house-style table (section 4), which the wave-Z3
   dossier predates.

## 5. Flags

1. **Same paper as `Zhang2024_extraction.md`**; do not cite the two as independent sources.
2. **Figure-derived targets in B17.** See "What is new" 2; the 8.7 % and 54.2 % dimer shares are
   figure read-offs.
3. **The five rate numbers have no units** and the control's are back-calculated from rounded
   ratios; the old dossier's warning stands: do not ingest as numbers.
4. **Run length.** 2.3 says "up to 60 min"; the regression has a break "at 90 min" and a segment
   after it; Fig. S1's duration is not printed. The kinetic run is longer than the method describes.
5. **Internal standard amount** "1 µL IPDS (0.0943 ng/mL)" cannot be right; the concentrations in
   ng/mL rest on external curves referenced to this IS, so the scale is only as good as the
   unstated true IS amount (ratios within the study are unaffected).
6. **Thiamine salt unstated** (44.5 or 45.8 mmol/L); **buffer concentration unstated**; **vessel,
   fill, headspace, atmosphere unstated** for every model. The 40-mL vial with 10 mL is the SPME
   vessel, not the reactor.
7. **Non-equimolar additives**: 15 mg/mL of each gives 124 (Cys), 49 (GSH), 62 (GCys) mmol/L; the
   "equal weight" comparison confounds identity with thiol concentration (Cys has 2.5 x GSH's
   thiol). The dose series (0.5-75 mg/mL) partly repairs this for Cys vs GSH.
8. **Methionine at 100 mmol/L (or 33.5 in the ternary) and pH 4.9**: a MeSH supply far above any
   food; the MeSH "up then down" with additive (3.2) is a competition statement about a saturated
   source.
9. **"Second-order thermal procedure"** means two temperature stages (115 C then 95 C), not a
   kinetic order; the optimum is an RSM output over a 27-run design with R2 0.87-0.91, not a
   mechanism.
10. **Real-food arm** (pork flavouring) is semi-quantitative and figure-only; the 112-day DMTS
    storage series is figure-only.
11. **SI absent** (Tables S1-S7, Figs. S1-S4): every concentration-time series of this paper is in
    the SI or in figures. Request the raw data (offered "upon reasonable request").
12. **Registry gaps**: MMFT is the species this paper is about and has no key; methionine,
    thiamine (molecule), xylose, cysteine, glutathione, cystine likewise.
