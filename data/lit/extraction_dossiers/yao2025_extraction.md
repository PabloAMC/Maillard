# Yao 2025 — EXTRACTION (methional 4300 µg/mL = 41 mmol/L alone, and methionine + glucose 100 + 100 mmol/L, in 40 mmol/L phosphate pH 6.24, sealed 20-mL vials, 100 C, 60 / 80 / 100 / 120 min, ± EGCG 50 mmol/L; sulfur volatiles by HS-SPME-GC×GC-TOF-MS — every model-system number is in Fig. 4 or the SI; green tea VSC levels in Table 1)
### The only pot on disk that heats methional by itself in buffer and follows methanethiol, DMDS, DMTS and the (methylthio)methyl sulfides against time — and prints none of the values.

**Source on disk:** `data/articles/Yao2025.pdf` (11 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Yao2025.txt`); Table 1 (tea infusions) came through clean and is
re-typed below. Fig. 4 (the two model systems' VSCs vs time, panels a-g), Fig. 1 (aroma profiles),
Fig. 3 (dose-response) were not read. The Supporting Information (Fig. S1 workflow, S2 adduct
spectra; Table S1 tea calibration, **S2 model calibration curves**, S3 UPLC gradient, **S4 VSCs
identified in the four models**) is NOT on disk. "Data will be made available on request."

## 0. Identity

| field | value |
|---|---|
| Title | "Characterization of cooked off-flavor volatile sulfur-containing compounds in green tea and their thermal inhibition via (−)-epigallocatechin gallate" |
| Authors | Xin Yao, Yangyang Li (equal), Jun Tang, Jieyao Yu, Yanyan Zhang, Xiaochun Wan, Guoyu Zhang, Xiaoting Zhai* (Anhui Agricultural University, Hefei; Hohenheim; Anhui University of Chinese Medicine) |
| Venue | Food Chemistry 463 (2025) 141143 (PII S0308814624027936). Received 16 May 2024, revised 26 August 2024, accepted 3 September 2024, online 5 September 2024 |
| DOI | 10.1016/j.foodchem.2024.141143 |
| Naming | VSCs = volatile sulfur-containing compounds; DMS / DMDS / DMTS / DEDS = dimethyl sulfide / disulfide / trisulfide, diethyl disulfide; Met-Glu = methionine + D-glucose (the authors write "Glu" for glucose); EGCG; MG-ARP; REF / OF1-3 = the four teas; GC-SCD = sulfur chemiluminescence detection |
| Companions | Zhang, Wang & Cao 2023 (`zhang2023_extraction.md`, cited for methional as the precursor of MeSH); Cheng 2020 (`cheng2020_extraction.md`); Pan 2021 (catechin-methional interaction, not on disk); Totlani & Peterson 2005 (catechin-sugar-fragment adducts, cited) |

## 1. Why it matters

The methionine chain needs the step methional -> methanethiol measured without a Maillard pot
supplying the methional. Section 2.7 of this paper does exactly that: methional alone in phosphate
buffer at 100 C for 60-120 min in a sealed vial, with the sulfur products identified and quantified
(HS-SPME-GC×GC-TOF-MS with an internal standard and calibration curves), and the same with EGCG
added. It also runs Met + glucose at 100 + 100 mmol/L under the same conditions. **But the paper
prints no concentration for either model**: the time courses are Fig. 4 (a-g), the identities are
Table S4, the calibration is Table S2. What the text does give is the product list (which sulfur
species a methional pot makes at 100 C), the direction of each species with time, and the
qualitative EGCG effect (EGCG "consumed the precursor methional completely", Table S4). So for the
repository this is a data-request target and a species list, not yet a rate source. The tea part
(Table 1) is a set of levels in infusions with a sensory threshold for DMTS (0.4 µg/L) — background
for off-note thresholds, not for the kinetic lane.

## 2. Methods as they matter to a model

- **Met-Glu model (2.6).** "Met (0.1 mmol) and Glu (0.1 mmol) were dissolved in a phosphate buffer
  solution with a pH value of 6.24 (0.04 mol/L, 1 mL) and then transferred into 20 mL headspace
  vial. The thermal reaction started when heating the substrate solution in a thermostatic metal
  bath at 100 C for 60, 80, 100, and 120 min, respectively. In addition, EGCG (0.05 mmol) was added
  to another thermal model solution of Met-Glu." So **[Met] = [Glc] = 100 mmol/L**, **[EGCG] = 50
  mmol/L** (22.9 g/L; solubility at that level is doubtful, Flags 6), phosphate 40 mmol/L pH 6.24,
  **1 mL liquid under 19 mL headspace, sealed** (the vial is the reactor and the SPME vessel: no
  transfer, no opening). Cooled to room temperature before extraction.
- **Methional model (2.7).** "a headspace vial contained a phosphate buffer solution with methional
  (4300 µg/mL), while another vial contained phosphate buffer solution with EGCG (0.05 mmol) and
  methional (4300 µg/mL). The reaction condition was set as the same with thermal model of Met-Glu."
  So **methional 4.3 g/L = 41.3 mmol/L** (M 104.17), presumably in 1 mL of the same buffer (volume
  not restated), 100 C, 60 / 80 / 100 / 120 min, ± EGCG 50 mmol/L.
- **Extraction and separation.** HS-SPME, PDMS/DVB 65 µm, 30 min equilibration with stirring at 35 C,
  30 min extraction at 35 C, desorption 250 C 5 min; GC×GC-TOF-MS (Agilent 8890 / 7250; BPX50 30 m x
  0.25 mm x 0.25 µm and BPX5 1 m x 0.1 mm x 0.1 µm; INSIGHT flow modulator 2.4 s; 40 C 2 min -> 260
  C at 6 C/min, 5 min; 50 Hz, m/z 50-500; RI with C7-C40). The same paragraph gives "injector 310 C,
  1 µL splitless", which belongs to liquid injection of the SAFE distillates, not the SPME runs.
- **Quantification in the models.** "VSCs in thermal reaction models were quantitated via the
  internal standard and 2-heptene-3-pentanone was selected as the internal standard. The
  calibration curves ... were performed by ... HS-SPME combining with comprehensive GC×GC-TOF-MS ...
  Table S2." The internal standard name as printed is not a chemical name; the reference-odorant
  list includes 2-methyl-3-heptanone (PubChem 25,611), a common internal standard, and that is
  probably what was used (Flags 4). Calibration curves not on disk. So the model numbers, when
  obtained, are IS-normalised, calibrated concentrations (unit not stated in the main text).
- **Tea infusions (for the summary).** 10 g leaf / 500 mL boiling water 4 min; DCM liquid-liquid
  extraction (1.2 L), SAFE at 42 C, concentration to ~200 µL; GC-SCD (DB-1 30 m x 0.32 mm x 1 µm)
  with **diallyl sulfide** as internal standard and response factors from five mass ratios (Table
  S1); DMS by GC-MS. OAV against Leibniz-LSB@TUM thresholds. Sensory panel of 25, 0-3 scale.
- **Replicates.** "More than three parallel experiments"; mean ± SD; Duncan p < 0.05; tea RSD < 20 %.

## 3. Tables re-typed

### Table 1. "Concentration and OAVs of VSCs detected in REF and OF samples" (tea infusions, µg/L)

Mean ± SD; letters compare the four teas in a row; "–" = not detected; DMS by GC-MS, the rest by
GC-SCD; OAVs as printed.

| odorant | off-odor note | REF concn. | OAV | OF1 concn. | OAV | OF2 concn. | OAV | OF3 concn. | OAV |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Dimethyl sulfide | cooked cabbage-like, pungent, rotten | 77.98 ± 4.88 c | 260 | 65.03 ± 3.96 c | 217 | 447.51 ± 49.76 a | 1492 | 383.47 ± 35.01 b | 1278 |
| Dimethyl trisulfide | fishy, cooked cabbage-like | – | – | 3.42 ± 0.475 a | 346 | 2.21 ± 0.385 b | 223 | 0.48 ± 0.014 c | 49 |
| Methanethiol | cooked cabbage-like | 6.51 ± 1.22 a | 11 | 5.08 ± 0.279 b | 8.6 | 3.47 ± 0.328 c | 5.9 | 5.38 ± 0.654 ab | 9.1 |
| Diethyl disulfide | roasted onion-like, rotten | 0.104 ± 0.0053 c | 5.2 | 0.32 ± 0.009 a | 16 | 0.18 ± 0.015 b | 9 | 0.15 ± 0.005 c | 7.5 |
| Dimethyl disulfide | cooked onion-like | 0.004 ± 0.0002 d | < 1 | 0.21 ± 0.014 a | < 1 | 0.09 ± 0.015 c | < 1 | 0.18 ± 0.029 b | < 1 |
| 2-(Methylsulfanyl)propane | cooked cabbage-like | 0.012 ± 0.0008 d | < 1 | 0.039 ± 0.0012 b | < 1 | 0.024 ± 0.002 c | < 1 | 0.053 ± 0.007 a | < 1 |
| 1-(Methylsulfanyl)propane | cooked cabbage-like | 0.004 ± 0.0003 b | < 1 | 0.031 ± 0.003 a | < 1 | – | – | – | – |
| Furfuryl methyl sulfide | cooked cabbage-like | 0.03 ± 0.005 b | < 1 | 0.095 ± 0.002 a | < 1 | 0.097 ± 0.009 a | < 1 | 0.10 ± 0.002 a | < 1 |
| 2-Acetylthiazole | roasty, popcorn-like | 0.005 ± 0.0007 b | < 1 | 0.123 ± 0.019 a | < 1 | 0.022 ± 0.005 b | < 1 | 0.025 ± 0.001 b | < 1 |
| Benzothiazole | rubber-like, cooked cabbage-like | 0.04 ± 0.003 c | < 1 | 0.11 ± 0.019 a | < 1 | 0.092 ± 0.006 b | < 1 | 0.10 ± 0.007 ab | < 1 |

Implied thresholds (concn./OAV, mine): DMS ≈ 0.30 µg/L, DMTS ≈ 0.0099 µg/L, methanethiol ≈ 0.59
µg/L, DEDS ≈ 0.020 µg/L.

### Model systems — everything the main text prints (Fig. 4 and Table S4 are FIGURE-ONLY / SI)

- **Species.** "a total of seven and five VSCs were identified in the Met-Glu and methional reaction
  models, respectively. They were DMDS, 2-(methylsulfanyl)propane, bis(methylthio)methane,
  methanethiol, DMTS, furfuryl methyl sulfide, and methyl (methylthio)methyl disulfide." The five of
  the methional model are not listed separately; the text names DMDS, DMTS, bis(methylthio)methane
  and methyl (methylthio)methyl disulfide as present in it; the fifth is not named in the main text
  (methanethiol is the obvious candidate, not asserted). "after adding EGCG, only four and two VSCs
  were detected" (Met-Glu + EGCG; methional + EGCG); in the methional + EGCG pot
  "bis(methylthio)methane and methyl (methylthio)methyl disulfide were totally inhibited".
- **Directions with time, Met-Glu, 60-120 min, 100 C.** DMDS, methional and furfuryl methyl sulfide
  "showed slight increases" (Fig. 4a, c, d); DMTS "a quite flat decrease" (4b); 2-(methylsulfanyl)-
  propane and bis(methylthio)methane "generated after heating 80 and 100 mins" and level at the end
  (4e, f). With EGCG: the last two absent; "dramatic amount declines of DMTS, methional, furfuryl
  methyl sulfide"; the DMTS inhibition "diminished slightly after heating 100 mins".
- **Methional model.** "prominent decreases of VSCs (DMDS, DMTS, and bis(methylthio)methane) in
  methional model with and without EGCG (Fig. 4 g)"; "additional EGCG consumed the precursor
  methional completely" (Table S4). An EGCG-methional adduct (product I or II, Fig. 5b) and an
  EGCG-glucose adduct were found by UPLC-Q-Exactive-Orbitrap-MS (Fig. S2).
- **Sensory of the models (0-3 scale).** Methional vs EGCG-methional: cooked potato-like 2.6 -> 1.0;
  cooked cabbage-like 1.4 vs 1.7, cooked onion-like 0.7 vs 1.0, rotten 1.5 vs 1.7 (the latter three
  higher with EGCG); "the same appearances" in Met-Glu vs EGCG-Met-Glu (values not printed).
- **Tea sensory.** Cooked note 1.4 / 1.3 / 1.2 (OF1-3) vs 0.6 (REF); addition tests 1.9 / 2.0 / 1.5;
  DMTS off-flavour threshold in the tea matrix **0.4 µg/L** ("hawthorn-like" below, "fishy / cooked
  cabbage-like" at 0.4-3.9 µg/L); DMS acceptable below 169.5 µg/L, plateau 169.5-352.5, pungent at
  447.5 µg/L.

## 4. Kinetic numbers the repository can use

Registry mapping: methional -> `methional`; methanethiol -> `methanethiol`; dimethyl disulfide ->
`dimethyl_disulfide`; dimethyl trisulfide -> `dimethyl_trisulfide`; dimethyl sulfide, diethyl
disulfide, bis(methylthio)methane, methyl (methylthio)methyl disulfide, 2-(methylsulfanyl)propane,
1-(methylsulfanyl)propane, furfuryl methyl sulfide, 2-acetylthiazole, benzothiazole, methionine,
glucose, EGCG, gallic acid -> not in registry.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| DMDS, DMTS, bis(methylthio)methane, methyl (methylthio)methyl disulfide (and one unnamed species) vs time from methional alone, ± EGCG | — | not stated (IS-calibrated) | methional 41.3 mmol/L, 40 mmol/L phosphate pH 6.24, 1 mL in sealed 20-mL vial, 100 C, 60 / 80 / 100 / 120 min; EGCG 50 mmol/L | — | Fig. 4g; Table S4 | figure_only (SI not on disk) |
| DMDS, methional, DMTS, furfuryl methyl sulfide, 2-(methylsulfanyl)propane, bis(methylthio)methane, methanethiol vs time from Met + Glc, ± EGCG | — | not stated | 100 + 100 mmol/L, same buffer, 100 C, 60-120 min | — | Fig. 4a-f; Table S4 | figure_only |
| methional remaining with EGCG | "consumed ... completely" | — | methional + EGCG pot | — | text 3.9 citing Table S4 | level_only (qualitative) |
| species list of a methional-alone pot at 100 C | DMDS, DMTS, bis(methylthio)methane, methyl (methylthio)methyl disulfide (+1) | — | as above | — | text 3.7 | qualitative (usable as a product-set constraint: the chain reaches the trisulfide and the CH2-bridged sulfides without any sugar) |
| direction with time, Met-Glu, 60-120 min | DMDS, methional, FMS up; DMTS down; 2-(MeS)propane and bis(MeS)methane appear after 80-100 min | — | 100 C | — | text 3.7 | qualitative |
| tea infusion VSC levels (10 teas x 4) | Table 1 | µg/L infusion (10 g / 500 mL) | brewed 4 min at 100 C | — | Table 1 | level_only (food matrix, not a model) |
| DMTS off-flavour threshold in green-tea matrix | 0.4 | µg/L | sensory dose-response, 25 panellists | — | text 3.5, Fig. 3d | threshold (sensory) |
| DMS acceptable ceiling in tea matrix | 169.5 | µg/L | same | — | text 3.5 | threshold (sensory) |

## 5. Flags

1. **No model-system concentration is printed.** Fig. 4 (seven panels) and Table S4 hold the whole
   kinetic content; Table S2 holds the calibration. The paper is a data-request target ("Data will
   be made available on request"): ask for the Fig. 4g series (methional alone, four times, ±
   EGCG) with the Table S2 curves and the unit. Until then the dossier carries no rate and no
   level for the models.
2. **Methional load 41 mmol/L (4.3 g/L)** — three to four orders of magnitude above any food or
   pot level in the corpus (Deng 2022: 2 µmol/L; Pan 2025: 6 µmol/L). Self-reactions of methional
   (aldol, oxidation to the acid, the CH2-bridged sulfides bis(methylthio)methane and methyl
   (methylthio)methyl disulfide, which need formaldehyde or a second methional-derived C1) are
   favoured at this load; the product spectrum and any rate from it do not scale linearly down.
3. **Sealed vial, 1 mL under 19 mL headspace, extracted at 35 C after cooling**: the reported amounts
   are headspace-derived totals of a closed system, good for mass balance but partition-weighted
   towards MeSH and DMS; the EGCG pots change the matrix (50 mmol/L of a polyphenol) and with it
   the partition — part of the "inhibition" may be retention, as the authors' own citation of Pan
   2021 concedes.
4. **Internal standard "2-heptene-3-pentanone"** is not a chemical name; 2-methyl-3-heptanone is in
   the reagent list and is the likely IS. The unit of the model concentrations is not stated in the
   main text.
5. **No time zero and no point before 60 min**; the fast early chemistry of methional at 100 C is
   not resolved; "slight increases" over 60-120 min may be the tail of a curve.
6. **EGCG at 50 mmol/L in 40 mmol/L phosphate** exceeds EGCG's usual aqueous solubility at room
   temperature (a few g/L); at 100 C EGCG hydrolyses to gallic acid + EGC (the authors invoke this
   for the adduct chemistry) and acidifies the pot — the pH of the EGCG pots after heating is not
   reported.
7. **Phosphate 40 mmol/L pH 6.24**: phosphate catalyses the Amadori and Strecker steps; the buffer
   is part of the rate.
8. **Methanethiol in the methional pot** is not named in the main text as one of the five species;
   do not assert its presence from this paper.
9. **The 1 µL splitless injector setting** in the GC×GC paragraph belongs to the liquid (SAFE)
   samples; the SPME desorption is 250 C / 5 min.
10. **Tea numbers**: the infusion is 10 g / 500 mL DCM-extracted and SAFE-concentrated, a different
    workflow from the models; RSD < 20 %; OAVs use water thresholds. Background only.
11. **Registry**: eight of the eleven sulfur species in this paper have no key; DMS is the largest
    gap if tea or juice off-notes are ever modelled.
