# Vauthey, Milo, Frossard, Garti, Leser & Watzke 2000 (J. Agric. Food Chem. 48:4808-4816) — EXTRACTION (file stem `vauthey2000`)
### *"Structured Fluids as Microreactors for Flavor Formation by the Maillard Reaction"*
### Two pot designs, each heated ONCE (**100 °C, 4 h**, Pyrex tube, oil bath) in water, an L2 microemulsion and / or a monoglyceride cubic phase: (1) L-cysteine + furfural 2 + 2 mmol in 10 mL 0.5 M phosphate pH 5.0; (2) **L-cysteine 6.67 mmol + D-ribose 20 mmol in 10 mL 0.5 M phosphate pH 5.0**. Volatiles as µg per sample (duplicates, "SD < 20%"), internal standards on the FID trace, **uncorrected for response factor and recovery; no stable-isotope dilution**. In the cysteine + ribose **aqueous control MFT is "nd" and FFT "tr"**; the authors attribute the missing MFT to losses in their own work-up. **One heating time, no time course anywhere in the paper.**

**Read 2026-11-05 (first reader) from `data/articles/vauthey2000.pdf`** (9 pp., ACS typeset version of record, clean
text layer). Prose and both tables from `pdftotext -layout`. Tables 1 and 2 were checked cell by cell against 150 dpi
page renders (pp. 4812, 4814). Fig. 7 (GC/FPD, ribose/cysteine) was re-read at 300 dpi for which peaks are present
(scratchpad `vauthey/p-5.png`, `p-7.png`, `f7hi-7.png`). Figs. 2, 6 and 7 are chromatograms; Figs. 1, 3, 4, 5 and 8 are a
phase diagram, mass spectra, a synthesis scheme, a proposed pathway and a structure sketch. **No figure carries a
quantity that can be digitised** (sec. 4.3). There is no Supporting Information.

**Reservation (M1c amendment A1).** This paper was declared in advance as a candidate source for M1c signature (d)
and for H0's refusal of absolute thiol levels (`docs/reference/cases/M1c_amendment_A1_vauthey2000.md`, on branch
`m1c-case`, commit 17710774). Every aqueous-control value below carries the mark
**RESERVED: unseen by the engine (M1c A1 roles 1–2); no engine run on these pots until the relevant design is frozen**
(abbreviated **[RESERVED]** in the tables after its first full statement in sec. 3). No generator, probe or
prediction was run on any pot of this paper while this dossier was written.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited by the authors · **[Z]** derived by me, never printed ·
**[!]** flag. Provenance: **first reader only, second reader pending.** Nothing here is in a ledger, fit table or
benchmark.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | S. Vauthey, Ch. Milo, Ph. Frossard, N. Garti, M. E. Leser, **H. J. Watzke*** |
| Affiliation | Nestlé Research Centre, Nestec Ltd., Vers-chez-les-Blanc, Lausanne (Garti's Hebrew University affiliation is not printed) |
| Venue | *J. Agric. Food Chem.* **48** (10) (2000) **4808-4816** |
| DOI | **10.1021/jf991254a** (printed on p. 4808), matches the brief |
| Dates | received 19 Nov 1999; revised 5 Jun 2000; accepted 5 Jun 2000; published on the web 7 Sep 2000 |
| Funding / SI | none printed / none |

**Laboratory.** Nestlé Research Centre, Lausanne, in a colloid / material-science group (Watzke, Leser) working with
Garti (structured fluids). Its flavour isolation is the DFA one: the vacuum transfer of Sen et al. 1991 [C, Grosch
group]. `data/keys/laboratories.yml` has `nestle_blank` (Nestlé Research Center flavour group: Blank, Davidek,
Kraehenbuehl). **Proposal, not applied:** for leave-one-lab-out the conservative choice is to merge this paper into
`nestle_blank` (same centre, same site). Whether "Ch. Milo" is the Milo of the Grosch-group papers is **not checked**
[!]. The amendment's "a laboratory outside the DFA network" holds for the corresponding author. The analytical chain is
DFA-derived. Owner call.

## 1. What this paper answers for the model

| model need | what Vauthey 2000 gives | tag |
|---|---|---|
| FFT / MFT time course in an aqueous Cys + ribose pot (M1c sig. (d)) | **none**: one heating time (4 h at 100 °C) per arm | [M] |
| Absolute aqueous MFT / FFT (H0 refusal) | **MFT nd, FFT tr** in the aqueous Cys + ribose control. No number, no SIDA, no recovery | [M] |
| Aqueous furfural / norfuraneol at 100 °C, 4 h, pH 5 | furfural **351 µg**, norfuraneol **291 µg** per ~2 g aqueous sample (one time, FID, uncorrected) | [M] |
| Aqueous FFT from Cys + furfural | **1.8 µg** per ~1 g sample (≈ 8.7 × 10⁻⁵ mol/mol [Z]), one time, uncorrected | [M] |
| Matrix (structured-fluid) effect | cubic / water ≈ 2.4-2.5× (furfural, norfuraneol), ≈ 6.4× (FFT from furfural); MFT only in cubic | [M] / [Z] |
| Rate constants, Ea, pH law | **none** | — |

## 2. Conditions [M] (pp. 4809-4810), aqueous control first

### 2.1 Cysteine + ribose, AQUEOUS CONTROL ("reference") — the pot named in M1c A1

| item | as printed | derived [Z] |
|---|---|---|
| Stock | "L-cysteine (6.67 mmol) and D-ribose (20 mmol) were dissolved in a phosphate buffer (10 mL; 0.5 M; pH 5.0)" | nominal **0.667 mol Cys and 2.0 mol ribose per L of buffer** (ratio 1 : 3, the Hofmann & Schieberle 1995 ratio) |
| Charge heated | "**2 g of buffer solution alone (reference)**" | stock mass ≈ 10.5-10.7 g buffer (density 1.05-1.07 assumed, not printed) + 0.808 g Cys + 3.003 g ribose = 14.3-14.5 g, so **2 g holds ≈ 0.92-0.93 mmol Cys and ≈ 2.76-2.80 mmol ribose**. In-solution molarity ≈ 0.54 M Cys / 1.6 M ribose if the solids add ≈ 2.3 mL (partial volumes assumed) [!] |
| Water / buffer | sodium or potassium phosphate **not stated**; 0.5 M; water "twice distilled" | — |
| pH and scale | **pH 5.0 of the buffer before the reactants were dissolved**. Scale, electrode and temperature of measurement **not stated**. **pH after dissolving 0.67 M cysteine + 2 M ribose not reported; pH after heating not reported** [!] | 0.5 M phosphate at pH 5 sits between pKa₁ ≈ 2.1 and pKa₂ ≈ 7.2, so its buffer capacity is weak there |
| Vessel | **Pyrex tube, i.d. 2 cm** (stated for the cubic samples; the reference is "treated following the same procedure"). Length, closure and headspace volume **not stated** | 2 g ≈ 1.7 mL fills ≈ 0.55 cm of a 2 cm tube; headspace unknown |
| Closure | not stated. 100 °C was chosen partly "to avoid the use of pressurizable reaction vessels" [V] | an aqueous sample at 100 °C in an unpressurised tube may lose water and volatiles [!] |
| Heating | "thermally treated following the same procedure as mentioned above", i.e. the cysteine/furfural protocol: [vortex, **oil bath 140 °C for 1 min**, stir], then **4 h at 100 °C**. The 140 °C / 1 min step forms the cubic phase. **Whether the aqueous reference also got it is ambiguous** [!]. The heating medium for the 4 h is not restated (oil bath implied) | — |
| Time points | **one: 4 h** | — |
| Work-up | cooled. Then **4 g melted Dimodan PV added** to turn the aqueous sample into a cubic phase "to ensure maximal similarity", extracted with diethyl ether (3 ×, 50 mL total) **spiked with the internal standards**, vacuum transfer (two liquid-N₂ traps, 100 mPa, 1.5 h), Na₂SO₄, Vigreux to 2 mL, N₂ to 0.5 mL | the reference work-up mass is 2 g + 4 g = 6 g, **not** the 80 : 20 ratio of the 10 g cubic sample [!] |
| Replicates | Table 2 footnote: "**Mean values of duplicates (SD < 20%)**" | a bound, not a per-value SD |

### 2.2 Cysteine + furfural, AQUEOUS CONTROL

| item | as printed | derived [Z] |
|---|---|---|
| Stock | L-cysteine 2 mmol + furfural 2 mmol in phosphate buffer 10 mL, 0.5 M, pH 5.0 | nominal 0.2 mol each per L of buffer |
| Charge heated | "**1 g of the buffer solution containing the reactants**" | stock ≈ 10.9-11.1 g, so 1 g holds **≈ 0.18 mmol Cys and ≈ 0.18 mmol furfural** |
| Vessel / heating | Pyrex tube (i.d. 1.5 cm for the cubic samples); "an identical thermal treatment": **4 h at 100 °C** (the 140 °C / 1 min step is again ambiguous for the reference) [!] | — |
| Work-up | cooled. **4 g melted Dimodan PV added**, then as 2.1. **Internal standard: benzyl mercaptan, 5 or 10 µg** in the ether | 1 g + 4 g = the same 80 : 20 mass as the cubic sample |
| Replicates | Table 1 footnote: duplicates, SD < 20 % | — |

### 2.3 Structured-fluid arms (matrix evidence only; sec. 5)

- **Cubic phase:** 2 g stock (ribose series) or 1 g (furfural series) + 8 g or 4 g melted (60 °C) **Dimodan PV**
  (saturated vegetable-oil monoglyceride, > 98 %), i.e. **monoglyceride : water-phase 80 : 20**. Vortex, 140 °C oil
  bath for 1 min, stir to a homogeneous cubic phase, then **4 h at 100 °C**. The Dimodan PV / water (80 : 20) cubic
  phase holds at 80-120 °C and reverts to L2 above 120 °C [C, Krog 1998].
- **L2 (reversed) microemulsion** (furfural series only): Dimodan PV replaced by **Dimodan LS** (sunflower-oil
  monoglyceride, C18:2 67 %, C18:1 21 %), same 80 : 20 and 100 °C. There is no microemulsion arm in the ribose series.
- Both structured fluids "exhibited a water activity of ∼1" [V] (method and number not printed).

## 3. The assay [M] (p. 4810, p. 4811)

| item | as printed |
|---|---|
| Ribose/cysteine internal standards (spiked into the extraction ether) | **350 µg 5-methylfurfural** for MFT, FFT and furfural; **350 µg abhexon** (5-ethyl-4-methyl-3-hydroxy-2(5H)-furanone) for norfuraneol; **4.7 µg 3-acetyl-2,5-dimethylthiophene** for MFT-MFT. No standard is named for 2-methyl-3(2H)-furanone [!] |
| Furfural/cysteine internal standard | benzyl mercaptan, 5 or 10 µg |
| Signal | "ratio of peak areas using the FID trace and was **not corrected for possible differences in FID responses and recovery yields** of internal standard and analyte" (p. 4810) |
| Second description [!] | Results (p. 4811): "**A multiple standard addition method** was used to quantify six compounds", again with recovery "not taken into account". The two descriptions do not match, and no standard-addition levels or curves are printed |
| Stable-isotope dilution | **none** (SIDA is mentioned only as Schieberle & Hofmann's method [C]) |
| Recovery | **not measured, not corrected**. The authors name "losses and poor recovery during isolation from the cubic phase" [Q] |
| LOD | "∼0.5 µg/5 g of cubic phase sample" (ribose/cysteine, p. 4812). **No LOQ.** "tr" is never defined: Table 2 cites a **footnote d that is not printed** [!] |
| GC | HP 5890, **DB-1701** 30 m × 0.32 mm × 0.25 µm, effluent split by a Valco valve to **FID and FPD** (250 °C), He 80 kPa, 1 µL splitless; 20 °C (2 min) → 70 °C/min → 50 °C (2 min) → 6 °C/min → 180 °C → 20 °C/min → 240 °C (10 min) |
| Identification | GC-MS on DB-5 (Finnigan SSQ 7000), EI 70 eV and CI (isobutane); RI on DB-5 and DB-1701 against reference compounds. 2-(2-furyl)thiazolidine confirmed by synthesis from cysteamine + furfural (Fig. 4); N-(2-mercaptovinyl)-2-(2-furyl)thiazolidine **tentative** |
| Replicate noise | duplicates, **SD < 20 %** (both tables); no per-value SD |
| Free cysteine, free ribose, H₂S, pH after heating | **not measured** |

## 4. Data

### 4.1 Table 2 (p. 4814): ribose + cysteine, µg per sample, 4 h at 100 °C, mean of duplicates (SD < 20 %) [M]

Every value in the **water** column is: **RESERVED: unseen by the engine (M1c A1 roles 1–2); no engine run on these
pots until the relevant design is frozen.**

| compound (as printed) | cubic phase | **water [RESERVED]** | water, molar [Z] (per ≈ 0.93 mmol Cys / ≈ 2.78 mmol ribose) |
|---|---:|---:|---|
| 2-methyl-3-furanthiol (MFT) | 18.4 | **nd** | — (not detected; see 4.4) |
| 2-methyl-3(2H)-furanone | 36.7 | **13.1** | 0.131 µmol if C₅H₈O₂ (2-methyltetrahydrofuran-3-one, MW 100.1) [!] |
| furfural | 874 | **351** | 3.65 µmol = **1.3 × 10⁻³ mol/mol ribose** |
| 2-furfurylthiol (FFT) | 12.0 | **tr** (footnote d not printed) | — (no number) |
| 3-mercapto-2-pentanone | tr | **nd** | — |
| norfuraneol | 698 | **291** | 2.55 µmol = **9.2 × 10⁻⁴ mol/mol ribose** |
| MFT-MFT (bis(2-methyl-3-furyl) disulfide) | 8.3 | **tr** | — |
| MFT-FFT (mixed disulfide) | tr | **tr** | — |

Cubic arm, molar [Z], matrix only: MFT 0.161 µmol (**1.7 × 10⁻⁴ mol/mol Cys**), FFT 0.105 µmol (**1.1 × 10⁻⁴**),
MFT-MFT 0.0367 µmol, furfural 9.10 µmol, norfuraneol 6.12 µmol.

### 4.2 Table 1 (p. 4812): furfural + cysteine, µg per sample, 4 h at 100 °C, mean of duplicates (SD < 20 %) [M]

| matrix | FFT | 2-(2-furyl)thiazolidine | N-(2-mercaptovinyl)-2-(2-furyl)thiazolidine |
|---|---:|---:|---:|
| **water [RESERVED]** | **1.8** | **nd** (text: "detectable in only one sample") | **nd** |
| microemulsion | 8.7 | 42 | 11.7 |
| cubic phase | 11.6 | 240 | 123 |

Text (p. 4811): from relative peak areas to benzyl mercaptan, "∼2 µg of FFT" in water, ∼9 µg in the microemulsion and
∼12 µg in the cubic phase, "a 6-fold increase" [M/Q]. Water FFT molar [Z]: 15.8 nmol per ≈ 0.18 mmol Cys =
**≈ 8.7 × 10⁻⁵ mol/mol** (0.0087 mol %) **[RESERVED]**. Ratios [Z]: microemulsion / water 4.8, cubic / water 6.4;
cubic / microemulsion 1.3 (FFT), 5.7 (thiazolidine), 10.5 (mercaptovinyl), as the text says (1.3, 6, 10).

### 4.3 Figures: nothing quantitative to digitise [F]

- **Fig. 7A (GC/FPD, ribose/cysteine in water) [RESERVED]**, read at 300 dpi, ±300 abundance units. Baseline
  ≈ 44,000. Peaks: FFT at ≈ 9.0 min, apex ≈ 44,700 (just above noise); MFT-MFT at ≈ 23.9 min, ≈ 45,400; **MFT-FFT at
  ≈ 26.6 min, ≈ 49,000**; unlabelled peaks at ≈ 18.2 (≈ 45,700) and ≈ 27.7 min (≈ 52,800). **There is no peak at
  MFT's retention (≈ 7.3 min).** These are FPD peak heights. The FPD response is non-linear in sulfur (≈ quadratic),
  and the paper prints no calibration, so they are presence / absence evidence, not amounts. The water trace shows
  the mixed disulfide as the largest labelled sulfur peak, although Table 2 calls it "tr".
- Fig. 7B (cubic): MFT and FFT off scale (> 230,000), MFT-MFT ≈ 112,000, MFT-FFT ≈ 56,000 over a baseline of ≈ 47,000.
- Figs. 2 (FPD, furfural series) and 6 (FID, ribose series): chromatograms, with no amounts beyond Tables 1 and 2.

### 4.4 The authors on the missing aqueous MFT [Q] (p. 4812)

MFT in the reference is "surprising, but not contradictory" to Hofmann & Schieberle 1995. It "was detectable in the
reference when the reaction system was directly extracted with diethyl ether", i.e. without forming a cubic phase
first, "and when the isolation of volatiles was done by vacuum distillation". **No amount is printed for those
extractions.** So the aqueous "nd" reports **a work-up loss, not formation below a floor**. It cannot serve as an
upper bound either. About 50 % of the cubic MFT is said to be dimerised; the dimer is "at least partly due to
artifacts" [Q] (see Conflicts).

### 4.5 Sensory [Q] (p. 4810)

Aqueous reaction (judged after conversion to a cubic phase): "weak in aroma, reminiscent of meat and lard with some
burnt character" **[RESERVED]**. Cubic: more intense, "rubber, eggy, roast chicken-like".

## 5. Structured-fluid arms: matrix evidence only

- The same reactant mass in 80 % monoglyceride gives **2.4-2.5× the furfural and norfuraneol** (874 vs 351; 698 vs 291
  µg), and MFT 18.4 µg against nd in water. FFT is 12.0 µg against tr. In the furfural series the cubic phase gives
  6.4× the FFT of water and large amounts of two thiazolidines that water barely forms.
- The authors' readings [Q]: macroscopic a_w ≈ 1 in both fluids, so a_w is not the cause. They propose interfacial
  concentration of the reactants, a cage effect, a changed cysteine environment ("cubic catalysis"), and MFT
  partitioning into the chain region, which protects it. None of this is tested quantitatively.
- **Confound [!]:** the aqueous arms are converted into a cubic phase only AFTER heating. Work-up recovery therefore
  differs between arms (the authors say MFT is lost in the reference), so the "enhancement" ratios mix formation and
  recovery. Nothing in the engine models a lipid mesophase; these arms are not engine pots.

## 6. Pre-declared roles (M1c amendment A1): verdicts

Criteria quoted from the amendment and not reinterpreted.

**Role 1: M1c signature (d), FFT shape. Verdict: NOT MET.**

| sub-criterion | finding | verdict |
|---|---|---|
| FFT and / or MFT at ≥ 3 heating times in one pot design | **one heating time** (4 h at 100 °C) in every arm; the aqueous Cys + ribose control has **FFT "tr", MFT "nd"** at that one time | **NOT MET** |
| stated Cys and ribose charges, temperature, time, pH / buffer on a stated scale | 6.67 mmol Cys + 20 mmol ribose in 10 mL 0.5 M phosphate pH 5.0, 2 g heated (≈ 0.93 mmol Cys, ≈ 2.78 mmol ribose [Z]); 100 °C; 4 h. **pH scale not stated, pH after mixing not measured** | **PARTLY** |
| replicate noise | duplicates, "SD < 20%", a bound only; **no SD can exist for nd / tr** | **PARTLY** |

Under the amendment's own fallback ("If not met (e.g. one time point), the source is evidence only and the case's
verdict stays development-grade"), Vauthey 2000 is **evidence only** for M1c.

**Role 2: second-laboratory absolute thiol data for H0's refusal. Verdict: NOT MET.**

| sub-criterion | finding | verdict |
|---|---|---|
| absolute MFT / FFT from the aqueous control | Cys + ribose: **MFT nd, FFT tr, no number** (and MFT-MFT tr, MFT-FFT tr). The Cys + furfural aqueous control prints **FFT 1.8 µg**, but it is not a Cys + ribose pot | **NOT MET** |
| quantification method stated (ideally SIDA; which standards) | 5-methylfurfural 350 µg (MFT, FFT), benzyl mercaptan 5-10 µg (furfural series); FID area ratio; **no SIDA**; a contradictory "multiple standard addition" sentence | stated, but not SIDA |
| recovery | **"not corrected for ... recovery yields"**; recovery never measured; the authors report MFT losses in the aqueous work-up | **NOT MET** |

**Consequence.** Neither role is met, so no pre-registration on these pots is triggered. The reservation marks stay
until the coordinator releases them. That is a call for the M1c / H0 owners, not for this dossier.

## 7. Conflicts and agreements with the repo's records

- **[!] The queue description overstates the paper.** The amendment quotes "MFT and / or FFT quantified" with an
  aqueous control. In the aqueous Cys + ribose control **neither is quantified** (nd / tr). `docs/guides/LITERATURE_SEARCH_2026-10_ROUND3.md`
  (line 250) says "FFT 5-7x higher in L2/cubic phases than water". Table 1 gives **4.8× (microemulsion) and 6.4×
  (cubic)**, and only for Cys + **furfural**. For Cys + ribose, water is "tr", so no ratio exists.
- **[!] In-paper arithmetic.** "8.3 µg of MFT-MFT indicates that ∼50% of the MFT was dimerized." By mass 8.3 / 18.4 =
  0.45 (dimer : free MFT). As the fraction of total MFT equivalents held in the dimer it is 2 × 0.0367 / (0.161 +
  0.073) = **0.31** [Z]. "~50 %" holds only as a mass ratio.
- **[!] Two quantification descriptions** (internal standard on FID vs "multiple standard addition"), and a missing
  footnote d for "tr". See sec. 3.
- **Agreement, `hofmann1995_extraction.md` / `hofmann1996_extraction.md`:** disulfides of MFT in solvent-extracted
  ribose/cysteine pots are largely work-up artefacts. Vauthey says the same of its MFT-MFT [Q].
- **No direct conflict with the DFA aqueous pot** (`schieberle1998_extraction.md`: MFT 19.8 µg, FFT 12.1 µg per 3.3
  mmol Cys at 145 °C / 20 min, by SIDA). The conditions differ: 100 °C / 4 h, 20× the concentration, ≈ 0.93 mmol
  Cys. Vauthey's aqueous "nd" is a work-up loss by the authors' own account, so it does not contradict a SIDA
  detection.
- **Consistent with `yaghmur2005_extraction.md`** (water arm FFT "< 1 %" of furfural; `YAGHMUR_FFT_SHARE_CEILING` = 0.012):
  Vauthey's aqueous Cys + furfural FFT ≈ 8.7 × 10⁻⁵ mol/mol [Z, RESERVED] is far below that ceiling. It is not a
  test of it (different T, time and work-up), and it is reserved.
- **Compound keys:** MFT (`2_methyl_3_furanthiol`), FFT (`2_furfurylthiol`), `furfural`, `norfuraneol`, bis(2-methyl-3-furyl)
  disulfide and `2_methyltetrahydrofuran_3_one` resolve in `data/keys/compounds.yml`. The printed
  "2-methyl-3(2H)-furanone" is assumed to be the latter (dihydro form) [!]. 2-(2-furyl)thiazolidine, its
  N-(2-mercaptovinyl) derivative and the MFT-FFT mixed disulfide have no key. Nothing was keyed (no data rows).

## 8. Held dossiers that already cite this paper

| dossier | what it says about Vauthey 2000 |
|---|---|
| `yaghmur2002_extraction.md` | Methods: products identified "tentatively" by library and literature, citing Vauthey 2000 (its ref. 26). Qualitative claims: 2-(2-furanyl)thiazolidine may react with mercaptoacetaldehyde (from cysteine Strecker degradation) to give the mercaptovinyl-thiazolidine, citing Vauthey 2000. **Matches** Vauthey's Fig. 5, which is itself labelled "hypothetical". Vauthey confirmed the thiazolidine by synthesis, but the mercaptovinyl compound is tentative there too. |
| `yaghmur2005_extraction.md` | Identity table, "Lineage": "the microreactor idea from [18] Vauthey et al. (Nestlé)". Methods and Flag 5: products "tentatively identified by comparing their mass spectra with ... previously published literature [18]", i.e. Vauthey 2000 is the reference library for the furfural-cysteine products, and no authentic standards were used. Correct as a description; Vauthey's own identification of the mercaptovinyl compound is tentative, so the chain of tentative identifications is two deep [!]. |
| `docs/guides/LITERATURE_SEARCH_2026-10_ROUND3.md` (not a dossier) | Row 250: "FFT 5-7x higher in L2/cubic phases than water", decision SKIP. See sec. 7. |

## 9. What it can and cannot be used for

- **Cannot:** any fit row, rate constant, time course, absolute aqueous thiol level, pH law or Ea. Cannot carry M1c
  signature (d) or lift H0's refusal (sec. 6).
- **Can (evidence only, after release of the reservation):** a qualitative matrix statement (a lipid mesophase raises
  recovered furfural, norfuraneol, FFT and MFT 2.5-6× at 100 °C, with recovery confounded). One-time aqueous furfural
  / norfuraneol / FFT levels at 100 °C, 4 h, pH 5.0 (start), uncorrected for recovery, are at most an order-of-magnitude
  plausibility check.
- **Proposed ledger entries: none.** No row here passes the units / quantification checks a second reader would apply.
