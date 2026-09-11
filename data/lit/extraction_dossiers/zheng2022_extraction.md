# Zheng 2022 — EXTRACTION (**COUNTER-EVIDENCE**) — GSH 0.5 mM + methylglyoxal / glyoxal / 3-deoxyglucosone 0.5 mM, cell-free, **25 mM potassium phosphate pH 7.4, 37 °C, 6 h**, GSH and GSSG both quantified by LC-TQ-MS; plus a U2OS Nrf2/Cytotox CALUX reporter series with GSH modulated by NAC and BSO

### THE HEADLINE, AND IT CUTS AGAINST THE PRE-REGISTRATION: in the one experiment in this batch that put a thiol and an α-dicarbonyl in a pot and looked for the **disulfide**, the disulfide barely appeared. "These reactions were monitored for 6 h during which **the formation of GSSG from GSH was limited as confirmed by LC-TQ-MS (data not shown)**" (§3.2, p. 7). The thiol went to an **adduct** — MGO-GSH at m/z 380, GO-GSH at 366, hydrous GO-GSH at 384, 3-DG-GSH at 470, all identified by LC-TOF-MS — and then **stopped**: GSH fell by ca. 18.1 % (MGO), 8.6 % (GO) and 1.7 % (3-DG) essentially instantaneously at equimolar loading and did not fall further over six hours. **But the word "limited" is never given a number**: the GSSG channel was monitored by four MRM transitions and the result is "data not shown". This is a real negative result whose magnitude the paper does not publish.

**Source on disk:** `data/articles/Zheng2022.pdf` (16 pp., Nutrients 2022, 14, 1364; open access CC BY).
Read from the `pdftotext -layout` text layer. The PDF carries a **doubled text layer** — a proof/"FOR PEER REVIEW" version is overlaid on the typeset version, so most passages appear twice with slightly different line breaks. Every quotation below was taken from the typeset layer (the one whose running head reads `Nutrients 2022, 14, 1364` with page `n of 16`) and checked against its proof twin; they agree on every number.
**This paper contains NO TABLES.** Everything is figures (Figs. 1–11) and running text. There is **no supplementary material** (no "Supplementary Materials" statement appears), and the Data Availability Statement says the data "are available on request from the corresponding author".
**Repo status before this dossier:** Zheng 2022 is cited nowhere in `src/kinetic_core/`, nowhere in `data/lit/reaction_rules.yml`, and has no extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "The Influence of Intracellular Glutathione Levels on the Induction of Nrf2-Mediated Gene Expression by α-Dicarbonyl Precursors of Advanced Glycation End Products" |
| Authors | Liang Zheng (corresponding, liang.zheng@wur.nl), Katja C. W. van Dongen, Wouter Bakker, Ignacio Miro Estruch, Ivonne M. C. M. Rietjens — Division of Toxicology, Wageningen University and Research, Stippeneng 4, 6708 WE Wageningen, The Netherlands |
| Venue | **Nutrients 2022, 14, 1364.** Academic Editor Leni Rose Rivera. Received 1 March 2022; accepted 21 March 2022; published 24 March 2022 |
| DOI | **10.3390/nu14071364** — printed on the title page in two places (the Citation block and the footer line) exactly as `https://doi.org/10.3390/nu14071364` |
| Funding / interest | China Scholarship Council scholarship No. 202008510115 (L. Zheng). "The authors declare no conflict of interest." |
| Computational content | **None.** No DFT, no molecular modelling, no docking anywhere in this paper, so the standing no-DFT policy excludes nothing here. |
| The α-dicarbonyls | **methylglyoxal (MGO)** 40 % in water, **glyoxal (GO)** 40 % in water (both Sigma-Aldrich), **3-deoxyglucosone (3-DG)** purity >99 % (Toronto Research Chemicals) |
| The thiol | **L-glutathione reduced (GSH)**, purity ≥98 %, Sigma-Aldrich. **L-glutathione oxidized (GSSG)**, purity ≥98 %, was also purchased — i.e. an authentic GSSG standard was in hand |
| The companion paper | Zheng 2023 (`zheng2023_extraction.md`) is the same first author and the same Wageningen group, on kaempferol vs GSH as MGO scavengers in SH-SY5Y cells |

## 1. Why it matters

**This is the paper in the batch that argues against `kinetic_core_b27_prereg.md`, and it should be
read before the other three.**

The pre-registration proposes a new step `ch_redox_dicarbonyl` — **the pot's own α-dicarbonyls,
reduced to hydroxyalkanones, oxidising thiols to disulfides** — taken from Whitfield & Mottram 1999,
which "prints no rate, no order and no barrier" for it (pre-reg §1). Zheng 2022 is the only paper in
this batch that ran exactly that experiment: a thiol and an α-dicarbonyl in the same buffer, with a
mass-spectrometric assay watching **both** the adduct channel and the disulfide channel, for six
hours.

**What it found, in order of importance to the wave:**

1. **The disulfide channel was minor.** §3.2 (p. 7): "These reactions were monitored for 6 h during
   which the formation of GSSG from GSH was **limited** as confirmed by LC-TQ-MS." The GSSG assay was
   not an afterthought — §2.5 (p. 4) prints **four MRM transitions** for GSSG (m/z 613.15 → 355.05,
   231.00, 484.15, 177.10) alongside the four for GSH, and an authentic GSSG standard was bought
   (§2.1). The channel was instrumented and it did not deliver.
2. **The adduct channel was the answer, and it was identified structurally.** Four GSH conjugates by
   LC-TOF-MS with the diagnostic neutral losses (§3.1, pp. 6–7; Fig. 11's scheme, p. 12): MGO-GSH
   *m/z* 380 → 308 (loss of 72, one MGO), GO-GSH *m/z* 366 → 308 (loss of 58, one GO), a **hydrous**
   GO-GSH *m/z* 384 → 366 (loss of H₂O), 3-DG-GSH *m/z* 470 → 308 (loss of 162, one 3-DG). The
   Discussion names the species: a **hemithioacetal** (p. 3 and §4, p. 13).
3. **The thiol loss was instantaneous and then flat — an equilibrium, not a consumption.** "GSH
   reacted instantaneously with the three α-dicarbonyl compounds, and the scavenging of the
   α-dicarbonyls by GSH **did not continue** during the subsequent 6 h of incubation, during which
   only a small further decrease in GSH was observed **likely due to some autoxidation of GSH**"
   (§3.2, p. 7). The Discussion restates it: "the reactions reach **equilibrium** rapidly" (§4, p. 11).
   **Note the attribution of the only slow GSH loss: autoxidation — i.e. O₂, not the dicarbonyl.**
4. **The extent of thiol capture at equimolar loading is small and ordered MGO > GO > 3-DG:**
   ca. **18.1 % / 8.6 % / 1.7 %** of the GSH (§3.2, p. 7–8). At 1:1 and 0.5 mM, more than 80 % of the
   thiol was still free after six hours with methylglyoxal.

**What this does refute, stated narrowly.** In a **cell-free aqueous phosphate pot at 37 °C and
pH 7.4, aerobic, over 6 h, at 0.5 mM equimolar**, with **glutathione** and with **methylglyoxal,
glyoxal or 3-deoxyglucosone**, the α-dicarbonyl does **not** function primarily as a thiol oxidant.
The dominant fate of the thiol–dicarbonyl encounter is a reversible adduct. Both branches were
assayed on the same instrument in the same samples, so this is a genuine **within-study branch
comparison**, which is exactly the comparison the pre-registration's §3(b) needs and does not have.

**What this does NOT refute, stated just as narrowly, because the pre-registration's claim lives
somewhere else:**

- **Temperature.** 37 °C against the wave's **140 °C**. A 103 K gap. If the redox route has a higher
  barrier than the adduct route — which is the ordinary expectation for a two-electron transfer
  against a nucleophilic addition — then a 37 °C experiment is systematically blind to it, and the
  branch ratio at 140 °C could be inverted. **Zheng 2022 measures the branch ratio at the one
  temperature where the redox route is least likely to show.**
- **pH.** 7.4 against the Whitfield pot's **4.5**. GSH's thiol pKa is ≈8.8–9.2; at pH 7.4 a few per
  cent is thiolate, at pH 4.5 ~1e-3 of that. Both branches are thiolate-gated, but not necessarily
  equally.
- **The thiol.** Glutathione is a **peptidic aliphatic thiol**. The disulfide share the wave must
  reproduce belongs to **2-methyl-3-furanthiol** — a heteroaromatic thiol on an electron-rich furan
  with a substantially lower thiol pKa and a much lower one-electron oxidation potential. Furan- and
  thiophene-thiols dimerise readily; glutathione, in a cell, is held reduced by an entire enzymatic
  apparatus. **These are not the same oxidation problem.**
- **The α-dicarbonyl.** MGO, GO and 3-DG are α-**oxoaldehydes** — each has an aldehyde carbon, which
  is precisely the carbon that makes the hemithioacetal. The three species the pre-registration §3(a)
  proposes to source from norfuraneol are **2,3-pentanedione, 2,4-pentanedione and 3,4-hexanedione**:
  **alkyl diketones with no aldehyde carbon at all**. The adduct branch that dominated in Zheng's pot
  is structurally unavailable to them in the same form. **This is the single strongest reason the
  refutation does not transfer**, and it cuts *in favour* of the pre-registration's chemistry.
- **The redox mechanism proposed is not the one tested.** Whitfield's proposal, as the pre-reg
  reports it, is that the α-dicarbonyl is **reduced to a hydroxyalkanone** — a two-electron reduction
  of the dicarbonyl to an α-hydroxyketone. Zheng 2022 never looks for a hydroxyalkanone, never
  measures the dicarbonyl (§3.1, p. 6: "Due to the poor ionization efficiency and stability of
  dicarbonyls in the ESI source, **these compounds were not detected**"), and never reports a mass
  balance. It reports the *absence of GSSG*, which is a necessary consequence of the mechanism but
  not the mechanism itself.
- **Concentration and matrix.** 0.5 mM in 25 mM phosphate against a food pot at tens of mmol/L, in
  0.5 M phosphate, with H₂S and a full Maillard inventory.

**And the negative result is weaker than it reads, for one specific reason.** "Limited" is
**unquantified**. There is no percentage, no GSSG concentration, no limit of detection, no
calibration curve for GSSG (§2.5 says the calibration curves were made "for quantification of GSH",
and GSSG is listed only as monitored), and the supporting data are "data not shown". A reader cannot
tell whether "limited" means 0.1 % of the GSH or 3 %. **If it meant 3 % over six hours at 37 °C, that
would be entirely compatible with a redox channel that dominates at 140 °C.** The result is
directionally clear and quantitatively empty.

**The honest summary for the pre-registration's ledger.** Zheng 2022 is the first piece of evidence
in the corpus that measures the adduct/redox branch ratio for a thiol and an α-dicarbonyl in the same
pot, and it comes out against the redox branch. It is a **cell-free physiological-condition
experiment**, and the wave concerns 140 °C, pH 4.5, a heteroaromatic thiol and alkyl diketones with
no aldehyde carbon. It should move the pre-registration's odds — see §5, Flags 1 — and it should not
be treated as decisive.

**What this paper does NOT give the repository:** any rate constant; any reaction order; any
activation energy; any equilibrium constant (the equilibrium is described, not measured); any
temperature other than 37 °C; any pH other than 7.4; any measurement of the α-dicarbonyl itself; any
mass balance; any GSSG number; any volatile, headspace or aroma measurement; any food matrix.

## 2. Methods as they matter to a model

**The cell-free pot (§2.3, p. 3) — the only pot in this paper that a kinetic model can read.**

- **Reactants and concentrations.** GSH at a **final 0.5 mM** starting concentration; MGO, GO or 3-DG
  added from stock to a **final 0.5 mM (the 1:1 kinetic runs, Fig. 4) or 5 mM (the 10:1 runs used for
  adduct identification, Figs. 2 and 3)**. One α-dicarbonyl per incubation; no mixtures.
- **Buffer and pH.** **25 mM potassium phosphate, pH 7.4.** All five stock solutions (MGO, GO, 3-DG,
  GSH, GSSG) were made up in the same buffer, so there is no co-solvent.
- **Temperature.** **37 °C in a water bath.**
- **Time.** **0, 0.5, 1, 2, 4 and 6 h** (the 1:1 kinetic series); **1 h** for the 10:1 adduct
  identification.
- **Atmosphere.** **Not controlled and not stated.** A water bath, no degassing, no inert gas, no
  headspace statement. The authors themselves attribute the slow residual GSH loss to
  "**autoxidation of GSH**" (§3.2, p. 7), which is a direct admission that **O₂ was present and
  active**. This matters for the pre-registration in both directions: the pot had an ambient oxidant
  and *still* made little GSSG, but it also means the small GSSG that did form has an alternative
  father.
- **Quench.** At each time point **4 µL of acetic acid** was added to **196 µL** of sample — i.e. a
  **2 % v/v acid quench (mine)** — "to stabilize GSH and its related adducts", then immediate storage
  at **−80 °C** until analysis. **The quench acidifies the sample, which will shift any
  hemithioacetal equilibrium**; the paper does not report what that does to the measured adduct
  (Flags 6).
- **How the thiol was quantified.** **LC-TQ-MS** (Shimadzu Nexera XR LC-20AD XR UHPLC + Shimadzu 8040
  triple quadrupole, ESI positive, **MRM**). Column Phenomenex Luna Omega Polar C18, 100 × 2.1 mm,
  1.6 µm; mobile phase A = water + 0.1 % formic acid, B = acetonitrile + 0.1 % formic acid;
  0.2 mL/min; gradient 0–1 min 100 % A, 1–5 min 100→35 % A, 5–7.5 min 35 % A, 7.5–7.6 min 35→100 % A,
  7.6–18 min 100 % A. Nebulising gas 3.0 L/min, drying and heating gas 10.0 L/min, interface 300 °C,
  heat block 400 °C.
  **GSH transitions:** *m/z* 307.90 → 179.05 (CE −12 eV), → 76.10 (−25 eV), → 162.05 (−16 eV),
  → 84.05 (−21 eV). **Quantified against calibration curves made with the reference compound.**
  Readout: `Remaining GSH (%) = detected amount of GSH in test samples / original amount of GSH in
  the samples × 100`.
- **How the DISULFIDE was quantified — this is the load-bearing method statement.**
  **GSSG transitions:** *m/z* 613.15 → 355.05 (CE −22 eV), → 231.00 (−22 eV), → 484.15 (−22 eV),
  → 177.10 (−30 eV), on the same MRM run. **Four transitions, an authentic ≥98 % GSSG standard in the
  chemicals list — and no calibration curve stated for GSSG and no number reported.** §2.5 says
  "Quantification of GSH in the samples was achieved via calibration curves"; GSSG is not named in
  that sentence. The entire disulfide result is the phrase "was limited as confirmed by LC-TQ-MS
  (data not shown)".
- **How the adducts were identified.** **LC-TOF-MS** (Agilent 1200 LC + Bruker micro-TOF), same
  column, 0.18 mL/min, 1 µL injection, ESI positive, *m/z* 100–1500, capillary −4500 V, nebuliser
  1.2 bar, drying gas 8 L/min at 200 °C. **Qualitative only** — the paper says so ("used to
  qualitatively detect the reaction products"). **No adduct was ever quantified.**
- **What was NOT measured in this pot.** The α-dicarbonyls themselves ("these compounds were not
  detected", §3.1, p. 6). Any hydroxyalkanone or α-hydroxy acid. Any mass balance. Any dissolved
  oxygen. Any second temperature. Any second pH.
- **Replication.** Mean ± **SEM** of **three independent replicates** (Fig. 4 caption).

**The cell pots (§2.2, §2.6–§2.11) — none of it is kinetic, and all of it is at 37 °C in medium.**

- **Cells.** Nrf2 CALUX and Cytotox CALUX — human osteosarcoma **U2OS** stably transfected with a
  luciferase reporter under four EpREs, and under a constitutive promoter respectively (BioDetection
  Systems). DMEM/F12 GlutaMAX + 7.5 % FCS + 1 % NEAA + 0.3 % pen/strep, 5 % CO₂, 37 °C; 200 µg/mL
  G418 weekly.
- **Exposure.** 2×10⁴ cells/well in white opaque 96-well plates; assay medium is DMEM/F12 without
  phenol red + 5 % dextran-coated charcoal-stripped FCS; **eight concentrations 100, 250, 500, 750,
  1000, 1250, 1500 and 1750 µM**; **24 h**; compounds dissolved in sterile ultrapure water, added
  from 200× stocks; **curcumin 25 µM** as positive control. Readout: luciferase in RLU on a
  GloMax-Multi, expressed as induction factor against medium control.
- **Viability.** WST-1, 1 h, absorbance 440 nm against a 620 nm reference, expressed as % of medium
  control.
- **ROS.** DCFDA 25 µM in HBSS + 0.4 % FCS, 45 min loading, then **6 h** exposure; ex 485 / em 535 nm;
  **TBHP 50 µM** positive control.
- **GSH modulation.** **NAC 10 mM for 4 h** (then washed off — the paper is explicit that NAC was not
  co-exposed "since NAC is known to scavenge α-dicarbonyl compounds"), and **BSO 100 µM** for 24 h
  pre-incubation **plus** 100 µM BSO maintained through the 24 h exposure. §3.5 describes this as
  "100 µM BSO for **48 h**" (Flags 8).
- **Intracellular GSH.** LC-TQ-MS after a **200 µL 2 % TCA** quench, 15 min on ice, ≥6 h at −80 °C,
  scraped, vortexed, centrifuged 12 000 rpm (13 523 × *g*) for 30 min; **normalised to protein by BCA**.
- **Statistics.** Mean ± SEM, ≥3 independent experiments; independent-samples *t* test or Mann-Whitney
  U after Shapiro-Wilk; p < 0.05; SPSS 25.0; GraphPad Prism 9; ChemDraw 20.0.

## 3. Tables re-typed

**THIS PAPER PRINTS NO TABLES.** The string "Table" does not occur anywhere in the text layer.
Everything quantitative is either in running text or inside Figures 4–10. The figures are:
Fig. 1 (detoxification schemes), Fig. 2 (LC-TOF-MS base peak chromatograms, 10:1, 1 h),
Fig. 3 (LC-TOF-MS spectra of GSH and the four adducts), **Fig. 4 (GSH content vs time, 0–6 h — the
counter-evidence figure)**, Fig. 5 (WST-1 viability), Fig. 6 (Nrf2 and Cytotox CALUX induction),
Fig. 7 (DCFDA ROS), Fig. 8 (intracellular GSH after NAC/BSO), Fig. 9 (viability with NAC/BSO),
Fig. 10 (Nrf2 induction with NAC/BSO), Fig. 11 (adduct formation scheme).

**Per house rule, values read off figure axes are not typed.** Fig. 4's six-point curves are
therefore recorded as figure-only; only the three "ca." percentages the authors state in the running
text are usable.

### Every number printed in the running text

| quantity | value | class | where |
|---|---|---|---|
| **GSH reduction, equimolar MGO, immediate** | **ca. 18.1 %** | `[M]` | §3.2, p. 7–8 |
| **GSH reduction, equimolar GO, immediate** | **ca. 8.6 %** | `[M]` | §3.2, p. 8 |
| **GSH reduction, equimolar 3-DG, immediate** | **ca. 1.7 %** | `[M]` | §3.2, p. 8 |
| **GSSG formation from GSH over 6 h** | **"limited"** — no number, no bound, no LOD; "as confirmed by LC-TQ-MS (**data not shown**)" | `[M]` (asserted, unquantified) | **§3.2, p. 7** |
| further GSH decrease over 0.5–6 h | "only a small further decrease … **likely due to some autoxidation of GSH**" — no number | `[M]` (asserted, unquantified) | §3.2, p. 7 |
| adduct, MGO-GSH | *m/z* **380** [M+H]⁺, fragment **308** [M−72+H]⁺, retention **5.7 min** | `[M]` | §3.1, p. 6 |
| adduct, GO-GSH | *m/z* **366** [M+H]⁺, fragment **308** [M−58+H]⁺, retention **4.2 min** | `[M]` | §3.1, p. 6 |
| adduct, hydrous GO-GSH | *m/z* **384** [M+H]⁺, fragment **366** [M−18+H]⁺, retention **2.7 min** | `[M]` | §3.1, pp. 6–7 |
| adduct, 3-DG-GSH | *m/z* **470** [M+H]⁺, fragment **308** [M−162+H]⁺, retention **2.9 min** | `[M]` | §3.1, p. 7 |
| GSH parent ion | *m/z* **308** [M+H]⁺ | `[M]` | §3.1, p. 6 |
| cell viability, MGO 1750 µM, 24 h | **66.6 %** remaining | `[M]` | §3.3, p. 8 |
| cell viability, GO and 3-DG, 100–1750 µM, 24 h | no cytotoxic effect | `[M]` | §3.3, p. 8 |
| Nrf2 induction, first significant concentration | **≥750 µM MGO, ≥750 µM GO, ≥500 µM 3-DG** (p < 0.05) | `[M]` | §3.3, p. 9 |
| Nrf2 fold-induction order at ≥1250 µM | **MGO > 3-DG > GO** | `[M]` | §3.3, p. 9 and §4, p. 12 |
| ROS induction, first significant concentration (6 h, DCFDA) | **≥500 µM MGO, ≥1250 µM GO, ≥1250 µM 3-DG** (p < 0.05); order **MGO > GO ≈ 3-DG** | `[M]` | §3.4, p. 10 |
| intracellular GSH after 10 mM NAC, 4 h | **1.4-fold increase** | `[M]` | §3.5, p. 10 |
| intracellular GSH after 100 µM BSO | **18.4-fold decrease** | `[M]` | §3.5, p. 10 |
| cell viability, BSO alone | **89.5 %** remaining | `[M]` | §3.5, p. 10 |
| luciferase induction, BSO alone | **4.0-fold** vs untreated control | `[M]` | §3.5, p. 11 |
| Nrf2 induction factor, MGO 1000 / 1250 / 1500 µM, ± NAC | **3.8 → 1.9 (50 % decrease); 9.1 → 5.0 (45 %); 23.1 → 9.1 (61 %)** | `[M]` | §3.5, p. 11 |
| Nrf2 induction factor, GO 1500 µM, ± NAC | **2.9 → 2.1 (28 % decrease)** | `[M]` | §3.5, p. 11 |
| Nrf2 induction factor, MGO 750 / 1000 µM, ± BSO | **1.9 → 33.7-fold; 3.8 → 51.5-fold** | `[M]` | §3.5, p. 11 |
| Nrf2 induction factor, GO 750 / 1000 µM, ± BSO | **1.9 → 15.4-fold; 1.8 → 11.5-fold** | `[M]` | §3.5, p. 11 |
| Nrf2 induction factor, 3-DG 500 µM, ± BSO | **2.7 → 8.9-fold** | `[M]` | §3.5, p. 11 |
| typical cellular concentrations of MGO, GO, 3-DG | **1–4 µM** | `[C]` (refs 7,8) | §1, p. 1 |
| MGO speciation in water | monohydrate **71 %**, dihydrate **28 %**, unhydrated **1 %** | `[C]` (refs 37,38) | §4, p. 11 |
| GO speciation in water | mainly dihydrate; dimers **1–2 %**, monohydrate **0.5 %**, unhydrated **0.005 %** | `[C]` (refs 37,38) | §4, p. 11 |
| 3-DG reactivity with arginine residues | **200-fold lower** than MGO and GO | `[C]` (ref 42) | §4, p. 12 |
| BSO cytotoxic synergy | the response with BSO is "considerably higher than the sum of the responses induced by each α-dicarbonyl compound and BSO separately" — **more than additive**, no number | `[M]` | §3.5, p. 11 |

### Arithmetic on the printed numbers (all mine)

**1. The scavenging ladder as within-study ratios.**

| pair | ratio (mine) |
|---|---:|
| MGO / GO on GSH consumed at 1:1 | 18.1 / 8.6 = **2.10×** |
| MGO / 3-DG | 18.1 / 1.7 = **10.6×** |
| GO / 3-DG | 8.6 / 1.7 = **5.06×** |

These reproduce the paper's stated order MGO > GO > 3-DG and put numbers on the gaps. They are
**consumption extents at a fixed equimolar loading and a fixed time**, not rate ratios and not
affinity ratios (Flags 3).

**2. An apparent equilibrium constant, entirely conditional and NOT to be shipped (mine).** If the
plateau in Fig. 4 is the 1:1 adduct equilibrium `dicarbonyl + GSH ⇌ adduct` reached by 0.5 h and held
to 6 h — which is what the authors say it is ("the reactions reach equilibrium rapidly", §4, p. 11) —
then from 0.5 mM each and the stated fractional GSH losses:

| α-dicarbonyl | GSH consumed | [adduct] | [free] each | K_app = [adduct]/([D][GSH]) (mine) |
|---|---:|---:|---:|---:|
| MGO | 18.1 % | 0.0905 mM | 0.4095 mM | **≈ 5.4e2 M⁻¹** |
| GO | 8.6 % | 0.0430 mM | 0.4570 mM | **≈ 2.1e2 M⁻¹** |
| 3-DG | 1.7 % | 0.0085 mM | 0.4915 mM | **≈ 3.5e1 M⁻¹** |

**Four assumptions, none of them the paper's:** 1:1 stoichiometry; the plateau is equilibrium and not
a kinetic stall; *all* GSH loss is adduct and none is GSSG or autoxidation; and the free-dicarbonyl
concentration is the total, which the paper's own speciation discussion says is wrong (only **1 %** of
MGO and **0.005 %** of GO is unhydrated). Correcting for hydration would raise these by 10²–10⁴.
**These three numbers are `derived_assumption` and belong in no fit row.** They are recorded only so a
later reader does not re-derive them and mistake them for measurements.

**3. What "limited" would have to mean to be compatible with the pre-registration (mine, a
sensitivity note, not a claim).** Two GSH are consumed per GSSG. If GSSG accounted for, say, 1 % of
the starting GSH over 6 h, that is 0.005 mM GSSG from 0.01 mM of GSH — about **5.5 % of the 18.1 %
that MGO consumed**, and it would plausibly be described as "limited". At 0.1 % it would be 0.55 %.
**The paper's phrasing spans at least an order of magnitude in the quantity that matters**, and
nothing in the text narrows it.

**4. The quench (mine).** 4 µL acetic acid into 196 µL is **2.0 % v/v**; glacial acetic acid is
≈17.4 M, so the quenched sample is roughly **0.35 M acetic acid (mine)** — a pH well below 3. The
paper says this "stabilizes GSH and its related adducts"; it is also a large pH jump applied to a
system the paper has just described as a rapid equilibrium (Flags 6).

## 4. Numbers the repository can use

Every row below shares: **cell-free, 0.5 mM GSH, 25 mM potassium phosphate, pH 7.4, 37 °C water bath,
atmosphere uncontrolled, LC-TQ-MS MRM against a GSH calibration curve, mean ± SEM of three
independent replicates** — except the cell rows, which are U2OS Nrf2/Cytotox CALUX in DMEM/F12,
5 % CO₂, 37 °C, 24 h (6 h for ROS).

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **GSSG (disulfide) formation from GSH in the presence of MGO, GO or 3-DG** | **"limited"** — no numeric value, no bound, no LOD published | — | 0.5 mM GSH + 0.5 mM dicarbonyl, pH 7.4, 37 °C, **6 h**, aerobic | **§3.2, p. 7** ("data not shown") | **`measured_bound`, qualitative only** — the direction is measured, the magnitude is not published |
| GSH consumed at equimolar MGO, essentially at t = 0 and flat to 6 h | **ca. 18.1** | % of starting GSH | as above | §3.2, pp. 7–8 | **`measured_ratio`** |
| GSH consumed at equimolar GO | **ca. 8.6** | % of starting GSH | as above | §3.2, p. 8 | `measured_ratio` |
| GSH consumed at equimolar 3-DG | **ca. 1.7** | % of starting GSH | as above | §3.2, p. 8 | `measured_ratio` |
| MGO / GO / 3-DG scavenging ladder | **2.10 / 10.6 / 5.06** (MGO:GO, MGO:3-DG, GO:3-DG) | ratio | as above | 18.1, 8.6, 1.7 (mine) | **`within_study_ratio`** |
| time to reach the plateau | **≤ 0.5 h** (the first sampled point after zero; "instantaneously") and **no further scavenging to 6 h** | h | as above | §3.2, p. 7; Fig. 4 | **`measured_bound`** — an upper bound on the equilibration time, not a rate |
| any rate constant, order or activation energy | **— none, anywhere in this paper** | — | — | — | **absent** |
| apparent 1:1 adduct equilibrium constant, MGO / GO / 3-DG | **≈ 5.4e2 / 2.1e2 / 3.5e1** | M⁻¹ | as above; ignores hydration speciation, assumes 1:1 and assumes all GSH loss is adduct | §3 item 2 (mine) | **`derived_assumption` — DO NOT SHIP**; hydration correction alone moves these 10²–10⁴ |
| residual slow GSH loss over 0.5–6 h, and its cause as the authors assign it | "small", **attributed to autoxidation of GSH** — no number | % | as above, **aerobic** | §3.2, p. 7 | `level_only` — an admission that O₂ was present |
| cell viability, MGO 1750 µM, 24 h, U2OS | **66.6** | % of medium control | DMEM/F12, 5 % CO₂, 37 °C, WST-1 | §3.3, p. 8 | `level_only` |
| Nrf2 induction thresholds | **≥750 / ≥750 / ≥500** | µM for MGO / GO / 3-DG | 24 h, CALUX, p < 0.05 | §3.3, p. 9 | `level_only` |
| ROS induction thresholds | **≥500 / ≥1250 / ≥1250** | µM for MGO / GO / 3-DG | 6 h, DCFDA | §3.4, p. 10 | `level_only` |
| intracellular GSH modulation achieved | **1.4-fold up (10 mM NAC, 4 h) / 18.4-fold down (100 µM BSO)** | fold | U2OS, LC-TQ-MS, BCA-normalised | §3.5, p. 10 | `within_study_ratio` |
| 3-DG vs MGO/GO reactivity with arginine | **200-fold lower** | ratio | not measured here | §4, p. 12, ref [42] | **`level_only`, CITED** — do not treat as this paper's measurement |
| MGO / GO speciation in water | MGO 71 / 28 / 1 %; GO dihydrate-dominant, 1–2 % dimers, 0.5 % monohydrate, 0.005 % unhydrated | % | aqueous, not stated | §4, p. 11, refs [37,38] | **`level_only`, CITED** — but structurally important, see Flags 4 |

### Adduct branch or redox branch?

**This paper is the only one in the batch that assays BOTH branches in the same samples, and its
verdict is: adduct branch dominant, redox branch "limited".** The assignment for each observation:

| observation | branch | how it was decided |
|---|---|---|
| the ca. 18.1 / 8.6 / 1.7 % GSH losses | **ADDUCT** | The corresponding conjugates were found by LC-TOF-MS in the same incubations, each with the diagnostic neutral loss of one intact dicarbonyl (72, 58, 162 Da) leaving the GSH ion at *m/z* 308 (§3.1). The Discussion names the species **hemithioacetal** (§4, p. 13, and §1, p. 3). Sulfur–carbon, not sulfur–sulfur. |
| the "limited" GSSG | **REDOX — measured and found minor** | Directly instrumented: four GSSG MRM transitions (§2.5, p. 4), authentic GSSG standard (§2.1, p. 3). **This is the only direct redox-branch measurement in the batch.** Its magnitude is unpublished. |
| the small slow GSH decline over 0.5–6 h | **REDOX, but attributed to O₂ and not to the dicarbonyl** | The authors' own words: "likely due to some **autoxidation** of GSH" (§3.2, p. 7). If any of it were dicarbonyl-driven oxidation, this experiment could not tell the two apart — there is **no dicarbonyl-free GSH control curve reported** (Flags 5). |
| the flat plateau after 0.5 h | **ADDUCT, and reversible** | "the reactions reach equilibrium rapidly" (§4, p. 11). A rapidly-established equilibrium with >80 % of the thiol still free is the signature of a reversible hemithioacetal, which is what the glyoxalase system is built to intercept. |
| everything in the cell assays | **neither** | Reporter-gene induction factors, WST-1 viability and DCFDA fluorescence. No branch is measured; ROS induction is a downstream cellular response, not a thiol-oxidation rate. |

**The one structural caveat that most limits the transfer of this verdict**, stated here because it
belongs with the branch assignment: **the adduct branch won because MGO, GO and 3-DG all carry an
aldehyde carbon**, and the hemithioacetal forms there. The pre-registration's proposed α-dicarbonyls
(2,3-pentanedione, 2,4-pentanedione, 3,4-hexanedione, from Whitfield's Table 1 rows 1–3) are
**symmetric or near-symmetric alkyl diketones with no aldehyde carbon**. Removing the aldehyde
removes the winning branch. Zheng 2022 therefore refutes the redox hypothesis for α-**oxoaldehydes**
and says considerably less about α-**diketones** — which is the class the pre-registration's §3(a)
source step would actually make.

## 5. Flags

1. **THE ODDS QUESTION, ANSWERED DIRECTLY.** Should this change the pre-registration's stated
   probabilities? **Yes, but modestly, and only two of the four.** Prediction 1 (T3 improves by >1
   decade on ≥2 of Zhou's shares, **60 %**) and prediction 2 (T6 holds, **65 %**) are about whether a
   channel that is currently structurally dead can carry gradient once it is fed — they are claims
   about the *objective's* geometry, not about the chemistry's direction, and Zheng 2022 touches
   neither. Prediction 3 (the wave does **not** ship, **50 %**) is the one that should move: the
   `ch_redox_dicarbonyl` constant is fitted with no anchor, and this paper is the first evidence in
   the corpus that the branch it represents is minor where it has been looked for. **A fitted
   constant with no anchor and now one contrary observation is a weaker structure than a fitted
   constant with no anchor and no evidence either way.** My own reading, offered as a reading and not
   as a rule: prediction 3 should go up by a few points, not by twenty — because the temperature, the
   pH, the thiol class and above all the **absence of an aldehyde carbon in the pot's own diketones**
   all cut the other way, and because "limited" is unquantified. **What the pre-registration should
   actually do is add a row to §5: report the fitted `ch_redox_dicarbonyl` constant against an
   explicit prior that says the redox branch is minor at 37 °C, and declare in advance what fitted
   value would be considered implausible.** That is a change to the record, not to the odds, and it is
   the more useful one.
2. **TEMPERATURE TRANSFER, IN BOTH DIRECTIONS, AND IT COSTS MORE HERE THAN ANYWHERE ELSE.** 37 °C
   against **140 °C**. Every other flag in this dossier is downstream of this one. The particular
   damage is that **a null result at low temperature is the weakest kind of null**: if the redox route
   has the higher barrier — the ordinary expectation for a two-electron transfer competing with a
   nucleophilic addition — then a 37 °C pot is precisely where it is least visible, and a 103 K rise
   would favour it by whatever the barrier difference implies. With no second temperature anywhere in
   this paper there is **no way to estimate that difference from this evidence**. The refutation is
   therefore strong about 37 °C and structurally uninformative about 140 °C. **Anyone who cites this
   paper as closing the question is over-reading it, and anyone who dismisses it because of the
   temperature is under-reading it** — it remains the only direct branch measurement in the corpus.
3. **"Limited" is never quantified and the data are not shown.** No percentage, no concentration, no
   LOD, no LOQ, no GSSG calibration curve, no figure. Four MRM transitions were monitored and the
   result is one adjective. **This is the single biggest weakness of the counter-evidence** and it is
   the first thing to request from the authors (Flags 11).
4. **The α-dicarbonyls are α-oxoaldehydes and the pot's are alkyl diketones.** MGO, GO and 3-DG each
   have an aldehyde carbon. 2,3-pentanedione, 2,4-pentanedione and 3,4-hexanedione do not. The
   winning branch in Zheng's pot is the one that needs that carbon. Additionally the paper's own
   cited speciation says only **1 % of MGO and 0.005 % of GO** is unhydrated in water — an alkyl
   diketone is far less hydrated, so its *reactive* fraction is orders of magnitude larger and the
   whole balance could sit elsewhere.
5. **There is no dicarbonyl-free GSH control curve in the reported data.** The residual slow GSH loss
   is attributed to autoxidation, but no blank incubation of GSH alone over 0–6 h is shown, so the
   attribution is an inference. Without that blank, the paper cannot separate "the dicarbonyl oxidises
   the thiol slowly" from "O₂ oxidises the thiol slowly" — and both would appear as a small,
   unquantified GSSG signal.
6. **The acid quench perturbs the very equilibrium being measured.** 2 % v/v acetic acid (≈0.35 M,
   mine) is added to a system the paper describes as a rapid, reversible equilibrium. Acid stabilises
   the thiol against autoxidation, which is why it is there, but it also shifts hemithioacetal
   formation and could partially reverse it in the seconds before freezing. **The measured
   "remaining GSH" is a post-quench quantity.** No recovery control is reported.
7. **Only GSH is quantified; the adducts and the dicarbonyls are not.** The LC-TOF-MS work is
   explicitly qualitative. So there is **no mass balance**: the 18.1 % of GSH that disappeared is
   never accounted for as adduct + GSSG + other. If some of it went somewhere unmeasured, nothing in
   this paper would catch it.
8. **A methods inconsistency on the BSO exposure.** §2.10 (p. 6) describes 24 h pre-incubation with
   100 µM BSO followed by 24 h co-exposure; §3.5 (p. 10) says "treated with 100 µM BSO for **48 h**".
   These are reconcilable (24 + 24 = 48) but the 18.4-fold GSH depletion is quoted against the 48 h
   figure while Fig. 8's legend describes exposure "to NAC or BSO only". Minor, but the depletion
   factor's exposure duration should be confirmed before it is quoted.
9. **The thiol is wrong for the target pot, in a way that matters mechanistically.** Glutathione is a
   tripeptidic aliphatic thiol maintained reduced *in vivo* by glutathione reductase; its
   autoxidation and its dimerisation are both slow by design. The repository's disulfide share
   belongs to **2-methyl-3-furanthiol** and **2-furfurylthiol**, heteroaromatic thiols whose radical
   and thiolate chemistry is quite different and which are known to dimerise readily in food systems.
   `parameters_sulfur.py` already treats those two as a separate class with their own barrier
   (`ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL`). **A null for GSH is not a null for MFT.**
10. **The concentrations are ~100× below the pot's, and the matrix is a 25 mM buffer.** No protein,
    no sugar, no amino acid other than the GSH itself, no H₂S, no melanoidin, no lipid. The pots the
    objective must reproduce carry all of these, several of which are themselves redox-active.
11. **What to request from the authors:** (i) **the GSSG numbers behind "limited"** — the percentage
    of starting GSH found as GSSG at each of 0, 0.5, 1, 2, 4 and 6 h for each dicarbonyl, plus the
    GSSG LOD and whether a GSSG calibration curve exists; (ii) the **GSH-alone blank** over the same
    six hours; (iii) the numeric values behind Fig. 4 (the paper's Data Availability Statement offers
    them on request); (iv) whether any incubation was ever run above 37 °C; (v) whether an
    α-hydroxyketone or α-hydroxy acid was ever searched for in the LC-TOF-MS traces.
12. **What this paper does NOT contain:** any table; any rate constant; any order; any activation
    energy; any equilibrium constant; any second temperature; any second pH; any measurement of the
    α-dicarbonyl itself; any mass balance; any GSSG quantity; any α-diketone (as opposed to
    α-oxoaldehyde); any heteroaromatic thiol; any food matrix; any volatile measurement; any
    supplementary material; any DFT or other computation (so nothing was excluded under the no-DFT
    policy).
