# Rizzi, Amba & Heineman 2010 — EXTRACTION (aldoses 0.10 M in 0.10 M K-phosphate pH 7.1, reflux ~100 °C, 10-300 min, NO amine: reducing species by 2,6-dichloroindophenol (DCIP) titration in ascorbic-acid equivalents; open-circuit Pt/AgCl redox potential after 0-4 h; one OPD trap of xylose at 2 h; cyclic voltammetry of an Amadori compound)

### PHOSPHATE ALONE (pH 7.1, 100 °C, no amine) TURNS PENTOSES INTO A **mM-SCALE POOL OF 2-ELECTRON REDUCTANTS**: after 300 min, **ribose 40.5, arabinose 22.7, xylose 15.3, glucose 3.3 µmol AAE per mmol sugar** (i.e. 4.1 / 2.3 / 1.5 / 0.33 mM from 0.1 M sugar), rising ~linearly after a 20-60 min lag at **≈ 15 / 9 / 6 / 1.3 µM min⁻¹**. Phosphate is required (earlier work); unheated sugars and an Amadori compound give zero DCIP titre. The one identified reductone (tentative) is **4,5-dihydroxy-2,3-pentanedione = 1-deoxy-2,3-pentosone = the engine's DPO**. This is the reductant class Kobayasi 1965's MAA → H2S + acetaldehyde (and furfural → furfuryl alcohol) needs, but the paper gives **no reduction rate** for any acceptor.

**Source on disk:** `data/articles/rizzi2010.pdf` (5 pp., ACS publisher PDF, clean text layer; download stamp 28 Oct
2026). Prose from `pdftotext -layout`; Figs. 2-4 digitised from 220-dpi renders (`scratch/lit8/riz_f23.png`,
`riz_f3.png`, `riz_f4.png`), tagged [D].

## 0. Identity

| field | value |
|---|---|
| Title | "Quantification of Chemically Reducing Species in the Phosphate Ion Catalyzed Degradation of Reducing Sugars" |
| Authors | **George P. Rizzi** (corresponding; Miami University Middletown), Eme E. Amba, William R. Heineman (University of Cincinnati) |
| Venue | J. Agric. Food Chem. **58** (17), 9739-9743, 2010. Received 25 Mar 2010, revised 2 Aug, accepted 4 Aug, web 16 Aug 2010 |
| DOI as printed | **10.1021/jf102433u** ✔ matches the brief. **Identity verified.** |

## 1. Why it matters

`kobayasi1965_extraction.md`: mercaptoacetaldehyde (MAA, the Strecker sulfur co-product of cysteine) → H2S +
acetaldehyde is a two-electron reduction and needs a reductant (enaminols / reductones). SULFUR_FLUX_MAP sec. 5's pH-7
arm is H2S-starved because the engine sends Strecker sulfur to an inert pool. The question here: how big and how fast is
a sugar-derived reductant pool under Maillard-like conditions, which sugars, which pH, which molecules — and does it bear
on furfural → furfuryl alcohol (Schirle-Keller 1992 finds furfuryl alcohol ≈ furfural in a Glc + Cys pot)?

## 2. Methods as they matter to a model (p. 9740)

- **Pots.** Sugar 0.10 M in 0.10 M potassium phosphate, pH 7.1 (KH2PO4 + K2HPO4), 30.0 mL, heated at reflux ("ca.
  100 °C"), open to air through the condenser (vessel otherwise not described). Sugars: ribose, arabinose, xylose, glucose
  (titration); + galactose, mannose, xylitol, sorbitol (EMF). **No amine** in any titration/EMF pot.
- **pH drift.** "Reaction pH remained nearly constant during 300 min, occasionally falling to as low as 6.8."
- **DCIP titration.** 1.00 mL aliquots, cooled < 1 min to ~23 °C, titrated under propane with 1 mM DCIP in the same buffer
  (standardised vs ascorbic acid; stoichiometry assumed equal to ascorbate's, i.e. **1 AAE = one 2-electron reductant
  equivalent**). Two independent experiments averaged; precision ±15 %. Titration beyond 300 min impossible (colour).
- **EMF.** Pt/AgCl combination electrode at 23 °C on cooled solutions; reported as Ef − Eb (mV, vs unheated buffer).
  Xylose/xylitol triplicate; others single (±15 %).
- **OPD trap.** 2 × ~465 mg xylose in 30 mL buffer, 2 h reflux; OPD 50 mg, 30 min, 50 °C; CH2Cl2 extraction → 10.5 mg
  residue → ~2 mg polar quinoxalines (λ 237, 318 nm); major, Rf 0.50, ESI-MS MH⁺ 205 → tentatively
  2-(1,2-dihydroxyethyl)-3-methylquinoxaline = OPD derivative of **4,5-dihydroxy-2,3-pentanedione (1-deoxy-2,3-xylosone)**.
- **CV.** 2 mM N-(1-deoxyfructos-1-yl)piperidine, 0.10 M phosphate pH 7.14, glassy carbon, ±1.5 V, 50 mV/s, in air.

## 3. Numbers

### Fig. 2 (p. 9741), reducing species, µmol AAE per mmol initial sugar, 0.10 M sugar, 0.10 M phosphate pH 7.1, ~100 °C [D, ±0.5]

**The x-axis is categorical (10, 20, 40, 60, 180, 300 min equally spaced)**: the visual "acceleration" in the figure is
largely an axis artefact.

| sugar | 10 | 20 | 40 | 60 | 180 | 300 min |
|---|---|---|---|---|---|---|
| ribose | 1.2 | 1.2 | 3.0 | 5.4 | 24.2 | **40.5** |
| arabinose | 0.6 | 0.6 | 1.2 | 2.0 | 12.9 | **22.7** |
| xylose | 0.5 | 0.5 | 0.9 | 1.4 | 9.0 | **15.3** |
| glucose | 0.2 | 0.2 | 0.2 | 0.2 | 1.8 | **3.3** |

### Fig. 3 (p. 9741), Ef − Eb (mV), 0.10 M in 0.10 M phosphate pH 7.1, 100 °C, mean ± 1 SD of 3 [D, ±5 mV] (x-axis categorical: 0, 1, 2, 4 h)

| | 0 h | 1 h | 2 h | 4 h |
|---|---|---|---|---|
| xylose | −42 ± 19 (text) | −272 ± 12 | −315 ± 14 | −346 ± 15 (text) |
| xylitol | −52 | −30 | −32 | −52 ± 19 (text) |

### Fig. 4 (p. 9742), Ef − Eb after 4 h, 0.10 M in 0.10 M phosphate pH 7.2, 100 °C, single runs (±15 %) [D, ±5 mV]

sorbitol −18; **glucose −242** (text); galactose −265; mannose −276; xylitol −44; arabinose −325; ribose −341; xylose −346.

### Arithmetic (mine)

1. **Concentrations.** ×0.1 M sugar: at 300 min ribose **4.05 mM AAE**, arabinose 2.27, xylose 1.53, glucose 0.33 mM.
2. **Production rates (60→300 min, ~linear):** ribose 0.146, arabinose 0.086, xylose 0.058, glucose 0.013 µmol AAE
   mmol⁻¹ min⁻¹ → **14.6 / 8.6 / 5.8 / 1.3 µM min⁻¹**, i.e. pseudo-first-order **1.5e-4 / 8.6e-5 / 5.8e-5 / 1.3e-5 min⁻¹**
   per sugar, as a *net* (formation − consumption) of titratable reductant. Ribose 60→180: 0.157; 180→300: 0.136 (slight
   slowing; "seemed to level off at 300 min").
3. **Lag:** ≤ 1.2 µmol/mmol for the first 20 min for all sugars — consistent with a consecutive scheme (sugar →
   enediol/deoxyosone → reductone).
4. **Ranking conflict:** DCIP ranks ribose > arabinose > xylose > glucose (Fig. 2), EMF ranks xylose ≥ ribose > arabinose ≫
   hexoses (Fig. 4). EMF is logarithmic in the redox couple and single-run; trust the titration for amounts.
5. **OPD yield (lower bound):** ~2 mg quinoxalines (M ≈ 204) ≈ 10 µmol from ~6.2 mmol xylose after 2 h ≈ **0.16 mol %**,
   after losses in extraction and repeated TLC.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| ribose → reductant pool (net) | rate | 14.6 (60-300 min) | µM AAE min⁻¹ | 0.1 M ribose, 0.1 M K-phosphate pH 7.1→≥6.8, ~100 °C, no amine | Fig. 2 [D] | derived from digitised |
| arabinose / xylose / glucose, same | rate | 8.6 / 5.8 / 1.3 | µM AAE min⁻¹ | same | Fig. 2 [D] | derived |
| reductant pool at 300 min | level | 4.05 / 2.27 / 1.53 / 0.33 | mM AAE | Rib / Ara / Xyl / Glc | Fig. 2 [D] | digitised |
| reductant identity (major, tentative) | — | 4,5-dihydroxy-2,3-pentanedione (= DPO) | — | xylose, 2 h | OPD + ESI-MS | tentative ID |
| phosphate requirement | — | required (ref. 10 = Rizzi 2004) | — | pH 7.1, 100 °C | text | cited |
| Hofmann 1999 comparison (cited) | yield | 60 µmol AAE/mmol (6.0 %) at 50 min | — | 1.5 M Glc + 1.5 M Ala, 0.5 M phosphate pH 7, reflux | text p. 9741 | secondary |
| reduction rates of any acceptor (MAA, furfural) | — | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** the size and growth rate of an amine-free, phosphate-catalysed reductant pool from pentoses at pH 7, 100 °C
(order: mM after hours; ribose > arabinose > xylose ≫ glucose); that the reductant is (at least partly) the 1-deoxy-2,3-
osone itself; that with an amine (Hofmann 1999, cited) the pool is larger and faster (6 % in 50 min at 1.5 M).
**Cannot:** no rate for reduction of MAA, furfural or any carbonyl; no temperature dependence (100 °C only); no pH
dependence (pH 7.1 only); no amine or thiol present (DCIP would also titrate thiols, so the method cannot be transferred
as-is to a cysteine pot); open-air reflux (O2 consumes reductant; the titre is net of autoxidation).

## 6. Flags

1. **"The maximum chemical yield of reducing species observed was 1.7 %, with ribose after 300 min"** (p. 9741) contradicts
   Fig. 2 (40.5 µmol/mmol = **4.05 %**). Possibly a different stoichiometric convention; unexplained. The Hofmann
   comparison in the same paragraph uses µmol/mmol × 0.1 = % (60 → 6.0 %), which supports 4.05 %.
2. Categorical time axes in Figs. 2 and 3 (see above).
3. Open-chain percentages quoted in the text (ribose 0.05 %, arabinose 0.03 %, glucose 0.002 %; Angyal 1984) differ from
   Dworkin & Miller 2000's kinetic values (ribose 0.040 %, arabinose 0.019 %, glucose 0.0014 % at 25 °C).
4. The CV section is qualitative (two anodic waves at +0.7 and +0.9 V vs Ag/AgCl for the piperidine Amadori compound; no
   cathodic return wave); no numbers for the model.
5. Propane as "inert" gas for titration; reflux under air for the reactions.

## 7. Use in the model

- **Reductant for MAA → H2S + acetaldehyde (Kobayasi 1965).** If B59/B-next writes that step as reductant-dependent, this
  paper supports using the **deoxyosone pool itself (DPO, and by analogy TDP)** as the reductant: the tentatively
  identified reductone is 4,5-dihydroxy-2,3-pentanedione, i.e. the engine's DPO. Scale check at pH 7, 100 °C, no amine:
  ~15 µM/min of 2-e⁻ equivalents from 0.1 M ribose (scale linearly with sugar concentration if first order, which is
  not tested here); more is expected at 145 °C, with 0.5 M phosphate and with an amine (all unmeasured in this paper).
  That is a **reductant supply ceiling, not a rate constant**; the MAA reduction constant remains unmeasured. Note the
  coupling: DPO used as reductant is DPO not available to NF (and the engine's pH-7 arm already has DPO in excess:
  SULFUR_FLUX_MAP sec. 5, 49 mM DPO, 13.5 mM NF unused).
- **Furfural → furfuryl alcohol.** Same pool; Schirle-Keller 1992 shows furfuryl alcohol ≈ furfural (by area) in Glc + Cys
  at 120-150 °C, so a reductive furfural sink is real and of order the furfural standing pool. This paper gives no rate
  for it; it only shows a pentose-derived reductant exists at mM level. A furfuryl-alcohol channel would be part of the
  "large unidentified" `k_fur_decay` sink rather than a new barrier datum.
- **Sugar ranking for reductant supply** (Rib 1 : Ara 0.56 : Xyl 0.38 : Glc 0.08 at 300 min) can be a within-study
  directional claim if the engine ever carries a reductant pool. Not a hold-out for any current engine observable.
