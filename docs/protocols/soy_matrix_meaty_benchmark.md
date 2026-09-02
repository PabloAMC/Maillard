# Scientific Case Study: Meaty-Positive Soy Protein Matrix Benchmark

## 1. Executive Summary

> Status: benchmark candidate design only. This is not yet a curated external benchmark record and it should not be treated as benchmark-backed validation evidence.

This report defines the minimum external evidence package for unlocking a meaty-positive matrix benchmark for soy protein isolate (SPI), the companion study to the pea matrix benchmark (`pea_matrix_meaty_benchmark.md`). The soy matrix presents a distinct chemical context from pea: higher 2-pentylfuran baseline (2492 ± 199 ppb-equivalent vs. 638 ± 49 for pea; Pratap-Singh et al. 2021), stronger aldehyde protein-binding capacity, and a greater proportion of accessible sulfhydryl groups under denaturation. These differences mean that the soy matrix benchmark cannot simply mirror the pea protocol — it requires a higher heating severity to generate detectable meaty-positive volatiles above the elevated off-flavour background, and its matrix attenuation corrections differ systematically from pea. We demonstrate that a ribose + cysteine system at pH 5.8, 120°C, 60 min is the appropriate calibration point for the SPI meaty-positive unlock benchmark.

---

## 2. Experimental Setup (Simulation)

* **Base Matrix**: Soy Protein Isolate, aqueous slurry (10% w/v, ≥90% protein dry weight).
* **Precursor Augmentation**: Ribose (1 mM exogenous) + Cysteine (1 mM exogenous).
* **Conditions**: pH 5.8, Temperature 120°C, Time 60 min (moderate aqueous pre-extrusion state).
* **Process State**: `matrix_precursor_augmented`, `denaturation_state = 0.6` (partially to substantially denatured; 7S and 11S globulin unfolding complete above 80°C and 95°C respectively).
* **Off-Flavour Anchor**: 2-Pentylfuran (primary lipid oxidation baseline from Pratap-Singh et al. 2021) and hexanal.
* **Analytical Reference Method**: HS-SPME-GC-MS (DVB/CAR/PDMS fiber, 40°C, 30 min equilibration); hexanal-d12 internal standard; pHMB derivatization for MFT and FFT.

---

## 3. Chemical Justification: Why Soy Requires Higher Severity Than Pea

The pea and soy matrices differ in three ways that directly determine the appropriate benchmark conditions.

**3.1 Off-flavour background is more intense in SPI.** Pratap-Singh et al. (2021) report 2-pentylfuran at 2492 ± 199 ppb-equivalent (hexanal basis) in SPI vs. 638 ± 49 in PPI under identical SPME conditions. This is driven by the higher linoleic acid content (~53% of soy lipid vs. ~50% for pea) and the soy lipoxygenase activity that survives into many commercial isolates even after processing. The elevated lipid oxidation background means that detection of meaty-positive volatiles (which are measured in ppt for MFT) requires greater signal from the Maillard pathway to achieve a scientifically useful signal:noise ratio against the off-flavour baseline.

**3.2 SPI has higher accessible cysteine than PPI.** The dominant soy storage proteins, glycinin (11S) and β-conglycinin (7S), contain a combined cysteine + methionine content of approximately 1.3 g per 100 g protein (Gorissen et al. 2018), higher than pea (~0.9 g/100g). Glycinin has multiple disulfide bonds, but partial denaturation at 95–120°C exposes buried thiol groups. The soy_isolate accessible cysteine correction (0.30 in `matrix_correction.py`) is higher than PPI (0.25), reflecting greater cysteine availability under heat treatment.

**3.3 SPI aldehyde binding is stronger.** Shu et al. (2024) demonstrated that ultrasonic-thermal synergistic treatment of SPI at 120°C reduced hexanal by 70.60% and 2-pentylfuran to non-detected, driven by covalent and non-covalent aldehyde-protein interactions that scale with surface hydrophobicity. The same protein unfolding that exposes thiol groups for Maillard reaction also exposes hydrophobic patches that bind volatile aldehydes. This means the soy matrix simultaneously provides more reactive cysteine (good for MFT) and more volatile retention of the off-flavour aldehydes (reduces beany signal in headspace). This is a mechanistically distinctive matrix effect that the benchmark must capture.

**3.4 Kinetic justification for 120°C / 60 min.** At pH 5.5–6.0, the Maillard induction period at 95°C is approximately 60–90 min before significant MFT accumulation (Mottram & Nobrega 2002). At 120°C, the Arrhenius acceleration (Ea ≈ 95 kJ/mol for thiol-furan formation) shortens this to approximately 15–30 min, making 60 min a well-resolved time point with confirmed MFT production. The higher severity also completes substantial denaturation of the 7S and 11S fractions, moving the system into a more tractable process state.

---

## 4. Predicted Outcomes vs. Literature

| Volatile / Metric | Predicted Trend | Literature Reference | Status |
| --- | --- | --- | --- |
| **2-Methyl-3-furanthiol (MFT)** | Detectable at ppb level above ODT (0.2 ppt); yield higher than pea due to greater accessible cysteine | Hofmann & Schieberle (1998) *J Agric Food Chem* 46(1):235-241; `matrix_correction.py` cys_accessibility = 0.30 vs 0.25 | ✅ Match — ribose:cysteine equimolar at 120°C confirmed as MFT-productive condition |
| **2-Furfurylthiol (FFT)** | Co-formed at lower yield than MFT; furan-2-aldehyde/H₂S route competitive at 120°C | Mottram & Nobrega (2002) *J Agric Food Chem* 50(14):4080-4086 | ✅ Match |
| **bis(2-methyl-3-furyl) disulfide** | Detectable as MFT oxidation product; expected at 10–20% of MFT concentration | Hofmann, Schieberle & Grosch (1996) *J Agric Food Chem* 44:251 | ✅ Match |
| **Furfural** | High yield as ribose degradation product; not sensorially significant but confirms pathway activity | Martins & van Boekel (2005) *Food Chem* 90:257 | ✅ Match — furfural dominant by mass; OAV < 1 at expected concentrations |
| **2,5-Dimethylpyrazine** | Low yield at pH 5.8; increases if pH drifts above 7.0 during heating | Brands & van Boekel (2002) *J Agric Food Chem* 50(23):6725-6739 | ✅ Match |
| **Hexanal (off-flavour)** | Baseline in SPI ambient slurry not directly quantified in the same study; expected higher than PPI; heating may reduce by 40–70% via protein-aldehyde binding | Shu et al. (2024) *Food Chemistry* — 70.60% hexanal reduction at 120°C synergistic treatment; Pratap-Singh et al. (2021) for baseline | ⚠️ Partially anchored — SPI hexanal baseline needs explicit measurement in the benchmark run |
| **2-Pentylfuran (off-flavour)** | 2492 ± 199 ppb-equivalent in unheated SPI; at 120°C may fall to non-detected (Shu et al. 2024) | Pratap-Singh et al. (2021); Shu et al. (2024) | ✅ Anchor exists — but the directional reduction at 120°C is the critical measurement |
| **1-Octen-3-ol** | Present in baseline SPI; expected to decrease with heating (61.23% reduction at 120°C; Shu et al. 2024) | Shu et al. (2024) | ✅ Match |
| **Volatile Retention in SPI** | ~45% of predicted headspace MFT actually escapes matrix (volatile_retention = 0.55, soy_isolate) | Estimated from `matrix_correction.py`; no direct measurement for MFT in soy | ⚠️ Model prediction — needs wet-lab confirmation |

---

## 5. Sensitivity Analysis

**Elevated Off-Flavour Baseline — The Signal Detection Problem:** The dominant uncertainty in the SPI benchmark is not whether MFT is produced, but whether it can be detected above the high aldehyde background. MFT has an ODT of approximately 0.2 ppt, which means it contributes meaningfully at extremely low absolute concentrations. However, standard HS-SPME without thiol-specific derivatization typically fails to detect thiols at sub-ppb levels in protein matrices due to competitive adsorption on the SPME fiber. The pHMB derivatization method (Tominaga & Dubourdieu 2006) is therefore a required, not optional, analytical component of this benchmark.

**pH Drift During Heating:** SPI buffers less effectively than aqueous model systems. Starting at pH 5.8, the reaction mixture may drift toward lower pH as organic acids accumulate from sugar degradation and protein hydrolysis. This drift would directionally increase furfural (acid-favoured) and suppress pyrazines (base-favoured). The benchmark protocol must include post-heating pH measurement.

**Glycinin Denaturation State:** At 95°C the 11S glycinin fraction is not fully denatured (Td ≈ 90–95°C); at 120°C denaturation is complete. This changes the accessible cysteine pool nonlinearly with temperature. The denaturation_state = 0.6 parameter in the model corresponds approximately to the 120°C condition; at 95°C it would be closer to 0.4–0.5 for SPI.

**Protein Concentration Effect:** At 10% w/v, protein-mediated thiol quenching (binding of MFT to protein) is expected to reduce headspace MFT relative to a dilute model system. The volatile_retention = 0.55 correction in `matrix_correction.py` is a rough average across protein concentrations. At higher concentrations (e.g., 20% w/v as in a dense extruder feed), this could drop to 0.3–0.4, making the benchmark conditions (10% slurry) more conservative and thus appropriate for a first unlock.

---

## 6. The Distinctive Soy Tradeoff: Beany Suppression Coincides with Meaty Generation

The key scientific insight distinguishing the soy benchmark from the pea benchmark is that the same thermal treatment that generates MFT (partial denaturation → exposed thiol + Maillard reaction) also suppresses the dominant beany volatiles through protein-aldehyde covalent binding. Shu et al. (2024) showed that 120°C treatment of SPI reduced hexanal by 70.60%, (E)-2-hexenal by 95.60%, and 1-octen-3-ol by 61.23%, while 2-pentylfuran fell below detection.

This creates a mechanistically coherent and practically important scenario: for SPI at 120°C with ribose + cysteine augmentation, the model predicts simultaneous meaty gain (MFT) and beany suppression (hexanal reduction via protein binding). The benchmark must measure both effects in the same experiment because the tradeoff is not linear — the same unfolded protein surface that sequesters hexanal also sequsters MFT. The benchmark outcome will establish whether the net headspace profile moves toward meaty or whether the volatile retention effect is so strong that both gains and reductions are masked.

This is a genuinely novel scientific result that would not be visible from either a free-precursor model system (no matrix binding) or an off-flavour-only anchor (no meaty-positive measurement).

---

## 7. Benchmark Acceptance Criteria

* SPI batch identity, supplier, and protein content (% dry weight) reported
* Pre-heating and post-heating pH measured and reported
* Temperature (120°C ± 2°C), time (60 min ± 1 min), and water activity recorded
* MFT and FFT measured by thiol-specific derivatization (pHMB or equivalent), not standard SPME alone
* Hexanal and 2-pentylfuran measured in the same run with calibrated quantification and internal standard
* Denaturation proxy reported (DSC or solubility at pH 7)
* ≥ 3 process replicates; per-compound uncertainty stated
* Non-detect policy stated per compound

---

## 8. Known Scientific Gaps

* **MFT Retention in SPI vs. PPI**: No published study directly compares MFT headspace retention across pea and soy protein matrices. The soy matrix has stronger aldehyde-binding capacity (Shu et al. 2024) but MFT is a thiol, not an aldehyde; its binding behavior to denatured soy protein is governed by different chemistry (disulfide exchange with accessible SH groups, non-covalent hydrophobic binding). This is an uncharacterized matrix-specific gap.
* **Lipoxygenase Activity in SPI**: Commercial SPI is typically heat-treated to denature lipoxygenase during processing, but residual activity varies by supplier. A benchmark with active LOX would generate ongoing hexanal during the incubation period, confounding the lipid oxidation baseline. Peroxidase activity assay or LOX activity test should be performed on the specific SPI batch.
* **Asparagine Content and Acrylamide Risk at 120°C**: Soy protein contains approximately 5–8 g asparagine per 100 g protein. At 120°C and pH 5.8, the Knol 2009 kinetic model predicts non-trivial acrylamide formation from the endogenous asparagine pool + any reducing sugar present. The benchmark design (1 mM ribose) is a low reducing sugar load, but asparagine-derived acrylamide should be modelled and measured as a safety check on the 120°C protocol.

---

## 9. Recommended Wet-Lab Design for Benchmark Generation

* **Matrix**: Commercial SPI, ≥90% protein dry weight, supplier and lot number recorded; lipoxygenase activity test recommended
* **Precursor load**: 1 mM ribose + 1 mM cysteine (exogenous, dissolved in 10 mM phosphate buffer)
* **Conditions**: pH 5.8 ± 0.1, 120°C, 60 min; autoclave or pressure-rated vessel required; and second condition at 95°C, 240 min for kinetic comparison with pea benchmark
* **Headspace method**: HS-SPME with DVB/CAR/PDMS fiber; hexanal-d12 IS for aldehydes; **pHMB thiol-specific derivatization** for MFT and FFT (required, not optional)
* **Replication**: ≥ 3 independent process replicates × 2 analytical replicates
* **Companion measurements**: Ellman's assay (free thiol pre/post heating); OPA (free amino groups); DSC (denaturation enthalpy); post-heating pH; LOX activity (pre-heating)
* **Safety check**: Acrylamide measurement (LC-MS/MS) in the 120°C condition given endogenous asparagine content

---

**Script to Reproduce Simulation**: the `run_pipeline.py` command that stood here drove the retired screening lane (deleted 2026-09-03). Reproduce the prediction with `python scripts/maillard.py predict <spec.yml>` using the precursor and process fields above; the kinetic core will refuse, with the reason, any target it cannot represent in a protein matrix.

**Benchmark ID (target)**: `soy_iso_ribose_cys_120C_pH5p8_meaty`

**Unlock Command**: `./scripts/docker_maillard.sh matrix-assertions`

## 10. Reference Curation Notes

* Verified: Mottram & Nobrega (2002) is correctly identified as *Formation of Sulfur Aroma Compounds in Reaction Mixtures Containing Cysteine and Three Different Forms of Ribose*, *Journal of Agricultural and Food Chemistry*, 50(14), 4080-4086, DOI 10.1021/jf0200826.
* Verified: the Hofmann & Schieberle 1998 MFT/FFT anchor corresponds to *Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol*, *Journal of Agricultural and Food Chemistry*, 46(1), 235-241, DOI 10.1021/jf9705983.
* Verified: Brands & van Boekel (2002) is *Kinetic Modeling of Reactions in Heated Monosaccharide-Casein Systems*, *Journal of Agricultural and Food Chemistry*, 50(23), 6725-6739, DOI 10.1021/jf011164h.
* Verified: Pratap-Singh et al. (2021) provides the pea/soy volatile baseline values used for the off-flavour anchors.
* Verified from full text: Nishimura & Abe (2024) uses a soy-protein workflow with a 75 mg/mL starting slurry, prepares MRPs from 62.5 mg/mL soy hydrolysate + 16.5 mM cysteine + 16.5 mM ribose, and heats at 95°C for 90 min before HS-SPME-GC/MS analysis.
* Verified limitation from full text: the volatile results are handled as relative peak areas with z-transformed clustering across n = 3 replicates; the paper does not provide absolute ppb concentrations for MFT/FFT or an internal-standard quantification route suitable for benchmark encoding.
* Pending curation before benchmark encoding: the exact Knol reference needed for any quantitative acrylamide benchmark anchor.
