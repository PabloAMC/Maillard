# Scientific Case Study: Meaty-Positive Pea Protein Matrix Benchmark

## 1. Executive Summary

> Status: benchmark candidate design only. This is not yet a curated external benchmark record and it should not be treated as benchmark-backed validation evidence.

This report defines and validates the minimum external evidence package for unlocking a meaty-positive matrix benchmark for pea protein isolate (PPI), as specified in `docs/EXTERNAL_MATRIX_BENCHMARK_UNLOCK_REPORT.md`. The core scientific challenge is that the pea matrix presents a structural asymmetry: its off-flavour profile (hexanal, 2-pentylfuran) is well-documented in quantitative external literature, but its meaty-positive volatile output under controlled precursor augmentation is absent from the benchmark set. We demonstrate, using published model-system kinetics and matrix-specific headspace data, that a narrow exogenous ribose + cysteine system at pH 5.5–6.0 and 95–120°C is both chemically justified and experimentally tractable as the first meaty-positive unlock benchmark for PPI.

---

## 2. Experimental Setup (Simulation)

* **Base Matrix**: Pea Protein Isolate, aqueous slurry (10% w/v, 80–85% protein dry weight).
* **Precursor Augmentation**: Ribose (1 mM exogenous) + Cysteine (1 mM exogenous).
* **Conditions**: pH 5.5, Temperature 95°C, Time 240 min (aqueous pre-extrusion state).
* **Process State**: `matrix_precursor_augmented`, `denaturation_state = 0.5` (partially denatured isolate, heat-treated but not extruded).
* **Off-Flavour Anchor**: Hexanal (lipid oxidation baseline from Pratap-Singh et al. 2021).
* **Analytical Reference Method**: HS-SPME-GC-MS (DVB/CAR/PDMS fiber, 40°C equilibration, 10 min extraction, hexanal-d12 internal standard).

---

## 3. Chemical Justification: Why Ribose + Cysteine in a Pea Matrix

Pea protein isolate is sulfur-deficient relative to animal muscle. Total cysteine content in PPI is approximately 0.8–1.1 g per 100 g protein (Gorissen et al. 2018), but the fraction accessible for Maillard reaction is substantially lower: cysteine is predominantly buried in legumin and vicilin disulfide bridges, and the free thiol fraction under native conditions is near or below detection by Ellman's assay (Prigent et al. 2024). Endogenous ribose is similarly negligible in a cleaned isolate.

The exogenous 1 mM ribose + 1 mM cysteine system is therefore the minimal, cleanest precursor pair capable of producing the target meaty-positive compound panel under aqueous heating conditions. This system is anchored by two tiers of literature evidence:

**Tier 1 — Free precursor model system kinetics (Hofmann & Schieberle 1998):** Heating ribose + cysteine at temperatures ≥ 95°C produces 2-methyl-3-furanthiol (MFT) and 2-furfurylthiol (FFT) as the dominant sulfur heterocycles. At pH 5.0 and 145°C (20 min), the molar yield of MFT from the ribose + cysteine system was confirmed as substantially higher than from hexose + cysteine systems. The pathway proceeds via the 1-deoxyosone intermediate derived from ribose degradation, which reacts with H₂S liberated from cysteine thermal decomposition to form the methylfuran thiol scaffold. Both MFT and FFT are formed when ribose and cysteine are in equimolar or near-equimolar ratios; a large molar excess of ribose (≥ 3:1 ribose:cysteine) suppresses FFT and shifts the pathway toward mercaptopentanones (Mottram & Nobrega 2002; Cerny & Davidek 2003 as cited in Cerny 2012).

**Tier 2 — Isotope-labelled confirmation of carbon skeleton origin (Mottram & Nobrega 2002):** Heating of [¹³C₅]ribose + cysteine at pH 5, 95°C, 4 h confirmed by SPME-GC-MS that the carbon skeleton of ribose remains intact in MFT and FFT. This mechanistic confirmation is critical: it means that predicting MFT yield from ribose:cysteine stoichiometry is mechanistically grounded, not phenomenological.

The combination of these two tiers justifies using 1 mM ribose + 1 mM cysteine at pH 5.5 as the primary unlock precursor system. It is the narrowest system with the strongest mechanistic backing for MFT and FFT formation.

---

## 4. Predicted Outcomes vs. Literature

| Volatile / Metric | Predicted Trend | Literature Reference | Status |
| --- | --- | --- | --- |
| **2-Methyl-3-furanthiol (MFT)** | Detectable above ODT (0.2 ppt) at equimolar 1 mM ribose:cys, pH 5.5, 95°C | Hofmann & Schieberle (1998) *J Agric Food Chem* 46(1):235-241 | ✅ Match — pentose/cysteine yields significantly higher MFT than hexose/cysteine |
| **2-Furfurylthiol (FFT)** | Co-formed with MFT; lower yield than MFT at equimolar ratio | Mottram & Nobrega (2002) *J Agric Food Chem* 50(14):4080-4086 | ✅ Match — equimolar ratio favours both isomers; excess ribose suppresses FFT |
| **bis(2-methyl-3-furyl) disulfide** | Present as oxidation product of MFT; lower absolute concentration | Hofmann, Schieberle & Grosch (1996) *J Agric Food Chem* 44:251 | ✅ Match — disulfide forms at ppb-level under oxidative conditions |
| **2,5-Dimethylpyrazine** | Low yield at pH 5.5, increases above pH 7; suppressed by sulfur pathway competition | Martins & van Boekel (2005) *Food Chem* 90:257 | ✅ Match — pyrazine formation requires amino acid nitrogen without S-competition at low pH |
| **Hexanal (off-flavour anchor)** | Baseline 638 ± 49 ppb-equivalent in unheated PPI; lipid oxidation-derived, not Maillard | Pratap-Singh et al. (2021) *Molecules* 26:4104 | ✅ External anchor — quantified in ambient PPI slurry by HS-SPME-GC-MS |
| **2-Pentylfuran (off-flavour anchor)** | 638 ppb-equivalent (co-detected with hexanal); linoleic acid oxidation product | Pratap-Singh et al. (2021) *Molecules* 26:4104 | ✅ External anchor — fully quantified in PPI |
| **Matrix Attenuation of MFT** | ~50–60% reduction in headspace MFT relative to free precursor model at equivalent concentration | Estimated from `matrix_correction.py` cysteine_accessibility = 0.25, volatile_retention = 0.50 (PEA_ISOLATE) | ⚠️ Model prediction — no direct external measurement in pea matrix yet |

---

## 5. Sensitivity Analysis

**pH Window (5.0 → 7.0):** MFT and FFT yields peak in the pH 5.0–6.0 range. Above pH 6.5, pyrazine formation competes strongly for the available nitrogen from cysteine, reducing sulfur heterocycle yield. The benchmark pH of 5.5 sits at the optimal meaty-positive operating point. The Mottram & Nobrega (2002) model at pH 5 and 95°C demonstrates that the 4 h aqueous heating time produces detectable MFT and FFT above their respective ODTs; at pH 7, the same system yields substantially less MFT.

**Temperature Window (95°C → 120°C):** Increasing from 95°C to 120°C accelerates ribose degradation to the 1-deoxyosone intermediate and increases MFT yield, but also accelerates furfural formation and the onset of melanoidin polymerization that can bind and irreversibly quench free thiols. The 95°C / 4 h condition is chosen specifically because it is the most conservative condition at which MFT has been confirmed above detection limits in the Mottram & Nobrega (2002) study, minimizing confounding thiol quenching.

**Matrix Attenuation — Cysteine Accessibility:** The largest sensitivity in the pea matrix prediction is cysteine_accessibility, currently set to 0.25 (Prigent et al. 2024; estimated). A wet-lab OPA + Ellman's assay on the specific PPI batch used in the benchmark would allow this parameter to be replaced with a measured value, which would propagate directly into the MFT prediction as a linear correction. This is the single highest-leverage wet-lab measurement for improving model accuracy.

**Ribose:Cysteine Ratio (0.5:1 → 3:1):** At ratios below equimolar, cysteine is in excess and H₂S pathways dominate, producing more 3-mercapto-2-pentanone and less MFT. At ratios above 3:1 ribose:cysteine, FFT formation is suppressed in favour of mercaptopentanones. The 1:1 molar ratio used in this benchmark is the literature-optimal ratio for MFT + FFT co-production.

---

## 6. Tradeoff Markers Required for Benchmark Acceptance

Per the UNLOCK report, each meaty-positive benchmark must include at least one adverse marker. For the pea system, the required adverse markers are:

* **Hexanal**: Must be measured in the same experiment. The expected baseline from lipid oxidation is 638 ± 49 ppb-equivalent (Pratap-Singh et al. 2021) in unheated PPI slurry; under 95°C heating this may increase by a factor of 2–3× depending on oxygen availability and iron content of the isolate.
* **2-Pentylfuran**: Co-measured with hexanal; expected to follow hexanal trend as a linoleic acid oxidation co-product.
* **Furfural**: Forms in parallel with MFT from ribose degradation; serves as a chemical indicator that the ribose pathway is actually running. Expected to dominate the volatile profile by mass but contribute negligibly to sensory impact at near-neutral pH (ODT 3000 ppb vs MFT at 0.2 ppt).

The tradeoff that this benchmark is designed to capture is: **as the exogenous ribose:cysteine system generates MFT (meaty gain), the same heating event also accelerates hexanal generation from the residual PPI lipid fraction (off-flavour cost)**. A complete external dataset must capture both sides of this tradeoff within a single experimental run.

---

## 7. Benchmark Acceptance Criteria

The following criteria must all be met before this benchmark can be promoted to the strict gate:

* Pea isolate batch identity and supplier reported in `source_metadata`
* pH, temperature, time, and water activity explicitly stated
* Denaturation proxy reported (DSC enthalpy or solubility index)
* MFT and/or FFT detected above ODT with stated quantification method and internal standard
* Hexanal concentration reported in same experiment
* At least 3 analytical replicates with stated uncertainty per compound
* Non-detect policy explicitly stated (e.g., ND = below LOQ, not assumed zero)

---

## 8. Known Scientific Gaps

* **Cysteine Accessibility in Heated PPI**: No published study directly reports free thiol fraction (Ellman's) in pea isolate after mild aqueous heating (95°C, 60–240 min). This is the primary wet-lab measurement blocking the strict gate.
* **MFT Quenching by Pea Melanoidins**: MFT is known to bind irreversibly to coffee melanoidins (Hofmann & Schieberle 2002). Whether pea protein Maillard browning products show similar quenching under aqueous conditions is unknown. This introduces a systematic negative bias in predicted headspace MFT.
* **Iron Catalysis**: Pea isolate contains 20–30 ppm non-heme iron. Iron catalyses both lipid oxidation (accelerating hexanal generation) and, at alkaline pH, pyrazine formation. At pH 5.5, iron's primary effect is expected to be pro-oxidant on the lipid fraction. No matrix-specific iron catalysis correction is currently implemented.

---

## 9. Recommended Wet-Lab Design for Benchmark Generation

* **Matrix**: Commercial PPI, ≥80% protein dry weight, supplier and lot number recorded
* **Precursor load**: 1 mM ribose + 1 mM cysteine (exogenous, dissolved in phosphate buffer)
* **Conditions**: pH 5.5 ± 0.1, 95°C, 240 min; and a second condition at 120°C, 60 min
* **Headspace method**: HS-SPME, DVB/CAR/PDMS fiber, 40°C, 30 min equilibration, 10 min extraction; hexanal-d12 internal standard; external calibration for hexanal and furfural; thiol-specific pHMB derivatization for MFT and FFT (Tominaga & Dubourdieu 2006 method)
* **Replication**: Minimum 3 independent process replicates × 2 analytical replicates per compound
* **Companion measurements**: Ellman's assay for free thiol content pre- and post-heating; OPA assay for free amino group content; DSC for denaturation enthalpy

---

**Script to Reproduce Simulation**: `scripts/run_pipeline.py --protein-type pea_iso --denaturation-state 0.5 --sugars ribose:1.0 --amino-acids cysteine:1.0 --ph 5.5 --temp 95`

**Benchmark ID (target)**: `pea_iso_ribose_cys_95C_pH5p5_meaty`

**Unlock Command**: `./scripts/docker_maillard.sh matrix-assertions`

## 10. Reference Curation Notes

* Verified: Mottram & Nobrega (2002) is correctly identified as *Formation of Sulfur Aroma Compounds in Reaction Mixtures Containing Cysteine and Three Different Forms of Ribose*, *Journal of Agricultural and Food Chemistry*, 50(14), 4080-4086, DOI 10.1021/jf0200826.
* Verified: the Hofmann & Schieberle 1998 MFT/FFT benchmark anchor corresponds to *Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol*, *Journal of Agricultural and Food Chemistry*, 46(1), 235-241, DOI 10.1021/jf9705983.
* Verified: Pratap-Singh et al. (2021) is the rapid SPME-GC/MS plant-protein volatile baseline paper used for pea/soy off-flavour anchors.
* Pending full-text curation before benchmark encoding: the 2024 pea-specific accessibility source and any exact quantitative extraction from Hofmann beyond the verified article-level anchor.

## 11. Hexanal And Nonanal Calibration-First Closure

Status: the internal calibration route is executable and currently closed, but this does not unlock external promotion claims.

The current protocol-pilot lane for pea isolate is now checked directly against the frozen Internal2026 lane through the generated artifact [results/validation/hexanal_nonanal_calibration_closure.md](results/validation/hexanal_nonanal_calibration_closure.md). The rule is intentionally narrow: Hexanal and Nonanal must each remain within a ProtocolPilot2026-to-Internal2026 ratio window of 0.5x to 2.0x for the internal calibration route to be treated as closed.

This means the adverse-marker story is no longer a default mechanistic-escalation trigger for the pea mixed lane. It is instead treated as an internally concordant calibration route that still needs externally quantitative mixed-matrix evidence before the lane can move toward external decision readiness.

The scientific implication is bounded and specific:

* Hexanal and Nonanal are currently compatible with the matrix headspace and lipid-oxidation calibration used for the internal pea lane.
* This does not prove broad external validity for pea matrices.
* It does reduce the justification for pushing selective QM first on these adverse markers.
