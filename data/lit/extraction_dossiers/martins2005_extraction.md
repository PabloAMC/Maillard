# Martins & van Boekel 2005 — A kinetic model for the glucose/glycine Maillard reaction pathways

### Per-paper extraction, 2026-09-04, from the full text on disk (`data/articles/martins2005.pdf`; it had been there since 2026-08-28 under a file name the citation census did not match). **This paper is the measured backbone of the trunk lane (`src/kinetic_core/parameters.py`, `MARTINS_M4`, ten steps) and the source of the sulfur module's acid-yield analogy. The extraction changed METADATA only: the ten steps' 95 % HPDs were transcribed and two Ea intervals corrected; no rate constant or activation energy moved.**

**Provenance codes:** **[M]** measured and printed · **[F]** fitted by the authors · **[D]** derived by this extraction · **[NEG]** verified negative.

---

## §0. IDENTITY — CONFIRMED, DOI ON THE PAGE

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/martins2005.pdf` (13 pages); `martins2005b.pdf` is the companion paper *Kinetics of the glucose/glycine Maillard reaction pathways: influences of pH and reactant initial concentrations* (12 pages, not read here) | read in full |
| title | "A kinetic model for the glucose/glycine Maillard reaction pathways" | p. 257 |
| authors | Sara I.F.S. Martins and Martinus A.J.S. Van Boekel, Wageningen University | p. 257 |
| DOI | `10.1016/j.foodchem.2004.04.006` | p. 257 footer |
| journal | *Food Chemistry* 90 (2005) 257–269 | p. 257 header |
| dates | received 23 January 2004; revised 19 April 2004; accepted 22 April 2004 | p. 257 |
| system | D-glucose 0.2 mol/L + glycine 0.2 mol/L, 0.1 mol/L phosphate, pH 6.8 initial (uncontrolled during heating, drops to ≈5.5 in 4 h at 100 °C), 10 mL screw-capped Schott tubes, oil bath, 80/90/100/110/120 °C | p. 259, Fig. 1c |
| method | multiresponse kinetic modelling (Athena Visual Workbench), determinant criterion; model discrimination by posterior probability and AIC; reparameterised Arrhenius k = X·exp(−Y·Ea) with X the rate at T_av | pp. 259, 264–266 |

## §1. TABLE 2 (SCHEME 4, MODEL M4) AGAINST `MARTINS_M4`

Table 2 **[F]**, "Estimated parameters (X and Ea) with temperature dependence included as found by kinetic modelling for the proposed model in Scheme 4"; X = rate at T_av ≈ 100 °C, Ea in kJ/mol, ± = 95 % highest-posterior-density interval.

| step (Scheme 4) | code key | X printed | X in code | Ea printed | Ea in code (before → after) | verdict |
|---|---|---|---|---|---|---|
| 1 Glu + Gly → E1 (→ DFG) | k_schiff | 1.6e-5 ± 3.3e-7 | 1.6e-5 | 96.8 ± 2.8 | 96.8 ± 2.8 | exact |
| 2 Glu → Fru | k_glc_fru | 1.6e-3 ± 1.0e-4 | 1.64e-3 | 122.6 ± 5.2 | 123.0 ± 5.0 → ± 5.2 | rounded; HPD added |
| 3 Fru → Glu | k_fru_glc | 9.2e-3 ± 1.9e-3 | 9.15e-3 | 93.4 ± 1.9 | 93.0 ± 3.0 → ± 1.9 | Ea HPD corrected; k HPD added |
| 4 E1 → 3-DG | k_ama_tdg | 1.1e-2 ± 4.0e-4 | 1.11e-2 | 97.1 ± 1.7 | 97.0 ± 2.0 → ± 1.7 | rounded; HPD added |
| 5 3-DG → FA | k_tdg_fa | 3.5e-2 ± 6.4e-3 | 3.45e-2 | 29.6 ± 8.5 | 30.0 ± 9.0 → ± 8.5 | rounded; HPD added |
| 6 DFG → MG | k_ama_mgo | 7.1e-3 ± 4.6e-4 | 7.08e-3 | 124.5 ± 4.7 | 125.0 ± 5.0 → ± 4.7 | rounded; HPD added |
| 7 DFG → 1-DG | k_ama_odg | 1.6e-2 ± 6.8e-4 | 1.57e-2 | 107.3 ± 7.3 | 107.0 ± 7.0 → ± 7.3 | rounded; HPD added |
| 8 1-DG → AA | k_odg_aa | 1.4 ± 6.8e-2 | 1.45 | 75.7 ± 3.8 | 76.0 ± 4.0 → ± 3.8 | rounded; HPD added |
| 9 3-DG + Gly → Mel | k_tdg_mel | 8.1e-4 ± 1.7e-5 | 8.12e-4 (± 1.7e-5) | 95.2 ± 2.3 | 95.2 ± 2.3 | exact (browning readout is the hold-out) |
| 10 Glu → FA + AA | k_fru_acids | 4.4e-5 ± 3.6e-5 | 4.41e-5 | 236.7 ± 63.4 | 237.0 ± 36.0 → **± 63.4** | **Ea HPD was mis-transcribed** (the digits of the k HPD); the paper: "this step was not estimated very precisely" |

The code's three-significant-digit X values come from the thesis Table 5.2 (the same fit printed with one more digit); every one rounds to the printed value. Units: X in 1/min for the first-order steps and L/(mmol·min) for the bimolecular steps 1 and 9, matching the code.

## §2. WHAT THE SULFUR MODULE TAKES FROM THIS PAPER

| code site | claim | what the paper prints | verdict |
|---|---|---|---|
| `sulfur.py` (ACID pool), `species_sulfur.py` (ACID), `ph_state.py` §1 and §4 | "Martins 2005 measures that PART of the hexose deoxyosone flux terminates as formic acid (step 5) and acetic acid (step 8) while another part terminates as browning polymer (step 9), so the acid yield of the pentose analogue is a fraction between 0 and 1" | Scheme 4 **[F]**: 3-DG → FA (step 5) and 3-DG + Gly → Mel (step 9) compete for 3-DG; 1-DG → AA (step 8); Glu → FA + AA (step 10). Fig. 1c **[M]** at 100 °C, 4 h: acetic acid ≈ 17 mmol/L, formic ≈ 3.5 mmol/L; text p. 260: "25 % of the degraded D-glucose was transformed into acetic acid whereas only 5 % was formic acid"; pH fell from 6.8 to ≈ 5.5 (Fig. 1c). Mass balance (Fig. 5) closes at 96 %. | **confirmed**: the acid yield of the hexose sink is a measured fraction (≈ 30 % of degraded glucose as acid, the rest as fructose, DFG, MG and melanoidin), which is exactly why the pentose analogue's yield is a free fraction. |
| `parameters_furanic.py` (`_ODG_REFERENCE_K_100C = 1.45`, `_ODG_REFERENCE_EA = 76.0`, "corroborated by Knol 2010 (75 ± 10)") | 1-DG → acetic acid at 100 °C | step 8: X = 1.4 ± 0.068 /min, Ea 75.7 ± 3.8 | **confirmed** |
| `parameters.py` k_tdg_fa note ("Ea 30 ± 9 CONFLICTS with Knol 2010's 84 ± 14") | | step 5: 29.6 ± 8.5; text p. 267: "formation of 3-DG and 1-DG showed much lower temperature dependence in formic (29.6 kJ/mol) and acetic acid (75.7 kJ/mol) formation" | the conflict is real and now both sides are on disk |

## §3. WHAT THE PAPER DOES NOT GIVE

- **[NEG]** No pentose data of any kind; the sulfur module's use is an analogy, as it says.
- **[NEG]** pH was not controlled ("the effect of a non-constant pH, if any, is concealed in the rate constants", p. 267); the pH-dependence paper is `martins2005b.pdf`, not read here.
- **[NEG]** HMF was 0–20 µmol/L at pH 6.8 and is not in the scheme.
- Model discrimination (Table 1): hypothesis A (E1 reversible with DFG, formation of E1 irreversible) has PPB −19.05 vs −24.11 / −23.91; the code's `SCHIFF_AMADORI_SPLIT` note that "the source refuses the split" is consistent with the paper's finding that E1 is not rate-determining.

## §4. VERDICT

All ten trunk constants are Table 2's values; nine were already exact or correctly rounded and one Ea interval was mis-transcribed (step 10, ± 36 → ± 63.4). The ten k-side HPDs are now in the code for the first time. Nothing that feeds a prediction changed: the B1 fit report, not these intervals, is what the engine and the envelope read.
