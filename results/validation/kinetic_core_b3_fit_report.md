# Kinetic core, Build Wave B3 -- the acrylamide / safety fit

Module 3 of `docs/reference/FIT_HOLDOUT_DECLARATION.md`.

- network: **51 species, 31 reactions** (15 inherited from B1's trunk, 16 new), carbon, nitrogen and sulfur balance enforced at import
- objective: **30 declared FIT rows**, **11 free parameters**, final cost **52.651**, reduced chi-square **5.54**, dof **19**
- acrylamide ELIMINATION channels: **5** (`a_acr_dp`, `a_acr_cys`, `a_acr_gln`, `a_acr_lys`, `a_acr_ala`). The old FAST lane had none, which is why it was ~40x under-responsive to time.
- competition multipliers: **0**. Competition is two named mass-action channels per amino acid, one of them measured.
- pH terms: **0**. a_w terms: **0**. Occurrences of the unsourceable 129 kJ/mol acrylamide barrier: **0**.

## Read this before reading anything else

**2 of 11 free parameters are individually identified.** Seven of the sixteen steps in this module carry a literature rate constant and are NOT fitted, which is why the panel can afford eleven free parameters at all -- but the row-by-row agreement below is still not evidence that the model is right. The rows worth reading as evidence are the ones with nothing fitted to them:

- `claeys_cys_kE_160C` and `dv2_cysteine_reduction`: **both of cysteine's channels are measured**, so these are parameter-free predictions of a competitor's effect.
- `claeys_*_EaF`: the formation barrier is a composite of two measured barriers and one fitted one, and the five systems share it.
- `knol2005_k6_*`: five temperatures from a third lab against a constant fitted to two others.

## Row-by-row

| row | kind | observed | predicted | fold / delta | w. residual | source |
|---|---|---:|---:|---:|---:|---|
| `claeys_control_kF_160C` | rate | 0.000451 | 0.0002685 | 0.60x | -0.75 | Claeys 2005 Biotechnol Prog 21:1525-1530 Table 2, control sy |
| `claeys_control_kE_160C` | rate | 0.1111 | 0.1265 | 1.14x | +0.19 | Claeys 2005 Biotechnol Prog 21:1525-1530 Table 2, control sy |
| `claeys_control_EaF` | ea | 168.2 | 164.3 | -4.0 kJ/mol | -0.16 | Claeys 2005 Biotechnol Prog 21:1525-1530 Table 2, control sy |
| `claeys_control_EaE` | ea | 167.2 | 136.1 | -31.1 kJ/mol | -1.24 | Claeys 2005 Biotechnol Prog 21:1525-1530 Table 2, control sy |
| `claeys_gln_kF_160C` | rate | 0.00164 | 0.0002685 | 0.16x | -2.62 | Claeys 2005 Table 2, glutamine competitor system, k_Fref |
| `claeys_gln_kE_160C` | rate | 0.2741 | 0.2538 | 0.93x | -0.11 | Claeys 2005 Table 2, glutamine competitor system, k_Eref |
| `claeys_gln_EaF` | ea | 166.8 | 164.3 | -2.5 kJ/mol | -0.10 | Claeys 2005 Table 2, glutamine competitor system, Ea_F |
| `claeys_gln_EaE` | ea | 103.9 | 101 | -2.9 kJ/mol | -0.11 | Claeys 2005 Table 2, glutamine competitor system, Ea_E |
| `claeys_cys_kF_160C` | rate | 0.000501 | 0.0002676 | 0.53x | -0.91 | Claeys 2005 Table 2, cysteine competitor system, k_Fref |
| `claeys_cys_kE_160C` | rate | 0.2687 | 0.1313 | 0.49x | -1.04 | Claeys 2005 Table 2, cysteine competitor system, k_Eref |
| `claeys_cys_EaF` | ea | 206.3 | 164.5 | -41.8 kJ/mol | -1.67 | Claeys 2005 Table 2, cysteine competitor system, Ea_F |
| `claeys_cys_EaE` | ea | 180 | 114.9 | -65.1 kJ/mol | -2.61 | Claeys 2005 Table 2, cysteine competitor system, Ea_E |
| `claeys_lys_kF_160C` | rate | 0.000587 | 0.0002513 | 0.43x | -1.23 | Claeys 2005 Table 2, lysine competitor system, k_Fref |
| `claeys_lys_kE_160C` | rate | 0.2802 | 0.1565 | 0.56x | -0.84 | Claeys 2005 Table 2, lysine competitor system, k_Eref |
| `claeys_lys_EaF` | ea | 179.3 | 168.8 | -10.5 kJ/mol | -0.42 | Claeys 2005 Table 2, lysine competitor system, Ea_F |
| `claeys_lys_EaE` | ea | 140 | 124.4 | -15.6 kJ/mol | -0.62 | Claeys 2005 Table 2, lysine competitor system, Ea_E |
| `claeys_ala_kF_160C` | rate | 0.000465 | 0.0002685 | 0.58x | -0.79 | Claeys 2005 Table 2, alanine competitor system, k_Fref |
| `claeys_ala_kE_160C` | rate | 0.1031 | 0.1265 | 1.23x | +0.30 | Claeys 2005 Table 2, alanine competitor system, k_Eref |
| `claeys_ala_EaF` | ea | 173.3 | 164.3 | -9.0 kJ/mol | -0.36 | Claeys 2005 Table 2, alanine competitor system, Ea_F |
| `claeys_ala_EaE` | ea | 169.7 | 136.1 | -33.6 kJ/mol | -1.34 | Claeys 2005 Table 2, alanine competitor system, Ea_E |
| `dv1_glucose_kE_160C` | rate | 0.1 | 0.1265 | 1.27x | +0.34 | De Vleeschouwer 2009 I Table 3, glucose k_Eref = 0.10 +/- 0. |
| `dv1_glucose_EaE` | ea | 113.2 | 136.1 | +22.9 kJ/mol | +0.92 | De Vleeschouwer 2009 I Table 3, glucose Ea_E = 113.2 +/- 32. |
| `dv2_cysteine_reduction` | ceiling | 0.01 | 0.8178 | 81.78x | +6.38 | De Vleeschouwer 2009 II abstract: '>99 % reduction' of acryl |
| `knol2005_k6_120C` | rate | 0.00796 | 0.002705 | 0.34x | -1.56 | Knol 2005 Table 1, k6 at 120 C |
| `knol2005_k6_140C` | rate | 0.0281 | 0.02031 | 0.72x | -0.47 | Knol 2005 Table 1, k6 at 140 C |
| `knol2005_k6_160C` | rate | 0.0881 | 0.1265 | 1.44x | +0.52 | Knol 2005 Table 1, k6 at 160 C |
| `knol2005_k6_180C` | rate | 0.25 | 0.6708 | 2.68x | +1.43 | Knol 2005 Table 1, k6 at 180 C |
| `knol2005_k6_200C` | rate | 0.65 | 3.089 | 4.75x | +2.26 | Knol 2005 Table 1, k6 at 200 C |
| `knol2005_EaF` | ea | 94.4 | 189.5 | +95.1 kJ/mol | +3.81 | Knol 2005 Table 1, Ea(k4) = 94.4 +/- 11 kJ/mol (formation) |
| `spi_extrusion_130C_acrylamide_ppb` | ppb | 150 | 0.03532 | 0.00x | -3.63 | data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json |

## Parameters

| parameter | value | unit | 95% half-width | identified? |
|---|---:|---|---:|---|
| `k_int1_mel` | -4.496 | log10 k(160 C) | 1316.322 | **UNIDENTIFIED** |
| `k_acr_dp` | -0.898 | log10 k(160 C) | 0.445 | yes |
| `k_gln_glc` | -7.836 | log10 k(160 C) | 2157390.472 | **UNIDENTIFIED** |
| `k_lys_glc` | -2.596 | log10 k(160 C) | 27.164 | **UNIDENTIFIED** |
| `k_ala_glc` | -6.872 | log10 k(160 C) | 446550.909 | **UNIDENTIFIED** |
| `k_acr_gln` | -1.895 | log10 k(160 C) | 3.126 | **UNIDENTIFIED** |
| `k_acr_lys` | -2.492 | log10 k(160 C) | 5.460 | **UNIDENTIFIED** |
| `k_acr_ala` | -8.825 | log10 k(160 C) | 2305878.147 | **UNIDENTIFIED** |
| `Ea_int1_mel` | 260.000 | kJ/mol | 122057.374 | **UNIDENTIFIED** |
| `Ea_acr_dp` | 136.109 | kJ/mol | 48.105 | yes |
| `Ea_competitor_sugar` | 20.000 | kJ/mol | 4383.492 | **UNIDENTIFIED** |

'Identified' means the 95 % HALF-WIDTH is inside 1 dex for a rate constant and 60 kJ/mol for a barrier -- the test is on the interval, not on the raw standard error, because a constant whose 95 % interval spans three orders of magnitude has not been determined by anything. `flat` means the direction is numerically flat and no interval exists at all.

The three acrylamide-scavenging channels' activation energy is NOT in that table because it is not fitted: it is held at the measured Ea_E2 = 51.3 kJ/mol.

## Deliberate under-fit, stated before the fit was run

- **Claeys 2005 T2, the GLUTAMINE competitor system's formation constant** -- UNDER-PREDICT it by roughly the measured promotion factor. Every competition channel in this module can only REMOVE precursor or REMOVE acrylamide, so the model's glutamine system can never form MORE acrylamide than its control.
  Why no term was added: Because the shape belongs to a HOLD-OUT. Inventory sec. B5.5: glutamine's acrylamide promotion GROWS with temperature in the liquid system and SHRINKS with temperature at a_w 0.92 -- same lab, same amino acid, and neither paper remarks on it. The liquid half is Claeys (FIT); the a_w 0.92 half is De Vleeschouwer II's glutamine system (HOLD-OUT). A promotion term fitted to the FIT half would be a term built toward a hold-out's shape, and the build brief for this wave forbids exactly that.

## Cross-lab conflicts carried, not averaged

- **acrylamide ELIMINATION rate constant at 160 C** (1/min): Claeys 2005 T2 control, aqueous pH 6 = 0.1111, De Vleeschouwer 2009 I glucose, a_w 0.92 = 0.1, Knol 2005 T1 k6, aqueous = 0.0881. Spread: 1.26x across three labs, two water activities and two decades. THE STRONGEST CROSS-LAB AGREEMENT IN THE MODULE, and it is on the very constant the inventory calls unidentifiable (sec. C.6). Both statements are true: the RATE is reproducible, the DEGRADATION MECHANISM behind it is not constrained by any of the three.
- **acrylamide ELIMINATION activation energy** (kJ/mol): Claeys 2005 T2 control = 167.21, De Vleeschouwer 2009 I glucose = 113.2, Knol 2005 T1 k6 = 85.1. Spread: 2.0x, intervals do not all overlap (167.21 +/- 4.30 vs 85.1 +/- 14). A STANDING OWNER QUESTION. The same three labs agree to 1.26x on the RATE at 160 C while spanning 2x on the BARRIER, which is the signature of a reference temperature sitting in the middle of every data set: k(T_ref) is well determined and the slope is not.
- **acrylamide FORMATION activation energy** (kJ/mol): Claeys 2005 T2 control, aqueous = 168.25, De Vleeschouwer 2009 I glucose, a_w 0.92 = 159.2, Knol 2005 T1 k4, aqueous = 94.4. Spread: 74 kJ/mol between Claeys and Knol; intervals (168.25 +/- 3.80 and 94.4 +/- 11) are nowhere near overlapping. Claeys and De Vleeschouwer (same lab, different a_w) agree to 9 kJ/mol; Knol is the outlier by 65-74. This is the largest unresolved numerical conflict in Module 3.
- **condensation barrier, Glc + Asn -> Schiff base** (kJ/mol): Knol 2005 T1 k1 = 57.6, Martins 2005 step 1 (glycine, B1's operative value) = 96.8, De Vleeschouwer 2009 I k_INTg (this module's operative value) = 117.5. Spread: 60 kJ/mol; no two intervals overlap. Owner question, already flagged by the inventory as Z2 #17.

## Diagnostics (computed, not scored)

- cysteine crossover k_E/k_E2: model **2.56 mmol/L** against the literature's 2.0 mmol/L.
- Claeys' own internal validation at 160 C: plateau 2886 ppb and t_max 49.8 min from the published constants; this model gives **875 ppb** at **32.6 min**.

### The binding constraint is B1's trunk, not a Module 3 parameter

Glucose remaining at the end of the 45 min window, Claeys control system, as a fraction of the charge:

| T | glucose left | fructose left | asparagine left |
|---|---:|---:|---:|
| 140C | 0.31 | 0.108 | 0.69 |
| 160C | 1.76e-06 | 3.46e-07 | 0.288 |
| 180C | 0 | 0 | 0.0127 |
| 200C | 0 | 0 | 7.11e-07 |

The asparagine initiation is SECOND ORDER, so it stops when the sugar does. Martins' sugar constants were measured over 80-120 C and this panel runs at 140-200 C; evaluated there, the trunk's own fructose fragmentation lane empties the sugar within minutes. **That, and not the Int1 partition, is what caps every Claeys formation row** -- and it is why the fitted partition sits at its lower bound, buying back acrylamide the sugar lane has already spent. It is a B1 trunk extrapolation, flagged on every run by the module's own warnings, and it is the first thing a later wave should look at.

### The cysteine reduction is real but TRANSIENT

| window | ratio to control | reduction |
|---:|---:|---:|
| 2 min | 0.0099 | 99.0 % |
| 5 min | 0.0121 | 98.8 % |
| 10 min | 0.0484 | 95.2 % |
| 20 min | 0.5504 | 45.0 % |
| 45 min | 0.8179 | 18.2 % |
| 120 min | 0.8179 | 18.2 % |

The mechanism: the same paper's k_Y = 0.35 min^-1 empties the cysteine pool with a 2 min half-life while its k_F = 3.57e-3 min^-1 keeps forming acrylamide for hours, so the scavenger is gone long before the acrylamide peak. The `dv2_cysteine_reduction` row is scored at the panel's declared 45 min window and FAILS there, while the same model reproduces the published >99 % at two minutes. This is an internal tension in De Vleeschouwer II's own parameter set -- a scavenger with a two-minute half-life cannot durably suppress a product that forms over hours -- exposed by writing the chemistry as mass action rather than as two lumped constants. It is reported, not resolved, and the row was NOT re-scored: the window was declared for every row before the fit; picking the window that passes would be fitting the objective to the answer.

## Repo defect fixed en route

data/lit/safety_reference_payloads.json entries[27] (knol_2005_acrylamide_kinetics) stated the Knol 2005 acrylamide barriers as 52.1 (formation) and 72.9 (degradation) kJ/mol. The true pair is 94.4 +/- 11 and 85.1 +/- 14 (inventory sec. A.2 lines 202-203, errata #3).

CORRECTED in place with a dated provenance note. No test asserted the old values; the schema is unchanged.

## Hold-out exposure disclosure

The build brief directed this wave to read k3_final_parameter_inventory.md sec. A.2 and k1_kinetic_parameters.md sec. 2, and BOTH print the declared HOLD-OUT columns beside the FIT ones: De Vleeschouwer 2009 Part I's FRUCTOSE and SUCROSE columns (sec. 2b prints all three sugars in one table, and the inventory prints the sucrose hydrolysis and isomerisation constants in its own sec. A.2), and De Vleeschouwer 2009 Part II's GLUTAMINE column (sec. 2c prints Gln and Cys side by side). Inventory sec. B5.5 additionally prints the glutamine promotion percentages that are the shape of the Part II hold-out.

Nothing held out entered a parameter, a bound, an initialisation or a fit row. tests/unit/test_kinetic_core_b3.py enforces it mechanically with a literal-grep firewall over the executable code of every runtime file and the fit script, in the pattern Wave B2 established.

## Out of scope for this wave

- **sucrose (and any non-reducing sugar)** -- strands: De Vleeschouwer 2009 I's sucrose system, whose hydrolysis and isomerisation constants are measured. It is a declared HOLD-OUT, so the module may not carry it -- but that means the hold-out is reported UNSCOREABLE rather than predicted. Reason: Adding a sucrose species would require its hydrolysis rate, and the only measurement of that rate sits in the hold-out column. Inventing one to be able to score the hold-out would defeat the hold-out.
- **the composed sulfur + acrylamide network** -- strands: the interaction between B2's four thiol sinks and B3's cysteine channels; any prediction of MFT/FFT and acrylamide in the same pot; the cysteine-rich extrusion systems where both matter. Reason: Two independently calibrated cysteine sink sets would consume the same pool. Resolving which sink applies at which temperature and a_w is a decision about the corpus, not a coding step.
- **glutamine promotion of acrylamide** -- strands: Claeys' glutamine competitor row, which the model under-predicts by roughly its measured promotion factor. Reason: The promotion's temperature SHAPE is a declared hold-out (inventory sec. B5.5). A term fitted to the FIT half would be a term built toward the hold-out.
- **asparagine's OTHER Maillard products** -- strands: everything between Int1 and the melanoidin pool: Strecker aldehydes, pyrazines, the 3-aminopropionamide route. Reason: The Int1 partition is measured only as an aggregate (and its literature constant is one the authors mark 'no physical meaning'), so resolving the non-acrylamide branch would be fitting structure to one number.
