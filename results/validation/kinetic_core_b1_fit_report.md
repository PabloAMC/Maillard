# Kinetic core B1 -- fit report

Generated 2026-08-28 on `audit-remediation` @ `05a5fae8` (dirty).

## What was fitted

- Fit corpus: `data/lit/timeseries/martins2005_glucose_glycine_80_100_120C_pH68.yml (nine responses)`
- Responses in the objective (9): glucose, glycine, fructose, DFG, 3-DG, 1-DG, methylglyoxal, formic_acid, acetic_acid
- **Never read**: melanoidins -- the Module 4 hold-out (docs/reference/FIT_HOLDOUT_DECLARATION.md D.6, Module 4)

- excluded: the melanoidin/browning response of the same file (Module 4 HOLD-OUT)
- excluded: epsilon = 0.64 L/(mmol*cm) (Martins) and 282 L/(mol*cm) (Knol) -- HOLD-OUT, absent from src/kinetic_core/ and from this script
- excluded: data/benchmarks/** (never read)
- excluded: data/lit/timeseries/brands_sugar_casein_120C_pH68.yml (directional C/N diagnostic only)

The ten Martins trunk constants are **fixed at their measured values** and are not fitted. Only four consumption steps, for which the corpus has no rate constant at all, are estimated.

## Variant A -- measured melanoidin sink

- observations 354, free parameters 8, dof 346
- **objective value (0.5 * SSR of weighted log residuals): 1989.9384**
- reduced chi-square 11.503; RMS weighted residual 3.353
- multistart: 1 of 12 random starts reached this optimum

| parameter | k(100 C) | 95% CI | Ea kJ/mol | 95% CI | at bound? | verdict (k) | verdict (Ea) |
|---|---:|---|---:|---|---|---|---|
| `k_glc_frag` | 1e-08 /min | unbounded | 180.7 | unbounded | **YES -- REJECTED BY THE DATA** | UNIDENTIFIED (flat direction; the data do not bound this parameter) | UNIDENTIFIED (flat direction) |
| `k_mgo_mel` | 0.02273 /min | 0.0136 - 0.038 | 20.0 | -42 - 82 | no | constrained (95% CI inside a factor of 2) | UNCONSTRAINED SIGN (95% CI reaches -42 kJ/mol; an Ea <= 0 is unphysical, so the data do not determine this step's temperature dependence) |
| `k_fa_frag` | 3.465e-08 /min | unbounded | 20.5 | unbounded | no | UNIDENTIFIED (flat direction; the data do not bound this parameter) | UNIDENTIFIED (flat direction) |
| `k_aa_frag` | 0.01181 /min | 0.00816 - 0.0171 | 20.0 | -21 - 61 | no | constrained (95% CI inside a factor of 2) | UNCONSTRAINED SIGN (95% CI reaches -21 kJ/mol; an Ea <= 0 is unphysical, so the data do not determine this step's temperature dependence) |

- **channels the data REQUIRE** (bounded away from zero): `k_mgo_mel`, `k_aa_frag`
- **channels the data REJECT** (pinned to the lower search bound, i.e. the fit prefers no such channel at all): `k_glc_frag`
- **channels the data cannot resolve** (flat direction; the estimate printed above is not information): `k_fa_frag`

### Per-response agreement (in the objective)

| species | T (C) | n | median fold error | max fold error | signed log10 bias |
|---|---:|---:|---:|---:|---:|
| 1-DG | 80 | 11 | 1.09x | 1.38x | -0.016 |
| 1-DG | 100 | 14 | 1.15x | 1.34x | +0.036 |
| 1-DG | 120 | 9 | 1.24x | 1.35x | +0.071 |
| 3-DG | 80 | 16 | 1.19x | 1.53x | +0.077 |
| 3-DG | 100 | 12 | 1.10x | 1.18x | +0.040 |
| 3-DG | 120 | 8 | 1.19x | 1.38x | +0.043 |
| DFG | 80 | 16 | 1.18x | 1.86x | -0.072 |
| DFG | 100 | 15 | 1.03x | 1.25x | +0.007 |
| DFG | 120 | 13 | 1.13x | 1.57x | -0.027 |
| acetic_acid | 80 | 12 | 1.09x | 1.31x | +0.039 |
| acetic_acid | 100 | 17 | 1.06x | 1.25x | -0.025 |
| acetic_acid | 120 | 13 | 1.11x | 1.20x | +0.015 |
| formic_acid | 80 | 8 | 1.12x | 1.27x | -0.018 |
| formic_acid | 100 | 12 | 1.15x | 1.70x | -0.059 |
| formic_acid | 120 | 11 | 2.07x | 2.49x | -0.316 |
| fructose | 80 | 13 | 1.47x | 1.91x | +0.168 |
| fructose | 100 | 16 | 1.45x | 1.69x | +0.161 |
| fructose | 120 | 13 | 1.40x | 2.60x | +0.146 |
| glucose | 80 | 17 | 1.02x | 1.04x | -0.007 |
| glucose | 100 | 15 | 1.02x | 1.08x | -0.008 |
| glucose | 120 | 12 | 1.07x | 1.15x | -0.028 |
| glycine | 80 | 19 | 1.02x | 1.05x | -0.011 |
| glycine | 100 | 13 | 1.01x | 1.01x | -0.002 |
| glycine | 120 | 12 | 1.01x | 1.04x | -0.004 |
| methylglyoxal | 80 | 10 | 1.26x | 1.87x | +0.055 |
| methylglyoxal | 100 | 17 | 1.16x | 2.25x | -0.041 |
| methylglyoxal | 120 | 10 | 1.16x | 1.38x | +0.060 |

## Variant B -- reactant-side melanoidin sink (the pre-registered predictor)

- observations 354, free parameters 9, dof 345
- **objective value (0.5 * SSR of weighted log residuals): 1984.2586**
- reduced chi-square 11.503; RMS weighted residual 3.348
- multistart: 2 of 12 random starts reached this optimum

| parameter | k(100 C) | 95% CI | Ea kJ/mol | 95% CI | at bound? | verdict (k) | verdict (Ea) |
|---|---:|---|---:|---|---|---|---|
| `k_glc_frag` | 1.074e-08 /min | unbounded | 20.6 | unbounded | no | UNIDENTIFIED (flat direction; the data do not bound this parameter) | UNIDENTIFIED (flat direction) |
| `k_mgo_mel` | 0.02276 /min | 0.0136 - 0.0381 | 20.5 | unbounded | no | constrained (95% CI inside a factor of 2) | UNIDENTIFIED (flat direction) |
| `k_fa_frag` | 1e-08 /min | unbounded | 236.9 | unbounded | **YES -- REJECTED BY THE DATA** | UNIDENTIFIED (flat direction; the data do not bound this parameter) | UNIDENTIFIED (flat direction) |
| `k_aa_frag` | 0.01178 /min | 0.00815 - 0.017 | 26.3 | unbounded | no | constrained (95% CI inside a factor of 2) | UNIDENTIFIED (flat direction) |

- **channels the data REQUIRE** (bounded away from zero): `k_mgo_mel`, `k_aa_frag`
- **channels the data REJECT** (pinned to the lower search bound, i.e. the fit prefers no such channel at all): `k_fa_frag`
- **channels the data cannot resolve** (flat direction; the estimate printed above is not information): `k_glc_frag`

**Melanoidin sink, estimated from the reactant side only:** k(100 C) = 0.0007528 L/(mmol*min) (Martins measured 0.000812; ratio fitted/measured = 0.927). well constrained (95% CI inside +/-26%).

### Per-response agreement (in the objective)

| species | T (C) | n | median fold error | max fold error | signed log10 bias |
|---|---:|---:|---:|---:|---:|
| 1-DG | 80 | 11 | 1.09x | 1.38x | -0.016 |
| 1-DG | 100 | 14 | 1.15x | 1.34x | +0.036 |
| 1-DG | 120 | 9 | 1.24x | 1.35x | +0.071 |
| 3-DG | 80 | 16 | 1.23x | 1.57x | +0.090 |
| 3-DG | 100 | 12 | 1.15x | 1.23x | +0.057 |
| 3-DG | 120 | 8 | 1.19x | 1.44x | +0.071 |
| DFG | 80 | 16 | 1.18x | 1.86x | -0.072 |
| DFG | 100 | 15 | 1.03x | 1.25x | +0.007 |
| DFG | 120 | 13 | 1.13x | 1.57x | -0.027 |
| acetic_acid | 80 | 12 | 1.11x | 1.32x | +0.046 |
| acetic_acid | 100 | 17 | 1.06x | 1.25x | -0.024 |
| acetic_acid | 120 | 13 | 1.10x | 1.20x | +0.012 |
| formic_acid | 80 | 8 | 1.10x | 1.25x | -0.016 |
| formic_acid | 100 | 12 | 1.13x | 1.69x | -0.042 |
| formic_acid | 120 | 11 | 2.02x | 2.44x | -0.306 |
| fructose | 80 | 13 | 1.47x | 1.91x | +0.168 |
| fructose | 100 | 16 | 1.45x | 1.69x | +0.161 |
| fructose | 120 | 13 | 1.40x | 2.60x | +0.146 |
| glucose | 80 | 17 | 1.02x | 1.04x | -0.007 |
| glucose | 100 | 15 | 1.02x | 1.08x | -0.008 |
| glucose | 120 | 12 | 1.07x | 1.15x | -0.028 |
| glycine | 80 | 19 | 1.02x | 1.05x | -0.011 |
| glycine | 100 | 13 | 1.01x | 1.01x | -0.002 |
| glycine | 120 | 12 | 1.01x | 1.04x | -0.004 |
| methylglyoxal | 80 | 10 | 1.26x | 1.87x | +0.055 |
| methylglyoxal | 100 | 17 | 1.16x | 2.25x | -0.041 |
| methylglyoxal | 120 | 10 | 1.16x | 1.38x | +0.059 |

## The five worst-fitting responses (variant A)

| species | T (C) | median fold error | direction |
|---|---:|---:|---|
| formic_acid | 120 | 2.07x | model low |
| fructose | 80 | 1.47x | model high |
| fructose | 100 | 1.45x | model high |
| fructose | 120 | 1.40x | model high |
| methylglyoxal | 80 | 1.26x | model high |

These are misfits of the MEASURED constants, not of anything fitted here: the ten Martins steps were held at their published values. The largest, formic acid at 120 C, sits in the direction that Knol 2010's conflicting formic-acid Ea (84 +/- 14 vs Martins' 30 +/- 9) would move it, which is recorded as a standing conflict rather than resolved by refitting.

## Where the carbon goes (100 C, 120 min, variant A)

| reaction | integrated flux mmol/L | share |
|---|---:|---:|
| `r_schiff` | 55.57 | 23.3% |
| `r_amadori` | 53.04 | 22.2% |
| `r_glc_fru` | 31.04 | 13.0% |
| `r_ama_odg` | 19.14 | 8.0% |
| `r_odg_aa` | 19.01 | 8.0% |
| `r_ama_tdg` | 13.53 | 5.7% |
| `r_fru_glc` | 12.94 | 5.4% |
| `r_tdg_mel` | 10.34 | 4.3% |
| `r_ama_mgo` | 8.63 | 3.6% |
| `r_aa_frag` | 7.827 | 3.3% |
| `r_mgo_mel` | 5.214 | 2.2% |
| `r_tdg_fa` | 2.442 | 1.0% |
| `r_fru_acids` | 0.06237 | 0.0% |
| `r_glc_frag` | 0.0001893 | 0.0% |
| `r_fa_frag` | 4.07e-06 | 0.0% |

## Melanoidin C/N -- directional diagnostic (not fitted, not scored)

- measured (Brands 2002, glucose/casein 120 C): {'10': 4.01, '30': 4.1, '60': 4.22, 'unheated_casein_reference': 3.97} -- direction **increasing**
- predicted: {'0': None, '10': 8.415719442617828, '30': 9.128435368486835, '60': 9.938160988411417} -- direction **increasing**
- Brands' melanoidin is protein-bound (casein backbone, unheated reference C/N 3.97), so its level is set by the protein, not by the Maillard carbon. Only the SIGN of the change transfers.

## Cross-lab comparators carried but NOT operative

| quantity | value | source | why not operative |
|---|---|---|---|
| Ea, sugar isomerisation | 61.0 kJ/mol | Knol et al. 2010, Table 2 | Martins measures the same transformation on the same system this module integrates, in the module's own 80-120 C window, and measures BOTH directions (123 +/- 5 forward, 93 +/- 3 reverse). Knol's 61 +/- 8 is a NET lumped isomerisation over 120-200 C. CONFLICT: 61 +/- 8 does not overlap either Martins direction. Reported, not averaged. |
| Ea, acetic acid formation | 75.0 kJ/mol | Knol et al. 2010, Table 2 | AGREES with the operative Martins step 8 value (76 +/- 4) to within 1 kJ/mol, across two labs, two amines and two temperature windows. This is the strongest cross-lab corroboration in Module 4 and it costs nothing, so it is reported as such rather than merged. |
| Ea, formic acid formation | 84.0 kJ/mol | Knol et al. 2010, Table 2 | DIRECT CONFLICT with the operative Martins step 5 value (30 +/- 9). The intervals are 54 kJ/mol apart and do not overlap. Martins is operative because it is the same system in the same window; the conflict is a standing owner question, not something this module resolves. |
| Ea, condensation (Glc + Asn -> Schiff base) | 57.6 kJ/mol | Knol et al. 2005, Table 1 | CONFLICT flagged by the inventory itself (Z2 #17): 57.6 +/- 8.0 versus Martins' 96.8 +/- 2.8 for the same condensation, intervals non-overlapping. Different amine (asparagine vs glycine) and a different window. Owner question before either is called 'the' condensation barrier. |
| Ea, browning A420 (glucose + casein / + BSA, alkaline) | [164.0, 120.0, 130.0, 126.0, 92.0, 95.0] | Ajandouz et al. 2008, Table 3 p. 1250 | THE ALKALINE-pH WALL. Ajandouz licenses Ea transfer to pH 5-7 only for glucose-loss and amino-loss, and EXPLICITLY NOT FOR BROWNING, whose Ea fall 15-29% between pH 8.0 and 9.7 (inventory sec. C.13). Carried as an unvalidated prior with its measurement pH attached, exactly as FIT_HOLDOUT_DECLARATION.md sec.5 decision 3 requires. It is a third-lab CONTEXT for Martins' melanoidin Ea 95.2, not a referee for it. |
| Ea, beta-elimination (cysteine) | 123.0 kJ/mol | Zheng & Ho 1994, Tables I/III/V | Cysteine is not in the glucose/glycine trunk; this constant belongs to the sulfur module that plugs into this core. It is registered here so that module inherits it WITH its pH 9 label rather than re-deriving it. It also condemns the repo's shipped beta_elimination_dha Ea of 79.5 (43.5 kJ/mol below the measured aqueous barrier, with a silent NaN prefactor). GAP REPORTED: the declaration assigns 'Zheng 1994 Tables I/III/V (36 k + 8 Ea)' to FIT, but this single value is the ONLY Zheng number transcribed anywhere in the repository's dossiers. The other 43 are not available to be used. |

## Pre-registration of the hold-out evaluation

- variant B (melanoidin sink estimated from the REACTANT side only) is the pre-registered out-of-sample prediction of the Martins step-9 browning hold-out.
- variant A (Martins' own measured step-9 constant) is reported alongside it as a REPRODUCIBILITY check, and is explicitly NOT out-of-sample: Martins fitted that constant to the held-out response.
- No parameter may change after the hold-out is read. The hold-out scorer is a separate script that reads the frozen parameters from the JSON companion of this report.
