# Yu, Seow, Ong & Zhou 2018 — EXTRACTION (thermal control arm: 0.1 M D-glucose + 0.1 M glycine, pH 10.0 unbuffered, 70/80/90 C, 15-90 min; multi-response fit with pyrazine steps)
### A glucose-glycine multi-response model whose Table 1 prints Ea ± SE and a rate constant at 80 C for four pyrazine-forming steps, second order in [1-deoxyglucosone] x [glycine]; the ultrasound arm is ignored here.

**Source on disk:** `data/articles/yu2018.pdf` (44 pp., owner's download, 2026-09-08). This is the
**accepted manuscript** (Elsevier "Accepted Manuscript" cover, PII S0308-8146(18)31192-0, DOI
10.1016/j.foodchem.2018.07.053 printed on p. 1); the published version is Food Chemistry 269 (2018)
628-637 (per Zhou 2024, ref 27) and may differ in detail. Read from the text layer
(`scratchpad/articles/yu2018.txt`), which is clean for prose but breaks superscripts in Table 1 and
drops Equations 17-23. **Manuscript pp. 13-14 (equations), p. 35 (Table 1), Figure 1 and Figure 4
were rasterised at 110-130 dpi** (`scratchpad/img/yu-14.png`, `yu-15.png`, `yu-36.png`, `yu-41.png`,
`yu-44.png`); every Table 1 cell below was read from the raster. Figure axes were read for units only;
no data value was read from any figure.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetic study of high-intensity ultrasound-assisted Maillard reaction in a model system of D-glucose and glycine" |
| Authors | Hang Yu, Yi-Xin Seow, Peter K. C. Ong, Weibiao Zhou* (National University of Singapore; NUS Suzhou Research Institute) |
| Venue | Food Chemistry 269 (2018) 628-637. Received 3 March 2018, revised 5 July 2018, accepted 8 July 2018 |
| DOI | 10.1016/j.foodchem.2018.07.053 |
| Naming | DFG = N-(1-deoxy-D-fructos-1-yl)-glycine; 1-DG = 1-deoxyglucosone; 3-DG = 3-deoxyglucosone; MG = methylglyoxal; "2,3,5-tetramethylpyrazine" in Eq. 14 is 2,3,5,6-tetramethylpyrazine (Fig. 1, Table 2); "3,5-dimethyl-vinylpyrazine" = 3,5-dimethyl-2-vinylpyrazine; "Timethylpyrazine" in Table 2 = 2,3,5-trimethylpyrazine |
| Companions | Yu, Seow, Ong & Zhou 2018, npj Sci. Food 2, 2 (the SPME-GC-MS extraction / identification / quantification method — NOT re-printed here); Yu et al. 2016 IFSET (the ultrasonic reactor); Martins & van Boekel 2005 (the trunk scheme and the reparameterised Arrhenius form; already in the repo as `martins2005_extraction.md`) |

## 1. Why it matters

The trunk of the engine is glucose/glycine -> Amadori -> deoxyosones. This is the only
glucose-glycine multi-response study in the corpus that hangs **pyrazine formation steps on the
1-deoxyglucosone node with fitted rate constants and activation energies (with standard errors)**,
in water, in a plain batch heating. Its thermal arm therefore offers a candidate pre-registration
for a pyrazine lane: r = k_i [1-DG][Gly] with Ea for 2,5-dimethylpyrazine (99.8 ± 6.7 kJ/mol),
2,3,5-trimethylpyrazine (117.8 ± 18.0), 3-ethyl-2,5-dimethylpyrazine (103.3 ± 2.7) and
tetramethylpyrazine (104.1 ± 5.1), and k at 80 C for each. Two things stop it being directly
usable: **the time unit of k is never printed** (s per Eq. 17 and the figure axes, min per the
residence-time text), and **every concentration is figure-only**. The Ea are transferable as
model-conditional barriers; the k values are transferable only once the unit is resolved from the
published version or the npj method paper. Methylpyrazine and the parent pyrazine are NOT among the
compounds detected in this system.

## 2. Methods as they matter to a model (thermal arm)

- **Reactants.** "Equal molar (0.1 mol) of both reactants, D-glucose (>= 99.0%) and glycine (>= 99.0%)
  were dissolved in 1 L deionized water. The sample solution's pH was adjusted to 10.0 by adding
  sodium hydroxide." So **100 mmol/L glucose + 100 mmol/L glycine, pH 10.0, no buffer**.
- **Heating.** "Thermal MR was conducted by placing a reaction tube filled with 10 mL of the prepared
  sample into a shaking water bath ... at a preset temperature", 70, 80, 90 C; "The temperature of
  prepared sample increased rapidly to the pre-set point in less than one min." Residence times "from
  15 to 90 min at an interval of 15 min" (six points), then flash-cooled in ice water, foil-wrapped.
  Tube type and closure not stated.
- **pH drift.** "the pH value was measured after each trial, and the change was found to be within
  the range of 0.3 - 0.6" (downward, from acid formation); the authors treat < 1 unit as negligible
  after Martins & van Boekel 2005.
- **Replicates.** "All experiments were conducted in triplicates."
- **Analytics.** Glucose, fructose, formic and acetic acid: HPLC, C-610H column, 0.1 % H3PO4, 0.5
  mL/min, 30 C, RI + DAD 210 nm. Glycine: Waters AccQ-Tag derivatisation, HPLC-DAD 210 nm. DFG, 3-DG,
  1-DG, MG: o-phenylenediamine derivatisation (1 mL sample + 1 mL water + 2 mL 1 mol/L OPD in
  methanol, 24 h room temperature), Sunfire C18, ammonium acetate 20 mmol/L pH 3.50 / acetonitrile
  gradient, DAD (Martins, Marcelis & van Boekel 2003 method). Melanoidins: A470 with epsilon = 0.64 ±
  0.03 L mmol-1 cm-1 (Martins & van Boekel 2003, same system). Volatiles: SPME, DVB/CX/PDMS 50/30 µm
  fibre, 4 mL sample, 20 min extraction (temperature not stated), Agilent 7890A GC-MS; "The
  extraction, identification, and quantification methods were the same as those described in (Yu,
  Seow, Ong, & Zhou, 2018)" — **internal standard, calibration and LOD are not given in this paper**.
  Figure 4 plots pyrazines in µmol/L, so an absolute calibration exists in the npj paper.
- **Scheme (Figure 1, re-typed).** k1 glucose -> fructose; k2 fructose -> glucose; k3 glucose +
  glycine -> DFG; k4 DFG -> 3-DG (+ glycine); k5 3-DG -> formic acid; k6 DFG -> MG (+ glycine); k7 DFG
  -> 1-DG (+ glycine); k8 1-DG -> acetic acid; k9 1-DG + glycine -> melanoidins; k10 1-DG + glycine ->
  2,5-dimethylpyrazine; k11 -> 2,3,5-trimethylpyrazine; k12 -> 3-ethyl-2,5-dimethylpyrazine; k13 ->
  2,3,5,6-tetramethylpyrazine; k14 -> 3,5-dimethyl-2-vinylpyrazine; k15 -> 3,5-diethyl-2-methylpyrazine
  (k14, k15 ultrasound only). Rationale: "the concentrations of 1-DG were more than 25 times higher
  than other dicarbonyl compounds ... the generation of melanoidins and flavor compounds was
  established to start from 1-DG due to its majority and high reactivity."
- **Rate laws (Eqs 1-16, verbatim structure).** r_glucose = -k3[Glu][Gly] + (-k1 + k2)[Fru] [sic;
  see Flags]; r_glycine = -k3[Glu][Gly] + (k4 + k6 + k7)[DFG]; r_fructose = k1[Glu] - k2[Fru]; r_DFG =
  k3[Glu][Gly] - (k4 + k6 + k7)[DFG]; r_3-DG = k4[DFG] - k5[3-DG]; r_MG = k6[DFG]; r_1-DG = k7[DFG] -
  k_sum[1-DG][Gly]; r_formic = k5[3-DG]; r_acetic = k8[1-DG]; r_melanoidins = k9[1-DG][Gly];
  **r_pyrazine_i = k_i[1-DG][Gly], i = 10..15**. Concentrations in mmol/L. "the depletion rate
  constant of 1-DG was written as k_sum, which was not numerically equivalent to the sum of
  generation rate constants of the quantified final MRPs ... (i.e. k_sum != k9 + ... + k15)". So the
  pyrazine steps are **second order (1-DG x glycine)**, with glycine at ~100 mmol/L and nearly
  constant — in practice pseudo-first order in 1-DG (my remark).
- **Fitting.** Batch reactor: MATLAB dsolve; k by nlinfit at each temperature; then Ea and A by
  nlinfit with the **reparameterised Arrhenius (Eqs 20-23, from the raster):** k = A exp(-Y Ea),
  A = A0 exp(-Ea/(R T_av)), Y = (1/R)(1/T - 1/T_av), T_av = sum(T)/n. With T = 70, 80, 90 C,
  **T_av = 80 C = 353.15 K, so the printed "A" is k at 80 C**, not a pre-exponential factor. Models
  "were fitted at all three temperatures simultaneously". Goodness of fit (all pyrazines, both arms):
  relative RMSE "0.0013% to 4.5036%", R2 "0.8994 to 0.9999" — not given per step.

## 3. Tables re-typed

### Table 1. "Kinetic model parameters of Ea (kJ mol-1) and its corresponding A for steps 1 to 15 (in Fig. 1) in ultrasonic and thermal MRs"

Footnote: "Within each reaction step, significant differences in Ea and A are indicated by different
superscript letters (p < 0.05). The unit of A is the same as that of the corresponding k." Values
are mean ± (standard error of the fit, presumably; the paper does not say which). Superscript
letters a/b shown after the Ea. Read from the raster.

| step | what it is | Ea ultrasonic | A ultrasonic | Ea thermal | A thermal |
|---:|---|---|---|---|---|
| 1 | glucose -> fructose | 100.8 ± 6.2 a | 6.1e-5 ± 7.5e-6 | 84.2 ± 5.7 b | 6.8e-5 ± 3.4e-6 |
| 2 | fructose -> glucose | 103.6 ± 3.6 a | 2.1e-4 ± 6.0e-5 | 95.2 ± 3.0 b | 8.5e-4 ± 6.9e-4 |
| 3 | glucose + glycine -> DFG | 72.9 ± 4.1 a | 8.5e-5 ± 4.8e-6 | 64.8 ± 9.0 a | 9.3e-5 ± 4.7e-6 |
| 4 | DFG -> 3-DG | 51.4 ± 6.9 b | 4.3e-7 ± 1.0e-7 | 89.5 ± 5.4 a | 8.1e-7 ± 6.7e-8 |
| 5 | 3-DG -> formic acid | 29.1 ± 5.2 a | 1.3e-2 ± 2.4e-3 | 43.9 ± 10.3 a | 1.7e-4 ± 6.7e-6 |
| 6 | DFG -> MG | 51.8 ± 12.9 a | 2.3e-6 ± 1.1e-7 | 30.2 ± 8.4 a | 1.7e-6 ± 1.6e-7 |
| 7 | DFG -> 1-DG | 60.9 ± 9.7 b | 8.3e-5 ± 4.2e-5 | 105.6 ± 9.9 a | 3.4e-5 ± 2.9e-5 |
| 8 | 1-DG -> acetic acid | 85.9 ± 18.0 a | 2.5e-3 ± 4.1e-4 | 56.1 ± 12.1 a | 5.3e-4 ± 1.3e-5 |
| 9 | 1-DG + Gly -> melanoidins | 87.7 ± 12.2 b | 8.1e-3 ± 2.3e-3 | 108.7 ± 1.5 a | 1.3e-4 ± 5.9e-6 |
| **10** | **1-DG + Gly -> 2,5-dimethylpyrazine** | 88.7 ± 3.4 b | 1.7e-6 ± 1.3e-7 | **99.8 ± 6.7 a** | **1.0e-5 ± 1.4e-6** |
| **11** | **1-DG + Gly -> 2,3,5-trimethylpyrazine** | 89.6 ± 13.4 b | 7.2e-7 ± 2.6e-8 | **117.8 ± 18.0 a** | **2.8e-6 ± 4.4e-7** |
| **12** | **1-DG + Gly -> 3-ethyl-2,5-dimethylpyrazine** | 76.6 ± 8.3 b | 1.1e-7 ± 1.7e-8 | **103.3 ± 2.7 a** | **3.9e-8 ± 2.1e-8** |
| **13** | **1-DG + Gly -> 2,3,5,6-tetramethylpyrazine** | 82.8 ± 8.0 b | 2.3e-7 ± 5.4e-8 | **104.1 ± 5.1 a** | **6.2e-8 ± 9.1e-9** |
| 14 | 1-DG + Gly -> 3,5-dimethyl-2-vinylpyrazine | 63.2 ± 28.9 | 4.3e-9 ± 5.2e-10 | -- | -- |
| 15 | 1-DG + Gly -> 3,5-diethyl-2-methylpyrazine | 101.90 ± 10.41 | 2.36e-9 ± 1.30e-10 | -- | -- |

Text cross-checks (thermal): step 1 "84.2 ± 5.7", step 2 "95.2 ± 3.0", step 3 "64.8 ± 9.0", step 4
"89.5 ± 5.4", step 6 "30.2 ± 8.4", step 7 "105.6 ± 9.9" (abstract prints "105.5 ± 9.9" — Flag 8),
step 9 "108.7 ± 1.5", step 10 "99.8 ± 6.7". The Martins & van Boekel 2005 comparators quoted by the
authors: 3-DG 97.1 ± 1.7, MG 124.5 ± 4.7, 1-DG 107.3 ± 7.3 kJ/mol.

**Units of A / k.** First-order steps (1, 2, 4-8): time-1. Second-order steps (3, 9-15): L mmol-1
time-1 (concentrations are in mmol/L per Eqs 1-16). The time unit is not printed anywhere. Eq. 17
defines residence time tau in **s**, and Figures 2-4 have time axes in s (0-6000 s); the Methods
give residence times in **min**. Both readings for the pyrazine steps:

| step | k(80 C) printed | if per second: L mol-1 s-1 | if per minute: L mol-1 s-1 |
|---:|---|---:|---:|
| 10 | 1.0e-5 L mmol-1 t-1 | 1.0e-2 | 1.67e-4 |
| 11 | 2.8e-6 | 2.8e-3 | 4.67e-5 |
| 12 | 3.9e-8 | 3.9e-5 | 6.5e-7 |
| 13 | 6.2e-8 | 6.2e-5 | 1.03e-6 |
| 9 (melanoidins) | 1.3e-4 | 1.3e-1 | 2.17e-3 |

(1 L mmol-1 = 1000 L mol-1.) The choice cannot be settled from this manuscript because the
concentrations that would allow a back-calculation are figure-only; the ratio between steps is
unit-independent: k10 : k11 : k12 : k13 = 1 : 0.28 : 0.0039 : 0.0062 at 80 C.

**k at the other two temperatures (mine, from the printed A and Ea via the reparameterised form,
same unresolved time unit):** step 10: 3.7e-6 (70 C), 1.0e-5 (80 C), 2.6e-5 (90 C); step 11: 8.7e-7,
2.8e-6, 8.5e-6; step 12: 1.4e-8, 3.9e-8, 1.0e-7; step 13: 2.2e-8, 6.2e-8, 1.7e-7 L mmol-1 t-1.
Conventional pre-exponential A0 = A exp(Ea/(R x 353.15)): step 10 5.8e9; step 11 7.4e11; step 12
7.4e7; step 13 1.6e8 (same units as k).

### Table 2. "Summary of all volatile MRPs in the D-glucose and glycine model system detected by GC-MS in ultrasonic and thermal MRs samples at residence time of 90 min under 90 C" (thermal column only; "--" = not detected)

| RT (min) | compound | thermal MR | class |
|---:|---|---|---|
| 8.969 | 2,5-dimethylpyrazine | present | pyrazine |
| 9.193 | 2,3-dimethylpyrazine | present | pyrazine |
| 9.517 | 2,6-dimethylpyrazine | -- | pyrazine |
| 11.460 | 2,3,5-trimethylpyrazine | present | pyrazine |
| 11.901 | 2,5-dimethyl-3-propylpyrazine | -- | pyrazine |
| 13.484 | 3-ethyl-2,5-dimethylpyrazine | present | pyrazine |
| 13.626 | 2-ethyl-3,5-dimethylpyrazine | -- | pyrazine |
| 13.676 | 2,3,5,6-tetramethylpyrazine | present | pyrazine |
| 14.001 | 3,5-dimethyl-2-vinylpyrazine | -- | pyrazine |
| 15.523 | 3,5-diethyl-2-methylpyrazine | -- | pyrazine |
| 11.247 | 2,4,6-trimethylpyridine | present (thermal only) | pyridine |
| 10.892 | butyl amide | present | amine |
| 12.733 | butyl amine | -- | amine |
| 20.514 | acryl amide | present | amine |
| 20.519 | octyl amine | present | amine |
| 5.135 | 2,2,5,5-tetramethyltetrahydrofuran | present (thermal only) | furan |
| 9.781 | 4-methyl-2-heptanone | present | ketone |
| 12.211 | 2-ethyl-1-hexanol | present | alcohol |
| 6.119 | ethyl butyrate | present (thermal only) | ester |
| 7.508 | ethyl 3-methylbutanoate | present (thermal only) | ester |
| 10.658 | methyl 2-methylpropanoate | present | ester |
| 20.372 | butyl butanoate | present | ester |

No concentrations are printed in Table 2. **Absent from the list altogether: methylpyrazine,
pyrazine, 2-ethylpyrazine, furfural, HMF, furaneol** (the octyl amine, 2-ethyl-1-hexanol and ester
entries look like fibre/lab background; my remark). 2,3-dimethylpyrazine is detected but has no
kinetic step.

### Figures (FIGURE-ONLY)

Figure 2: glucose, glycine, fructose, DFG vs time (thermal panels B, D, F, H). Figure 3: 3-DG, 1-DG,
MG, formic acid, acetic acid, melanoidins (thermal B, D, F, H, J, L). Figure 4: 2,5-dimethylpyrazine,
2,3,5-trimethylpyrazine, 3-ethyl-2,5-dimethylpyrazine, tetramethylpyrazine (thermal B, D, F, H),
y-axis "Concentration (µmol/L)", x-axis "Reaction time (s)" 0-6000, three temperatures, error bars,
fitted curves. No data table exists in the manuscript; the published version's supplementary
material (if any) is not on disk.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): 2,5-dimethylpyrazine -> `2_5_dimethylpyrazine`;
2,3-dimethylpyrazine -> `2_3_dimethylpyrazine`; 2,6-dimethylpyrazine -> `2_6_dimethylpyrazine`;
2,3,5-trimethylpyrazine -> `trimethylpyrazine` (registry SMILES Cc1cnc(C)c(C)n1 is the 2,3,5 isomer);
2,3,5,6-tetramethylpyrazine -> `tetramethylpyrazine`; **3-ethyl-2,5-dimethylpyrazine -> not in
registry** (the registry's `2_ethyl_3_5_dimethylpyrazine`, CCc1nc(C)cnc1C, is a different isomer,
which this paper lists separately and finds only under ultrasound — do not conflate); methylpyrazine
(`methylpyrazine`) and the parent pyrazine: not detected here; glucose, fructose, glycine, DFG, 1-DG,
3-DG, methylglyoxal -> not in `compounds.yml` (`reaction_rules.yml` uses short species names such as
Glc, Gly, MGO).

Shared conditions: 100 mmol/L glucose + 100 mmol/L glycine, water, initial pH 10.0 (NaOH,
unbuffered; -0.3 to -0.6 by the end), 10 mL batch in a shaking water bath, 70/80/90 C, 15-90 min,
triplicates; multi-response fit with the Figure 1 scheme.

| quantity | value | unit | conditions | reaction order (model) | source location | evidence class |
|---|---|---|---|---|---|---|
| Ea, 1-DG + Gly -> 2,5-dimethylpyrazine (step 10) | 99.8 ± 6.7 | kJ/mol | 70-90 C, pH 10 | second order, [1-DG][Gly] | Table 1, thermal | measured_barrier (model-conditional) |
| Ea, -> 2,3,5-trimethylpyrazine (step 11) | 117.8 ± 18.0 | kJ/mol | same | same | Table 1 | measured_barrier |
| Ea, -> 3-ethyl-2,5-dimethylpyrazine (step 12) | 103.3 ± 2.7 | kJ/mol | same | same | Table 1 | measured_barrier |
| Ea, -> 2,3,5,6-tetramethylpyrazine (step 13) | 104.1 ± 5.1 | kJ/mol | same | same | Table 1 | measured_barrier |
| k(80 C), step 10 | 1.0e-5 ± 1.4e-6 | L mmol-1 per (s or min — unresolved) | same | second order | Table 1 ("A" = k at T_av = 80 C) | measured_rate, time unit unresolved |
| k(80 C), step 11 | 2.8e-6 ± 4.4e-7 | same | same | same | Table 1 | measured_rate, time unit unresolved |
| k(80 C), step 12 | 3.9e-8 ± 2.1e-8 | same | same | same | Table 1 | measured_rate, time unit unresolved (SE ~ 54 % of value) |
| k(80 C), step 13 | 6.2e-8 ± 9.1e-9 | same | same | same | Table 1 | measured_rate, time unit unresolved |
| k10 : k11 : k12 : k13 at 80 C | 1 : 0.28 : 0.0039 : 0.0062 | — | same | — | derived from Table 1 | within-study ratio (unit-free) |
| k9 (1-DG + Gly -> melanoidins) vs k10 at 80 C | 1.3e-4 / 1.0e-5 = 13 | — | same | — | derived | within-study ratio |
| Ea, trunk steps (thermal): glucose->fructose 84.2 ± 5.7; fructose->glucose 95.2 ± 3.0; Glc+Gly->DFG 64.8 ± 9.0; DFG->3-DG 89.5 ± 5.4; 3-DG->formic 43.9 ± 10.3; DFG->MG 30.2 ± 8.4; DFG->1-DG 105.6 ± 9.9; 1-DG->acetic 56.1 ± 12.1; 1-DG+Gly->melanoidins 108.7 ± 1.5 | as listed | kJ/mol | same | Fig. 1 orders | Table 1 | measured_barrier (context for the trunk; compare `martins2005_extraction.md`) |
| k(80 C), trunk steps (thermal) | Table 1 "A thermal" column | time-1 or L mmol-1 time-1, unresolved | same | — | Table 1 | measured_rate, time unit unresolved |
| presence/absence of volatiles at 90 C / 90 min | Table 2 | — | same | — | Table 2 | level_only (qualitative) |
| pH change over a run | -0.3 to -0.6 | pH units | all runs | — | text | level_only |
| all concentration-time data (glucose, glycine, fructose, DFG, 3-DG, 1-DG, MG, formic, acetic, melanoidins, four pyrazines) | — | mmol/L; pyrazines µmol/L | 70/80/90 C | — | Figures 2-4 | figure_only |

## 5. Flags

1. **Time unit of every k / A is not printed.** Eq. 17 and the figure axes say seconds; the Methods
   say minutes. A factor of 60 hangs on this. Until the published article or the npj method paper
   settles it, the k values are ratios only.
2. **"A" is k(T_av = 80 C)**, not the Arrhenius pre-exponential (Eqs 20-23 on the raster; the text
   layer drops them). Anyone reading Table 1 from the text layer would take 1.0e-5 as a prefactor and
   be off by ~15 orders of magnitude.
3. **Accepted manuscript.** Numbers could have changed in proof. The published pagination and any
   supplementary tables are not on disk.
4. **All concentrations are figure-only**; no raw table, no per-step R2/RMSE; the global ranges
   (RMSE 0.0013-4.5 %, R2 0.899-0.9999) cover both arms and all pyrazines.
5. **The second-order law is a modelling choice, not a determination.** Glycine is ~100 mmol/L and
   changes little over 90 min at 70-90 C, so [1-DG][Gly] and [1-DG] alone are not distinguishable;
   with the bimolecular form, k absorbs 1/[Gly] ~ 0.01 L/mmol. The 1-DG balance uses a separate
   k_sum (explicitly not the sum of k9-k15), so the scheme is not mass-conserving at the 1-DG node.
6. **pH 10, unbuffered, -0.3 to -0.6 drift**; far above the engine's usual 5-8 range; NaOH only.
7. **Quantification method not in this paper** (delegated to the npj Sci. Food companion): internal
   standard, calibration range, LOD and SPME temperature unknown. Figure 4 values are absolute
   (µmol/L) on the strength of that companion.
8. **Small internal inconsistencies.** Abstract: 1-DG thermal Ea "105.5 ± 9.9" vs Table 1 and text
   "105.6 ± 9.9". Eq. 1 prints "(-k1 + k2) x [Fru]" where the scheme requires -k1[Glu] + k2[Fru]
   (Eq. 3 has it right). Figure references "Fig. 3G, 4I and Fig. 3H, 4J" and "Fig. 3K and 4L" should
   read Fig. 3 throughout. Eq. 14 names "2,3,5-tetramethylpyrazine". The sentence "However,
   significant difference was observed between the two MRs regarding the concentration of MG" sits
   next to "no significant difference in Ea value for the generation of MG" — one of the two lacks a
   "no".
9. **Not detected / not modelled:** methylpyrazine and pyrazine (absent from Table 2, so this paper
   gives nothing for the methylpyrazine channel the brief asked about); 2,3-dimethylpyrazine
   detected but unmodelled; 2,6-dimethylpyrazine and 2-ethyl-3,5-dimethylpyrazine absent in the
   thermal arm at 90 C / 90 min.
10. **Temperature window 70-90 C** — below the engine's 100-145 C; the Ea are three-point fits (with
    SE from the simultaneous fit, not from replicate Arrhenius plots).
11. Step 12's k has SE 54 % of its value (3.9e-8 ± 2.1e-8); step 11's Ea has SE 15 %. Treat those two
    as loosely constrained.
