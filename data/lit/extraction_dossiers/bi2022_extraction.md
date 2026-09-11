# Bi 2022 — EXTRACTION (pea protein isolate 10 g/L in 0.01 M potassium phosphate pH 7.6, 37 C, 2 h; NON-COVALENT binding of (Z)-2-penten-1-ol, hexanal and (E)-2-octenal by static headspace GC/MS with a Klotz fit, fluorescence quenching at 25/30/37 C, CD, surface hydrophobicity, phase-ratio-variation partition coefficients and docking)

### THE PEA ROW THE BINDING TABLE DOES NOT HAVE: `REVERSIBLE_BINDING` in `src/kinetic_core/parameters_matrix.py` carries dairy, soy and beta-lactoglobulin constants and not one pea constant, while `data/species/protein_matrices.yml` charges a full `pea_isolate` site table — and this paper prints headspace-derived Klotz binding constants for pea protein on hexanal (684.46 ± 109.43 M^-1) and on an alkenal ((E)-2-octenal, 8207.72 ± 2223.82 M^-1), plus a complete phase-ratio-variation partition pair (matrix/gas and buffer/gas) from which the repository's own per-gram form can be computed.

**Source on disk:** `data/articles/1-s2.0-S0308814622010068-main.pdf` (10 pp., Food Chemistry 389 (2022) 133044).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/1-s2.0-S0308814622010068-main.txt`); **Tables 1 and 2 came through clean** and
are re-typed in full below. **Tables S1, S2, S3 and S4 are supplementary and are NOT on disk** —
they hold the hydrophobicity constants and the real-pea concentrations (S1), the surface
hydrophobicity values (S2), the PRV partition coefficients (S3) and the bond-disrupting-agent key
(S4). Several of the S3 and S1 numbers are repeated verbatim in the running text of sections 3.1.1,
3.5 and 3.8 and those repeats ARE typed here; anything that appears only in a supplementary table is
recorded as absent. Figures 1A (retention vs concentration), 1B-D (Klotz plots), 2A-C
(Stern-Volmer plots, which carry the Ksv values inside the panels), 3A-C (fluorescence decay),
4A-F (docking) and S1/S2 are images: **every Stern-Volmer constant Ksv and every effective
quenching constant Ka in this paper is figure-only or unprinted.** Repo status before this dossier:
Bi 2022 is **not** cited anywhere in `src/kinetic_core/parameters_matrix.py`,
`src/kinetic_core/matrix_sites.py` or `data/species/protein_matrices.yml`, and has no extraction
dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Non-covalent interactions of selected flavors with pea protein: Role of molecular structure of flavor compounds" |
| Authors | Shuang Bi (BTBU), Xin Pan, Wentao Zhang, Zhuo Ma, Fei Lao, Qun Shen, Jihong Wu (corresponding, wjhcau@hotmail.com) — College of Food Science and Nutritional Engineering, China Agricultural University, Beijing 100083 |
| Venue | Food Chemistry 389 (2022) 133044. Received 20 January 2022, revised 27 March 2022, accepted 20 April 2022, online 22 April 2022 |
| DOI | 10.1016/j.foodchem.2022.133044 |
| Funding | National Key R&D Plan, 13th Five-Year Plan of China (2017YFD0401202) |
| The three ligands | **(Z)-2-penten-1-ol** (C5 alcohol, PubChem CID 5364919), **hexanal** (C6 n-alkanal, CID 6184), **(E)-2-octenal** (C8 alpha,beta-unsaturated aldehyde, CID 5283324). All three were identified as key pea aroma compounds by the same group (Bi et al. 2020, JAFC 68:2718-2727) |
| Protein | pea protein **extracted in the authors' own laboratory** from fresh peas (Dongfangliang, Shanxi) by hexane defatting -> alkaline extraction at pH 9.0 -> isoelectric precipitation at pH 4.5 -> resuspension to pH 7.0 -> freeze-drying. **Protein content is never stated.** The docking used the 11S crystal structure, PDB 3KSC |
| Naming | "binding percentage" = headspace depletion by Eq. 1; "K" = the Klotz/Scatchard binding constant in M^-1; "nK" = sites x constant; "Ksv" = Stern-Volmer quenching constant; "Ka" = effective quenching constant from the modified Stern-Volmer equation; "R (%)" = PRV percentage retention by Eq. 13 |
| Companions on disk | `guo2020_extraction.md` (soy isolate, the same Klotz/headspace family — cited by this paper as Guo et al. 2019/2020), `barallatperez2024_extraction.md` (lupin), `bornhorst2017_extraction.md` / `bornhorst2017b_extraction.md` (whey) |

## 1. Why it matters

**These numbers are NON-COVALENT throughout, and the paper says so in its title, its abstract and
its conclusion ("the interactions of pea protein with the three flavor compounds were reversible
and involved hydrophobic interactions, as well as hydrogen bonds").** They therefore belong to the
**matrix-retention layer, `src/kinetic_core/parameters_matrix.py`**, in the `REVERSIBLE_BINDING`
tuple — NOT to the covalent binding layer in `src/kinetic_core/matrix_sites.py`, whose
`BINDING_CLASSES` table is an aldehyde-to-lysine-amine and thiol-to-disulfide channel with rate
constants in M^-1 s^-1. Nothing in this paper is a rate: every constant here is an equilibrium
constant or a partition coefficient.

**The gap it fills.** `REVERSIBLE_BINDING` currently holds fifteen rows over four matrices —
`skim_milk` (Meynier 2002), `caseinate_1pct` (Leksrisompong 2010), `soy_protein` (Damodaran 1981 and
Arai 1970) and `beta_lactoglobulin` (Andriot 2000) — and **no pea row at all**, while
`data/species/protein_matrices.yml` carries a fully populated `pea_isolate` entry (free thiol
0.0159, disulfide 0.0257, amine 0.47 mmol per gram of protein, with bands). So the repository can
today charge covalent pea sites and cannot state a pea non-covalent constant. This paper is the
matching half.

**Where it lands in the existing method taxonomy.** The registry's `method` field is first-class
because k2 sec. B.3 measured a 35x aldehyde gap between headspace-depletion and dialysis
determinations, and `binding_constant_for` refuses to cross that boundary for aldehydes. **Every
constant in this paper is on the headspace side of that boundary** — the Klotz K is fitted to GC/MS
headspace peak areas (Eqs. 1-4), and the PRV K is a static-headspace partition coefficient (Eqs.
9-13). Its correct `method` values are therefore `headspace_depletion` (the Klotz K, matching
Andriot 2000's beta-lactoglobulin rows) and `static_headspace_partition` (the PRV K, matching
Meynier and Leksrisompong). It must NOT be pooled with the Damodaran soy dialysis aldehyde value
`kg_nonanal_soy` (4.38e-2 L/g), which was measured with 2-mercaptoethanol precisely to suppress the
cysteine-aldehyde chemistry a headspace determination counts.

**The structural trend the title promises, and what it does to the unsaturation penalty.** The
paper's whole point is the ladder **(E)-2-octenal > hexanal > (Z)-2-penten-1-ol**, reproduced by
five independent measurements (binding percentage, Klotz n, Klotz K, Ksv, PRV retention). The
alkenal/alkanal contrast on the Klotz K is **8207.72 / 684.46 = 12.0x (mine)**, and on the PRV
matrix/gas partition **2203.85 / 116.37 = 18.9x (mine)**. The registry's
`ALPHA_BETA_UNSATURATION` term is currently fitted on two FIT-row observations, one of which —
Meynier's t-2-hexenal / skim-milk row, `kg_t_2_hexenal_dairy` — is explicitly **quarantined as a
binding constant** because Amendment 6 ruling 4 sizes 22-33 % of it as irreversible chemistry. Bi's
(E)-2-octenal row has the **same structural hazard and no chemical control**: it is a Michael
acceptor, held 2 h at 37 C against a protein carrying 0.47 mmol lysine amine per gram, with the
bound fraction inferred *by disappearance from the headspace*. The paper never demonstrates
reversibility for (E)-2-octenal by an independent route (see Flags 2). **It should enter as a FIT
row for the unsaturation contrast and be quarantined as a reversible binding constant on exactly
the Meynier precedent.** The hexanal and (Z)-2-penten-1-ol rows carry no such hazard.

**A caution the registry already owns.** The parameters module refuses "any shipped matrix term
that is a monotone function of log P" (k4b hold-out guard #4). This paper argues its ladder *is*
hydrophobicity-ordered ("the retention of the three flavor compounds by pea protein was consistent
with their hydrophobicities"), but the hydrophobicity constants live in Table S1, which is not on
disk, so the correlation cannot be checked here and **nothing in this dossier licenses a log P
term**. The ladder is also confounded: chain length, functional group and unsaturation all move
together across the three compounds (Flags 3).

What this paper does NOT give the repository: any temperature above 37 C; any protein content for
the isolate it made; any molar mass for pea protein (so the Klotz K in M^-1 cannot be put on the
registry's per-gram basis without an assumption — Flags 1); any covalent adduct measurement; any
sensory measurement; any pH other than 7.6; any measurement on a commercial isolate.

## 2. Methods as they matter to a model

- **The protein.** Made in-house, not bought. Fresh peas milled to flour, defatted with 3 volumes
  hexane (2 h stirring, repeated twice, centrifuged 5000 x g / 4 C / 20 min); defatted flour
  dispersed in deionised water 1:10 w/v, **pH raised to 9.0 with 2.0 M NaOH**, stirred 1.5 h with pH
  monitored every 10 min; centrifuged 8000 x g / 4 C / 25 min; supernatant taken to **pH 4.5 with
  2 M HCl**; centrifuged 5000 x g / 4 C / 20 min; washed twice with deionised water; **resuspended
  to pH 7.0**; frozen at -80 C for 24 h and freeze-dried at -40 C for 48 h. The absence of the three
  target compounds in the finished protein was confirmed by GC/MS. **No protein content (N x 6.25),
  no purity, no molar mass and no thiol or amine assay is reported anywhere in the paper.** This is
  an alkaline-extracted, isoelectric pea isolate — the same preparation family as Gao 2020's, which
  is what `protein_matrices.yml` uses for the pea site densities.
- **Buffer and pH.** 0.01 M (10 mM) potassium phosphate, **pH 7.6**, throughout. A 2 % (w/v) pea
  protein stock was stirred 1 h at room temperature. **Methanol is present**: the flavour stock
  solutions were made at 10 mM in buffer with methanol added for dissolution at a **methanol:buffer
  ratio of 1:19**, and 0.5 mL of that stock goes into a 2 mL sample, so the assay carries roughly
  **1.25 % v/v methanol (mine)** — an unremarked co-solvent that competes for hydrophobic sites.
- **Loading.** Headspace samples: 1 mL of 2 % (w/v) protein + 0.5 mL buffer + 0.5 mL diluted flavour
  stock = 2 mL, so the **final protein concentration is 1 % w/v = 10 g/L**, stated in the text and
  repeated in the Fig. 1A caption as "10 mg/mL". Final flavour concentrations **0.05, 0.1, 0.15,
  0.2, 0.25, 0.5, 1, 1.5, 2 and 2.5 mM**. Loaded into a **20 mL vial**, sealed with a PTFE/silicone
  septum and magnetic crimp cap, held at **37 C for 2 h**.
- **Method 1 — static headspace GC/MS, binding percentage (what it measures).** Binding % =
  (1 - peak area with protein / peak area without protein) x 100, Eq. 1, after Wang & Arntfield
  2015a. The control is 1.5 mL buffer + 0.5 mL flavour stock, i.e. protein replaced by buffer at the
  same total volume. **This measures the DISAPPEARANCE OF THE COMPOUND FROM THE HEADSPACE and
  nothing else.** It cannot distinguish reversible partitioning into the protein phase from
  irreversible covalent capture; that distinction is made in this paper only by the bond-disrupting
  agents of section 3.6, and only qualitatively. GC/MS: Agilent 5975 MSD, HP-WAX 30 m x 0.25 mm x
  0.25 um, injector 250 C splitless, 1 mL headspace pumped by syringe, oven 45 C -> 240 C at
  15 C/min, helium 1.0 mL/min, selected-ion monitoring, EI 70 eV, ion source 230 C, auxiliary heater
  250 C, quadrupole 150 C.
- **Method 2 — the Klotz/Scatchard fit (what the constant is).** v/[L] = nK - vK (Eq. 2, Scatchard),
  in double-reciprocal Klotz form 1/v = 1/n + 1/(nK[L]) (Eq. 3), with v = ([L]t - [L])/Cp (Eq. 4),
  where v is **moles of flavour bound per mole of pea protein**, [L] is the free flavour
  concentration in mol/L taken from the GC peak areas, [L]t the total, Cp the pea protein
  concentration, n the number of binding sites and K the binding constant in M^-1. **K is therefore
  a per-MOLE-OF-PROTEIN constant whose numerical value depends entirely on the molar mass used for
  Cp, and that molar mass is never printed** (Flags 1). Standard errors of n, K and nK were
  propagated from the SEs of the Klotz intercept and slope by Eqs. 5-7. Fitted at **37 C**, on the
  **high-concentration branch (0.5-2.5 mM)** — the n values in Table 1 are stated in the text as
  "at high flavor compound concentrations".
- **Method 3 — fluorescence quenching (a different object entirely).** FS5 spectrofluorimeter
  (Edinburgh), excitation fixed at **290 nm**, emission 300-450 nm, slits 2 nm. 500 uL of 2 mg/mL
  pea protein + 250 uL flavour + 4250 uL buffer, so the **final protein is 0.2 mg/mL** and the
  flavour runs 0 to 2.5 mM in 0.5 mM steps. Water bath **25, 30 and 37 C for 2 h**. Stern-Volmer
  F0/F = 1 + Ksv[Q] (Eq. 8) gives Ksv; the modified form F0/dF = 1/(fa Ka [Q]) + 1/fa (Eq. 14) gives
  the effective quenching constant Ka, and the van 't Hoff pair (Eqs. 15-16) turns Ka's temperature
  dependence into dH, dS, dG. **A quenching constant is not a headspace binding constant**: it
  measures the perturbation of tryptophan emission at a 50x lower protein loading (0.2 mg/mL against
  10 mg/mL), and its magnitude and its temperature sign are properties of the fluorophore's
  environment. This paper's own two answers differ by two orders of magnitude (Flags 4).
- **Method 4 — time-resolved fluorescence.** Fluorolog-3-2ultrafast (HORIBA) with a DeltaDiode
  295 nm laser, excitation 290 nm, emission 330 nm, 37 C, at 0 / 1 / 2 mM flavour. **Lifetimes are
  reported only as curves (Fig. 3); no lifetime in ns is printed.**
- **Method 5 — circular dichroism.** Chirascan, far-UV 200-240 nm, 1 mm path, 200 nm/min, 0.1 nm
  resolution, five scans averaged; protein fixed at **0.2 mg/mL**, flavour 0-2.5 mM in 0.5 mM steps,
  37 C, samples centrifuged 10 000 x g / 20 min / 20 C first. Secondary structure by CDPro
  (CONTIN/LL).
- **Method 6 — surface hydrophobicity.** Kato & Nakai 1980 as modified: protein 0.1-0.4 mg/mL,
  4 mL supernatant + 40 uL ANS (8.0 mM in the same buffer), excitation 390 nm, emission 470 nm,
  slits 1 nm; H0 is the initial slope of fluorescence against protein concentration. Reported as a
  percentage of the flavour-free control (= 100 %).
- **Method 7 — bond-disrupting agents (the reversibility test, such as it is).** 4 M guanidine
  hydrochloride (weakens hydrophobic interactions AND inhibits hydrogen bonding) and 80 % propylene
  glycol (disrupts hydrophobic interactions but PROMOTES hydrogen bonding), after Ustunol 1992 and
  Wang & Arntfield 2016b. 10 mL of 2 % protein + 5 mL agent, stirred 1 h (liquid A); the buffer
  control is liquid B; 1.5 mL of A or B + 0.5 mL flavour, vortexed. Read as headspace peak area.
  **This is the only evidence in the paper that anything released again.**
- **Method 8 — phase ratio variation (PRV), after Tehrany 2007.** 1/A = a + b*beta (Eq. 9) with
  a = K/(fi Cin) (Eq. 10) and b = 1/(fi Cin) (Eq. 11), so **K = a/b** (Eq. 12) is the **matrix/gas
  partition coefficient** — note the direction: high K means the compound prefers the liquid.
  Retention R = (1 - K1/K2) x 100 (Eq. 13) with K1 the buffer/gas coefficient and K2 the
  matrix/gas coefficient; a positive R is retention, a negative R is release. Volumes of 50, 75,
  100, 200, 500, 1000, 2000, 3000, 4000 and 5000 uL in a 20 mL vial give beta = 399, 265.67, 199,
  99, 39, 19, 9, 5.7, 4 and 3. All other conditions as section 2.3.1 (so 37 C, 2 h, pH 7.6). **The
  flavour concentrations here are the real-pea levels, not the assay ladder**: (Z)-2-penten-1-ol
  901 ± 119 ug/kg, hexanal 1260 ± 114 ug/kg, (E)-2-octenal 22.9 ± 3.54 ug/kg (Bi et al. 2020, via
  Table S1, quoted in section 3.8).
- **Method 9 — docking.** AutoDock Vina 1.1.2 against pea 11S, PDB 3KSC, waters and ligands
  stripped, grid box over the whole protein, visualised in PyMOL. **In silico; no binding energy in
  kcal/mol is printed, only contact distances in angstrom.**
- **Replication and statistics.** All experiments in triplicate; mean ± SD; one-way ANOVA with
  Duncan's test at p < 0.05; Origin 9.1.

## 3. Tables re-typed

### Table 1 (p. 5). "Binding parameters characterizing binding pea proteins to different selected flavors at 37 °C: (Z)-2-penten-1-ol, hexanal, and (E)-2-octenal."

| Flavor compounds | n | K (M^-1) | nK (M^-1) |
|---|---|---|---|
| (Z)-2-penten-1-ol | 2.55 ± 1.34 | 438.68 ± 259.27 | 1119.90 ± 273.72 |
| hexanal | 4.84 ± 0.36 | 684.46 ± 109.43 | 3313.23 ± 308.45 |
| (E)-2-octenal | 26.95 ± 2.92 | 8207.72 ± 2223.82 | 221238.94 ± 39964.95 |

Header note: the nK column is printed with the unit `(M^-1)`, which is the unit of K, not of nK;
nK carries the same unit only because n is dimensionless (moles per mole). The table came through
the text layer clean; the (E)-2-octenal row wrapped across two lines in the extraction and has been
rejoined, and the values are unambiguous.

**Internal check (mine): the columns do not close.** n x K should equal nK. It does not:

| compound | n x K (mine) | printed nK | ratio |
|---|---|---|---|
| (Z)-2-penten-1-ol | 1118.6 | 1119.90 | 1.001 |
| hexanal | 3312.8 | 3313.23 | 1.000 |
| (E)-2-octenal | 221 197 | 221 238.94 | 1.000 |

They close to 4 significant figures for all three. (Recorded because it is the only arithmetic check
available on this table, and it passes.)

### Table 2 (p. 6). "Thermodynamic parameters of the pea protein with different selected flavors at 37 °C: (Z)-2-penten-1-ol, hexanal, and (E)-2-octenal."

| Flavor compounds | T (K) | ΔH (kJ mol^-1) | ΔS (kJ mol^-1 K^-1) | ΔG (kJ mol^-1) |
|---|---:|---:|---:|---:|
| (Z)-2-penten-1-ol | 298 | 17.944 | 0.090 | −8.876 |
| | 303 | 17.944 | 0.096 | −11.144 |
| | 310 | 17.944 | 0.095 | −11.506 |
| hexanal | 298 | −150.916 | −0.468 | −11.452 |
| | 303 | −150.916 | −0.466 | −9.718 |
| | 310 | −150.916 | −0.466 | −6.456 |
| (E)-2-octenal | 298 | 105.728 | 0.383 | −8.356 |
| | 303 | 105.728 | 0.396 | −14.210 |
| | 310 | 105.728 | 0.396 | −16.982 |

The table title says "at 37 °C" but the table itself spans 298, 303 and 310 K; 37 C = 310.15 K is
the last row of each block. ΔH is constant down each block by construction (the van 't Hoff
assumption of Eq. 15). Note that ΔS is **not** constant down the (Z)-2-penten-1-ol and (E)-2-octenal
blocks even though ΔH is, which is internally odd (Flags 5).

**Internal check (mine): ΔG = ΔH − TΔS.** (Z)-2-penten-1-ol: 17.944 − 298(0.090) = −8.876 ✓;
17.944 − 303(0.096) = −11.144 ✓; 17.944 − 310(0.095) = −11.506 ✓. hexanal:
−150.916 + 298(0.468) = −11.452 ✓; −150.916 + 303(0.466) = −9.718 ✓; −150.916 + 310(0.466) =
−6.456 ✓. (E)-2-octenal: 105.728 − 298(0.383) = −8.406 against a printed −8.356;
105.728 − 303(0.396) = −14.260 against −14.210; 105.728 − 310(0.396) = −17.032 against −16.982 —
**a constant 0.050 kJ/mol offset in all three (E)-2-octenal rows**, consistent with ΔS being rounded
for display (ΔS ≈ 0.38283 and 0.39584 reproduce the printed ΔG). No number is wrong; the display is
rounded.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| retention span across all three compounds and all concentrations | **0.70 % to 97.35 %** | §3.1.1 |
| (E)-2-octenal retention span | **80.52 % to 97.35 %** | §3.1.1 |
| (Z)-2-penten-1-ol retention span | **0.70 % to 11.69 %** | §3.1.1 |
| hexanal retention span | **not printed** (intermediate; Fig. 1A only) | §3.1.1 |
| concentration branches | retention RISES over 0.05-0.25 mM, FALLS over 0.5-2.5 mM, and the high branch sits below the low branch | §3.1.1 |
| the reversal explained | at 2.5 mM the flavour "exceeded the binding capacity of pea protein" | §3.1.1 |
| alpha-helix, no flavour | **17.9 %** | §3.4 |
| alpha-helix with 1 mM (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **15.9 % / 15.4 % / 14.4 %** | §3.4 |
| beta-sheet and random coil | "increased" (no value printed) | §3.4 |
| surface hydrophobicity, 1 mM (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **98.38 % / 98.13 % / 89.89 %** of the flavour-free control | §3.5 |
| surface hydrophobicity, 2 mM (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **94.59 % / 96.92 % / 89.26 %** | §3.5 |
| (E)-2-octenal peak area, propylene glycol effect | rose **1.29 x 10^4 -> 1.55 x 10^4** (hydrophobic interaction confirmed) | §3.6, Fig. S1C |
| hexanal peak area, propylene glycol effect | fell **3.64 x 10^7 -> 2.93 x 10^7** (hydrogen bonding confirmed) | §3.6, Fig. S1B |
| (Z)-2-penten-1-ol, propylene glycol effect | "did not change significantly" -> both bond types present | §3.6, Fig. S1A |
| guanidine hydrochloride | peak area "increased significantly for all three systems" | §3.6 |
| **PRV: K matrix/gas, (Z)-2-penten-1-ol, pea protein / buffer** | **2065.44 / 1460.49** | §3.8 (from Table S3) |
| **PRV: K matrix/gas, hexanal, pea protein / buffer** | **116.37 / 32.90** | §3.8 |
| **PRV: K matrix/gas, (E)-2-octenal, pea protein / buffer** | **2203.85 / 455.93** | §3.8 |
| **PRV retention R, (Z)-2-penten-1-ol / hexanal / (E)-2-octenal** | **29.29 % / 71.72 % / 79.31 %** | §3.8 |
| concentration in peas, (Z)-2-penten-1-ol | 901 ± 119 ug kg^-1 | §3.8 (from Table S1, originally Bi 2020) |
| concentration in peas, hexanal | 1260 ± 114 ug kg^-1 | §3.8 |
| concentration in peas, (E)-2-octenal | 22.9 ± 3.54 ug kg^-1 | §3.8 |
| docking, (Z)-2-penten-1-ol | 2 H-bonds to **Arg147** at 2.98 and 2.87 Å; hydrophobic contact with **Lys357** at 3.65 and 3.30 Å | §3.7 |
| docking, hexanal | hydrophobic contacts **Leu12** 3.59 Å and **Lys39** 3.29 Å; 1 H-bond to **Gln40** at 2.89 Å | §3.7 |
| docking, (E)-2-octenal | hydrophobic **Lys482** at 3.68 / 3.77 / 3.69 Å and **Thr449** at 3.64 Å; H-bonds to **Asn336, Lys482, Phe483, Leu484** at 3.52 / 3.56 / 3.45 / 3.20 Å | §3.7 |
| quenching mechanism | Ksv RISES with T for (Z)-2-penten-1-ol and (E)-2-octenal (**dynamic**); Ksv FALLS with T for hexanal (**static**) | §3.2 |
| Ka temperature trend | Ka falls with T for hexanal; rises for (Z)-2-penten-1-ol and (E)-2-octenal | §3.2 |
| interaction assignment | hydrophobic for (Z)-2-penten-1-ol and (E)-2-octenal (ΔH > 0, ΔS > 0); van der Waals / hydrogen bonds for hexanal (ΔH < 0, ΔS < 0), by Ross & Subramanian 1981 | §3.2 |

**Ksv is FIGURE-ONLY.** The Stern-Volmer constants are printed inside the panels of Fig. 2A-C
(three compounds x three temperatures = nine values) and appear nowhere in the text or a table. Per
house rule they are not typed as numbers. The same holds for every Ka: the text describes only their
temperature trend. **Every retention percentage in Fig. 1A other than the four span endpoints quoted
above is figure-only**, as are the fluorescence decay curves (Fig. 3), the CD spectra (Fig. S2), the
Klotz plots themselves (Fig. 1B-D) and the bond-disrupting peak areas other than the four quoted.

### Arithmetic on the printed constants (all mine)

**1. The structural ladder, as ratios.** On the Klotz K: (E)-2-octenal / hexanal = 8207.72/684.46 =
**12.0x**; hexanal / (Z)-2-penten-1-ol = 684.46/438.68 = **1.56x**; (E)-2-octenal /
(Z)-2-penten-1-ol = **18.7x**. On nK (which is the quantity that actually scales a bound amount):
221238.94/3313.23 = **66.8x** and 3313.23/1119.90 = **2.96x**. On the number of sites n:
26.95/4.84 = **5.57x** and 4.84/2.55 = **1.90x**. **The alkenal's advantage is carried more by n
than by K** — it recruits 5.6x more sites at 12x higher affinity — which is what an irreversible
Michael addition to a protein carrying many lysines would also look like (Flags 2).

**2. The PRV pair converted to the registry's per-gram form.** The registry stores
K_g = (K_aw,water / K_aw,matrix − 1) / protein_g_per_L, where K_aw is an AIR/matrix coefficient.
This paper's K is MATRIX/gas, the reciprocal, so K_aw,water/K_aw,matrix = K2/K1 here (matrix over
buffer). Taking the protein loading as **10 g/L** (Flags 6 — it is stated for the headspace assay
and the PRV section says "all the conditions were the same", but it is not restated in §2.10):

| compound | K2/K1 (mine) | K_g = (ratio − 1)/10, L/g (mine) | the registry's nearest existing row |
|---|---:|---:|---|
| (Z)-2-penten-1-ol | 2065.44/1460.49 = **1.414** | **4.14e-2** | none (no alcohol in the table) |
| hexanal | 116.37/32.90 = **3.537** | **2.54e-1** | `kg_hexanal_dairy` 1.151e-2 L/g (skim milk, static headspace, 30 C) |
| (E)-2-octenal | 2203.85/455.93 = **4.834** | **3.83e-1** | `kg_t_2_hexenal_dairy` 1.734e-1 L/g (quarantined) |

**Pea protein's hexanal constant computed this way is 22x Meynier's skim-milk value** and
**173x Damodaran's dialysis-derived denatured-soy value** (`kg_hexanal_soy_denatured`, 1.47e-3
L/g). The second gap is the method gap k2 sec. B.3 already names and it is here five times wider
than the 35x it recorded; the first is protein-to-protein and matrix-to-matrix at the same method
family. **Neither ratio should be read as a pea-vs-dairy fact until the loading in §2.10 is
confirmed** (Flags 6) — the whole K_g scale is inversely proportional to it.

**3. Cross-check: the two headspace methods do not agree, and the paper says they nearly do.** The
PRV retention (29.29 / 71.72 / 79.31 %) is compared in §3.8 against the §3.1.1 binding percentages
and called a match "with a slight difference". For (E)-2-octenal the §3.1.1 span is 80.52-97.35 %
against a PRV 79.31 % — the PRV number is **below the entire §3.1.1 span**. For (Z)-2-penten-1-ol
the §3.1.1 span is 0.70-11.69 % against a PRV **29.29 %** — the PRV number is **2.5x above the top
of the entire §3.1.1 span**. The authors attribute this to the different (real-pea, far lower)
concentrations used in the PRV run. That attribution is plausible and untested, and it means **the
two headspace determinations in this one paper differ by a factor of 2.5 on the same compound in the
same buffer at the same temperature** (mine). Any single number taken from this paper inherits that
spread.

**4. The Klotz K in M^-1 cannot be put on a per-gram basis from this paper.** The registry's
soy and beta-lactoglobulin rows were built as n*K/MW with the molar mass either stated by the source
(Damodaran, 100 000 g/mol) or recovered by arithmetic (Andriot, 36 800 g/mol dimer). **Bi prints no
molar mass.** The docking used pea 11S (legumin), whose hexamer is commonly taken near 360 000
g/mol and whose acidic+basic subunit pair is near 60 000. The arithmetic, entirely conditional:

| compound | nK, M^-1 | K_g at MW 360 000 (mine) | K_g at MW 60 000 (mine) |
|---|---:|---:|---:|
| (Z)-2-penten-1-ol | 1119.90 | 3.11e-3 L/g | 1.87e-2 L/g |
| hexanal | 3313.23 | 9.20e-3 L/g | 5.52e-2 L/g |
| (E)-2-octenal | 221238.94 | 6.15e-1 L/g | 3.69 L/g |

The 360 kDa reading puts hexanal at 9.20e-3 L/g, within **1.25x** of Meynier's skim-milk
1.151e-2 L/g — a striking agreement, and **it is an artefact of a molar mass I chose**. It is
recorded so that a later reader does not rediscover it and mistake it for evidence. **Do not ship
any of these six numbers.** The defensible route to a pea per-gram constant from this paper is the
PRV pair in item 2, which needs no molar mass at all.

**5. The quenching constants recovered from ΔG (mine, and they disagree with the Klotz K by 10-60x).**
ΔG = −RT ln K (Eq. 16), so K = exp(−ΔG/RT) with R = 8.314 J mol^-1 K^-1. This recovers the effective
quenching constant Ka that the paper computed but never printed:

| compound | Ka at 298 K | Ka at 303 K | Ka at 310 K | Klotz K at 310 K (Table 1) | Klotz / Ka at 310 K |
|---|---:|---:|---:|---:|---:|
| (Z)-2-penten-1-ol | 36.0 M^-1 | 83.4 M^-1 | **86.9 M^-1** | 438.68 | **5.0x** |
| hexanal | 101.7 M^-1 | 47.4 M^-1 | **12.2 M^-1** | 684.46 | **56x** |
| (E)-2-octenal | 29.2 M^-1 | 281.6 M^-1 | **728.5 M^-1** | 8207.72 | **11x** |

**A headspace-derived constant and a fluorescence-derived one are not the same object, and this
table is the corpus's cleanest demonstration of it**: same laboratory, same protein preparation,
same buffer, same pH, same 37 C, same three ligands — and the two determinations differ by 5x, 56x
and 11x, in a direction that is not even consistent across the ladder. They also disagree on the
ORDER: on Ka at 310 K the ladder is (E)-2-octenal > (Z)-2-penten-1-ol > hexanal, i.e. **hexanal and
the alcohol swap places** against every headspace measurement in the paper. `binding_constant_for`'s
refusal to cross the method boundary is vindicated here on a single paper's own internal data.

**6. Ka is also 50x lower in protein loading than K.** The quenching runs are at 0.2 mg/mL; the
headspace runs at 10 mg/mL. Any comparison between the two carries a 50x aggregation-state
difference on top of the method difference.

**7. What the structure data says about the mechanism.** Alpha-helix falls 17.9 % -> 14.4 % with
1 mM (E)-2-octenal, an absolute loss of **3.5 percentage points and a relative loss of 19.6 %
(mine)**; hexanal and the alcohol lose 2.5 and 2.0 points. Surface hydrophobicity falls only
**10.1 %** at 1 mM (E)-2-octenal and **1.6-1.9 %** for the other two. Both trends order with the
binding ladder. Neither is a binding constant and neither should be modelled.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `hexanal` is keyed (id `hexanal`).
`(E)-2-octenal` is keyed (id `e_2_octenal`). **`(Z)-2-penten-1-ol` is NOT in the registry** and has
no `COMPOUND_STRUCTURE` entry in `parameters_matrix.py` either — it would be the first alcohol in
that structural table (the nearest members are `1_hexanol` and `1_octen_3_ol` in `compounds.yml`,
neither of which is in `COMPOUND_STRUCTURE`). `COMPOUND_STRUCTURE` already carries `hexanal`
(n_alkanal) and `t_2_octenal` (alkenal, alpha,beta-unsaturated carbonyl = True) — note the key is
the *trans* isomer and this paper's is (E), which is the same isomer under a different name, so the
existing key applies.

Every row below shares: in-house alkaline-extracted isoelectric pea protein isolate (protein content
unstated), **0.01 M potassium phosphate, pH 7.6, ~1.25 % v/v methanol, 37 C, 2 h, sealed 20 mL vial,
triplicate**. The headspace rows are at **10 g/L protein**; the fluorescence rows at 0.2 mg/mL.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| Klotz binding constant K, pea protein x hexanal | **684.46 ± 109.43** | M^-1 (per mole of pea protein; **molar mass not stated**) | 37 C, pH 7.6, 10 g/L protein, 0.5-2.5 mM hexanal, **static headspace GC/MS depletion** | Table 1, p. 5 | **binding_constant** |
| Klotz binding constant K, pea protein x (E)-2-octenal | **8207.72 ± 2223.82** | M^-1, as above | as above | Table 1, p. 5 | **binding_constant** — QUARANTINE as reversible (Michael acceptor, no reversibility control; Flags 2) |
| Klotz binding constant K, pea protein x (Z)-2-penten-1-ol | **438.68 ± 259.27** | M^-1, as above | as above | Table 1, p. 5 | **binding_constant** (SE is 59 % of the value) |
| binding sites n, hexanal / (E)-2-octenal / (Z)-2-penten-1-ol | **4.84 ± 0.36 / 26.95 ± 2.92 / 2.55 ± 1.34** | mol per mol protein | as above | Table 1, p. 5 | binding_constant (dimensionless companion) |
| nK, hexanal / (E)-2-octenal / (Z)-2-penten-1-ol | **3313.23 ± 308.45 / 221238.94 ± 39964.95 / 1119.90 ± 273.72** | M^-1 | as above | Table 1, p. 5 | binding_constant |
| **matrix/gas partition coefficient, hexanal, pea protein solution** | **116.37** | dimensionless (matrix/gas) | 37 C, pH 7.6, real-pea level 1260 ug/kg, **PRV, static headspace** | §3.8, p. 8 (Table S3) | **binding_constant** |
| **matrix/gas partition coefficient, hexanal, buffer control** | **32.90** | dimensionless | as above | §3.8, p. 8 | **binding_constant** (the control leg) |
| matrix/gas partition coefficient, (E)-2-octenal, pea / buffer | **2203.85 / 455.93** | dimensionless | 37 C, pH 7.6, real-pea level 22.9 ug/kg, PRV | §3.8, p. 8 | binding_constant |
| matrix/gas partition coefficient, (Z)-2-penten-1-ol, pea / buffer | **2065.44 / 1460.49** | dimensionless | 37 C, pH 7.6, real-pea level 901 ug/kg, PRV | §3.8, p. 8 | binding_constant |
| PRV retention R, (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **29.29 / 71.72 / 79.31** | % retained by the pea protein solution | as above | §3.8, p. 8 | **retention_percent** |
| headspace binding percentage, (E)-2-octenal | **80.52 to 97.35** | % | 37 C, pH 7.6, 10 g/L, 0.05-2.5 mM | §3.1.1, p. 4 | **retention_percent** (span endpoints only; the ten-point curve is Fig. 1A) |
| headspace binding percentage, (Z)-2-penten-1-ol | **0.70 to 11.69** | % | as above | §3.1.1, p. 4 | **retention_percent** (span endpoints only) |
| headspace binding percentage, hexanal | — | % | as above | Fig. 1A | **figure_only** (no hexanal span is printed) |
| ΔH / ΔS / ΔG, hexanal | **−150.916 / −0.466 / −6.456** | kJ mol^-1, kJ mol^-1 K^-1, kJ mol^-1 (at 310 K) | 0.2 mg/mL protein, **fluorescence quenching**, van 't Hoff over 298-310 K | Table 2, p. 6 | binding_constant (thermodynamic companion; **fluorescence-derived, not headspace**) |
| ΔH / ΔS / ΔG, (E)-2-octenal | **105.728 / 0.396 / −16.982** | as above (at 310 K) | as above | Table 2, p. 6 | binding_constant, fluorescence-derived |
| ΔH / ΔS / ΔG, (Z)-2-penten-1-ol | **17.944 / 0.095 / −11.506** | as above (at 310 K) | as above | Table 2, p. 6 | binding_constant, fluorescence-derived |
| effective quenching constant Ka at 310 K, (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **86.9 / 12.2 / 728.5** | M^-1 | 0.2 mg/mL, fluorescence, 37 C | recovered from Table 2 ΔG by K = exp(−ΔG/RT) (mine) | **derived_assumption** (arithmetic only; the paper computed these and did not print them) |
| Stern-Volmer constant Ksv, all three compounds, 25/30/37 C | — | M^-1 | 0.2 mg/mL, fluorescence | Fig. 2A-C (printed inside the panels) | **figure_only** |
| per-gram binding constant from PRV, hexanal | **2.54e-1** | L/g protein | 37 C, pH 7.6, static headspace partition, **assumes 10 g/L** | (K2/K1 − 1)/10 from §3.8 (mine) | **derived_assumption** — the registry's own K_g form; loading assumed (Flags 6) |
| per-gram binding constant from PRV, (E)-2-octenal | **3.83e-1** | L/g protein | as above | (mine) | **derived_assumption**, and quarantined as reversible (Flags 2) |
| per-gram binding constant from PRV, (Z)-2-penten-1-ol | **4.14e-2** | L/g protein | as above | (mine) | **derived_assumption** |
| per-gram binding constant from the Klotz nK | — | L/g protein | needs a pea protein molar mass | §3 item 4 (mine) | **derived_assumption, DO NOT SHIP** — the molar mass is not printed and the answer moves 6x between 60 and 360 kDa |
| alkenal / alkanal contrast, Klotz K | **12.0** | — | 37 C, pH 7.6, 10 g/L, headspace | 8207.72/684.46 (mine) | **within_study_ratio** — FIT candidate for `ALPHA_BETA_UNSATURATION` |
| alkenal / alkanal contrast, PRV matrix/gas K | **18.9** | — | 37 C, pH 7.6, PRV | 2203.85/116.37 (mine) | **within_study_ratio** |
| alkenal / alkanal contrast, PRV per-gram K_g | **1.51** | — | as above | 3.83e-1 / 2.54e-1 (mine) | **within_study_ratio** — note it COLLAPSES from 18.9 to 1.51 once the water baseline is divided out; the registry's term operates on the per-gram scale |
| aldehyde / alcohol contrast, Klotz K | **1.56** | — | 37 C, headspace | 684.46/438.68 (mine) | within_study_ratio |
| method disagreement, headspace Klotz K vs fluorescence Ka at 310 K | **5.0x / 56x / 11x** for (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | — | same lab, same protein, same buffer, same temperature | (mine) | **within_study_ratio** — the strongest single-paper evidence for the `method` field being first-class |
| disagreement between the paper's own two headspace routes, (Z)-2-penten-1-ol | PRV 29.29 % against a §3.1.1 span of 0.70-11.69 % (**2.5x above the top**) | — | 37 C, pH 7.6 | (mine) | within_study_ratio |
| alpha-helix content, no flavour / +1 mM (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **17.9 / 15.9 / 15.4 / 14.4** | % | 0.2 mg/mL, 37 C, CD 200-240 nm, CDPro CONTIN/LL | §3.4, p. 7 | level_only (a conformational response, not a binding term) |
| surface hydrophobicity at 1 mM / 2 mM, (Z)-2-penten-1-ol | **98.38 / 94.59** | % of flavour-free control | ANS, 0.1-0.4 mg/mL, ex 390 / em 470 nm | §3.5, p. 8 | level_only |
| surface hydrophobicity at 1 mM / 2 mM, hexanal | **98.13 / 96.92** | % of control | as above | §3.5, p. 8 | level_only |
| surface hydrophobicity at 1 mM / 2 mM, (E)-2-octenal | **89.89 / 89.26** | % of control | as above | §3.5, p. 8 | level_only |
| concentration in peas, (Z)-2-penten-1-ol / hexanal / (E)-2-octenal | **901 ± 119 / 1260 ± 114 / 22.9 ± 3.54** | ug kg^-1 | raw peas | §3.8, p. 8, from Table S1 via Bi 2020 | level_only (**not measured in this paper**; carried from the group's 2020 paper) |
| docking contact distances | 2.87-3.77 | Å | AutoDock Vina against PDB 3KSC | §3.7, p. 8 | **derived_assumption** (in silico; no energy printed) — not a measurement |

### Can these be put on the same basis as the shipped binding constants?

**(a) The PRV pair can, and it is the row to take.** `REVERSIBLE_BINDING`'s Meynier and
Leksrisompong rows were all built from a within-study partition RATIO divided by a protein loading,
precisely because the absolute static-headspace scale is suspect (6.24x low against literature in
Meynier's case, 6-17x low in Leksrisompong's). Bi's PRV gives the ratio directly and in the same
form — a matrix leg and a water leg measured in the same run on the same instrument — so **the same
offset cancels the same way**. `method` = `static_headspace_partition`, `ph_of_measurement` = 7.6,
`temperature_c` = 37.0, `medium` = a new `pea_isolate` matrix. The one thing needed is the loading
(Flags 6).

**(b) The Klotz K cannot be, without a molar mass.** It is a good number in its own units and it is
useless in L/g until the pea protein molar mass used for Cp is known. Ask, or ship only the PRV
form.

**(c) The temperature is 37 C, which is inside the existing spread but is a body temperature, not a
process one.** The shipped rows sit at 25 C (Damodaran), 30 C (Meynier, Andriot) and 40 C
(Leksrisompong); 37 C needs no transport. But this is an **in-mouth / consumption temperature**
study: it says nothing about what a pea matrix does to hexanal at 90 or 140 C, and the covalent
layer's own note in `matrix_sites.py` is that the sites are charged once at the start of the cook
and their change with heating is not modelled. **Nothing here licenses extrapolation above 37 C.**

**(d) The chain-length slope does not apply.** `CHAIN_LENGTH_SLOPE_PER_CH2` = 2.81x per CH2, from
two ketone series. Bi's ladder crosses functional groups (alcohol -> alkanal -> alkenal) as well as
chain length, so it cannot be used to check or to extend that slope. The C5 alcohol to C6 alkanal
step is 1.56x on K and the C6 alkanal to C8 alkenal step is 12.0x over two carbons (3.46x per CH2 if
one pretends the double bond is not there, which one must not).

**(e) Nothing here goes to `matrix_sites.py`.** That module wants second-order rate constants in
M^-1 s^-1 and activation energies in kJ/mol. This paper has neither. Its ΔH values are van 't Hoff
enthalpies of an equilibrium, not activation energies, and reading −150.916 kJ/mol as an Ea would be
a category error.

## 5. Flags

1. **The Klotz K's molar basis is not printed and the paper is unusable in L/g without it.** Eq. 4
   divides by Cp, "the concentration of pea protein", which must be molar for v to be
   "moles bound per mole of pea protein" and for K to be in M^-1. **No molar mass appears anywhere
   in the paper.** For a mixed isolate of legumin (11S), vicilin (7S) and convicilin there is no
   single right answer. The registry's precedent is Damodaran's stated 100 000 g/mol and Andriot's
   recovered 36 800; here neither route is open. **Request the molar mass used.** Until then the
   Klotz numbers travel only as within-study ratios.
2. **(E)-2-octenal is a Michael acceptor held for 2 h at 37 C against a protein carrying ~0.47 mmol
   lysine amine per gram, and its "binding" is measured by disappearance.** The paper's only
   reversibility evidence is that guanidine hydrochloride and propylene glycol raise the headspace
   peak area (§3.6, Fig. S1C: 1.29e4 -> 1.55e4, a **20 % recovery (mine)**). A 20 % recovery is not
   a demonstration that 80 % was reversible. `matrix_sites.py` charges the
   `unsaturated_aldehyde_amine` class at 5.3-7.9e-5 M^-1 s^-1 at 20 C, and Amendment 6 ruling 4
   sized 22-33 % of Meynier's t-2-hexenal skim-milk "partition" as irreversible chemistry. **This
   row must enter quarantined on the Meynier precedent**: FIT for the unsaturation contrast, refused
   as a reversible binding constant. The hexanal row's propylene-glycol response is in the OPPOSITE
   direction (3.64e7 -> 2.93e7, a **20 % further loss (mine)**), which the authors read as
   hydrogen bonding strengthened by the glycol; it is not a reversibility test either, but hexanal
   carries no Michael hazard.
3. **The three compounds confound four variables at once.** Chain length (C5 / C6 / C8), functional
   group (alcohol / aldehyde / aldehyde), unsaturation (none / none / alpha,beta) and hydrophobicity
   all move together. The paper's own attribution of the ladder to hydrophobicity rests on Table S1,
   which is **not on disk**. **No single-variable conclusion can be drawn from these three
   compounds**, and in particular the 12.0x alkenal/alkanal contrast is also a two-carbon contrast.
   The registry's log-P refusal (k4b guard #4) means the hydrophobicity reading must not be shipped
   in any case.
4. **The two methods disagree by 5x, 56x and 11x within this one paper** (section 3, item 5) and
   they disagree on the ORDER of hexanal against the alcohol. Do not present a Bi number without its
   method. This is the cleanest same-lab, same-day demonstration of the method boundary in the
   corpus and it is worth citing in `parameters_matrix.py`'s §2 header comment beside the k2 sec.
   B.3 35x aldehyde gap.
5. **Table 2's ΔS varies down a block in which ΔH is held constant.** The van 't Hoff treatment of
   Eq. 15 assumes ΔH is temperature-independent; if it is, ΔS = (ΔH − ΔG)/T inherits whatever
   temperature dependence ΔG has, so a varying ΔS is arithmetically consistent but physically means
   the van 't Hoff assumption is failing. For hexanal ΔS moves only −0.468 -> −0.466; for
   (E)-2-octenal 0.383 -> 0.396 (3.4 %); for (Z)-2-penten-1-ol 0.090 -> 0.096 -> 0.095 (**not even
   monotone**). The (Z)-2-penten-1-ol thermodynamics are the least trustworthy in the table, which
   matches its Klotz K carrying a 59 % standard error.
6. **The protein loading in the PRV experiment (§2.10) is not restated.** §2.10 says only "pea
   protein solution containing flavor compound" and "all the conditions were the same as those
   described in Section 2.3.1", where the final protein is 1 % w/v = 10 g/L; the Fig. 1A caption
   independently says 10 mg/mL. **The 10 g/L reading is an inference and every per-gram number in
   section 3 item 2 scales inversely with it.** Confirm before shipping.
7. **Roughly 1.25 % v/v methanol is present in every headspace assay (mine)** — the flavour stocks
   were made up with methanol at a 1:19 methanol:buffer ratio and contribute a quarter of the final
   volume. Methanol competes for hydrophobic sites and changes the air/water partition of all three
   compounds. It is never mentioned again after §2.3.1, is absent from the registry's existing rows,
   and is not corrected for.
8. **The protein is characterised by nothing.** No protein content, no purity, no SDS-PAGE, no
   thiol, no free amine, no solubility. `protein_matrices.yml`'s `pea_isolate` densities come from
   Gao 2020, Xiao 2024, Chen 2022, Shen 2022 and Chihi 2016 on *other* preparations. Pairing Bi's
   binding constant with those site densities is a cross-preparation pairing and should be labelled
   as one.
9. **The retention-vs-concentration curve is non-monotone and the paper's explanation is
   incomplete.** Retention rises over 0.05-0.25 mM and falls over 0.5-2.5 mM, and the high branch
   sits below the low branch. Saturation explains the fall within a branch; it does not explain a
   discontinuity between the two branches, which were run as separate ladders. **The Klotz constants
   in Table 1 are fitted on the HIGH branch only** ("at high flavor compound concentrations"), i.e.
   on the saturating, lower-retention regime. A constant fitted on the low branch would be a
   different number and the paper does not report it.
10. **All supplementary material is off disk.** Table S1 (hydrophobicity constants and the pea
    concentrations), Table S2 (surface hydrophobicity values), Table S3 (the PRV partition
    coefficients as a table), Table S4 (the bond-disrupting-agent key), Fig. S1 (bond-disrupting
    peak areas) and Fig. S2 (CD spectra). Enough of S1 and S3 is quoted verbatim in §3.5 and §3.8
    that the PRV pairs and the pea concentrations are recovered; **S2 and S4 are not recoverable**.
11. **What this paper does NOT contain**: any temperature above 37 C; any pH other than 7.6; any
    ionic-strength or salt series; any heat-treated protein (contrast Guo 2020, which is a preheat
    study); any covalent adduct measurement or mass-spectrometric adduct search; any sensory panel;
    any odour threshold; any ketone (so no point on the `CHAIN_LENGTH_SLOPE_PER_CH2` ladder); any
    commercial isolate; any error bar on the PRV coefficients; any printed Ksv or Ka.
12. **What to request from the authors**: (i) the molar mass used for Cp in Eq. 4; (ii) the protein
    content and purity of the in-house isolate; (iii) the protein loading in the PRV runs; (iv) the
    nine Ksv values and the nine Ka values as numbers rather than figure labels; (v) the ten-point
    binding-percentage table behind Fig. 1A, especially the missing hexanal span; (vi) the
    hydrophobicity constants of Table S1; (vii) a reversibility control on (E)-2-octenal —
    exhaustive dialysis or a borohydride-reduction check — that would settle how much of the 12x
    alkenal advantage is non-covalent at all.
13. **Registry gaps against `data/keys/compounds.yml`**: `hexanal` present; `e_2_octenal` present;
    **`(Z)-2-penten-1-ol` absent** and it would be the first simple aliphatic alcohol carried in
    `COMPOUND_STRUCTURE` (`parameters_matrix.py` §0), which has no `alcohol` binding class — adding
    it means adding a class, not just a key. Separately, `parameters_matrix.py`'s `MATRIX_LOADING`
    table holds exactly five keys — `water`, `skim_milk`, `caseinate_1pct`, `gelatin_3pct`,
    `soy_paste_hong` — and **no pea entry**, so a pea binding row needs a matrix loading created as
    well as a compound key. (Note that `soy_protein` and `beta_lactoglobulin` appear as `medium`
    strings on existing `REVERSIBLE_BINDING` rows without having `MATRIX_LOADING` entries, so a
    `pea_isolate` medium string at 10 g/L would follow that precedent; the loading itself is
    printed by this paper, which is more than either of those two has.)
