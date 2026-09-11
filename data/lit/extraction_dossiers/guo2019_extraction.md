# Guo 2019 — EXTRACTION (soy protein isolate made in-house, preheated at 80, 90 or 100 C, then held at 1 % w/v = 9.21 g/L protein in water at pH 7.0-7.2 and 37 C for 24 h with two esters, two terpene esters, two terpene alcohols and two terpene hydrocarbons; Klotz binding constants from HS-SPME/GC-MS at 0.012-0.08 mM)

### THE PREHEAT ANSWER IS BIDIRECTIONAL, AND IT IS AQUEOUS. `matrix_sites.py` charges its binding sites once at the start of the cook and does not change them with heating. This paper heats soy protein isolate to **80, 90 and 100 C**, then measures binding **in water at 37 C** against an unheated control from the same lot — and the effect **reverses sign between compound classes**: the two simple esters lose **18-21x of their nK** (hexyl acetate 18000 → 850 M^-1; heptyl acetate 44000 → 2400 M^-1, Table 2 p. 22) while the two terpene esters, the two terpene alcohols **gain 1.15-1.71x**, and the two terpene hydrocarbons — which are *salted out* by native SPI — are salted out **less** after preheating (Fig. 4, p. 22). **A single scalar "preheat multiplier" on a site pool is refuted by this one table.** The mechanism the paper proposes is the same one Crowther 1980 proposed on dry soy 39 years earlier and in the opposite direction: high-affinity **primary** sites (buried hydrophobic cavities at the 11S subunit interface) are **destroyed** by heating, while low-affinity **secondary** sites on the newly exposed hydrophobic surface **multiply** — surface hydrophobicity rises **3.6x** from 119.6 to 428.7 at 80 C (Table 1, p. 18). **And the binding constants here CAN be put on the registry's per-gram basis**, because the paper prints both a protein loading (1 % w/v) and a molar mass (220 000 Da) — but that molar mass is **2.2x** the 100 000 Da Damodaran stated for the same protein and on which three shipped soy rows rest (Flags 4).

**Source on disk:** `data/articles/guo2019.pdf` (8 pp., Food Chemistry 290 (2019) 16-23).
Read from the `pdftotext -layout` text layer; **Tables 1 (p. 18) and 2 (p. 22) came through clean and were both
verified against the rendered page images.** **Figures 1, 2 and 3 are the Klotz double-reciprocal plots
themselves** (eight panels each, native / 80 / 90 / 100 C × two compounds) and carry no printed numbers — every
n, K and nK derived from them is in Table 2, so nothing is lost, but **the raw 1/v against 1/[L] points, the
regression lines, and the breakpoint between the primary and secondary branches are figure-only**. **Figure 4
(p. 22) is a two-panel bar chart** of limonene and myrcene release with a-b-c-c significance letters; **its bar
heights are figure-only and are NOT typed as numbers here** — only the fact that every bar exceeds 100 % and the
ordering native > 80 C > 90 C ≈ 100 C. **Supplementary Fig. 1** (the molecular structures of the terpenoids) is
**not on disk**. Repo status before this dossier: Guo 2019 has **no extraction dossier** and is **not cited** in
`src/kinetic_core/parameters_matrix.py`, `src/kinetic_core/matrix_sites.py` or
`data/species/protein_matrices.yml`. **Note the near-name collision on disk**: `guo2020_extraction.md` already
exists and is a *different* paper by the same group (Bi 2022 cites "Guo et al. 2019/2020"); this dossier is for
`guo2019.pdf`, Food Chemistry 290:16-23.

## 0. Identity

| field | value |
|---|---|
| Title | "Binding of aroma compounds with soy protein isolate in aqueous model: Effect of preheat treatment of soy protein isolate" |
| Authors | **Jun Guo** (a,b,c), **Zhiyong He** (a,c), **Shengfang Wu** (a,c), **Maomao Zeng** (a,c, corresponding, mmzeng@jiangnan.edu.cn), **Jie Chen** (a,c, corresponding, chenjie@jiangnan.edu.cn). (a) State Key Laboratory of Food Science and Technology, and School of Food Science and Technology, Jiangnan University, Wuxi 214122, Jiangsu, PR China; (b) Department of Chemistry and Material Engineering, Chizhou University, Chizhou 247100, PR China; (c) International Joint Laboratory on Food Safety, Jiangnan University |
| Venue | Food Chemistry **290** (2019) **16-23** |
| DOI | **`https://doi.org/10.1016/j.foodchem.2019.03.126`** — printed at the foot of p. 16, exactly as given here |
| Dates | Received 22 October 2018; received in revised form 22 March 2019; accepted 24 March 2019; available online 25 March 2019. 0308-8146 / © 2019 Elsevier Ltd |
| Funding | Natural Science Foundation of China (Grant No. 31471583); National First-class discipline program of Food Science and Technology (Grant No. JUFSTR20180201). Conflict of interest: none declared |
| Protein | **Soy protein isolate made in-house** from soybean cv. **Suyun 626** (Fengyuan Seed Co., Lianyungang) by alkaline pH extraction / isoelectric precipitation after Guo et al. 2015. **Protein purity 92.1 % dry weight** by micro-Kjeldahl with a **nitrogen conversion factor of 5.71** (Morr 1985). **Average molecular weight taken as 220 000 Da** — cited from **Beyeler & Solms 1974**, not measured here |
| The 8 ligands | **elongated esters**: hexyl acetate (**HxAc**), heptyl acetate (**HpAc**); **terpene esters**: linalyl formate (**LiFo**), linalyl acetate (**LiAc**); **terpene alcohols**: geraniol, linalool; **terpene hydrocarbons**: limonene, myrcene. All 98 %, J&K Chemical Ltd (Shanghai) |
| Naming | "NSPI"/"nature SPI" = native (unheated); "PSPI" = preheated; "n" = number of binding sites per mole of protein; "K" = Klotz/Scatchard binding constant, M^-1; "nK" = n × K, M^-1; "binding ratio (%)" = the bound fraction; "R_s" = percentage of release (Eq. 6) |
| Companions on disk | `crowther1980_extraction.md` and `aspelund1983_extraction.md` (**the same question — preheat, and functional group — asked on DRY soy 39-40 years earlier, with the opposite sign for the esters**), `damodaran1981_extraction.md` (**the source of three shipped soy rows, and the paper Guo argues against by name**), `bi2022_extraction.md` (the same Klotz/headspace family on pea at the same 37 C), `guo2020_extraction.md` (**a different paper by the same group — do not confuse**), `Xu2022_extraction.md` (the pea-protein preheat analogue) |

## 1. Why it matters

**1. It is the aqueous preheat experiment, and it says the effect is not a scalar.**
`src/kinetic_core/matrix_sites.py` charges free thiol, disulfide and amine pools in mmol per gram once at the
start of the thermal programme and integrates binding against those fixed pools; **how those pools change as
the protein is heated is not modelled**. Crowther 1980, on dry soy, gave one direction and one magnitude
(binding falls 35-49 % after autoclaving). **Guo gives the direction as a function of the ligand**, on soy, in
water, at 37 C, against a native control:

| class | compound | nK native → nK preheated | direction |
|---|---|---|---|
| simple ester | hexyl acetate | 18 000 → 850 M^-1 at 100 C | **falls 21x** |
| simple ester | heptyl acetate | 44 000 → 2 400 M^-1 at 100 C | **falls 18x** |
| terpene ester | linalyl formate | 780 → 1 200 M^-1 at 90 C | **rises 1.54x** |
| terpene ester | linalyl acetate | 2 160 → 3 400 M^-1 at 80 C | **rises 1.57x** |
| terpene alcohol | linalool | 900 → 1 300 M^-1 at 80 C | **rises 1.44x** |
| terpene alcohol | geraniol | 700 → 1 200 M^-1 at 90/100 C | **rises 1.71x** |
| terpene hydrocarbon | limonene, myrcene | release **above 100 %** (salting out), falling towards 100 % with preheat | **an ENHANCEMENT, reduced by heating** |

**Any preheat correction that is a single number applied to a site pool is refuted by this table**, in the same
way and for the same reason that `parameters_matrix.py` refuses "a general matrix correction factor" (k2 sec.
D.1: a uniform 33x misplaces the two extreme compounds by 10x and 28x in *opposite* directions). The two
extreme compounds here move by 21x down and 1.7x up.

**2. Its constants convert to the registry's per-gram form, and by the registry's own construction.** The
shipped soy rows `kg_2_heptanone_soy` (4.40e-3 L/g), `kg_2_octanone_soy` (1.24e-2) and `kg_2_nonanone_soy`
(3.72e-2) were built as **n·K / MW** with Damodaran's molar mass **stated by source**, and the module's
provenance field records exactly that (`"molar_basis": "stated_by_source"`). Guo prints nK directly and prints
a molar mass, so the same arithmetic runs: **K_g = nK / 220 000**, giving 8.18e-2 L/g for hexyl acetate on
native SPI down to 3.86e-3 after 100 C preheat, and 2.00e-1 down to 1.09e-2 for heptyl acetate (section 3).
The paper also prints the protein loading (1 % w/v) and, better, **reveals its own effective protein basis**:
the figure captions read "SPI (1 %, **0.0418 mM**)", and 0.0418 mM × 220 000 g/mol = **9.20 g/L**, which is
1 % w/v **times the 92.1 % purity** — so the molarity was computed on protein, not on solids (mine, and it
closes to 0.2 %). **The effective loading is 9.21 g/L of protein, not 10 g/L.**

**3. It is a direct, aqueous, same-protein check on `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81` — and it passes on the
native protein and fails after preheat.** Hexyl acetate and heptyl acetate differ by exactly one methylene.
On the **native SPI primary sites** the constant ratio is **40 000 / 14 000 = 2.86x per CH2 (mine)** — within
**2 %** of the shipped 2.81, which was itself the geometric mean of Andriot's 2.72 (beta-lactoglobulin,
headspace, 30 C) and Damodaran's 2.9 (soy, dialysis, 25 C). The paper states its own version in words:
*"the binding constant increased by 3-fold for each methylene group increment"* (p. 19), and cites Damodaran &
Kinsella's 2-3-fold prediction. **This is the closest independent confirmation of the shipped slope in the
corpus** — a third protein preparation, a third method, and an ester class rather than a ketone class. But it
holds only on the native primary sites: on the native **secondary** sites the ratio is 1.88x, at 80 C 2.07x,
at 90 C 1.90x and at 100 C 3.00x (mine, and the 100 C K's carry 80-160 % relative errors). **A chain-length
transfer rule licensed on native protein is not licensed on preheated protein by this paper's own data**, and
`CHAIN_LENGTH_SLOPE_PER_CH2` carries no protein-state qualifier today.

**4. It does NOT lift the Wave B26 temperature limit, and the distinction must be stated precisely.** B26's
record says "37 C is an in-mouth temperature, not a process one: nothing here licenses a pea binding constant
at 90 or 140 C." **Guo's binding measurements are also at 37 C.** The 80, 90 and 100 C are the temperatures at
which the *protein was cooked beforehand*; every K, n and nK in Table 2 was then measured in a 37 C water bath
over 24 h. **So this paper supplies a binding constant for a protein that HAS BEEN to 100 C, measured at 37 C.
It does not supply a binding constant AT 100 C.** That distinction matters for the model, because the two
questions have different answers: what a matrix does to a volatile *while it is hot* (unmeasured, everywhere in
this corpus) versus what a *heat-damaged* matrix does at mouth temperature (this paper, and Crowther, and
Xu 2022). **The honest use of Guo is for the second question.** It is nonetheless the most process-relevant
soy binding measurement in the registry's reach, because the protein in a real cooked product *has* been to
90-100 C.

**5. It measures the structural change alongside the binding change, on the same samples.** Table 1 (p. 18)
gives, for the same four protein states, the fluorescence emission maximum (**332.0 → 335.5 → 336.0 → 338.0 nm**,
a red shift = tryptophan moving to a more polar environment), the volume-average particle size D(4,3)
(**598.3 → 422.4 → 355.1 → 345.0 nm**, dissociation not aggregation), and the ANS surface hydrophobicity S0
(**119.6 → 428.7 → 419.3 → 417.2**, a **3.58x jump at 80 C then a plateau**, mine). **The surface hydrophobicity
saturates at 80 C while the ester binding keeps collapsing through 90 and 100 C** — so the fall in ester binding
is *not* tracking exposed hydrophobic surface, which is a direct measured argument against a hydrophobicity- or
log-P-shaped explanation. `parameters_matrix.py` already refuses any log-P term under k4b hold-out guard #4;
**this is corroborating evidence on soy, and it is worth citing there.**

**6. A measured enhancement, which the layer is short of.** Limonene and myrcene are **salted out** by SPI — the
headspace release exceeds 100 % of the protein-free control in every bar of Fig. 4. `REVERSIBLE_BINDING` carries
exactly two negative rows (`kg_delta_decalactone_caseinate` −2.30e-2 and `kg_furaneol_caseinate` −6.20e-2, both
Leksrisompong), and the module note says they are "the only thing in the FIT column that lets this layer emit a
shift below 1 at all". **Guo's terpene hydrocarbons are a third and fourth instance of a protein making an
odourant MORE volatile, on soy rather than on caseinate** — but the magnitudes are figure-only (Flags 8), so
they cannot be shipped as values, only as an existence claim.

**What this paper does NOT give the repository**: any binding measurement above 37 C; any preheating protocol
(**time, protein concentration during heating and cooling regime are all absent — Flags 1**); any aldehyde or
ketone (so no direct overlap with any shipped soy row and **no contact with `matrix_sites.py`'s
aldehyde-amine channel**); any 2-alkenal (so **nothing touches `ALPHA_BETA_UNSATURATION`**); any rate constant
or activation energy; any covalent-adduct measurement; any thermodynamic quantity (no ΔH, ΔS or ΔG anywhere);
any pH other than 7.0-7.2; any thiol, disulfide or free-amine assay of its own isolate; any numeric value for
the limonene and myrcene release.

## 2. Methods as they matter to a model

- **The protein.** Made in-house from soybean **Suyun 626** by alkaline pH extraction and isoelectric
  precipitation, "detailed by Guo et al. (2015)" — the protocol itself is by reference and is not restated. After
  neutralisation to **pH 7.0**, protein content by **micro-Kjeldahl** with **N × 5.71** and total solids by oven
  drying (105 C overnight): **protein purity 92.1 % on a dry-weight basis**. An **8 % (w/v)** SPI suspension in
  deionised water was **centrifuged 10 000 × g for 20 min** to remove particulates. **Ionic strength, expressed
  as NaCl, 0.03-0.04 M**, by conductivity meter (Mettler Toledo S30 SevenEasy). **Average molecular weight taken
  as 220 000 Da**, cited to Beyeler & Solms 1974 and apportioned "based on molecular weight and percentage of
  each protein fraction" (Fukushima 1991; O'Keefe 1991). **No thiol, disulfide or amine assay.**
- **THE PREHEAT TREATMENT IS NEVER DESCRIBED.** Section 2.2 makes the isolate; section 2.3 is headed
  *"Evaluation of characteristics of preheated SPI"* and goes straight into fluorescence, surface hydrophobicity
  and particle size. **Nowhere in the paper is there a heating time, a heating vessel, a protein concentration
  during heating, a heating rate, a hold duration or a cooling protocol** — only the three temperatures 80, 90
  and 100 C, which appear first as column headers in Table 1. This was checked against the rendered page image
  of p. 17 and the omission is real, not a text-layer artefact (Flags 1).
- **The pot.** **SPI at 1 % (w/v) in water, pH 7.0-7.2**, sealed in **18 mL flasks**, **10 mL sample volume**
  with SPI and **10 mL control volume** without. **Effective protein loading 9.21 g/L** (1 % w/v × 92.1 %
  purity), confirmed by the figure captions' "0.0418 mM" against the stated 220 000 Da (mine, section 3).
  **Note: water, not a buffer** — pH 7.0-7.2 is the pH of the neutralised protein suspension itself, unbuffered.
- **The ligands and their dosing.** A **1 mM stock of each flavour in PROPYLENE GLYCOL**, made by gradient
  dilution. Final concentrations **0.012 to 0.08 mM** for HxAc, HpAc, LiFo, LiAc, geraniol and linalool;
  **limonene and myrcene fixed at 0.024 mM**. The choice of propylene glycol is justified explicitly:
  *"Propylene glycol is a suitable solvent for flavor compounds because of its low vapor pressure, avoiding
  competition between flavors and the solvent in the fiber coating"* — i.e. it was chosen to protect the SPME
  fibre, not to be inert towards the protein (Flags 6). **Table 2's own C(flavor) column prints 0.12-0.24,
  0.24-0.8 and 0.12-0.8 mM — exactly 10x the Methods and the figure captions (Flags 2).**
- **Equilibration.** *"Sample and control were made and shaken for 24 h at 37 C to reach equilibration."* Then,
  during extraction, stirred at **250 rpm** in a water bath at **37 C for 5 min**.
- **Measurement family: HEADSPACE SOLID-PHASE MICROEXTRACTION, HS-SPME/GC-MS.** Fibre **50/30 µm
  DVB/CAR/PDMS** (AnPu, Shanghai), conditioned per the manufacturer; fibre type, extraction time and agitation
  all optimised beforehand; **5 min extraction at 37 C**, then **5 min desorption** in a 250 C splitless
  injector. GC/MS: Bruker SCION SQ456, **Supelcowax 10 fused silica 30 m × 0.25 mm (Agilent DB-WAX)**, helium at
  **0.8 mL/min**, injector and detector 250 C, oven **40 C for 3 min → 90 C at 5 C/min → 250 C at 10 C/min**.
  EI 70 eV, detector 350 V, **m/z 33-450**, 3.00 scans/s.
  **Where this sits in the registry's `method` taxonomy.** The constant is fitted to **headspace peak areas
  with and without protein** by a Klotz plot — structurally identical to Bi 2022's construction, which the
  registry classes `headspace_depletion` (matching Andriot 2000's beta-lactoglobulin rows). **So `method` =
  `headspace_depletion`, with a provenance note that the sampling is SPME rather than static headspace.** It is
  emphatically **not** on the dialysis side of the k2 sec. B.3 boundary, and must not be pooled with
  `kg_nonanal_soy` (equilibrium dialysis with 2-mercaptoethanol) — though the compound classes do not overlap in
  any case. The paper's own reason for choosing SPME is sensitivity: it cites Roberts 2000 for **"1800-fold
  more sensitivity for non-polar flavors than top space concentration"**, which is what allows the 0.012-0.08 mM
  *unsaturated* dosing range that is the paper's methodological claim against Damodaran's saturation-level
  dialysis.
- **The arithmetic, Eqs. 1-6 (pp. 17-18).** Scatchard v/[L] = nK − vK (Eq. 1); the isotherm
  v = n[L]/(1 + K[L]) (Eq. 2); the Klotz double reciprocal **1/v = 1/n + 1/(nK[L])** (Eq. 3), so **1/n is the
  intercept and 1/(nK) the slope**. The free concentration and the bound moles come from peak areas:
  **[L] = ([HS]_P / [HS]_C) × O** (Eq. 4) and **v = {([HS]_C − [HS]_P)/[HS]_C × O} / C_P** (Eq. 5), where O is
  the flavour concentration of the control, C_P the protein concentration of the sample, [HS]_C and [HS]_P the
  control and sample peak areas. **R_s = ([HS]_P/[HS]_C) × 100** (Eq. 6) is the percentage of release — the
  quantity plotted in Fig. 4, and **values above 100 mean the protein DRIVES the compound into the headspace**.
  Standard errors of K, n and nK were propagated from the SEs of the intercept and slope using SPSS 19.0, with
  the variance expressions of Kühn 2007.
- **The structural assays.** **Fluorescence**: Hitachi F-2700, excitation **280 nm**, emission 300-450 nm, 5 nm
  slits, SPI at **0.02 % w/v in 10 mM phosphate buffer pH 7.2** — note this is a *different* medium and a
  **500x lower protein concentration** than the binding assay. **Surface hydrophobicity**: ANS probe after
  Wagner 2000, H0 as the initial slope of fluorescence index against protein concentration. **Particle size**:
  Microtrac S3500, distilled water dispersant, relative refractive index **1.095**, volume-average diameter
  **d(4,3)** recorded.
- **Replication and statistics.** Table 1 values are "means of three determinations". ANOVA was performed on
  fluorescence, surface hydrophobicity, particle size and R_s; letters denote p < 0.05. **Table 2 carries
  standard errors but the paper never states the number of replicates behind them** — they are propagated from
  a single Klotz regression per condition (Flags 5).

## 3. Tables re-typed

Both tables were verified against the rendered page images. Evidence marks: `[M]` measured in this study,
`[C]` cited from elsewhere, `[F]` fitted.

### Table 1 (p. 18). "Structure and surface property changes of NSPI and heat-denatured SPI at different temperatures."

| parameters | native SPI | preheated SPI, 80 (°C) | 90 (°C) | 100 (°C) |
|---|---|---|---|---|
| fluorescence wave_max (nm) | 332^c ± 0.6 `[M]` | 335.5^b ± 1.5 `[M]` | 336.0^b ± 1.7 `[M]` | 338.0^a ± 1.5 `[M]` |
| particle size, D(4,3) (nm) | 598.3^a ± 12.6 `[M]` | 422.4^b ± 16.9 `[M]` | 355.1^c ± 21.4 `[M]` | 345.0^c ± 18.5 `[M]` |
| surface hydrophobicity (S0) | 119.6^c ± 2.6 `[M]` | 428.7^a ± 6.9 `[M]` | 419.3^b ± 11.4 `[M]` | 417.2^b ± 8.5 `[M]` |

Footnote exactly as printed: *Values are means of three determinations; different lower case letters (a-c) in
the same row indicate significant difference among the values at the 95 % confidence level (p < 0.05).*

**S0 is dimensionless as printed** (no unit is given for surface hydrophobicity anywhere in the paper).
**Read the letters: for surface hydrophobicity the 80 C sample is `a` and the 90 and 100 C samples are `b`, so
S0 PEAKS at 80 C and then falls slightly** — it does not keep rising, contrary to the running text on p. 18
("the surface hydrophobicity of PSPI continued to increase"), which the table's own letters contradict
(Flags 3).

### Table 2 (p. 22). "Binding parameters of SPI with HxAc, HpAc, LiFo, LiAc, linalool and geraniol."

Header exactly as printed, including the typographical slip "binding rato (%)".

| flavor | T (°C) | C(flavor) (mM) | binding rato (%) | n | K (M^-1) | nk (M^-1) |
|---|---|---|---|---|---|---|
| HxAc | native | 0.12-0.24 | 46-75 `[M]` | 1.3 ± 0.5 `[F]` | 14000 ± 6000 `[F]` | 18000^a ± 4000 `[F]` |
| HxAc | native | 0.24-0.8 | | 6.0 ± 3.8 `[F]` | 800 ± 300 `[F]` | 4800 ± 240 `[F]` |
| HxAc | 80 | 0.12-0.24 | 33-60 `[M]` | 0.8 ± 0.5 `[F]` | 14000 ± 4800 `[F]` | 11000^b ± 3000 `[F]` |
| HxAc | 80 | 0.24-0.8 | | 11.0 ± 5.4 `[F]` | 290 ± 240 `[F]` | 3200 ± 300 `[F]` |
| HxAc | 90 | 0.12-0.8 | 17-37 `[M]` | 5.0 ± 3.4 `[F]` | 200 ± 120 `[F]` | 1000^c ± 180 `[F]` |
| HxAc | 100 | 0.12-0.8 | 20-42 `[M]` | 8.5 ± 5.1 `[F]` | 100 ± 80 `[F]` | 850^d ± 160 `[F]` |
| | | | | | | |
| HpAc | native | 0.12-0.24 | 39-88 `[M]` | 1.4 ± 0.8 `[F]` | 40000 ± 9000 `[F]` | 44000^a ± 10000 `[F]` |
| HpAc | native | 0.24-0.8 | | 6.5 ± 2.5 `[F]` | 1500 ± 600 `[F]` | 10000 ± 2700 `[F]` |
| HpAc | 80 | 0.12-0.24 | 25-66 `[M]` | 0.8 ± 0.3 `[F]` | 29000 ± 9000 `[F]` | 24000^b ± 4100 `[F]` |
| HpAc | 80 | 0.24-0.8 | | 10.0 ± 4.0 `[F]` | 510 ± 220 `[F]` | 5100 ± 350 `[F]` |
| HpAc | 90 | 0.12-0.8 | 16-38 `[M]` | 5.2 ± 3.8 `[F]` | 380 ± 290 `[F]` | 1900^c ± 280 `[F]` |
| HpAc | 100 | 0.12-0.8 | 24-44 `[M]` | 8.0 ± 4.1 `[F]` | 300 ± 160 `[F]` | 2400^c ± 360 `[F]` |
| | | | | | | |
| LiFo | native | 0.12-0.8 | 13-29 `[M]` | 3.9 ± 2.2 `[F]` | 200 ± 80 `[F]` | 780^b ± 210 `[F]` |
| LiFo | 80 | 0.12-0.8 | 20-34 `[M]` | 7.5 ± 4.2 `[F]` | 120 ± 100 `[F]` | 900^a ± 280 `[F]` |
| LiFo | 90 | 0.12-0.8 | 17-28 `[M]` | 4.0 ± 2.3 `[F]` | 300 ± 260 `[F]` | 1200^a ± 220 `[F]` |
| LiFo | 100 | 0.12-0.8 | 18-30 `[M]` | 3.5 ± 2.0 `[F]` | 300 ± 180 `[F]` | 1050^a ± 200 `[F]` |
| | | | | | | |
| LiAc | native | 0.12-0.8 | 17-29 `[M]` | 3.6 ± 1.8 `[F]` | 600 ± 400 `[F]` | 2160^c ± 200 `[F]` |
| LiAc | 80 | 0.12-0.8 | 23-32 `[M]` | 6.8 ± 4.0 `[F]` | 530 ± 320 `[F]` | 3400^a ± 570 `[F]` |
| LiAc | 90 | 0.12-0.8 | 18-28 `[M]` | 4.0 ± 2.2 `[F]` | 620 ± 400 `[F]` | 2500^b ± 480 `[F]` |
| LiAc | 100 | 0.12-0.8 | 20-32 `[M]` | 4.3 ± 2.7 `[F]` | 570 ± 320 `[F]` | 2300^b ± 410 `[F]` |
| | | | | | | |
| Linalool | native | 0.12-0.8 | 7-19 `[M]` | 1.0 ± 0.8 `[F]` | 900 ± 390 `[F]` | 900^b ± 210 `[F]` |
| Linalool | 80 | 0.12-0.8 | 13-26 `[M]` | 3.6 ± 1.8 `[F]` | 380 ± 190 `[F]` | 1300^a ± 230 `[F]` |
| Linalool | 90 | 0.12-0.8 | 11-22 `[M]` | 2.2 ± 1.6 `[F]` | 500 ± 300 `[F]` | 1100^a ± 260 `[F]` |
| Linalool | 100 | 0.12-0.8 | 9-19 `[M]` | 1.9 ± 1.1 `[F]` | 570 ± 340 `[F]` | 1000^b ± 170 `[F]` |
| | | | | | | |
| Geraniol | native | 0.12-0.8 | 9-18 `[M]` | 1.4 ± 0.8 `[F]` | 500 ± 260 `[F]` | 700^b ± 110 `[F]` |
| Geraniol | 80 | 0.12-0.8 | 12-24 `[M]` | 3.6 ± 2.2 `[F]` | 280 ± 140 `[F]` | 1000^a ± 220 `[F]` |
| Geraniol | 90 | 0.12-0.8 | 14-26 `[M]` | 1.5 ± 0.9 `[F]` | 810 ± 460 `[F]` | 1200^a ± 290 `[F]` |
| Geraniol | 100 | 0.12-0.8 | 13-25 `[M]` | 1.2 ± 0.8 `[F]` | 1000 ± 450 `[F]` | 1200^a ± 300 `[F]` |

Footnote exactly as printed: *Different lower case letters (a-c) in the same row indicate significant
differences among values at the 95 % confidence level (p < 0.05).* **The footnote is wrong: the letters vary
down the nK COLUMN within each compound block, not across a row** (a row has only one lettered value). And
HxAc uses four letters a-d while the footnote says a-c (Flags 3).

**Read the structure of this table before using it.** For **HxAc and HpAc on native and 80 C protein there are
TWO ROWS PER CONDITION** — a high-affinity **primary** class fitted on the low-concentration branch and a
low-affinity **secondary** class fitted on the high-concentration branch. At **90 and 100 C the primary class
is GONE** and a single fit spans the whole range. For **LiFo, LiAc, linalool and geraniol there is only ever
one class**, at every temperature: *"Primary binding sites for LiFo, LiAc, linalool, and geraniol were not
observed"* (p. 20). **So the native-versus-preheated comparison for the two esters is a comparison between a
two-class and a one-class fit, and there is no single number on either side of it** (Flags 7).

**Internal check (mine): does n × K equal the printed nK?** It should, since nK is the product.

| row | n × K (mine) | printed nK | agreement |
|---|---:|---:|---|
| HxAc native primary | 18 200 | 18 000 | 1 % ✓ |
| HxAc native secondary | 4 800 | 4 800 | exact ✓ |
| HxAc 80 primary | 11 200 | 11 000 | 2 % ✓ |
| HxAc 80 secondary | 3 190 | 3 200 | 0.3 % ✓ |
| HxAc 90 | 1 000 | 1 000 | exact ✓ |
| HxAc 100 | 850 | 850 | exact ✓ |
| **HpAc native primary** | **56 000** | **44 000** | **27 % — DOES NOT CLOSE** |
| HpAc native secondary | 9 750 | 10 000 | 2.5 % ✓ |
| HpAc 80 primary | 23 200 | 24 000 | 3 % ✓ |
| HpAc 80 secondary | 5 100 | 5 100 | exact ✓ |
| HpAc 90 | 1 976 | 1 900 | 4 % ✓ |
| HpAc 100 | 2 400 | 2 400 | exact ✓ |
| LiFo, all four | 780 / 900 / 1 200 / 1 050 | 780 / 900 / 1 200 / 1 050 | exact ✓ |
| LiAc, all four | 2 160 / 3 604 / 2 480 / 2 451 | 2 160 / 3 400 / 2 500 / 2 300 | ≤ 6 % ✓ |
| Linalool, all four | 900 / 1 368 / 1 100 / 1 083 | 900 / 1 300 / 1 100 / 1 000 | ≤ 8 % ✓ |
| Geraniol, all four | 700 / 1 008 / 1 215 / 1 200 | 700 / 1 000 / 1 200 / 1 200 | ≤ 1 % ✓ |

**Twenty-seven of twenty-eight rows close to 8 % or better. One does not: the HpAc native primary row, where
1.4 × 40 000 = 56 000 against a printed 44 000.** That is the **single highest-affinity number in the paper**
and the headline of the whole preheat comparison for heptyl acetate. One of n, K or nK in that row is
mis-transcribed and there is no way to tell which from the printed data (Flags 5).

### Figure 4 (p. 22). "Release of limonene (A), myrcene (B) in SPI (1 %, 0.0418 mM) solution, at 37 °C (0.024 mM)."

**A two-panel bar chart with a y-axis "Release of flavor (%)" running 0-160 and four bars per panel — Native,
80 °C, 90 °C, 100 °C — carrying significance letters a, b, c, c in both panels.** The bar heights are
**figure-only and are not typed here** per house rule. What is unambiguous and load-bearing:
- **Every bar in both panels lies ABOVE 100 %.** Since R_s = ([HS]_P/[HS]_C) × 100 (Eq. 6), a value above 100
  means the protein solution puts **more** of the compound into the headspace than the protein-free control.
  This is the **salting-out effect** the paper names in its abstract and discusses at length.
- **The ordering is Native > 80 °C > 90 °C ≈ 100 °C** in both panels, with the letters `a`, `b`, `c`, `c`
  confirming that 90 and 100 C are statistically indistinguishable and both differ from native and from 80 C.
- **Preheating REDUCES the salting-out**, i.e. moves the release back towards 100 %. The conclusion states it:
  *"The salting out effect of myrcene and limonene was also observed, and PSPI reduced the salting out effect."*
- **No numeric release value for limonene or myrcene appears anywhere in the text or in a table.**

### Numbers printed in the running text

| quantity | value | where | class |
|---|---|---|---|
| SPI protein purity | **92.1 %** dry weight, micro-Kjeldahl, **N × 5.71** | p. 17 §2.2 | `[M]` |
| SPI molecular weight used for molarity | **220 000 Dalton** | p. 17 §2.2 | **`[C]` — Beyeler & Solms 1974, NOT measured here** |
| ionic strength as NaCl | **0.03-0.04 M** | p. 17 §2.2 | `[M]` |
| protein loading in the binding assay | **1 % (w/v)** in water, **pH 7.0-7.2**, 10 mL in an 18 mL flask | p. 17 §2.4.1 | `[M]` |
| the same, as a molarity | **0.0418 mM** | Fig. 1, 2, 3 and 4 captions | `[M]` — **implies 9.21 g/L of protein (mine)** |
| flavour concentration range | **0.012 to 0.08 mM**; limonene and myrcene fixed at **0.024 mM** | p. 17 §2.4.1; Fig. 1-3 captions | `[M]` — **and Table 2 prints 10x these (Flags 2)** |
| equilibration | **24 h at 37 °C**, shaken | p. 17 §2.4.1 | `[M]` |
| SPME extraction | **5 min at 37 °C**, 250 rpm stirring; 5 min desorption at 250 °C | p. 17 §2.4.2 | `[M]` |
| SPME sensitivity advantage | **1800-fold** more sensitive for non-polar flavours than headspace concentration | p. 16 | `[C]` (Roberts 2000) |
| fluorescence assay protein | **0.02 % w/v in 10 mM phosphate buffer pH 7.2** (a 50x lower loading and a different medium) | p. 17 §2.3.1 | `[M]` |
| particle-size refractive index | 1.095 | p. 17 §2.3.3 | `[C]` (Tang 2011) |
| **the chain-length claim** | *"the binding constant increased by 3-fold for each methylene group increment in the chain"* | p. 19 | **`[M]`/interpretation — 2.86x on the printed numbers (mine)** |
| the prediction it is compared with | binding affinity increases **2-3 fold** per additional methylene for hydrophobic interactions | p. 19 | `[C]` (Damodaran & Kinsella 1981b) |
| primary binding sites, HxAc and HpAc, native and 80 C | **1-2**, "very high affinity ... at low concentration" | p. 19 | `[F]` |
| surface (secondary) binding sites, same | **5-8**, "low affinity at high concentration" | p. 19 | `[F]` |
| denaturation temperature of the 7S subunit | **below 80 °C** | pp. 18, 20 | `[C]` (Sorgentini 1995) |
| denaturation temperature of the 11S subunit | **approximately 95 °C** | pp. 18, 20 | `[C]` (Sorgentini 1995) |
| the 90 C result | *"the primary binding site of SPI disappeared, and secondary binding sites increased ... binding constants and the value of nK decreased"* | p. 19 | `[M]` |
| the mechanism proposed | 11S dissociating into smaller subunits destroys the **hydrophobic cavities**; the dissociated subunits have more hydrophobic surface, so **secondary sites increase** | p. 19 | interpretation |
| why the terpenoids behave differently | **steric hindrance** limits access of LiFo and LiAc to the hydrophobic cavities, so they bind on the **hydrophobic surface** instead | p. 20 | interpretation |
| LiAc > LiFo explained | the acetate group has "an additional methylene group" over the formate | p. 20 | interpretation |
| **the contrary literature the paper argues against** | Damodaran & Kinsella 1981b: binding of **2-nonanone INCREASED** on preheated SPI (90 °C, 1 h) by equilibrium dialysis, with the number of sites unchanged and the constant changed | pp. 16, 20 | **`[C]` — the opposite sign to this paper, on the same protein** |
| the conditions Guo blames for the difference | Damodaran's **pH 8.0, 30 mM Tris-HCl, 0.02 % sodium azide and 10 mM 2-mercaptoethanol**, and saturation-level dosing | p. 20 | `[C]` |
| a third literature answer | Semenova 2002, HxAc on heat-treated **broad-bean 11S** by ultrafiltration in 50 mM phosphate pH 7.2: heating **increased the number of sites but decreased the binding constant** | p. 20 | `[C]` |
| a fourth | Heng 2004: binding capacity **decreased** with heat treatment of **vicilin** at 90 °C for 30 min at 0.006-0.04 mM | p. 21 | `[C]` |
| a fifth | McNeill & Schmidt 1993: vanillin binding capacity **increased** with heat treatment of milk proteins at 85 °C above 70 ppm | p. 21 | `[C]` |
| the authors' own verdict on the literature | *"It is difficult to compare results and reach conclusions based on thermally denatured SPI because of differences in approaches and experimental conditions used in different studies."* | p. 21 | — |

**Figure-only quantities.** Every Klotz plot point, regression line and breakpoint in Figs. 1-3 (24 panels
total); every bar height in Fig. 4; the molecular structures in Supplementary Fig. 1 (not on disk). Per house
rule none is typed as a number here.

### Arithmetic on the printed constants (all mine)

**1. The effective protein loading, recovered from the figure captions (mine).** The captions read
"SPI (1 %, 0.0418 mM)". At the stated 220 000 g/mol, **0.0418 mM × 220 000 = 9.196 g/L**. Against a 1 % w/v
solids loading of 10 g/L and a stated purity of 92.1 %, **10 × 0.921 = 9.21 g/L → 9.21/220 000 = 0.0419 mM**,
which reproduces the printed 0.0418 to **0.2 %**. **The molarity was computed on protein, not on total solids,
so the protein loading is 9.21 g/L.** (For contrast, using 10 g/L would give 0.0455 mM, 8.7 % off the printed
value.) This matters because every per-gram conversion divides by it.

**2. The constants on the registry's per-gram basis, K_g = nK / MW (mine).** This is the construction
`kg_2_heptanone_soy` and its two siblings were built with, using a molar mass **stated by the source**.
With MW = 220 000 g/mol:

| compound | protein state | nK, M^-1 | **K_g = nK/220 000, L/g (mine)** |
|---|---|---:|---:|
| hexyl acetate | native, primary | 18 000 | **8.18e-2** |
| hexyl acetate | native, secondary | 4 800 | 2.18e-2 |
| hexyl acetate | 80 C, primary | 11 000 | 5.00e-2 |
| hexyl acetate | 80 C, secondary | 3 200 | 1.45e-2 |
| hexyl acetate | 90 C | 1 000 | 4.55e-3 |
| hexyl acetate | 100 C | 850 | **3.86e-3** |
| heptyl acetate | native, primary | 44 000 | **2.00e-1** |
| heptyl acetate | native, secondary | 10 000 | 4.55e-2 |
| heptyl acetate | 80 C, primary | 24 000 | 1.09e-1 |
| heptyl acetate | 80 C, secondary | 5 100 | 2.32e-2 |
| heptyl acetate | 90 C | 1 900 | 8.64e-3 |
| heptyl acetate | 100 C | 2 400 | **1.09e-2** |
| linalyl formate | native / 80 / 90 / 100 C | 780 / 900 / 1 200 / 1 050 | 3.55e-3 / 4.09e-3 / 5.45e-3 / 4.77e-3 |
| linalyl acetate | native / 80 / 90 / 100 C | 2 160 / 3 400 / 2 500 / 2 300 | 9.82e-3 / 1.55e-2 / 1.14e-2 / 1.05e-2 |
| linalool | native / 80 / 90 / 100 C | 900 / 1 300 / 1 100 / 1 000 | 4.09e-3 / 5.91e-3 / 5.00e-3 / 4.55e-3 |
| geraniol | native / 80 / 90 / 100 C | 700 / 1 000 / 1 200 / 1 200 | 3.18e-3 / 4.55e-3 / 5.45e-3 / 5.45e-3 |

**For scale**: the shipped soy rows are 4.40e-3 (2-heptanone), 1.24e-2 (2-octanone), 3.72e-2 (2-nonanone) and
4.38e-2 L/g (nonanal), all by dialysis at 25 C. **Native SPI's heptyl acetate primary constant, 2.00e-1 L/g, is
4.6x the largest shipped soy value; after a 100 C preheat it is 1.09e-2 L/g, below the 2-octanone row.** The
whole span of this one paper, 3.18e-3 to 2.00e-1 L/g, is **63x**, and it is spanned by protein state and ligand
class alone at one temperature and one pH. **These numbers are all conditional on the 220 000 Da basis
(Flags 4) and none should be shipped without resolving that.**

**3. The preheat factor, as a ratio on nK (mine).** Preheated nK divided by native nK, per compound:

| compound | 80 C / native | 90 C / native | 100 C / native |
|---|---:|---:|---:|
| hexyl acetate (primary class; 90 and 100 are single-class fits) | **0.61** | **0.056** | **0.047** |
| hexyl acetate (secondary class, native vs 80 C only) | **0.67** | — | — |
| heptyl acetate (primary class) | **0.55** | **0.043** | **0.055** |
| heptyl acetate (secondary class, native vs 80 C only) | **0.51** | — | — |
| linalyl formate | **1.15** | **1.54** | **1.35** |
| linalyl acetate | **1.57** | **1.16** | **1.07** |
| linalool | **1.44** | **1.22** | **1.11** |
| geraniol | **1.43** | **1.71** | **1.71** |

**The spread across ligands at a single preheat temperature is 0.047 to 1.71 — a factor of 36 (mine).** Note
also that the terpenoids are **non-monotone in preheat temperature**: LiAc and linalool peak at 80 C and fall
back, LiFo peaks at 90 C, geraniol plateaus at 90-100 C. **There is no monotone preheat function here for any
class**, and for the esters the 90 → 100 C step even reverses for HpAc (1900 → 2400).

**Comparing like with like on the secondary sites only**, which is the one class that exists at every
temperature for the esters: HxAc 4 800 (native) → 3 200 (80 C) → 1 000 (90 C) → 850 (100 C), a **5.6x fall**;
HpAc 10 000 → 5 100 → 1 900 → 2 400, a **4.2x fall**. **The 18-21x figure quoted in the headline is
primary-against-single-class and overstates the like-for-like collapse by roughly 4x** (Flags 7). The
like-for-like number, **4-6x**, is the defensible one.

**4. The chain-length slope, HpAc over HxAc (one CH2), per protein state (mine).**

| protein state and class | on K | on nK |
|---|---:|---:|
| native, primary | **2.86x** | 2.44x |
| native, secondary | 1.88x | 2.08x |
| 80 C, primary | 2.07x | 2.18x |
| 80 C, secondary | 1.76x | 1.59x |
| 90 C, single class | 1.90x | 1.90x |
| 100 C, single class | 3.00x | 2.82x |

**The native primary value, 2.86x, lands within 2 % of the shipped `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81`.** That
is the strongest independent confirmation of that constant in the corpus: a different protein preparation, a
different method (SPME headspace vs dialysis and static headspace), a different compound class (esters vs
2-alkanones) and a different temperature (37 C vs 25-30 C). **But it is ONE methylene step on ONE pair of
compounds**, whereas Andriot's 2.72 and Damodaran's 2.9 each came from a three-member series, and the other five
determinations in this same table range 1.59-3.00. **The arithmetic mean of all six K-based values is 2.24x and their
geometric mean 2.19x (mine)** — which is where Aspelund's 2.23-2.27x and Crowther's 1.90-2.15x also sit.
**Read together, the corpus's chain-length evidence clusters near 2.2-2.3x with the 2.81 sitting at the top of
the range, and the shipped value is defensible only for the high-affinity, native-protein case.**

**5. Surface hydrophobicity does not track ester binding (mine).** S0 rises **3.58x** from native (119.6) to
80 C (428.7), then falls slightly to 419.3 (90 C) and 417.2 (100 C) — a **2.7 % decline** across the 80-100 C
range. Over exactly that range, HxAc's nK falls from 11 000 to 850, a **12.9x collapse**. **A 2.7 % change in
exposed hydrophobic surface accompanies a 13x change in ester binding.** The two are not proportional, not
monotone together, and not plausibly the same mechanism. This is a clean measured argument against a
hydrophobicity-driven account of the ester result, and it supports `parameters_matrix.py`'s standing refusal of
any log-P-shaped matrix term.

**6. Particle size and binding (mine).** D(4,3) falls **1.42x** native → 80 C, **1.19x** 80 → 90 C, and
**1.03x** 90 → 100 C — i.e. the dissociation is **essentially complete by 90 C**. The ester binding, by
contrast, does most of its collapsing **between 80 and 90 C** (HxAc 11 000 → 1 000, an 11x fall in the step
where particle size moves only 19 %). **The binding collapse coincides with the 11S denaturation temperature
the paper quotes (≈ 95 C), not with the particle-size change**, which is consistent with the paper's own
cavity-destruction mechanism and inconsistent with a simple aggregation-state explanation.

**7. Binding ratio as an independent read of the same effect (mine).** The printed binding-ratio ranges give
HxAc 46-75 % (native) → 33-60 % (80 C) → 17-37 % (90 C) → 20-42 % (100 C), and HpAc 39-88 % → 25-66 % →
16-38 % → 24-44 %. **On the upper end of the range the fall is 75 → 42 % for HxAc (1.8x) and 88 → 44 % for HpAc
(2.0x)** — an order of magnitude gentler than the 18-21x fall in nK, and roughly consistent with the 4-6x
like-for-like nK fall once one accounts for a bound fraction being a saturating function of nK. **A model
consumer that cares about how much odourant is actually held should read the binding ratio, not nK.**

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **None of this paper's eight compounds is in the registry** —
`hexyl_acetate`, `heptyl_acetate`, `linalyl_formate`, `linalyl_acetate`, `linalool`, `geraniol`, `limonene` and
`myrcene` all return nothing on a key search. Nor is any of them in `COMPOUND_STRUCTURE` in
`parameters_matrix.py`, which carries no ester, no terpene and no terpene hydrocarbon class at all. (Note the
same is true of Meynier's `isoamyl_acetate`, `amyl_acetate` and `ethyl_pentanoate`, which appear as `compound`
strings on `REVERSIBLE_BINDING` rows without having entries in `compounds.yml` — so a Guo row would follow that
precedent rather than break new ground.) **No compound here is a panel target**, so nothing in this paper
changes a prediction directly; its value to the repository is entirely in the **ratios** and in the
**preheat-mechanism finding**.

Every row below shares: **in-house alkaline-extracted isoelectric soy protein isolate (92.1 % protein), preheated
at the stated temperature by an UNDESCRIBED protocol, then 1 % w/v (= 9.21 g/L protein) in unbuffered water at
pH 7.0-7.2, 10 mL in an 18 mL sealed flask, flavour dosed from a 1 mM propylene-glycol stock to 0.012-0.08 mM,
shaken 24 h at 37 C, HS-SPME 5 min at 37 C on a 50/30 DVB/CAR/PDMS fibre, GC-MS, Klotz double-reciprocal fit.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **Klotz K, hexyl acetate × NATIVE SPI, primary sites** | **14 000 ± 6 000** | M^-1 (per mole of protein, MW 220 000 stated) | 37 C, pH 7.0-7.2, 9.21 g/L, 0.012-0.024 mM | Table 2, p. 22 | **`binding_constant`** |
| Klotz K, hexyl acetate × native SPI, secondary sites | **800 ± 300** | M^-1 | as above, 0.024-0.08 mM | Table 2, p. 22 | `binding_constant` |
| **Klotz K, hexyl acetate × 100 C PREHEATED SPI** | **100 ± 80** | M^-1 | as above, single class over the whole range | Table 2, p. 22 | `binding_constant` (SE is 80 % of the value) |
| **Klotz K, heptyl acetate × NATIVE SPI, primary sites** | **40 000 ± 9 000** | M^-1 | 37 C, pH 7.0-7.2, 9.21 g/L | Table 2, p. 22 | **`binding_constant`** — the largest constant in the paper, **and its row fails the n×K check (Flags 5)** |
| Klotz K, heptyl acetate × 100 C preheated SPI | **300 ± 160** | M^-1 | as above | Table 2, p. 22 | `binding_constant` |
| Klotz K and nK, all 6 compounds × 4 protein states | see Table 2 | M^-1 | as above | Table 2, p. 22 | `binding_constant` |
| binding sites n, all rows | 0.8 to 11.0 | mol per mol protein | as above | Table 2, p. 22 | `binding_constant` (dimensionless companion) |
| **per-gram binding constant, hexyl acetate, native primary** | **8.18e-2** | L/g protein | 37 C, pH 7.0-7.2 | nK/220 000 (mine) | **`derived_assumption`** — the registry's own K_g form; MW is stated but cited (Flags 4) |
| per-gram binding constant, hexyl acetate, 100 C preheat | **3.86e-3** | L/g protein | as above | (mine) | `derived_assumption` |
| per-gram binding constant, heptyl acetate, native primary / 100 C | **2.00e-1 / 1.09e-2** | L/g protein | as above | (mine) | `derived_assumption` |
| per-gram binding constants, all 20 rows | see section 3 item 2 | L/g protein | as above | (mine) | `derived_assumption` |
| **PREHEAT FACTOR on nK, simple esters, like-for-like (secondary class only)** | **0.18-0.24** (a 4.2-5.6x fall) | × | 80-100 C preheat vs native | Table 2 (mine) | **`within_study_ratio`** — the defensible ester number |
| preheat factor on nK, simple esters, primary-vs-single-class | **0.043-0.056** (an 18-21x fall) | × | 100 C preheat vs native | Table 2 (mine) | `within_study_ratio` — **overstates the like-for-like fall ~4x (Flags 7)** |
| **PREHEAT FACTOR on nK, terpene esters and terpene alcohols** | **1.07 to 1.71** (a RISE) | × | 80-100 C preheat vs native | Table 2 (mine) | **`within_study_ratio`** — the opposite sign |
| spread of the preheat factor across ligands at one temperature | **0.047 to 1.71, a factor of 36** | × | 100 C preheat | Table 2 (mine) | **`within_study_ratio`** — the refutation of a scalar preheat term |
| preheat factor on the binding ratio (bound fraction), esters | **0.56-0.50** (upper end of range: 75 → 42 %, 88 → 44 %) | × | 100 C preheat vs native | Table 2 (mine) | `within_study_ratio` — **the quantity a model consumer should use** |
| binding ratio, all rows | 7 % to 88 % | % bound | 37 C, 9.21 g/L, 0.012-0.08 mM | Table 2, p. 22 | **`measured_ratio`** |
| **chain-length slope, HpAc / HxAc on native primary K** | **2.86** | × per CH2 | 37 C, aqueous, soy, esters | 40 000/14 000 (mine) | **`within_study_ratio`** — within 2 % of the shipped 2.81 |
| the same, across all six protein-state/class determinations | **1.76 to 3.00**, geometric mean **2.19** | × per CH2 | as above | Table 2 (mine) | `within_study_ratio` |
| the paper's own statement of it | "increased by **3-fold** for each methylene group increment" | × per CH2 | native SPI | p. 19 | `within_study_ratio` |
| fluorescence emission maximum, native / 80 / 90 / 100 C | **332.0 / 335.5 / 336.0 / 338.0** | nm | 0.02 % w/v SPI, 10 mM phosphate pH 7.2, ex 280 nm | Table 1, p. 18 | `level_only` |
| particle size D(4,3), native / 80 / 90 / 100 C | **598.3 / 422.4 / 355.1 / 345.0** | nm | Microtrac S3500, water | Table 1, p. 18 | `level_only` |
| **surface hydrophobicity S0, native / 80 / 90 / 100 C** | **119.6 / 428.7 / 419.3 / 417.2** | dimensionless (no unit printed) | ANS probe, initial slope | Table 1, p. 18 | **`level_only`** — a **3.58x jump at 80 C then a plateau (mine)** |
| **S0 does not track ester binding** | S0 moves **2.7 %** over 80-100 C while HxAc nK falls **12.9x** | — | as above | Tables 1 and 2 (mine) | **`within_study_ratio`** — evidence against a log-P/hydrophobicity term |
| effective protein loading | **9.21** | g/L protein | 1 % w/v × 92.1 % purity, confirmed by the printed 0.0418 mM at 220 000 Da | §2.2, §2.4.1, Fig. captions (mine) | `derived_assumption` (arithmetic on printed values) |
| SPI molecular weight basis | **220 000** | Da | — | §2.2, p. 17 | **`level_only` — `[C]`, Beyeler & Solms 1974; 2.2x Damodaran's 100 000 for the same protein (Flags 4)** |
| **limonene and myrcene are SALTED OUT by SPI** | release **above 100 %** in all four protein states; ordering native > 80 C > 90 C ≈ 100 C, letters a/b/c/c | % of protein-free control | 37 C, 0.024 mM, 9.21 g/L | Fig. 4, p. 22 | **`level_only`** — direction and ordering only; **every magnitude is figure-only (Flags 8)** |
| primary binding sites exist only for HxAc and HpAc, and only on native and 80 C protein | a structural finding, no value | — | 37 C | pp. 19-20, Figs. 1-3 | **`structural_gate`** |
| 7S denaturation temperature / 11S denaturation temperature | **below 80 °C / approximately 95 °C** | °C | — | pp. 18, 20 | `level_only` — **`[C]`, Sorgentini 1995** |

### Can these be put on the same basis as the shipped binding constants, i.e. converted to K_g in L/g?

**Yes — this is the first of the five papers in this batch for which the answer is yes, and by the registry's
own construction. But one input is cited rather than measured, and it is the input that sets the whole scale.**

- **The protein loading is printed**, twice and consistently: 1 % (w/v) in §2.4.1, and 0.0418 mM in every figure
  caption, which back-solves to **9.21 g/L of protein** at the stated molar mass and purity (section 3 item 1).
  This is better provenance than Bi 2022's pea loading, which had to be inherited from a cross-reference.
- **The molar mass is printed**: 220 000 Da. That closes the gap that made Bi 2022's Klotz K unusable in L/g.
  **So K_g = nK / 220 000 runs, exactly as `kg_2_heptanone_soy` = n·K/100 000 does for Damodaran.**
- **And that is precisely where the problem is.** Damodaran's molar mass for **soy protein** was **100 000
  g/mol, stated by source**, and three shipped rows carry `"molar_basis": "stated_by_source"` on that footing.
  Guo's molar mass for **soy protein** is **220 000 Da, cited from a 1974 paper**. **The two differ by 2.2x for
  the same protein species, and a K_g computed from a Klotz K is inversely proportional to it.** Neither figure
  is wrong in its own terms — soy protein isolate is a mixture of 7S (trimer, ~150-180 kDa), 11S (hexamer,
  ~320-360 kDa) and their dissociation products, and any single number is a weighted convention. **But it means
  a Guo K_g and a Damodaran K_g are not on a common scale, and the difference is a bookkeeping choice, not a
  measurement.** If a Guo row is ever shipped it must carry `"molar_basis": "cited_not_measured (220 000 Da,
  Beyeler & Solms 1974)"` and a note that the registry's other soy rows use 100 000.
- **A further wrinkle specific to this paper: the molar mass is a moving target across its own conditions.**
  The paper's central finding is that **preheating DISSOCIATES the protein** — D(4,3) falls from 598 to 345 nm
  and 11S "subunits broken", and the discussion turns on 11S dissociating into "smaller dissociation subunits".
  A fixed 220 000 Da was nevertheless used for the native and all three preheated states. **So the preheated
  rows' K in M^-1 is computed against a molar mass that the paper itself argues no longer applies**, and the
  per-gram conversion inherits that. **The per-gram form is actually MORE defensible than the M^-1 form here**,
  because a per-gram constant does not care how the protein is subdivided — but it is reached by dividing by a
  molar mass, so the error does not cancel. **The clean route would be to refit from the raw peak areas on a
  per-gram basis, which needs data this paper does not print.**
- **What travels without any of this: the RATIOS.** Preheated nK over native nK cancels the molar mass exactly,
  since both legs use the same 220 000. **The preheat factors (0.18-0.24 for the esters like-for-like, 1.07-1.71
  for the terpenoids) and the chain-length ratio (2.86x) are molar-mass-free and are the safest objects in this
  paper.** They are also the objects the repository actually lacks.

**On the B26 temperature limit: this paper does not lift it either, and the reason is worth stating in the
registry.** Every constant here is measured at **37 C**. What varies at 80-100 C is the *history of the
protein*, not the temperature of the measurement. **The correct field to carry that in is not `temperature_c`
— which would be 37.0 — but a new provenance key**, e.g. `preheat_c`, so that a caller can see the difference
between "a pea constant measured at 37 C on native protein" (the Bi 2022 rows) and "a soy constant measured at
37 C on protein that has been to 100 C" (these rows). **Conflating the two would be a category error of the
same family as reading a van 't Hoff enthalpy as an activation energy.**

**Nothing here goes to `matrix_sites.py` as a rate.** No rate constant, no time series (only a single 24 h
equilibration), no activation energy and no thermodynamic quantity of any kind appears in this paper — there is
**no ΔH, ΔS or ΔG anywhere in it**, so the enthalpy/activation-energy confusion cannot arise here. What it
offers `matrix_sites.py` is the *shape* of a preheat correction on a site pool, and the finding that the shape
is ligand-dependent and bidirectional.

## 5. Flags

1. **THE PREHEAT PROTOCOL IS NOT REPORTED.** This is the most serious defect in a paper whose entire subject is
   preheating. Section 2.2 prepares the isolate; section 2.3, headed "Evaluation of characteristics of preheated
   SPI", moves straight to fluorescence, hydrophobicity and particle size. **There is no heating time, no
   protein concentration during heating, no vessel, no heating or cooling rate, and no statement of whether the
   protein was heated as a solution or as a powder** — and the temperatures 80/90/100 C appear first as bare
   column headers in Table 1. Verified against the rendered image of p. 17: the omission is real. **Without a
   time-at-temperature, none of the preheat factors in this paper can be placed on a thermal-history axis, and
   they cannot be compared with Crowther's 20-min autoclave or with Damodaran's 90 C / 1 h.** This alone blocks
   shipping any preheat correction derived here as a parameter.
2. **Table 2's flavour-concentration column is 10x the Methods and the figure captions.** §2.4.1 (p. 17) says
   "0.012 to 0.08 mM"; every figure caption says "0.012-0.080 mM"; **Table 2 prints 0.12-0.24, 0.24-0.8 and
   0.12-0.8 mM.** The Methods and the four figure captions agree with each other against the one table, and a
   dropped decimal point is the obvious explanation, but **the paper never resolves it and the concentration
   axis is what the Klotz fit is against**. If Table 2's range were the true one, every K in M^-1 would be
   wrong by a factor of 10 (Klotz slope 1/(nK[L]) scales with [L]). **Confirm before using any constant.**
3. **Both table footnotes are wrong about their own letters.** Table 2's footnote says "different lower case
   letters (a-c) in the same **row**", but each row carries at most one letter and the letters vary **down the
   nK column within a compound block**; it also says a-c while the HxAc block uses **four** letters, a-d.
   Table 1's footnote says "in the same row", which is correct there, but the **running text on p. 18 contradicts
   Table 1's own letters**: it says "the surface hydrophobicity of PSPI continued to increase" at 90 and 100 C,
   while the table letters mark 80 C as `a` and 90/100 C as `b` — i.e. S0 **peaks at 80 C and then falls
   significantly**. Trust the table.
4. **The molar mass is 220 000 Da, cited from 1974, and it is 2.2x the basis of three shipped soy rows.**
   Detailed in section 4. Damodaran stated 100 000 g/mol for soy protein and `parameters_matrix.py` records
   `"molar_basis": "stated_by_source"` on `kg_2_heptanone_soy`, `kg_2_octanone_soy` and `kg_2_nonanone_soy`.
   Guo takes 220 000 from Beyeler & Solms 1974, apportioned by fraction after Fukushima 1991 and O'Keefe 1991,
   and does not measure it. **Any per-gram constant from this paper is inversely proportional to a number chosen
   by convention, and the convention differs from the one already in the registry by more than a factor of two.**
   Worse, the paper's own central claim is that preheating **dissociates** the protein, yet the same 220 000 is
   used for all four states.
5. **One row does not close arithmetically, and it is the most important one.** HpAc × native SPI, primary
   sites: n = 1.4, K = 40 000, **n × K = 56 000 against a printed nK of 44 000 — a 27 % discrepancy.** The other
   twenty-seven rows close to 8 % or better and eight of them exactly. **This is the highest-affinity value in
   the paper and the numerator of the headline 18x preheat collapse for heptyl acetate.** One of the three
   numbers in that row is mis-transcribed and the printed data cannot say which. Separately, **the paper never
   states the number of replicates behind Table 2's standard errors** — they are propagated from the intercept
   and slope of a single Klotz regression per condition via the Kühn 2007 expressions, so they are regression
   uncertainties, not experimental reproducibility.
6. **Propylene glycol is present in every sample and is never accounted for.** Every flavour was dosed from a
   1 mM stock in **propylene glycol**, chosen explicitly to protect the SPME fibre from solvent competition —
   i.e. chosen for the *instrument*, not for the *protein*. Propylene glycol is a small amphiphilic diol that
   competes for hydrophobic surface and changes the air/water partition of every one of these eight compounds.
   The paper does not state the final co-solvent fraction, does not run a solvent-only control, and never
   mentions it again. (Bi 2022 carries the same defect with methanol at ~1.25 % v/v.)
7. **The native-versus-preheated comparison for the esters is not like-for-like, and the headline number
   inflates the effect ~4x.** Native and 80 C SPI give **two** binding classes for HxAc and HpAc; 90 and 100 C
   give **one**. Comparing the native *primary* nK (18 000) with the 100 C *single-class* nK (850) is comparing
   a high-affinity subpopulation against a whole-population fit and yields 21x. Comparing the two **secondary**
   classes, the only class that plausibly persists, gives **5.6x for HxAc and 4.2x for HpAc**. **The 18-21x
   figure should never be quoted without the class qualification.** The registry's own precedent applies here:
   a constant fitted on a different concentration branch is a different number (the same problem Bi 2022's
   high-branch-only Klotz fit carries).
8. **The salting-out result is real, important and unquantified.** Limonene and myrcene are the only compounds
   in the paper for which the protein **increases** headspace concentration, and this is one of very few
   measured enhancements in the whole corpus (the others being Leksrisompong's two negative caseinate rows).
   But **no numeric release value appears anywhere in the text or in a table** — the entire result lives in
   Fig. 4's bar heights. Per house rule no magnitude is typed here. **The existence and direction of the effect
   can be carried; the size cannot.** Note also that limonene and myrcene were dosed at a **single fixed
   concentration (0.024 mM)**, so no isotherm and no constant exists for them at all.
9. **The medium is unbuffered water.** "1 % (W/V) in water (pH 7.0-7.2)" — the pH is a range, not a set point,
   and there is no buffer to hold it. The residual ionic strength is stated as **0.03-0.04 M as NaCl**, carried
   over from the isolate preparation. Every shipped `REVERSIBLE_BINDING` row with a `ph_of_measurement` was
   measured in a real buffer (30 mM Tris for Damodaran, 10 mM potassium phosphate for Bi). **A `pH` field of
   7.1 on a Guo row would be an average of a drifting quantity**, and the fluorescence assay was run in a
   *different* medium again (10 mM phosphate pH 7.2 at a 50x lower protein loading), so the structural
   characterisation and the binding measurement were not made under the same conditions.
10. **The literature this paper sits in disagrees with it, and the paper says so.** Damodaran & Kinsella 1981b
    found 2-nonanone binding **INCREASED** on preheated SPI (90 C, 1 h) by equilibrium dialysis — the opposite
    sign to Guo's esters, on the same protein species. Semenova 2002 found, on broad-bean 11S, that heating
    **raised the number of sites and lowered the constant**. Heng 2004 found vicilin binding **fell** on heating.
    McNeill & Schmidt 1993 found milk-protein vanillin binding **rose**. Guo's own summing-up (p. 21):
    *"It is difficult to compare results and reach conclusions based on thermally denatured SPI because of
    differences in approaches and experimental conditions used in different studies."* **The corpus therefore
    contains at least four signs for the preheat effect, and Guo's own table contains two of them.** No single
    preheat direction is licensed by the literature.
11. **The relative errors on the preheated constants are enormous.** HxAc at 100 C: K = 100 ± 80 (**80 %**);
    n = 8.5 ± 5.1 (**60 %**). HpAc at 90 C: K = 380 ± 290 (**76 %**). LiFo at 80 C: K = 120 ± 100 (**83 %**).
    Geraniol at 90 C: K = 810 ± 460 (**57 %**). **Every terpenoid K in the table has a standard error above
    40 % of its value**, and the terpenoid preheat "rises" of 1.07-1.71x are mostly inside those errors. The nK
    errors are smaller (11-30 %) and the lettering is applied to nK, which is why the paper argues on nK — but a
    ratio of two nK values with 20 % errors each carries a ~28 % error, so **a 1.07x rise (LiAc at 100 C) is not
    distinguishable from no change.** Only the geraniol (1.71x) and linalool/LiAc-at-80 C (1.44/1.57x) rises
    clear their own error bars.
12. **The 90 → 100 C step is non-monotone for three of six compounds.** HpAc nK rises 1 900 → 2 400; LiFo falls
    1 200 → 1 050 after rising; LiAc falls 3 400 → 2 500 → 2 300 after its 80 C peak. **There is no monotone
    function of preheat temperature here for any ligand class**, and any interpolation between the three
    temperatures is unsupported.
13. **The protein is characterised by three bulk assays and nothing chemical.** Fluorescence maximum, particle
    size and ANS hydrophobicity — no SDS-PAGE, no DSC, no thiol, no disulfide, no free amine, no solubility.
    `data/species/protein_matrices.yml` charges soy site densities from other preparations, so pairing a Guo
    preheat ratio with those densities is a cross-preparation pairing and should be labelled as one. The 7S and
    11S denaturation temperatures the whole mechanism rests on (below 80 C and ≈ 95 C) are **cited from
    Sorgentini 1995, not measured on this isolate**.
14. **English-language and typographical defects are frequent enough to warrant checking any quoted sentence
    against the page.** Examples: "nature SPI" for native SPI throughout the abstract; "NPSI" and "NSPI" both
    used for the native isolate; "binding rato" in Table 2's header; "masa spectrometry" on p. 17; "It guess
    that" on p. 19; "Bluker SCION" for Bruker. None affects a number, but the abstract's own summary sentence —
    "for LiFo, LiAc, geraniol, and linalool, nature < 80 °C < 90 °C < 100 °C PSPI" — **describes a monotone
    increase that Table 2 does not show** (see Flags 12).
15. **What this paper does NOT contain**: any binding measurement above 37 C; any preheat protocol; any
    aldehyde, ketone, 2-alkenal, pyrazine, pyridine, furan or sulfur compound; any thermodynamic quantity
    (no ΔH, ΔS, ΔG); any rate constant or activation energy; any covalent-adduct measurement; any pH other than
    7.0-7.2; any buffer; any time series; any sensory measurement or odour threshold; any numeric limonene or
    myrcene release; any raw Klotz data behind Figs. 1-3; any replicate count for Table 2.
16. **What to request from the authors**: (i) **the preheating protocol — time, concentration, vessel, cooling**
    — without which nothing here is transferable to a thermal programme; (ii) whether Table 2's concentration
    column or the Methods' range is correct; (iii) the corrected HpAc native-primary row; (iv) the numeric
    limonene and myrcene release values behind Fig. 4; (v) the replicate count behind Table 2's standard errors;
    (vi) the final propylene-glycol fraction in the assay; (vii) the raw peak areas, which would allow a
    per-gram refit free of the 220 000 Da assumption.
