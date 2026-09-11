# Damodaran & Kinsella 1981 — EXTRACTION (four methyl ketones and one aldehyde bound to 1 % native whole soy protein by equilibrium dialysis, 30 mM Tris-HCl pH 8.0 + 10 mM 2-mercaptoethanol, at 5 / 25 / 45 C; seven intrinsic binding constants and free energies on a stated 100 000 g/mol basis)

### THE MATRIX LAYER'S OWN SUBJECT, MEASURED ON SOY: seven (n, K, ΔG) triples with the 100 kDa basis PRINTED rather than inferred — the strongest provenance in the binding batch, and the paper the repository has been citing at second hand for four of its five soy constants while leaving the 5 C row, the 5-nonanone row and the succinylated row on the floor.

**Source on disk:** `data/articles/damodaran1981.pdf` — the article itself is 5 pp. (J. Agric. Food
Chem. **1981**, 29 (6), 1249-1253), and the file continues into **page 1253 of the immediately
following paper by the same two authors** ("Interaction of Carbonyls with Soy Protein:
**Conformational** Effects", pp. 1253-1257, of which only the first page and abstract are present).
**Do not confuse the two**: the 11S/7S fractionation, the urea work and the fluorescence are the
*second* paper's; only the succinylated row of Table I here overlaps.
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/damodaran1981.txt`). The text layer is a clean two-column OCR;
**Table I came through complete and legible** and is re-typed below without reconstruction. The four
figures — 1A/1B (binding isotherms and double-reciprocal plots for the three methyl ketones at
25 C), 2 (double-reciprocal plots for nonanal / 2-nonanone / 5-nonanone), 3A/3B (the temperature
series on 2-nonanone at 5, 25 and 45 C) and 4 (native vs partially denatured soy) — are images:
**every binding isotherm and every point of ν versus free ligand in this paper is figure-only**, and
the constants in Table I are the fitted intercepts and slopes of those plots. There is no
supplementary material. Repo status before this dossier: Damodaran 1981 is cited by
`src/kinetic_core/parameters_matrix.py` in five `REVERSIBLE_BINDING` rows (lines 378-418) as
"damodaran1981 via k2_matrix_and_thresholds.md sec. (b)", and `k2_matrix_and_thresholds.md` states
in its own header that "the soy binding constants and the 100 kDa basis are NOT re-derived here —
they are cited and extended", pointing at a `damodaran1981_extraction.md` that **did not exist on
disk**. This file closes that loop.

## 0. Identity

| field | value |
|---|---|
| Title | "Interaction of Carbonyls with Soy Protein: Thermodynamic Effects" |
| Authors | Srinivasan Damodaran (corresponding) and John E. Kinsella — Institute of Food Science, Stocking Hall, Cornell University, Ithaca, New York 14853 |
| Venue | J. Agric. Food Chem. **1981**, **29** (6), 1249-1253. Received for review 14 May 1981; accepted 20 August 1981. Supported in part by NSF grant CPE 80-18394 |
| DOI / article ID | none printed; the ACS code line reads `0021-8561/81/1429-1249$01.25/0` |
| Protein | **whole soy protein**, prepared in-house from **defatted, low-heat-treated soy flour** (Central Soya, Chicago; lot 878, code 3040) by isoelectric precipitation. **Not a commercial isolate** |
| Method | **equilibrium dialysis** (Spectrapor-2 membrane, matched acrylic half-cells), free ligand extracted into isooctane and quantified by **GC-FID** |
| Basis | **"the number of moles of ligand bound per 100 000 g of soy protein"** — the 100 kDa basis is stated in the Methods, in the Figure 1 legend and again in the Discussion. **It is printed, not inferred** |
| Naming | `ν` (rendered variously as `v`, `P`, `i,` and `t` by the OCR) = moles bound per mole of protein; `n` = total binding sites; `K` = intrinsic binding constant; `ΔG` in **kcal/mol** |
| Companion paper | Damodaran & Kinsella 1981b, "Conformational Effects", the immediately following article — the 11S fraction binds 2-nonanone almost not at all while the 7S fraction binds like whole soy; urea and succinylation both reduce binding |
| Earlier work by the same pair | Damodaran & Kinsella 1980, J. Agric. Food Chem. **28**, 567 — the same experiment on **bovine serum albumin**, quoted here for K(2-nonanone) ≈ 1 800 M^-1 and six binding sites |
| Second-hand content | Arai et al. 1970 (hexanal and 1-hexanol on partially denatured soy by gel filtration); Beyeler & Solms 1974 (2-butanone on a soy isolate suspension); Abraham 1980 (-540 cal/mol per CH2 in model solvents) |

## 1. Why it matters

The protein-matrix layer has two halves and this paper feeds both.

**(a) The reversible-binding half.** `REVERSIBLE_BINDING` in
`src/kinetic_core/parameters_matrix.py` is the table of per-gram constants `k_g` in L/g that the
matrix layer uses to size how much of an odourant a protein pot holds back. Five of its rows are
this paper, and **the arithmetic that produced them is `k_g = n·K / MW` with MW = 100 000 g/mol**:

| repository row | value in the code | this paper's numbers | check (mine) |
|---|---|---|---|
| `kg_2_heptanone_soy` | 4.40e-3 L/g | n = 4, K = 110 M^-1 | 4 x 110 / 100 000 = **4.40e-3** ✓ |
| `kg_2_octanone_soy` | 1.24e-2 L/g | n = 4, K = 310 M^-1 | 4 x 310 / 100 000 = **1.24e-2** ✓ |
| `kg_2_nonanone_soy` | 3.72e-2 L/g | n = 4, K = 930 M^-1 | 4 x 930 / 100 000 = **3.72e-2** ✓ |
| `kg_nonanal_soy` | 4.38e-2 L/g | n = 4, K = 1 094 M^-1 | 4 x 1 094 / 100 000 = **4.376e-2** ✓ |
| `kg_hexanal_soy_denatured` | 1.47e-3 L/g | **Arai 1970, quoted here**: K = 173.4 M^-1, 0.847 mg hexanal bound per g protein at saturation | (0.847e-3 / 100.16) x 173.4 = **1.466e-3** ✓ — the route is via the printed **mg/g** saturation figure, which is why the code's provenance says `per_gram_no_molar_mass_needed` |

**All five reproduce exactly.** The `k_g` values in the code are correct arithmetic on this paper's
printed table, and the 100 kDa basis is the paper's own, not an assumption bolted on downstream.
That matters because the corresponding beta-lactoglobulin rows (`kg_*_blg`, from Andriot 2000) carry
`"molar_basis": "recovered_by_arithmetic (36 800 g/mol dimer)"` — a recovered basis — whereas these
carry `"molar_basis": "stated_by_source"`, and now that statement can be pointed at a sentence.

**(b) The covalent-binding half, by contrast.** `BINDING_CLASSES` in
`src/kinetic_core/matrix_sites.py` is a different quantity: **second-order rate constants** in
M^-1 s^-1 for aldehydes and HMF consuming amine and thiol **sites**, sourced from the adduct-kinetics
synthesis. Nothing in this paper is a rate. More sharply: **this paper's chemistry is deliberately
the complement of that table's.** The dialysis buffer carries **10 mM 2-mercaptoethanol** throughout,
which keeps the protein's cysteines reduced and, more to the point for the aldehyde row,
**out-competes the protein for any aldehyde that would otherwise form a thiol adduct**. So
`kg_nonanal_soy` is a purely reversible, hydrophobic constant that **excludes** exactly the covalent
chemistry `matrix_sites.py` prices — which is what the code's own note on that row says
("DIALYSIS + 2-mercaptoethanol: this constant EXCLUDES the cysteine-aldehyde chemistry a headspace
determination would count. It must never be pooled with a headspace aldehyde value"). **Confirmed
from the Methods; the note is right and it is load-bearing.**

**What this dossier adds that the code does not have.** Three of Table I's seven rows are on the
floor: the **5 C** measurement (K = 2 000 M^-1 with n dropping to **2**), the **5-nonanone** row
(K = 541 M^-1, the internal-keto control), and the **succinylated** row (K = 850 M^-1 with n = 2).
The partially denatured row (K = 1 240) is used only through Arai's hexanal, not in its own right.
And the paper's **temperature behaviour** — flat between 25 and 45 C, then a discontinuity at 5 C —
is the only direct statement in the corpus about whether a soy binding constant may be carried across
the temperature range a cook spans. It says: **between 25 and 45 C, yes; below 25 C, no.**
`REVERSIBLE_BINDING` rows carry a single `temperature_c` field (25.0 for all four soy rows) and no
temperature dependence at all; this paper is the evidence that, over the warm part of the range, none
is needed.

## 2. Methods as they matter to a model

- **The protein, and it is not a commercial isolate.** Whole soy protein prepared from
  **defatted, low-heat-treated soy flour** by isoelectric precipitation: extracted with **30 mM
  Tris-HCl pH 8.0 containing 10 mM 2-mercaptoethanol**, meal:buffer **1:20**, centrifuged, the
  supernatant adjusted to **pH 4.8 with 2 N HCl**, the precipitate redissolved in the Tris buffer,
  **dialysed against water at pH 8.0 for 24 h**, and lyophilised. The authors are explicit that the
  point is to work on **native** protein, because "studies with denatured proteins may not be
  indicative of the true molecular nature of interaction". This is a different material from the
  commercial isolates behind `data/species/protein_matrices.yml`'s `soy_isolate` entry (Ruan 2014,
  Shimada 1988, Xiao 2024, Jaeger 2023, Gorissen 2018), all of which are partially denatured.
- **The buffer, in full, because two of its components change the chemistry.** **30 mM Tris-HCl,
  pH 8.0, 10 mM 2-mercaptoethanol, 0.02 % sodium azide, protein at 1 % (w/v) in every experiment.**
  Protein concentration by A280 with **an absorptivity of 8.02 for a 1 % solution** (Thanh &
  Shibasaki 1976). **pH 8.0 is high** — above the `PH_ADDUCT_GATE` thresholds in
  `parameters_matrix.py`, and high enough that lysine epsilon-amines are substantially unprotonated,
  which would *favour* Schiff-base chemistry were the thiol reductant not there.
- **Equilibrium dialysis, exactly.** Matched acrylic half-cells separated by a **Spectrapor-2**
  membrane and clamped. **3 mL of protein solution** one side, **3 mL of buffer containing a known
  amount of ligand** the other. Shaken **at the required temperature for at least 18 h** to reach
  equilibrium. Then **1 mL from each side** into a vial with **1 mL isooctane**, shaken to extract.
  The authors justify quantitative extraction two ways: the aqueous/isooctane partition coefficient
  of the carbonyls is very large (the printed value is lost in the scan — "of the order of ___"), and
  **a second extraction of the protein-side aqueous phase recovered no ligand**. **The difference in
  ligand concentration across the membrane is the bound amount; the buffer side is the free
  concentration [L].**
- **The fit.** Double-reciprocal, **1/ν = 1/n + 1/(nK[L])**. The intercept gives **n**, the total
  number of binding sites; the slope gives **1/(nK)**. So **n and K are not independently measured —
  they are an intercept and a slope from the same straight line**, and a set of ligands that share an
  intercept share an n by construction of the reading, not by separate assay.
- **Quantification.** Perkin-Elmer Model 900 GC with FID; a stainless-steel column (diameter and
  length garbled in the scan: "___ in. diameter and 10___%   length") packed with **10 % Apiezon on
  Chromosorb**; hydrogen 15 mL/min, oxygen 300 mL/min, nitrogen 40 mL/min. **No internal standard is
  mentioned and no calibration is described.**
- **Ligands and purities.** 2-Nonanone (99+ %), 2-octanone (98 %), nonanal (98 %) and 5-nonanone
  from Aldrich; **2-heptanone's supplier and purity are not stated** although it is in Table I.
  Purity checked by GC. Spectral-grade isooctane (Fisher). Distilled deionised water throughout.
- **The three perturbations.**
  - **Temperature**: the same 2-nonanone experiment at **5, 25 and 45 C** (Figs. 3A/3B).
  - **Partial denaturation**: "**heating a 1 % solution at 90 C for 1 h**" (Fig. 4 legend), then the
    same dialysis at 25 C.
  - **Succinylation**: no method is given in this paper; the succinylated row of Table I is the only
    trace of it, and the procedure belongs to the following paper.
- **What is never measured here.** No site chemistry of any kind — **no thiol assay, no amine assay,
  no lysine count, no free-amino-group determination**. The "four binding sites per 100 kDa" is a
  curve intercept, not a chemical inventory. No kinetics: every number is an equilibrium constant,
  and the paper never asks how long binding takes beyond the 18 h it allows for it. No irreversible
  or covalent binding is looked for; the whole treatment assumes reversibility and the Discussion
  closes on that assumption ("since the binding process is reversible, it should be possible to
  remove the off-flavors under appropriate conditions").

## 3. Tables re-typed

### Table I. "Thermodynamic Constants for the Binding of Carbonyls to Soy Protein at 25 °C"

Column headings exactly as printed: `ligand | type of soy preparation | n | K, ~-5 | ΔG, kcal/mol`.
**The K column header is damaged in the scan** — it reads `K, ,` over `~-5` where the journal
almost certainly prints `K, M^-1` (see Flags 2; the running text quotes these same values as
"930 M^-1", "1240 M^-1", "1800 M^-1", so **the unit is M^-1 and the values are as printed, with no
power-of-ten multiplier**).

| ligand | type of soy preparation | n | K (M^-1) | ΔG (kcal/mol) |
|---|---|---:|---:|---:|
| 2-heptanone | native | 4 | 110 | -2.781 |
| 2-octanone | native | 4 | 310 | -3.395 |
| 2-nonanone | native | 4 | 930 | -4.045 |
| 2-nonanone | part. denatured | 4 | 1240 | -4.215 |
| 2-nonanone | succinylated | 2 | 850 | -3.992 |
| 5-nonanone | native | 4 | 541 | -3.725 |
| nonanal | native | 4 | 1094 | -4.141 |

**All at 25 °C.** The table title says so and no other temperature appears in it.

### Numbers printed in the running text but NOT in Table I

| quantity | value | conditions | where | whose measurement |
|---|---|---|---|---|
| **2-nonanone K at 5 °C** | **2 000 M^-1** | 1 % native soy, 30 mM Tris pH 8.0, 10 mM 2-ME | Results, "Effect of Temperature" | **this paper** (the plot is Fig. 3B) |
| **2-nonanone n at 5 °C** | **2** binding sites (against 4 at 25 and 45 C) | as above | same | this paper |
| 2-nonanone ΔG at 5 °C | **-4.221 kcal/mol** | as above | Results, "Effect of Temperature" | this paper |
| 2-nonanone K at 45 °C | **930 M^-1** — "at 25 and 45 °C it is only 930" | as above | same | this paper (**identical to 25 C; the slope and intercept are stated to be identical**) |
| binding-site count, general | "about **four to five** binding sites for methyl ketones in the native soy protein (on the basis of 100 000 molecular weight)" | 25 C | Results, first paragraph | this paper (**the table says 4; the prose says 4-5**, Flags 3) |
| K increment per CH2 | "the binding constant increases **3-fold** for each increment in the chain length" | 25 C | Results | this paper |
| ΔG increment per CH2, this study | "about **-600 cal/CH2 residue**" | 25 C | Results | this paper (**my arithmetic on Table I gives -632**, §3 below) |
| ΔG increment per CH2, BSA | "about **-550 cal/CH2 residue**" | — | Results | Damodaran & Kinsella **1980**, quoted |
| ΔG increment per CH2, model solvents | "about **-540 cal/mol of CH2**" | water -> apolar solvent | Results | **Abraham 1980**, quoted |
| keto-position penalty | "the hydrophobic free energy of association becomes more positive by **105 cal/mol**" for each shift of the keto group from position 1 toward the centre | 25 C | Results | this paper |
| 2-nonanone on **bovine serum albumin** | **≈1 800 M^-1**, **six** binding sites, "highly hydrophobic" | — | Discussion | Damodaran & Kinsella 1980, quoted |
| denaturation effect | the binding affinity "increases by about **30 %**" while the number of sites does not change | 90 C / 1 h on a 1 % solution, measured at 25 C | Results | this paper (1240/930 = **1.33x**, mine) |
| **hexanal** on partially denatured soy | **K = 173.4 M^-1**; **0.847 mg bound per g protein** at saturation | gel filtration, partially denatured soy | Discussion | **Arai et al. 1970**, quoted |
| **1-hexanol** on partially denatured soy | **K = 80.3 M^-1**; **0.889 mg/g** at saturation | as above | Discussion | Arai et al. 1970, quoted |
| sites implied by Arai's numbers | "about **one**" per 100 000 g | Damodaran's own recalculation of Arai | Discussion | derived by **this paper** from Arai (**and rejected by it** as "very low") |
| 2-butanone on a soy isolate | **5 174 M^-1** at 20 °C, pH 7.0; molal ratios "up to **1 200**" on an assumed **50 000** g/mol | equilibrium dialysis on insoluble suspensions | Discussion | **Beyeler & Solms 1974**, quoted **and refuted** (Flags 6) |
| hexanal on **native** soy | "it may be **speculated** that the binding constant for hexanal would be about **40 M^-1**" | 25 C, by extrapolating the 3-fold-per-CH2 rule | Discussion | **this paper, explicitly as speculation** — Flags 5 |

**Binding isotherms: FIGURE-ONLY.** Figures 1A, 1B, 2, 3A, 3B and 4 hold every measured (ν, [L])
pair in the study — three methyl ketones at 25 C, three nonanone isomers/nonanal at 25 C, 2-nonanone
at three temperatures, and native vs denatured soy. Per house rule they are not typed as numbers.
Table I and the running-text values above are the entire numeric content of the paper.

### Arithmetic on the printed constants (all mine)

**1. The chain-length rule, checked against the table.** K: 110 -> 310 -> 930 for C7 -> C8 -> C9
methyl ketones, i.e. **x2.82 then x3.00 per CH2**, geometric mean **x2.91**. ΔΔG: -3.395 - (-2.781)
= **-614 cal/mol**; -4.045 - (-3.395) = **-650 cal/mol**; mean **-632 cal/mol per CH2**. The paper's
prose rounds this to "-600"; **its own abstract quotes 550, which is the bovine-serum-albumin number
from the authors' 1980 paper, not this one** (Flags 1).

**2. Internal consistency of ΔG against K.** ΔG = -RT ln K with R = 1.987e-3 kcal/mol/K and
T = 298.15 K gives: 2-heptanone -2.784 (printed -2.781); 2-octanone -3.399 (-3.395); 2-nonanone
-4.049 (-4.045); part. denatured -4.220 (-4.215); succinylated -3.996 (-3.992); 5-nonanone -3.728
(-3.725); nonanal -4.146 (-4.141). **Every row reproduces to within 5 cal/mol, a constant offset
consistent with the authors using T = 298 K exactly. Table I is internally consistent and the ΔG
column carries no information beyond the K column** — it is not a second measurement.

**3. The keto-position penalty, checked.** ΔG: nonanal (carbonyl at C1) **-4.141**, 2-nonanone
**-4.045**, 5-nonanone **-3.725**. From position 1 to 2: **+96 cal/mol**. From 2 to 5:
+320 cal/mol over three positions = **+107 cal/mol per position**. Over the whole span 1 -> 5:
416/4 = **+104 cal/mol per position**. **The printed 105 cal/mol is reproduced.** In K terms,
nonanal binds **1.18x** more tightly than 2-nonanone and **2.02x** more tightly than 5-nonanone.

**4. The aldehyde-vs-ketone comparison, which is what the odour layer wants.** At equal chain length
(C9) and equal n, **nonanal K = 1 094 against 2-nonanone K = 930: a ratio of 1.18x.** Per gram,
4.38e-2 vs 3.72e-2 L/g. **The aldehyde is barely a better binder than the ketone under these
conditions** — and that is the signature of a purely hydrophobic interaction with the carbonyl
chemistry suppressed. Any determination that finds an aldehyde binding a soy protein *far* harder
than the matched ketone is measuring something this experiment was designed to exclude.

**5. Denaturation, and the direction it runs.** 1 240 / 930 = **1.33x** stronger binding after
90 C / 1 h, with **n unchanged at 4**. So heating a soy pot makes it hold odourants **more** tightly,
by about a third, without creating new sites. The corresponding per-gram constant, not currently in
the code, would be **4 x 1 240 / 100 000 = 4.96e-2 L/g (mine)** for 2-nonanone on heat-denatured soy.

**6. The 5 °C discontinuity, in per-gram terms.** At 5 C, n = 2 and K = 2 000, so
**n·K = 4 000 against 3 720 at 25 C — a per-gram constant of 4.0e-2 L/g, only 1.08x the 25 C value
(mine)**, even though the *intrinsic* constant more than doubles. **This is the number to quote if
anyone asks whether the cold row matters**: it is a large change in the microscopic picture (half as
many sites, twice as tight) and almost no change in the macroscopic capacity. The paper reports the
K and the n separately and never multiplies them; doing so is the honest way to compare the cold row
with the warm ones.

**7. Temperature transport, and its limit.** K(2-nonanone) is **930 at both 25 and 45 C** — the
authors state the slope and intercept are identical. From ΔG = ΔH - TΔS, a K that does not move over
20 K implies **ΔH ≈ 0**, which is exactly what the paper argues (hydrophobic association is entropy-
driven). **Consequence for the engine: within 25-45 C, the four soy `k_g` values need no temperature
correction at all, and that is a measured statement, not an assumption.** Above 45 C nothing is
measured, and the 90 C denaturation experiment shows that heating far enough changes the protein
rather than the equilibrium — so extrapolating these constants to a cook temperature is not licensed
by this paper.

**8. Soy against BSA, the paper's own cross-protein comparison.** 2-nonanone: **930 M^-1 on soy with
4 sites; ≈1 800 M^-1 on BSA with 6 sites.** Ratio of intrinsic constants **1.94x**; ratio of
capacities n·K **2.90x**. The authors attribute it to BSA's sites being more hydrophobic. Note for
the repository: `kg_2_nonanone_blg` from Andriot 2000 is **6.63e-2 L/g** against this paper's
**3.72e-2 L/g** for the same ligand on soy — a factor of **1.78x** across two proteins, two methods
and two laboratories nineteen years apart, which is the "agrees ... to 1.8x" the code's note on
`kg_2_nonanone_soy` claims. **Reproduced (mine): 6.63e-2 / 3.72e-2 = 1.78.**

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** **None of this paper's five ligands is
keyed.** `2_heptanone`, `2_octanone`, `2_nonanone`, `5_nonanone` and `nonanal`— of these only
**`nonanal` is present** (id `nonanal`); the four ketones are absent, as is `2_butanone` and
`1_hexanol` from the quoted work. `hexanal` is keyed and is the compound of the Arai row. Every row
below shares: **1 % (w/v) whole native soy protein prepared by isoelectric precipitation from
defatted low-heat soy flour, in 30 mM Tris-HCl pH 8.0 with 10 mM 2-mercaptoethanol and 0.02 % sodium
azide; equilibrium dialysis across a Spectrapor-2 membrane, 3 + 3 mL, shaken ≥18 h; free ligand
extracted into isooctane and read by GC-FID; double-reciprocal fit; basis 100 000 g protein per
mole.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| 2-heptanone on native soy, intrinsic constant | 110 | M^-1 | **25 C**, pH 8.0 | Table I p. 1250 | **binding_constant** |
| 2-heptanone, binding sites | 4 | per 100 000 g protein | 25 C | Table I | **binding_constant** (an intercept, Flags 3) |
| 2-heptanone, free energy | -2.781 | kcal/mol | 25 C | Table I | binding_constant (derived from K by the authors, §3.2) |
| 2-octanone on native soy | 310 | M^-1 | 25 C | Table I | **binding_constant** |
| 2-octanone, sites / ΔG | 4 / -3.395 | — / kcal/mol | 25 C | Table I | binding_constant |
| 2-nonanone on native soy | 930 | M^-1 | **25 C** | Table I | **binding_constant** |
| 2-nonanone on native soy | 930 | M^-1 | **45 C** — "at 25 and 45 °C it is only 930"; slope and intercept identical | Results, "Effect of Temperature" | **binding_constant** — the temperature-invariance row |
| 2-nonanone on native soy | **2 000** | M^-1 | **5 C**, with **n = 2** | Results, "Effect of Temperature"; plot Fig. 3B | **binding_constant** — **NOT currently in the repository** |
| 2-nonanone at 5 C, free energy | -4.221 | kcal/mol | 5 C | Results | binding_constant |
| 2-nonanone, sites / ΔG at 25 C | 4 / -4.045 | — / kcal/mol | 25 C | Table I | binding_constant |
| 2-nonanone on **partially denatured** soy | 1 240 | M^-1 | 25 C; protein heated **90 C / 1 h at 1 %** beforehand | Table I; Fig. 4 | **binding_constant** — **not currently in the repository** |
| 2-nonanone, part. denatured: sites / ΔG | 4 / -4.215 | — / kcal/mol | 25 C | Table I | binding_constant |
| 2-nonanone on **succinylated** soy | 850 | M^-1, with **n = 2** | 25 C; succinylation method not given in this paper | Table I | **binding_constant** — **not currently in the repository**; the method belongs to the companion paper |
| 2-nonanone, succinylated: ΔG | -3.992 | kcal/mol | 25 C | Table I | binding_constant |
| **5-nonanone** on native soy | 541 | M^-1 | 25 C | Table I | **binding_constant** — **not currently in the repository**; the internal-keto control |
| 5-nonanone, sites / ΔG | 4 / -3.725 | — / kcal/mol | 25 C | Table I | binding_constant |
| **nonanal** on native soy | 1 094 | M^-1 | 25 C | Table I | **binding_constant** — and see the 2-ME caveat below |
| nonanal, sites / ΔG | 4 / -4.141 | — / kcal/mol | 25 C | Table I | binding_constant |
| per-gram constant, 2-heptanone / 2-octanone / 2-nonanone / nonanal | 4.40e-3 / 1.24e-2 / 3.72e-2 / 4.38e-2 | L/g | 25 C; n·K / 100 000 | derived from Table I (mine); **already in `REVERSIBLE_BINDING`** | binding_constant (derived arithmetic on printed values) |
| per-gram constant, 5-nonanone | **2.16e-2** | L/g | 25 C | derived (mine) | binding_constant — **available and unused** |
| per-gram constant, 2-nonanone on heat-denatured soy | **4.96e-2** | L/g | 25 C after 90 C / 1 h | derived (mine) | binding_constant — **available and unused** |
| per-gram capacity, 2-nonanone at 5 C | **4.0e-2** | L/g | 5 C (n·K = 2 x 2 000) | derived (mine) | binding_constant — **available and unused** |
| chain-length effect | 2.91 (K), -632 (ΔG) | x per CH2, cal/mol per CH2 | 25 C, C7-C9 methyl ketones | derived from Table I (mine) | within_study_ratio |
| carbonyl-position effect | +105 | cal/mol per position from C1 toward the centre | 25 C, C9 series | Results; reproduced as +104 (mine) | within_study_ratio |
| aldehyde vs matched methyl ketone | 1.18 | x (nonanal / 2-nonanone) | 25 C, C9 | derived (mine) | within_study_ratio |
| heat denaturation effect | 1.33 | x on K, with n unchanged | 90 C / 1 h at 1 %, measured at 25 C | derived (mine); the paper says "about 30 %" | within_study_ratio |
| succinylation effect | 0.91 x on K but **0.46 x on n·K** | — | 25 C | derived (mine) | within_study_ratio |
| temperature invariance, 25 -> 45 C | 1.00 | x on K and on n | 2-nonanone | Results | within_study_ratio — **the licence for carrying the 25 C rows to 45 C, and only to 45 C** |
| soy vs bovine serum albumin | 1.94 x on K, 2.90 x on n·K | — | 2-nonanone, 25 C | Discussion (BSA value from Damodaran & Kinsella 1980) | within_study_ratio (**cross-study**; the BSA number is not measured here) |
| hexanal on **partially denatured** soy | 173.4 M^-1; 0.847 mg/g at saturation | M^-1; mg/g | gel filtration, not dialysis; temperature not stated in the quotation | Discussion (**Arai et al. 1970**) | **binding_constant** — second-hand; the source of `kg_hexanal_soy_denatured` = 1.47e-3 L/g |
| 1-hexanol on partially denatured soy | 80.3 M^-1; 0.889 mg/g | M^-1; mg/g | as above | Discussion (Arai 1970) | **binding_constant** — second-hand, **not in the repository** |
| 2-butanone on a soy isolate suspension | 5 174 M^-1 at 20 C, pH 7.0 | M^-1 | insoluble suspension, ν/C slope, MW assumed 50 000 | Discussion (**Beyeler & Solms 1974**) | **REFUSE** — second-hand **and refuted in place** (Flags 6) |
| hexanal on **native** soy | ~40 | M^-1 | 25 C, by extrapolating 3-fold-per-CH2 from 2-heptanone | Discussion, explicitly "may be speculated" | **derived_assumption** — the authors' own speculation, never measured (Flags 5) |
| binding isotherms and double-reciprocal plots, all ligands and all temperatures | — | ν vs [L] | 5 / 25 / 45 C | Figs. 1A, 1B, 2, 3A, 3B, 4 | **figure_only** |

### What can and cannot be put next to these

**(a) They are soy, native, at pH 8.0, with a thiol reductant present — three qualifications, each of
which matters.** *Native*: the repository's `soy_isolate` site densities come from commercial
isolates that are partially denatured, and this paper measures that the denatured form binds **1.33x
harder**. *pH 8.0*: higher than any food pot the engine cooks, and above the `PH_ADDUCT_GATE`
thresholds. *2-mercaptoethanol at 10 mM*: a large excess of a competing thiol over the protein's own
free cysteines (soy isolate carries ~0.0067 mmol free SH per gram, so a 1 % solution has ~0.067 mM
protein thiol against 10 mM 2-ME — **a 150-fold excess of the competitor, mine**). Any aldehyde-thiol
adduct chemistry is comprehensively suppressed.

**(b) The four per-gram constants already in the code are correct and now have a primary source.**
No value changes. What changes is that `"molar_basis": "stated_by_source"` can be pointed at three
separate sentences, and that the exclusion note on `kg_nonanal_soy` can be pointed at the buffer
composition.

**(c) The temperature question the code leaves open is answered over part of its range.** All four
soy rows carry `temperature_c = 25.0` and no temperature model. This paper measures that
K(2-nonanone) is unchanged at 45 C. **Carry the 25 C constants to 45 C on this evidence; do not carry
them above it** — the 90 C experiment shows the protein itself changes, and nothing between 45 and
90 C was measured.

**(d) What cannot be transported.** Nothing here reaches beta-lactoglobulin, pea, or any real food;
nothing here is a rate; nothing here bears on covalent adduction, on thiols, or on Schiff bases;
and the four-site count is a property of a double-reciprocal intercept on a 100 kDa basis, not a
count of chemical groups — **it must never be used as a site density in
`data/species/protein_matrices.yml`**, whose entries are lysine and cysteine counts in mmol per gram
and are a completely different quantity. (For scale: 4 sites per 100 000 g is **0.04 mmol/g**,
against a lysine density of **0.36 mmol/g** for soy isolate — the hydrophobic binding sites are about
a ninth as numerous as the amines, mine.)

## 5. Flags

1. **The abstract contradicts the body twice, and one of the contradictions is a factor of 300.**
   (i) The abstract says "The binding constant increased with the chain length of the ligand by
   **3 orders of magnitude** for each methylene group increase in the chain." The body says
   "**3-fold**", Table I shows 110 -> 310 -> 930, and the whole hydrophobic argument depends on it
   being 3-fold (Wishnia's prediction is "factors of 2-3"). **The abstract is wrong; "orders of
   magnitude" should read "fold".** (ii) The abstract says "The favorable change in the hydrophobic
   free energy was **550 cal/mol** of CH2 residue", but 550 is the **bovine serum albumin** value
   from the authors' 1980 paper; this study's own value is stated as "about -600" in the body and
   computes to **-632** from Table I. **Quote the body and the table, never the abstract.**
2. **The unit on Table I's K column is damaged in the scan.** The header renders as `K, ,` over
   `~-5`, which is not readable as printed. It is resolved unambiguously by the running text, which
   quotes the very same values with units — "the binding constant for native soy is 930 **M^-1**",
   "for partially heat denatured soy is 1240 **M^-1**", "2-nonanone to bovine serum albumin is about
   1800 **M^-1**". **K is in M^-1 and the tabulated numbers carry no multiplier.** The `~-5` fragment
   is most likely the remains of a "M^-1" that the OCR mangled; it is **not** read as "x 10^-5", which
   would make every ΔG in the table wrong by ~7 kcal/mol. The ΔG cross-check in §3.2 settles it: with
   K in plain M^-1 every row reproduces to 5 cal/mol.
3. **"n = 4" is a curve intercept, and the paper is inconsistent about it by one site.** Table I
   prints 4 for six of seven rows; the Results text says "about **four to five** binding sites";
   the Discussion says "about four". More importantly, **n and K come from the intercept and slope of
   the same double-reciprocal line**, so they are correlated, and the statement that several ligands
   "share the same intercept" is made by eye on a figure, not by a statistical test. **No error bar,
   confidence interval or replicate count is given for any n or any K in this paper.** That is the
   single largest gap: seven binding constants, no uncertainty on any of them.
4. **10 mM 2-mercaptoethanol is in every buffer, and it changes what the aldehyde row means.** It is
   there to keep the soy protein reduced during extraction and dialysis, but it is also a 10 mM
   nucleophilic thiol sitting in the same cell as the nonanal. The `kg_nonanal_soy` note in
   `parameters_matrix.py` already says this constant "EXCLUDES the cysteine-aldehyde chemistry"; it
   is confirmed. **A corollary the code does not say: the same buffer means the four-site count is a
   count of *hydrophobic* sites on a *reduced* protein, and gives no information about the covalent
   capacity the `matrix_sites.py` pools represent.**
5. **The "~40 M^-1 for hexanal" figure is the authors' speculation and must not be ingested.** It
   appears once, in the Discussion, prefaced by "it may be **speculated** that", and is obtained by
   running the 3-fold-per-CH2 rule down from 2-heptanone. **It is not a measurement, has no error,
   and is a ketone rule applied to an aldehyde** — which this paper's own nonanal/2-nonanone
   comparison shows is worth another 1.18x. Class it `derived_assumption` and prefer the Arai value
   (173.4 M^-1, denatured soy, gel filtration) if a hexanal-on-soy number is needed, with the caveat
   that Arai's is a different protein state and a different method.
6. **Beyeler & Solms 1974's 5 174 M^-1 is quoted here only to be refuted, and the refutation is
   sound.** Damodaran's objections: (i) they took K as the **slope of the binding isotherm** rather
   than the reciprocal of the free concentration at half-saturation, which is not the definition of
   an equilibrium binding constant; (ii) their molal ratios reach **1 200 mol ligand per mole of
   protein** on an assumed 50 kDa protein, i.e. about **3 molecules of 2-butanone per amino acid
   residue**, which is not binding but adsorption or entrapment; (iii) they worked on **insoluble
   protein suspensions**. **Do not ingest the 5 174 number from any source.** The general lesson
   transfers to any headspace-depletion determination on a suspension.
7. **The succinylated row is orphaned.** Table I gives K = 850 and n = 2 for succinylated soy, but
   **this paper describes no succinylation procedure at all** — the method is in the following
   article. If the succinylated row is ever used, it must be sourced to the companion paper
   (Damodaran & Kinsella 1981b, J. Agric. Food Chem. 29:1253-1257), whose first page is in this PDF
   and whose abstract says succinylation "profoundly affected both the binding affinity and the
   binding capacity".
8. **The 5 C result is a structural change, not a binding measurement, and the authors say so.** n
   falls from 4 to 2 and K rises from 930 to 2 000; the explanation offered is that the 11S fraction
   precipitates at low temperature and the subunits reorganise. **It is evidence that this protein is
   not the same object at 5 C, and it should be carried as such rather than as a cold-temperature
   value of the same constant.** The per-gram capacity barely moves (§3.6), which is the safe way to
   summarise it.
9. **Sample sizes and replication are nowhere stated.** No number of replicate dialysis cells, no
   number of ligand concentrations per isotherm, no error bars on any figure that the text describes.
   "Within the experimental error, the intercept in Figure 1B ... is the same" is the only reference
   to error in the paper, and it is qualitative.
10. **What this paper does not contain**: any rate constant; any covalent or irreversible binding;
    any thiol, amine or lysine assay; any pH series (one pH, 8.0); any protein concentration series
    (one, 1 %); any temperature above 45 C except the 90 C denaturation pre-treatment; any measurement
    on a commercial isolate; any 11S/7S fractionation (that is the companion paper); any
    uncertainty on any constant; any supplementary material.
11. **What to request**: (i) the isotherm data behind Figures 1-4, which would allow the seven
    constants to be refitted with intervals and would test whether the shared intercepts really are
    shared; (ii) the number of replicates; (iii) the succinylation procedure and degree of
    modification; (iv) a measurement between 45 and 90 C, which is the range a cooking model actually
    needs and which nobody in this corpus has measured for soy.
12. **Registry gaps against `data/keys/compounds.yml`**: `nonanal` and `hexanal` are keyed;
    **`2_heptanone`, `2_octanone`, `2_nonanone` and `5_nonanone` are all absent**, and three of those
    four are used as `compound` keys in `REVERSIBLE_BINDING` (`kg_2_heptanone_soy`,
    `kg_2_octanone_soy`, `kg_2_nonanone_soy`, and again in the `_blg` rows). As with the aldehyde
    thresholds, **the matrix layer is keyed on compound names the registry does not know**. Note also
    that `2_heptanone` appears in `matrix_oav.py`'s `WATER_THRESHOLDS` with a back-solved 140 ug/L
    value, so the same unkeyed name is used in two independent tables.
