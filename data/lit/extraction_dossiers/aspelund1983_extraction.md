# Aspelund 1983 — EXTRACTION (spray-dried isoelectric soy protein isolate Edi-Pro A packed DRY into a GC column, 0.2019 m^2/g, inverse gas chromatography at 80, 90 and 100 C, five homologous series of alcohols, aldehydes, ketones, n-alkanes and methyl esters, n = 9 per point)

### THE FIRST SOY FLAVOUR-AFFINITY MEASUREMENT IN THE CORPUS ABOVE 60 C — AND IT IS DRY. Wave B26's record says "37 C is an in-mouth temperature, not a process one: nothing here licenses a pea binding constant at 90 or 140 C." This paper measures soy protein flavour affinity at exactly **80, 90 and 100 C**, and its answer is a measured **direction** — adsorption weakens monotonically as the protein is heated (Table III, p. 544) — plus an independent **chain-length slope of 578 cal per CH2 at 90 C (= 2.23x/CH2, mine)** against the registry's shipped `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81`, on the same 2-alkanone series Damodaran measured. But there is **no water anywhere in this experiment**: it is gas-solid chromatography on a dry powder, the authors themselves write that "a reordering of the binding affinities will occur" in aqueous solution (p. 544), and **no number here can be shipped as a `REVERSIBLE_BINDING` row**.

**Source on disk:** `data/articles/aspelund1983.pdf` (7 pp., J. Agric. Food Chem. 1983, 31, 539-545).
Read from the `pdftotext -layout` text layer, with **Tables I (p. 542) and II/III (p. 544) verified cell-by-cell against
the rendered page images** — the text layer scrambled Table I's column order and mis-OCR'd two cells of Table III
(2-hexanone at 100 C as "131" for the printed **731**; methyl pentanoate at 100 C as "713" for the printed **773**).
The re-typed tables below are from the page images and are correct. **Figures 1-9 are images and carry no printed
numbers**: Fig. 1 is four scanning electron micrographs, Figs. 2-4 are the ln(t_cor) vs 1/T regression lines for
2-pentanone, 1-heptanal and 2-octanone (three columns each), and Figs. 5-9 are the -ΔH-vs-carbon-number plots for
hydrocarbons, aldehydes, ketones, methyl esters and alcohols with 3-standard-error bars. Every value plotted in
Figs. 5-9 is also printed in Table II, so nothing is lost there; the **standard errors themselves are figure-only**.
Repo status before this dossier: Aspelund 1983 has **no extraction dossier**, and is **not cited** in
`src/kinetic_core/parameters_matrix.py`, `src/kinetic_core/matrix_sites.py` or `data/species/protein_matrices.yml`.

## 0. Identity

| field | value |
|---|---|
| Title | "Adsorption of Off-flavor Compounds onto Soy Protein: A Thermodynamic Study" |
| Authors | Thomas G. Aspelund* and Lester A. Wilson — Department of Food Technology, Iowa State University, Ames, Iowa 50011 |
| Venue | J. Agric. Food Chem. 1983, **31** (3), **539-545** |
| DOI | **NO DOI IS PRINTED IN THE PDF.** This is a 1983 ACS article; the footer carries only the CODEN-era line `0021-8561/83/1431-0539$01.50/0` and "© 1983 American Chemical Society". Nothing DOI-shaped appears anywhere in the seven pages. |
| Dates | Received for review 14 July 1982; revised manuscript received 30 December 1982; accepted 10 January 1983 |
| Funding / provenance | Iowa State University Research Foundation and the Iowa Agriculture and Home Economics Station (Journal Paper No. J-10721, Project 2164, contributing to North Central Regional Research Project NC-136). Portions presented at the World Soybean Research Conference-II, 26-29 March 1979, Raleigh NC (p. 545) |
| Protein | **Edi-Pro A**, Ralston Purina — a **spray-dried isoelectric soy protein isolate**, used as **dry GC column packing**. Surface area **0.2019 ± 0.0011 m^2/g** by Micromeritics, "no difference in surface area for samples degassed at 30 or 100 C" (p. 540) |
| The 17 significant ligands | n-nonane, n-decane; 2-hexanone, 2-heptanone, 2-octanone; 1-hexanal, 1-heptanal, 1-octanal; methyl pentanoate, methyl hexanoate, methyl heptanoate, methyl octanoate; 1-butanol, 1-pentanol, 1-hexanol, 1-heptanol, 1-octanol. CAS numbers for all 17 are printed in the Registry No. paragraph (p. 544) |
| The 7 NON-significant ligands | 2-butanone, 2-pentanone, 1-butanal, 1-pentanal, n-hexane, n-heptane, n-octane (p. 542) |
| Naming | "-ΔH" = heat of adsorption, kcal/mol, from the slope of ln(t_cor) vs 1/T times R; "V_S" = specific retained volume, mL/m^2; "-ΔG" = RT ln V_S, cal/mol; "-ΔS" = back-computed from ΔG = ΔH − TΔS, cal/(mol K) |
| Companions on disk | `crowther1980_extraction.md` (the same laboratory, the same Edi-Pro A, the processing/heat-treatment arm — Aspelund cites it as "Crowther, A.; Wilson, L. A.; Glatz, C. E. *J. Food Process Eng.* **1981**, *4*, 99"), `damodaran1981_extraction.md` (the aqueous-dialysis soy ketone series this paper compares itself to by name), `guo2019_extraction.md` and `Xu2022_extraction.md` (the preheat arm on soy and pea), `bi2022_extraction.md` (pea, 37 C) |

## 1. Why it matters

**The temperature limit is the whole point.** `REVERSIBLE_BINDING` in `src/kinetic_core/parameters_matrix.py`
carries 21 rows, and their measurement temperatures are 25 C (Damodaran soy dialysis), 30 C (Meynier skim milk,
Andriot beta-lactoglobulin), 37 C (the Wave B26 Bi 2022 pea rows) and 40 C (Leksrisompong caseinate). The B26
record states the limit in plain words: **"37 C is an in-mouth temperature, not a process one: nothing here
licenses a pea binding constant at 90 or 140 C."** This paper is the first source in the corpus that measures
flavour affinity to a food protein **at 80, 90 and 100 C** — squarely inside the process band the model actually
runs in. That is why it was read.

**And it is the wrong phase.** Every shipped `REVERSIBLE_BINDING` row is an aqueous constant: a compound
partitioning between air and a *protein solution*, stored as `K_g = (K_water/K_matrix − 1) / protein_g_per_L`
in L per gram of protein. Aspelund's experiment has **no aqueous phase at all**. Edi-Pro A is packed dry into a
3-ft glass GC column, nitrogen is the mobile phase, and what is measured is a **gas-solid** distribution — a
retention volume per unit *surface area* of dry powder. The two quantities are not the same object and the
authors say so themselves on p. 544: *"it is to be expected that, when these compounds are in an aqueous
solution, a reordering of the binding affinities will occur (that is, the alcohols may not be bound to the same
extent due to their interactions with water)."* **No row in this paper may be shipped into `REVERSIBLE_BINDING`,
and the 80-100 C temperature range does not by itself lift the B26 limit** — lifting that limit needs an
*aqueous* constant above 60 C, which this paper does not contain.

**What it does give, and it is not nothing.**

1. **A measured temperature dependence, with a sign, over 80-100 C.** Table III (p. 544) runs four
   six-carbon compounds — 1-hexanal, 2-hexanone, 1-hexanol, methyl pentanoate — at 80, 90 and 100 C, and −ΔG
   falls monotonically for all four (1-hexanal 989 → 825 → 540 cal/mol; 2-hexanone 1014 → 857 → 731;
   1-hexanol 2750 → 2470 → 2110; methyl pentanoate 1100 → 915 → 773). The paper's reading: *"higher
   temperatures are less favorable to adsorption, which is characteristic of a physical adsorption process."*
   The registry has **no measured temperature coefficient of a binding constant of any kind**, and this is one,
   on soy, over the process band. Its transfer to an aqueous constant is a separate question (Flags 3).
2. **An independent chain-length slope at 90 C.** The registry's `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81` is the
   geometric mean of Andriot's 2.72x/CH2 (beta-lactoglobulin, headspace, 30 C) and Damodaran's 2.9x/CH2 (soy,
   dialysis, 25 C), and the module comment calls it "the ONLY licence in this layer to move a constant from one
   chain length to another". This paper prints a third determination: **578 cal per methylene group** (p. 544),
   averaged over five homologous series on soy protein at 90 C — which is **2.23x per CH2 (mine)** at 363.15 K.
   Aspelund makes the comparison himself, to the same source the registry leans on: *"remarkably similar to the
   value of 600 cal/CH2 obtained by Damodaran and Kinsella (1981) using an aqueous equilibrium dialysis system."*
   And the raw V_S column of Table I gives the same answer independently: over all twelve consecutive
   chain-length steps in the five series, the geometric mean is **2.27x per CH2 (mine)**. So the shipped 2.81 is
   **21-26 % high** against a dry-phase soy determination at 90 C. That is a check, not a replacement (Flags 4).
3. **The one homologous series that overlaps the registry exactly.** 2-hexanone / 2-heptanone / 2-octanone at
   90 C here; 2-heptanone / 2-octanone / 2-nonanone at 25 C in Damodaran (`kg_2_heptanone_soy` 4.40e-3,
   `kg_2_octanone_soy` 1.24e-2, `kg_2_nonanone_soy` 3.72e-2 L/g). Same protein species, same compound class,
   two carbons of overlap, 65 C apart, two methods. The V_S step ratios here are **2.32x and 2.46x (mine)**
   against Damodaran's 2.82x and 3.00x (mine, from the shipped values).
4. **A functional-group ordering measured on soy, and it is NOT the aqueous one.** On dry soy at 90 C the
   ladder is **alcohols >> aldehydes > methyl esters ≈ ketones >> n-alkanes** (Table II, p. 544). The alcohols
   sit 5-7 kcal/mol above the carbonyls at equal carbon number (9-10 kcal/mol for 1-heptanol), which the paper
   attributes to a *second* hydrogen bond. This is the reverse of what an aqueous hydrophobic-partition picture
   would predict, and it is a direct, measured caution against any log-P-shaped matrix term — which
   `parameters_matrix.py` already refuses under k4b hold-out guard #4. **This paper is corroborating evidence
   for that refusal, on soy, and should be cited beside it.**
5. **The preheat question, by reference only.** p. 544 reports the companion result: *"Crowther et al. (1981)
   observed that the adsorption coefficient (K) decreased with heat treatment of Edi-Pro A. They postulated that
   as the protein denatured more nonpolar regions were exposed, decreasing both the solubility and availability
   of polar binding sites."* `matrix_sites.py` charges its binding sites once at the start of the cook and does
   not model their change with heating; that is the gap Guo 2019 and Xu 2022 speak to directly, and this is a
   third, older voice on it — **but it is a citation, not a measurement made here** (`[C]`).

**One thing this paper is NOT, and the registry has been burned on this before.** The −ΔH values in Tables I
and II are **heats of adsorption obtained from a van 't Hoff-type plot** — ln(corrected retention time) against
1/T, slope × R, after Gale & Beebe 1964 (p. 540). They are **THERMODYNAMIC (isosteric) enthalpies of an
adsorption equilibrium. They are NOT activation energies.** Reading 8.89 kcal/mol for 1-hexanal as an E_a and
feeding it to `matrix_sites.py`'s `ea_band_kj_mol` (which carries 15-20 kJ/mol for the aldehyde-amine channel,
from a genuine rate measurement) would be exactly the category error this repository has flagged before. The
sign convention makes it worse: these are printed as **−ΔH**, i.e. adsorption is **exothermic**, and an
activation energy cannot be negative. Nothing in this paper is a rate; there is no time axis anywhere in it.

**What this paper does NOT give the repository**: any aqueous measurement; any binding constant in M^-1 or L/g;
any pH (there is no solvent to have one); any protein concentration in g/L (the protein is the stationary
phase); any rate constant; any activation energy; any covalent-adduct evidence; any measurement below 80 C;
any measurement of a pyrazine, pyridine, furan, thiol or sulfur compound; any 2-alkenal or other
alpha,beta-unsaturated carbonyl (so **nothing here touches `ALPHA_BETA_UNSATURATION`**); any error bar as a
number (the 3-SE bars are figure-only).

## 2. Methods as they matter to a model

- **The pot, such as it is.** There is no pot. The "matrix" is a **3-ft silanized glass GC column, 2 mm i.d.,
  dry-packed with Edi-Pro A** using a Millipore vacuum pump and hand-tapping (a vibrator packed it too densely
  and cut the flow rate). **Exact packing mass: 1.4 g (column 1), 1.4 g (column 2), 1.5 g (column 3)**, glass
  wool plugs at both ends, conditioned overnight at 80 C under 20 mL/min nitrogen (p. 540). Three separate
  columns were built and the whole study repeated on each.
- **Loading, expressed the only way it can be.** The registry's `protein_g_per_L` has no analogue. What this
  paper gives instead is **surface area, 0.2019 ± 0.0011 m^2/g**, and a specific retained volume in **mL per
  m^2**. Multiplying gives a per-gram retention volume (section 4), which is a *gas-solid* distribution constant
  and not the registry's K_g.
- **Temperature and its programme.** **80 C**, then up to **90 C**; the next day repeated at **90 C**, then up
  to **100 C** (p. 540). So 90 C is measured twice, on different days, and the 80 → 90 and 90 → 100 legs come
  from different runs. Nitrogen flow "frequently monitored and kept constant (a prerequisite for the heat of
  adsorption determination) at **20 mL/min at all three temperatures**". Injector and detector held at 150 C.
- **Dosing.** Compounds stored in **5-mL vials** sealed with Teflon-coated septa and aluminium seals; **5 µL of
  the headspace above the neat compound at room temperature** drawn with a 50-µL gas-tight syringe and injected
  onto the soy column. **Three replications per compound per column temperature, on three columns → n = 9**
  (Table I footnote a, p. 542). Quantities per 5-µL injection, from vapour-pressure data (p. 540): **~10^-6 g**
  for the lower-molecular-weight ketones, aldehydes and hydrocarbons; **~10^-9 g** for the larger-molecular-weight
  alcohols and aldehydes; **10^-7-10^-8 g** for the methyl esters. The paper stresses that such minute quantities
  "was a prerequisite for this type of study" — i.e. this is the **infinite-dilution (Henry's-law) regime**, one
  molecule at a time on a bare surface, with no competition and no saturation. That is a different regime from
  every aqueous binding study in the corpus, all of which dose in the 0.05-2.5 mM range.
- **Measurement family: INVERSE GAS CHROMATOGRAPHY (gas-solid).** The protein is the stationary phase and the
  measured quantity is a retention time. This is a **sixth** method family beside the five the registry's
  `method` field already distinguishes (`headspace_depletion`, `equilibrium_dialysis`, `gel_filtration`,
  `static_headspace_partition`, `sensory_BET`). It is not a member of any of them, and given that k2 sec. B.3
  found a 35x gap between just two of the existing five on aldehydes, **it must not be pooled with any of them**
  (Flags 1).
- **Dead-time correction.** Every retention time corrected by subtracting the retention time of **methane**, a
  non-adsorbed gas (p. 540). This matters: the seven compounds that failed significance did so because their
  retention times "were very close to that of methane", so slight injection errors gave large standard errors
  (p. 542). The dynamic range of the method is bounded from below by the methane peak.
- **The heat of adsorption.** **ln T_cor = −ΔH/(RT) + C at constant flow rate**, Gale & Beebe 1964 (p. 540).
  Slope of ln(corrected average retention time) against 1/T, times R. **Three points only** (80, 90, 100 C =
  1/T × 10^3 of 2.83, 2.75, 2.68 — the x-axis of Figs. 2-4). A three-point van 't Hoff line over a 20 C span.
- **The free energy.** **V_S = (T_cor)(flow rate)/[(surface area)(mass of the packing material)]**, Sawyer &
  Brookman 1968, and **−ΔG = RT ln V_S** (p. 544). Note what this means: **ΔG is not independent of V_S; it IS
  V_S, on a log scale.** The −ΔG column of Table I carries no information the V_S column does not.
- **The entropy.** −ΔS is **not measured**: it is back-computed from the Gibbs-Helmholtz relation
  ΔG = ΔH − TΔS (Table I footnote b, p. 542), i.e. from the other two columns. It is arithmetic, not a third
  observation.
- **Statistics.** Linear regression on the average corrected retention time against 1/T, per column; the
  analysis returned (1) the significance level of the fit, (2) the estimated slopes and standard errors, (3)
  whether the slope exceeded zero significantly, (4) temperature effects, (5) within-column and within-replicate
  variation (p. 541). Table I is footnoted "Statistically significant compounds with a **99 % confidence limit**
  (n = 9)"; Table II is headed "Statistically Significant (**95 %**)"; Table III is footnoted "**99 % confidence
  level** ... Average from three columns and nine replications".
- **The one compound with a broken fit.** **2-octanone** was significant but showed real between-column
  variation (F-test); a mass-temperature interaction ruled out the packing-mass differences as the cause. **For
  2-octanone the average slope is declared INVALID** and the three individual column slopes must be used:
  **column 1 = 10.05, column 2 = 10.42, column 3 = 11.37 kcal/mol** (Table I footnote c, p. 542). The 10.61
  printed in Table II is the average of those three and is flagged as such.
- **The blank.** Qualitative headspace analysis of Edi-Pro A itself in 100-mL vials (14.5-18.0 g of isolate) at
  21, 45 and 80 C, sampled at 1 h, 24 h and 5 days on Carbowax 20M and SE-30 columns. At 21 and 45 C: virtually
  nothing. **At 80-100 C: 7-13 compounds appeared, tentatively n-hexane and n-pentane** — i.e. the protein
  itself emits volatiles at the study's own temperatures. Conditioning the packed column 6-8 h at 80 C removed
  them below the FID's most sensitive setting of 1 × 10^-12 A/s (p. 540). **The blank was cleared, but only
  after the protein had been held at 80 C for 6-8 h, which is itself a heat treatment** (Flags 6).
- **Structural control.** Scanning electron micrographs of Edi-Pro A before and after the adsorption studies
  (Fig. 1A-D, 104× and 1040×): *"no noticeable heating effects on the protein's structure were observed"* (p. 540).
  This is a **morphological** observation at 104-1040× magnification. It is not a denaturation assay, not a DSC
  trace, not an SDS-PAGE, and it says nothing about tertiary structure (Flags 6).

## 3. Tables re-typed

Every cell below was read off the rendered page image, not the text layer. Evidence marks: `[M]` measured in
this study, `[C]` cited from elsewhere, `[F]` fitted/extrapolated by the authors.

### Table I (p. 542). "Thermodynamic Quantities Determined for Soy Isolate at 90 °C"

| compounds^a | V_S, mL/m^2 | −ΔH, kcal/mol | −ΔG, cal/mol | −ΔS,^b cal/(mol K) |
|---|---:|---:|---:|---:|
| *n*-nonane | 3.75 `[M]` | 6.52 `[M]` | 952 `[M]` | 15.32 `[F]` |
| *n*-decane | 8.39 `[M]` | 8.36 `[M]` | 1530 `[M]` | 18.79 `[F]` |
| 2-hexanone | 3.29 `[M]` | 6.04 `[M]` | 857 `[M]` | 14.26 `[F]` |
| 2-heptanone | 7.62 `[M]` | 8.11 `[M]` | 1470 `[M]` | 18.29 `[F]` |
| 2-octanone^c | 18.73 `[M]` | 10.61 `[M]` | 2110 `[M]` | 23.38 `[F]` |
| 1-hexanal | 3.14 `[M]` | 8.89 `[M]` | 825 `[M]` | 22.18 `[F]` |
| 1-heptanal | 6.17 `[M]` | 9.61 `[M]` | 1310 `[M]` | 22.81 `[F]` |
| 1-octanal | 12.54 `[M]` | 13.52 `[M]` | 1820 `[M]` | 32.44 `[F]` |
| methyl pentanoate | 3.58 `[M]` | 6.71 `[M]` | 915 `[M]` | 15.95 `[F]` |
| methyl hexanoate | 8.16 `[M]` | 8.19 `[M]` | 1510 `[M]` | 18.35 `[F]` |
| methyl heptanoate | 19.81 `[M]` | 10.48 `[M]` | 2160 `[M]` | 22.91 `[F]` |
| methyl octanoate | 46.52 `[M]` | 12.56 `[M]` | 2760 `[M]` | 26.93 `[F]` |
| 1-butanol | 5.77 `[M]` | 10.45 `[M]` | 1260 `[M]` | 25.29 `[F]` |
| 1-pentanol | 10.99 `[M]` | 11.25 `[M]` | 1730 `[M]` | 26.18 `[F]` |
| 1-hexanol | 31.24 `[M]` | 13.89 `[M]` | 2470 `[M]` | 31.40 `[F]` |
| 1-heptanol | 70.95 `[M]` | 18.06 `[M]` | 3070 `[M]` | 41.22 `[F]` |
| 1-octanol | 165.52 `[M]` | 16.39 `[M]` | 3690 `[M]` | 24.94 `[F]` |

Footnotes exactly as printed: *^a Statistically significant compounds with a 99 % confidence limit (n = 9).*
*^b ΔS values were calculated from the Gibbs-Helmholtz equation (ΔG = ΔH − TΔS).* *^c For 2-octanone, the
average slope cannot be used for ΔH determinations. The individual slopes to be used are as follows:
column 1 = 10.05 kcal/mol; column 2 = 10.42 kcal/mol, column 3 = 11.37 kcal/mol.*

**Note the columns are not independent.** −ΔG = RT ln V_S by construction (p. 544), and −ΔS = (ΔH − ΔG)/T by
construction (footnote b). **V_S and −ΔH are the only two measured quantities in this table.**

**Internal check (mine): −ΔG against RT ln V_S at 363.15 K.** With R = 1.987 cal/(mol K), RT = 721.6 cal/mol.
1-hexanal: 721.6 × ln(3.14) = 825.4 against a printed **825** ✓. 1-octanol: 721.6 × ln(165.52) = 3684 against
**3690** ✓. 2-heptanone: 721.6 × ln(7.62) = 1466 against **1470** ✓. methyl octanoate: 721.6 × ln(46.52) = 2769
against **2760** ✓. The relation closes to better than 0.4 % on every row tested. **The −ΔG column carries no
information beyond V_S.**

**Internal check (mine): the alcohol series is NOT monotone in −ΔH.** 1-heptanol is printed at **18.06** and
1-octanol at **16.39** kcal/mol — the eight-carbon alcohol adsorbs *less* strongly than the seven-carbon one,
by 1.67 kcal/mol, while its V_S is 2.33x *larger* (70.95 → 165.52). **ΔH and V_S disagree on the direction of
the last step of the alcohol series** (Flags 5). Every other series is monotone in both.

### Table II (p. 544). "Statistically Significant (95%) and Extrapolated Heats of Adsorption on Edi-Pro A"

Values are −ΔH in kcal/mol; the column header is "carbon no.".

| compounds | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ketones | 2.3^a `[F]` | 4.3^a `[F]` | 6.04 `[M]` | 8.11 `[M]` | 10.61^b `[M]` | 12.1^a `[F]` | 14.1^a `[F]` |
| methyl esters | 2.1^a `[F]` | 4.6^a `[F]` | 6.71 `[M]` | 8.19 `[M]` | 10.48 `[M]` | 12.56 `[M]` | 14.7^a `[F]` |
| aldehydes | 4.7^a `[F]` | 6.6^a `[F]` | 8.89 `[M]` | 9.61 `[M]` | 13.52 `[M]` | 14.3^a `[F]` | 16.3^a `[F]` |
| alcohols | 10.45 `[M]` | 11.25 `[M]` | 13.89 `[M]` | 18.06 `[M]` | 16.39 `[M]` | 19.6^a `[F]` | 21.6^a `[F]` |
| hydrocarbons | 0^a `[F]` | 0^a `[F]` | 0.4^a `[F]` | 2.4^a `[F]` | 4.3^a `[F]` | 6.52 `[M]` | 8.36 `[M]` |

Footnotes exactly as printed: *^a Extrapolated values.* *^b Average of three slopes for 2-octanone.*

**Which cells are real.** Only 17 of the 35 cells are measured; **18 are extrapolated off the carbon-number
plots of Figs. 5-9 and are marked `[F]`**. In particular *every* hydrocarbon cell from C4 to C8 is extrapolated,
two of them to a value of **exactly zero** — which is a fitted intercept artefact, not a measurement that
n-butane does not adsorb. The paper is explicit about why the extrapolations exist: the seven non-significant
compounds eluted too close to methane, so "extrapolated values would be the most accurate ΔH values for the
nonsignificant compounds" (p. 542). **Do not carry any `[F]` cell as an observation.**

**Note on carbon numbering for the esters.** Table II indexes methyl pentanoate under carbon no. **6**, methyl
hexanoate under **7**, methyl heptanoate under **8** and methyl octanoate under **9** — i.e. the count includes
the methyl of the ester group, not just the acyl chain. Confirmed against Table I (methyl pentanoate −ΔH 6.71 =
the C6 ester cell).

### Table III (p. 544). "Temperature Dependence of Gibb's Free Energy and Entropy of Binding for Six Carbon Flavor Compounds on Soy Protein Isolate^a"

| compound | temp, °C | −ΔG, cal/mol^b | −ΔS, cal/(mol K)^b |
|---|---:|---:|---:|
| 1-hexanal | 80 | 989 `[M]` | 22.36 `[F]` |
| | 90 | 825 `[M]` | 22.18 `[F]` |
| | 100 | 540 `[M]` | 22.38 `[F]` |
| 2-hexanone | 80 | 1014 `[M]` | 14.23 `[F]` |
| | 90 | 857 `[M]` | 14.26 `[F]` |
| | 100 | **731** `[M]` | 14.23 `[F]` |
| 1-hexanol | 80 | 2750 `[M]` | 31.54 `[F]` |
| | 90 | 2470 `[M]` | 31.40 `[F]` |
| | 100 | 2110 `[M]` | 31.56 `[F]` |
| methyl pentanoate | 80 | 1100 `[M]` | 15.91 `[F]` |
| | 90 | 915 `[M]` | 15.95 `[F]` |
| | 100 | **773** `[M]` | 15.92 `[F]` |

Footnotes exactly as printed: *^a 99 % confidence level.* *^b Average from three columns and nine replications.*

The two bolded cells are the ones the text layer mis-OCR'd (as 131 and 713); the page image reads **731** and
**773**, and both are confirmed by the arithmetic below.

**Internal check (mine): −ΔG = −ΔH + T(−ΔS), with −ΔH taken from Table I.** 2-hexanone (−ΔH 6.04 kcal/mol,
−ΔS 14.23): at 353.15 K, −6040 + 353.15 × 14.23 = **−1015** against a printed 1014 ✓; at 373.15 K,
−6040 + 373.15 × 14.23 = **−730** against a printed 731 ✓. 1-hexanal (8.89, 22.38 at 100 C):
−8890 + 373.15 × 22.38 = **−539** against 540 ✓. 1-hexanol (13.89, 31.56 at 100 C):
−13890 + 373.15 × 31.56 = **−2113** against 2110 ✓. methyl pentanoate (6.71, 15.92 at 100 C):
−6710 + 373.15 × 15.92 = **−769** against 773 ✓ (0.5 %). **The table closes.** This also confirms that −ΔS is
pure arithmetic: it is whatever makes the identity hold, which is why it is near-constant down each block by
construction rather than by observation.

### Numbers printed in the running text

| quantity | value | where | class |
|---|---|---|---|
| surface area of Edi-Pro A | **0.2019 ± 0.0011 m^2/g**, "no difference ... degassed at 30 or 100 C" | p. 540 | `[M]` |
| packing mass per column | **1.4 / 1.4 / 1.5 g** (columns 1 / 2 / 3) | p. 540 | `[M]` |
| nitrogen flow | **20 mL/min** at all three temperatures | p. 540 | `[M]` |
| injector / detector temperature | **150 C** | p. 540 | `[M]` |
| detector gas flows | H2 30 mL/min, air 300 mL/min | p. 540 | `[M]` |
| injection quantity, low-MW ketones/aldehydes/hydrocarbons | **~10^-6 g** per 5-µL injection | p. 540 | `[M]` |
| injection quantity, high-MW alcohols and aldehydes | **~10^-9 g** | p. 540 | `[M]` |
| injection quantity, methyl esters | **10^-7-10^-8 g** | p. 540 | `[M]` |
| FID sensitivity floor | 1 × 10^-12 A/s | p. 540 | `[M]` |
| protein's own headspace, 21 and 45 C | "virtually no compounds present" | p. 540 | `[M]` |
| protein's own headspace, 80-100 C | **7-13 compounds**, tentatively n-hexane and n-pentane; removed by 6-8 h conditioning at 80 C | p. 540 | `[M]` |
| **mean free-energy increment per methylene group** | **578 cal/CH2** | p. 544 | `[M]` — the headline number |
| the comparison the author draws | "remarkably similar to the value of **600 cal/CH2** obtained by Damodaran and Kinsella (1981) using an aqueous equilibrium dialysis system" | p. 544 | `[C]` |
| group B (ketone / aldehyde / ester) advantage over hydrocarbons at equal carbon number | **6-8 kcal/mol** higher −ΔH | p. 543 | `[M]` |
| alcohol advantage over group B at equal carbon number | **5-7 kcal/mol** (9-10 kcal/mol for 1-heptanol) | p. 544 | `[M]` |
| hydrogen-bond energy of formation, for interpretation | 3-10 kcal/mol, generally 5 | p. 543 | `[C]` (MacKenzie 1962) |
| −ΔG range, alcohols | **−4 to −1 kcal/mol** | p. 544 | `[M]` |
| −ΔG range, all other compounds | **−2.5 to −0.5 kcal/mol** | p. 544 | `[M]` |
| −ΔS range, alkanes/esters/ketones/aldehydes | **−30 to −15 cal/(mol K)** | p. 544 | `[F]` (from the computed ΔS column) |
| −ΔS range, alcohols | **−40 to −25 cal/(mol K)** | p. 544 | `[F]` |
| 2-octanone individual column slopes | **10.05 / 10.42 / 11.37 kcal/mol** | Table I footnote c, p. 542 | `[M]` |
| **the Crowther heat-treatment result** | "the adsorption coefficient (K) **decreased with heat treatment** of Edi-Pro A ... as the protein denatured more nonpolar regions were exposed, decreasing both the solubility and availability of polar binding sites" | p. 544 | `[C]` — Crowther et al. 1981, not measured here |
| the authors' own aqueous caveat | "when these compounds are in an aqueous solution, a **reordering of the binding affinities will occur** (that is, the alcohols may not be bound to the same extent due to their interactions with water)" | p. 544 | `[C]`/interpretation |
| electron-microscope verdict on heating | "no noticeable heating effects on the protein's structure were observed" (104× and 1040×) | p. 540 | `[M]` (morphological only) |

**Figure-only quantities.** The standard errors of every −ΔH (plotted as 3-SE bars in Figs. 5-9) are **not
printed as numbers anywhere**. The per-column regression slopes are figure-only except for 2-octanone's three,
which are printed in Table I's footnote. The regression lines of Figs. 2-4 carry no printed intercepts. Per
house rule none of these is typed as a number here.

### Arithmetic on the printed constants (all mine)

**1. The 578 cal/CH2 as a multiplicative slope (mine).** The registry stores a chain-length effect as a *ratio*
per CH2, not a free-energy increment. Converting at the temperature the number was measured, T = 363.15 K,
R = 1.987 cal/(mol K): **exp(578 / (1.987 × 363.15)) = 2.228x per CH2 (mine)**. For reference, converting
Damodaran's cited 600 cal/CH2 at his 25 C gives **exp(600/(1.987 × 298.15)) = 2.753x (mine)**, which is within
5 % of the 2.9x/CH2 the registry actually carries for Damodaran — so the two determinations are being compared
on a defensible common basis. **Note the temperature trap**: 578 cal/CH2 evaluated at 298 K instead would give
2.65x, nearly the same as Damodaran's. **The gap between 2.23 and 2.75 is almost entirely the 65 C temperature
difference, not a difference in the underlying free-energy increment (578 vs 600 cal/CH2 is a 3.7 % gap).**
That is itself the finding: a chain-length *ratio* is temperature-dependent even when the free-energy increment
behind it is not, and the registry stores the ratio.

**2. The same slope from the raw V_S column (mine).** Consecutive-step ratios of V_S within each series at 90 C:

| series | consecutive V_S ratios (mine) | geometric mean per CH2 (mine) |
|---|---|---:|
| ketones (2-hexanone → 2-heptanone → 2-octanone) | 2.316, 2.458 | **2.386** |
| aldehydes (1-hexanal → 1-heptanal → 1-octanal) | 1.965, 2.032 | **1.998** |
| methyl esters (pentanoate → hexanoate → heptanoate → octanoate) | 2.279, 2.428, 2.348 | **2.351** |
| alcohols (1-butanol → 1-pentanol → 1-hexanol → 1-heptanol → 1-octanol) | 1.905, 2.843, 2.271, 2.333 | **2.314** |
| hydrocarbons (n-nonane → n-decane) | 2.237 | **2.237** |
| **all five series, 12 steps pooled** | | **2.272** |

The pooled 2.272x agrees with the 2.228x recovered from the printed 578 cal/CH2 to **2 %**, which is the
expected consistency (the 578 is that pooled average, expressed as an energy). **Against the shipped
`CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81`, this is 19-21 % low.**

**3. The direct overlap with the shipped Damodaran soy ketones (mine).** Damodaran's shipped per-gram constants
are 4.40e-3 (2-heptanone), 1.24e-2 (2-octanone), 3.72e-2 (2-nonanone) L/g, giving step ratios **2.818x and
3.000x (mine)**. Aspelund's dry-phase V_S on the overlapping window gives **2.316x (2-hexanone → 2-heptanone)
and 2.458x (2-heptanone → 2-octanone)**. Same protein species, same compound class, **the ketone chain-length
slope is 18-22 % shallower on dry soy at 90 C than on aqueous soy at 25 C.** This is one measured contrast, not
a temperature coefficient: phase and temperature move together between the two studies and cannot be separated.

**4. V_S put on a per-gram basis (mine), and why it is still NOT K_g.** V_S is printed per m^2; the surface area
is printed as 0.2019 m^2/g; so V_S × 0.2019 is a retention volume **per gram of protein**:

| compound | V_S, mL/m^2 | per-gram retention volume (mine) |
|---|---:|---:|
| 2-hexanone | 3.29 | 6.64e-4 L/g |
| 2-heptanone | 7.62 | **1.54e-3 L/g** |
| 2-octanone | 18.73 | **3.78e-3 L/g** |
| 1-hexanal | 3.14 | 6.34e-4 L/g |
| 1-heptanal | 6.17 | 1.25e-3 L/g |
| 1-octanal | 12.54 | 2.53e-3 L/g |
| methyl pentanoate | 3.58 | 7.23e-4 L/g |
| methyl hexanoate | 8.16 | 1.65e-3 L/g |
| methyl heptanoate | 19.81 | 4.00e-3 L/g |
| methyl octanoate | 46.52 | 9.39e-3 L/g |
| 1-butanol | 5.77 | 1.17e-3 L/g |
| 1-pentanol | 10.99 | 2.22e-3 L/g |
| 1-hexanol | 31.24 | 6.31e-3 L/g |
| 1-heptanol | 70.95 | 1.43e-2 L/g |
| 1-octanol | 165.52 | 3.34e-2 L/g |
| n-nonane | 3.75 | 7.57e-4 L/g |
| n-decane | 8.39 | 1.69e-3 L/g |

**READ THE NEXT SENTENCE BEFORE USING ANY OF THESE.** These numbers are in L/g and the registry's K_g is in
L/g, and **they are not the same quantity**. K_g is defined as (K_water/K_matrix − 1) / protein_g_per_L: it is
the *excess* affinity of a protein solution over plain water, per gram, and it goes to zero when the protein
does nothing. A gas-solid retention volume has **no water leg to divide out** and cannot go to zero. The
apparent agreement between 2-heptanone here (1.54e-3 L/g) and Damodaran's `kg_2_heptanone_soy` (4.40e-3 L/g) —
a factor of 2.9 — is **numerological**, and it is recorded here only so that a later reader does not rediscover
it and mistake it for cross-validation. **Do not ship any of these seventeen numbers into `REVERSIBLE_BINDING`.**

**5. The temperature coefficient of −ΔG, as the registry would need it (mine).** From Table III, the drop in
−ΔG per 10 C, averaged over the two legs: 1-hexanal **−225 cal/mol per 10 C**; 2-hexanone **−142**;
1-hexanol **−320**; methyl pentanoate **−164**. Converted to a multiplicative fall in V_S per 10 C at these
temperatures (V_S = exp(−ΔG/RT), so the ratio also picks up the RT change): 1-hexanal
**3.14 → 2.07, i.e. 1.52x lower at 100 C than at 90 C (mine)**; 2-hexanone 3.29 → 2.68 (**1.23x**);
1-hexanol 31.24 → 17.22 (**1.81x**); methyl pentanoate 3.58 → 2.84 (**1.26x**).
**So adsorption on dry soy falls by roughly 1.2-1.8x per 10 C over 90-100 C, compound-dependent.** This is the
only measured temperature coefficient of a protein-flavour affinity in the corpus, and its transfer to an
aqueous constant is unlicensed (Flags 3).

**6. What the group contrasts are worth as ratios (mine).** At six carbons and 90 C, on V_S:
1-hexanol / 2-hexanone = **9.5x**; 1-hexanol / 1-hexanal = **9.9x**; 1-hexanal / 2-hexanone = **0.95x**
(indistinguishable); methyl pentanoate / 2-hexanone = **1.09x**. At nine carbons, methyl octanoate / n-nonane =
**12.4x**. **The alcohol/carbonyl contrast is ~10x on dry protein.** The corresponding aqueous contrast is not
measured here and, per the authors' own caveat, is expected to be different in sign as well as size.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Of the 17 significant ligands, the ones already keyed and
relevant: `hexanal` (this paper's "1-hexanal"), `2_heptanone`, `2_octanone`, `nonanal` is *not* studied here
(the aldehyde series stops at 1-octanal), `1_hexanol`. **`2_hexanone`, `1_octanal`, `1_heptanal`, `1_butanol`,
`1_pentanol`, `1_heptanol`, `1_octanol`, `n_nonane`, `n_decane` and the four methyl esters are not carried**;
none of them is a panel target and none needs adding on this paper's account, since no row here is shippable.
`COMPOUND_STRUCTURE` in `parameters_matrix.py` carries `hexanal`, `2_heptanone`, `2_octanone`, `2_nonanone`,
`nonanal` and `t_2_octenal`; it gained its first alcohol (`z_2_penten_1_ol`) only in Wave B26 and has no
`n_alkane` or `methyl_ester` class at all.

Every row below shares: **Edi-Pro A spray-dried isoelectric soy protein isolate, DRY, packed as a GC stationary
phase, 0.2019 m^2/g, nitrogen carrier at 20 mL/min, infinite-dilution injections of 10^-9 to 10^-6 g,
n = 9 (three columns × three replicates), inverse gas chromatography.** There is **no solvent, no pH and no
protein loading in g/L** on any of them.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| specific retained volume V_S, 2-heptanone | **7.62** | mL/m^2 | dry soy isolate, 90 C, inverse GC | Table I, p. 542 | `binding_constant` (gas-solid; NOT an aqueous K) |
| specific retained volume V_S, 2-octanone | **18.73** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, 2-hexanone | **3.29** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, 1-hexanal | **3.14** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, 1-heptanal / 1-octanal | **6.17 / 12.54** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, 1-butanol … 1-octanol | **5.77 / 10.99 / 31.24 / 70.95 / 165.52** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, methyl C5…C8 esters | **3.58 / 8.16 / 19.81 / 46.52** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| specific retained volume V_S, n-nonane / n-decane | **3.75 / 8.39** | mL/m^2 | as above | Table I, p. 542 | `binding_constant` (gas-solid) |
| **heat of adsorption −ΔH, 1-hexanal** | **8.89** | kcal/mol | dry soy, van 't Hoff slope over 80-100 C | Table I, p. 542; Table II, p. 544 | `binding_constant` — **THERMODYNAMIC, NOT AN ACTIVATION ENERGY** (Flags 2) |
| heat of adsorption −ΔH, 2-heptanone / 2-octanone | **8.11 / 10.61** | kcal/mol | as above (2-octanone's average slope declared invalid: use 10.05 / 10.42 / 11.37 per column) | Table I, p. 542 | `binding_constant` — thermodynamic |
| heat of adsorption −ΔH, 1-hexanol | **13.89** | kcal/mol | as above | Table I, p. 542 | `binding_constant` — thermodynamic |
| heat of adsorption −ΔH, all 17 significant compounds | see Table I | kcal/mol | as above | Table I, p. 542 | `binding_constant` — thermodynamic |
| heats of adsorption, the 18 extrapolated cells of Table II | 0 to 21.6 | kcal/mol | read off the carbon-number regression of Figs. 5-9 | Table II, p. 544, footnote a | `derived_assumption` — **the authors' extrapolation, not a measurement; two cells are exactly zero** |
| **−ΔG at 80 / 90 / 100 C, 1-hexanal** | **989 / 825 / 540** | cal/mol | dry soy, three temperatures | Table III, p. 544 | `within_study_ratio` — the corpus's only measured temperature series of a protein-flavour affinity |
| −ΔG at 80 / 90 / 100 C, 2-hexanone | **1014 / 857 / 731** | cal/mol | as above | Table III, p. 544 | `within_study_ratio` |
| −ΔG at 80 / 90 / 100 C, 1-hexanol | **2750 / 2470 / 2110** | cal/mol | as above | Table III, p. 544 | `within_study_ratio` |
| −ΔG at 80 / 90 / 100 C, methyl pentanoate | **1100 / 915 / 773** | cal/mol | as above | Table III, p. 544 | `within_study_ratio` |
| **free-energy increment per methylene group** | **578** | cal/CH2 | dry soy, 90 C, pooled over five homologous series | p. 544 | **`within_study_ratio`** — the FIT-candidate check on `CHAIN_LENGTH_SLOPE_PER_CH2` |
| the same, as a multiplicative slope | **2.23** | × per CH2 | at 363.15 K | exp(578/RT) (mine) | `derived_assumption` (arithmetic only) |
| the same, recovered from V_S directly | **2.27** | × per CH2 | 12 pooled steps, 90 C | Table I ratios (mine) | `derived_assumption` (arithmetic only) |
| ketone chain-length slope, dry soy 90 C | **2.32 and 2.46** | × per CH2 | 2-hexanone → 2-heptanone → 2-octanone | Table I (mine) | `within_study_ratio` |
| Damodaran's cited increment, for comparison | **600** | cal/CH2 | aqueous equilibrium dialysis, soy | p. 544 | `level_only` — **`[C]`, cited, not measured here** |
| alcohol advantage over ketones/aldehydes/esters at equal carbon number | **5-7** (9-10 for 1-heptanol) | kcal/mol on −ΔH | dry soy, 90 C | p. 544 | `within_study_ratio` |
| group-B (carbonyl/ester) advantage over n-alkanes at equal carbon number | **6-8** | kcal/mol on −ΔH | dry soy, 90 C | p. 543 | `within_study_ratio` |
| alcohol / ketone contrast at C6 on V_S | **9.5** | × | dry soy, 90 C | 31.24/3.29 (mine) | `within_study_ratio` |
| fall in adsorption per 10 C, 90 → 100 C | **1.23x (2-hexanone) to 1.81x (1-hexanol)** | × | dry soy | from Table III via V_S = exp(−ΔG/RT) (mine) | `derived_assumption` — **not transferable to an aqueous constant** (Flags 3) |
| surface area of the isolate | **0.2019 ± 0.0011** | m^2/g | Micromeritics; unchanged degassing at 30 vs 100 C | p. 540 | `measured_bound` (a physical property of the powder) |
| per-gram retention volumes | see section 3 item 4 | L/g | dry soy, 90 C | V_S × 0.2019 (mine) | **`derived_assumption` — DO NOT SHIP.** Same units as K_g, different quantity |
| entropy of adsorption −ΔS, all rows | 14.23 to 41.22 | cal/(mol K) | dry soy | Tables I and III | `derived_assumption` — **back-computed from ΔH and ΔG, not an independent observation** |
| adsorption coefficient K falls with heat treatment of Edi-Pro A | direction only, no value | — | — | p. 544, citing Crowther et al. 1981 | `level_only` — **`[C]`, a citation; go to `crowther1980_extraction.md` for the number** |
| the protein's own volatile emission at 80-100 C | **7-13 compounds**, tentatively n-hexane and n-pentane | count | Edi-Pro A held at 80-100 C in a sealed vial | p. 540 | `level_only` |

### Can these be put on the same basis as the shipped binding constants, i.e. converted to K_g in L/g?

**No, and the obstruction is structural rather than a missing number.** The registry's form is
`K_g = (K_water/K_matrix − 1) / protein_g_per_L`. It needs three inputs: an air/water partition coefficient, an
air/matrix partition coefficient measured on the **same** system, and a protein loading in **g/L of solution**.

- **The water leg does not exist.** There is no aqueous control in this paper, and there cannot be: the protein
  is the chromatographic stationary phase and the mobile phase is nitrogen. A gas-solid retention volume has
  nothing to be normalised against, so the "− 1" that makes K_g an *excess* over water cannot be formed. This
  is not a reporting omission the authors could fix; it is what inverse GC is.
- **The protein loading in g/L does not exist either**, and its analogue is not interchangeable. What the paper
  prints instead — 1.4-1.5 g of protein per column, and 0.2019 m^2/g — is a **mass and a surface area**, not a
  concentration. The correct normaliser for a gas-solid constant is surface area, which is why V_S is printed
  per m^2, and section 3 item 4 shows what happens if one converts anyway: numbers that look like K_g, are in
  the same units as K_g, and mean something else.
- **A molar mass is not needed** (unlike Bi 2022's Klotz K), because V_S is already an extensive-per-gram
  quantity. That obstruction, at least, is absent.
- **What the paper *can* be put on a common basis for is RATIOS.** The chain-length slope (578 cal/CH2 →
  2.23x/CH2), the functional-group contrasts, and the 80/90/100 C temperature series are all within-study
  ratios in which the gas-solid scale cancels, exactly as Amendment 4 lets Meynier's and Leksrisompong's
  suspect absolute static-headspace scales cancel in a ratio. **Those are the transferable objects here, and
  the chain-length slope is the one with a live consumer in the code.**

**On lifting the B26 temperature limit: this paper does not lift it.** It shows that soy protein still adsorbs
flavour compounds at 100 C, that the adsorption weakens as temperature rises, and roughly by how much on a dry
surface. It does **not** show what an aqueous soy or pea binding constant is at 90 C, because it never puts the
protein in water. The honest statement to carry forward is: *the direction of the temperature effect above 60 C
is measured (adsorption weakens), the magnitude is measured only in the gas-solid phase, and the aqueous
magnitude remains unmeasured.*

**Nothing here goes to `matrix_sites.py`.** That module wants second-order rate constants in M^-1 s^-1 and
activation energies in kJ/mol, and this paper contains neither. See Flags 2 for why its enthalpies must not be
put in the `ea_band_kj_mol` field.

## 5. Flags

1. **A sixth method family, and the corpus already knows what method mixing costs.** Inverse gas chromatography
   on a dry protein is not `headspace_depletion`, not `equilibrium_dialysis`, not `gel_filtration`, not
   `static_headspace_partition` and not `sensory_BET`. k2 sec. B.3 found a **35x** gap between dialysis and
   headspace on aldehydes alone, and Bi 2022 showed a same-lab, same-day **5-56x** gap between headspace and
   fluorescence. The gas-solid/aqueous gap has never been measured in this corpus and there is no reason to
   expect it to be smaller. **If any Aspelund number is ever carried, it needs its own `method` value, e.g.
   `inverse_gc_dry`, and `binding_constant_for` must refuse to cross into it.**
2. **The enthalpies are van 't Hoff quantities, not activation energies — say it out loud.** −ΔH comes from
   ln T_cor = −ΔH/(RT) + C (p. 540, after Gale & Beebe 1964): a three-point van 't Hoff plot of an equilibrium
   retention, over a 20 C span. It is the **isosteric heat of adsorption** of a reversible physical adsorption,
   and the paper's own reading of the negative sign is that the process "is characteristic of a physical
   adsorption process" (p. 544). It is **not** an E_a, it does not belong anywhere near
   `matrix_sites.py`'s `ea_band_kj_mol` (15-20 kJ/mol for the aldehyde-amine channel, from a rate measurement),
   and note the magnitudes would invite exactly that error: 8.89 kcal/mol = 37.2 kJ/mol (mine), which would
   look like a plausible E_a and is not one. There is **no rate and no time axis anywhere in this paper.**
3. **Dry is not wet, and the authors say so.** p. 544: *"it is to be expected that, when these compounds are in
   an aqueous solution, a reordering of the binding affinities will occur (that is, the alcohols may not be
   bound to the same extent due to their interactions with water)."* The mechanism they propose for the
   alcohols' 10x advantage is a **second hydrogen bond to the protein surface** — a bond that in water would be
   competing against bulk solvent for both the alcohol's proton and its oxygen. **The functional-group ordering
   measured here is a dry-phase ordering and must not be transferred.** The chain-length slope is on safer
   ground (a CH2 increment is a hydrophobic-surface effect in both phases) but is still a dry-phase number.
4. **The chain-length disagreement with the shipped 2.81 is mostly temperature, and that is the interesting
   part.** 578 cal/CH2 here against Damodaran's cited 600 is a **3.7 % gap in the free energy**; the same two
   numbers become **2.23x and 2.75x per CH2** once each is evaluated at its own temperature, a **19 % gap in the
   ratio**. The registry stores the ratio and applies it without a temperature argument. **`CHAIN_LENGTH_SLOPE_PER_CH2`
   is being used at process temperatures with a value measured at 25-30 C, and this paper says that costs
   roughly 20 % per application at 90 C** — in the direction of over-crediting long chains. That is a real,
   sizeable, previously-unmeasured bias, and it is the single most actionable thing in this dossier. It is also
   a dry-phase measurement of it (Flags 3), so it is a **caution to record, not a value to swap in**.
5. **The alcohol series breaks monotonicity in ΔH and the paper does not remark on it.** 1-heptanol's
   −ΔH is 18.06 kcal/mol and 1-octanol's is 16.39 — the longer chain adsorbs *less* strongly — while V_S moves
   the other way (70.95 → 165.52) and −ΔG rises as expected (3070 → 3690). Because −ΔS is back-computed, the
   inconsistency is absorbed into the entropy column, where 1-heptanol's 41.22 cal/(mol K) is the largest value
   in the table and 1-octanol's 24.94 is smaller than 1-hexanol's 31.40. **Either the 1-heptanol slope or the
   1-octanol slope is wrong, and there is no way to tell which from the printed data.** Fig. 9 plots both with
   3-SE bars but the SEs are not printed. Treat the whole alcohol ΔH series as the least reliable in the paper.
6. **The protein was heat-treated before every measurement, and the "no structural change" evidence is weak.**
   The column was conditioned **overnight at 80 C**, then a further **6-8 h at 80 C** to clear the isolate's own
   volatiles (p. 540), before any retention time was taken; and the 90 C and 100 C data were taken after that
   protein had already been held at 80-90 C for a day or more. **This is a preheated soy protein throughout**,
   which is directly relevant to the `matrix_sites.py` question of whether binding sites change with heating —
   and this paper cannot answer it, because it has no unheated arm. The evidence offered against a heating
   effect is a scanning electron micrograph before and after at 104× and 1040× magnification (Fig. 1): a
   **morphological** observation at particle scale. It rules out gross sintering or fusion of the spray-dried
   spheres. It says nothing whatever about unfolding, about surface hydrophobicity, or about the availability of
   polar binding sites — which is precisely what Crowther, in the same laboratory, reports *does* change
   (p. 544). **The paper's own citation contradicts the reassurance its micrograph offers.**
7. **The blank was not clean at the study's own temperatures.** Edi-Pro A at 80-100 C emitted 7-13 compounds
   including, tentatively, n-hexane and n-pentane (p. 540). Conditioning removed them below the FID floor
   *before* the study. But n-hexane, n-heptane and n-octane are three of the seven compounds that failed
   significance in the study itself. The paper attributes that failure to short retention near methane, which is
   plausible; **it does not test whether residual endogenous alkane emission at 80-100 C contributed to the
   scatter**, and no re-blanking at 90 or 100 C is reported.
8. **Three points, twenty degrees, two of them from different days.** Every −ΔH in this paper is the slope of a
   line through three temperatures spanning 80-100 C, and the design (p. 540) took 80 and 90 C on one day, then
   90 and 100 C the next after re-equilibrating. The 90 C point is therefore measured twice and the two
   measurements are averaged into one; **no test of day-to-day agreement at 90 C is reported.** Column repacking
   between the three columns is a separate variance source, and it was large enough to invalidate one compound
   outright (2-octanone).
9. **Eighteen of the thirty-five cells in Table II are extrapolations, including two hard zeroes.** The
   hydrocarbon row is extrapolated from C4 to C8 and reads 0, 0, 0.4, 2.4, 4.3 kcal/mol; only C9 and C10 are
   measured. A −ΔH of exactly 0 for n-butane is a regression intercept, not a finding that butane does not
   adsorb on soy protein. Similarly the ketone C9 and C10 cells (12.1, 14.1) are extrapolations that happen to
   sit on the same carbon numbers as Damodaran's measured 2-nonanone — **do not compare an extrapolated cell
   with a measured constant.**
10. **The −ΔG and −ΔS columns are not independent observations.** −ΔG = RT ln V_S exactly (p. 544) and
    −ΔS = (ΔH − ΔG)/T exactly (Table I footnote b). Table I looks like four measurements per compound and is
    two. Any statistical treatment that counts −ΔG or −ΔS as separate evidence is double-counting. This also
    explains the "AS is not a function of temperature" claim on p. 544: **ΔS was constructed to be so**, given
    a ΔH that the van 't Hoff fit holds constant by assumption.
11. **The infinite-dilution regime is not the food regime.** Doses are 10^-9 to 10^-6 g per injection onto
    1.4-1.5 g of protein — a surface coverage so low that the paper calls it "a prerequisite for this type of
    study". Every shipped `REVERSIBLE_BINDING` row comes from a study dosing at 0.05-2.5 mM into a protein
    solution, and Bi 2022 showed retention on pea *reversing direction* between its low and high concentration
    branches. **Nothing here speaks to saturation, competition between odourants, or the number of sites** —
    there is no n, no Klotz plot, no Scatchard plot in this paper.
12. **The isolate is a 1980s commercial product that no longer exists.** Edi-Pro A (Ralston Purina) is a
    spray-dried isoelectric isolate; its protein content, nitrogen factor, thiol, disulfide, amine density,
    residual lipid and moisture are **not reported anywhere in the paper**. `data/species/protein_matrices.yml`
    charges soy site densities from other, modern preparations. Any pairing is cross-preparation and 40 years
    apart, and should be labelled as one.
13. **No DOI is printed** (section 0). Cite by volume/page.
14. **What this paper does NOT contain**: any aqueous phase, pH, ionic strength or buffer; any binding constant
    in M^-1 or L/g; any protein concentration in g/L; any n or number of sites; any rate constant or activation
    energy; any covalent-adduct evidence; any 2-alkenal, pyrazine, pyridine, furan or sulfur compound; any
    measurement below 80 C; any unheated control arm; any printed standard error; any sensory measurement; any
    odour threshold; any supplementary material.
15. **What to request, if the authors or their successors could be reached** (Iowa State, 1983 — realistically
    unobtainable, recorded for completeness): (i) the per-column regression slopes and their standard errors for
    all 17 compounds, not just 2-octanone's three; (ii) the two independent 90 C determinations separately;
    (iii) the protein content and moisture of the Edi-Pro A lot; (iv) whether the isolate was re-blanked at 90
    and 100 C after the initial conditioning.
