# Solina 2007 — EXTRACTION (wheat starch + 1 % soy protein isolate, with 1 % glucose and/or 1 % acid-hydrolysed vegetable protein, twin-screw extruded at 150 C / 20 % moisture and 180 C / 16 % moisture, ~60 s residence; 94 volatiles by dynamic headspace GC-MS in triplicate with SD, ng/10 g; GCO at three dilutions on two extrudates; NO rate constant, NO barrier, NO time series)

### A four-feedstock x two-condition volatile inventory in a real extruder, with standard deviations on every number: it prints the Strecker aldehydes rising 25-fold when free amino acids are added, prints 1 % added glucose doing nothing at all, and prints one 3-methylbutanal level that is HIGHER at 150 C than at 180 C — but temperature and moisture are changed together, so no barrier can be read from any of it.

**Source on disk:** `data/articles/solina2007.pdf` (17 pp., Food Chemistry 104 (2007) 1522-1538).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/solina2007.txt`, 955 lines). Born-digital Elsevier PDF; the running text is
clean. **Tables 1 and 2 are landscape-set six-column tables and the layout extractor packs the
value and its parenthesised SD into ragged columns**; they came through legibly and are re-typed in
full below, but each row had to be assigned to its column by counting value/SD pairs left to right
against the printed column heads (RM, 150/20, 180/16 for each of the two feedstocks). **Where a
row's pairs did not resolve unambiguously it is marked, not guessed** (Flag 1). Tables 3 and 4
(GCO) came through clean and are re-typed verbatim. The Elsevier ligature substitutions are
pervasive ("eﬀects", "Whitﬁeld", "diﬀerent") and are silently normalised. The paper has **no
figures at all** — four tables and running text only. No supplementary material (2007). Repo status
before this dossier: Solina 2007 has **no extraction dossier**; the closely related
`whitfield1988_extraction.md`, `whitfield1999_extraction.md` and `whitfield2001_extraction.md` are
on disk (same corresponding author), and the direct companion Solina, Johnson & Whitfield 2007,
Food Chem. 100:678-692, on which a great deal of this paper's comparison rests, is **not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of soy protein isolate, acid-hydrolysed vegetable protein and glucose on the volatile components of extruded wheat starch" |
| Authors | Marica Solina (The University of Western Sydney, Centre for Advanced Food Research, Locked bag 1797, South Penrith DC, NSW 1797, Australia); Robert L. Johnson; **Frank B. Whitfield** (corresponding; Food Science Australia, P.O. Box 52, North Ryde, NSW 1670, Australia) |
| Journal | Food Chemistry **104** (2007) **1522-1538**; received 20 July 2006, revised 21 January 2007, accepted 22 February 2007 |
| DOI | **10.1016/j.foodchem.2007.02.031** (printed) |
| Copyright | "Crown Copyright © 2007 Published by Elsevier Ltd" |
| Naming | SPI = soy protein isolate; HVP = hydrolysed vegetable protein; **aHVP = acid-hydrolysed vegetable protein**; RM = raw material (the unextruded feedstock); "150/20" and "180/16" = extrusion temperature (C) / feed moisture content (%); LRI = linear retention index; OTC = odour threshold concentration **in water**; SME = specific mechanical energy; GCO = gas chromatography olfactometry; "tr" = trace, < 0.5 ng/10 g; "(-)" = not detected, LOD 0.1 ng/10 g |
| **Same-year companion, heavily relied on, NOT on disk** | Solina M., Johnson R. L., Whitfield F. B., "Effects of glucose and acid-hydrolysed vegetable protein on the volatile components of extruded wheat starch", **Food Chem. 100:678-692 (2007)** — supplies the starch-alone, starch/glucose, starch/aHVP and starch/glucose/aHVP comparators, the starch composition and the reducing-sugar figures. Cited ~20 times. |
| **Ingredient-characterisation paper, NOT on disk** | Solina M., Baumgartner P., Johnson R. L., Whitfield F. B., "Volatile aroma components of soy protein isolate and acid-hydrolysed vegetable protein", **Food Chem. 90:861-873 (2005)** — supplies the SPI and aHVP compositions. Cited ~10 times. |
| Related dossiers on disk | `whitfield1988_extraction.md` (cited here as Whitfield, Mottram, Brock, Puckey & Salter 1988, on phospholipid interference with heterocycle formation), `whitfield1999_extraction.md`, `whitfield2001_extraction.md`, `hofmann1998b_extraction.md` (cited here as Hofmann & Schieberle 1998, for the proline route to 2-propionyl-1-pyrroline), `frankel1989_extraction.md` (the same group's 2-alkylfuran chemistry), `schieberle1989_extraction.md` (the proline / Acp route this paper's pH argument turns on) |

## 1. Why it matters

This is a **real-matrix, real-processing** paper of exactly the kind the repository's scorecard
keeps asking for — a twin-screw extruder, a starch base, protein rather than free amino acid as the
nitrogen source, 16-20 % moisture rather than dilute aqueous — and it is quantitative, in
triplicate, with a standard deviation on every entry. But it contains **no kinetics at all**: two
temperatures, each with its own moisture content, one residence time, no time series. Nothing in it
can be put on an Arrhenius axis and nothing can enter a rate registry.

What it *can* do is test three structural assumptions the trunk makes, and it contradicts one of
them.

**(a) It confirms that protein alone supplies the Strecker nitrogen.** The starch/SPI feedstock
contains **no free amino acids and no reducing sugars** (both measured, in the 2005 companion), yet
extruding it at 180 C produces 3-methylbutanal at 57 ng/10 g, 2-methylbutanal at 16, phenyl-
acetaldehyde at 14 and 2-methylpropanal at 5. The authors' reading is that thermal and mechanical
disruption of the 82.5 % protein liberates enough valine, leucine, isoleucine and phenylalanine,
and that thermal cleavage of the starch supplies the reducing sugar. This is the same regime
`balagiannis2010_extraction.md` and `parker2013_extraction.md` record for a protein-rich meat
extract — free amino acid is not the limiting pool — and it is the regime the isolate programme
will meet.

**(b) It says added glucose at 1 % does essentially nothing, and the authors call this
"perplexing".** At 180 C, adding 1 % glucose moves 3-methylbutanal from **57 to 46 ng/10 g**
(starch/SPI -> starch/glucose/SPI) and from **1420 to 1200 ng/10 g** (starch/SPI/aHVP ->
starch/glucose/SPI/aHVP). **Both changes are downward.** The authors write that "the addition of
glucose at the 1 % level had little effect on the principal chemical reactions occurring in the
starch/SPI system" and that "the apparent absence of glucose decomposition products and the absence
of increased levels of Strecker aldehydes is somewhat perplexing". The resolution the paper itself
offers is that starch cleavage already delivers **75-85 mg reducing sugar per 10 g extrudate**
(second-hand, from the 2007 companion) — i.e. 7.5-8.5 mg/g, against 1 % added glucose = 10 mg/g. So
the added glucose roughly doubles a pool that was already non-limiting. **A model that makes
Strecker aldehyde first order in added reducing sugar will over-predict this experiment.**

**(c) It contradicts a monotonic temperature response in one of the four feedstocks.** For
starch/glucose/SPI/aHVP, 3-methylbutanal is **2960 ng/10 g at 150 C and 1200 ng/10 g at 180 C** —
a factor 2.5 the wrong way — and the total Strecker aldehydes go 3183 -> 1680 ng/10 g, and the
total volatiles 4062 -> 2840 ng/10 g. For starch/SPI/aHVP the same comparison runs the expected way
(3-MB 44 -> 1420). **The two feedstocks that differ only by 1 % glucose disagree in the sign of
their temperature response.** This cannot be read as a barrier of any sign, because moisture falls
from 20 % to 16 % at the same time (Flag 2), and because a headspace measurement of a volatile
aldehyde in a hotter, drier extrudate confounds formation with loss and with matrix binding. It is
nonetheless a real, replicated (n = 3, SD 310 and 320) datum that any predictive model of a
starch/protein extrudate has to survive.

What this paper does **not** give: any rate constant, any activation energy, any time course, any
residence-time variation, any temperature at constant moisture, any pH other than the single
statement "our extrudates have a pH of < 6", any water activity, and any absolute concentration in
the extrudate (every number is a dynamic-headspace recovery — Flag 3).

## 2. Methods as they matter to a model

- **Feedstocks (four).** starch/1 % SPI; starch/1 % glucose/1 % SPI; starch/1 % SPI/1 % aHVP;
  starch/1 % glucose/1 % SPI/1 % aHVP. Mixed 10 min, sifted through a 2 mm sieve. **Calcium
  triphosphate 0.05 % w/w** added to every feedstock "in order to improve the flow properties of the
  starch" — a phosphate salt in a Maillard pot, at 500 mg/kg, and its catalytic effect is nowhere
  discussed (Flag 6). 15 kg of each feedstock processed.
- **Wheat starch** (The Manildra Group, Auburn NSW), "free of all noticeable odours": particle size
  ~75 um, moisture 12.5 %, protein 0.3 % w/w, unbound lipids including free fatty acids 0.2 %.
  **Reducing sugars not detected in the raw starch**; extrusion under mild and extreme conditions
  gave a reducing sugar content in the extrudates of **75 and 85 mg/10 g** (both figures
  second-hand, from Solina et al. 2007, not measured here).
- **SPI** (ADM Protein Specialists, Decatur IL): moisture 10.3 %, **pH 6.6**, protein **82.5 % w/w**,
  **no detectable free amino acids or reducing sugars**, bound lipid 3 % w/w, unbound lipid
  including free fatty acids 0.5 % w/w (all second-hand, from Solina et al. 2005). The discussion
  later quotes the SPI lipid content as "3.5 %" (= 3 + 0.5) and its lipid contribution as
  "300 ng/10 g"; the latter unit is a volatile amount, not a lipid content, and reads as a slip.
- **aHVP** (Halcyon Protein, Dandenong VIC), derived from soy protein: moisture 6.4 % w/w,
  **free amino acids 18.4 % w/w**, unbound lipid 0.4 % w/w, **free of reducing sugars** (second-hand,
  Solina et al. 2005). **So 1 % aHVP delivers 0.184 % free amino acid to the feedstock = 1840 mg/kg
  = 1.84 mg/g.** Proline is named as a component.
- **Extruder.** APV Baker MPF 40 **co-rotating twin-screw**. Screw speed **225 rpm**, held. Overall
  screw length/diameter **20D**. Screw configuration from feed section to die: 3D feed screws, 1D
  lead screws, 3 forward paddles at 60 deg, 2D lead screws, 4 forward paddles at 60 deg, 2D lead
  screws, **4 reverse paddles at 60 deg**, 1D lead screws. **6 mm die.** Mass temperature read from
  **six probes in contact with the fluid mass** along the barrel, and barrel sections heated or
  cooled to hold a constant mass temperature in each — so the quoted 150 and 180 C are *mass*
  temperatures, not barrel set-points, which is better than most extrusion papers.
- **Residence time.** "The standard screw configuration gave a **median retention time of 60 s** at
  a feed rate of 15 kg/h, as obtained from the residence time distribution of the marker compound,
  Erythrosin B." **One residence time only.** The distribution itself is not printed.
- **The two conditions.** Mild = **150 C and 20 % moisture**; extreme = **180 C and 16 % moisture**.
  **These are changed together and never separated** (Flag 2). SME varied between **26 % and 27 %**
  across these feedstocks. Die pressure highest at **240 psi** (wheat starch/SPI, mild) and lowest
  at **40-100 psi** (starch/glucose/SPI, extreme). Extrudates collected over 3-5 min while the
  variables were steadiest, cooled, mixed, milled to a homogeneous powder, sealed in
  polyethylene-foil-polyester laminate bags, **stored at -20 C** (printed as "20 °C"; read as -20 C,
  Flag 1) until analysis.
- **Volatile collection — dynamic headspace, and this is the key limitation.** 10 g of powdered
  sample + **80 ml water** in a 250 ml conical flask, stirred slowly, **held at 37 C in a water
  bath**, purged with **oxygen-free nitrogen at 40 ml/min for 1 h** onto **10 mg Tenax TA** in a
  glass-lined stainless-steel tube (115 mm x 0.75 mm i.d.); then 5 min of dry nitrogen to remove
  residual water. **Chlorododecane (100 ng in 100 ul ethanol) added to the flask** as a recovery
  internal standard; **chlorotetradecane (100 ng in 1 ul pentane) added to the front of the trap
  just before analysis** as the quantification standard. **The reported recovery from the
  chlorododecane is never given anywhere in the paper** (Flag 3).
- **Quantification.** Peak area against the 100 ng chlorotetradecane, **"assuming all response
  factors were 1"**, reported as **ng/10 g sample**. Footnote f of both tables: "The GC-MS response
  factors for each component are assumed to be 1:1. Consequently, the reported quantities are
  considered as approximate values." **LOD 0.1 ng/10 g; "tr" < 0.5 ng/10 g. Averages of three
  separate isolations, with SD printed in parentheses.**
- **GC-MS.** HP 5890 series II Plus + HP 5972 MSD, G 1701 BA ChemStation; **HP5-Trace Analysis
  25 m x 0.2 mm i.d., 1 um film** plus a 30 cm x 0.32 mm uncoated deactivated retention gap. Tenax
  desorbed **10 min at 280 C** with the pre-column cooled to **-78 C** with solid carbon dioxide
  (printed as "78 °C"; read as -78 C). Oven 40 C, then **5 C/min to 280 C**, hold 5 min. EI 70 eV,
  emission 50 uA, source 250 C, scan 35-400 amu at 1 s/decade. **LRI from n-alkanes C5-C24 run under
  the same conditions.** Identification: mass spectrum against NIST/EPA/NIH and Wiley, then LRI
  against authentic compounds or published values; the tables mark each row **"MS + LRI"** or
  **"MS"** (mass spectrum only — a weaker identification, and 10 rows in Table 1 and 8 in Table 2
  carry it).
- **GCO.** HP 5890 series II Plus with a CHIS injection port and a humidified odour port; same
  pre-column, column and oven as the GC-MS; effluent 0.6 ml/min **split 1:8 (v/v) between FID and
  odour port**; detector and transfer line 250 C; humidified air 40 ml/min at the port. **Five
  pre-screened assessors**, each evaluating extracts from **three sample sizes: 10, 1 and 0.1 g**
  (i.e. undiluted, one tenth, one hundredth) of the two 180 C extrudates only. Assessors described
  odours in their own words; **only odours reported by at least three of the five assessors are
  listed**. Intensity scale: "low", "moderate", "strong", "very strong" (tabulated as W, M, S, VS,
  with an S/M and an M/W category also appearing).
- **Qualitative sensory.** 10 g of powder in a screw-capped jar, wetted with 30 ml water just before
  assessment; **four panellists**, room temperature, agreed descriptive terms.
- **Extrudate pH.** Stated once, in the discussion: **"our extrudates have a pH of < 6"**. No value,
  no method, no per-feedstock figure.

## 3. Tables re-typed

Units for Tables 1 and 2 throughout: **ng/10 g sample**, with the SD of three isolations in
parentheses. "-" = not detected (LOD 0.1 ng/10 g); "tr" = < 0.5 ng/10 g. Column heads: **RM** = raw
material (unextruded feedstock); **150/20** and **180/16** = temperature (C) / moisture (%).

### Table 1 (pp. 1526-1528). "Relative concentrations of headspace volatiles of starch/SPI and starch/glucose/SPI feedstocks extruded under different conditions of temperature and moisture content"

| Identity | LRI | S/SPI RM | S/SPI 150/20 | S/SPI 180/16 | S/Glc/SPI RM | S/Glc/SPI 150/20 | S/Glc/SPI 180/16 | ID method |
|---|---:|---|---|---|---|---|---|---|
| **Lipid-derived — Aldehydes** | | | | | | | | |
| Pentanal | 722 | 8 (1.7) | 2 (2.7) | 17 (4.9) | 2 (1.1) | 39 (8.1) | 11 (2.0) | MS + LRI |
| Hexanal | 817 | 220 (44) | 170 (40) | 230 (5.5) | 150 (36.4) | 280 (21) | 250 (34) | MS + LRI |
| Heptanal | 907 | 7 (1.7) | 21 (1.0) | 18 (4.2) | 4 (2.9) | 14 (2.9) | 23 (5.4) | MS + LRI |
| Octanal | 1008 | 3 (0.1) | 2 (0.2) | 5 (0.1) | 3 (0.5) | 3 (<0.1) | 4 (1.2) | MS + LRI |
| Nonanal | 1105 | 10 (2.6) | 5 (0.1) | 9 (<0.1) | 12 (2.9) | - | - | MS + LRI |
| Decanal | 1204 | 2 (0.7) | 3 (0.6) | 2 (0.3) | 3 (0.7) | 1 (0.2) | 2 (0.6) | MS + LRI |
| (E)-2-Hexenal | 858 | - | 1 (0.1) | 1 (<0.1) | - | 1 (0.2) | 1 (0.2) | MS + LRI |
| (E)-2-Heptenal | 961 | 3 (1.3) | 5 (0.8) | 3 (0.4) | 2 (0.5) | 5 (0.1) | 6 (1.5) | MS + LRI |
| (Z)-2-Octenal | 1062 | 6 (1.6) | 10 (1.5) | 15 (1.4) | 6 (0.8) | 9 (0.4) | 11 (1.8) | MS + LRI |
| (E)-2-Nonenal | 1162 | 10 (1.2) | 12 (1.3) | 10 (2.9) | 12 (2.5) | 8 (0.6) | 16 (2.2) | MS + LRI |
| (E,E)-2,4-Octadienal | 1120 | - | 1 (0.1) | 1 (<0.1) | tr (0.3) | - | - | MS + LRI |
| (E,E)-2,4-Nonadienal | 1217 | 1 (0.3) | 1 (0.1) | - | 1 (0.2) | - | - | MS + LRI |
| (E,Z or Z,E)-2,4-Decadienal | 1295 | - | 22 (3.0) | 13 (1.0) | tr (0.1) | 30 (5.4) | 52 (14) | MS |
| (E,E)-2,4-Decadienal | 1317 | 2 (0.5) | 36 (4.9) | 29 (3.3) | 3 (0.3) | 59 (8.1) | 60 (15) | MS + LRI |
| **Ketones** | | | | | | | | |
| 2-Pentanone | 699 | - | - | 4 (<0.1) | - | - | - | MS + LRI |
| 2-Hexanone | 802 | - | - | 2 (0.1) | - | - | 1 (0.2) | MS + LRI |
| 3-Heptanone | 886 | - | - | 1 (<0.1) | - | - | tr (<0.1) | MS + LRI |
| 2-Heptanone | 898 | 12 (2.6) | 24 (3.3) | 48 (2.0) | 6 (2.8) | 22 (0.4) | 35 (3.8) | MS + LRI |
| 3-Octen-2-one | 1045 | - | 3 (0.5) | 4 (0.4) | - | 2 (0.4) | 4 (0.3) | MS + LRI |
| 2-Nonanone | 1093 | - | 2 (<0.1) | 3 (0.1) | 2 (1.5) | - | - | MS + LRI |
| (E,E)-3,5-Octadien-2-one | 1100 | - | - | 1 (0.1) | - | - | - | MS + LRI |
| **Alcohols** | | | | | | | | |
| 1-Pentanol | 783 | 2 (1.8) | 3 (0.6) | 4 (0.9) | 3 (1.3) | 2 (0.5) | 3 (2.2) | MS + LRI |
| 1-Hexanol | 880 | - | 7 (0.2) | 11 (1.2) | 9 (4.5) | 13 (2.2) | 11 (3.1) | MS + LRI |
| 1-Heptanol | 974 | 2 (0.4) | - | - | 2 (0.7) | - | - | MS + LRI |
| 1-Octen-3-ol | 984 | 5 (<0.1) | 3 (0.3) | 3 (0.7) | 4 (0.8) | 2 (0.3) | 3 (0.6) | MS + LRI |
| 1-Octanol | 1071 | 2 (0.2) | - | - | 1 (0.5) | - | - | MS + LRI |
| 1-Nonanol | 1160 | - | - | - | 1 (0.6) | - | - | MS + LRI |
| **Furans** | | | | | | | | |
| 2-Methylfuran | <650 | - | - | 1 (0.9) | - | - | - | MS + LRI |
| 2-Ethylfuran | 699 | - | 3 (2.5) | 21 (6.0) | - | - | 25 (5.7) | MS + LRI |
| 2,5-Dimethylfuran | 706 | - | - | - | - | - | 6 (1.2) | MS + LRI |
| 2-Propylfuran | 791 | - | - | 3 (0.7) | - | 1 (0.9) | 5 (1.3) | MS + LRI |
| 2-Pentylfuran | 995 | 110 (8.3) | 64 (9.3) | 110 (27) | 74 (9.2) | 92 (9.8) | 150 (22) | MS + LRI |
| **Sugar-derived — Ketones** | | | | | | | | |
| 2,3-Octanedione | 991 | - | - | 7 (0.8) | - | - | 5 (0.1) | MS + LRI |
| **Furans** | | | | | | | | |
| 2-Furfural | 841 | - | - | - | - | - | 2 (0.4) | MS + LRI |
| **Amino acid-derived — Aldehydes** | | | | | | | | |
| 2-Methylpropanal | <650 | - | - | **5 (1.2)** | - | - | **7 (1.3)** | MS + LRI |
| 3-Methylbutanal | 669 | - | tr (0.5) | **57 (5.0)** | - | **17 (0.9)** | **46 (13)** | MS + LRI |
| 2-Methylbutanal | 677 | - | - | **16 (0.2)** | - | **3 (0.2)** | **33 (5.3)** | MS + LRI |
| Benzaldehyde ^f | 979 | 5 (0.6) | 7 (1.1) | 15 (2.0) | 5 (0.7) | 10 (0.3) | 16 (2.1) | MS + LRI |
| Phenylacetaldehyde | 1063 | - | - | **14 (0.7)** | - | - | **-** | MS + LRI |
| **Maillard reaction-derived — Pyrroles** | | | | | | | | |
| 1-H-Pyrrole | 771 | - | - | - | - | - | 1 (0.2) | MS + LRI |
| **Pyridines** | | | | | | | | |
| Pyridine | 765 | 48 (17) | 7 (1.5) | 88 (14) | 18 (7.9) | 41 (3.6) | 9 (0.5) | MS + LRI |
| 2-Propylpyridine | 1201 | - | - | 1 (<0.1) | - | - | - | MS |
| **Pyrazines** | | | | | | | | |
| Pyrazine | 736 | - | - | 5 (1.1) | - | - | 1 (0.1) | MS + LRI |
| Methylpyrazine | 839 | - | - | 3 (0.6) | - | - | 8 (2.6) | MS + LRI |
| **Sulphur-containing — Aliphatic** | | | | | | | | |
| Dimethyl disulphide | 756 | 3 (1.1) | 2 (0.4) | 23 (3.3) | - | 4 (1.2) | 30 (3.1) | MS + LRI |
| Dimethyl trisulphide | 975 | 1 (<0.1) | 1 (0.3) | 3 (0.1) | - | - | 1 (0.1) | MS + LRI |
| Methyl pentyl disulphide | 1137 | - | 1 (<0.1) | 2 (<0.1) | - | - | - | MS |
| **Thiazoles** | | | | | | | | |
| Benzothiazole | 1257 | - | - | - | 1 (0.1) | - | - | MS + LRI |
| **Thiophenes** | | | | | | | | |
| Thiophene | <650 | - | - | 4 (0.6) | - | - | - | MS + LRI |
| 3-Methylthiophene | 805 | - | - | tr (<0.1) | - | - | tr (0.1) | MS |
| 2-Pentylthiophene | 1164 | - | - | 4 (0.8) | - | - | - | MS + LRI |
| **Thiapyrans** | | | | | | | | |
| 2-Pentylthiapyran | 1324 | - | - | 1 (<0.1) | - | - | - | MS + LRI |
| **Derived from other sources — Aldehydes** | | | | | | | | |
| 2-Butylacrolein | 847 | - | - | 5 (0.2) | - | - | 3 (0.7) | MS |
| 2-Butyl-2-octenal | 1367 | 4 (0.2) | 2 (0.3) | 2 (0.4) | 4 (0.5) | 1 (0.3) | - | MS + LRI |
| 4-Pentylbenzaldehyde | 1505 | - | - | 1 (<0.1) | - | - | 2 (0.3) | MS |
| **Ketones** | | | | | | | | |
| 3-Cyclohepten-1-one | 823 | - | - | 2 (0.5) | - | - | 2 (0.6) | MS |
| 5-Decanone | 1179 | 1 (<0.1) | - | - | - | - | - | MS |
| (E,E)-6,10-Dimethyl-5,9-undecadien-2-one | 1443 | 1 (0.2) | 1 (0.1) | 1 (0.9) | 2 (0.1) | - | - | MS |
| 2,6-bis(1,1-Dimethyl)-2,5-cyclohexadiene-1,4-dione | 1461 | 1 (0.1) | 1 (1.1) | 1 (0.2) | 1 (0.1) | - | - | MS |
| 2-Tridecanone | 1475 | - | 1 (0.1) | - | - | - | - | MS |
| **Alcohols** | | | | | | | | |
| 3-Methylbutanol | 722 | - | - | - | tr (0.1) | - | - | MS + LRI |
| **Furans** | | | | | | | | |
| 2,2,4,4-Tetramethyltetrahydrofuran | 783 | - | - | 10 (0.1) | - | - | - | MS |
| **Phenols** | | | | | | | | |
| Phenol | 756 | - | - | 1 (<0.1) | - | - | 2 (0.2) | MS + LRI |
| **Hydrocarbons** | | | | | | | | |
| Toluene | 784 | 11 (11) | 4 (0.6) | 15 (0.9) | 4 (3.3) | 7 (0.3) | 17 (2.1) | MS + LRI |
| Octane | 800 | 4 (2.9) | 1 (0.3) | 3 (2.1) | tr (0.2) | - | 7 (3.5) | MS + LRI |
| 1,4-Dimethylbenzene | 878 | 18 (6.2) | 8 (2.8) | 15 (2.1) | 7 (4.9) | 8 (3.2) | 11 (6.2) | MS + LRI |
| alpha-Pinene | 937 | tr (0.2) | - | 1 (0.1) | - | - | - | MS + LRI |
| Limonene | 1034 | 4 (2.4) | 9 (1.7) | 26 (1.1) | 4 (1.0) | 3 (0.3) | 10 (1.3) | MS + LRI |
| **Miscellaneous** | | | | | | | | |
| Ethyl acetate | <650 | 36 (53) | - | 13 (2.1) | 5 (4.6) | - | - | MS + LRI |
| Propyl acetate | 707 | - | - | 2 (<0.1) | - | - | - | MS + LRI |
| 1,1-Diethoxyethane | 738 | - | - | - | 1 (0.5) | - | - | MS + LRI |
| Butyl acetate | 816 | - | - | 1 (0.2) | - | - | - | MS + LRI |
| 1-tert-Butoxy-2-methoxyethane | 844 | 22 (15) | 14 (1) | 120 (27) | 3 (0.3) | 64 (16.3) | 54 (12) | MS |
| Hexanenitrile | 882 | - | - | 2 (0.3) | - | 2 (<0.1) | 2 (0.5) | MS |
| 1-tert-Butoxy-2-ethoxyethane | 904 | - | - | 40 (17) | - | 16 (6.3) | 20 (5.8) | MS |
| 1-Nitropentane | 947 | - | - | 2 (0.2) | - | - | 2 (0.5) | MS |
| 1-Nitrohexane | 1050 | 7 (0.9) | 7 (0.3) | - | 8 (0.4) | 9 (2.2) | 8 (1.4) | MS |

Footnotes as printed: a) Linear retention index. b) MS + LRI: identified by comparison of mass
spectra and LRI with those of an authentic compound or previously published data; MS: mass spectrum
agrees with the reference spectrum from the NIST/EPA/NIH Mass Spectral Database. c) Raw material.
d) Extrusion variables, temperature (°C)/moisture content (%). e) Concentration (ng/10 g) obtained
by comparing GC-MS peak area with that from 100 ng chlorotetradecane internal standard added to the
Tenax trap after volatile collection; the averages of triplicate analyses are shown; (-) not
detected (limit of detection 0.1 ng/10 g sample); (tr) volatiles in concentrations of < 0.5 ng/10 g.
f) **Benzaldehyde may be either lipid- or amino acid-derived.** (The GC-MS response factors for each
component are assumed to be 1:1. Consequently, the reported quantities are considered as
approximate values).

**Note on 2-ethylfuran, starch/glucose/SPI RM.** The extracted line reads
`2-Ethylfuran 699 – – 3 (2.5) 21 (6.0) – – 25 (5.7)` with one dash-pair short of the six columns;
the RM and 150/20 entries of the starch/glucose/SPI block cannot both be resolved from the text
layer. The 180/16 value 25 (5.7) is unambiguous (it is the last pair). Recorded as **not resolved**
for the two earlier columns rather than guessed (Flag 1).

### Table 2 (pp. 1530-1532). "Relative concentrations of headspace volatiles of starch/SPI/aHVP and starch/glucose/SPI/aHVP feedstocks extruded under different conditions of temperature and moisture content"

| Identity | LRI | S/SPI/aHVP RM | S/SPI/aHVP 150/20 | S/SPI/aHVP 180/16 | S/Glc/SPI/aHVP RM | S/Glc/SPI/aHVP 150/20 | S/Glc/SPI/aHVP 180/16 | ID method |
|---|---:|---|---|---|---|---|---|---|
| **Lipid-derived — Aldehydes** | | | | | | | | |
| Pentanal | 722 | 3 (1.5) | 17 (2.4) | 12 (1.0) | 13 (4.4) | 38 (6.9) | 6 (0.5) | MS + LRI |
| Hexanal | 817 | 210 (12) | 180 (29) | 250 (35) | 210 (10) | 310 (6.7) | 310 (39) | MS + LRI |
| Heptanal | 907 | 7 (0.3) | 17 (0.8) | 48 (11) | 7 (1.4) | 24 (0.6) | 41 (7.7) | MS + LRI |
| Octanal | 1008 | 3 (0.2) | 2 (0.3) | 17 (4.8) | 2 (1.2) | 4 (0.2) | 14 (3.8) | MS + LRI |
| Nonanal | 1105 | 8 (2.2) | 3 (0.4 | 14 (3.5) | 10 (2.5) | 6 (0.8) | 15 (1.9) | MS + LRI |
| Decanal | 1204 | 1 (0.2) | 1 (0.1) | 3 (0.3) | 2 (0.4) | 2 (0.3) | - | MS + LRI |
| (E)-2-Hexenal | 858 | - | 1 (0.1) | 4 (0.3) | - | 1 (<0.1) | - | MS + LRI |
| (E)-2-Heptenal | 961 | 3 (0.4) | 4 (0.2) | 3 (0.8) | 3 (0.4) | 5 (0.2) | 4 (0.9) | MS + LRI |
| (E)-2-Octenal | 1062 | 7 (0.3) | 12 (1.8) | 23 (4.2) | 6 (0.8) | 11 (0.8) | 32 (5.9) | MS + LRI |
| (E)-2-Nonenal | 1162 | 14 (0.6) | 12 (1.3) | 21 (4.9) | 13 (2.1) | 16 (3.2) | 24 (5.1) | MS + LRI |
| (E)-2-Decenal | 1269 | - | tr (0.1) | - | - | - | - | MS + LRI |
| (E,E)-2,4-Octadienal | 1120 | - | 1 (<0.1) | - | - | - | - | MS + LRI |
| (E,E)-2,4-Nonadienal | 1217 | 1 (<0.1) | 1 (0.1) | 2 (0.3) | 1 (0.1) | - | - | MS + LRI |
| (E,Z or Z,E)-2,4-Decadienal | 1295 | 1 (0.2) | - | 14 (6.1) | 1 (0.1) | 32 (3.0) | 31 (9.0) | MS |
| (E,E)-2,4-Decadienal | 1317 | 4 (0.3) | 79 (3.3) | 59 (9.2) | 4 (0.5) | 63 (5.4) | 69 (16) | MS + LRI |
| **Ketones** | | | | | | | | |
| 3-Heptanone | 886 | - | - | 2 (0.5) | - | - | 1 (0.3) | MS + LRI |
| 2-Heptanone | 898 | 12 (1.8) | 21 (1.3) | 68 (12) | 15 (2.4) | 27 (0.3) | 41 (7.3) | MS + LRI |
| 3-Octen-2-one | 1045 | - | 3 (0.2) | - | - | 3 (0.2) | - | MS + LRI |
| 2-Nonanone | 1093 | 2 (0.1) | - | 4 (0.6) | 3 (0.3) | - | 4 (0.7) | MS + LRI |
| 3-Nonen-2-one | 1140 | - | 1 (0.1) | - | - | - | - | MS + LRI |
| **Alcohols** | | | | | | | | |
| 1-Pentanol | 783 | 3 (0.5) | - | 2 (1.8) | 4 (0.7) | 7 (1.6) | 5 (1.0) | MS + LRI |
| 1-Hexanol | 880 | 11 (2.2) | 7 (1.5) | 20 (5.7) | 15 (2.7) | 6 (1.3) | 19 (4.9) | MS + LRI |
| 1-Heptanol | 974 | 2 (<0.1) | - | - | 2 (0.4) | - | - | MS + LRI |
| 1-Octen-3-ol | 984 | 6 (0.3) | 2 (0.5) | 3 (0.6) | 6 (0.7) | 2 (0.5) | - | MS + LRI |
| 1-Octanol | 1071 | 1 (0.1) | - | - | 2 (1.6) | - | - | MS + LRI |
| 1-Nonanol | 1160 | 1 (<0.1) | - | - | 1 (0.1) | - | - | MS + LRI |
| **Furans** | | | | | | | | |
| 2-Ethylfuran | 699 | - | - | - | - | 41 (10) | - | MS + LRI |
| 2-Propylfuran | 791 | - | - | - | - | 5 (1.2) | 4 (0.8) | MS + LRI |
| 2-Pentylfuran | 995 | 100 (2.2) | 80 (13) | 140 (25) | 86 (14) | 120 (10) | 170 (40) | MS + LRI |
| **Sugar-derived — Ketones** | | | | | | | | |
| 2,3-Octanedione | 991 | - | - | 5 (1.1) | - | - | 4 (1.0) | MS + LRI |
| **Furans** | | | | | | | | |
| 1-(2-Furanyl)-ethanone [= 2-acetylfuran] | 916 | - | - | 3 (0.6) | - | - | 3 (0.9) | MS |
| **Amino acid-derived — Aldehydes** | | | | | | | | |
| 2-Methylpropanal | <650 | - | - | **120 (8.2)** | - | - | **150 (40)** | MS + LRI |
| 3-Methylbutanal | 669 | 2 (1.8) | **44 (11)** | **1420 (210)** | 5 (0.2) | **2960 (310)** | **1200 (320)** | MS + LRI |
| 2-Methylbutanal | 677 | - | **17 (4.1)** | **99 (1.4)** | - | **147 (40)** | **110 (30)** | MS + LRI |
| Benzaldehyde ^f | 979 | 6 (1.1) | 8 (1.8) | 29 (6.0) | 7 (0.8) | 14 (2.1) | 26 (4.2) | MS + LRI |
| Phenylacetaldehyde | 1063 | - | **16 (0.8)** | **240 (63)** | - | **76 (5.9)** | **220 (31)** | MS + LRI |
| **Maillard reaction-derived — Pyrroles** | | | | | | | | |
| 1-H-Pyrrole | 771 | - | - | 9 (2.5) | - | - | 5 (0.4) | MS + LRI |
| 1-Ethyl-1H-pyrrole | 829 | - | - | 4 (0.5) | - | - | 4 (0.6) | MS + LRI |
| **Pyridines** | | | | | | | | |
| Pyridine | 765 | 9 (4.6) | 8 (0.7) | 4 (7.2) | 21 (5.8) | - | - | MS + LRI |
| **Pyrazines** | | | | | | | | |
| Methylpyrazine | 839 | - | - | 28 (6.1) | - | 3 (SD not printed) | 34 (10) | MS + LRI |
| 2,5- or 2,6-Dimethylpyrazine | 916 | - | - | 6 (1.5) | - | - | 23 (4.3) | MS + LRI |
| Ethylpyrazine | 919 | - | - | 19 (4.3) | - | - | 22 (4.9) | MS + LRI |
| Ethenylpyrazine | 934 | - | - | 3 (0.7) | - | - | 3 (0.6) | MS + LRI |
| 2-Vinyl-6-methylpyrazine | 1016 | - | - | 6 (1.2) | - | - | 5 (1.4) | MS + LRI |
| 3-Ethyl-2,5-dimethylpyrazine | 1081 | - | - | 11 (2.7) | - | - | - | MS + LRI |
| 2-(3-Methylbutyl)-6-methylpyrazine | 1248 | - | - | 5 (1.2) | - | - | - | MS |
| **Oxazoles** | | | | | | | | |
| 2,4,5-Trimethyloxazole | 848 | - | - | 3 (0.5) | - | - | - | MS + LRI |
| **Sulphur-containing — Aliphatic** | | | | | | | | |
| Dimethyl disulphide | 756 | 1 (0.8) | 12 (1.3) | - | 3 (1.3) | **100 (15)** | **47 (2.6)** | MS + LRI |
| Dimethyl trisulphide | 975 | 1 (0.2) | 2 (0.2) | 16 (4.2) | 2 (0.7) | - | 16 (3.4) | MS + LRI |
| Methyl pentyl disulphide | 1137 | - | 1 (<0.1) | - | - | - | - | MS |
| **Thiazoles** | | | | | | | | |
| Benzothiazole | 1257 | 1 (0.1) | 1 (0.1) | - | 1 (<0.1) | - | - | MS + LRI |
| **Derived from other sources — Aldehydes** | | | | | | | | |
| 3-Methyl-2-butenal | 806 | - | - | 8 (1.5) | - | - | - | MS + LRI |
| 2-Butyl-2-octenal | 1367 | 4 (0.8) | 2 (0.2) | 3 (0.7) | 5 (0.6) | - | 5 (1.2) | MS + LRI |
| **Ketones** | | | | | | | | |
| 5-Methyl-2-hexanone | 857 | - | - | - | - | 1 (0.1) | 11 (2.7) | MS + LRI |
| (E,E)-6,10-Dimethyl-5,9-undecadien-2-one | 1443 | 1 (0.1) | 1 (0.3) | - | 1 (0.2) | - | - | MS |
| 2,6-bis(1,1-Dimethyl)-2,5-cyclohexadiene-1,4-dione | 1461 | 1 (0.2) | - | - | 1 (0.3) | tr (0.2) | - | MS |
| **Furans** | | | | | | | | |
| 2,5-Dihydrofuran | 918 | - | - | - | - | - | 10 (1.0) | MS |
| 3-Phenylfuran | 1228 | - | 1 (0.1) | 15 (4.1) | - | 8 (0.2) | 40 (7.5) | MS |
| **Hydrocarbons** | | | | | | | | |
| Toluene | 784 | 5 (3.1) | 3 (1.7) | 12 (2.4) | 11 (5.5) | 5 (0.8) | 26 (5.0) | MS + LRI |
| Octane | 800 | 1 (1.3) | - | 1 (1.9) | 1 (0.2) | - | - | MS + LRI |
| 1,4-Dimethylbenzene | 878 | 13 (11) | 5 (2.9) | 36 (15) | 9 (2.5) | 10 (5.0) | 40 (12) | MS + LRI |
| alpha-Pinene | 937 | - | - | - | tr (0.1) | - | - | MS + LRI |
| Limonene | 1034 | 4 (0.2) | 4 (0.4) | 5 (0.7) | 2 (0.5) | 5 (0.4) | 5 (1.7) | MS + LRI |
| **Miscellaneous** | | | | | | | | |
| Ethyl acetate | <650 | 72 (23) | - | - | 14 (13) | - | - | MS + LRI |
| Propyl acetate | 707 | - | - | 4 (0.8) | - | - | 2 (1.8) | MS + LRI |
| 1,1-Diethoxyethane | 738 | tr (0.3) | - | - | - | - | - | MS + LRI |
| 1-tert-Butoxy-2-methoxyethane | 844 | 8 (1.5) | 5 (2.0) | 54 (7.0) | 12 (1.7) | 8 (2.2) | 42 (8.4) | MS |
| Hexanenitrile | 882 | - | - | 4 (0.7) | 2 (2.7) | 2 (0.3) | 3 (0.1) | MS |
| 1-Nitropentane | 947 | 3 (0.4) | 2 (0.4) | - | 3 (0.3) | - | 2 (0.4) | MS |
| 1-Nitrohexane | 1050 | 9 (0.5) | - | - | 7 (0.3) | - | - | MS |
| 2,3-Dihydro-1H-indole | 1114 | - | tr (<0.1) | 4 (0.9) | - | - | - | MS |

Footnotes identical to Table 1's a-f.

**Two typographic defects in Table 2, as printed.** (i) Nonanal, starch/SPI/aHVP 150/20 prints
`3 (0.4` — the closing parenthesis is missing; the value 3 and the SD 0.4 are both legible.
(ii) Methylpyrazine, starch/glucose/SPI/aHVP 150/20 prints `3 –`, i.e. a value of 3 with a dash in
the SD position; **no SD is printed for that one cell.** Neither is guessed.

### Table 3 (p. 1534). "GCO analysis of volatile components at different concentrations of starch/glucose/SPI extrudates processed under extreme conditions"

Intensity key as printed: VS = "very strong", S = "strong", M = "moderate", W = "weak", "-" = "not
present". Footnote c: "Odour description need not necessarily relate to compound(s) identified in
this region of the chromatogram."

| LRI | Odour description | 10 g | 1 g | 0.1 g | Major compound in region of odour |
|---:|---|---|---|---|---|
| 585 | Cooked meat, sweet | M/W | W | - | Unknown |
| 625 | Caramel, fruity | M | - | - | Unknown |
| 645 | Over-ripe apple, musty on dilution | W | W | - | 2-Methylpropanal |
| 669 | Burnt toffee, caramel on dilution | M | W | - | 3-Methylbutanal |
| 677 | Toffee apple, toffee, fruity on dilution | S | M | W | 2-Methylbutanal |
| 682 | Crushed insect | **VS** | S | M | **Unknown** |
| 690 | Caramel, stale biscuit on dilution | M | W | - | Unknown |
| 699 | Apple, pear, caramel, toffee on dilution | M | W | - | 2-Ethylfuran |
| 756 | Sulphury, rubber-like on dilution | S/M | W | - | Dimethyl disulphide |
| 771 | Solvent-like, apple, sweet and piercing on dilution | S | S | - | 1-H-Pyrrole |
| 817 | Unpleasant, acrid crushed ant-like grassy on dilution | M | W | - | Hexanal |
| 825 | Unpleasant, sulphur-like, sweet oil | M | - | - | Unknown |
| 877 | Acrid, crushed ants, sweet and berry-like on dilution | S | W | - | Unknown |
| 898 | Acrid, egg-like, sweet, wine, cheese-like on dilution | S/M | W | - | 2-Heptanone |
| 907 | Woody, leather | S | M | - | Heptanal |
| 961 | Putty | W | W | - | (E)-2-Heptenal |
| 984 | Mushroom-like | S/M | S | M | 1-Octen-3-ol |
| 995 | Metallic, green | S | W | - | 2-Pentylfuran |
| 1008 | Sweet, citrus, honey-like | M | M | - | Octanal |
| 1062 | Crushed ants, green and capsicum-like on dilution | S | W | - | (E)-2-Octenal |
| 1103 | Green, ant-like, slightly rancid | M | - | - | Unknown |
| 1162 | Fried vegetables | W | - | - | (E)-2-Nonenal |
| 1204 | Sweaty rubber, sweet on dilution | S | W | - | Decanal |
| 1295 | Fatty | **VS** | M | W | (E,Z or Z,E)-2,4-Decadienal |
| 1317 | Sweet, citrus | S | M | - | (E,E)-2,4-Decadienal |

### Table 4 (p. 1535). "GCO analysis of volatile components at different concentrations of starch/glucose/SPI/aHVP extrudates processed under extreme conditions"

| LRI | Odour description | 10 g | 1 g | 0.1 g | Major compound in region of odour |
|---:|---|---|---|---|---|
| 560 | Sweet cabbage | M | - | - | Unknown |
| 620 | Tomato | S/M | M | W | Unknown |
| 645 | Vegetable-like, tomato | M | W | W | 2-Methylpropanal |
| 669 | Sweet, caramel, condensed milk, golden syrup-like on dilution | S/M | S | **S** | **3-Methylbutanal** |
| 677 | Sweet, nail varnish, pear-like on dilution | S | S | W | 2-Methylbutanal |
| 681 | Crushed ants | S/M | M | M | Unknown |
| 707 | Vegetable, rubber, onion-like | W | W | W | Propyl acetate |
| 722 | Sweet caramel | S | W | W | Pentanal |
| 753 | Fried onions, sulphury | S/M | - | - | Unknown |
| 756 | Pungent, flue gas, sweaty, vegetable on dilution | M | M | W | Dimethyl disulphide |
| 771 | Fruity, apple, cinnamon, biscuits | M | - | - | 1-H-Pyrrole |
| 800 | Sweet green, green apples | M | M | W | Unknown |
| 817 | Green, rancid | S | W | - | Hexanal |
| 903 | Curried apple, onion-like | W | - | - | Unknown |
| 907 | Sweet biscuits | **VS** | S | M | Heptanal |
| 916 | Gravy, meat stew | S | S | - | 2,5 or 2,6-Dimethylpyrazine |
| 919 | Oasts [sic], biscuits, tobacco-like on dilution | S | M | W | Ethylpyrazine |
| 961 | Cut grass | S | - | - | (E)-2-Heptenal |
| 975 | Roast beef | **VS** | S | W | Dimethyl trisulphide |
| 984 | Mushrooms | M | M | W | 1-Octen-3-ol |
| 995 | Metallic, vegetable | **VS** | M | - | 2-Pentylfuran |
| 1008 | Stink bug, citrus-like | S | M | - | Octanal |
| 1016 | Biscuit | W | - | - | 2-Vinyl-6-methylpyrazine |
| 1062 | Savoury, pizza, cheese cracker | S | S | W | (E)-2-Octenal |
| 1105 | Soft, floral, citrus, fatty | S | M | - | Nonanal |
| 1115 | Toasted cheese, burger ring, onion- and garlic-like on dilution | S | S | W | **Unknown** |
| 1162 | Cracker cheese biscuit | M | M | - | (E)-2-Nonenal |
| 1317 | Fried vegetables | S | M | - | (E,E)-2,4-Decadienal |

### Totals and shares printed in the running text

| feedstock | state | compounds identified | total volatiles (ng/10 g) | lipid-derived | Strecker aldehydes | Maillard heterocycles |
|---|---|---:|---:|---|---|---|
| starch/SPI | RM | 33 | — (lipid fraction 571 = 71 %) | 17 compounds, 571 ng/10 g, 71 % | none | none |
| starch/SPI | 150/20 | 38 | **466** | 22 compounds, 86 % | tr only | none |
| starch/SPI | 180/16 | **64** (the largest of the eight) | **1088** | 27 compounds, 52 % | **9 %** (= 92 ng/10 g, mine) | pyrazine 5, methylpyrazine 3 |
| starch/glucose/SPI | RM | 37 | — (lipid fraction 363 = 83 %) | 22 compounds, 363 ng/10 g, 83 % | none | none |
| starch/glucose/SPI | 150/20 | **31** (the smallest of the eight) | **768** | 18 compounds, 76 % | **3 %** | none |
| starch/glucose/SPI | 180/16 | 48 | **978** | 22 compounds, 69 % | **9 %** | pyrrole 1, pyrazine 1, methylpyrazine 8 |
| starch/SPI/aHVP | RM | 39 | — (lipid fraction 541 = 74 %) | 21 compounds, 541 ng/10 g, 74 % | none detected | **none detected** |
| starch/SPI/aHVP | 150/20 | 38 | **576** | 20 compounds, **444 ng/10 g**, 77 % | 4 compounds, **15 %** | **none** |
| starch/SPI/aHVP | 180/16 | 51 | **2885** | 20 compounds, **709 ng/10 g**, 25 % | 5 compounds, **1879 ng/10 g, 65 %** | 8 compounds, **94 ng/10 g, 3 %** |
| starch/glucose/SPI/aHVP | RM | 40 | **513** | 21 compounds, **406 ng/10 g**, 79 % | none | none |
| starch/glucose/SPI/aHVP | 150/20 | 34 | **4062** (the largest of the eight) | 17 compounds, **722 ng/10 g**, 18 % | 4 compounds, **3183 ng/10 g, 78 %** | 1 compound |
| starch/glucose/SPI/aHVP | 180/16 | 44 | **2840** | 20 compounds, **790 ng/10 g**, 28 % | 5 compounds, **1680 ng/10 g, 59 %** | 5 compounds, **87 ng/10 g, 3 %** |

Also printed: **94 compounds in total across the eight extrudates**; the top-four share of the
starch/SPI feedstock's volatiles is 73 % (hexanal 220, 2-pentylfuran 110, pyridine 48, ethyl acetate
36 ng/10 g), and of the starch/glucose/SPI feedstock 67 % (hexanal 150, 2-pentylfuran 74, pyridine
18 ng/10 g).

**Comparators quoted from the companion paper (Solina, Johnson & Whitfield 2007, Food Chem. 100:678-692) — ALL SECOND-HAND, primary not on disk:**

| quantity | value | attributed to |
|---|---|---|
| reducing sugar in starch extrudates, mild and extreme | **75 and 85 mg/10 g** | Solina et al. 2007 |
| starch alone, 180 C | 33 compounds, **534 ng/10 g** | Solina et al. 2007 |
| starch/glucose/aHVP | 60 compounds at 150 C, **67 at 180 C**; **54 odour points** at 180 C, at least 10 of them biscuit-like and unidentified | Solina et al. 2007 |
| odour of starch alone at 180 C | "wet paper" | Solina et al. 2007 |
| odour of starch/glucose at 180 C | "rice crispbread" | Solina et al. 2007 |
| odour of starch/aHVP at 180 C | "savoury", "malty", moderate intensity | Solina et al. 2007 |
| odour of starch/glucose/aHVP at 180 C | "bakery", "cheese cracker", moderate to strong | Solina et al. 2007 |
| aHVP's own volatiles: Maillard products | 11 % of the volatile content of the ingredient | Solina et al. 2005 |

**Sensory verdicts on the four 180 C extrudates of THIS paper** (four panellists, qualitative):
starch/SPI = "bread dough", **low**; starch/glucose/SPI = "malty and fatty", **low**;
starch/SPI/aHVP = "cheese-like", **low**; **starch/glucose/SPI/aHVP = "savoury", "melted cheese on
bread", moderate** — the only one above "low", and the paper's headline result.

### Arithmetic on the printed numbers (all mine)

**1. The Strecker totals close exactly, and they EXCLUDE benzaldehyde.**
starch/SPI/aHVP 180/16: 120 + 1420 + 99 + 240 = **1879** ✓ (printed 1879; adding benzaldehyde 29
would give 1908).
starch/glucose/SPI/aHVP 150/20: 2960 + 147 + 76 = **3183** ✓ (printed 3183).
starch/glucose/SPI/aHVP 180/16: 150 + 1200 + 110 + 220 = **1680** ✓ (printed 1680).
starch/SPI 180/16: 57 + 16 + 14 + 5 = 92; 92/1088 = **8.5 %**, printed as 9 % ✓.
starch/glucose/SPI 180/16: 46 + 33 + 7 = 86; 86/978 = **8.8 %**, printed as 9 % ✓.
starch/glucose/SPI 150/20: 17 + 3 = 20; 20/768 = **2.6 %**, printed as 3 % ✓.
**Tables 1 and 2 are internally consistent with the running text to better than a percentage
point**, with one exception noted in Flag 5 (the "15 %" for starch/SPI/aHVP at 150 C computes as
13.4 % excluding benzaldehyde or 14.8 % including it — the only place benzaldehyde appears to be
counted in a ng total).
Shares: 1879/2885 = 65.1 % ✓; 3183/4062 = 78.4 % ✓; 1680/2840 = 59.2 % ✓. Lipid shares:
444/576 = 77.1 % ✓; 709/2885 = 24.6 % ✓; 722/4062 = 17.8 % ✓; 790/2840 = 27.8 % ✓.

**2. What 1 % aHVP buys, at 180 C (the free-amino-acid effect).** Adding 1 % aHVP (= 1.84 mg free
amino acid per g of feedstock) multiplies:

| compound | starch/SPI 180/16 | starch/SPI/aHVP 180/16 | factor |
|---|---:|---:|---:|
| 3-methylbutanal | 57 | 1420 | **24.9x** |
| 2-methylbutanal | 16 | 99 | 6.2x |
| phenylacetaldehyde | 14 | 240 | **17.1x** |
| 2-methylpropanal | 5 | 120 | **24.0x** |
| total Strecker aldehydes | 92 | 1879 | **20.4x** |
| total volatiles | 1088 | 2885 | 2.7x |
| Maillard heterocycles | 8 (pyrazine + methylpyrazine) | 94 | **11.8x** |

The same comparison with glucose present (starch/glucose/SPI vs starch/glucose/SPI/aHVP at 180 C):
3-MB 46 -> 1200 = **26.1x**; total Strecker 86 -> 1680 = **19.5x**. **The aHVP effect is a factor
of ~20 on the Strecker aldehydes and is the same size with and without added glucose.**

**3. What 1 % added glucose buys: nothing, or slightly less than nothing.**

| pair, at 180 C | without glucose | with glucose | factor |
|---|---:|---:|---:|
| 3-methylbutanal, no aHVP | 57 | 46 | **0.81x** |
| 3-methylbutanal, with aHVP | 1420 | 1200 | **0.85x** |
| total Strecker, no aHVP | 92 | 86 | 0.93x |
| total Strecker, with aHVP | 1879 | 1680 | 0.89x |
| total volatiles, no aHVP | 1088 | 978 | 0.90x |
| total volatiles, with aHVP | 2885 | 2840 | 0.98x |

**Every one of these is at or below unity.** The paper's own explanation is the starch-derived
reducing sugar pool (75-85 mg/10 g, i.e. 7.5-8.5 mg/g against 10 mg/g of added glucose). Two of
the four Strecker pairs are within their combined SDs and none is a large effect; the honest
statement is **"1 % added glucose has no detectable effect on Strecker aldehyde yield in this
system"**, which is a null result of real value to a model that would otherwise make the step first
order in added sugar.

**4. The temperature inversion, and how far it goes.** Comparing 180/16 against 150/20:

| feedstock | 3-MB 150 -> 180 | total Strecker 150 -> 180 | total volatiles 150 -> 180 |
|---|---|---|---|
| starch/SPI | tr -> 57 | ~0 -> 92 | 466 -> 1088 (**2.3x up**) |
| starch/glucose/SPI | 17 -> 46 (2.7x up) | 20 -> 86 (4.3x up) | 768 -> 978 (1.3x up) |
| starch/SPI/aHVP | 44 -> 1420 (**32x up**) | ~77 -> 1879 (24x up) | 576 -> 2885 (5.0x up) |
| **starch/glucose/SPI/aHVP** | **2960 -> 1200 (2.5x DOWN)** | **3183 -> 1680 (1.9x DOWN)** | **4062 -> 2840 (1.4x DOWN)** |

Three of four go up; the fourth goes down on all three measures. The two aHVP feedstocks, which
differ only by 1 % glucose, disagree in the sign. At **150 C the glucose-containing feedstock makes
67x more 3-methylbutanal than the one without** (2960 vs 44); at **180 C it makes 0.85x** (1200 vs
1420). **So the effect of added glucose reverses with the processing condition**, from enormously
positive at 150 C / 20 % moisture to slightly negative at 180 C / 16 %. The paper reports the
numbers and does not comment on this reversal at all. It is the single most model-relevant thing in
the paper and it is confounded (Flag 2).

**5. Dimethyl disulphide does the same thing.** starch/glucose/SPI/aHVP: **100 ng/10 g at 150 C
against 47 at 180 C** (Table 2) — down, in the same feedstock and the same direction as the
Strecker aldehydes. In starch/SPI/aHVP it runs 12 -> not detected. In the two aHVP-free feedstocks
it runs up (2 -> 23 and 4 -> 30). **The same inversion, in a sulphur compound from a different
precursor, in the same two feedstocks.** That is a coherence check: whatever suppresses the volatile
yield at 180 C / 16 % in the aHVP feedstocks is not specific to the Strecker aldehydes.

**6. Molar scale of the largest yield.** 3-methylbutanal 2960 ng/10 g = 296 ng/g; at 86.13 g/mol
that is **3.4 nmol/g of extrudate**. The free amino acid delivered by 1 % aHVP is 1.84 mg/g, or
roughly **14 umol/g** at an average residue mass near 130. So the headspace-recovered 3-MB is
**~0.02 % of the total free amino acid** — and leucine is only one of the residues, so the
conversion on leucine is higher by whatever leucine's share is. This is a **recovery**, not a
yield: it is what one hour of nitrogen purge at 37 C stripped out of a 10 g slurry, and the true
extrudate content is larger by an unreported factor (Flag 3).

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Present and directly usable:
`3_methylbutanal`, `2_methylbutanal`, `2_methylpropanal`, `phenylacetaldehyde`, `benzaldehyde`,
`hexanal`, `heptanal`, `nonanal`, `e_2_octenal`, `1_hexanol`, `1_octen_3_ol`, `2_pentylfuran`,
`furfural` (Table 1's "2-Furfural"), `2_acetylfuran` (Table 2's "1-(2-Furanyl)-ethanone" is the same
compound), `dimethyl_disulfide`, `dimethyl_trisulfide`, `methylpyrazine`, `2_ethylpyrazine`,
`2_ethyl_3_5_dimethylpyrazine` (Table 2's "3-ethyl-2,5-dimethylpyrazine"), and the pair
`2_5_dimethylpyrazine` / `2_6_dimethylpyrazine`.

**Absent from the registry**, among the compounds this paper quantifies: pentanal, octanal, decanal,
(E)-2-hexenal, (E)-2-heptenal, (Z)-2-octenal, (E)-2-decenal, (E)-2-nonenal, the 2,4-octadienals,
2,4-nonadienals and both 2,4-decadienal isomers, 2-pentanone, 2-hexanone, 3-heptanone, 2-heptanone,
3-octen-2-one, 2-nonanone, 3-nonen-2-one, 3,5-octadien-2-one, 1-pentanol, 1-heptanol, 1-octanol,
1-nonanol, 2-methylfuran, 2-ethylfuran, 2,5-dimethylfuran, 2-propylfuran, 3-phenylfuran,
2,3-octanedione, **pyridine**, 2-propylpyridine, **pyrazine (the parent — the registry has the group
id `pyrazines` but no bare `pyrazine`)**, 1-H-pyrrole, 1-ethyl-1H-pyrrole, ethenylpyrazine,
2-vinyl-6-methylpyrazine, 2-(3-methylbutyl)-6-methylpyrazine, 2,4,5-trimethyloxazole, thiophene,
3-methylthiophene, 2-pentylthiophene, 2-pentylthiapyran, methyl pentyl disulphide, benzothiazole,
and the reactants (wheat starch, glucose, soy protein, the free amino acids).

**Mapping trap.** Table 2 reports **"2,5- or 2,6-Dimethylpyrazine"** unresolved at LRI 916, while
the registry carries `2_5_dimethylpyrazine` and `2_6_dimethylpyrazine` as **two separate ids**. This
row cannot be assigned to either without splitting the paper's own measurement; carry it as
ambiguous or not at all.

Every row below shares: wheat starch base + 0.05 % calcium triphosphate, APV Baker MPF 40
co-rotating twin-screw, 225 rpm, 20D, 6 mm die, **median residence time 60 s**, SME 26-27 %,
extrudate pH < 6, milled and stored at -20 C; analysed by **1 h dynamic headspace at 37 C from a
10 g / 80 mL water slurry** onto Tenax, GC-MS against 100 ng chlorotetradecane with **all response
factors assumed 1**; n = 3, SD printed.

| step | quantity | value | unit | conditions | order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| Strecker of leucine (protein-bound source) | 3-methylbutanal | **57 (5.0)** | ng/10 g | starch + 1 % SPI, 180 C / 16 % | **none — one residence time** | Table 1 | measured_level (headspace recovery) |
| " | " | tr (0.5) | ng/10 g | same, 150 C / 20 % | — | Table 1 | level_only (< 0.5) |
| " | " | **1420 (210)** | ng/10 g | starch + 1 % SPI + 1 % aHVP, 180 C / 16 % | — | Table 2 | measured_level |
| " | " | **44 (11)** | ng/10 g | same, 150 C / 20 % | — | Table 2 | measured_level |
| " | " | **1200 (320)** | ng/10 g | starch + 1 % glucose + 1 % SPI + 1 % aHVP, 180 C / 16 % | — | Table 2 | measured_level |
| " | " | **2960 (310)** | ng/10 g | same, 150 C / 20 % — **the largest single value in the paper, and it is at the LOWER temperature** | — | Table 2 | measured_level |
| " | " | 46 (13) / 17 (0.9) | ng/10 g | starch + 1 % glucose + 1 % SPI, 180/16 and 150/20 | — | Table 1 | measured_level |
| Strecker of isoleucine | 2-methylbutanal | 16 (0.2) / 33 (5.3) / 99 (1.4) / 17 (4.1) / 110 (30) / 147 (40) | ng/10 g | the six extrudates in the order S/SPI 180, S/Glc/SPI 180, S/SPI/aHVP 180, S/SPI/aHVP 150, S/Glc/SPI/aHVP 180, S/Glc/SPI/aHVP 150 | — | Tables 1, 2 | measured_level |
| Strecker of valine | 2-methylpropanal | 5 (1.2) / 7 (1.3) / 120 (8.2) / 150 (40); **not detected at 150 C in any feedstock** | ng/10 g | as above | — | Tables 1, 2 | measured_level + null_result at 150 C |
| Strecker of phenylalanine | phenylacetaldehyde | 14 (0.7) / **not detected** / 240 (63) / 16 (0.8) / 220 (31) / 76 (5.9) | ng/10 g | S/SPI 180, S/Glc/SPI 180, S/SPI/aHVP 180, S/SPI/aHVP 150, S/Glc/SPI/aHVP 180, S/Glc/SPI/aHVP 150 | — | Tables 1, 2 | measured_level; **the S/Glc/SPI 180 non-detection against 14 ng/10 g without glucose is a printed anomaly** |
| **free-amino-acid effect** | Strecker aldehyde multiplier on adding 1 % aHVP (= 1.84 mg free AA/g) | **20.4x** (no glucose) / **19.5x** (with glucose), at 180 C | — | as above | — | derived (mine) from Tables 1, 2 | within_study_ratio (**the most transportable number here**) |
| **added-glucose effect** | Strecker aldehyde multiplier on adding 1 % glucose at 180 C | **0.93x** (no aHVP) / **0.89x** (with aHVP) | — | as above | — | derived (mine) | **null_result — added reducing sugar at 1 % does not raise the Strecker yield** |
| " | the same multiplier at 150 C / 20 %, with aHVP | **67x** (3-MB: 44 -> 2960) | — | as above | — | derived (mine) | within_study_ratio — **the glucose effect reverses between the two conditions** |
| **temperature/moisture effect** | 3-MB, 150/20 -> 180/16 | **32x up** (S/SPI/aHVP) but **2.5x DOWN** (S/Glc/SPI/aHVP) | — | as above | — | derived (mine) | within_study_ratio — **NOT a barrier; T and moisture move together (Flag 2)** |
| pyrazine formation | total Maillard heterocycles | **94 (S/SPI/aHVP) and 87 (S/Glc/SPI/aHVP)** ng/10 g, both **3 %** of volatiles, both **only at 180 C** | ng/10 g | as above; **zero at 150 C in both aHVP feedstocks** | — | text pp. 1533, and Table 2 | measured_level + **null_result at 150 C** |
| " | methylpyrazine | 3 (0.6) / 8 (2.6) / 28 (6.1) / 34 (10) | ng/10 g | S/SPI 180, S/Glc/SPI 180, S/SPI/aHVP 180, S/Glc/SPI/aHVP 180 | — | Tables 1, 2 | measured_level |
| " | ethylpyrazine | 19 (4.3) / 22 (4.9) | ng/10 g | the two aHVP feedstocks at 180 C | — | Table 2 | measured_level |
| " | 3-ethyl-2,5-dimethylpyrazine | 11 (2.7); **not detected with glucose** | ng/10 g | S/SPI/aHVP 180 | — | Table 2 | measured_level |
| sulphur route | dimethyl disulphide | 23 (3.3) / 30 (3.1) / **100 (15) at 150 C vs 47 (2.6) at 180 C** | ng/10 g | S/SPI 180, S/Glc/SPI 180, S/Glc/SPI/aHVP at the two conditions | — | Tables 1, 2 | measured_level (**the same temperature inversion, section 3 arithmetic 5**) |
| " | dimethyl trisulphide | 3 (0.1) / 1 (0.1) / 16 (4.2) / 16 (3.4) | ng/10 g | S/SPI 180, S/Glc/SPI 180, and both aHVP feedstocks at 180 C | — | Tables 1, 2 | measured_level |
| sugar-fragment marker | 2,3-octanedione; 2-furfural; 2-acetylfuran | 7 / 5; 2; 3 / 3 | ng/10 g | 180 C only, various | — | Tables 1, 2 | measured_level (**the only sugar-derived rows in either table — the sugar side of this system is almost invisible to the method**) |
| lipid background | hexanal | 170-310 across all eight extrudates and all four feedstocks | ng/10 g | all | — | Tables 1, 2 | measured_level (**effectively constant — a useful internal normaliser**) |
| reducing sugar produced by starch cleavage | — | **75 and 85** | mg/10 g | mild and extreme extrusion of wheat starch | — | text p. 1523 | **second-hand from Solina et al. 2007 (Food Chem. 100:678-692), NOT on disk** |
| extrudate pH | — | **< 6** | — | all extrudates | — | text p. 1537 | level_only (a bound, no value, no method) |
| odour thresholds | which compounds exceed their OTC in water | 5-14 compounds per extrudate, listed by name | — | **thresholds themselves never printed** | — | text, citing Badings 1970 and Fors 1983 | **pointer only — no threshold value appears anywhere in this paper** |
| odour intensity ranking at three dilutions | — | Tables 3 and 4 | VS/S/M/W | the two 180 C extrudates, 5 assessors | — | Tables 3, 4 | ordinal_only |
| any rate constant, activation energy, time course, or residence-time series | — | **not printed** | — | — | — | whole paper | — |

### Can anything here be put on the trunk's axes? Step by step.

**(a) No rate, no barrier, and the two temperatures cannot be made into one.** Two conditions, each
with its own moisture; one residence time; no intermediate sampling. Even if moisture had been held
constant, a two-point Arrhenius on a headspace recovery in a shearing extruder would not be a step
barrier. **Nothing from this paper enters a rate registry.**

**(b) The amino-acid supply ratio is the transportable number.** The 20x rise in Strecker aldehydes
for 1.84 mg/g of added free amino acid, reproduced twice (with and without glucose) at the same
temperature, is a within-study ratio measured on one method in one matrix. It supports the
per-amino-acid supply term the identity wave drafts in
`results/validation/kinetic_core_b19_prereg_draft.md`, and it does so in a low-moisture protein
matrix rather than a dilute aqueous pot. It is a **ratio of levels, not of rates**, and it is net of
every sink and of the headspace partition.

**(c) The added-glucose null is the second transportable result, and it is a falsifier.** Four
independent comparisons at 180 C all give a multiplier at or below 1.0 for 1 % added glucose. Any
model that makes the Strecker step first order in added reducing sugar predicts a rise here. The
escape is that starch cleavage already supplies 7.5-8.5 mg/g — but that figure is second-hand from a
paper not on disk, so **the escape itself is unverified in this corpus** (Flag 8).

**(d) Nothing can be used as a benchmark concentration.** Every number is ng recovered per 10 g of
sample by a one-hour nitrogen purge at 37 C, with response factors assumed to be 1 and with the
chlorododecane recovery never reported. These are **not concentrations in the extrudate** and the
factor between them and the true content is unreported and compound-dependent (it depends on each
compound's volatility and on its binding to the SPI protein, which the discussion itself raises as a
mechanism). Class every level as a headspace recovery.

**(e) The Maillard heterocycles are a clean on/off result.** In both aHVP feedstocks, pyrroles,
pyrazines and oxazoles are **entirely absent at 150 C / 20 % and present at 3 % of volatiles at
180 C / 16 %** — 8 compounds and 5 compounds respectively. Against them, the Strecker aldehydes are
already at 15 % and 78 % of volatiles at 150 C. **The heterocycle branch has a much sharper
condition threshold than the Strecker branch in the same pot**, which is the qualitative shape a
higher-barrier condensation step would produce and is consistent with the pyrazine lane's measured
barriers (103.1 and 114.9 kJ/mol, `FROZEN_B18`) sitting above the Strecker-cascade barriers
recorded in `chan1994b_extraction.md` (80-90 kJ/mol) and `cremer2000_extraction.md` (115-124
kJ/mol at aw 0.52). It is a sign check, not a measurement.

**(f) The SPI suppression result is real but unexplained.** Comparing against the companion paper,
adding 1 % SPI to starch/glucose/aHVP cuts the compound count from 67 to 44 at 180 C and the odour
points from 54 to 29, with **eight of the missing compounds being Maillard heterocycles** and at
least ten of the missing odour points "biscuit-like". The authors offer three candidate mechanisms
and choose none: volatile binding to the added protein; interaction of SPI with Maillard precursors
or products; and carbonyls from oxidation of the SPI's 3.5 % lipid interfering with Maillard steps
(citing `whitfield1988_extraction.md`). **All of the comparator numbers in this paragraph are
second-hand.**

## 5. Flags

1. **Table extraction and printing defects.** The two large tables are landscape six-column layouts
   and the layout extractor packs value/SD pairs raggedly; each row was assigned by counting pairs
   against the printed column heads. **One row did not resolve: 2-ethylfuran's RM and 150/20 entries
   in the starch/glucose/SPI block, marked "not resolved" and not guessed.** Two defects are in the
   print itself, not the OCR: nonanal (S/SPI/aHVP 150/20) prints `3 (0.4` with no closing
   parenthesis, and methylpyrazine (S/Glc/SPI/aHVP 150/20) prints a value of 3 with **a dash where
   the SD should be**. The Methods print "stored at 20 °C" and "cooling the pre-column to 78 °C";
   both are minus signs lost in typesetting and are read as **-20 C** and **-78 C** (a Tenax
   pre-column focus at +78 C and frozen storage at +20 C are both impossible).
2. **TEMPERATURE AND MOISTURE ARE CONFOUNDED, AND SO IS EVERYTHING DOWNSTREAM OF THEM.** The only
   two conditions are 150 C / 20 % and 180 C / 16 %. There is no 150 C / 16 % and no 180 C / 20 %.
   **No temperature coefficient, no activation energy, and no moisture coefficient can be extracted
   from this paper, in either direction.** Die pressure also varies by a factor of 2.4-6 across the
   runs (240 psi down to 40-100 psi) and is a third uncontrolled variable. The striking inversion in
   the starch/glucose/SPI/aHVP feedstock (3-MB 2960 -> 1200) is therefore a **condition** effect,
   not a temperature effect, and must never be quoted as evidence about temperature alone.
3. **Every level is a dynamic-headspace recovery, and the recovery is never reported.** 1 h of N2 at
   40 ml/min through a stirred 10 g / 80 mL slurry at 37 C, onto 10 mg of Tenax. Chlorododecane was
   added to the flask expressly "to estimate the recovery of the volatile compounds" and **that
   estimate appears nowhere in the paper**. Response factors are assumed to be 1 by the authors' own
   footnote, which they flag as making the quantities "approximate". So the ng/10 g values are
   **neither extrudate concentrations nor calibrated amounts**; they are internal-standard-normalised
   trap loadings. Comparisons within a compound across feedstocks are sound; comparisons between
   compounds are not.
4. **Both key ingredient characterisations and every cross-feedstock comparator are second-hand.**
   The starch composition and the 75/85 mg/10 g reducing sugar figures come from Solina, Johnson &
   Whitfield 2007 (Food Chem. 100:678-692); the SPI and aHVP compositions from Solina, Baumgartner,
   Johnson & Whitfield 2005 (Food Chem. 90:861-873). **Neither is on disk.** So the free-amino-acid
   loading (18.4 %), the protein content (82.5 %), the absence of free amino acids in SPI, the
   absence of reducing sugars in all three ingredients, and the starch-derived sugar pool that
   explains the added-glucose null are **all** second-hand within this corpus. Fetching both is the
   first requirement before any number here is used quantitatively.
5. **The treatment of benzaldehyde in the "Strecker aldehyde" totals is inconsistent.** Footnote f
   of both tables warns that benzaldehyde "may be either lipid- or amino acid-derived", and the text
   says it "has been included among these aldehydes" for the **compound counts** (4 and 5). But the
   printed **ng totals** (1879, 3183, 1680) all reproduce exactly *without* benzaldehyde, while the
   one printed **percentage** for starch/SPI/aHVP at 150 C ("15 %") reproduces only *with* it
   (14.8 % vs 13.4 %). Recompute any Strecker total from the table rows rather than relying on the
   text's aggregates, and state the benzaldehyde convention explicitly.
6. **Calcium triphosphate at 0.05 % w/w is in every feedstock and is never discussed.** It is added
   as a flow aid, but phosphate is a well-documented Maillard catalyst (the same claim appears in
   `chan1994b_extraction.md` and in Potman & van Wijk, cited there). At 500 mg/kg it is a constant
   across the eight extrudates, so it does not confound the between-feedstock comparisons, but it
   means the absolute levels are those of a phosphate-catalysed system and should not be compared
   with unbuffered pots.
7. **A printed value in the text disagrees with the table.** Section 3.4 lists the greatest
   components of the starch/glucose/SPI/aHVP extrudates as including "dimethyl disulphide (**160**
   and 47 ng/10 g)", while **Table 2 prints 100 (15) and 47 (2.6)** for the same cells. The 47 agrees;
   the 150 C value does not. **Use the table (100), and record the discrepancy.** This is the only
   text-vs-table conflict found; the Strecker aggregates all check out (section 3 arithmetic 1).
8. **The explanation for the added-glucose null rests on an off-disk number.** The whole resolution
   of the authors' own "perplexing" observation is that starch cleavage already delivers 75-85 mg
   reducing sugar per 10 g, i.e. roughly the same as the 10 mg/g added. That figure is measured in
   the companion paper, not here, and the companion is not on disk. **Until it is fetched, the null
   is an observation without a verified mechanism.**
9. **No pH value, no water activity, no temperature profile through the barrel.** pH is given once
   as "< 6" with no method and no per-feedstock number, in a discussion aside about whether the
   proline route to roast-smelling pyrrolines (citing Hofmann & Schieberle 1998,
   `hofmann1998b_extraction.md`, which works at pH 7) can operate here. Water activity is never
   measured — only feed moisture, which is not the same thing. The six barrel probes are described
   but no profile is printed.
10. **Identification strength is uneven and is marked in the tables.** 10 rows in Table 1 and 8 in
    Table 2 carry **"MS" only** (mass spectrum against a library, no retention-index confirmation),
    including both 2,4-decadienal isomer rows, 3-phenylfuran, 2-acetylfuran and several unknowns.
    Table 2's LRI-916 pyrazine is printed as **"2,5- or 2,6-dimethylpyrazine"** — unresolved between
    two compounds the registry keys separately. The GCO tables carry an explicit footnote that the
    "major compound in region of odour" **need not be the compound responsible for the odour**, and
    the text says so directly for LRI 771 ("the relatively weak odorant 1-H-pyrrole could not be
    responsible for the 'strong' odour at LRI 771") and for the very-strong unknown at LRI 682 and
    the strong unknown at LRI 1115. **Do not read Tables 3 and 4 as compound-to-odour assignments.**
11. **Odour thresholds are used throughout and never printed.** The phrase "exceeded their OTC in
    water" carries every sensory argument in the paper, citing Badings 1970 and Fors 1983, and **not
    one threshold value appears anywhere.** The authors themselves warn that "the OTC values in the
    extrudate may be very different from those in water, due to the different nature of the matrix".
    No threshold can be taken from this paper.
12. **What this paper does not contain**: any rate constant; any activation energy; any time course
    or residence-time series; any second residence time; any temperature at constant moisture; any
    water activity; any pH value; any threshold; any recovery figure; any quantification of the
    proteins' released free amino acids (the central premise of the whole interpretation is that
    extrusion liberates them, and **it is never measured**); any figure of any kind; any
    supplementary material.
13. **What to fetch next, in priority order**: (i) **Solina, Johnson & Whitfield 2007, Food Chem.
    100:678-692** — the companion, which carries the starch-alone, starch/glucose, starch/aHVP and
    starch/glucose/aHVP arms of the same experiment plus the reducing-sugar measurements, and
    without which half of this paper's comparisons cannot be checked; (ii) **Solina, Baumgartner,
    Johnson & Whitfield 2005, Food Chem. 90:861-873** — the SPI and aHVP compositions, including the
    18.4 % free amino acid figure and the amino acid profile that would let the leucine-specific
    conversion be computed; (iii) Fors 1983, ACS Symp. Ser. 215, pp. 185-286, and Badings 1970 —
    the two threshold compilations behind every OTC claim here; (iv) Bruechert et al. 1988, J. Food
    Sci. 53:1444-1447, on lipid contribution to volatiles in extruded corn model systems.
14. **Registry gaps against `data/keys/compounds.yml`**: the Strecker aldehydes, the pyrazines this
    paper quantifies, both dimethyl sulphides, 2-pentylfuran, 2-acetylfuran, furfural, hexanal,
    heptanal, nonanal, (E)-2-octenal, 1-hexanol and 1-octen-3-ol are all keyed. **Absent and needed
    before the lipid background or the heterocycle rows could be written**: pentanal, octanal,
    decanal, the four 2-alkenals other than (E)-2-octenal, both 2,4-decadienal isomers, the
    2-alkanones, **pyridine** (which reaches 88 ng/10 g in one extrudate — the largest unkeyed
    Maillard-relevant value in the paper), **pyrazine (the parent; only the group id `pyrazines`
    exists)**, 1-H-pyrrole, 1-ethyl-1H-pyrrole, ethenylpyrazine, 2-vinyl-6-methylpyrazine,
    2-(3-methylbutyl)-6-methylpyrazine, 2,4,5-trimethyloxazole, 2,3-octanedione, and the thiophene
    and thiapyran series. **And the LRI-916 row cannot be mapped at all** while the registry splits
    2,5- and 2,6-dimethylpyrazine that the paper does not resolve.
