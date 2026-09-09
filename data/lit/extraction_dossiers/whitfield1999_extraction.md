# Whitfield & Mottram 1999 — EXTRACTION (fed norfuraneol 50 mmol/L with cysteine 50 mmol/L or H2S ~97 mmol/L, 0.5 M phosphate pH 4.5, flame-sealed 5 mL ampules, 140 C for 60 min; 69 compounds by dynamic headspace GC/MS against methyl decanoate, duplicates printed, in µg per 10 mg of norfuraneol fed)

### THE PARENT OF THE pH-6.5 PAPER, and the compound-level table its dossier could not verify — now verified row by row. What is new here and bears on the thiol sink: **26 disulfides are measured in this pot, they carry about 7.5 % of the whole thiol pool and about 35 % of the 2-methyl-3-furanthiol specifically (both mine, from Table 1)** — a third laboratory, a third temperature and a third pH agreeing with Zhou 2023's 6.5-9.6 % and Zhang 2024's 8.7 %, against the 0.03-0.93 % the refused B17 runs produced. And the paper **names the oxidant**: not air, not the analysis, but the pot's own α-dicarbonyls, in an explicit redox couple (Figure 6) that is the same reduction step which makes the mercaptoketones.

**Source on disk:** `data/articles/whitfield1999.pdf` (9 pp., J. Agric. Food Chem. 1999, 47 (4),
1626-1634). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/whitfield1999.txt`, 610 lines). **Table 1 — the paper's only table, 69 rows,
four numeric columns — came through clean** and is re-typed in full below. The four numeric columns
are cysteine duplicate 1, cysteine duplicate 2, H2S duplicate 1, H2S duplicate 2, and **that reading
is verified three independent ways** against the class totals printed in the companion paper's
Table 2 (section 3, "column check"). Figures 1-12 are reaction schemes and carry no data. There is
no supplementary material. Repo status before this dossier: Whitfield 1999 is a declared **FIT**
source (`docs/reference/FIT_HOLDOUT_DECLARATION.md` line 41, "the only independent replication of
the Hofmann NF chemistry; needed to pin the NF channel"), supplies **three scored rows** in every
B2.x fit generator (`whitfield_nf_cys_MFT` 0.150 mol %, `whitfield_nf_h2s_MFT` 0.120 mol %,
`whitfield_mercaptoketone_over_MFT` 16.3), and is quoted in the `r_nf_mp3p` note in
`src/kinetic_core/sulfur.py`. **It has never had a dossier of its own**, and the companion dossier
`whitfield2001_extraction.md` says so in its own caveat 7: the pH-4.5 compound levels there were
"taken from the repo's k3 inventory, not from this PDF". This dossier closes that.

## 0. Identity

| field | value |
|---|---|
| Title | "Investigation of the Reaction between 4-Hydroxy-5-methyl-3(2H)-furanone and Cysteine or Hydrogen Sulfide at pH 4.5" |
| Authors | Frank B. Whitfield (Food Science Australia, North Ryde NSW) and Donald S. Mottram (corresponding; Food Science and Technology, University of Reading) |
| Venue | J. Agric. Food Chem. 1999, 47 (4), 1626-1634. Received 1 September 1998, accepted 12 January 1999, web 5 March 1999 |
| DOI / article ID | 10.1021/jf980980v (printed as `JF980980V`) |
| Naming | the paper's "HMF" is **4-hydroxy-5-methyl-3(2H)-furanone = norfuraneol = the repository's `NF`**, NOT hydroxymethylfurfural |
| Why pH 4.5 | "this pH has been shown to favor the production of sulfur compounds with meatlike aromas (Madruga and Mottram, 1995)" and it is "slightly on the acidic side of the isoelectric point of cysteine (pH 5.1), a condition that is known to favor the formation of sulfur-substituted furans" |
| Companions on disk | `whitfield2001_extraction.md` — **the same experiment at pH 6.5**, same laboratory, same charges, same instrument, same reporting unit; `whitfield1988_extraction.md`; this cluster's `cerny2003_extraction.md` (which uses fed norfuraneol as a labelling probe) and `hofmann1995_extraction.md` |

## 1. Why it matters

**What it contributes to the thiol-sink question. A number, and it is a big one.**

The three sink refusals turned on a single quantitative failure that all three shared: **T3, the
disulfide share.** The model makes 0.03-0.93 % of its thiol as disulfide; Zhou 2023 measures 6.5 to
9.6 % at 120 C and pH 6-8, and Zhang 2024 measures 8.7 % at 115 C. B17 variant (b)'s outcome note
concluded the channel is "oxidant-limited" and that a larger reservoir or a second oxidant would be
needed. **This paper is a third measurement of the same share, at 140 C and pH 4.5, in a pot that is
already a FIT source, and it agrees with the other two.** From Table 1 (my arithmetic, section 3
item 2):

- **Disulfide-bound thiol is ~7.5 % of the total thiol in the cysteine system** (24 µg of disulfides
  against 294 µg of free thiols and mercaptoketones, per the class totals the companion paper prints
  for exactly these data), and **~29 % in the H2S system** (18 against 44).
- **For 2-methyl-3-furanthiol specifically the share is far higher: ~35 %.** Free MFT is 15 µg per
  10 mg of norfuraneol; the six MFT-bearing disulfides in the cysteine column (compounds 47, 48, 54,
  55, 62, 64) carry a further **~8.0 µg of MFT equivalent (mine)**. So the MFT actually made in this
  pot is **0.230 mol % of the norfuraneol fed, not 0.150** — and **the `whitfield_nf_cys_MFT` FIT
  row, which targets 0.150 mol %, is scoring free MFT against a model species that in the network is
  the whole MFT pool.** That is a 1.5x systematic, in the direction of making the model's MFT look
  too high, on one of the three rows this paper supplies.

**And the paper names the oxidant, which is precisely what the B17 (b) refusal said was missing.**
p. 1631, in full: "the conditions in the reaction systems were not conducive to aerial oxidation,
because of the relatively high concentrations of hydrogen sulfide. Furthermore, recent observations
on the stability of thiols and disulfides during headspace collection and thermal desorption have
shown that thiols are not converted to disulfides during the analysis procedures employed in this
work (Mottram et al., 1998). **Therefore, another redox system must be involved in the formation of
the disulfides.** In the formation of mercaptoalkanones, it has been proposed that α-dicarbonyls are
reduced to hydroxyalkanones before substitution by hydrogen sulfide (Figures 1 and 2). This could
provide the redox system required for the formation of the disulfides (Figure 6)." Figure 6 is
captioned "Proposed redox reaction between thiols and α-dicarbonyl compounds."

**This is an internal, sugar-derived, continuously regenerated oxidant that the module does not
have.** `THIOL_CHANNELS`'s dimerisation entry explains its own gating by Ngamchuea's result that
metal-free thiol autoxidation in water is negligible and says the channel is therefore first order
in an explicit `OX` pool — a pool the B11 wave ships **inert**, sized as an ambient reservoir. The
refusal's own diagnosis was that this reservoir "runs out". Whitfield and Mottram's proposal makes
the oxidant **stoichiometrically coupled to the mercaptoketone-forming flux**: every α-dicarbonyl
that is reduced to a hydroxyalkanone on the way to a mercaptoketone is an oxidising equivalent
delivered to a thiol. In this pot the mercaptoketone flux is enormous (244 µg per 10 mg NF against
15 µg of free MFT), so the oxidant supply scales with the very chemistry the lane already models.
**A variant (c) or (d) that sources `OX` from the α-dicarbonyl reduction flux rather than from a
fixed ambient reservoir is the structure this paper argues for, and it was never pre-registered.**
It is a claim about mechanism, made by the authors, not a measured rate — that boundary is kept.

**What else this paper adds over its pH-6.5 companion.** (a) The **compound-level pH-4.5 column
itself**, previously reaching the code only through the k3 inventory; every number the generators use
is now checked against the print (section 3 item 1) and all four verify. (b) A **26-disulfide
inventory**, 20 of them not previously reported in any food or model system, with the specific
finding that 3-thienyl-containing disulfides occur **only** in the cysteine system while
1-mercaptobutan-2-one-derived disulfides occur **only** in the H2S system. (c) The **explicit charge
of the H2S ampule**, which resolves the printed defect the companion dossier flags: the 2001 paper's
Methods sentence omits norfuraneol from the H2S ampule, whereas this paper prints "mixing equal
quantities (1 mL each) of the solutions of HMF (containing 11.4 mg) with those of either cysteine
(12.1 mg) **or the hydrogen sulfide (~6.6 mg)**". Since the 2001 paper says its method was "as
described previously (8)", the 2001 H2S charge is now grounded in a printed parent method rather
than in an analogy. (d) A **negative that constrains the mechanism**: none of the di- or
tetrahydro derivatives of MFT or of 2-methyl-3-thiophenethiol was detected, against van den Ouweland
& Peer's report of them under a 100-fold excess of H2S. (e) The observation that **MFT and
2-methyl-3-thiophenethiol come out at similar levels from cysteine and from H2S alone** (15/15
against 16/8; 10/5 against 9/-), "which suggests that the pathways involved require only hydrogen
sulfide and not other cysteine degradation products", while **3-thiophenethiol is found only in the
cysteine system** — the same conclusion Cerny 2003 reaches by isotope (95 % unlabelled, i.e. from
cysteine).

What this paper does NOT give the repository: any rate; any time course (one time point, 60 min);
any barrier; any second temperature; any pH other than 4.5 (that is the companion's job); any
absolute concentration in the solution (everything is a headspace-recoverable amount with response
factors assumed 1); any total (the class totals come from the companion's Table 2); and any
quantification of norfuraneol consumed.

## 2. Methods as they matter to a model

- **Charges, exactly as printed.** Separate 0.1 M solutions of **HMF (norfuraneol)** and of
  **cysteine** in 0.5 M phosphate buffer at pH 4.5. A **saturated H2S solution (~0.2 M)** made by
  passing H2S through the same buffer cooled in ice at 0 C. "Reaction mixtures were prepared by
  mixing equal quantities (**1 mL each**) of the solutions of HMF (containing **11.4 mg**) with
  those of either cysteine (**12.1 mg**) or the hydrogen sulfide (**~6.6 mg**) in **5 mL glass
  ampules**." So, in the final 2 mL: **norfuraneol 50 mmol/L** (11.4 mg / 114.10 = 0.0999 mmol),
  **cysteine 50 mmol/L** (12.1 mg / 121.16 = 0.0999 mmol), or **H2S ~97 mmol/L** (6.6 mg / 34.08 =
  0.194 mmol) — an **H2S : NF ratio of 1.94 : 1 (mine)**.
- **Heating.** Ampules **flame sealed**, then **140 C for 60 min** in an oven. No stirring during
  reaction. The headspace is the sealed ampule's own air over 2 mL in a 5 mL ampule.
- **A second, larger run.** "The reaction of HMF with cysteine was also carried out at a higher
  concentration (**2-fold**) and in **greater quantities (10 mL each)** to obtain larger amounts of
  some trace volatile components." Table 1's footnote a marks the entries found only in that run as
  **tr\*** and describes it as "reaction mixtures with **10 times higher quantities** of cysteine and
  HMF" — the two descriptions (2-fold concentration in 10 mL rather than 1 mL, i.e. 20x the
  material) are not the same statement (Flags 4).
- **Isolation.** After cooling, transferred to a 250 mL conical flask with **20 mL of 0.5 M
  phosphate pH 4.5**, tube rinsed twice with 2 mL; **methyl decanoate, 100 µg in 0.1 mL ethanol**,
  added as internal standard. **Dynamic headspace**: the diluted, stirred solution held at **60 C**,
  volatiles swept onto a room-temperature Tenax GC trap (115 mm x 0.75 mm i.d.) with **oxygen-free
  nitrogen at 60 mL/min for 1 h**, then 5 min of dry nitrogen.
- **GC/MS.** Varian 1440 with a Unijector in concentrator-headspace mode, **BP5 50 m x 0.32 mm**,
  helium 1 mL/min. Thermal desorption **260 C for 5 min** with the oven cryofocused at **0 C** under
  liquid nitrogen, then 60 C for 5 min and up to 250 C at **4 C/min**. Varian-MAT 311A
  double-focusing MS, **EI 70 eV**, source 250 C, continuous scan **34-340 amu at 2 s/decade**.
  C6-C20 n-alkanes as external LRI standards.
- **Quantification — RELATIVE, and the paper says so.** "the approximate concentrations of selected
  compounds were determined by comparing their GC/MS chromatogram peak areas with the area of the
  internal standard, methyl decanoate, which was taken as 100 µg, and **assuming all response factors
  were 1**. The concentrations of these compounds are reported as **micrograms per 10 mg of HMF
  used** in the reaction." **Detection limit 0.1 µg/10 mg HMF** (3x background noise); "trace" (tr)
  = **< 1 µg/10 mg HMF**. Note this differs from the companion paper, which defines tr as **0.1 to
  1** µg/10 mg — the same band, stated as an interval there and as an upper bound here.
- **Replication.** "duplicate analyses are shown" — both values are printed for every entry, which
  is a great deal more than most papers in this corpus give, and the spread is often wide (Flags 3).
- **The conversion the repository uses.** 10 mg norfuraneol = **87.64 µmol**, so
  mol % = µg / (MW x 0.8764). For MFT (MW 114.17) the divisor is **100.06**, i.e. **µg per 10 mg NF
  is numerically the mol % x 100** to within 0.1 % — 15 µg = **0.150 mol %**.
- **Sensory.** Three assessors from the laboratory. Cysteine system: "sulfurous, rubbery, and boiled
  meat", and "when diluted, the mixture had increased meatlike characteristics". H2S system:
  "dominated by the smell of hydrogen sulfide when the reaction vial was first opened, but as this
  odor disappeared, a caramel, meatlike aroma developed."
- **What the buffer does.** 0.5 M phosphate at 50 mM reactants is a genuine buffer here (10:1), so
  unlike most pots in this corpus the pH 4.5 is likely to have held. The paper does not re-measure it.

## 3. Tables re-typed

### Table 1. "Volatile Compounds Obtained from Reactions between 4-Hydroxy-5-methyl-3(2H)-furanone and Cysteine or Hydrogen Sulfide"

Column head as printed: **"approx concn^a"** over **cysteine** and **H2S**, each with its two
duplicate values; then **"method of ID^b"** and **"LRI^c"**. The mass-spectral column (footnote d) is
omitted here except where noted. Footnote a: "Concentrations (µg/10 mg of HMF) obtained by comparing
GC/MS peak area with that from 100 µg of methyl decanoate internal standard added to the HMF
solution before volatile collection; duplicate analyses are shown; **-, not detected (limit of
detection ~0.1 µg/10 mg of HMF); tr, < 1 µg/10 mg of HMF; tr\*, found only in reaction mixtures with
10 times higher quantities of cysteine and HMF.**" Footnote b: "MS + LRI, identified by comparison of
mass spectrum and LRI with those of authentic compound; MS, tentative identification by comparison
with mass spectrum reported in the literature; ms, tentative identification by interpretation of mass
spectrum." Footnote e: "Not reported previously in meat or meatlike model systems." Footnote f:
"Pairs of diastereoisomers."

| no. | compound | cys 1 | cys 2 | H2S 1 | H2S 2 | ID | LRI |
|---:|---|---:|---:|---:|---:|---|---:|
| 1 | 2,3-pentanedione | 4 | 15 | - | - | MS+LRI | 680 |
| 2 | 2,4-pentanedione | 12 | 16 | 47 | 31 | MS+LRI | 787 |
| 3 | 3,4-hexanedione | 4 | 4 | - | - | MS+LRI | 800 |
| 4 | 4,5-dihydro-5-methyl-3(2H)-furanone | 2 | 8 | 6 | 6 | MS | 808 |
| 5 | **3-mercaptobutan-2-one**^e | **84** | **100** | **6** | **5** | MS+LRI | 817 |
| 6 | **2-methyl-3-furanthiol** | **15** | **15** | **16** | **8** | MS+LRI | 867 |
| 7 | 1-mercaptobutan-2-one | tr | - | 8 | 6 | ms | 886 |
| 8 | **3-mercaptopentan-2-one** | **81** | **68** | **9** | **6** | MS+LRI | 902 |
| 9 | **2-mercaptopentan-3-one** | **78** | **77** | **8** | **6** | MS+LRI | 908 |
| 10 | 3-thiophenethiol | 26 | 14 | - | - | MS | 972 |
| 11 | 2-acetyl-5-methylfuran | tr | - | 6 | 6 | MS+LRI | 977 |
| 12 | dihydro-5-methylthiophen-3(2H)-one | 3 | 2 | 5 | 3 | MS | 982 |
| 13 | dihydro-2-methylthiophen-3(2H)-one | 66 | 64 | 26 | 33 | MS+LRI | 990 |
| 14 | 5-mercaptohexan-2-one^e | 10 | 5 | - | - | ms | 993 |
| 15 | 4,5-dihydro-2,4-dimethylthiophen-3(2H)-one (E or Z) | 22 | 17 | 2 | 3 | MS | 1016 |
| 16 | 4,5-dihydro-2,4-dimethylthiophen-3(2H)-one (E or Z) | 5 | 5 | - | - | MS | 1027 |
| 17 | 2-methyl-3-thiophenethiol | 10 | 5 | 9 | - | MS | 1060 |
| 18 | **3-methyl-1,2-dithiolan-4-one** | 2 | 1 | **204** | **196** | MS | 1071 |
| 19 | 2-ethyl-4,5-dihydrothiophen-3(2H)-one | 3 | tr | 6 | - | MS | 1082 |
| 20 | 2-acetylthiophene | tr | tr | - | tr | MS+LRI | 1092 |
| 21 | 3,5-dimethyl-1,2-dithiolan-4-one (E or Z) | 16 | 10 | 112 | 116 | MS | 1098 |
| 22 | 3,5-dimethyl-1,2-dithiolan-4-one (E or Z) | 12 | 8 | 77 | 86 | MS | 1105 |
| 23 | 2-formyl-5-methylthiophene | 7 | tr | - | tr | MS+LRI | 1124 |
| 24 | (3-thienyl)-2-propanone | 4 | tr | - | - | MS | 1134 |
| 25 | 2-methyl-(3-methylthio)thiophene | tr | tr | - | - | MS | 1141 |
| 26 | 2-acetyl-5-methylthiophene | 7 | 10 | 4 | 3 | MS+LRI | 1157 |
| 27 | 3-ethyl-1,2-dithiolan-4-one | 15 | 6 | 48 | 30 | MS | 1167 |
| 28 | 1-(3-thienyl)-1-propanone | 6 | 6 | 12 | 1 | MS+LRI | 1183 |
| 29 | 3-ethyl-5-methyl-1,2-dithiolan-4-one (E or Z) | tr | tr | 7 | tr | MS | 1193 |
| 30 | 3-ethyl-5-methyl-1,2-dithiolan-4-one (E or Z) | tr | tr | 3 | 2 | MS | 1197 |
| 31 | 2,3-dihydro-6-methylthieno[2,3-c]furan | tr | - | tr | 4 | MS+LRI | 1199 |
| 32 | 3-ethyl-2-formylthiophene | 57 | 17 | - | - | MS | 1206 |
| 33 | thieno[3,2-b]thiophene | tr | - | - | - | MS | 1213 |
| 34 | 3-methyl-1,2-dithian-4-one | 3 | 3 | - | 5 | MS | 1220 |
| 35 | diformylthiophene^e | 8 | 3 | - | - | MS | 1245 |
| 36 | 3,5-dimethyl-1,2-dithian-4-one (E or Z) | 3 | 1 | - | - | MS | 1251 |
| 37 | 3,5-dimethyl-1,2-dithian-4-one (E or Z) | 1 | 0 | - | - | MS | 1261 |
| 38 | a dihydrothienothiophene | tr | - | - | - | MS | 1314 |
| 39 | a methylthienothiophene | tr | 0 | 5 | 1 | MS | 1317 |
| 40 | a methylthienothiophene | 16 | 9 | 3 | 1 | MS | 1355 |
| 41 | a dihydromethylthienothiophene | tr | 1 | 5 | 3 | MS | 1378 |
| 42 | a dihydromethylthienothiophene | 18 | 14 | tr | - | MS | 1409 |
| 43 | a dihydromethylthienothiophene | 3 | 2 | tr | - | MS | 1418 |
| 44 | bis(1-methyl-2-oxopropyl) disulfide^e,f | 2 | 1 | - | - | MS+LRI | 1469 |
| 45 | bis(1-methyl-2-oxopropyl) disulfide^e,f | 2 | 1 | - | - | MS+LRI | 1474 |
| 46 | 1-[2-methyl-(3-furyldithio)]propan-2-one^e | - | - | 3 | tr | ms | 1476 |
| 47 | **3-[2-methyl-(3-furyldithio)]butan-2-one** | **6** | **4** | **5** | **1** | MS+LRI | 1501 |
| 48 | **bis(2-methyl-3-furyl) disulfide** | **2** | **3** | **7** | **2** | MS+LRI | 1537 |
| 49 | 3-(1-methyl-2-oxopropyldithio)pentan-2-one^e | tr | tr | - | - | MS+LRI | 1539 |
| 50 | 3-(2-oxobutyldithio)butan-2-one^e | - | - | 1 | - | ms | 1547 |
| 51 | 2-(1-methyl-2-oxopropyldithio)pentan-3-one^e,f | 2 | 2 | - | - | MS+LRI | 1555 |
| 52 | 2-(1-methyl-2-oxopropyldithio)pentan-3-one^e,f | 2 | 1 | tr | - | MS+LRI | 1561 |
| 53 | 1-[2-methyl-(3-furyldithio)]butan-2-one^e | - | - | 4 | 1 | ms | 1572 |
| 54 | **3-[2-methyl-(3-furyldithio)]pentan-2-one** | **1** | **1** | **3** | **tr** | MS+LRI | 1574 |
| 55 | **2-[2-methyl-(3-furyldithio)]pentan-3-one** | **5** | **3** | **4** | **1** | MS+LRI | 1589 |
| 56 | bis(2-oxobutyl) disulfide^e | - | - | tr | - | ms | 1616 |
| 57 | 3-(1-methyl-2-oxobutyldithio)pentan-2-one^e,f | 1 | - | - | - | MS+LRI | 1624 |
| 58 | 3-(1-methyl-2-oxobutyldithio)pentan-2-one^e,f | 1 | - | - | - | MS+LRI | 1630 |
| 59 | bis(1-methyl-2-oxobutyl) disulfide^e,f | tr\* | | - | - | MS+LRI | 1647 |
| 60 | bis(1-methyl-2-oxobutyl) disulfide^e,f | tr\* | | - | - | MS+LRI | 1651 |
| 61 | 3-(3-thienyldithio)butan-2-one^e | 3 | 1 | - | - | ms | 1657 |
| 62 | **2-methyl-3-(3-thienyldithio)furan^e** | **1** | **-** | **-** | **-** | ms | 1697 |
| 63 | 3-[2-methyl-(3-thienyldithio)]butan-2-one^e | 1 | tr | 1 | - | ms | 1711 |
| 64 | **2-methyl-3-[2-methyl-(3-thienyldithio)]furan** | **tr** | **tr** | **2** | **tr** | MS | 1744 |
| 65 | 2-(3-thienyldithio)pentan-3-one^e | 2 | tr | - | - | ms | 1747 |
| 66 | 3-[2-methyl-(3-thienyldithio)]pentan-2-one^e | tr | - | tr | - | ms | 1780 |
| 67 | 1-[2-methyl-(3-thienyldithio)]butan-2-one^e | tr\* | | 1 | - | ms | 1787 |
| 68 | 2-[2-methyl-(3-thienyldithio)]pentan-3-one^e | tr | - | - | - | ms | 1795 |
| 69 | bis(2-methyl-3-thienyl) disulfide | tr\* | | - | - | ms | 1955 |

Rows 59, 60, 67 and 69 print a single **tr\*** spanning the two cysteine columns (found only in the
larger run); the layout is reproduced above with an empty second cell rather than invented.

**Mass-spectral data printed in Table 1 for the tentatively identified compounds** (footnote d says
that where no spectrum or reference is given the reference spectrum is in the NIST/EPA/NIH database):
7, m/z 57 (100), 47 (28), 104 (19), 45 (10), 42 (9); 14, 43 (100), 71 (50), 55 (22), 61 (20), 41
(15), 47 (13), 90 (7), 132 (5); 46, 43 (100), 113 (65), 114 (32), 45 (26), 202 (25), 85 (12), 51
(11), 81 (10); 50, 43 (100), 57 (63), 59 (18), 104 (14), 103 (12), 61 (10), 206 (7); 53, 57 (100),
113 (67), 114 (45), 216 (38), 43 (26), 45 (20), 81 (14), 85 (13), 51 (12); 56, 57 (100), 43 (16),
104 (10), 45 (9), 71 (4), 206 (3); 61, 43 (100), 218 (46), 115 (39), 71 (36), 116 (32), 45 (28), 59
(18), 141 (17), 111 (17), 57 (13); 62, 113 (100), 228 (54), 43 (43), 45 (43), 71 (40), 164 (24), 114
(21), 115 (20), 116 (15), 51 (14); 63, 129 (100), 43 (72), 97 (54), 34 (53), 130 (46), 232 (42), 59
(37), 57 (23), 85 (19), 125 (15); 65, 57 (100), 71 (25), 232 (24), 115 (22), 116 (20), 45 (17), 59
(13), 141 (12); 66, 43 (100), 129 (74), 97 (71), 45 (60), 130 (39), 246 (36), 85 (21), 161 (18), 73
(17), 59 (14), 41 (14); 67, 57 (100), 45 (84), 129 (67), 130 (60), 97 (49), 43 (42), 59 (32), 85
(19), 53 (10), 71 (9), 111 (8), 232 (8); 68, 57 (100), 123 (94), 130 (50), 45 (50), 246 (45), 97
(45), 59 (37), 125 (19), 85 (19); 69, 129 (100), 45 (53), 258 (39), 130 (30), 85 (22), 97 (18), 59
(14), 131 (9).

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| compounds found | **69 total; 65 in the cysteine system, 43 in the H2S system** | abstract and p. 1627 |
| class counts | disulfides **26**, thiols **7**, dithiolanones **6**, thiophenones/dihydrothiophenones **4**, dithianones **3**, alkanediones **3**, thienothiophenes **6** | abstract and p. 1627 |
| total volatile output | "Both systems produced approximately the same total quantity of volatile compounds" — **no number** | abstract |
| mercaptopentanones in the H2S system | "the quantities were only **~10 %** of those found in the cysteine-containing reaction" | p. 1629 |
| disulfides new to the literature | "**20 of these compounds had not been reported previously in any food or model system**"; 46, 50, 53, 56, 61-63, 65-69 reported for the first time | p. 1627, p. 1631 |
| MFT and 2-methyl-3-thiophenethiol | "formed in similar quantities in the hydrogen sulfide and cysteine systems, which suggests that the pathways involved require only hydrogen sulfide and not other cysteine degradation products" | p. 1630 |
| 3-thiophenethiol | "found only in the cysteine-containing system"; a known thermal degradation product of cysteine (Shu 1985) | p. 1630 |
| hydrogenated derivatives | "In the present work, **none of these hydrogenated derivatives were detected**" (di- and tetrahydro MFT / 2-methyl-3-thiophenethiol) | p. 1630 |
| why not | earlier work used higher concentrations and "a **100-fold excess of hydrogen sulfide**, whereas the present work involved dilute aqueous solutions with similar concentrations of each reactant" | p. 1630 |
| the oxidant argument | "the conditions ... were not conducive to aerial oxidation, because of the relatively high concentrations of hydrogen sulfide ... thiols are not converted to disulfides during the analysis procedures employed in this work ... **Therefore, another redox system must be involved**" | p. 1631 |
| the proposed oxidant | α-dicarbonyls reduced to hydroxyalkanones, "This could provide the redox system required for the formation of the disulfides (Figure 6)" | p. 1631 |
| disulfide occurrence pattern | 3-thienyl-containing disulfides **only** in the cysteine systems; disulfides from 1-mercapto-2-butanone **only** in the H2S system; mercaptoalkanone-derived disulfides, with one exception, **only** in the cysteine system | p. 1631 |
| odour threshold, bis(2-methyl-3-furyl) disulfide | **2 x 10^-5 µg/kg** (cited to Buttery 1984, **not measured here**) | p. 1631 |
| odour threshold, 3-thiophenethiol | **5-10 µg/kg** (cited, not measured here) | p. 1630 |
| 2,3-pentanedione in the H2S system | not isolated, although the mercaptopentanones were formed | p. 1629 |

**No figure in this paper carries a datum.** Figures 1-12 are all proposed pathways or a proposed
redox couple.

### Arithmetic on the printed table (all mine)

**Column check — the reading of the four numeric columns is verified three ways.** The companion
paper (Whitfield & Mottram 2001) prints, in its Table 2, the **class totals of exactly these data**
under the heading "pH 4.5 (from ref 8)". Summing this paper's Table 1 under the reading "cys-1,
cys-2, H2S-1, H2S-2" and taking duplicate means with tr ≈ 0:

| class | my sum, cysteine | companion's printed value | my sum, H2S | companion's printed value |
|---|---:|---:|---:|---:|
| dithiolanones + dithianones (18, 21, 22, 27, 29, 30, 34, 36, 37) | **40.5** | **41** | **443** | **443** |
| thiols + mercaptoketones (5, 6, 7, 8, 9, 10, 14, 17) | **294** | **294** | **43.5** | **44** |
| disulfides (44-69) | **~26** | **24** | **~18** | **18** |

Three independent classes reproduce, one of them (the H2S dithiolanones, 443) exactly. **The column
assignment is correct** and the two papers' data are the same data.

**1. The three shipped FIT rows, checked against the print.**

| generator row | target | Table 1 | verdict |
|---|---|---|---|
| `whitfield_nf_cys_MFT` | 0.150 mol % | row 6 cysteine, 15 and 15 µg/10 mg → mean 15 → **0.150 mol %** | **verified** |
| `whitfield_nf_h2s_MFT` | 0.120 mol % | row 6 H2S, 16 and 8 → mean 12 → **0.120 mol %** | **verified**, but the duplicates differ 2-fold (Flags 3) |
| `whitfield_mercaptoketone_over_MFT` | 16.3 | rows 5, 8, 9 cysteine means 92 + 74.5 + 77.5 = **244**, against MFT 15 → **16.27** | **verified**; note it counts the three classical mercaptoketones and excludes rows 7 and 14 |
| the `r_nf_mp3p` note, "MFT only 2.6 % of everything fed NF produces" | 2.6 % | 15 against the companion's printed total of **580** → **2.59 %** | **verified**, but the 580 is printed in the *companion*, not here |
| `cerny_isomer_split`'s supporting statement, "the two isomers at 74.5 : 77.5" | ~1:1 | rows 8 and 9 cysteine means | **verified** |

**2. The disulfide share — the number this cluster was sent for.**
Free thiols and mercaptoketones, cysteine system: **294 µg/10 mg NF**. Disulfides: **24 µg/10 mg NF**
(the companion's printed class total; my own sum of rows 44-69 gives ~26 with tr ≈ 0). A disulfide's
mass is very nearly the sum of its two thiol moieties less 2 amu, so µg of disulfide is µg of
thiol-equivalent to within about 1 %. **Disulfide-bound share of the thiol pool = 24/(294+24) =
7.5 %** in the cysteine system and **18/(44+18) = 29 %** in the H2S system.

For MFT alone, converting each MFT-bearing disulfide to MFT equivalents at MFT MW 114.17:

| no. | compound | cys mean, µg | MW | MFT moieties | MFT-eq, µg |
|---:|---|---:|---:|---:|---:|
| 47 | MFT–S–S–(3-mercaptobutan-2-one) | 5 | 216.3 | 1 | 2.64 |
| 48 | bis(2-methyl-3-furyl) disulfide | 2.5 | 226.3 | 2 | 2.52 |
| 54 | MFT–S–S–(3-mercaptopentan-2-one) | 1 | 230.4 | 1 | 0.50 |
| 55 | MFT–S–S–(2-mercaptopentan-3-one) | 4 | 230.4 | 1 | 1.98 |
| 62 | MFT–S–S–(3-thienyl) | 0.5 | 228.4 | 1 | 0.25 |
| 64 | MFT–S–S–(2-methyl-3-thienyl) | ~0.25 (tr) | 242.4 | 1 | 0.12 |
| | **total MFT-equivalent bound** | | | | **~8.0** |

Against **15 µg of free MFT**, that is **8.0/(15 + 8.0) = 35 % of the MFT in disulfide form**, and a
**true MFT yield of 23.0 µg/10 mg NF = 0.230 mol %** rather than the 0.150 the FIT row targets.
*Caveats carried:* response factors are assumed 1 for the disulfides as for everything else, and a
disulfide of MW 226 is much less volatile than a thiol of MW 114, so a dynamic-headspace method at
60 C will under-recover the disulfides — **35 % is therefore a lower bound if anything.** The
duplicate spread on the disulfide rows is also wide.

**3. The same three shares side by side (mine).**

| source | temperature | pH | disulfide share of the thiol |
|---|---:|---:|---:|
| **this paper, cysteine system** | **140 C** | **4.5** | **~7.5 % of all thiol; ~35 % of MFT** |
| this paper, H2S system | 140 C | 4.5 | ~29 % of all thiol |
| Zhou 2023 | 120 C | 6 / 7 / 8 | 8.6 / 6.5 / 9.6 % |
| Zhang 2024, cysteine arm | 115 C | 4.9 | 8.7 % (figure-read; the printed text supports only an ordering) |
| **B16/B17 model** | 100-145 C | 6-8 | **0.03 to 0.93 %** |

**4. The cysteine-versus-H2S ratios (mine), which are what a fitted model must reproduce.**
Mercaptoketones: 3-mercaptobutan-2-one **16.7x** more in the cysteine system (92 vs 5.5),
3-mercaptopentan-2-one **9.9x** (74.5 vs 7.5), 2-mercaptopentan-3-one **11.1x** (77.5 vs 7) — the
paper's own "~10 %" summary. **MFT: 1.25x (15 vs 12) and 2-methyl-3-thiophenethiol 1.7x (7.5 vs
4.5)** — essentially flat, which is the paper's mechanistic point. Dithiolanones run the other way:
3-methyl-1,2-dithiolan-4-one is **133x** higher with H2S (200 vs 1.5) and the class **10.9x** (443 vs
40.5). 2,4-Pentanedione is **2.8x** higher with H2S (39 vs 14). **The sulfur in this pot is
partitioned completely differently by the two donors even though the total output is similar** — a
strong, cheap structural test that costs the objective nothing to score, since both systems are
already configured as fit rows.

**5. The two mercaptopentanone isomers are 1:1 here (74.5 : 77.5 = 0.96)** and the module uses that
against Cerny 2007's isotope split. Worth noting alongside it that **2-mercapto-3-pentanone, which
Cerny 2003 could not detect at all from ribose + cysteine at 95 C, is one of the two largest
mercaptoketones here at 140 C from fed norfuraneol** — consistent with Cerny's own reading that the
norfuraneol route makes it (96 % unlabelled in his pot D) and the ribose route does not.

**6. A charge mismatch in the fit generators (mine, and it is a defect to check).** Every B2.x
generator configures the two Whitfield systems as `{"NF": 20.0, "Cys": 20.0}` and
`{"NF": 20.0, "H2S": 40.0}` in mmol/L. **The paper's charges are 50 mmol/L norfuraneol with
50 mmol/L cysteine, or 50 mmol/L norfuraneol with ~97 mmol/L H2S.** The 1:1 and ~1:2 ratios are
preserved but the absolute concentration is **2.5x low**, and the H2S:NF ratio is set to exactly 2.0
where the paper's is 1.94. For a mol %-of-fed target this is harmless only if every step on the path
is first order in norfuraneol; the H2S-addition steps are written second order, so a 2.5x
concentration error propagates into those rows. Two of the three shipped rows from this paper are
mol % targets and one is a ratio.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Keyed: `norfuraneol` (the fed reactant),
`2_methyl_3_furanthiol`, `bis_2_methyl_3_furyl_disulfide`, `hydrogen_sulfide`,
`2_3_butanedione` (inferred but not measured here). **Not keyed and central to this paper:**
3-mercaptobutan-2-one, 3-mercaptopentan-2-one, 2-mercaptopentan-3-one, 1-mercaptobutan-2-one,
5-mercaptohexan-2-one, 3-thiophenethiol, 2-methyl-3-thiophenethiol, all 26 disulfides except
compound 48, the dithiolanones and dithianones, the dihydrothiophenones, 2,3-pentanedione,
2,4-pentanedione, 3,4-hexanedione, 2-acetyl-5-methylthiophene, 3-ethyl-2-formylthiophene, and the
thienothiophenes. See Flags 8.

Every row below shares: **norfuraneol 50 mmol/L, 0.5 M phosphate at pH 4.5, flame-sealed 5 mL glass
ampule, 140 C for 60 min, duplicate, dynamic headspace at 60 C for 1 h onto Tenax, GC/MS against
100 µg of methyl decanoate with all response factors assumed 1, detection limit 0.1 µg per 10 mg of
norfuraneol.** The partner is **cysteine 50 mmol/L** or **H2S ~97 mmol/L**. Conversion:
mol % = µg / (MW x 0.8764).

| step / quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **NF + cysteine -> MFT** | 15, 15 (mean **15**) = **0.150** | µg/10 mg NF; mol % | cysteine system | none — a single-point yield | Table 1 row 6 | **fed_intermediate_yield** — the shipped `whitfield_nf_cys_MFT` row |
| **NF + H2S -> MFT** | 16, 8 (mean **12**) = **0.120** | µg/10 mg NF; mol % | H2S system | none | Table 1 row 6 | fed_intermediate_yield — the shipped `whitfield_nf_h2s_MFT` row |
| NF + cysteine -> 3-mercaptobutan-2-one | 84, 100 (mean 92) = **0.883 mol %** | as above | cysteine | none | Table 1 row 5 | fed_intermediate_yield |
| NF + cysteine -> 3-mercaptopentan-2-one | 81, 68 (mean 74.5) = **0.719 mol %** | as above | cysteine | none | Table 1 row 8 | fed_intermediate_yield |
| NF + cysteine -> 2-mercaptopentan-3-one | 78, 77 (mean 77.5) = **0.748 mol %** | as above | cysteine | none | Table 1 row 9 | fed_intermediate_yield |
| the same three with H2S | 5.5 / 7.5 / 7 | µg/10 mg NF | H2S | none | Table 1 rows 5, 8, 9 | fed_intermediate_yield |
| **mercaptoketones : MFT** | **16.3** | — | cysteine | — | derived from Table 1 (mine) | **within_study_ratio** — the shipped ratio row |
| MFT isomer partner, 2-methyl-3-thiophenethiol | 10, 5 (cys) / 9, - (H2S) | µg/10 mg NF | both | none | Table 1 row 17 | fed_intermediate_yield |
| 3-thiophenethiol | 26, 14 (cys) / not detected (H2S) | µg/10 mg NF | both | none | Table 1 row 10 | fed_intermediate_yield — **cysteine-derived**, cf. Cerny 2003's 95 % unlabelled |
| 3-methyl-1,2-dithiolan-4-one | 2, 1 (cys) / **204, 196** (H2S) | µg/10 mg NF | both | none | Table 1 row 18 | fed_intermediate_yield — the largest single product in the paper |
| dihydro-2-methylthiophen-3(2H)-one | 66, 64 (cys) / 26, 33 (H2S) | µg/10 mg NF | both | none | Table 1 row 13 | fed_intermediate_yield |
| **bis(2-methyl-3-furyl) disulfide** | 2, 3 (cys) / 7, 2 (H2S) | µg/10 mg NF | both | none | Table 1 row 48 | fed_intermediate_yield |
| **total MFT-bearing disulfide** | **~8.0** MFT-equivalent | µg/10 mg NF | cysteine | — | derived from Table 1 rows 47, 48, 54, 55, 62, 64 (mine) | **derived_assumption** (response factors 1; a lower bound because disulfides are less volatile) |
| **MFT disulfide-bound share** | **~35** | % of the MFT made | cysteine, 140 C, pH 4.5 | — | derived (mine) | **within_study_ratio (mine)** |
| **total MFT yield including bound** | **0.230** | mol % of NF fed | cysteine | — | derived (mine) | derived_assumption |
| **whole-pool disulfide share** | **7.5** (cysteine) / **29** (H2S) | % of thiol as disulfide | 140 C, pH 4.5 | — | derived from the companion's printed class totals (mine) | **within_study_ratio (mine)** |
| cysteine : H2S ratios | mercaptoketones 10-17x; **MFT 1.25x**; 2-methyl-3-thiophenethiol 1.7x; dithiolanones 0.09x | — | 140 C, pH 4.5 | — | derived (mine) | within_study_ratio |
| mercaptopentanone isomer ratio | **74.5 : 77.5 = 0.96** | — | cysteine | — | Table 1 rows 8, 9 | within_study_ratio — used by `cerny_isomer_split` |
| di- and tetrahydro MFT / thiophenethiol | **not detected** | — | both systems | — | p. 1630 | **threshold** |
| 2,3-pentanedione in the H2S system | not detected | — | H2S | — | Table 1 row 1 | threshold |
| bis(2-methyl-3-furyl) disulfide odour threshold | 2 x 10^-5 µg/kg | µg/kg in water | — | — | cited to Buttery 1984 | **threshold** (**not measured here** — do not attribute) |
| 3-thiophenethiol odour threshold | 5-10 µg/kg | µg/kg | — | — | cited | threshold (not measured here) |
| the oxidant behind the disulfides | α-dicarbonyl reduction, coupled to the mercaptoketone-forming flux; aerial oxidation and analytical artefact both argued out | — | 140 C, pH 4.5 | — | p. 1631 + Figure 6 | **level_only** — a mechanism proposal by the authors, with two supporting negatives |

### Can these be put on the same basis as the shipped rows? Yes for the yields; no for the disulfides — yet.

- **The three shipped rows are correct** and now rest on the primary table rather than on an
  inventory summary. Two corrections travel with them: the **charge mismatch** (section 3 item 6) and
  the **free-versus-total MFT question** (item 2). The second is the more consequential: if the
  model's `MFT` is the whole pool and the target is free MFT, the row is systematically 1.5x low, and
  the honest fix is either to target 0.230 mol % with a stated derivation or to score free MFT
  against a model that carries the disulfide explicitly — which is exactly what the `ch_dimer_*` and
  `ch_dimer_release_*` steps exist for.
- **The disulfide share is a candidate FIT or hold-out row that does not exist yet.** It is a
  *ratio* within one analysis, so it is immune to the response-factor caveat that limits the levels,
  and it sits at a temperature (140 C) between the fed panel's 145 C and Zhou's 120 C. Making it a
  row would give the objective its **only** disulfide constraint that is not at 115-120 C, and would
  do it inside a pot the objective already simulates. Whether to score it as FIT or hold-out is the
  owner's call; the declaration currently lists this paper as FIT and Zhou/Zhang's dimer fractions as
  hold-outs, so a disulfide row from this paper would be the first fitted one.
- **Nothing here transports as a rate.** One time point, one temperature.

## 5. Flags

1. **Relative quantification with all response factors assumed 1, by dynamic headspace at 60 C.**
   Every µg in Table 1 is a methyl-decanoate-equivalent of a headspace-recoverable amount, not a
   solution concentration. This biases *systematically by volatility*, which matters most for the
   disulfides (MW ~220-260 against the thiols' ~114) — so the disulfide share computed in section 3
   is a **lower** bound, and the mercaptoketone : MFT ratio, comparing molecules of similar
   volatility, is the safest number in the paper.
2. **The measurement is of the sealed ampule's chemistry but the analysis is of a diluted, aerated,
   stirred flask.** The reaction mixture is transferred into 20 mL of buffer in a 250 mL conical flask
   and purged with nitrogen for an hour at 60 C. The authors address the obvious worry directly —
   Mottram et al. 1998 is cited for the finding that thiols are not converted to disulfides during
   this procedure — and that citation is load-bearing for the whole disulfide argument. It is not
   re-measured here.
3. **Duplicate spread is wide and is not summarised.** Both values are printed, which is admirable,
   but the pairs include 4/15 (row 1), 16/8 (row 6 H2S — the very row a shipped FIT target is the
   mean of), 26/14 (row 10), 57/17 (row 32), 12/1 (row 28 H2S), 48/30 (row 27 H2S). **The H2S MFT
   row's duplicates differ 2-fold**, so the 0.120 mol % target carries at least a factor-of-1.4
   uncertainty either way and the generator's `sigma_log=0.5` is not generous.
4. **The "higher concentration" run is described two incompatible ways.** Methods: "at a higher
   concentration (2-fold) and in greater quantities (10 mL each)". Table 1 footnote a: "found only in
   reaction mixtures with **10 times higher quantities** of cysteine and HMF". 2-fold concentration in
   10 mL rather than 1 mL is 20x the material, not 10x. The eight **tr\*** entries (59, 60, 67, 69)
   depend on which is meant, but no quantity is attached to any of them, so nothing numerical hangs
   on it.
5. **Three of the seven thiols and most of the disulfides are tentative identifications.**
   3-Thiophenethiol (10) and 2-methyl-3-thiophenethiol (17) are "MS" only — by comparison with a
   literature spectrum, not an authentic standard — and 1-mercaptobutan-2-one (7) and
   5-mercaptohexan-2-one (14) are "ms", i.e. interpreted from the spectrum with no reference at all.
   **MFT itself (6) is MS + LRI against an authentic compound**, so the shipped rows are on the firm
   side of this line; the disulfide sums in section 3 are not — compounds 62, 63, 65-68 are all "ms".
6. **The paper prints no total**, so the "MFT is 2.6 % of everything the fed NF produces" claim in
   `sulfur.py` is arithmetic on **the companion paper's** Table 2 (580 µg), not on this one. That is
   legitimate — the companion states its column is "from ref 8" and three class totals reproduce
   against Table 1 — but the citation in the code should name both papers.
7. **What to request from the authors**: (i) a direct quantification of the disulfides against
   authentic standards, which would turn section 3's 7.5 % and 35 % from response-factor-assumed
   estimates into measurements; (ii) the α-dicarbonyl concentrations, which would test the Figure 6
   redox proposal stoichiometrically — the paper says 2,3-butanedione and 2-oxobutanal were "too
   volatile" to be entrained, which is exactly the pool the oxidant argument depends on; (iii) a
   time course at 140 C, the absence of which means nothing here bounds a rate; (iv) the norfuraneol
   remaining at 60 min, without which no yield in this paper can be put on a converted-substrate
   basis; (v) confirmation of whether the larger run was 10x or 20x the material.
8. **Registry gaps against `data/keys/compounds.yml`.** The two that matter most are
   **3-mercaptopentan-2-one** and **2-mercaptopentan-3-one** — they carry the isomer-split
   diagnostic the module scores and are absent from the registry, as this cluster's `cerny2003`
   dossier also records. After those, **3-mercaptobutan-2-one** (the largest thiol in the paper) and
   the mixed disulfide **3-[2-methyl-(3-furyldithio)]butan-2-one** (compound 47, the largest single
   MFT-bearing disulfide, and the same species the pH-6.5 companion finds as its only MFT-containing
   product). No dithiolanone, dithianone or thienothiophene has a key, and the dithiolanones are the
   dominant product class of the whole H2S system.
9. **A defect to check in the generators, not in the paper.** All B2.x generators configure
   `whitfield_nf_cys` and `whitfield_nf_h2s` at **NF 20 mmol/L** against the paper's **50 mmol/L**,
   and H2S at 40 against ~97 (section 3 item 6). Ratios are preserved; absolute concentrations are
   2.5x low. This dossier does not change any code.
10. **What this paper does not contain**: any time course; any second temperature; any pH other than
    4.5; any rate or barrier; any norfuraneol mass balance; any measurement of H2S consumed; any
    absolute solution concentration; any oxygen or headspace-composition variation; any statement of
    the ampule's headspace volume or atmosphere beyond "flame sealed"; and any supplementary
    material.
