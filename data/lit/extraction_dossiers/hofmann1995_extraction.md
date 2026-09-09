# Hofmann & Schieberle 1995 — EXTRACTION (ribose 100 mmol/L + cysteine 33 mmol/L in 0.5 mol/L phosphate pH 5.0, 100 mL, ramped from 20 to 145 C over 20 min in an autoclave; aroma extract dilution analysis of an acidic and a neutral/basic fraction, static-headspace olfactometry, and fourteen odour thresholds in air — thirty odorants ranked, and NOT ONE concentration)

### THE PAPER THAT DEFINES THE POT the whole sulfur lane is built on, and it contains no quantity: every number in it is a flavour dilution factor or an odour threshold. On the thiol-sink question it carries **no rate, no yield, no binding plateau and no mass balance** — but it carries the one caution that most bears on the sink evidence in this cluster: **the disulfides in a solvent-extracted ribose/cysteine pot are largely workup artefacts** ("no disulfide was detected in the headspace analysis"; more than 50 % of MFT stored in diethyl ether converts to its disulfide), so the disulfide share of the reference pot at 145 C has never been measured by a method that could see it.

**Source on disk:** `data/articles/hofmann1995.pdf` (8 pp., J. Agric. Food Chem. 1995, 43 (8),
2187-2194). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/hofmann1995.txt`, 569 lines). **The text layer of this 1995 scan is OCR'd and
is the worst in this cluster**: the body text is legible throughout, but Table 1 carries a dozen
character-level corruptions — `<1` variously rendered `41`, `'1`, `(1`, `<l`; `1431` rendered
`14.31`; `2145` rendered `2 145`; `<800` rendered `e800`; the fraction label `AF I/NB I` rendered
`AFJf-NBI`; and several compound names run together. **Every such cell is flagged in the re-typed
table below and each correction is justified against the paper's own running text**, which quotes
most of the important values independently. Tables 2 and 3 came through clean, and Table 2 is
internally self-checking (section 3). Figures 1, 2 and 3 are mass spectra and Scheme 1 is a proposed
pathway — **no figure in this paper carries a datum**. There is no supplementary material. Repo
status before this dossier: Hofmann & Schieberle 1995 is the pot that `van Seeventer 2001`
reproduces, the parent of the `hofmann1998` papers already dossiered
(`hofmann1998b_extraction.md`, `hofmann1998_reconciliation.md`), and it is cited by Cerny 2003 (ref
8) and by Whitfield 1999. **It has never had a dossier of its own.**

## 0. Identity

| field | value |
|---|---|
| Title | "Evaluation of the Key Odorants in a Thermally Treated Solution of Ribose and Cysteine by Aroma Extract Dilution Techniques" |
| Authors | Thomas Hofmann and Peter Schieberle (corresponding) — Deutsche Forschungsanstalt für Lebensmittelchemie and Institut für Lebensmittelchemie der Technischen Universität München, Lichtenbergstrasse 4, D-85748 Garching |
| Venue | J. Agric. Food Chem. 1995, 43 (8), 2187-2194. Received 20 October 1994, revised 4 April 1995, accepted 26 May 1995; Advance ACS Abstracts 1 July 1995 |
| DOI / article ID | printed as `JF9405968` |
| Abbreviations | AEDA = aroma extract dilution analysis; FD = flavour dilution factor; SHO = static headspace analysis/olfactometry; AF = acidic fraction; NB = neutral/basic fraction; ADT = 5-acetyl-2,3-dihydro-1,4-thiazine; HDF = 4-hydroxy-2,5-dimethyl-3(2H)-furanone (the registry's `hdmf`); norfuraneol = 4-hydroxy-5-methyl-3(2H)-furanone (the repository's `NF`) |
| Novelty | first report of **5-acetyl-2,3-dihydro-1,4-thiazine** among the volatiles of any Maillard model system or food; its structure paper is the companion Hofmann, Hassner & Schieberle 1995, JAFC 43:2195-2198 |
| Companions on disk | `hofmann1996_extraction.md` (= "Hofmann et al. 1995b, submitted" in this paper's own reference list — the ether-storage thiol oxidation study whose > 50 %/10 days figure this paper quotes in advance), `hofmann1998b_extraction.md` and `hofmann1998_reconciliation.md` (the quantitative sequels), this cluster's `vanseeventer2001_extraction.md` (which reproduces this pot at 130 C and stores it), `cerny2003_extraction.md` (which cites this paper for the 145 C mercaptopentanone pair) and `whitfield1999_extraction.md` |

## 1. Why it matters

**What it contributes to the thiol-sink question. No number — and one caution that changes how a
number in this same cluster should be read.**

This paper prints thirty flavour dilution factors, fifteen static-headspace dilution factors and
fourteen odour thresholds. **It prints no concentration, no yield, no rate and no mass balance.** An
FD factor is a dilution step on a factor-of-two ladder; an odour threshold is a threshold. Neither
is a quantity of substance in the pot. So on the plain question — does this paper carry a thiol loss
rate, a binding plateau, a mass balance that does not close, or a measured sink partner that the
sink objective has never seen? — **the answer is no, and nothing in it could have changed the
objective that refused the three sink structures.**

But it carries a methodological finding that bears directly on the disulfide evidence, and it is
worth more to this cluster than most numbers would be. p. 2191, in full:

> "Unexpectedly, the FD factors and, therefore, the concentrations of these neutral compounds were
> higher in the acidic than in the neutral/basic fraction (cf. columns I and II, Table 1),
> **suggesting that the disulfides were formed as artifacts during the workup procedure.** This
> assumption was very recently corroborated by our findings (Hofmann et al., 1995b) that after 10
> days storage of, e.g., 2-methyl-3-furanthiol in diethyl ether at 6 °C, **more than 50 % of the
> thiol had been converted into the corresponding disulfide**."

and

> "Interestingly, **no disulfide was detected in the headspace analysis**, corroborating the
> assumption that a significant proportion of these flavor compounds detected by AEDA result as
> artifacts during the workup procedure."

Three consequences, and they need keeping apart.

1. **Eight disulfides are seen in this 145 C ribose/cysteine pot** (compounds 18, 20, 22, 27, 28, 29,
   30 and the mixed pair), and one of them — bis(2-methyl-3-furyl) disulfide — is among the six
   highest-FD odorants in the paper. But their authors do not believe them, on an internally
   consistent argument (a *neutral* compound should not enrich into the *acidic* fraction).
2. **The static-headspace check is weak in the same direction and the paper says so.** No disulfide
   was detected in the headspace — but the paper's own explanation for why most of Table 1's
   compounds are missing from Table 2 is "possibly due to their lower vapor pressure". A disulfide of
   MW 226 against a thiol of MW 114 is exactly that case. **Absence from a static headspace at 40 C
   is not evidence of absence in the pot.**
3. **This does NOT undercut the disulfide share measured in this cluster's `whitfield1999`
   dossier**, and the reason is method. Hofmann's artefact mechanism is *storage of a solvent
   extract in diethyl ether*; Whitfield & Mottram used **dynamic headspace onto Tenax with no solvent
   contact at any stage** and cite a dedicated study (Mottram et al. 1998) for the finding that their
   procedure does not convert thiols to disulfides. Cerny 2003 makes the same point about his own
   SPME work, citing this laboratory's ether result (his ref 40) and noting that "solvent extraction
   was avoided in the present study". **So the three papers are consistent: solvent extraction makes
   disulfides, headspace methods do not, and the disulfides Whitfield measures at 140 C are real.**
   What follows for the lane is a gap, not a contradiction: **the disulfide share of the reference
   ribose/cysteine pot at 145 C has never been measured by a method capable of seeing it**, and the
   experiment that would fix that is the one the B25 outcome note already calls for.

**What this paper actually gives the repository, which is not small.** It **defines the pot**. The
sulfur module's reference temperature `T_REF_S_K = 418.15` (145 C) is described in
`parameters_sulfur.py` as "the temperature of the Hofmann 1998 SIDA anchors and of the Cerny isotope
experiments" — and this is the paper in which that pot is specified: **ribose 100 mmol/L, cysteine
33 mmol/L, 0.5 mol/L phosphate at pH 5.0, 100 mL, in a laboratory autoclave, taken from 20 to
145 C over 20 minutes.** Note what that last clause says: **the "20 min" is the ramp, not a hold**
(Flags 2). It also supplies:

- **The ranking the whole lane's species list is chosen from.** The six highest-FD odorants are
  2-furfurylthiol (1024), 3-mercapto-2-pentanone (512), 2-methyl-3-furanthiol (256),
  5-acetyl-2,3-dihydro-1,4-thiazine (256 in the neutral/basic fraction), 3-mercapto-2-butanone (128)
  and bis(2-methyl-3-furyl) disulfide (128). Every one except the thiazine is a species the network
  carries.
- **Fourteen odour thresholds in air, measured here by HRGC-olfactometry** (Table 3), four of them
  adapted from the literature and ten of them this paper's own. These are the kind of number
  `matrix_oav` and `desirable_targets.yml` need and mostly lack — they are **thresholds in air**, not
  in water or in a matrix (Flags 5).
- **A pH and a stoichiometry screen, reported as sensory description only, data not shown**: pH 5.0
  gives "sulfury, meatlike" and pH 7.0 "caramel, burnt"; a 1:1 ribose:cysteine mixture is "rubbery,
  sulfury, and H2S-like" and 10:1 "burnt and caramel-like"; **3:1 (100 : 33 mmol/L) gives "the most
  prominent roasty, meatlike odor"**, and that is why the pot is 3:1.
- **The observation that 2-thenyl mercaptan and ethyl mercaptan are major headspace odorants that
  the solvent extract badly under-represents** (2-thenyl mercaptan FD 8 in the extract against 125 in
  the headspace; ethyl mercaptan not detected at all in the extract, FD 125 in the headspace) —
  "Presumably, the last thiol is relatively labile, which would lead to significant losses of this
  odorant during concentration of the extract." **A thiol that is lost during a room-temperature
  solvent concentration is a thiol whose measured levels anywhere are suspect**, which is a general
  warning this corpus should carry.

What this paper does NOT give the repository: any concentration of anything; any yield; any rate;
any time course; any second temperature; any pH series with data attached; any mass balance; any
statement of how much ribose or cysteine reacted; and any quantification of the disulfides it
believes are artefacts.

## 2. Methods as they matter to a model

- **The pot.** "L-Cysteine (**3.3 mmol**) and D-ribose (**10 mmol**) in phosphate buffer (**100 mL,
  0.5 mol/L, pH 5.0**) were allowed to react in a laboratory autoclave (Model II; Roth, Karlsruhe)
  **by raising the temperature within 20 min from 20 to 145 C**." So **ribose 100 mmol/L, cysteine
  33 mmol/L, molar ratio 3:1**, and Table 1's footnote a repeats "reacted for **20 min**". **No hold
  at 145 C is described** (Flags 2). Headspace of the autoclave: not stated. Stirring: not stated.
- **Why this pot.** Preliminary sensory work (data not shown) on pH 5.0 versus 7.0 and on
  ribose:cysteine from 1:1 to 10:1, at a fixed ribose of 100 mmol/L with cysteine from 100 down to
  10 mmol/L. Verdicts quoted in section 1.
- **Extraction and fractionation.** After cooling, the "dark yellow" mixture was extracted **five
  times with diethyl ether (total 300 mL)**; the combined organic phase was washed **three times with
  0.5 mol/L aqueous sodium bicarbonate (total 100 mL)** to remove acids. The ether layer is the
  **neutral/basic fraction NB**; the bicarbonate layer, re-acidified to **pH 3.0** with 1 mol/L HCl
  and back-extracted with ether (150 mL), is the **acidic fraction AF**. Each fraction was then
  distilled from non-volatile material (SAFE-precursor apparatus, Sen et al. 1991).
- **Column chromatography of NB.** Concentrated to 1 mL on a Vigreux column, then flash silica
  (30-60 µm) at 5 mL/min under nitrogen, water-cooled to 10-12 C, eluted in six cuts:
  **NB I** n-pentane 150 mL; **NB II** pentane/ether 95+5; **NB III** 8+2; **NB IV** 1+1;
  **NB V** ether/methanol 99.5+0.5; **NB VI** ether/methanol 9+1.
- **Thiol enrichment of AF by affinity chromatography** (Full & Schreier 1994): Bio-Rad Affi-Gel 501
  (**phenylmercury chloride**) in 2-propanol; unbound compounds eluted with pentane/dichloromethane
  2+1 (60 mL) = **AF I**; the mercury-complexed thiols released with the same solvent containing
  **1,4-dithiothreitol (10 mmol)** = **AF II**. The paper notes with interest that this thiol-affinity
  step also enriched **enol-oxo compounds** — norfuraneol, HDF, 3-hydroxy-2H-pyran-2-one and sotolon
  all appear in AF II. **A dithiothreitol elution is a reducing elution**, which is worth noticing
  for a paper arguing about disulfides (Flags 3).
- **AEDA.** The original 100 µL extract of NB or AF from **100 mL of reaction mixture** was diluted
  stepwise **1+1** with diethyl ether; 1.0 µL aliquots by HRGC-olfactometry, **FFAP for the acidic
  fraction and DB-5 for the neutral/basic**. **Two assessors, in duplicate; "the data differed to not
  more than two FD factors"** (Table 1 footnote d) — i.e. up to a factor of four.
- **Static headspace olfactometry (SHO).** **100 mL of the reaction mixture** in a septum-sealed
  **240 mL** thermostatted vessel, **equilibrated 1 h at 40 C**; decreasing headspace volumes from
  **20 mL down to 0.02 mL** taken by gastight syringe onto the FFAP column. **FD = 20 mL divided by
  the smallest volume in which the odorant was still detectable**, so FD = 1 by definition at 20 mL.
- **Odour thresholds.** By HRGC-olfactometry (Ullrich & Grosch 1987), **with (E)-2-decenal (threshold
  2.7 ng/L of air) as the reference odorant instead of hexanal**. Reported as **ng per litre of
  air**.
- **Identification.** Retention index on **two** capillaries (FFAP and SE-54), **MS in both EI (70 eV)
  and CI (isobutane, 115 eV)**, and odour quality at the sniffing port — a four-criterion standard,
  and unusually strict for 1995. Reference compounds were bought or synthesised in-house; the
  syntheses of 2-methyl-3-thiophenethiol, 2-thenyl mercaptan, 3-mercapto-2-butanone,
  1-mercapto-2-propanone and five **mixed disulfides** (by CuSO4 oxidation of thiol pairs, then
  RP-18 HPLC) are given in full.
- **Scale for identification work.** For the AF structure work, **25 runs were combined (2.5 L)**;
  for the NB work, the neutral/basic volatiles from **2.5 L** were flash-chromatographed.

## 3. Tables re-typed

### Table 1. "Intense Odorants (FD >= 4) Formed by Heating of an Aqueous Ribose/Cysteine Solution^a"

Footnote a exactly as printed: "A mixture of ribose (10 mmol) and cysteine (3.3 mmol) was reacted
for 20 min in phosphate buffer (100 mL; 0.5 mol/L; pH 5.0)." Footnote b: identification by retention
index on two capillaries, MS (EI) and MS (CI), and odour quality at the sniffing port. Footnote c:
"Fraction in which most of the compound appeared after chromatography on silica gel or Affi-Gel,
respectively." Footnote d: "Flavor dilution (FD) factor determined in extracts containing the acidic
(I) or the neutral/basic volatiles (II). **Analyses were performed by two assessors in duplicates.
The data differed to not more than two FD factors.**" Footnote e: "Reported in the literature as
volatile reaction product of cysteine, cystine, or cysteamine reacted in the presence of
carbohydrates or carbohydrate degradation products", with the numbered key given below. Footnote f:
"**3-Mercapto-2-pentanone contained small amounts of its isomer 2-mercapto-3-pentanone.**"
Footnote g: "The MS signals were too weak for an unequivocal interpretation. The compound was
identified on the basis of the remaining criteria given in footnote b."

**⚠ OCR NOTE.** Cells marked † below were corrupted in the text layer and are given here as the
only chemically and typographically possible reading; the corruption and the justification are
listed under the table. No value has been guessed where the reading was not forced.

| no. | odorant | fraction | RI FFAP | RI SE-54 | FD I (acidic) | FD II (neutral/basic) | reported in |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | 3-mercapto-2-butanone | AF II | 1282 | 820 | **128** | <1 | 5, 7 |
| 2 | 2-methyl-3-furanthiol | AF II | 1300 | 869 | **256** | <1 | 4, 6 |
| 3 | 3-mercapto-2-pentanone^f | AF II | 1347 | 907 | **512** | <1 | 6, 7 |
| 4 | 1-mercapto-2-propanone | AF II | 1357 | <800 † | 16 | <1 | 5 |
| 5 | 2-furfurylthiol | AF II | 1431 † | 909 | **1024** | <1 † | 3-7 |
| 6 | 2-methyltetrahydrothiophen-3-one | AF I / NB III | 1512 | 1017 | 8 | 16 | 3-7 |
| 7 | 2-methyl-3-thiophenethiol | AF II | 1559 | 1059 | 16 | <1 | 2, 6 |
| 8 | 2-acetylthiazole | NB III | 1615 | 1142 | <1 | 8 | 1, 3-6, 11 |
| 9 | 2-formylthiophene | NB III | 1669 | 1000 | <1 | **64** | 4-6 |
| 10 | 2-thenyl mercaptan | AF II | 1682 | 1092 | 8 | <1 | 4, 5, 12 |
| 11 | 2-acetyl-2-thiazoline | NB IV | 1720 | 1111 | <1 | **64** | 8 |
| 12 | acetylthiophene | NB IV | 1740 | 1090 | <1 † | 16 | 4-6, 11 |
| 13 | unknown | — | 1784 | — | <1 | 16 | |
| 14 | unknown | — | 1800 | — | <1 | 8 | |
| 15 | 2-(hydroxymethyl)thiophene | AF I / NB IV | 1917 | 1030 | 2 | <1 | 5 |
| 16 | 3-hydroxy-2H-pyran-2-one | AF II | 1965 | 989 | 16 | <1 | 5 |
| 17 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone (HDF) | AF II | 2016 | 1100 | 32 | <1 | 5, 11 |
| 18 | bis(2-methyl-3-furyl) disulfide | AF I / NB I † | 2100 | 1526 | **128** | 16 | 7, 8, 10 |
| 19 | 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol) | AF II | 2105 | 1044 | **64** | <1 | 5, 6 |
| 20 | 2-methyl-3-furyl 2-oxo-3-butyl disulfide | AF I | 2127 | 1491 | 4 | <1 | 9 |
| 21 | 2-mercaptopropanoic acid | AF II | 2130 | 1057 | 4 | <1 | 5 |
| 22 | 2-methyl-3-furyl 2-oxo-3-pentyl disulfide | AF I | 2145 † | 1561 | 8 | <1 | 9 |
| 23 | 3-hydroxy-4,5-dimethyl-2(5H)-furanone (sotolon)^g | AF II | 2153 | 1112 | 16 | <1 | |
| 24 | **5-acetyl-2,3-dihydro-1,4-thiazine** | NB VI | 2177 | 1374 | <1 † | **256** | |
| 25 | 3-mercaptopropanoic acid | AF II | 2193 | 1057 | 16 | <1 | 5 |
| 26 | 5-propionyl-2,3-dihydro-1,4-thiazine | NB VI | 2235 | 1456 | <1 | 8 | |
| 27 | 2-furfuryl 2-methyl-3-furyl disulfide | AF I / NB I | 2323 | 1624 | 16 | 4 | 10 |
| 28 | 2-furfuryl 2-oxo-3-butyl disulfide | AF I | 2352 | 1577 | 4 | <1 | |
| 29 | 2-furfuryl 2-oxo-3-pentyl disulfide | AF I | 2385 | 1649 | 4 | <1 † | |
| 30 | bis(2-furfuryl) disulfide | AF I / NB I | 2465 | 1673 | 16 | 4 | 10 |

**The † cells, and why each reading is forced.** (4) SE-54 index printed `e800`; the compound is
1-mercapto-2-propanone (MW 90) eluting before every other entry, and the column has no value below
820 elsewhere — read **`<800`**. (5) FFAP index printed `14.31`, which is not a retention index; the
sequence runs 1357 → ? → 1512 and Table 2 independently prints **1431** for 2-furfurylthiol on FFAP
— read **1431**. (5, 29) FD II printed `<l` and `<I` — read **`<1`**, the only value in that column
below 2. (12) FD I printed `(1` — read **`<1`**. (18) fraction printed `AFJf-NBI`; every other
two-fraction entry is written `AF I/NB I` and this compound is quantified in both columns — read
**AF I / NB I**. (22) FFAP index printed `2 145` with a space — read **2145** (it must lie between
2130 and 2153). (24) FD I printed `41`; the running text says compound 24 "showed the highest odor
activity in the neutral basic fraction" and lists it among the compounds identified **only** in NB,
and 41 is not a power of two — read **`<1`**. **The abstract independently confirms the six highest
FD factors**: 2-furfurylthiol, 3-mercapto-2-pentanone, 2-methyl-3-furanthiol,
5-acetyl-2,3-dihydro-1,4-thiazine, 3-mercapto-2-butanone and bis(2-methyl-3-furyl) disulfide, which
is exactly the set {1024, 512, 256, 256, 128, 128} in the table as re-typed.

**Footnote e key, as printed:** (1) Mulders 1973; (2) van den Ouweland & Peer 1975; (3) Mussinan &
Katz 1973; (4) Whitfield et al. 1988; (5) Martin 1988; (6) Farmer et al. 1989; (7) Güntert et al.
1990; (8) Sakaguchi & Shibamoto 1978; (9) Whitfield et al. 1993; (10) Farmer & Patterson 1991;
(11) Silwar 1992; (12) Güntert et al. 1993. **Compounds 13, 14, 23, 24, 26, 28 and 29 carry no
reference — they had not been reported in such systems before**, and 24 and 26 not anywhere in
Maillard chemistry.

### Table 2. "Results of Static Headspace Analysis/Olfactometry (SHO) of the Cysteine/Ribose Model Solution"

Footnote a: identification by mass spectrum, retention index, odour quality and odour activity
against the reference compound (cf. Table 1). Footnote b: "Lowest headspace volume required to
perceive the odorant at the sniffing port." Footnote c: "Calculated by dividing the largest volume
analyzed (20 mL) by the lowest headspace volume required to perceive the odorant at the sniffing
port." Footnote d: "3-Mercapto-2-pentanone contained a smaller proportion of its isomer
2-mercapto-3-pentanone." Footnote e: "MS signals were too weak for an unequivocal interpretation."

| odorant^a | RI on FFAP | vol^b (mL) | FD factor^c |
|---|---:|---:|---:|
| hydrogen sulfide^d | <900 | 5 | 4 |
| methyl mercaptan | <900 | 2.5 | 8 |
| **ethyl mercaptan** | <900 | **0.16** | **125** |
| 2,3-butanedione (diacetyl) | 968 | 20 | 1 |
| 3-mercapto-2-butanone | 1282 | 2.5 | 8 |
| **2-methyl-3-furanthiol** | 1300 | **0.08** | **250** |
| 3-mercapto-2-pentanone^d | 1347 | 0.64 | 32 |
| dimethyl trisulfide | 1356 | 20 | 1 |
| **2-furfurylthiol** | 1431 | **0.02** | **1000** |
| unknown (sulfury) | 1438 | 10 | 2 |
| unknown (sulfury) | 1511 | 5 | 4 |
| 2-methyl-3-thiophenethiol | 1559 | 0.64 | 32 |
| **2-thenyl mercaptan** | 1682 | **0.16** | **125** |
| 2-acetyl-2-thiazoline^e | 1720 | 20 | 1 |
| 5-acetyl-2,3-dihydro-1,4-thiazine^e | 2177 | 20 | 1 |

**Self-check (mine): every FD reproduces as 20/vol** — 20/5 = 4, 20/2.5 = 8, 20/0.16 = 125,
20/20 = 1, 20/0.08 = 250, 20/0.64 = 31.25 → 32, 20/0.02 = 1000, 20/10 = 2. **Table 2 is internally
consistent to the last row**, which is a strong indication that this table's text layer is
uncorrupted. Note the footnote-d marker is printed against **hydrogen sulfide** in the OCR but the
footnote text is about 3-mercapto-2-pentanone; the marker plainly belongs to the latter (Flags 6).

**Fifteen odorants detectable in a 20 mL headspace sample, twelve of them identified.** Five —
hydrogen sulfide, methyl mercaptan, ethyl mercaptan, diacetyl, dimethyl trisulfide — were **not**
found by AEDA of the solvent extract, "undoubtedly due to the high volatility of these four
odorants causing losses during concentration of the solvent extracts" (the paper writes "four"
where it has just listed five). Conversely compounds 4, 6, 8, 9, 12, 15-23 and 25-30 of Table 1
were not detectable in the headspace, "possibly due to their lower vapor pressure".

### Table 3. "Odor Thresholds of Selected Sulfur-Containing Odorants in the Cysteine/Ribose System"

Footnote a: "Determined by HRGC-O/olfactometry as recently described (Ullrich and Grosch, 1987).
**(E)-2-Decenal (odor threshold = 2.7 ng/L of air) instead of hexanal was used as the reference
odorant.**" Footnote b: "Adapted from the literature (Gasser and Grosch, 1990b)." Footnote c:
"Threshold in water = **0.6 µg/L (retronasally)** and **1.7 µg/L (nasally)** (**unpublished
results**)."

| odorant | odour threshold (ng/L of air) |
|---|---|
| 2-methyl-3-thiophenethiol | 0.0032-0.0128 |
| 2-thenyl mercaptan | 0.003-0.012 |
| **5-acetyl-2,3-dihydro-1,4-thiazine**^c | **0.02-0.08** |
| 2-acetyl-2-thiazoline | 0.02-0.08 |
| 3-mercapto-2-pentanone | 0.05-0.2 |
| 3-mercapto-2-butanone | 0.2-0.8 |
| 1-mercapto-2-propanone | 1.7-6.8 |
| 2-furfurylthiol | 0.0025-0.01^b |
| 2-methyl-3-furanthiol | 0.0025-0.01^b |
| bis(2-methyl-3-furyl) disulfide | 0.0006-0.0024^b |
| bis(2-furfuryl) disulfide | 0.00015-0.0006^b |
| 2-furfuryl 2-methyl-3-furyl disulfide | 0.0004-0.0016 |
| 2-methyl-3-furyl 3-oxo-2-butyl disulfide | 0.01-0.04 |
| 2-furfuryl 3-oxo-2-butyl disulfide | 0.004-0.016 |

Every entry is a **factor-of-four band** (the HRGC-O bracketing method gives the last dilution
perceived and the first not perceived), except the two adapted values which are factor-of-four
bands as well. The abstract quotes the thiazine as "**0.06 ng/L of air**", a point value inside its
own band (Flags 7).

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| odorants in the acidic fraction | **21**, FD range **4-1024** | p. 2190 |
| odorants in the neutral/basic fraction | **12** | p. 2190 |
| total intense odorants tabulated | **29** (abstract) / 30 rows in Table 1, two of them "unknown" | abstract, Table 1 |
| odorants detected in a 20 mL headspace | **15**, of which **12 identified** | p. 2191 |
| disulfides identified | **8** (18, 20, 22, 27-30 and the mixed pair) | p. 2191 |
| **MFT in diethyl ether at 6 C** | "after **10 days** storage ... **more than 50 %** of the thiol had been converted into the corresponding disulfide" | p. 2191 — **cited forward to Hofmann et al. 1995b, i.e. `hofmann1996` on disk; NOT measured here** |
| disulfides in the headspace | **none detected** | p. 2191 |
| pH screen | pH 5.0 "sulfury, meatlike"; pH 7.0 "caramel, burnt" — **data not shown** | p. 2189 |
| ratio screen | 1:1 "rubbery, sulfury, H2S-like"; 10:1 "burnt and caramel-like"; **3:1 gives "the most prominent roasty, meatlike odor"** — data not shown | p. 2189 |
| synthesis yields (reference compounds, not the model reaction) | 3-bromo-2-methylthiophene 69 %; 2-methyl-3-thiophenethiol 58 %; 2-thenyl mercaptan 39 %; 3-mercapto-2-butanone 10 %; 1-mercapto-2-propanone 71 % | Experimental |
| 2-thenyl mercaptan | not detected when cysteine was omitted (Güntert 1993), "suggesting that cysteine is the main precursor of this odorant" | p. 2192 |
| thiamine/cysteine vs ribose/cysteine | thiamine/cysteine is "a much more effective precursor mixture in the formation of 2-methyl-3-furanthiol" (Grosch & Zeiler-Hilgart 1993) — **cited, not measured** | p. 2192 |
| mass spectra printed | 2-thenyl mercaptan (Fig. 1), five mixed disulfides (Fig. 2), 5-acetyl-2,3-dihydro-1,4-thiazine (Fig. 3, m/z 143 base), 5-propionyl-2,3-dihydro-1,4-thiazine (m/z 157 (100), 100 (64), 72 (42), 45 (40), 57 (34), 101 (22)), plus the synthesised standards | throughout |

**Everything else in this paper is a proposed pathway (Scheme 1) or a mass spectrum.** There is no
plot with a numeric axis anywhere.

### Arithmetic on the printed numbers (all mine)

**1. The two dilution ladders disagree about which thiol matters, and by how much.**

| odorant | AEDA FD (acidic extract) | SHO FD (headspace) | ratio SHO/AEDA |
|---|---:|---:|---:|
| 2-furfurylthiol | 1024 | 1000 | 1.0 |
| 2-methyl-3-furanthiol | 256 | 250 | 1.0 |
| 3-mercapto-2-pentanone | 512 | 32 | **0.06** |
| 3-mercapto-2-butanone | 128 | 8 | **0.06** |
| 2-methyl-3-thiophenethiol | 16 | 32 | 2.0 |
| **2-thenyl mercaptan** | 8 | **125** | **15.6** |
| **ethyl mercaptan** | not detected | **125** | **>= 125** |
| 2-acetyl-2-thiazoline | 64 (in NB) | 1 | 0.016 |
| 5-acetyl-2,3-dihydro-1,4-thiazine | 256 (in NB) | 1 | 0.004 |

The two furan/furfuryl thiols agree to within the method's own resolution; **everything else moves
by one to two orders of magnitude**, in both directions. The two ladders measure different things
(an extract's dilution against a headspace's dilution, i.e. a concentration ladder against a
partial-pressure ladder), so this table is a **volatility map, not a disagreement** — but it is the
sharpest available warning that **an FD factor from a solvent extract cannot be read as a
concentration ranking**, which is precisely how such tables are usually used.

**2. What the FD factors do and do not bound.** An FD factor is the last dilution at which the
odorant is still smelled, so **FD ≈ concentration / odour threshold** in whatever medium was
sniffed. Combining Table 1's FD with Table 3's threshold gives an *odour activity*, not a
concentration, and the two are only convertible if the extraction recovery is known — which it is
not, for any compound, in this paper. **No concentration can be recovered from this paper by any
arithmetic**, and this dossier does not attempt one.

**3. The assessor spread.** Table 1's footnote d says the duplicate analyses by two assessors
"differed to not more than two FD factors", i.e. up to **a factor of 4** on a base-2 ladder. So the
gap between 2-furfurylthiol (1024) and 2-methyl-3-furanthiol (256) — a factor of 4 — is at the edge
of the method's resolution, and the abstract's ordering of the six top odorants should not be read
as a strict ranking below the top pair.

**4. The pot, in the module's units.** Ribose **100 mmol/L**, cysteine **33 mmol/L**, phosphate
**500 mmol/L**, pH 5.0, ramped 20 → 145 C over 20 min. **This is the same charge van Seeventer 2001
reproduces** (their footnote a: "According to the procedure of Hofmann and Schieberle (2), a
phosphate-buffered (0.5 M, pH 5.0) solution of D-ribose (100 mM) and L-cysteine (33 mM)") — but van
Seeventer changed the thermal history to **a 10 min ramp to 130 C plus a 20 min hold**. The two pots
are therefore *the same chemistry with different cooks*, and the 55 % / 75 % precursor conversions
that this cluster's van Seeventer dossier records belong to the 130 C version, not to this one.

**5. Which of this paper's odorants the network carries.** Of the six highest-FD compounds, five are
in the sulfur network (`MFT`, `FFT`, `MP3P`, the mercaptobutanone family, and the MFT dimer) and one
— 5-acetyl-2,3-dihydro-1,4-thiazine — is **not modelled at all** and has the highest FD factor of
anything in the neutral/basic fraction. `sulfur.py`'s OUT_OF_SCOPE block names the thiazole family
beyond 2-acetylthiazole as out of scope for this wave; **the 1,4-thiazines are a different ring
system and are not named in that block at all** (Flags 8).

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Keyed: `2_furfurylthiol`, `2_methyl_3_furanthiol`,
`bis_2_methyl_3_furyl_disulfide`, `mercapto_2_propanone` (= 1-mercapto-2-propanone), `norfuraneol`,
`hdmf`, `hydrogen_sulfide`, `methanethiol` (= methyl mercaptan), `dimethyl_trisulfide`,
`2_3_butanedione`. **Not keyed, and appearing in Table 1, Table 2 or Table 3:** 3-mercapto-2-butanone,
3-mercapto-2-pentanone, 2-mercapto-3-pentanone, 2-methyl-3-thiophenethiol,
**2-acetylthiazole** (which the B2 objective already scores as a species — see Flags 8),
2-acetyl-2-thiazoline, 2-formylthiophene, acetylthiophene, 2-methyltetrahydrothiophen-3-one,
2-(hydroxymethyl)thiophene, 2-thenyl mercaptan, ethyl mercaptan, 3-hydroxy-2H-pyran-2-one, sotolon,
2- and 3-mercaptopropanoic acid, **5-acetyl-2,3-dihydro-1,4-thiazine**,
5-propionyl-2,3-dihydro-1,4-thiazine, bis(2-furfuryl) disulfide and the four mixed disulfides.

**There is no rate, no yield and no concentration in this paper.** Every row below is a dilution
factor or a threshold. Shared conditions for Tables 1 and 2: **ribose 100 mmol/L + cysteine
33 mmol/L in 100 mL of 0.5 mol/L phosphate at pH 5.0, autoclave, 20 → 145 C over 20 min**; Table 1
by AEDA on a diethyl-ether extract split into acidic and neutral/basic fractions, two assessors in
duplicate; Table 2 by static headspace at 40 C from a 240 mL vessel.

| quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| FD, 2-furfurylthiol | **1024** (acidic extract) | dilution steps | AEDA | none | Table 1 no. 5 | **threshold** |
| FD, 3-mercapto-2-pentanone | 512 | " | AEDA | none | Table 1 no. 3 | threshold |
| FD, 2-methyl-3-furanthiol | 256 | " | AEDA | none | Table 1 no. 2 | threshold |
| FD, 5-acetyl-2,3-dihydro-1,4-thiazine | 256 (neutral/basic) | " | AEDA | none | Table 1 no. 24 | threshold |
| FD, 3-mercapto-2-butanone | 128 | " | AEDA | none | Table 1 no. 1 | threshold |
| FD, bis(2-methyl-3-furyl) disulfide | 128 (acidic) / 16 (neutral/basic) | " | AEDA | none | Table 1 no. 18 | threshold — **and see Flags 3: the authors call it an artefact** |
| FD, norfuraneol | 64 | " | AEDA | none | Table 1 no. 19 | threshold |
| FD, 2-formylthiophene / 2-acetyl-2-thiazoline | 64 each (neutral/basic) | " | AEDA | none | Table 1 nos. 9, 11 | threshold |
| FD, HDF | 32 | " | AEDA | none | Table 1 no. 17 | threshold |
| the other 21 rows of Table 1 | 2 to 16 | " | AEDA | none | Table 1 | threshold |
| SHO FD, 2-furfurylthiol | **1000** (0.02 mL) | dilution steps | static headspace, 40 C | none | Table 2 | threshold |
| SHO FD, 2-methyl-3-furanthiol | 250 (0.08 mL) | " | " | none | Table 2 | threshold |
| SHO FD, ethyl mercaptan | **125** (0.16 mL), **absent from the extract** | " | " | none | Table 2 | threshold |
| SHO FD, 2-thenyl mercaptan | **125** (0.16 mL) against FD 8 in the extract | " | " | none | Table 2 | threshold |
| SHO FD, hydrogen sulfide | 4 (5 mL) | " | " | none | Table 2 | threshold |
| SHO, disulfides | **none detected** | — | " | — | p. 2191 | threshold (an absence, weakened by volatility — Flags 3) |
| odour threshold in air, 2-furfurylthiol | 0.0025-0.01 | ng/L air | HRGC-O, (E)-2-decenal reference | — | Table 3 | **threshold** — adapted from Gasser & Grosch 1990b, **not measured here** |
| odour threshold in air, 2-methyl-3-furanthiol | 0.0025-0.01 | ng/L air | " | — | Table 3 | threshold — adapted, not measured here |
| odour threshold in air, bis(2-methyl-3-furyl) disulfide | 0.0006-0.0024 | ng/L air | " | — | Table 3 | threshold — adapted |
| odour threshold in air, bis(2-furfuryl) disulfide | **0.00015-0.0006** — the lowest in the paper | ng/L air | " | — | Table 3 | threshold — adapted |
| odour threshold in air, 2-methyl-3-thiophenethiol | 0.0032-0.0128 | ng/L air | " | — | Table 3 | **threshold — measured here** |
| odour threshold in air, 2-thenyl mercaptan | 0.003-0.012 | ng/L air | " | — | Table 3 | threshold — measured here |
| odour threshold in air, 3-mercapto-2-pentanone | 0.05-0.2 | ng/L air | " | — | Table 3 | threshold — measured here |
| odour threshold in air, 3-mercapto-2-butanone | 0.2-0.8 | ng/L air | " | — | Table 3 | threshold — measured here |
| odour threshold in air, 1-mercapto-2-propanone | 1.7-6.8 | ng/L air | " | — | Table 3 | threshold — measured here, **the highest (weakest) in the paper** |
| odour threshold in air, 2-acetyl-2-thiazoline | 0.02-0.08 | ng/L air | " | — | Table 3 | threshold |
| odour threshold in air, 5-acetyl-2,3-dihydro-1,4-thiazine | 0.02-0.08 (abstract: 0.06) | ng/L air | " | — | Table 3 | threshold |
| odour threshold in **water**, 5-acetyl-2,3-dihydro-1,4-thiazine | **0.6 retronasal / 1.7 nasal** | µg/L water | — | — | Table 3 footnote c | threshold — **"unpublished results", no method given** |
| the pot | ribose 100, cysteine 33, phosphate 500 | mmol/L, pH 5.0 | 20 → 145 C over 20 min, 100 mL, autoclave | — | Experimental + Table 1 footnote a | **level_only** — the definition of the reference system |
| MFT → disulfide in diethyl ether | > 50 % in 10 days at 6 C | % | ether, 6 C, dark not stated | — | p. 2191, **cited to Hofmann et al. 1995b** | **level_only** — do not attribute to this paper; it is `hofmann1996` |

### Can any of this be put on the same basis as a shipped constant? No, and the reason is categorical.

An FD factor is the number of two-fold dilutions an odorant survives at a sniffing port. To turn one
into a concentration you need the odour threshold **in the sniffed medium** and the **recovery of
the isolation**, and this paper reports the threshold only in air and the recovery not at all.
An odour threshold, in turn, is exactly what the house rules call a `threshold` — usable for
odour-activity work, never as an amount. The correct use of this paper is:

- **as the specification of the reference pot** (charge, buffer, pH, thermal history), which is what
  `vanseeventer2001` and both `hofmann1998` papers build on and what the module's 145 C reference
  descends from;
- **as a species-selection argument** — which odorants a ribose/cysteine pot at 145 C actually
  presents to a nose, and in what order;
- **as fourteen thresholds in air**, ten of them this laboratory's own, which is a genuine
  contribution to `desirable_targets.yml` and to any OAV work, **as long as they are labelled as air
  thresholds** (Flags 5);
- **as a methodological warning** about disulfides in solvent extracts, and about the loss of
  ethyl mercaptan and 2-thenyl mercaptan during extract concentration.

## 5. Flags

1. **The 1995 scan's text layer is corrupted in Table 1 and this dossier lists every repair.** Eight
   cells were unreadable as printed; all eight are re-typed above with the reading forced by the
   running text, by the abstract's independent list of the six top FD factors, by Table 2's
   independent printing of the FFAP retention index 1431, and by the fact that every FD factor must
   be a power of two. **No value was guessed.** Two rows carry an unresolved oddity: compounds **21
   and 25** (2- and 3-mercaptopropanoic acid) both print SE-54 retention index **1057** although
   their FFAP indices differ by 63 units — one of the two is probably an OCR duplication, and the
   page image should be checked before either index is used.
2. **The "20 min" is the RAMP, not a hold.** "raising the temperature within 20 min from 20 to
   145 °C" is the whole thermal description; no isothermal period at 145 C is stated, and Table 1's
   footnote a compresses this to "reacted for 20 min". **The pot may never dwell at 145 C at all.**
   The module's reference temperature and every "145 C / 20 min" condition string in the corpus that
   traces to this paper should carry that caveat; the sequels (`hofmann1998`) should be checked for
   whether they hold or ramp, and this dossier does not assume they match.
3. **The paper's own disulfides are declared artefacts by their authors, and the counter-evidence is
   weak in the same direction.** The argument is internal (a neutral compound enriched into the
   acidic fraction) plus a forward citation to a not-yet-published ether-storage experiment, plus a
   headspace non-detect that the paper elsewhere explains away by low vapour pressure for exactly
   this molecular-weight class. **A further complication nobody has raised**: the acidic-fraction
   thiols were eluted from the phenylmercury affinity column with **10 mmol of dithiothreitol**, a
   reducing agent, which would *reduce* disulfides rather than create them — yet the disulfides come
   out higher in that fraction. Whatever is happening in AF I (where the disulfides sit, i.e. the
   **unbound** cut, not the DTT cut), the fractionation chemistry here is not simple, and the
   artefact claim, while plausible, is not demonstrated in this paper.
4. **No recovery, no internal standard, no quantification of any kind.** There is no isotope
   dilution here (contrast Hofmann's own 1998 and 2002 papers, which use SIDA), no external standard
   and no calibration. That is appropriate for an AEDA screen and it means the paper supports no
   level.
5. **The thresholds are in AIR, by a bracketing HRGC-O method, against an unusual reference.**
   (E)-2-Decenal at 2.7 ng/L replaces the customary hexanal, which makes these values a scale of
   their own; four of the fourteen are **adapted from Gasser & Grosch 1990b and were not measured
   here** (2-furfurylthiol, 2-methyl-3-furanthiol and the two symmetrical disulfides — footnote b).
   The only threshold **in water** in the paper is the thiazine's, and it is marked "unpublished
   results" with no method, no panel size and no dilution scheme. `matrix_oav` needs matrix
   thresholds; **these are air thresholds and are not that.**
6. **A misplaced footnote marker in Table 2.** Footnote d, "3-Mercapto-2-pentanone contained a
   smaller proportion of its isomer 2-mercapto-3-pentanone", is attached in the OCR to the **hydrogen
   sulfide** row; it plainly belongs to the 3-mercapto-2-pentanone row, where the identical footnote
   appears in Table 1 (footnote f). More substantively: **the paper's 3-mercapto-2-pentanone is not
   pure** — it contains its isomer, in both tables, in an unstated proportion. Since the isomer
   split is the diagnostic the module uses to separate the norfuraneol route from the ribose route
   (`cerny_isomer_split`), **this paper's FD 512 and its 0.05-0.2 ng/L threshold are both for a
   mixture**, and Cerny 2003's finding that 2-mercapto-3-pentanone is absent from a ribose/cysteine
   pot at 95 C sits directly against it. Cerny explicitly cites this paper's 145 C conditions as the
   contrast case ("Higher temperatures and shorter reaction times seem to favor a different reaction
   mechanism leading to both isomers").
7. **The abstract's 0.06 ng/L for the thiazine is not in the table.** Table 3 gives the band
   0.02-0.08; 0.06 is a point value inside it that appears nowhere else. Quote the band.
8. **Two registry gaps that are already biting.** (i) **`2_acetylthiazole` has no id in
   `data/keys/compounds.yml`**, although the B2.x objective scores a row named `zhou_pH7_ACTZ`
   against a species `ACTZ` and this paper lists 2-acetylthiazole at FD 8 in the neutral/basic
   fraction. (ii) **5-acetyl-2,3-dihydro-1,4-thiazine and its propionyl homologue are absent from
   the registry and from the network**, and the acetyl compound has the highest FD factor of any
   neutral/basic odorant in the reference pot together with an odour threshold of 0.02-0.08 ng/L in
   air. `sulfur.py`'s OUT_OF_SCOPE block names 2-acetyl-1-pyrroline and "the thiazole family beyond
   2-acetylthiazole"; **the 1,4-thiazines are a distinct ring system and are named nowhere**, so
   this is an unrecorded omission rather than a declared one. Also unkeyed and prominent here:
   3-mercapto-2-pentanone, 2-mercapto-3-pentanone, 3-mercapto-2-butanone, 2-methyl-3-thiophenethiol,
   2-thenyl mercaptan, ethyl mercaptan, sotolon and bis(2-furfuryl) disulfide.
9. **What to request from the authors**: (i) any quantification at all of this pot — the paper's own
   sequels supply some, but not for the thiazines or the mercaptopropanoic acids; (ii) the pH and
   ratio screen data, which are described only as odour impressions ("data not shown") and which
   would be a genuine pH row for a lane that has few; (iii) the thermal profile of the autoclave ramp
   and whether a hold followed; (iv) a headspace or SPME measurement of the disulfides, which would
   settle the artefact question in the reference pot and is the one number this paper's chemistry
   most needs; (v) the method behind the thiazine's water thresholds.
10. **What this paper does not contain**: any concentration; any yield; any rate; any barrier; any
    time course; any second temperature; any pH series with data; any mass balance; any recovery
    figure; any replicate dispersion beyond "not more than two FD factors"; any measurement of ribose
    or cysteine consumed; and any supplementary material.
