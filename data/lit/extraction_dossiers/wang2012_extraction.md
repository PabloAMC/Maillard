# Wang, Yang & Song 2012 — EXTRACTION
### *"Key meat flavour compounds formation mechanism in a glutathione–xylose Maillard reaction"*
### Five unbuffered pots (pH 4.8, 132 °C, 90 min), all carrying 1.2 mM thiamine and a 1:1 [13C5]/[12C] xylose mix. Their sulfur compounds are IS-normalised semi-quantities. The isotope table says MFT is >99 % non-xylose and furfural is intact xylose. FFT reads as only ~40 % xylose-derived, which contradicts both Cerny papers. There is no H2S measurement, no time course and no kinetics. Read 2026-10-27.

**Full read of `data/articles/wang2012.pdf`** (6 pages, publisher PDF with a clean text layer). Every
table was checked against a 150-dpi raster of the page. Figures 2 and 3 are reaction schemes with no
data. Related dossiers: `cerny2007b_extraction.md` (the same CAMOLA-type question at 100 mM, which
this paper cites), `cerny2007_extraction.md` (the pH ladder), `zheng1994_extraction.md` (the GSH → H2S
rate constants the model uses).

**Provenance:** **[M]** measured, printed as a number · **[Q]** qualitative · **[Z]** derived by me
from the paper's own numbers, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Ran Wang, Chao Yang, Huanlu Song* |
| Title | "Key meat flavour compounds formation mechanism in a glutathione–xylose Maillard reaction" |
| Venue | *Food Chemistry* **131** (2012) **280–285** |
| DOI | 10.1016/j.foodchem.2011.08.079 |
| Affiliation | Beijing Key Laboratory of Flavor Chemistry, Beijing Technology and Business University, China |
| Dates | received 16 May 2011, revised 14 June 2011, accepted 26 Aug 2011, online 1 Sep 2011 |
| Funding | NSFC 31071510; Beijing PHR200906110 |
| PDF character | Table 1 (recipes), Table 2 (30 volatiles × 5 pots), Table 3 (6 isotope splits), Table 4 (DKP, 5 values); Figs. 2–3 are schemes only |

---

## 1. THE ONE-PARAGRAPH ANSWER

The pots contain xylose (26 mM, half of it [13C5]) and thiamine HCl (1.19 mM), plus one sulfur
co-reactant: glutathione (A, 7.4 mM), Glu + Cys + Gly (B, Cys 9.6 mM), Cys + Gly (C, Cys 12.1 mM),
Glu + Cys-Gly (D) or Cys-Gly (E). Each is 25 mL of distilled water at "pH 4.8", unbuffered, heated
in a Parr steel reactor at 132 °C for 90 min. The free-cysteine pots give the most sulfur
volatiles. At 7.4 mM GSH, pot A gives MFT 0.021, FFT 0.049, thiophene 0 and DMDS 0.279 (Table 2
units, read as mg/L). At 9.6 mM free Cys, pot B gives 0.082 / 0.066 / 0.312 / 0.157. The
Cys-Gly pots give **no FFT and no thiophene**, but the most furfural (D 0.575) and the most MFT
(D 0.179). The CAMOLA-type label (Table 3, pot not stated; see [!] in §4) gives four results.
Furfural is 52:48 unlabelled:labelled, i.e. **intact xylose**. **MFT is >99:<1, so <~2 %
xylose-derived and essentially all thiamine** at these dilute loadings. **FFT is 80:20**, which
would mean only ~40 % intact xylose. That conflicts with the >99 % xylose found by Cerny 2007b
and Cerny & Briffod 2007 under every condition, and I treat it as unreliable. Thiophene is 88:12,
attributed to cysteine C2 units. For GSH, the paper proposes hydrolysis of the γ-glutamyl bond to
pyroglutamic acid + Cys-Gly, then cyclo(Cys-Gly) (DKP), then desulfurisation to
3-methyl-2,5-piperazinedione. The evidence is LC-TOF ions whose printed masses do not match the
assigned formulae, plus an HPLC-UV DKP number that equals **0.035 mol % of the GSH** in pot A.
**Nothing here is a rate constant, a yield-vs-time curve or an H2S measurement.**

---

## 2. METHODS [M]

| item | as printed | note |
|---|---|---|
| Recipes | Table 1 (§3); mg in **25 mL distilled water** | [Z] molarities below |
| pH | "dissolved in distilled water (25 ml, pH 4.8)" | **no buffer**; how pH 4.8 was set, and whether it applies to every pot, is not stated; final pH not reported [!] |
| Vessel | "High-pressure Stainless Reactor (Parr instrument Co.)"; reactor volume, headspace, atmosphere and stirring not stated | [!] headspace unknown |
| Heating | **132 °C, 90 min** ("the reaction parameters had been optimised" after Lv & Song 2010); come-up time not stated | single time point |
| Quench / storage | cooled rapidly under running tap water, then **stored at 4 °C for 24 h** "to afford the after reactions" | [!] thiol/disulfide chemistry continues during storage |
| Volatile isolation | purge & trap (Tekmar Stratum): **10 mL** reaction product + **1 µL internal standard (2-methyl-3-heptanone, 1.68 µg/µL in methanol)**, sparge 50 °C, N2 40 mL/min, purge 10 min, desorb 250 °C 2 min, bake 260 °C 10 min | IS in sample = 1.68 µg / 10 mL = **0.168 mg/L** [Z] |
| GC-O-MS | Agilent 7890A + 7000B triple quad, Gerstel Sniffer 9000; DB-5MS UI 30 m × 0.25 mm × 0.25 µm; He 1.2 mL/min; 40 °C (3 min) → 5 °C/min → 200 °C → 15 °C/min → 260 °C (5 min); EI 70 eV, m/z 35–550; three GC-O panellists | |
| Identification | MS + RI + odour vs authentic references (Table 2 fn a) and NIST 2.0 | [!] some identities implausible (§8) |
| Replicates | Table 2: "means of at least three assays; maximum SD + 10%" | Table 3 and 4 replication not stated |
| DKP | HPLC-UV (ZORBAX SB-C18, 1–30 % MeCN/0.1 % TFA over 20 min, 220/254 nm), external standard 0.5–10.0 mg/L | identity by retention time only |
| Non-volatiles | LC-ESI-TOF-MS (Agilent 6210), C18, 15 % MeOH, m/z 100–1000, lock masses 121.0508 / 149.0233 / 922.0097 | qualitative |
| Isotope | [13C5]xylose **98 % enrichment** (Sigma), mixed **1:1 by mass** with unlabelled xylose in every pot | [Z] mole fraction labelled = 0.492 |

**Molar charges [Z]** (25 mL; MW thiamine·HCl 337.27, GSH 307.32, Cys 121.16, Gly 75.07, Glu 147.13,
Cys-Gly 178.21, xylose 150.13 / [13C5] 155.09):

| pot | xylose total (mM) | thiamine (mM) | sulfur co-reactant (mM) | other (mM) |
|---|---:|---:|---|---|
| A | 26.2 | 1.19 | **GSH 7.42** | — |
| B | 26.2 | 1.19 | **Cys 9.64** | Glu 7.42, Gly 7.62 |
| C | 26.2 | 1.19 | **Cys 12.15** | Gly 9.80 |
| D | 26.2 | 1.19 | **Cys-Gly 6.67** | Glu 7.42 |
| E | 26.2 | 1.19 | **Cys-Gly 12.79** | — |

[!] The pots are **not equimolar in sulfur**. B carries 1.30× the S of A, C 1.64×, D 0.90×,
E 1.72×. Glu in B and D is equimolar with A's GSH, so B and D look like attempts at "GSH
constituents", but the Cys and Gly masses do not match that. Cys-Gly is listed as "(90%
enrichment)", probably meaning purity, in which case D and E are ~10 % lower still.

---

## 3. TABLES AS PRINTED

### Table 1. "Model reactions." Amount (mg)

| | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|
| Xyl | 50 | 50 | 50 | 50 | 50 |
| [13C5] Xyl | 50 | 50 | 50 | 50 | 50 |
| Thiamine HCl | 10 | 10 | 10 | 10 | 10 |
| GSH | 57 | | | | |
| Glu | | 27.3 | | 27.3 | |
| Cys | | 29.2 | 36.8 | | |
| Gly | | 14.3 | 18.4 | | |
| Cys-Gly | | | | 29.7 | 57 |

### Table 2. "Volatiles of the model reaction." Relative concentration (unit printed as "µg µl", see [!])

Footnotes: a, identified vs authentic references by MS, RI and odour; b, RI from C7–C22 n-alkanes;
c, odour at the Sniffer 9000; d, "Data are means of at least three assays; maximum SD + 10%".
"–" = not detected.

| No. | compound | RI | odour | A | B | C | D | E |
|---:|---|---:|---|---:|---:|---:|---:|---:|
| | *Sulphur-containing compounds* | | | | | | | |
| 1 | Dimethyldisulphide | 751 | Sulphury | 0.279 | 0.157 | 0.242 | 0.515 | 0.227 |
| 2 | **2-Methyl-3-furanthiol** | 865 | Meaty, vitamin | **0.021** | **0.082** | **0.118** | **0.179** | **0.034** |
| 3 | **2-Furfurylthiol** | 908 | Roasted, coffee | **0.049** | **0.066** | **0.170** | – | – |
| 4 | Thiophene | 660 | Scallion | – | 0.312 | 0.339 | – | – |
| 5 | 2-Methyl-thiophene (printed "2-Me thy1-thiophene") | 774 | Onion | – | – | 0.111 | – | – |
| 6 | Dimethyl trisulphide | 973 | Mint | – | – | 0.052 | 0.392 | 0.052 |
| 7 | 2-Pentyl-thiophene | 1168 | Fatty | – | 0.005 | 0.006 | – | – |
| | *Aldehydes, Ketones* | | | | | | | |
| 8 | 2-Methyl-3-pentanone | 712 | – | 0.094 | 0.103 | – | – | 0.121 |
| 9 | 6-Hydroxy-2-hexanone | 760 | – | 0.161 | 0.114 | – | – | – |
| 10 | Hexanal | 799 | Green, grass | 0.089 | 0.063 | 0.051 | 0.112 | 0.073 |
| 11 | **Furfural** | 831 | Almond, pungent | **0.152** | **0.066** | **0.089** | **0.575** | **0.169** |
| 12 | Heptanal | 890 | Cheesy, fatty | 0.015 | – | – | – | 0.022 |
| 13 | 3-Octanone | 982 | Rusty | 0.025 | – | – | – | – |
| 14 | 3-Octanol | 990 | Earthy | 0.094 | – | – | – | – |
| 15 | Octanal | 1002 | Fruity, sour | – | – | – | – | 0.011 |
| 16 | Nonanal | 1078 | Green, stale | – | – | – | – | 0.006 |
| | *Heterocyclic compounds* | | | | | | | |
| 17 | 3-Methyl-2,5-piperazinedione | 953 | – | 0.019 | – | – | – | – |
| 18 | 3,6-Dimethylpiperazine-2,5-dione | 1612 | – | 0.661 | – | – | – | – |
| 19 | 3-Ethyl-1H-pyrrole | 811 | – | – | 0.075 | 0.054 | – | – |
| 20 | 2-Pentylfuran | 987 | Fruity | – | 0.021 | 0.020 | – | – |
| 21 | D-Limonene | 994 | Sweet | – | 0.245 | 0.188 | 0.277 | – |
| 22 | 4-Methyl-5-thiazoleethanol | 1250 | Roasted | – | – | – | – | 0.025 |
| | *Others* | | | | | | | |
| 23 | 1,2-Ethanediol monoformate | <600 | – | 0.102 | – | 0.104 | 0.135 | – |
| 24 | 1,5-Pentanediol dinitrate | <600 | – | 0.606 | – | – | – | 1.038 |
| 25 | 3-Methyl-1-butanol acetate | 878 | Solvent | 0.046 | 0.099 | – | – | – |
| 26 | Butylbenzene | 1062 | – | 0.030 | 0.022 | 0.023 | – | 0.012 |
| 27 | 1-Methyl-2-(1-methylethyl)-benzene | 1031 | – | – | 0.070 | 0.040 | 0.060 | 0.008 |
| 28 | 1-Methyl-4-(1-methylethenyl)-benzene | 1096 | – | – | 0.140 | 0.070 | – | 0.009 |
| 29 | Pentylbenzene | 1163 | – | – | 0.003 | 0.003 | 0.023 | – |
| 30 | 2-Methylbutylbenzene | 1163 | – | – | – | – | – | 0.002 |

**[!] The unit.** The raster prints "Relative concentration (µg µl)" with no exponent. Read
literally as µg/µL, pot D's MFT would be 179 mg/L = 1.57 mM. That exceeds the 1.19 mM thiamine
that, per Table 3, supplies its carbon, so the literal reading is impossible. The internal
standard sits at 0.168 mg/L in the sparged sample, and the values span 0.002–1.04. That fits
**IS-equivalent mg/L (= µg/mL), response factor 1**. I use mg/L below and mark every derived
number [Z]. These are semi-quantities: there are no response factors, and the purge efficiency
of a thiol at 50 °C is not that of 2-methyl-3-heptanone.

**Not reported:** bis(2-methyl-3-furyl) disulfide, bis(2-furfuryl) disulfide, H2S, 2-acetylthiazole,
mercaptoketones (3-mercapto-2-pentanone etc.), 4,5-dihydro-2-methyl-3-furanthiol. The purge &
trap at 50 °C is poor for the high-boiling disulfides. **Their absence is not evidence of
absence** [!].

### Table 3. "Proportion of isotopomers formed from xylose/[13C5] xylose (1:1) model reaction."

Footnote a: "Values represent the proportion between labelled and the indicated unlabelled isotopomers
and are based on the abundance of the respective analysed ion signals." Which of pots A–E is not
stated.

| No. | compound | m/z (M+) | m/z (analysed) | unlabeled carbon atoms (%) | labeled carbon atoms (%) | number of labelled carbon atoms |
|---:|---|---:|---|---:|---:|---:|
| 1 | 2-Methyl-3-furanthiol | 114 | 114; 119 | >99 | <1 | 0 |
| 2 | 2-Furfurylthiol | 114 | 114; 119 | 80 | 20 | 5 |
| 3 | Thiophene | 84 | 84; 88 | 88 | 12 | 5 |
| 4 | 2-Methyl-thiophene | 98 | 98; 103 | >99 | <1 | 0 |
| 5 | 2-Pentyl-thiophene | 154 | 154; 159 | >99 | <1 | 0 |
| 6 | Furfural | 96 | 96; 101 | 52 | 48 | 5 |

[!] Thiophene is C4H4S. The analysed ions 84/88 are M/M+4, yet the last column prints "5".
[!] Only M and M+5 (M+4 for thiophene) were read, so any scrambled isotopomers (M+1…M+3 from
fragment recombination) are invisible.

### Table 4. "The amount of DKP producing in the model reaction." (DKP = cyclo(Cys-Gly), 3-mercaptomethyl-2,5-diketopiperazine)

| model reaction | amount of DKP (mg/l) |
|---|---:|
| A | 0.420 |
| B | 0.231 |
| C | 0.379 |
| D | 11.322 |
| E | 17.960 |

### LC-TOF-MS ions (text, p. 283–284) [M] as printed, [!] my mass check

| printed [M+H]+ | assigned to | found in | [Z] theoretical [M+H]+ | verdict |
|---|---|---|---:|---|
| 122.0689 | cysteine | all systems | 122.0270 (C3H8NO2S+) | **[!] off by 0.042 Da (~340 ppm), so not cysteine at TOF accuracy** |
| 130.0475 | pyroglutamic acid (PCA) | A, B, D (the Glu-containing pots) | 130.0499 (C5H8NO3+) | −18 ppm, plausible |
| 140.0792 | 1-(2-thienyl)-1-propanone, unlabelled | — | 141.0369 (C7H9OS+) | **[!] wrong nominal mass** |
| 144.0458 | the same, "4-times-labelled" | — | 145.05 (or 144.09 if 140.0792 + 4 × 1.00335) | **[!] inconsistent either way** |

---

## 4. ISOTOPE LABELLING: WHAT THE 1:1 DESIGN ACTUALLY MEASURES

**The design.** Every pot carries 0.333 mmol [12C]xylose + 0.322 mmol [13C5]xylose, a labelled mole
fraction of **0.492** [Z]. With 98 % atom enrichment, the fully [13C5] fraction of the labelled
xylose is ~0.98^5 = 0.90. A product made **entirely from one intact xylose** should therefore show
M+5 / (M + M+5) ≈ **0.47–0.49**, not 1.00. The paper is a mixed-label CAMOLA design, and its
"labeled carbon atoms (%)" column is **not** the xylose share. **Xylose share ≈ labelled % / 0.48**
[Z]:

| compound | printed unlabelled:labelled | [Z] intact-xylose share | authors' reading | my reading |
|---|---|---:|---|---|
| furfural | 52 : 48 | **~100 %** | "carbon skeleton of xylose chain remained intact" | correct, a built-in positive control that validates the method |
| 2-methyl-3-furanthiol | >99 : <1 | **< ~2 %** | from thiamine and/or cysteine; "supposedly … from thiamine" | **thiamine**. A cysteine-carbon MFT route without sugar "seems unlikely" (Cerny 2007b, quoted here by the authors) |
| 2-furfurylthiol | 80 : 20 | **~42 %** | "partly from xylose"; rest "may be from thiamine and/or cysteine" | **[!] contradicts Cerny 2007b (>99 % xylose in the ternary, in BOTH controls) and Cerny & Briffod 2007 (>99 % at pH 4–6)**. No carbon route from thiamine or cysteine to FFT is known. Most likely an integration or interference artefact at m/z 114 (MFT shares M+ 114; RI 865 vs 908 is resolvable, but the trace signals are small). **Do not use** |
| thiophene | 88 : 12 | **~25 %** (as M+4) | partly from furfural; the rest via cysteine-derived mercaptoacetaldehyde dimerisation (Fig. 2) | the authors' mechanism is a proposal with no labelled-cysteine test |
| 2-methylthiophene | >99 : <1 | < ~2 % | thiamine and/or cysteine | thiamine is the known source (Cerny 2007b Table 5: >99 % unlabelled with no cysteine) |
| 2-pentylthiophene | >99 : <1 | < ~2 % | same | a lipid-type product with no lipid in the pot, so identity is doubtful |

**[!] Which pot is Table 3?** Not stated; the caption says "model reaction", singular. In Table 2
**thiophene appears only in B and C, 2-methylthiophene only in C**, and FFT in A–C. **Only pot C
(free Cys + Gly, no GSH) contains all six Table 3 compounds** [Z]. The isotope evidence in a paper
titled "glutathione–xylose" therefore most likely comes **from the free-cysteine pot, not the GSH
pot**. If it is pooled, it is still not pot A alone.

---

## 5. GLUTATHIONE: HOW IT RELEASES H2S OR CYSTEINE (the question the model asks)

| claim | evidence in the paper | strength |
|---|---|---|
| GSH is cleaved at the γ-glutamyl–cysteinyl amide to **PCA (5-oxoproline) + Cys-Gly** (Fig. 3, "+H2O") | an LC-TOF ion at 130.0475 assigned to PCA, seen only in the Glu-bearing pots A, B, D | **[Q] moderate.** The mass fits PCA, but B and D contain free Glu, which cyclises to PCA by itself, so PCA does not prove GSH cleavage |
| GSH also gives **free cysteine** (abstract; Fig. 2 draws GSH → cysteine directly) | an ion "122.0689" in all systems, assigned to cysteine | **[!] weak.** The mass is not cysteine (§3), and it also appears in D and E, which contain no free Cys |
| Cys-Gly cyclises to **cyclo(Cys-Gly)** (3-mercaptomethyl-DKP) | HPLC-UV vs a standard: A 0.420, D 11.32, E 17.96 mg/L | **[M]** but identified by retention time only. B and C (no peptide) also read 0.23–0.38 mg/L, i.e. co-elution or genuine free-amino-acid condensation [!] |
| The DKP loses its SH to give **3-methyl-2,5-piperazinedione** (cyclo(Ala-Gly)); 3,6-dimethylpiperazine-2,5-dione is also "only found in GSH system" | Table 2 nos. 17–18, pot A only (0.019 and 0.661) | **[Q]**. Desulfurising the DKP would release H2S from **peptide-bound** cysteine, but the paper does not say so or measure it. Cyclo(Ala-Ala) (no. 18) has no Gly and cannot come from cyclo(Cys-Gly) by desulfurisation alone [!] |
| H2S release from cysteine is faster than from GSH, giving 4× the volatiles | **cited**: Zhang, Chien & Ho 1988; not measured here | [Q] cited. Note that Zheng & Ho 1994 (the model's source) reports K(GSH)/K(Cys) ≈ 0.9–1.9 at pH 5, 90–110 °C, i.e. the opposite direction |
| H2S is the vehicle | text: "This may be related to the release of hydrogen sulphide" | [Q] speculation; **no H2S was measured** |

**[Z] Scale of the DKP channel.** Pot A makes 0.420 mg/L DKP = 2.6 µM from 7,419 µM GSH, i.e.
**0.035 mol %**. Pot E converts **0.88 mol %** of its Cys-Gly and pot D 1.06 %. Cyclisation is
therefore a minor fate of GSH at pH 4.8 / 132 °C / 90 min. That is consistent with the model's
choice to carry only GSH → H2S, but it says nothing about the H2S rate.

---

## 6. DERIVED QUANTITIES [Z] (mine, assuming Table 2 = mg/L)

**Molar concentrations, µM:**

| compound (MW) | A (GSH) | B (Glu+Cys+Gly) | C (Cys+Gly) | D (Glu+Cys-Gly) | E (Cys-Gly) |
|---|---:|---:|---:|---:|---:|
| MFT (114.17) | 0.18 | 0.72 | 1.03 | 1.57 | 0.30 |
| FFT (114.17) | 0.43 | 0.58 | 1.49 | 0 | 0 |
| furfural (96.08) | 1.58 | 0.69 | 0.93 | 5.99 | 1.76 |
| thiophene (84.14) | 0 | 3.71 | 4.03 | 0 | 0 |
| 2-methylthiophene (98.17) | 0 | 0 | 1.13 | 0 | 0 |
| DMDS (94.20) | 2.96 | 1.67 | 2.57 | 5.47 | 2.41 |
| DMTS (126.26) | 0 | 0 | 0.41 | 3.11 | 0.41 |

**Yields:**
- MFT per thiamine charged: A 0.016 %, B 0.061 %, C 0.087 %, **D 0.13 %**, E 0.025 % (thiamine
  1.19 mM). If MFT is thiamine-derived in every pot (Table 3 tests only one), the **thiamine →
  MFT yield varies 8.5× with the co-reactant**, and GSH gives the lowest.
- FFT per xylose: 0.002–0.006 %. FFT per sulfur co-reactant: A 0.0058 %, B 0.0060 %, C 0.012 %.
- **A vs B (GSH vs free Cys):** FFT ratio 0.74, against an S-charge ratio of 0.77, so **per mole
  of cysteinyl sulfur, GSH and free Cys give the same FFT (0.96×)**. MFT ratio 0.26 and thiophene
  0 vs 3.7 µM. The two FFT-forming pots agree with Zheng & Ho's K(GSH)/K(Cys) ≈ 0.9–1.0 at pH 5,
  ≥ 110 °C. MFT, which is thiamine-carbon, does not track the S charge at all.
- **Cys-Gly pots D and E make no FFT and no thiophene** while making the most furfural (D). An
  N-terminal-cysteinyl dipeptide that cyclises (§5) apparently does not release H2S at the rate
  free Cys or GSH does. This is a qualitative contrast, with no LOD.
- Total sulfur volatiles ≈ 3–11 µM against 7–13 mM S charged, i.e. **≤ 0.1 % of the sulfur**.

---

## 7. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**The model's relevant parts:** `r_gsh_h2s` (GSH → H2S, Zheng & Ho 1994 Table III, pH-resolved,
first order, residue unassigned); `r_cys_h2s`; the thiamine route `r_thi_hmp → HMP → MFT / MP2P`
(first order, sulfide-independent) and `r_thi_mesh`; B50 `r_thi_h2s` (zero below pH 6.25); the
sugar routes `r_nf_mft`, `r_fur_fft`.

**Cannot anchor (a) any rate constant or Ea.** This is one temperature, one time, no time course,
unbuffered, and gives semi-quantities.

**Cannot anchor (b) a GSH → H2S or GSH → Cys yield.** No H2S was measured. The "cysteine" ion is
mis-assigned. The DKP data bound only a minor side channel (0.035 mol %).

**Usable as directional HOLD-OUT tests only** (none should become FIT rows):
1. **MFT is thiamine-carbon (<~2 % xylose) at dilute loadings** (thiamine 1.2 mM, xylose 26 mM,
   Cys ~10 mM, pH ~4.8, 132 °C, 90 min). This is a third point on the concentration → route-mix
   trend: Cerny 2007b 46 % xylose at 100/300 mM, Cerny & Briffod 15 % at 50/150 mM, here <2 %.
   The model's first-order-thiamine vs second-order-sugar structure predicts the direction. A
   prediction at Wang's pot C conditions should give a thiamine share ≥ ~0.98. The xylose:thiamine
   ratio (22:1 vs 3:1), temperature, time and unbuffered pH all differ, so the test is a sign,
   not a number.
2. **Furfural intact-C5 from xylose** (52:48 ≈ the label mole fraction). This is consistent with
   `TDP → FUR`; nothing new.
3. **GSH vs free Cys, FFT per mole cysteinyl S ≈ 1 (A vs B).** This is consistent with
   `k_gsh_h2s ≈ k_cys_h2s` at pH 5 and ≥ 110 °C, a weak directional check of the Zheng & Ho
   transfer to a Maillard pot. Confounders: B also carries Glu/Gly, the FFT label is suspect
   (§4), and the pH is unbuffered.
4. **The thiamine-route MFT changes 8.5× with the co-reactant** (A 0.18 → D 1.57 µM) at fixed
   thiamine. The model's `THI → HMP → MFT` is independent of cysteine, GSH and H2S, so, sinks
   aside, it would predict roughly constant thiamine-MFT across A–E. If the effect is real, it
   says the thiamine route is co-reactant-sensitive. The sign for A < B < C (more MFT with more
   free Cys) is the **opposite** of what Dwivedi & Arnold 1971's cysteine protection of thiamine
   would give. The pot-to-pot pH is unknown, so this is a flag for the backlog, not a test.

**Not usable:** the FFT isotope split (80:20, contradicted by two Firmenich papers); the thiophene
mechanism (no species in the model and no labelled-cysteine test); the LC-TOF assignments; the
disulfides (not reported; the method is blind to them); DMDS/DMTS (sources unassigned, and in a
pot with thiamine the model's MeSH is thiamine-only anyway).

---

## 8. CAVEATS

1. **No buffer, "pH 4.8"** with no statement of how it was set or whether it held. Xylose
   degradation acidifies an unbuffered pot during 90 min at 132 °C, so the effective pH is
   unknown and probably drifts below 4.8.
2. **Every pot contains thiamine (1.19 mM)**, so there is no thiamine-free control and no pot
   without a sulfur amino acid. Any "GSH effect" is GSH + thiamine.
3. **Semi-quantities**: one IS (2-methyl-3-heptanone), no response factors, purge & trap at 50 °C;
   the unit in Table 2 is misprinted (§3). Ratios between pots for the same compound are the most
   robust use.
4. **24 h at 4 °C before analysis.** Thiols oxidise and add during storage, so the thiol numbers
   are post-storage.
5. **Identifications are questionable in places.** 1,5-pentanediol dinitrate (no nitrogen-oxide
   source), D-limonene and the alkylbenzenes (contaminants), 2-pentylthiophene (no lipid) and
   6-hydroxy-2-hexanone are all claimed against "authentic reference compounds". The claim cannot
   hold for all of them, so treat footnote a as boilerplate.
6. **Table 3's pot is not stated**; by detection pattern it is pot C (§4).
7. The **"(90% enrichment)" Cys-Gly** is probably 90 % purity, which would make D and E ~10 % lower
   in peptide than the [Z] molarities.
8. Reactor volume and headspace are not stated, so no gas-phase partition or H2S loss can be
   computed.
