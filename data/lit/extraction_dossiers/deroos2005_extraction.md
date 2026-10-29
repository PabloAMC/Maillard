# de Roos, Wolswinkel & Sipma 2005 — EXTRACTION
### *"Amadori Compounds of Cysteine and Their Role in the Development of Meat Flavor"*
### Industrial (Givaudan) chapter: cysteine + aldose condense first to the thiazolidine (TTCA-type) in a pre-equilibrium reached in 1 h at 50 °C; on heating (85-95 °C, pH 4-5, 20-40 % water) the thiazolidine drains slowly into the cysteine Amadori compound (ketose-L-cysteine, which exists as the cyclic tetrahydro-1,4-thiazine), which peaks at ~2 h at 90 °C and never exceeds **35 % yield**. The Amadori compound, not the thiazolidine, is the meat-flavour precursor. Numbers are sparse: one yield table (mg thiol per mol cysteine), one relative-peak-area table, one pH-stability curve, one browning chart. **No rate constant, no equilibrium constant, no MFT yield from the purified Amadori compound.** Read 2026-10-28.

**Source on disk:** `data/articles/deroos2005.pdf` (13 pp., ABBYY OCR of the ACS Symposium Series chapter, downloaded via pubs.acs.org). Prose and Table I/II read from the `pdftotext -layout` layer and **checked against 200-dpi page renders** (the text layer garbles Figure 2 and the Table II row structure); **Figures 3 and 4 digitised from 400-dpi renders** (scratchpad `img/drfig3-07.png`, `img/drfig4crop.png`); Figure 1 (strip-chart chromatograms) read qualitatively only.

**Provenance tags:** **[M]** printed number · **[Q]** qualitative statement · **[D]** digitised by me (read error stated) · **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Kris B. de Roos, Kees Wolswinkel, Gerben Sipma |
| Affiliation | Givaudan Nederland BV, P.O. Box 414, 3770 AK Barneveld, The Netherlands |
| Title | "Amadori Compounds of Cysteine and Their Role in the Development of Meat Flavor" |
| Venue | Chapter 9 in *Process and Reaction Flavors*, Weerasinghe, D.K. & Sucan, M.K., Eds.; **ACS Symposium Series 905**; American Chemical Society: Washington, DC, **2005**; pp. **117-129** |
| DOI as printed (page footer) | **10.1021/bk-2005-0905.ch009** ✔ matches the brief; "Publication Date: July 28, 2005" |
| Character | Industrial review-plus-data chapter; no replicates, no SD, no statistics anywhere. Most numbers are single values; several are attributed to ref. 5 (de Roos 1992, ACS Symp. Ser. 490, pp. 203-216; **not on disk**) and ref. 13 (Y. Laarhoven, Givaudan, **unpublished**) |

---

## 1. The one-paragraph answer (for the TTCA-routing decision)

The chapter's mechanism (its Figure 2, "on the analogy of Kallen and Hodge") is
**Cys + aldose ⇌ Schiff base (6) ⇌ thiazolidine-4-carboxylic acid (7a/7b)** and
**Schiff base (6) → Amadori compound (8, open) ⇌ tetrahydro-1,4-thiazine (9a/9b, cyclic hemithioacetal)**.
The thiazolidine is therefore **not** on the direct line to the Amadori compound: it is a
reservoir in equilibrium with the Schiff base, and the Schiff base is the common intermediate
for both. "At pH < 7, the route to Amadori compounds is preferred"; "at pH > 9 … the opposite
occurs" (equilibria shift to the thiazolidinecarboxylic acid, the strongest acid) [Q]. The
evidence offered is a sugar ranking (Figure 3): the sugars that condense most to thiazolidine
also brown fastest at pH ~3-4 (RIB > XYL > GLC > FRU), and the ranking collapses at pH ~9 [D].
The Amadori compound is "the key intermediate to brown pigments and flavor compounds"; its
decomposition to flavour is **first order** and "not influenced by dilution", and its formation
"is not so strongly influenced by dilution" (Table II) [Q]. The thiazolidine → Amadori route is
"general-acid-catalyzed by the other amino acids" (glycine), while the Amadori route of glycine
"is inhibited by cysteine" [Q, inferred from Table II]. Quantitatively, the chapter gives only:
the Amadori compound **peaks at ~2 h at 90 °C, pH 5.0, 20 % water**; **maximum yield 35 %**
(batch, optimum 1:1, 30-40 % water, pH 4-5; 22.5 % in another batch comparison); **~42 % of the
purified xylose-cysteine Amadori decomposes in 2 h in boiling 0.5 M phosphate at pH 5, ~82 % at
pH 6** [D]; and thiol yields (MFT 36 / FFT 3.3 mg per mol cysteine from Cys + xylose) **from the
cysteine-sugar mixture, not from the purified Amadori compound**. The "10 times" superiority of
the Amadori compound over the cysteine-sugar mixture is a **sensory intensity** statement for canned
soups and meats, with no data shown.

---

## 2. Methods as they matter to a model [M]

| item | as printed |
|---|---|
| Reagents | L-cysteine, L-cysteine·HCl monohydrate, D-xylose (Diamalt AG, Munich); other chemicals Merck; ion-exchange resins Bio-Rad |
| HPLC | LiChrosorb NH₂, 300 × 4.6 mm, 10 µm; **refractive-index detector**; 0.02 M ammonium formate in ACN : water 60 : 40, pH 6 (formic acid); 2 mL/min; external standards of the pure compounds |
| General reaction | **L-cysteine·HCl·H₂O 17.6 g (0.1 mol) in 50 mL demineralised water** in a 150 mL round-bottom flask; optionally **Na₂HPO₄ 14.2 g (0.10 mol) or 90 % lactic acid 10 g (0.10 mol)**; pH adjusted to **5.0** with 50 % NaOH; **sugar 0.11 mol** added; **stirred 1 h at 50 °C "to obtain an equilibrium mixture of thiazolidine, cysteine and sugar"**; cooled; pH re-adjusted to the desired value (50 % NaOH or conc. HCl); water content adjusted (vacuum evaporation or water addition); then **2-8 h at 85-95 °C**, reflux condenser, mechanical stirring, thermostated water bath; HPLC samples "at regular time intervals" |
| Isolation | dilute to ~20 % solids; Bio-Rad **AG 50W-X4 (H⁺)**, 200-400 mesh; wash 2 column volumes water; displace with **0.3 M ammonia**; elution order **thiazolidine, then Amadori compound, then cysteine**; Amadori fraction dried, redissolved, ethanol to cloudiness, slow crystallisation (colourless) |
| Isolated compounds (all as crystalline ketose-L-cysteines) | **D-xylulose-L-cysteine** (from D-xylose; decomp. > 170 °C); **D-ribulose-L-cysteine** (from D-ribose; m.p. 132-135 °C dec.); **L-ribulose-L-cysteine** (from L-arabinose; 171-173 °C dec.); **D-fructose-L-cysteine** (from D-glucose; 175-177 °C dec.); **6-deoxy-L-fructose-L-cysteine** (from L-rhamnose; dec. > 183 °C). ¹³C NMR (D₂O, 25.0 MHz) given for major/minor isomer of each |
| Elemental analysis, D-xylulose-L-cysteine | Found C 37.87, H 5.96, N 5.51, S 12.65, O 37.80 %; calc. for C₈H₁₅NO₆S (**253.27**) C 37.94, H 5.97, N 5.53, S 12.66, O 37.90 %; dried in vacuo over P₂O₅ 18 h at 69 °C; no crystal water (Karl Fischer) |
| Structure | "2-hydroxy-2-(1′(S),2′(R),3′-trihydroxypropyl)-tetrahydro-2H-1,4-thiazine-5(S)-carboxylic acid" = cyclic hemithioacetal of D-xylulose-L-cysteine; two isomers in water (O at C2 cis or trans to the C5 carboxyl); "predominantly in the monocyclic form 9a" (NMR, elemental analysis, acetylation) |
| Continuous reactor | 45 × 2.5 cm ID jacketed glass column, **Bio-Rex 9 anion exchange resin (OH⁻)**, 100-200 mesh, **90 °C**; 500 × 0.3 mm Teflon restriction at outlet; feed **50 % thiazolidine solution** (equimolar Cys + Xyl in water, **1 h at 50 °C**) at **0.5 mL/min** |

[Z] Concentration of the general procedure: 0.1 mol cysteine in ~50 mL water plus ~15-16 g solids
(Cys·HCl·H₂O 17.6 g, xylose 16.5 g, buffer 10-14 g) → **~1.2-1.5 mol/L cysteine** before any
evaporation; the heated pots are then concentrated further ("moisture 20 %", "30-40 %").

---

## 3. Results — every number in the chapter

### 3.1 Figure 1 — time course of the thiazolidine and the Amadori compound, 90 °C [Q]

Conditions (caption): **cysteine·HCl + xylose, 90 °C, pH 5.0 (lactate buffer), moisture content
20 %**; chromatograms at **0.5, 1, 2, 3 h**. Peaks: a = solvent + xylose + cysteine (co-eluting);
? = unknown, later identified as D-xylulose-L-cysteine; b = chloride; c = lactate; d =
2-xylulosylthiazolidine-4-carboxylic acid (2 isomers).

Text: "the initially high concentrations of the two 2-xylulosylthiazolidine-4-carboxylic acid
isomers (peak d) decreased slowly, whereas at the same time a new unknown compound was formed …
After a reaction time of **2 hours** the concentration of this compound achieved a **maximum** and
then decreased again."

**[!] Not digitisable.** Strip-chart recorder traces on chart paper, peaks off-scale at 2 h and 3 h,
attenuation unknown, RI detector (response factors differ by compound), and the thiazolidine
peak (d) is labelled only in the 0.5 h trace. From the render [Q]: the "?" peak at 0.5 h is small;
at 1 h it is among the tallest peaks; at 2 h and 3 h it runs off the top of the chart. **No
concentration can be read.** Free cysteine and xylose co-elute with the solvent (peak a), so the
figure carries **no free-cysteine time course**.

### 3.2 Figure 3 — browning of cysteine + sugar, dilute, 90 °C [D]

Caption: "Browning in cysteine-sugar systems in dilute aqueous medium (**74-79 % water**) at
**90 °C**. a: equimolar mixtures of cysteine and sugar, b: equimolar mixtures of cysteine, sugar
and NaOH" (**reprinted from ref. 5**, de Roos 1992). y-axis "absorbance / M sugar", **log scale
0.001-10**; wavelength **not stated** [!]. Bars at 1 h (white), 2 h (checked), 4 h (black). A number
is printed above each 4 h bar (a: 3.9 / 3.7 / 3.4 / 3.3; b: 9.2 / 9.6 / 9.2 / 9.0) — **not
explained in the caption; almost certainly the final pH** [!]. Read error **± 0.05 decade (± 12 %)**.

| panel | sugar | 1 h | 2 h | 4 h | printed number (pH?) |
|---|---|---|---|---|---|
| a (Cys + sugar) | FRU | 0.0047 | 0.0090 | 0.013 | 3.9 |
| a | GLC | 0.030 | 0.065 | 0.13 | 3.7 |
| a | XYL | 0.054 | 0.11 | 0.30 | 3.4 |
| a | RIB | 0.11 | 0.32 | 0.96 | 3.3 |
| b (Cys + sugar + NaOH) | FRU | 0.155 | 0.60 | 1.6 | 9.2 |
| b | GLC | 0.155 | 0.53 | 1.3 | 9.6 |
| b | XYL | 0.115 | 0.23 | 0.56 | 9.2 |
| b | RIB | 0.10 | 0.16 | 0.36 | 9.0 |

[Z] At acidic pH the browning order RIB : XYL : GLC : FRU at 4 h is ≈ 74 : 23 : 10 : 1; with NaOH
the order inverts (FRU ≈ GLC > XYL > RIB) and the pentoses brown **less** than at acidic pH (XYL
0.56 vs 0.30, i.e. only ×1.9; RIB 0.36 vs 0.96, i.e. ×0.4), which is the chapter's "highest
inhibition for the most reactive sugars at pH > 9".

### 3.3 Table I (verbatim) — Likens-Nickerson extraction of mixtures of cysteine with carbonyl compounds under optimum conditions for Amadori compound formation [M]

**Quantity produced (mg/mole of cysteine)**

| compound | Ribose | Xylose | Glucose | Norfuraneol |
|---|---|---|---|---|
| 2-Methylfuran-3-thiol | 43.5 | 36 | < 0.05 | 0.9 |
| 2-Furfurylthiol | 2.8 | 3.3 | < 0.05 | 10.4 |
| 3-Mercapto-2-pentanone | 8.0 | 5.6 | 0.5 | < 0.05 |
| 2-Mercapto-3-pentanone | 16.4 | 10.5 | < 0.05 | 24.3 |

Footnote a (verbatim): "Cysteine : carbonyl compound : lactate = 1 : 1 : 1.67; water content: 37 %
(w/w); pH 5.0; boiling for 1.5 h." Source attribution in text: **ref. 13 (Laarhoven, Givaudan,
unpublished)**.

[Z] Molar yields (MW: MFT and FFT 114.17; mercaptopentanones 118.20):

| compound | Ribose | Xylose | Norfuraneol |
|---|---|---|---|
| MFT, mol/mol Cys | **3.8 × 10⁻⁴ (0.038 %)** | **3.2 × 10⁻⁴ (0.032 %)** | 7.9 × 10⁻⁶ |
| FFT, mol/mol Cys | 2.5 × 10⁻⁵ | **2.9 × 10⁻⁵ (0.0029 %)** | **9.1 × 10⁻⁵** |
| 3-mercapto-2-pentanone | 6.8 × 10⁻⁵ | 4.7 × 10⁻⁵ | < 4 × 10⁻⁷ |
| 2-mercapto-3-pentanone | 1.39 × 10⁻⁴ | 8.9 × 10⁻⁵ | 2.06 × 10⁻⁴ |
| MFT : FFT | 15.5 | **10.9** | 0.09 |

[Z] Concentration: per mol cysteine, 121 g Cys + 150 g pentose + 150 g lactic acid ≈ 420 g solids
at 37 % water → **~250 g water per mol Cys, i.e. ~4 mol/kg water** (more if NaOH for pH 5 is
counted as solids). These are **extremely concentrated** pots.

Text attached to Table I [Q]: norfuraneol "does not produce higher yields of 2-methylfuran-3-thiol
and 3-mercapto-2-pentanone during boiling than ribose", so norfuraneol "is not playing a
significant role" in their formation during boiling (with Cerny & Davidek 2003, ref. 19); high
MFT from norfuraneol + Cys/H₂S needs "140 °C or higher in a closed system". The glucose column
shows ketohexose-L-cysteines do not make the boiled-meat thiols.

**[!] What Table I is and is not.** It is **cysteine + free sugar**, run under conditions the
authors judge optimal for Amadori formation (1:1, pH 5, 37 % water). It is **not** a yield from
the purified ketose-L-cysteine, and the chapter prints **no** MFT/FFT yield from the purified
Amadori compound anywhere. Whether "boiling for 1.5 h" is the reaction time or the duration of
the Likens-Nickerson simultaneous distillation-extraction is ambiguous; in an SDE the pot is
diluted and boiled, so the two are probably the same operation.

### 3.4 Thiamine comparison [M, text p. 126]

Cysteine + thiamine "when carried out under the conditions of Table I, generates **about 25 times
more** 2-methylfuran-3-thiol than the reaction of cysteine with ribose" [Z: ≈ 1.1 g MFT per mol
cysteine ≈ 0.95 mol %]. Yet Grosch et al. (ref. 24) found that adding cysteine + thiamine to meat
does not enhance MFT during boiling; the authors attribute this to the reaction being **second
order** and too slow at boiled-meat dilution, versus the **first-order** decomposition of
ketose-cysteines, "not influenced by dilution".

### 3.5 Table II (re-typed from the page render) — Formation of Amadori compounds in cysteine-amino-acid-glucose model systems at 90 °C [M]

| CysSH (moles) | Glc (moles) | amino acid (moles) | water (w %) | reaction time (h) | D-Fru-L-Cys | D-Fru-AA ᵃ |
|---|---|---|---|---|---|---|
| 1 | 1.5 | Glycine (1) | 15 | 1 | 26 | 29 |
| | | | | 2 | 33 | 29 |
| 1 | 1.5 | L-Proline (1) | 15 | 1 | 40 | 6 |
| | | | | 2.5 | 46 | 16 |
| 1 | 1 | Glycine (1) | 79 | 2 | 1.5 | 0.3 |
| | | | | 4 | 4 | 1 |
| 1 | 1 | L-Proline (1) | 79 | 2 | 3.8 | 0.1 |
| | | | | 4 | 7.5 | 0.25 |

Columns 6-7 are **"Relative HPLC peak areas"** (RI detector). ᵃ "D-Fru-AA = Amadori compound of
glucose with amino acid other than cysteine." Text: the 15 %-water pots had "no pH adjustment,
initial pH 5-6"; the 79 %-water pots are called "dilute neutral media" [!: pH not printed].

[Z] Ratios (peak-area units, response factors not given, so mol ratios are uncertain by the
RI response ratio, likely within ×1.5):
- Fru-Cys / Fru-Gly: **0.9-1.1** at 15 % water; **4-5** at 79 % water.
- Fru-Cys / Fru-Pro: 2.9-6.7 at 15 % water; **30-38** at 79 % water.
- Dilution penalty (15 % → 79 % water, glycine arm, 2 h): Fru-Cys **33 → 1.5 (÷22)**; Fru-Gly
  **29 → 0.3 (÷97)**. Proline arm (2-2.5 h): Fru-Cys 46 → 3.8 (÷12); Fru-Pro 16 → 0.1 (÷160).
- The chapter's claims: Cys-Amadori formation is "not so strongly influenced by dilution" and
  "cysteine becomes clearly the most reactive amino acid in dilute neutral media"; glycine is
  "about 100 times more reactive than cysteine" in the Maillard reaction (ref. 8, Ames 1986),
  so equal Fru-Cys and Fru-Gly rates at 15 % water are read as **glycine general-acid catalysis of
  the cysteine route and cysteine inhibition of the glycine route**.

**[!] Table II is glucose, not a pentose.** It gives the hexose Amadori (D-fructose-L-cysteine)
only; no pentose analogue is tabulated.

### 3.6 Figure 4 — pH stability of purified D-xylulose-L-cysteine, boiling 0.5 M phosphate [D]

Caption: "Effect of pH on the decomposition of D-xylulose-L-cysteine in a boiling 0.5 M aqueous
phosphate buffer." y = "% left after 2 h" (0-100), x = pH 2-7. A single dashed curve, no markers,
spanning pH 3-6. **Amadori concentration not stated; products not stated.** Read error **± 3 %
(absolute), ± 0.05 pH**.

| pH | 3.0 | 3.5 | 4.0 | 4.5 | 5.0 | 5.25 | 5.5 | 6.0 |
|---|---|---|---|---|---|---|---|---|
| % left after 2 h | 51 | 53 | 54.5 | 56 | **58** (max) | 54.5 | 43 | **18** |
| [Z] k₁ (min⁻¹), -ln(f)/120 | 5.6 × 10⁻³ | 5.3 × 10⁻³ | 5.1 × 10⁻³ | 4.8 × 10⁻³ | **4.5 × 10⁻³** | 5.1 × 10⁻³ | 7.0 × 10⁻³ | **1.43 × 10⁻²** |

[Z] Half-life at pH 5 ≈ 2.5 h; at pH 6 ≈ 0.8 h, at **~100 °C** (boiling, 0.5 M phosphate; the
buffer is dilute in reactant, so this is the "dilute" first-order decay). Stability maximum at
**pH ~5**; d log k / d pH ≈ **+0.5 per unit between pH 5 and 6** (+0.6 between 5.5 and 6). Text: the
"relatively low reactivity in products under kitchen conditions" frustrates in-situ use.

### 3.7 Yields and process numbers [M]

| item | value | conditions |
|---|---|---|
| Thiazolidine pre-equilibrium | "stirred for **1 h at 50 °C** to obtain an **equilibrium mixture** of thiazolidine, cysteine and sugar" | Cys·HCl 0.1 mol + sugar 0.11 mol in 50 mL water, pH 5.0 |
| Max D-xylulose-L-cysteine yield, batch | **35 %** ("only 35 %, which is due to its instability under the conditions of its formation") | optimum: **1:1 molar, water 30-40 %, buffer pH 4-5** |
| Batch yield quoted beside the continuous run | **22.5 %** | "the corresponding batch process" (conditions not further specified) |
| Continuous reactor conversion of thiazolidine → Amadori | **12 % after 2 h → 35 % after 6 h** (time on stream) | 90 °C, Bio-Rex 9 (OH⁻), 50 % thiazolidine feed, 0.5 mL/min |
| Continuous vs batch | "allowed **50 % higher conversions**" (35 vs 22.5 %) but "the maximum yields were again not higher than 35 %" | — |
| Off-flavour cut-off | stop the Cys + Xyl reaction "before **75 %** of the maximum Amadori compound concentration has been achieved" | crude reaction product as precursor |
| Flavour intensity, Amadori vs Cys + sugar | "about **10 times** higher" | canned soups and meats; sensory, **no data shown** |
| Two pathways | Cys-Amadori route and Cys + other-amino-acid Amadori route "probably of about equal importance" for boiled meat and process flavourings | inferred from Table II |

### 3.8 Qualitative claims worth carrying [Q]

1. Sugars that convert most to thiazolidines also brown and flavour most (pH < 7); the Schiff base
   is the shared intermediate (Fig. 2).
2. Ketopentose-L-cysteines make MFT, FFT and 3-mercapto-2-pentanone; ketohexose-L-cysteines do not
   (Table I glucose column).
3. Thiol instability: thiols react with melanoidins and with each other; second-order in thiol, so
   high-yield pots go "rubbery" (refs 6, 25: van Seeventer 2001, Hofmann & Schieberle 2002).
4. Unsaturated aldehydes from fat bind cysteine as thiazolidines and bind thiols as
   (hemi)thioacetals and sulfides; proposed to explain why chicken (more unsaturated fat) yields
   fewer thiols than beef at similar free cysteine (refs 20-23).
5. In meat, free cysteine ≪ the sum of other free amino acids, so much meat flavour may form via
   other amino acids' Amadori compounds + cysteine, especially in concentrated crusts.

---

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| Cys + pentose ⇌ thiazolidine | time to equilibrium | **≤ 1 h at 50 °C** | — | ~1.5-2 mol/L, pH 5 | §2 general procedure | qualitative_bound (no K printed) |
| thiazolidine → Cys-Amadori (xylose) | time of Amadori maximum | **~2 h** | h | 90 °C, pH 5.0 lactate, 20 % water | Fig. 1 text | within_study_timing |
| same | max yield | **35 %** (22.5 % in one batch) | mol % of Cys | 1:1, 30-40 % water, pH 4-5, 85-95 °C | §3.7 | level_only |
| same, continuous | conversion | 12 % (2 h) → 35 % (6 h) | % | 90 °C, anion exchanger | §3.7 | process-specific; not transferable |
| D-xylulose-L-cysteine decomposition | % left after 2 h | 51-58 % (pH 3-5), 43 % (pH 5.5), 18 % (pH 6) | % | boiling 0.5 M phosphate | Fig. 4 [D] | digitised (± 3 %) |
| same | apparent first-order k | **4.5 × 10⁻³ (pH 5) - 1.4 × 10⁻² (pH 6)** | min⁻¹ | ~100 °C, dilute | §3.6 [Z] | derived_from_digitised (assumes first order) |
| Cys + Xyl → MFT | yield | **3.2 × 10⁻⁴** | mol/mol Cys | 37 % water, pH 5, lactate, boiling 1.5 h | Table I [Z] | level_only (single run, unpublished source) |
| Cys + Xyl → FFT | yield | **2.9 × 10⁻⁵** | mol/mol Cys | same | Table I [Z] | level_only |
| Cys + Rib → MFT / FFT | yield | 3.8 × 10⁻⁴ / 2.5 × 10⁻⁵ | mol/mol Cys | same | Table I [Z] | level_only |
| Cys + norfuraneol → MFT / FFT | yield | 7.9 × 10⁻⁶ / 9.1 × 10⁻⁵ | mol/mol Cys | same | Table I [Z] | level_only |
| Cys + thiamine → MFT | ratio to Cys + Rib | **~25×** | — | Table I conditions | §3.4 | within_study_ratio |
| Fru-Cys vs Fru-Gly formation | relative peak area | 0.9-1.1 (15 % water), 4-5 (79 % water) | — | 90 °C, glucose | Table II [Z] | within_study_ratio (RI response unknown) |

**[Z] A rough two-step reading of the formation curve (derived_assumption, not a prior).**
Treat the cysteine-bearing pool (Cys + thiazolidine) as A, the Amadori compound as B, and
everything downstream as C, with first-order A → B → C. A maximum at t = 120 min with peak yield
0.35 requires k_B/k_A ≈ 1.1, i.e. **k_A ≈ k_B ≈ 8 × 10⁻³ min⁻¹** (90 °C, 20-40 % water, pH 4-5); with
peak yield 0.225, **k_A ≈ 5 × 10⁻³, k_B ≈ 1.2 × 10⁻² min⁻¹**. Two conclusions survive the
assumption: (i) under formation conditions the Amadori compound decays **at least as fast as it
forms**, which is why the yield never exceeds ~35 %; (ii) that decay (~1 × 10⁻² min⁻¹ at 90 °C,
concentrated) is **faster** than the Fig. 4 value at pH 5 in dilute boiling buffer (4.5 × 10⁻³ min⁻¹
at ~100 °C), i.e. the concentrated medium (or the free sugar/cysteine in it) accelerates Amadori
loss by ≳ 3-5× after correcting for the 10 K. Note also that the Fig. 1 maximum at 2 h would, with
the Fig. 4 decay rate scaled to 90 °C (~2 × 10⁻³ min⁻¹), imply a peak yield near 80 %, which the
chapter says is never reached; the simple dilute-decay number cannot be transplanted into the
concentrated pot.

---

## 5. What this can and cannot anchor in the model

**Model context.** `src/kinetic_core/sulfur.py` (B58 candidate, uncommitted) adds Cys + PENT ⇌ TTCA
with K from Cao 2017 (≈ 60 L/mol, 90 °C, pH 4.2-4.7) and considers routing the pentose's
cysteine-catalysed loss **through** TTCA → Cys-Amadori → deoxyosones. `parameters_sulfur.py`
already carries Zhang 2025b's Ea for Cys → Cys-Amadori (100.9 kJ/mol) and Cys-Amadori → α-DC
(85.7 kJ/mol).

**Can:**
1. **Topology.** The chapter is the primary source behind the Cao/Hou/Zhai group's statement that
   TTCA → Cys-Amadori is the meat-flavour gate (Cao 2017 cites it for the standards and the
   mechanism). Its Figure 2 makes the **Schiff base**, not TTCA, the branch point: TTCA is an
   off-path reservoir. In a mass-action model with a fast pre-equilibrium, "TTCA → CysARP at
   k_iso[TTCA]" and "Schiff → CysARP at k[Schiff] with Schiff ∝ [TTCA]" are kinetically
   indistinguishable, so routing the loss as k_iso·[TTCA] (= k_iso·K·[Cys][PENT]) is **consistent**
   with de Roos; what is *not* consistent is a TTCA step that degrades to deoxyosones without
   passing through the Amadori compound.
2. **Pre-equilibrium is fast.** "1 h at 50 °C to obtain an equilibrium mixture" at ~1.5 mol/L
   supports Cao's inference that TTCA is equilibrated within 30 min at 90 °C.
3. **Cys-Amadori is short-lived at 90-100 °C.** Fig. 4 gives a first-order decay of the purified
   pentose Cys-Amadori of **4.5 × 10⁻³ min⁻¹ at pH 5 rising to 1.4 × 10⁻² min⁻¹ at pH 6** (~100 °C,
   dilute), with a stability maximum at pH ~5. This is a usable **hold-out band** on any
   `CysARP → DPO/TDP` step at pH 3-6 and a pH slope (+0.5 log units per pH unit above 5). Pair
   with Zhang 2025b's Ea 85.7 kJ/mol to transfer to other temperatures (assumption).
4. **Yield ceiling for the intermediate.** Cys-Amadori peaks at ≤ 35 % of charged cysteine (90 °C,
   concentrated) — a falsifiable ceiling on the modelled CysARP pool, and evidence that its decay
   rate under formation conditions is comparable to its formation rate.
5. **MFT/FFT yields from Cys + pentose at pH 5 (Table I):** MFT 3.2-3.8 × 10⁻⁴ and FFT 2.5-2.9 ×
   10⁻⁵ mol/mol Cys, **MFT : FFT ≈ 11-16**. Usable as an order-of-magnitude hold-out for very
   concentrated (≈ 4 mol/kg) open-boiling pots only; note the SDE steam-distils the thiols out as
   they form, which protects them from the second-order losses the chapter itself describes.
6. **Norfuraneol is a poor MFT source at 100 °C** (MFT 0.9 vs 36-43.5 mg/mol from the pentoses) —
   supports keeping the pentose → MFT flux off a norfuraneol-only route at boiling temperatures.

**Cannot:**
- **No equilibrium constant for Cys + sugar ⇌ thiazolidine** and no thiazolidine concentration
  at any time; Cao 2017 remains the only K.
- **No rate constant for TTCA → Cys-Amadori.** Figure 1 is an un-digitisable strip chart; the
  continuous-reactor conversion is resin-catalysed (OH⁻ anion exchanger) and not transferable.
- **No MFT/FFT yield from the purified Amadori compound**, so the chapter cannot decide whether
  per-mole flavour yields via the Amadori route exceed those via free Cys + sugar. The "10 ×"
  is sensory, in canned products, with no data.
- No temperature dependence (90 °C, ~100 °C only), no H₂S, no deoxyosone measurement, no
  replicates.
- Table II is glucose, relative RI peak areas only.

---

## 6. Flags

1. **Industrial chapter, mostly secondary/unpublished data.** Table I is ref. 13 (unpublished);
   Figure 3 is reprinted from ref. 5 (1992). No replicates, no error bars, no SDs.
2. **Figure 3 axis and annotations.** Absorbance wavelength not stated; the numbers above the 4 h
   bars are unexplained (read as final pH). Absorbance per M sugar on a log scale.
3. **"Optimum conditions" are stated three ways:** 1:1, 37 % water, pH 5.0 (Table I); 1:1, 30-40 %
   water, pH 4-5 (p. 127); 90 °C, 20 % water, pH 5.0 lactate (Fig. 1). Yields (35 %, 22.5 %) are not
   tied to one of them explicitly.
4. **Figure 4 lacks the starting concentration, the number of points, and the products.** The
   line is drawn as a dashed curve with no markers; treat as a smoothed sketch (± 3 % at best).
5. **Ref. 6 is miscited in the chapter** ("Van Seeventer … J. Agric. Food Chem. 2002, 49,
   4292-4295"); volume 49 is 2001 (`vanseeventer2001_extraction.md`).
6. **Table I "boiling for 1.5 h"** is probably the Likens-Nickerson SDE itself, i.e. volatiles are
   continuously stripped; yields are therefore *formation* yields with losses suppressed, not
   standing concentrations in a closed pot.
7. **Registry gaps** (`data/keys/compounds.yml`): D-xylulose-L-cysteine (Cys-Amadori), the
   thiazolidine (TTCA) and 2-mercapto-3-pentanone should be checked before any row is written;
   paper not checked against `data/keys/papers.yml`.

---

## 7. Relations to existing dossiers

| dossier | relation |
|---|---|
| `cao2017_extraction.md` | Uses this chapter's isolation procedure (AG 50W-X4, 0.3 M NH₃, TTCA from 50 °C / 1 h) for its TTCA and Cys-Amadori standards and adopts its mechanism ("acid-catalysed by other amino acids"; Cys-Amadori, not TTCA, is the flavour precursor). Cao supplies what de Roos lacks: K ≈ 60 L/mol and time-resolved TTCA / Cys-Amadori. Cao's blank Cys-Amadori/TTCA keeps rising to 5-7 h at 90 °C and 0.4 M, whereas de Roos's concentrated pot peaks at 2 h — consistent with de Roos's claim that concentration and other amino acids accelerate the route. Cao's glycine acceleration of TTCA drain is the same effect as Table II's "general acid catalysis by glycine". |
| `zhai2019_extraction.md` | Zhai's ≤ 6 % Cys-ARP/TTCA at 90 °C, pH 7.4 vs de Roos's ≤ 35 % at pH 4-5 fits the de Roos pH picture (Amadori route preferred at pH < 7; Fig. 4 shows the Amadori compound itself is 3× less stable at pH 6 than at pH 5). |
| `zhang2025b_extraction.md` | Supplies Ea for Cys → Cys-Amadori (100.9 kJ/mol) and Cys-Amadori → α-DC (85.7 kJ/mol) at pH 6.5, 120 °C; de Roos's Fig. 4 decay band at ~100 °C, pH 5-6 is an independent cross-check on Zhang's absolute rates. Zhang's lysine effect is the same "other amino acids catalyse the cysteine route" family as de Roos's glycine. |
| `hou2017_extraction.md` | Hou ranks reactivity Gly-Amadori > Cys-Amadori > TTCA (fed, 140 °C); de Roos's "Amadori superior to thiazolidine or equivalent cysteine-sugar mixture" is the same ordering (sensory). |
| `wang2026_extraction.md` | Wang's Cys-Amadori / Glu-Amadori thiol ladder is the quantitative version of de Roos's second pathway (cysteine + other amino acids' Amadori compounds). |
| `zhou2022_extraction.md`, `zhou2023c_extraction.md` | Ala-xylose Amadori (the engine's `ARP`) + cysteine: these are measured instances of de Roos's "cysteine reacts with (decomposition products of) Amadori compounds of other amino acids" pathway; de Roos says it is of comparable importance to the Cys-Amadori route. |
| `vanseeventer2001_extraction.md` | ref. 6, thiol instability. |
