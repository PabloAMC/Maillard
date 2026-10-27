# Hidalgo, León & Zamora 2016 — EXTRACTION
### *"Amino acid decarboxylations produced by lipid-derived reactive carbonyls in amino acid mixtures"*
### Phenylalanine + 2,4-decadienal + one of 19 second amino acids, 140 °C, pH 3. Only phenylalanine's products (2-phenylethylamine, phenylacetaldehyde) are measured. Leucine, isoleucine, valine and methionine were added, but their own Strecker aldehydes were NOT measured. Read 2026-10-27.

**Full read of `data/articles/hidalgo2016.pdf`** (6 pages, `pdftotext -layout`). Tables 1–2 came
through cell by cell. Figures 1–4 were rendered to PNG and inspected. The Supplementary data
(Figs. S-1 to S-5: the time courses behind every Ea) are NOT on disk. Companions:
`zamora2015_extraction.md` (the same protocol with phenylalanine alone: the 12 × 12 amine/aldehyde
table), `zamora2013_extraction.md` (Ea ladder; 2,4-decadienal under nitrogen 78.0 kJ/mol),
`zamora2007_extraction.md`, `zamora2008_extraction.md`.

**Provenance:** **[M]** measured, printed · **[F]** figure only, not digitised · **[Q]** qualitative ·
**[Z]** derived here · **[!]** flag.

---

## 0. IDENTITY (matches)

| field | value as printed |
|---|---|
| Authors | Francisco J. Hidalgo, M. Mercedes León, Rosario Zamora* |
| Venue | *Food Chemistry* 209 (2016) 256–261 |
| DOI | 10.1016/j.foodchem.2016.04.032 (PMID 27173560) |
| Dates | received 10 Dec 2015; revised 8 Apr 2016; accepted 12 Apr 2016; online 13 Apr 2016 |
| Affiliation | Instituto de la Grasa, CSIC, Seville; corresponding rzamora@ig.csic.es |
| Funding | FEDER + Plan Nacional I+D, AGL2012-35627, AGL2015-68186-R |
| PDF character | 6 pages; Tables 1–2; Figures 1–4 (Fig. 4 = Arrhenius plots); Supplementary Figs S-1 to S-5 not on disk |
| Registry | DOI not in `data/keys/papers.yml` at time of reading |

---

## 1. THE ONE-PARAGRAPH ANSWER

Equimolar phenylalanine + lipid carbonyl + a second amino acid (10 µmol each on 600 mg sand,
110 µL pH 3 citrate/water, closed tube, 140 °C, 1 h) was used to see how a second amino acid shifts
phenylalanine's amine/aldehyde split. **Only 2-phenylethylamine (PEA) and phenylacetaldehyde (PAC)
were quantified** (GC-MS ions m/z 91 and 122). No Strecker aldehyde, amine or α-keto acid of the
second amino acid was measured or mentioned. Most second amino acids **raise PEA** (arginine most:
+309 % under N2) and **lower PAC** (tryptophan, cysteine and proline most: −75 to −90 %). With
2,4-decadienal in air, equimolar leucine, isoleucine, valine or methionine cut PAC from 129.7 to
42.7–57.0 µmol/mmol Phe (−56 to −67 %) [Z]. The authors attribute this to the second amino acid's
α-amino group and its side-chain amino/guanidino groups reacting with the decadienal. **Ea is
unchanged by the second amino acid**: PEA 47–54 kJ/mol, PAC 71–81 kJ/mol (DD, N2, pH 3, ~100–180 °C;
measured for none, Arg, Lys and Cys only). **For the model's question:** this paper gives **no
number for 3-methylbutanal, 2-methylbutanal, 2-methylpropanal or methional**. What it does give is
measured evidence that Leu, Ile, Val and Met compete with phenylalanine for 2,4-decadienal at
comparable strength, which is indirect support for treating them as kinetically similar substrates.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reactants | phenylalanine + lipid derivative + second amino acid, **10 µmol each** (equimolar) |
| Matrix | homogenised with **600 mg sand** (50–70 mesh), **30 µL 0.3 M sodium citrate pH 3 + 80 µL water**. Nominal **90.9 mmol/L** each if all is dissolved in 110 µL [Z]; lipids are poorly water-soluble, so the effective phase is undefined |
| pH | **3 only** (citrate). No pH 6 run |
| Atmosphere | "controlled": nitrogen or air (Table 1); Table 2 and all kinetics **nitrogen** ("to avoid oxidations") |
| Heat | "most experiments" **1 h at 140 °C**; kinetics at several temperatures (Fig. 4 abscissa 0.0022–0.0027 K⁻¹; points at ≈179, 162, 140, 121, 100 °C [F, read off ticks, approximate]) |
| Cooling | 5 min room temperature, 15 min at −30 °C (**[!]** the minus sign is an unmapped glyph in the text layer, which reads "30 °C"; −30 °C as in Zamora 2015) |
| Work-up | + 20 µL internal standard (24.09 mg 2-ethylpyridine / 50 mL methanol) + 1 mL methanol–water 80:20; centrifuge 10 min at 2000g; 700 µL supernatant + 1 mg NaBH4, 30 min (PAC → 2-phenylethanol); + 500 µL acetone; GC-MS "as described previously" (Zamora 2015) |
| Ions | IS m/z 106; **PEA m/z 91; PAC (as 2-phenyl-1-ethanol) m/z 122**. No other analyte monitored |
| Calibration | "described previously" (Zamora 2015: ten-level curves in sand/pH 6 phosphate) |
| Second amino acids (Table 1) | **19**: Gly, Ala, **Val, Leu, Ile**, Pro, Trp, **Met**, Ser, Thr, Cys, Tyr, Asn, Gln, Asp, Glu, Lys, Arg, His. Plus methylguanidine and butylamine (Figs 2–3) |
| Lipids (Table 2) | 13-HPODE (acid), its methyl ester, methyl 13-hydroperoxy-linolenate; 2-pentenal, 2-octenal; 2,4-hexadienal, 2,4-heptadienal, 2,4-decadienal; 4,5-epoxy-2-heptenal, 4,5-epoxy-2-decenal; 4-oxo-2-hexenal, 4-oxo-2-nonenal; 4-hydroxy-2-nonenal |
| Kinetics | [PEA] = [PEA]0 + kt and [PAC] = [PAC]0 + kt (zero-order, linear in time); k → Arrhenius |
| Statistics | mean ± SD of ≥ 3 independent experiments; Tukey, p < 0.05 (Origin 7.0) |

Units throughout: **µmol per mmol phenylalanine** (= 0.1 % molar yield on Phe per unit).

---

## 3. RESULTS

### 3a. Table 1 — Phe + 2,4-decadienal + second amino acid (140 °C, 1 h, pH 3) [M]

Asterisks are significance against the "None" control (* p < 0.05, ** p < 0.01, *** p < 0.001).
The four columns on the right are [Z] (change versus the control, from the printed means).

| second amino acid | PEA N2 | PEA air | PAC N2 | PAC air | ΔPEA N2 | ΔPEA air | ΔPAC N2 | ΔPAC air |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| None | 21.9 ± 3.6 | 10.6 ± 2.4 | 21.1 ± 3.4 | 129.7 ± 4.6 | — | — | — | — |
| Glycine | 50.6 ± 2.0*** | 26.1 ± 2.1*** | 10.7 ± 0.5*** | 119.3 ± 3.0** | +131 % | +146 % | −49 % | −8 % |
| Alanine | 30.5 ± 0.7** | 23.8 ± 6.9*** | 14.0 ± 0.3*** | 103.9 ± 9.9** | +39 % | +125 % | −34 % | −20 % |
| **Valine** | 40.4 ± 5.3*** | 16.7 ± 0.6** | 11.6 ± 3.2*** | 57.0 ± 7.8*** | +84 % | +58 % | −45 % | −56 % |
| **Leucine** | 35.4 ± 8.1*** | 19.5 ± 1.5*** | 7.4 ± 3.9*** | 45.2 ± 13.0*** | +62 % | +84 % | −65 % | −65 % |
| **Isoleucine** | 48.8 ± 1.4*** | 21.2 ± 1.6*** | 6.3 ± 1.8*** | 42.7 ± 9.9*** | +123 % | +100 % | −70 % | −67 % |
| Proline | 45.1 ± 10.4*** | 18.1 ± 2.0*** | 5.2 ± 2.5*** | 32.3 ± 1.9*** | +106 % | +71 % | −75 % | −75 % |
| Tryptophan | 33.2 ± 0.8*** | 17.0 ± 2.6** | 2.4 ± 0.7*** | 13.3 ± 1.3*** | +52 % | +60 % | −89 % | −90 % |
| **Methionine** | 24.6 ± 1.2 | 12.9 ± 0.2 | 9.4 ± 1.7*** | 46.9 ± 7.1*** | +12 % (n.s.) | +22 % (n.s.) | −55 % | −64 % |
| Serine | 43.5 ± 1.7*** | 16.6 ± 5.0* | 24.7 ± 3.1 | 123.4 ± 5.1 | +99 % | +57 % | +17 % (n.s.) | −5 % (n.s.) |
| Threonine | 24.5 ± 1.7 | 15.1 ± 1.0* | 15.4 ± 0.1* | 79.4 ± 10.7*** | +12 % (n.s.) | +42 % | −27 % | −39 % |
| Cysteine | 63.9 ± 5.8*** | 24.3 ± 11.8** | 3.4 ± 0.7*** | 25.5 ± 10.9*** | +192 % | +129 % | −84 % | −80 % |
| Tyrosine | 31.9 ± 6.4** | 13.7 ± 0.4 | 15.8 ± 1.1** | 65.9 ± 6.6*** | +46 % | +29 % (n.s.) | −25 % | −49 % |
| Asparagine | 35.4 ± 5.3*** | 11.3 ± 2.4 | 14.5 ± 0.6** | 76.2 ± 14.5*** | +62 % | +7 % (n.s.) | −31 % | −41 % |
| Glutamine | 31.4 ± 4.4** | 9.2 ± 2.1 | 14.3 ± 1.2** | 73.9 ± 9.3*** | +43 % | −13 % (n.s.) | −32 % | −43 % |
| Aspartic acid | 38.3 ± 0.8*** | 16.9 ± 6.8* | 18.1 ± 4.8 | 77.7 ± 14.5*** | +75 % | +59 % | −14 % (n.s.) | −40 % |
| Glutamic acid | 39.6 ± 4.3*** | 16.6 ± 3.6* | 23.2 ± 0.6 | 111.3 ± 21.8** | +81 % | +57 % | +10 % (n.s.) | −14 % |
| Lysine | 62.3 ± 11.6*** | 26.6 ± 7.7*** | 15.4 ± 1.1*** | 61.4 ± 9.8*** | +184 % | +151 % | −27 % | −53 % |
| Arginine | 89.6 ± 7.2*** | 24.1 ± 3.3*** | 18.6 ± 3.3 | 99.5 ± 23.2** | +309 % | +127 % | −12 % (n.s.) | −23 % |
| Histidine | 44.9 ± 1.0*** | 15.0 ± 1.5* | 14.8 ± 2.1** | 35.4 ± 7.1*** | +105 % | +42 % | −30 % | −73 % |

**Text checks [Z]:** Arg +309 % PEA under N2 (89.6/21.9 = 4.09) ✓; Cys and Lys "about 200 %" (+192,
+184) ✓; in air Lys, Gly, Cys, Arg, Ala "124–151 %" (+125 to +151) ✓; Trp, Cys, Pro PAC decrease
"75–88 % under nitrogen and 75–90 % under air" (75–89 %, 75–90 %) ✓; Arg raises PEA + PAC by "152 %"
((89.6 + 18.6)/(21.9 + 21.1) = 2.52) ✓. **Control vs Zamora 2015 [Z]:** the "None" row matches
Zamora 2015 Table 1, DD pH 3 (PEA 30.99 / 11.23; PAC 21.07 / 126.82 under N2 / air) within the SDs:
the same protocol, reproduced.

**Total phenylalanine degraded (PEA + PAC) with the meat-aroma amino acids [Z]:**

| second AA | N2 | air |
|---|---:|---:|
| None | 43.0 | 140.3 |
| Valine | 52.0 | 73.7 |
| Leucine | 42.8 | 64.7 |
| Isoleucine | 55.1 | 63.9 |
| Methionine | 34.0 | 59.8 |

Under N2 the second amino acid mostly moves phenylalanine from aldehyde to amine at a roughly
constant total. **In air, equimolar Leu/Ile/Val/Met halves phenylalanine's total degradation**
(−47 to −57 %). That is what competition for the oxidised carbonyl would produce. It is an
inference: the paper does not measure what the second amino acid becomes.

### 3b. Table 2 — Phe + arginine + lipid oxidation product (N2, 140 °C, 1 h, pH 3) [M]

Letters: means in a column that share no letter differ (p < 0.05).

| lipid oxidation product | PEA | PAC |
|---|---:|---:|
| None | 3.9 ± 0.8 kl | 0.1 ± 0.2 hi |
| 13-Hydroperoxyoctadeca-9,11-dienoic acid | 64.1 ± 0.3 cd | 9.2 ± 0.3 cdgh |
| Methyl 13-hydroperoxyoctadeca-9,11-dienoate | 70.3 ± 4.6 bd | 17.6 ± 2.2 bdej |
| Methyl 13-hydroperoxyoctadeca-9,11,15-trienoate | 84.4 ± 5.6 ab | 46.9 ± 4.3 a |
| 2-Pentenal | 13.5 ± 3.9 jl | 7.3 ± 2.6 cgij |
| 2-Octenal | 60.8 ± 6.2 cde | 4.9 ± 1.3 dik |
| 2,4-Hexadienal | 47.0 ± 4.3 efg | 16.5 ± 1.9 cef |
| 2,4-Heptadienal | 72.8 ± 3.5 bc | 19.4 ± 3.4 bc |
| 2,4-Decadienal | 89.6 ± 7.2 a | 18.6 ± 3.3 bc |
| 4-Oxo-2-hexenal | 20.2 ± 0.5 hij | 15.9 ± 1.9 bfgk |
| 4-Oxo-2-nonenal | 31.7 ± 4.3 gi | 29.4 ± 1.2 b |
| 4,5-Epoxy-2-heptenal | 31.9 ± 4.9 gh | 20.7 ± 0.5 bc |
| 4,5-Epoxy-2-decenal | 57.9 ± 4.8 df | 46.7 ± 9.2 a |
| 4-Hydroxy-2-nonenal | 14.3 ± 0.8 jk | 0.3 ± 0.1 hi |

- The 2,4-decadienal row is **identical** to Table 1's arginine/N2 cells (89.6 ± 7.2; 18.6 ± 3.3):
  the same experiment printed twice [Q].
- Text: DD and LnOOMe raise PEA "more than twenty times" (89.6/3.9 = 23; 84.4/3.9 = 22 [Z]) ✓;
  13-HPODE, 2-pentenal, 2-octenal and HNE do not raise PAC significantly.
- **[!]** Text: LnOOMe and 4,5-epoxy-2-decenal raise PAC "more than seven hundred times". From the
  printed means: 46.9/0.1 = 469 and 46.7/0.1 = 467 [Z]. The control value 0.1 ± 0.2 is rounded,
  so the ">700" cannot be checked. Do not quote a fold-change from this table.

### 3c. Dose series (Figs 1–3; Phe + DD, N2, 140 °C, 1 h) [F] with printed fits [M]

| series | printed statement |
|---|---|
| Arginine 0–10 µmol (Fig. 1) | PEA linear 0–6 µmol (r = 0.9993, p = 0.0007), then rises more slowly; PAC falls linearly over 0–10 µmol (r = 0.997, p < 0.0001) |
| Lysine | PEA linear 0–4 µmol (r = 0.9991, p = 0.023); PAC falls linearly 0–4 µmol (r = 0.993, p = 0.07) |
| Cysteine | PEA linear 0–10 µmol (r = 0.998, p < 0.0001); PAC falls linearly 0–8 µmol (r = 0.969, p = 0.001) |
| Methylguanidine 0–2 µmol (Fig. 2) | PEA and PAC both rise linearly (r = 0.998 and 0.9998, p < 0.0001); more active than Arg for PEA; **raises** PAC |
| Butylamine 0–2 µmol (Fig. 2) | PEA rises (r = 0.998), PAC falls (r = −0.999), both p < 0.0001 |
| Additivity (Fig. 2) | methylguanidine + butylamine effects summed give a line with slope "quite similar" to arginine's |
| Fig. 3 | butylamine slope for PEA ≈ cysteine's, below lysine's; 2× butylamine ≈ lysine (two amino groups). For PAC, Lys and Cys slopes similar and smaller than butylamine's; the authors suggest the second amino acid also reacts with the formed PAC (e.g. Cys thiazolidines) |

Fig. 1 axes [F]: PEA 20 → ~90 µmol/mmol Phe at 10 µmol Arg; PAC ~21 → ~4 at 8–10 µmol Cys. Not
digitised.

### 3d. Activation energies (Fig. 4; Phe + DD ± second AA, N2, pH 3) [M]

| system | Ea PEA (kJ/mol) | Ea PAC (kJ/mol) | location |
|---|---:|---:|---|
| Phe + DD (no second AA) | **51** | **78** | §3.4, Figs S-1A/B |
| + arginine | **51** | **74** | §3.4, Figs S-2A/B |
| + lysine | **54** | **71** | §3.4, Figs S-3A/B |
| + cysteine | **47** | **81** | §3.4, Figs S-4A/B |
| summary (Discussion) | 47–54 | 71–81 | p. 261 |

- Temperatures: ≈100–180 °C, five points, read off the Fig. 4 abscissa [F]. The text gives no
  temperature list. The k values are only in Supplementary figures, which are not on disk; the
  ln k unit is not stated.
- Rough check against the plotted end-points of the control line [Z]: PEA Δln k ≈ 2.9 over
  Δ(1/T) ≈ 0.00047 K⁻¹ → ~51 kJ/mol ✓; PAC Δln k ≈ 4.7 → ~83 kJ/mol (printed 78; consistent within
  reading error).
- **[!] Copy error:** the cysteine paragraph says "in the presence of 2,4-decadienal and **lysine**";
  the context (and Fig. 4 squares = cysteine) shows it is cysteine.
- **[!] Cross-paper:** Zamora 2015 (same lab and protocol, DD + Phe, pH 3, N2) printed PEA 51 / PAC
  **74** kJ/mol. Here the no-second-AA control is 51 / **78**, and **74** is the arginine system.
  Zamora 2013 (buffer, not sand) printed DD + Phe under N2 = 78.0. Either the 2015 PAC value is the
  arginine series, or the 74 vs 78 difference is replicate noise. Treat PAC Ea under N2 as **74–78**,
  not two independent confirmations.
- **No Ea was measured for any system with Val, Leu, Ile or Met.**

---

## 4. WHAT THE MODEL CAN TAKE FROM IT

**The question: does it quantify Strecker aldehydes or amines of leucine, isoleucine, valine or
methionine?** **No.** The four meat-aroma amino acids were each run as the *second* amino acid in
Table 1 only (DD, 140 °C, 1 h, pH 3, N2 and air). The analytes were PEA and PAC alone (ions 91 and
122). There is **no 3-methylbutanal, 2-methylbutanal, 2-methylpropanal or methional value, no
isoamylamine, 2-methylbutylamine, isobutylamine or 3-(methylthio)propylamine value, and no Ea**
for any of them. The Ea work used Arg, Lys and Cys only.

What it does give:

1. **Competition evidence [M → Z].** Equimolar Leu, Ile, Val or Met removes 56–67 % of
   phenylalanine's PAC in air and ~50 % of its total degradation. So these amino acids take a
   phenylalanine-sized share of the oxidised decadienal, which supports (does not prove) the
   assumption that a lipid-driven Strecker rate constant for Phe can be applied to Leu/Ile/Val/Met
   **within a band**. A row built this way is an **assumption with a band**, not a measured value.
   Closing it needs a primary paper that measures 3-methylbutanal or methional from a lipid
   carbonyl. See `hidalgo2015_extraction.md` §4: the review cites none.
2. **Matrix effect for a multi-amino-acid matrix [M].** In a protein hydrolysate or a plant-based
   patty, free amino acids are a mixture. In this model system every amino acid tested lowers Phe's
   Strecker aldehyde significantly, by up to 90 %, except Ser, Asp, Glu and Arg under N2 and Ser in
   air (the paper's own list). Any lipid-driven Strecker step
   calibrated on single-amino-acid systems (Zamora 2007/2008/2013/2015) **over-predicts per-amino-acid
   aldehyde in a mixture**. That is a directional claim the panel can carry: "adding a second amino
   acid lowers the lipid-driven Strecker aldehyde of the first".
3. **Ea robustness [M].** The PEA (47–54) and PAC (71–81 kJ/mol, N2) barriers do not move when
   Arg, Lys or Cys is added. Together with Zamora 2013 (DD in air 27.6, N2 78.0), this says the
   **oxygen state, not the amino-acid mix, sets the lipid-Strecker barrier**. A row built on these
   values needs an atmosphere flag.
4. **Limits:** pH 3 only (a plant-based patty is ~pH 6, where Zamora 2015 shows PAC 5–10× lower);
   sand slurry with nominal 91 mM reactants; 1 h at 140 °C; nitrogen for every Ea.
5. Nothing here is currently scoreable: the kinetic core has no phenylalanine species and no
   lipid-driven Strecker step.
