# Tai & Ho 1998 (JAFC) — EXTRACTION (glutathione or glutathionesulfonic acid 0.0012 mol + D-glucose (amount not printed), 50 mL water, initial pH 3.0 / 6.0 / 8.0, sealed 120 mL cylinder, 160 °C oven, 2 h; two end-point tables of solvent-extracted volatiles in mg per mol against a tridecane internal standard)

### THE ONE pH LADDER FOR SULFUR VOLATILES FROM A GLUTATHIONE + SUGAR POT IN THE CORPUS, and a clean negative control: the sulfonic-acid (irreversibly oxidised) glutathione gives NO sulfur volatile at any pH. Reduced GSH + glucose: total sulfur volatiles at pH 6 ≈ pH 8 ≈ 2.1× pH 3; thiofuraneol 240-259 mg/mol dominates; FFT 9-15 mg/mol, nearly flat in pH; **no MFT, no 2-mercapto-3-pentanone, no open-chain disulfide, no H2S, no free GSH or cysteine, no time course.** It does NOT study the disulfide (GSSG) — the title's "glutathione oxidation" is the sulfonic acid. Read 2026-10-31.

**Source on disk:** `data/articles/tai1998.pdf` (6 pp., born-digital ACS PDF, clean text layer;
publisher download stamp in the margin). Prose and both tables read from `pdftotext -layout`;
**Tables 1 and 2 re-checked cell by cell against 120 dpi renders of pp. 2262 and 2264** (the
polysulfide rows, whose pH-3 column is blank, are placed by the render and by the printed totals).
There are no figures. No Supporting Information exists for this article.

**Provenance codes:** **[M]** measured and printed · **[V]** verbatim quotation (short) ·
**[C]** cited by the authors from another paper · **[Z]** derived by this reading, never printed ·
**[Q]** qualitative statement. Nothing in this paper is digitised or fitted.

| field | value |
|---|---|
| Title | "Influence of Glutathione Oxidation and pH on Thermal Formation of Maillard-Type Volatile Compounds" |
| Authors | Chao-Ying Tai, **Chi-Tang Ho** (corresponding, ho@aesop.rutgers.edu) |
| Affiliation | Department of Food Science, Cook College, Rutgers, The State University, New Brunswick NJ |
| Venue | J. Agric. Food Chem. **1998, 46 (6), 2260-2265**. Received 30 Dec 1997, revised 2 Apr 1998, accepted 3 Apr 1998, published on the Web 15 May 1998 |
| DOI as printed | **10.1021/jf971111t** — no "DOI" line printed (1998 ACS style); manuscript code **S0021-8561(97)01111-4**, article id **JF971111T** (p. 2265), publisher stamp `jf971111t.pdf`. ✔ matches the brief |
| Companion | `tai1997_extraction.md` (same lab / vessel / oven / 160 °C / 2 h / extraction; cysteine, cysteinesulfinic and cysteic acid + glucose at pH 6) |
| Repo history | `results/literature/deep_research_backlog.json` carried this DOI under the confabulated string "Lee, Jo & Kim (1998), JAFC 46:1388"; its 2026-08-26 repair to "Tai & Ho (1998), JAFC 46:2260-2265" is **confirmed** by this reading (Flag 9) |

## 1. What this paper answers for the model

**The brief asked for** yields / time courses of thiols (MFT, FFT, 2-mercapto-3-pentanone),
thiophenes, thiazoles, H2S and disulfides, their dependence on cysteine vs cystine and on pH, with
conditions and quantitation. What the paper holds:

| asked for | in this paper? | what there is |
|---|---|---|
| **MFT** (2-methyl-3-furanthiol) | **NO** | not reported at any pH (glucose is a poor MFT parent; cf. `hofmann1998_reconciliation.md`, `zhang1991b_extraction.md`) |
| **FFT** (2-furfurylthiol) | **yes** | **11.9 / 14.7 / 9.1 mg/mol** at pH 8 / 6 / 3 [M] = 1.0 / 1.3 / 0.8 × 10⁻⁴ mol per mol GSH [Z] |
| other thiols | **yes** | 5-methyl-2-furfurylthiol 47.8 / 47.3 / 27.8; 2-mercapto-3-butanone 7.2 / 8.1 / 3.3; 2-thiophenethiol 1.6 / 3.5 / 8.7; 5-methyl-2-thienylmethanethiol 2.6 / 3.6 / nd; 3-mercaptopropionic acid nd / nd / **26.8** mg/mol [M] |
| 2-mercapto-3-pentanone | **NO** | — |
| thiophenes | **yes** | total **326.8 / 311.6 / 93.1** mg/mol [M]; thiofuraneol (2,5-dimethyl-4-hydroxy-3(2H)-thiophenone) **259.3 / 240.4 / 40.3** |
| thiazoles | **yes** | total 11.6 / 11.7 / 4.7 mg/mol [M] |
| disulfides | **cyclic polysulfides only** | total 8.1 / 9.0 / **nd** mg/mol [M]; no bis(2-methyl-3-furyl) or bis(2-furfuryl) disulfide, no GSSG |
| H2S | **NO** | argued from Zheng & Ho 1994's activation energies, mis-cited (Flag 5) [C] |
| GSH vs **GSSG / cystine** | **NO** | the oxidised arm is **glutathionesulfonic acid** (Cys residue → cysteic acid). The disulfide is mentioned only in the introduction [Q] |
| pH dependence | **yes — initial pH only** | 3.0 / 6.0 / 8.0 set at the start, unbuffered, reagent not named, never re-measured |
| time course | **NO** | one end point (2 h, oven setpoint 160 °C) |
| free GSH / cysteine / residual sugar | **NO** | — |
| quantitation | **semi-quantitative** | FID area vs **one** internal standard (tridecane, 0.25 mL of 1000 ppm), response factors not stated (implicitly 1); unit "milligrams per mole" with the denominator **not stated** (Flag 1); no replicates, no SD; identities by MS (Wiley 138) **plus DB-1 retention index** |

**Sulfur recovery** (assuming mg per mol of GSH): Σ (mg/mol ÷ MW × n_S) over the 26 sulfur compounds
= **3.38 / 3.36 / 1.57 mmol S per mol GSH = 0.34 / 0.34 / 0.16 %** of the GSH sulfur at pH 8 / 6 / 3
[Z]. As in Tai 1997 (≈ 0.18 %), the volatile fraction is **< 0.5 % of the charged sulfur**.

**Implications for the two open questions:**

- **Free cysteine 2-40× too high in fed-TTCA pots (missing sink or oxidation?).** Nothing here
  measures a thiol pool, so the paper **cannot decide it** (confidence 0.95). Two things carry over.
  (i) Sulfur volatiles are a negligible sink (< 0.5 % of S); whatever removes the cysteine is a
  non-volatile product (cystine, cysteine-carbonyl adducts, melanoidin-bound S, H2S lost to the
  headspace), which is the content of `r_cys_thermal`'s lump. (ii) The sulfur volatiles **rise** from
  pH 3 to pH 6-8 (×2.1). Thiol autoxidation runs through the thiolate and accelerates with pH, so if
  in-pot oxidation were competing strongly for the sulfur donor, the high-pH arms would be the
  depleted ones; they are not. That is weak evidence (one end point, an O2 supply that was not
  controlled, H2S release itself pH-dependent) against oxidation being a dominant thiol loss in
  **sealed** Maillard pots (my weight: likelihood ratio ≈ 1.5 against, not more).
- **Thiol oxidation / disulfide formation as a loss channel** (`kinetic_core_thiol_submodel_prereg.md`
  sec. 14 gives "oxidation / measurement artefacts" a prior of 0.04 for the MFT lag). The paper shows
  that **irreversible** S-oxidation (sulfonic acid) abolishes the sulfur chemistry completely at every
  pH [M]: zero sulfur compounds out of 31-38 identified. It does not show that such oxidation happens
  in situ, and it gives no disulfide data. It neither raises nor lowers that 0.04.
- **Headspace O2 arithmetic (§3, item 6)**: this vessel holds enough air to take **all** of its 24 mM
  GSH to the disulfide (capacity ≈ 48 mM thiol), whereas the 0.48 M cysteine pots of Tai 1997 could lose
  at most 10 %. Any oxidation rival in the engine is therefore a **charge-dependent** effect (large in
  dilute pots, negligible at 0.4 M), which is a usable signature: Cao 2017 / Zhao 2017 (0.4 M) fit while
  the 5-10 mM Zhai pots do not. Against it: the same prereg's Zhai 2021 Fig. 6 excess is **largest at
  pH 5.5 (+1.55 dex) and smallest at pH 8 (+0.72)**, the opposite of a thiolate-driven oxidation. My
  probability that thiol → disulfide oxidation is the main cause of the fed-TTCA cysteine excess:
  **≈ 0.2** (the vessel headspaces of the Zhai pots are not stated, which is what would move it).

**Fit rows: none.** **Hold-out / directional candidates: Table 1's sulfur rows at three pH** (§4) —
worth a cheap two-reader keying only if a glucose + GSH pot is inside the engine's scope.

## 2. Methods as they matter to a model [M]

- **Charge.** "Glutathione or glutathionesulfonic acid (**0.0012 mol**) ... with α-D-glucose ... was
  dissolved in **50 mL** of distilled water" → GSH **24 mmol/L** [Z]. **The glucose amount is not
  printed** (Flag 2); the same lab's Zhang & Ho 1991b used equimolar charges.
- **pH.** "adjusted to pH, **3.0, 6.0, or 8.0**" [V] — reagent not named, **unbuffered**, final pH not
  measured. GSH carries two carboxyls and an amine; glutathionesulfonic acid adds a strong acid, so the
  base load at a given setpoint differs between the tables [Z].
- **Vessel and heating.** "sealed in a **120 mL cylinder reactor**", "heated in a **160 °C oven for
  2 h**", ice-water quench. Headspace ≈ 70 mL air [Z]; material and come-up time not stated.
- **Isolation.** Whole mixture + **0.25 mL tridecane, 1000 ppm** (≈ 250 µg if the solvent density is
  ~1; solvent not named); CH2Cl2 **4 × 50 mL**; Na2SO4; Kuderna-Danish to 5 mL; N2 to 1 mL.
- **GC.** Varian 3400, FID, **DB-1 60 m × 0.32 mm × 1 µm**, split 25 : 1, 40 → 280 °C at 2 °C/min,
  He 1 mL/min. **GC-MS** HP 5790A + HP 5970A, EI 70 eV. Identification by **Wiley 138 library + retention
  index vs authentic compounds or literature**.
- **Replication.** None stated.

## 3. Tables re-typed (verbatim) [M]

### Table 1. "Volatiles (Milligrams per Mole) Obtained from the Reaction of Glutathione and Glucose"

| class | compound | RI (DB-1) | pH 8.0 | pH 6.0 | pH 3.0 |
|---|---|---|---|---|---|
| carbonyls | 2,3-butanedione | 625 | 3.0 | 2.2 | 0.9 |
| | 2-butanone | 631 | 61.3 | 36.0 | 7.3 |
| | 1-hydroxy-2-propanone | 668 | 31.1 | 25.8 | 21.6 |
| | 5-hydroxy-2-pentanone | 687 | 2.3 | | |
| | 2-pentanone | 693 | 4.8 | 3.6 | |
| | 3-pentanone | 701 | 2.7 | 2.4 | |
| | acetoin | 705 | 16.6 | 13.0 | 4.8 |
| | 1-hydroxy-2-butanone | 748 | 5.7 | 5.3 | |
| | 3-hexanone | 771 | 1.8 | | |
| | 2-mercapto-3-butanone | 789 | 7.2 | 8.1 | 3.3 |
| | 2,4-hexanedione | 862 | 19.6 | 24.8 | 3.1 |
| | 2,5-hexanedione | 888 | | | 14.8 |
| | 3-methyl-3-penten-2-one | 902 | | | 1.8 |
| | 2-methyl-1,3-cyclohexanedione | 1006 | 19.6 | 22.4 | 12.5 |
| | cyclotene | 1006 | 19.0 | 4.3 | 10.7 |
| | **total carbonyls** | | 194.7 | 147.9 | 80.8 |
| furans | 2,5-dimethylfuran | 715 | 11.5 | 11.5 | 3.6 |
| | dihydro-2-methyl-3(2H)-furanone | 781 | 3.9 | 3.7 | 1.4 |
| | 2-furfuryl alcohol | 835 | 30.8 | 32.3 | 5.1 |
| | 2-acetylfuran | 884 | 29.2 | 44.9 | 43.4 |
| | **2-furfurylthiol** | 888 | **11.9** | **14.7** | **9.1** |
| | 2,5-dimethyl-3(2H)-furanone | 917 | 3.1 | | |
| | 5-methyl-2-furfural | 934 | 9.5 | 17.7 | 94.7 |
| | **5-methyl-2-furfurylthiol** | 993 | **47.8** | **47.3** | **27.8** |
| | 1-(2-furyl)-1,2-propanedione | 1023 | 3.5 | 4.3 | |
| | 2,5-dimethyl-4-hydroxy-3(2H)-furanone | 1046 | 10.9 | 14.9 | 9.5 |
| | 5-hexyldihydro-2(3H)-furanone | 1144 | 8.8 | 7.9 | 21.0 |
| | 2-furancarbothioic acid, S-methyl ester | 1182 | | | 5.7 |
| | 5-(hydroxymethyl)-2-furfural | 1196 | | | 30.5 |
| | 5-acetyl-2-(2-furfuryl)furan | 1596 | 7.6 | 3.7 | |
| | **total furans** | | 178.5 | 202.9 | 242.3 |
| thiophenes | thiophene | 685 | 4.1 | 0.4 | 6.2 |
| | 2-methylthiophene | 767 | 11.2 | 12.4 | 9.5 |
| | 3-ethylthiophene | 854 | 2.5 | | |
| | 2,5-dimethylthiophene | 859 | 2.0 | 2.3 | 2.5 |
| | tetrahydrothiophene-3-one | 914 | 7.6 | 6.6 | 3.0 |
| | 2-thiophenethiol | 953 | 1.6 | 3.5 | 8.7 |
| | 2-methyltetrahydrothiophen-3-one | 958 | 13.9 | 23.5 | 12.2 |
| | 2-thiophenecarboxaldehyde | 960 | | | 5.6 |
| | 2-acetylthiophene | 1067 | 8.9 | 8.9 | 5.1 |
| | 2-formyl-5-methylthiophene | 1095 | | 6.0 | |
| | **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** | 1158 | **259.3** | **240.4** | **40.3** |
| | 5-methyl-2-thienylmethanethiol | 1183 | 2.6 | 3.6 | |
| | 3-thiophenecarboxylic acid | 1199 | 9.2 | | |
| | 5-methyl-2-thiophenecarboxylic acid | 1255 | 3.9 | 4.0 | |
| | **total thiophenes** | | 326.8 | 311.6 | 93.1 |
| thiazoles | thiazole | 726 | 6.0 | 5.9 | 4.7 |
| | 3-methylisothiazole | 796 | 2.0 | 2.4 | |
| | 2,4,5-trimethylthiazole | 981 | 3.6 | 3.4 | |
| | **total thiazoles** | | 11.6 | 11.7 | 4.7 |
| polysulfides | 3,5-dimethyl-1,2,4-trithiolane | 1111 | 5.5 | 5.2 | |
| | 3,6-dimethyl-1,2,4,5-tetrathiane | 1401 | 0.2 | 2.7 | |
| | 4,7-dimethyl-1,2,3,5,6-pentathiepane | 1667 | 2.4 | 1.1 | |
| | **total cyclic polysulfides** | | 8.1 | 9.0 | |
| pyrans | 3,4-dihydro-6-methyl-2H-pyran | 769 | 4.8 | 12.4 | |
| | 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | 1124 | 6.2 | 16.9 | 12.3 |
| | 3,5-dihydroxy-2-methyl-4H-pyran-4-one | 1178 | 6.0 | 7.0 | 10.0 |
| | **total pyrans** | | 17.0 | 30.3 | 10.5 |
| N-heterocyclics | pyridine | 738 | 2.0 | 2.6 | |
| | methylpyrazine | 803 | 3.2 | 2.1 | |
| | 2-acetylpyrrole | 1043 | 28.2 | 25.6 | 10.5 |
| | **total N-containing heterocyclics** | | 33.4 | 30.3 | 10.5 |
| others | acetic acid | 650 | 1.3 | 7.6 | 5.3 |
| | ethanethioic acid, S-methyl ester | 758 | 10.6 | 12.5 | 3.9 |
| | 3-mercaptopropionic acid | 973 | | | 26.8 |

### Table 2. "Volatiles (Milligrams per Mole) Obtained from the Reaction of Glutathionesulfonic Acid and Glucose"

| class | compound | RI (DB-1) | pH 8.0 | pH 6.0 | pH 3.0 |
|---|---|---|---|---|---|
| carbonyls | 3-buten-2-one | 621 | 6.3 | 11.7 | 3.9 |
| | 2-butanone | 631 | 5.1 | 13.8 | 14.9 |
| | 1-hydroxy-2-propanone | 668 | 165.1 | 58.3 | 8.9 |
| | 2-pentanone | 680 | 6.7 | 10.6 | 13.9 |
| | 3-methyl-3-buten-2-one | 693 | 22.0 | | 2.0 |
| | acetoin | 705 | 31.7 | 13.9 | 2.6 |
| | 4-methyl-2-pentanone | 732 | 2.4 | 3.7 | 6.8 |
| | 1-hydroxy-2-butanone | 748 | 1.0 | 2.8 | |
| | 4-methyl-3-penten-2-one | 780 | 3.4 | 7.9 | 11.4 |
| | 5-methyl-5-hexen-2-one | 792 | | | 2.2 |
| | 2-cyclopenten-1-one | 799 | | | 2.8 |
| | 4-hydroxy-4-methyl-2-pentanone | 814 | | 4.4 | 10.0 |
| | 1-(acetyloxy)-2-propanone | 826 | 4.5 | | |
| | cyclopent-2-ene-1,4-dione | 846 | 1.1 | | |
| | cyclohexanone | 900 | | | 1.5 |
| | 3,4-dimethyl-3-penten-2-one | 914 | 0.7 | | |
| | 2-methylcyclopentanone | 970 | 2.8 | | |
| | 4-propyl-1,3-cyclohexandione | 981 | | | 5.9 |
| | cyclotene | 1006 | 7.7 | 10.1 | 5.5 |
| | **total carbonyls** | | 260.5 | 137.2 | 92.3 |
| furans | 2-methylfuran | 737 | 2.6 | | 5.2 |
| | dihydro-2-methyl-3(2H)-furanone | 781 | 3.7 | 1.6 | |
| | 2-furfural | 803 | 42.7 | 140.0 | 98.2 |
| | 2-furfuryl alcohol | 835 | 23.0 | 52.3 | |
| | 2(5H)-furanone | 863 | 0.7 | 12.6 | 7.9 |
| | 2-acetylfuran | 884 | 2.3 | 2.5 | 4.6 |
| | 2-(1,1-dimethylethyl)-4-methylfuran | 915 | | | 16.5 |
| | 5-methyl-2-furfural | 934 | 5.6 | 1.3 | 3.0 |
| | 5-ethyl-2-furfural | 1032 | 18.6 | 9.8 | |
| | 1-(2-furyl)-2-hydroxyethanone | 1056 | 7.2 | 9.6 | 4.1 |
| | 5-(hydroxymethyl)-2-furfural | 1195 | 159.4 | 1244.9 | 895.9 |
| | **total furans** | | 265.8 | 1474.6 | 1035.4 |
| pyrans | 3,6-dihydro-2H-pyran-2-one | 883 | | | 4.1 |
| | 3-hydroxy-2-methyl-4H-pyran-4-one | 1085 | | 5.1 | 5.1 |
| | 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | 1124 | | 46.5 | 14.5 |
| | 3,5-hydroxy-2-methyl-4H-pyran-4-one [sic] | 1178 | 2.1 | 6.0 | 5.2 |
| | **total pyrans** | | 2.1 | 57.6 | 28.9 |
| pyrazines | pyrazine | 722 | 7.0 | 1.3 | |
| | methylpyrazine | 803 | 16.2 | 5.0 | |
| | 2,5-dimethylpyrazine | 886 | 64.5 | 13.4 | |
| | ethylpyrazine | 889 | 1.5 | | |
| | 2,3-dimethylpyrazine | 893 | 2.5 | | |
| | trimethylpyrazine | 980 | 13.1 | 6.0 | |
| | **total pyrazines** | | 104.8 | 25.7 | |
| pyrroles | 1H-pyrrole-2-carboxaldehyde | 964 | 2.1 | 7.9 | |
| | 2,5-pyrrolidinedione | 1084 | 15.6 | 80.4 | 71.3 |
| | 5-methyl-1H-pyrrole-2-carboxaldehyde | 1099 | | 3.7 | 3.7 |
| | **total pyrroles** | | 17.7 | 92.0 | 75.0 |
| others | methylbutenol | 640 | 40.6 | | 1.3 |
| | acetic acid | 650 | 7.8 | | |
| | 2-methyl-2-butanol | 658 | 10.5 | 15.2 | 10.2 |
| | 3-penten-2-ol | 694 | 39.8 | 2.2 | 4.3 |
| | ethyl acetate | 766 | 1.4 | | |

**Table 2 has no sulfur-containing compound at any pH** [M] (the text: "could not produce
sulfur-containing compounds").

### Values printed in the text [M]

| quantity | value | where |
|---|---|---|
| compounds identified, GSH + Glc | **37 / 47 / 51** at pH 3.0 / 6.0 / 8.0 | p. 2261 |
| compounds identified, GSH-SO3H + Glc | **31 / 31 / 38** at pH 3.0 / 6.0 / 8.0 | p. 2263 |
| HMF share of total furans, GSH-SO3H | **87 / 84 / 60 %** at pH 3.0 / 6.0 / 8.0 | p. 2263 |
| total thiophenes | 93.1 (pH 3), 311.6 (pH 6), 326.8 (pH 8) mg/mol | p. 2261 |
| 5-methyl-2-furfural highest at pH 3 | 94.7 mg/mol | p. 2261 |
| odour thresholds of furaneol and thiofuraneol | "identical, at the value of 1 ng/L in air" [C Blank & Schieberle 1993] | p. 2263 |
| Ea of H2S release from GSH (as cited) | pH 9.0 17.4 < pH 3.0 18.5 < pH 7.0 22.1 < pH 5.0 30.9 kcal/mol [C Zheng & Ho 1994] — **mis-cited, Flag 5** | p. 2261 |

### Arithmetic [Z] (`scratchpad tai_scripts/tai_arith.py`, run in the maillard env; not shipped)

1. **Counts reproduce:** Table 1 has 51 / 47 / 37 entries at pH 8 / 6 / 3, Table 2 38 / 31 / 31 — all
   six match the text.
2. **Printed totals vs the rows.** All carbonyl, thiophene, thiazole, polysulfide and N-heterocycle
   totals reproduce, and all Table 2 totals reproduce. **Three Table 1 totals do not** (Flag 3):
   pyrans pH 6 rows sum **36.3** vs printed 30.3; pyrans pH 3 rows sum **22.3** vs printed 10.5 (= the
   N-heterocycle total of the same column); furans pH 3 rows sum **251.8** vs printed 242.3 (the
   difference, 9.5, is exactly the pH-3 furaneol cell).
3. **Molar yields per mol GSH** (MW: FFT 114.17, 5-methyl-FFT 128.19, thiofuraneol 144.19,
   2-mercapto-3-butanone 104.17, 2-thiophenethiol 116.20, 3-mercaptopropionic acid 106.14):

   | compound | pH 8.0 | pH 6.0 | pH 3.0 |
   |---|---|---|---|
   | FFT | 1.04e-4 | 1.29e-4 | 0.80e-4 |
   | 5-methyl-2-furfurylthiol | 3.73e-4 | 3.69e-4 | 2.17e-4 |
   | thiofuraneol | 1.80e-3 | 1.67e-3 | 2.79e-4 |
   | 2-mercapto-3-butanone | 6.9e-5 | 7.8e-5 | 3.2e-5 |
   | 2-thiophenethiol | 1.4e-5 | 3.0e-5 | 7.5e-5 |
   | 3-mercaptopropionic acid | nd | nd | 2.5e-4 |

4. **Sulfur in volatiles:** 3.38 / 3.36 / 1.57 mmol S per mol GSH (0.34 / 0.34 / 0.16 %) at
   pH 8 / 6 / 3, counting n_S per molecule (polysulfides 3-5 S, thiophenethiol and thienylmethanethiol
   2 S).
5. **Within-study pH ratios** (the only robust use of a single-IS table): FFT pH 6 / pH 3 = **1.62**,
   pH 8 / pH 3 = **1.31**; total thiophenes 3.35 / 3.51; thiofuraneol 5.97 / 6.43; 5-methyl-FFT 1.70 /
   1.72; 2-mercapto-3-butanone 2.45 / 2.18; 2-thiophenethiol 0.40 / 0.18 (the only sulfur compound
   that **falls** with pH); polysulfides present at pH 6 and 8, absent at pH 3.
6. **Thiofuraneol vs furaneol** (molar): 21 / 14 / 3.8 at pH 8 / 6 / 3 — the sulfur swap outruns its
   oxygen analogue by an order of magnitude at pH ≥ 6.
7. **Headspace oxygen budget.** 70 mL air → **0.60 mmol O2** = 12 mmol per L of liquid; at 4 thiol per
   O2 (→ disulfide) that oxidises up to **48 mM thiol**, i.e. all 24 mM GSH; to the sulfonic acid
   (1.5 O2 per thiol) it could take ~33 % of the GSH. In Tai 1997's 0.48 M cysteine pots the same
   headspace caps oxidation at 10 % of the cysteine.
8. **Oxidised vs reduced GSH at pH 3:** HMF 895.9 vs 30.5 mg/mol (**29×**); at pH 6, 1244.9 vs nd.

## 4. What the repository can use

| item | value | conditions | source | evidence class | use |
|---|---|---|---|---|---|
| FFT from GSH + Glc vs pH | 11.9 / 14.7 / 9.1 mg/mol; ratios pH6/3 = 1.62, pH8/3 = 1.31 | GSH 24 mM, Glc not stated, initial pH 8 / 6 / 3, 160 °C oven 2 h, sealed | Table 1 [M]; arithmetic 5 [Z] | semi_quantitative_single_IS, within_study_ratio | **hold-out candidate (ratio)** if a glucose + GSH pot is in scope |
| sulfur classes vs pH | thiophenes 326.8 / 311.6 / 93.1; thiazoles 11.6 / 11.7 / 4.7; polysulfides 8.1 / 9.0 / nd | same | Table 1 [M] | within_study_ratio | **directional claims**: S volatiles pH 6 ≈ pH 8 > pH 3; polysulfides absent at pH 3 |
| thiofuraneol (DMHFS) | 259.3 / 240.4 / 40.3 mg/mol | same | Table 1 [M] | semi_quantitative | context for `DMHFS` (rate fixed at zero by declaration); the only pH ladder for it in a sugar + thiol pot — **not** a rate |
| oxidised GSH (sulfonic acid) gives no S volatiles | zero of 31-38 compounds at every pH | GSH-SO3H 24 mM + Glc | Table 2 [M] | directional (absence) | negative control; nothing in the engine to score (no sulfonic species) |
| HMF from oxidised vs reduced GSH | 29× at pH 3 | same | Tables 1-2 [M] | within_study_ratio | context: a thiol / H2S suppresses the 3-deoxyosone → HMF route |
| S recovered in volatiles | 0.34 / 0.34 / 0.16 % | assumes mg per mol GSH | arithmetic 4 [Z] | derived_assumption | bound: volatiles are not the thiol sink |
| headspace O2 capacity | 48 mM thiol → disulfide | 120 mL vessel, 50 mL liquid | arithmetic 7 [Z] | derived | design input for the oxidation rival: charge-dependent |

**Two-reader keying.** Table 1's 26 sulfur rows × 3 pH are cheap to key twice and are the only
pH-resolved sulfur-volatile data from a GSH + sugar pot in the corpus. Key them **only as hold-out
ratios / directional claims, never as fit rows**: one lab, one end point, single internal standard
with unstated response factors, glucose charge unstated, initial pH only, an oven-heated vessel of
unknown thermal history, and 160 °C, above the sulfur-lane pots the engine is fitted on. Before
keying, check that the engine can simulate a glucose + GSH pot at all (GSH is carried only as an H2S
source; the hexose → furfural → FFT path must be active), otherwise these belong in the directional-
claims register and nowhere else. **Table 2: nothing to key** beyond the zero-sulfur claim.

## 5. Flags

1. **"Milligrams per mole" — of what is not stated.** Read here as per mole of GSH (the same lab's
   Zhang & Ho 1991 / 1991b state "per mol of Cys or GSH"). Per mole of glucose would give the same
   numbers only if the charge was equimolar (Flag 2).
2. **Glucose amount not printed** ("with α-D-glucose").
3. **Three printed totals in Table 1 disagree with their rows** (arithmetic 2): pyrans pH 6 (30.3 vs
   36.3), pyrans pH 3 (10.5 vs 22.3), furans pH 3 (242.3 vs 251.8). Use the rows, not the totals.
4. **Possible co-elution of FFT.** FFT is assigned RI **888** on DB-1, the same as 2,5-hexanedione
   (888, pH 3 only) and 4 units from 2-acetylfuran (884, 29-45 mg/mol in the same arms). FID
   quantitation of a 9-15 mg/mol peak next to a 30-45 mg/mol one is at risk; the FFT cells are the least
   secure of the thiol numbers.
5. **Zheng & Ho 1994 is mis-cited.** Tai & Ho quote GSH H2S-release Ea of 17.4 (pH 9), 18.5 (pH 3),
   22.1 (pH 7), 30.9 (pH 5) kcal/mol; Zheng & Ho's Table V prints **19.9, 18.8, 22.8, 30.8**
   (`zheng1994_extraction.md` §5). The ranking survives, the numbers do not. The inference drawn from
   them — H2S release "likely in the order of pH 3.0 > pH 8.0 > pH 6.0" [V] — is invalid as stated: an
   activation energy alone does not order rates (the pre-exponentials differ by ~10³ across pH in
   Zheng's own data). Take H2S-from-GSH kinetics from `zheng1994_extraction.md`, never from here.
6. **Oxidation-state naming reversed in the introduction:** "Glutathionesulfonic acid and
   glutathionesulfinic acid ... in which the cysteine residue is oxidized to cysteinesulfinic acid and
   cysteic acid, respectively" [V] — the sulfonic acid carries the cysteic-acid residue, the sulfinic the
   cysteinesulfinic. The experiments use the sulfonic acid only.
7. **2-mercapto-3-butanone (table) vs "3-mercapto-2-butanone" (text)** — the same molecule
   (CH3-CO-CH(SH)-CH3), named two ways.
8. **Same-lab conflict on FFT from glucose.** Zhang & Ho 1991b (glucose + GSH, 180 °C, 1 h, pH 7.5,
   SDE) report **no FFT and no MFT**; this paper (160 °C, 2 h, solvent extraction) reports FFT 9-15 and
   5-methyl-FFT 28-48 mg/mol. Temperature, time and isolation (SDE vs solvent extraction) differ, and
   both papers print DB-1 retention indices. Not registered as a formal conflict (different conditions),
   but any FFT hold-out from this table should be read with it.
9. **Repo citation history.** `results/literature/deep_research_backlog.json` repaired "Lee, Jo & Kim
   (1998), JAFC 46:1388" to this paper on 2026-08-26 with the anchor "zero sulfur compounds from
   GSH-SO3H" — **confirmed** (Table 2) — and noted that the corpus's "pH optimum 7.0 for sulfur
   heterocycles" is not what the paper says — **confirmed**: pH 6.0 and 8.0 are both favourable and
   statistically indistinguishable on this evidence; no pH 7 arm exists.
10. **No final pH, no replicates, no vessel material or come-up time, no residual GSH / glucose.**
11. **Registry gaps:** thiofuraneol, 5-methyl-2-furfurylthiol, 2-mercapto-3-butanone,
    glutathionesulfonic acid and the cyclic polysulfides are absent from `data/keys/compounds.yml`
    (FFT, MFT and glutathione are present).
