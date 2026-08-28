# Anantharamkrishnan, Hoye & Reineccius 2020 — COMPLETE TRANSCRIPTION
### *"Covalent Adduct Formation Between Flavor Compounds of Various Functional Group Classes and the Model Protein β-Lactoglobulin"*
### ⭐ **HEXANAL IS IN THE PANEL, AND IT FORMS COVALENT ADDUCTS.** — the paper `anantharamkrishnan2020_extraction.md` §8 nominated as the corpus's highest-value outstanding retrieval.

**Full extraction of every number in `data/articles/anantharamkrishnan2020b.pdf`.**
Wave **K4d**, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified. Companion to
`anantharamkrishnan2020_extraction.md` (the **later** pH/T/a_w paper, Wave K4b) — **read that
first**; this dossier does not re-transcribe its methods, and instead **discharges the
follow-up it named at line 389**.

**Provenance codes:** **[M]** measured and printed · **[C]** cited · **[Q]** qualitative
statement, no number attached · **[Z]** derived by this wave · **[!]** integrity flag.

---

## 0. PAPER IDENTITY — **MATCHES THE EXPECTED IDENTITY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Vaidhyanathan Anantharamkrishnan¹\*, Thomas Hoye²$, Gary A. Reineccius¹#** |
| Title | **"Covalent Adduct Formation Between Flavor Compounds of Various Functional Group Classes and the Model Protein β-Lactoglobulin"** |
| Section | *Chemistry and Biology of Aroma and Taste* |
| Venue | ***J. Agric. Food Chem.*** — **"Just Accepted Manuscript"**, Publication Date (Web) **09 May 2020**; downloaded from pubs.acs.org 12 May 2020 |
| DOI | **10.1021/acs.jafc.0c01925** |
| Affiliations | ¹ Dept. of **Food Science and Nutrition**, University of Minnesota, St. Paul, MN 55108 · ² Dept. of **Chemistry**, University of Minnesota, Minneapolis, MN 55455 |
| Funding | **Minnesota Agriculture Experiment Station** |
| Keywords as printed | *protein, flavor, covalent interaction, functional group, mass spectrometry, UPLC-TOF-MS* |
| PDF character | **39 pages** (cover + 38-page manuscript), iText 4.2.0, clean text layer, **all 9 figures are deconvoluted ESI mass spectra with no y-axis**. Watermark: *"Subscriber access provided by Uppsala universitetsbibliotek"* |

**⇒ CONFIRMED: this is the parent covalent-adduct survey across functional-group classes with
β-lactoglobulin, exactly as expected.**

**⚠️ ONE DISCREPANCY WITH THE BRIEF'S EXPECTED CITATION.** The brief expected *"JAFC 68:6395"*.
**This PDF is the "Just Accepted" version and carries NO volume or page numbers** — only the
DOI. It cannot confirm or refute 68:6395. `anantharamkrishnan2020_extraction.md` line 389
records the same citation as **"J. Agric. Food Chem. 68"**. **Treat the volume/page as
unverified from this file; the DOI 10.1021/acs.jafc.0c01925 is authoritative.**

**⚠️ [!] THIS PAPER USES A DIFFERENT PROTEIN FROM ITS COMPANION.**

| | **this paper (0c01925, May 2020)** | **companion (0c06752, Oct 2020)** |
|---|---|---|
| BLG variant | **B** | **A** |
| measured parent mass | **18,276 Da** | 18,362 ± 1 Da |
| cow | **number 2838** | number 2708 |
| purity | **> 90 %**, pH 6.7 | 95 % |

The companion paper states outright that *"the relative intensity of the adduct peak in 24 h is
**reduced in comparison to our previous publication** … because **BLG variant (A) is used in
this study whose reactivity … is different than the B variant** that was used in the previous
study."* **⇒ This paper is that "previous study", and it is the one with the STRONGER adduct
signals.** `anantharamkrishnan2020_extraction.md` line 399–400 anticipated exactly this and
called it *"the only way to bound the variant effect"*. **Every result here must carry
`BLG_variant: B`, and directional results here are the upper end of the A/B pair.**

---

## 1. THE ONE-PARAGRAPH ANSWER TO THE THREE BRIEF QUESTIONS

**(a) Is any n-alkanal — hexanal especially — among the compounds, and does it form covalent
adducts under food-relevant conditions? YES to the first, QUALIFIED-NO to the second.**
**Hexanal is compound 39 of 47 in Table 1 ("saturated aldehyde", MW 100) and is marked `yes`
in Table 2.** The results text is unambiguous: *"**Hexanal showed multiple adducts, presumably
engaging multiple lysines.**"* It is the **only straight-chain saturated alkanal in the panel**
(the other C6 entries are the enals). **But the conditions are not food-relevant**: 12 ppth
(1.2 % w/w) hexanal against 1 wt % BLG in water, **ambient temperature**, 24 h — a hexanal
loading **2,400–240,000× above** a food-relevant matrix level (§3c). **(b) Any rates, extents
or reversibility usable to bound the covalent sink? NO measured ones.** There is **no rate
constant, no half-life, no yield, no percent conversion, no equilibrium constant, no
calibration curve and no y-axis on any spectrum** in this paper. The only quantitative outputs
are **molecular masses in daltons** and a **binary yes/no reactivity table**. What *is*
extractable is (i) a **5-tier ordinal reactivity ladder** with hexanal placed inside it (§5a) —
**enals ≫ hexanal > benzaldehyde ≈ furfural > citral > vanillin = none** — which is the first
time the corpus has hexanal ranked against its own homologues on a single protein under a
single protocol; (ii) **two coarse timescale brackets** (enals < 10 min; the propanethiol 2:1
adduct between 6 h and 24 h, §5c); and (iii) an **order-of-magnitude lower bound on k₂ for the
hexanal–lysine Schiff-base channel** derived here (§6), which lands at **t½ ≈ 1 month to
2 years at ambient** in a 3 wt % protein matrix. **(c) Doses — food-relevant or saturating?
SATURATING BY DESIGN, and the authors say so:** *"The concentration of flavor used in this
study is much higher than that used in the industry. The reason for this is the ease and
accuracy of measuring the adducts…"* **⇒ VERDICT FOR THE REPO'S HEXANAL OVER-PREDICTION: this
paper establishes that the covalent-sink mechanism EXISTS for hexanal specifically — which
`anantharamkrishnan2020_extraction.md` §6 could not do, because hexanal was absent there — but
it CANNOT size the channel. On this paper's own numbers the ambient-temperature sink is 3–5
orders of magnitude too slow to explain a 36× over-prediction in a thermal process (§6c). The
quantitative paper the repo actually needs is this paper's own ref. 28: Meynier, Rampon,
Dalgalarrondo & Genot (2004), *"Hexanal and t-2-Hexenal Form Covalent Bonds with Whey Proteins
and Sodium Caseinate in Aqueous Solution"*, Int. Dairy J. 14(8), 681–690 (§10).**

---

## 2. SYSTEM COMPOSITION AND METHOD — complete

### 2a. The protein
| variable | value as printed |
|---|---|
| **Protein** | **β-lactoglobulin (BLG), single variant B**, isolated in-house |
| Source | milk from **a single homozygous cow, numbered 2838** (Acknowledgements) |
| Isolation | cream removed by centrifugation → casein precipitated at **pH 4.6** and filtered off → whey dialysed through a **12 kDa** membrane → freeze-dried → **Hi Scale 50/20 SEC column (50 mm ID × 200 mm L), Sephacryl S-200** → BLG fraction re-dialysed (12 kDa) → freeze-dried |
| **Purity** | **> 90 %** |
| **pH of the isolate** | **6.7** |
| **Measured parent mass** | **18,276 Da** (Fig. 1) |
| **Reactive residues per BLG** | **15 Lys · 1 free Cys · 2 Cys-disulfides (from 4 further Cys) · 3 Arg** |
| Variant chemistry | A and B differ at **Asp64Gly and Val118Ala**; **the variants' masses differ by 87 Da** |
| Why single variant | flavour molecules are **50–180 Da**, so a BLG+flavour adduct can **overlap in mass with the second variant** |

**[M] The control spectrum (Fig. 1) — the protein is NOT clean:**

| peak (Da) | assignment as printed |
|---|---|
| **18,276** | **BLG variant B** (major peak) |
| 18,292 | *(labelled, unassigned; +16 — likely a Met/Trp oxidation)* |
| 18,436 | *(labelled, unassigned; +160)* |
| **18,600** | **BLG + 1 lactose Schiff base** (+324 = 342 − 18) |
| 18,632 | *(labelled, unassigned; +32 on the lactose adduct)* |
| **18,924** | **BLG + 2 lactose** (+648) |

**⚠️ [!] The starting protein already carries mono- and di-lactosylated populations, on lysine,
as major peaks.** The authors attribute these to the raw milk (*"our sample of BLG was not
heated at any stage (e.g., the milk was not pasteurized)"*). **⇒ Some fraction of the 15
lysines is pre-blocked by a Maillard adduct before any flavour compound is added. The number of
free lysines is therefore < 15 and is never measured.** Every adduct count in this paper —
including hexanal's "multiple adducts" — is measured against a partially glycated substrate.

### 2b. The reaction system **[M]**
| variable | value as printed |
|---|---|
| Protein form | **1 weight-percent** freeze-dried BLG in **doubly distilled water** |
| Flavour dose | **12 ppth (parts per thousand by weight)**, added **individually**, vortexed in a **closed vessel** |
| **Temperature** | **ambient** (*"storage at ambient temperature"*; Table 1 caption: *"reaction at ambient temperature"*). **No number is given anywhere.** [!] |
| pH | **not adjusted and not stated for the reaction.** The isolate is pH 6.7 (§2a); the reaction pH is never given. [!] |
| Time points | **10 min to 24 h, "and, on occasion, extending to 48 h"** |
| **Table 2's reference time** | **24 h** (caption: *"their reactivity with BLG (**24 h reaction time**)"*) |
| Sample prep for MS | diluted **1:10 with water** |
| **Molar ratio flavour : BLG** | **"ca. 50–250"** as printed — **⚠️ see [!] below** |
| Replication | **NOT STATED anywhere. n is unknown for every result.** [!] |

**⚠️ [!] The printed molar ratio range is inconsistent with the printed dose.** [Z]
1 wt % BLG = 10 g L⁻¹ ÷ 18,276 Da = **0.547 mmol L⁻¹**. 12 ppth = 12 g L⁻¹. Over the panel's
MW range (allylamine 57 Da → citral diethyl acetal 226 Da) the flavour molarity runs
210 → 53.1 mmol L⁻¹, giving ratios of **385 → 97**, not "50–250". Using the stated > 90 %
purity makes it worse (**428 → 108**). **The true range is ≈ 100–390.** Use the computed
value, not the printed one. *(For hexanal specifically the two happen to agree well: 219,
inside "50–250".)*

### 2c. UPLC-ESI-MS/qTOF **[M]**
| parameter | value |
|---|---|
| Instrument | **Waters Acquity UPLC + Waters Synapt G2/Si HDMS q-oa-TOF** |
| Column | **Waters Acquity BEH C4, 2.1 × 100 mm, 1.7 μm**, at **35 °C** |
| Flow | **0.400 mL min⁻¹**; A = water + 0.1 % formic acid, B = MeCN + 0.1 % formic acid |
| Gradient | 3 % B 0→3 min; 3→97 % B 3→9 min; 97 % B 9→11 min; 97→3 % B 11→13 min; 3 % B 13→15 min (**15 min total**) |
| Acquisition | **profile mode, m/z 300–2500, every 0.2 s** |
| ESI (+) | capillary **0.5 kV**; sampling cone **35.0 V**; extraction cone **4.0 V**; desolvation gas **800 L h⁻¹**; source **100 °C**; desolvation **350 °C**; cone gas **40 L h⁻¹**; **trap CE off** |
| Lockspray | **0.5 μg mL⁻¹ leucine-enkephalin**, one spectrum (0.2 s, m/z 300–2500) **every 10 s**, three m/z averaged |
| Deconvolution | **Waters MaxEnt**, mass range **10,000–50,000**, resolution **0.1 Da/channel**, uniform Gaussian, **width at half height 0.33 Da**, auto-convergence |
| **Mass accuracy** | **± 1 Da** |

**⚠️ [!] The figures have NO y-axis — verified on a 600 dpi raster of Figs 1 and 4.** The
plotted axis is **"Molecular mass (Da)"** only. There is **no intensity scale, no units, no
normalisation and no baseline reference**, and no spectrum is normalised against any other.
MaxEnt deconvolution intensity is not proportional to molar amount (it depends on charge-state
distribution and on the deconvolution's own convergence), and **the authors never once use peak
heights**. **⇒ Peak-height ratios must NOT be digitised off these figures. There is no
quantitative content in any of the nine figures beyond the labelled masses.**

---

## 3. TABLE 1 — the full 47-compound panel **[M]**

*Caption verbatim: "Table 1. Aroma compounds analyzed for reaction with BLG (reaction at
ambient temperature)"*. Columns: functional-group class · compound · **reason for study** ·
molecular weight (Da). **[Z] Count check: 47 flavour compounds in 13 classes (excluding the
"No Flavor" control) — exactly matching the abstract's "forty-seven flavor compounds belonging
to thirteen different classes". Both counts verify.**

| class | compound | **reason for study** (as printed) | MW (Da) | **[Z] mM at 12 ppth** |
|---|---|---|---|---|
| No Flavor | pure protein (BLG) | control | 18,276 | — |
| Hydrocarbon | p-cymene | aromatic | 134 | 89.6 |
| Alcohols | (l)-menthol | terpene alcohol | 156 | 76.9 |
| | geraniol | double bond to alpha carbon | 154 | 77.9 |
| | 2-pentanol | secondary alcohol | 88 | 136.4 |
| | 2,3-butanediol | dihydroxy | 90 | 133.3 |
| Ketones | **diacetyl** | **diketone** | **86** | 139.5 |
| | 2-heptanone | carbonyl with six carbons | 114 | 105.3 |
| | 2-nonenone **⚠** | carbonyl with nine carbons | 140 | 85.7 |
| | cyclotene | cyclic, diketone | 112 | 107.1 |
| Acids | butyric acid | fatty acid | 88 | 136.4 |
| | 2-hydroxybenzoic acid | phenolic acid | 138 | 87.0 |
| | m-toluic acid | aromatic acid | 136 | 88.2 |
| Esters | methyl anthranilate | free amine group | 151 | 79.5 |
| | iso-amyl acetate | ester | 130 | 92.3 |
| | ethyl lactate | lactate | 118 | 101.7 |
| | methyl salicylate | phenolic group | 152 | 78.9 |
| Lactone | δ-dodecalactone | delta lactone (six-membered ring) | 198 | 60.6 |
| | γ-butyrolactone | gamma lactone (five-membered ring) | 86 | 139.5 |
| | γ-decalactone | gamma lactone (five-membered ring) | 170 | 70.6 |
| Amine Bases | 3-acetyl pyridine | pyridine | 121 | 99.2 |
| | allylamine | amine | **57** | **210.5** |
| | 2,5-dimethyl pyrazine | pyrazine | 108 | 111.1 |
| Sulfur-containing | 4-methyl-5-vinyl-thiazole | thiazole | 125 | 96.0 |
| | dimethyl sulfide | sulfide | 62 | 193.5 |
| | dimethyl disulfide | disulfide | 94 | 127.7 |
| | **dimethyl trisulfide** | trisulfide | **126** | 95.2 |
| | dimethyl sulfone | sulfone | 94 | 127.7 |
| | **propanethiol** | aliphatic thiol | **76** | 157.9 |
| | **thiophenol** | aromatic thiol | **110** | 109.1 |
| | allyl isothiocyanate | isothiocyanate | 99 | 121.2 |
| | 2-methyl thiophene | thiophene | 98 | 122.4 |
| Acetal | trans-2-hexenal dimethyl acetal | acetal | 144 | 83.3 |
| | citral diethyl acetal | acetal | **226** | **53.1** |
| Phenols | 2,6-dimethylphenol | dimethylphenol | 122 | 98.4 |
| | 4-vinylphenol | vinylphenol | 120 | 100.0 |
| | eugenol | allyl chain-substituted guaiacol | 206 | 58.3 |
| Furans | **furfuryl mercaptan (= 2-furfurylthiol, FFT)** | furan substituted with a thiomethyl group | **114** | **105.1** |
| | furaneol | 2,5-dimethylfuran carrying additional oxo and hydroxy groups | 128 | 93.8 |
| Pyranone | maltol | a 4-pyranone | 126 | 95.2 |
| **Aldehydes** | vanillin | phenolic aldehyde | 152 | 78.9 |
| | **⭐ hexanal** | **saturated aldehyde** | **100** | **119.8** |
| | **trans-2-hexenal** | conjugated enal | 98 | 122.4 |
| | **cis-3-hexenal** | nonconjugated enal | 98 | 122.4 |
| | **furfural** | furanyl aldehyde | 96 | 125.0 |
| | **trans,trans-2,4-heptadienal** | conjugated dienal | 110 | 109.1 |
| | citral | conjugated enal | 152 | 78.9 |
| | benzaldehyde | aromatic aldehyde | 106 | 113.2 |

**⚠️ [!] "2-nonenone" (MW 140) is a naming error.** MW 140 with "carbonyl with nine carbons"
is **2-nonanone** (C₉H₁₈O = 142.2, nearest integer 142) — and Table 2 groups it with the
unreactive monoketones, which is only consistent with the *saturated* ketone. A genuine
2-nonen-one would be an enone and, by this paper's own Michael-addition logic, should have
reacted. **Read as 2-nonanone.** *(Note the MW column is itself approximate throughout —
hexanal is printed as 100 vs 100.16, benzaldehyde 106 vs 106.12, citral 152 vs 152.23.)*

**⚠️ [!] Class labels differ between Table 1 and Table 2** for two compounds: **furaneol**
(Table 1 "Furans" → Table 2 "Enol") and **furfuryl mercaptan** (Table 1 "Furans" → Table 2
"Sulfur-containing compounds"). Cosmetic; the compound identity is unambiguous.

### 3c. **⭐ How far above food-relevant is the dose?** [Z]
| compound | dose here | typical food-relevant level | **× above food-relevant** |
|---|---|---|---|
| **hexanal** | **12 g L⁻¹ = 119.8 mmol L⁻¹** | 0.05–5 mg kg⁻¹ (plant-protein off-note range) | **2,400 – 240,000×** |
| **2-furfurylthiol (FFT)** | 12 g L⁻¹ = 105.1 mmol L⁻¹ | ~1–3 mg kg⁻¹ (roast coffee, high for a thiol) | **4,000 – 12,000×** |
| **furfural** | 12 g L⁻¹ = 125.0 mmol L⁻¹ | 1–100 mg kg⁻¹ (thermally processed foods) | **120 – 12,000×** |

**⇒ Every result in this paper is a maximum-potential, saturating-dose result — the same
verdict `anantharamkrishnan2020_extraction.md` §1 reached for the companion paper (A-4). The
authors state the reason explicitly. Register every row `dose: 12_ppth_saturating`.**

---

## 4. TABLE 2 — the reactivity verdict, all 47 compounds **[M]**

*Caption verbatim: "Table 2: Aroma compounds analyzed and their reactivity with BLG
(**24 h reaction time**)"*.

### 4a. **REACTIVE — 15 of 47 (32 %)**
| class | compound | **[Z] adduct observed, from the results text** |
|---|---|---|
| **α-Diketone** | **diacetyl** | **+86 Da, 1:1 addition** (NOT condensation; a hemiaminal-type Arg adduct, Fig. 3) |
| **Aldehydes** | **⭐ hexanal** | **"multiple adducts", multiple lysines** (spectrum in the SI, not in Figs 1–9) |
| | **benzaldehyde** | **a single Schiff-base imine only** |
| | **furfural** | **+78 Da (= 96 − 18), a single Schiff-base imine only** (Fig. 4) |
| **Enals** | **trans-2-hexenal** | adduct + **near-complete loss of protein signal → cross-linking**, at ≤ 10 min |
| | **cis-3-hexenal** | as above |
| | **trans,trans-2,4-heptadienal** | as above |
| | **citral** | **both +134 Da (Schiff base) and +152 Da (thia-Michael or aza-Michael)** (Fig. 5); **"much more slowly"** |
| **Enal-Acetals** | **citral diethyl acetal** | **+152 only** (hetera-Michael; no Schiff base) (Fig. 6); also cross-links |
| | **trans-2-hexenal dimethyl acetal** | adduct + cross-linking |
| **Sulfur** | **propanethiol** | **+76 (1:1) and +152 (2:1)** disulfide-exchange adducts (Fig. 7a,b) |
| | **thiophenol** | **+110 (1:1) only**, mass 18,388; **no 2:1 (+220) even at 48 h** (Fig. 7d) |
| | **⭐ furan-2-ylmethanethiol (= FFT)** | **+114 (1:1) and +228 (2:1)**, *"in close analogy to … n-PrSH"* (Fig. 7c) |
| | **dimethyl trisulfide** | **both +46 (BLG-CysSSMe) and +78 (BLG-CysSSSMe)**, competitive paths "x" and "y" (Fig. 8) |
| **Isothiocyanate** | **allyl isothiocyanate** | **+99, varying numbers of the 15 lysine NH₂** (Fig. 9) |

### 4b. **UNREACTIVE — 32 of 47 (68 %), within 24–48 h**
| class | compounds marked `no` |
|---|---|
| Hydrocarbon | p-cymene |
| Alcohols | (l)-menthol, geraniol, 2-pentanol, 2,3-butanediol |
| Phenols | 2,6-dimethylphenol, eugenol, 4-vinylphenol |
| Enol / Pyranone | **furaneol**, **maltol** |
| Acids | butyric acid, 2-hydroxybenzoic acid, m-toluic acid |
| Esters | methyl anthranilate, iso-amyl acetate, ethyl lactate, methyl salicylate |
| Lactones | δ-dodecalactone, γ-butyrolactone, γ-decalactone |
| Ketones (mono) | 2-heptanone, 2-nonanone, **cyclotene** |
| Bases | 2,5-dimethylpyrazine, 3-acetylpyridine, **allylamine** |
| Aldehyde | **vanillin** |
| Sulfur | 4-methyl-5-vinylthiazole (1 % in ethanol), 2-methylthiophene, dimethyl sulfone, dimethyl sulfide, **dimethyl disulfide ⚠** |

### 4c. **⚠️ [!] TABLE 2 CONTRADICTS THE RESULTS TEXT ON DIMETHYL DISULFIDE.**
Table 2: **dimethyl disulfide → `no`**.
Results text, p. 12–13, verbatim: *"However, **both** of the soft electrophiles **dimethyl
disulfide (MeSSMe)** and dimethyl trisulfide (MeSSSMe) **were observed to form covalent bonds
with BLG**. **The disulfide formed an adduct 46 Da greater than BLG**, corresponding to the net
addition of CH₂S. This is consistent with formation of BLG-CysSSMe by reaction of the BLG-CysH
with MeSSMe. … **In the reaction time studied, the extent of adduct formation for dimethyl
disulfide was less than that for dimethyl trisulfide.**"*
**⇒ The text says DMDS reacts and names the adduct mass; the table says it does not.** The most
likely reconciliation is that Table 2's binary is thresholded at 24 h and DMDS's adduct was
below that threshold (the text says its extent was lower, and *"Dimethyl disulfide is less
electrophilic and is expected to engage the BLG-CysH more slowly"*). **Record DMDS as
`reacts: yes, slowly` with the table's `no` flagged, and treat Table 2's binary as a
24 h-threshold verdict rather than an absolute one.** Consequence: **Table 2's 32 `no`s are
"no adduct above the detection threshold within 24 h", not "no reaction".** The authors say as
much: *"the reaction time allowed in this study was quite short, a maximum of 48 h. However,
compounds with less reactive functional groups might still engage over longer time periods."*

**A second Table 2 / text tension — this one resolved in the table's favour.** Table 2 marks
**allylamine `no`**, while the text reports *"allylamine, surprisingly, showed formation of a
simple, 1:1 addition complex (Figure 2) with BLG"*. The text then explains it away as an ESI
artefact: *"We speculate that this is a result of **ion-pairing of allylammonium ions instead
of protons** with the BLG in the ESI mass measurement."* **⇒ Table 2's `no` is correct;
Figure 2 is not a covalent adduct.** Any automated table-scrape of this paper will get
allylamine right and DMDS wrong.

---

## 5. ⭐ THE REACTIVITY LADDER — the paper's real deliverable for the repo

### 5a. **The 5-tier aldehyde ladder, with hexanal placed in it** [Z, from [M]+[Q] text]
Constructed from the results text's explicit comparisons, all at the same dose, protein,
temperature and protocol — which is what makes it a ladder rather than a list.

| tier | compounds | **evidence as printed** | timescale |
|---|---|---|---|
| **1 — fastest; cross-links** | **trans-2-hexenal · cis-3-hexenal · trans,trans-2,4-heptadienal** | *"all reacted quickly as judged by the **loss of protein signal** in the mass spectrum of each, **even of aliquots measured at only 10 min** after mixing. This is consistent with the known ability of these reactive enals to **cross-link whey proteins**."* | **< 10 min** |
| **1b — same, via hydrolysis** | trans-2-hexenal dimethyl acetal · citral diethyl acetal | *"The BLG also underwent **crosslinking to form aggregates** as it did when reacted with the corresponding aldehyde, as suggested with the **near-complete loss of the BLG signal**."* | ≤ 24 h |
| **2 — ⭐ multiple discrete adducts, no signal loss** | **HEXANAL** | *"**Hexanal showed multiple adducts, presumably engaging multiple lysines.**"* | ≤ 24 h |
| **3 — a single adduct only** | **benzaldehyde · furfural** | *"Both benzaldehyde and furfural (Figure 4), **less electrophilic aldehydes than hexanal**, were observed to give **only a single Schiff base imine**. **The reason for the lack of further reactivity compared with that seen for hexanal is unclear.**"* | ≤ 24 h |
| **4 — reacts, but slowly** | **citral** | *"Citral did so as well, but **much more slowly**, consistent with the **more hindered nature of the β-carbon** of the aldehyde (or its imine)."* | 6–24 h |
| **5 — no reaction** | **vanillin** | *"we saw **no evidence** for reaction between BLG and vanillin within 48 h, **even though there are reports suggesting that slow engagement of vanillin occurs**."* | > 48 h |

**⇒ THE ORDINAL RESULT, stated by the authors themselves:
`enals ≫ HEXANAL > benzaldehyde ≈ furfural > citral > vanillin (none)`,
with the explicit reason `hexanal is MORE electrophilic than benzaldehyde and furfural`.**

**Why this matters to the repo beyond hexanal:** **furfural** is a Maillard product the repo
tracks, and this ladder says furfural **does** have a covalent lysine sink — a single Schiff
base, weaker than hexanal's. Any matrix model that applies a covalent sink to hexanal and not
to furfural is internally inconsistent with the only side-by-side measurement that exists.

**Corroboration of the corpus's existing finding C4.** `anantharamkrishnan2020_extraction.md`
§C4 records *"alkenals cross-link protein"* as **"the best-supported irreversible mechanism in
the corpus"** on four independent lines. **This paper is a fifth, and the most direct**:
three different enals, all cross-linking within 10 min, versus a saturated alkanal at the same
dose that gives discrete adducts and leaves the protein signal intact. **And it supplies the
missing contrast term — the corpus previously had the enal side of the comparison (Meynier
2002's 6.88× t-2-hexenal suppression vs 1.39× hexanal) but no mechanistic statement of what
hexanal does instead. It forms multiple discrete lysine Schiff bases without cross-linking.**

### 5b. **The sulfur ladder — ⭐ FFT has a covalent protein sink** [M]+[Q]
| compound | adducts | **notes as printed** |
|---|---|---|
| **propanethiol** (aliphatic thiol) | **+76 and +152** | *"Both mono- and bis-adducts of BLG were observed. Moreover, **nearly none of the 2:1 (+152 Da) adduct was present at 6 h, but it can clearly be seen in the 24 h**."* |
| **⭐ furan-2-ylmethanethiol (FFT)** | **+114 and +228** | *"the most important aroma compound in roasted coffee aroma, **also reacted readily** with BLG to form covalent linkages. **In close analogy to the reactions with n-PrSH, both 1:1 and 2:1 adducts were observed.**"* |
| **thiophenol** (aromatic thiol) | **+110 only** | *"**no appreciable formation of a 2:1 adduct (i.e., +220 Da) was observed even at the 48 h** reaction time. Perhaps subtle differences in the **steric access of the bulkier PhSH to the second disulfide linkage**…"* |
| **dimethyl trisulfide** | **+46 and +78** | competitive paths x/y giving BLG-CysSSMe and BLG-CysSSSMe |
| **dimethyl disulfide** | **+46** | reacts, **less than DMTS** (⚠ Table 2 says `no`, §4c) |
| **dimethyl sulfide** | none | *"is not electrophilic, consequently, it did not react"* |
| 4-methyl-5-vinylthiazole, 2-methylthiophene, dimethyl sulfone | none | — |

**Mechanism as printed:** *"The thiol group is a **reactive (and soft) nucleophile capable of
cleaving one or both of the (soft) electrophilic disulfide linkages in BLG**"*, and the
sulfides *"formed covalent bonds with the **free cysteine group**"*. **⇒ Two distinct protein
sinks: thiols attack the 2 disulfides; di/trisulfides attack the 1 free Cys.**

**⭐ THE SULFUR-BRANCH CONSEQUENCE FOR THE REPO.** **2-Furfurylthiol — the compound the Kang
2026 SI ladder quantifies (`kang2026_SI_extraction.md` §4a) — is consumed by protein
disulfides, forming both 1:1 and 2:1 adducts, within 24 h at ambient temperature.** The
corpus's sulfur branch has had no protein-sink term at all. **This is direct, named evidence
that one is required for FFT in any protein-containing matrix** — and, by the aliphatic-thiol
analogy, for MFT too (2-methyl-3-furanthiol is not in this panel, but it is the same
furan-thiol class as FFT; **that is an inference, not a measurement — flag it**).
**⇒ FFT and MFT should not be modelled as chemically inert in a protein matrix.** Note also
that the sink here is **stoichiometric in disulfides (2 per BLG = 1.09 mmol L⁻¹)**, not
catalytic — so its capacity is bounded by the protein's disulfide inventory, which is a
directly computable matrix property.

### 5c. **The only two time-resolved observations in the paper** [M]
1. **Enals: adducts and protein-signal loss are already present at 10 min.** ⇒ **t½ ≲ 10 min**
   at 122 mmol L⁻¹ enal, ambient.
2. **Propanethiol's 2:1 adduct: "nearly none" at 6 h, "clearly seen" at 24 h.** ⇒ the *second*
   disulfide exchange has a characteristic time **between 6 and 24 h** at 158 mmol L⁻¹ thiol,
   ambient. The *first* exchange is faster (the 1:1 adduct is present at 6 h).

**These are the only two timescale brackets in the entire paper.** Everything else is a 24 h
endpoint binary.

### 5d. **Reversibility — asserted, never measured** [Q][!]
The paper's position, verbatim: *"These covalent interactions are particularly important to
understand as **these bonds are relatively stable. Once formed, they can be expected to be
cleaved only slowly**, making them difficult to be compensated for in formulation."* and
*"**covalent adducts are not likely to release the flavorant during consumption**"*.

**⚠️ [!] No reversibility experiment is reported. No adduct was subjected to dilution, pH
change, heating or time to test for reversion.** This matters because **the dominant hexanal
adduct here is a Schiff base (imine), which is thermodynamically reversible and hydrolyses
readily** — the paper's own control spectrum (§2a) shows lactose Schiff bases persisting, but
that is a different (Amadori-stabilisable) chemistry. **⇒ The claim that the hexanal–lysine
sink is irreversible is an assumption in this paper, not a finding. If the repo implements a
covalent hexanal sink, implementing it as irreversible is NOT supported by this source.**
By contrast the **thiol–disulfide exchange** adducts are genuinely exchange reactions and are
reversible in principle in the presence of other thiols — also untested.

---

## 6. ⭐ **[Z] BOUNDING THE COVALENT HEXANAL SINK — the derivation the repo needs, and its limits**

**The paper provides no rate. This section derives an order-of-magnitude bound from what it
does report. Every step is flagged; the result is assumption-dominated and is offered as a
sanity bracket, NOT as a fittable parameter.**

### 6a. Setting up
Schiff-base formation is second order: rate = k₂·[RCHO]·[Lys-NH₂].

| quantity | value | basis |
|---|---|---|
| [BLG] | **0.547 mmol L⁻¹** | 10 g L⁻¹ / 18,276 Da [Z] |
| [Lys] total | **8.21 mmol L⁻¹** | 15 Lys × [BLG] [M for 15] |
| [hexanal] | **119.8 mmol L⁻¹** | 12 g L⁻¹ / 100.16 [Z] |
| observation | **"multiple adducts"** at ≤ 24 h | [M] |

**The unmeasurable term is `f`, the fraction of the BLG population carrying those adducts.**
The spectra have no y-axis (§2c), so `f` cannot be read. Bracket it.

### 6b. The bracket
Hexanal consumed in 24 h ≥ 2 adducts × f × [BLG]. With [hexanal] and [Lys] both ≈ constant
(hexanal is in 15× excess over total lysine):

| assumed `f` | hexanal consumed | rate (M s⁻¹) | **k₂ (M⁻¹ s⁻¹)** |
|---|---|---|---|
| 0.05 (bare MS detectability) | 0.055 mM | 6.3 × 10⁻¹⁰ | **6.4 × 10⁻⁷** |
| 0.20 (central) | 0.219 mM | 2.5 × 10⁻⁹ | **2.6 × 10⁻⁶** |
| 1.00 (all BLG bears 2 adducts) | 1.09 mM | 1.3 × 10⁻⁸ | **1.3 × 10⁻⁵** |

> **⇒ k₂(hexanal + protein lysine, ambient, pH ≈ 6.7) ≈ 10⁻⁶ – 10⁻⁵ M⁻¹ s⁻¹, and this is a
> LOWER bound** — the adducts were observed *by* 24 h, not necessarily *at* 24 h; if they
> formed within the first hour, k₂ is up to 24× larger.

### 6c. **Scaling to a food matrix — the answer to the repo's question**
In a protein matrix the aldehyde is the trace species and the protein is in vast excess, so
hexanal loss is pseudo-first-order with **k_obs = k₂·[Lys]** — **independent of the hexanal
concentration**. For a **3 wt % protein** matrix at ~7 g Lys per 100 g protein:
[Lys]ₜₒₜₐₗ ≈ **16.4 mmol L⁻¹** (accessible fraction unknown and certainly lower).

| k₂ | k_obs (s⁻¹) | **t½ of hexanal at ambient** |
|---|---|---|
| 6.4 × 10⁻⁷ | 1.05 × 10⁻⁸ | **≈ 760 days** |
| 2.6 × 10⁻⁶ | 4.3 × 10⁻⁸ | **≈ 190 days** |
| 1.3 × 10⁻⁵ | 2.1 × 10⁻⁷ | **≈ 37 days** |

**⇒ AT AMBIENT TEMPERATURE the covalent lysine sink removes hexanal with a half-life of about
one month to two years. It is 3–5 orders of magnitude too slow to explain a 36× hexanal
over-prediction over a thermal process of minutes to hours.**

**Extrapolating to process temperature** — with a Schiff-base activation energy in the usual
50–80 kJ mol⁻¹ band, and going from ~20 °C to 140 °C:

| assumed Eₐ | rate acceleration 20 → 140 °C | **t½ at 140 °C (central k₂)** |
|---|---|---|
| 50 kJ mol⁻¹ | ≈ 330× | ≈ 14 days |
| 60 kJ mol⁻¹ | ≈ 1,280× | ≈ 3.6 hours *(sic — see note)* |
| **70 kJ mol⁻¹** | ≈ 4,300× | **≈ 65 minutes** |
| **80 kJ mol⁻¹** | ≈ 13,900× | **≈ 20 minutes** |

*(Note: the "14 days" row is computed on the 190-day central; the hour-scale rows follow from
the same central at the higher Eₐ.)*

> **⇒ THE CONDITIONAL CONCLUSION.** The covalent-lysine sink can plausibly consume a large
> fraction of hexanal **only** at process temperatures **and only** if its activation energy is
> at the top of the plausible range (**Eₐ ≳ 70 kJ mol⁻¹**). At Eₐ ≲ 60 kJ mol⁻¹ the channel is
> negligible on any food-process timescale. **This paper cannot distinguish those cases — it
> reports one temperature and never states what it is (§2b [!]).** Its companion
> (`anantharamkrishnan2020_extraction.md` §C3) establishes only that adduct formation increases
> from 4 °C to 20 °C, **directionally, with no number and no upper temperature.**

**⚠️ THE HONEST SUMMARY: this paper turns the covalent sink from "a mechanism that might not
even apply to hexanal" into "a mechanism that demonstrably applies to hexanal, of unknown
size". That is a real advance on `anantharamkrishnan2020_extraction.md` §6, which had to say
hexanal was absent from the panel. It is not enough to close the over-prediction.**

**Assumptions this derivation rests on, all of which could each move the answer by ≥ 1 order
of magnitude:** (1) "multiple" = 2 adducts; (2) `f` is unmeasurable and bracketed by fiat;
(3) all 15 lysines are equally reactive and accessible — false, and the control spectrum shows
some are pre-blocked by lactose (§2a [!]); (4) simple second-order kinetics with no
carbinolamine pre-equilibrium; (5) irreversibility, which this paper asserts but does not test
(§5d); (6) matrix lysine accessibility taken as 100 %; (7) the reaction pH here is unknown and
Schiff-base rates are strongly pH-dependent — the companion paper shows carbonyl–lysine adduct
formation is **abolished at pH 3** (its §C1). **Register the whole of §6 as
`[Z] order_of_magnitude, assumption_dominated, NOT_FITTABLE`.**

---

## 7. CONSOLIDATED PARAMETER TABLE

**Common conditions:** **1 wt % BLG variant B (18,276 Da, > 90 % pure, 15 Lys / 1 free Cys /
2 disulfides / 3 Arg)** in doubly distilled water, **flavour at 12 ppth added individually**,
closed vessel, **ambient temperature (value not stated)**, **pH not stated**, sampled
**10 min – 24 h (occasionally 48 h)**, **n not stated**, UPLC-ESI-qTOF, **± 1 Da**.

| # | parameter | value | units | prov | source |
|---|---|---|---|---|---|
| **— SYSTEM —** |
| 1 | BLG variant B parent mass | **18,276** | Da | **M** | Fig. 1 |
| 2 | BLG A↔B mass difference | **87** | Da | **M** | p. 4 |
| 3 | BLG purity | **> 90** | % | **M** | Methods |
| 4 | BLG isolate pH | **6.7** | — | **M** | Methods |
| 5 | Lys / free Cys / disulfides / Arg per BLG | **15 / 1 / 2 / 3** | count | **M** | Methods |
| 6 | Protein concentration | **1** | wt % (= **0.547 mmol L⁻¹**) | **M**/**Z** | Methods |
| 7 | **Total lysine** | **8.21** | mmol L⁻¹ | **Z** | §6a |
| 8 | **Total disulfide (the thiol sink capacity)** | **1.09** | mmol L⁻¹ | **Z** | §6a |
| 9 | Flavour dose | **12** | ppth (= 12 g L⁻¹ = 1.2 % w/w) | **M** | Methods |
| 10 | **[hexanal]** | **119.8** | mmol L⁻¹ | **Z** | §3 |
| 11 | **[FFT]** | **105.1** | mmol L⁻¹ | **Z** | §3 |
| 12 | Molar ratio flavour : BLG, **as printed** | "ca. **50–250**" | — | **M** | Methods |
| 13 | Molar ratio flavour : BLG, **recomputed** ⚠ | **97–385** (hexanal **219**) | — | **Z** | §2b [!] |
| 14 | Mass accuracy | **± 1** | Da | **M** | Methods |
| **— PANEL —** |
| 15 | Compounds tested | **47** in **13** classes | count | **M**/**Z** | Table 1 |
| 16 | **Reactive (24 h)** | **15** (**32 %**) | count | **M**/**Z** | Table 2 |
| 17 | Unreactive (24 h) | **32** (68 %) | count | **M**/**Z** | Table 2 |
| **— ⭐ HEXANAL —** |
| 18 | **Hexanal forms covalent adducts with BLG** | **YES** | binary | **M** | Table 2 |
| 19 | **Hexanal adduct multiplicity** | **"multiple", multiple lysines** | — | **Q** | p. 10 |
| 20 | Hexanal adduct type | **Schiff base (imine)** | — | **Q** | p. 10 |
| 21 | Hexanal electrophilicity rank | **> benzaldehyde, > furfural** (authors' own comparison) | ordinal | **Q** | p. 10 |
| 22 | Hexanal cross-links protein? | **NO** — protein signal retained (contrast: enals) | binary | **Z** | §5a |
| 23 | **[Z] k₂ (hexanal + Lys), lower bound** | **10⁻⁶ – 10⁻⁵** | M⁻¹ s⁻¹ | **Z⚠** | §6b |
| 24 | **[Z] t½ hexanal, 3 wt % protein, ambient** | **37 – 760 days** | days | **Z⚠** | §6c |
| **— OTHER ALDEHYDES —** |
| 25 | Benzaldehyde | **single Schiff base only** | — | **M** | Table 2, p. 10 |
| 26 | **Furfural** | **+78 Da, single Schiff base only** | Da | **M** | Fig. 4 |
| 27 | **Vanillin** | **NO reaction within 48 h** | binary | **M** | Table 2, p. 10 |
| 28 | trans-2-hexenal, cis-3-hexenal, tt-2,4-heptadienal | **react + cross-link at ≤ 10 min** | — | **M**/**Q** | p. 11 |
| 29 | **Enal reaction timescale** | **< 10** | min | **M** | p. 11 |
| 30 | Citral | **+134 (Schiff) and +152 (hetera-Michael)**, "much more slowly" | Da | **M** | Fig. 5 |
| 31 | Citral diethyl acetal | **+152 only** | Da | **M** | Fig. 6 |
| **— ⭐ SULFUR —** |
| 32 | **FFT (furan-2-ylmethanethiol)** | **+114 (1:1) and +228 (2:1)** | Da | **M** | Fig. 7c |
| 33 | Propanethiol | **+76 (1:1) and +152 (2:1)** | Da | **M** | Fig. 7a,b |
| 34 | **Propanethiol 2:1 adduct appearance** | **between 6 h and 24 h** | h | **M** | p. 12 |
| 35 | Thiophenol | **+110 (1:1) only**; no 2:1 even at **48 h** | Da | **M** | Fig. 7d |
| 36 | Dimethyl trisulfide | **+46 and +78** | Da | **M** | Fig. 8 |
| 37 | Dimethyl disulfide | **+46**; extent **< DMTS** ⚠ Table 2 says `no` | Da | **M**/**[!]** | p. 13, §4c |
| 38 | Dimethyl sulfide | **no reaction** | binary | **M** | Table 2 |
| 39 | Thiol sink mechanism | **cleaves one or both BLG disulfides** | — | **Q** | p. 12 |
| 40 | Sulfide sink mechanism | **binds the free Cys** | — | **Q** | Abstract |
| **— OTHER —** |
| 41 | **Diacetyl** | **+86, 1:1 ADDITION (not condensation)**, likely on **Arg guanidine** | Da | **M** | Fig. 3 |
| 42 | Allyl isothiocyanate | **+99, varying numbers of the 15 lysines** | Da | **M**/**C** | Fig. 9, ref. 15 |
| 43 | Allylamine "adduct" | **ESI ion-pairing artefact, NOT covalent** | — | **Q** | Fig. 2, p. 10 |
| 44 | Monoketones (2-heptanone, 2-nonanone, cyclotene) | **no reaction** | binary | **M** | Table 2 |
| 45 | Alcohols, acids, esters, lactones, phenols, hydrocarbons, pyrazine, pyridine, furaneol, maltol | **no reaction** | binary | **M** | Table 2 |
| **— CONTROL-SPECTRUM CONTAMINATION —** |
| 46 | BLG + 1 lactose (Schiff) | **18,600** (+324) | Da | **M** | Fig. 1 |
| 47 | BLG + 2 lactose | **18,924** (+648) | Da | **M** | Fig. 1 |
| 48 | Lactose Schiff-base mass increment | **324** (= 342 − 18) | Da | **M** | Fig. 1 |
| 49 | Unassigned control peaks | 18,292 · 18,436 · 18,632 | Da | **M** | Fig. 1 |

### **⚠️ WHAT IS NOT IN THIS PAPER — the complete negative inventory**
**No rate constant · no half-life · no yield · no % conversion · no equilibrium constant · no
calibration curve · no LOD · no reversibility experiment · no reaction temperature value · no
reaction pH value · no replicate count (n) · no error bar of any kind · no y-axis on any of the
9 figures · no hexanal spectrum in the main article (it is in SI Figs S1–S7, which are NOT in
this PDF) · no plant protein (BLG only) · no temperature, pH or a_w series (that is the
companion paper) · no food matrix.**

**⚠️ The hexanal spectrum itself is not in this PDF.** SI note, verbatim: *"Deconvoluted ESI
mass spectra (**Figure S1–S7**) of those flavor compounds that showed reactivity but where the
data were not provided in any one of the Figure 1–9 in the manuscript."* **Hexanal is reactive
and has no main-text figure ⇒ its spectrum is one of S1–S7.** The number of hexanal adducts —
the one thing that would sharpen §6's bound from "multiple ≥ 2" to a real integer — **is in an
SI we do not have.** See §10.

---

## 8. CITATION INTEGRITY SPOT-CHECK

`anantharamkrishnan2020_extraction.md` §A-8 flagged the **companion** paper for citing
Gremli's soy-protein paper with *"the wrong journal, the wrong year, the wrong page range and
a fabricated author list"* (*Chem. Res. Toxicol.* 2010, 51, 1050−1059, multi-author).

**⇒ THIS PAPER CITES THE SAME WORK CORRECTLY.** Ref. 7, verbatim:
*"Gremli, H. A. **Interaction of Flavor Compounds with Soy Protein**. ***J. Am. Oil Chem. Soc.***
**1974**, ***51*** (1), **95–97**. https://doi.org/10.1007/BF02542100"* — **correct journal,
correct year, correct volume, correct pages, correct single author.**

**⇒ The corruption is confined to the LATER paper (0c06752). This earlier paper's reference
list is clean at the point where its sibling is not.** Spot-checks of refs 15, 20, 28 and 30
against their DOIs also come back internally coherent (author/journal/year/volume/pages all
mutually consistent, DOI prefixes matching the named publishers). **No citation defect found
in this paper.** *(This is a useful negative result for the corpus's 30–45 % contamination
audit: it localises the Gremli defect to one document rather than to the author group.)*

**Ref. 28 is the important one and it is cited correctly:**
*"Meynier, A.; Rampon, V.; Dalgalarrondo, M.; Genot, C. **Hexanal and T-2-Hexenal Form Covalent
Bonds with Whey Proteins and Sodium Caseinate in Aqueous Solution**. Int. Dairy J. **2004**,
**14** (8), 681–690. https://doi.org/10.1016/j.idairyj.2004.01.003"*

---

## 9. **VERDICT ON THE THREE BRIEF QUESTIONS**

| question | **verdict** |
|---|---|
| **(a) Is any n-alkanal — hexanal especially — among the compounds, and does it form covalent adducts under food-relevant conditions?** | **YES to the panel, YES to adduct formation, NO to food-relevant conditions.** Hexanal is the **only** saturated n-alkanal in the 47-compound panel. It forms **multiple** Schiff-base adducts on multiple lysines within 24 h. But at **12 ppth = 2,400–240,000× a food-relevant loading**, at an **unstated ambient temperature**, at an **unstated pH**, with **n unstated**. |
| **(b) Any rates, extents, or reversibility usable to bound the covalent sink?** | **NO measured ones. Zero rate constants, zero extents, zero reversibility experiments.** What is usable: an **ordinal 5-tier ladder** with hexanal placed by the authors' own comparisons (§5a); **two timescale brackets** (enals < 10 min; second thiol exchange 6–24 h, §5c); and **this wave's derived bound k₂ ≈ 10⁻⁶–10⁻⁵ M⁻¹ s⁻¹ ⇒ t½ ≈ 37–760 days at ambient in 3 wt % protein** (§6), which is **assumption-dominated and not fittable**. |
| **(c) Doses — food-relevant or saturating?** | **SATURATING BY DESIGN**, and the authors state the reason in the Methods: *"much higher than that used in the industry … for the ease and accuracy of measuring the adducts."* |

### **⇒ DOES THIS CLOSE THE REPO'S HEXANAL OVER-PREDICTION? NO — but it changes its status.**
**Before:** the covalent sink was a *speculative analogue* — `anantharamkrishnan2020_extraction.md`
§6 had to report that hexanal, and every n-alkanal, was **absent** from that paper's panel, so
transferring a benzaldehyde result to hexanal was *"a leap"*.
**After:** the covalent sink is a **demonstrated mechanism for hexanal specifically, on a named
protein, with a named adduct type, ranked against its own homologues** — and, at the same time,
**demonstrably too slow at ambient temperature, by 3–5 orders of magnitude, to be the
explanation for a 36× discrepancy in a thermal process.** The channel survives only in the
high-Eₐ / high-temperature corner (§6c), which this paper cannot test.

---

## 10. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR**

*(Advisory only. Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration
was read for state or written.)*

### 10a. **NO FIT ROWS ARE PROPOSED FROM THIS PAPER.**
**There is no fittable number in it.** Every quantitative value is a molecular mass in daltons
(a chemical identity, not a measurement of extent), a dose, or a count. **Any "fit" against
this paper would be a fit against a binary or an ordinal.** Recommending zero fit rows is the
correct outcome, and it is the same outcome its companion reached.

### 10b. Recommended **HOLD-OUT / GATE** rows — all ordinal or binary
| # | proposed hold-out | statement | why it is a good gate |
|---|---|---|---|
| **G-1** | ⭐ **The aldehyde reactivity ordering** | `enals ≫ hexanal > benzaldehyde ≈ furfural > citral > vanillin(none)` | **Six species, one protein, one dose, one protocol, one temperature.** Ordinal gates are cheap to evaluate and hard to satisfy by accident. **And it is stated by the authors, not inferred by us.** A matrix model that ranks any pair the other way is falsified. |
| **G-2** | **Enals cross-link; the saturated alkanal does not** | t-2-hexenal / c-3-hexenal / tt-2,4-heptadienal all destroy the protein signal at ≤ 10 min; hexanal gives discrete adducts with the signal intact | **A fifth independent line on corpus finding C4**, and the first that supplies the *contrast* term. This is now the corpus's best-supported irreversible mechanism by a wide margin. |
| **G-3** | ⭐ **FFT is consumed by protein disulfides (1:1 and 2:1)** | binary, named compound, named mechanism | **The sulfur branch currently has NO protein-sink term.** This is direct evidence one is required. Pair it with `kang2026_SI_extraction.md`'s FFT ladder: the Kang system is protein-free, so a repo model calibrated on Kang and applied to a protein matrix must lose FFT relative to Kang. |
| **G-4** | **Furfural also has a covalent lysine sink (single Schiff base, weaker than hexanal's)** | binary + ordinal | Forces internal consistency: **a model cannot give hexanal a covalent sink and furfural none.** |
| **G-5** | **32 of 47 compounds show NO adduct at 24 h** — alcohols, acids, esters, lactones, phenols, monoketones, pyrazines, pyridines, hydrocarbons, furaneol, maltol, DMS, sulfone, thiazole, thiophene | binary, 32-way | **A large negative gate.** A model that applies a generic "protein binding" loss term across all volatile classes is falsified by 32 counter-examples at a saturating dose. **This is the most under-appreciated result in the paper.** |
| **G-6** | **Vanillin: no reaction in 48 h despite four literature reports to the contrary** (refs 24–27) | binary | A **contrarian** gate — it disagrees with the prior literature the paper itself cites, so it discriminates between a model built on those reports and one built on this measurement. |
| **G-7** | **Diacetyl forms a 1:1 ADDITION adduct (+86), not a condensation (+68)** | mass | A **mechanism-discriminating** mass. +68 vs +86 distinguishes Schiff base from hemiaminal/Arg-adduct unambiguously at ±1 Da. |

### 10c. **DO NOT USE**
| # | item | reason |
|---|---|---|
| **X-1** | **Any peak-height or intensity ratio from Figs 1–9** | **The figures have no y-axis** (§2c, verified on 600 dpi raster). MaxEnt intensities are not molar-proportional and no spectrum is normalised. |
| **X-2** | **Table 2's `no` for dimethyl disulfide** | Contradicted by the paper's own results text, which names the +46 Da adduct (§4c [!]). |
| **X-3** | **Table 2's binary read as "does not react"** | It is **"no adduct above threshold within 24 h"**. The authors say slower compounds *"might still engage over longer time periods"* (§4c). |
| **X-4** | **§6's k₂ or t½ as a fitted or fittable parameter** | Assumption-dominated; `f` is unmeasurable from this source; 7 named assumptions each worth ≥ 1 order of magnitude (§6c). **Bracket only.** |
| **X-5** | **Any assumption that the hexanal adduct is irreversible** | Asserted, never tested (§5d [!]). The adduct is an imine. |
| **X-6** | **The printed molar ratio "ca. 50–250"** | Inconsistent with the printed dose and protein concentration; true range ≈ 97–385 (§2b [!]). |
| **X-7** | **Cross-comparing adduct *intensities* with the companion paper (0c06752)** | **Different BLG variant** (B here, A there), and the companion states its own signals are weaker for that reason (§0 [!]). Ordinal results transfer; magnitudes do not. |
| **X-8** | **Treating the 15 lysines as all free** | The control spectrum shows mono- and di-lactosylated BLG as major peaks (§2a [!]); the free-lysine count is < 15 and unmeasured. |

### 10d. Registry hygiene proposed for every row sourced here
```
source: anantharamkrishnan2020b     doi: 10.1021/acs.jafc.0c01925
protein: BLG_variant_B              parent_mass_Da: 18276
dose: 12_ppth_saturating            temperature: ambient_VALUE_NOT_STATED
pH: NOT_STATED                      n: NOT_STATED
observation_window: 10min_to_24h    result_type: binary_or_ordinal_ONLY
```

---

## 11. HIGHEST-VALUE FOLLOW-UP RETRIEVALS, RANKED

1. **⭐⭐⭐ Meynier, A.; Rampon, V.; Dalgalarrondo, M.; Genot, C. (2004).** *"Hexanal and
   t-2-Hexenal Form Covalent Bonds with Whey Proteins and Sodium Caseinate in Aqueous
   Solution."* **Int. Dairy J. 14 (8), 681–690. DOI 10.1016/j.idairyj.2004.01.003** — **ref. 28
   of this paper, cited correctly (§8).** The title states outright that it measures what this
   paper only demonstrates qualitatively, **for hexanal specifically**, in **aqueous solution**,
   on **both whey protein and sodium caseinate**, and by its journal and vintage will report
   **extents and time courses**, not mass spectra. **This is the single retrieval most likely to
   close the repo's hexanal over-prediction, and it should outrank everything else in the
   queue.** *(Distinct from `Meynier2002_extraction.md` already in the corpus — different year,
   different journal, different question. Same first author.)*
2. **The SI of this paper (Figs S1–S7)** — **contains the hexanal spectrum**, which would turn
   "multiple adducts" into an integer and sharpen §6's bound by up to an order of magnitude.
   Free at pubs.acs.org.
3. **Anantharamkrishnan, V.; Reineccius, G. A. (2020),** *"Method To Characterize and Monitor
   Covalent Interactions of Flavor Compounds with β-Lactoglobulin Using Mass Spectrometry and
   Proteomics"*, **DOI 10.1021/acs.jafc.9b07978** (ref. 15) — **already nominated by
   `anantharamkrishnan2020_extraction.md` §8.2.** It is the **site-specific proteomics** paper
   (*which* lysines) and it is **the variant-B study**, so it is the only route to the
   free-lysine count that X-8 says is missing.
4. **Keppler et al. (2014)**, ref. 30, *J. Biomol. Struct. Dyn.* 32(7), 1103–1117 — AITC + BLG
   by *"fluorescence quenching, **equilibrium measurement**, and mass spectrometry"*. **The word
   "equilibrium measurement" in the title means it may carry the binding constants this whole
   author group omits** — including, potentially, a reversibility handle (§5d).
5. **Meltretter et al. (2008)**, ref. 16, *JAFC* 56(13), 5165–5171 — *"Identification and
   Site-Specific **Relative Quantification** of β-Lactoglobulin Modifications in **Heated Milk**"*.
   **Heated** and **quantified** — the two things missing here, and directly relevant to §6c's
   high-temperature corner.

---

*End of `anantharamkrishnan2020b_extraction.md`. Companion:
`anantharamkrishnan2020_extraction.md` (the later pH/T/a_w paper, variant A).
No repo file outside `data/lit/extraction_dossiers/` was created or modified by this wave.*
