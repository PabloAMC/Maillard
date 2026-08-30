# Zamora, Alcón & Hidalgo 2013 — COMPLETE TRANSCRIPTION
### *"Strecker-Type Degradation of Phenylalanine Initiated by 4-Oxo-2-alkenals in Comparison to That Initiated by 2,4-Alkadienals, 4,5-Epoxy-2-alkenals, or 4-Hydroxy-2-nonenal"*
### ⭐⭐ **THE Eₐ LADDER ACROSS CARBONYL CLASSES — eight measured activation energies in one table, one lab, one method.** Wave **K6b**, 2026-08-29.

**Full extraction of every number in `data/articles/zamora2013.pdf`.**
Read-only wave. No repo file outside `data/lit/extraction_dossiers/` created or modified.

**Read first:** `research_round4_nulls.md` §A.1 (which found this paper by abstract and
recorded `Ea = 55–64 · 58–67 · 28–38 [ABS]` — **§3 below corrects and completes that record**).

**Provenance:** **[M]** measured · **[C]** cited · **[Q]** qualitative · **[Z]** derived here · **[!]** flag.

---

## 0. IDENTITY — **MATCHES EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Rosario Zamora, Esmeralda Alcón, Francisco J. Hidalgo\*** |
| Venue | ***J. Agric. Food Chem.* 2013, 61 (43), 10231−10237** |
| DOI | **dx.doi.org/10.1021/jf305007y** ✔ **exactly the expected DOI** |
| Affiliation | **Instituto de la Grasa, CSIC, Avenida Padre García Tejero 4, 41012 Seville, Spain** |
| Special issue | ***ISMR11 — 100 Years of the Maillard Reaction*** |
| Dates | Received Nov 22 2012 · Revised Jan 15 2013 · Accepted/Published Jan 29 2013 |
| Corresponding | fhidalgo@ig.csic.es · +34 954 611 550 |
| Funding | EU FEDER + Plan Nacional de I+D, MINECO Spain, **AGL2009-07638, AGL2012-35627**; José L. Navarro thanked |
| PDF character | 7 pages, clean text layer, **Tables 1–2**, Figures 1–7 (Fig. 4 = three Arrhenius plots) |

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER

**Table 2 is the object the corpus came for: eight activation energies for carbonyl–amine
chemistry, all measured by the same lab with the same method (linear zero-order fit of product
formation vs time at 4–7 temperatures, then an Arrhenius plot), spanning four structural
classes of lipid-derived carbonyl.** It is the **densest single-table Eₐ resource in the
corpus** and it is what turns `research_round4_nulls.md` §A from a null into a ladder.

**Three things must be said immediately about what these numbers ARE:**
1. **They are Eₐ of a DOWNSTREAM PRODUCT (phenylacetaldehyde), not of adduct formation.** The
   reaction is `Phe + carbonyl → imine → decarboxylation → hydrolysis → phenylacetaldehyde`,
   and the fitted rate is the appearance of the last species. The Eₐ is that of the
   rate-determining step of the whole cascade, most plausibly the **decarboxylation**, not the
   carbonyl–amine condensation.
2. **The amine is a free amino acid, not a protein.** No protein appears in this paper.
3. **⭐ The temperature range is 120–190 °C and the pH is 3.** This is a **process**
   measurement, in the range the repo's Maillard path actually runs — a genuine strength,
   and a genuine mismatch with the ambient-storage covalent-sink question.

**⇒ These eight numbers bound the class-level Eₐ landscape for lipid-derived carbonyl–amine
chemistry. They do NOT supply the aldehyde–protein adduct Eₐ; `shepelev2024_extraction.md`
does that.** Both are needed and they answer different questions (`k6b_adduct_kinetics_synthesis.md` §2).

---

## 2. SYSTEM AND METHOD — complete transcription **[M]**

| variable | value as printed |
|---|---|
| **Amine** | **L-phenylalanine, 10 µmol** |
| **Carbonyl** | **0–10 µmol** of the lipid derivative; Table 2 uses **equimolecular amounts** |
| **Solvent** | **0.5 mL of buffer** |
| **Vessel** | **Schott Duran test tubes, 16 × 1.5 cm**, closed |
| **Heating** | **heater block, 120−190 °C, 0−1 h** |
| **Atmosphere** | **air**, unless otherwise indicated (one series under **nitrogen**) |
| **Buffers** | **0.3 M sodium citrate, pH 2.15−5.0** · **0.3 M sodium phosphate, pH 6.0−8.0** · **0.3 M sodium borate, pH 9.0−10** |
| **pH used for all kinetics** | **pH 3** — *"previous studies … were carried out at pH 3, this pH was also employed for the rest of the experiments"* |
| Workup | cooled, + **1 mL acetonitrile** + **50 µL internal standard (54.8 mg methyl heptanoate in 25 mL methanol)** |
| **Internal standard** | **methyl heptanoate** |
| GC-MS | HP 6890 GC Plus + Agilent 5973 MSD; **HP5-MS 30 m × 0.25 mm i.d., 0.25 µm**; **1 µL pulsed splitless**; He 1 mL/min constant flow; injector **250 °C**; oven **40 °C (1 min) → 240 °C at 5 °C/min → 300 °C at 10 °C/min**; transfer line **280 °C**; **EI 70 eV**; ion source **230 °C**; **m/z 28−550** |
| Quantitation | **eight-point standard curve** of phenylacetaldehyde prepared in the same 1.55 mL GC-MS solution; **r = 0.999, p < 0.0001**; **coefficients of variation < 10 %** |
| Statistics | *"mean values of **several** independent determinations"* ⚠️ **n never given [!]**; ANOVA + **Tukey test**; Origin v. 7.0; **p < 0.05** |

### 2a. **The four carbonyl classes and how each was made [M]**
| compound | preparation |
|---|---|
| **4-oxo-2-nonenal** | from **2-pentylfuran** with **N-bromosuccinimide** (refs. 13, 14) |
| **4-oxo-2-hexenal** | from **2-ethylfuran** with NBS |
| **4,5-epoxy-2-heptenal** | epoxidation of **2,4-heptadienal** with **3-chloroperoxybenzoic acid** (refs. 15, 16) |
| **4,5-epoxy-2-decenal** | epoxidation of **2,4-decadienal** with mCPBA |
| **4-hydroxy-2-nonenal** | synthesized per **Gardner et al. 1992** (ref. 17) |
| 2,4-heptadienal, 2,4-decadienal | purchased (Aldrich/Sigma/Fluka/Merck, analytical grade) |

### 2b. **The rate law — verbatim, and note that it is ZERO order [M]**
> *"the amount of phenylacetaldehyde **increased linearly** (r > 0.972, p < 0.003) as a function
> of reaction time between 120 and 180 °C … Reaction rates at the different assayed
> temperatures were calculated by using the equation
> **[phenylacetaldehyde] = [phenylacetaldehyde]₀ + kt**
> where [phenylacetaldehyde]₀ represents the intercept, k is the rate constant, and t is the
> reaction time. These rate constants were used in an Arrhenius plot."*

⚠️ **[!] THE FIT IS ZERO-ORDER IN A SYSTEM WHERE BOTH REACTANTS ARE BEING CONSUMED.**
`k` has units of concentration/time and is **not** a rate constant in the Arrhenius sense
unless the reactant concentrations are effectively constant over the fitted window — which the
authors ensure only by staying in the initial-rate region (linear, r > 0.96–0.99). The
resulting Eₐ is an **apparent** Eₐ of an initial rate. **This is the same estimator class as
`hidalgo1993_extraction.md` (zero-order in colour) and `hidalgo2010_extraction.md`
(zero-order in acrylamide loss); the whole Instituto de la Grasa Eₐ corpus shares it.
Register `estimator: zero_order_initial_rate_Arrhenius` on every row.**

### 2c. **The temperature sets actually used, per figure [M]** — they are NOT all the same
| figure | system | **temperatures (°C)** | n points |
|---|---|---|---|
| **Fig. 3** | oxoalkenals (4-oxo-2-hexenal, 4-oxo-2-nonenal) | **180, 160, 140, 120** | **4** |
| **Fig. 5** | epoxyalkenals + 4-HNE | **190, 180, 170, 160, 150, 140, 120** | **7** |
| **Fig. 6** | 2,4-alkadienals (air and nitrogen) | **180, 160, 150, 140, 130, 120** | **6** |

**⇒ The oxoalkenal Eₐ values rest on 4 points; the epoxyalkenal/HNE values on 7; the alkadienal
values on 6. Weight accordingly.** *(Text says "linearly … between 120 and 180 °C" for Figs. 3
and 6 and "between 120 and 190 °C" for Fig. 5.)*

---

## 3. ⭐⭐ TABLE 2 — THE Eₐ LADDER — TRANSCRIBED IN FULL **[M]**

*Caption verbatim: "**Table 2.** Activation Energies (Eₐ) and Amount of Phenylacetaldehyde
Produced in the Degradation of Phenylalanine by Lipid-Derived Reactive Carbonyls."
Footnote a: "Equimolecular amounts of phenylalanine and the oxidized lipid were heated at the
indicated temperature for the specified reaction time. **Means in the same column with the same
letters are not significantly different.**" Footnote b: "ND, not determined."*

| carbonyl compound | **Eₐ (kJ/mol)** | **1 h at 180 °C** | **2 h at 80 °C** | **24 h at 37 °C** |
|---|---:|---:|---:|---:|
| **none** (control) | **ND** | 8.5 ± 2.0 c | 4.5 ± 0.4 c | 4.8 ± 0.1 c |
| **4-oxo-2-hexenal** | **55.3** | 140.9 ± 8.5 d | 13.6 ± 1.1 de | 9.5 ± 3.3 cd |
| **4-oxo-2-nonenal** | **63.5** | **168.0 ± 9.3 e** | **35.5 ± 3.0 f** | 14.5 ± 2.7 dfg |
| **4,5-epoxy-2-heptenal** | **57.9** | 99.8 ± 6.3 f | 10.9 ± 1.1 d | 10.6 ± 1.6 def |
| **4,5-epoxy-2-decenal** | **58.9** | 92.4 ± 6.1 f | 16.6 ± 2.6 e | 15.8 ± 1.6 eg |
| **4-hydroxy-2-nonenal** | **67.1** | 21.1 ± 0.9 c | 5.1 ± 0.4 c | 5.6 ± 0.2 cf |
| **2,4-heptadienal** | **37.5** | 88.6 ± 7.3 f | 18.1 ± 0.5 e | 18.0 ± 2.2 g |
| **2,4-decadienal** | **27.6** | 133.5 ± 2.5 d | 29.3 ± 2.3 g | **22.9 ± 1.6 h** |
| **2,4-decadienal (nitrogen)** | **78.0** | ND | ND | ND |

*Phenylacetaldehyde units: **µmol per mmol of phenylalanine**. All ± values are as printed
(SD, n unstated).*

### 3a. ⭐ **CORRECTIONS TO THE ROUND-4 ABSTRACT-ONLY RECORD**
`research_round4_nulls.md` §A.1 recorded, from the abstract: *"55–64 (4-oxo-2-alkenals) ·
58–67 (4,5-epoxy-2-alkenals + 4-HNE) · 28–38 (2,4-alkadienals)"*. The full-text values are:

| the round-4 [ABS] record | **the printed [M] value** | verdict |
|---|---|---|
| oxoalkenals 55–64 | **55.3 and 63.5** | ✔ **correct** |
| epoxyalkenals + HNE 58–67 | **57.9, 58.9, 67.1** | ⚠️ **range should be 57.9–67.1** (the abstract rounded 57.9 up to 58) |
| alkadienals 28–38 | **27.6 and 37.5** | ⚠️ **range should be 27.6–37.5** |
| — | ⭐ **2,4-decadienal under NITROGEN = 78.0 kJ/mol** | ⭐ **ENTIRELY ABSENT from the abstract — a ninth value, and the largest in the table** |

**⇒ The abstract-only record was directionally right and missed the single most informative
row. Register the corrected ladder from this dossier, not from `research_round4_nulls.md`.**

### 3b. ⭐⭐ **THE OXYGEN RESULT — 2.8× ON THE ACTIVATION ENERGY [M]**
> *"in the absence of air, the Eₐ for the reaction with decadienal was **2.8 times higher** than
> when the reaction was carried out in the presence of air."*
**27.6 (air) → 78.0 (N₂). [Z] Ratio = 2.83 ✔ the arithmetic checks.**

The authors' interpretation [Q]:
> *"These results suggest that, **in the presence of oxygen, alkadienals should be oxidized and
> the product(s) of this oxidation should be the responsible for the degradation of
> phenylalanine** … the conversion of alkadienals into epoxyalkenals was proposed. The results
> obtained in the present study suggest that **conversion of alkadienals into other more
> reactive intermediates should also be produced because the Eₐ determined for phenylalanine
> degradation in the presence of alkadienals was LOWER than the Eₐ determined … in the presence
> of epoxyalkenals.**"*

**⇒ ⭐ A MEASURED, 2.8×, ATMOSPHERE-DEPENDENT ACTIVATION ENERGY FOR THE SAME NOMINAL REACTANT
PAIR. This is the strongest single caution in the Eₐ ladder: an Eₐ is a property of the
mechanism actually operating, and the mechanism here is set by an unlisted variable (O₂).
Any repo row that carries an Eₐ for an alkadienal must also carry the atmosphere.
Register `2,4-decadienal + Phe: Ea = 27.6 (air) | 78.0 (N2)`, never a single value.**
*(Note this is the same lab's finding as `hidalgo2010_extraction.md` §, where oxygen changes
the acrylamide-mercaptan pathway entirely — a consistent methodological signature.)*

### 3c. ⭐⭐ **THE Eₐ↔YIELD DECOUPLING — the paper's own most important caveat [M]**
> *"this **difference in Eₐ only correlated with the amount of phenylacetaldehyde produced at
> 37 °C**. At higher temperatures, 4-oxo-2-nonenal was the lipid-derived carbonyl compound that
> produced the highest amount of the Strecker aldehyde."*

Quantified in the Discussion, verbatim [M]:
> *"when phenylalanine/lipid-derived carbonyl reaction mixtures were heated at **37 °C for 24 h,
> there was a correlation (r = −0.83205, p = 0.02)** between the Eₐ determined for the different
> reactions and the amount of phenylacetaldehyde produced. On the contrary, at high temperature
> other reactions should be competing … and **there was not any correlation** between Eₐ and the
> amount of phenylacetaldehyde produced when reactions were heated either at **80 °C for 2 h
> (r = −0.37, p = 0.4)** or at **180 °C for 1 h (r = −0.24, p = 0.6)**."*

| condition | r (Eₐ vs yield) | p | significant? |
|---|---:|---:|---|
| **37 °C, 24 h** | **−0.832** | **0.02** | ✅ |
| 80 °C, 2 h | −0.37 | 0.4 | ❌ |
| 180 °C, 1 h | −0.24 | 0.6 | ❌ |

**⇒ ⭐⭐ A LOW ACTIVATION ENERGY PREDICTS HIGH YIELD ONLY AT LOW TEMPERATURE. At process
temperature the ranking is set by competing channels, not by Eₐ.** *"At low temperatures, only
reactions having very low Eₐ are favored."* **This is the cleanest statement in the entire
corpus of why an Eₐ table cannot be turned into a yield prediction, and it is measured, with a
p-value, three times.** Register as a **standing caution on every Eₐ row in the repo.**

### 3d. **The reactivity orderings at the three conditions, verbatim [M]**
- **1 h at 180 °C:** *oxononenal > oxohexenal ∼ decadienal > epoxyheptenal ∼ epoxydecenal ∼ heptadienal > hydroxynonenal ∼ control*
- **2 h at 80 °C:** *oxononenal > decadienal > heptadienal ∼ epoxydecenal ∼ oxohexenal ∼ epoxyheptenal > hydroxynonenal ∼ control*
- **24 h at 37 °C:** *decadienal > heptadienal ∼ epoxydecenal ∼ oxononenal ∼ epoxyheptenal ∼ oxohexenal ∼ hydroxynonenal ∼ control*

**⇒ The ordering INVERTS between 37 °C and 180 °C** (decadienal first at 37 °C, oxononenal first
at 180 °C). Quantified [M]: *"the amount of phenylacetaldehyde produced by oxononenal was
**25 % higher** than … decadienal (and **90 % higher** than … heptadienal) when the reaction was
heated for 1 h at 180 °C. This difference **decreased to 21 %** compared to decadienal when
samples were heated for 2 h at 80 °C. Different from these results, **decadienal produced more
phenylacetaldehyde than oxononenal** when … heated at 37 °C for 24 h."*
**⭐ At 37 °C, 4-HNE is statistically indistinguishable from the no-carbonyl control
(5.6 ± 0.2 vs 4.8 ± 0.1, both letter c/cf). The most-cited biological lipid-peroxidation
electrophile is inert for Strecker chemistry at body temperature.**

---

## 4. TABLE 1 — REACTION PRODUCTS AND RETENTION INDICES **[M]**

*Caption: "Retention Indices of the Products Formed in Phenylalanine/Oxoalkenal Reaction Mixtures"*

| compound | retention index | origin |
|---|---:|---|
| styrene | **893** | amine elimination from Phe |
| **2-ethylpyrrole** ᵃ | **915** | lipid (4-oxo-2-hexenal only) |
| benzaldehyde | **960** | Phe, one carbon shorter than PAC |
| **phenylacetaldehyde** ⭐ | **1047** | **the Strecker aldehyde — the quantified species** |
| phenylethylamine | **1097** | Phe decarboxylation |
| **2-pentylpyrrole** ᵇ | **1198** | lipid (4-oxo-2-nonenal only) |
| benzenepropionic acid | **1332** | ⚠️ *"identification … should be considered only **tentative**"* |
| cinnamic acid | **1467** | |
| **2-ethyl-1-phenethyl-1H-pyrrole** ᵃ | **1604** | carbonyl–amine adduct (4-oxo-2-hexenal) |
| **2-pentyl-1-phenethyl-1H-pyrrole** ᵇ | **1885** | carbonyl–amine adduct (4-oxo-2-nonenal) |

ᵃ *only in the 4-oxo-2-hexenal reaction* · ᵇ *only in the 4-oxo-2-nonenal reaction.*
All identified by RI **and** mass spectra against pure reference compounds, **except**
benzenepropionic acid.

**⇒ The 2-alkylpyrrole / 2-alkyl-1-phenethylpyrrole pair is the mechanistic proof: the
alkyl chain tracks the parent oxoalkenal (ethyl from C6, pentyl from C9). The lipid is not a
catalyst — it is consumed and incorporated.**

---

## 5. THE OTHER MEASURED DEPENDENCES **[M]**

### 5a. pH — Figure 1
> *"The amount of phenylacetaldehyde produced by phenylalanine degradation in the presence of
> oxoalkenals **decreased linearly (r < −0.983, p < 0.0001) as a function of reaction pH**"*
(4-oxo-2-hexenal and 4-oxo-2-nonenal, **180 °C, 1 h**, pH 2.15–10).
**⭐⭐ [!] THE pH DEPENDENCE IS THE OPPOSITE SIGN TO EVERY OTHER CARBONYL–AMINE RESULT IN THE
CORPUS.** `hidalgo1993_extraction.md` Table 1 measures **base catalysis** (kB rises 116× from
pH 3 to pH 12) for epoxyheptenal + lysine browning; `yuan2023_extraction.md` §10 cites
*"lower at pH 3 and higher at pH 7 or pH 8"* for flavor–BLG adduction;
`meynier2004_extraction.md` §6e is consistent with amine free-base availability governing.
**Here the Strecker YIELD falls monotonically with pH.** The most likely reconciliation [Z]:
these are **different steps** — free-base amine availability (base-favoured) governs
condensation, while the **decarboxylation/hydrolysis** that liberates phenylacetaldehyde is
**acid-favoured**, and it is the latter that this paper's observable reports. **Do not merge
the two pH rules. Register `pH_rule: Strecker_aldehyde_yield_ACID_favoured (this paper) ≠
adduct_formation_BASE_favoured (hidalgo1993, yuan2023)`.**
*Reaction yields at optimum: **"∼20 %"** [M].*

### 5b. Carbonyl concentration — Figure 2
> *"the amount of the produced aldehyde **increased as a function of oxoalkenal concentration**
> … **This increase was not linear** and **higher increases were observed at low concentrations
> of oxoalkenals** than when high concentrations … were employed."*
(180 °C, 1 h, pH 3.) **⇒ Saturating in carbonyl. Register `order_in_carbonyl: <1, saturating`.**

### 5c. The proposed pathway — Figure 7 **[Q]**
> *"The first step … is the formation of the corresponding **imine** between the amino group of
> the amino acid and the carbonyl group of the oxoalkenal. **Because the oxoalkenal has two
> carbonyls, two imines can be produced.** The second step is the **loss of carbon dioxide**.
> This loss induces an electronic rearrangement and the formation of **a new imine, which, after
> hydrolysis, is the origin of phenylacetaldehyde**. As a consequence … the oxoalkenal is
> transformed into a **hydroxylamino derivative**, which can be later converted into a **pyrrole
> derivative by dehydration** or be polymerized."*

And the class-level unification [Q]:
> *"This similarity among the proposed reaction pathways for all tertiary lipid oxidation
> products is **in agreement with the similarity found for the Eₐ** of the different reactions
> involving tertiary lipid oxidation products (**55.3−67.1 kJ/mol**, Table 2)."*

**⇒ The paper's own claim is that a shared Eₐ band implies a shared mechanism. Its own
alkadienal rows (27.6 air / 78.0 N₂) are the counter-example that proves the converse:
a DIFFERENT Eₐ implies a different rate-determining step — which is exactly why
`shepelev2024`'s 15–20 kJ/mol for direct aldehyde–protein adduction sits outside this band
and should.**

---

## 6. WHAT THIS PAPER DOES AND DOES NOT SUPPLY THE REPO

| question | answer |
|---|---|
| Eₐ for lipid-carbonyl → Strecker aldehyde, per class, at process T | ✅ **eight values, §3** |
| Eₐ for aldehyde + PROTEIN covalent adduct formation | ❌ **no protein in this paper** |
| Eₐ for a saturated n-alkanal | ❌ **every carbonyl here is α,β-unsaturated and activated** |
| Whether Eₐ predicts yield | ✅ **only at 37 °C; measured, §3c** |
| Atmosphere dependence of Eₐ | ✅ **2.8×, measured, §3b** |
| pH rule for Strecker aldehyde liberation | ✅ **acid-favoured, linear, §5a** — and it conflicts in sign with the adduct-formation rule |
| Reversibility | ❌ not tested |
| n / replicate count | ❌ *"several independent determinations"*, never stated **[!]** |

**SI: this paper has no Supporting Information. Nothing further to retrieve.**

---

## 7. INTEGRITY FLAGS

| # | flag | severity |
|---|---|---|
| **1** | **The Eₐ values are for a DOWNSTREAM Strecker product, not for adduct formation. Never substitute them into an adduct sink term.** | ⭐ HIGH |
| **2** | **2,4-decadienal has TWO Eₐ values (27.6 air / 78.0 N₂) differing by 2.8×. A single-valued row is wrong.** | ⭐ HIGH |
| **3** | Eₐ predicts yield only at 37 °C (r = −0.83, p = 0.02); not at 80 °C or 180 °C (p = 0.4, 0.6) | ⭐ HIGH |
| 4 | Zero-order initial-rate estimator; `k` has concentration/time units | MEDIUM |
| 5 | n never stated for any value; SDs printed without a sample size | MEDIUM |
| 6 | Eₐ rests on 4 points (oxoalkenals), 6 (alkadienals) or 7 (epoxyalkenals/HNE) — not stated in the table | MEDIUM |
| 7 | pH dependence has the **opposite sign** to the adduct-formation pH rule elsewhere in the corpus | ⭐ HIGH (do not merge) |
| 8 | Benzenepropionic acid identification is *"only tentative"* | LOW |
| 9 | `research_round4_nulls.md` §A.1's abstract-derived ranges omit the 78.0 kJ/mol N₂ row entirely and mis-round two endpoints | MEDIUM (superseded by §3a here) |
