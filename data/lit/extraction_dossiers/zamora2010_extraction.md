# Zamora, Delgado & Hidalgo 2010 — COMPLETE TRANSCRIPTION
### *"Model Reactions of Acrylamide with Selected Amino Compounds"*
### ⭐⭐⭐ **THE CORPUS'S FIRST — AND ONLY — MEASURED FORWARD/REVERSE ACTIVATION-ENERGY PAIR FOR A CARBONYL-CHEMISTRY ADDUCT: 44 kJ/mol forward, 52 kJ/mol reverse.** Wave **K6b**, 2026-08-30.

**Full extraction of every number in `data/articles/zamora2010.pdf`.**
Read-only wave. No repo file outside `data/lit/extraction_dossiers/` created or modified.

**Read first:** `hidalgo2010_extraction.md` (its thiol counterpart — same lab, same year, same
electrophile, **opposite reversibility verdict**) and `research_round4_nulls.md` §A.1, which
recorded of this paper: *"**no numbers in the abstract**; the abstract states only: 'acrylamide
seems to be converted into its Michael adduct with a **lower** activation energy than the
elimination reaction of the Michael adduct'"*. **§3 below supplies the numbers.**

**Provenance:** **[M]** measured · **[C]** cited · **[Q]** qualitative · **[Z]** derived here · **[!]** flag.

---

## 0. IDENTITY — **MATCHES EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Rosario Zamora, Rosa M. Delgado, Francisco J. Hidalgo\*** |
| Venue | ***J. Agric. Food Chem.* 2010, 58 (3), 1708–1713** |
| DOI | **10.1021/jf903378x** ✔ **exactly the expected DOI** |
| Affiliation | **Instituto de la Grasa, CSIC, Avenida García Tejero 4, 41012 Seville, Spain** |
| Published on Web | **01/15/2010**; received Sep 25 2009, revised Dec 16 2009, accepted Dec 19 2009 |
| Corresponding | fhidalgo@ig.csic.es · fax +34 954 616 790 |
| Keywords | *Acrylamide; 3-(alkylamino)propionamides; amines; amino acids; **carbonyl-amine reactions**; **Hoffman elimination**; Maillard reaction; **Michael addition*** |
| Funding | EU FEDER · **Junta de Andalucía P07-AGR-2846** · Plan Nacional I+D **AGL2006-01092, AGL2009-07638**; José L. Navarro thanked |
| PDF character | **6 pages**, clean text layer, **Tables 1–2**, Figures 1–8 |

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER

**This paper measures both directions of one reaction and reports both activation energies.**

| direction | reaction | **Eₐ (kJ/mol)** |
|---|---|---:|
| **forward** | acrylamide + **butylamine** → 3-(butylamino)propanamide (aza-Michael) | **44** |
| **reverse** | 3-(butylamino)propanamide → acrylamide + butylamine (elimination) | **52** |
| forward | acrylamide + **glycine** → Michael adduct | **52** |

**⇒ `Eₐ(reverse) − Eₐ(forward) = +8 kJ/mol`. [Z] Therefore ΔH_rxn ≈ −8 kJ/mol: the adduct is
favoured, but by only 8 kJ/mol, and the equilibrium constant FALLS with temperature by
`exp(8000/R · (1/T₁ − 1/T₂))`. From 25 °C to 180 °C, `K_eq` drops 3.0×.** (§4.)

**⭐⭐ THAT IS THE STRUCTURAL RESULT THE COVALENT-SINK QUESTION NEEDS, AND IT IS THE ONLY
MEASUREMENT OF ITS KIND IN THE CORPUS.** A conjugate-addition adduct of a soft electrophile at
an amine nitrogen — the same bond class as an aza-Michael aldehyde adduct, one step from the
same class as a Schiff base — has a **higher barrier out than in**, so **heating does not merely
fail to build the adduct pool faster than it destroys it; heating actively shifts the
equilibrium back toward free electrophile.** The paper demonstrates this experimentally and
dramatically in §5: acrylamide that had **completely disappeared** during 14 days at 60 °C was
**recovered by reheating at 180 °C for 20 min**.

**And the paper names the food-safety consequence itself [Q]:** *"they also point to
**3-(alkylamino)propionamides as possible compounds in which acrylamide might be potentially
hidden**"* — an analyte that is absent by measurement and present by chemistry.

---

## 2. SYSTEM AND METHOD — complete transcription **[M]**

| variable | value as printed |
|---|---|
| **Electrophile** | **acrylamide, 0.2 µmol** |
| **Nucleophiles** | **butylamine · glycine · lysine · two poly-L-lysine hydrobromides**, at **0–50 µmol** |
| **Poly-L-lysine MWs** | **4 200** and **68 300** (both *"determined by viscosity"*, Sigma certificate of analysis) |
| **Support** | **0.063–0.200 mm silica gel 60 (300 mg)**, Macherey-Nagel |
| **Buffer** | **30 µL of 0.3 M** — sodium citrate pH 3–6, sodium phosphate pH 6–8 — **+ 180 µL water** |
| **Water activity** | **a_w = 0.95**, measured with a **Pawkit Decagon analyzer** |
| **Standard condition** | **heater block, under nitrogen, 180 °C, 10 min**, unless otherwise indicated |
| **pH maintenance** | *"The reaction pH was **maintained upon heating**."* |
| **Cooling** | ⚠️ **"15 min at −20 °C"** *(as printed; the companion `hidalgo2010` prints "15 min at 20 °C" for the identical step — a sign discrepancy between the two papers [!])* |
| Internal standard | **10 µL of 1 mg/mL [1,2,3-¹³C₃]acrylamide in methanol** |
| Extraction | + **2 mL of 0.3 M sodium citrate pH 2.2**, stirred **1 min**, filtered |
| Derivatization / GC-MS | **identical to `hidalgo2010_extraction.md` §2** (KBr/Br₂ → 2-bromopropenamide; m/z 149 analyte, 154 IS; HP5-MS; 50 → 240 → 325 °C) |
| Calibration | **15 levels, 0–200 µg**, prepared in the silica gel and carried through; **r = 0.999, p < 0.0001**; **CV < 10 %** |
| Replication | *"mean ± SD for, **at least, two independent experiments**"* (Table 2) |

### 2a. The two synthesized reference adducts **[M]** — full characterization

**Compound 3 = 3-(butylamino)propanamide.** Acrylamide (**14.5 mmol**) in methanol (**29 mL**)
+ butylamine (**14.5 mmol**), **60 °C for 24 h**, solvent evaporated, crystallized and
recrystallized from **toluene/hexane (1:1)**.
- **¹H NMR (CD₃OD):** δ 0.95 t (3H, J = 7.3, CH₃), 1.37 sx (2H, J = 7.3), 1.50 qu (2H, J = 7.3), 2.43 t (2H, J = 6.9, CH₂CONH₂), 2.61 t (2H, J = 7.3), 2.84 t (2H, J = 6.9, CH₂CH₂NH)
- **¹³C NMR (CD₃OD):** δ 14.28, 21.46, 32.44, 35.33 (CH₂CONH₂), 46.38 (CH₂CH₂NH), 50.06, **177.22 (CONH₂)**
- **MS:** m/z **144 (1, M⁺)**, 101 (55), 84 (65), 56 (22), **44 (100)**, 42 (77)

**Compound 4 = 2-(bis(3-amino-2-oxopropyl)amino)acetic acid** *(the bis-adduct of glycine)*.
Acrylamide (**3.75 mmol**) + glycine (**7.5 mmol**) in water (**14.5 mL**), **2 M KOH** to basic
pH, **60 °C for 24 h**, then **2 M acetic acid** to acid pH; recrystallized from water.
- **¹H NMR (D₂O):** δ 2.65 t (4H), 3.37 t (4H), 3.63 s (2H, CH₂COOH)
- **¹³C NMR (D₂O):** δ 31.16, 53.72, 58.82, **172.65 (CONH₂)**, **176.81 (COOH)**
- **MS of the TMS derivative:** m/z **433 (3, M⁺)**, 416 (1), 343 (5), 316 (28), 305 (15), **185 (100)**, 173 (74), 144 (62), 128 (15), 73 (67)

⚠️ **[!] THE MONO-ADDUCT OF GLYCINE COULD NOT BE MADE.** *"Although this Michael adduct (formed
by reaction of **one** molecule of glycine and **one** molecule of acrylamide) **could not be
prepared**, the thermal heating of an **analogous** Michael adduct (formed by reaction of one
molecule of glycine and **two** molecules of acrylamide) produced acrylamide very rapidly to a
high extent (Figure 5)."* **⇒ The glycine reverse reaction is demonstrated on the BIS-adduct
only, and no Eₐ is reported for it. Only the butylamine pair (§3) is a matched forward/reverse
Eₐ measurement.**

---

## 3. ⭐⭐⭐ THE ACTIVATION ENERGIES **[M]**

### 3a. Temperature set and rate law
> *"Acrylamide disappearance **increased linearly as a function of time** at the four assayed
> temperatures … By employing the disappearance rates obtained from the **slopes of the lines in
> Figure 2A**, it was possible to determine the activation energy (Eₐ) … by means of an
> **Arrhenius plot (Figure 3A)**."*

**Temperatures (Figs. 2 and 4 captions): 180 (□), 160 (○), 140 (△) and 120 °C (▽) — FOUR
temperatures, at pH 7.** *(Fewer than the eight in `hidalgo2010`; weight accordingly.)*

### 3b. **THE THREE VALUES — the corpus's reversibility pair**

| # | reaction | figure | **Eₐ (kJ/mol)** |
|---|---|---|---:|
| **F1** | **acrylamide + butylamine → adduct 3** (0.2 : 20 µmol) | Fig. 2A → Fig. 3A ○ | **44** |
| **R1** | **adduct 3 → acrylamide** (0.2 µmol, alone) | Fig. 2B → Fig. 3A △ | **52** |
| **F2** | **acrylamide + glycine → adduct** (0.2 : 20 **and** 0.2 : 50 µmol) | Fig. 4A,B → Fig. 3B | **52** |

Verbatim: *"The activation energy (Eₐ) for the disappearance of acrylamide in the presence of
**butylamine was 44 kJ/mol**."* · *"The activation energy of **compound 3 decomposition was
52 kJ/mol**."* · *"The activation energy (Eₐ) found for the disappearance of acrylamide in the
presence of **glycine was 52 kJ/mol**, which is **slightly higher** than that above-described for
butylamine/acrylamide mixtures."*

**⭐ A NOTEWORTHY INTERNAL CONTROL ON F2:** *"When **50 µmol** of glycine was employed, reaction
rates were **slightly higher** (Figure 4B). However, **the same activation energies for both
reactions were determined** by means of an Arrhenius plot (Figure 3B). As observed in this last
figure, **two lines were obtained, but the slopes were identical**."*
**⇒ Eₐ is invariant to a 2.5× change in nucleophile concentration while the rate is not.
This is a real, printed demonstration that the estimator separates the pre-exponential from
the barrier — exactly the property `shepelev2024_extraction.md` §6f relies on.**

### 3c. **The paper's own statement of the asymmetry, verbatim [M]**
> Abstract: *"the activation energies (Eₐ) of both reactions **are not the same**. In fact,
> acrylamide seems to be converted into its Michael adduct **with a lower activation energy
> than the elimination reaction** of the Michael adduct."*
> Discussion: *"Thus, when butylamine was employed, the activation energy of the first reaction
> was **44 kJ/mol** and the activation energy of its reverse was **52 kJ/mol**."*

⚠️ **[!] ALL THREE VALUES ARE BARE INTEGERS: no R², no SE, no CI, four temperatures each.**
The **difference** of 8 kJ/mol is the load-bearing quantity and it is **smaller than the
plausible uncertainty on either value taken alone**. §4 therefore states the consequence in
sign-and-order terms and gives the arithmetic conditional on the difference being real.
**Register `Ea_forward = 44 ± ~4`, `Ea_reverse = 52 ± ~4`, `ΔEa = +8, sign asserted by the
authors, magnitude uncertain`.**

---

## 4. ⭐⭐ [Z] WHAT THE 44/52 PAIR IMPLIES — the arithmetic

`K_eq(T) = k_f/k_r ∝ exp(−(Eₐ,f − Eₐ,r)/RT) = exp(+8000/RT)` × (A_f/A_r).

| T | exp(8000/RT) | **relative K_eq** (normalized to 25 °C) |
|---|---:|---:|
| **25 °C** | 25.3 | **1.00** |
| **60 °C** | 17.8 | **0.70** |
| **100 °C** | 12.9 | **0.51** |
| **140 °C** | 9.94 | **0.39** |
| **180 °C** | 8.42 | **0.33** |

**⇒ The adduct/free-electrophile equilibrium constant falls 3.0× from 25 °C to 180 °C.**
Equivalently: **ΔH ≈ −8 kJ/mol — the adduct is only marginally enthalpically favoured, and
the reaction is entropically penalized (two molecules → one), so heating releases it.**

**⭐⭐ THE TRANSFERABLE PRINCIPLE, AND IT IS THE ONE THAT MATTERS FOR THE COVALENT SINK:**
> **For a reversible conjugate addition at an amine nitrogen, the reverse activation energy
> exceeds the forward one. A covalent sink built on such a bond is therefore MOST effective at
> LOW temperature and is partly unwound by processing.**

**This is the mechanistic explanation for the otherwise puzzling temperature inversion that
`shepelev2024_extraction.md` §4d and §7d measure directly — decanal bound to BLG peaks at
25 °C, not 65 °C, and declines at 65 °C after day 7.** Two laboratories, two chemistries, two
decades apart, one conclusion. `k6b_adduct_kinetics_synthesis.md` §4 develops the convergence.

⚠️ **The transfer's limit, stated honestly:** the bond here is a **C–N σ bond from aza-Michael
addition to an activated alkene**, not a **C=N Schiff base** from an aldehyde. Schiff-base
hydrolysis has its own (generally low) barrier and its equilibrium is water-activity-dependent
in a way this system, fixed at a_w 0.95 on silica, cannot report. **Register as
`analogue_bound`, not as a Schiff-base measurement.**

---

## 5. ⭐⭐ THE STORAGE-AND-REHEAT EXPERIMENT — TABLE 2, TRANSCRIBED IN FULL **[M]**

*Caption verbatim: "**Table 2.** Effect of Storage at 37 or 60 °C and Reheating at 180 °C on
Acrylamide Disappearance ᵃ"
Footnote a: "Values are **acrylamide disappearance (%)** in a model system containing
**0.2 µmol of acrylamide and 20 µmol of glycine**. Values are **mean ± SD for, at least, two
independent experiments**. Means with different letters in the same column are significantly
different (p < 0.05). **nd, not determined**."*

| time (days) | **storage 37 °C, non-reheated** | **37 °C, reheated** | **storage 60 °C, non-reheated** | **60 °C, reheated** |
|---:|---:|---:|---:|---:|
| **0** | **0.0 ± 2.8** b | **77.8 ± 1.1** b | **0.0 ± 2.8** b | **77.8 ± 1.1** b |
| **7** | **51.5 ± 3.2** c | nd | **85.3 ± 1.6** c | nd |
| **14** | **77.1 ± 4.5** d | nd | **97.6 ± 0.1** d | nd |
| **21** | **81.6 ± 3.4** d | nd | **99.0 ± 0.1** d | nd |
| **28** | **87.6 ± 3.8** d | **90.8 ± 3.5** c | **99.3 ± 0.1** d | **89.4 ± 3.2** c |

*"Reheated" = the stored sample then heated **20 min at 180 °C**.*

### 5a. **How to read this table — and it is counter-intuitive**
The column is **disappearance**, so **recovery = 100 − disappearance**:

| sample | disappearance | **acrylamide RECOVERED** |
|---|---:|---:|
| unstored, then reheated 180 °C/20 min | 77.8 % | **22.2 %** |
| stored 28 d at 37 °C (no reheat) | 87.6 % | 12.4 % |
| **stored 28 d at 37 °C, then reheated** | **90.8 %** | **9.2 %** |
| stored 28 d at 60 °C (no reheat) | 99.3 % | 0.7 % |
| **stored 28 d at 60 °C, then reheated** | **89.4 %** | **10.6 %** |

**⭐⭐ THE 60 °C ROW IS THE RESULT.** After 28 days at 60 °C the acrylamide is **99.3 % gone
(0.7 % left)**. Twenty minutes at 180 °C brings it back to **10.6 %** — a **15-fold recovery of
an analyte that had, to all measurement, disappeared.**

Verbatim: *"when acrylamide was stored in the presence of glycine at 60 °C (Table 2),
acrylamide **disappeared after 14 days**. However, when samples heated for 28 days at 60 °C were
**heated again for 20 min at 180 °C, the equilibrium was reestablished and a significant amount
of acrylamide was detected in samples in which acrylamide was previously absent**. In fact,
when both unheated samples and samples heated for 28 days at 37 or 60 °C were heated for
20 min at 180 °C, **all assayed samples showed very similar amounts of acrylamide contents**,
therefore suggesting that **these last conditions are sufficient to reestablish the
equilibrium**."*

**⇒ THE 180 °C/20 min TREATMENT IS AN EQUILIBRATION, AND ALL THREE HISTORIES CONVERGE ON IT
(22.2 %, 9.2 %, 10.6 % recovered — the last two statistically indistinguishable, both letter c).
This is a direct experimental demonstration that a covalent-adduct pool built at storage
temperature is substantially released by a process-temperature step.**

### 5b. **How much of the loss is real — the authors quantify it [M]**
> *"this content was **not exactly the same**. In fact, the content of previously stored samples
> was **lower**. This difference … might be a consequence of the **destruction of acrylamide
> during storage**, which should take place by other reaction mechanisms, such as **oxidation or
> polymerization**. Different from the formation of the Michael adduct, **these acrylamide losses
> are real** because acrylamide cannot be produced again … Nevertheless, under the reaction
> conditions employed in this study, **these other reactions contributed to only ∼10 % of
> acrylamide losses** observed in samples stored for 1 month at either 37 or 60 °C."*

**⇒ [Z] REVERSIBLE : IRREVERSIBLE ≈ 90 : 10 over one month of storage.** ⭐ **The corpus has
never had this split measured for any carbonyl-chemistry sink. Register
`storage_sink_reversible_fraction: ~90 % (acrylamide + glycine, 37–60 °C, 28 d) [M]`.**

---

## 6. TABLE 1 — THE pH SERIES — TRANSCRIBED IN FULL **[M]**

*Caption: "**Table 1.** Effect of pH on Acrylamide Determined after Heating either an
Acrylamide/Amino Compound Model System or Compound 3 ᵃ"
Footnote a: "Values are **acrylamide disappearance (%)** in acrylamide/amino compound
(**0.2 µmol / 20 µmol**) model systems and **acrylamide formation (%)** in the thermal
decomposition of 3-(butylamino)propanamide (**3**). Samples were heated for **10 min at 180 °C**."
ᵇ sodium citrate buffers · ᶜ sodium phosphate buffers · ᵈ not determined.*

| **pH** | buffer | **acrylamide/butylamine** *(disappearance %)* | **compound 3** *(formation %)* | **acrylamide/glycine** *(disappearance %)* |
|---:|---|---:|---:|---:|
| **2.15** | citrate | **64.8** | **8.8** | **30.5** |
| **3** | citrate | 69.1 | 14.7 | 28.0 |
| **4** | citrate | 71.1 | 37.1 | 40.2 |
| **5** | citrate | 73.1 | 50.0 | 52.7 |
| **6** | citrate | **74.5** | **60.0** | **52.9** |
| **6** | phosphate | **70.0** | **27.3** | **53.5** |
| **6.5** | phosphate | 63.0 | 24.3 | nd |
| **7** | phosphate | 61.9 | 24.4 | 49.4 |
| **8** | phosphate | 64.4 | 30.0 | 55.1 |

### 6a. ⚠️⚠️ **[!] THE BUFFER CHANGES THE RESULT BY 2.2× AT CONSTANT pH — AND THIS TABLE PROVES IT**
**At pH 6, in citrate, compound 3 releases 60.0 % of its acrylamide. At pH 6, in phosphate,
it releases 27.3 %.** Same pH, same compound, same time, same temperature; **2.20× apart.**
The forward reaction moves too (74.5 → 70.0, a 6 % drop).

**⇒ THIS IS THE CLEANEST BUFFER-CONFOUND DEMONSTRATION IN THE CORPUS, because the paper
deliberately measured the same pH in both buffers.** It retrospectively validates the flags
raised in `hidalgo1993_extraction.md` §2 (the pH 7.5 → 8 k_B dip at the phosphate→borate
boundary) and `hidalgo2010_extraction.md` §6b. **Register a corpus-wide rule:
`pH_series_across_buffer_boundaries: NEVER fit a single curve; the buffer is a covariate,
and its measured effect on a reverse reaction here is 2.2×.`**

### 6b. The pH readings themselves [M][Q]
- **Butylamine forward: essentially flat.** *"the pH **did not seem to play a major role** in acrylamide/butylamine reactions, and **similar amounts of acrylamide were recovered in the range of pH 2.15−8**."* (61.9–74.5 %, a 1.2× span.)
- **Compound 3 reverse: strongly pH-dependent in citrate.** *"the acrylamide produced **increased with pH in the range of 2.15−6**, when using citrate buffer, but it was **less pH dependent with phosphate buffer** in the range of 6−8."* (**8.8 → 60.0 %, a 6.8× span in citrate.**)
- **Glycine forward: rises then flattens.** *"the amount of disappeared acrylamide **increased as a function of pH** when using citrate buffer, and **no clear effect** was observed at pH 6−8."* (28.0 → 55.1 %.)

**⇒ [Z] The REVERSE reaction is far more pH-sensitive than the forward one (6.8× vs 1.2×).
Combined with §4's Eₐ asymmetry, the equilibrium position of an aza-Michael adduct is governed
by pH and temperature mainly through the ELIMINATION step. A repo sink term that models only
the forward rate will get the pH and temperature responses of the net sink wrong.**

---

## 7. THE NUCLEOPHILE COMPARISONS **[M]**

### 7a. Concentration dependence — Figure 1 (N₂, 180 °C, 10 min)
- **Butylamine:** *"an amount of **∼12 µmol of butylamine was needed to reduce acrylamide content to half**"* (against 0.2 µmol acrylamide — **a 60 : 1 excess for 50 % conversion**).
  The authors' explanation [Q]: *"Although acrylamide and butylamine react equimolecularly, the much higher concentration … needed may be a consequence of both the **low concentration of acrylamide** employed … and a **low reaction yield**."*
- **Glycine:** *"Acrylamide content **decreased linearly (r = −0.993, p = 0.00076)** in the presence of glycine in the range of **0−15 µmol** of the amino acid, and **higher amounts of glycine did not produce further losses**."* ⇒ **saturating at 15 µmol.**
- Control: *"acrylamide **remained unchanged when the amine was not present**"* — a clean blank.

### 7b. Butylamine vs glycine vs lysine — Figure 6 (0.2 : 20 µmol, pH 7, 180 °C)
> *"only **slight differences** … **Glycine exhibited the lowest reactivity, and lysine was the
> most reactive**. The reactivity of butylamine was somewhat **intermediate** … At **short heating
> times**, butylamine reactivity was **quite similar to that of lysine**. However, at **long
> heating times**, it was close to that of **glycine**."*
**⇒ Ordering: lysine > butylamine > glycine, and butylamine's rank is TIME-DEPENDENT.**
Explanation offered [Q]: *"a higher reactivity of butylamine and **the presence of two amino
groups in the molecule of lysine**."*

### 7c. ⭐ Lysine vs polylysine — Figure 7 (0.2 µmol acrylamide + **2 mg** of each, pH 7, 180 °C)
> *"reactivity was **inversely proportional to the molecular weight**. The **most reactive was
> lysine**, followed by the polylysine with the **lower molecular weight (4 200)**, and the
> polylysine with the **higher molecular weight (68 300) was the least reactive**. This is likely
> related to **the number of amino groups** (the higher molecular weight, the smaller number of
> terminal amino groups) **and also to the accessibility of amino groups** (the higher molecular
> weight of the protein, the **lesser accessibility** of amino groups to acrylamide)."*
> *"there was a **linear relationship (r = 0.989, p = 0.095)** between the amount of acrylamide
> recovered at the end of the heating period and the **natural logarithm of the molecular
> weight** of the amino compound added."*

**⭐⭐ THIS IS THE MOST DIRECTLY REPO-RELEVANT RESULT IN THE PAPER AND IT IS NOT IN THE
ABSTRACT.** A **free amino acid, a 4.2 kDa polymer and a 68.3 kDa polymer of the SAME residue,
at the SAME mass loading**, ranked by reactivity — i.e. **a direct measurement of the
free-amine → polymer-bound-amine reactivity penalty.**
**⇒ Register `polymerization_penalty: reactivity ∝ −ln(MW), free Lys > polyLys-4200 >
polyLys-68300 [M]`. It is the quantitative counterpart of
`shepelev2024_extraction.md` §7e's qualitative principle (*"the number of reactive sites on the
SURFACE … plays a greater role than the total amount of reactive amino acids"*) — and it means
every rate constant measured on a free amino acid OVERSTATES the protein case.**

⚠️ **[!] `p = 0.095` IS NOT SIGNIFICANT** at the paper's own declared α = 0.05, on **three**
points. **The ORDERING is a solid observation; the log-linear FORM is not established.**

### 7d. The bis-adduct question — Figures 5 and 8 **[M][Q]**
*"Although the Michael addition of this adduct to a **second molecule of acrylamide** is
possible, the formation of the corresponding 3,3′-(butylazanediyl)dipropanamide **was not
observed** under the conditions employed."* And: *"it is also **not expected to be produced in
foods** because acrylamide will usually be produced to a much lower extent than the content of
amino compounds."*
**⇒ A measured null for the 2:1 adduct at food-relevant stoichiometry.** Figure 8 is the
summary scheme (compounds 1–4). Figure 5 shows compound **4** (the glycine bis-adduct)
releasing acrylamide *"very rapidly to a high extent"* at pH 7, 180 °C.

---

## 8. WHAT THIS PAPER DOES AND DOES NOT SUPPLY THE REPO

| question | answer |
|---|---|
| **A matched forward/reverse Eₐ pair for an adduct** | ✅ ⭐⭐⭐ **44 / 52 kJ/mol — the corpus's only one** |
| **A measured reversible : irreversible split for a storage sink** | ✅ ⭐ **≈ 90 : 10 over 28 d, §5b** |
| **A demonstration that heating releases a storage-formed adduct** | ✅ ⭐⭐ **15× recovery from 0.7 % to 10.6 %, §5a** |
| **A free-amine vs polymer-bound-amine reactivity penalty** | ✅ ⭐ **lysine > polyLys 4.2k > polyLys 68.3k, §7c** |
| **A quantified buffer confound at constant pH** | ✅ ⭐ **2.2× on the reverse reaction, §6a** |
| An **aldehyde** | ❌ the electrophile is acrylamide |
| A **protein** | ❌ polylysine is the closest surrogate (deliberately: *"to avoid the competition of amino groups with other functional groups"*) |
| Eₐ for the reverse of the **glycine** reaction | ❌ only the bis-adduct was available; no Eₐ reported |
| Statistical support for the Eₐ values | ⚠️ four temperatures, bare integers, no R²/SE/CI |

**SI: none listed. Nothing further to retrieve.**
**One lead [C]: Zamora, Delgado & Hidalgo (2009), *Mol. Nutr. Food Res.* 53, 1512–1520,
"Conversion of 3-aminopropionamide and 3-alkylaminopropionamide into acrylamide in model
systems" (ref. 20) — the mechanism paper for the elimination. LOW priority; the Eₐ is here.**

---

## 9. INTEGRITY FLAGS

| # | flag | severity |
|---|---|---|
| **1** | **At pH 6, changing citrate → phosphate changes the reverse reaction by 2.20× (60.0 % → 27.3 %). Buffer is a covariate, not a control.** | ⭐⭐ HIGH — and it retro-flags `hidalgo1993` and `hidalgo2010` |
| **2** | **The 8 kJ/mol forward/reverse gap is smaller than the plausible uncertainty on either bare-integer value. The SIGN is the paper's claim; the magnitude is soft.** | ⭐ HIGH |
| **3** | The glycine mono-adduct **could not be synthesized**; the reverse reaction is shown only for the **bis**-adduct, and no glycine reverse Eₐ exists | MEDIUM |
| 4 | Eₐ from **four** temperatures; no R², no SE, no CI | MEDIUM |
| 5 | `p = 0.095` on three points is presented as a *"linear relationship"* (§7c) | MEDIUM |
| 6 | Cooling step printed as **"15 min at −20 °C"** here vs **"15 min at 20 °C"** in `hidalgo2010` for the identical protocol | LOW |
| 7 | Butylamine's reactivity **rank changes with heating time** (§7b) — a single-number ranking is unsafe | MEDIUM |
| 8 | The observable is acrylamide **disappearance/formation**, not adduct concentration | MEDIUM |
| 9 | Silica gel at a_w 0.95, 120–180 °C — a low-moisture surface system, not a solution | MEDIUM |
| 10 | **The 44/52 pair is for an aza-Michael C–N σ bond, NOT a Schiff base. Transfer to aldehyde–amine chemistry as an `analogue_bound` only.** | ⭐ HIGH |
