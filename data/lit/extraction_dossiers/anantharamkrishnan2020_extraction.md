# Anantharamkrishnan & Reineccius 2020 — COMPLETE TRANSCRIPTION
### Covalent flavour–β-lactoglobulin adduct formation vs pH, temperature and water activity — **and why it yields no fittable number.**

**Full extraction of every number in `data/articles/anantharamkrishnan2020.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **Vaidhyanathan Anantharamkrishnan\*, Gary A. Reineccius\*** |
| Title | **"Influence of pH, Temperature, and Water Activity on Covalent Adduct Formation between Selected Flavor Compounds and Model Protein β-Lactoglobulin"** |
| Venue | ***J. Agric. Food Chem.*** — ASAP article, **2020**; the PDF footer reads *"J. Agric. Food Chem. XXXX, XXX, XXX−XXX"*, i.e. **no volume or page numbers were assigned at the time of this PDF** |
| DOI | **10.1021/acs.jafc.0c06752** · file ID **jf0c06752** |
| Dates | Received 25 Oct 2020; Revised 29 Oct 2020; **Accepted 30 Oct 2020**; PDF created 10 Nov 2020 |
| Affiliation | **Department of Food Science and Nutrition, University of Minnesota, St. Paul, MN 55108** |
| Funding | **Minnesota Agriculture Experiment Station** |
| Keywords as printed | *protein, flavor, covalent interaction, pH, temperature, water activity, mass spectrometry* |
| PDF character | 11 pages, clean text layer. **All 16 figures are deconvoluted ESI mass spectra.** |

**⇒ The wave brief's expectation — "the covalent flavor-compound + beta-lactoglobulin adduct
study vs pH/T/a_w, JAFC 2020" — is confirmed exactly.** The brief's further expectation,
**"adduct yields/rates vs conditions"**, is **NOT met by the source**: see §1.

**Provenance codes:** **[M]** measured and printed · **[C]** cited · **[Q]** qualitative
statement, no number attached · **[Z]** derived by this wave.

---

## 1. THE ONE-PARAGRAPH ANSWER — **⚠️ READ THIS BEFORE PLANNING ANY FIT**

**This paper contains no table, no rate constant, no yield, no percentage conversion, no
activation energy, no equilibrium constant, and no calibration curve.** Every result is a
**deconvoluted ESI mass spectrum reproduced as a figure**, and every conclusion is stated as a
direction (*"more slowly at low pH"*, *"increased with temperature"*, *"led to an increase"*).
**The only numbers in the entire results section are molecular masses in daltons**, plus the
protein concentration (**1 wt %**), the flavour dose (**12 ppth**), the three pH values, the
three temperatures, the seven water activities and the three time points. **There is nothing
here to fit and nothing here to hold out as a numeric target.** What the paper *does* deliver
— and it is genuinely valuable — is **a complete, mass-verified inventory of which covalent
adducts form on which residue class under which conditions**, and a **hard directional result
that carbonyl–lysine adduct formation is abolished at pH 3 and proceeds at pH 7–8** (§3).
**⚠️ For the repo's hexanal over-prediction question specifically: the four flavour compounds
studied are benzaldehyde, citral, allyl isothiocyanate and dimethyl tri-/disulfide.
HEXANAL IS NOT AMONG THEM, and neither is any n-alkanal.** The paper is therefore a
**mechanism analogue**, not a mechanism measurement, for that problem (§6).

---

## 2. SYSTEM COMPOSITION AND METHOD — complete

### 2.1 The protein

| variable | value as printed |
|---|---|
| **Protein** | **β-lactoglobulin (BLG), single variant A**, isolated in-house |
| Source | **milk from a single homozygous cow (BLG Variant A), cow number 2708**, Dairy Barn, St. Paul campus, University of Minnesota |
| Isolation | defatted by centrifugation **16 000 rpm, 1 h, 4 °C**; casein precipitated at **pH 4.6 with 1 M HCl**, removed on Whatman No. 1; whey adjusted to **pH 7.0**, passed through a **10 kDa size-exclusion membrane**; concentrate adjusted to **pH 3.8**, **heated and held at 45 °C for 30 min**; precipitated α-lactalbumin, BSA and immunoglobulin removed by centrifugation (**16 000 rpm, 1 h, 4 °C**); supernatant passed through a **10 kDa membrane** and **freeze-dried** |
| **Purity** | *"a protein fraction comprising **95 % β-lactoglobulin**"* |
| **Measured mass of the unmodified protein** | **18 362 Da** — *"corresponds to the A Variant of BLG found in nature"* [C] Bewley et al. 1997 |
| **Lysine residues per BLG molecule** | **15** (stated as the basis for the dosing choice) |
| **Cysteine residues per BLG molecule** | **five, four of which are engaged in disulfide bonds with each other** ⇒ **one free cysteine** |
| **Lactose adduct** | **18 686 Da = BLG + 324 Da**, present in the control spectrum |
| ⚠️ variant change vs the authors' own prior work | *"The relative intensity of the adduct peak in 24 h is **reduced in comparison to our previous publication**. This could be because **BLG variant (A) is used in this study whose reactivity and native structure is different than the B variant** that was used in the previous study. It could also be because of the **variation in the glycation level between the two variants**"* ⇒ **the same experiment on variant B gave different adduct intensities.** [Q] |

### 2.2 The flavour compounds

| compound | purity as printed | used for pH & T study | used for a_w study | reactive target |
|---|---|:--:|:--:|---|
| **Benzaldehyde** | **> 99.5 %** | ✅ | ✅ | **lysine ε-amine (Schiff base)** |
| **Citral** (mix of cis and trans) | **> 96 %** | ✅ | ✅ | **lysine (Schiff base + aza-Michael)** |
| **Allyl isothiocyanate (AITC)** | **> 95 %** | ✅ | ✅ | **lysine (addition)** |
| **Dimethyl trisulfide (DMTS)** | **> 98 %** | ✅ | ❌ | **cysteine / disulfide exchange** |
| **Dimethyl disulfide (DMDS)** | **> 99 %** | ❌ | ✅ | **cysteine** |

All from **Sigma-Aldrich (St. Louis, MO)**, analytical grade; **HCl ACS reagent 37 %**;
**all water doubly distilled.**
⚠️ **DMTS and DMDS are not the same compound and are not interchangeable.** The pH and
temperature arms used **DMTS**; the water-activity arm used **DMDS**. The abstract mixes them:
it says *"dimethyl trisulfide"* for the pH result and *"dimethyl disulfide"* for the a_w
result — correctly, but a reader will not notice the substitution. **Record the compound per
arm.**
Compound selection basis: *"These flavor compounds were chosen because they have been shown to
react with BLG protein in our previous study"* [C] ref. 10.

### 2.3 The three experimental arms

| | **pH arm** | **temperature arm** | **water-activity arm** |
|---|---|---|---|
| **protein form** | **1 wt % aqueous BLG solution** | **1 wt % aqueous BLG solution** | **freeze-dried BLG powder** |
| **variable** | **pH 3, 7, 8** — adjusted with **0.5 M HCl** or **0.5 M NaOH**, measured on a VWR SympHony meter | **4 °C, ambient (ca. 20 °C), 45 °C** — at **pH 7** | **a_w 0.11, 0.23, 0.33, 0.44, 0.53, 0.64, 0.75** |
| **equilibration** | — | equilibrated to each temperature before dosing | **6 WEEKS in a desiccator** over the respective saturated salt |
| **salts** | — | — | **LiCl 0.11 · CH₃CO₂K 0.23 · MgCl₂·6H₂O 0.33 · K₂CO₃ 0.44 · Mg(NO₃)₂ 0.53 · NaNO₃ 0.64 · NaCl 0.75** |
| **flavour dose** | **12 parts per thousand (ppth)**, added individually and vortexed | **12 ppth** | **12 ppth**, added individually and mixed thoroughly |
| **vessels** | — | — | **amber vials, closed tightly, stored in the dark** |
| **time points** | **1 h, 6 h, 24 h** | **1 h, 6 h, 24 h** | **24 h and 48 h** |
| **sample prep for MS** | **diluted 1:10 with water** | **diluted 1:10 with water** | **powder mixed with water to a 0.1 % solution** |
| **temperature of the a_w arm** | — | — | **room temperature** (stated in the figure captions: *"incubated at room temperature"*) |
| ⚠️ **background lactose adducts** | *"low levels (closer to noise) since they were done shortly after protein isolation and the samples were maintained at low a_w (**0.18**) until use"* | same | *"conducted **6 months after protein extraction**, during which time, **lactose adducts spontaneously formed** and thus are observed in all of the mass spectra"* |

**⚠️ THE DOSE IS THE DESIGN'S DEFINING CHOICE, and the authors state it plainly:**
> *"The high concentrations of model flavor compounds used in this study were chosen since
> **this high level would provide adequate flavor compound to occupy all lysine residues (15
> per BLG molecule)**. This would provide **a measure of the maximum potential for reaction**."*

**⇒ 12 ppth = 12 g/kg = 1.2 % w/w. Every result in this paper is a maximum-potential,
saturating-dose result. It bounds what CAN happen, not what DOES happen at a food-relevant
flavour loading (typically ppm to ppb).** [Z] For scale: a hexanal off-note in a plant-protein
matrix is a **µg/kg** problem; 12 ppth is **10⁶–10⁷×** larger.

### 2.4 The instrument

| variable | value as printed |
|---|---|
| System | **Waters Acquity UPLC** + **Waters Synapt G2/Si HDMS quadrupole orthogonal-acceleration TOF** (UPLC/qTOF/MS) |
| Column | **Waters Acquity BEH C4, 2.1 mm × 100 mm, 1.7 µm, at 35 °C** |
| Gradient | **15 min at 0.400 mL/min**; A = water + 0.1 % formic acid, B = acetonitrile + 0.1 % formic acid; **3 % B 0–3 min; 3→97 % B 3–9 min; 97 % B 9–11 min; 97→3 % B 11–13 min; 3 % B 13–15 min** |
| MS acquisition | **profile mode, m/z 300–2500, every 0.2 s** |
| ESI+ parameters | **capillary 0.5 kV; sampling cone 35.0 V; extraction cone 4.0 V; desolvation gas 800 L/h; source 100 °C; desolvation 350 °C; cone gas 40 L/h; trap CE off** |
| Lockspray | **0.5 µg/mL leucine-enkephalin**, one spectrum (0.2 s, m/z 300–2500) **every 10 s**, three measurements averaged |
| Deconvolution | **MaxEnt (Waters)**, mass range **10 000 : 50 000 Da**, resolution **0.1 Da/channel**, **uniform Gaussian, width at half height 0.33 Da**, iterations to autoconvergence |
| **Stated mass accuracy** | **"The results of the deconvoluted spectra reported are accurate to ± 1 Da."** |
| ⚠️ **quantification** | **NONE.** No internal standard, no calibration curve, no response factor, no normalisation. All comparisons are of **relative peak intensity within a single deconvoluted spectrum**. |

---

## 3. THE COMPLETE ADDUCT INVENTORY — every mass shift the paper reports [M]

**This is the paper's real deliverable and it is transcribed exhaustively.** All masses in Da;
all are ± 1 Da per §2.4. Base protein **18 362 Da**.

| compound | **mass shift** | **assigned structure as printed** | residue | conditions where observed |
|---|---:|---|---|---|
| **Benzaldehyde** | **+88** | **Schiff base** | lysine ε-NH₂ | pH 7, pH 8; **4 °C (reduced), 20 °C, 45 °C**; a_w **0.75 only** (as 18 450 Da) |
| **Citral** | **+152** | **aza-Michael addition adduct** | lysine | pH 7, pH 8; all temperatures; all a_w |
| **Citral** | **+268** | *"could be **two Schiff base adducts** formed with multiple lysine groups **OR** the **pyridinium adduct** formed with a single lysine group and two citral molecules"* ⚠️ **ambiguous assignment, stated as such** | lysine | pH 7, pH 8 |
| **Citral** | **+402** | *"additional Schiff base and Michael adducts at the other reactive sites (e.g. **lysine and arginine**)"* | lysine / arginine | pH 7, pH 8 |
| **Citral** | **+554** | " | " | pH 7, pH 8 |
| **Citral** | **+88** (as **18 450 Da**) | the adduct whose intensity rises from 4 → 20 °C | lysine | temperature arm |
| **AITC** | **+99** | **addition to lysine** | lysine | all pH; all temperatures; all a_w |
| **AITC** | multiple, to **18 659 Da** = **+297** (three adducts) at a_w 0.75/48 h; **18 560 Da** = **+198** (two) and **18 461 Da** = **+99** (one) at pH 7/8 and in the temperature arm | *"up to **six adducts** can be observed"* at a_w 0.75 | lysine (15 available) | see §4 |
| **DMTS** | **+46** | first cysteine adduct | **free cysteine** | pH 7, pH 8; 20 °C, 45 °C |
| **DMTS** | **+78** | second cysteine adduct | free cysteine | pH 7, pH 8 |
| **DMTS** | **increments of +32** | *"corresponding to the mass of a **sulfur atom** ... This suggests that the adducts are formed **after cleaving the Cys−Cys disulfide linkage**"* | **Cys–Cys disulfide** | **pH 8 only** |
| **DMDS** | **+47** (as **18 409 Da**) | not assigned | cysteine | **a_w 0.75 only** |
| **dehydration product** | **−18** (**18 344 Da**) | *"18 362 (BLG) − 18 (water molecule)"* | protein | AITC at pH 7 and 8 |
| **dehydration + adduct** | **18 443, 18 452, 18 641 Da** | *"adducts of this dehydrated compound"* | protein | AITC at pH 7 and 8 |
| **lactose (background)** | **+324** (**18 686 Da**) | Amadori lactose adduct | lysine | **all spectra in the a_w arm**; near noise in the pH and T arms |

**⇒ 14 distinct mass shifts, all mass-verified to ±1 Da, spanning four reaction chemistries:
Schiff base, aza-Michael addition, isothiocyanate addition, and disulfide exchange.** This is
a genuinely complete mechanism inventory and it is the reason to keep the paper.

---

## 4. THE THREE DIRECTIONAL RESULTS — verbatim, with every qualification the authors attach

### 4.1 pH: **3 < 7 < 8**, and pH 3 abolishes carbonyl chemistry entirely

| compound | **pH 3** | **pH 7** | **pH 8** |
|---|---|---|---|
| **Benzaldehyde** | **"no measurable adducts formed"** — *"In acidic conditions, the intermediate step for forming the covalent adduct does not form since the molecule would be protonated"* | **+88 Schiff base observed** | **+88 observed** |
| **Citral** | **"there were no adducts observed at the lowest pH"** — *"Both mechanisms ... require the participation of the **free base, primary amines**"* | **+152, +268, +402, +554** | same, but *"The intensities of the unreacted protein and adduct peaks were **reduced at the basic pH** ... due to the **cross-linking of the protein by citral**, thereby making the mass of the cross-linked molecule **above the detection limits of the MS**"* [C] ref. 18 |
| **AITC** | *"the adducts **form but at a slower rate** than in comparison to neutral and basic conditions"* | multiple adducts; major peak **18 560 Da** | *"the intensity of the protein peak reduced and **ultimately disappeared**"*; major peak **18 659 Da** |
| **DMTS** | **"there were no appreciable adducts formed"** | **+46 and +78** | **multiple adducts in +32 increments** ⇒ disulfide cleavage |

**Summary as printed:** *"for the types of covalent adduct formation studied (i.e. **Schiff
base, Michael addition, and disulfide exchange reaction**), **the rate of covalent adduct
formation increases at higher pH, consistent with the general Maillard reaction**"* [C] refs.
20–25. *"At a basic pH, **multiple adducts were observed that were not observed at the lower pH
values, suggesting that the structure of the protein could be changing and reactive sites that
were stable at acidic pH are exposed in basic pH values.**"*

**⚠️ AITC is the exception that matters: it is the ONLY compound that forms adducts at pH 3.**
Isothiocyanate addition does not require a free base amine in the same way a Schiff base does.
**⇒ pH gating is mechanism-specific, not compound-generic.**

### 4.2 Temperature: **4 °C < 20 °C ≤ 45 °C**, and cross-linking destroys the signal at 45 °C

| compound | **4 °C** | **20 °C** | **45 °C** |
|---|---|---|---|
| **Benzaldehyde** | *"the +88 Da Schiff adduct could be observed **but to a lesser extent**"* | — | ⚠️ *"There were **no appreciable differences** in the adducts formed for the samples incubated at 20 and 45 °C ... consistent with the previous studies where no difference in reactivity was observed for samples stored at 25 and 45 °C"* [C] refs. 26–28 |
| **Citral** | *"an **increase in the intensity of the adduct peak (mass 18 450 Da)** for samples incubated at 4 and 20 °C"* | " | ⚠️ **peaks at NOISE LEVEL** — *"likely due to **extensive cross-linking of the protein**, making them too large in molecular weight to be measured by the MS"*; *"The enals are known to cross-link proteins"* [C] refs. 9, 18, 29 |
| **AITC** | major peak **18 362 Da** (unreacted) | major peak **18 461 Da** (+99) | ⚠️ **peaks in the noise level** — *"The increase in temperature appears to have **induced cross-linking**, resulting in the molecular weights of the adducts being **above the mass range of the ESI/MS measurement**"* |
| **DMTS** | *"only a small amount of adduct ... a small spike at +46 Da, but it was **not significant**"* | **+46 Da present at 24 h** | *"**multiple adducts** were observed. The major unreacted protein peak had **completely disappeared** ... At 45 °C, the **Cys−Cys disulfide bond would be cleaved**"* |

**⚠️ FOR THREE OF FOUR COMPOUNDS, THE 45 °C SPECTRUM IS UNINTERPRETABLE.** The signal vanishes
into the noise because the products leave the **10 000–50 000 Da** MaxEnt window. **The method
loses sensitivity exactly where the reaction is fastest.** The authors state this repeatedly
and honestly. **⇒ The "adducts increase with temperature" conclusion is inferred from the
DISAPPEARANCE of the monomer peak, not from the appearance of a quantified product.**

**Proposed AITC cross-linking mechanism (Figure 10):** *"the addition reaction between allyl
isothiocyanate and BLG is **partially reversible**. When the reaction reverses, **the
isothiocyanate group could be formed with the protein instead of the allyl group**, (b) product
in Figure 10. This would further react to **induce cross-linking or form higher-molecular-
weight compounds**"* [C] ref. 28.
**Proposed DMTS mechanism (Figures 11–12):** *"Once the free cysteine group is reacted, the
disulfide bond present in the Cys−Cys would be cleaved and the adduct with the DMTS would be
formed ... It can be a **mix and match of the combination between the +46 and +78 Da ions**."*

### 4.3 Water activity: **only AITC responds; three of four compounds show nothing below a_w 0.75**

| compound | result as printed |
|---|---|
| **AITC** | *"An increase in water activity led to an increase in the rate of formation of adducts ... due to the **increase in mobility of the flavor molecule to reach the reactive amino acid moiety**"* [C] ref. 30. **At a_w 0.11 the major peak is 18 362 Da (unreacted) with +99 adducts; at a_w 0.75 up to SIX adducts, major peak 18 659 Da.** Time dependence: at a_w 0.75, **major peak 18 461 Da at 24 h → 18 659 Da at 48 h**. *"For lower water activities, the rate of adduct formation is slower and therefore **no significant change is observed between time points 24 and 48 h**."* |
| **Citral** | *"citral forms covalent adducts with BLG at **all water activities**. However, there was **no observable difference** between the samples at different water activities. This could be because of the **shorter reaction period studied** — 24 h and 48 h. As time progresses, there could be larger differences."* |
| **Benzaldehyde** | *"**no adduct peaks formed for water activities between 0.11 and 0.64 in 48 h**. However ... There is a peak with **relatively low intensity, corresponding to a molecular weight of 18 450 ... at 0.75 water activity**. There is **no observable increase in the intensity of the adducts between 24 and 48 h**."* |
| **DMDS** | same as benzaldehyde: nothing to a_w 0.64; **a low-intensity 18 409 Da peak at a_w 0.75 only** |

Mechanistic framing as printed: *"In general, one would expect that an increase in water
activity (**to some point, some studies have shown it to be 0.75**) would increase reactant
solubility and mobility and thereby affect the rate of chemical reactions, especially the
Maillard reaction. **The reaction rate would decrease with the further increase in a_w, as the
moisture content would dilute the concentration of reactants.**"* [C] refs. 30–32.
Background lactose adducts as an internal a_w probe: *"Lactose adducts (18 686 Da, BLG +324 Da)
can be observed in Figure 13a,b at 0.11 a_w and they **generally increased in intensity with
increases in storage water activity**. Previous studies have shown that **lysine-lactose adduct
formation increases with a_w to a maximum around 0.6−0.7**"* [C] refs. 31, 32.

**⇒ THE LACTOSE ADDUCT IS THE PAPER'S ONLY INTERNAL POSITIVE CONTROL FOR THE a_w ARM, and it
behaves as the literature predicts. That is worth recording: it validates the a_w
equilibration and shows the powder was chemically live across the whole range**, which makes
the three null results (benzaldehyde, citral, DMDS) genuine nulls rather than a failed
experiment.

The authors' own closing qualification: *"**this lack of observed effect may be due to the
rate of reaction being too slow to be detected in the timeframe of this study**"* (abstract),
and *"**Further investigations on the effect of the environmental conditions on reaction rates
and extent are needed over longer reaction times**."*

---

## 5. THE FINAL SUMMARY, VERBATIM — the paper's own statement of what it established

> *"Results from this study have demonstrated that the covalent binding of flavor molecules with
> BLG **varies with pH, temperature, and water activity. The rate (and, therefore, extent) of
> covalent bond formation increased with pH (from pH 3, 7, 8), temperature (4−20 °C), and,
> somewhat, water activity (0.11−0.75).** The results also show that **the effect of these
> variables on the rate and extent of the reaction depends on the reactive group (functional
> group) on the flavor compound, which, in turn, dictates the nature of the adduct-forming
> reaction.**"*

⚠️ Note the temperature range in that sentence is **"4−20 °C"**, not 4–45 °C — the authors
themselves exclude the 45 °C point from the claim, because at 45 °C the signal is lost to
cross-linking (§4.2).

---

## 6. ⚠️ RELEVANCE TO THE REPO'S HEXANAL OVER-PREDICTION — an honest assessment

The wave brief nominated this paper as *"a candidate mechanism for the repo's known hexanal
over-prediction in protein matrices."* **The assessment is: the mechanism is right, the
measurement is absent, and the compound is missing.**

**What supports the hypothesis:**
1. **Benzaldehyde — an aldehyde — forms a Schiff base with BLG lysine at pH 7 and 8, and forms
   NOTHING at pH 3.** That is the exact chemistry `k2_matrix_and_thresholds.md` §B.3 named as
   the reason headspace-derived and dialysis-derived **aldehyde** binding constants differ by
   7–35×, and it is now demonstrated at the level of an intact-protein mass shift.
2. **Citral — an α,β-unsaturated aldehyde — does more: Schiff base AND aza-Michael addition,
   AND protein cross-linking at 45 °C.** That is the mechanism behind
   `k2_matrix_and_thresholds.md` §D.3's ≈2–3× α,β-unsaturation threshold penalty and behind
   Meynier 2002's **6.88× t-2-hexenal suppression in skim milk** (this wave). **Three
   independent lines now converge on 2-alkenals being irreversibly consumed by protein.**
3. **The pH gate is sharp and mechanism-specific**, and it is corroborated in this same wave by
   `leksrisompong2010_extraction.md` §7, where diacetyl binds caseinate at **pH 7.0 (P < 0.05)
   and not at pH 5.5 (N.S.)** — a different protein, a different carbonyl, the same direction.

**What does NOT support using it quantitatively:**
1. **Hexanal is not studied.** Neither is any straight-chain alkanal. Benzaldehyde is
   **aromatic** — its carbonyl is conjugated to a ring, which changes both its electrophilicity
   and the stability of the resulting imine. **Transferring a benzaldehyde result to hexanal is
   an extrapolation across compound class.**
2. **No yield, no rate, no half-life, no equilibrium constant.** There is **no number** to put
   into a loss term.
3. **The dose is saturating by design (12 ppth ≈ 1.2 % w/w, 10⁶–10⁷× a food-relevant hexanal
   loading).** The paper measures the **maximum potential** for reaction. A repo that needs to
   know *what fraction of 50 µg/kg hexanal is lost to lysine in 2 hours* gets no help.
4. **The protein is 1 wt % purified BLG in water, or freeze-dried BLG powder — not a food
   matrix.**
5. **The 45 °C data are uninterpretable** for three of four compounds (§4.2), and 45 °C is at
   the bottom of the temperature range the repo's thermal lane cares about.

**⇒ VERDICT: ingest as a MECHANISM CONSTRAINT with directional force and zero numeric content.
Specifically usable as three qualitative gates — `carbonyl_lysine_adduct: {pH3: none, pH7:
yes, pH8: yes}`, `alkenal_crosslinks_protein: yes, above ~45 °C`, and `low_a_w_suppresses_
Schiff_base: yes, threshold near a_w 0.75 for benzaldehyde`. It cannot close the hexanal
over-prediction gap by itself, and no wave should plan to fit it.** The paper that *would*
close it is **Fares et al. 1998** (quantitative covalent diacetyl–caseinate binding with two
site classes) or the authors' own **ref. 9/10** — see §9.

---

## 7. NAMED LAUNDERING HAZARDS

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| A-1 | **"The rate (and, therefore, extent) of covalent bond formation increased with pH, temperature and water activity"** | **No rate was measured.** There is no rate constant, no yield, no calibration and no internal standard anywhere in the paper. Every conclusion is read from the **relative height of peaks in a deconvoluted spectrum**, which MaxEnt does not preserve quantitatively across samples. The word "rate" in this sentence means "apparent order of appearance". | Conclusion, p. I |
| A-2 | **"An increase in temperature led to an increase in the rate of formation of adducts"** for AITC and citral | For both compounds the **45 °C spectra are at NOISE LEVEL** — the products left the 10–50 kDa MaxEnt window through cross-linking. The increase is inferred from the **disappearance of the monomer**, not from any observed product. The authors say so; a citing paper will not. | §4.2 |
| A-3 | The abstract's **"dimethyl trisulfide"** (pH/T arm) and **"dimethyl disulfide"** (a_w arm) | **Two different compounds were used in different arms** — DMTS for pH and temperature, DMDS for water activity. Correct as written, but trivially conflated. **The a_w result does not apply to DMTS and the pH result does not apply to DMDS.** | §2.2 |
| A-4 | **12 ppth** flavour dose, presented without comment in the results | **1.2 % w/w — chosen explicitly to saturate all 15 lysines and measure "the maximum potential for reaction."** It is **10⁶–10⁷×** a food-relevant aroma loading. Every directional result is a saturating-dose result. | p. B |
| A-5 | **The +268 Da citral adduct**, listed alongside the unambiguous ones | The paper offers **two mutually exclusive assignments** — *"could be two Schiff base adducts ... or the pyridinium adduct"* — and does not resolve them. The +402 and +554 assignments are likewise hedged (*"e.g. lysine and arginine"*). **Three of five citral adducts are structurally unassigned.** | §3 |
| A-6 | Adduct intensities compared against the authors' **previous publication** | *"The relative intensity of the adduct peak in 24 h is reduced in comparison to our previous publication ... because **BLG variant (A) is used in this study** whose reactivity and native structure is different than the **B variant**."* **The same experiment on the other natural variant of the same protein gives different results.** Any repo record must carry `BLG_variant: A`. | p. B |
| A-7 | **Lactose adducts** present in every a_w-arm spectrum | The a_w arm was run **6 months after protein isolation**, during which lactose adducts *"spontaneously formed"*. **The starting protein in the a_w arm is not the same material as in the pH and temperature arms** (which were run at a_w 0.18 shortly after isolation, near-noise lactose). The three arms are **not on a common protein baseline.** | §4.3 |
| A-8 | Citation of **"Gremli, H. A.; Dimitrova, N.; ... Interaction of Flavor Compounds with Soy Protein. *Chem. Res. Toxicol.* 2010, 51, 1050−1059"** (ref. 2) | Gremli's soy-protein flavour-interaction paper is ***J. Amer. Oil Chem. Soc.* 51 (1974) 95A–97A**, single-authored. The reference as printed has the wrong journal, the wrong year, the wrong page range and a fabricated author list. ⚠️ **The same Gremli paper is cited correctly by Meynier 2002 (this wave, `Meynier2002_extraction.md` §9) as the source of soy hexanal/t-2-hexenal retention values.** Do not follow this citation. | Ref. 2, p. J |

---

## 8. CONSOLIDATED NEW-PARAMETER TABLE

**⚠️ There are no kinetic or thermodynamic parameters in this paper.** What follows is the
complete set of measured quantities plus the directional constraints, which is what the source
actually contains.

**Common conditions:** **β-lactoglobulin variant A, 95 % pure, measured mass 18 362 ± 1 Da,
15 lysines, 1 free cysteine + 2 disulfides**; flavour at **12 ppth**; **UPLC-ESI-qTOF, ±1 Da,
MaxEnt 10 000–50 000 Da**; **no quantification of any kind.**

| # | parameter | value | units | conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1 | **BLG variant A molecular mass** | **18 362 ± 1** | Da | measured, control spectrum, 24 h | M | Fig. 1 |
| 2 | **BLG + lactose (Amadori) adduct** | **18 686 (= +324)** | Da | background; rises with a_w | M | Fig. 1, §4.3 |
| 3 | **BLG dehydration product** | **18 344 (= −18)** | Da | AITC, pH 7 and 8 | M | §3 |
| 4 | **benzaldehyde–lysine Schiff base** | **+88** | Da | pH 7, 8; 4/20/45 °C; a_w 0.75 only | M | §3 |
| 5 | **citral aza-Michael adduct** | **+152** | Da | pH 7, 8 | M | §3 |
| 6 | **citral +268 adduct** | **+268** ⚠️ **two competing assignments** | Da | pH 7, 8 | M (structure Q) | §3 |
| 7–8 | **citral higher adducts** | **+402, +554** | Da | pH 7, 8; lysine and/or arginine | M (structure Q) | §3 |
| 9 | **citral adduct at 18 450 Da** | **+88** | Da | temperature arm, rises 4 → 20 °C | M | §4.2 |
| 10 | **AITC–lysine adduct, single** | **+99** (18 461 Da) | Da | all pH, all T, all a_w | M | §3 |
| 11 | **AITC, two adducts** | **+198** (18 560 Da) | Da | pH 7 | M | §3 |
| 12 | **AITC, three adducts** | **+297** (18 659 Da) | Da | pH 8; a_w 0.75 at 48 h | M | §3, §4.3 |
| 13 | **AITC, maximum adducts observed** | **six** | count | a_w 0.75, 48 h | M | §4.3 |
| 14–15 | **DMTS–cysteine adducts** | **+46, +78** | Da | pH 7, 8; 20 and 45 °C | M | §3 |
| 16 | **DMTS disulfide-cleavage series** | **increments of +32 (sulfur)** | Da | **pH 8 only** | M | §3 |
| 17 | **DMDS–cysteine adduct** | **+47** (18 409 Da) | Da | **a_w 0.75 only** | M | §4.3 |
| 18 | **BLG lysine count** | **15** | residues/molecule | — | M | p. B |
| 19 | **BLG cysteine count** | **5, of which 4 in disulfides ⇒ 1 free** | residues/molecule | — | M | §4.1 |
| 20 | **flavour dose (saturating by design)** | **12** | **ppth (= 12 g/kg = 1.2 % w/w)** | chosen to occupy all 15 lysines | M | p. B |
| 21 | **protein concentration** | **1** | wt % aqueous | pH and temperature arms | M | p. B |
| 22 | **a_w equilibration time** | **6** | weeks | desiccator over saturated salt | M | p. B |
| 23 | **a_w levels** | **0.11, 0.23, 0.33, 0.44, 0.53, 0.64, 0.75** | a_w | LiCl/CH₃CO₂K/MgCl₂/K₂CO₃/Mg(NO₃)₂/NaNO₃/NaCl | M | p. B |
| 24 | **pH levels** | **3, 7, 8** | pH | 0.5 M HCl / 0.5 M NaOH | M | p. B |
| 25 | **temperatures** | **4, ~20, 45** | °C | pH 7 | M | p. B |
| 26 | **time points** | **1, 6, 24 h** (pH, T); **24, 48 h** (a_w) | h | — | M | p. B |
| **C1** | **carbonyl–lysine adduct formation is ABOLISHED at pH 3** | benzaldehyde **none**; citral **none**; DMTS **none appreciable**; **AITC forms adducts, more slowly** | **directional gate** | 1 wt % BLG, 12 ppth, 24 h, ~20 °C | **Q** | §4.1 |
| **C2** | **adduct formation increases pH 3 → 7 → 8** | direction only; at pH 8, adducts appear that do not exist at pH 7 | **directional gate** | " | **Q** | §4.1 |
| **C3** | **adduct formation increases 4 °C → 20 °C** | direction only; ⚠️ **the authors' own claim stops at 20 °C** | **directional gate** | pH 7, 12 ppth | **Q** | §4.2, §5 |
| **C4** | **α,β-unsaturated aldehyde (citral) CROSS-LINKS BLG at 45 °C** | qualitative — the monomer peak vanishes into noise | **mechanism** | pH 7, 12 ppth, 24 h | **Q** | §4.2 |
| **C5** | **AITC cross-links BLG at 45 °C via a partially reversible addition** | qualitative; mechanism drawn in Fig. 10 | **mechanism** | " | **Q** | §4.2 |
| **C6** | **DMTS cleaves the Cys–Cys disulfide at pH 8 and at 45 °C** | qualitative; mechanism drawn in Figs. 11–12 | **mechanism** | " | **Q** | §4.1, §4.2 |
| **C7** | **water activity gates Schiff-base formation** | benzaldehyde and DMDS: **nothing from a_w 0.11 to 0.64; a trace at a_w 0.75.** AITC: monotone increase 0.11 → 0.75. Citral: **no a_w effect in 48 h** | **directional gate** | powder, 12 ppth, room temperature | **Q** | §4.3 |
| **C8** | **lysine–lactose adduct formation increases with a_w** | observed here as an internal positive control; **literature maximum at a_w 0.6–0.7** | **directional + cited** | powder | Q + C | §4.3 |
| **C9** | **the response depends on the flavour compound's functional group, not on the compound** | the paper's own closing conclusion | **scoping rule** | — | **Q** | §5 |
| — | **rate constants, yields, half-lives, activation energies, equilibrium constants** | **NONE EXIST IN THIS PAPER** | — | — | — | §1 |

---

## 9. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Anantharamkrishnan & Reineccius 2020 is a **new source**; a declaration
> amendment would be required before any wave fits it. **This dossier's recommendation is that
> no amendment is needed, because there is nothing fittable.**

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **The 14 adduct mass shifts (§3)** | 14 | **INGEST as a MECHANISM INVENTORY — not a fit target, not a hold-out** | Mass-verified to ±1 Da. Their value is that they identify **which** reaction occurs on **which** residue for **which** functional-group class. That is structural metadata, not a parameter. |
| **C1 — no carbonyl–lysine adduct at pH 3** | 1 | **HOLD-OUT — the wave's strongest mechanism gate** | Corroborated independently in this same wave by `leksrisompong2010_extraction.md` §7 (diacetyl binds caseinate at pH 7, not pH 5.5). **Proposed guard: any repo aldehyde-loss-to-protein term must go to zero at acid pH.** |
| **C2, C3, C7 — the pH, temperature and a_w directions** | 3 | **HOLD-OUT — directional gates only** | No magnitudes exist. Usable to reject a model that predicts the wrong sign; useless to calibrate one. ⚠️ C3 must be recorded as **4 → 20 °C only** (§5). |
| **C4 — alkenals cross-link protein** | 1 | **HOLD-OUT — supports, and is supported by, three other corpus results** | `k2_matrix_and_thresholds.md` §D.3 (≈2–3× α,β-unsaturation threshold penalty, two matrices), `Meynier2002_extraction.md` §5.1 (6.88× t-2-hexenal skim-milk suppression vs 1.39× for hexanal), and §B.3's aldehyde method-split. **Four independent lines. This is now the best-supported irreversible mechanism in the corpus.** |
| **C5, C6 — AITC and DMTS mechanisms** | 2 | **METADATA** | Outside the repo's compound scope. |
| **C8 — lactose adducts rise with a_w** | 1 | **METADATA / positive control** | Validates that the a_w-arm nulls are real nulls. |
| **Everything else** | — | **NOT INGESTIBLE** | No rates, no yields, no calibration, no internal standard (A-1). |
| **The 45 °C spectra for citral, AITC and benzaldehyde** | 3 arms | **REJECT** | At noise level; products outside the MaxEnt window (A-2). |
| **All results** | — | **must carry `dose: 12_ppth_saturating` and `BLG_variant: A`** | A-4 and A-6. |

**⇒ Recommended orchestrator action: record Anantharamkrishnan 2020 as a
`mechanism_reference` in the citation graph, NOT as a data source, and do NOT amend the FIT
declaration for it.** If the repo needs a *quantitative* covalent aldehyde-loss term, the
retrievals in §10 are the route.

---

## 10. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Anantharamkrishnan, V., Hoye, T. & Reineccius, G. (2020)**, *J. Agric. Food Chem.* **68**,
   6395–6402 — *"Covalent Adduct Formation Between Flavor Compounds of Various Functional Group
   Classes and the Model Protein β-Lactoglobulin"* (ref. 10). **The parent survey**, covering
   **"various functional group classes"** — which is the only plausible place an **n-alkanal**
   (and possibly hexanal itself) would appear. **Highest-value retrieval from this paper by a
   wide margin**, and directly targeted at the repo's hexanal question.
2. **Anantharamkrishnan, V. & Reineccius, G. A. (2020)**, *J. Agric. Food Chem.*,
   DOI 10.1021/acs.jafc.9b07978 (ref. 9) — *"Method To Characterize and Monitor Covalent
   Interactions of Flavor Compounds with β-Lactoglobulin Using Mass Spectrometry and
   Proteomics."* The **bottom-up (digest + MS/MS) companion**, which identifies **which specific
   lysine residues** are modified. It is also the study that used **variant B** and got
   **stronger adduct signals** (A-6) — so it is the only way to bound the variant effect.
3. **Kühn, J., Considine, T. & Singh, H. (2008)**, *J. Agric. Food Chem.* **56**, 10218–10224 —
   *"Binding of Flavor Compounds and Whey Protein Isolate as Affected by Heat and High Pressure
   Treatments"* (ref. 27). A **quantitative** whey-protein binding study with a **heat**
   treatment axis — i.e. the quantitative complement to this paper's qualitative temperature
   arm.
4. **Adams, R. L., Mottram, D. S., Parker, J. K. & Brown, H. M. (2001)**, *J. Agric. Food Chem.*
   **49**, 4333–4336 — *"Flavor-Protein Binding: Disulfide Interchange Reactions between
   Ovalbumin and Volatile Disulfides"* (ref. 25). **The quantitative disulfide-exchange paper**
   this one's DMTS/DMDS arm is qualitative about — directly relevant to the repo's sulfur lane.
5. ⚠️ **Do NOT use this paper's ref. 2 to locate Gremli.** The correct citation is
   **Gremli, H. A. (1974), *J. Amer. Oil Chem. Soc.* 51, 95A–97A** (see A-8 and
   `Meynier2002_extraction.md` §13).
