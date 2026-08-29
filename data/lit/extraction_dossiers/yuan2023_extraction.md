# Yuan, Anantharamkrishnan, Hoye & Reineccius 2023 — COMPLETE TRANSCRIPTION
### *"Covalent Adduct Formation between β-Lactoglobulin and Flavor Compounds under Thermal Treatments That Mimic Food Pasteurization or Sterilization"*
### ⭐⭐ **THE COMPOUND-CLASS × TEMPERATURE REACTIVITY MAP — 46 compounds × 4 thermal treatments. And the paper that supplies hexanal's own row.** Wave **K6b**, 2026-08-29.

**Full extraction of every number and every table cell in `data/articles/Yuan2023.pdf`.**
Read-only wave: no repo file outside `data/lit/extraction_dossiers/` was created or modified.

**Read first:** `anantharamkrishnan2020b_extraction.md` (this paper's direct ambient-temperature
predecessor — ref. 17 — whose 46-compound panel is re-run here) and
`shepelev2024_extraction.md` (its quantitative sister from the same lab).

**Provenance codes:** **[M]** measured and printed · **[C]** cited · **[Q]** qualitative ·
**[Z]** derived by this wave · **[!]** integrity flag.

---

## 0. PAPER IDENTITY — **MATCHES THE EXPECTED IDENTITY EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Jieyao Yuan, Vaidhyanathan Anantharamkrishnan, Thomas R. Hoye, Gary A. Reineccius\*** |
| Title | as above |
| Venue | ***J. Agric. Food Chem.* 2023, 71 (24), 9481−9489** |
| DOI | **10.1021/acs.jafc.3c01220** ✔ **exactly the expected DOI** |
| Affiliations | Dept. of Food Science and Nutrition, **U. Minnesota, St. Paul** (Yuan, Anantharamkrishnan, Reineccius); Dept. of Chemistry, **U. Minnesota, Minneapolis** (Hoye) |
| Present addresses | Yuan → Changsha Univ. of Science and Technology, Hunan, China · Anantharamkrishnan → Lux Flavours, Chennai |
| Dates | Received **Feb 25, 2023** · Revised **May 14, 2023** · Accepted **May 17, 2023** · Published **Jun 6, 2023** |
| Keywords | *flavor, protein, thermal processing, covalent reactions, UPLC−ESI−QTOF−MS* |
| Acknowledgments | Dr. Chi Chen and Dana Yao (UPLC−ESI−QTOF−MS); **Dr. Brian Crooker, for milk from a single cow (numbered 2708)** |
| PDF character | **9 pages**, clean text layer, **Figures 1–8**, **Tables 1–4** (all four came out of the text layer intact and are transcribed in full below), SI PDF |

**⇒ CONFIRMED: expected authors, journal, volume, DOI, year, protein (BLG), compound count
(46), and the expected temperature set. The brief's prediction of "22/63/72/110/130 °C" is
correct as a list of the temperatures *attempted*; §3b records that **110 °C yielded no data at
all** and is absent from every table.**

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER

**This paper is an ORDINAL map, not a kinetic one. Its entire result set is a 0/1/2 score,
defined by whether an unreacted-BLG peak remains in the deconvoluted mass spectrum.** There
is no rate, no extent, no percentage, no adduct stoichiometry, no time course and no
activation energy anywhere in it. What it delivers that nothing else in the corpus does is
**the same 46-compound panel scored at four thermal severities, including two that mimic real
UHT and pasteurization processes**, plus **three compounds newly found reactive on heating**.

**The four headline results.**
**(a) ⭐ HEXANAL'S ROW IS FLAT: `1, 1, 1, 1` — the identical partial-reaction score at ambient
22 °C/30 min, HTST 72 °C/15 s, IC 63 °C/30 min AND UHT 130 °C/30 s** (§4). **Hexanal's covalent
reactivity does not visibly increase across a 108 K span. On the paper's own scale, hexanal
carries ZERO temperature information — which is exactly what a small activation energy looks
like on an ordinal scale, and exactly what `shepelev2024_extraction.md` §6 measures
quantitatively for a homologue.**
**(b) The aldehyde class is uniformly reactive and mostly saturating.** All seven non-vanillin
aldehydes score ≥ 1 at every condition; **five of them coagulate the protein at ≥ 1 condition**
(the `*` cells), which is a *stronger* outcome than a score of 2, not a missing value (§3c [!]).
**(c) ⭐ THE SEVERITY ORDERING IS NON-MONOTONE IN TEMPERATURE AND THE AUTHORS SAY SO:
`HTST (72 °C, 15 s) < IC (63 °C, 30 min) ≈ UHT (130 °C, 30 s)`.** A **lower**-temperature,
**longer** treatment matches the highest-temperature one. §6 turns that statement into an
Arrhenius bracket **and then explains why that bracket must not be used as an aldehyde Eₐ**.
**(d) Three compounds cross from inert to reactive on heating: eugenol, 4-vinylphenol,
3-nonen-2-one — and for eugenol the authors show the reactive species is an autoxidation
CONTAMINANT, not eugenol itself** (§5a). A cleanly reported, self-falsified positive.

---

## 2. ⚠️ WHAT THIS PAPER IS NOT — the negative inventory, up front

**Absent entirely:** any rate constant · any half-life · any activation energy · any time
course (each compound × condition is **one time point**) · any replicate count (**n is never
stated for any cell in any table**) · any error bar · any statistical test · any quantification
of adduct amount · any pH measurement or control (**1 wt % BLG in doubly distilled water,
unbuffered, pH never stated**) · any water-activity axis (explicitly deferred to ref. 18) ·
any protein other than BLG · any 110 °C data (§3b) · any measurement between 72 and 130 °C.

**⚠️ [!] THE SCORE IS NOT AN EXTENT.** §Results defines it: *"we have rated the extent of
reaction (R) with a scoring system of **R = 0−2**. An R score of **0** indicates that there was
no observed reaction … while an R score of **2** indicates **complete reaction (i.e., no visible
MS peak for the unreacted protein)** under, **at least, the most forcing UHT conditions**. …
The value of **1** indicates **partial reaction** (i.e., there is visible MS peak for remaining
unreacted protein)."* **⇒ A score of 1 spans everything from a barely-detectable adduct to
99 % conversion. The dynamic range inside "1" is unbounded, and it is where hexanal lives at
every temperature.** Register every row `scale: ordinal_0_1_2_MS_peak_presence`.

**⚠️ [!] AND THE SCORE IS PARTLY A TURBIDITY JUDGEMENT.** Verbatim: *"Based on a combination of
the MS profiles **and turbidity observations**, we have rated the extent of reaction."*
An unaided visual turbidity call is folded into the same integer. **Never treat a 1 → 2 step
as a factor.**

---

## 3. SYSTEM AND METHOD — complete transcription

### 3a. The protein **[M]**
| variable | value as printed |
|---|---|
| Source | milk from **a single homozygous cow with BLG Variant A**, UMN herd |
| **Molecular weight** | **18,363 Da** (stated in Methods); **18,362 Da** used throughout Results; **"18.2 kDa"** in the Introduction — ⚠️ **[!] three different values in one paper; use 18,363** |
| Isolation | fat removed by centrifugation **16,000 rpm, 1 h, 4 °C** → casein precipitated at **pH 4.6**, filtered (Whatman no. 1) → whey neutralized to **pH 7.0** → salts/lactose removed on a **10 kDa** size-exclusion membrane → filtrate to **pH 3.8**, heated **45 °C, 30 min** to precipitate α-lactalbumin, BSA and immunoglobulin → centrifuged **16,000 rpm, 1 h, 4 °C** → salts removed on a second **10 kDa** membrane → freeze-dried |
| **Purity** | **95 % (Dumas analysis)** *(cf. `shepelev2024_extraction.md`'s 90 % for the same lab's protocol)* |

### 3b. The reaction system and the four thermal treatments **[M]**
| variable | value as printed |
|---|---|
| **Protein solution** | **1 weight-percent BLG** in **doubly distilled water** |
| **Flavor dose** | **12 parts per thousand (ppth), weight/weight**, *"**Excess amounts of flavor compounds**"*, added individually, vortexed; solids pre-dissolved in **ethanol** |
| **pH** | ⚠️ **never stated, unbuffered** |
| Vessel (63 °C) | **amber glass vial**, directly in a temperature-controlled **water bath**, **1 min thermal equilibration** allowed |
| Vessel (72 / 110 / 130 °C) | **100 µL** into a **flame-sealed glass capillary (1.5 mm × 150 mm)** on a metal rack in a temperature-controlled **oil bath**, **6 s equilibration** (justified from ref. 9, Datta et al. 2005) |
| Quench | capillary **immediately into ice water**; surface deoiled with **dichloromethane**; sample stored at **−80 °C** |
| LC dilution | **10× with LC−MS grade water** |
| **Control** | pure BLG solution, no flavor, **held at the same time/temperature as each treatment** |

**THE FOUR TREATMENTS, and their real-process meanings [M]:**

| code | condition | mimics | data obtained? |
|---|---|---|---|
| **ambient** | **~22 °C**, held **30 min / 15 s / 30 s** to match each partner | storage | ✅ |
| **HTST** | **72 °C for 15 s** | high-temperature-short-time pasteurization | ✅ |
| **IC** | **63 °C for 30 min** | in-container ("standard") pasteurization | ✅ |
| **UHT** | **130 °C for 30 s** | ultra-high-temperature sterilization | ✅ |
| ~~sterilization~~ | ~~**110 °C for 30 min**~~ | in-container sterilization | ❌ **NO DATA — see below** |

⚠️ **[!] THE MOST SEVERE CONDITION PRODUCED NO DATA AT ALL, FOR ANY COMPOUND.** Verbatim:
*"Unfortunately, our methodology **did not permit us to obtain meaningful data using the most
aggressive standard sterilization thermal conditions (110 °C for 30 min)** because **extensive
aggregation/coagulation removed essentially all of the BLG protein from the reaction mixtures
prior to MS analysis**."* And: *"there are **no useful spectra** for samples heated at 110 °C for
30 min since the BLG **always precipitated**."*
**⇒ The paper's temperature axis is 22 → 130 °C with a hole at 110 °C, and the hole is caused
by the observable disappearing, not by the chemistry stopping. Register
`110C_30min: no_data_protein_lost`.**

### 3c. ⭐⭐ **[!] THE ASTERISK IS NOT A MISSING VALUE — IT IS THE STRONGEST POSSIBLE OUTCOME**
Table 2 footnote *a*, verbatim: *"**\* = unable to measure the extent of reaction due to
complete protein coagulation, leaving no measurable soluble protein.**"*

And the Results text: *"BLG reactions with **some aldehydic flavor molecules** led to
**agglomeration/coagulation of the protein** to the extent that the derivatized protein in
solution **could no longer be observed by LC−ESI−MS** (see \* in Table 2)."*

**⇒ A `*` means the flavor compound cross-linked and precipitated the entire protein
population. Chemically that is MORE reaction than a score of 2, not less.** An automated ingest
that treats `*` as `null` will systematically under-rate exactly the aldehydes the repo cares
about, and will invert the temperature trend for **benzaldehyde (1, \*, \*, \*)**,
**citral (1, 1, \*, \*)** and **E-2-hexenal (2, \*, \*, \*)**.
**Register `asterisk_means: reaction_so_complete_protein_precipitated`, and treat `*` as
`≥ 2` in any ordinal analysis, never as missing.** *(This is the single most important ingest
rule this paper carries.)*

### 3d. Analytical method, complete **[M]**
| item | parameters as printed |
|---|---|
| LC | **Waters Acquity UPLC**, **BEH C4, 2.1 × 100 mm, 1.7 µm**, **35 °C** |
| Mobile phases | A = water + **0.1 % formic acid**; B = acetonitrile + **0.1 % formic acid** |
| Gradient | **15 min**, **0.40 mL/min**: 3 % B 0–3 min; **3→97 % B 3–9 min**; 97 % B 9–11 min; 97→3 % B 11–13 min; 3 % B 13–15 min |
| MS | **Waters Synapt G2/Si**, QTOF, **positive ESI**, **centroid**, **m/z 300–2500 every 0.2 s** |
| MS params | capillary **0.5 kV**; sampling cone **35.0 V**; extraction cone **4.0 V**; desolvation gas **800 L/h**; source **100 °C**; desolvation **350 °C**; cone gas **40 L/h**; **trap collision energy OFF** |
| Lock mass | **leucine enkephalin, [M+H]⁺ = 556.2771**, injected intermittently in real time |
| Deconvolution | **MaxEnt (Waters)**, mass range **10,000–20,000 Da**, **0.1 Da/channel**, **uniform Gaussian, width at half height 0.33 Da**, minimum left/right intensity ratios **33 %**, auto-convergence, **±1 Da accuracy** |

### 3e. The control spectrum — the masses that are always present **[M]**
From Fig. 1 (BLG alone at all four conditions):

| mass (Da) | assignment as printed |
|---:|---|
| **18,362** | **BLG A variant** (the reference peak) |
| **18,686** | **BLG + lactose Schiff base** — *"342 Da (lactose) − 18 (water) = **324 Da**"* |
| **18,416** | *"additional peak … common to all samples"*, **unidentified** (Δ = +54) |
| **18,462** | *"additional peak … common to all samples"*, **unidentified** (Δ = +100) |

*"No effort was devoted to determining their identities because these peaks had very low
intensities."* And critically: *"the control spectra are **essentially identical (qualitatively
not quantitatively) across sample heating conditions**."*

**⭐ [Z] THAT LACTOSE ADDUCT IS ITSELF A RESULT THE CORPUS SHOULD REGISTER.** A **+324 Da
lactose Schiff base is present on BLG in the untreated control at 22 °C** — i.e. **the protein
arrives already partially glycated**, before any flavor compound is added, and the authors
report **no qualitative change in it from 22 °C to 130 °C**. **A lactosylation channel that is
already saturated at ambient and does not visibly grow through UHT is the same signature as
hexanal's flat row (§4) — see `k6b_adduct_kinetics_synthesis.md` §4.**

---

## 4. ⭐⭐ TABLE 2 — CARBONYL COMPOUNDS — TRANSCRIBED IN FULL **[M]**

*Caption verbatim: "**Table 2.** Reaction of Carbonyl Compounds with BLG at Room Temperature
(30 min), 72 °C (15 s; HTST Pasteurization), 63 °C (30 min, Standard Pasteurization), and
130 °C (30, UHT Sterilization)"* *(the caption's "130 °C (30," is missing its "s" — as printed).*
*Column header printed as "**reaction rate**" in every column; footnote b defines it as the
0–2 extent score, so **"rate" here means extent. [!]*

| carbonyl compound | **ambient (30 min)** | **72 °C (15 s)** | **63 °C (30 min)** | **130 °C (30 s)** |
|---|:---:|:---:|:---:|:---:|
| **— aldehydes —** | | | | |
| benzaldehyde | **1** | **\*** | **\*** | **\*** |
| citral (geranial) | **1** | **1** | **\*** | **\*** |
| furfural | **1** | **1** | **1** | **1** |
| *E,E*-2,4-heptadienal | **\*** | **1** | **\*** | **\*** |
| ⭐ **hexanal** | **1** | **1** | **1** | **1** |
| *E*-2-hexenal | **2** | **\*** | **\*** | **\*** |
| *Z*-3-hexenal | **\*** | **\*** | **\*** | **\*** |
| vanillin | **0** | **0** | **0** | **0** |
| **— acetals —** | | | | |
| citral diethyl acetal | **1** | **1** | **1** | **1** |
| *E*-2-hexenal dimethyl acetal | **1** | **1** | **1** | **1** |
| **— ketones —** | | | | |
| diacetyl (butane-2,3-dione) | **1** | **1** | **1** | **1** |
| furaneol | **0** | **0** | **0** | **0** |
| 2-heptanone | **0** | **0** | **0** | **0** |
| 2-hydroxy-3-methyl-2-cyclopenten-1-one (**cyclotene**) | **0** | **1** | **1** | **1** |
| *E*-3-nonene-2-one | **1** | **1** | **1** | **1** |

### 4a. ⭐⭐ **HEXANAL'S ROW, READ CAREFULLY**
**`1 / 1 / 1 / 1`.** Hexanal forms a detectable partial adduct at ambient temperature in
30 minutes and **the score never moves** through HTST, IC pasteurization or UHT sterilization.

**What this does and does not license [Z]:**
- ✅ **It licenses:** *hexanal reacts covalently with BLG at ambient temperature at a saturating
  dose in ≤ 30 min, and no thermal treatment up to 130 °C/30 s changes that qualitative
  outcome.* Both the floor and the ceiling of the ordinal scale are untouched.
- ✅ **It licenses, weakly, an Eₐ inference:** compounds whose covalent chemistry accelerates
  strongly with T show up in this table as **0 → 1 crossings (cyclotene, 4-vinylphenol) or
  1 → \* crossings (benzaldehyde, citral, E-2-hexenal)**. **Hexanal makes neither crossing.**
  It is one of only **six of the 46 compounds** (hexanal, furfural, both acetals, diacetyl,
  3-nonen-2-one, plus the three thiols and DMDS in Table 3) that are **flat at `1` across all
  four conditions.** **Flatness across 108 K is the ordinal signature of a small activation
  energy** — and `shepelev2024_extraction.md` §6 supplies the number (15–20 kJ/mol) for a
  homologous n-alkanal on the same protein in the same lab.
- ❌ **It does NOT license any quantitative Eₐ.** The `1` bin is unbounded (§2 [!]). A tenfold
  rate change is invisible inside it.

**⇒ Register `hexanal_BLG_covalent_T_dependence: NO_ORDINAL_CHANGE_22_to_130C
[M, ordinal, non-quantitative]`.**

### 4b. **The aldehyde-class reading, compound by compound [Z]**

| tier | compounds | evidence |
|---|---|---|
| **Protein-precipitating at every severity** | ***Z*-3-hexenal** | `*` at **all four** conditions — the only compound in the paper that coagulates BLG even at 22 °C/30 min |
| **Precipitating on any heating** | **benzaldehyde**, ***E*-2-hexenal** | score at ambient (1 and **2**) then `*` at all three heated conditions |
| **Precipitating on the longer/hotter treatments** | **citral (geranial)**, ***E,E*-2,4-heptadienal** | `*` at 63 and 130 °C; heptadienal anomalously `*` at ambient but `1` at HTST — ⚠️ **[!] a non-monotone cell, see §7** |
| **⭐ Flat partial reaction, no T dependence** | **hexanal**, **furfural** | `1,1,1,1` |
| **Inert at every severity** | **vanillin** | `0,0,0,0` |

**⭐ The α,β-unsaturation ladder is visible and it matches the rest of the corpus.**
Hexanal `1,1,1,1` vs *E*-2-hexenal `2,*,*,*` vs *Z*-3-hexenal `*,*,*,*`. **Same carbon number,
same lab, same dose, same protein — the only variable is the position/presence of the C=C.**
This is the ordinal counterpart of `meynier2004_extraction.md` §5d (t-2-hexenal consumes
2.0–2.2× more His than Lys and 7× more Lys than hexanal) and of
`k2_matrix_and_thresholds.md` §D.3 (the 2–3× alkenal matrix penalty). **Three independent
methods, one conclusion: Michael-acceptor unsaturation, not chain length or logP, is what
makes an aldehyde a large covalent sink.**

### 4c. **Vanillin — the paper's most interesting negative, and the authors flag it themselves [M][Q]**
> *"The inability to measure any adduct formation of vanillin with BLG is **surprising** based on
> previous authors having proposed this reaction²⁸⁻³² and that **a color change (assumed to be
> the Maillard reaction) is often noted on longer-term storage of dry, vanillin-containing
> flavorings (personal observation)**. Presumably, the **reduced electrophilic character of the
> vanillin carboxaldehyde arising from electron donation from the neutral phenol (and,
> especially so, in the conjugate base: **vanillin pKₐ = 7.4**) is responsible for this lack of
> reactivity."*

**⇒ A named, measured, four-condition NULL for vanillin + BLG that contradicts five cited
papers (refs. 28–32: Hansen & Heinis 1991; O'Neill 1996; Li, Grun & Fernando 2000;
McNeill & Schmidt 1993; Ng, Hoehn & Bushuk 1989). Register
`vanillin_BLG_covalent: NULL_at_all_four_conditions [M]` with the contradicted citation list
attached.** **[Z] Note the pKₐ mechanism makes the null pH-dependent by construction: at
pH < 7.4 the neutral phenol dominates and the null should hold; above it the conjugate base
deactivates the carbonyl further. Since this paper's pH is unstated (§3b), the null is not
transferable to an acidic matrix without care.**

### 4d. **Acetals — an unremarked result worth registering [Z]**
**Both acetals score `1` at all four conditions**, identical to their parent aldehydes'
best-behaved case. A dimethyl/diethyl acetal is a **protected** carbonyl; for it to adduct at
all, it must hydrolyse back to the aldehyde first. **⇒ The table records, without comment, that
acetal hydrolysis is fast enough at 22 °C in 30 min in unbuffered water to feed the covalent
channel.** Relevant to any repo module treating acetals as inert carbonyl reservoirs.

---

## 5. THE THREE NEWLY REACTIVE COMPOUNDS — and one is self-falsified

### 5a. ⭐ **Eugenol — the paper disproves its own positive [M][Q]**
Table 4 gives eugenol `1, 1, 1, 1` with footnote *a*. The mass shift observed was
**+178 Da**, not the +164 Da eugenol's molecular weight demands.
> *"We were puzzled by the fact that the new adduct appearing in the incubation with eugenol
> reaction was **178 rather than 164 amu** higher than that of BLG. It occurred to us that this
> might indicate that **a compound other than eugenol itself was responsible**, especially given
> that **there is no obvious functionality in eugenol itself that should be reactive with
> protein functional groups**. One compound that came to mind was **a vinyl ketone that can
> arise from autoxidation of eugenol** … **GCMS analysis of our eugenol revealed the presence of
> a trace amount of that enone** … the **limited amount of adductation with the BLG is most
> likely a reflection of simply exhausting the small supply of the enone** rather than an
> indication of its inherent reactivity with the protein."*

**⇒ EUGENOL IS NOT REACTIVE. A trace autoxidation enone contaminant is. The paper says so and
supplies the GC-MS evidence (Figure S14).** The adduct mass is +178, the enone's.
**Register `eugenol_BLG: FALSE_POSITIVE_contaminant_enone, NOT_eugenol [M]`.**
**⭐ This is a first-class methodological warning for the whole corpus: at a 12 ppth saturating
dose, a sub-percent reactive impurity in a nominally inert compound produces a full positive
on a presence/absence assay. Any `0 → 1` transition in this table is only as clean as the
reagent.** *(Mass check [Z]: eugenol C₁₀H₁₂O₂ = 164.20; the vinyl-ketone autoxidation product
is 178, i.e. +14 = one more oxygen minus two hydrogens — consistent with an allylic oxidation
to the conjugated enone. The arithmetic supports the authors' assignment.)*

### 5b. **4-Vinylphenol (4VL, 120 Da) [M][Q]** — Table 4: `0, 0, 1, 1`
*"Slow, progressive appearance of increased levels of a **1:1 addition adduct**."* Mechanism
proposed (Figure 5): the **para-quinone methide tautomer**, *"although thermodynamically less
stable than the aromatic 4VL isomer, this quinone methide is sufficiently kinetically stable
(it has been isolated and structurally characterized²⁴) that, once formed, even in low
equilibrium concentration, it can be expected to have a sufficiently long lifetime to engage
BLG through a **conjugate addition of nucleophilic lysine (amine) and/or cysteine (thiol)**."*
Precedent cited: 4VL adds methanol to a benzylic methyl ether (ref. 25).
**⭐ Note the pattern: `0` at HTST (72 °C, 15 s) but `1` at IC (63 °C, 30 min) — the LOWER
temperature with the LONGER time. Time, not temperature. §6.**

### 5c. **3-Nonen-2-one (140 Da) [M][Q]** — Table 2: `1, 1, 1, 1`
*"observed to form a covalent **1:1 addition adduct** having 140 amu higher than BLG itself …
by **conjugate addition of a lysine ε-amino and/or cysteine thiol group to the β-carbon (C4)
of the conjugated enone** (i.e., via **aza-Michael²⁶ or thia-Michael²⁷** reactions, Figure 7).
3-Nonen-2-one is sufficiently reactive that, although **unreacted BLG … is found in the
ambient and HTST-heated samples**, **no unreacted BLG remains in either of the more forcefully
heated samples** (i.e., IC pasteurization and UHT)."*
⚠️ **[!] The prose says the IC and UHT samples show NO unreacted BLG — which by the paper's own
definition is a score of 2. Table 2 nevertheless prints `1` for both.** A direct
text-vs-table contradiction; **prefer the prose, and register `3-nonen-2-one: 1,1,2,2 (prose)`
alongside the printed row.** *(This also shows again that the enone, not the parent, is the
Michael acceptor — the same motif as §5a.)*

---

## 6. ⭐⭐ THE SEVERITY ORDERING — the paper's one quasi-quantitative statement, and what may be done with it

### 6a. The statement, verbatim from the abstract **[M][Q]**
> *"An overall view of the data shows that the **HTST heat treatment (72 °C for 15 s) had the
> least effect on the extent of reaction** while **in-container pasteurization conditions
> (63 °C for 30 min) produced a similar extent of reaction as the UHT (130 °C 30 s)** heat
> treatment. These varying extents of adductation are **in reasonable accord with what one might
> expect, given that the rates of most classes of chemical reactions occurring near ambient
> temperature increase by a factor of 2−4 for each increase of 10 K in temperature.**"*

And in Results, for eugenol specifically: *"the mass spectra of the IC pasteurization (63 °C
for 30 min) and the UHT (130 °C for 30 s) sterilization processes (Figure 2b,d) **look very
similar** while the HTST sample (Figure 2c) **shows less eugenol adductation** than either."*

⚠️ **[!] THE Q₁₀ = 2–4 IS A TEXTBOOK GENERALITY, NOT A MEASUREMENT.** *"what one might expect,
given that the rates of most classes of chemical reactions…"* — no rate was measured, no Q₁₀
was computed, and the sentence is doing rhetorical, not quantitative, work.
**Do not ingest `Q10 = 2-4` as a finding of this paper.** *(For reference [Z], Q₁₀ = 2–4 near
300 K corresponds to Eₐ = 53–99 kJ/mol; `shepelev2024_extraction.md` measures **Q₁₀ ≈ 1.2** for
decanal + BLG, i.e. **outside and far below** the generality this sentence invokes.)*

### 6b. **[Z] Converting "IC ≈ UHT" into an Arrhenius bracket — the arithmetic, and then the reasons not to trust it**
If two treatments produce equal extent, their **time-integrated rates are equal**:
`k(130 °C) × 30 s = k(63 °C) × 1800 s` ⇒ `k(130)/k(63) = 60`.

> `Eₐ = R ln(60) / (1/336.15 − 1/403.15) = 8.314 × 4.094 / 4.944 × 10⁻⁴` = **68.9 kJ/mol** [Z]

**⚠️ FOUR REASONS THIS NUMBER IS NOT AN ALDEHYDE ACTIVATION ENERGY, and the corpus must not
oppose it to `shepelev2024`'s 15–20:**

1. **⭐ It is a statement about the WHOLE 46-compound panel, and the panel's response is
   dominated by the compounds that actually move** — isothiocyanates, thiols, activated enones,
   quinone methides. **The compounds that are flat (hexanal, furfural, diacetyl, the acetals,
   the thiols) contribute nothing to a "similar extent" impression, because they look the same
   everywhere.** An eyeball comparison of mixed spectra is weighted by the movers.
   *(Consistently: `shepelev2024_extraction.md` §7c measures **Eₐ = 36–43 kJ/mol** for an
   isothiocyanate on the same protein — an order-of-magnitude-appropriate value for a panel
   whose most reactive members are isothiocyanates.)*
2. **"Similar" is an eyeball comparison of two mass spectra** (Fig. 2b vs 2d), not a
   measurement, and it is asserted for **eugenol** — the one compound the paper itself shows is
   reacting via a **contaminant present in limited supply** (§5a). **The comparison compound is
   the worst possible choice: its extent is limited by exhaustion of the reactant, not by rate.**
3. **The 130 °C extent is biased DOWNWARD by coagulation.** Five aldehydes and the two most
   reactive sulfur compounds lose their protein at the forcing conditions (§3c). Removing
   adducted protein from the soluble phase makes 130 °C look *less* reactive than it is, which
   makes 68.9 kJ/mol a **lower** bound on whatever the panel's effective Eₐ is — pushing the
   number the *wrong* way for reconciliation with Shepelev, and confirming that it is measuring
   something else entirely.
4. **The 130 °C sample was in a 1.5 mm capillary in an oil bath with 6 s equilibration
   (i.e. 20 % of its 30 s hold), the 63 °C sample in a bulk vial in a water bath with 60 s
   equilibration (3 % of its 30 min hold).** The effective thermal histories are not the clean
   step functions the arithmetic assumes.

**⇒ REGISTER: `yuan2023_panel_Ea_inference: 68.9 kJ/mol [Z], class: mixed_panel_ordinal_eyeball,
DO NOT USE AS AN ALDEHYDE Ea, DO NOT OPPOSE TO shepelev2024`. Its legitimate use is as an
order-of-magnitude sanity check on the ISOTHIOCYANATE/thiol/enone end of the ladder, where it
sits within 1.6–1.9× of Shepelev's measured PEITC value.**

### 6c. **[Z] The one place the table genuinely does bound an Eₐ from below — and it is not an aldehyde**
**Cyclotene** goes `0 (22 °C, 30 min) → 1 (72 °C, 15 s)`. A 0 → 1 crossing means the extent
rose above the detection floor despite the exposure time falling 120×:
`k(72)/k(295 K) ≥ 120` ⇒ **Eₐ(cyclotene) ≥ 81 kJ/mol** [Z].
**4-Vinylphenol** goes `0 (72 °C, 15 s) → 1 (63 °C, 30 min)`, i.e. time beat temperature:
this bounds Eₐ from **above**: `k(72)/k(63) ≤ 120` ⇒ **Eₐ(4VL) ≤ 513 kJ/mol** — vacuous.
**⇒ The table yields exactly one non-trivial Eₐ bound, it is a lower bound, it belongs to a
cyclic enolone, and it says nothing about alkanals.** *(It does, usefully, show the ordinal
scale is capable of registering a large Eₐ when one is present — which makes hexanal's
flatness informative rather than merely uninformative.)*

---

## 7. TABLE 3 — SULFUR COMPOUNDS — TRANSCRIBED IN FULL **[M]**

| sulfur compound | **ambient (30 min)** | **72 °C (15 s)** | **63 °C (30 min)** | **130 °C (30 s)** |
|---|:---:|:---:|:---:|:---:|
| 4-methyl-5-vinylthiazole | 0 | 0 | 0 | 0 |
| 2-methylthiophene | 0 | 0 | 0 | 0 |
| dimethyl sulfone | 0 | 0 | 0 | 0 |
| **thiophenol** | **1** | **1** | **1** | **1** |
| **propanethiol** | **1** | **1** | **1** | **1** |
| **2-furfurylmercaptan** | **1** | **1** | **1** | **1** |
| dimethyl sulfide | 0 | 0 | 0 | 0 |
| **dimethyl disulfide** | **1** | **1** | **1** | **1** |
| **dimethyl trisulfide** | **1** | **1** | **1** | **2** |
| ⭐ **allyl isothiocyanate** | **1** | **1** | **2** | **2** |

**[Z] Two observations the paper does not draw:**
1. **The three free thiols and DMDS are as flat as hexanal (`1,1,1,1`).** Combined with
   `hidalgo2010_extraction.md`'s **measured Eₐ = 28–30 kJ/mol for thiol-Michael addition**, the
   flat-row/small-Eₐ association now has an independent quantitative anchor in a second
   laboratory and a second chemistry. **Flat row ⇔ small Eₐ is not a guess.**
2. **Only two compounds in the entire paper reach a score of 2 on heating: dimethyl trisulfide
   (at 130 °C only) and allyl isothiocyanate (at 63 and 130 °C).** These are the panel's movers,
   and §6b's 68.9 kJ/mol is essentially their Eₐ, not the aldehydes'.
   ⭐ **Note allyl isothiocyanate scores 2 at IC (63 °C, 30 min) as well as UHT — the direct
   evidence behind the "IC ≈ UHT" claim, and it comes from an isothiocyanate.**

Author commentary [Q][C]: *"**allyl (iso)thiocyanate and dimethyl trisulfide were the most
reactive of all flavor compounds studied.** … The thiocyanates have the ability to react with
**both amines and sulfur-groups** in a protein."* Keppler et al. (refs. 10, 11, 12, 33) cited as
the reference work on AITC + BLG/WPI, with the note that *"the focus of Keppler's studies has
been on the **effect of this reaction on the functionality of the protein** rather than the
effect of the reaction on the **flavor/sensory character**."*

---

## 8. TABLE 4 — NITROGEN AND MISCELLANEOUS COMPOUNDS — TRANSCRIBED IN FULL **[M]**

| compound | **ambient (30 min)** | **72 °C (15 s)** | **63 °C (30 min)** | **130 °C (30 s)** |
|---|:---:|:---:|:---:|:---:|
| **2,5-dimethylpyrazine** | **0** | **0** | **0** | **0** |
| **3-acetylpyridine** | **0** | **0** | **0** | **0** |
| **4-vinylphenol** | **0** | **0** | **1** | **1** |
| **eugenol** ᵃ | **1** | **1** | **1** | **1** |

ᵃ *"See the text (and Figure 3) for hypothesis of what may limit the further consumption and
adductation of BLG by eugenol"* — i.e. the contaminant story, §5a.

**⭐ A MEASURED FOUR-CONDITION NULL FOR PYRAZINES, THROUGH UHT.** *"As one would expect, the
pyridine and pyrazine were unreactive."* **Register `pyrazine_protein_covalent_sink: NULL
through 130 °C/30 s [M]`. Any repo module that puts a protein sink on the pyrazine lane is
contradicted by direct measurement.**

---

## 9. TABLE 1 — THE 46-COMPOUND PANEL BY FUNCTIONAL-GROUP CLASS — TRANSCRIBED IN FULL **[M]**

| class | compounds as printed |
|---|---|
| hydrocarbons | p-cymene |
| alcohols | L-menthol · geraniol · 2-pentanol · 2,3-butanediol |
| phenols | **eugenol ᵃ** · 2,6-dimethylphenol · **4-vinylphenol ᵃ** |
| enols | furaneol |
| acids | butyric acid · **4-hydroxybenzoic acid ᵇ** · m-toluic acid |
| esters | methyl anthranilate · iso-amyl acetate · ethyl lactate · methyl salicylate |
| lactones | δ-dodecalactone · γ-butyrolactone · γ-decalactone |
| ketones | diacetyl · 2-heptanone · **3-nonen-2-one ᵃ** · cyclotene |
| **aldehydes** | trans-2-hexenal · cis-3-hexenal · trans,trans-2,4-heptadienal · citral (geranial) · vanillin · **hexanal** · furfural · benzaldehyde |
| acetals | citral diethyl acetal · trans-2-hexenal dimethyl acetal |
| sulfur-containing | 4-methyl-5-vinyl-thiazole · 2-methyl thiophene · dimethyl sulfone · propanethiol · thiophenol · 2-furfuryl mercaptan · dimethyl sulfide · dimethyl disulfide · dimethyl trisulfide · allyl isothiocyanate |
| pyridines | 3-acetylpyridine · **2,5-dimethylpyrazine** ⚠️ *(printed under "pyridines" — a pyrazine is not a pyridine; a table-layout error [!])* |
| pyranone | maltol |

ᵃ *"Newly found to be reactive on heating."* ᵇ *"Substituted for 2-hydroxybenzoic acid
(salicylic acid) used in the previous study."*

**⚠️ [!] COUNT CHECK [Z]: Table 1 as printed lists **44** compounds; the abstract and title
claim **46** across **13 functional-group classes**. Table 1 shows **13 class headings** ✔ but
two compounds short. The Results text additionally names **ethanol** among the unreactive
alcohols, which is absent from Table 1 — accounting for one. **The 46th is unaccounted for.**
Register `panel_size: 44_printed_46_claimed`.**

### 9a. **The comprehensive unreactive list, verbatim from Results [M]** — a large measured null set
> *"Many classes of compounds showed **no measurable reactivity with BLG, irrespective of the
> method of heat treatment**: **alcohols (ethanol, menthol, and 2,3-butanediol)**, **phenols
> (2,6-dimethylphenol)**, **enols (maltol)**, **acids (butyric, hydroxybenzoic, and toluic)**,
> **alkenes (cymene and geraniol)**, **esters (isoamyl acetate, ethyl lactate, and methyl
> salicylate)**, **lactones (δ-dodecalactone, γ-butyrolactone, and γ-decalactone)**, **ketones
> (2-heptanone and 3-acetylpyridine)**, and several miscellaneous compounds (**methyl
> anthranilate, vanillin, 2-methylthiophene, and dimethyl sulfone**)."*

⚠️ **[!] Two mis-classifications inside that sentence:** **geraniol is an alcohol, not an
alkene**, and **3-acetylpyridine is a pyridine, not a ketone** (it does have a ketone group, but
it is listed under pyridines in Table 1 and Table 4). Cosmetic, but they will corrupt a
class-level ingest.

**⇒ Register the whole list as `covalent_protein_sink: NULL through UHT [M]` — 19 named
compounds spanning 9 classes. For a repo that must decide which volatiles need a protein sink
term, this negative inventory is worth as much as the positives.**

---

## 10. WHAT THE PAPER CLAIMS IT ESTABLISHES — verbatim, for the record **[Q]**

> *"As previously shown, **all aldehydes, unsaturated ketones or diketones, and most sulfur
> compounds (thiols, disulfides, trisulfides, and isothiocyanates) readily reacted with BLG
> during thermal processing.** The three newly noted reactants, eugenol (**or, most likely, a
> contaminant therein**), 4-vinyl phenol, and 3-nonen-2one were also reactive. It is of interest
> that the **levels of reaction under conditions mimicking IC pasteurization (63 °C for 30 min)
> were similar to the UHT process (130 °C for 30 s)**. This work clearly shows that several
> classes of flavor molecules, previously unreactive at ambient temperatures, will react with
> BLG (and likely many other proteins) under more forcing thermal conditions."*

And the environmental caveats it carries forward from ref. 18 [C]:
> *"reaction rates have been observed to be **lower at pH 3 and higher at pH 7 or pH 8**. There
> were **no observable differences in the rate of adduct formation as a function of water
> activity** for three of the compounds studied that were found to be very reactive, i.e.,
> **benzaldehyde, citral, and dimethyl disulfide**. The observation that water activity had no
> effect on the extent of reaction was considered to be due to the low levels of reaction due
> to **short reaction times (generally less than 30 min) and/or low reaction temperatures
> (≤ 45 °C)**."*
**⇒ The a_w null is explicitly self-limited by the authors and must not be transferred to a
low-moisture process. The pH direction (3 ≪ 7 ≈ 8) is consistent with
`meynier2004_extraction.md` §6e and with `hidalgo1993_extraction.md`'s measured base catalysis.**

---

## 11. SUPPORTING INFORMATION — what it holds and whether it is needed

Contents as listed on p. 9489:
1. *"Deconvoluted LC−ESI−MS spectra of the products formed between **one representative compound
   from each class of functional group**"* — **13 spectra, one per class.** **⚠️ Hexanal is
   almost certainly NOT the aldehyde representative** (the paper's Figures 2/4/6 feature the
   three newly reactive compounds; the aldehyde exemplar is unstated).
2. *"GC−EI−MS analysis of eugenol showing a trace amount of the enone"* — **Figure S14**, the
   evidence behind §5a.
3. *"the **reason for the choice of each specific flavor compound** within any one class"* —
   **Table S1**, reproduced from ref. 17 (`anantharamkrishnan2020b_extraction.md` §, which
   already transcribes it).

**⇒ SI PRIORITY: LOW.** It contains **no additional quantitative data, no time courses, no
replicates and no per-compound extents.** The only item of independent value is **Figure S14**,
and only if the eugenol false-positive is ever disputed. **Register `yuan2023_SI: not_required`.**

---

## 12. INTEGRITY FLAGS RAISED BY THIS WAVE — consolidated

| # | flag | where | severity |
|---|---|---|---|
| **1** | **`*` means the reaction went so far the protein precipitated — it is `≥ 2`, NOT a missing value** | §3c | ⭐ **HIGH — inverts the T trend for 5 aldehydes if mishandled** |
| **2** | **Column header reads "reaction rate"; the footnote defines an EXTENT score. No rate is measured anywhere.** | §4 | ⭐ **HIGH — a laundering hazard by wording alone** |
| **3** | **Eugenol's positive is a contaminant enone (+178 Da ≠ +164 Da), as the authors themselves demonstrate** | §5a | ⭐ HIGH |
| 4 | The `1` bin is unbounded; a 10× rate change is invisible inside it | §2 | HIGH |
| 5 | 110 °C/30 min produced no data for any compound (protein lost) | §3b | MEDIUM |
| 6 | Text says no unreacted BLG remains for 3-nonen-2-one at IC and UHT (⇒ score 2); Table 2 prints 1 | §5c | MEDIUM |
| 7 | BLG mass given as 18,363 / 18,362 / "18.2 kDa" in three places | §3a | LOW |
| 8 | Table 1 lists 44 compounds; abstract claims 46 | §9 | MEDIUM |
| 9 | 2,5-dimethylpyrazine filed under "pyridines"; geraniol called an alkene; 3-acetylpyridine called a ketone | §9, §9a | LOW (class-ingest corruption) |
| 10 | E,E-2,4-heptadienal is `*` at ambient but `1` at HTST — a non-monotone cell with no explanation | §4b | MEDIUM |
| 11 | pH never stated, unbuffered; a_w null explicitly self-limited by the authors | §3b, §10 | MEDIUM |
| 12 | "rates increase by a factor of 2−4 per 10 K" is a textbook generality, **not** a measurement of this paper | §6a | ⭐ HIGH (would be laundered as a measured Q₁₀) |
| 13 | Aldehyde-reactivity claim here (aliphatic > aromatic, measured) contradicts `shepelev2024`'s cited claim (aromatic > aliphatic, from 1977 haemoglobin work) | `shepelev2024_extraction.md` §8c | LOW–MEDIUM |
