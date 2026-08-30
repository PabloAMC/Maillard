# Wave K6b — ADDUCT KINETICS SYNTHESIS
### ⭐⭐⭐ **THE COVALENT-SINK VERDICT: the bracket `meynier2004` left open is CLOSED, and it closes on the "does not matter at process temperature" side.**
### Plus: the milk-unit question, SETTLED by an exact seven-row numeric anchor.

**Wave K6b, 2026-08-30.** Eight papers extracted; seven new per-paper dossiers written; this
synthesis consolidates them. **Read-only outside `data/lit/extraction_dossiers/`. No repo file
was modified, no declaration was edited, nothing was committed.**

**Dossiers written by this wave** — all new files, none pre-existing:
`shepelev2024_extraction.md` · `yuan2023_extraction.md` · `zamora2013_extraction.md` ·
`hidalgo1993_extraction.md` · `hidalgo2010_extraction.md` · `zamora2010_extraction.md` ·
`tian2019_extraction.md` · `tian2020b_extraction.md` · **this file.**

**Cross-referenced, not re-extracted:** `hamzalioglu2018_extraction.md` (K5a),
`meynier2004_extraction.md` (K4e), `anantharamkrishnan2020b_extraction.md` (K4d),
`k2_matrix_and_thresholds.md`, `k3_final_parameter_inventory.md`, `research_round4_nulls.md`.

---

# §0. ⭐ FILE IDENTITY — **ALL EIGHT PDFs MATCH THEIR EXPECTED IDENTITY. ZERO WRONG-FILE INCIDENTS.**

*Reported first and prominently, as the brief requires.*

| file | expected DOI | **found** | verdict |
|---|---|---|---|
| `Shepelev2024.pdf` | 10.1021/acs.jafc.4c00134 | **10.1021/acs.jafc.4c00134** — Shepelev & Reineccius, *JAFC* 72:10579–10583 | ✅ **exact** |
| `Yuan2023.pdf` | 10.1021/acs.jafc.3c01220 | **10.1021/acs.jafc.3c01220** — Yuan, Anantharamkrishnan, Hoye & Reineccius, *JAFC* 71:9481–9489 | ✅ **exact** |
| `zamora2013.pdf` | 10.1021/jf305007y | **10.1021/jf305007y** — Zamora, Alcón & Hidalgo, *JAFC* 61:10231–10237 | ✅ **exact** |
| `hidalgo1993.pdf` | 10.1111/j.1365-2621.1993.tb04352.x | **that DOI's paper** — Hidalgo & Zamora, *J. Food Sci.* 58(3):667–670 | ✅ **exact** |
| `hidalgo2010.pdf` | 10.1016/j.foodchem.2010.03.016 | **10.1016/j.foodchem.2010.03.016** — Hidalgo, Delgado & Zamora, *Food Chem.* 122:596–601 | ✅ **exact** |
| `zamora2010.pdf` | 10.1021/jf903378x | **10.1021/jf903378x** — Zamora, Delgado & Hidalgo, *JAFC* 58:1708–1713 | ✅ **exact** |
| `tian2020b.pdf` | 10.3168/jds.2019-17495 | **10.3168/jds.2019-17495** — Tian, Yu, Yu & Chen, *J. Dairy Sci.* 103:7957–7967 | ✅ **exact** |
| `tian2019.pdf` | 10.3168/jds.2019-16796 | **10.3168/jds.2019-16796** — Tian, Xu, Chen & Yu, *J. Dairy Sci.* 102 | ✅ **exact** |

**⇒ NO MIS-FILES. NO WRONG IDENTITIES. NOTHING TO REPORT UNDER THIS HEADING.**

**Four naming/labelling hazards were found INSIDE correctly-filed papers, and they are worse
than a mis-file would have been because they survive an identity check:**

| # | hazard | consequence if unflagged |
|---|---|---|
| **H1 ⭐⭐** | **`Shepelev2024`: the Abbreviations list defines `DEC = decane` (the inert blank), but Figure 3, the Results heading and the Summary all use `DEC` to mean DECANAL.** | **The entire aldehyde dataset gets assigned to an inert alkane control.** |
| **H2 ⭐⭐** | **`Yuan2023`: a `*` cell means the reaction went so far the protein COAGULATED — chemically `≥ 2`, not a missing value.** | **Inverts the temperature trend for five aldehydes.** |
| **H3 ⭐** | **`Yuan2023`: every table column is headed "reaction rate"; the footnote defines a 0–2 EXTENT score. No rate exists in the paper.** | A pure laundering hazard by wording. |
| **H4 ⭐** | **`hidalgo1993`: the OCR text layer of the 1993 scan renders Table 1's pH-6.5 rate constant as `6.12` where the paper prints `0.12` — a 51× error that inverts the pH trend**; and the abstract's `50 kJ/mol` as `SOKJ/mol`. | **Any automated ingest of that PDF's text layer is unsafe.** Wave K6b re-read the table at 400 dpi. |

---

# §1. ⭐⭐⭐ THE CONSOLIDATED PARAMETER TABLE

## 1a. THE ACTIVATION-ENERGY LADDER — every measured Eₐ now in the corpus, sorted

*Provenance: **[M]** printed by the paper · **[Z]** derived by wave K6b · **[F]** printed, refit
by an earlier wave · **[C]** cited by a corpus paper from a source not held.*

| **Eₐ (kJ/mol)** | reaction | nucleophile | electrophile | T range | n_T | source | prov |
|---:|---|---|---|---|---:|---|---|
| **3.3** | HMF + arginine, coffee (low moisture) | Arg | aldehyde | 5–50 °C | 3 | `hamzalioglu2018` | **[F] ⛔ REFUSED — chemically impossible barrier** |
| **10.0** | ⭐ **HMF + LYSINE, aqueous pH 3.5** | **Lys** | **aldehyde** | **5–50 °C** | 3 | `hamzalioglu2018` | **[F]** |
| **11.6** | HMF + cysteine, coffee | Cys | aldehyde | 5–50 °C | 3 | `hamzalioglu2018` | [F] |
| **12.3** | HMF + lysine, coffee | Lys | aldehyde | 5–50 °C | 3 | `hamzalioglu2018` | [F] |
| **15–20** | ⭐⭐⭐ **DECANAL + β-LACTOGLOBULIN covalent adduct** | **protein Lys** | **n-alkanal** | **25–65 °C** | **3** | **`shepelev2024` Fig. 3** | **[Z] K6b** |
| **24.4** | acrylamide + cysteine | Cys | activated alkene | — | — | Hidalgo 2011, via `hamzalioglu2018` | [C] |
| **25.5** | HMF + arginine, aqueous | Arg | aldehyde | 5–50 °C | 3 | `hamzalioglu2018` | [F] |
| **27.6** | 2,4-decadienal + Phe → phenylacetaldehyde, **in AIR** | Phe | 2,4-alkadienal | 120–180 °C | 6 | `zamora2013` | [M] |
| **28** | acrylamide + **benzyl mercaptan** (thia-Michael) | thiol | activated alkene | **80–180 °C** | **8** | `hidalgo2010` | [M] |
| **29.6** | HMF + cysteine, aqueous | Cys | aldehyde | 5–50 °C | 3 | `hamzalioglu2018` | [F] |
| **30** | acrylamide + **N-acetylcysteine** | thiol | activated alkene | 80–180 °C | 8 | `hidalgo2010` | [M] |
| **36–43** | ⭐ **PEITC + β-lactoglobulin covalent adduct** | protein Lys/Cys | isothiocyanate | 25–65 °C | 3 | `shepelev2024` Fig. 4 | **[Z] K6b** |
| **37.5** | 2,4-heptadienal + Phe → PAC | Phe | 2,4-alkadienal | 120–180 °C | 6 | `zamora2013` | [M] |
| **44** | ⭐ acrylamide + butylamine, **FORWARD** (aza-Michael) | amine | activated alkene | 120–180 °C | 4 | `zamora2010` | [M] |
| **50** | epoxyheptenal + lysine → **FLUORESCENCE** | Lys | epoxyalkenal | **4–80 °C** | 4 | `hidalgo1993` | [M] |
| **52** | acrylamide + glycine, forward | amino acid | activated alkene | 120–180 °C | 4 | `zamora2010` | [M] |
| **52** | ⭐ 3-(butylamino)propanamide → acrylamide, **REVERSE** | — | — | 120–180 °C | 4 | `zamora2010` | [M] |
| **55.3** | 4-oxo-2-hexenal + Phe → PAC | Phe | 4-oxo-2-alkenal | 120–180 °C | 4 | `zamora2013` | [M] |
| **57.9** | 4,5-epoxy-2-heptenal + Phe → PAC | Phe | epoxyalkenal | 120–190 °C | 7 | `zamora2013` | [M] |
| **58.9** | 4,5-epoxy-2-decenal + Phe → PAC | Phe | epoxyalkenal | 120–190 °C | 7 | `zamora2013` | [M] |
| **63.5** | 4-oxo-2-nonenal + Phe → PAC | Phe | 4-oxo-2-alkenal | 120–180 °C | 4 | `zamora2013` | [M] |
| **66.5** | epoxyheptenal + lysine → **COLOUR (ΔE)** | Lys | epoxyalkenal | 4–80 °C | 4 | `hidalgo1993` | [M] |
| **67.1** | 4-hydroxy-2-nonenal + Phe → PAC | Phe | 4-hydroxyalkenal | 120–190 °C | 7 | `zamora2013` | [M] |
| **≈69** | ⚠️ Yuan's mixed-panel ordinal "IC ≈ UHT" inference | 46 compounds | mixed | 63 vs 130 °C | 2 | `yuan2023` | **[Z] ⛔ NOT AN ALDEHYDE Eₐ — see §2c** |
| **78.0** | ⭐ 2,4-decadienal + Phe → PAC, **under NITROGEN** | Phe | 2,4-alkadienal | 120–180 °C | 6 | `zamora2013` | [M] |
| **≥ 81** | ⚠️ cyclotene + BLG, ordinal **lower bound only** | protein | cyclic enolone | 22 vs 72 °C | 2 | `yuan2023` | [Z] |
| **81.0** | asparagine decarboxylation by 2,4-decadienal | Asn | 2,4-alkadienal | — | — | Hidalgo 2010 *JAFC* | [C] |
| *(100–133)* | *milk colour change on heating* | — | — | — | — | Pagliarini 1990, via `hidalgo1993` | [C] *implied* |

## 1b. ⭐⭐ THE LADDER HAS THREE BANDS, AND THE BAND BOUNDARY IS THE WHOLE ARGUMENT

| band | **Eₐ** | what the reaction IS | members |
|---|---:|---|---|
| **A — DIRECT ADDUCT FORMATION** | **10–30** | a nucleophile adds to a carbonyl or an activated alkene; the observable is the **adduct** or the **electrophile's disappearance** | **HMF + Lys 10.0** · **decanal + BLG 15–20** · HMF + Arg 25.5 · HMF + Cys 29.6 · acrylamide + thiol 28–30 · acrylamide + Cys 24.4 |
| **B — HINDERED ADDITION AND ITS REVERSE** | **36–52** | addition at a less accessible or less activated site, and the elimination back out | PEITC + BLG 36–43 · acrylamide + butylamine **44 fwd** · acrylamide + glycine 52 · **adduct → acrylamide 52 rev** |
| **C — DOWNSTREAM CASCADES** | **50–81** | the observable is a **product several steps past the adduct**: a Strecker aldehyde, a colour, a fluorophore, a biogenic amine | Strecker 27.6–67.1 · **browning 66.5** · **fluorescence 50** · decadienal/N₂ 78.0 · Asn decarboxylation 81.0 |

**⇒ ⭐⭐⭐ THE FINDING THAT REFRAMES THE WHOLE QUESTION.
Adduct FORMATION sits in band A, at 10–30 kJ/mol. Every intuition the corpus had about
"carbonyl–amine activation energies" was built on band C — Strecker chemistry, browning,
fluorescence — at 50–81 kJ/mol. The `≥ 70 kJ/mol` threshold that
`anantharamkrishnan2020b_extraction.md` set for the covalent sink to matter at process
temperature is a BAND C number applied to a BAND A reaction.** That is the single sentence
this wave adds to the corpus.

**And band A is now populated by two independent laboratories, two different aldehydes, two
different nucleophile presentations and two overlapping-but-distinct temperature ranges:**

| | `hamzalioglu2018` (Hacettepe) | `shepelev2024` (Minnesota) |
|---|---|---|
| electrophile | **HMF** (an aromatic aldehyde) | **decanal** (an aliphatic n-alkanal) |
| nucleophile | **free L-lysine** | **lysine ε-NH₂ on β-lactoglobulin** |
| temperatures | **5, 25, 50 °C** | **25, 45, 65 °C** |
| method | HPLC-UV loss of HMF, pseudo-first-order | **¹⁴C scintillation counting of protein-bound label** |
| matrix | aqueous, pH 3.5, buffered | aqueous, unbuffered |
| **Eₐ** | **10.0 kJ/mol** | **15–20 kJ/mol** |

**Two labs that have never cited each other, two entirely unrelated analytical methods, one
answer: the barrier to aldehyde–amine adduct formation is 10–20 kJ/mol.**

## 1c. THE COVALENT-SINK RATE CONSTANTS — consolidated

| system | **k₂ (M⁻¹ s⁻¹), ambient** | **t½ in a 3 wt % protein matrix** | source |
|---|---:|---|---|
| **decanal + BLG, 25 °C** | **≈ 2.9 × 10⁻⁴** | *(this experiment is protein-limited; see below)* | **[Z] K6b, `shepelev2024` §4f** |
| t-2-hexenal + Na-CN/WP Lys, 20 °C | 5.3–7.9 × 10⁻⁵ (Lys) · 1.5–2.4 × 10⁻⁴ (His) | 2.8–5.0 d | `meynier2004` §8 |
| **hexanal + Na-CN/WP Lys, 20 °C** | **≤ 2.5 × 10⁻⁵ (upper bound)** | **≥ 18 d** | `meynier2004` §1(b) |
| hexanal + BLG, ambient (MS adduct counting) | 10⁻⁶ – 10⁻⁵ | 37–760 d | `anantharamkrishnan2020b` §6 |
| **⇒ the corpus's hexanal bracket** | **6 × 10⁻⁶ – 2.5 × 10⁻⁵** | **⭐ 18–74 d (overlap region 37–74 d)** | `meynier2004` §1(c) |

**⭐ [Z] K6b's chain-length correction, and it matters.** Decanal's k₂ is **12–15× larger** than
hexanal's upper bound. The corpus's two independently measured hydrophobic chain-length slopes
(**Andriot 2.72×/CH₂**, **Damodaran 2.9×/CH₂**, `k2` §B.6) predict an enhancement of
**3× to 61×** for four extra CH₂ groups if adduction proceeds from a hydrophobically bound
state. **The observed 12–15× is inside that bracket.**
**⇒ `k₂(decanal)` IS AN UPPER BOUND ON `k₂(hexanal)`, NOT AN ESTIMATE OF IT. `shepelev2024`
does NOT move the corpus's hexanal rate bracket; it supplies the missing Eₐ and leaves the rate
where `meynier2004` and `anantharamkrishnan2020b` put it.**

## 1d. THE NON-Eₐ PARAMETERS THIS WAVE ADDS

| parameter | value | source | note |
|---|---|---|---|
| **Non-covalent carry-over floor, ¹⁴C method** | **< 0.1 % of added dose, all T, all t** | `shepelev2024` §3d | the control `meynier2004` and `anantharamkrishnan2020b` both lacked |
| **Mono-ketone vs n-alkanal covalent selectivity** | **16–18×** (2-undecanone vs decanal, matched dose/T) | `shepelev2024` §5 [Z] | two independent routes agree |
| **Isothiocyanate vs n-alkanal extent** | **4.4×** (71.8 % vs 16.4 % of dose) | `shepelev2024` [M] | |
| **Reversible : irreversible split, storage sink** | **≈ 90 : 10** over 28 d at 37–60 °C | `zamora2010` §5b [M] | **the corpus's only measured split** |
| **K_eq temperature dependence, aza-Michael adduct** | **falls 3.0× from 25 → 180 °C** | `zamora2010` §4 [Z] | from Eₐ,rev − Eₐ,fwd = +8 |
| **Free-amine → polymer-bound reactivity penalty** | **Lys > polyLys 4.2 kDa > polyLys 68.3 kDa**, ∝ −ln(MW) | `zamora2010` §7c [M] | ⚠️ *p* = 0.095 on 3 points |
| **Oxygen effect on an alkadienal Eₐ** | **2.83×** (27.6 air → 78.0 N₂) | `zamora2013` §3b [M] | atmosphere is a required field |
| **Buffer effect at constant pH** | **2.20×** on a reverse reaction at pH 6 (citrate vs phosphate) | `zamora2010` §6a [M] | ⭐ corpus-wide caution |
| **Eₐ ↔ yield correlation** | *r* = **−0.83, p = 0.02** at 37 °C; **n.s.** at 80 °C (p = 0.4) and 180 °C (p = 0.6) | `zamora2013` §3c [M] | ⭐ an Eₐ table is not a yield table |
| **Binary carbonyl synergy in a dairy matrix** | **1.81–2.31×** threshold reduction | `tian2020b` §7a [M] | near-threshold only (5/72 mixtures) |
| **Matrix/water threshold ratios, yogurt** | diacetyl **1,810×** · acetaldehyde **140×** · acetoin **3.6×** | `tian2020b` §2 [M] | **500× spread within one class, one matrix, one panel** |
| **Fifty literature odour thresholds, μg/kg** | 0.027 – 740 | `tian2019` §4d [C] | ⚠️ **all unsourced, all mislabelled "in air"** |
| **Milk-fan matrix thresholds** | **7 rows, 160 – 52,000 μg/kg** | `tian2020` Table 1, unit resolved here | **§3** |

---

# §2. ⭐⭐⭐ THE COVALENT-SINK VERDICT

## 2a. The bracket, restated exactly as `meynier2004` left it

`meynier2004_extraction.md` §1(d) and §9c:
> *"**25 % + < 1 % ≠ 100 %.** The covalent channel only becomes material at **process
> temperature**, and its activation energy is **unmeasured in this paper and in every other
> paper in the corpus**."*

`anantharamkrishnan2020b_extraction.md`, via `research_round4_nulls.md` §A.0, made the
dependence explicit:
> *"| assumed Eₐ | rate acceleration 20 → 140 °C | t½ at 140 °C | … | **70 kJ mol⁻¹** | ≈ 4,300× |
> **≈ 65 minutes** | **80 kJ mol⁻¹** | ≈ 13,900× | **≈ 20 minutes** |"*

**⇒ The question was binary and precisely posed: is Eₐ ≳ 70 kJ/mol?**

## 2b. **THE ANSWER: NO. Eₐ = 15–20 kJ/mol, measured. The threshold is missed by 3.5–4.7×.**

**Primary measurement** — `shepelev2024_extraction.md` §6, derived by this wave from Figure 3
(the paper computes no kinetics):

| estimator | Eₐ | R² |
|---|---:|---:|
| raw day-1 extents (3.6 / 5.0 / 7.5 mg/g at 25 / 45 / 65 °C) | **15.2 kJ/mol** | 0.998 |
| **saturation-corrected first-order approach to plateau** | **20.0 kJ/mol** | 0.9997 |
| full sensitivity envelope over every digitization choice | **14 – 23 kJ/mol** | — |
| propagated from the printed SDs (n = 4) | **20 ± 4 (1 σ)** → 95 % CI **12 – 28** | — |

**Four independent corroborations, none of which shares a method with the primary:**

| # | corroboration | value | why it is independent |
|---|---|---|---|
| **1 ⭐** | **`hamzalioglu2018` HMF + lysine, aqueous, 5–50 °C** | **10.0 kJ/mol** | different lab, different aldehyde, different nucleophile presentation, different analytical method (HPLC-UV vs ¹⁴C) |
| **2 ⭐** | **`shepelev2024`'s own PEITC control**, same figures, same apparatus | **36–43 kJ/mol** | proves the method **resolves** Eₐ — a 2.0–2.5× separation within one experiment. The low decanal number is not an artefact of the estimator |
| **3** | **`yuan2023` hexanal's ordinal row: `1, 1, 1, 1` across 22 → 130 °C** | flat | ordinal, so weak — but the same table registers 0 → 1 and 1 → `*` crossings for compounds that DO accelerate (cyclotene ≥ 81 kJ/mol), so the flatness is informative, not merely uninformative |
| **4** | **`hidalgo2010` thiol-Michael, 8 temperatures, 80–180 °C** | **28–30 kJ/mol** | a second lab confirming that *direct addition to an electrophile* lives in band A while *the same lab's* Strecker and browning cascades live at 50–67 |

**The most generous defensible upper bound is 45 kJ/mol** (`shepelev2024` §6e: adding back a
hypothetical hydrophobic pre-equilibrium ΔH_bind of −25 kJ/mol, which the paper does not
measure and which the repo would not want anyway, since its sink term acts on total rather
than bound aldehyde). **Even there, 45 < 70.**

## 2c. ⚠️ THE ONE COMPETING NUMBER, AND WHY IT IS NOT A COUNTER-EXAMPLE

`yuan2023`'s abstract states that **IC pasteurization (63 °C, 30 min) produced a similar extent
of reaction to UHT (130 °C, 30 s)**. Converting equal extent to equal time-integrated rate
gives `k(130)/k(63) = 1800/30 = 60`, hence **Eₐ ≈ 68.9 kJ/mol [Z]** — right at the threshold.
**It must not be opposed to the 15–20, for four reasons developed in
`yuan2023_extraction.md` §6b:**

1. **It describes a 46-compound panel whose response is dominated by the compounds that
   actually move** — isothiocyanates, trisulfides, quinone methides, activated enones. The flat
   compounds (hexanal, furfural, diacetyl, both acetals, three thiols, DMDS) contribute nothing
   to a "similar extent" impression because they look identical everywhere.
   **Consistently, `shepelev2024` measures 36–43 kJ/mol for an isothiocyanate on the same
   protein — within 1.6–1.9× of the panel inference.**
2. **"Similar" is an eyeball comparison of two mass spectra, asserted for EUGENOL** — the one
   compound the paper itself proves is reacting via a **trace autoxidation contaminant in
   limited supply**, i.e. whose extent is exhaustion-limited, not rate-limited.
3. **The 130 °C extent is biased downward by coagulation** (the `*` cells), so 68.9 is a lower
   bound on whatever the panel's effective Eₐ is — pushing away from reconciliation and
   confirming it measures a different thing.
4. **For hexanal specifically, `yuan2023` carries ZERO temperature information** (`1,1,1,1`).

**⇒ Register `yuan2023_panel_Ea: 68.9 [Z], class: mixed_panel_ordinal_eyeball, DO NOT USE AS
AN ALDEHYDE Eₐ`.**

## 2d. ⭐⭐⭐ THE ARITHMETIC — what Eₐ = 15–20 does to the process-temperature sink

**Inputs.** The corpus's hexanal covalent-sink bracket: **t½ = 37–74 d at ambient in a 3 wt %
protein matrix** (`meynier2004` §1(c); the overlap of two independent methods). This wave uses
**t½(25 °C) = 37 d**, the FAST end, throughout — every number below is therefore an
**over**-estimate of the sink.

**Acceleration factor** `A(T) = exp[(Eₐ/R)(1/298.15 − 1/T)]`:

| Eₐ (kJ/mol) | **→ 100 °C** | **→ 130 °C** | **→ 140 °C** | **→ 180 °C** |
|---:|---:|---:|---:|---:|
| **15** | 3.4× | 4.5× | 5.4× | 7.9× |
| **20 ⭐ (measured)** | **5.1×** | **8.2×** | **9.5×** | **15.8×** |
| 45 *(generous bound)* | 38× | 100× | 156× | 498× |
| **70** *(the assumption)* | **291×** | **1,560×** | **2,596×** | **15,700×** |

**⇒ FRACTION OF THE HEXANAL DOSE CONSUMED BY THE COVALENT CHANNEL DURING A REAL PROCESS:**

| process | Eₐ = 15 | **Eₐ = 20 (measured)** | Eₐ = 45 (bound) | **Eₐ = 70 (the assumption)** |
|---|---:|---:|---:|---:|
| **UHT, 130 °C, 30 s** | 0.003 % | **0.005 %** | 0.07 % | **1.0 %** |
| **Extrusion, 140 °C, 60 s** | 0.007 % | **0.012 %** | 0.20 % | **3.3 %** |
| **Retort hold, 100 °C, 30 min** | 0.13 % | **0.20 %** | 1.5 % | **4.6 %** |
| **Bake / fry, 180 °C, 10 min** | 0.10 % | **0.21 %** | 6.3 % | **90.7 %** |

**⇒ ⭐⭐⭐ THE VERDICT.
At the measured activation energy, the covalent aldehyde–protein channel removes
0.005 % – 0.21 % of the hexanal during any realistic thermal process. Under the 70 kJ/mol
assumption it would have removed 1 % – 91 %. THE ASSUMPTION WAS WRONG BY TWO TO THREE ORDERS
OF MAGNITUDE IN THE QUANTITY THAT MATTERS.**

**Against the gap it was invoked to explain.** `k2_matrix_and_thresholds.md` §A.1 records the
hexanal cooked-beef/water threshold shift as **1,304× = 3.115 decades**. Removing a fraction
*f* of the hexanal contributes `−log₁₀(1 − f)` decades. At **f = 0.21 %** that is
**0.0009 decades — 0.03 % of the observed shift** [Z]. `meynier2004` §9b computed **0.06 % over
a 3 h ambient panel and 14 % over a full month**; **this wave's contribution is to show that
going to process temperature does NOT rescue the channel, because the acceleration is 9.5×
rather than 2,596×.**

## 2e. ⭐ TWO INDEPENDENT ARGUMENTS THAT DO NOT DEPEND ON THE ARITHMETIC AT ALL

**The rate argument above could in principle be attacked on the Eₐ value. These two cannot.**

**(i) THE CAPACITY OF THE SINK FALLS WITH TEMPERATURE — measured directly.**
`shepelev2024` Fig. 3 at day 28: **25 °C = 25.6 mg/g > 45 °C ≈ 20.1 > 65 °C ≈ 17.4.**
The 65 °C series **peaks at day 7 and then declines 26 %**. The authors state the mechanism
themselves: *"storage of a protein at elevated temperatures may alter protein structure such
that it becomes **less available/reactive**. **We favor this explanation.**"* The same inversion
recurs in the PEITC extrapolation. And `shepelev2024` §7e gives the governing principle:
*"the **number of reactive sites on the surface** of a protein plays a greater role in
determining reaction rate/extent than the **total amount of reactive amino acids**."*
**⇒ Heating does not merely fail to speed the sink up; it shrinks the sink's ceiling.**

**(ii) THE EQUILIBRIUM ITSELF SHIFTS BACK ON HEATING — measured directly.**
`zamora2010` measures **Eₐ,forward = 44** and **Eₐ,reverse = 52 kJ/mol** for a conjugate
addition at an amine nitrogen — **the corpus's only matched forward/reverse pair.**
`Eₐ,rev > Eₐ,fwd` ⇒ `ΔH ≈ −8 kJ/mol` ⇒ **K_eq falls 3.0× from 25 °C to 180 °C [Z]**.
And the paper demonstrates it experimentally and dramatically: acrylamide that was **99.3 %
gone after 28 days at 60 °C was restored to 10.6 % by 20 minutes at 180 °C** — a **15×
recovery of an analyte that had, to all measurement, disappeared** — with the authors
quantifying that only **~10 % of the storage loss was irreversible**.

**⇒ THREE INDEPENDENT MECHANISMS, THREE LABORATORIES, ALL POINTING THE SAME WAY:
the rate barely rises (Eₐ 15–20), the capacity falls (Shepelev), and the equilibrium unwinds
(Zamora). The covalent aldehyde–protein channel is an AMBIENT-STORAGE channel. It is not a
process-temperature channel, and no amount of heat makes it one.**

## 2f. **What the channel IS good for, stated positively**

| regime | verdict |
|---|---|
| **Ambient storage, weeks to months, high-protein matrix** | ⭐ **REAL and SIZEABLE.** Decanal loses 14.5–20 % of dose over 28 d at 25 °C (`shepelev2024`); hexanal's t½ is 18–74 d. **This is where the sink belongs in the repo.** |
| **A 3 h sensory panel at ambient** | negligible — **0.06 %** of the hexanal log-shift (`meynier2004` §9b) |
| **Any thermal process, seconds to tens of minutes** | ⛔ **NEGLIGIBLE — 0.005 % to 0.21 % (§2d).** |
| **Explaining the 1,304× hexanal matrix threshold shift** | ⛔ **NO — 0.03 % of it.** `k2` §D.2's diagnosis stands: the shift is criterion mismatch + unverified pre-cook dosing, not chemistry. |

## 2g. ⚠️ WHAT THE VERDICT DOES **NOT** SETTLE — the honest residue

1. **⭐ Reversibility of an aldehyde–protein adduct is STILL never measured.**
   `meynier2004` §4a, `anantharamkrishnan2020b` §5d and `shepelev2024` §3d all assert
   irreversibility and none tests it. **`hidalgo2010` §4 shows exactly what the test looks
   like** (synthesize the adduct → heat it alone → quantify the released electrophile against a
   labelled IS), and `zamora2010` shows the amine analogue of the same electrophile behaves
   **oppositely** to the thiol one. **The null survives K6b. It is now a well-specified
   experiment rather than an open-ended gap.**
2. **Hexanal itself is still not measured with a temperature axis.** §5 argues this no longer
   blocks anything.
3. **The Eₐ transfer from decanal to hexanal is an argument, not a measurement.** It rests on
   Eₐ being a property of the reaction coordinate while chain length moves the pre-exponential
   (`shepelev2024` §9). The corroboration from `hamzalioglu2018`'s HMF + lysine at 10.0 kJ/mol
   — a completely different aldehyde — is what makes it more than an assertion.
4. **pH is uncontrolled in `shepelev2024` and unstated in `yuan2023`.** Both papers'
   aqueous systems are unbuffered. `zamora2010` §6a shows a **2.2× buffer effect at constant
   pH** for a reverse reaction. **The Eₐ is a ratio and is largely protected; the extents are not.**

---

# §3. ⭐⭐⭐ THE MILK-UNIT RESOLUTION — SETTLED

## 3a. The question
`k3_final_parameter_inventory.md` §C.18: *"Tian et al. 2019 (`10.3168/jds.2019-16796`) — the
only way to settle the `?/kg` unit; without it **seven milk thresholds carry a factor-of-1000
basis risk**. **OPEN**."*
`k2_matrix_and_thresholds.md` §A.4: *"The units cell prints a literal question mark, verified
at 900 dpi. Contextually µg/kg; **not asserted**. **Blocks any OAV**."*

## 3b. ⭐ THE ANCHOR — an exact numeric identity, not an inference
The `?/kg` table is **`tian2020.pdf` = Tian, Xu, Sun, Chen & Yu (2020),
*J. Dairy Sci.* 103:5863–5873, `10.3168/jds.2019-17880`, Table 1.** Its own footnote reads:
*"Key aroma compounds detected in the milk fan sample **Y6** in our previous study
(**Tian et al., 2019**)."*
**Its `Concentration (?/kg)` column is that Y6 column, and `tian2019` Table 1 is headed
`Concentration (μg/kg; RSD in % in parentheses)` in plain type.**

| compound | `tian2020` `(?/kg)` conc. | `tian2019` Y6, **μg/kg** | match |
|---|---:|---:|:---:|
| Propanoic acid | **347** | **347** | ✅ exact |
| Butanoic acid | **7,030** | **7,030** | ✅ exact |
| Octanoic acid | **1,719** | **1,720** | ✅ |
| Octanal | **29** | **29.3** | ✅ |
| Nonanal | **198** | **199** | ✅ |
| 2-Nonanone | **244** | **244** | ✅ exact |
| Ethyl hexanoate | **1,001** | **1,000** | ✅ |

**SEVEN OF SEVEN.** In `tian2020` Table 1 **the same `(?/kg)` notation heads both the
concentration and the threshold column.** The concentration column is now identified as a
verbatim reproduction of a column its source labels **μg/kg**.

> ## ⭐⭐⭐ **`?/kg` = `μg/kg`. The factor-of-1000 basis risk on the seven milk rows is CLOSED. `k3` §C.18's entry can be marked RESOLVED.**

## 3c. The seven rows, now unit-resolved and OAV-enabled

| compound | **matrix threshold (μg/kg)** | SD | `tian2019` reference threshold (μg/kg) | **[Z] matrix / reference** |
|---|---:|---:|---:|---:|
| Propanoic acid | **51,200** | 84.07 | 3.00 | **17,067×** |
| Butanoic acid | **7,500** | 24.64 | 13.0 | **577×** |
| Octanoic acid | **25,600** | 120.88 | 5.10 | **5,020×** |
| Octanal | **160** | 3.29 | 0.880 | **182×** |
| Nonanal | **1,600** | 29.41 | 1.10 | **1,455×** |
| 2-Nonanone | **52,000** | 116.04 | 32.0 | **1,625×** |
| Ethyl hexanoate | **1,024** | 34.28 | 3.00 | **341×** |

**Two independent confirmations of the resolution:**
- **`k2` §A.4 already recovered nonanal's aqueous threshold as ≈1.1 μg/kg from `Xin2026b` and
  computed ≈1,450×. This wave gets 1,455× from Tian's own cited 1.10. Agreement to 0.3 %.**
- **`tian2020b` (yogurt, same laboratory) measures dairy-matrix thresholds of
  5.43–29.0 mg/L = 5,430–29,000 μg/kg for small carbonyls. The seven milk rows read as μg/kg
  span 160–52,000 — the same range. Read as ng/kg they would be 30–300× LOWER than the yogurt
  set; read as mg/kg, 1,000× higher and physically absurd.**

## 3d. ⚠️ THREE CAVEATS THAT MUST TRAVEL WITH THE ROWS

1. **⭐ "Milk fan" is a CHEESE. The numerator and denominator matrices differ.** Concentrations
   come from milk-fan cheese; thresholds were measured by dosing into *"a fresh milk solution
   matrix"* whose fat, protein and pH are never stated. **`same_matrix: FALSE`.**
2. **The reference thresholds are UNSOURCED and mislabelled "in air".** `tian2019` states only
   *"The threshold values were taken from the literature"* and computes OAV against
   *"its detection threshold **in air**"*. **They are almost certainly aqueous values — which is
   why nonanal's 1.10 matches the corpus's independently recovered aqueous figure. Ingest as
   `unsourced_literature_threshold`; the ratio column is `cross_study_cross_method` per
   `k2` §D.2.**
3. **The SDs are NOT threshold uncertainty.** 0.16–2.1 % RSD on a 15-panellist sensory
   threshold is impossible, and **five of the seven values are exact power-of-two dilution-series
   steps** (1,024 = 2¹⁰; 25,600 = 25 × 1,024; 51,200 = 50 × 1,024). `k2` §A.4's flag stands and
   is now explained.

**Bonus: `tian2020b` adds three MEASURED dairy-matrix thresholds in an unambiguous `mg/L`**
(diacetyl 5.43, acetaldehyde 15.4, acetoin 29.0) **in the corpus's most fully specified matrix**
(0.170 % fat, 2.34 % protein, 6.75 % NFMS, pH 4.6, pasteurized **before** dosing ⇒
`thermal_step_after_dosing: FALSE`, the clean tier). **And a 500× matrix-penalty spread across
three small carbonyls in ONE matrix with ONE panel — a tighter refutation of a uniform matrix
factor than `k2` §D.1's cross-compound one.**

---

# §4. HOW THE PIECES FIT — three convergences worth naming

**C1 — THE FLAT-ROW / SMALL-Eₐ ASSOCIATION IS NOW CALIBRATED.**
`yuan2023`'s ordinal table has exactly nine compounds flat at `1` across 22 → 130 °C: hexanal,
furfural, diacetyl, both acetals, thiophenol, propanethiol, furfuryl mercaptan, DMDS.
**Four of those nine now have an independent measured Eₐ in band A** — the thiols at 28–30
(`hidalgo2010`) and hexanal's homologue decanal at 15–20 (`shepelev2024`). Meanwhile the
compounds that *do* cross bins (cyclotene 0→1) bound Eₐ **≥ 81**. **The ordinal scale is
therefore weakly but genuinely informative about Eₐ, and hexanal's flatness is evidence, not
absence of evidence.**

**C2 — THREE MEASURES OF THE SAME PROTEIN-ACCESSIBILITY EFFECT.**
`shepelev2024` §7e (qualitative: surface sites, not total residues) · `zamora2010` §7c
(quantitative: Lys > polyLys 4.2 kDa > polyLys 68.3 kDa) · `shepelev2024` §4d (thermal: the
plateau extent falls as the protein denatures). **⇒ Every rate constant measured on a FREE
amino acid overstates the protein case, and the overstatement grows with polymer size and with
thermal history. Register as a standing correction on every free-amino-acid rate row in the
repo.**

**C3 — THE OBSERVABLE-SUBSTITUTION HAZARD, NOW WITH THREE INDEPENDENT FALSIFICATIONS.**
`meynier2004` §6b (fluorophore intensity moves opposite to residue loss by ~7× each in one
experiment) · `hidalgo1993` §5b (colour and fluorescence are "always correlated" yet have
**different Eₐ**, 66.5 vs 50, a k_F/k_B ratio varying 7.4× across pH, and different 20-day pH
optima) · `zamora2013` §3c (Eₐ predicts yield at 37 °C, r = −0.83, p = 0.02, but **not** at
80 or 180 °C, p = 0.4 and 0.6). **⇒ Browning, fluorescence, Strecker-product formation and
adduct formation are four different observables with four different activation energies. The
repo must never substitute one lane's Eₐ into another's sink term. This is the mechanism by
which the 70 kJ/mol assumption arose in the first place.**

---

# §5. DO HEXANAL-SPECIFIC KINETICS STILL NEED THE TWO UMN THESES?

**⚠️ First, an honest note on the premise.** No file in
`data/lit/extraction_dossiers/` names "two UMN theses". The only thesis lead in the corpus is
`research_round4_nulls.md` §L4 — *"a University of Minnesota thesis by **F. Chan**"*, recorded
with *"⚠️ I did not find the record. Search the Conservancy before assuming it exists."*
**The natural referents for "the two UMN theses" are the dissertations underlying the two
papers this wave extracted — Igor Shepelev's and Jieyao Yuan's, both Reineccius-lab, both
UMN.** This section answers on that reading and flags the ambiguity.

**VERDICT: NO. Neither thesis is needed, and neither would close the remaining gap.**

| what a thesis would add | value | why it does not change anything |
|---|---|---|
| **Shepelev**: the raw DPM tables behind Figs. 2–4 | replaces this wave's ±0.15 mg/g digitization with exact values and per-point SDs | **The Eₐ is stable at 14–23 kJ/mol across a ±25 % perturbation of every single point (`shepelev2024` §6d). Exact values would tighten a number that is already 3.5× clear of the decision threshold.** |
| **Shepelev**: possibly a fourth temperature or a fourth compound | would improve the Arrhenius fit | R² is already 0.998–0.9997 on three points, and the PEITC control independently validates the estimator |
| **Yuan**: per-compound MS intensities behind the 0/1/2 scores | would convert an ordinal map into a semi-quantitative one | **would not produce a rate**, because each compound × condition is **one time point** with **no replicate count**. An extent at one time is not a kinetic. |
| **either**: hexanal | ⛔ | **Neither paper ran hexanal quantitatively. Shepelev's panel is decane / decanal / 2-undecanone / PEITC; Yuan's hexanal datum is the ordinal `1,1,1,1`. A thesis cannot contain data its paper's design never generated.** |

**What WOULD close the residual gap, in priority order:**

| # | experiment | what it settles | why it is now LOWER priority than before this wave |
|---:|---|---|---|
| **1** | **A reversibility test on an aldehyde–protein adduct** — `hidalgo2010`'s protocol (synthesize/isolate the adduct, heat it alone, quantify released aldehyde against a labelled IS) | **the corpus's only surviving hard null in this area** | Still the top item. But §2e(ii) already supplies the *sign* from a measured analogue. |
| **2** | Hexanal + a food protein at ≥ 3 temperatures, quantitative | would replace the decanal→hexanal Eₐ transfer with a direct measurement | **The verdict is 3.5–4.7× clear of its decision threshold and is corroborated by a second aldehyde (HMF) in a second lab. A hexanal measurement would refine a number, not reverse a conclusion.** |
| **3** | A same-method matrix-vs-water threshold pair | `k2` §D.2(i)'s standing null | untouched by this wave |

**⇒ RECOMMENDATION: do not order either thesis. The Eₐ question is answered; the reversibility
question needs a new experiment, not an old thesis.**

---

# §6. WHICH SUPPORTING-INFORMATION FILES ARE ACTUALLY NEEDED

| paper | SI exists? | contents | **needed?** |
|---|---|---|---|
| **`Shepelev2024`** | ✅ one PDF | *"Description of preliminary study"* — extraction-protocol development, decane recovery vs solvent and salt | ⛔ **NO.** Method development only; **no time-course, no temperature, no extent data.** Would speak only to the residual doubt in §2g(1), and weakly. |
| **`Yuan2023`** | ✅ one PDF | (i) deconvoluted spectra for **one representative compound per functional-group class** (13 spectra) · (ii) **Fig. S14**, GC-EI-MS of eugenol showing the trace enone · (iii) **Table S1**, "reason for study" per compound | ⛔ **NO.** No extents, no time courses, no replicates. Table S1 is already transcribed in `anantharamkrishnan2020b_extraction.md`. **Fig. S14 only if the eugenol false-positive is ever disputed.** |
| `zamora2013` | ❌ none | — | n/a |
| `hidalgo1993` | ❌ none | — | n/a |
| `hidalgo2010` | ❌ none | — | n/a |
| `zamora2010` | ❌ none | — | n/a |
| `tian2019` | ❌ none | — | n/a |
| `tian2020b` | ❌ none | — | n/a |

> ## ⇒ **NO SUPPORTING-INFORMATION FILE FROM ANY OF THE EIGHT PAPERS IS REQUIRED. Nothing in this wave's retrieval queue remains open.**

**The one adjacent retrieval this wave newly justifies is NOT an SI:**
**Zamora, Delgado & Hidalgo (2009), *Mol. Nutr. Food Res.* 53:1512–1520**, *"Conversion of
3-aminopropionamide and 3-alkylaminopropionamide into acrylamide in model systems"* — the
mechanism paper behind the elimination whose Eₐ = 52 kJ/mol is the corpus's only reverse-reaction
number. **Priority: LOW** — the Eₐ itself is already in hand from `zamora2010`.

---

# §7. DRAFT FIT vs HOLD-OUT ROLES

> ## ⚠️ **STATUS: DRAFT ONLY. NOT APPLIED. `docs/reference/FIT_HOLDOUT_DECLARATION.md` WAS NOT OPENED AND NOT EDITED. This section is a proposal for the orchestrator to accept, amend or reject, and it obeys the three rules of `k3` §D.0.**

## 7a. Rule compliance
1. **No dataset appears in both columns** — verified by construction below.
2. **The 21 frozen hold-outs of `k3` §D.1 are untouched** — none is moved, renamed, re-scoped
   or promoted. This wave proposes no change to any existing bundle.
3. **Every module touched holds out at least one dataset** — checked in §7d.

## 7b. Proposed roles

| dataset | module | **proposed role** | rationale |
|---|---|---|---|
| **`shepelev2024` Fig. 3, decanal, 25 / 45 / 65 °C** | covalent aldehyde–protein sink | ⭐ **FIT — the Eₐ only** | It is the only temperature-resolved quantitative aldehyde–protein measurement in existence. Holding it out would leave the sink term with **no** Eₐ at all, which is the state that produced the 70 kJ/mol error. **Fit `Ea = 20 kJ/mol`; do NOT fit the pre-exponential from it (chain-length mismatch, §1c).** |
| **`shepelev2024` Fig. 4, PEITC** | covalent sink | **HOLD-OUT** | An independent electrophile class measured in the same experiment. If the sink term's Eₐ machinery is right, it should predict 36–43 kJ/mol for PEITC without being shown it. **A free, clean, same-apparatus hold-out.** |
| **`shepelev2024` Fig. 2, 2-undecanone** | covalent sink | **HOLD-OUT** | A near-null (0.9 % of dose). Tests that the sink term does **not** fire on a mono-ketone. |
| **`hamzalioglu2018` HMF + Lys, Eₐ = 10.0** | HMF sink / covalent sink | ⚠️ **ALREADY ASSIGNED by K5a — do not re-scope.** Use here as a **cross-module consistency check** only, not as a new role. | Its dossier's USE-Q label and its `Ea(Coffee-Arg) = 3.275` refusal both stand. |
| **`yuan2023` Tables 2–4 (46 compounds × 4 conditions)** | covalent sink — **reactivity gate, not kinetics** | ⭐ **HOLD-OUT** | An **ordinal** classification with ~19 measured nulls (§9a of its dossier). Ideal for a **binary gate test**: does the repo predict a covalent sink for the compounds Yuan scores ≥ 1 and none for those it scores 0? **It cannot be a fit target** — no rate, no extent, no n. |
| **`zamora2010` 44 / 52 kJ/mol pair** | reversibility of adduct sinks | **FIT — the SIGN and ORDER only** | The corpus's only forward/reverse pair. Fit the qualitative structure (`Ea_rev > Ea_fwd` ⇒ K_eq falls with T). **The 8 kJ/mol magnitude is softer than the uncertainty (§9 of its dossier) — do not fit to it numerically.** |
| **`zamora2010` Table 2 storage/reheat** | reversibility | **HOLD-OUT** | The 15× recovery from 0.7 % → 10.6 % is a sharp, falsifiable prediction target for any reversible-sink implementation. |
| **`hidalgo2010` 28 / 30 kJ/mol** | thiol/sulfur sink | **FIT** | 8 temperatures, the densest series in the set; a clean band-A anchor for thiol chemistry. |
| **`hidalgo2010` §4 irreversibility of the thioether** | sulfur sink | **HOLD-OUT** | A measured negative: the repo must predict no thermal release from a thioether. |
| **`zamora2013` Table 2, eight Strecker Eₐ** | ⚠️ **lipid-oxidation / Strecker lane — NOT the covalent-sink lane** | **split: FIT the four tertiary-carbonyl rows (55.3 / 57.9 / 58.9 / 63.5); HOLD OUT 4-HNE (67.1) and both alkadienal rows (27.6 air, 78.0 N₂)** | The held-out rows are the ones that test structure rather than interpolation: 4-HNE is the class outlier and the alkadienal pair is the atmosphere test. |
| **`zamora2013` §3c Eₐ↔yield correlations** | any lane | **HOLD-OUT (a meta-test)** | The repo should reproduce that Eₐ ranks yields at 37 °C and fails to at 80 and 180 °C. |
| **`hidalgo1993` Eₐ = 66.5 (colour) / 50 (fluorescence)** | ⚠️ **browning lane ONLY** | **FIT the browning lane** | ⛔ **MUST NOT enter the covalent-sink lane** — the observables are optical and `meynier2004` §6b falsifies fluorescence as a covalent-extent proxy. |
| **`hidalgo1993` Table 1, ten-point pH series** | browning lane | **HOLD-OUT** | The corpus's densest carbonyl–amine pH curve; a strong pH-response test. ⚠️ carries a buffer confound at pH 7.5→8. |
| **`tian2020` seven milk-fan thresholds (unit now resolved)** | matrix threshold lookup | **HOLD-OUT** | They were previously unusable. **Newly usable ⇒ newly valuable as a hold-out, and they should NOT be spent as fit data.** Class: `cross_study_cross_method`. |
| **`tian2019` fifty literature thresholds** | matrix threshold lookup | ⛔ **NEITHER — reference data only** | Unsourced and mislabelled "in air". Usable as a **denominator of last resort** with an explicit provenance flag; never as a fit target or a hold-out. |
| **`tian2020b` three yogurt thresholds** | matrix threshold lookup | **HOLD-OUT** | Fully specified matrix, largest panel in the corpus (n = 45), clean thermal ordering. **The single best matrix-threshold hold-out available.** |
| **`tian2020b` §7a binary synergy 1.81–2.31×** | OAV composition | **HOLD-OUT** | Tests whether the repo's OAV summation is additive where it should not be. |

## 7c. What must NOT be done — three explicit refusals
1. ⛔ **Do not fit a hexanal covalent rate constant to `shepelev2024`.** Decanal is 12–15×
   faster (§1c). Fit the Eₐ; keep the rate at the `meynier2004` + `anantharamkrishnan2020b`
   bracket.
2. ⛔ **Do not use `yuan2023`'s 68.9 kJ/mol anywhere.** §2c.
3. ⛔ **Do not move `hidalgo1993`'s 66.5/50 into the covalent-sink lane.** §7b.

## 7d. Rule-3 check
| module | holds out at least one dataset? |
|---|---|
| covalent aldehyde–protein sink | ✅ PEITC · 2-undecanone · Yuan's 46-compound gate |
| reversibility | ✅ `zamora2010` Table 2 · `hidalgo2010` §4 |
| sulfur / thiol sink | ✅ `hidalgo2010` §4 |
| lipid-oxidation / Strecker | ✅ 4-HNE row · both alkadienal rows · the Eₐ↔yield meta-test |
| browning | ✅ `hidalgo1993` Table 1 pH series |
| matrix thresholds | ✅ `tian2020` seven rows · `tian2020b` three rows · the synergy factor |

---

# §8. CONSOLIDATED INTEGRITY FLAGS FROM THIS WAVE

| # | flag | paper | severity |
|---|---|---|---|
| **1** | `DEC` = decanal in Fig. 3 / Results / Summary but decane in the Abbreviations list | `shepelev2024` | ⭐⭐ CRITICAL |
| **2** | The text's "highest temperature/longest storage" maximum is contradicted by its own Fig. 3, whose maximum is the **25 °C** day-28 point | `shepelev2024` | ⭐⭐ HIGH — **and the contradiction is the result** |
| **3** | `*` = protein coagulated = `≥ 2`, not a missing value | `yuan2023` | ⭐⭐ HIGH |
| **4** | OCR text layer renders Table 1's pH-6.5 constant as `6.12` vs the printed `0.12` (51×) | `hidalgo1993` | ⭐⭐ CRITICAL |
| **5** | Columns headed "reaction rate" contain an ordinal **extent** score; no rate exists | `yuan2023` | ⭐ HIGH |
| **6** | "rates increase 2–4× per 10 K" is a **textbook generality**, not a measurement | `yuan2023` | ⭐ HIGH |
| **7** | Eugenol's positive is a **contaminant enone** (+178 ≠ +164 Da), shown by the authors | `yuan2023` | ⭐ HIGH |
| **8** | Buffer changes a reverse reaction **2.20×** at constant pH 6 (citrate vs phosphate) | `zamora2010` | ⭐⭐ HIGH — retro-flags `hidalgo1993`, `hidalgo2010` |
| **9** | Diacetyl's water threshold is `0.003` (Table 2) vs `0.0003` (Methods) — 10× internal | `tian2020b` | ⭐ HIGH |
| **10** | Part of the 1,810× matrix ratio is attributed by the authors to a **criterion change** (E1432 vs E679-04), not the matrix | `tian2020b` | ⭐ HIGH |
| **11** | The σ-τ "optimal settings" all fall below the Table 3 optima, and diacetyl's below its own threshold | `tian2020b` | ⭐ HIGH |
| **12** | 2,4-decadienal has **two** Eₐ (27.6 air / 78.0 N₂) — a single-valued row is wrong | `zamora2013` | ⭐ HIGH |
| **13** | Eₐ predicts yield only at 37 °C (r = −0.83, p = 0.02); not at 80 or 180 °C | `zamora2013` | ⭐ HIGH |
| **14** | Strecker pH rule is **acid-favoured**, the **opposite sign** to the adduct-formation rule elsewhere in the corpus — do not merge | `zamora2013` vs `hidalgo1993`/`yuan2023` | ⭐ HIGH |
| **15** | Eₐ 66.5 / 50 are **browning/fluorescence**, not adduct, activation energies | `hidalgo1993` | ⭐ HIGH |
| **16** | Thiol-adduct irreversibility must not transfer to a Schiff base or aza-Michael adduct | `hidalgo2010` | ⭐ HIGH |
| **17** | The 44/52 gap (8 kJ/mol) is smaller than the plausible uncertainty on either bare integer | `zamora2010` | ⭐ HIGH |
| **18** | 50 reference thresholds are **unsourced** and declared "in air" while divided into μg/kg solids | `tian2019` | ⭐ HIGH |
| **19** | "Milk fan" is a **cheese**; the threshold matrix is a separate unspecified "fresh milk solution" | `tian2019`/`tian2020` | ⭐ HIGH |
| **20** | `tian2020` threshold SDs (0.16–2.1 % RSD) are not threshold uncertainty; 5 of 7 values are 2ⁿ series steps | `tian2020` | ⭐ HIGH — carried from `k2` §A.4, now explained |
| 21 | "0.1 **mmol** of UDO per mol of BLG" is wrong by ~1000× | `shepelev2024` | MEDIUM |
| 22 | Dose recipe implies 1.53 : 1 flavour : lysine, not the stated 1 : 1 (1.38× molar-basis inconsistency) | `shepelev2024` | MEDIUM |
| 23 | pH never stated and unbuffered | `shepelev2024`, `yuan2023` | MEDIUM |
| 24 | 110 °C/30 min yielded no data for any compound (protein lost) | `yuan2023` | MEDIUM |
| 25 | Table 1 lists 44 compounds; abstract claims 46 | `yuan2023` | MEDIUM |
| 26 | Pentanoic acid Y6 prints `2` where the arithmetic requires ≈248 μg/kg | `tian2019` | MEDIUM |
| 27 | Several printed OAV cells disagree with concentration ÷ threshold by 1–2 digits | `tian2019` | MEDIUM |
| 28 | `research_round4_nulls.md` §A.1's abstract-derived Eₐ ranges **omit the 78.0 kJ/mol N₂ row** and mis-round two endpoints | `zamora2013` | MEDIUM — superseded by §1a here |
| 29 | k_FI at pH 7/25 °C is 3.86 (Table 1) vs 4.2 (text) — 9 % internal | `hidalgo1993` | MEDIUM |
| 30 | Matrix is an **acidified, unfermented skim-milk simulant** labelled "yogurt matrix" | `tian2020b` | MEDIUM |

---

# §9. THE THREE SENTENCES

1. **The activation energy of aldehyde–protein covalent adduct formation is 15–20 kJ/mol, not
   the ≥ 70 the corpus assumed — measured by ¹⁴C on β-lactoglobulin at three temperatures and
   corroborated at 10.0 kJ/mol by an unrelated laboratory on a different aldehyde — so the
   covalent channel removes 0.005 % to 0.21 % of the hexanal during any real thermal process
   instead of the 1 % to 91 % the assumption implied.**
2. **The channel is further disabled from two independent directions that need no Arrhenius
   arithmetic: its capacity FALLS with temperature (the coldest sample binds the most aldehyde,
   measured), and the equilibrium itself unwinds on heating (the reverse activation energy
   exceeds the forward one, measured, with a 15× recovery of "vanished" analyte to prove it).**
3. **`?/kg` = `μg/kg`: seven of seven concentration values in the ambiguous table are a
   digit-for-digit reproduction of a column its own source paper labels μg/kg, so the
   factor-of-1000 basis risk on the seven milk threshold rows is closed by arithmetic rather
   than by inference.**
