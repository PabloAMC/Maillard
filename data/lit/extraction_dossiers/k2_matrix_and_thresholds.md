# Wave K2 — matrix thresholds, protein binding, and precursor-dosed extrusion
### Ten-paper synthesis, 2026-08-28. Read-only wave: **no repo file was written, staged or modified.**

Ten dossiers accompany this file: `brewer1995_`, `vega1994_`, `tian2020_`,
`starkenmann2008_`, `andriot2000_`, `Barallat-Perez2024_`, `guo2020_`, `Conti2025_`,
`Conti2025b_`, `Xin2026_extraction.md`. Every paper read in full; Brewer Table 1, Vega
Table 1, Tian Table 1, Barallat-Pérez Tables 1–3 and Figure 4, Andriot Figures 1–7, Guo
Tables 1–2 and Xin Figure 7 were re-read from **200–900 dpi renders** because the text
layer was rotated, scrambled or (Tian) because the published PDF prints a literal `?`
where its unit should be.

Builds on `damodaran1981_extraction.md` and `z0_index.md` (the sibling Z0 wave). **The soy
binding constants and the 100 kDa basis are NOT re-derived here** — they are cited and
extended.

---

# THE FOUR-LINE ANSWER

1. **The matrix/water threshold ratios are NOT consistent enough to support a general
   correction factor. They span 3.4× to 6 714× — a factor of 2 000 — and the 1-σ band on
   any single factor is 27–41× wide.** Detail in §D.
2. **They are not even the same kind of quantity as the aqueous values they are compared
   to.** No same-method matrix-vs-water threshold pair exists anywhere in these ten papers,
   and **Brewer's beef "thresholds" are nominal doses added BEFORE a 70 °C cook** — a large
   part of the 100–6 700× is thermal loss, not perception. §D.2 is the load-bearing section
   of this file.
3. **The binding model cannot deliver matrix odour activity, and this is now quantified.**
   The measured ceiling on reversible hydrophobic binding is **1.3–3.7× at 4 % protein**
   (Andriot) and **≤7.6× extrapolated to 10 %**. Against threshold shifts of 100–6 700×,
   binding is a single-digit term. §B.4.
4. **What the batch does deliver is nine per-gram binding constants that need no molar
   mass, one recovered molar basis (β-lactoglobulin = the dimer), three fully specified
   precursor-dosed extrusion runs, and a sign-crossing result that inverts pentose ≫ hexose
   in a real extruder.** §B, §C.

---

# (a) THE THRESHOLD TABLE

**Reading conventions.** `M` = measured in the cited paper. `C` = cited by it from another
paper. Ratios in **bold** are computed by this wave and are **all `cross_study_cross_method`**
— see §D.2 before using any of them.

## A.1 Cooked lean ground beef (Brewer & Vega 1995) — the protein-matrix set

Matrix: lean beef ground and **mixed 1:1 with distilled water** (≈100 g protein/L), 15 g in
a sealed 25 mL vial, **cooked to 70 °C internal, held at 45 °C**, room 22 °C / 60 % RH.
Method: **triangle test, ascending series, 50 % correct above chance, ASTM E-1432-91**,
10 panellists, n/geometric mean = 180, group geometric mean.

| compound | matrix | threshold | aqueous threshold | **matrix/water** | method | source |
|---|---|---:|---:|---:|---|---|
| pentanal | cooked ground beef 1:1, 45 °C | **2.67 ppm** (M) | 12 ppb (C, Guadagni 1963) | **223×** | triangle / ASTM E-1432-91 | Brewer 1995 T1 |
| **hexanal** | " | **5.87 ppm** (M) | 4.5 ppb (C, Guadagni 1972) | **1 304×** | " | Brewer 1995 T1 |
| heptanal | " | **0.23 ppm** (M) | 3 ppb (C) | **77×** | " | Brewer 1995 T1 |
| t-2-hexenal | " | **7.87 ppm** (M) | 3 ppb (C) | **2 623×** | " | Brewer 1995 T1 |
| t-2-octenal | " | **4.20 ppm** (M) | 3 ppb (C) | **1 400×** | " | Brewer 1995 T1 |
| t,t-2,4-decadienal | " | **0.47 ppm** (M) | 0.07 ppb (C) | **6 714×** | " | Brewer 1995 T1 |

Individual-panellist ranges (printed): pentanal 1.41–5.55, hexanal **2.37–36.81**, heptanal
0.004–0.484, t-2-hexenal 2.04–10.62, t-2-octenal 1.12–26.30, decadienal 0.013–3.08 ppm.
**Hexanal's inter-panellist range alone is 15.5×.**
⚠️ Brewer's own "(M)" molar column is wrong by 10³ on four rows and 10⁴ on two — do not
ingest it (`brewer1995_extraction.md` §1).
⚠️ Brewer's prose claim of "three orders of magnitude" vs gelatin is refuted by its own
tables (true: 2.9–101×) — a named laundering hazard (`brewer1995_extraction.md` §3).

## A.2 3 % gelatin gel (Vega & Brewer 1994) — the clean protein-gel set

Matrix: **3 % w/v gelatin in distilled water (30 g protein/L), lipid-free**, 100 mL in a
250 mL sniffing flask, dosed at 22 °C and held 18 h at 4 °C. Method: **single-sample
sniffing vs a gelatin control, ascending series, LINEAR fit, 75 % detection (uncorrected)**,
16 panellists, 3 replicates.

| compound | 4 °C | **22 °C** | 37 °C | 60 °C | aqueous (C) | **gelatin/water @22 °C** |
|---|---:|---:|---:|---:|---:|---:|
| pentanal | 47 | **41** | 34 | 22 | 12 (Guadagni 1963) | **3.4×** |
| hexanal | 90 | **58** | 34 | 38 | 4.5 (Guadagni 1972) | **12.9×** |
| heptanal | 108 | **79** | 62 | 50 | 3 (Guadagni 1972) | **26.3×** |
| t-2-hexenal | 170 | **109** | 79 | 60 | 3 (Guadagni 1963) | **36.3×** |
| t-2-octenal | 140 | **109** | 105 | 81 | 3 (Guadagni 1972) | **36.3×** |
| t,t-2,4-decadienal | 112 | **64** | 89 | 64 | 0.07 (Guadagni 1972) | **914×** |
*(all values ppb; DOT fit `CD` ranges 0.49–0.99 — the hexanal 22 °C and 37 °C rows are
fitted at CD = 0.51 and 0.49.)*

**Temperature effect, within-study: 1.7–2.8× over 4 → 60 °C, non-monotone in 2 of 6.**
Far smaller than the compound-to-compound spread.

## A.3 Paraffin oil (Guadagni et al. 1972, **quoted second-hand** by Vega 1994 p. 234–235)

| compound | paraffin oil @22 °C (C) | aqueous (C) | **oil/water** |
|---|---:|---:|---:|
| hexanal | 120 ppb | 4.5 | **26.7×** |
| heptanal | 250 ppb | 3 | **83×** |
| t,t-2,4-decadienal | 135 ppb | 0.07 | **1 929×** |

A third medium for free, same qualitative shape as gelatin (polyunsaturated ≫ saturated),
larger in magnitude. `provenance: quoted_second_hand_by_vega1994`.

## A.4 Fresh milk (Tian et al. 2020) — ⚠️ UNITS UNREADABLE IN THE PUBLISHED PDF

Matrix: "**a fresh milk solution matrix**" — **fat, protein and pH not stated**. 10 mL in a
covered brown bottle at **20 °C**. Method: **3-AFC (ISO 13301), descending 2-fold series,
BET = geometric mean of last-missed and adjacent higher (ASTM E679 cited), 50 % detection**,
15 trained panellists.

| compound | milk threshold | printed unit | aqueous | **milk/water** |
|---|---:|---|---:|---:|
| propanoic acid | 51 200 | **`?/kg`** | — | — |
| butanoic acid | 7 500 | `?/kg` | — | — |
| octanoic acid | 25 600 | `?/kg` | — | — |
| **octanal** | **160** | `?/kg` | — | — |
| **nonanal** | **1 600** | `?/kg` | ≈1.1 µg/kg (recovered from Xin 2026, §A.6) | **≈1 450×** |
| **2-nonanone** | **52 000** | `?/kg` | — | — |
| ethyl hexanoate | 1 024 | `?/kg` | — | — |

**The units cell prints a literal question mark**, verified at 900 dpi. Contextually µg/kg;
**not asserted**. Blocks any OAV until Tian et al. 2019 (`10.3168/jds.2019-16796`) is
retrieved. The paper prints **no** aqueous comparison value for any compound.
⚠️ Do not ingest the SDs as threshold uncertainty (0.16–0.47 % RSD on a 15-panellist
sensory threshold; five of seven values are exact binary series steps).

## A.5 Meat slurry (Wick et al. 1967) — **second-hand in TWO of these papers, measured in neither**

| compound | matrix | threshold | source | how it reaches this batch |
|---|---|---:|---|---|
| methional | "beef" / "meat slurries" | 6.1 ppm | Wick 1967 | cited by Brewer 1995 p. 593 **and** Vega 1994 p. 234 |
| phenylacetaldehyde | " | 0.94 ppm | Wick 1967 | " |
| **1-nonanal** | " | **7.6 ppm** | Wick 1967 | " |

⚠️ **Described as "in beef" by Brewer and as "in meat slurries" by Vega — the same three
numbers, two different matrix labels.** Neither paper measured them. Attributing them to
Brewer or Vega would be laundering.
Cross-check worth having: **Wick's meat nonanal 7 600 ppb vs Tian's milk nonanal 1 600 —
a meat/milk ratio of 4.75×** for the same compound.

## A.6 Water — the measured aqueous anchors, and the ones recovered by arithmetic

| compound | matrix | threshold | method | source | class |
|---|---|---:|---|---|---|
| **(R/S)-3-sulfanylhexan-1-ol** | mineral water | **22 ng/L** | **3-AFC**, ascending 4–500 ng/L, 30 trained panellists, geometric mean of 2 sessions, **RETRONASAL (30 mL held in mouth 5 s)** | Starkenmann 2008 p. 9577 | **M** |
| same | water | 1–60 ng/L | not stated | Vermeulen et al. 2005, via Starkenmann | **C** |
| **hexanal** | water | **5.0 µg/kg** | not stated | **recovered by this wave from Xin 2026's own OAV arithmetic** | derived |
| **nonanal** | water | **≈1.1 µg/kg** | not stated | " | derived |
| **2-pentylfuran** | water | **5.800 µg/kg** | not stated | " (exact to 4 s.f. at both ends of the printed range) | derived |

**⚠️ The Starkenmann 22 ng/L is already an IN-MOUTH number** — the panellist holds the
sample for 5 s before judging. **A saliva correction applied on top of it would
double-count.**

## A.7 Saliva — the "quenching factor", corrected

| quantity | value | modality | note |
|---|---:|---|---|
| analytical free-thiol LOD, **water** | 0.001 mg/L | HPLC-fluorescence, Acrylodan | Starkenmann p. 9578 |
| analytical free-thiol LOD, **crude saliva** | **60 ± 10 mg/L** | same method, same instrument | " |
| **⇒ same-method saliva/water quench** | **6 × 10⁴** | **analytical** | **derived by this wave** |
| the paper's own printed figure | "3 × 10⁶" | **analytical LOD ÷ SENSORY threshold** | **cross-modality — do not use as a threshold ratio** |
| sensory saliva quench | **no number exists** | qualitative only | spat-out 0.1 mg/L solution "no longer detectable" |

Salivary protein measured: **0.59 ± 0.2 mg/mL, pH 7.6 ± 0.1**, pooled from 4 donors.
⚠️ The paper quotes a literature range "0.6 to 1.2 mg/mL **or from 230 to 250 mg/mL**" —
the second half is a transcription error in the source. Do not ingest.

## A.8 THE CROSS-MATRIX RATIO MATRIX — everything computable, in one place

| compound | **beef/water** | **gelatin/water** | **oil/water** | **beef/gelatin** |
|---|---:|---:|---:|---:|
| pentanal | 223× | 3.4× | — | 65× |
| hexanal | **1 304×** | **12.9×** | 26.7× | **101×** |
| heptanal | 77× | 26.3× | 83× | **2.9×** |
| t-2-hexenal | 2 623× | 36.3× | — | 72× |
| t-2-octenal | 1 400× | 36.3× | — | 38.5× |
| t,t-2,4-decadienal | 6 714× | 914× | 1 929× | 7.3× |
| **geometric mean** | **906×** | **33×** | 244× | **27×** |
| **max / min** | **87×** | **269×** | 72× | **35×** |
| **1-σ band (log10 SD)** | **0.71 dec ⇒ 27× wide** | **0.80 dec ⇒ 41× wide** | — | **0.63 dec ⇒ 18× wide** |

---

# (b) THE BINDING TABLE

`per-gram` in the basis column means **no protein molar mass is needed**; that is the
strongest provenance available.

| compound | protein | constant + units | **molar basis** | conditions | source |
|---|---|---:|---|---|---|
| **2-heptanone** | β-lactoglobulin | K_b = **330 M⁻¹** (global) ⇒ **8.97 × 10⁻³ L/g** | **36 800 g/mol — RECOVERED BY ARITHMETIC from Fig. 3 (band 32 500–42 800); dimer. NOT stated by source** | pH 3, 50 mM NaCl, 30 °C, 50 µL/L aroma, 0–4 % protein | Andriot 2000 |
| **2-octanone** | β-lactoglobulin | K_b = **950 M⁻¹** ⇒ **2.58 × 10⁻² L/g** | " | " | Andriot 2000 |
| **2-nonanone** | β-lactoglobulin | K_b = **2 440 M⁻¹** ⇒ **6.63 × 10⁻² L/g** | " | " | Andriot 2000 |
| **hexanal** | lupin protein isolate (ProLupin 10600, 91 % protein) | **5.15 × 10⁻² L/g** (from 34 % binding at 10 g/L) | **per-gram — none needed** | pH 7.0, 30 °C, 3 h, 125 rpm, 5 mg/L aroma, static HS-GC-MS | Barallat-Pérez 2024 Fig. 4 |
| **nonanal** | lupin protein isolate | **3.17 × 10⁻¹ L/g** (76 %) | **per-gram** | " | Barallat-Pérez 2024 Fig. 4 |
| **2-nonanone** | lupin protein isolate | **4.71 × 10⁻² L/g** (32 %) | **per-gram** | " | Barallat-Pérez 2024 Fig. 4 |
| **hexanal** | pig gastric mucin | **6.4 × 10⁻¹ L/g** (6 % at 0.1 g/L) | **per-gram** | " | Barallat-Pérez 2024 Fig. 4 |
| **nonanal** | pig gastric mucin | **2.35 L/g** (19 %) | **per-gram** | " | Barallat-Pérez 2024 Fig. 4 |
| **2-nonanone** | pig gastric mucin | **1.36 L/g** (12 %) | **per-gram** | " | Barallat-Pérez 2024 Fig. 4 |
| **3-sulfanylhexan-1-ol** | whole human saliva | **STRANDED — no basis, no site stoichiometry, mechanism unresolved** (see below) | — | pH 7.6, 0.59 g/L protein, 22 °C, 45 min | Starkenmann 2008 |
| *2-heptanone* | *soy (Damodaran, sibling wave)* | *n·K/MW ⇒ 4.40 × 10⁻³ L/g* | *100 000 g/mol, **stated by source*** | *pH 8.0, 30 mM Tris, 10 mM 2-ME, 25 °C* | *Damodaran 1981* |
| *2-octanone* | *soy* | *1.24 × 10⁻² L/g* | *stated* | *"* | *Damodaran 1981* |
| *2-nonanone* | *soy* | *3.72 × 10⁻² L/g* | *stated* | *"* | *Damodaran 1981* |
| *nonanal* | *soy* | *4.38 × 10⁻² L/g* | *stated* | *"* | *Damodaran 1981* |
| *hexanal* | *soy, part. denatured* | *1.47 × 10⁻³ L/g* | *per-gram (Arai capacity)* | *gel filtration* | *Arai 1970 via Damodaran* |
| *1-hexanol* | *soy, part. denatured* | *6.99 × 10⁻⁴ L/g* | *per-gram* | *"* | *Arai 1970 via Damodaran* |

**Non-constant supporting quantities from Andriot 2000** (all digitised from figures,
±5–8 %): gas–liquid partition coefficients `K_ga` = 7.2 × 10⁻³ / 1.03 × 10⁻² / 2.4 × 10⁻²
(C7/C8/C9, pH 3, 30 °C); mass-transfer coefficients `h_D` ≈ 2.0 / 2.5 / 5.0 × 10⁻⁵ m/s
(**apparatus-specific — the authors attribute the fit deviation to unstirred layers; do not
port**).

## B.1 The recovered molar basis for β-lactoglobulin — how, and how confident

Andriot never states the MW behind `c_b` (M). Inverting eq 3, `c_b = (K_ga/K_ga^eff − 1)/K_b`,
at 4 % w/w with the printed K_b and the digitised Figure 3:

| ligand | K_ga^eff/K_ga | c_b (M) | MW at 40 g/L | MW at 36 g/L (>90 % purity) |
|---|---:|---:|---:|---:|
| 2-heptanone | 0.75 | 1.010 × 10⁻³ | 39 600 | 35 600 |
| 2-octanone | 0.53 | 9.34 × 10⁻⁴ | 42 800 | 38 500 |
| 2-nonanone | 0.27 | 1.108 × 10⁻³ | 36 100 | 32 500 |

**Monomer 18 400; dimer 36 800. Three ligands agree at 32 500–42 800 — the dimer, and not
the monomer, by a factor of 2.** Corroborated by the paper's own cited stoichiometry:
"**one binding site of 2-nonanone per β-lactoglobulin dimer at pH 3**" (Charles et al. 1996,
quoted p. 4246). **Weaker provenance than Damodaran's printed 100 kDa — label it
`recovered_by_arithmetic`, not `stated_by_source`.**

## B.2 The three-protein, three-lab agreement on 2-nonanone

| protein | K_eff (L/g) | method | lab / decade |
|---|---:|---|---|
| soy (Damodaran) | 3.72 × 10⁻² | equilibrium dialysis + Klotz | Cornell, 1981 |
| lupin (Barallat-Pérez) | 4.71 × 10⁻² | static headspace depletion | Wageningen, 2024 |
| β-lactoglobulin (Andriot) | 6.63 × 10⁻² | headspace + model fit | INRA Dijon, 2000 |

**Three proteins, three methods, three labs, forty-three years — agreeing within 1.8×.**
That is the strongest cross-validation in the batch and it says the **ketone** per-gram
constants are transferable.

## B.3 The aldehyde constants are NOT transferable, and the mechanism is named

| compound | lupin (headspace) | soy (dialysis / gel filtration) | ratio |
|---|---:|---:|---:|
| hexanal | 5.15 × 10⁻² | 1.47 × 10⁻³ | **35×** |
| nonanal | 3.17 × 10⁻¹ | 4.38 × 10⁻² | **7.2×** |
| 2-nonanone | 4.71 × 10⁻² | 3.72 × 10⁻² | **1.3×** |

**A 35× gap on an aldehyde where the ketone gap is 1.3× is not a species difference.** The
mechanism is printed in Barallat-Pérez p. 8735: "**Aldehydes can bind to proteins through
reversible or irreversible mechanisms, such as cysteine-aldehyde condensation reactions and
Schiff base formation under certain conditions (e.g., pH 6−10), forming strong amide
linkages**", against "**ketones predominantly bind through weaker hydrophobic
interactions**". Barallat-Pérez is at **pH 7.0** (inside that window) and measures by
**headspace depletion**, which counts irreversible loss; Damodaran is at pH 8.0 but with
**10 mM 2-mercaptoethanol**, which suppresses cysteine–aldehyde chemistry, and measures by
**dialysis**, which does not count it.

**⇒ `method` must be a first-class field on every aldehyde binding record, and headspace-
derived and dialysis-derived aldehyde constants must never be pooled.** This is a new,
quantified constraint on `binding_constants.yml`.

## B.4 THE CEILING ON REVERSIBLE BINDING — the number that closes the mission's question

**Measured**, Andriot 2000 Figure 1, 0 → 4 % (40 g/L) β-lactoglobulin:
**headspace suppression 1.25× (heptanone), 2.4× (octanone), 3.7× (nonanone).**
**Extrapolated** by eq 3 to 100 g/L: **1.9× / 3.6× / 7.6×**.
**Barallat-Pérez's strongest constant** (nonanal, lupin, headspace) at 100 g/L: **32.7×**
— and that one already includes irreversible Schiff-base loss (§B.3).
**Damodaran's soy hexanal** at 100 g/L: **1.15×**.

| observed threshold shift | best-case binding explanation | fraction of the log-shift explained |
|---|---:|---:|
| hexanal, gelatin/water **12.9×** (30 g/L collagen) | ≤2.5× (using the *lupin* constant, an over-estimate for gelatin) | ≤36 % |
| hexanal, beef/water **1 304×** (100 g/L) | **6.2×** | **25 %** |
| t,t-2,4-decadienal, beef/water **6 714×** | no constant exists; even 33× would give | ≤40 % |

**⇒ Reversible hydrophobic binding is a single-digit-to-low-tens factor. Matrix threshold
shifts are two to four orders. The repo cannot obtain matrix-relevant odour activity from
`f_free`, and the gap is not closable by tuning `a_p`.**

## B.5 Two independent measurements say headspace suppression ≠ perceived suppression

| paper | protein | **headspace / nose-space suppression** | **perceived-intensity suppression** | ratio |
|---|---|---:|---:|---:|
| Andriot 2000 | 1 % β-lactoglobulin | 10 / 20 / 39 % (C7/C8/C9) | 15–24 / 11–23 / 4–27 % | rank order **does not** agree |
| Barallat-Pérez 2024 | 1 wv% lupin isolate | **31.4 / 72.4 / 36.8 %** (hexanal / nonanal / 2-nonanone) | **19.4 / 15.5 / 11.4 %** | **1.6× / 4.7× / 3.2×** |

Andriot states it outright (p. 4250): "**this effect is not well correlated with the
retention of the aromas for the protein**"; and his ANOVA finds **A × BLG not significant**,
i.e. **the protein effect on perceived intensity is an additive offset on the intensity
scale, not a multiplicative factor** — the opposite of what `f_free` produces.

**⇒ This pushes the same way as §B.4 but harder: perception is LESS sensitive to matrix
binding than headspace is. So a headspace-calibrated `f_free` OVER-predicts the threshold
shift — and the measured shifts are still 100–1 000× larger than any `f_free` can produce.
Whatever drives the matrix threshold effect, it is neither binding nor partition.** §D.2
identifies what it mostly is.

## B.6 The additional binding facts worth having

- **Mucin is 7–29× a stronger binder per gram than the food protein** (§b table). At the
  0.1 g/L used it removes 6–19 % of the free aroma; at the 2.16 g/L in the paper's own
  artificial-saliva recipe the same constants predict **58 / 84 / 75 %**.
- **Protein + mucin is strongly super-additive**: observed combined binding exceeds
  independent-Langmuir prediction by **+40 pp (hexanal)** and **+35 pp (2-nonanone)**;
  ≈0 for nonanal (saturated). **Never compose saliva and matrix binding as independent
  terms.**
- ⚠️ **Named laundering hazard**: Barallat-Pérez's abstract, "**In vitro mucin addition
  increased aroma binding four to 12-fold**", is the ratio to **mucin alone**. The effect of
  adding mucin **to the protein system** is **1.07–2.2×**. Ingesting the abstract's number
  would be wrong by 2–11×.
- **Chain-length slope, two independent proteins:** Andriot β-lg **2.72×/CH₂**, Damodaran
  soy **2.9×/CH₂** — agreeing to 6 %, against the shipped `a_p·Pow` implied **≈3.6×/CH₂**.
  Second lab, second protein, same conclusion as the S4 [P] item.
- **Saliva thiol binding is STRANDED and probably not binding at all.** At the 1 mg/L
  quench point the thiol is **7.45 µM** and the salivary protein **59 mg/L**; a 1:1
  stoichiometry would require an effective protein MW below **≈8 000 g/mol**, smaller than
  any major salivary protein. The apparent per-gram avidity is **10³–10⁴ L/g**, three to
  five orders above every hydrophobic constant in the table. The authors concede: "**Whether
  this absorption results from physicochemical interactions, chemical transformations, or
  covalent linkage to glycoproteins remains unclear.**"

---

# (c) THE EXTRUSION TABLE

| study | protein base | **added precursor + concentration** | process conditions | quantified products |
|---|---|---|---|---|
| **Guo et al. 2020**<br>*Food Hydrocolloids* 105, 105752 | **SPI 94.2 % protein** + wheat gluten at **0 / 10 / 20 / 30 / 40 %** of SPI (dry mass)<br>⚠️ abstract says "concentrate", Materials says isolate | **NOT a Maillard precursor.** "Natural flavor powder (spice flavoring extract, **β-cyclodextrin**, 11 % fat)" at **1 % dry mass of SPI**. ⚠️ **encapsulated delivery** | **twin-screw** FT-36, 8 zones **20/50/80/150/140/100/80/60 °C**, **L/D 26:1**, compression **4.6:1**, **die 10.0 mm**, material moisture **50/60/70/80 %**, feed **6 kg/h *or* 30 g/min (contradictory, 3.3×)**, screw speed not stated | **retention %** of 16 terpenes/phenylpropenes vs the raw blend, HS-GC-MS vs 2,4,6-collidine.<br>**Total: 52.46 → 44.07 → 32.24 → 23.04 %** over 50 → 80 % moisture (**2.28× loss, every step significant, 14/14 compounds agree**).<br>Wheat gluten: **35.46 → 42.04 → 44.07 → 43.78 → 33.39 %** (**optimum at 20–30 %**). |
| **Conti et al. 2025 I + II**<br>*Food Res. Int.* 208, 116169 and 218, 116938 | **Soy protein CONCENTRATE**, Arcon SM (ADM), **≥70 g/100 g protein, dry basis** | **THIAMINE hydrochloride, purity > 99 %, at 1.5 % w/w**, added **2 h before extrusion**. **One dose level only.** | **single-screw** RXPQ Labor 24, 5 zones, zones 1–3 **~40/60/80 °C**, zone 4 = zone 5 − 15 °C, **zone 5 = 180 / 160 / 140 °C**, SPC moisture **30 / 34 / 38 % DRY BASIS (= 23–28 % wet — LOW-moisture TVP, not HME)**, **L/D 15.5:1**, compression **3.3:1**, predie 5.8 mm, **die 3.6 mm**, feed **170 g/min**, **216 rpm**. ⚠️ **moisture confounded with temperature by design (diagonal, not factorial)**.<br>Post-extrusion: rehydrate **1:4 in ~100 °C water, 15 min** (mass gain **435 / 411 / 406 %**), + NaCl 1.0 g/100 g, + MSG 0.4 g/100 g, ± oil 7.0 g/100 g, freeze −18 °C, serve at 50 °C. | **98 and 71 volatiles** (30 % M/180 °C and 38 % M/140 °C, both no oil), **µg/g, SEMI-QUANT vs hexanal-D12**, + GC-O with detection frequencies.<br>**pyrazines 49.74 vs 3.38 µg/g (14.7×)** · **all sulfur 44.81 vs 19.78 (2.27×)** · **lipid-oxidation aldehydes 3.80 vs 7.65 (0.50× — MILD wins)** · **hexanal 1.95 vs 4.70 µg/g (the ONLY absolute values; the IS is hexanal-D12)** · largest single volatile **4-methyl-5-(2-hydroxyethyl)thiazole 20.25 µg/g** (the first-formed thiamine product) · **methional at "trace" but detected by 50–75 % of GC-O assessors in both.** |
| **Xin et al. 2026**<br>*Food Hydrocolloids* 182, 113124 | **NPI : SPI : WG = 20 : 70 : 10** (w/w, dry basis) | **EIGHT CARBOHYDRATES, each at 6 % w/w dry basis, ADDITIVE not substitutive, plus a true no-carbohydrate CONTROL**: glucose, **xylose**, **ribose**, fructose, β-glucan, wheat starch, maltodextrin, sucrose. **One dose level only.** | **twin-screw**, screw **704 mm × 22 mm**, **L/D 32:1**, **cooling die at 70 °C**, **240 rpm**, dry feed **2 kg/h**, water **3.7 kg/h**, **final moisture ≈65 % — TRUE HIGH-MOISTURE EXTRUSION**, zones **60/80/110/135/140/140/140 °C**, **3 independent runs per formulation** | **45 volatiles**, µg/kg **semi-quant vs cyclohexanone**, 16 with OAV > 1 (⚠️ OAVs computed against **aqueous** thresholds, mislabelled "in air").<br>**Total pyrazines: FR 6 621.64 > SU 2 200.79 > RI 1 389.14 > GL 1 149.61 > {BG, WS, MA, control} ≈700–780 > XY ≈600 µg/kg.**<br>2-ethyl-3,6-dimethylpyrazine in FR **5 638.91** · 2,3-dimethyl-5-propylpyrazine in FR **357.76**.<br>2-pentylfuran **BG 100 259.10 / WS 97 794.41 / control ≈97 000 / RI ≈52 000** µg/kg.<br>Aldehyde totals **GL 3 506.49 / FR 3 483.59 / MA 3 381.41 / control 2 433.19 / SU 2 346.88**; benzaldehyde **control 922.04 → GL 1 698.56**.<br>Ketone totals **BG 7 364.61 / WS 7 003.79 / control 6 778.15 / FR 6 536.71**. Alkane+alkene **FR 5 557.19 / RI 5 320.48 / BG 4 862.06**. |

## C.1 The three extrusion papers do not contradict each other — they measure different things

- **Guo measures RETENTION only** (survival of a pre-dosed, encapsulated flavour) and finds
  it falls with process severity.
- **Conti measures the NET** (formation minus retention, with a dosed precursor) and finds
  the severe condition wins **3.15× overall and 14.7× on pyrazines**.
- **Conti's own text reconciles them** (p. 8): shear-driven loss is real, "**Nevertheless,
  the greater number and quantity of volatile compounds in the FTSPs obtained at 30 % M/
  180 °C may be attributed to thiamine degradation, which results in the generation of more
  compounds.**"

**⇒ FORMATION BEATS RETENTION under precursor dosing. The repo's extrusion lane needs both
terms and must not fit one to data generated by the other.**

## C.2 ⚠️ THE PENTOSE ≫ HEXOSE ORDERING INVERTS IN A REAL EXTRUDER

| system | ordering | source |
|---|---|---|
| aqueous, 145 °C, 0.5 M phosphate, cysteine + sugar | **ribose ≫ glucose, MFT 10.4× / FFT 4.3×** | Hofmann & Schieberle 1998 (`hofmann1998_extraction.md`) |
| **HME, 140 °C, 65 % moisture, sheared, plant protein** | **fructose 4.77× ribose; xylose ≈0.86× the no-sugar control; sucrose (non-reducing) 1.58× ribose; ribose only 1.21× glucose** | **Xin 2026** |

**Corroborated independently by colour.** Xin's Table 1: **xylose L\* 58.49 and ribose 58.80
vs control 57.60; ΔE 42.19 and 41.96 vs control 43.41** — **both pentoses brown LESS than
adding no sugar at all.** A colorimeter and a GC-MS agreeing makes this very hard to write
off as an SPME artefact.

**This does not refute Hofmann** — those are sulfur volatiles from a cysteine/sugar system in
0.5 M buffer at 145 °C. **It bounds the transferability of the ordering to a sheared,
65 %-moisture, protein-rich, 140 °C matrix, and it is the only direct evidence on that
question the repo has.** Owner call: the directional panel's pentose > hexose row should
probably be **re-scoped by matrix**, exactly as Z0's finding 10 re-scoped PH-05/PH-07.

⚠️ Xin's own explanation (ribose diverted to furans) is **contradicted by its own data**:
**RI has the LOWEST 2-pentylfuran and 2-butylfuran of all nine samples**. The mechanism is
open.
⚠️ **Unremarked by the paper, and worth an owner's attention: sucrose — non-reducing — is
second on pyrazines AND is the only addition that LOWERS hexanal (−39 %).** The parsimonious
reading is that sucrose hydrolyses under HME shear and the active species is fructose,
which is also the top performer. That would make **fructose the single most effective
Maillard precursor in a high-moisture extruder**, in two independent arms.

## C.3 A clean composition-vs-process dissociation on the beany marker

| paper | axis varied | 2-pentylfuran response |
|---|---|---|
| Conti 2025b | **process severity** (30 %/180 °C vs 38 %/140 °C) | **0.26 vs 0.24 µg/g — NOT SIGNIFICANT** |
| Xin 2026 | **composition** (8 carbohydrates at 6 %) | **RI −46 %, GL −30 %, FR −21 % vs control; BG +3 %, WS +1 %** |

**Composition moves 2-pentylfuran; process severity does not.** Two papers, one dissociation,
directly relevant to the repo's beany off-note lane.

## C.4 What the extrusion batch does NOT contain

- **No dose–response for any precursor.** Conti uses one thiamine level (1.5 %); Xin uses one
  carbohydrate level (6 %); Guo uses one flavour level (1 %). **No sensitivity of yield to
  precursor loading is recoverable from any of the three.**
- **No zero-precursor control in Conti.** Every "thiamine-derived" assignment is
  literature-attributed (Dreher et al. 2003, a model *orange juice* system), not
  demonstrated. And the largest class in Conti's table, the pyrazines, is **not** a thiamine
  product at all — it is sugar/amino-acid Maillard chemistry from the SPC itself.
- **No absolute quantification, except Conti's hexanal.** Both Conti and Xin are
  single-internal-standard semi-quantifications with no response factors, and Xin says
  "relative concentrations" itself. Xin's 2-pentylfuran at **50–100 ppm** is a DVB/CAR/PDMS
  partition artefact, not a concentration.
- **Xin's Tables S2 and S3 (all 45 volatiles × 9 formulations, and the OAVs and thresholds)
  are online supplementary material and are NOT in the PDF.** Retrieving them would turn
  ~12 quoted compounds into **405 quantified cells**.

---

# (d) THE VERDICT ON A GENERAL MATRIX CORRECTION FACTOR

## D.1 The direct answer: NO. The ratios are compound- and matrix-specific.

| statistic | beef/water | gelatin/water | beef/gelatin |
|---|---:|---:|---:|
| geometric mean | **906×** | **33×** | **27×** |
| minimum | 77× (heptanal) | 3.4× (pentanal) | 2.9× (heptanal) |
| maximum | 6 714× (decadienal) | 914× (decadienal) | 101× (hexanal) |
| **max / min** | **87×** | **269×** | **35×** |
| **log10 SD** | **0.71 decades** | **0.80 decades** | **0.63 decades** |
| **1-σ band width** | **27×** | **41×** | **18×** |

**What a single factor would cost.** Applying the gelatin/water geometric mean (33×)
uniformly to the six compounds Vega measured:

| compound | true ratio | error of a 33× uniform factor |
|---|---:|---|
| pentanal | 3.4× | **9.8× too high** |
| hexanal | 12.9× | 2.6× too high |
| heptanal | 26.3× | 1.27× too high |
| t-2-hexenal | 36.3× | 0.92× (about right) |
| t-2-octenal | 36.3× | 0.92× |
| t,t-2,4-decadienal | 914× | **27.7× too LOW** |

**A single factor would misplace the OAVs of the two extreme compounds by 10× and 28× in
opposite directions.** That is not an improvement over doing nothing; it would create
confident wrong answers where the repo currently has an acknowledged gap.

## D.2 ⚠️ AND THE RATIOS ARE NOT MEASURING WHAT THEY APPEAR TO MEASURE

**This is the most important finding of the wave and it changes how the whole (a) table
should be read.**

**(i) NO SAME-METHOD MATRIX-VS-WATER PAIR EXISTS ANYWHERE IN THESE TEN PAPERS.** Every
ratio in §A.8 crosses labs, decades and criteria:

| dataset | threshold criterion | task |
|---|---|---|
| Brewer 1995, beef | **50 % correct above chance** | triangle (forced-choice), ASTM E-1432-91 |
| Vega 1994, gelatin | **75 % detection, UNCORRECTED for chance** | single sample vs a named control |
| Tian 2020, milk | **50 % detection probability** | 3-AFC, ISO 13301, descending series |
| Guadagni 1963/72, water | not stated in either quoting paper | classical dilution |
| Starkenmann 2008, water | 50 % | 3-AFC, ascending, **retronasal** |

Vega asserts its water references were obtained "**using similar sensory techniques**"; that
is **not verifiable from the paper** and is contradicted by its own methods (a 75 %
uncorrected difference criterion is not a forced-choice criterion). **A 75 % uncorrected
criterion yields a systematically HIGHER threshold than a 50 % forced-choice one.** A
plausible size for that bias is ~2×, sign known, magnitude not.

**(ii) BREWER'S BEEF NUMBERS ARE NOMINAL DOSES ADDED *BEFORE* A 70 °C COOK.**

Brewer, p. 593: the aldehyde is mixed into **raw** ground beef, sealed, then "**cooked in a
circulating water bath to 70 °C (internal temperature) and maintained at 45 °C**", 15–20 min
before evaluation. Vega, p. 231: the gelatin is microwaved **first**, cooled to 22 °C, and
**then** dosed, then held 18 h at 4 °C. **Neither paper analytically verifies the
concentration at the moment of sniffing.**

**⇒ Brewer's "5.87 ppm hexanal in beef" is the dose you must ADD TO RAW BEEF BEFORE COOKING
to make hexanal detectable AFTER cooking. It is not the concentration present in the matrix
at the moment of perception.** Every thermally driven loss — evaporation into the vial
headspace, Schiff-base condensation with lysine on a denaturing protein, further oxidation —
comes straight off the numerator. Brewer's own discussion says the mechanism is active:
"**During heating, protein unfolding increased, enhancing the binding of aldehydes with
protein.**"

**This single protocol difference plausibly accounts for most of the beef/gelatin ratio for
the saturated and mono-unsaturated aldehydes (65×, 101×, 72×, 38.5×)** — and it explains the
otherwise inexplicable exceptions: **heptanal (2.9×) and decadienal (7.3×)** are the two
compounds Brewer reports at *sub-ppm* thresholds, i.e. the two dosed at concentrations low
enough that a fixed fractional loss matters least in the ratio.

**⇒ RECLASSIFY. Brewer 1995 is a `dose_added_pre_cook` dataset, not a
`threshold_in_matrix` dataset. Vega 1994 is much closer to a true matrix threshold (no
thermal step after dosing), and Tian 2020 and Starkenmann 2008 are cleanest of all (dosed
into a liquid at the serving temperature).** The (a) table should carry a
`concentration_verified: false` and a `thermal_step_after_dosing: true/false` field on every
row, and the beef/water column should be labelled **`added_dose_ratio`, not
`threshold_ratio`.**

**(iii) THE PHYSICS CANNOT PRODUCE THE OBSERVED NUMBERS.** §B.4 caps reversible binding at
6.2× for hexanal in a 100 g/L matrix; §B.5 shows perception is *less* sensitive to matrix
binding than headspace is, so the true perceptual factor is smaller still. Gelatin is a
collagen hydrolysate with essentially no hydrophobic binding pockets, and it still shows
3.4–914×. **No partition or binding mechanism in this batch can generate two to four orders
of magnitude.** What remains: (1) analytical non-verification of the numerator (ii);
(2) criterion mismatch (i); (3) **irreversible chemistry**, which is the only mechanism that
scales with the right variable — see D.3.

## D.3 THE ONE STRUCTURE THAT DOES TRANSFER: α,β-UNSATURATION, ≈2–3×

Two independent matrices, two independent panels, same compound pair:

| matrix | hexanal | t-2-hexenal | **alkenal penalty** |
|---|---:|---:|---:|
| 3 % gelatin, 22 °C | 12.9× (vs water) | 36.3× | **2.8×** |
| cooked beef, 45 °C | 1 304× | 2 623× | **2.0×** |

**And it is chemically motivated**: 2-alkenals are Michael acceptors for protein
lysine/cysteine, so a matrix penalty that *grows* with unsaturation while *ignoring*
hydrophobicity is exactly the signature of irreversible adduction rather than partition.
t-2-Hexenal has a **lower** logP than hexanal and a **2.8× larger** matrix penalty.
And the polyunsaturated **t,t-2,4-decadienal is the extreme outlier in all three media** —
914× (gelatin), 1 929× (paraffin oil), 6 714× (beef).

**The chain-length structure, by contrast, does NOT transfer.** In gelatin the saturated
n-alkanals are monotone — **3.4 → 12.9 → 26.3 for C5 → C6 → C7, i.e. 2.78× per CH₂**, and
strikingly close to the two independently measured binding slopes (**Andriot 2.72×/CH₂,
Damodaran 2.9×/CH₂**). **In beef it collapses: 223 → 1 304 → 77.** Heptanal breaks the
monotone by 17×.

## D.4 THE RECOMMENDATION

**Do NOT implement a general matrix correction factor.** Instead:

1. **Ship a LOOKUP TABLE, not a formula.** Per compound × per matrix × per temperature,
   populated only where a measurement exists (§A), with an explicit
   **`no_factor_available`** state everywhere else. Twenty-four gelatin values, six beef
   values, seven milk values and three paraffin-oil values is the whole world; anything
   beyond it is extrapolation.
2. **Carry the provenance fields the ratios actually need**: `criterion`
   (50 %-forced-choice / 75 %-uncorrected / BET), `thermal_step_after_dosing`,
   `concentration_verified`, `aqueous_reference_source`, and
   `cross_study_cross_method: true` (which is **every** ratio in this batch).
3. **Relabel Brewer 1995 as an added-dose dataset** (§D.2 ii). It remains valuable — it is
   an engineering answer to "how much can I let form before someone smells it" — but it is
   not a threshold in the repo's sense.
4. **If any parametric term is wanted, the only defensible one is the α,β-unsaturation
   penalty of ≈2–3×** (§D.3), which reproduces in two matrices and has a mechanism.
   The chain-length term reproduces in gelatin and fails in beef; do not ship it.
5. **Stop deriving matrix odour activity from `f_free`.** §B.4 and §B.5 close that route
   quantitatively: binding gives single digits where the data show two to four orders, and
   perception responds *less* than headspace does. The binding model's proper job is
   predicting **headspace concentration**, which is what it was calibrated on.
6. **Fix the OAV denominator before anything else.** Xin 2026 (§A.6, §5 of that dossier) is
   a live 2026 example of the repo's exact defect: `matrix concentration ÷ aqueous
   threshold`, with the denominator mislabelled "in air". The recovered thresholds it
   actually used — **hexanal 5.0, nonanal ≈1.1, 2-pentylfuran 5.800 µg/kg** — are the
   classical water values. Citing this when the correction is argued for shows the practice
   is standard and standardly wrong.

---

# NAMED LAUNDERING HAZARDS — the "342/200 list" for this batch

Ranked by how likely a future wave is to import them.

| # | claim, as printed | reality | source |
|---:|---|---|---|
| **1** | **"OAV ... calculated by dividing its concentration in the extrudate (μg/kg) by its odor threshold in air (μg/kg)"** | The thresholds used are the classical **AQUEOUS** ones (recovered exactly: 2-pentylfuran **5.800**, hexanal **5.0**, nonanal **≈1.1** µg/kg). Both the operation and the label are wrong. | Xin 2026 p. 3 |
| **2** | **"esters, alkanes, alkenes, phenols, aldehydes, and alcohols"** retention ranking | **Zero alcohols are measured**; alkanes have no class total; esters are one compound; "maltol" is a pyranone, not a phenol. The tables cannot produce the ranking. | Guo 2020 abstract + conclusions |
| **3** | **"three orders of magnitude for pentanal, hexanal, t-2-hexenal, t-2-octenal, and two orders for heptanal and t,t-2,4-decadienal"** vs gelatin | True ratios from the same authors' own two tables: **65×, 101×, 72×, 38.5×, 2.9×, 7.3×.** Heptanal is wrong by 35×. | Brewer 1995 p. 593 |
| **4** | **"In vitro mucin addition increased aroma binding four to 12-fold"** | That is the ratio to **mucin alone**. Adding mucin **to the protein system** gives **1.07–2.2×**. | Barallat-Pérez 2024 abstract |
| **5** | **"the detection limit ... in crude saliva was 60 mg/L, which is 3 × 10⁶ above its odor threshold"** | An **analytical LOD** divided by a **sensory threshold**. The same-method figure is **6 × 10⁴**. | Starkenmann 2008 p. 9578 |
| **6** | **Table 1's molar column** (3.0 × 10⁻⁸ M for a 2.67 ppm pentanal threshold) | Off by 10³ on four rows and 10⁴ on two. Only the ppm column is self-consistent. | Brewer 1995 T1 |
| **7** | **Table 1's three class-total rows** (alkenes, phenols, aldehydes) | **Shifted one column right**, proved against the paper's own shared replicate (sample 3 ≡ sample 7). Fourteen individual rows and the grand total are correct. | Guo 2020 T1 |
| **8** | **"Imax_S decreased by 25.92 %"** (hexanal, protein addition) | Table 2 gives 67 → 54 = **−19.4 %**. Not rounding. | Barallat-Pérez 2024 p. 8736 |
| **9** | **"hexanal was 46.94 % less persistent than nonanal"** and **"surpassing 2-nonanone by 15.93 %"** | Wrong-base percentages. Raw seconds are 60 / 74 / 88; the true figures are 31.8 % and 18.9 %. | Barallat-Pérez 2024 p. 8736 |
| **10** | **methional 6.1 ppm / phenylacetaldehyde 0.94 ppm / nonanal 7.6 ppm "in beef"** | **Wick et al. 1967**, cited by **both** Brewer and Vega, measured by neither, and described as "in beef" by one and "in meat slurries" by the other. | Brewer 1995 & Vega 1994 |
| **11** | **"the retention rates ... is not affected by the wheat gluten and moisture content"** | Table 2's aldehyde total falls **2.7×** and its alkene total **2.5×** across the moisture series, at different rates. | Guo 2020 conclusions |
| **12** | **"2,5-dimethyl-3-(3-methylbutyl) pyrazine"** listed among the compounds higher at 38 % M/140 °C | **That compound appears nowhere in the paper's Table 2.** Nine of the ten named exceptions do; that one does not. | Conti 2025b p. 8 |
| **13** | Milk-fan threshold **SDs** of 84–121 on values of 25 600–52 000 | **0.16–0.47 % RSD** on a 15-panellist sensory threshold, with five of seven values being exact binary series steps. Not threshold uncertainty. | Tian 2020 T1 |
| **14** | Barallat-Pérez **Table 1** physicochemical values | All **cited from ref 31**, none measured. (And the "nonanal" structure drawn is a methyl ketone.) | Barallat-Pérez 2024 T1 |

**Plus four unresolved internal contradictions**, recorded but not resolvable from the PDFs:
Guo's feed rate (**6 kg/h vs 30 g/min, 3.3×**); Guo's LF-NMR shared replicate (3 of 7 rows
disagree where all 15 flavour rows agree); Conti 2025b's meat-aroma direction (**stated
three ways in one paper**: 38 % is meatier p. 3, no effect p. 5, 30 % is meatier p. 11);
and Tian's σ-τ ladders (**pairs (a) and (b) seeded with the wrong compound's threshold**).

---

# WHAT THE OWNER SHOULD DECIDE — five items, in order

1. **Adopt the lookup-table decision (§D.4.1) rather than a correction factor**, and add the
   provenance fields in §D.4.2 to `odour_threshold` records. This is the mission's headline
   question and the answer is unambiguous.
2. **Reclassify Brewer 1995 as `dose_added_pre_cook`** (§D.2 ii). It is the single change
   that most improves the honesty of the threshold lane, and it downgrades the most
   dramatic-looking numbers in this batch.
3. **Decide whether the pentose ≫ hexose directional row should be re-scoped by matrix**
   (§C.2). Xin 2026 inverts it in a real 65 %-moisture extruder, with independent
   corroboration from colour. This is a panel-membership change and therefore an owner call,
   directly analogous to Z0's PH-05/PH-07 item.
4. **Split aldehyde binding constants by method** (§B.3). Headspace-derived and
   dialysis-derived aldehyde constants differ by 7–35× for mechanistic reasons the source
   states; pooling them is currently possible in `binding_constants.yml` and should not be.
5. **Promote the β-lactoglobulin constants with `recovered_by_arithmetic` provenance, or
   don't promote them** (§B.1). The dimer basis is well supported (three ligands, plus the
   paper's own cited stoichiometry) but it is an inference from a digitised figure, not a
   printed number, and it must not be recorded at the same confidence as Damodaran's stated
   100 kDa.

# THE THREE RETRIEVALS WORTH REQUESTING NEXT

1. **Xin et al. 2026 Supporting Information** (`10.1016/j.foodhyd.2026.113124`): **Tables S2
   and S3** — all 45 volatiles × 9 formulations in µg/kg, and the OAVs **with the threshold
   list**. Turns ~12 quoted compounds into **405 quantified cells** and lets §A.6's
   threshold recovery be audited directly. **Highest yield of any single retrieval across
   all ten papers.**
2. **Tian et al. 2019**, *J. Dairy Sci.* 102:9639–9650 (`10.3168/jds.2019-16796`). The only
   way to settle the `?/kg` unit and the milk matrix composition. Without it, seven milk
   thresholds carry a **factor-of-1000 basis risk** and are unusable.
3. **Harrison & Hills 1997**, *J. Agric. Food Chem.* 45, 1883–1890. Primary source of
   Andriot's eqs 1–6; settles the eq-4 bracket ambiguity and states the `c_b` convention
   Andriot leaves implicit — which would upgrade §B.1 from `recovered_by_arithmetic` to
   `stated_by_source`.

*Runners-up:* **Milani, Menis-Henrique & Conti 2022** (the only plausible source of a
thiamine **dose–response**, which neither Conti paper provides) and **Guadagni, Buttery &
Turnbaugh 1972**, *J. Sci. Food Agric.* 23, 1435–1444 — the primary source of **four of the
six aqueous thresholds** that every ratio in §A.8 divides by, currently held only
second-hand through Vega 1994.
