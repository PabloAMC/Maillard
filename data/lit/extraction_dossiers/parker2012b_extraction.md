# Parker, Balagiannis, Higley, Smith, Wedzicha & Mottram 2012 — EXTRACTION (file stem `parker2012b`)
DOI: 10.1021/jf302415n (printed in every page footer). J. Agric. Food Chem. **60**(36):9321-9331; received 4 Jun 2012,
published 28 Aug 2012. Crossref 2026-11-03: "Kinetic Model for the Formation of Acrylamide during the Finish-Frying of
Commercial French Fries", JAFC 60(36) 9321-9331, issued 2012-09-04, authors Parker, Balagiannis, Higley, Smith,
Wedzicha, Mottram. **Matches the PDF and the brief.**
**Naming caveat.** The stem is `parker2012b` because `data/articles/parker2012.pdf` is a different paper: Parker JK,
"The kinetics of thermal generation of flavour", J. Sci. Food Agric. 93:197-208 (2013, online Nov 2012,
10.1002/jsfa.5943), aliased to `parker2013_extraction.md` in `data/keys/pdf_dossier_aliases.yml` l. 42. That review
is single-author, has no primary data, and cites **this** paper as its ref. 42 (its dossier l. 105 lists it as not on
disk; it now is, under this stem). **Follow-up paper:** `balagiannis2019_extraction.md` refits these data with maltose.

### Commercial-style shoestring fries (Ranger Russet, 7.5 × 7.5 mm), pre-processed on a pilot line (blanch, glucose OR fructose dip 0-2 %, dry, par-fry 190 °C 50 s, freeze), finish-fried in canola oil at 165 °C (10 dip arms × 10 times, 0-5 min) and 175 / 185 °C (one 0.5 % glucose arm). Glucose, fructose, total free amino acids, acrylamide (13C3-IDA LC-MS/MS) on a defatted dry-weight basis. One multiresponse Wedzicha-type model (Glc + AAs → Int1 → Int2 → acrylamide; Fru → Int2) in a **1-D slab with fitted T and moisture profiles**; **one universal Ea fixed at 100 kJ/mol, not fitted**. Table 2: six constants with 95 % HPD. **No Ea fitted, no acrylamide elimination (returned as zero), no Asn time course, no measured T or moisture curves shown.**

**Read 2026-11-03 from `data/articles/parker2012b.pdf`** (11 pp., ACS born-digital). Prose and Table 1 from
`pdftotext -layout` (Table 1 cross-checked by eye on the 150 dpi render of p. 9323, not separately rendered; values
coherent with Fig. 1 t = 0). **Pp. 9324-9330 rendered at 150 dpi** (`ingest/parker2012b_work/p-04.png` … `p-10.png`):
Table 2 and Fig. 7 (p. 9329), Figs. 5-6 (p. 9328), Figs. 1, 2, 8 digitised (t = 0 and t = 300 s only), read error
±2 mmol/kg (sugars, 0-80 axis), ±10 mmol/kg (amino acids, 0-250 axis), ±3 µmol/kg (acrylamide, 0-100 axis). Eqs. 1-2
read from a 300 dpi crop of p. 9327 (`eq-07.png`). Derived numbers [Z] by hand.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food and Nutritional Sciences, **University of Reading** (Parker, Balagiannis, Mottram),
with ConAgra Foods, Omaha (Higley, Smith; pilot-line processing and analyses) and Leeds (Wedzicha). Corresponding
author **D. S. Mottram**. Existing lab id **`reading_mottram`** ("Reading (Mottram / Parker / Elmore)"). Same-lab
dossiers: `balagiannis2009`, `balagiannis2010`, `balagiannis2015`, `parker2013`, `parker2014`, `goncuoglu2026`,
`kocadagli2021`, and `balagiannis2019` (same team, same pilot line, refits these data). Analytical methods are those
of Higley et al. 2012 JAFC 60:8763 [C] (ConAgra; not on disk).

## 1. What this paper answers for the model

| model need | what Parker 2012 gives | tag |
|---|---|---|
| acrylamide formation rate constants | Table 2: k1 (Glc + AAs → Int1), k2 (Int1 → Int2), k6 (Fru → Int2), k8 / k9 (Glc ⇄ Fru), F_Asn; 95 % HPD | [M] |
| Ea | **one universal Ea = 100 kJ/mol, set, not fitted** ("This was set to 100 kJ/mol" [V]) | [M] assumption |
| T_ref of the k | **not stated**; "kref = rate constant at the frying temperature" [V] | [!] |
| acrylamide elimination | rate constants "consistently returned as zero" [Q]; step removed | [M]/[Q] |
| measured data: Glc, Fru, total AAs, acrylamide vs time | Figs. 1, 2 (165 °C, 10 arms × 11 points 0-300 s), Fig. 8 (175, 185 °C, 1 arm) | [F] |
| asparagine | Table 1 only (stages of pre-processing); R_Asn = Asn/AAs = 0.37 fixed | [M] |
| maltose, sucrose | sucrose Table 1 only, "no decrease" in 5 min [Q]; no maltose | [M]/[Q] |
| surface / core T | thermocouples at centre and just above surface, 165/175/185 °C, every 2 s — **curves not shown**; only the fitted Eq. 1 | [!] |
| moisture vs time | bulk moisture measured, **not shown**; Eq. 2 fitted (A = 57.7 ± 1.4, B = 11.2) | [!] |
| colour | Fig. 9 Agtron vs acrylamide (wet weight), figure only | [F] |
| pH, a_w | **not measured** | [!] |

**Answers.**
1. **Reaction scheme** (Fig. 5, final model; confidence 99 %): glucose + AAs →(k1) Int1 →(k2) Int2; fructose →(k6)
   Int2 (unimolecular, no amino acid); glucose ⇄ fructose (k8 Glc→Fru, k9 Fru→Glc); Int2 is a fast transient: it
   yields "Maillard products" consuming amino acids and acrylamide at rate (k2[Int1] + k6[Fru])·R_Asn·F_Asn (Fig. 6).
   Concentrations in mmol/kg defatted dry weight. Sucrose hydrolysis and acrylamide degradation were tested and
   dropped; a 2-step specific-Asn (Schiff base) model "did not fit" [Q]; a two-parallel-pathway model fitted but the
   generic/specific split was indeterminate [Q].
2. **Fitted constants** (Table 2, 95 % CI): k1 = 0.76 × 10⁻³ mmol⁻¹ kg s⁻¹ (12 %), k2 = 0.40 s⁻¹ (58 %), k6 = 18 × 10⁻³
   s⁻¹ (11 %), k8 = 7.0 × 10⁻³ s⁻¹ (25 %), k9 = 1.1 × 10⁻³ s⁻¹ (152 %, CI spans zero), F_Asn = 5.9 × 10⁻³ (7 %).
   Fit on 165 + 175 + 185 °C data combined. **They are local (surface-referenced) constants inside a slab model with
   (1 − M)² weighting and a T(x, t) integral; they are not bulk rate constants** (confidence ~90 %). Bulk fructose
   falls only ~30 % in 300 s in F5 while k6 alone would leave < 1 % at the surface (e^{−5.4}) [Z].
3. **T range / Ea**: oil 165-185 °C; strip interior < 100 °C [Q]. **The paper gives no information on Ea**: 100 kJ/mol
   for every step is an assumption the authors call validated because the 175/185 °C arm fits (Fig. 8) [Q]. The
   175/185 °C data were in the Table 2 fit, and the "165 °C-only iteration" used for the Fig. 8 prediction has no
   printed parameters (confidence ~85 % that Fig. 8 is not an independent test of Table 2).
4. **Does the lane use any of this? No** (confidence ~97 %). `grep` of `src/kinetic_core/` (all acrylamide modules,
   `parameters_acrylamide.py`) for the DOI, Parker, Balagiannis, Mottram, fries, F_Asn, R_Asn, Int2-of-Wedzicha,
   maltose: no hit tied to this paper (Int1/Int2 in the lane are De Vleeschouwer's). The lane runs mmol/L, min,
   T_ref 160 °C, Arrhenius on Claeys / De Vleeschouwer / Knol; no frying geometry.
5. **Transferable to the lane? Not as constants** (confidence ~85 %): per-kg-dry concentration basis, (1 − M) factors
   whose meaning is unclear (sec. 5), fixed Ea, unstated T_ref, and a 1-D slab transfer model the lane lacks. Usable as
   (a) a **directional check**: no measurable acrylamide elimination during ≤ 5 min finish-frying at 165-185 °C
   (lane includes Knol k6 / Claeys k_E elimination); (b) a future **fry-geometry hold-out** (initial composition +
   time courses are complete for Glc, Fru, AAs, acrylamide); (c) F_Asn = 0.6 % of Asn consumed → acrylamide [M].
6. **Quarantined Parker2012 benchmark** (sec. 6): the README's "nearest real paper" description is **correct**
   (confidence ~97 %); none of the benchmark's conditions or its 1500 ppb value is in this paper.

## 2. Conditions

- **Potato** [M]: Ranger Russet, 2011 harvest, stored 7-9 °C; 6-7 kg batches → 3-4 kg par-fried frozen fries.
- **Pre-processing** [M]: wash, peel, preheat 54 °C water 30 min; cut 7.5 × 7.5 mm (7.5-12.5 cm long); blanch 74 °C,
  10 min (~70 L, fresh water every two batches); dip 30 s at 70-75 °C in glucose (batch G) or fructose (batch F) at 0,
  0.5, 1.0, 1.5, 2.0 % w/v; dry at 60 °C to ~10 % weight loss; par-fry 190 °C, 50 s, ~225 g lots, 20.5 kW Frymaster,
  ~30 L canola oil; freeze −23 °C.
- **Finish-fry** [M]: ~225 g frozen, same fryer, canola oil 165, 175 or 185 °C, 30 s to 5 min (30 s steps), hold-down
  grid; preliminary runs to 15 min (overcooked > 5 min [Q]); quench 40 s in liquid N2; frozen, homogenised.
- **Temperature** [M]: thermocouples in 3 thawed-then-refrozen strips (tip ~centre, 10-12 mm insertion) plus one wire
  "just above the surface"; 165/175/185 °C, 5 min, logged every 2 s. **Data not shown.**
- **175/185 °C arm** [M]: caption says "fries dipped in 0.5 % glucose" but t = 0 composition (Glc ~19, Fru ~6, AAs ~200
  mmol/kg [F]) differs from G2 (27 / 14 / 175) — **a third, undescribed potato lot** [!].
- **Replicates**: not stated for any analyte [!]. pH, a_w: not measured.

## 3. Analytical methods

- **Acrylamide** [M]: 1.0 g + 9.9 mL water + 13C3-acrylamide IS; SPE (Oasis HLB, Bond Elut Accucat); Agilent 1200 /
  AB Q-Trap 4000, ESI+, m/z 72→55 and 75→58. **Isotope dilution (sida).** LOD not given ("nd" in Table 1).
- **Free amino acids** [M]: OPA/Fmoc pre-column derivatisation, Agilent 1100 FLD; norvaline + sarcosine IS; 18 amino
  acids. "Total amino acids" = sum (not stated which; presumably the 18).
- **Sugars** [M]: AOAC 977.20 modified; HPLC-ELSD, Prevail Carbohydrate ES; external standards.
- **Moisture** oven 125 °C to constant weight; **fat** AOAC 996.06 acid hydrolysis, FAME GC vs C11:0 IS. All
  concentrations converted to **defatted dry weight**. **Colour** Agtron E-30 (NIR/green ratio).
- **Modelling** [M]: Athena Visual Studio; Simpson integration over 0-3.75 mm (n = 51); 95 % HPD intervals. The ±
  in Table 2 is the HPD half-width.

## 4. Data

### 4a. Table 1 (p. 9323): pre-processing, glucose-dip batch G, defatted dry weight [M]

| analyte (unit) | stage | raw | 0 % | 0.5 % | 1 % | 1.5 % | 2 % |
|---|---|---|---|---|---|---|---|
| moisture (%) | raw / pre-par-fry / post-par-fry | 78 | 78 / 63 | 77 / 63 | 77 / 62 | 76 / 61 | 77 / 63 |
| glucose (mmol/kg) | raw / pre / post | 30 | 14 / 16 | 30 / 27 | 39 / 36 | 48 / 48 | 66 / 60 |
| fructose (mmol/kg) | raw / pre / post | 23 | 11 / 14 | 13 / 13 | 13 / 14 | 12 / 13 | 16 / 16 |
| sucrose (mmol/kg) | raw / pre / post | 33 | 18 / 24 | 19 / 21 | 18 / 21 | 18 / 21 | 20 / 21 |
| asparagine (mmol/kg) | raw / pre / post | 89 | 64 / 62 | 60 / 60 | 60 / 59 | 54 / 57 | 60 / 58 |
| total amino acids (mmol/kg) | raw / pre / post | 239 | 177 / 166 | 171 / 168 | 172 / 162 | 146 / 155 | 163 / 156 |
| acrylamide (µmol/kg) | raw / pre / post | nd | nd / 1 | nd / 4 | nd / 5 | nd / 8 | nd / 6 |
| fat (g/100 g dry) | raw / pre / post | 0.2 | 0.6 / 21 | 0.6 / 21 | 0.8 / 22 | 0.8 / 20 | 0.6 / 20 |

**No equivalent table for the fructose-dip batch F** [!]. R_Asn from Table 1 post-par-fry: 0.37 / 0.36 / 0.36 / 0.37 /
0.37 [Z] (printed fit values 0.37 G, 0.38 F, combined 0.37, R² 0.96 [M]).

### 4b. Table 2 (p. 9329): final model, all G + F arms, 165/175/185 °C [M]

| parameter | units | estimate × 10³ | 95 % CI × 10³ | lane units [Z] |
|---|---|---|---|---|
| k1 | mmol⁻¹ kg s⁻¹ | 0.76 | ±0.1 (12 %) | 0.0456 kg mmol⁻¹ min⁻¹ |
| k2 | s⁻¹ | 400 | ±230 (58 %) | 24 min⁻¹ |
| k6 | s⁻¹ | 18 | ±2 (11 %) | 1.08 min⁻¹ |
| k8 | s⁻¹ | 7.0 | ±1.7 (25 %) | 0.42 min⁻¹ |
| k9 | s⁻¹ | 1.1 | ±1.7 (152 %) | 0.066 min⁻¹ |
| F_Asn | — | 5.9 | ±0.4 (7 %) | — |

Fixed, not fitted: Ea = 100 kJ/mol (all steps), R_Asn = 0.37. With Ea = 100 kJ/mol a rate at the surface rises ×1.85
from 165 to 175 °C and ×3.31 from 165 to 185 °C; moving T_ref 165 → 160 °C multiplies every k by 0.728 [Z] (only if
T_ref = 165 °C, which is not printed).
Transfer-model constants [M]: Eq. 1 T = T_ref − 780x / (−4.96 × 10⁻⁵ t² + 9.81 × 10⁻² t + 1.1657) (t s, x mm from
surface, T_ref = oil T just above surface); Eq. 2 M = M_ref − 100t / (A e^{Bx}), A = 57.7 ± 1.4, B = 11.2.

### 4c. Figs. 1 and 2, 165 °C, t = 0 → t = 300 s, defatted dry weight [F]

| arm (dip) | Glc mmol/kg | Fru mmol/kg | AAs mmol/kg | Acr µmol/kg | model Acr at 300 s |
|---|---|---|---|---|---|
| G1 (0 %) | 16 → 9 | 14 → 11 | 173 → 157 | 1 → 21 | ~27 |
| G2 (0.5 % Glc) | 27 → 14 | 14 → 10 | 175 → 137 | 3 → 34 | ~40 |
| G3 (1 %) | 36 → 21 | 14 → 13 | 166 → 141 | 4 → 49 | ~50 |
| G4 (1.5 %) | 50 → 24 | 14 → 10 | 163 → 133 | 8 → 60 | ~66 |
| G5 (2 %) | 63 → 27 | 16 → 16 | 162 → 115 | 6 → 78 | ~73 |
| F1 (0 %) | 8 → 6 | 12 → 9 | 190 → 222 | 0 → 8 | ~19 |
| F2 (0.5 % Fru) | 9 → 6 | 29 → 22 | 163 → 161 | 2 → 29 | ~32 |
| F3 (1 %) | 8 → 6 | 36 → 27 | 170 → 187 | 8 → 43 | ~43 |
| F4 (1.5 %) | 12 → 6 | 54 → 38 | 182 → 156 | 9 → 68 | ~62 |
| F5 (2 %) | 12 → 6 | 72 → 50 | 186 → 157 | 11 → 86 | ~78 |

11 time points each (0, 30 … 300 s). Acrylamide shows a ~60 s lag in every arm [F/Q]. AA data scatter ±20-40
mmol/kg point to point (single determinations?) [F]. Text [M]: G5 final acrylamide 4 × G1; final acrylamide linear in
dip concentration (R² 0.995 for both sugars).

### 4d. Fig. 8, 0.5 % glucose arm (third lot), t = 0 → 300 s [F]

| oil T | Glc | Fru | AAs | Acr (µmol/kg) |
|---|---|---|---|---|
| 175 °C | 19 → 8 | 6 → 4 | 197 → 150 | 2 → 29 |
| 185 °C | 19 → ~7 (270 s) | 6 → 3 | 200 → 130 | 2 → 54 |

### 4e. Other figures
Fig. 7 observed vs predicted (all G/F arms): slopes "1 ± 0.02" [M] for Glc, Fru, acrylamide; AAs under-predict loss
(predicted ~150-200 vs observed 110-240 [F]). Fig. 9 Agtron score (100-20) vs acrylamide 0-~4500 µg/kg **wet weight**,
fructose-dipped above glucose-dipped at equal colour [F/Q]; not digitised.

## 5. Flags and internal consistency checks

1. [!] **k1 unit**: Table 2 prints mmol⁻¹ kg s⁻¹; text prints "7.6 × 10⁻⁴ mol⁻¹ kg s⁻¹" [V]. With concentrations in
   mmol/kg the table is right; the text's "mol⁻¹" is a ×1000 typo.
2. [!] **T_ref undefined**: Eq. 4 uses e^{−Ea/R(1/T − 1/Tref)} with T_ref "the frying temperature" (Eq. 1 defines T_ref
   as the oil T just above the surface). Read literally, one k_ref across 165/175/185 °C would make the surface rate
   oil-T-independent; more likely (~60 %) the Arrhenius T_ref was 165 °C, as Balagiannis 2019 states for its refit.
3. [!] **Eq. 1 is unphysical at short times / depth**: at x = 3.75 mm, t = 60 s it gives T ≈ −260 °C; at t = 300 s the
   centre is 53 °C and the 100 °C isotherm sits 2.2 mm deep (165 °C oil) [Z]. Eq. 2 (t in s, M in g/100 g) makes the
   surface dry at ~36 s and the dry front only ~0.19 mm deep at 300 s [Z] — so T > 100 °C extends ~2 mm into still-wet
   tissue, inconsistent with evaporative clamping (~70 % confidence, units of x and t not printed in Eq. 2).
4. [!] "(1 − M)" is described as converting dry-basis to aqueous concentrations [Q]; multiplying mmol/kg dry by the
   dry-matter fraction gives mmol/kg of whole tissue, not an aqueous concentration [Z]. M is in g/100 g, so the factor
   must be (1 − M/100) (presentation slip).
5. [!] **Model misses the zero-dip arms**: predicted acrylamide at 300 s exceeds observed ×1.3 (G1) and ×2.4 (F1) [F];
   the same low-sugar over-prediction recurs in Balagiannis 2019 (M0, ×3.5).
6. [!] Fig. 6 d[AAs]/dt = −(k2[Int1] + k6[Fru]) consumes one AA per Int2, while the text calls amino acids under-
   consumed by the model; Fig. 7 agrees (AA slope not 1).
7. [!] 175/185 °C lot undescribed (sec. 2); Fig. 8 is not independent of Table 2 (sec. 1, ans. 3).
8. [!] Fig. 1 G2 and Fig. 2 F2 acrylamide axes are labelled "mmol/kg dry wt" (all others µmol/kg): label typo.
9. [!] No replicate count or error bars anywhere; batch F pre-processing composition not tabulated.
10. Checks that pass [Z]: Table 1 post-par-fry Glc (16-60) and acrylamide (1-8) agree with Fig. 1 t = 0; 2 % w/v
    glucose = 111 mM, as printed; moisture loss in par-fry 77 → 63 % ("15 %" [M]); Table 2 % CIs recompute (k9 155 %
    vs 152 % printed, rounding).

## 6. Existing records in the repo that cite this paper

`grep -rIn` (2026-11-03) of `data docs src scripts tasks` for `jf302415n`, `parker2012`, `Parker.{0,25}2012`; also
git history of the deleted benchmark (`git show ab9c5cbf^:data/benchmarks/quarantined/acrylamide_asparagine_glucose_Parker2012.json`).

### 6a. `data/benchmarks/quarantined/README.md` (l. 3, 30, 98, 295-320), "nearest real paper" row

| claim | paper | verdict |
|---|---|---|
| Parker et al. (2012), `10.1021/jf302415n` | same | **matches** |
| "a kinetic model of acrylamide in commercial French fries during finish frying" | title and content | **matches** |
| "Not a 10 mM asparagine/glucose aqueous model system" | fries, mmol/kg defatted dry wt, Asn ~60, Glc 8-63 | **matches** |
| "different matrix, different precursor basis, different measurement" | yes on all three | **matches** |
| claimed DOI `10.1021/jf3032779` dead | not this paper's DOI | consistent (not re-checked) |

### 6b. Deleted benchmark `acrylamide_asparagine_glucose_Parker2012.json` (git history; file gone since 2026-10-27)

| field | benchmark | paper | verdict |
|---|---|---|---|
| temp_C | 180 | oil 165 / 175 / 185 °C | **not in paper** |
| time_min | 20 | 0.5-5 min (preliminary to 15) | **not in paper** |
| ph | 6.8 | not measured | **not in paper** |
| water_activity | 0.30 | not measured | **not in paper** |
| L-Asparagine / D-Glucose | 10.0 mM each | Asn 57-62, Glc 8-63 mmol/kg defatted dry wt | **not in paper** |
| acrylamide | 1500 ppb ± 20 % | max ~86 µmol/kg dry (≈ 6.1 mg/kg dry [Z]); Fig. 9 0-4500 µg/kg wet, no condition at 1500 | **not in paper** |
| tier PRIMARY / contract max_ratio 1.05 | — | — | unsupported |
Pattern guess [Z, ~40 %]: pH 6.8 at 180 °C looks like a Knol 2005 phosphate pot, not anything Reading published here.

### 6c. Other records

| record | claim | verdict |
|---|---|---|
| `scripts/generators/generate_prefactor_audit.py` l. 418 | "Parker 2012 (comment: 'Ea_e typically ~90-110 kJ/mol')" for the acrylamide **elimination** Ea | **not in paper**: elimination constants returned zero and the step was removed; the only Ea is a universal 100 kJ/mol set, not fitted, for formation steps. [!] Side note: the generator reads `src/safety.py`, which no longer exists |
| `scripts/generators/generate_gap_heatmap.py` l. 313 | label "Parker 2012 (Asparagine + Glucose)" | **mis-describes the paper** (fries, glucose and fructose dips, total AAs); label for a deleted benchmark |
| `parker2013_extraction.md` l. 105 | ref 42 = Parker, Balagiannis, Higley, Smith, Wedzicha, Mottram, JAFC 60:9321-9331 (2012); on disk "NO" | citation **matches**; on-disk status **stale** (now `parker2012b.pdf`) |
| `docs/history/tasks/audit_remediation.md` l. 152-153 | "real Parker 2012 is a different system (finish-fried …)" | **matches** |
| `data/keys/pdf_dossier_aliases.yml` l. 108 (pending) | Parker et al. 2012 JAFC, acrylamide kinetic model in commercial French fries (Reading) | **matches** |
| `src/kinetic_core/` (acrylamide lane) | — | **no use** of any constant or datum from this paper |

## 7. Proposed ledger entries (`data/lit/measurement_ledger/acrylamide.yml`; NOT written)

```yaml
  - id: parker2012b_tab2_finish_fry_multiresponse_constants
    dossier: parker2012b_extraction.md
    paper: "Parker, Balagiannis, Higley, Smith, Wedzicha & Mottram 2012 J Agric Food Chem 60:9321"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried Ranger Russet shoestring fries (7.5 x 7.5 mm), glucose- or fructose-dipped 0-2 %; Glc 8-63, Fru 9-72, total free AA 160-200, Asn ~60 mmol/kg defatted dry wt"
    medium: "finish-fried in canola oil (pilot Frymaster, ~225 g lots), frozen par-fried strips"
    temperature_c: [165, 175, 185]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [acrylamide, glucose, fructose, total_free_amino_acids]
    kind: rate_constant
    location: "Table 2 (k1, k2, k6, k8, k9, F_Asn with 95 % HPD); scheme Fig. 5, rate laws Fig. 6, slab T/moisture Eqs. 1-2 (dossier sec. 4b)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "One universal value for the activation energy (Ea) could be applied to all kinetic steps"
    flags: [matrix_food, low_moisture_no_aw_term, ph_not_measured, "local_constants_inside_slab_transfer_model", "ea_fixed_not_fitted", no_replicates, "complements: parker2012b_figs1_2_finish_fry_timecourses_165c"]
    note: "Local surface-referenced constants inside a 1-D slab model with (1-M)^2 weighting; Ea FIXED at 100 kJ/mol for every step (not fitted); T_ref not printed (165 C likely); k1 unit typo in text (mol vs mmol). Acrylamide elimination returned zero and dropped. Not bulk kinetics; not transferable to the mmol/L lane without the transfer model."
  - id: parker2012b_figs1_2_finish_fry_timecourses_165c
    dossier: parker2012b_extraction.md
    paper: "Parker, Balagiannis, Higley, Smith, Wedzicha & Mottram 2012 J Agric Food Chem 60:9321"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried fries, 5 glucose-dip (G1-G5) + 5 fructose-dip (F1-F5) arms, 0-2 % w/v; initial composition per arm at t = 0"
    medium: "finish-fried in canola oil at 165 C, quenched in liquid N2"
    temperature_c: [165]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [acrylamide, glucose, fructose, total_free_amino_acids]
    kind: time_series
    location: "Fig. 1 (G arms) and Fig. 2 (F arms), 11 times 0-300 s; t = 0 and 300 s digitised (dossier sec. 4c); Table 1 pre-processing composition (G only)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "were finish-fried at 165 C for different times (30 s intervals up to 5 min)"
    flags: [matrix_food, low_moisture_no_aw_term, heating_up_included_in_time, no_replicates, ph_not_measured]
    note: "Defatted dry-weight basis. Only endpoints digitised; full 11-point series figure-only. Acrylamide ~60 s lag; no acrylamide loss. Batch F pre-processing composition not tabulated. Needs a fry heat/mass-transfer model to simulate."
  - id: parker2012b_fig8_finish_fry_175_185c
    dossier: parker2012b_extraction.md
    paper: "Parker, Balagiannis, Higley, Smith, Wedzicha & Mottram 2012 J Agric Food Chem 60:9321"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "par-fried fries, '0.5 % glucose' dip, undescribed third lot (t = 0: Glc ~19, Fru ~6, AA ~200 mmol/kg dry)"
    medium: "finish-fried in canola oil at 175 and 185 C"
    temperature_c: [175, 185]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0, 5]
    species: [acrylamide, glucose, fructose, total_free_amino_acids]
    kind: time_series
    location: "Fig. 8 (2 T x 11 times); endpoints digitised (dossier sec. 4d)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "fries dipped in 0.5% glucose and cooked at 175 and 185 C"
    flags: [matrix_food, low_moisture_no_aw_term, condition_uncertain, no_replicates]
    note: "Presented as a prediction from a 165 C-only fit, but these data are also in the Table 2 fit; the 165 C-only parameters are not printed. Only T series in the paper (with fixed Ea it cannot identify an Ea)."
  - id: parker2012b_tab1_preprocessing_precursors
    dossier: parker2012b_extraction.md
    paper: "Parker, Balagiannis, Higley, Smith, Wedzicha & Mottram 2012 J Agric Food Chem 60:9321"
    lab: "Reading (Mottram / Parker / Elmore)"
    system: "Ranger Russet strips, raw / pre-par-fry / post-par-fry, glucose dip 0-2 %"
    medium: "pilot line: blanch 74 C 10 min, dip 70-75 C 30 s, dry 60 C, par-fry 190 C 50 s"
    temperature_c: [190]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [0.83]
    species: [acrylamide, asparagine, glucose, fructose, sucrose, total_free_amino_acids, water content]
    kind: endpoint
    location: "Table 1 (dossier sec. 4a)"
    status: digitised_in_dossier
    quant_class: sida
    simulable: false
    condition_quote: "They were dried at 60 C to achieve about a 10% weight loss before par-frying"
    flags: [matrix_food, low_moisture_no_aw_term, no_replicates]
    note: "Only Asn data in the paper (R_Asn = Asn/AA 0.36-0.37). Par-fry forms 1-8 umol/kg acrylamide. Batch F not tabulated."
```

## 8. Supplementary material

None referenced (no Supporting Information paragraph). The measured temperature and bulk-moisture curves exist
[Q] but are published nowhere in this paper; Higley et al. 2012 (JAFC 60:8763) [C] is the companion analytical paper.
