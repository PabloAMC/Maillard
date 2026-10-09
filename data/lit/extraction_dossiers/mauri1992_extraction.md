# Mauri, Alzamora & Tomio 1992 — EXTRACTION (first-order thiamine loss at aw 0.95, pH 4.0 and 5.5, 80-100 °C, with NaCl, KCl or Na2SO4 as humectant)

**Source on disk:** `data/articles/mauri1992.pdf` (a 2003 scan, 5 pages, with a poor OCR text layer).
Every number below was read from the page images and checked by eye on 2026-10-09; Table 1 was
additionally re-rendered at 300 dpi and read from that crop. **The OCR text layer is unreliable for this
paper and was NOT used as a source**: it garbles most of Table 1 (e.g. it renders the pH 5.5 Ea column
"27·4, 27·4, 27·7" as "21.4, 21.4, 21.1", and "0·047" as "OW7") and the abstract's "aw = 0·95" as "0.93".
Where the text layer is legible it agrees with the images (e.g. 0.093, 0.035, 0.014, 0.0262, 28.1). The
journal prints decimals with a raised dot (0·093); they are written 0.093 below. Page numbers are journal
pages (PDF page = journal page − 18). Written for the thiamine route of the sulfur lane
(`src/kinetic_core/sulfur.py`, `r_thi_hmp` / `r_thi_mesh`), whose activation energy is borrowed.

| field | value |
|---|---|
| Title | "Effect of electrolytes on the kinetics of thiamine loss in model systems of high water activity" |
| Authors | L. M. Mauri, S. M. Alzamora, J. M. Tomio (Universidad de Buenos Aires) |
| Venue | Food Chemistry 1992, 45, 19-23 |
| DOI | 10.1016/0308-8146(92)90006-N |

## 1. Methods (pp. 19-20)

- **Thiamine form:** thiamine **hydrochloride** (Sigma), **0.1 mg/100 ml** = 1 mg/L = **2.97 µM**
  (derived here: 1 mg/L / 337.27 g/mol). Mononitrate not used.
- **Buffer:** phosphate (KH₂PO₄ and Na₂HPO₄·2H₂O, **0.07 M**) "to maintain a constant pH (4.0 or 5.5)
  during the runs".
- **Water activity:** adjusted to **aw 0.95** with NaCl, KCl or Na₂SO₄ ("electrolyte concentration between
  7.7 and 17.2% w/w", p. 19); concentrations computed with Pitzer & Mayorga, aw estimated by the Ross
  equation, checked with a Novasina Thermoconstanter hygrometer. The **control ("Buffer") has no
  humectant**; its aw is not printed in Table 1 (the text on p. 22 calls the humectant-free system
  "aw ≅ 1"), although Table 1's title says "(aw = 0·95)" for all rows.
- **Heating:** screw-cap glass tubes at 100, 90 or 80 °C in a shaken water-glycerol bath; some runs at
  lower temperature ("mainly at 55°C") in a forced-convection oven. Oxygen/headspace not described.
- **Assay:** every sample by the **thiochrome fluorometric method** (HgCl₂ oxidant; Ryan & Ingle 1980;
  ex 365 / em 435 nm). For some systems (pH 5.5, aw 0.95, 100 °C, NaCl or KCl) also by HPLC (carbohydrate
  column, acetonitrile-NH₄HCO₃ 80:20, UV 254 nm); rate constants from the two methods agreed within
  "10-12%" relative error. Both measure intact thiamine: the quantity is **total loss of thiamine by all
  routes**; no product is measured.
- **Reaction order:** first order ("Good straight lines were obtained in all cases", p. 20; Fig. 1 for
  pH 5.5; pH 4.0 "results not shown"). k and Ea "and corresponding errors, were calculated by regression
  analysis" (p. 20). **The paper does not say whether Δk and ΔEa are standard errors, confidence
  half-widths or something else** (it cites a statistics text, Boquet 1984, without detail). No R².

## 2. Findings that matter

### Table 1 (p. 21), "Kinetic parameters for thiamine degradation in model systems (aw = 0·95)"

Units exactly as printed: k and Δk in **h⁻¹**; Ea and ΔEa in **kcal mole⁻¹**. The kJ/mol column is
**derived here** (× 4.184).

| pH | system | k 100 °C (h⁻¹) | Δk | k 90 °C (h⁻¹) | Δk | k 80 °C (h⁻¹) | Δk | Ea (kcal mole⁻¹) | ΔEa | Ea (kJ/mol), derived here |
|---|---|---|---|---|---|---|---|---|---|---|
| 5.5 | Buffer | 0.093 | 0.004 | 0.035 | 0.003 | 0.014 | 0.001 | 27.4 | 2 | 114.6 ± 8.4 |
| 5.5 | NaCl | 0.065 | 0.004 | 0.027 | 0.003 | 0.0087 | 0.0005 | 27.4 | 0.5 | 114.6 ± 2.1 |
| 5.5 | KCl | 0.067 | 0.004 | 0.021 | 0.001 | 0.0080 | 0.0008 | 27.7 | 0.4 | 115.9 ± 1.7 |
| 5.5 | Na₂SO₄ | 0.069 | 0.006 | 0.0262 | 0.0009 | 0.0081 | 0.0004 | 28.1 | 0.3 | 117.6 ± 1.3 |
| 4.0 | Buffer | 0.047 | 0.002 | 0.018 | 0.003 | 0.0065 | 0.0005 | 29 | 3 | 121.3 ± 12.6 |
| 4.0 | NaCl | 0.024 | 0.001 | 0.0107 | 0.0003 | 0.0031 | 0.0001 | 29 | 1 | 121.3 ± 4.2 |
| 4.0 | KCl | 0.0281 | 0.0009 | 0.0098 | 0.0005 | 0.0030 | 0.0001 | 29 | 2 | 121.3 ± 8.4 |
| 4.0 | Na₂SO₄ | 0.044 | 0.005 | 0.0170 | 0.0008 | 0.0052 | 0.0002 | 29 | 3 | 121.3 ± 12.6 |

Abstract (p. 19) and text (p. 22): Ea "ranged from 27 to 29 kcal mole⁻¹" (= 113-121 kJ/mol, derived here),
"independent of the humectant ... as well as of pH". Text (p. 22): at pH 5.5 the electrolyte systems' k
is "about 30-40% lower" than the buffer; at pH 4.0 NaCl/KCl are "40-50% lower" and Na₂SO₄ "from 6 to 9%"
lower. Lower pH gives greater retention (p. 20).

**Which points enter the Ea regression is not stated.** Fig. 2 (p. 21) plots ln(k × 10⁴) vs 1/T × 10⁴
including k values at 55 °C for the pH 5.5 systems (and, for Na₂SO₄, literature points at 45-65 °C from
Fernández 1984), but these 55 °C k values are **not tabulated**. A three-point refit of the tabulated
80/90/100 °C k (derived here, unweighted least squares on ln k vs 1/T) gives 103.7 (Buffer), 110.3
(NaCl), 116.3 (KCl), 117.4 (Na₂SO₄) kJ/mol at pH 5.5, and 108.4, 112.3, 122.6, 117.1 kJ/mol at pH 4.0.
The KCl and Na₂SO₄ refits reproduce the printed Ea; the Buffer and NaCl refits fall 4-11 kJ/mol below it,
consistent with the printed values including extra (e.g. 55 °C) points, but the paper does not say.

### Derived here: k at 100 °C (measured) and 140 °C (extrapolated)

k(100 °C) is measured (Table 1). k(140 °C) = k(100 °C) · exp[−(Ea/R)(1/413.15 − 1/373.15)] with the
printed Ea × 4.184 kJ/mol, R = 8.314 J mol⁻¹ K⁻¹. **140 °C is 40 °C above the highest measured
temperature, a long extrapolation.** The paper prints no pre-exponential factor.

| pH | system | k(100 °C), measured | k(100 °C) in min⁻¹ | t½ at 100 °C | k(140 °C), derived | k(140 °C) in min⁻¹ |
|---|---|---|---|---|---|---|
| 5.5 | Buffer | 0.093 h⁻¹ | 0.00155 | 7.5 h | 3.33 h⁻¹ | 0.0555 |
| 5.5 | NaCl | 0.065 h⁻¹ | 0.00108 | 10.7 h | 2.33 h⁻¹ | 0.0388 |
| 5.5 | KCl | 0.067 h⁻¹ | 0.00112 | 10.3 h | 2.49 h⁻¹ | 0.0416 |
| 5.5 | Na₂SO₄ | 0.069 h⁻¹ | 0.00115 | 10.0 h | 2.71 h⁻¹ | 0.0451 |
| 4.0 | Buffer | 0.047 h⁻¹ | 0.00078 | 14.7 h | 2.07 h⁻¹ | 0.0346 |
| 4.0 | NaCl | 0.024 h⁻¹ | 0.00040 | 28.9 h | 1.06 h⁻¹ | 0.0176 |
| 4.0 | KCl | 0.0281 h⁻¹ | 0.00047 | 24.7 h | 1.24 h⁻¹ | 0.0207 |
| 4.0 | Na₂SO₄ | 0.044 h⁻¹ | 0.00073 | 15.8 h | 1.94 h⁻¹ | 0.0323 |

Worked example (pH 5.5 Buffer): Ea = 27.4 × 4.184 = 114.6 kJ/mol; exp[−(114642/8.314)(1/413.15 −
1/373.15)] = exp(3.578) = 35.8; 0.093 × 35.8 = 3.33 h⁻¹ = 0.0555 min⁻¹.

## 3. What it means for the model

- **Total-thiamine-loss Ea, phosphate-buffered, 80-100 °C: 27.4-29 kcal/mol = 114.6-121.3 kJ/mol**
  (derived conversion), at pH 4.0 and 5.5, with or without salt. That is within the 55-145 kJ/mol band
  and ~15-20 kJ/mol above its ~100 centre; it is ~51-57 kJ/mol above the 64.08 kJ/mol the live engine
  actually gives `k_thi_hmp` (frozen lumped formation Ea, `kinetic_core_b9_fit_report.json`, rate anchored
  at 145 °C). It agrees with Ramaswamy 1990's water series (102.6-118.0
  kJ/mol, 110-150 °C) but not with Ramaswamy's glucose/glycine/ascorbate mixture (71-73 kJ/mol).
- pH matters more than Ea: k(100 °C) in buffer roughly doubles from pH 4.0 (0.047 h⁻¹) to pH 5.5
  (0.093 h⁻¹). Meat pH (~5.5-6) is at the faster end of this paper's range; nothing above 5.5 is measured.
- Salt at aw 0.95 slows loss by 30-50 % (except Na₂SO₄ at pH 4.0); it does not change Ea within the
  printed errors.
- The measured quantity is **total thiamine loss** (thiochrome, HPLC cross-check). The model's thiamine
  sinks are branches of it; their summed k can be capped by these numbers, and branch Ea = total Ea only
  if the branching ratio does not depend on temperature, which is not tested here.
- Cross-paper check (derived here): Mauri pH 5.5 buffer k(100 °C) 0.00155 min⁻¹ vs Ramaswamy water
  k(100 °C) extrapolated 0.00078-0.00135 min⁻¹ (unbuffered, pH not printed); same order of magnitude.

## What it does not give

- Nothing above 100 °C (except by extrapolation) and nothing above pH 5.5; no meat or tissue matrix.
- No product measurement (no HMP, MFT, H₂S or volatiles), no branching.
- Thiamine at 3 µM only; no concentration dependence, no mononitrate, no phosphorylated forms.
- The statistic behind Δk and ΔEa is not defined; no R²; no individual retention data tabulated (Fig. 1
  graph only, pH 5.5).
- The 55 °C k values shown in Fig. 2 are not tabulated, and it is not stated whether they enter the
  printed Ea.
- The buffer control's actual aw is not printed (Table 1's title applies "aw = 0·95" to all rows; the text
  implies the control is ≅ 1).
- Fig. 3 (p. 22: pH 4.0, 100 °C, non-electrolyte humectants glycerol, sorbitol, propylene glycol, sucrose,
  glucose) is cited to "Mauri et al. (1990) ... (submitted to Food Chem.)", whose reference-list title
  is identical to this paper's own title; the non-electrolyte k values are not tabulated here.
