# Xin, Liu, Cui, Nie, Huang & Qin 2026 — EXTRACTION (eight carbohydrates at 6 % w/w in a high-moisture extruded plant-based meat analogue: texture, structure and volatiles)

**Source on disk:** `data/articles/Xin2026.pdf` (10.0 MB; downloaded before 2026-09-07). Read 2026-09-11
via `pdftotext -layout` and `pypdf` (86 k characters, full text layer). Dossier written because the
2026-09-11 audit found this PDF was the ONLY one on disk with no dossier of its own that was not
already declared "no dossier" by another (its companion `Xin2026b_extraction.md` refers to an
`Xin2026_extraction.md` that did not exist; `k2_matrix_and_thresholds.md` and
`k3_final_parameter_inventory.md` cite "Xin 2026" for the HS-SPME same-sample dispersion fact).

| field | value |
|---|---|
| Title | "Effects of different carbohydrates on the textural, structural, and flavor properties of high-moisture extruded plant-based meat analogs" |
| Venue | Food Hydrocolloids 182 (2027) 113124; received 12 May 2026, accepted 9 July 2026, online 9 July 2026 |
| DOI | 10.1016/j.foodhyd.2026.113124 |
| Group | Dalian Polytechnic University (Lei Qin, Xuhui Huang) — the same laboratory as `Xin2026b` |
| Systems | HME-PBMA base (protein blend; see sec. 2.1) with **6 % (w/w)** of one of eight carbohydrates: glucose (GL), xylose (XY), ribose (RI), fructose (FR), β-glucan (BG), wheat starch (WS), maltodextrin (MA), sucrose (SU); plus a control with none |
| What is measured | rheology of the pre-extrusion blend, texture, fibrous degree, colour (Table 1), cooking yield, LF-NMR water distribution, macro/micro structure (SEM), FTIR secondary structure, DSC, E-nose, HS-SPME-GC-MS volatiles (relative concentrations against an internal standard, µg/kg), OAV against air thresholds, sensory |

## 1. What this paper is and is not, for this model

It is a **matrix/formulation study**, not a kinetic one: one extrusion condition, one time point,
no rate constant, barrier, order or time course anywhere. Its value to the repository is
**directional** (which sugar makes which volatile class in a real extrudate) and, at most, as
external-matrix levels for a benchmark family this repository already refuses to score on absolute
terms (HME extrudates, `external_validation_li_2026_spi_wg_hme_control`).

**Only ONE table in the paper prints numbers: Table 1, the colour coordinates.** Every volatile
concentration is in Figures 6–7 (a PLS-DA score plot, a biplot, a VIP list, a 45-compound heatmap,
class-sum bars and "key compound" bars) and in running prose. The prose values below are the only
volatile numbers that can be transcribed without digitising a figure.

## 2. Table 1 — colour (verbatim, n = 3, letters = Duncan groups)

| sample | L* | a* | b* | ΔE |
|---|---:|---:|---:|---:|
| Control | 57.60 ± 0.10 | −0.39 ± 0.09 | 8.81 ± 0.07 | 43.41 ± 0.10 |
| GL | 55.70 ± 0.08 | 0.08 ± 0.02 | 10.43 ± 0.14 | 45.62 ± 0.22 |
| XY | 58.49 ± 0.05 | −0.80 ± 0.07 | 6.91 ± 0.09 | 42.19 ± 0.03 |
| RI | 58.80 ± 0.17 | −0.67 ± 0.02 | 7.37 ± 0.22 | 41.96 ± 0.12 |
| FR | 56.53 ± 0.08 | 0.55 ± 0.03 | 10.46 ± 0.06 | 44.82 ± 0.06 |
| BG | 60.06 ± 0.03 | −0.07 ± 0.04 | 8.32 ± 0.06 | 40.90 ± 0.50 |
| WS | 59.64 ± 0.05 | −0.35 ± 0.04 | 8.89 ± 0.07 | 41.43 ± 0.04 |
| MA | 58.82 ± 0.09 | −0.42 ± 0.02 | 7.37 ± 0.13 | 41.94 ± 0.06 |
| SU | 57.98 ± 0.10 | 0.02 ± 0.00 | 8.42 ± 0.09 | 42.97 ± 0.07 |

The browning ordering the model would be asked about: **FR ≈ GL > SU > control > pentoses (XY, RI)
> polysaccharides** by ΔE. Note the pentoses BROWN LESS than the hexoses here, the opposite of the
aqueous model-system ordering the trunk was fitted on — an extrudate at high moisture and a few
seconds' residence is not a 120 °C aqueous hold, and the sugars are 6 % of a protein matrix.

## 3. Volatile numbers printed in the prose (everything else is figure-only)

| quantity | value | where |
|---|---|---|
| Maillard-derived pyrazine total, FR extrudate (the maximum) | **6621.64 µg/kg** | abstract |
| 2-pentylfuran, polysaccharide-added extrudates | **100 259.10 and 97 794.41 µg/kg** (two of BG/WS/MA) | sec. 3.11 |
| 2-pentylfuran, monosaccharide-added extrudates (GL, RI, FR) | **50 230.59–77 137.05 µg/kg** | sec. 3.11 |
| aldehyde class total, highest (one polysaccharide extrudate) | 3506.49 µg/kg; FR 3483.59; MA 3381.41; two others 2433.19 and 2346.88 | sec. 3.11 |
| hexanal, control extrudate (the minimum) | 922.04 µg/kg; GL 1698.56; FR (number cut at page break) | sec. 3.11 |
| compounds with OAV > 1 | 16: three ketones, six aldehydes (hexanal, benzaldehyde, nonanal, heptanal, octanal, decanal), two furans (2-butylfuran, 2-pentylfuran), three alcohols, two pyrazines (2-ethyl-3,6-dimethylpyrazine, 2,3-dimethyl-5-propylpyrazine) | sec. 3.11 |
| pyrazines identified | 2,5-dimethylpyrazine, 2-ethylpyrazine, 2-methyl-5-(1-methylethyl)pyrazine, 2-ethyl-3,6-dimethylpyrazine, 2,3-dimethyl-5-propylpyrazine | sec. 3.11 |

**Quantification caveat, verbatim in kind:** concentrations are "relative concentrations (µg/kg)"
against one internal standard with no response factors — the same class of number as Frankel 1981's
peak areas and Whitfield 1999's unit-response-factor assumption. They carry the HS-SPME dispersion
band this repository already applies (10–23× same-sample, `k2_matrix_and_thresholds.md`).

## 4. What the repository takes

- **Two directional claims** (sec. 3.11): (i) monosaccharides raise the pyrazine total over the
  control and polysaccharides, with fructose highest; (ii) monosaccharide addition LOWERS
  2-pentylfuran against polysaccharide addition (50–77 vs 98–100 thousand µg/kg). The second is the
  more useful one for this model: it is a lipid-lane observable moved by a sugar, which the model
  cannot represent (the lipid lane charges a hydroperoxide pool and knows nothing about sugar), so it
  goes into the directional panel as a NOT EVALUABLE claim with its reason, not as a miss.
- **No rate, no fit row, no hold-out row.** A single extrudate at one condition with relative
  concentrations is not a level this model should be scored on; the HME family is already refused on
  absolute terms.
- The colour table is the only thing here a browning lane could ever be checked against, and this
  repository has no browning observable on the panel.

## 5. Flags

1. Volatile concentrations are **figure-only** except the prose values above; Figure 7B's "key
   compound" bars would need digitising.
2. The base protein blend and the extrusion profile are in sec. 2.1–2.2 and were not transcribed
   here; if this paper is ever used as a benchmark, they are the conditions block.
3. The same group's `Xin2026b` (companion, already extracted) shares methods; do not double-count
   any claim across the two.
