# Zamora, Gallardo & Hidalgo 2007 — EXTRACTION
### *"Strecker Degradation of Phenylalanine Initiated by 2,4-Decadienal or Methyl 13-Oxooctadeca-9,11-dienoate in Model Systems"*
### At 180 °C in air, 2,4-decadienal turns ~8 % of phenylalanine into phenylacetaldehyde in 1 h and 9–14 % of itself into hexanal. The first measured number for the decadienal → hexanal edge the lipid lane declares as a gap. Read 2026-10-27.

**Full read of `data/articles/zamora2007.pdf`** (text layer, 7 pages). Companions:
`zamora2008_extraction.md` (hydroperoxides, same protocol), `zamora2013_extraction.md` (Ea ladder,
including 2,4-decadienal + Phe 27.6 kJ/mol in air), `hidalgo2004_extraction.md` (epoxyalkenals, 37 °C).

**Provenance:** **[M]** measured · **[Q]** qualitative · **[F]** figure only, not digitised · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Rosario Zamora, Emerenciana Gallardo, Francisco J. Hidalgo* |
| Venue | *J. Agric. Food Chem.* 2007, 55 (4), 1308–1314 |
| DOI | 10.1021/jf062838x (PMID 17300152) |
| Affiliation | Instituto de la Grasa, CSIC, Seville; corresponding fhidalgo@ig.csic.es |
| PDF character | 7 pages; Table 1 (retention indices), Schemes 1–2, Figures 1–6 |

---

## 1. THE ONE-PARAGRAPH ANSWER

2,4-Decadienal and the linoleate ketodiene (methyl 13-oxooctadeca-9,11-dienoate) degrade
phenylalanine to **phenylacetaldehyde at 180 °C** but only "very slowly" at 37 °C. Both first oxidise
to the epoxy compound (4,5-epoxy-2-decenal; the epoxyketo ester), detected by GC-MS at acid pH, and
under an inert atmosphere almost nothing forms: **oxygen is required**. Yield of phenylacetaldehyde:
**~8 % after 1 h** from an equimolar decadienal/Phe mixture, **~6 % after 3 h** from the ketodiene,
both at pH 3. Phenylacetaldehyde falls **exponentially with pH** (pH 3–12, r² = 0.96). In the
decadienal system **9–14 % of the decadienal becomes hexanal**, with no clear pH dependence, and
hexanal and phenylacetaldehyde both rise **linearly for the first hour** (zero-order). Co-products:
2-pentylpyridine (from decadienal), methyl 8-(6-pentylpyridin-2-yl)octanoate (from the ketodiene),
benzaldehyde, 2-pentylfuran, pyrroles in trace amounts.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reactants | 0–0.05 mmol 2,4-decadienal (93 %, Aldrich) or methyl 13-oxooctadeca-9,11-dienoate + **0.05 mmol L-phenylalanine** in 1 mL (equimolar = **50 mM** each) |
| Medium | **[!]** Methods: "1 mL of buffer", buffer pre-bubbled 15 min with O2. Every figure caption: "**acetonitrile/buffer (2:1)**". The captions are more specific; treat the medium as 2:1 acetonitrile/buffer until the authors say otherwise |
| Buffers | 0.3 M citrate pH 2.15–6; 0.3 M phosphate pH 6–8 and 11–12; 0.3 M borate pH 8–10; **pH 3** for the dose and time series |
| Vessel / atmosphere | closed Schott Duran tubes 16 × 1.5 cm, **air**; inert-atmosphere runs mentioned, data not shown |
| Temperature / time | **180 °C**; 1 h for the pH and dose series; time course in Figs 5–6 |
| Quantification | phenylacetaldehyde and **hexanal** by GC-MS against 8-level external curves in the injection matrix, internal standard 2-ethylpyridine (r > 0.992); **CV < 10 %** |
| Other measures | yellowness index (YI = 142.86 b*/L*), fluorescence (ex 350 / em 430 nm, quinine sulfate standard), pyrrolization (Ehrlich reagent, ~570 nm) |

---

## 3. RESULTS

| finding | value | provenance |
|---|---|---|
| Phenylacetaldehyde yield, decadienal + Phe, equimolar, pH 3, 180 °C, **1 h** | **~8 %** | [M] Discussion |
| Phenylacetaldehyde yield, ketodiene + Phe, equimolar, pH 3, 180 °C, **3 h** | **~6 %** | [M] Discussion |
| **Hexanal from decadienal** | **9–14 % of the decadienal** "under the assayed conditions" | [M] abstract and Discussion |
| Phenylacetaldehyde vs pH (both lipids) | exponential decrease, pH 3–12 (r² = 0.96) | [F] Fig. 2A |
| Hexanal vs pH | produced "to a high extent" at every pH, no clear correlation | [F] Fig. 2B |
| Dose response (pH 3, 1 h) | phenylacetaldehyde and hexanal linear in lipid amount; browning, fluorescence, pyrrolization linear (r > 0.95) | [F] Figs 3–4 |
| Time course, decadienal (pH 3, 180 °C) | phenylacetaldehyde and hexanal **linear for the first hour** (r > 0.997, p < 0.003), zero-order | [F] Fig. 6 |
| Time course, ketodiene | phenylacetaldehyde linear over the whole period (r = 0.990) | [F] Fig. 6A |
| Without lipid | Strecker degradation "almost absent" | [Q] |
| At 37 °C | decadienal + Phe "occurred very slowly" | [Q] data not shown |
| Inert atmosphere | much less browning and fluorescence; pyrroles and phenylacetaldehyde "only formed in very small amounts" | [Q] data not shown |
| Browning vs pH | slight linear rise (r = 0.81); fluorescence Gaussian, max pH ~8; pyrrolization linear to pH 9 | [F] Fig. 1 |

---

## 4. WHAT THE MODEL CAN TAKE FROM IT

- **`DECADIENAL_TO_HEXANAL_GAP`** (`src/kinetic_core/parameters_lipid.py:600`) says "No rate for it
  exists in the corpus". This paper gives the first measured magnitude: **9–14 % of decadienal to
  hexanal in 1 h at 180 °C, in air, with phenylalanine present, pH 3**, rising linearly over that hour
  (Fig. 6B). Two limits before it becomes a constant: the amino acid is present (part of the hexanal
  comes from pyrrole formation, Scheme 1), and the medium is probably 2:1 acetonitrile/buffer. It does
  show the edge is live at searing temperature, so the gap's warning about the 13-OOH/9-OOH split is
  not hypothetical. Any constant is a new wave, not an edit.
- **Lipid-driven Strecker step:** one yield point per lipid (8 %/1 h, 6 %/3 h) plus the Ea ladder in
  Zamora 2013 would be enough to declare a decadienal + amino acid → Strecker aldehyde step with a
  wide band. The core has no phenylalanine species; for meat aroma the relevant amino acids are
  methionine (methional), leucine and isoleucine (3- and 2-methylbutanal), and whether the same
  yields carry over to them is untested here.
- **Directional claims:** acid pH favours the Strecker yield (both lipids); oxygen is required;
  37 °C is too slow for dienals, unlike the epoxyalkenals of Hidalgo 2004.
