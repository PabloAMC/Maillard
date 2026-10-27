# Hidalgo & Zamora 2004 — EXTRACTION
### *"Strecker-type Degradation Produced by the Lipid Oxidation Products 4,5-Epoxy-2-Alkenals"*
### Lipid oxidation products make the Strecker aldehyde at 37 °C. Qualitative and relative only: no concentration, yield or rate constant is printed. Read 2026-10-27.

**Full read of `data/articles/hidalgo2004.pdf`** (text layer, 6 pages). Read because the deleted prior
`hidalgo_zamora_2004_4hne_pentylpyrrole_v1` cited this paper (see §4).

**Provenance:** **[M]** measured · **[Q]** qualitative · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Francisco J. Hidalgo, Rosario Zamora* |
| Venue | *J. Agric. Food Chem.* 2004, 52 (23), 7126–7131 |
| DOI | 10.1021/jf048883r (PMID 15537327) |
| Affiliation | Instituto de la Grasa, CSIC, Seville |
| PDF character | 6 pages; Schemes 1–3, Figures 1–4; no tables |

---

## 1. THE ONE-PARAGRAPH ANSWER

4,5-Epoxy-2-alkenals, secondary lipid oxidation products (4,5(E)-epoxy-2(E)-heptenal from n−3 and
4,5(E)-epoxy-2(E)-decenal from n−6 fatty acids), react with phenylalanine at **37 °C** to give the
Strecker aldehyde **phenylacetaldehyde**, alongside the main products, N-substituted
2-(1-hydroxyalkyl)pyrroles and N-substituted pyrroles. Phenylacetaldehyde is **absent with the
methyl ester** (the free carboxyl is needed), is confirmed as formed at 37 °C (not in the GC injector)
by cysteamine trapping, and peaks at **pH ~6.9** (Gaussian, r² = 0.92). The same reaction gives
2-ethylpyridine / 2-pentylpyridine. **[Q]** The paper prints no absolute amounts: Figure 3 is a
phenylacetaldehyde/internal-standard GC-MS area ratio against pH.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reactants | 0.1 mmol epoxyalkenal (5a or 5b) + 0.1 mmol L-phenylalanine (6) or its methyl ester (6M) |
| Medium | 1 mL acetonitrile–water (2:1), or 1 mL buffer: 0.3 M citrate pH 4, 5, 5.5, 6; 0.3 M phosphate pH 6, 6.5, 7, 7.4, 8; 0.3 M borate pH 8, 9, 10 |
| Concentration | 0.1 M each (0.1 mmol in 1 mL) |
| Temperature / time | **37 °C**, "different incubation times"; figures are **overnight** incubations |
| Analysis | GC-MS (HP 6890 / 5973 MSD, HP5-MS); internal standard 3(Z)-nonenol (1.35 mg/mL, 25 µL); aldehydes also as cysteamine thiazolidines (Yasuhara et al. 1998) |
| Quantification | **none absolute**: analyte/IS area ratios only (Figure 3) |
| Browning | ΔE from CIELAB, aqueous and chloroform extracts (Figure 4) |

---

## 3. RESULTS

| finding | value | provenance |
|---|---|---|
| Phenylacetaldehyde vs pH | Gaussian, centre **pH 6.9**, r² = 0.92; "mainly produced at pH 6-7" | [M] Fig. 3 (area ratios) |
| Phenylacetaldehyde with Phe methyl ester | **absent** (both epoxyalkenals) | [Q] Fig. 1B, 1D |
| Formed at 37 °C, not in the injector | yes: thiazolidine 14T detected after cysteamine trapping | [Q] Fig. 2 |
| Co-products | propanal (5a) / hexanal (5b), 2-ethyl- / 2-pentylpyridine, 2-pentylfuran (5b), 2,4-heptadienal / 2,4-decadienal, 2-octenal (5b) | [Q] Figs 1–2 |
| Browning, aqueous extract | linear in pH 5–10 (r = 0.994, p < 0.0001) | [M] Fig. 4A |
| Browning, chloroform extract | Gaussian, centre pH 6.4 (r² = 0.96); correlates with phenylacetaldehyde (r = 0.92, p = 0.00013) | [M] Fig. 4B |
| Mechanism | imine → decarboxylation → hydrolysis (Scheme 3), analogous to α-dicarbonyl Strecker degradation; the hydroxyl-amino intermediate 13 gives the 2-alkylpyridines | [Q] proposed |

---

## 4. WHAT THE REPOSITORY CLAIMED, AND WHAT THE PAPER SAYS [!]

The prior `hidalgo_zamora_2004_4hne_pentylpyrrole_v1` (deleted 2026-10-27 with 36 other unanchored
priors) described this citation as **"4-HNE + phenylalanine"**, effect "flag_4hne_derived_pyrrole_risk",
with `absolute_concentrations_reported: true` and `loq_reported: true`. All three are wrong for this
paper: the reactant is a **4,5-epoxy-2-alkenal, not 4-hydroxynonenal**; the headline product is the
**Strecker aldehyde**, not a pentylpyrrole risk (a 2-pentyl-1-phenethylpyrrole is only "tentatively
identified" among side products); and **no absolute concentration or LOQ is reported**. The prior
held no number, so nothing numeric was lost. Its intake entry
(`hidalgo_zamora_2004_4hne_pentylpyrrole`, which cited Globisch, Kaden & Henle 2015) was deleted too.

---

## 5. WHAT THE MODEL CAN TAKE FROM IT

- **Nothing quantitative.** No rate constant, no yield, no absolute concentration; 37 °C only; organic
  co-solvent for most runs.
- **One qualitative claim for the directional panel, if wanted:** a lipid oxidation product
  (epoxyalkenal) Strecker-degrades phenylalanine to phenylacetaldehyde, maximal near pH 7, and the
  free carboxyl is required. The kinetic core has no epoxyalkenal species and no lipid-driven
  Strecker step, so the claim would be refused today, not scored.
- **Relevance to plant-based meat:** a route by which oil oxidation could supply Strecker aldehydes
  (malty, honey, cocoa notes) without sugar-derived dicarbonyls. Strength at cooking temperatures is
  untested here.
