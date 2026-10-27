# Zamora, Gallardo & Hidalgo 2008 — EXTRACTION
### *"Model Studies on the Degradation of Phenylalanine Initiated by Lipid Hydroperoxides and Their Secondary and Tertiary Oxidation Products"*
### At 180 °C, linoleate and linolenate hydroperoxides turn phenylalanine into phenylacetaldehyde: up to 12 % at pH 2.15, 0.25 % at pH 12. Two yields printed; the rest is in figures. Read 2026-10-27.

**Full read of `data/articles/zamora2008.pdf`** (text layer, 6 pages). Companion to
`zamora2013_extraction.md` (Ea ladder across lipid carbonyl classes) and `hidalgo2004_extraction.md`
(epoxyalkenals at 37 °C).

**Provenance:** **[M]** measured · **[Q]** qualitative · **[F]** figure only, not digitised · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Rosario Zamora, Emerenciana Gallardo, Francisco J. Hidalgo* |
| Venue | *J. Agric. Food Chem.* 2008, 56 (17), 7970–7975 |
| DOI | 10.1021/jf801409w (PMID 18707112) |
| Affiliation | Instituto de la Grasa, CSIC, Seville; corresponding fhidalgo@ig.csic.es |
| PDF character | 6 pages; Table 1 (retention indices), Figures 1–5 |

---

## 1. THE ONE-PARAGRAPH ANSWER

Primary (hydroperoxides MeLOOH, MeLnOOH), secondary (hydroxydiene MeLOH, ketodiene MeLCO) and
tertiary (epoxyhydroxy MeLEPOH, epoxyketo MeLEPCO) linoleate-derived oxidation products all degrade
phenylalanine to **phenylacetaldehyde at 180 °C**. The yield **falls exponentially with pH, from 12 %
at pH 2.15 to 0.25 % at pH 12** (1 h, equimolar). It rises **linearly for the first 30 min** and then
plateaus. Reactivity in air after 1 h: **MeLOOH < MeLOH ≤ MeLEPOH < MeLCO ≤ MeLEPCO**, a small
spread; under nitrogen every product except MeLEPCO loses most of its activity. n−6 and n−3
hydroperoxides behave alike. The proposed route is oxidative interconversion towards the epoxyketo
ester MeLEPCO, plus fragmentation to short-chain aldehydes.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reactants | 0–25 µmol lipid derivative + **25 µmol phenylalanine** in 0.5 mL buffer (equimolar = **50 mM** each) |
| Lipids | methyl esters, prepared by lipoxygenase (hydroperoxides), NaBH4 reduction (MeLOH), CrO3 oxidation (MeLCO), mCPBA epoxidation (MeLEPCO), NaBH4 (MeLEPOH); purity by TLC, NMR, GC |
| Vessel | closed Schott Duran tubes 16 × 1.5 cm, **air** headspace unless "nitrogen" |
| Temperature / time | **180 °C**; 1 h for the pH, dose and reactivity series; time course in Fig. 3 |
| Buffers | 0.3 M citrate pH 2.15–6.0; 0.3 M phosphate pH 6.0–8.0 and 11.0–12.0; 0.3 M borate pH 8.0–10.0. **pH 3 for everything except Fig. 1**, chosen to match the group's earlier work |
| Quantification | phenylacetaldehyde by GC-MS against an **8-level external calibration** in the injection matrix, internal standard methyl heptanoate; response linear (r = 0.999); **CV < 10 %**; means of **at least two** independent experiments |

---

## 3. RESULTS

| finding | value | provenance |
|---|---|---|
| Yield vs pH (MeLOOH or MeLnOOH + Phe, equimolar, 180 °C, 1 h) | exponential decrease (r² > 0.96): **12 % at pH 2.15 → 0.25 % at pH 12** | [M] text; curve [F] Fig. 1 |
| n−6 vs n−3 hydroperoxide | "very similar behaviors" | [Q] Fig. 1, 3A |
| Dose response, hydroperoxide | non-linear, Boltzmann fit (r² > 0.97); largest gains at small doses | [F] Fig. 2A |
| Dose response, MeLOH and MeLCO | linear (r > 0.97, p < 0.001) | [F] Fig. 2B |
| Time course, hydroperoxides | linear for the first 30 min (r > 0.98, p < 0.02), then a slowly rising plateau | [F] Fig. 3A |
| Early products | at 5–10 min the main product is MeLOH; later, phenylacetaldehyde | [Q] |
| Reactivity, air, 1 h | MeLOOH < MeLOH ≤ MeLEPOH < MeLCO ≤ MeLEPCO ("a small increase") | [F] Fig. 4 |
| Reactivity, nitrogen, 1 h | MeLEPCO ≫ MeLEPOH > MeLOH ≈ MeLCO > MeLOOH; MeLEPCO unchanged, the rest "much reduced" | [F] Fig. 4 |
| Co-products identified (Table 1) | hexanal, 2-hexenal (n−3), benzaldehyde, 2-pentylfuran, 2-octenal, 2,4-decadienal, methyl octanoate, methyl 8-oxooctanoate, methyl 9-oxononanoate, 1-phenethylpyrrole, 2-pentyl-1-phenethylpyrrole, MeLOH, MeLCO | [Q] retention indices |

**[!] Wording mismatch:** the abstract orders reactivity "MeLEPCO ≥ MeLCO > MeLEPOH ≥ MeLOH > MeLOOH ≈
MeLnOOH"; the Results order MeLEPOH below MeLCO with "≤" between MeLOH and MeLEPOH. Same ranking,
different inequality signs; Figure 4 is the source.

---

## 4. WHAT THE MODEL CAN TAKE FROM IT

- **Two printed yields** (12 % at pH 2.15; 0.25 % at pH 12, 180 °C, 1 h, 50 mM equimolar). Everything
  else is in Figures 1–4 and would need digitising before any row is declared.
- **No rate constant and no activation energy.** The temperature dependence of lipid-driven Strecker
  degradation is in `zamora2013_extraction.md` (carbonyl classes, not hydroperoxides).
- **Directional claims it supports:** acid pH favours the reaction across the whole lipid oxidation
  cascade; oxygen is needed for every product except the epoxyketo ester; n−6 and n−3 hydroperoxides
  behave alike.
- **Relevance to plant-based meat:** linoleate dominates pea, soy and sunflower oils, so this is the
  route by which their hydroperoxides could make Strecker aldehydes during searing. Caveat: at the
  pH of a plant-based patty (~6) the yield is far below the pH 3 condition most of the paper uses
  (Fig. 1); read the pH 6 point off the figure before quoting a number.
- The kinetic core has no lipid-driven Strecker step and no phenylalanine species today, so none of
  this can be scored yet.
