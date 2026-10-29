# Barallat-Perez 2023 — EXTRACTION (commercial pea, soy, lupin and whey protein isolates 1 % w/v in 50 mM sodium phosphate pH 7.0, ten carbonyl/alcohol flavours at 5 mg/L each, held 3 h at 30 / 70 / 90 °C then stored at 5 °C and measured weekly; % bound by static-headspace GC-MS against a protein-free control; surface hydrophobicity H0 and sulfhydryl content of the four isolates)

### THE "SEALED LUPIN" PARENT STUDY, READ IN FULL — and it is not the parent of the sealed keys either: it prints **no binding constant (only % bound, eq 1)**, **never measures nonanal** and **has no mucin arm**, while `HOLDOUT_SEALED_BINDING` registers `kg_*` (L/g) constants for hexanal, nonanal and 2-nonanone on lupin and on mucin. Its lupin numbers are bars in Figs. 3–4; its printed numbers are pooled PPI/SPI/LPI ranges, fat contents, H0 (Table 1) and –SH (Table 2).

**Source on disk:** `data/articles/barallatperez2023.pdf` (11 pp., J. Agric. Food Chem. 2023, 71, 20274−20284, CC-BY 4.0;
DOI 10.1021/acs.jafc.3c05991; PMCID PMC10739987). Read from the `pdftotext -layout` text layer on 2026-10-29. Tables 1 and 2
came through the text layer intact. Figures 1–4 are bar charts (images): **every per-protein binding percentage is
figure-only** and is not digitised here. Supporting Information (Tables S1–S4, Figure S1) is **not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Unraveling the Role of Flavor Structure and Physicochemical Properties in the Binding Phenomenon with Commercial Food Protein Isolates" |
| Authors | Cristina Barallat-Pérez (corresponding), Hans-Gerd Janssen (also Unilever Foods Innovation Centre), Sara Martins (also AFB International EU), Vincenzo Fogliano, Teresa Oliviero (corresponding) — Agrotechnology and Food Science, Wageningen |
| Journal | J. Agric. Food Chem. 2023, 71(50), 20274−20284; received Aug 25, revised Nov 9, accepted Nov 21, published Dec 7, 2023 |
| Repo keys | `data/lit/binding_constants.yml` source `barallat_perez_2023_jafc`; `data/keys/papers.yml` paper_id `acs_jafc_3c05991_ppi_spi_partitioning` (dossier: null until now) |
| Companions | Barallat-Pérez 2024 JAFC 72:8731 (lupin in-mouth, mucin; `barallatperez2024_extraction.md`); Barallat-Pérez 2025 npj Sci Food 9:174 (mucin, Klotz constants for SPI/WPI/LPI; `binding_constants.yml` source `barallat_perez_2025_npjsf`) |

## 1. Why it matters

The matrix layer cites this paper for lupin binding and for the sign of the heat effect on hexanal binding. What it
actually contains: (i) **pooled** (PPI + SPI + LPI) % bound ranges for heptanal, trans-2-heptenal, cis-4-heptenal, 2-octanone,
octanal and 2-octanol at 1 % w/v isolate; (ii) a residual-fat (defatting) test on 2-heptanone … 2-decanone with per-isolate
differences; (iii) hexanal binding over weeks of storage after 5 / 70 / 90 °C, figure-only, with a qualitative ordering
SPI > LPI > PPI ≈ WPI; (iv) **measured surface hydrophobicity and sulfhydryl contents** of the four commercial isolates,
non-heated and heated — the only numeric protein-property data in the paper.

## 2. Methods as they matter to a model

- **Flavours** (Sigma-Aldrich, ≥ 95 %): hexanal, heptanal, trans-2-heptenal, cis-4-heptenal, octanal, 2-octanol, 2-heptanone,
  2-octanone, 2-nonanone, 2-decanone. **Ten compounds; nonanal is not among them.** Stock 10 mg/L in sodium phosphate
  (pH 7.0, 50 mM), 1 h ultrasound at 20 °C.
- **Proteins:** SPI Supro XT219D IP (Solae); PPI FYPP-85-C-EU (AGT); **LPI 10 600 (ProLupin, Lupinus angustifolius)**; WPI BiPro
  (Davisco). N × 6.25. Stored 10–15 °C, dark, dry.
- **Samples:** 2 % w/v protein stock in the same buffer; *"1 mL was added into a 20 mL GC-MS vial followed by the addition of
  1 mL of flavor stock solution. Thus, a final protein solution of 1 wv% and 5 mg/L flavor concentration was obtained."*
  **Protein basis: 1 % w/v of isolate powder (10 g/L isolate), not of protein.** Reference: buffered protein, no flavour.
- **Hold:** *"kept in a water bath ... at 30, 70, or 90 °C, respectively, shaking at 125 rpm for 3 h before headspace analysis.
  Samples were prepared in triplicate. After the preparation of these samples, they were stored at 5 °C and measured with a
  weekly frequency."* Which hold temperature applies to Figs. 1–3 is **not stated** (30 °C is the natural reading).
- **Headspace:** static HS GC-MS (Agilent 7890A/5975C), **incubated and shaken 14 min at 40 °C**, 1 mL headspace, split 1:10,
  DB-WAX 20 m. Flavours analysed one at a time.
- **Binding (eq 1):** binding(%) = (1 − (HS1 − HS2)/HS3) × 100, HS1 = protein + flavour, HS2 = protein without flavour, HS3 =
  flavour without protein. **A % headspace depletion at one protein and one flavour level — not a binding constant.** No
  calibration to concentration; no partition coefficient; no Klotz/Scatchard analysis.
- **Defatting:** chloroform/methanol 1:2 (Bligh & Dyer), 1:9 product:solvent; fat by NMR.
- **H0:** 8-ANS fluorescence (ex 390 / em 470 nm), protein 0.025–0.4 mg/mL in 10 mM phosphate pH 7.0; H0 = slope of FI vs
  concentration, 95 % CI.
- **–SH:** Ellman's reagent (DTNB) in Tris-glycine-EDTA pH 8.0 ± 8 M urea (total vs exposed), 1 h at 95 °C, A412;
  SH (μmol/g) = 73.53 × A412 × D / C.
- **Statistics:** ANOVA + Tukey, p < 0.05.

## 3. Tables re-typed

### Table 1. "Experimental Values of the Surface Hydrophobicity (H0) of Pea (PPI), Soy (SPI), Lupin (LPI), and Whey (WPI) Protein Isolates at a 95% Confidence Interval"

| isolate | H0 5 °C | H0 70 °C | H0 90 °C |
|---|---|---|---|
| PPI | 953 [649, 1257] | 745 [505, 985] | 564 [463, 666] |
| SPI | 2010 [1605, 2415] | 1916 [1541, 2291] | 1385 [1188, 1582] |
| LPI | 834 [752, 916] | 1291 [1056, 1526] | 566 [368, 765] |
| WPI | 1692 [1585, 1797] | 1579 [1245, 1914] | 2374 [2303, 2446] |

### Table 2. "Experimental Values of Sulfhydryl Content (−SH (μmol/g)) of Pea (PPI), Soy (SPI), Lupin (LPI), and Whey (WPI) Protein Isolates Non-Heated (NH) and Heated (H)"

| isolate | exposed –SH (μmol/g) | total –SH (μmol/g) |
|---|---|---|
| PPI_NH | 3.98 ± 0.03 | 9.14 ± 0.24 |
| PPI_H | 2.88 ± 0.04 | 5.34 ± 0.086 |
| SPI_NH | 9.16 ± 0.19 | 15.62 ± 0.25 |
| SPI_H | 5.92 ± 0.15 | 11.54 ± 0.36 |
| LPI_NH | 7.13 ± 0.63 | 11.42 ± 0.64 |
| LPI_H | 4.19 ± 0.18 | 10.26 ± 0.24 |
| WPI_NH | 4.53 ± 0.02 | 8.78 ± 0.04 |
| WPI_H | 2.86 ± 0.08 | 6.39 ± 0.02 |

"Heated" = *"heated (>90 °C)"* (text); the Ellman step itself incubates at 95 °C for 1 h. Per gram of isolate (C = protein
concentration in mg/mL, but the isolates are 83–93 % protein; basis not further specified).

### Numbers printed in the running text

| quantity | value | anchor |
|---|---|---|
| heptanal → trans-2-heptenal, PPI/SPI/LPI | *"significant binding increase from 13.73 to 54.60% ± 3.95 for PPI, SPI, and LPI"*; WPI not significant | Results, Fig. 1 |
| cis-4-heptenal → trans-2-heptenal | *"significantly increased the protein binding from 18.19 to 54.60% ± 2.98"* (all four isolates) | Results, Fig. 1 |
| 2-octanol, all isolates | *"in the range of 28.85−57.7% ± 5.77"* | Results, Fig. 2 |
| 2-octanone, all isolates | *"17.53−23.30% ± 3.36"* | Results, Fig. 2 |
| 2-octanone → octanal, PPI/SPI/LPI | *"a significant binding increase from 14.73 to 52.76% ± 4.65"*; WPI not significant | Results, Fig. 2 |
| fat, manufacturer / defatted | PPI 8.3 wt % → DPPI 2.46; SPI 3.1 → DSPI 1.10; LPI 3 → DLPI 0.11; WPI < 0.05 % (not defatted) | Results, Fig. 3 |
| fat in the 1 % w/v solution | PPI 0.083 wt %; DPPI 0.0246 wt % | Results |
| PPI vs DPPI excess binding | 2-octanone 15.38 % ± 2.92; 2-nonanone 8.9 % ± 3.87; 2-decanone 6.07 % ± 1.53 (not significant) | Results, Fig. 3 |
| SPI vs DSPI (2-heptanone); LPI vs DLPI (2-decanone) | 10.29 % ± 7.95; 8.34 % ± 10.45 higher in the non-defatted isolate | Results, Fig. 3 |
| H0 change at 70 °C | *"decrease by 21.82, 6.62, and 4.70% ± 9.38 for PPI, WPI, and SPI"* | Results, Table 1 |
| H0 change at 90 °C | *"decreased further up to 56, 27.67, and 24.26% ± 17.47 for LPI, SPI, and PPI"* | Results, Table 1 |
| hexanal over storage | SPI, LPI, WPI rise over the first 3 weeks then plateau; SPI > LPI; PPI and WPI lowest; "often not statistically significant"; higher temperature → "slight increase" | Results, Fig. 4 (no numbers) |
| protein content | "PI are mainly composed of protein (>83%, Table S2)" | Results |

### Arithmetic on the printed numbers (all mine) `[D]`

- **H0 at 70 °C** vs 5 °C: PPI −21.8 %, WPI −6.7 %, SPI −4.7 % ✔ (the text's three); **LPI +54.8 %** — an *increase* the text
  does not mention. **H0 at 90 °C** — the text's 56 / 27.67 / 24.26 % are decreases **relative to 70 °C** (LPI 1291 → 566
  = −56.2 %; SPI 1916 → 1385 = −27.7 %; PPI 745 → 564 = −24.3 %), not relative to 5 °C (−32 %, −31 %, −41 %).
- **–SH loss on heating:** exposed −28 % (PPI), −35 % (SPI), −41 % (LPI), −37 % (WPI); total −42 %, −26 %, −10 %, −27 %.
- **2-octanone, internal inconsistency:** the low end for PPI/SPI/LPI is printed as **14.73 %** in one sentence, while another
  sentence gives the 2-octanone range over all four isolates as **17.53−23.30 %** — a pooled range that includes the plant
  isolates cannot start above the plant low end. One of the two is misprinted.
- **Heptanal ↔ octanal:** 13.73 % (heptanal low end) vs 52.76 % (octanal high end) are opposite endpoints of different pooled
  ranges; their ratio is not a chain-length effect.

## 4. Numbers the repository can use

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| pooled % bound, 6 flavours | as §3 | % headspace depletion | 1 % w/v isolate, 5 mg/L, pH 7.0 phosphate 50 mM, HS 40 °C/14 min; hold T not stated (30 °C likely) | text | **level_only, pooled** — endpoints of ranges, not means |
| H0, 4 isolates × 3 T | Table 1 | fluorescence slope, arbitrary | 10 mM phosphate pH 7.0 | Table 1 | measured_property |
| exposed / total –SH, 4 isolates × NH/H | Table 2 | μmol/g isolate | Ellman, ± 8 M urea | Table 2 | **measured_property — a protein thiol density for lupin, pea, soy, whey** (`also_relevant: sulfur`) |
| fat contents | §3 | wt % | NMR / manufacturer | text | measured_level |
| per-isolate lupin binding (2-heptanone…2-decanone; hexanal over storage) | — | % | Figs. 3, 4 | figures | **figure_only** |

**Not in the paper:** any binding constant (K, K_g, Ka, n), any partition coefficient, any nonanal measurement, any mucin arm,
any per-protein printed % bound, any log P model (log P values are in the off-disk Table S1).

## 5. ★ The code's "sealed lupin" values against this paper (checked 2026-10-29)

`grep -rn Barallat src/ data/species/ data/` (py, yml) and the JSON records that cite the same DOI:

| # | file:line | what the repo says | what the paper has | mismatch |
|---|---|---|---|---|
| 1 | `src/kinetic_core/parameters_matrix.py:609-612` (`HOLDOUT_SEALED_BINDING`) | `kg_hexanal_lupin`, `kg_nonanal_lupin`, `kg_2_nonanone_lupin` = "Barallat-Perez 2023 (JAFC 71:20274) lupin" | % bound only (eq 1), lupin values only as bars (2-nonanone in Fig. 3; hexanal in Fig. 4) | **Quantity**: a `kg_` per-gram constant (L/g) does not exist in this paper; any such value would be a derived conversion of a figure-read % at 10 g/L isolate. **Basis**: g isolate, not g protein |
| 2 | `src/kinetic_core/parameters_matrix.py:611` | `kg_nonanal_lupin` → this paper | nonanal is never measured (it appears only in a sentence citing Kühn et al.) | **Compound absent from the paper** |
| 3 | `src/kinetic_core/parameters_matrix.py:613-615` | `kg_*_mucin` → "Barallat-Perez 2024 pig gastric mucin" | 2023 has no mucin; 2024 prints no binding value (the file's own comment at :529-531 says so) | **Mucin keys cite a paper that prints no value**; the comment at :529-534 moved only the lupin citation |
| 4 | `src/kinetic_core/parameters_matrix.py:529-534` | "the in vitro parent study is Barallat-Perez ... 2023 ... The keys stay sealed ... only the citation changes" | see 1–3 | **The re-attribution is itself wrong** for nonanal and for every mucin key. The compound trio hexanal / nonanal / 2-nonanone × lupin / mucin is the design of Barallat-Pérez **2024** (Fig. 4, figure-only) and is contained in the **2025 npj Sci Food** Klotz table (LPI columns, mucin conditions; `binding_constants.yml:247-262, 938-972`), which reports K in M⁻¹ per site, not L/g |
| 5 | `src/kinetic_core/parameters_matrix.py:441-442` | soy-denatured hexanal note: "35x BELOW the headspace-derived lupin value for the same compound" | no lupin hexanal number is printed anywhere in this paper | **A lupin hexanal value (≈ 5 × 10⁻² L/g implied) is referenced but has no printed source here** |
| 6 | `data/lit/binding_constants.yml:620-647` (`barallat2023_2octanone_plant_isolates_pooled`, 14.73 %) | quote correct | the paper also prints 2-octanone 17.53−23.30 % over all four isolates | **Internal inconsistency in the paper, not flagged in the record** |
| 7 | `data/lit/binding_constants.yml:600` (`sample_hold_celsius: 30.0`, octanal record) | 30 °C | methods list 30, 70 or 90 °C without assigning one to Figs. 1–3 | **Inferred, not stated** |
| 8 | `data/lit/process_state_calibrations.json:139, 144, 147-150` (`acs_jafc_3c05991_ppi_spi_partitioning`) | `logp_model_available: true`, `ka_values_reported: true`, target observables hexanal, 2-pentylfuran, **MFT, pyrazines** | no Ka, no log P model, no 2-pentylfuran / MFT / pyrazine | **Three false claims about the paper's content** |
| 9 | `data/lit/benchmark_intake_registry.json:1152, 1171-1172, 1182` (same id) | `chemistry_family: lipid_maillard_crosstalk`, `benchmark_eligible_parameter`, `ready_for_calibration_encoding`, `ka_values_reported: true` | a binding-% paper with no Ka and no Maillard chemistry | **Mis-described; not calibration-ready** |

What **matches**: `binding_constants.yml` records `barallat2023_octanal_plant_isolates_pooled` (52.76 ± 4.65), `..._heptanal_`
(13.73), `..._trans2heptenal_` (54.60 ± 3.95), `..._cis4heptenal_` (18.19), `barallat2023_hexanal_ordering` (qualitative) and
`barallat2023_temperature_effect` (qualitative, ~32 % second-hand canola) — values, units (% bound), 10 g/L isolate basis,
5 mg/L, pH 7.0, 50 mM phosphate, HS 40 °C and method all agree with the paper. `data/species/` has no Barallat entry and no
lupin matrix.

**Net:** no sealed lupin VALUE can be checked because none is carried; but the citation the code gives for them cannot
supply them (wrong quantity for all three, wrong compound for nonanal, no mucin at all). The most likely numeric parent is the
2025 npj paper (Klotz K, M⁻¹ — not convertible to L/g without a protein molar mass, as `binding_constants.yml:41-46` already
rules), or a figure-read of 2024 Fig. 4.

## 5b. Flags

1. Binding is a single-point % depletion (eq 1), not a constant; per-protein values are figure-only.
2. Pooled PPI/SPI/LPI ranges; endpoints, not means.
3. 2-octanone 14.73 % vs 17.53−23.30 % — internal inconsistency.
4. Hold temperature for Figs. 1–3 not stated; Fig. 4 samples are stored at 5 °C and measured weekly after a 3 h hold.
5. H0 90 °C percentages are relative to 70 °C; LPI H0 rises 55 % at 70 °C (unremarked).
6. Fig. 3 caption typo "DDSPI".
7. SI (Tables S1–S4, Fig. S1) not on disk: flavour log P / volatility, manufacturer composition, amino-acid content.

## Full text (read 2026-10-29)

New dossier, written from the full paper `data/articles/barallatperez2023.pdf` on 2026-10-29. Sections 0–5b above are the reading;
the figures were not digitised.
