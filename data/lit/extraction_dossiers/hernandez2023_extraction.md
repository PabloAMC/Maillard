# Hernandez, Woerner, Brooks & Legako 2023 — EXTRACTION (cooked commercial plant-based burgers against ground beef; the source of this repository's PBMA identity row)

**Source on disk:** `data/articles/hernandez2023.pdf` (1.2 MB; downloaded 2026-09-11). Read
2026-09-11 via `pdftotext -layout` with Tables 3 and 4 verified against page renders. Wave B37.

| field | value |
|---|---|
| Title | "Descriptive Sensory Attributes and Volatile Flavor Compounds of Plant-Based Meat Alternatives and Ground Beef" |
| Venue | Molecules 2023, 28, 3151 |
| DOI | 10.3390/molecules28073151 |
| Group | Texas Tech University, Department of Animal and Food Sciences |
| Bundle | `resconi_2023_pbma_beef_identity_benchmark` — whose vessel note said this paper was **NOT ON DISK**; it is, as of 2026-09-11, and B37 corrects the note |
| Products | Beyond Burger, Impossible Burger, a third retail plant-based brand, plus lean (≥93 %) and regular (80–85 %) ground beef, collected from six US cities |
| Cook | verbatim: "cooked on an enamel-lined cast-iron skillet heated to a surface temperature of **200 ± 10 °C** … cooked to an **internal temperature of 71 °C** and were flipped at 35 °C"; 150 g patty; **no cook time is stated** |
| Quantification | SPME (85 µm Carboxen/PDMS), 5 g cooked homogenate, 65 °C incubation, 25 min extraction; 1,2-dichlorobenzene internal standard; **five-level calibration curve against authentic standards**; **ng per gram of sample** |

## 1. The one number the repository already uses

Table 3, furfural (ng/g), least-squares means:

| lean ground beef | regular ground beef | Beyond Meat | Impossible Burger | third retail brand | SEM | p |
|---:|---:|---:|---:|---:|---:|---:|
| 18.67 ᵇ | 20.99 ᵇ | **987.41 ᵃ** | 64.71 ᵇ | **1093.54 ᵃ** | 138.580 | < 0.001 |

The bundle's three plant-based furfural values are confirmed against the print.

## 2. Table 3 — the other Maillard rows the trunk carries (ng/g)

| compound | lean GB | regular GB | Beyond | Impossible | third brand | SEM | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2-Furan methanol | 23.19 | 16.25 | 175.73 | 39.18 | 136.85 | 19.454 | <0.001 |
| 2-Methyl-3-furanthiol | 21.33 | 1.58 | 11.11 | 24.04 | 41.09 | 9.246 | <0.001 |
| 5-Methylfurfural | 6.24 | 5.76 | 24.33 | 6.99 | 21.74 | 2.194 | <0.001 |
| 2,3-Butanedione | 30.15 | 42.01 | 12.79 | 11.84 | 13.95 | 6.348 | 0.002 |
| Methylpyrazine | 12.73 | 10.26 | 30.88 | 38.94 | 24.18 | 5.303 | <0.001 |
| **2,5-Dimethylpyrazine** | 14.46 | 16.97 | 32.37 | 30.92 | 32.09 | 5.203 | 0.017 |
| Trimethylpyrazine | 2.98 | 3.91 | 7.82 | 9.67 | 15.56 | 1.707 | <0.001 |
| 2-Ethyl-3,5/6-dimethylpyrazine | 6.73 | 7.35 | 10.33 | 12.48 | 36.20 | 6.338 | 0.001 |
| **Benzaldehyde** | 36.67 | 23.65 | 62.07 | 164.42 | 89.35 | 17.477 | <0.001 |
| **Methional** | 2.28 | 4.66 | 8.38 | 4.89 | 3.77 | 1.088 | 0.001 |
| 2-Methylbutanal | 18.10 | 23.51 | 12.74 | 20.74 | 15.56 | 1.915 | <0.001 |
| 3-Methylbutanal | 15.69 | 16.91 | 16.61 | 16.22 | 16.26 | 1.071 | 0.918 |
| 2-Acetylpyrrole | 169.69 | 141.75 | 559.78 | 306.20 | 725.88 | 132.360 | <0.001 |

## 3. Table 4 — the lipid rows (ng/g)

| compound | lean GB | regular GB | Beyond | Impossible | third brand | SEM | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Hexanal** | 39.98 | 75.63 | 197.61 | 103.87 | 136.41 | 23.278 | <0.001 |
| **Nonanal** | 30.23 | 48.12 | 69.56 | 35.45 | 56.09 | 8.879 | 0.011 |
| Heptanal | 9.30 | 20.42 | 26.36 | 7.53 | 13.50 | 2.591 | <0.001 |
| Octanal | 11.81 | 22.63 | 26.85 | 14.89 | 17.31 | 3.145 | 0.005 |
| Pentanal | 3.85 | 6.89 | 12.90 | 5.53 | 15.49 | 1.775 | <0.001 |
| Decanal | 9.89 | 10.82 | 19.03 | 20.54 | 13.11 | 3.850 | 0.132 |
| (E,E)-2,4-Decadienal | 15.57 | 17.22 | 174.19 | 17.73 | 95.07 | 34.280 | 0.002 |
| **2-Pentylfuran** | 3.24 | 4.04 | 24.17 | 6.85 | 25.46 | 3.195 | <0.001 |

## 4. Two defects in the published paper, recorded rather than corrected

1. **The dispersion column is not a standard deviation.** It is "largest standard error of the
   least squares means" — one pooled number per compound, not a per-column uncertainty. A bundle
   scoring these values may not read the SEM as that column's SD.
2. **The lean/regular footnotes are transposed and mutually inconsistent.** Table 3 footnotes read
   "¹ Includes 15% and 20% fat ground beef. ² Includes less than or equal to 7% fat ground beef",
   which reverses the Methods (lean = ≥93 % lean = ≤7 % fat); Table 4 prints "² Includes greater than
   or equal to 7% fat ground beef". The Methods are authoritative. The beef columns are not used by
   this repository, so nothing downstream turns on it, but it is recorded.

## 5. What the repository takes

The provenance correction, and nothing else. The bundle's own note already says what this pot is: a
**commercial plant-based product cooked as a reference**, whose formulation and thermal history
before the skillet are unknown (`process_metadata.extrusion_history = commercial_pbma_unknown`), and
whose "150 °C / 60 min" in the bundle is a proxy for an unspecified commercial process rather than
this paper's skillet cook. The print does not change that: it gives a skillet surface temperature,
an internal endpoint and no time. **The paper being on disk makes the citation first-hand; it does
not make the pot chargeable.**
