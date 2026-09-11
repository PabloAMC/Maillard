# Cai, Zhu, Ma, Thakur, Zhang & Wei 2021 — EXTRACTION (soybeans roasted at four temperatures against an unroasted control, with pyrazines and furfural quantified)

**Source on disk:** `data/articles/cai2021.pdf` (3.2 MB; downloaded 2026-09-11). Read 2026-09-11 via
`pdftotext -layout`, Table 4 checked against rendered page images. Wave B37.

| field | value |
|---|---|
| Title | "Effects of roasting level on physicochemical, sensory, and volatile profiles of soybeans using electronic nose and HS-SPME-GC–MS" |
| Venue | Food Chemistry 340 (2021), article 127880 |
| DOI | 10.1016/j.foodchem.2020.127880 |
| System | whole soybeans, cultivar Suinong 25, **lipoxygenase-positive**, roasted in a lab oven |
| Treatments | **140, 170, 200, 230 °C**; the design says 10, 20 and 30 min, and **the volatile table's caption does not state which time it is** (see §4) |
| Control | **NR, unroasted beans** — verbatim, "Unroasted beans were used as control (NR)" |
| Quantification | HS-SPME (50/30 µm DVB/CAR/PDMS), 2.0 g, 60 °C, 40 min; **external standard 1,2-dichlorobenzene run in a devolatilised soybean-flour matrix**, reported as µg/kg sample, n = 3 |

## 1. Why this paper matters

Of the sixteen papers read in B37, this is **the only one that prints a Maillard product series with an
unheated control, as concentrations**. Eleven pyrazines, furfural, 5-methylfurfural, 2-furanmethanol
and maltol are quantified from an unroasted zero to a 230 °C endpoint. The repository's trunk makes
pyrazines and furfural and has almost no external anchor for either in a real plant matrix.

## 2. Table 4 — the rows this model carries (µg/kg sample, mean ± SD, n = 3; `nd` = not detected)

| compound | NR (unroasted) | 140 °C | 170 °C | 200 °C | 230 °C |
|---|---:|---:|---:|---:|---:|
| **Hexanal** | 2.33 ± 0.14 | 0.17 ± 0.02 | nd | nd | nd |
| **1-Hexanol** | 13.42 ± 0.95 | 5.19 ± 0.67 | nd | nd | nd |
| **1-Octen-3-ol** | 2.61 ± 0.22 | 1.11 ± 0.34 | nd | nd | nd |
| **Benzaldehyde** | 1.50 ± 0.08 | 3.44 ± 0.13 | 4.70 ± 0.21 | 2.86 ± 0.11 | nd |
| **2-Pentylfuran** | nd | nd | 2.14 ± 0.12 | 6.53 ± 0.74 | 9.06 ± 0.60 |
| **Furfural** | nd | nd | 1.82 ± 0.56 | 9.43 ± 0.81 | 17.76 ± 1.34 |
| 5-Methylfurfural | nd | nd | nd | nd | 7.20 ± 0.67 |
| 2-Furanmethanol | nd | nd | 1.32 ± 0.09 | 3.56 ± 0.06 | 2.22 ± 0.15 |
| Maltol | nd | 28.63 ± 2.28 | 79.13 ± 2.76 | 91.01 ± 3.15 | 98.67 ± 2.97 |
| **2,5-Dimethylpyrazine** | nd | 16.58 ± 1.34 | 82.43 ± 7.61 | 278.91 ± 6.83 | 411.18 ± 9.55 |
| 2-Methylpyrazine | nd | 0.30 ± 0.05 | 4.38 ± 0.19 | 18.18 ± 0.47 | 54.04 ± 2.82 |
| 2,3,5-Trimethylpyrazine | nd | 0.66 ± 0.27 | 3.62 ± 0.44 | 17.51 ± 0.95 | 72.68 ± 5.33 |
| 3-Ethyl-2,5-dimethylpyrazine | nd | 27.55 ± 1.23 | 63.50 ± 4.06 | 73.28 ± 3.85 | 123.59 ± 9.46 |
| 2-Ethyl-5-methylpyrazine | nd | 35.38 ± 4.04 | 77.09 ± 8.65 | 192.96 ± 4.71 | 43.67 ± 3.32 |
| 2,3-Diethyl-5-methylpyrazine | nd | 7.59 ± 0.83 | 12.57 ± 0.69 | 12.31 ± 1.77 | 21.67 ± 1.23 |
| 2,5-Diethylpyrazine | nd | 2.83 ± 0.17 | 3.53 ± 0.14 | 7.11 ± 0.32 | 18.12 ± 0.95 |
| 2-Methyl-5-vinylpyrazine | nd | 2.39 ± 0.31 | 5.10 ± 0.27 | 6.75 ± 0.74 | 8.96 ± 1.55 |
| 2-Acetylpyrazine | nd | nd | 4.25 ± 1.19 | 8.43 ± 0.90 | 26.19 ± 3.10 |
| 2-Acetyl-3-methylpyrazine | nd | nd | nd | 6.24 ± 1.03 | 26.79 ± 2.20 |
| 2-Acetyl-3-ethylpyrazine | nd | nd | 2.29 ± 0.41 | 18.57 ± 2.49 | 17.73 ± 1.36 |
| 2-Acetylpyrrole | nd | nd | nd | 4.01 ± 0.36 | 12.62 ± 0.78 |
| **Total** | **26.91** | **160.09** | **421.98** | **893.10** | **1140.51** |

**Nonanal, pentanal, heptanal, octanal and decanal are not measured in this paper at all.**

## 3. The shape of the result, which is the model's own story

The lipid-derived volatiles (hexanal, 1-hexanol, 1-octen-3-ol) **fall to nothing by 170 °C** while the
Maillard products (pyrazines, furfural, 2-pentylfuran) **rise from an unroasted zero**. The unroasted
control has 26.91 µg/kg of total volatiles; the 230 °C bean has 1140.51. The crossover the model
would have to reproduce — lipid down, Maillard up, both over one temperature axis — is printed here
with a genuine zero at the cold end.

## 4. Four limits, all of which must travel with any use of this table

1. **The roasting TIME for Table 4 is not printed.** The caption gives no time; §2.2 describes 10,
   20 and 30 min; the electronic-nose section says the 20-min group was used "as representatives",
   and the discussion speaks of 20 min. **20 minutes is the probable but unstated reading.** A
   benchmark built on this table must record the time as unstated, or not be built.
2. **The quantification is single-external-standard.** All 41 analytes are 1,2-dichlorobenzene
   equivalents; there is no compound-specific response factor. It is one tier below a calibrated
   concentration and must carry that class.
3. **`nd` is not zero.** No detection limit is printed. Hexanal at 170 °C and above, and every
   Maillard product in the NR column, are "not detected", which is a bound, not a measurement.
4. **A whole soybean is not a pot.** This is a dry bean roasted in air, not a charged aqueous
   system: precursor concentrations are not stated in any usable form, and the model would have to
   declare them.

## 5. Verdict

The best external Maillard anchor of the sixteen, and the only one with an unheated zero. It is a
candidate for a new external-validation family — a **dry roast temperature ladder** — and that is a
pre-registered decision for a later wave, not something wave B37 installs. Fig. 1 (moisture, protein,
reducing sugar, fat, lipoxygenase, peroxidase against temperature and time) is figure-only.
