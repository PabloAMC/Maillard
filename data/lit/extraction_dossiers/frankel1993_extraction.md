# Frankel 1993 — EXTRACTION (apparent activation energies for hexanal formation by thermal decomposition of oxidized oils)

**Source on disk:** `data/articles/frankel1993.pdf` (0.6 MB, scan with a text layer; downloaded
2026-09-11 at this repository's request). Read 2026-09-11 via `pdftotext -layout`. Wave B37.

| field | value |
|---|---|
| Title | "Formation of Headspace Volatiles by Thermal Decomposition of Oxidized Fish Oils vs. Oxidized Vegetable Oils" |
| Venue | Journal of the American Oil Chemists' Society 70(8):767–772 (August 1993) |
| DOI | 10.1007/BF02542598 — **NOT PRINTED IN THE PAPER.** This is a 1993 scan that carries no DOI; the identifier comes from the Springer landing page for this article, supplied with the 2026-09-11 reading list. Recorded as second-hand rather than read |
| Author | E. N. Frankel, Department of Food Science and Technology, University of California, Davis |
| Systems | menhaden, sardine, soybean, safflower, canola and high-oleic sunflower oils, **pre-autoxidized at 40 or 50 °C to a stated peroxide value**, then thermally decomposed |
| What is measured | headspace volatiles (propanal, pentane, **hexanal**) by static equilibrium headspace GC (Tekmar 7000 autosampler, 9-mL headspace vials), at decomposition temperatures from 40 to 180 °C in 10 °C intervals; Arrhenius plots per 30 °C window |

## 1. Why this paper matters to this model

`src/kinetic_core/parameters_lipid.py` carries the lipid lane's rate as a single anchor measured at
**25 °C** (Schroen & Berton-Carabin 2022, k4 = 6e-3 /h, hand-fitted by visual agreement) with the
flag `temperature_dependence_UNMEASURED`, and bridges the gap to cooking temperatures with
`Q10_ASSUMPTION`, a constant 2–3 taken from the same authors' prose. Every lipid row in the
scorecard therefore prints an extrapolation warning. **This paper measures the temperature
dependence of hexanal formation directly**, in the decomposition of an existing hydroperoxide pool,
which is the same step the model calls `k_LOOH_decomp` → hexanal.

## 2. What was actually done (and the caveat that governs everything below)

The oils were **first autoxidized at 40–50 °C to a stated peroxide value, then heated** and the
headspace analysed. The measured activation energy is therefore the activation energy of
**decomposition of an already-formed peroxide pool**, not of the whole autoxidation chain. That is
the right quantity for this model. But it is measured in **BULK OIL**, not in an aqueous emulsion
and not in a protein matrix, and the repository's standing ask was for a protein-containing food
matrix. The transfer is from oil to matrix and must be declared as such.

## 3. Table 2 — apparent activation energy of thermal decomposition, TOTAL volatiles (kcal/mol)

Verbatim, all 29 runs as printed ("Apparent Activation Energy of Thermal Decomposition of Fish and
Vegetable Oils Autoxidized at 50 °C").

| run | oil | PV (meq/kg) | temperature (°C) | Ea (kcal/mol) | group average |
|---:|---|---:|---|---:|---|
| 1 | Menhaden | 3.1 | 40–70 | 6.7 | |
| 2 | Menhaden | 2.8 | 40–70 | 6.3 | |
| 3 | Menhaden | 4.6 | 40–70 | 6.5 | |
| 4 | Sardine | 57.0 | 40–70 | 6.7 | **6.6 ± 0.2** |
| 5 | Menhaden | 5.8 | 70–100 | 11.5 | |
| 6 | Menhaden | 5.4 | 70–100 | 11.7 | |
| 7 | Sardine | 7.3 | 70–100 | 11.9 | **11.7 ± 0.2** |
| 8 | Menhaden | 5.9 | 100–130 | 18.2 | |
| 9 | Sardine | 3.6 | 100–130 | 18.6 | **18.4 ± 0.2** |
| 10 | Menhaden | 2.0 | 130–160 | 18.4 | |
| 11 | Sardine | 19.2 | 130–160 | 18.9 | **18.7 ± 0.3** |
| 12 | Soybean | 4.9 | 100–130 | 21.4 | |
| 13 | Safflower | 19.4 | 100–130 | 22.1 | |
| 14 | Canola | 7.3 | 100–130 | 20.2 | **21.2 ± 0.8** |
| 15 | Soybean | 4.9 | 130–160 | 26.2 | |
| 16 | Soybean | 11.9 | 130–160 | 25.1 | |
| 17 | Safflower | 8.3 | 130–160 | 26.3 | |
| 18 | Safflower | 16.0 | 130–160 | 27.7 | |
| 19 | Canola | 7.3 | 130–160 | 27.9 | |
| 20 | Canola | 10.1 | 130–160 | 26.3 | **26.6 ± 1.0** |
| 21 | Soybean | 3.1 | 150–180 | 45.5 | |
| 22 | Safflower | 3.0 | 150–180 | 46.1 | |
| 23 | Canola | 1.3 | 150–180 | 48.6 | |
| 24 | Canola | 10.5 | 150–180 | 48.7 | **47.2 ± 1.4** |
| 25 | High-oleic sunflower | 6.2 | 130–160 | 51.1 | |
| 26 | High-oleic sunflower | 7.1 | 130–160 | 45.3 | |
| 27 | High-oleic safflower | 11.2 | 130–160 | 53.6 | **50.0 ± 3.5** |
| 28 | High-oleic sunflower | 0.8 | 150–180 | 56.9 | |
| 29 | High-oleic sunflower | 5.1 | 150–180 | 62.4 | **59.7 ± 2.8** |

## 4. Table 3 — apparent activation energy of formation of INDIVIDUAL volatiles (kcal/mol)

This is the table this repository needs: a **hexanal-specific** activation energy.

| run | oil | temperature window (from Table 2) | pentane | propanal | **hexanal** |
|---:|---|---|---:|---:|---:|
| 5 | Menhaden | 70–100 °C | 11.8 | 11.5 | – |
| 6 | Menhaden | 70–100 °C | 12.0 | 11.7 | – |
| 8 | Menhaden | 100–130 °C | 18.5 | 18.2 | **19.2** |
| 9 | Sardine | 100–130 °C | 22.8 | 22.2 | **23.7** |
| 15 | Soybean | 130–160 °C | 26.1 | 25.0 | **27.2** |
| 18 | Safflower | 130–160 °C | 28.8 | 28.1 | **29.2** |
| 19 | Canola | 130–160 °C | 26.0 | 25.3 | **27.3** |
| 24 | Canola | 150–180 °C | 47.4 | 46.8 | **49.0** |

In SI, the rows that matter for a linoleate-carrying plant matrix (the n-6 vegetable oils):

| oil, window | Ea (kcal/mol) | Ea (kJ/mol) |
|---|---:|---:|
| Soybean, 130–160 °C | 27.2 | **113.8** |
| Safflower, 130–160 °C | 29.2 | **122.2** |
| Canola, 130–160 °C | 27.3 | **114.2** |
| Canola, 150–180 °C | 49.0 | **205.0** |

## 5. The finding that constrains how the number may be used

**The apparent activation energy RISES with the decomposition temperature window**, in every oil:
fish oils 6.6 → 11.7 → 18.4 → 18.7 kcal/mol across 40–70, 70–100, 100–130 and 130–160 °C;
polyunsaturated vegetable oils 21.2 → 26.6 → 47.2 across 100–130, 130–160 and 150–180 °C. The
authors say why, verbatim: *"Because activation energies increased with temperature, volatile
formation is apparently more difficult by thermal decomposition of secondary products than from the
corresponding hydroperoxide precursors."* And: *"The changes observed in apparent activation energy
with temperature of decomposition indicate that mechanistic changes occurred."*

**A single Arrhenius line therefore does not describe this system across the model's 60–140 °C
window, and neither does a single Q10.** The paper's own Arrhenius plots are fitted only inside
30 °C windows (r = 0.945–0.999). Any use of these numbers must name its window.

Peroxide value had **no** significant effect on the activation energy (runs 15/16, 17/18, 23/24,
25/27 span 1.8- to 8-fold differences in PV) — which is a genuine convenience: the extrapolation
does not have to know the pot's starting peroxide value to use the barrier.

## 6. Verdict

The first measured, **hexanal-specific** activation energy for hydroperoxide decomposition in the
corpus, at cooking temperatures. It is FIT-class evidence for the lipid lane's missing temperature
term, with three declared limits: bulk oil rather than an aqueous or protein matrix; the barrier is
window-dependent by the authors' own analysis; and the oils are n-6 seed oils, which is the right
fatty-acid class for pea and soy but not the right physical phase. Not acted on in B37.
