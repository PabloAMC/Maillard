# Liu, Deng, Yang, Zhang, Guo, Li, Chen, He, Chen & Zeng 2025 — EXTRACTION (the lipophilic-protein fraction of soy heated 120–180 °C: a clean temperature series whose per-compound data did not arrive)

**Source on disk:** `data/articles/Liu2025.pdf` (downloaded 2026-09-11). Read 2026-09-11 via
`pdftotext -layout`, methods re-checked against a page render. Wave B37.

| field | value |
|---|---|
| Title | "Effect of Lipophilic Proteins on the Bean Flavor Formation: Lipids Oxidation and Protein Modification" |
| Venue | Food Biophysics (2025) 20:190 |
| DOI | 10.1007/s11483-025-10070-z |
| Material | **lipophilic protein (LP)**, a fraction of soy protein isolate — verbatim, "Lipophilic proteins (LPs) constitute approximately 31% of soy protein isolate and contain 11–13% lipids". **Not SPI.** |
| Treatment | verbatim: "2 g of the sample … deionized water … at a powder-to-water weight ratio of 1:4 … placed in a forced convection drying oven that has been heated to 120 °C for 10 min. This procedure is to be repeated to prepare samples at 140 °C, 160 °C, and 180 °C, **with unheated samples serving as the control**." |

## 1. The design is exactly right and the data are not in the paper

A four-point **isochronal temperature series at 10 minutes and 20 % solids — 120, 140, 160, 180 °C —
with an unheated control**: the cleanest moist-heat temperature ladder among the sixteen papers.
**But there is no table of individual volatile compounds anywhere in the main text.** Every numbered
display item is a figure. The per-compound, per-temperature values are cited to **Supplementary
Table 1, which is not on disk** and for which the paper states "No datasets were generated or
analysed during the current study".

## 2. What IS printed

**Group totals of the 22 bean-flavour volatiles, in µg/L, without standard deviations:**

| 120 °C | 140 °C | 160 °C | 180 °C | LP + curcumin 180 °C | LP + tea polyphenols 180 °C |
|---:|---:|---:|---:|---:|---:|
| 22.72 | 34.94 | 35.95 | 67.70 | 48.54 | 29.04 |

**The unheated control's total is not printed.** The abstract confirms 22.72 µg/L is the 120 °C point.

**Malondialdehyde, by an external-standard method (µmol/g):** unheated **0.19**, then 0.11 (120 °C),
0.13 (140 °C), 0.17 (160 °C), 0.18 (180 °C). Note this runs **below** the unheated control at
120–160 °C; the authors attribute it to malondialdehyde being consumed into Schiff bases.

**Schiff base:** 7.77 unheated, 14.91 at 120 °C, 10.99 at 180 °C (140 and 160 °C figure-only; units
not given). **Protein carbonyl at 180 °C:** 1.83 mmol/g without additive, 1.66 with curcumin, 1.53
with tea polyphenols.

**Compound counts** (Fig. 2, an UpSet plot): 59 unheated, 62 at 120 °C, 79 at 140 °C, 81 at 160 °C,
102 at 180 °C.

**Fatty-acid losses at 180 °C, as percentages:** C18:1 −8.40 %, C18:2 −20.17 %, C18:3n3 −31.90 %.

## 3. And the quantification basis is never stated

There is **no internal standard, no calibration curve, no response factor and no library statement**
in the volatiles methodology — only the HS-SPME and GC-MS instrument conditions — yet totals are
reported in µg/L, a volume basis that is itself undefined for 2.0 g of a 20 %-solids paste weighed
into a vial. By contrast the paper *does* state a basis where it has one ("The malondialdehyde
content was calculated by the external standard method"). The absence is therefore real, not a
text-extraction artefact.

## 4. Verdict

**Data unavailable.** The material is a sub-fraction rather than the isolate, the per-compound table
is in a missing supplement, the control's total is not printed, and the quantification basis is not
declared. Recorded in full so that the design is on file: if Supplementary Table 1 is ever obtained,
this becomes a four-temperature moist-heat ladder with a blank, which is close to what
`docs/guides/EXPERIMENTS.md` asks for. Added to the reading list as a supplement request.
