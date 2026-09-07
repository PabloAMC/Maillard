# Hofmann, Münch & Schieberle 2000 — EXTRACTION (Strecker aldehyde AND acid yields from glucose and from four α-dicarbonyls; air vs argon; pH; Cu²⁺)

**Source on disk:** `data/articles/hofmann2000.pdf` (owner's download, 2026-08-28). Read-only extraction,
2026-09-07, from `pdftotext -layout`; every table's text layer was clean (numbers re-typed verbatim below).
This is the paper `tasks/audit_remediation.md` §(b) called "**the** Strecker paper — nothing open substitutes
for it", and the "Hofmann 2000" that `data/lit/arrhenius_params.yml` (`strecker.A_value` comment) credits
without a citation. ⚠ It prints **no rate constant, no Arrhenius parameter and no temperature series**
(single temperature, 98 °C / reflux); the `A_value = 1e10` and `Ea = 80.33 kJ/mol` cannot come from here.

**Provenance codes:** **[M]** measured and printed · **[D]** derived by this extraction · **[FIG]** figure only,
no digitised numbers · **[NEG]** verified negative.

## 0. Identity

| field | value |
|---|---|
| Title | "Quantitative Model Studies on the Formation of Aroma-Active Aldehydes and Acids by Strecker-Type Reactions" |
| Authors | Thomas Hofmann, Petra Münch, Peter Schieberle (Deutsche Forschungsanstalt für Lebensmittelchemie, Garching) |
| Venue | J. Agric. Food Chem. 2000, 48 (2), 434–440 |
| DOI | 10.1021/jf990954c (received 24 Aug 1999, accepted 4 Nov 1999, web 4 Jan 2000) |
| Cited by the companion | Hofmann & Schieberle 2000 (JAFC 48:4301, `hofmann2000b_extraction.md`) cites this as "Hofmann et al., 1999/2000" |

## 1. Why it matters

The corpus's only SIDA-quantified Strecker aldehyde **and acid** yields from a free sugar + amino acid and
from each of the four α-dicarbonyls the trunk lane carries (2-oxopropanal = methylglyoxal, glyoxal,
3-deoxy-2-hexosulose = 3-DG, 2-hexosulose = glucosone), each measured **with and without oxygen** in the same
pot. It (i) ranks the dicarbonyls as Strecker donors, (ii) shows the aldehyde is oxygen-INdependent while the
acid is oxygen-dependent (×4–5.5 under air) and Cu²⁺-catalysed, (iii) shows by ¹³C₂ labelling that the acid
is NOT formed by oxidising the aldehyde (a parallel branch from the hemiaminal/enaminol, Fig. 6), and
(iv) gives a 500-min time course of PA and PAA from MGO + Phe with the fraction of Phe consumed. The model's
sugar path carries Strecker degradation with Martins 2005 constants and an "oxidant" pool nothing pins;
this paper is the first in the corpus to measure an oxygen effect on a Strecker product.

## 2. Methods as they matter to a model

| item | value | where |
|---|---|---|
| Glucose/amino acid pot (Tables 1, 2, Fig. 2, Fig. 3, Table 5) | glucose 1.0 mmol + L-Phe (or L-Ala, L-Leu) 1.0 mmol in **10 mL phosphate buffer pH 7.0**, i.e. 0.1 mol/L each | p. 435, Table 1–2 footnotes |
| Buffer strength | ⚠ **inconsistent in the paper**: Experimental says "0.5 mmol/L" (p. 435, twice); every table footnote says **0.5 mol/L** (Tables 1, 2, 3, 5) or **0.1 mol/L** (Table 4). 0.5 mol/L is the Schieberle-group standard (cf. Hofmann 1998) and the likely truth; "mmol" is a typo | p. 435 vs footnotes |
| Temperature | "heated at 98 °C" (Experimental) = "refluxed" (footnotes); aqueous reflux, so ≈ 98–100 °C | p. 435 |
| Time | 30 min (Tables 1, 2, 5); 10–500 min (Table 3, Fig. 2, Fig. 3); 60 min (Table 4) | footnotes |
| Vessel | "closed vials" for the dicarbonyl pots and the MGO time course; reflux for the glucose pots | p. 435, Table 3 fn |
| Atmosphere | **Tables 1, 2, 3, 5, Figs 2–5: NOT controlled** ("oxygen, the concentration of which was not controlled in the experiments displayed in Table 3", p. 437) — i.e. air in a closed vial. **Table 4 only:** model A = argon bubbled ≥ 20 min before heating, closed vial ("oxygen was absent"); model B = "air oxygen" ("oxygen was present"); expt 2 = air + CuSO₄ 0.05 mmol (5 mmol/L) | p. 435, 437 |
| Dicarbonyl pots (Table 4, Fig. 4–5) | α-dicarbonyl 1.0 mmol + L-Phe 1.0 mmol, 10 mL phosphate pH 7.0 (0.1 mol/L per Table 4 fn), 98 °C, 60 min, closed vial | Table 4 fn |
| pH series (Figs 4, 5) | glucose or each dicarbonyl + Phe at pH 3.0 / 5.0 / 7.0 / 9.0 (figure only; buffer for pH 3 and 9 not stated) | p. 438–439 |
| Analytes and method | PA, PAA, 3-methylbutanal, 3-methylbutanoic acid, acetic acid: **stable isotope dilution assay**, GC-MS/CI, [¹³C₂] or [²H₂] internal standards, ether extraction at pH 3.0; acetaldehyde: headspace GC-MS/CI with [¹³C₂]acetaldehyde. Dicarbonyls: 1,2-diaminobenzene trapping, 3 h 30 °C → quinoxalines; 3-DG and glucosone by HPLC-DAD 320 nm, glyoxal and MGO by GC-MS/CI vs [¹³C₄]butane-2,3-dione | p. 435–436 |
| Quantification basis | **µmol product per mmol amino acid charged** (= 0.1 mol %); the text converts 0.42 and 0.29 µmol/mmol to "0.042 and 0.029 mol %" (p. 436). Table 3 prints absolute µmol in the 10 mL / 1 mmol pot, numerically the same scale | p. 436 |
| Replicates | **[NEG]** not stated anywhere; no error bars or SDs printed |

## 3. Every table, verbatim

### Table 1 — "Key Odorants (FD ≥ 16) Generated from L-Phenylalanine by Thermal Treatment in the Presence of Glucose"
Footnote a: "A solution of glucose (1.0 mmol) and L-phenylalanine (1.0 mmol), dissolved in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0), was refluxed for 30 min."

| odorant | odor quality | FD factor |
|---|---|---:|
| phenylacetaldehyde | flowery | 1024 |
| phenylacetic acid | honey-like | 512 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone | caramel-like | 64 |

### Table 2 — "Amounts of Strecker Aldehydes and Acids Generated from Three Different α-Amino Acids" **[M]**
Footnote a: "A solution of glucose (1.0 mmol) and the respective amino acid (1.0 mmol), dissolved in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0), was refluxed for 30 min. b na, not analyzed." Units: **amount (µmol/mmol)** generated from the amino acid.

| odorant | phenylalanine | alanine | leucine |
|---|---:|---:|---:|
| phenylacetaldehyde | 0.42 | na | na |
| phenylacetic acid | 0.29 | na | na |
| acetaldehyde | na | 0.48 | na |
| acetic acid | na | 0.30 | na |
| 3-methylbutanal | na | na | 0.89 |
| 3-methylbutanoic acid | na | na | 0.49 |

**[D]** In mol % of amino acid: PA 0.042, PAA 0.029, acetaldehyde 0.048, acetic acid 0.030, 3-methylbutanal 0.089,
3-methylbutanoic acid 0.049. Aldehyde/acid ratio at 30 min: Phe 1.45, Ala 1.60, Leu 1.82. Leu/Phe aldehyde 2.1×.
⚠ The alanine "acetic acid 0.30 µmol/mmol" is the Strecker acid only in the authors' reading; acetic acid is
also a sugar-degradation product (Martins 2005 step 8, 1-DG → AA), and the SIDA cannot tell the two apart —
do not read 0.30 as a Strecker-only yield.

### Table 3 — "Time Course of the Formation of Phenylacetaldehyde (PA) and Phenylacetic Acid (PAA) from L-Phenylalanine and 2-Oxopropanal" **[M]**
Footnote a: "2-Oxopropanal (1 mmol) and L-phenylalanine (1 mmol), dissolved in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0), were heated (98 °C) in a closed vial." Atmosphere NOT controlled (air).

| reaction time (min) | PA (µmol) | PAA (µmol) | phenylalanine degraded (%) |
|---:|---:|---:|---:|
| 10 | 6.2 | 10.1 | 9 |
| 30 | 9.5 | 17.2 | 18 |
| 60 | 10.6 | 20.6 | 36 |
| 120 | 12.3 | 23.1 | 49 |
| 240 | 13.5 | 24.8 | 58 |
| 300 | 13.3 | 25.8 | 61 |
| 500 | 12.2 | 27.6 | 69 |

**[D]** PA + PAA as a share of the Phe consumed (1000 µmol charged): 10 min 16.3/90 = **18 %**; 30 min 26.7/180 = 15 %;
60 min 31.2/360 = 8.7 %; 120 min 35.4/490 = 7.2 %; 500 min 39.8/690 = **5.8 %**. PA plateaus at ~120–240 min
and falls 10 % by 500 min while Phe keeps being consumed; PAA/PA rises from 1.6 (10 min) to 2.3 (500 min).

### Table 4 — "Amounts of Phenylacetaldehyde (PA) and Phenylacetic Acid (PAA) Generated from L-Phenylalanine in the Presence of Different α-Dicarbonyls: Influence of Oxygen and Copper Ions" **[M]** — THE OXYGEN TABLE
Footnotes: "a A solution of the α-dicarbonyl (1.0 mmol) and L-phenylalanine (1.0 mmol) in phosphate buffer (10 mL; 0.1 mol/L, pH 7.0) was refluxed for 60 min in a closed vial. b Oxygen was absent. c Oxygen was present. d The reaction was performed in the presence of oxygen and copper ions (0.05 mmol). e na, not analyzed." Text (p. 437): model A = "atmosphere of ... argon", model B = "air oxygen".

| expt | dicarbonyl | PA model A (argon) | PA model B (air) | PAA model A (argon) | PAA model B (air) |
|---:|---|---:|---:|---:|---:|
| 1 | 2-oxopropanal | 11.1 | 10.2 | 4.7 | 20.4 |
| 2 | 2-oxopropanal + Cu²⁺ (air) | na | 10.0 | na | 27.1 |
| 3 | glyoxal | 9.2 | 8.3 | 2.2 | 12.2 |
| 4 | 3-deoxy-2-hexosulose | 8.4 | 7.8 | 1.5 | 3.3 |
| 5 | 2-hexosulose | 8.1 | 6.9 | 1.3 | 2.1 |

Units µmol/mmol Phe (÷10 for mol %). **[D]** derived ratios:

| quantity | MGO | glyoxal | 3-DG | glucosone |
|---|---:|---:|---:|---:|
| PA air / argon | 0.92 | 0.90 | 0.93 | 0.85 |
| PAA air / argon | **4.34** ("factor of nearly 4") | **5.5** | 2.2 | 1.6 |
| PA / PAA under argon | 2.4 | 4.2 | 5.6 | 6.2 |
| PA / PAA under air | 0.50 | 0.68 | 2.4 | 3.3 |
| PA (argon) relative to MGO | 1.00 | 0.83 | 0.76 | 0.73 |
| PAA (air) relative to MGO | 1.00 | 0.60 | 0.16 | 0.10 |
| PA + PAA under air, mol % | 3.06 | 2.05 | 1.11 | 0.90 |

Cu²⁺ (5 mmol/L) under air: PAA 27.1 vs 20.4 = **×1.33**; PA unchanged (10.0 vs 10.2).

### Table 5 — "Amounts of ¹³C₂-Labeled Phenylacetic Acid Generated upon Heating of [¹³C₂]Phenylacetaldehyde in the Presence of L-Phenylalanine and Glucose" **[M]**
Footnotes: "a A solution of glucose (1.0 mmol), L-phenylalanine (1.0 mmol), and [¹³C₂]phenylacetaldehyde (50.1 µg) in phosphate buffer (10 mL; 0.5 mol/L, pH 7.0) was refluxed for 30 min. b The acid was quantified by using furan-2-carboxylic acid as the internal standard."

| acid isotopomer | yield (µg/mmol of phenylalanine) |
|---|---:|
| phenylacetic acid | 40.1 |
| [¹³C₂]phenylacetic acid | < 0.5 |

**[D]** 40.1 µg/mmol ÷ 136.15 g/mol = 0.295 µmol/mmol, identical to Table 2's 0.29 (a consistency check that
passes). Of the 50.1 µg labelled aldehyde added (≈ 0.41 µmol, i.e. the same order as the PA the pot itself makes),
< 0.5 µg (< 1 %) became labelled acid ⇒ **aldehyde → acid oxidation is < 1 % per 30 min at reflux under air**.

### Figures (no numbers in the text layer) **[FIG]**
- Fig. 2 (glucose/Phe, reflux, air): PA slightly ahead of PAA for the first 60 min, PAA then dominates; after ~8 h
  PAA ≈ 2 × PA (text p. 436).
- Fig. 3 (dicarbonyls in glucose/Phe, reflux): 3-DG and glyoxal rise fastest and peak at ~60 min; glucosone maximum
  **0.009 mol %** of glucose; MGO keeps rising and after 500 min is **3 ×** glyoxal (text p. 436–437). 1-DG was not
  measured.
- Fig. 4 (PA vs pH 3/5/7/9, each donor + glucose): PA maximal at **pH 5.0** for every donor; MGO gives **2 ×** more
  PA at pH 5 than at pH 9; all dicarbonyls beat glucose at every pH.
- Fig. 5 (PAA vs pH): PAA rises monotonically pH 3 → 9; MGO gives **7 ×** more PAA at pH 9 than at pH 3. Abstract:
  PA:PAA = **3:1 at pH 3.0**, **1:5 at pH 9.0** (MGO donor).

## 4. What the repo could take

Under the owner's rule (within-study ratios are primary evidence and may FIT; end-to-end levels VALIDATE).

**Fed-intermediate yields (levels — VALIDATE rows; each is one pot, 98 °C, pH 7, 60 min unless stated):**

| row | donor + amino acid | atmosphere | product | mol % of amino acid |
|---|---|---|---|---:|
| H00-L1 | MGO 0.1 M + Phe 0.1 M | argon | PA | 1.11 |
| H00-L2 | same | air | PA | 1.02 |
| H00-L3 | same | argon | PAA | 0.47 |
| H00-L4 | same | air | PAA | 2.04 |
| H00-L5 | glyoxal + Phe | argon / air | PA | 0.92 / 0.83 |
| H00-L6 | glyoxal + Phe | argon / air | PAA | 0.22 / 1.22 |
| H00-L7 | 3-DG + Phe | argon / air | PA | 0.84 / 0.78 |
| H00-L8 | 3-DG + Phe | argon / air | PAA | 0.15 / 0.33 |
| H00-L9 | glucosone + Phe | argon / air | PA | 0.81 / 0.69 |
| H00-L10 | glucosone + Phe | argon / air | PAA | 0.13 / 0.21 |
| H00-L11 | glucose + Phe, 30 min reflux, air | — | PA / PAA | 0.042 / 0.029 |
| H00-L12 | glucose + Leu, 30 min | — | 3-methylbutanal / -oic acid | 0.089 / 0.049 |
| H00-L13 | glucose + Ala, 30 min | — | acetaldehyde / acetic acid | 0.048 / 0.030 (acetic acid not Strecker-only) |
| H00-L14 | MGO + Phe time course (Table 3, air, 0.5 M buffer) | — | PA, PAA at 7 times | see Table 3 |

**Within-study ratios (FIT candidates):**

| row | ratio | value | note |
|---|---|---:|---|
| H00-R1 | PAA(air)/PAA(argon), MGO donor | 4.34 | the oxygen effect on the Strecker ACID |
| H00-R2 | PAA(air)/PAA(argon), glyoxal | 5.5 | |
| H00-R3 | PAA(air)/PAA(argon), 3-DG / glucosone | 2.2 / 1.6 | C6 donors barely oxygen-sensitive |
| H00-R4 | PA(air)/PA(argon), all four donors | 0.85–0.93 | the ALDEHYDE from a dicarbonyl is oxygen-independent (slightly lower under air) |
| H00-R5 | PA per donor under argon, MGO : GO : 3-DG : glucosone | 1 : 0.83 : 0.76 : 0.73 | donor identity moves the aldehyde < 1.4× |
| H00-R6 | PAA per donor under air, MGO : GO : 3-DG : glucosone | 1 : 0.60 : 0.16 : 0.10 | donor identity moves the acid 10× |
| H00-R7 | PAA(air + 5 mM Cu²⁺)/PAA(air), MGO | 1.33 | |
| H00-R8 | PA from MGO donor (30 min, Table 3) / PA from glucose (30 min, Table 2) | 9.5 / 0.42 = 23 | same buffer, same T, same time |
| H00-R9 | PA at pH 5 / PA at pH 9, MGO donor | 2 | Fig. 4, text |
| H00-R10 | PAA at pH 9 / PAA at pH 3, MGO donor | 7 | Fig. 5, text |
| H00-R11 | labelled acid / labelled aldehyde added, 30 min | < 1 % | Table 5 |
| H00-R12 | (PA + PAA) / Phe consumed, MGO donor | 18 % at 10 min → 5.8 % at 500 min | Table 3 |
| H00-R13 | MGO / glyoxal at 500 min in glucose/Phe | 3 | Fig. 3, text |

**Directional claims a model must reproduce:**
1. The Strecker ACID is a parallel product of the hemiaminal/enaminol intermediate (Fig. 6), oxygen- and
   metal-dependent, **not** a downstream oxidation of the aldehyde (Table 5). A model that makes PAA by oxidising PA
   is refuted (R11).
2. The Strecker ALDEHYDE from an α-dicarbonyl does not need oxygen (R4). In the model's language: the
   dicarbonyl-initiated Strecker step should not draw on the oxidant pool; only the acid branch should.
3. C6 dicarbonyls (3-DG, glucosone) give aldehyde but almost no acid; C2/C3 fragments (glyoxal, MGO) give both
   (R6, Fig. 7 argument: cyclic hemiacetal forms lack the enaminol to oxidise).
4. Aldehyde favoured at pH 5, acid favoured at pH 9 (R9, R10) — a pH switch on the branch ratio of 15× between
   pH 3 and 9.
5. In glucose/Phe under air the aldehyde leads for the first hour (3-DG era, Fig. 3) and the acid dominates after
   (MGO era) — the donor time course explains the product time course.
6. Strecker products account for at most ~18 % of the amino acid consumed by MGO (R12); most Phe goes elsewhere
   (condensation, imines). The model's "amine returns" bookkeeping should not assume the amino acid is conserved
   through the Strecker step.

## 5. Caveats

- **No kinetics.** One temperature; no rate constant, Ea or pre-exponential anywhere. `arrhenius_params.yml`'s
  "Estimated from Hofmann 2000" for `strecker.A_value` is unsupported by this paper (and by the 2000b companion).
- **Atmosphere is controlled only in Table 4.** Tables 2, 3, 5 and all figures are "oxygen not controlled" (air in
  a closed vial or reflux); treat them as air but not as a defined oxygen charge.
- **Buffer-strength typos.** "0.5 mmol/L" in Experimental vs "0.5 mol/L" in footnotes; Table 4 says 0.1 mol/L.
  The 60-min MGO numbers in Table 3 (0.5 M, air) and Table 4 model B (0.1 M, air) agree within 4 % (PA 10.6 vs
  10.2; PAA 20.6 vs 20.4), so either buffer strength has no effect or the footnotes differ only on paper.
- **No replicates or uncertainties printed.** Treat ratios below ~1.3× (R4, R5, R7) as directional, not as
  numbers to fit tightly.
- **Table 3's "phenylalanine degraded (%)" column coincides almost exactly with Table 3 of the 2000b companion**
  (Glc/Phe, reflux): 9/18/36/49/61 % here vs 8.0/17.5/36.0/49.0/60.2 % there at 10/30/60/120/300 min, from two
  different pots (MGO vs glucose donor). Either the Phe series was reused, or the coincidence is real; do not use
  both as independent evidence of amino-acid consumption.
- **Cu²⁺ dose is single** (0.05 mmol in 10 mL = 5 mmol/L) and far above food levels; R7 is a direction, not a dose
  response.
- Figures 2–5 are not digitised; the pH numbers (2×, 7×, 3:1, 1:5) are the text's own readings of them.
- PA under air is consistently 7–15 % below argon for every donor (R4). The authors do not comment; it may be
  aldehyde loss to other oxidative sinks. Record, do not fit.
