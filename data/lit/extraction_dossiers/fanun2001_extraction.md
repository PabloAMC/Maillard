# Fanun, Leser, Aserin & Garti 2001 (Colloids Surf. A 194:175-187) — EXTRACTION (file stem `fanun2001`)
### *"Sucrose ester microemulsions as microreactors for model Maillard reaction"*
### L-cysteine + furfural (0.12 + 0.06 mmol) in W/O sucrose-ester microemulsions (sucrose stearate S1570 or laurate L1695 / phosphate buffer pH 4-8 / dodecane + C2-C6 alcohol), sealed tubes, **30-60 °C, up to 8 h**; an **aqueous reference ("1 g of aqueous furfural/cysteine mixture")** heated the same way. Quantified: **furfural conversion (%)** and ratios of three sulfur products; furfural and FFT by a **mesitylene internal standard that "does not entirely compensate for losses"**. **No sugar, no MFT. In the aqueous reference FFT is the main sulfur product but no amount is printed**; its furfural conversion is "Very slow" at every pH. Microemulsion furfural loss: pseudo-first-order k 0.003-0.012 min⁻¹ (30-60 °C), Ea 38.8 kJ/mol.
DOI: 10.1016/s0927-7757(01)00786-5 (printed on p. 175 as PII S0927-7757(01)00786-5)

**Read 2026-11-05 (first reader) from `data/articles/fanun2001.pdf`** (13 pp., Elsevier typeset PDF, clean text layer).
Prose and Tables 1-5 from `pdftotext -layout`; the tables are short and came through intact. Figs. 2-15 are a phase
diagram, schemes, an FPD chromatogram, mass spectra and conversion-vs-time curves; **the curves (Figs. 8-14) were not
digitised**: Tables 3 and 5 print the values the text uses. Arithmetic: scratchpad `tools/m4_arith.py` (maillard env,
`-I`).

**Role (fixed in the M4 case sweep, sec. 2.9: "not M4: a Vauthey sibling ..., passed to H0's question").** This dossier
judges it against H0's criterion only (sec. 8). Every aqueous-reference thiol statement below carries
**RESERVED (H0): unseen by the engine; no engine run on these pots until H0's design is frozen** (**[RESERVED H0]**).
No engine run on any pot of this paper.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[Q]** qualitative statement
by the authors · **[C]** cited · **[Z]** derived by me · **[!]** flag. Provenance: **first reader only, second reader
pending.** Nothing here is in a ledger, fit table or benchmark.

## 0. Identity

| field | value as printed (p. 175) |
|---|---|
| Authors | **M. Fanun**ᵃ, **M. Leser**ᵇ, **A. Aserin**ᵃ, **N. Garti**ᵃ* |
| Affiliations | ᵃCasali Institute of Applied Chemistry, The Hebrew University of Jerusalem; ᵇNestlé Research Center, Nestec Ltd., Lausanne |
| Venue | *Colloids and Surfaces A: Physicochemical and Engineering Aspects* **194** (2001) **175-187** |
| DOI | **10.1016/S0927-7757(01)00786-5** (from the printed PII), matches the brief |
| Dates | received 19 Mar 2001; accepted 14 May 2001 |

**Laboratory.** Garti group, Hebrew University (corresponding), with Leser (Nestlé, the colloid group of
`vauthey2000_extraction.md`). Not DFA. Not in `data/keys/laboratories.yml`. **Proposal, not applied:** group it with
Vauthey 2000 as one structured-fluid collaboration (Garti / Leser); unlike Vauthey, no DFA-derived isolation is used
(plain ether extraction, mesitylene standard).

## 1. What this paper answers for the model

| model need | what Fanun 2001 gives | tag |
|---|---|---|
| Absolute aqueous MFT / FFT with a stated method and recovery (H0) | **none**: no sugar pot, no MFT; aqueous Cys + furfural reference: FFT "mainly found", **no amount** | [Q] |
| Aqueous furfural loss with Cys | "Very slow" at pH 4, 5, 6, 8 (Table 4); no number | [M] |
| Microemulsion furfural loss vs T | k 0.003 / 0.004 / 0.0058 / 0.008 / 0.012 min⁻¹ at 30 / 37 / 43 / 50 / 60 °C; Ea 38.8 kJ/mol | [M] |
| Microemulsion furfural conversion vs time, ratio, pH | Table 3 (60 °C) | [M] |
| Product ratios in microemulsion | Table 2 (60 °C, 1 and 8 h) | [M] |

## 2. Conditions [M] (sec. 2.3, p. 178)

| item | microemulsion pots | aqueous reference |
|---|---|---|
| Charge | furfural **0.06 mmol** in alcohol-dodecane-sucrose ester (1 : 1 : 1.5, "43 dilution line"), diluted with phosphate buffer (pH 4-8) containing L-cysteine **0.12 mmol** | "**1 g of aqueous furfural/cysteine mixture**" [V]: concentrations, buffer and pH **not stated** [!] |
| Molar ratio | furfural : Cys 1 : 2 (also 2 : 1 and 1 : 1 in Table 3) | not stated |
| Water content | 8, 16, 30 (or 32) wt% of the microemulsion | water |
| Buffer / pH | phosphate (molarity not stated), **pH 4-8, label (as prepared)**; not re-measured | not stated |
| Vessel | sealed tube (culture tubes with viton-lined screw caps are named for phase work) | "the same heating procedure" |
| Heating | water bath, **30-60 °C, 8 h** (chosen to stay in the one-phase region, ~30-70 °C) | same |
| Headspace / O2 | not stated | not stated |
| Replicates | not stated | not stated |

## 3. The assay [M] (secs. 2.4-2.5)

| item | as printed |
|---|---|
| Isolation | diethyl ether extraction (sucrose ester co-extracted, decanted), MgSO4 |
| Quantification | furfural and FFT by **mesitylene (20 mg) added after extraction** as internal standard; "this method does not entirely compensate for losses during the work-up procedure" [V, paraphrase < 15 words] |
| SIDA | **none** |
| Recovery | **not measured** (the authors state losses are not compensated) |
| Response factors / calibration | not stated |
| LOD | **not stated** |
| GC | HP 5890 GC-MS, HP5 30 m; FPD signal shown for the microemulsion (Fig. 4); RI by van den Dool & Kratz |
| Identification | RI + EI mass spectra; 2-(2-furanyl)-thiazolidine (RI 1521) and N-(2-mercaptovinyl)-2-(2-furanyl)-thiazolidine (RI 1997, tentative) |
| Replicate noise | not stated |

## 4. Data

### 4.1 Aqueous reference [Q] / [M] [RESERVED H0]

- Text (p. 179): in the aqueous phase (reference) "one sulfur compound (FFT) was mainly found" [V].
- Table 4 (p. 183), water: furfural conversion rate "Very slow" at pH 4, 5, 6 and 8; products "A" (= FFT) at each pH.
- Conclusions: furfural + Cys in water at 60 °C gives one main sulfur compound, FFT; the microemulsion gives "higher
  amounts" of FFT [Q]. **No FFT amount, no yield, no conversion number for water is printed anywhere.**

### 4.2 Table 2 (p. 181): product ratios in microemulsions, 60 °C [M]

| water (wt%) | L1695 A/B (1 h) | L1695 B/C (1 h) | L1695 B/C (8 h) | S1570 A/B (1 h) | S1570 B/C (1 h) | S1570 B/C (8 h) |
|---|---:|---:|---:|---:|---:|---:|
| 8 | 0.02 | 2.00 | 0.07 | 0.02 | 2.00 | 0.15 |
| 16 | 0.05 | 1.64 | 0.10 | 0.03 | 1.50 | 0.35 |
| 32 | 0.15 | 1.20 | 0.30 | 0.07 | 1.00 | 0.45 |

A = furfurylthiol (FFT); B = "2-(2-furanyl)-thiazole" in the footnote; C = N-(2-mercaptovinyl)-2-(2-furanyl)-
thiazolidine. **[!]** The text calls the quantified pair 2-(2-furanyl)-thiazolidine and the mercaptovinyl compound, at
water contents **8, 16 and 30 wt%**; the table says thiazole and **32** wt%. FFT is ≤ 15 % of B in every microemulsion
cell, i.e. the thiazolidine route dominates there.

### 4.3 Table 3 (p. 182): furfural conversion (%) in the S1570 microemulsion, 60 °C [M]

| time (h) | 2/1, pH 5 | 1/1, pH 5 | 1/2, pH 5 | 2/1, pH 8 | 1/1, pH 8 | 1/2, pH 8 |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 48 | 60 | 65 | 70 | 67 | 67 |
| 4 | 58 | 72 | 77 | 85 | 83 | 81 |
| 8 | 69 | 82 | 87 | 93 | 92 | 90 |

(furfural / Cys molar ratio). **[!] Abstract: rates "decrease ... with the increase of pH from 4 to 8"; Table 3 shows
higher conversion at pH 8 than at pH 5 at every time, and sec. 3.1 says the rate "increased with increasing pH"**
(Fig. 9, "not significant").

### 4.4 Table 5 (p. 186): pseudo-first-order k for furfural loss, S1570 microemulsion, 16 wt% water, pH 5, furfural / Cys 1 / 2 [M]

| T (°C) | 30 | 37 | 43 | 50 | 60 |
|---|---:|---:|---:|---:|---:|
| k (min⁻¹) | 0.003 | 0.004 | 0.0058 | 0.008 | 0.012 |

Printed: **Ea 38.8 kJ/mol, k0 0.0029 min⁻¹** (rate law d[F]/dt = −k[F], with [H+] taken constant).

### 4.5 Arithmetic [Z]

1. Least-squares Arrhenius on Table 5: **Ea 39.7 ± 1.4 kJ/mol**; two-point 30 / 60 °C: 38.8 (the printed value).
   **[!] The printed k0 = 0.0029 min⁻¹ is wrong**: with Ea 38.8 it gives k(60 °C) = 1.7 × 10⁻⁹ min⁻¹; the fit's
   intercept is **k0 ≈ 2.0 × 10⁴ min⁻¹**. The k values are microemulsion-only (furfural partitioned into butanol).
2. Half-life of furfural at 60 °C in the microemulsion: ln 2 / 0.012 = **58 min**, i.e. 51 % converted at 1 h. Table 3
   at 60 °C (1/2, pH 5) reads 65 % at 1 h and is not first order over 1-8 h (65 → 77 → 87 %); its water content is not
   stated, so the two tables need not be the same pots.

## 5. What it can and cannot test in the model

**Cannot:** anything on MFT; any absolute aqueous FFT (none printed); any sugar-derived thiol; any rate in water (only
"very slow"); pH effects in water.

**Can (matrix evidence only, outside the validated domain by judgement: 30-60 °C W/O microemulsions are not an engine
matrix):** the microemulsion furfural loss with Cys at pH 5 (k 0.012 min⁻¹ at 60 °C, Ea ~39 kJ/mol), and that the
interface favours the thiazolidine route over FFT. Not an engine pot (no lane carries a microemulsion phase).

## 6. Conflicts and agreements with existing records

- **Agreement with `vauthey2000_extraction.md`** (same collaboration, Nestlé, 100 °C, 4 h): structured fluids raise FFT
  and make 2-(2-furyl)thiazolidine and its N-(2-mercaptovinyl) derivative; in water FFT is the main thiol. Vauthey
  prints an aqueous FFT amount (1.8 µg, uncorrected, reserved); Fanun prints none.
- **`yaghmur2002_extraction.md`** cites this paper ("sucrose-laurate W/O (Fanun 2001) favours D", in that dossier's
  own lettering, not checked here). Table 2 agrees that in sucrose-laurate W/O the thiazolidine products dominate FFT.
- **[!] Internal:** thiazole vs thiazolidine and 30 vs 32 wt% (sec. 4.2); the abstract's pH direction vs Table 3 (sec.
  4.3); k0 (sec. 4.5).
- No conflict with any ledger value.

## 7. Provenance tiers (first reader only)

Tables 2, 3, 5: `direct_measurement`, quant class relative (internal standard, losses uncorrected), confidence low;
matrix not an engine domain. Aqueous reference: `qualitative`. Derived Ea / k0: `derived`. **Proposed ledger entries:
none.**

## 8. H0 role

**H0's criterion (as given in the brief):** absolute MFT / FFT from an aqueous control, quantified with a stated method
(ideally SIDA) and recovery, from a non-DFA laboratory.

| element of the criterion | what the paper has | verdict |
|---|---|---|
| an aqueous control | yes: "1 g of aqueous furfural/cysteine mixture", 30-60 °C, 8 h; **no sugar**; concentrations and pH not stated | met (Cys + furfural, not a sugar pot) |
| absolute MFT or FFT in it | **no**: MFT not studied; FFT "mainly found", **no amount printed** | **not met** |
| stated quantification method | mesitylene internal standard after extraction (no SIDA) | stated, not SIDA |
| recovery | **not measured**; the authors say losses are not compensated | **not met** |
| non-DFA laboratory | yes: Hebrew University (Garti) with Nestlé (Leser) | met |

**Verdict: NOT MET.** There is no aqueous thiol number to reserve; the qualitative aqueous statements are marked
[RESERVED H0] all the same. Not an M4 source (the M4 sweep passed it to H0).
