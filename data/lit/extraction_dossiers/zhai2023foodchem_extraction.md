# Zhai et al. 2023 (Food Chemistry) — TTCA + extra xylose: the ORIGIN of the corpus's sulfur temperature ladder

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★★★ THIS IS THE PRIMARY SOURCE. Kang 2026's Table S4 — the repo's only sulfur temperature ladder — is **this paper's Table 1, TTCA columns**, published **3.5 years earlier**. And this paper states a DIFFERENT quantification basis for the same numbers. See §2 and §8.

---

# §0. ★ WRONG-FILE IDENTITY — REPORT FIRST

**THE TWO ZHAI FILES ARE SWAPPED RELATIVE TO THE K6a BRIEF.** `Zhai2023.pdf` (1 567 301 B) is
**this** paper — *Food Chemistry* `10.1016/j.foodchem.2022.134420` — which the brief expected under
the stem `Zhai2023b`. The JAFC paper `10.1021/acs.jafc.3c04166` is in `Zhai2023b.pdf`
(dossier `zhai2023jafc_extraction.md`). Dossiers here are named by paper identity, per the K5a
precedent for the two Kocadağlı files.

⚠️ **And the year is ambiguous by design.** Received 3 Jul 2022, accepted 25 Sep 2022, **online
28 Sep 2022**, but assigned to ***Food Chemistry* volume 404 (2023)**. Cite as **Zhai et al. 2023**
(the journal's own assignment, and the form Kang 2026 would have used) but **date the data to
September 2022** when reasoning about precedence.

---

## §1. IDENTITY `[M]`

| item | value | anchor |
|---|---|---|
| file | `data/articles/Zhai2023.pdf`, 1 567 301 B, 14 pages, Elsevier/Distiller, born-digital | `pdfinfo` |
| title | **"Reduced asynchronism between regenerative cysteine and fragments of deoxyosones promoting formation of sulfur-containing compounds through extra-added xylose and elevated temperature during thermal processing of 2-threityl-thiazolidine-4-carboxylic acid"** | p. 1 |
| authors | **Yun Zhai**, Xue Xia, Shibin Deng, Heping Cui, Khizar Hayat, **Xiaoming Zhang\***, **Chi-Tang Ho\*** | p. 1 |
| affiliations | State Key Laboratory of Food Science and Technology, **Jiangnan University**, Wuxi (Zhai, Xia, Deng, Cui, Zhang); Miami University Ohio (Hayat); Rutgers (Ho) | p. 1 |
| **DOI** | ★ **`10.1016/j.foodchem.2022.134420`** ✔ matches the brief's expectation for `Zhai2023b` | p. 1 footer |
| journal | ***Food Chemistry* 404 (2023) 134420** | running head |
| dates | Received **3 Jul 2022** · Revised **19 Sep 2022** · Accepted **25 Sep 2022** · Online **28 Sep 2022** | p. 1 |
| version | version of record, 0308-8146/© 2022 Elsevier Ltd | p. 1 |
| competing interests | none declared | p. 13 |
| SI | **NOT ON DISK.** Referenced: **Fig. S1** (intermediates/flavour map of Cys–Xyl MRPs), **Fig. S2** (Cys-alone degradation), **Table S1** (ROAV values), **Table S2** (MGO–GO–Cys validation model) | throughout |
| ⚠ funding | **no funding statement is printed** — unusual for this group, which acknowledges NSFC 32172330 in the 2023 JAFC companion | — |

**Note the institutional difference from the JAFC companion:** this work is **Jiangnan University**
(Zhang's lab); the 2023 JAFC paper's lead moved to **Zhejiang Gongshang**, and Kang 2026 is
Zhejiang Gongshang. The data travel with the first author, Yun Zhai.

---

## §2. ★★ QUANTIFICATION BASIS — **AND IT CONTRADICTS THE JAFC COMPANION FOR THE SAME NUMBERS**

### 2.1 What this paper says, verbatim (p. 3, §2.7)

> *"1,2-Dichlorobenzene (3 μL, 0.018 μg/μL in methanol) was selected as an internal standard and
> added into the mixture of 3 g sample and 2 g saturated sodium chloride solution. … **The detected
> compounds were quantified by the comparison of peak areas and the correction factor was 1.**"*

with the explicit formula

> **`Wi = f′ × (Ai × ms) / (As × V)`** — *"where Wi represented the obtained concentration (μg/L)
> of the detected sample, Ai meant the peak area of i, As was the peak value of
> 1,2-dichlorobenzene, ms was the corresponding mass of 1,2-dichlorobenzene, V and f′ were the
> sample volume and correction factor, respectively."*

### 2.2 ★ THE VERDICT

> **THIS IS SINGLE-INTERNAL-STANDARD SEMI-QUANTIFICATION WITH THE RESPONSE FACTOR SET TO 1.
> There is no calibration curve, no authentic-standard response factor, and no recovery
> correction anywhere in this paper.**

**⇒ EVERY μg/L VALUE IN THIS PAPER IS A PEAK-AREA RATIO DRESSED IN CONCENTRATION UNITS.**

**★ THE WAVE RULE, APPLIED.** A constant unknown response factor cancels in a ratio and therefore
in an Arrhenius slope. So:
- ✅ **Every Ea, fold-change and time-course SHAPE in §9 is legitimate**, compound by compound.
- ❌ **Absolute μg/L values are not**, and neither is any comparison of one compound's magnitude to
  another's (their true response factors differ by up to an order of magnitude on a DB-Wax/EI
  system), nor any comparison to another paper's absolute concentration.
- ❌ **The class subtotals are sums over compounds with different unknown response factors** — they
  are *not* physical totals and should be treated as an index, not a mass.

### 2.3 ★★ AND THIS IS A DIRECT CONTRADICTION WITH THE COMPANION PAPER

The 2023 JAFC companion reports the **identical numbers** (§8) and says of them:

> *"The flavor compounds were quantified by the corresponding **calibration curves (Table S2)
> constructed by the standards**."* (JAFC p. 14302)

**The same measurements are described as `f′ = 1` semi-quant in 2022 and as calibration-curve
absolute quantification in 2023.** Both cannot be true of numbers that agree to three decimals.
**The earlier, more specific statement — `f′ = 1` — is the one that describes how the numbers were
produced**, because the formula is printed and the correction factor is named.

**★ CONSEQUENCE FOR THE REPO — a correction to wave K5's Tier A / Tier B split.**
`kang2026_SI_extraction.md` §4d classified Kang's Table S4 rows as **Tier A** ("external standard
curve in Table S3") or **Tier B**, and concluded that ~93 % of the 140 °C sulfur subtotal was
Tier A, assigning it **± 15 %** uncertainty on that basis. **Since Kang's Table S4 numbers ARE these
numbers, produced with `f′ = 1`, the Tier A designation cannot apply to the temperature ladder.**
Kang's Table S3 calibration curves may well be real, but they were not used to generate S4 —
identical values cannot come from two different quantification routes.
**⇒ Recommended: treat the whole Kang/Zhai temperature ladder as Tier B — semi-quant, response
factor assumed 1 — and widen its uncertainty accordingly (× ÷ 2 to × ÷ 3 on absolute magnitude;
Ea and shape unaffected).**

### 2.4 The rest of the analytical method `[M]`
| parameter | value | note vs the JAFC companion |
|---|---|---|
| sample | **3 g sample + 2 g saturated NaCl** | ★ **the JAFC paper omits the salt** — a salting-out difference that changes headspace partition |
| IS | 1,2-dichlorobenzene, 3 μL of 0.018 μg/μL in methanol | identical |
| ⚠ pyrazine samples | ★ ***"the solutions for pyrazines detection were adjusted to the pH of 8"*** | ★ **a separate, pH-adjusted extraction for the N-heterocycles only** — so pyrazine values are NOT on the same headspace basis as the sulfur values. Not mentioned in the JAFC paper |
| fibre | CAR/PDMS/DVB 75 μm, headspace, **60 °C, 20 min, with agitation** | identical |
| desorption | 250 °C, 10 min, splitless | identical |
| column | **DB-WAX 213**, 30 m × 0.250 mm × 0.25 μm | JAFC says "DB-Wax" |
| oven | 40 °C (2 min) → 3 °C/min → 80 °C (3 min) → 2 °C/min → 120 °C (2 min) → 8 °C/min → 230 °C (2 min) | ★ **printed in full here; the JAFC paper defers to a reference.** Total run ≈ 51 min |
| carrier | He 99.9 %, 1 mL/min; ion source 230 °C; **full scan 50–500 m/z**, EI 70 eV | — |
| ID | NIST 17 + RI vs C₇–C₃₀ n-alkanes | identical |
| **non-volatiles** | TTCA, Cys, Xyl by **HPLC-ELSD with pure TTCA / Xyl / Cys as external standards** (absolute); α-dicarbonyls by **OPD derivatisation → HPLC-DAD, quinoxaline external standards** (absolute) | ★ **the non-volatile data ARE absolute; only the volatiles are semi-quant** |
| ⚠ | *"The specific parameters were not shown in this study"* (HPLC-ELSD) and *"The specific derivation and detection steps … were not shown in this study"* (dicarbonyls) | **neither non-volatile method is reproducible from this paper** `[NEG]` |
| replication | **n = 3**, mean ± SD, **Duncan's multiple comparison**, p < 0.05, SPSS 19.0 | ★ **and unlike the JAFC companion, the SDs ARE printed** (Table 1) |

### 2.5 μ / unit raster check
Elsevier PDF, clean text layer. Raster-verified at 900 dpi: Table 1 caption **`(μg/L)`**; Fig. 1a
axis and all heat-map cells (figure, no text layer — **the entire §5 grid is raster-read**);
Fig. 1b y-axis **`Concentration (μg/L)`**; Figs. 2–3 y-axes **`(mmol/L)`**; IS **`0.018 μg/μL`**;
fibre **`75 μm`**; column **`0.25 μm`**. **No μ→m corruption.** ✅

---

## §3. SYSTEM AND CONDITIONS `[M]`

| system | composition | temperatures | times |
|---|---|---|---|
| **TTCA (T)** | TTCA **10 mmol/L**, water, **pH 7 ± 0.01** set with NaOH (2 and 6 mol/L) | **100, 120, 140 °C** | **20, 40, 60, 80, 100, 120 min** (Table 1 at 120 min; Fig. 1a extends to **140 min** at 100 °C) |
| **TTCA-Xyl (T-X)** | TTCA 10 mmol/L **+ Xyl 10 mmol/L** | 100, 120, 140 °C | as above |
| **TTCA-[¹³C₅-Xyl]** | TTCA + **[¹³C₅]-D-xylose (99 %)** 10 mmol/L | 100, 120, 140 °C | **120 min only** |
| **Cys alone** | Cys **10 mmol/L**, pH 7 | 100, 120, 140 °C | 20–120 min |
| **fresh MRP control** | equimolar Cys + Xyl **0.0827 mol/L** | **120 °C** | **120 min** |
| **MGO-GO-Cys** | MGO 10 + GO 10 + Cys 10 mmol/L, pH 7.0 ± 0.1 | 100, 120, 140 °C | 20–120 min (results in **Table S2, not on disk**) |

- Vessel: **pressure-resistant bottles**, collector-type magnetic stirrer **DF-101S**; ice-bath quench.
- **TTCA synthesis:** Cys + Xyl 0.0827 mol/L, pH 7.4 ± 0.01, **90 °C / 40 min**, Dowex 50WX4 (H⁺)
  then semi-prep RP-HPLC on XBridge Amide; UPLC-ESI-MS + 400 MHz NMR confirmation.
- ⚠️ **UNBUFFERED.** pH set with NaOH only. **`[NEG]` no buffer anywhere.**
- ★ **BUT THE pH DRIFT IS REPORTED HERE, unlike in the JAFC companion:**
  > *"3-DX … was the dominant compound in the reaction system compared with 1-DX demonstrated by
  > the **rapid decrease of pH from 7.0 to 4.9**, which was beneficial to the 1,2-enolization and
  > generating the 3-DX."* (p. 5)
  ★★ **This is the origin of the "measured pH 4.9" that `kinetic_core_b2_2_diagnosis.md` §3 uses
  as its out-of-sample comparator, attributed there to "Kang p. 3242". It is Zhai's number, from
  2022, and it is stated without a temperature or a time.** The B2.2 diagnosis's caveat
  ("the temperature of that pH series is `100_or_120_C_UNRESOLVED`") is if anything understated:
  **in the primary source it is not attached to any single condition at all.**

---

## §4. ★★ TABLE 1 — the temperature ladder, complete, WITH SDs and Duncan letters `[M]`

**Caption, verbatim:** *"Volatile flavor compounds (μg/L) identified in the TTCA and TTCA-Xyl model
reactions under elevated temperatures of 100, 120 and 140 °C at an initial pH value of 7.0."*
*"Results were presented as means ± standard deviation, data within a row with different letters
are significantly different (p < 0.05) using Duncan's multiple comparison test (n = 3)."*
`–` = not detected. `RIᵃ` = measured on DB-WAX; `RIᵇ` = literature (flavornet / NIST WebBook).
**Time = 120 min throughout.**

### Sulfur-containing compounds (25 rows)

| compound | RIᵃ | RIᵇ | **100 °C TTCA** | **100 °C TTCA+Xyl** | **120 °C TTCA** | **120 °C TTCA+Xyl** | **140 °C TTCA** | **140 °C TTCA+Xyl** |
|---|---|---|---|---|---|---|---|---|
| Thiophene | 1016 | 1022 | – | – | – | 0.036 ± 0.009ᵇᶜ | 0.118 ± 0.072ᵃᵇ | 0.219 ± 0.089ᵃ |
| 2-Methylthiophene | 1135 | 1095 | – | – | ⚠ **1.128 ± 0.126ᵇ** | ⚠ **1.128 ± 0.126ᵇ** | 0.892 ± 0.115ᶜ | 2.519 ± 0.09ᵃ |
| 2-Methylthiazole | 1272 | 1278 | 0.142 ± 0.008ᶜ | 0.173 ± 0.027ᶜ | 0.194 ± 0.014ᶜ | 1.347 ± 0.129ᵇ | 0.254 ± 0.109ᶜ | 2.426 ± 0.082ᵃ |
| Thiazole | 1240 | 1262 | 0.896 ± 0.053ᶠ | 1.214 ± 0.086ᵉ | 1.796 ± 0.12ᵈ | 2.798 ± 0.087ᵇ | 2.038 ± 0.139ᶜ | 4.408 ± 0.071ᵃ |
| 3-Mercapto-2-butanone | 1273 | 1283 | 1.303 ± 0.153ᵈ | 1.673 ± 0.115ᶜ | 1.798 ± 0.127ᶜ | 3.769 ± 0.274ᵃ | 2.398 ± 0.141ᵇ | 3.429 ± 0.099ᵃ |
| 2,5-Dimethyl-thiazole | 1301 | 1326 | 0.079 ± 0.013ᵉ | 0.245 ± 0.013ᵉ | 2.203 ± 0.142ᵈ | 3.038 ± 0.076ᵇ | 2.553 ± 0.096ᶜ | 4.165 ± 0.188ᵃ |
| **2-Methyl-3-furanthiol (MFT)** | 1302 | 1305 | **1.237 ± 0.142ᶜ** | **1.729 ± 0.152ᶜ** | **1.388 ± 0.256ᶜ** | **2.498 ± 0.13ᶜ** | **5.907 ± 0.085ᵃᵇ** | **3.238 ± 0.735ᵇ** |
| 2-Ethyl-thiazole | 1319 | 1304 | 0.245 ± 0.035ᵇ | 0.318 ± 0.063ᵇ | 0.273 ± 0.006ᵇ | 0.113 ± 0.011ᶜ | 0.289 ± 0.025ᵇ | 0.432 ± 0.088ᵃ |
| 3-Mercapto-2-pentanone | 1352 | 1343 | 0.973 ± 0.099ᶜ | 1.129 ± 0.007ᶜ | 1.298 ± 0.253ᶜ | 2.286 ± 0.142ᵃ | 1.736 ± 0.234ᵇ | 2.103 ± 0.208ᵃᵇ |
| 2,4,5-Trimethylthiazole | 1373 | 1390 | 0.419 ± 0.052ᵈ | 0.537 ± 0.061ᵈ | 3.199 ± 0.075ᶜ | 4.367 ± 0.095ᵇ | 3.843 ± 0.188ᵇ | 4.999 ± 0.601ᵃ |
| **2-Furfurylthiol (FFT)** | 1426 | 1402 | **3.734 ± 0.085ᵈ** | **4.736 ± 0.639ᶜ** | **4.107 ± 0.137ᶜᵈ** | **6.398 ± 0.266ᶜᵈ** | **11.439 ± 0.265ᵃ** | **6.123 ± 0.235ᵇ** |
| 5-Ethyl-2,4-dimethyl-thiazole | 1437 | – | – | 0.129 ± 0.013ᵉ | 0.587 ± 0.078ᵈ | 2.179 ± 0.064ᵇ | 0.889 ± 0.085ᶜ | 3.038 ± 0.151ᵃ |
| 2-Ethyl-4-methyl-thiazole | 1449 | 1410 | 0.117 ± 0.005ᵉ | 0.437 ± 0.063ᵉ | 1.278 ± 0.117ᵈ | 3.389 ± 0.271ᵇ | 2.112 ± 0.370ᶜ | 5.473 ± 0.521ᵃ |
| **3-Thiophenethiol** | 1564 | 1530 | – | – | 0.178 ± 0.027ᶜ | ★ **7.389 ± 0.289ᵇ** | 0.268 ± 0.062ᶜ | ★ **11.497 ± 0.649ᵃ** |
| 2-Acetylthiazole | 1646 | 1643 | 3.079 ± 0.215ᶜ | 4.024 ± 0.146ᶜ | 8.795 ± 0.094ᵇ | 23.397 ± 0.405ᵇ | 9.858 ± 1.215ᵇ | ★ **29.478 ± 0.415ᵃ** |
| 2-Thiophenecarboxaldehyde | 1674 | 1679 | 0.09 ± 0.013ᶜ | 0.179 ± 0.012ᶜ | 2.196 ± 0.25ᵇ | 7.774 ± 0.251ᵇ | 2.28 ± 0.144ᵇ | **18.541 ± 0.823ᵃ** |
| 5-Methyl-2-thiophenecarboxaldehyde | 1701 | 1785 | 0.936 ± 0.06ᶜ | 1.178 ± 0.294ᶜ | 3.649 ± 0.099ᵇ | 5.237 ± 0.087ᵇ | 5.706 ± 0.644ᵇ | **16.309 ± 2.416ᵃ** |
| 2-Thiophenemethanethiol | 1702 | 1713 | 0.119 ± 0.074ᶜ | 0.146 ± 0.023ᶜ | 0.196 ± 0.015ᶜ | 0.779 ± 0.152ᵇ | 0.332 ± 0.21ᶜ | 2.479 ± 0.182ᵃ |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio)furan | 1732 | – | 0.021 ± 0.004ᶜ | 0.178 ± 0.039ᵃᵇ | 0.098 ± 0.033ᵇᶜ | 0.178 ± 0.031ᵃᵇ | 0.137 ± 0.036ᵃᵇ | 0.217 ± 0.080ᵃ |
| 1,2,3-Trithiolane | 1794 | – | 0.208 ± 0.029ᵇ | 1.158 ± 0.018ᵃ | 0.378 ± 0.071ᵇ | 1.278 ± 0.233ᵃ | 0.497 ± 0.185ᵇ | 1.479 ± 0.289ᵃ |
| 3,3′-Dithiobisthiophene | 1845 | – | – | – | – | 0.037 ± 0.011ᵇ | – | 0.178 ± 0.025ᵃ |
| Thieno[3,2-b]thiophene | 1868 | 1843 | – | – | 1.302 ± 0.079ᶜ | 5.598 ± 0.163ᶜ | 3.369 ± 0.454ᵇ | ★ **23.979 ± 0.766ᵃ** |
| 2,5-Thiophenedicarboxaldehyde | 1911 | – | 0.346 ± 0.039ᵃ ⚠ | 0.793 ± 0.107ᵈ | 0.724 ± 0.139ᵈᵉ | 2.178 ± 0.059ᵇ | 1.283 ± 0.122ᶜ | 4.728 ± 0.405ᵃ |
| 2-Methylthieno[2,3-b]thiophene | 1947 | – | – | – | 0.19 ± 0.03ᶜ | 4.379 ± 0.113ᶜ | 2.123 ± 0.221ᵇ | **11.635 ± 0.255ᵃ** |
| Bis(2-furfuryl)sulfide | 2419 | – | 0.034 ± 0.007ᵇ | 0.215 ± 0.083ᵃᵇ | 0.042 ± 0.004ᵇ | 0.237 ± 0.113ᵃᵇ | 0.078 ± 0.018ᵇ | 0.339 ± 0.158ᵃ |
| **Subtotal (printed)** | | | **13.978 ± 0.568ᶠ** | **20.19 ± 0.514ᵉ** | ⚠ **35.866 ± 1.382ᵈ** | **91.807 ± 1.67ᵇ** | **60.4 ± 1.62ᶜ** | **163.432 ± 3.102ᵃ** |
| **Subtotal recomputed `[D]`** | | | **13.978** ✔ | **20.191** ✔ | ⚠ **36.997** | ⚠ **91.807** ✔ | **60.399** ✔ | **163.431** ✔ |

### Nitrogen-containing heterocycles (4 rows)
| compound | RIᵃ | RIᵇ | 100 T | 100 T+X | 120 T | 120 T+X | 140 T | 140 T+X |
|---|---|---|---|---|---|---|---|---|
| Pyrazine | 1216 | 1209 | – | – | – | 0.348 ± 0.061ᵇ | – | 1.389 ± 0.071ᵃ |
| Methylpyrazine | 1214 | 1263 | – | – | – | 1.039 ± 0.059ᵇ | – | 2.315 ± 0.129ᵃ |
| 3-Methyl-pyridine | 1307 | 1346 | – | – | – | – | – | 0.116 ± 0.072ᵃ |
| 2,5-Dimethyl-pyrazine | 1320 | 1328 | – | – | – | 2.578 ± 0.106ᵇ | – | 6.475 ± 0.741ᵃ |
| **Subtotal** | | | **0.000** | **0.000** | **0.000** | **3.965 ± 0.132ᵇ** | **0.000** | **10.296 ± 0.755ᵃ** |

★★ **The N-heterocycle block is EXACTLY ZERO in the TTCA-alone arm at all three temperatures, and
non-zero only when exogenous xylose is present.** This is the sharpest structural statement in the
paper and it is a hard constraint: **pyrazine formation in a TTCA pot is limited by α-dicarbonyl
supply, not by amine supply** — the Cys is there either way.

### Oxygen-containing heterocycles (5 rows)
| compound | RIᵃ | RIᵇ | 100 T | 100 T+X | 120 T | 120 T+X | 140 T | 140 T+X |
|---|---|---|---|---|---|---|---|---|
| 2-Methyl-furan | 851 | 829 | – | – | 0.019 ± 0.007ᶜ | 1.182 ± 0.007ᵇ | 0.133 ± 0.088ᶜ | 2.578 ± 0.249ᵃ |
| Furan | 797 | 798 | – | – | 0.138 ± 0.017ᵈ | 0.298 ± 0.037ᶜ | 1.108 ± 0.092ᵇ | 2.319 ± 0.133ᵃ |
| **Furfural** | 1457 | 1460 | **3.381 ± 0.089ᵉ** | 4.128 ± 0.139ᵉ | **5.793 ± 0.422ᵈ** | 8.497 ± 0.514ᶜ | **11.039 ± 0.302ᵇ** | **26 ± 0.725ᵃ** |
| 1-(2-Furanyl)-ethanone | 1497 | 1501 | – | – | 0.136 ± 0.017ᶜ | 0.178 ± 0.011ᶜ | 0.379 ± 0.087ᵇ | 0.867 ± 0.173ᵃ |
| 2(5H)-Furanone | 1748 | 1767 | – | – | 0.072 ± 0.012ᵇ | 0.189 ± 0.016ᵇ | 0.098 ± 0.015ᵇ | 0.462 ± 0.261ᵃ |
| **Subtotal** | | | **3.381 ± 0.089ᵉ** | **4.128 ± 0.139ᵉ** | **6.158 ± 0.425ᵈ** | **10.344 ± 0.48ᶜ** | **12.757 ± 0.469ᵇ** | **32.226 ± 0.812ᵃ** |
| **TOTAL (printed)** | | | **17.359 ± 0.527ᶠ** | **24.318 ± 0.517ᵉ** | ⚠ **42.024 ± 1.737ᵈ** | **106.116 ± 2.239ᵇ** | **73.157 ± 1.557ᶜ** | **205.954 ± 2.4ᵃ** |
| **TOTAL recomputed `[D]`** | | | **17.359** ✔ | **24.319** ✔ | ⚠ **43.155** | **106.116** ✔ | **73.156** ✔ | **205.952** ✔ |

### 4.1 ★ THE ONE ARITHMETIC DEFECT, DIAGNOSED — and it is the seed of a corpus-wide error

**Recomputed vs printed, all 6 subtotals + 6 totals: 10 of 12 reproduce exactly. The two failures
are the same failure**, the 120 °C TTCA sulfur subtotal (36.997 vs 35.866) and total
(43.155 vs 42.024), both short by **1.131 ≈ 1.128 + rounding**.

★ **The 2-Methylthiophene 120 °C row reads `1.128 ± 0.126ᵇ` in BOTH the TTCA and the TTCA+Xyl
column — identical mean, identical SD, and identical Duncan letter.** Two independent triplicate
measurements do not produce identical SDs. **This is a copy-paste of the TTCA+Xyl cell into the
TTCA cell**, and the printed subtotal — computed before the slip, from a TTCA value of `–` —
is the *correct* one.

**⇒ Diagnosis: the true value of 2-methylthiophene, TTCA arm, 120 °C is `–` (not detected). The
printed `1.128` is spurious.**

**⇒ AND THIS RESOLVES THE ONE DISAGREEMENT BETWEEN THIS PAPER AND KANG 2026** (§8): Kang's
Table S4 prints `0.000` in that cell. **Kang corrected the slip; the JAFC 2023 companion's
Fig. 3b propagated it.** `zhai2023jafc_extraction.md` §7.3 point 3 flagged the cell as "one of the
two has been altered" — this dossier identifies which, and in which direction: **Kang is right,
the two Zhai papers carry the artefact.** No hold-out row is affected (the cell is not used in any
declared row), but the corpus should record the correct value as `not detected`.

⚠️ **One further typographic defect:** 2,5-thiophenedicarboxaldehyde's 100 °C TTCA cell carries
Duncan letter **ᵃ** while its value (0.346) is the smallest in the row — every other row runs
`f < e < d < c < b < a` upward. **Letter error only.**

⚠️ **And a text/table mismatch:** the Conclusion states the TTCA+Xyl pyrazine total at 140 °C
*"reached … 10.020 μg/L"*; **Table 1's own N subtotal is 10.296.** Use the table.

---

## §5. ★★ FIGURE 1a — the 100 °C TIME COURSE GRID, 41 compounds × 15 columns `[M]`

A second numeric heat map: **41 compounds × (7 time points × 2 systems + 1 MRP control) = 615
values**, all at **100 °C**. Entirely raster-read at 900 dpi.
**Column key** (printed under the panel): `A` T-20 · `B` T-X-20 · `C` T-40 · `D` T-X-40 ·
`E` T-60 · `F` T-X-60 · `G` T-80 · `H` T-X-80 · `I` T-100 · `J` T-X-100 · `K` T-120 ·
`L` T-X-120 · `M` T-140 · `N` T-X-140 · **`O` MRPs** (fresh Cys+Xyl, 120 °C/120 min).

⚠️ **This figure's compound list is NOT the Table 1 list.** It has **41 rows**, adds seven species
(2,5-dimethyl-thiophene, dihydro-2-methyl-3(2H)-thiophenone, dihydro-3-(2H)-thiophenone,
2-acetyl-3-methylthiophene, 2,3′-bithiophene, pyrrole, benzofuran) and orders them differently.
**Do not index the two by row number.**

| # | compound | A | B | C | D | E | F | G | H | I | J | **K** | L | M | N | **O MRP** |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.079 |
| 2 | 2,3-Dihydro-5-methyl-thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.035 |
| 3 | 2,5-Dimethyl-thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.065 |
| 4 | 2-Methylthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.141 |
| 5 | Thiazole | 0 | 0.139 | 0.097 | 0.224 | 0.186 | 0.424 | 0.289 | 0.709 | 0.555 | 0.999 | **0.896** | 1.214 | 0.698 | 1.309 | 0.063 |
| 6 | 2-Methylthiazole | 0 | 0 | 0.001 | 0.005 | 0.007 | 0.010 | 0.022 | 0.043 | 0.038 | 0.749 | **0.142** | 0.173 | 0.119 | 0.204 | 0.045 |
| 7 | 3-Mercapto-2-butanone | 0 | 0 | 0 | 0.028 | 0.079 | 0.129 | 0.198 | 0.349 | 0.239 | 0.386 | **1.303** | 1.673 | 1.529 | 2.018 | 0.728 |
| 8 | 2,5-Dimethyl-thiazole | 0 | 0 | 0 | 0.079 | 0.009 | 0.119 | 0.028 | 0.186 | 0.063 | 0.217 | **0.079** | 0.245 | 0.084 | 0.379 | 0.083 |
| **9** | **2-Methyl-3-furanthiol (MFT)** | 0 | 0.083 | 0.078 | 0.213 | 0.100 | 0.334 | 0.509 | 0.939 | 0.980 | 1.537 | **1.237** | 1.729 | 1.179 | 1.928 | 0.477 |
| 10 | 2-Ethyl-thiazole | 0 | 0 | 0 | 0 | 0.128 | 0.197 | 0.163 | 0.238 | 0.211 | 0.297 | **0.245** | 0.318 | 0.218 | 0.379 | 0.023 |
| 11 | 3-Mercapto-2-pentanone | 0.079 | 0.129 | 0.127 | 0.308 | 0.308 | 0.598 | 0.554 | 0.887 | 1.098 | 1.102 | **0.973** | 1.129 | 0.997 | 1.428 | 0.123 |
| 12 | 2,4,5-Trimethylthiazole | 0 | 0.004 | 0.004 | 0.133 | 0.100 | 0.833 | 0.392 | 2.204 | 1.718 | 0.179 | **0.419** | 0.537 | 0.399 | 1.012 | 0.050 |
| **13** | **2-Furfurylthiol (FFT)** | 0.253 | 0.372 | 0.708 | 1.218 | 0.842 | 1.397 | 1.206 | 2.101 | 3.022 | 3.499 | **3.734** | 4.736 | 3.029 | 4.875 | 0.064 |
| 14 | 5-Ethyl-2,4-dimethyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.018 | 0 | 0.079 | 0 | 0.129 | 0 | 0.179 | 0.113 |
| 15 | 2-Ethyl-4-methyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0.078 | 0.012 | 0.218 | 0.097 | 0.329 | **0.117** | 0.437 | 0.129 | 0.608 | 0.057 |
| 16 | Dihydro-2-methyl-3(2H)-thiophenone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.144 |
| 17 | Dihydro-3-(2H)-thiophenone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.095 |
| 18 | 3-Thiophenethiol | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.013 |
| 19 | 2-Acetylthiazole | 0.477 | 0.728 | 0.634 | 0.998 | 0.708 | 1.794 | 0.814 | 3.084 | 2.488 | 3.477 | **3.079** | 4.024 | 3.398 | 4.798 | 1.299 |
| 20 | 2-Thiophenecarboxaldehyde | 0 | 0 | 0 | 0 | 0 | 0.028 | 0 | 0.089 | 0.019 | 0.124 | **0.090** | 0.179 | 0.079 | 0.217 | 0.105 |
| 21 | 5-Methyl-2-thiophenecarboxaldehyde | 0 | 0 | 0 | 0 | 0 | 0.489 | 0.235 | 1.130 | 0.726 | 0.714 | **0.936** | 1.178 | 0.674 | 1.380 | 0.177 |
| 22 | 2-Thiophenemethanethiol | 0 | 0 | 0 | 0.004 | 0 | 0.013 | 0 | 0.098 | 0 | 0.129 | **0.119** | 0.146 | 0.114 | 0.217 | 0.079 |
| 23 | 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.079 | **0.021** | 0.178 | 0.035 | 0.218 | 0.015 |
| 24 | 2-Acetyl-3-methylthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.113 |
| 25 | 1,2,3-Trithiolane | 0 | 0 | 0 | 0 | 0 | 0.336 | 0 | 0.599 | 0.156 | 0.819 | **0.208** | 1.158 | 0.179 | 1.298 | 0.200 |
| 26 | 3,3′-Dithiobisthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.009 |
| 27 | Thieno[3,2-b]thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.012 | 0.238 |
| 28 | 2,5-Thiophenedicarboxaldehyde | 0 | 0 | 0 | 0 | 0 | 0.124 | 0.095 | 0.238 | 0.173 | 0.424 | **0.346** | 0.793 | 0.129 | 0.129 | 0.095 |
| 29 | 2-Methylthieno[2,3-b]thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.023 | 0.117 |
| 30 | 2,3′-Bithiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.017 |
| 31 | Bis(2-furfuryl)sulfide | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.114 | **0.034** | 0.215 | 0.179 | 0.449 | 0.072 |
| 32 | Pyrazine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.017 |
| 33 | Methylpyrazine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.213 |
| 34 | 3-Methyl-pyridine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.129 |
| 35 | 2,5-Dimethyl-pyrazine | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.033 |
| 36 | Pyrrole | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.048 |
| 37 | Furan | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.079 |
| 38 | 2-Methyl-furan | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.128 |
| **39** | **Furfural** | 0.214 | 0.428 | 0.307 | 1.106 | 0.351 | 2.976 | 2.891 | 4.028 | 3.315 | 4.894 | **3.381** | 4.128 | 3.381 | 4.683 | **5.566** |
| 40 | Benzofuran | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.124 | 0 | 0 | 0 | 0 | 0.248 |
| 41 | 4-Hydroxy-5-methyl-3(2H)-furanone | 0 | 0.214 | 0.103 | 0 | 0.142 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| — | **COLUMN TOTAL `[D]`** | **1.023** | 2.097 | 2.059 | 4.316 | 2.960 | 9.879 | 7.408 | 17.158 | 14.898 | 20.271 | **17.359** | 24.319 | **16.549** | 27.743 | **11.395** |

### 5.1 Internal-consistency audit of the grid `[Z]` — **passes on every quotable number**
| check | recomputed | paper's text | Δ |
|---|---:|---:|---|
| TTCA total at 20 min | **1.023** | *"1.023 μg/L at 20 min"* (p. 3) | **0.000** ✔ |
| TTCA total at 120 min | **17.359** | *"17.359 μg/L at 120 min"* | **0.000** ✔ |
| TTCA sulfur at 20 min | **0.809** | *"the sulfur-containing compounds increased from 0.809 μg/L"* | **0.000** ✔ |
| TTCA sulfur at 120 min | **13.978** | *"to 13.978 μg/L"* | **0.000** ✔ |
| Column K vs **Table 1's 100 °C TTCA column** | identical, all 34 shared rows | — | **0.000** ✔ |
| Column L (T-X-120) vs Table 1's 100 °C TTCA+Xyl total | **24.319** | Table 1 prints 24.318 | 0.001 (rounding) ✔ |
| MRP total | **11.395** | *"the total concentration … from TTCA model at 120 min was **1.56 times** higher than that generated from MRPs"* | 17.359/11.395 = **1.523**; the paper says 1.56 — ⚠ a **2 % discrepancy**, small enough to be a rounding or a different total basis, recorded but not material |

★ **Note the MRP column is at 120 °C** while the rest of the panel is at 100 °C — the caption does
not say so; only §2.3 of the Methods does. **A reader comparing column O to column K is comparing
different temperatures.** `[!]`

### 5.2 ★ What the grid shows that Table 1 cannot
| finding | numbers |
|---|---|
| ★ **The TTCA arm PEAKS at 120 min at 100 °C and declines by 140 min** | totals **14.898 → 17.359 → 16.549** at 100/120/140 min; MFT **0.980 → 1.237 → 1.179**; FFT **3.022 → 3.734 → 3.029** |
| ★ **The TTCA+Xyl arm does NOT peak in the window** | totals **20.271 → 24.319 → 27.743** — still rising at 140 min |
| ★ **so the turnover is a substrate-limitation effect, and extra sugar postpones it** | this is the mechanism the JAFC companion shows at 140 °C, seen here at 100 °C |
| **Extra Xyl roughly doubles output at every time point** | T-X/T = 2.05 / 2.10 / 3.34 / 2.32 / 1.36 / 1.40 / 1.68 at 20–140 min |
| ★ **The fresh-MRP control has a completely different profile** | 40 species vs 19; MRP's largest peak is **furfural 5.566** (49 % of its total) whereas TTCA's is FFT; MRP has **detectable pyrazine, methylpyrazine, 3-methylpyridine, 2,5-dimethylpyrazine, pyrrole** — all **zero** in every TTCA and TTCA-Xyl column at 100 °C |
| **MRP total (11.395) < TTCA total (17.359)** | the paper's headline: the intermediate out-yields the fresh reaction **1.5×** on intensity while losing **half** the species count |

---

## §6. ★ TABLE 2 — the ¹³C₅-xylose isotope tracing, complete `[M]`

**Caption:** *"Isotope distribution patterns of the flavor compounds … derived from thermal
reaction of TTCA-[¹³C₅]-Xylose model under temperatures of 100, 120 and 140 °C, pH value of 7."*
Values are **% of the isotopomer population**; `M+` is the unlabelled molecular ion m/z.
**Reaction time 120 min.**

| compound | M⁺ (m/z) | T (°C) | **0** | 1 | **2** | **3** | 4 | **5** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **2-Acetylthiazole** | 127 | 100 | **56** | 0 | 0 | **41** | 0 | <1 |
| | | 120 | 76 | 0 | 3 | 16 | 0 | 5 |
| | | 140 | **54** | 0 | **11** | **21** | 0 | **14** |
| **2-Thiophenecarboxaldehyde** | 112 | 100 | 76 | 0 | 0 | 8 | 0 | 16 |
| | | 120 | 59 | 0 | 0 | 17 | 0 | 24 |
| | | 140 | **48** | 0 | 0 | **32** | 0 | **20** |
| **5-Methyl-2-thiophenecarboxaldehyde** | 126 | 100 | 81 | 0 | 14 | 0 | 5 | 0 |
| | | 120 | 63 | 0 | 16 | 0 | 21 | 0 |
| | | 140 | **48** | 0 | **25** | 0 | **27** | 0 |
| **3-Thiophenethiol** | 116 | 100 | **100** | 0 | 0 | 0 | 0 | 0 |
| | | 120 | 81 | 0 | **19** | 0 | 0 | 0 |
| | | 140 | **66** | 0 | **34** | 0 | 0 | 0 |
| **Thieno[3,2-b]thiophene** | 140 | 100 | **100** | 0 | 0 | 0 | 0 | 0 |
| | | 120 | 79 | 0 | **21** | 0 | 0 | 0 |
| | | 140 | **55** | 0 | **45** | 0 | 0 | 0 |
| **2-Methylthieno[2,3-b]thiophene** | 154 | 100 | 99 | 0 | 0 | 0 | 0 | <1 |
| | | 120 | 79 | 0 | 0 | 0 | 0 | **21** |
| | | 140 | **45** | 0 | 0 | 0 | 0 | **45** |

### 6.1 ★★ What this table is, and why it is the most transferable thing in the paper `[D]`

**It is a direct, temperature-resolved measurement of how much of each product's carbon skeleton
comes from the EXOGENOUS sugar rather than from the intermediate's own.** The "0-label" column is
the fraction built entirely from TTCA's own (unlabelled) threityl chain.

**Exogenous-carbon incorporation, 1 − (0-label fraction) `[D]`:**

| compound | 100 °C | 120 °C | 140 °C | **the pattern** |
|---|---:|---:|---:|---|
| 2-Acetylthiazole | **44 %** | 24 % | **46 %** | non-monotone ⚠ |
| 2-Thiophenecarboxaldehyde | 24 % | 41 % | **52 %** | monotone ↑ |
| 5-Methyl-2-thiophenecarboxaldehyde | 19 % | 37 % | **52 %** | monotone ↑ |
| **3-Thiophenethiol** | **0 %** | 19 % | **34 %** | ★ switched on between 100 and 120 °C |
| **Thieno[3,2-b]thiophene** | **0 %** | 21 % | **45 %** | ★ switched on between 100 and 120 °C |
| **2-Methylthieno[2,3-b]thiophene** | **<1 %** | 21 % | **45 %** | ★ switched on between 100 and 120 °C |

★★ **Three compounds go from ZERO exogenous carbon at 100 °C to 19–21 % at 120 °C and 34–45 % at
140 °C. All three are the compounds that Table 1 shows are absent at 100 °C and huge at 140 °C in
the +Xyl arm.** ⇒ **There is a genuine, isotopically-verified temperature threshold between 100 and
120 °C for the entry of exogenous sugar carbon into the Strecker/retro-aldol branch.**

**⇒ THIS IS A REAL SWITCH-ON, AND IT IS NOT THE ONE KANG CLAIMED.** It sits at **100→120 °C**, not
120→140 °C; it applies to **thiophene-ring products of the retro-aldol/Strecker route**, not to the
free furanthiols MFT and FFT (which are **not in this table at all** — the authors did not trace
them); and it is measured by **isotope ratio**, a quantity that is immune to response factors, hold
time and headspace partitioning — i.e. immune to every artefact that §9 shows contaminates the
yield ladders. **It is the single most robust temperature-structure observation in the K6a corpus.**

⚠️ **Caveats:** (i) one hold time only (120 min), so a time-slice artefact cannot be excluded —
but an isotope *ratio* is far less time-sensitive than a yield; (ii) no SDs, no n stated for this
table; (iii) 2-acetylthiazole's non-monotone 44/24/46 is unexplained and should be treated as
noise-dominated; (iv) `<1 %` and `0` are not distinguished from below-detection.

---

## §7. FIGURES 2 AND 3 — the non-volatile time courses `[M]`, `[D]` (digitised, coarse)

**Digitisation basis:** pages 5 and 8 at 200–600 dpi; **peak values and peak times only**, ±10–15 %.
Six series per panel (T and T+X at 100/120/140 °C); at the available raster **the marker glyphs
separating 120 °C from 140 °C were not reliably resolvable**, so series are assigned by the
peak-time ordering argument (higher T ⇒ earlier, sharper peak) and flagged as such.
**All four are ABSOLUTE (external-standard) measurements, unlike the volatiles.**

| panel | species | scale | **peak, TTCA arm** | **peak, TTCA+Xyl arm** | peak time |
|---|---|---|---|---|---|
| 2a | **3-DX** | 0–1.5 mmol/L | **≈ 0.82** | ★ **≈ 1.13** | T+X peaks **60 min**; T peaks **80–100 min** |
| 2b | **1-DX** | 0–0.25 mmol/L | ≈ 0.11 | ≈ 0.175 | 80–100 min |
| 3a | **regenerative Cys** | 0–1.5 mmol/L | ≈ 1.0 (plateau from ~90 min) | ≈ 1.0, **reached earlier and falling faster** | — |
| 3b | **diacetyl (DA)** | 0–0.004+ mmol/L | ≈ 0.0029 | ★ **≈ 0.0040, peaking at 40 min then collapsing to 0 by 80 min** | ★ **the sharpest turnover in the paper** |
| 3c | **MGO** | 0–0.10 mmol/L | ≈ 0.072 | ≈ 0.082 | 40–70 min |
| 3d | **GO** | 0–0.010 mmol/L | ≈ 0.0057 | ≈ 0.0098 | 40–60 min |

★ **The α-dicarbonyl hierarchy, same panel, same lab, absolute units:**
**3-DX (0.82–1.13) ≫ 1-DX (0.11–0.18) > MGO (0.072–0.082) ≫ GO (0.0057–0.0098) ≫ DA (0.003–0.004)
mmol/L at peak** — a **~200-fold span from 3-DX to DA**, and it is **the same ordering and nearly
the same magnitudes as the JAFC companion's Fig. 5** (3-DX 0.81, 1-DX 0.11, MGO 0.072, GO 0.0057),
which is a genuine internal replication across two papers and two sugar arms. `USE-Q`.

★ **The DA story is mechanistically load-bearing** and the authors state it explicitly (p. 6):
> *"the temperature-dependent DA content curves (Fig. 3b) proved that an elevated temperature could
> significantly promote the generation of DA, however, **a sharp consumption of DA appeared when it
> rapidly reached to its maximum.** Meanwhile, 3-mercapto-2-butanone derived from a reaction
> between DA and H₂S correspondingly exhibited an elevation trend when the temperature increased
> from 100 to 140 °C."*

**⇒ A named, measured precursor→product pair with opposite temperature responses: DA is consumed
faster than it forms above ~120 °C, and its H₂S adduct (3-mercapto-2-butanone) rises to compensate.
That is a directly testable structural constraint on the α-dicarbonyl → mercaptoketone edge.**

---

## §8. ★★★ THE PROVENANCE CHAIN — the corpus's sulfur ladder has ONE source, published FOUR times

| # | publication | date | where the ladder appears | what it adds |
|---|---|---|---|---|
| **1** | ★ **Zhai et al., *Food Chem.* 404:134420 — THIS PAPER** | **Sep 2022** | **Table 1, TTCA columns** — 34 compounds × 3 temperatures, **with SDs and Duncan letters** | ★ **the primary source.** Also the ONLY place the +Xyl arm, the SDs, the ¹³C tracing and the pH-4.9 drift appear |
| 2 | Zhai et al., *JAFC* 71:14300 | Sep 2023 | **Fig. 3, column E** (T-120 min) | 5 extra compounds, 3 extra hold times, the +Cys arm; **inherits the 2-methylthiophene slip**; SDs dropped |
| 3 | Kang et al., *Sust. Food Technol.* 4:3239, **Table S4** | Mar 2026 | the ladder as an SI table | ★ **corrects the 2-methylthiophene slip to `0.000`**; SDs re-attached (see below) |
| 4 | Kang et al., **Table S5, pH-7 column** | Mar 2026 | the 100 °C column again, as a pH point | (already identified as a duplicate in `kang2026_SI_extraction.md` §4c) |

**Verification:** the 34-row × 3-temperature block agrees **cell for cell** across sources 1, 2 and 3
(101 of 102 cells; the single exception is the 2-methylthiophene 120 °C cell, diagnosed in §4.1),
and all six class subtotals plus all three grand totals agree exactly, **including the arithmetic
error at 120 °C, which is reproduced identically in sources 1 and 3.**

### 8.1 ★ AND THE SDs TRAVEL TOO — which settles a question the K5 wave left open
`kang2026_SI_extraction.md` §4c found that **"Table S5's pH-7 column IS Table S4's 100 °C column"**
and, separately, that the printed SDs *"average ≈ 5 % relative — implausibly tight for HS-SPME at
sub-μg/L levels, and demonstrably not reproducible between S4 and S5."*
**The SDs printed in Kang's Table S4 are the SDs printed in this paper's Table 1** — e.g.
MFT `1.237 ± 0.142` / `1.388 ± 0.256` / `5.907 ± 0.085`; FFT `3.734 ± 0.085` / `4.107 ± 0.137` /
`11.439 ± 0.265`; furfural `3.381 ± 0.089` / `5.793 ± 0.422` / `11.039 ± 0.302`; sulfur subtotal
`13.978 ± 0.568` / `35.866 ± 1.382` / `60.4 ± 1.62`. **Identical to the last digit.**
**⇒ The SDs are not fabricated and not re-derived — they are this paper's original triplicate SDs,
carried through. The K5 objection should be restated: they are ORIGINAL, but they are tight for
the method, and they are internally inconsistent with Table S5 because Table S5's column is a
duplicate that was re-typeset, not re-measured.**

### 8.2 ★★ THE SEVEN CONSEQUENCES FOR THE REPO
1. **The ladder is n = 1.** Four publications, one experiment. `research_round4_nulls.md` §C.2's
   null ("no second sulfur temperature ladder") is not merely unfalsified — **the appearance of a
   second source was an artefact of the duplication.**
2. ★ **The declared provenance of every `kang_*` row should be re-pointed to
   `zhai2023foodchem_table1`**, with Kang recorded as a re-publication.
3. ★ **Any weighting that counts Kang S4 and Kang S5 as two observations is wrong by 2×**; adding
   either Zhai paper without noticing makes it 3× or 4×.
4. ★ **The Tier A / Tier B uncertainty split for the ladder must be revised to Tier B**, because
   the primary source used `f′ = 1` (§2.3).
5. ★ **The "measured pH 4.9" comparator used in `kinetic_core_b2_2_diagnosis.md` §3 originates
   here**, unattached to a temperature or a time (§3). Its use as a diagnostic is fine; it must
   never become a scored row.
6. **The 2-methylthiophene 120 °C cell is `not detected`**, not 1.128 (§4.1).
7. ★ **The +Xyl arm, the ¹³C table, the 100 °C 41 × 15 grid, the MRP control, the Cys-alone
   ladder and all four non-volatile time courses are NEW to the corpus** — they have never been
   republished and never been ingested. **That is the material this wave actually adds.**

---

## §9. ★ DERIVED ACTIVATION ENERGIES `[D]`

`Eₐ = R·ln(y₂/y₁)/(1/T₁ − 1/T₂)`, kJ mol⁻¹, from Table 1 (120 min hold), **both arms**.

### 9.1 Class level
| class | **TTCA arm** low / high / overall | **TTCA+Xyl arm** low / high / overall |
|---|---|---|
| **Sulfur (25)** | 59.4 / 33.1 / **46.9** | ★ **92.4 / 38.9 / 67.0** |
| **Nitrogen (4)** | — (zero at every T) | — / 64.4 / — |
| **Oxygen (5)** | 36.6 / 49.2 / **42.6** | 56.0 / 76.7 / **65.9** |
| **TOTAL** | 55.5 / 35.6 / **46.1** | ★ **89.9 / 44.8 / 68.5** |

★★ **Adding an equimolar sugar raises the apparent activation energy of the whole system by
~22 kJ mol⁻¹ (46.1 → 68.5).** Both arms **decelerate** on the top leg, which is the normal
precursor-depletion signature. **The +Xyl arm's low-leg Eₐ of ~90 kJ mol⁻¹ is the largest
class-level value in the K6a corpus outside the extrusion lane.**

### 9.2 Compound level — the compounds that matter
| compound | **TTCA** low / high / overall | **TTCA+Xyl** low / high / overall | note |
|---|---|---|---|
| **MFT** | ★ **7.0 / 97.8 / 50.1** | **22.4 / 17.5 / 20.1** | ★★ **the "switch-on" appears in the TTCA arm and VANISHES when sugar is added** |
| **FFT** | ★ **5.8 / 69.2 / 35.9** | **18.3 / −3.0 / 8.2** | ★★ same — and the +Xyl high leg goes **negative** |
| 2-Acetylthiazole | 64.0 / 7.7 / 37.3 | 107.4 / 15.6 / 63.8 | decelerating in both arms |
| Furfural | 32.8 / 43.5 / 37.9 | 44.0 / 75.5 / 59.0 | the only cleanly Arrhenius species in either arm |
| 3-Mercapto-2-butanone | 19.6 / 19.4 / **19.5** | 49.5 / −6.4 / 23.0 | ★ **the TTCA arm is almost perfectly Arrhenius (two legs within 1 %)** |
| 3-Mercapto-2-pentanone | 17.6 / 19.6 / **18.6** | 43.0 / −5.6 / 19.9 | ★ same |
| Thiazole | 42.4 / 8.5 / 26.3 | 50.9 / 30.7 / 41.3 | |
| 2-Thiophenecarboxaldehyde | 194.8 / 2.5 / 103.6 | 230.0 / 58.7 / 148.7 | extreme low leg, zero high leg |
| 5-Methyl-2-thiophenecarboxaldehyde | 83.0 / 30.2 / 57.9 | 91.0 / 76.7 / **84.2** | |
| 3-Thiophenethiol | — / 27.6 / — | — / 29.9 / — | zero at 100 °C in both arms |
| Thieno[3,2-b]thiophene | — / 64.2 / — | — / 98.2 / — | " |
| 2-Methylthieno[2,3-b]thiophene | — / 163.0 / — | — / 66.0 / — | " |

### 9.3 ★★ THE KANG SWITCH-ON CROSS-CHECK — **it does not survive the co-substrate axis either**

The `kang_switch_on_*` feature (MFT ×1.12 then ×4.26; FFT ×1.10 then ×2.79; Eₐ climbing from
~6 to 70–98 kJ mol⁻¹) is, in this primary source, **a property of the TTCA-alone arm at 120 min**.
In the **+Xyl arm of the same experiment, at the same three temperatures and the same hold time**:

| | TTCA arm ×(100→120) | ×(120→140) | **TTCA+Xyl** ×(100→120) | ×(120→140) |
|---|---:|---:|---:|---:|
| **MFT** | 1.12× | **4.26×** | **1.44×** | **1.30×** |
| **FFT** | 1.10× | **2.79×** | **1.35×** | **0.96×** ★ *falls* |

**⇒ Both thiols become nearly flat, and FFT actually DECREASES from 120 to 140 °C, once an
equimolar co-substrate is present.** Combined with `zhai2023jafc_extraction.md` §8.3 (the feature
is absent at 40, 80 and 180 min in the same system), the switch-on now fails on **two independent
axes — hold time and co-substrate — within its own source laboratory.**

★ **A THIRD, INDEPENDENT PIECE OF EVIDENCE POINTS THE OTHER WAY.** The ¹³C table (§6.1) DOES show a
genuine temperature threshold, but at **100→120 °C**, on **thiophene-ring products**, measured as
an **isotope ratio**. **The corpus should replace "MFT/FFT switch on between 120 and 140 °C" with
"exogenous-sugar carbon enters the retro-aldol/Strecker branch between 100 and 120 °C" — the second
statement is supported by a response-factor-free measurement and the first is not.**

### 9.4 ★ THE CYSTEINE-DEGRADATION LADDER — and the origin of the repo's 55.1 kJ mol⁻¹
`[M]` p. 6, verbatim: *"The Cys alone model (Fig. S2) revealed that its consumption increased from
**16.2 %** at 100 °C to **38.8 %** and **62.2 %** at 120 and 140 °C, respectively."*
(Cys 10 mmol/L, pH 7, **120 min**.)

`[D]` As pseudo-first order, `k = −ln(1 − x)/120 min`:

| T | conversion | **k (min⁻¹)** |
|---:|---:|---:|
| 100 °C | 16.2 % | **1.473 × 10⁻³** |
| 120 °C | 38.8 % | **4.092 × 10⁻³** |
| 140 °C | 62.2 % | **8.107 × 10⁻³** |

**Eₐ: 62.3 (low leg) / 46.2 (high leg) / ★ 54.7 kJ mol⁻¹ overall, R² = 0.996.**

★★ **This is where the repo's `MEASURED_EA_OVERRIDES` value of 55.1 kJ mol⁻¹ for cysteine depletion
comes from.** `kang2026_SI_extraction.md` §6b derived it by digitising Kang's Fig. S4 (62.6 % and
38.7 % conversions) — **the same experiment, one publication downstream, read off a figure instead
of a sentence.** The agreement (54.7 vs 55.1, 0.7 %) is a validation of that digitisation, **not
an independent confirmation.** `kinetic_core_b2_2_diagnosis.md`'s claim K-4, that the 140 °C
cysteine row is *"governed by Kang's own immovable 55.1 kJ/mol"*, should read **"Zhai's"**.

⚠️ **And note the same non-Arrhenius shape**: 62.3 then 46.2 — the cysteine sink **decelerates**
across the very interval where the thiols were said to accelerate. **A single 55.1 kJ mol⁻¹ is a
40 K average of a mildly curved line, and the curvature has the opposite sign to the switch-on.**

---

## §10. VERIFIED NEGATIVES `[NEG]`

| # | sought | verdict |
|---|---|---|
| 1 | Any rate constant, half-life, Eₐ, or kinetic model printed by the authors | **ABSENT.** No kinetics anywhere; every number in §9 is `[D]` |
| 2 | H₂S measured directly | **ABSENT** — inferred throughout from Cys consumption |
| 3 | TTCA or Xyl depletion curves | **ABSENT from this paper** — the methods describe TTCA/Cys/Xyl quantification and then *"the specific parameters were not shown"*; only **regenerative Cys** is plotted (Fig. 3a). TTCA depletion is deferred to Zhai et al. 2021 |
| 4 | Buffer | **ABSENT** — unbuffered, NaOH-adjusted |
| 5 | pH as a function of time or temperature | **ABSENT** — a single unattached "7.0 → 4.9" |
| 6 | LOD / LOQ for the volatiles | **ABSENT** (they exist for TTCA/Cys/Xyl but are deferred to Zhai 2021) |
| 7 | Response factors / calibration curves for the volatiles | ★ **ABSENT BY CONSTRUCTION — `f′ = 1`** (§2) |
| 8 | Any temperature outside 100/120/140 °C | **ABSENT** |
| 9 | Any hold time beyond 140 min | **ABSENT** (the JAFC companion goes to 180 min) |
| 10 | Sensory panel | **ABSENT** — ROAV only, computed from literature thresholds (Table S1, not on disk) |
| 11 | MFT/FFT in the ¹³C table | ★ **ABSENT — the two flavour-critical thiols were NOT isotope-traced**, which is exactly the gap that would have settled the switch-on |
| 12 | Melanoidin / browning measurement | **ABSENT** — no A₄₂₀ in this paper at all (the JAFC companion has it) |

---

## §11. CONSOLIDATED PARAMETER TABLE

**Common conditions:** TTCA 10 mmol/L (± Xyl 10 mmol/L), water, **unbuffered**, pH 7.0 set with
NaOH, pressure bottle with stirring, ice quench, HS-SPME (CAR/PDMS/DVB 75 μm, 60 °C/20 min, 3 g
sample **+ 2 g sat. NaCl**) GC-MS on DB-WAX, **IS = 1,2-dichlorobenzene with `f′ = 1`**, n = 3,
**μg/L (semi-quant)**.

| # | parameter | value | units | condition | prov | anchor |
|---:|---|---|---|---|---|---|
| 1 | **Table 1, complete** | 34 compounds × 6 columns, with SDs | μg/L | 100/120/140 °C × ±Xyl, 120 min | **M** | §4 |
| 2 | **Fig. 1a grid, complete** | 41 compounds × 15 columns | μg/L | 100 °C, 7 times × ±Xyl + MRP | **M** | §5 |
| 3 | **Table 2, complete** | 6 compounds × 3 T × 6 isotopomers | % | ¹³C₅-Xyl, 120 min | **M** | §6 |
| 4 | ⚠ **2-methylthiophene, TTCA, 120 °C** | ★ **not detected** (the printed 1.128 is a copy-paste) | — | | **D** | §4.1 |
| 5 | **MFT ladder, TTCA** | 1.237 / 1.388 / 5.907 | μg/L | 100/120/140 °C | **M** | Table 1 |
| 6 | **MFT ladder, TTCA+Xyl** | ★ **1.729 / 2.498 / 3.238** | μg/L | " | **M** | ★ new |
| 7 | **FFT ladder, TTCA** | 3.734 / 4.107 / 11.439 | μg/L | " | **M** | Table 1 |
| 8 | **FFT ladder, TTCA+Xyl** | ★ **4.736 / 6.398 / 6.123** | μg/L | " | **M** | ★ new |
| 9 | **Sulfur-class Eₐ, TTCA** | 59.4 / 33.1 / **46.9** | kJ mol⁻¹ | low/high/overall | **D** | §9.1 |
| 10 | **Sulfur-class Eₐ, TTCA+Xyl** | ★ **92.4 / 38.9 / 67.0** | kJ mol⁻¹ | " | **D** | §9.1 |
| 11 | **Total-volatile Eₐ, TTCA / +Xyl** | **46.1** / **68.5** | kJ mol⁻¹ | 100→140 °C | **D** | §9.1 |
| 12 | ★ **Cys-alone conversion** | **16.2 / 38.8 / 62.2** | % | 100/120/140 °C, 10 mM, pH 7, 120 min | **M** | p. 6 |
| 13 | ★ **Cys-alone k** | 1.473 / 4.092 / 8.107 × 10⁻³ | min⁻¹ | " (1st-order assumption) | **D** | §9.4 |
| 14 | ★★ **Cys-degradation Eₐ** | ★ **54.7** (legs 62.3 / 46.2), R² 0.996 | kJ mol⁻¹ | 100–140 °C | **D** | §9.4 — ★ **the true origin of the repo's 55.1** |
| 15 | **Exogenous-carbon incorporation, 3-thiophenethiol** | **0 / 19 / 34** | % | 100/120/140 °C | **M** | Table 2 |
| 16 | **… thieno[3,2-b]thiophene** | **0 / 21 / 45** | % | " | **M** | Table 2 |
| 17 | **… 2-methylthieno[2,3-b]thiophene** | **<1 / 21 / 45** | % | " | **M** | Table 2 |
| 18 | **… 2-thiophenecarboxaldehyde** | 24 / 41 / 52 | % | " | **D** from Table 2 | §6.1 |
| 19 | ★ **3-DX peak** | **≈ 0.82** (TTCA) / **≈ 1.13** (+Xyl) | mmol/L | 80–100 / 60 min | **D** | §7 |
| 20 | 1-DX peak | ≈ 0.11 / ≈ 0.175 | mmol/L | | **D** | §7 |
| 21 | MGO peak | ≈ 0.072 / ≈ 0.082 | mmol/L | | **D** | §7 |
| 22 | GO peak | ≈ 0.0057 / ≈ 0.0098 | mmol/L | | **D** | §7 |
| 23 | ★ **DA peak** | ≈ 0.0029 / **≈ 0.0040, collapsing to 0 by 80 min** | mmol/L | | **D** | §7 |
| 24 | **regenerative Cys plateau** | ≈ 1.0 of 10 loaded, i.e. **≈ 10 % of TTCA's N released as free Cys** | mmol/L | | **D** | §7 |
| 25 | ★ **pH drift** | **7.0 → 4.9** | — | unattached to T or t | **M** | p. 5 |
| 26 | **TTCA vs fresh MRP intensity** | TTCA total = **1.52×** MRP (paper says 1.56×) | — | 100 °C/120 min vs 120 °C/120 min ⚠ | **D** | §5.1 |
| 27 | **species count** | TTCA **19**, MRP **40**, TTCA+Xyl at 140 °C ~34 | — | | **M** | p. 4 |
| 28 | ★ **ROAV, 140 °C TTCA+Xyl** | 3-thiophenethiol **195**; 2-acetylthiazole **100**; 2-thiophenecarboxaldehyde **0.04** | — | Table S1 (not on disk) | **M** (text) | p. 6 |

---

## §12. USABILITY VERDICTS

### `USE`
- **Table 2 (isotope tracing)** — a dimensionless ratio, immune to the `f′ = 1` problem, with SDs
  absent but the quantity robust. ★ **the strongest artefact in this paper.**
- **Cys-alone conversions 16.2 / 38.8 / 62.2 %** — absolute, external-standard HPLC.
- **The α-dicarbonyl peak values** (absolute quantification), as magnitudes with ±15 %.

### `USE-Q`
| item | qualification |
|---|---|
| **Cys-degradation Eₐ 54.7 kJ mol⁻¹** | first-order assumed from a single conversion per temperature; ⚠ **not independent of the repo's existing 55.1** |
| **All volatile Eₐ (§9)** | semi-quant basis cancels in the slope, but they are **yield-response coefficients at a 120-min slice**, not barriers; quote both legs |
| **Sulfur-class Eₐ 46.9 (TTCA) / 67.0 (+Xyl)** | the class is a sum over unequal response factors — an index, not a mass |
| **α-dicarbonyl hierarchy** | replicated across both this paper and the JAFC companion |

### `RATIO-ONLY`
- Every **+Xyl / TTCA** ratio, and every **time / time** ratio within a column.
- **TTCA vs fresh-MRP** comparisons ⚠ **at different temperatures** (§5.1).

### `STRUCTURAL`
1. ★ **N-heterocycles are exactly zero in a TTCA-only pot at 100, 120 and 140 °C, and appear as
   soon as free xylose is present.** Pyrazine formation is dicarbonyl-limited, not amine-limited.
2. ★ **Exogenous sugar carbon enters the thiophene-ring branch only above 100 °C** (§6.1).
3. ★ **The TTCA arm's volatile total peaks at 120 min at 100 °C; the +Xyl arm does not peak at
   all within 140 min.** Extra substrate postpones turnover.
4. **Diacetyl peaks early and is consumed to zero while its H₂S adduct rises.**
5. **3-DX ≫ 1-DX > MGO ≫ GO ≫ DA at peak, over ~200×.**
6. **The intermediate out-yields the fresh Maillard reaction on intensity (1.5×) and loses half
   the species** — the paper's thesis, and a useful sanity constraint on any MRI-vs-MRP model.

### `PRIOR-ONLY`
- The ROAV figures (195 / 100 / 0.04) — computed from literature thresholds in a table not on disk.

### ★ `REFUSE`
| item | reason |
|---|---|
| **Absolute μg/L volatile magnitudes** | `f′ = 1` semi-quant (§2.2) |
| **Cross-compound volatile comparisons** | different unknown response factors |
| **`1.128` for 2-methylthiophene, TTCA, 120 °C** | copy-paste artefact (§4.1) |
| **`42.024` and `35.866` as the 120 °C TTCA total/subtotal** | ⚠ they are the *correct* sums **given the corrected cell**, so keep them — but never alongside the 1.128 cell |
| **`10.020 μg/L`** (Conclusion's pyrazine total) | contradicts Table 1's 10.296 |
| **Any Kang-attributed use of this ladder as independent** | §8 |
| **Cross-arm pyrazine comparisons at face value** | pyrazine samples were extracted at **pH 8**, the others at pH 7 (§2.4) |

---

## §13. DRAFT FIT / HOLD-OUT ROLES — **FOR THE ORCHESTRATOR, NOT A DECLARATION EDIT**

*(No declaration file was opened or edited.)*

### 13a. Corrections to existing rows
- ★ Re-point the provenance of **every `kang_*` sulfur row and the 55.1 kJ mol⁻¹ cysteine
  override** to `zhai2023foodchem_table1` / `zhai2023foodchem_p6`, recording Kang as a
  re-publication (§8.2).
- ★ **Retire `kang_switch_on_*`** — the feature fails on hold time (JAFC companion) and on
  co-substrate (§9.3) inside its own source lab.

### 13b. Recommended **FIT**
| candidate | why |
|---|---|
| **Cys-alone conversion ladder (16.2 / 38.8 / 62.2 %)** | absolute, external-standard, three rungs, already implicitly in the repo — but should be attributed correctly |
| **α-dicarbonyl peak magnitudes and peak times, TTCA arm** | absolute, and the trunk needs them |

### 13c. Recommended **HOLD-OUT** — all genuinely unseen
| candidate | what it tests |
|---|---|
| ★ **The +Xyl arm at all three temperatures** (Table 1, 34 compounds) | a co-substrate perturbation the model has never seen; sulfur class must go 20.19 → 91.81 → 163.43 |
| ★ **N-heterocycles = 0 in TTCA-only, non-zero with Xyl** | a hard on/off structural test |
| ★ **The ¹³C exogenous-carbon fractions (0 / 19–21 / 34–45 %)** | ★ **the best hold-out in the wave** — response-factor-free, and it tests carbon routing rather than yield |
| **The 100 °C 41 × 15 time grid** | the time axis at the temperature the Yiltirak family fails on |
| **The MRP control column** | tests whether the model can distinguish an intermediate pot from a fresh Maillard pot |

### 13d. **DO NOT USE**
- Table 1's TTCA columns as new anchors (already in the repo via Kang).
- Absolute volatile μg/L for any calibration.

---

## §14. DECLARED GAPS

| # | gap | closure |
|---|---|---|
| G1 | **Table S2 (MGO−GO−Cys validation model)** is the only quantitative link from the dicarbonyl pool to the pyrazine channel here, and it is not on disk | order the SI of `10.1016/j.foodchem.2022.134420` |
| G2 | **Fig. S2 (Cys-alone degradation curves)** — only three summary percentages are in the text; the curves would give shape, not just endpoint | same SI |
| G3 | **Table S1 (ROAV / thresholds)** | same SI |
| G4 | **MFT and FFT were not isotope-traced**, which is the measurement that would settle the switch-on | would require new experiments |
| G5 | **TTCA depletion is deferred to Zhai et al. 2021** (*JAFC* 67:8632 / the 2021 degradation paper), which is **not on disk** | order `Zhai et al. 2021, Food Chem./JAFC` — ★ **it is also the source of the LODs and the HPLC methods this paper omits** |
| G6 | **No SDs on Table 2 or on any figure series** | — |
| G7 | ★ **Whether Kang 2026 cites this paper as the source of Table S4** — the repo's independence bookkeeping depends on it | read Kang's reference list |
