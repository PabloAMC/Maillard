# Kocadağlı & Gökmen 2016 (Food Chemistry) — glucose / wheat flour, low moisture

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

---

## §0. IDENTITY — ⚠️ AND A TEXT-LAYER DEFECT THAT MUST BE RECORDED

| item | value | how verified |
|---|---|---|
| **file on disk** | `data/articles/Kocadagli2016.pdf` (613,181 bytes) | `ls` |
| **DOI** | **`10.1016/j.foodchem.2016.05.150`** | Elsevier cover page of the PDF, line 9 |
| PII | `S0308-8146(16)30839-1` · Reference `FOCH 19299` | cover page |
| title | *"Multiresponse Kinetic Modelling of Maillard Reaction and Caramelisation in a Heated Glucose/Wheat Flour System"* | cover page |
| authors | Tolgahan Kocadağlı, Vural Gökmen | cover page |
| journal | *Food Chemistry* **211**, 892–902 (2016) — the volume/pages are **not** in this PDF; they are taken from the citation the JAFC companion paper gives (`Kocadağlı, T., Gökmen, V., Food Chem. 2016, 211, 892–902`, ref. 24 of `Kocada2016.pdf`) | cross-document |
| dates | Received 26 Feb 2016 · Revised 19 May 2016 · Accepted 23 May 2016 | cover page |
| version | **Elsevier "Accepted Manuscript"** — *"a PDF file of an unedited manuscript … The manuscript will undergo copyediting, typesetting, and review of the resulting proof before it is published in its final form."* | cover page, verbatim |

**⚠️ WRONG-FILE WARNING — see `kocadagli2016jafc_extraction.md` §0.** The two Kocadağlı 2016 files
are mapped the *opposite* way round from the naïve guess: `Kocadagli2016.pdf` (longer stem) is
**this** Food Chemistry wheat-flour paper, and `Kocada2016.pdf` (shorter stem) is the JAFC NaCl
caramelization paper.

### 0.1 ★ THE PDF TEXT LAYER OF THIS FILE IS CIPHER-GARBLED. IT WAS DECODED BY THIS WAVE.

`pdftotext` on this file returns mojibake for **876 of its 2 035 lines** — every line of the
manuscript body (the ASCII Elsevier cover page is unaffected). Example, raw:

```
DƵůƚŝƌĞƐƉŽŶƐĞ <ŝŶĞƚŝĐ DŽĚĞůůŝŶŐ ŽĨ DĂŝůůĂƌĚ ZĞĂĐƚŝŽŶ
```

The embedded subset font carries no usable `ToUnicode` CMap, so the extractor falls back to glyph
indices. The map is a **strict order-preserving per-glyph substitution**, which I solved by
aligning the garbled title against the plaintext title on the ASCII cover page and then closing the
remaining letters from context (`Kocadağlı`, `van Boekel`, `Yaylayan`, `Belitz`, `Prof. Dr.`).

**The full recovered cipher, for any future wave that has to re-read this file:**

| class | mapping |
|---|---|
| digits 0–9 | `U+03EC + d` (i.e. `Ϭ`=0 … `ϵ`=9) |
| lowercase a–z | `0x102 0x10F 0x110 0x11A 0x11E 0x128 0x150 0x15A 0x15D 0x169 0x16C 0x16F 0x175 0x176 0x17D 0x189 0x18B 0x18C 0x190 0x19A 0x1B5 0x1C0 0x1C1 0x1C6 0x1C7 0x1CC` |
| uppercase (observed) | A`0x04` B`0x11` C`0x12` D`0x18` E`0x1C` F`0x26` G`0x27` H`0x2C` I`0x2F` K`0x3C` L`0x3E` M`0x44` N`0x45` O`0x4B` P`0x57` Q`0x59` R`0x5A` S`0x5E` T`0x64` U`0x68` V`0x73` W`0x74` Y`0x7A` |
| punctuation | `.`=`0x358` `,`=`0x355` `-`=`0x372` `(`=`0x37E` `)`=`0x37F` `&`=`0x398` `:`=`0x357` `;`=`0x356` `–`=`0x374` `[`=`0x123E` `]`=`0x123F` space=`0x03` `/`=`0x02CD` |
| symbols | `°C`=`0x3A3` (single glyph) `α`=`0x272` `%`=`0x439` `→`=`0x13A` `µ`=`0x285` `±`=`0x446` `+`=`0x43D` `×`=`0x43F` `ö`=`0x182` `ğ`=`0x152` `ı`=`0x166` `ü`=`0x1BA` `*`=`0x38E` |
| Appendix-A math font | separate encoding, additionally **right-to-left reversed within each token** (`[݈ܿܩ]` reads `Glc`). `k`=`0x747` `d`=`0x740` `t`=`0x750` `A`=`0x723` `G`=`0x729` `D`=`0x726` `=`=`0xD4C` `−`=`0xD46` `+`=`0xD45`; digit subscripts in the `0x0B34–0x0B3B` range |

Decoded working copy written to the wave scratchpad only. **Every quotation below was read from the
decoded text and spot-verified against the raw codepoints.** ⚠️ Any earlier grep of this file for
`HMF`, `hydroxymethyl`, `Ea`, or any English word **returned zero hits regardless of content** —
so any prior "verified negative" on this file that was established by grep is **void**.
(Note: the corpus-wide HMF sweep recorded in `research_round3_channels.md` §A.1 lists six files
with >3 `HMF` hits and this file is not among them — **because its text layer is unreadable, not
because it lacks HMF.** It contains an HMF node and eleven HMF-bearing rate constants.)

---

## §1. SYSTEM AND CONDITIONS — verbatim, `[M]`

| item | value | anchor |
|---|---|---|
| matrix | **10 g all-purpose wheat flour + 1 g glucose + 6 mL deionised water → dough**, rolled thin, dried **50 °C for 2 h** (air circulation), finely ground | §2.2.1 |
| moisture | **7 % humidity**, by AACC Method 44-15.02 — *"a homogenous low moisture system (7 % humidity) mimicking the crust of bakery products"* | §2.2.1 |
| **buffer** | **NONE, deliberately** — *"This food system eliminated the necessity of a buffer, which would affect the fate of the reaction."* | §1, verbatim |
| **pH** | ⚠️ **NOT MEASURED AND NOT REPORTED ANYWHERE IN THIS PAPER.** | verified negative |
| sample | **100 mg** ground mixture per glass tube, PTFE-sealed screw caps | §2.2.2 |
| temperatures | **160, 180 and 200 °C** | §2.2.2 |
| times | **up to 20 min** | §2.2.2 |
| replication | **four repetitions** per temperature–time (2 extracted for solutes, 2 for protein hydrolysis) | §2.2.2–2.2.3 |
| control | a **glucose-only** freeze-dried system (0.5 mL solution containing 9 mg glucose, frozen −80 °C 2 h, freeze-dried 6 h at 0.1 Pa, ice condenser −76 °C), **sugar analysis only** | §2.2.2 |
| units of the model | **µmol/g** | throughout §3 |
| initial glucose | **432.9 ± 11.6 µmol/g** | §3.1 |
| initial fructose | **3.4 ± 0.3 µmol/g** | §3.1 |
| initial total free amino acids | **7.7 ± 0.1 µmol/g** | §3.1 |
| initial Lys residues | **14.2 ± 0.3 µmol/g** | §3.1 |
| initial Arg residues | **15.8 ± 0.3 µmol/g** | §3.1 |
| peak fructose | **54.7 ± 4.5** (160 °C × 20 min), **73.5 ± 0.7** (180 °C × 10 min), **73.4 ± 1.7 µmol/g** (200 °C × 5 min) | §3.1 |

⚠️ **"AA" IS A LUMP.** Verbatim, §3.1: *"The results were expressed as 'total amino acids (AA)'
comprising the mole quantities of **all individual free amino acids and Lys and Arg residues** for
modelling purposes."* So `[AA]₀ = 7.7 + 14.2 + 15.8 = 37.7 µmol/g`, mixing free amino acids with
**protein-bound** lysine and arginine side chains at a single reactivity. The authors name this as
a defect: *"The peak concentrations of 1-deoxyglucosone and 3-deoxyglucosone were underestimated by
the proposed model at 180 °C and 200 °C. This could be due to the **different reaction rates of free
amino acids, lysine and arginine residues**."*

### 1.1 Analytical quality `[M]`

Same instruments and derivatisation as the JAFC companion (o-phenylenediamine → quinoxalines,
Agilent 1200 / 6130 SQ MS, SIM). **The same semi-quantitation caveat applies and the authors state
it explicitly**, §3.4 verbatim:

> "It should be also considered that **1-deoxyglucosone, glucosone and 3,4-dideoxyglucosone were
> semi-quantitated** and this can be a source of error for the multi-response model."

⇒ **`absolute_concentration: false` on 1-DG, G and 3,4-DG.** 3,4-DG is the immediate precursor of
HMF on the 3-DG limb, so **k₁₂ and k₁₃ both inherit an unknown response-factor scale.**

Also: **Amadori products were never measured** — *"Amadori products (APs) corresponding to each free
amino acid and also glycated lysine and arginine residues were **not quantitated, due to the lack of
commercially available standard compounds**"* (§3.3). AP is therefore a second unmeasured node
alongside `Int` (§2).

---

## §2. THE REACTION NETWORK — 26 elementary steps (Figure 2), and the HMF ODE

Abbreviations verbatim (Table 1 caption): `Glc` glucose; `Fru` fructose; `1,2-ED` 1,2-enediol;
`AA` total amino acids; `AP` Amadori product; `1-DG`; `3-DG`; `3,4-DG` 3,4-dideoxyglucosone;
`G` glucosone; `GO` glyoxal; `MG` methylglyoxal; `DA` diacetyl; `HMF`; `Int` undetermined
intermediate; `P` products.

**Full edge list, from Table 1 (steps marked `*` are bimolecular):**

```
 1  Glc → 1,2-ED             14  Fru → Int
 2  1,2-ED → Glc             15  Int → HMF          ★ terminal, FRUCTOSE LIMB
 3  1,2-ED → Fru             16  G → GO
 4  Fru → 1,2-ED             17  1-DG → MGO
 5  Glc → G                  18  1-DG → DA
 6  Glc → 3-DG               19  HMF → P8           ★ THE HMF SINK (unimolecular)
 7* Glc + AA → AP            20* 3,4-DG + AA → P7
 8  AP → 1-DG                21  3-DG → P6
 9  AP → G                   22* GO + AA → P5
10  AP → GO                  23  G → P4
11  AP → 3-DG                24* DA + AA → P3
12  3-DG → 3,4-DG            25* MG + AA → P2
13  3,4-DG → HMF   ★ terminal, 3-DG LIMB      26* 1-DG + AA → P1
```

Two further steps are **named in the text and excluded by model discrimination**:
`k₂₇` Fru → 1-DG (*"estimated as zero or very low … Therefore, it was excluded"*, §3.4) and
`k₂₈` 3-DG → MGO (*"was zero at all temperatures studied"*, §3.5).

**★ THE HMF NODE ODE (Appendix A, decoded; the `Ĩ` glyph and the RTL math font are resolved):**

> **`d[HMF]/dt = k₁₃·[3,4-DG] + k₁₅·[Int] − k₁₉·[HMF]`**

**Structurally identical to the JAFC companion**: two parallel first-order sources, one
first-order sink. Note what is *absent*: **the HMF + AA bimolecular sink was TESTED AND REJECTED
here** — see S4 below.

Other decoded ODEs of interest:
```
d[Glc]/dt   = −k7[Glc][AA] − (k1+k5+k6)[Glc] + k2[1,2-ED]
d[Fru]/dt   =  k3[1,2-ED] − (k4+k14)[Fru]
d[AA]/dt    = −k7[Glc][AA] + (k8+k9+k10+k11)[AP] − k26[1-DG][AA] − k24[DA][AA]
               − k22[GO][AA] − k25[MGO][AA] − k20[3,4-DG][AA]
d[3-DG]/dt  =  k6[Glc] + k11[AP] − (k12+k21)[3-DG]
d[3,4-DG]/dt=  k12[3-DG] − k13[3,4-DG] − k20[3,4-DG][AA]
d[Int]/dt   =  k14[Fru] − k15[Int]
d[AP]/dt    =  k7[Glc][AA] − (k8+k9+k10+k11)[AP]
d[1,2-ED]/dt=  k1[Glc] + k4[Fru] − (k2+k3)[1,2-ED]
```

⚠️ **`Int` AND `AP` ARE BOTH UNMEASURED.** As in the JAFC paper, `[Int]` is never quantified, so
**k₁₄ and k₁₅ are identified only up to a common scale on the `Int` pool.** `AP` is likewise
unmeasured, so k₇–k₁₁ share a second unidentified scale. **Five of the six numbers on the HMF
fructose limb and the whole Amadori block are scale-free-only quantities.**

---

## §3. TABLE 1 — the complete 26-step × 3-temperature constant set

**Transcribed verbatim from the decoded Table 1 (PDF p.51 of the accepted manuscript).**
Caption verbatim: *"Estimated reaction rate constants (k, min⁻¹×10³ or (\*) g mixture·µmol⁻¹·min⁻¹×10³)
with 95 % highest posterior density (HPD) intervals."*
Footnote: **`ᵃ Not applicable`** — the cells rendering as a bare `0` followed by the footnote marker.
Provenance: **`[F]`** throughout.

| # | elementary step | k(160) | HPD | k(180) | HPD | k(200) | HPD |
|---|---|---|---|---|---|---|---|
| 1 | Glc → 1,2-ED | 48.3 | ±26.1 | 81.8 | ±31.1 | 204 | ±43.2 |
| 2 | 1,2-ED → Glc | 618 | ±445 | 219 | ±302 | 255 | ±105 |
| 3 | 1,2-ED → Fru | 170 | ±54 | 511 | ±367 | 946 | ±294 |
| 4 | Fru → 1,2-ED | **0** | ᵃ | 168 | ±243 | 289 | ±107 |
| 5 | Glc → G | 0.038 | ±0.01 | 1.69 | ±0.44 | 1.38 | ±0.24 |
| 6 | Glc → 3-DG | 4.50 | ±0.9 | 12.2 | ±4.1 | 27.6 | ±7.8 |
| 7\* | Glc + AA → AP | 0.84 | ±0.09 | 1.16 | ±0.26 | 1.65 | ±0.15 |
| 8 | AP → 1-DG | 290 | ±57 | 507 | ±178 | 994 | ±177 |
| 9 | AP → G | 0.92 | ±0.50 | **0** | ᵃ | **0** | ᵃ |
| 10 | AP → GO | 0.67 | ±0.19 | **0** | ᵃ | **0** | ᵃ |
| 11 | AP → 3-DG | 239 | ±45.7 | 46.7 | ±37.5 | **0** | ᵃ |
| 12 | 3-DG → 3,4-DG | 46.5 | ±10.0 | 117 | ±48 | 139 | ±68 |
| **13** | **3,4-DG → HMF** | **101** | **±29** | **410** | **±108** | **303** | **±96** |
| 14 | Fru → Int | 18.9 | ±13.9 | 90.3 | ±49.6 | 368 | ±37 |
| **15** | **Int → HMF** | **775** | **±559** | **10.1** | **±6.7** | **71.4** | **±8.9** |
| 16 | G → GO | 153 | ±174 | 184 | ±64 | 289 | ±67 |
| 17 | 1-DG → MGO | 1606 | ±455 | 1922 | ±234 | 2420 | ±290 |
| 18 | 1-DG → DA | 47.9 | ±12.1 | 84.9 | ±23.4 | 175 | ±26 |
| **19** | **HMF → P8** | **694** | **±491** | **0** | ᵃ | **233** | **±31** |
| 20\* | 3,4-DG + AA → P7 | 11.8 | ±7.4 | 18.2 | ±25.0 | 30.5 | ±39 |
| 21 | 3-DG → P6 | 970 | ±193 | 908 | ±353 | 1897 | ±660 |
| 22\* | GO + AA → P5 | 616 | ±260 | 2341 | ±737 | 1227 | ±289 |
| 23 | G → P4 | 14.9 | ±10.0 | 22.1 | ±14.5 | 20.8 | ±13 |
| 24\* | DA + AA → P3 | 15.1 | ±6.1 | 3.47 | ±6.0 | **0** | ᵃ |
| 25\* | MG + AA → P2 | 115 | ±34 | 71.0 | ±10.3 | 40.7 | ±7.0 |
| 26\* | 1-DG + AA → P1 | 373 | ±46 | 575 | ±99 | 1061 | ±194 |

### 3.1 ★ THERE IS **NO** TABLE 2 AND **NO** ACTIVATION ENERGY IN THIS PAPER — verified negative `[NEG]`

I grepped the decoded full text for `Arrhenius`, `activation energ`, `Ea`, `Table 2`. **One hit
only**, in the Conclusion, and it is a *disclaimer*, verbatim:

> "Reaction rate constants found in certain steps **do not necessarily follow the Arrhenius
> equation** because parallel and consecutive reactions happening during heating at elevated
> temperatures **may switch one to another**. Water was proposed to play an important role in
> switching of the cascade of these reactions and the origin of α-dicarbonyl compounds formation
> according to the estimated rate constants."

**This paper fits each temperature separately and publishes no Ea, no pre-exponential, and no
reparameterised Arrhenius.** (Contrast the JAFC companion, which does — see §5 of that dossier.)
Anything attributing an activation energy to this DOI is wrong.

### 3.2 Cells where the HPD equals or exceeds the estimate ⇒ **REFUSE**

Step 2 @180 (219 ± 302), step 4 @180 (168 ± 243), step 16 @160 (153 ± 174),
**step 15 @160 (775 ± 559, 72 %)**, **step 19 @160 (694 ± 491, 71 %)**,
step 20 @180 (18.2 ± 25.0) and @200 (30.5 ± 39), step 24 @180 (3.47 ± 6.0),
step 23 @180 (22.1 ± 14.5) and @200 (20.8 ± 13, 63 %).
**10 of 78 cells (12.8 %)** — the same rate as the JAFC companion.

---

## §4. ★ WAVE K5a CROSS-CHECK — a 3-point Arrhenius refit that the authors never performed `[D]`

Because the paper publishes no Ea, this refit is **the only temperature-dependence information that
exists for this system**, and it is derived, never printed. OLS `ln k` vs `1/T` over 433.15 /
453.15 / 473.15 K.

| step | Ea refit (kJ/mol) `[D]` | R² | flag |
|---|---|---|---|
| 1 Glc→1,2-ED | 56.8 | 0.969 | ✔ |
| 3 1,2-ED→Fru | 68.1 | 0.981 | ✔ |
| 5 Glc→G | 144.2 | 0.728 | **non-monotonic** (0.038→1.69→1.38) |
| 6 Glc→3-DG | **71.8** | **0.999** | ★ the cleanest Arrhenius line in the paper |
| 7 Glc+AA→AP | 26.7 | 0.998 | ✔ but suspiciously low (see below) |
| 8 AP→1-DG | 48.7 | 0.994 | ✔ |
| 11 AP→3-DG | **NOT COMPUTABLE** | — | **k = 0 at 200 °C**; k *decreases* 239→46.7→0 |
| 12 3-DG→3,4-DG | 43.7 | 0.880 | ✔ |
| **13 3,4-DG→HMF** | **44.3** | **0.578** | ⚠️ **non-monotonic (101→410→303)** |
| 14 Fru→Int | **117.5** | **1.000** | ★ perfectly log-linear |
| **15 Int→HMF** | **−97.6** | **0.322** | ⚠️⚠️ **catastrophically non-monotonic (775→10.1→71.4): a 77× DROP from 160 to 180 °C** |
| 16 G→GO | 25.0 | 0.934 | low |
| 17 1-DG→MGO | 16.2 | 0.991 | ⚠️ **chemically implausible barrier — the Ma-2022 plateau signature: k moves only 1.5× over 40 °C** |
| 18 1-DG→DA | 51.2 | 0.992 | ✔ |
| 19 HMF→P8 | **NOT COMPUTABLE** | — | **k = 0 at 180 °C**, non-monotonic 694→0→233 |
| 21 3-DG→P6 | 26.1 | 0.655 | non-monotonic |

### 4.1 What the cross-check establishes — and it is the opposite of the JAFC result

**★ FINDING 1 — IN THE WHEAT-FLOUR MATRIX, *NEITHER* HMF LIMB HAS A RECOVERABLE TEMPERATURE
DEPENDENCE, AND THE SINK HAS NONE EITHER.** All three of the paper's HMF-adjacent constants fail:
- `k₁₅` (Int→HMF) falls **77×** between 160 and 180 °C, then rises 7×. R² = 0.322, implied
  Ea = **−97.6 kJ/mol**.
- `k₁₃` (3,4-DG→HMF) is non-monotonic, R² = 0.578.
- `k₁₉` (the sink) is **exactly zero at 180 °C** with 694 ± 491 at 160 °C and 233 ± 31 at 200 °C.

**This is the direct opposite of the amine-free JAFC system**, where `k₆` and `k₇` were the *only*
two internally reproducible steps (R² = 1.000 in both systems). **Same lab, same instruments, same
three temperatures, same 18-vs-26-step formalism, same modelling software — and the fructose-limb
constants go from the most reliable numbers in the paper to the least.** The difference between the
two systems is the presence of a solid flour matrix and 37.7 µmol/g of amine. ⇒ **The fructose-limb
rate constants are matrix-labile, not transferable.**

**★ FINDING 2 — the authors themselves flag a rate constant that DECREASES with temperature and
explain it as a mechanism switch, not a barrier.** Verbatim, §3.4:

> "The rate constant of 3-deoxyglucosone formation from AP (k₁₁) was **decreased by increase in
> temperature** … It was estimated to be 239 × 10⁻³, 0.47 × 10⁻³ and 0 min⁻¹ at 160, 180 and 200 °C.
> Since the amount of water becomes more limited by increase in temperature, the rate of hydrolysis
> is expected to decrease. **However, the rate constant should in principle not decrease with
> increase in temperature for an elementary reaction with a fixed order kinetic and a
> temperature-independent energy of activation.** The observed rate constants indicate a more
> complex mechanism in a place where physical conditions of the system change with temperature."

⚠️ **Transcription note:** the printed value at 180 °C reads `0 47 × 10⁻³` in the running text
(missing decimal point) against **`46.7`** in Table 1 (units ×10³ min⁻¹, i.e. 46.7 × 10⁻³ min⁻¹).
**The text's "0.47" is a typo; Table 1's 46.7 is the value.** Recorded so a later wave does not
ingest 0.47.

**★ FINDING 3 — k₁₇ (1-DG→MGO) is the Ma-2022 plateau artefact, in this corpus, in a Gökmen paper.**
k = 1606 / 1922 / 2420 ×10⁻³ min⁻¹ — a **1.5× move over 40 °C**, giving Ea = 16.2 kJ/mol with
R² = 0.991. Round-3 §A.3.2 refused Ma's Ea(MGO) = 1.84 kJ/mol on exactly this signature (k moved
1.04×). **A high R² on three nearly-flat points is not evidence of a small barrier; it is evidence
the response has saturated.** The paper itself calls these estimates *"very high"* and consistent
with 1-DG's reactivity — i.e. the step is fast enough to be at its ceiling.

**★ FINDING 4 — the "3,4-DG→HMF is faster than 3-DG→3,4-DG" claim holds only at 180 and 200 °C.**
Verbatim, §3.4: *"dehydration of 3-deoxyglucosone to 3,4-dideoxyglucosone is a **rate-limiting
step**, as the rate constant of dehydration of 3,4-dideoxyglucosone to HMF was found to be
faster."* Check `[D]`: k₁₃/k₁₂ = **2.17** (160 °C), **3.50** (180 °C), **2.18** (200 °C). ✔ True at
all three, though the margin is 2–3.5×, not the "almost 5×" the hazelnut paper reports. The authors
then immediately flag their own anomaly, verbatim:

> "The rate constant of 3,4-dideoxyglucosone dehydration to HMF (k₁₃) was **lower at 200 °C than
> 180 °C**. The rate constant **could be expected to increase**, as it did from 160 °C to 180 °C
> … **This discrepancy leads us to think that the reactions of 3,4-dideoxyglucosone should be more
> complex than predicted by the model.** Recent studies demonstrated that 3,4-dideoxyglucosone is an
> intermediate compound in the isomerisation of 3-deoxyglucosone to 3-deoxygalactosone … **However,
> the resolution of chromatography was not enough to separate 3-deoxyglucosone and
> 3-deoxygalactosone in the present study.**"

★ **That last sentence is an author-declared co-elution: the "3-DG" channel (SIM m/z 235) may
contain 3-deoxygalactosone.** It applies to every 3-DG number in every Gökmen-school paper using
this method — i.e. to all five multiresponse papers in this K5a cluster.

---

## §5. THE HMF SECTION — §3.7, transcribed in full because it is the paper's contribution

> "Consecutive removal of 3 water molecules from a hexose sugar ends with the formation of HMF. For
> the formation of the 5-membered ring of HMF, a ring opening is necessary from glucose. On the
> contrary, it has been shown that fructose can also dehydrate **without ring opening, i.e.
> fructofuranose ring intact** (Antal, Mok, & Richards, 1990; Locas & Yaylayan, 2008). Dehydration
> of fructose via cyclic intermediates was considered with an **undetermined intermediate (Int)** in
> the model. **Although fructose was lower in the reaction system, a kinetic model without the
> formation of HMF from fructose did not fit the experimentally observed data (Fig. 5b). If we omit
> HMF formation from fructose through an undetermined intermediate, the amount of HMF was remarkably
> underestimated through the formation only via 3-deoxyglucosone.** It should be noted here again
> that the rate constant of 1-deoxyglucosone formation from fructose was estimated to be zero since
> its formation was predominantly from AP degradation. **However, HMF formation from fructose was a
> kinetically important pathway.**"

> "Nguyen et al (2016) indicated that the origin of HMF formed from glucose and fructose in biscuits
> **depended on the ratio of fructose to glucose**. In the biscuits containing only glucose as a
> sugar source, the rate constant of HMF formation from glucose was reported to be zero and HMF
> formation from fructose was shown to be the predominant pathway. In their study, HMF was
> considered to form **only via caramelisation** of glucose or fructose through the same
> undetermined intermediate. The results of the present paper reveal more detailed investigation of
> HMF formation through the measurement of 3-deoxyglucosone and 3,4-dideoxyglucosone as important
> intermediates. **Intermediates of the HMF formation from fructose in foods are needed to be
> clarified in future studies.**"

**★ AND THE STEP THE REPO MOST NEEDS — the HMF + amine sink, TESTED AND REJECTED here:**

> "Similar to the observations for the reaction of amino acids with 3-deoxyglucosone, **HMF did not
> show a good fit when its degradation included amino acids in a bimolecular reaction. Therefore,
> this step was excluded from the reaction network** (Fig. 2). On the contrary, the involvement of
> HMF in Maillard reaction was reported at elevated temperatures (Nikolov & Yaylayan, 2011a, 2011b).
> **Since the carbonyl group of HMF is conjugated to the aromatic ring its reactivity largely
> depends on temperature.** … According to the estimated reaction rate constants, **amino acids were
> primarily consumed by other more reactive carbonyl compounds** under the investigated conditions
> … **It could be proposed that the reaction of amino compounds with HMF is of minor importance in a
> quantitative way.** However, since the Schiff base adduct of HMF and primary amino acids does not
> involve Amadori rearrangement, **HMF could affect the early stage fate of the Maillard reaction**
> as in the case of predominant acrylamide formation from asparagine/HMF Schiff base."

⚠️ **This is a DIRECT CONFLICT with Şen & Gökmen 2022, which uses `− k₁₉[HMF][AA]` as its only HMF
sink, and with Hamzalıoğlu & Gökmen 2018, which measures the HMF + amino-acid reaction and finds
cysteine consumes 97 % of HMF in 7 days at 50 °C.** Same lab, three papers, three different answers
about whether the HMF + amine edge exists. See the synthesis §4. **The reconciliation the authors
offer is competitive kinetics, not absence: at 160–200 °C in flour, GO / 1-DG / MGO out-compete HMF
for the amine pool.**

---

## §6. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL`

| # | constraint | anchor |
|---|---|---|
| S1 | **HMF formation from fructose is a key step; a 3-DG-only model remarkably underestimates HMF** | §3.7, Fig. 5b |
| S2 | **1,2-enolisation IS rate-limiting when amines are present, and is NOT when they are absent** — the authors run the omit-1,2-ED test in both systems | §3.2: *"This kinetic evaluation and model discrimination showed that the enolisation in the wheat flour system is a **rate-limiting step**… It was not sensible to omit that step"*; contrast the JAFC paper, lines 280–282, where 1,2-ED is omitted because interconversion is *"obviously fast"*. ★ **A same-lab demonstration that a network's topology is amine-dependent.** |
| S3 | **Fru → 1-DG is zero here; 1-DG comes from AP.** In the amine-free JAFC system 1-DG comes *only* from fructose. | §3.4 vs JAFC lines 371–374 — ★ **the same node changes parent between the two matrices** |
| S4 | **HMF + AA is rejected by model discrimination** | §3.7, quoted above |
| S5 | **3-DG + AA is also rejected**, because 3-DG is *"a thermodynamically more stable intermediate"* | §3.6 |
| S6 | **MGO comes from 1-DG here (k₂₈ = 3-DG→MGO is zero at all T)** — the reverse of the JAFC system | §3.5 |
| S7 | Reactivity order of dicarbonyls toward amino acids: **glyoxal > 1-DG > methylglyoxal > diacetyl > 3,4-dideoxyglucosone** | §3.6, verbatim |
| S8 | **k for GO+AA, MGO+AA and DA+AA all *decrease* with temperature — a volatility artefact, not chemistry** | §3.6: *"This can be due the high volatility of these compounds leaving the system to the headspace of the tube. The rate constants of 3,4-dideoxyglucosone and 1-deoxyglucosone increased, as they are not volatile."* ⚠️ same transport contamination as the JAFC paper's S8 |
| S9 | **No acetic or formic acid was observed**, unlike every aqueous model | §3.5 |
| S10 | 3-DG : 1-DG ratio = **10.6 ± 2.8 / 7.8 ± 2.7 / 6.9 ± 1.5** at 160/180/200 °C; 3-DG : glucosone ≈ **26** | §3.1 `[M]` |
| S11 | Fructose's higher Maillard reactivity is a **melting-point** effect, not a chemical one, and **does not apply here** because the fructose is generated *after* glucose melts | §3.2, citing Robert 2004 |
| S12 | Glucose melts at **146 °C**; fructose at **104 °C** | §3.2 `[C]` |

---

## §7. VERIFIED NEGATIVES `[NEG]`

- **No Ea, no Arrhenius fit, no pre-exponential.** (§3.1)
- **No pH measurement.**
- **No water activity**; only "7 % humidity" by AACC 44-15.02.
- **No browning / colour / melanoidin measurement.** The Conclusion asks for them to be added.
- **No furfural.**
- **No absolute HMF concentration table** — HMF vs time exists only in Fig. 1 (raster).
- **Amadori product never quantified** (no standard available).
- **3-DG and 3-deoxygalactosone were not chromatographically resolved** (author-declared, §3.4).
- Supplementary material (mass balance; Figs. S1–S10 equivalents) is **not in the PDF on disk**.

---

## §8. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| Table 1, all 78 k ± HPD cells | `[F]` | **USE-Q** — per-temperature constants for **this matrix only**, `min⁻¹×10³` (or `g·µmol⁻¹·min⁻¹×10³` for the 6 starred steps); drop the 10 cells of §3.2 |
| **k₁₄, k₁₅ absolute (fructose limb)** | `[F]` | **REFUSE for k₁₅** (77× non-monotonic drop, HPD 72 % at 160 °C); `k₁₄` **RATIO-ONLY** (`[Int]` unmeasured) |
| **k₁₃ (3,4-DG→HMF)** | `[F]` | **USE-Q, per-temperature only** — non-monotonic; also carries the semi-quant response factor of 3,4-DG |
| **k₁₉ (HMF sink)** | `[F]` | **REFUSE** — zero at 180 °C, 71 % HPD at 160 °C |
| k₆ (Glc→3-DG) | `[F]` | ★ **USE** — the best-behaved row in the paper (R² of refit 0.999) |
| k₇ (Glc + AA → AP) | `[F]` | **USE-Q** — but `[AA]` is a free+bound lump and `[AP]` is unmeasured |
| **any activation energy** | — | ★ **DOES NOT EXIST IN THIS PAPER. REFUSE any attribution.** |
| Ea refits of §4 | `[D]` | **PRIOR-ONLY**, and only for steps 1, 3, 6, 8, 12, 18 (R² ≥ 0.88 and monotonic). **Never for 5, 11, 13, 15, 17, 19, 21.** |
| initial concentrations, peak fructose, 3-DG:1-DG ratios (S10) | `[M]` | **USE** |
| S1–S12 | `[M]`/`[F]` | **STRUCTURAL** — S2, S3, S6 are the highest-value rows: three same-lab topology switches driven by matrix |
| 1-DG, glucosone, 3,4-DG concentrations | `[M]` | **`absolute_concentration: false`** (author-declared semi-quantitation) |
| every 3-DG number | `[M]` | **flag `possible_3-deoxygalactosone_coelution: true`** (author-declared) |
