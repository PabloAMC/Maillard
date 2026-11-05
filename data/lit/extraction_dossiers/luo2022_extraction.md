# Luo, Xu, Li & Ho 2022 (J. Food Bioact. 20:90-97) — EXTRACTION (file stem `luo2022`)
**DOI:** 10.31665/JFB.2022.18332 (printed on p. 90; Crossref 2026-11-05: 10.31665/jfb.2022.18332, authors Yue Luo, Baojun Xu, Shiming Li, Chi-Tang Ho, J. Food Bioact. pp. 90-97, issued 2022-12) — matches the brief.
### *"Reactivity comparison of Amadori rearrangement product formation of glycine, diglycine, and triglycine in a competing Maillard system"*
### **Glucose + glycine + diglycine (R1) or glucose + glycine + triglycine (R2), 20 mM each (1:1:1), in Milli-Q water, initial pH 3 / 5 / 7 / 9**, sealed pressure vessels in an oil bath (600 rpm), **80 / 100 / 130 / 160 °C, 60 and 120 min**; Gly-ARP (GA), Digly-ARP (DA), Trigly-ARP (TA) by HILIC-LC-MS/MS (MRM) against purified standards; glycine's own ARP separated from peptide-derived glycine by a ¹⁵N-glycine twin batch. **In one pot at fixed glucose, DA/GA and TA/GA are 0.43-2.33 at 80-100 °C and pH 3-7 (near 1), 0.08-0.17 at pH 9, and 0.02-1.15 at 130-160 °C.** No replicates or SDs printed; pH after heating measured but not printed.

**Read 2026-11-05 (first reader) from `data/articles/luo2022.pdf`** (8 pp., journal PDF, clean text layer). Prose and
Tables 1-3 from `pdftotext -raw`; Table 3 ratios recomputed from Tables 1-2 [Z] (they agree to the printed 2-3 digits).
Figs. 1-3 are MS spectra; Fig. 4 (A420) not digitised. Arithmetic: scratchpad `m4tools/m4_arith.py`.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Yue Luo, **Baojun Xu** (corresponding), Shiming Li, **Chi-Tang Ho** (corresponding) |
| Affiliation | Dept. of Food Science, Rutgers University, New Brunswick NJ (Luo, Li, Ho); Food Science and Technology Program, BNU-HKBU United International College, Zhuhai (Xu) |
| Venue | *Journal of Food Bioactives* **20** (2022) 90-97 |
| DOI | 10.31665/JFB.2022.18332 |
| Dates | received 3 Dec 2022; revised and accepted 28 Dec 2022 |

**Laboratory.** Rutgers (Ho) with UIC Zhuhai (Xu). The case calls it "Ho network, as Xia 2022"; `xia2022` is the Jiangnan
group (Zhang X, with Ho; `data/keys/laboratories.yml` `jiangnan_zhang_x`). The experiments here were run at Rutgers /
UIC (acknowledgement thanks Heping Cui of Jiangnan for advice only). First-reader view: a different bench from Xia 2022
sharing a senior author; whether Ho co-authorship merges them for leave-one-laboratory-out is an owner call. Same group
as `luo2023`.

## 1. What this paper answers for the model

| model need | what Luo 2022 gives | tag |
|---|---|---|
| amino acid vs peptide competing for one sugar in one pot | GA, DA, TA (mM) at 4 pH × 4 T × 2 times (Tables 1-2) and their ratios (Table 3) | [M] |
| peptide-ARP degradation to a shorter-peptide ARP | DA formed from TA in R2 at 130-160 °C: 0.0024-0.0073 mM | [M] |
| ARP yields | GA up to 0.38 mM (1.9 % of glucose) at 130 °C / 1 h / pH 7 (R1) | [M] / [Z] |
| amine consumption, glucose loss, pH after heating | **not printed** | — |
| thiol, MFT, NF | **none** | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots | R1: glucose / glycine / diglycine 1:1:1; R2: glucose / glycine / triglycine 1:1:1; **0.02 mol/L each**, Milli-Q water; each in a normal-glycine batch (NGB) and a ¹⁵N-glycine batch (15NGB) |
| pH | adjusted to 3.0 / 5.0 / 7.0 / 9.0 before heating (NaOH / HCl); label pH; **measured after heating ("pH and UV-Vis analysis") but not reported** [!] |
| Vessel / heating | "high-temperature high-pressure resistant reaction vessels" (volume and fill not stated), sealed, magnetic stirring 600 rpm, oil bath, 80 / 100 / 130 / 160 °C, 60 and 120 min; ice-bath cooling |
| Headspace | not stated; air presumed |
| Replicates | **not stated; no SDs** |

## 3. Assays [M]

- **ARPs:** HILIC (Waters, 5 µm, 4.6 × 250 mm) MS/MS, MRM 238.1 → 220.1 / 202.0 (GA), 295.2 → 277.1 / 211.1 (DA),
  352.1 → 334.1 / 268.0 (TA); external standards = purified GA, DA, TA (two-step purification, ≥ 95 % by LC-MS/MS and
  NMR, data not shown). **No LOD / LOQ, recovery or matrix-effect check printed.**
- **Glycine's own ARP:** GA(actual) = GA(NGB) − GA(15NGB), i.e. a difference of two batches; ¹⁵N-GA is invisible in the
  GA MRM channel, so the 15NGB channel counts only glycine released from the peptides [Q]. [!] The difference inherits
  both batches' noise and is not checked against a ¹⁵N-GA standard.
- A420 by plate reader (Fig. 4).

## 4. Data

### 4.1 Table 1 (verbatim) — R1, glucose / glycine / diglycine, mmol/L

| | GA pH 3 | GA pH 5 | GA pH 7 | GA pH 9 | DA pH 3 | DA pH 5 | DA pH 7 | DA pH 9 |
|---|---|---|---|---|---|---|---|---|
| 80 °C 1 h | 0.0021 | 0.0033 | 0.0035 | 0.0695 | 0.0025 | 0.0035 | 0.0068 | 0.0059 |
| 80 °C 2 h | 0.0021 | 0.0043 | 0.0058 | 0.1099 | 0.0021 | 0.0043 | 0.0081 | 0.0092 |
| 100 °C 1 h | 0.0049 | 0.0061 | 0.0063 | 0.1905 | 0.0047 | 0.0088 | 0.0146 | 0.0143 |
| 100 °C 2 h | 0.0080 | 0.0172+ | 0.0218 | 0.3049 | 0.0077 | 0.0159 | 0.0282 | 0.0254 |
| 130 °C 1 h | 0.0563 | 0.2000 | 0.3825 | 0.2724 | 0.0465 | 0.1029 | 0.1235 | 0.1469 |
| 130 °C 2 h | 0.0637 | 0.3514 | 0.3238 | 0.2520 | 0.0479 | 0.0981 | 0.1121 | 0.1440 |
| 160 °C 1 h | 0.0543 | 0.1832 | 0.0793 | 0.0260 | 0.0523 | 0.0690 | 0.0739 | 0.0462 |
| 160 °C 2 h | 0.0326 | 0.2963 | 0.2223 | 0.0421 | 0.0315 | 0.0397 | 0.0439 | 0.0136 |

[!] "0.0172+" is printed with a stray "+".

### 4.2 Table 2 (verbatim) — R2, glucose / glycine / triglycine, mmol/L

| | GA pH 3 | GA pH 5 | GA pH 7 | GA pH 9 | TA pH 3 | TA pH 5 | TA pH 7 | TA pH 9 | DA (from TA) pH 3 / 5 / 7 / 9 |
|---|---|---|---|---|---|---|---|---|---|
| 80 °C 1 h | 0.0053 | 0.0128 | 0.0105 | 0.0876 | 0.0042 | 0.0123 | 0.0101 | 0.0129 | — |
| 80 °C 2 h | 0.0042 | 0.0065 | 0.0091 | 0.1166 | 0.0055 | 0.0065 | 0.0080 | 0.0198 | — |
| 100 °C 1 h | 0.0065 | 0.0161 | 0.0251 | 0.1604 | 0.0053 | 0.0129 | 0.0107 | 0.0213 | — |
| 100 °C 2 h | 0.0085 | 0.0211 | 0.0262 | 0.2502 | 0.0074 | 0.0189 | 0.0139 | 0.0222 | — |
| 130 °C 1 h | 0.0854 | 0.1662 | 0.1773 | 0.1319 | 0.0758 | 0.0988 | 0.1711 | 0.1299 | 0.0027 / 0.0029 / 0.0031 / 0.0073 |
| 130 °C 2 h | 0.1205 | 0.1231 | 0.1164 | 0.0788 | 0.0623 | 0.0996 | 0.1108 | 0.0896 | 0.0026 / 0.0025 / 0.0034 / 0.0069 |
| 160 °C 1 h | 0.0326 | 0.2506 | 0.0642 | 0.0520 | 0.0132 | 0.0111 | 0.0107 | 0.0209 | 0.0034 / 0.0048 / 0.0029 / 0.0056 |
| 160 °C 2 h | 0.0156 | 0.0593 | 0.0568 | 0.0487 | 0.0143 | 0.0018 | 0.0010 | 0.0024 | 0.0048 / 0.0030 / 0.0031 / 0.0035 |

### 4.3 Table 3 (verbatim) — ratios

| | DA/GA pH 3 | TA/GA pH 3 | DA/GA pH 5 | TA/GA pH 5 | DA/GA pH 7 | TA/GA pH 7 | DA/GA pH 9 | TA/GA pH 9 |
|---|---|---|---|---|---|---|---|---|
| 80 °C 1 h | 1.21 | 0.78 | 1.06 | 0.96 | 1.96 | 0.97 | 0.08 | 0.15 |
| 80 °C 2 h | 1.01 | 1.32 | 0.98 | 1.00 | 1.40 | 0.88 | 0.08 | 0.17 |
| 100 °C 1 h | 0.96 | 0.82 | 1.44 | 0.80 | 2.33 | 0.43 | 0.08 | 0.13 |
| 100 °C 2 h | 0.96 | 0.87 | 0.92 | 0.90 | 1.30 | 0.53 | 0.08 | 0.09 |
| 130 °C 1 h | 0.83 | 0.89 | 0.51 | 0.59 | 0.32 | 0.96 | 0.54 | 0.98 |
| 130 °C 2 h | 0.75 | 0.52 | 0.28 | 0.81 | 0.35 | 0.95 | 0.57 | 1.14 |
| 160 °C 1 h | 0.96 | 0.40 | 0.38 | 0.04 | 0.93 | 0.17 | 1.15 | 0.40 |
| 160 °C 2 h | 0.97 | 0.92 | 0.13 | 0.03 | 0.20 | 0.02 | 0.32 | 0.05 |

DA/TA (Table 3 lower block), pH 3 / 5 / 7 / 9: 80 °C 1 h 1.6 / 1.1 / 2.0 / 0.6; 2 h 0.8 / 1.0 / 1.6 / 0.5; 100 °C 1 h 1.2 /
1.8 / 5.4 / 0.6; 2 h 1.1 / 1.0 / 2.4 / 0.9; 130 °C 1 h 0.9 / 0.9 / 0.3 / 0.5; 2 h 1.5 / 0.3 / 0.4 / 0.5; 160 °C 1 h 2.4 /
8.5 / 5.6 / 4.4; 2 h 1.1 / 4.4 / 11.8 / 6.6. [!] DA and TA come from different pots (R1, R2), so DA/TA is a between-pot
ratio, unlike DA/GA and TA/GA.

### 4.4 Arithmetic [Z] (`m4_arith.py`)

1. **Within-study reproducibility.** GA is measured in both R1 and R2 under the same nominal glucose / glycine charge (the
   third amine differs). log10(GA_R2 / GA_R1) over the 32 cells: mean +0.01, **rms 0.34 dex** (max 0.70); if all of it is
   noise, ~0.24 dex per value. The partner peptide may also change GA, so this is an upper bound on noise only if the
   partner has no effect.
2. Yields: GA ≤ 0.38 mM = 1.9 mol % of glucose; DA ≤ 0.15 mM; TA ≤ 0.17 mM. At 80-100 °C and pH 3-7 all ARPs are
   ≤ 0.03 mM (≤ 0.15 %).
3. pH 9 at 80-100 °C: GA 0.07-0.30 mM vs DA / TA 0.006-0.025: glycine's ARP ×6-13 the peptides'.

## 5. What it can and cannot test in the model

**Can:**
- **Amine competition for one sugar inside one pot** (free α-amino acid vs di- / tripeptide), at four pH and four
  temperatures; a scale-free ratio (no absolute calibration needed if response factors are right).
- The temperature at which peptide ARPs stop accumulating (≥ 130 °C) relative to Gly-ARP.

**Cannot:**
- Formation constants: the ratios are of standing ARP levels at two times, mixing formation and degradation (the
  authors say so [Q]).
- Glucose and amine consumption, pH drift (measured, not printed), replicate noise (none printed).
- Thiol chemistry, cysteine, MFT, NF.

## 6. Conflicts and agreements with the repo's records

1. **[!] Against `xia2022`'s amine-class ratio (D7 of case M4).** Xia 2022 (separate pots, 200 + 200 mM, water, pH 7.5
   initial): TriG-ARP forms ×6.9 (100 °C) and ×5.0 (130 °C) faster than G-ARP (case sec. 2.4). Luo 2022 (competing in one
   pot, 20 mM each, initial pH 7): TA/GA **0.43-0.97** at 80-100 °C and 0.95-0.96 at 130 °C; DA/GA
   1.30-2.33 at 80-100 °C. Ratios of standing levels, not rates, and a 10× lower concentration — but the sign of
   the triglycine advantage is not reproduced in competition. Both papers carry Ho as an author; they are not two
   laboratories for this contrast in any conservative count.
2. [!] At pH 9 glycine's ARP dominates (×6-13) although glycine's amino group has the higher pKa (9.60, quoted) — the
   opposite of the authors' own unprotonated-amine argument; unexplained.
3. [!] No replicates, no SDs, stray "+" in Table 1, pH after heating not printed; the "actual GA" is a difference of two
   batches (sec. 3).
4. Agrees with `luo2023` (same group) that peptide ARPs are less heat-stable than Gly-ARP above 100 °C.

## 7. Provenance tiers (first reader only; second reader pending)

- Tables 1-3: `direct_measurement` (printed), LC-MS/MS with purified external standards, confidence **low-medium** (no
  replicates or SDs; GA is a two-batch difference; 0.34 dex rms between R1 and R2 for the same nominal GA).
- Ratios and yields in sec. 4.4: [Z], `derived`.

## M4 role (case 8660c005)

**Role fixed before reading:** part A evidence — peptide vs amino-acid competition for one sugar inside one study (D7;
Ho network, as Xia 2022); measurement that would count: "ARP ratios at fixed sugar across pH and T"; P ~0.5.

- **The named measurement is met:** DA/GA and TA/GA in one pot at fixed glucose, 4 pH × 4 T × 2 times.
- **What it supplies (D7, part A).** Within-pot ratios of peptide-ARP to amino-acid-ARP (glucose, 20 mM each, water,
  initial pH 3-9, 80-160 °C, 1-2 h): **0.43-2.33 at 80-100 °C and pH 3-7 (geometric centre ≈ 1); 0.08-0.17 at pH 9;
  0.02-1.15 at 130-160 °C** (peptide ARPs degrade). Transfer: glycine oligopeptides only (no Cys, no protein lysine),
  water, label pH, standing levels; noise ≈ 0.24-0.34 dex per value (sec. 4.4) with no printed SDs.
- **Unexpected, bearing on D7 (not re-roled).** D7 reads A-R3 ("equal per alpha-amine") as contradicted by Xia's ×5-7 at
  equal molarity. In competition at fixed sugar, the same network's later study finds peptide / amino-acid ARP ratios
  near 1 at 80-100 °C and pH 3-7, i.e. what A-R3 predicts there, and far below 1 at pH 9 and ≥ 130 °C. On standing
  levels this does not contradict A-R3 at the low-temperature, acid-to-neutral corner and contradicts A-R1's peptide
  advantage there. D7 is one-laboratory, development-grade evidence in the case; this paper adds a within-pot contrast
  that points the other way from Xia 2022, so the D7 row cannot be read as one-directional.
- **Decisive criterion:** criterion (d) ("an amine-class contrast at fixed sugar") is the one this paper was chosen for,
  and it holds that measurement; it does not touch the part G decisive rows D1-D3 (no thiol, no MFT).
