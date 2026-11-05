# Du, Wang, Ma, Li, Bai, Cui, Fan & Wang 2024 (LWT 192:115741) — EXTRACTION (file stem `du2024b`)
**DOI:** 10.1016/j.lwt.2024.115741 (printed on p. 1; Crossref 2026-11-05: same title, authors Wenbin Du, Yutang Wang, Qianli Ma, Xuejie Li, Shuang Bai, Weiye Cui, Bei Fan, Fengzhong Wang, LWT 192, article 115741, issued 2024-01) — matches the brief.
### *"Effect of 2-alkenals on the Maillard reaction of cysteine-glucose: Initial stage intermediate formation and flavor compound generation pathways"*
### **L-Cysteine 1 mmol + D-glucose 1 mmol (+ (E)-2-hexenal, -heptenal or -octenal 0.5 mmol) in 5 mL 0.2 M phosphate pH 6.50**, 15 mL screw-cap glass vials, oil bath **145 °C**: non-volatiles (Cys, Cys-Amadori, three alkenal thiazolidines) by UPLC-MS/MS at 10-180 min (Fig. 2, g/L, figure only); volatiles at 120 min by SPME-GC-MS(-O) (Table 1, µg/L, one internal standard) and a 1:1 [13C6]/[12C6]-glucose CAMOLA (Table 2). **The cysteine + glucose control (CG): Cys-Amadori peaks at ~35 g/L (≈ 124 mM, 62 % of the cysteine) at 10 min and falls to ~5.5 g/L by 30 min — but cysteine has dropped only to 158 mM by then, so Cys + Cys-Amadori exceeds the cysteine charged by 40 % [!].** No volatile data for the control pot; MFT only in the 2-heptenal pot (6.30 µg/L); no NF; no CAMOLA on MFT.

**Read 2026-11-05 (first reader) from `data/articles/du2024b.pdf`** (10 pp., Elsevier open-access CC BY PDF, clean
text layer). Prose and Tables 1-2 from `pdftotext -raw`. **Fig. 2 (p. 4) read** from a 400 dpi render (scratchpad
`m4img/du2024b_400-04.png`, crops `m4img/du_amad.png`, `m4img/du_cys.png`; dark-pixel runs for the axis ticks,
`m4tools/runs.py`): **read error ±1.5 g/L for Cys-Amadori and cysteine** (four overlapping line styles; the solid CG line
is separable at 0-15 min only), ±0.5 g/L for the thiazolidines. Arithmetic: scratchpad `m4tools/m4_arith.py`. SI
(Table S1 MRM, Table S2 lipid-derived volatiles, Fig. S1 NMR) not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Wenbin Du, Yutang Wang, Qianli Ma, Xuejie Li, Shuang Bai, Weiye Cui, **Bei Fan** (corresponding), **Fengzhong Wang** (corresponding) |
| Affiliation | Institute of Food Science and Technology, Chinese Academy of Agricultural Sciences (CAAS), Beijing (Du, Wang, Cui, Fan, Wang); First Hospital of Lanzhou University (Ma); BTBU Beijing (Li, Bai) |
| Venue | *LWT - Food Science and Technology* **192** (2024) 115741 |
| DOI | 10.1016/j.lwt.2024.115741 |
| Dates | received 6 Oct 2023; revised 31 Dec 2023; accepted 7 Jan 2024; online 9 Jan 2024 |

**Laboratory.** Fan / Wang group, CAAS Beijing — **the same laboratory as `du2023`, `du2024` and `du2025`** (Du 2023 is
case M4's D3 source: Cys + 1:1 [13C6]-glucose, pH 6.5, 150 °C, MFT M+0 : M+5). Case M4 sec. 2.9 calls this paper
"BTBU"; the corresponding authors are CAAS, with two BTBU co-authors. Either way it is D3's laboratory.

[!] **Possible duplicate control.** `du2024` (Food Chem.; same group) shows a Cys-Amadori curve for its cysteine +
glucose control that peaks at ~34 g/L at 10 min and is gone by 60 min (its dossier, Fig. 2 "M"); this paper's CG curve
has the same peak, time and decay. The two papers may report one control pot twice. `du2024` did not state that pot's
charge or temperature; this paper does (200 + 200 mM, 145 °C).

## 1. What this paper answers for the model

| model need | what Du 2024b gives | tag |
|---|---|---|
| Cys-hexose Amadori time course, Cys + glucose, stated charge and T | Fig. 2 CG: ~0 / 35 / 16 / 5.5 / 3 / ~0.8 / ≤ 0.5 g/L at 0 / 10 / 15 / 30 / 45 / 60 / 120-180 min | [F] |
| cysteine loss in the same pot | Fig. 2 CG: 24.3 / 19.1 / 13.3 / 7.4 / 3.3 / ~1.6 g/L at 0 / 10 / 15 / 30 / 45 / 180 min | [F] |
| alkenal competition for cysteine | Cys-Amadori at 10 min CG > OCT > HEP > HEX [Q]; thiazolidines PTCA / HTCA 1 / HTCA 2 peak at 10 / 30 / 30-45 min | [Q] / [F] |
| MFT | only in HEP (2-heptenal pot): **6.30 ± 0.28 µg/L** at 120 min; "–" (not detected) in HEX and OCT; **no CG volatiles** | [M] |
| NF, FFT, furfural, H2S | none reported (2-furanmethanol 43.67 µg/L in HEP; 5-methyl-2-furancarboxaldehyde 10.5-99 µg/L) | [M] |
| CAMOLA (1:1 [13C6]/[12C6]-glucose) | Table 2: twelve compound × system rows; MFT not among them | [M] |
| replicate noise | n = 3; Fig. 2 bars: Cys-Amadori ±3.5 g/L at the CG peak; Table 1 SDs 1-22 % | [F] / [M] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots | HEX / HEP / OCT: Cys 1 mmol + glucose 1 mmol + alkenal 0.5 mmol; **CG control: Cys 1 mmol + glucose 1 mmol**; in **5 mL 0.2 M phosphate pH 6.50** [Z: 200 + 200 (+ 100) mM; buffer : cysteine 1 : 1] |
| Vessel / heating | 15 mL glass vial with screw cap (10 mL headspace, air presumed), magnetic stirring, oil bath **145 °C**; non-volatiles sampled at 10, 15, 30, 45, 60, 120, 180 min; volatiles at 120 min |
| pH | label pH 6.50 (buffer as purchased); **not re-measured**; 0.2 M phosphate against 200 mM cysteine and glucose-derived acids may drift |
| CAMOLA pots | as above with glucose 0.5 mmol [13C6] + 0.5 mmol [12C6] |
| Replicates | "averages of three replicates" [V] |
| Standards | Cys-Amadori synthesised (Cys + glucose, "0.2 mol" each in 5 mL [!], 55 °C 1 h then 85 °C 6 h, AG 50W-X4 column, ≥ 95 % by HPLC-UV); thiazolidines from Cys + alkenal at 90 °C 30-40 min, polyamide column, ≥ 95 % |

[!] "0.2 mol" of cysteine and glucose in 5 mL (40 M) is impossible; `du2024` prints 0.02 mol for the same preparation.

## 3. Assays [M]

- **Non-volatiles:** UPLC-MS/MS (Xevo TQ-XS, ESI+, MRM in Table S1, not on disk), **external standard = the synthesised
  intermediates**. No LOD, LOQ, recovery or matrix check printed.
- [!] **The Cys-Amadori standard's identity is weakly documented:** printed [M+H]+ 284.1208 vs "calcd" 284.0804 (142
  ppm; the true monoisotopic [C9H18NO7S]+ is 284.0799 [Z]); and the text states "ARP was not found in the NMR data of the
  ARP sample" [V], reporting two six-membered-ring isomers instead. If the standard's purity or identity is off, every
  Cys-Amadori concentration scales with it.
- **Volatiles:** SPME (CAR/PDMS/DVB, 50 °C 30 + 30 min), GC-MS DB-Wax, 1 µL 1,2-dichlorobenzene internal standard
  (concentration printed "100 m/L" [!]), semi-quantified with RF = 1; GC-O by three sniffers. No LODs.
- **CAMOLA:** isotopomer shares from M+ to M+n, corrected for natural 13C, 33S, 34S and M−1.

## 4. Data

### 4.1 Fig. 2 — non-volatiles (g/L) [F]

| t (min) | 0 | 10 | 15 | 30 | 45 | 60 | 120 | 180 |
|---|---|---|---|---|---|---|---|---|
| **CG Cys-Amadori** | ~0 | **35** (SD ~3.5) | 16 | 5.5 | 3 | 0.8 | ≤ 0.5 | ≤ 0.5 |
| HEX / HEP / OCT Cys-Amadori | ~0 | 20 / 22 / 24.5 | 10-17 | 4-7 | 2-4 | ~1 | ≤ 0.5 | ≤ 0.5 |
| **CG cysteine** | 24.3 | 19.1 | 13.3 | 7.4 | 3.3 | — | — | ~1.6 |
| HEX / HEP / OCT cysteine | 24 | 13-17 | 7.5-8 | 2.5-3 | 2-2.5 | — | — | 0.6-1.3 |
| PTCA (HEX) | 0 | 10.3 | 5.8 | — | — | 3.8 | 2.6 | 2.1 |
| HTCA 1 (HEP) | 0 | 3.5 | — | 6.6 | 3.0 | 2.7 | 1.3 | 1.3 |
| HTCA 2 (OCT) | 0 | 2.0 | — | 4.2 | 4.4 | 2.7 | 2.6 | 2.4 |

### 4.2 Table 1 — volatiles at 120 min (µg/L, mean ± SD, n = 3) [M]

87 compounds in HEX / HEP / OCT (no CG column). Sums: thiophenes 1000.35 / 1047.96 / 1209.25; thiazoles 305.53 /
333.32 / 434.61; furans 57.74 / 373.99 / 120.66; furanones 0 / 134.27 / 122.70; pyridines 150.84 / 0 / 247.93; others
532.91 / 220.80 / 675.60; total 2047.37 / 2110.33 / 2810.75. Rows bearing on M4: **2-methyl-3-furanthiol (No. 76)** – /
**6.30 ± 0.28** / –; 2-furanmethanol – / 43.67 ± 5.39 / –; 5-methyl-2-furancarboxaldehyde – / 99.38 ± 8.23 / 10.52 ± 0.59;
3-thiophenethiol – / 2.49 ± 0.28 / –; 2-acetylthiazole 17.06 / 17.70 / 17.04; 4,5-dimethylthiazole 234.71 / – / –. The
remaining rows are lipid-alkenal products (M5's chemistry) and are not transcribed.

[!] Table 1's column headers read "HEX HEP OTC"; Table 2 labels the heptenal pot "HPT".

### 4.3 Table 2 — CAMOLA isotopomer shares (%, M+0 … M+n) [M]

| No. | compound (M+) | system | shares M+0, M+1, M+2, … |
|---|---|---|---|
| 2 | 2-ethylthiophene (112) | HEX | 99, <1 … |
| 4 | 2-propylthiophene (126) | HEX / HPT | 95, 4, <1 … / 99 … |
| 8 | 2-butylthiophene ("196" [!], M = 140) | HEX / HPT / OCT | 85, <1, 15 … / 96, 4 … / 99 … |
| 21 | 1-(2-thienyl)-ethanone (126) | HEX / HPT / OCT | 42, <1, 13, <1, 40, <1, 5 / 41, <1, 14, <1, 40, <1, 5 / 41, <1, 14, <1, 41, <1, 4 |
| 20 | 2-methyl-5-thiophenecarboxaldehyde (126) | HEX / HPT / OCT | 45, <1, 10, <1, 40, <1, 5 / 44, <1, 9, <1, 41, <1, 6 / 45, <1, 9, <1, 40, <1, 6 |
| 26 | 2-methylthieno[2,3-b]thiophene (154) | OCT | 47, <1, 5, <1, 39, <1, 6, <1 |
| 32 | 2-butanoylthiophene (154) | OCT | 99, <1 … |
| 45 | 2-acetylthiazole (127) | HEX / HPT / OCT | 61, <1, <1, 39 … / 60 … 40 / 60 … 40 |
| 40 | 4,5-dimethylthiazole (113) | HEX | 65, 1, 2, 1, 28, 3 |
| 82 | 2,3-dihydro-1H-indole (119) | HEX | 44, <1, 7, <1, <1, <1, 42, <1, 9 |

[!] The authors describe Nos. 20 and 26 as "approximately 1:1" unlabelled vs labelled with an intact glucose C5 [Q], but
the printed labelled share sits in the M+4 column (40 %), not M+5; either a column shift in the table or a C4 + C1 / C2
assembly. The text-layer column order matches the printed table.

### 4.4 Arithmetic [Z] (`m4_arith.py`)

1. **CG Cys-Amadori in molar terms** (283.3 g/mol): 124 mM at 10 min (61.8 % of the 200 mM cysteine), 56 mM at 15,
   19 mM at 30, 11 mM at 45, 3 mM at 60 min.
2. **Mass balance [!]:** at 10 min cysteine 158 mM + Cys-Amadori 124 mM = 281 mM > 200 mM charged (×1.4); Cys lost by
   10 min (43 mM) is a third of the Amadori formed (124 mM). Not explainable by read error (±1.5 g/L = ±5 mM Amadori,
   ±12 mM Cys).
3. **Net Cys-Amadori loss after the peak** (formation continuing, so a lower bound on the decay constant): 0.157 min⁻¹
   (10 → 15), 0.093 (10 → 30), 0.056 min⁻¹ (15 → 45).
4. **Comparison band:** `deng2021`'s Met-glucose ARP decay (0.0048 min⁻¹ at 100 °C, label pH 5.5, slope 0.2 dex per pH
   unit) moved to 145 °C with Ea 85.7-107 kJ/mol and to pH 6.5: 0.15-0.31 min⁻¹. Du's net 0.06-0.16 min⁻¹ lies at and
   below the low end — consistent, as a lower bound, with the band's centre within ×0.3-3.
5. Formation: ≥ 124 mM in ≤ 10 min from 200 + 200 mM (if the scale holds): ≥ 12 mM min⁻¹ average.

## 5. What it can and cannot test in the model

**Can:**
- The existence, size and lifetime of the cysteine-glucose Amadori pool in a buffered Cys + glucose pot at 145 °C and pH
  6.5 — close to Hofmann's T1 temperature (145 °C) and to Du 2023 / Zhao 2019b's pots (pH 6.5 phosphate).
- Amadori competition: alkenals lower the 10 min Cys-Amadori peak by 30-43 % and consume cysteine faster.

**Cannot:**
- MFT or NF in a cysteine + glucose pot (no CG volatiles); hexose MFT's carbon origin (no MFT CAMOLA).
- A clean Cys-Amadori constant: two time points before the peak (0, 10 min), net loss only, and the quantification fails
  the cysteine mass balance (sec. 4.4) with a weakly documented standard (sec. 3).
- pH drift (not measured); headspace O2 (not stated).

## 6. Conflicts and agreements with the repo's records

1. **[!] Internal mass balance** (sec. 4.4 item 2): either the Cys-Amadori response is ~×2-3 too high (standard purity /
   identity), or cysteine is under-recovered at 0 min and over-recovered later; the paper does not discuss it.
2. **[!] Probable duplicate of `du2024`'s control curve** (sec. 0); not independent evidence.
3. Agrees with `zhang2026` (GSH + glucose, 130-150 °C, BTBU) and `deroos2005` Table II that a cysteine-hexose Amadori
   compound forms readily; here it is gone within ~60 min at 145 °C.
4. Agrees with `du2023` (same laboratory) on the system and pH; MFT in Du 2023's control pot is the D3 datum, and this
   paper adds none.
5. [!] Smaller defects: "0.2 mol" in 5 mL; 284.1208 vs 284.0804; "Cystine-amadori" in the abstract and keywords for a
   cysteine compound; "196" as the M+ of 2-butylthiophene.

## 7. Provenance tiers (first reader only; second reader pending)

- Fig. 2 values: `direct_measurement`, figure-read, external standard; confidence **low** for Cys-Amadori (mass balance,
  standard), low-medium for cysteine.
- Table 1 / 2: `direct_measurement` (printed); semi-quantitative (RF = 1) for Table 1; confidence medium for CAMOLA shares.
- Sec. 4.4: [Z], `derived`.

## M4 role (case 8660c005)

**Role fixed before reading:** band source for the Cys-hexose Amadori steps (its no-aldehyde control), by amendment;
never decisive (BTBU's pots are D3's laboratory); measurement that would count: "a Cys-glucose Amadori time course with
charges and T stated"; P ~0.3.

- **The named measurement is met, with a defect:** the CG control's Cys-Amadori (and cysteine) time course with the
  charge (200 + 200 mM, 0.2 M phosphate pH 6.50) and temperature (145 °C) stated, triplicates, figure-only.
- **What it supplies.** Coordinates: "Glc + Cys → Cys-hexose Amadori, formation" and "Cys-hexose Amadori → 1-DG + Cys,
  decay" (sec. 5.2). Formation ≥ 62 % of cysteine to the Amadori compound within 10 min at 145 °C, pH 6.5; net post-peak
  loss 0.06-0.16 min⁻¹ (a lower bound on the decay constant), consistent with the band's 145 °C / pH 6.5 projection of
  Deng 2021 (0.15-0.31 min⁻¹) inside the ×0.3-3 width. Transfer: the same compound class (the cysteine-glucose Amadori
  compound itself), Hofmann's temperature, Du 2023's buffer and pH; **but** the cysteine mass balance fails by 40 % and
  the standard's identity is weakly documented, so it can at most centre a band, with its width carried as required
  agreement. It is D3's laboratory (CAAS), as the case anticipated, and probably the same control pot as `du2024`.
- **Decisive criterion:** **none met.** No NF and no MFT in the cysteine + glucose pot (criterion (a)); MFT appears
  only in the 2-heptenal pot (6.30 µg/L) without a label (criterion (c)); one pH (criterion (b)).
