# Zhao et al. 2019 (Food Chem) — EXTRACTION ([¹³C₅]-TTCA 30 mmol/L + glycine 30 + cysteine 30 + unlabelled Gly-Amadori 3 / 30 / 120 mmol/L, 0.2 M phosphate pH 5.5, 140 °C, 30 min; END-POINT volatiles as relative peak areas, ¹³C isotopomer split per compound)

### A carbon-tracing paper, not a kinetics paper: no TTCA, ARP or cysteine concentration at any time. Its one model-relevant number is that **97 % of the 2-methyl-3-furanthiol carries the intact TTCA xylose skeleton** when TTCA dominates the intermediate pool — confirming that TTCA's C5 is an MFT source, as `r_ttca_deg` → DPO → MFT assumes.

**Source on disk:** `data/articles/zhao2019.pdf` (42 pp., Elsevier **accepted manuscript**, FOCH 23438).
Read from the `pdftotext -layout` text layer; **Table 1 and Table 2 re-read from 130-170 dpi page
renders** because the manuscript typesets each number across two lines (e.g. "12.\n01ᶜ"), which the
text layer splits. Figures 1-4 are **pathway schemes only** (no data). **Table S1** (pH, A420 and
sensory scores of the six pots) is in the Supporting Information, **not on disk**.

## 0. Identity

| field | value |
|---|---|
| Title | "Meat flavor generation from different composition patterns of initial Maillard stage intermediates formed in heated cysteine-xylose-glycine reaction systems" |
| Authors | Jian Zhao, Tianze Wang (equal contribution), **Jianchun Xie** (corresponding, xjchun@th.btbu.edu.cn), Qunfei Xiao, Wenbin Du, Yaxin Wang, Jie Cheng, Shi Wang |
| Affiliations | BTBU Beijing (Xie group); Institute of Quality Standard and Testing Technology for Agro-products, CAAS |
| Venue | Food Chemistry (2019) — accepted manuscript; received 2 Apr 2018, revised 20 Aug 2018, accepted 21 Aug 2018 |
| DOI as printed | **10.1016/j.foodchem.2018.08.096** ✔ matches the brief |
| Funding | NSFC 31671895; Beijing Municipal NSF 6172004; National Key R&D 2017YFD0400106 |
| Same lab | `hou2017_extraction.md` (protocol parent: same buffer, temperature, time, vessel), `cao2017_extraction.md` (the 1 : 12 Cys : Gly optimum this paper re-tests) |

## 1. Why it matters (and why it mostly does not)

It re-tests Cao 2017's finding (meaty volatiles peak when Cys-Amadori ≈ Gly-Amadori) by charging
**purified intermediates** instead of free amino acids, and adds **[¹³C₅]-TTCA** so every volatile
splits into a TTCA-derived (labelled) and a Gly-Amadori-derived (unlabelled) fraction.

For this engine the paper has **no kinetic content**: one temperature, one time, no intermediate
measured after heating, volatiles as relative peak areas. What it does give:

1. **Carbon origin of MFT and FFT from TTCA.** In system (a), where TTCA outnumbers Gly-Amadori 10 : 1,
   **97.00 % of MFT, 93.00 % of FFT and 94.99 % of bis(2-methyl-3-furyl) disulfide are ¹³C-labelled**
   (Table 1), and in system (b) MFT is **65 % M+5 / 31 % M+0** (Table 2) — an intact C5 module. This
   is consistent with the engine's routing of TTCA's sugar through `r_ttca_deg` → DPO → NF/DDP → MFT
   and through PENT → TDP → FUR → FFT, and it **contradicts Xu 2013**'s claim that MFT and FFT carry
   no xylose label (`xu2013_extraction.md`).
2. **Evidence that TTCA releases xylose in the pot.** "[¹³C₅]-Gly-Amadori compounds were also found in
   the reaction solutions by LC-MS (data not shown)", which the authors attribute to "[¹³C₅]-xylose
   [that] was reversibly dissociated from [¹³C₅]-TTCA in the heated reaction systems" and then
   captured by glycine. That is a **direct tracer demonstration of the ring-opening `r_ttca_cys`
   (TTCA → Cys + pentose) at 140 °C, pH 5.5** — qualitative only, the amounts are not shown.
3. **Cysteine is sequestered by Gly-Amadori into "relatively stable thiazolidine compounds"** (after
   Ota 2006) — a cysteine sink the engine does not carry. Asserted, not measured.

### How it relates to the companion dossiers

| | `hou2017_extraction.md` | **this paper** | `cao2017_extraction.md` |
|---|---|---|---|
| charge | one intermediate + one amino acid, 30 mM each | TTCA + Gly + Cys 30 mM each, Gly-Amadori 3/30/120 mM | free Cys + Xyl 0.4 M, Gly 0-9.6 M |
| T, t | 140 °C, 5-30 min | 140 °C, 30 min | 90 °C, 0.5-5 h |
| tracer | none | **[¹³C₅]-TTCA** | none |
| intermediates measured | no | **no** (only a qualitative LC-MS sighting) | **yes, time-resolved** |
| output | peak areas, FD | peak areas, FD, isotopomer % | peak areas + concentrations |

## 2. Methods as they matter to a model

- **Intermediates** prepared "following the procedure described by Hou et al. (2017)": Gly-Amadori
  (N-(1-deoxy-D-xylulos-1-yl)-glycine), [¹²C₅]-TTCA and [¹³C₅]-TTCA (from [¹³C₅]-D-xylose 99 %,
  Sigma-Aldrich). **No purity stated in this paper.**
- **Charge.** "The dosage of [¹³C₅]-D-TTCA, [¹²C₅]-D-TTCA, glycine, or cysteine was **0.15 mmol**,
  while that of Gly-Amadori was varied as follows: **0.015, 0.15, and 0.60 mmol**", in **5 mL** sodium
  phosphate buffer (**pH 5.5, 0.2 M**) → TTCA **30 mmol/L**, glycine 30, cysteine 30, Gly-Amadori
  **3 / 30 / 120 mmol/L** in systems a / b / c. Controls a0, b0, c0 use [¹²C₅]-TTCA.
- **Heating.** **140 °C, 30 min, stirred**, Parallel Synthesis Poly-block 4 (HEL). Vessel by reference
  to Hou 2017: 15 mL pressure vial, 5 mL fill, sealed under air (headspace 10 mL).
- **Replicates.** "Two replicate experiments were performed for each reaction system"; RSD ≤ 13.61 %.
- **After heating:** pH and A420 measured (Table S1, not on disk).
- **Volatiles.** HS-SPME, 50/35 µm CAR/PDMS/DVB, **20 min at 50 °C** with stirring; "Immediately
  before desorption of the fiber, **0.1 µL of an internal standard (65 ng of 1,2-dichlorobenzene in
  1 µL of methanol) was injected**" — the IS goes onto the fibre, not into the sample (Flag 3). GC-MS
  on DB-Wax 30 m × 0.25 mm × 0.25 µm; AEDA/GC-O by three panellists, split 2 : 1 to 128 : 1; FD ≥ 128
  = potent.
- **Isotopomers.** Peak areas of [M] to [M+n] normalised, corrected for natural ¹³C 1.10 %, ³³S 0.76 %,
  ³⁴S 4.20 %, and for [M−1]⁺ hydrogen loss by the (M+1)/M⁺ ratio of the unlabelled control.
  "Sub-amount" = relative peak area × labelled (or unlabelled) fraction.

## 3. Tables re-typed

### Table 1 (verbatim). Identified volatile flavor compounds with amounts and FD values for the [¹³C₅]-D-TTCA/glycine/glycine-Amadori/cysteine reaction systems by GC-MS and GC-O

All amounts are **relative peak areas to 1,2-dichlorobenzene (dimensionless)**, mean of two
replicates, RSD ≤ 13.61 %; superscripts = Duncan groups within a row for one item; "-" not detected.
Systems a / b / c: [¹³C₅]-TTCA : Gly-Amadori = 10 : 1 / 1 : 1 / 1 : 4 by the Methods doses (the table
footnote prints "0.1:1, 1:1, and 1:4" — Flag 1). "% ¹³C" = total proportion of labelled
isotopomers. FD from AEDA on DB-Wax. ID: RI, MS, O (odour), S (authentic standard). **Bold** in the
source = potent meaty compounds used as EAS.

| RI | Compound | amount a | b | c | ¹³C sub a | b | c | unlab. sub a | b | c | % ¹³C a | b | c | odour | FD a | b | c | ID |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1015 | thiophene | 0.96ᶜ | 51.81ᵃ | 22.30ᵇ | 0.35 | 15.54 | - | 0.62ᶜ | 36.27ᵃ | 22.30ᵇ | 36.00 | 30.00 | - | garlic | 1 | 64 | 32 | RI,MS,O,S |
| 1093 | **2-methylthiophene** | 13.46ᶜ | 355.26ᵃ | 119.81ᵇ | 12.92ᶜ | 269.99ᵃ | 28.75ᵇ | 0.54ᵇ | 85.26ᵃ | 91.06ᵃ | 96.00 | 76.00 | 24.00 | sulfury | 64 | ≥128 | ≥128 | RI,MS,O,S |
| 1160 | 2,5-dimethylthiophene | 0.20 | 4.51 | - | 0.16 | 0.72 | - | 0.04 | 3.79 | - | 80.84 | 16.01 | - | sulfury | 1 | 16 | 1 | RI,MS,O,S |
| 1174 | 2-ethylthiophene | 0.02ᶜ | 4.06ᵃ | 0.64ᵇ | 0.01 | 0.77 | - | 0.01ᶜ | 3.29ᵃ | 0.64ᵇ | 58.33 | 19.00 | - | meaty | - | 4 | - | RI,MS,O,S |
| 1254 | thiazole | 0.30ᶜ | 7.26ᵃ | 4.69ᵇ | 0.23ᶜ | 2.11ᵃ | 0.47ᵇ | 0.06ᶜ | 5.16ᵃ | 4.22ᵇ | 79.01 | 29.00 | 10.01 | meaty | 1 | 16 | 16 | RI,MS,O,S |
| 1304 | **2-methyl-3-furanthiol** | 12.01ᶜ | 161.50ᵃ | 33.98ᵇ | 11.65ᵇ | 111.43ᵃ | 5.78ᶜ | 0.36ᶜ | 50.06ᵃ | 28.20ᵇ | 97.00 | 69.00 | 17.00 | meaty | ≥128 | ≥128 | ≥128 | RI,MS,O,S |
| 1347 | **3-mercapto-2-pentanone** | 6.83ᵇ | 14.57ᵃ | 2.55ᶜ | 6.15ᵇ | 12.82ᵃ | 1.17ᶜ | 0.68ᶜ | 1.75ᵃ | 1.38ᵇ | 90.01 | 88.00 | 45.99 | meaty | ≥128 | ≥128 | ≥128 | RI,MS,O,S |
| 1368 | 2,4,5-trimethylthiazole | 0.18ᶜ | 10.08ᵃ | 2.75ᵇ | 0.13ᵇ | 2.01ᵃ | 0.11ᵇ | 0.04ᶜ | 8.06ᵃ | 2.64ᵇ | 74.87 | 20.00 | 4.01 | sulfury | - | 16 | 4 | RI,MS,O,S |
| 1428 | **2-furfurylthiol** | 1.54ᶜ | 25.02ᵃ | 9.32ᵇ | 1.44ᶜ | 19.26ᵃ | 3.07ᵇ | 0.11ᵇ | 5.75ᵃ | 6.24ᵃ | 93.00 | 77.00 | 33.00 | coffee | ≥128 | ≥128 | ≥128 | RI,MS,O,S |
| 1438 | 5-ethyl-2,4-dimethylthiazole | 0.07ᶜ | 5.51ᵃ | 1.81ᵇ | 0.07 | 2.10 | - | - | 3.42ᵃ | 1.81ᵇ | 94.94 | 38.00 | - | - | - | - | - | RI,MS |
| 1512 | 2-methyltetrahydrothiophen-3-one | 0.19ᶜ | 34.16ᵇ | 69.88ᵃ | 0.17ᶜ | 9.22ᵃ | 6.29ᵇ | 0.02ᶜ | 24.94ᵇ | 63.59ᵃ | 88.04 | 27.00 | 9.00 | meaty | 1 | 64 | 64 | RI,MS,O,S |
| 1568 | **3-thiophenethiol** | 9.70ᶜ | 29.53ᵃ | 14.31ᵇ | 3.10ᵇ | 6.79ᵃ | 0.43ᶜ | 6.60ᶜ | 22.74ᵃ | 13.88ᵇ | 32.00 | 23.00 | 3.00 | meaty | ≥128 | ≥128 | ≥128 | RI,MS,O,S |
| 1573 | anti-3,5-dimethyl-1,2,4-trithiolane | 0.04ᶜ | 66.65ᵃ | 42.05ᵇ | - | - | - | 0.04ᶜ | 66.65ᵃ | 42.05ᵇ | - | - | - | meaty | - | ≥128 | ≥128 | RI,MS,O |
| 1594 | syn-3,5-dimethyl-1,2,4-trithiolane | 0.04ᶜ | 17.37ᵃ | 10.59ᵇ | - | - | - | 0.04ᶜ | 17.37ᵃ | 10.59ᵇ | - | - | - | meaty | - | 64 | 64 | RI,MS,O |
| 1631 | 2-acetylthiazole | 0.16ᶜ | 15.61ᵃ | 9.92ᵇ | 0.13ᶜ | 7.18ᵃ | 1.79ᵇ | 0.03ᵇ | 8.43ᵃ | 8.14ᵃ | 82.86 | 46.00 | 18.00 | roasted | 1 | 64 | 64 | RI,MS,O,S |
| 1675 | thiophene-2-carboxaldehyde | 1.19ᶜ | 9.55ᵇ | 27.48ᵃ | 0.97ᶜ | 4.78ᵇ | 6.05ᵃ | 0.21ᶜ | 4.78ᵇ | 21.43ᵃ | 81.99 | 50.00 | 22.00 | sulfury | 4 | 16 | 32 | RI,MS,O,S |
| 1754 | 1,2,3-trithiolane | 0.12ᶜ | 8.69ᵃ | 1.15ᵇ | - | - | - | 0.12ᶜ | 8.69ᵃ | 1.15ᵇ | - | - | - | sulfury | - | 16 | - | RI,MS,O |
| 1765 | 5-methylthiophene-2-carboxaldehyde | 0.02ᶜ | 4.13ᵃ | 1.31ᵇ | 0.01ᶜ | 1.16ᵃ | 0.16ᵇ | 0.01ᶜ | 2.97ᵃ | 1.15ᵇ | 38.46 | 27.99 | 11.97 | meaty | 1 | 16 | 8 | RI,MS,O,S |
| 1795 | 3-methylthiophene-2-carboxaldehyde | 0.53ᶜ | 28.73ᵇ | 43.04ᵃ | 0.45ᶜ | 10.34ᵃ | 6.03ᵇ | 0.08ᶜ | 18.39ᵇ | 37.01ᵃ | 85.03 | 36.00 | 14.00 | meaty | 2 | 64 | 64 | RI,MS,O,S |
| 1842 | thieno[3,2-b]thiophene | 3.06ᶜ | 137.80ᵃ | 37.21ᵇ | 1.16ᶜ | 27.56ᵃ | 3.72ᵇ | 1.90ᶜ | 110.24ᵃ | 33.49ᵇ | 38.01 | 20.00 | 10.00 | meaty | 16 | ≥128 | 64 | RI,MS,O |
| 1898 | thiophene-2,5-dicarboxaldehyde | 0.26ᶜ | 11.29ᵃ | 1.37ᵇ | 0.24ᵇ | 5.98ᵃ | 0.21ᵇ | 0.02ᶜ | 5.31ᵃ | 1.16ᵇ | 93.99 | 53.00 | 15.03 | sulfury | - | 16 | - | RI,MS,O |
| 1911 | 2-methylthieno[3,2-b]thiophene | 0.09ᶜ | 12.67ᵇ | 26.60ᵃ | 0.07ᵇ | 3.29ᵃ | 3.19ᵃ | 0.02ᶜ | 9.37ᵇ | 23.41ᵃ | 75.51 | 26.00 | 12.00 | meaty | 1 | 64 | 64 | RI,MS,O |
| 1932 | 2-methylthieno[2,3-b]thiophene | 0.02ᶜ | 4.10ᵇ | 10.00ᵃ | 0.02ᶜ | 1.64ᵃ | 1.20ᵇ | - | 2.46ᵇ | 8.80ᵃ | 80.00 | 40.00 | 12.00 | meaty | 1 | 64 | 64 | RI,MS,O |
| 1964 | thieno[2,3-b]thiophene | 0.09ᶜ | 5.71ᵇ | 9.22ᵃ | 0.07ᶜ | 2.63ᵃ | 2.12ᵇ | 0.01ᶜ | 3.09ᵇ | 7.10ᵃ | 84.38 | 45.99 | 23.00 | meaty | - | 64 | 64 | RI,MS,O |
| 2032 | 2-methylthieno[3,4-b]thiophene | 0.15ᶜ | 8.61ᵃ | 7.76ᵇ | 0.12ᶜ | 2.50ᵃ | 0.78ᵇ | 0.03ᶜ | 6.12ᵇ | 6.98ᵃ | 80.12 | 29.00 | 10.00 | - | - | - | - | RI,MS |
| 2407 | **bis(2-methyl-3-furyl) disulfide** | 0.64ᶜ | 26.22ᵇ | 44.64ᵃ | 0.61ᶜ | 22.54ᵃ | 9.37ᵇ | 0.03ᶜ | 3.68ᵇ | 35.26ᵃ | 94.99 | 86.00 | 21.00 | meaty | ≥128 | ≥128 | ≥128 | RI,MS,O,S |
| 2453 | 2-methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0.07 | - | - | 0.07 | - | - | - | - | - | 98.63 | - | - | - | - | - | - | RI,MS |
| | **Subtotal (TAS)** | 51.94 | 1060.42 | 554.35 | 40.31 | 542.38 | 80.68 | 11.62 | 518.03 | 473.68 | 77.61 | 51.15 | 14.55 | | | | | |
| 801 | furan | 1.57ᵇ | 68.66ᵃ | 72.95ᵃ | 1.42ᶜ | 33.65ᵃ | 9.48ᵇ | 0.16ᶜ | 35.02ᵇ | 63.47ᵃ | 89.98 | 49.00 | 13.00 | - | - | - | - | RI,MS |
| 873 | 2-methylfuran | 0.75ᶜ | 101.05ᵇ | 270.39ᵃ | 0.63ᶜ | 45.47ᵃ | 29.74ᵇ | 0.12ᶜ | 55.58ᵇ | 240.65ᵃ | 84.04 | 45.00 | 11.00 | - | - | - | - | RI,MS |
| 1456 | furfural | 0.31ᶜ | 28.76ᵇ | 492.13ᵃ | 0.26ᶜ | 19.27ᵇ | 118.11ᵃ | 0.05ᶜ | 9.49ᵇ | 374.02ᵃ | 83.93 | 67.00 | 24.00 | caramel | - | - | 64 | RI,MS,O,S |
| 1495 | 2-acetylfuran | - | 0.81 | 8.24 | - | 0.19 | 0.91 | - | 0.62 | 7.33 | - | 24.00 | 11.00 | - | - | - | - | RI,MS |
| 1602 | 2,2'-methylenedifuran | - | 17.60 | 47.34 | - | 9.15 | 3.79 | - | 8.45 | 43.56 | - | 52.00 | 8.00 | - | - | - | - | RI,MS |
| 1897 | 1-(2-furanylmethyl)-1H-pyrrole | - | 0.85 | 4.42 | - | 0.17 | 0.13 | - | 0.68 | 4.29 | - | 20.00 | 2.99 | - | - | - | - | RI,MS |
| | **Subtotal (other)** | 2.63 | 217.73 | 895.48 | 2.31 | 107.90 | 162.16 | 0.33 | 109.84 | 733.31 | 87.58 | 49.56 | 18.11 | | | | | |

Read notes: superscripts on a few sub-amount cells are illegible in the manuscript typesetting and
are omitted where not printed; the ID column "S" marks authentic-standard confirmation. Values were
cross-checked against Table 2 (e.g. MFT 69 % labelled in Table 1 vs 31 % M+0 / 65 % M+5 / 4 %
minor isotopomers in Table 2; furfural 67 % vs 33 % M+0).

### Table 2 (verbatim). Isotope distribution patterns of the identified flavor compounds in the reaction mixture of system (b)

Columns = number of ¹³C atoms in the isotopomer (%); corrected for natural ¹³C / ³³S / ³⁴S and for
[M−1]⁺; blank = not reported.

| No. | Compound | M⁺ (m/z) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2-methylfuran | 82 | 55 | <1 | <1 | <1 | 2 | 43 | | | | | |
| 2 | 2-methylthiophene | 98 | 24 | 1 | 2 | 3 | 2 | 67 | | | | | |
| 3 | 2-methyl-3-furanthiol | 114 | 31 | <1 | <1 | 1 | 3 | 65 | | | | | |
| 4 | 3-mercapto-2-pentanone | 118 | 12 | <1 | <1 | 1 | 5 | 82 | | | | | |
| 5 | 2-furfurylthiol | 114 | 23 | <1 | 2 | <1 | 3 | 72 | | | | | |
| 6 | furfural | 96 | 33 | <1 | <1 | 2 | <1 | 65 | | | | | |
| 7 | bis(2-methyl-3-furyl)-disulfide | 226 | 13 | 1 | <1 | <1 | 2 | 35 | 3 | 4 | <1 | 1 | 40 |
| 8 | furan | 68 | 51 | 3 | <1 | 4 | 41 | | | | | | |
| 9 | thiophene | 84 | 70 | 3 | 7 | 3 | 16 | | | | | | |
| 10 | 2-methyltetrahydrothiophen-3-one | 116 | 73 | 2 | 2 | 2 | 1 | 19 | | | | | |
| 11 | 3-methylthiophene-2-carboxaldehyde | 126 | 64 | 3 | 2 | 2 | 27 | <1 | 1 | | | | |
| 12 | 2-acetylthiazole | 127 | 54 | 2 | 2 | 39 | 2 | 1 | | | | | |
| 13 | 3-thiophenethiol | 116 | 77 | 3 | 13 | 3 | 2 | 2 | | | | | |
| 14 | 2-ethylthiophene | 112 | 81 | <1 | 3 | <1 | 15 | <1 | <1 | | | | |
| 15 | 2,5-dimethylthiophene | 112 | 84 | 2 | 1 | <1 | 13 | <1 | <1 | | | | |
| 16 | thiophene-2,5-dicarboxaldehyde | 140 | 47 | <1 | <1 | 3 | 4 | 46 | <1 | | | | |
| 17 | 2-methylthieno[2,3-b]thiophene | 154 | 60 | 3 | 2 | 4 | 2 | 26 | 2 | 1 | | | |
| 18 | 2-methylthieno[3,4-b]thiophene | 154 | 71 | 1 | <1 | 4 | 3 | 18 | 1 | <1 | | | |
| 19 | thieno[2,3-b]thiophene | 140 | 54 | <1 | 4 | 3 | 32 | 3 | 4 | | | | |
| 20 | thieno[3,2-b]thiophene | 140 | 80 | 3 | 14 | 1 | 2 | <1 | <1 | | | | |
| 21 | 5-methylthiophene-2-carboxaldehyde | 126 | 72 | 1 | 14 | 1 | 11 | <1 | <1 | | | | |
| 22 | 2-acetylfuran | 110 | 76 | <1 | <1 | 24 | <1 | <1 | <1 | | | | |
| 23 | thiophene-2-carboxaldehyde | 112 | 50 | 2 | 1 | 16 | 2 | 29 | | | | | |
| 24 | 2-methylthieno[3,2-b]thiophene | 154 | 74 | 2 | 2 | 12 | 1 | 8 | 1 | <1 | | | |
| 25 | 2,4,5-trimethylthiazole | 127 | 80 | 2 | 1 | <1 | 14 | 1 | 2 | | | | |
| 26 | 5-ethyl-2,4-dimethylthiazole | 141 | 62 | 3 | <1 | <1 | <1 | 35 | <1 | <1 | | | |
| 27 | thiazole | 85 | 71 | 26 | 4 | <1 | | | | | | | |
| 28 | 1-(2-furanylmethyl)-1H-pyrrole | 147 | 80 | <1 | <1 | <1 | 20 | <1 | <1 | <1 | <1 | <1 | |
| 29 | 2,2'-methylenedifuran | 148 | 48 | 4 | 4 | <1 | 2 | 42 | <1 | <1 | <1 | <1 | |

**Internal consistency checks (pass):** Table 2 (system b) has thieno[3,2-b]thiophene 80 % unlabelled and
2-methylthieno[3,2-b]thiophene 74 % unlabelled, i.e. 20 % / 26 % labelled — matching Table 1's
20.00 / 26.00 % ✔. Table 2's MFT 65 % M+5 + 3 % M+4 + 1 % M+3 ≈ 69 % matches Table 1's 69.00 % ✔.

### Arithmetic (mine)

**1. TTCA-derived MFT does not scale with Gly-Amadori the way a single-pathway model would.** The
labelled MFT sub-amount is 11.65 / 111.43 / 5.78 in a / b / c with TTCA, Cys and Gly held at 30 mM.
A 9.6× rise from a to b with only Gly-Amadori increased (3 → 30 mM), then a 19× fall from b to c
(30 → 120 mM). If TTCA → MFT ran only through its own C5 (TTCA → Cys-Amadori → 1-DX → MFT), the
labelled MFT should not depend on Gly-Amadori at all. The authors explain it by [¹³C₅]-xylose
released from TTCA forming [¹³C₅]-Gly-Amadori (so the label also flows through glycine), and by
Gly-Amadori consuming cysteine/H2S at high dose. Either way the **sulfur supply, not the TTCA carbon,
limits labelled MFT in system c** — a direction the engine can reproduce (Cys competition in
`r_dpo_ddp`, H2S as the MFT co-reactant) but not test quantitatively here.

**2. Ratios, not levels.** Relative peak areas with the IS added to the fibre are not concentrations.
Only within-row ratios (b/a, c/b) and labelled fractions are transferable.

## 4. Kinetic numbers the repository can use

**There is no rate constant, no activation energy and no time series in this paper.**

| step / observable | quantity | value | unit | conditions | order | source | evidence class |
|---|---|---|---|---|---|---|---|
| carbon origin of MFT from TTCA | labelled fraction | **97.00 / 69.00 / 17.00** | % | a / b / c, 140 °C 30 min pH 5.5 | — | Table 1 | **within_study_ratio** (isotope fraction) |
| carbon origin of FFT | labelled fraction | **93.00 / 77.00 / 33.00** | % | same | — | Table 1 | within_study_ratio |
| carbon origin of furfural | labelled fraction | 83.93 / 67.00 / 24.00 | % | same | — | Table 1 | within_study_ratio |
| MFT intact-C5 module | M+5 / M+0 | **65 / 31** | % | system b | — | Table 2 | within_study_ratio |
| 3-mercapto-2-pentanone intact-C5 | M+5 / M+0 | 82 / 12 | % | system b | — | Table 2 | within_study_ratio |
| 2-acetylthiazole C3 module | M+3 / M+0 | 39 / 54 | % | system b | — | Table 2 | within_study_ratio — consistent with MGO (C3 from xylose) + Cys, i.e. `r_cys_actz` |
| TTCA ring-opening releases xylose at 140 °C | [¹³C₅]-Gly-Amadori detected | present | — | all systems | — | text, "data not shown" | **threshold** (presence only) |
| every volatile amount | relative peak areas | Table 1 | dimensionless | — | — | Table 1 | **level_only** |
| pH, A420, sensory after heating | — | — | — | — | — | Table S1 | **not on disk** |

## 5. What this can and cannot anchor in the model

**Can:** (i) a **qualitative topology check** that TTCA's intact C5 reaches MFT, FFT,
3-mercapto-2-pentanone and furfural (all ≥ 65 % M+5 in system b, ≥ 90 % labelled in system a);
(ii) a **qualitative check that `r_cys_actz` takes a C3 dicarbonyl** (2-acetylthiazole 39 % M+3);
(iii) a qualitative tracer demonstration that TTCA releases free xylose at 140 °C / pH 5.5 (the
`r_ttca_cys` direction).

**Cannot:** anything kinetic or absolute. No TTCA, Cys-Amadori, cysteine or H2S is quantified; one
time point; relative peak areas; IS on the fibre; n = 2. Glycine and Gly-Amadori have no engine
species, so systems b and c are outside the module's scope; system a (Gly-Amadori 3 mM against
TTCA 30 mM) is the closest to a TTCA + Cys + Gly pot but still carries glycine.

## 6. Flags

1. **The ratio label is inverted in the Table 1 footnote.** Methods: Gly-Amadori 0.015 / 0.15 /
   0.60 mmol against 0.15 mmol TTCA → TTCA : Gly-Amadori = **10 : 1 / 1 : 1 / 1 : 4**. Footnote: "the
   molar ratios between [¹³C₅]-D-TTCA and Gly-Amadori were **0.1:1**, 1:1, and 1:4". The text
   ("Gly-Amadori relative to that of TTCA was varied at the ratio of 0.1:1, 1:1, or 4:1") states it
   as Gly-Amadori : TTCA. The doses are unambiguous; the footnote mixes the two orientations.
2. **No intermediate quantified.** The title is about "composition patterns of initial Maillard stage
   intermediates", but those patterns are *imposed* by the charge, not measured after heating. The
   only post-heating intermediate observation ([¹³C₅]-Gly-Amadori by LC-MS) is "data not shown".
3. **Internal standard added to the SPME fibre**, not to the sample: it normalises injection and MS
   response but not headspace partitioning or matrix effects. Relative areas across systems with
   very different browning and pH are therefore not proportional to concentrations.
4. **n = 2** with Duncan's multiple-range test on two replicates — the significance letters carry
   little weight.
5. **Accepted manuscript**; volume/pages not printed on disk.
6. **Purity of the [¹³C₅]-TTCA and Gly-Amadori not stated** in this paper (deferred to Hou 2017).
   Unlabelled carbon in system a's MFT (3 %) is within the plausible range of label impurity + the
   3 mM Gly-Amadori, so system a cannot distinguish a small non-TTCA MFT route.
7. **Registry gaps:** TTCA, Gly-Amadori, glycine, cysteine absent; `2_methyl_3_furanthiol`,
   `2_furfurylthiol`, `2_methylthiophene`, `furfural` present. Paper not in `papers.yml`.
8. **What to request:** Table S1 (pH and A420 after heating — would give a measured pH for each
   pot), the LC-MS data for [¹³C₅]-Gly-Amadori and residual TTCA, and the "0.90 mmol" and "8 : 1"
   arms mentioned as "data not shown".

---

## Supplementary material (read 2026-10-29)

**No supplementary material exists to read.** The owner searched the publisher's page for Zhao et al. 2019 Food Chem
(10.1016/j.foodchem.2018.08.096) on 2026-10-29 and found no supplementary file. The main text nevertheless cites a
**"Table S1"** four times (sec. 2: the six systems' compositions; sec. 3.1: compositions, A420, pH; sec. 3.3: sensory) — and
this dossier (§2, §4, §8) and the sulfur ledger (`zhao2019_tables1_ph_a420`) recorded it as an on-line table "not on disk".
It should now be treated as **not available at all**.

What the main text supports instead (sec. 3.1, verbatim): *"Table S1 illustrates that from systems (a) to (c), greater dosage
of Gly-Amadori compounds used resulted in greater absorbance at 420 nm and drop in the values of pH of the reaction
solution"*, and the sensory score *"first increased and then went down with the increase in the Gly-Amadori concentration"*.
**Direction only, no values**: the ledger status moves from `supplement_missing` to `text_only` with the flag
`si_not_found_at_publisher_2026-10-29` (change sent as a patch, since sulfur.yml was being edited by another agent).
The system compositions are fully given in the main-text methods (§2 here), so nothing else depends on Table S1.
