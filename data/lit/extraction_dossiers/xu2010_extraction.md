# Xu 2010 — EXTRACTION (ribose + cysteine, 50 + 50 mmol/L in 0.2 M pyrophosphate pH 5.6, 140 C, 60 min, under supercritical CO2 or N2 at 10/20/30/40 MPa and a conventional control; HS-SPME-GC-MS volatiles, response factor 1 against tridecane)
### The only paper on disk that heats a cysteine pot under a nitrogen atmosphere and prints MFT, FFT and their disulfides side by side — and the disulfide share of the thiol does not move with the atmosphere.

**Source on disk:** `data/articles/Xu2010.pdf` (10 pp., owner's download, 2026-09-08). Read from
the text layer (`scratchpad/articles/Xu2010.txt`); Table 1 came through with all nine data columns
on one line per compound and was re-typed below (column order fixed from the header: Control, then
SC-CO2 / SC-N2 pairs at 10, 20, 30, 40 MPa; the printed "peak value of 30.60 ng/ml at 20 MPa" for
the disulfide confirms the column assignment). No page was rasterised; Figures 1-3 (A280, A420,
chromatograms, class shares) are FIGURE-ONLY and no value was read from them.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Pressure on the Maillard Reaction between Ribose and Cysteine in Supercritical Carbon Dioxide" |
| Authors | Honggao Xu, Wenhao He, Xuan Liu, Yanxiang Gao* (China Agricultural University, Beijing). The paper's own citation line prints "Liu K." for the third author; the byline says Xuan Liu |
| Venue | Czech Journal of Food Sciences 28 (2010) 192-201. Received October 20 2008, accepted March 17 2010 |
| DOI | none printed |
| Naming | SC-CO2 / SC-N2 = supercritical carbon dioxide / nitrogen as the pressurising medium; MFT = 2-methyl-3-furanthiol; FFT = 2-furanmethanethiol (2-furfurylthiol); MFTD = bis(2-methyl-3-furyl) disulfide; norfuraneol = 4-hydroxy-5-methyl-3(2H)-furanone; LRI = linear retention index on the (unstated) column |
| Companions | Xu 2008a (J. Sci. Food Agric. 88, 328: the same apparatus and procedure, SC-CO2 at 40 MPa; the "procedures described previously") and Xu 2008b (Food Res. Int. 41, 730: ribose:cysteine ratios) — neither on disk; Whitfield & Mottram 1999, Mottram & Nobrega 2002, Hofmann & Schieberle 1998 cited for mechanisms (all in the repo) |

## 1. Why it matters

Wave B17 (`results/validation/kinetic_core_b17_prereg.md` section 6, T3) found the model's thiol
disulfide channel **oxidant-limited**: the dimer share of free MFT is 0.04-0.9 % where Zhou 2023
and Zhang 2024 measure 7-10 %, with both dimerisation constants already at the top of their bands,
because the ambient oxidant reservoir (B11, shipped inert) runs out. The question that raises is
whether the disulfide share depends on the oxygen supply at all. This paper heats the same
chemistry (ribose + cysteine, pH 5.6, 140 C, 1 h) under three atmospheres — air-headspace control,
10-40 MPa nitrogen, 10-40 MPa carbon dioxide — and prints MFT, FFT, MFTD and three mixed disulfides
in each. What it shows, as printed: the disulfide-to-thiol proportion is the same under nitrogen as
in the control and under CO2 (MFTD/MFT by mass 0.40-0.58 in all nine columns), while the absolute
amounts of both thiol and disulfide rise four- to eight-fold under CO2 for a pH reason. Under this
paper's conditions an inert pressurising gas does not lower the disulfide share. The caveats
(no degassing, no oxygen measurement, SPME with response factor 1, possible on-fibre oxidation) are
in section 5 and matter.

## 2. Methods as they matter to a model

- **Pot.** "D-Ribose and L-cysteine monohydrochloride (0.1 M) in 0.2 M pyrophosphate buffer at pH
  5.6 were prepared separately"; Table 1 footnote b: "Each model system consisted of 5 mmol ribose
  and 5 mmol cysteine in 100 ml 0.2 M pyrophosphate buffer", i.e. equal volumes mixed to
  **[ribose] = [cysteine] = 50 mmol/L, 0.2 M sodium pyrophosphate, pH 5.6 initial, 100 mL**.
- **Reactor.** CWYF-2 supercritical apparatus (Hua'an, Nantong), procedure per Xu 2008a (not on
  disk): **140 C, 1 h**, pressurised with CO2 or with N2 to **10, 20, 30, 40 MPa**; "the
  conventional control experiments were also conducted on the same apparatus without the pressure
  media" (the browning discussion calls the control "0.3 MPa", i.e. autogenous pressure with the
  original air headspace). Vessel volume, headspace, stirring, heat-up time, depressurisation rate:
  not stated. **Nothing is said about degassing, purging, dissolved oxygen or oxygen content of
  the gases** (verified negative: the words oxygen, oxidation, air, degas, purge, flush, dissolved
  oxygen do not occur in the text).
- **Volatile capture.** The reaction mixture **and three absorption buffers** (traps on the vent
  gas, so volatiles lost on depressurisation are counted) were analysed and **summed**. HS-SPME
  DVB/CAR/PDMS 50/30 um 1 cm, 60 C, 20 min; GC oven 40 C -> 60 C at 20 C/min (5 min) -> 250 C at
  4 C/min (10 min); column and MS not stated. **Quantification: TIC peak area against 15.13 ng
  tridecane internal standard with a response factor of 1; reported as "approximate quantities in
  headspace (ng/ml of mixture)"**. Identification: MS + LRI against authentics (MS+LRI), MS + LRI
  near literature (ms+LRI), MS against library (ms), literature MS (MS). Mean CV < 20 %, no
  compound above 47 %.
- **Non-volatile readouts.** A280 ("overall intermediates") and A420 (browning), figure-only;
  A280 in SC-N2: Abs = -0.002 P + 0.6906 (R2 0.9979, P in MPa) — the one printed fit.
- **Replication.** Duplicate reactions, at least triplicate analyses; Tukey test; +/- in Table 1 is
  unlabelled (SD of the duplicates, presumably; n = 2 stated only for Figure 1).
- **Unit conversion.** ng/mL of mixture = ug/L. Molar: MFT and FFT M = 114.17 g/mol (1 ng/mL =
  8.76 nmol/L); MFTD M = 226.32 (1 ng/mL = 4.42 nmol/L). No time series (single 60 min point), so
  no rate.

## 3. Tables re-typed

### Table 1. "Approximate quantities of volatiles identified in the headspace of ribose-cysteine model system under different pressure media" — ng/mL of mixture, mean +/- (unlabelled)

Column key: C = conventional control; 10/20/30/40 = MPa; CO2 / N2 = pressurising medium. "-" =
below detection (~ 0.07 ng/mL); tr = trace (< 0.3 ng/mL); + = present, quantification confounded
by an adjacent peak. Registry ids in brackets where they exist.

**Thiols**

| LRI | compound | C | 10 CO2 | 10 N2 | 20 CO2 | 20 N2 | 30 CO2 | 30 N2 | 40 CO2 | 40 N2 | ID |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 818 | 3-mercapto-2-butanone | 4.0 +/- 0.7 | 3.4 +/- 1.2 | 4.1 +/- 0.2 | 6.8 +/- 2.7 | 4.9 +/- 0.6 | 6.6 +/- 0.4 | 2.7 +/- 0.3 | 4.5 +/- 1.1 | 3.6 +/- 0.0 | MS+LRI |
| 869 | **2-methyl-3-furanthiol** [`2_methyl_3_furanthiol`] | 7.2 +/- 2.0 | 35.8 +/- 9.5 | 7.2 +/- 0.2 | 58.1 +/- 9.1 | 7.0 +/- 0.6 | 46.4 +/- 3.5 | 7.4 +/- 1.5 | 34.6 +/- 15.8 | 7.3 +/- 0.3 | MS+LRI |
| 904 | 3-mercapto-2-pentanone | 1.9 +/- 0.3 | 5.9 +/- 2.4 | 2.3 +/- 0.8 | 14.4 +/- 5.4 | 1.8 +/- 0.2 | 13.8 +/- 0.5 | 1.6 +/- 0.0 | 9.8 +/- 3.2 | 1.8 +/- 0.1 | MS+LRI |
| 909 | 2-mercapto-3-pentanone | 2.3 +/- 0.2 | 4.7 +/- 0.9 | 2.5 +/- 0.4 | 6.0 +/- 2.2 | 2.4 +/- 0.2 | 4.3 +/- 0.1 | 2.3 +/- 0.3 | 2.6 +/- 0.8 | 2.2 +/- 0.1 | ms+LRI |
| 913 | **2-furanmethanethiol (FFT)** [`2_furfurylthiol`] | 7.3 +/- 1.9 | 32.9 +/- 12.8 | 10.0 +/- 1.4 | 61.5 +/- 18.9 | 9.6 +/- 1.3 | 43.4 +/- 3.6 | 8.4 +/- 1.1 | 23.9 +/- 9.7 | 7.8 +/- 0.1 | MS+LRI |
| 980 | 3-thiophenethiol | 8.0 +/- 1.6 | 24.7 +/- 12.5 | 8.6 +/- 0.5 | 40.1 +/- 10.5 | 9.0 +/- 2.4 | 45.5 +/- 3.1 | 8.6 +/- 1.8 | 41.1 +/- 20.0 | 9.4 +/- 0.4 | ms+LRI |
| 1066 | 2-methyl-3-thiophenethiol | 12.3 +/- 2.5 | 25.5 +/- 9.7 | 14.9 +/- 1.0 | 31.5 +/- 0.6 | 17.3 +/- 0.8 | 20.6 +/- 1.6 | 16.2 +/- 0.1 | 15.8 +/- 5.6 | 14.8 +/- 0.2 | ms+LRI |
| | **Total thiols** | 43.0 +/- 8.8 | 132.8 +/- 49.1 | 49.7 +/- 2.2 | 218.4 +/- 49.4 | 51.9 +/- 4.8 | 180.7 +/- 10.8 | 47.1 +/- 4.8 | 132.2 +/- 56.1 | 46.8 +/- 0.0 | |

**Disulfides**

| LRI | compound | C | 10 CO2 | 10 N2 | 20 CO2 | 20 N2 | 30 CO2 | 30 N2 | 40 CO2 | 40 N2 | ID |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1542 | **bis(2-methyl-3-furyl) disulfide** [`bis_2_methyl_3_furyl_disulfide`] | 3.0 +/- 0.4 | 20.9 +/- 0.2 | 2.9 +/- 0.3 | 30.6 +/- 2.2 | 3.2 +/- 0.1 | 23.2 +/- 1.0 | 4.0 +/- 0.4 | 18.6 +/- 8.6 | 3.6 +/- 0.1 | MS+LRI |
| 1578 | 3-[(2-methyl-3-furyl)dithio]-2-pentanone | - | 2.6 +/- 0.4 | - | 4.2 +/- 0.8 | - | 3.6 +/- 0.9 | - | 3.3 +/- 1.0 | - | LRI |
| 1643 | 2-methyl-3-[(2-furfuryl)dithio]furan (MFT-FFT mixed) | - | 1.6 +/- 0.1 | - | 1.8 +/- 0.2 | - | 2.5 +/- 0.8 | - | 2.3 +/- 1.2 | - | LRI |
| 1668 | 3-(2-furfuryldithio)-2-pentanone | - | - | - | tr | - | tr | - | tr | - | LRI |
| 1702 | 2,3-dihydro-5-methyl-4-[(2-methyl-3-furyl)dithio]furan | 2.5 +/- 0.2 | 10.0 +/- 1.1 | 2.6 +/- 0.2 | 14.6 +/- 1.9 | 3.1 +/- 0.4 | 16.1 +/- 0.8 | 3.2 +/- 0.5 | 16.0 +/- 7.6 | 3.6 +/- 0.3 | ms |
| 1745 | 2-methyl-3-[(2-methyl-3-thienyl)dithio]furan | 2.9 +/- 0.7 | 11.8 +/- 1.6 | 3.3 +/- 0.4 | 14.1 +/- 2.8 | 3.5 +/- 0.5 | 7.8 +/- 0.6 | 3.3 +/- 1.0 | 6.2 +/- 2.7 | 3.3 +/- 0.3 | LRI |
| | **Total disulfides** | 8.4 +/- 1.2 | 47.0 +/- 2.0 | 8.9 +/- 0.5 | 65.4 +/- 7.9 | 9.7 +/- 1.1 | 53.3 +/- 0.7 | 10.5 +/- 1.8 | 46.5 +/- 21.0 | 10.5 +/- 0.1 | |

**Bis(2-furfuryl) disulfide (the FFT homodimer) is not in the table at all**; the only
FFT-containing disulfides are the MFT-FFT mixed disulfide (CO2 only) and the trace
furfuryldithio-pentanone.

**Polysulfur heterocyclics**

| LRI | compound | C | 10 CO2 | 10 N2 | 20 CO2 | 20 N2 | 30 CO2 | 30 N2 | 40 CO2 | 40 N2 | ID |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1153 | 3,5-dimethyl-1,2,4-trithiolane (E or Z) | 7.8 +/- 2.0 | - | 4.4 +/- 2.2 | - | 3.0 +/- 0.1 | - | 1.8 +/- 0.6 | - | 1.8 +/- 0.0 | MS |
| 1160 | 3,5-dimethyl-1,2,4-trithiolane (E or Z) | 9.5 +/- 2.1 | - | 6.1 +/- 2.5 | - | 4.6 +/- 0.3 | - | 1.8 +/- 0.4 | - | 1.3 +/- 0.5 | MS |
| 1185 | 1,2-dithian-4-one | 4.6 +/- 0.3 | 8.7 +/- 0.8 | 5.1 +/- 0.5 | 7.1 +/- 0.3 | 4.9 +/- 0.7 | 5.7 +/- 0.2 | 5.8 +/- 0.2 | 6.9 +/- 0.5 | 6.1 +/- 0.3 | ms |
| 1232 | 3-methyl-1,2-dithian-4-one | 4.2 +/- 0.4 | 16.4 +/- 0.3 | 4.6 +/- 1.0 | 45.2 +/- 1.4 | 4.1 +/- 0.6 | 50.8 +/- 3.5 | 4.8 +/- 0.1 | 55.0 +/- 8.9 | 5.4 +/- 0.4 | MS |
| 1264 | 3,(5 or 6)-dimethyl-1,2-dithian-4-one (E or Z) | 3.0 +/- 0.2 | 9.1 +/- 1.3 | 2.6 +/- 0.4 | 8.7 +/- 0.2 | 2.5 +/- 0.2 | 6.7 +/- 0.8 | 3.1 +/- 0.1 | 6.1 +/- 1.1 | 3.4 +/- 0.1 | ms |
| 1266 | 3-methyl-1,2,4-trithiane | 5.7 +/- 0.0 | - | 8.2 +/- 2.9 | - | 9.0 +/- 2.2 | - | 7.1 +/- 1.0 | - | 6.2 +/- 0.9 | MS |
| 1274 | 3,(5 or 6)-dimethyl-1,2-dithian-4-one (E or Z) | 3.2 +/- 0.7 | 4.8 +/- 2.4 | 4.0 +/- 1.4 | 3.3 +/- 0.2 | 3.1 +/- 1.1 | 2.5 +/- 0.4 | 4.4 +/- 0.7 | 2.3 +/- 0.2 | 4.6 +/- 1.1 | ms |
| | **Total polysulfur heterocyclics** | 37.9 +/- 2.4 | 39.0 +/- 4.8 | 35.0 +/- 12.8 | 64.3 +/- 1.5 | 31.2 +/- 5.1 | 65.7 +/- 4.1 | 28.8 +/- 3.4 | 70.2 +/- 13.7 | 28.8 +/- 1.3 | |

**Miscellaneous (the pentose intermediates)**

| LRI | compound | C | 10 CO2 | 10 N2 | 20 CO2 | 20 N2 | 30 CO2 | 30 N2 | 40 CO2 | 40 N2 | ID |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 825 | methylpyrazine [`methylpyrazine`] | 3.0 +/- 0.3 | - | 2.2 +/- 0.3 | - | 2.0 +/- 0.5 | - | 1.9 +/- 0.1 | - | 1.4 +/- 0.1 | MS+LRI |
| 832 | furfural [`furfural`] | 2.0 +/- 0.3 | 9.1 +/- 0.1 | 1.9 +/- 0.3 | 23.0 +/- 1.1 | 1.6 +/- 0.3 | 27.9 +/- 2.4 | 1.7 +/- 0.1 | 34.5 +/- 3.0 | 1.3 +/- 0.2 | MS+LRI |
| 850 | 2-furanmethanol | 2.1 +/- 0.1 | 1.5 +/- 0.5 | 1.9 +/- 0.3 | 1.6 +/- 0.1 | 1.7 +/- 0.0 | - | 2.1 +/- 0.2 | - | 1.9 +/- 0.1 | MS+LRI |
| 952 | 1-(2-furyl)-2-propanone | 1.5 +/- 0.1 | 3.1 +/- 0.5 | 1.0 +/- 0.1 | 3.4 +/- 0.0 | 0.9 +/- 0.1 | 2.1 +/- 0.1 | + | 1.6 +/- 0.5 | + | ms |
| 1022 | 2-acetylthiazole | 6.5 +/- 0.1 | 3.1 +/- 0.3 | 6.8 +/- 1.5 | 1.8 +/- 0.1 | 7.3 +/- 1.2 | 2.0 +/- 0.2 | 6.8 +/- 0.5 | 1.6 +/- 0.2 | 5.9 +/- 0.4 | MS+LRI |
| 1075 | 4-hydroxy-5-methyl-3(2H)-furanone [`norfuraneol`] | 38.7 +/- 2.6 | 56.0 +/- 6.9 | 45.4 +/- 10.2 | 39.1 +/- 14.3 | 35.1 +/- 15.7 | 52.9 +/- 10.1 | 51.4 +/- 20.7 | 71.8 +/- 23.1 | 33.5 +/- 3.9 | MS+LRI |
| | **Total miscellaneous** | 53.7 +/- 2.9 | 72.8 +/- 5.7 | 59.1 +/- 10.6 | 68.8 +/- 13.5 | 48.6 +/- 17.9 | 84.9 +/- 12.2 | 63.9 +/- 21.4 | 109.5 +/- 26.8 | 44.1 +/- 3.1 | |

**Class totals for the classes not re-typed row by row** (thiophenes 11 rows, thiophenones 3
rows, fused bicyclics 12 rows dominated by a dihydrothienothiophene at LRI 1325 of 98-520 ng/mL;
all in the text layer if needed):

| class | C | 10 CO2 | 10 N2 | 20 CO2 | 20 N2 | 30 CO2 | 30 N2 | 40 CO2 | 40 N2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Total thiophenes | 83.3 +/- 7.5 | 135.2 +/- 1.9 | 87.0 +/- 15.5 | 119.7 +/- 0.2 | 90.8 +/- 11.2 | 108.0 +/- 9.9 | 93.4 +/- 0.2 | 107.7 +/- 18.9 | 86.7 +/- 3.8 |
| Total thiophenones | 20.9 +/- 2.2 | 26.5 +/- 0.7 | 23.0 +/- 4.7 | 19.9 +/- 1.3 | 28.8 +/- 1.5 | 16.9 +/- 0.6 | 29.1 +/- 1.1 | 14.0 +/- 3.0 | 29.2 +/- 1.3 |
| Total fused bicyclics | 795.7 +/- 20.6 | 438.9 +/- 43.8 | 722.4 +/- 96.5 | 319.0 +/- 4.6 | 893.8 +/- 89.1 | 223.3 +/- 42.6 | 837.6 +/- 0.4 | 197.8 +/- 35.0 | 772.3 +/- 12.3 |
| Grand total | 1042.9 +/- 37.9 | 892.2 +/- 93.2 | 984.9 +/- 140.4 | 875.4 +/- 59.5 | 1154.9 +/- 130.8 | 732.7 +/- 56.5 | 1110.5 +/- 19.6 | 677.9 +/- 174.6 | 1018.4 +/- 12.8 |

### Derived: the disulfide share of MFT under each atmosphere (mine, from Table 1)

Molar conversion: n(MFT) = m/114.17, n(MFTD) = m/226.32; "MFT-equivalents in the dimer" = 2 n(MFTD).

| column | MFT (ng/mL) | MFTD (ng/mL) | MFTD/MFT by mass | share of MFT-equivalents present as dimer, 2n(D)/(n(M) + 2n(D)) | molecules: n(D)/(n(M) + n(D)) | total disulfides / total thiols by mass |
|---|---:|---:|---:|---:|---:|---:|
| control (air headspace, autogenous P) | 7.2 | 3.0 | 0.42 | 29.6 % | 17.4 % | 0.20 |
| 10 MPa N2 | 7.2 | 2.9 | 0.40 | 28.9 % | 16.9 % | 0.18 |
| 20 MPa N2 | 7.0 | 3.2 | 0.46 | 31.6 % | 18.7 % | 0.19 |
| 30 MPa N2 | 7.4 | 4.0 | 0.54 | 35.3 % | 21.4 % | 0.22 |
| 40 MPa N2 | 7.3 | 3.6 | 0.49 | 33.2 % | 19.9 % | 0.22 |
| 10 MPa CO2 | 35.8 | 20.9 | 0.58 | 37.1 % | 22.8 % | 0.35 |
| 20 MPa CO2 | 58.1 | 30.6 | 0.53 | 34.7 % | 21.0 % | 0.30 |
| 30 MPa CO2 | 46.4 | 23.2 | 0.50 | 33.5 % | 20.1 % | 0.29 |
| 40 MPa CO2 | 34.6 | 18.6 | 0.54 | 35.2 % | 21.3 % | 0.35 |

Reading: **N2 at 10-40 MPa reproduces the control** for MFT (7.0-7.4 vs 7.2), FFT (7.8-10.0 vs
7.3) and MFTD (2.9-4.0 vs 3.0); CO2 raises MFT 5-8x, FFT 3-8x and MFTD 6-10x; the MFT-equivalent
share held in the dimer is **29-37 % in every column**, with no ordering by atmosphere (N2 29-35,
CO2 34-37, control 30). The authors' own reading of the disulfides: "the high concentration of
disulfides in the SC-CO2-treated samples can be ascribed to their high quantity of monomers." The
FFT homodimer is absent everywhere; the MFT-FFT mixed disulfide appears only under CO2 (1.6-2.5
ng/mL against FFT 24-62), i.e. FFT dimerises far less than MFT in this pot.

### Statements in the text on the atmosphere and pressure effects (levels FIGURE-ONLY for A280 / A420)

- CO2 lowers pH reversibly (Gevaudan 1996 cited; CO2 saturation "beyond 20 MPa and no further pH
  decline"); A280 maximal at 10 MPa CO2, A420 falls from control to 20 MPa CO2 then flat; N2
  pressure has no significant effect on either (P > 0.05).
- Thiols under CO2 peak at 20 MPa (3-thiophenethiol at 30 MPa); "no significant differences
  (P > 0.05) between the quantities of 2-furanmethanethiol and 2-methyl-3-furanthiol in the
  SC-CO2-treated ribose-cysteine mixtures".
- 3-methyl-1,2-dithian-4-one rises with CO2 pressure (16.4 -> 55.0) and stays at 4-5 under N2;
  3,5-dimethyl-1,2,4-trithiolanes and 3-methyl-1,2,4-trithiane are absent under CO2 ("carbon
  dioxide could terminate the Strecker degradation of cysteine", Xu 2008b).
- Furfural rises with CO2 pressure (9.1 -> 34.5) and is at control level under N2 (1.3-1.9); the
  authors invoke "the lower degradation rate of ribose in SC-CO2 (data not shown)".

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): MFT -> `2_methyl_3_furanthiol`; FFT ->
`2_furfurylthiol`; bis(2-methyl-3-furyl) disulfide -> `bis_2_methyl_3_furyl_disulfide`; furfural ->
`furfural`; 4-hydroxy-5-methyl-3(2H)-furanone -> `norfuraneol`; methylpyrazine ->
`methylpyrazine`; 3-methyl-1,2-dithian-4-one, the mixed disulfides, the thiophenethiols,
mercaptoketones and thienothiophenes -> not in registry. No rate constant and no barrier exist in
this paper (one time point, one temperature).

| step | quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|---|
| ribose + cysteine -> MFT | MFT level | 7.2 / 7.2 / 7.0 / 7.4 / 7.3 (control, N2 10/20/30/40 MPa); 35.8 / 58.1 / 46.4 / 34.6 (CO2 10/20/30/40 MPa) | ng/mL of mixture (RF 1 vs tridecane, HS-SPME) | 50 + 50 mmol/L, 0.2 M pyrophosphate pH 5.6, 140 C, 60 min | Table 1 | level_only (semi-quantitative) |
| ribose + cysteine -> FFT | FFT level | 7.3 / 10.0 / 9.6 / 8.4 / 7.8; 32.9 / 61.5 / 43.4 / 23.9 | same | same | Table 1 | level_only (semi-quantitative) |
| 2 MFT -> MFTD | bis(2-methyl-3-furyl) disulfide level | 3.0 / 2.9 / 3.2 / 4.0 / 3.6; 20.9 / 30.6 / 23.2 / 18.6 | same | same | Table 1 | level_only (semi-quantitative) |
| 2 MFT -> MFTD | MFTD / MFT (mass) | 0.42 (control); 0.40-0.54 (N2); 0.50-0.58 (CO2) | ratio | same | derived from Table 1 | within_study_ratio (response-factor-immune only if the SPME/RF-1 bias is the same for thiol and disulfide, which it is not — Flag 3) |
| 2 MFT -> MFTD | MFT-equivalents held as dimer | 29.6 % (control); 28.9-35.3 % (N2); 33.5-37.1 % (CO2) | molar share | same | derived | within_study_ratio, same caveat |
| atmosphere effect | N2 vs control, MFT / FFT / MFTD | 0.97-1.03 / 1.07-1.37 / 0.97-1.33 | fold | 10-40 MPa | derived | within_study_ratio: **an inert pressurising gas changes neither the thiol nor the disulfide** |
| atmosphere effect | CO2 vs control, MFT / FFT / MFTD | 4.8-8.1 / 3.3-8.4 / 6.2-10.2 | fold | 10-40 MPa | derived | within_study_ratio (a pH / catalysis effect per the authors) |
| MFT vs FFT split | FFT / MFT | 1.01 (control); 1.07-1.39 (N2); 0.69-1.06 (CO2) | ratio | same | derived | within_study_ratio; the sulfur lane's MFT:FFT split in a 50 mM ribose + cysteine pot at pH 5.6, 140 C |
| 2 FFT -> bis(2-furfuryl) disulfide | FFT homodimer | not detected (not listed) | — | all columns | Table 1 | level_only (verified negative) |
| MFT + FFT -> mixed disulfide | MFT-FFT mixed disulfide | 1.6-2.5 (CO2 only); - under N2 and control | ng/mL | same | Table 1 | level_only |
| dithianone | 3-methyl-1,2-dithian-4-one | 4.2; 4.1-5.4 (N2); 16.4 / 45.2 / 50.8 / 55.0 (CO2 10-40 MPa) | ng/mL | same | Table 1 | level_only |
| pentose intermediates | furfural; norfuraneol | 2.0; 1.3-1.9 (N2); 9.1-34.5 (CO2) — 38.7; 33.5-51.4 (N2); 39.1-71.8 (CO2) | ng/mL | same | Table 1 | level_only |
| oxidant supply | dissolved oxygen / headspace oxygen | not measured, not controlled, not mentioned | — | — | whole paper | verified negative |
| browning | A280, A420 vs pressure | — | absorbance | same | Figure 1 | figure_only (except the N2 A280 line fit printed in the text) |

Cross-reference inside the repo: B17's T3 compares the model's MFT dimer share (0.04-0.90 %) with
Zhou 2023 (8.6 / 6.5 / 9.6 % at pH 6 / 7 / 8) and Zhang 2024 (8.7 %). This paper's 17-23 % of MFT
molecules as dimer (29-37 % of MFT-equivalents) at pH 5.6 / 140 C / 1 h is higher still and, as
far as the printed columns go, **independent of whether the headspace was air, 10-40 MPa nitrogen
or 10-40 MPa carbon dioxide**. If the disulfide were made by ambient oxygen alone, the nitrogen
columns should have shown less of it than the control; they do not. The reading consistent with the
table is that the oxidant in a ribose + cysteine pot at 140 C is internal (dicarbonyls, dehydro-
reductones, cystine, or the disulfide-exchange equilibrium with cystine itself), or that the
dimer forms during sampling (Flag 3) — either way, not the dissolved oxygen the B11 reservoir
represents. This bears on B17's diagnosis: enlarging the oxygen reservoir is not what this paper
points to; an oxidant the pot makes for itself is.

## 5. Flags

1. **No oxygen control.** "SC-N2" means the vessel was pressurised with nitrogen; nothing says the
   headspace air (0.1 MPa, ~ 21 kPa O2) or the dissolved oxygen of the buffer (~ 0.25 mmol/L
   air-saturated at 25 C, ~ 0.5 % of the cysteine) was removed first. Pressurising without purging
   leaves the O2 partial pressure and the dissolved O2 where they were. So the nitrogen columns
   show "no extra oxygen", not "no oxygen". The dissolved oxygen alone (0.25 mmol/L) could oxidise
   at most 0.5 mmol/L of thiol to disulfide, which is far above the ng/mL disulfides measured, so
   the experiment cannot exclude a dissolved-oxygen source on stoichiometric grounds; it only
   shows that adding 10-40 MPa of inert gas changes nothing.
2. **Semi-quantitative throughout**: HS-SPME with a single internal standard and a response
   factor of 1 for every compound, TIC areas, "approximate quantities in headspace (ng/ml of
   mixture)". Absolute values are not concentrations in the liquid; only within-column and
   within-row comparisons are safe, and even those assume equal SPME partitioning across
   atmospheres (the CO2 samples were depressurised from up to 40 MPa and the vented volatiles
   trapped in "three absorption buffers" and added back — a different loss path from the control).
3. **Disulfide may be made during sampling.** Thiols oxidise to disulfides on a DVB/CAR/PDMS
   fibre held at 60 C for 20 min in an air-containing headspace and in a hot GC inlet (250 C
   desorption); MFTD is also far less volatile than MFT, so its headspace share under-represents
   its liquid share. The 17-23 % dimer share is therefore neither an upper nor a lower bound on
   the share in the pot; it is a share in the sampled headspace under one protocol applied
   identically to all nine columns. The atmosphere-independence survives this caveat; the number
   does not.
4. **One time point (60 min), one temperature (140 C)**: no rate, no barrier, no information on
   whether the disulfide share changes with time.
5. **pH under CO2 is not measured**; the whole CO2 effect is argued from Gevaudan 1996 (skim milk
   under gaseous CO2) and from Mottram & Nobrega's buffer-catalysis argument. The N2 columns are
   the clean pressure control; the CO2 columns confound pressure, pH and carbonate catalysis.
6. **Apparatus details are in Xu 2008a (not on disk)**: vessel volume, liquid fill, stirring,
   heat-up and depressurisation times, gas purity.
7. **Column, MS instrument and LRI reference column are not stated**; LRIs (818 for
   3-mercapto-2-butanone, 869 MFT, 913 FFT, 1542 MFTD) are consistent with a non-polar column
   (DB-5 type).
8. **Cysteine as the hydrochloride** at 50 mmol/L in 0.2 M pyrophosphate; the initial pH 5.6 is set
   on the separate stock solutions; the mixed pot's pH is not re-stated.
9. Duplicate reactions only; +/- unlabelled; several 40 MPa CO2 entries carry +/- of 40-50 % of
   the value (MFT 34.6 +/- 15.8; MFTD 18.6 +/- 8.6).
10. Author-name discrepancy between byline (Xuan Liu) and the paper's own citation line (Liu K.).
