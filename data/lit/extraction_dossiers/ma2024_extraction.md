# Ma, Fu, Cheng & Liu 2024 — EXTRACTION (extrusion parameters and CML, CEL and acrylamide in a soy-protein-isolate meat analogue)

**Source on disk:** `data/articles/Ma2024.pdf` (1.6 MB; downloaded 2026-09-11 at the 2026-09-11
reading-list row's request). Read 2026-09-11 via `pdftotext`; pages 6–7 (sec. 2.3 and Figure 2)
rendered and read as images because the scored value is a bar in Figure 2D. Wave B36.

| field | value |
|---|---|
| Title | "Impact of Extrusion Parameters on the Formation of Nε-(Carboxymethyl)lysine, Nε-(Carboxyethyl)lysine and Acrylamide in Plant-Based Meat Analogues" |
| Venue | Int. J. Mol. Sci. 2024, 25, 8668 |
| DOI | 10.3390/ijms25168668 |
| Authors | Yurong Ma, Shuang Fu, Ka-Wing Cheng, Bin Liu |
| Systems | soy protein isolate : corn starch 9:1 (w/w), pilot twin-screw extruder with ten heating zones; single-factor series over moisture (20/30/40/50/60 %), screw speed (120/150/180 rpm), feed rate (4/6/8 kg/h) and barrel temperature (three zone sets, see §2); an unextruded control |
| What is measured | water, protein and total amino acids (Table 1); CML, CEL, GO, MGO and acrylamide by UHPLC-MS/MS with isotope-labelled internal standards (acrylamide-d3); **acrylamide only in Figure 2, no table** |
| Bundle | `acrylamide_spi_extrusion_130C_ACSRef3` (a FIT row; "ACSRef3" is a name older than this DOI, as the bundle's correction note records) |

## 1. What this paper is and is not, for this model

A **process-parameter survey** with one acrylamide value per condition, all in bars. It is the only
source behind the model's single extrusion acrylamide row, and that row scores 4 247× low for
reasons the 2026-09-11 WAVES row already names (a process declared as one isothermal hold). What
the print adds is the process itself, and one thing the bundle did not know.

## 2. The extrusion, from sec. 3.2 (page 12–13), verbatim where quoted

"the temperature of the first 5 heating zones from the feed to the die was maintained at 80, 80, 85,
90 and 100 °C, respectively, while the temperatures of the last 5 heating zones was set at 120, 140,
150, 150 and 150 °C. The screw speed was 150 rpm, with a feed rate of 6.0 kg/h for raw material and
2.57 kg/h for water. The moisture content of the mixed raw material was 30%." In the temperature
series "the temperature of the last 5 heating zones was set as follows: (1) 110, 120, 130, 130 and
130 °C; (2) 120, 140, 150, 150 and 150 °C; (3) 120, 140, 170, 170 and 170 °C."

So the bundle's **"130 °C" is arm (1)**: zones 80/80/85/90/100/110/120/130/130/130 °C at 150 rpm,
6.0 + 2.57 kg/h, 30 % moisture. **The paper prints no residence time.** The bundle's 25 s is not in
the source; B36 labels it an assumption in the vessel note and leaves it, because a benchmark bundle
is an isothermal hold and the row is a fit row whose report is frozen.

## 3. Figure 2D, read from the image (µg/kg, mean ± SD, n = 3; 0–160 axis)

| condition | bar |
|---|---:|
| 130 °C arm | **≈150** (letter a) |
| 150 °C arm | ≈120 (b) |
| 170 °C arm | ≈82 (c) |
| unextruded control | **≈38** (d) |

Prose (sec. 2.3): "The lowest and highest acrylamide content were found in PBMAs extruded by 20%
moisture and 130 °C, respectively. These values fall within the range reported for commercially sold
PBMAs (32–187 µg/kg)". Figure 2A–C: moisture bell-shaped with a maximum at 40 %; screw speed a weak
effect; feed rate rising 4 → 8 kg/h.

## 4. The bundle, checked against the print (wave B36)

| item | bundle | print | verdict |
|---|---|---|---|
| acrylamide | 150 µg/kg ± 20 % (figure read-off) | Figure 2D 130 °C bar ≈150 | matches |
| process | 130 °C / 25 s isothermal | ten zones, arm (1); no residence time | 130 °C is the last three zones; 25 s unsupported, now labelled |
| starting level | none declared | control ≈38 µg/kg | **named, not acted on** (below) |

**The control bar.** The raw blend already carries ≈38 µg/kg of acrylamide before extrusion, so the
amount *formed* in the 130 °C arm is ≈112 µg/kg, not 150. Declaring the 38 as a carried level is the
cure Trikusuma took in B35, but here it is a figure read on a FIT row with a frozen B3 report, and a
25 % change in the formed amount on a row that is 4 247× low moves no verdict. Recorded here and in
`docs/guides/EXPERIMENTS.md`; not applied.

## 5. M0 data corrections, batch 1 (2026-11-06; CF-131)

The bundle's charge {SPI 50, L-Asn 42, D-Glc 42 mM}, pH 6.0 and a_w 0.35 are **not in this paper** (first reader:
`m0_integrity_batch1_other_lanes`, MF-39; second reader 2026-11-06 against the PDF). Secs. 3.1-3.2 name soy protein
isolate and corn starch (9:1 w/w) as the only extrusion materials -- no asparagine or glucose is added -- and Table 1
prints water, protein and total amino acids only. "Water activity" occurs nowhere; the only pH values are an amino-acid
hydrolysate's (6-8) and an assay buffer's (7). No specific mechanical energy is printed either (the bundle's 250 kJ/kg).
Applied: the charge and a_w withdrawn (nulls with notes, `value_withdrawal`), the acrylamide row REFUSED
(`scoring_refused`: "precursor charge not printed"). The pH stays as a number only because the schema requires one
(`ph_note`). Fig. 2D re-read: 130 C bar ~150 ug/kg, kept as the datum. The frozen B3 report still carries the row
(CF-131, open).

**Added the same day (coordinator review of PR #88).** Also withdrawn as not in the paper: the bundle's
`sme_kj_per_kg` 250 (no kJ/kg value is printed; SME is discussed only qualitatively); the intake record
`acs_ref3_spi_acrylamide_fast_kinetics`'s 22.36 -> 62.62 ug/kg and 25 s (no such numbers; acrylamide appears only in
Fig. 2), with its three citation aliases ("ACS Food Sci. Tech. 2022 (Ref. 3)", "ACS ResearchGate ...") closed; and its
twin in `safety_reference_payloads.json` (112 ug/kg "in heated SPI at 150 C"; the 150 C bar reads ~120 and no number is
printed), deleted. None of them was read by a scored row.
