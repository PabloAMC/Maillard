# Kim, Jung & Jo 2024 — EXTRACTION (3D-cultured pig muscle stem cells: free amino acids, nucleotides)

**Source on disk:** `data/articles/Kim2024.pdf` (the publisher's PDF). First read 2026-09-14 from the
open-access full text at PMC (PMC11811360) through an automated fetch; Tables 2 to 4, the culture
protocol and the moisture statement then checked against the PDF by eye the same day. Every number in
section 2 matched the PDF. One correction: the construct's ~90 % moisture is stated in the text and
Fig. 3A, not in Table 3 (Table 3 is total amino acids per gram of protein). The stem is `kim2024b`
because `kim1998` already exists and a `kim2024` may follow. Written for the cultivated-tissue
composition box, not for the core fit.

| field | value |
|---|---|
| Title | "Fundamental study on structural formation, amino acids and nucleotide-related compounds of cultivated meat from 3D-cultured pig muscle stem cells" |
| Venue | Food Science and Biotechnology, 34(2), 457-469 (online 2024) |
| DOI | 10.1007/s10068-024-01793-9 |

## 1. Methods

Porcine muscle stem cells at 1e7 cells/mL in a cross-linked gelatin hydrogel, 44-day culture,
analysed at day 16. Control: commercial pig leg meat. Scaffold-only arm analysed alongside. Free
amino acids and nucleotide-related compounds by HPLC, reported in mg/kg (Table 4 for free amino
acids). Moisture of the cultivated construct about 90 % (text and Fig. 3A).

## 2. Findings that matter

Free amino acids, mg/kg (Table 4; letters are the paper's significance groups):

| compound | pork | cultivated | scaffold only |
|---|---|---|---|
| methionine | 61.9 | 10.5 | 12.4 |
| leucine | 115.5 | 35.8 | 49.6 |
| glycine | 72.8 | 26.3 | 18.4 |
| alanine | 138.9 | 82.8 | 4.6 |
| glutamic acid | 115.2 | 69.7 | 5.5 |
| total free amino acids | 1 463 | 550 | 609 |

Cysteine is not reported in Table 4.

Nucleotide-related compounds, mg/kg:

| compound | pork | cultivated | scaffold only |
|---|---|---|---|
| IMP | 792 | 0.11 | not detected |
| inosine | 630 | 1.22 | not detected |
| hypoxanthine | 121 | 1.48 | not detected |

## 3. What the repo takes

Cultivated-side points for the composition box, in mM in tissue water at 90 % moisture
(mg/kg ÷ molar mass ÷ 0.90):

- **IMP**: 0.11 mg/kg ÷ 348.2 g/mol ÷ 0.90 = 0.00035 mM. Three and a half decades below the
  pork control (792 mg/kg = 3.2 mM) and three and a half decades below the one bovine cultured
  measurement (Joo 2022, 1.98 mmol/kg). The two published cultivated IMP values disagree by
  four orders of magnitude; the box carries both corners.
- **leucine (free)**: 35.8 mg/kg ÷ 131.2 ÷ 0.90 = 0.30 mM. CONFOUNDED: the scaffold-only arm
  carries 49.6 mg/kg of free leucine, more than the construct, so the gelatin contributes an
  unknown share of the construct's value. Carried as one point of a wide range, pig not bovine.

Nothing here for ribose, glucose, thiamine, cysteine or ribose 5-phosphate in cultured tissue. The PDF was
searched for each of them on 2026-09-14: glucose appears only in a cited C2C12 fasting experiment
(no medium glucose, alanine rose); ribose, thiamine and cysteine do not appear at all.
