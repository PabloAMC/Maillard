# Joo et al. 2022 — EXTRACTION (bovine and chicken satellite-cell cultured tissue: IMP, amino acid composition)

**Source on disk:** `data/articles/Joo2022.pdf` (the journal's open-access PDF). First read 2026-09-14
from the HTML through an automated fetch; Table 1, the IMP paragraph and the methods then checked
against the PDF by eye later the same day: every number in section 2 matched. One open question was
closed by the methods (see section 1). Written for the cultivated-tissue composition box, not for
the core fit.

| field | value |
|---|---|
| Title | "A Comparative Study on the Taste Characteristics of Satellite Cell Cultured Meat Derived from Chicken and Cattle Muscles" |
| Venue | Food Science of Animal Resources 2022, 42(1), 175-185 |
| DOI | 10.5851/kosfa.2021.e72 |

## 1. Methods

Satellite cells from chicken pectoralis major and from Hanwoo biceps femoris (24-27-month steers).
Proliferation in glucose-free DMEM with 30 % FBS, 1 % Glutamax, 5 ng/mL bFGF; differentiation in
DMEM with 2 % horse serum; three passages, two weeks total. Cultured muscle tissue (CMT) compared
with the corresponding traditional meat (TM). Amino acids by Biochrom 30+ analyser AFTER ACID
HYDROLYSIS (1 g tissue, 6 M HCl, 110 °C, 24 h; confirmed at the PDF check), reported as PERCENT OF
TOTAL amino acids (Table 1): these are protein-bound totals, not a free pool. Nucleotide-related
compounds by HPLC against standards, in mmol/kg of the sample (Figure 2; wet or dry basis not stated);
electronic tongue.

## 2. Findings that matter

- Table 1 is a composition in percent, not a free pool in absolute units. Cattle CMT cysteine
  1.47 % vs TM 0.54 %; leucine 7.35 % vs 9.02 %; glycine 7.68 % vs 4.98 %. They are
  HYDROLYSED (total) amino acids, so the cysteine figure is protein cysteine plus any free cysteine,
  as a share of total amino acids; no absolute total is printed. Nothing converts to a free pool in mM.
- Figure 2: IMP in cattle CMT about 1.98 mmol/kg; cattle TM higher (roughly 2.5 to 3, read from
  the figure, not printed). AMP, inosine and hypoxanthine under 0.2 mmol/kg in CMT. Wet or dry
  basis not stated.

## 3. What the repo takes

- **IMP, cultivated bovine**: 1.98 mmol/kg, taken as wet basis at 75 % moisture → 2.6 mM. This
  is the upper corner of the cultivated IMP range; Kim 2024b's pig construct is the lower corner,
  four decades down.

The cysteine percentage is NOT taken: it is a percent of total amino acids after acid hydrolysis, so
it is mostly protein cysteine and cannot become a free-cysteine concentration. Free cysteine in
cultured bovine muscle remains unmeasured in the public record as of this read.

Note the proliferation medium was glucose-free, which is unusual and bears directly on what a
cultured cell's free glucose pool could be. Not measured in this paper.
