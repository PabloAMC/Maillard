# CF-041 (IMP in cultured / cultivated meat): pre-registration

Written 2026-11-03 after reading the two PDFs' methods, units and tables (no test yet run). No earlier cluster
worked this entry (`tests: []`).

## What the pages print (read by eye)

- **Joo 2022** (data/articles/Joo2022.pdf, p. 5 methods, p. 7 Fig. 2): cell sheet scraped, washed 3x in PBS and
  centrifuged (500 x g, 10 min each); "muscle tissues (1 g)" homogenised in 0.5 M perchloric acid; results "mmol/kg
  of the samples" (wet basis not written, but a 1 g as-is pellet, no drying step). Text: IMP 0.81 (chicken CMT)
  and 1.98 (cattle CMT) mmol/kg. Fig. 2B (read by eye, not printed): chicken TM ~24, cattle TM ~14 mmol/kg.
  The dossier's "cattle TM roughly 2.5 to 3" is a misreading of Fig. 2B.
- **Kim 2024b** (data/articles/Kim2024.pdf, p. 4 methods, p. 10 Table 4): "Contents (mg/kg)"; IMP pork 792.19,
  cultivated 0.11, scaffold ND; 5 g sample extracted with perchloric acid, freeze-dried EXTRACT, 1H NMR (Chenomx).
  Basis: as-is construct (cells in gelatin hydrogel, ~90 % moisture), not a cell pellet, not dry tissue.

## Hypotheses

- **H1 (unit / basis slip in one source).** H1a: Kim's 0.11 is a slip (e.g. g/kg or per-mg). H1b: Joo's absolute
  scale is inflated (all values, TM included).
- **H2 (real difference, two causes):** (i) Kim's construct is mostly hydrogel (1e7 cells/mL), so cell matter is
  diluted ~10-50x; (ii) neither cultured tissue runs the post-mortem ATP -> IMP pathway to completion, so IMP sits
  near the living-cell level, while Joo's pellet handling (scrape, three washes, spins) may let some ATP
  degrade to IMP.
- H1 and H2 are not exclusive.

## Tests (each uses evidence outside the conflict)

- **T1 (Kim's scale):** compare Kim's own pork control (792 mg/kg = 2.28 umol/g wet) with the reference meats:
  chicken breast 2.7-5.3 umol/g (aliani2013), beef 2.5-6.3 mmol/kg (koutsidis2008a/b), beef D1-D14 3.4-7.6
  umol/g (muroya2019). In range -> Kim's mg/kg scale is right, H1a rejected.
- **T2 (Joo's scale):** compare Joo's TM (~24 chicken, ~14 cattle mmol/kg) with the same references and with the
  ceiling of the total adenine nucleotide pool (ATP at slaughter ~6.5 umol/g in muroya2019 D0; IMP cannot exceed
  the ATP + ADP + AMP pool at death). TM above the pool by > 2x -> Joo's absolute scale inflated (H1b supported);
  within the references -> H1b rejected.
- **T3 (living-state IMP):** muroya2019 D0 (pre-rigor, ATP still 6.5 umol/g) IMP = 78 nmol/g. If the
  living-state IMP x a plausible cell fraction of Kim's construct (2-10 %) brackets Kim's 0.32 nmol/g within a
  decade, H2 explains Kim's corner without any slip.
- **T4 (box effect, no regeneration):** the frozen artifacts results/cultivated_tissue_invariance/*.json hold every
  draw. IMP is in UNRANKABLE, so it can only enter `naive_ranking_all_declared`. Remap each draw's cultivated IMP
  to a proposed range through the same uniform quantile (the draw count per entry does not change, so this is
  what a regeneration with the same seed would give) and count changed fields. Prediction: zero change in the
  engine ranking, naive candidate ranking and every aggregate; only the per-draw all-declared ranking moves.

## Decision rule

- T1 in range and T2 above pool -> status `explained` for the 4-decade gap as H1b (~0.5 decade) + H2 (rest);
  Joo's absolute values carried as ratios to its own TM, not as absolute mM.
- T1 out of range -> H1a; re-read Kim's table for a slip.
- T2 inside the references -> the gap is all H2; keep Joo's absolute value.
- Box change only if T2 supports H1b and T4 shows no effect on any aggregate; otherwise note-only edit.

## What would falsify

- A Joo TM inside 2-8 umol/g on re-read (then H1b falls). A Kim pork value outside the reference range (then
  Kim's scale is suspect). A cell fraction in Kim's construct near 100 % (then dilution fails and the gap
  needs a slip).
