# Pre-registration: wave B35, the benchmark provenance audit (written 2026-09-11, BEFORE any prediction was probed)

## 1. Why

B34 found that a hold-out bundle scored **one** of the six species its own paper measures, cited an
author who is not on the paper, and declared the source missing while the PDF sat on disk. That was
found by accident. **This wave asks whether it was the only one**, by cross-referencing all 48
benchmark bundles against the 245 PDFs and 255 dossiers on disk.

## 2. What the audit found, before anything is changed

**A. Three provenance defects.**

1. `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021` **still** says `SOURCE NOT ON DISK` in its
   **vessel** note. B34 corrected the buffer note and the citation and missed this one — my own fix,
   incomplete, found by the audit that followed it.
2. `external_validation_liu_2023_ppi_offnote_baseline` **contradicts itself**: its buffer note says
   "SOURCE PAPER NOT ON DISK; second-hand" and its vessel note says "the PDF is on disk
   (`data/articles/liu2023.pdf`, read 2026-09-04)". The PDF exists, 717 KB, dated 4 September. The
   vessel note was updated when the paper arrived and the buffer note was not.
3. The other **eleven** bundles claiming a missing source are telling the truth: no dossier, no PDF.
   That is worth recording as a clean result, not just the failures.

**B. Three bundles score fewer compounds than their on-disk papers print, with SDs.**

| bundle | scores | the same table also prints |
|---|---|---|
| `external_validation_bi_2020_roasted_pea_hexanal` | hexanal | **furaneol 2780 ± 125**, **2,5-dimethylpyrazine 5960 ± 77.6**, **furfural 327 ± 6.82** µg/kg |
| `external_validation_bi_2020_raw_pea_hexanal` | hexanal | **nonanal 69.8 ± 7.39** µg/kg |
| `pea_isolate_uht_140C_Trikusuma2019` | hexanal, 2-pentylfuran, nonanal | **2,5-dimethylpyrazine**, **methional**, **2-acetyl-1-pyrroline**, **(E,E)-2,4-decadienal** — each with a *control* column as well as a UHT one |

## 3. What is built

**No constant moves and nothing is fitted.** Every value below is an end-of-cook level in an external
or hold-out bundle: they VALIDATE. Three changes.

1. The two stale provenance notes are corrected, with the false claims retained and labelled
   superseded, as B34 did.
2. The eight measurements above are added as targets, each with its printed value, SD and a verbatim
   quote.
3. Where the source prints a control column (Trikusuma), the control value is ALSO declared under
   B31's `conditions.carried_volatiles`, because that is what that field is for and declaring it for
   three compounds and not the other four would be arbitrary.

## 4. What counts as success, declared before the run

- **T1 no existing row moves.** Every currently scored row keeps its prediction to the last bit.
- **T2 every added row is either answered or REFUSED BY NAME** — none may return a silent zero or a
  number in the wrong unit. B34 found the mmol/L fallback; this is the test that it stays fixed.
- **T3 the refusals are the ones the record predicts.** Methional is B22's product and 2-acetyl-1-
  pyrroline is B24's, and both waves were refused with their steps left inert, so the engine must
  refuse those two by name rather than answer them.
- **T4 the headline is reported in both directions**, as B34's was.

Ship rule: **SHIP if T1 and T2 hold.** T3 and T4 are reported.

## 5. Predictions, written before probing any of them

1. T1 holds. **95 %.**
2. **Methional and 2-acetyl-1-pyrroline are refused by name.** **85 %.** Their waves were refused and
   their steps ship inert, so an answer would mean a structural zero dressed as a prediction — which
   is the failure B28 and B34 both turned out to be.
3. **Furaneol and 2,5-dimethylpyrazine in the roasted pea are answered and miss by more than 10×.**
   **70 %.** Both are trunk products of a real roasted food with no declared sugar or amine charge;
   the bundle's charge was built to answer hexanal.
4. **At least one added row lands within 3×.** **35 %.** Four of the eight are on lanes that have
   never been tested on a whole food.
5. The panel's within-3× *rate* falls. **65 %.** Eight rows added, few expected inside the band.
6. **No further bundle has this defect.** **60 %** — the audit covered every bundle with a DOI and a
   dossier, but a paper on disk with no dossier, or a bundle whose DOI is absent, could still hide one.

---

# Outcome (2026-09-11)

## Ship: T1 and T2 hold

| test | result | pass |
|---|---|---|
| T1 no existing row moves | 44 rows compared, **not one prediction changed** | **yes** |
| T2 every added row answered or refused BY NAME, no silent zero | **failed on the first run and is the wave's main finding** — fixed, then held | **yes** |
| T3 the refusals are the ones the record predicts | methional (B22), 2-acetyl-1-pyrroline (B24), 2,5-dimethylpyrazine (B18), all refused by name | reported |
| T4 both directions | within-3× **9/44 → 10/45**; refused **25 → 32**; median fold **10.62× → 9.31×** | reported |

**SHIPPED.**

## T2 failed first, and that is the third instance of one bug family

On the first run `furaneol` and `furfural` were **answered with predicted = 0.0** against measurements
of 2780 and 327 µg/kg. Not a miss — an absence of a prediction reported as one.

The cause is general and was latent in every matrix-only bundle. Such a pot charges a protein isolate
and nothing else; the isolate is a LIPID CARRIER and is deliberately kept out of `mapped_precursors`,
so the trunk, sulfur and acrylamide networks integrate from an all-zero state and **every species in
them is zero by construction**. `5-HMF` and `acrylamide` would have done the same on any of the panel's
seven matrix-only bundles the moment anyone asked.

This is the same family as B28's 2-pentylfuran (reported in the wrong unit, diagnosed for a day as a
routing problem) and B34's silent mmol/L fallback (still live three days after B28). B28's own record
named the principle and the code kept violating it one level up: *"A near-zero is the absence of a
prediction dressed as one."*

`declare_envelope` now refuses a non-lipid target on a pot with no charged precursor, names the
compounds, and names the cure. Two tests hold it: the refusal itself, and a panel-wide assertion that
**no scored row anywhere has a prediction of exactly zero**.

## What the eight audited measurements did

| bundle | outcome |
|---|---|
| Trikusuma | **(E,E)-2,4-decadienal answered at 1.62×, inside the band** — a fourth hit in that pot. 2,5-dimethylpyrazine, methional and 2-acetyl-1-pyrroline **refused by name**, their waves (B18, B22, B24) having been refused with their steps inert. |
| Bi 2020 roasted | furaneol and furfural **refused** (no precursor charged); 2,5-dimethylpyrazine **refused** (B18, lane). |
| Bi 2020 raw | nonanal **refused** — the pot was never cooked, by B31's rule, exactly as its hexanal already was. |

Every one of the seven Trikusuma compounds the source prints a control column for now has a carried
level declared under B31, instead of three of seven.

## Predictions, scored

1. T1 holds (95 %) — **right**.
2. Methional and 2-acetyl-1-pyrroline refused by name (85 %) — **right**.
3. Furaneol and 2,5-dimethylpyrazine answered and missing by >10× (70 %) — **wrong, and the outcome
   is better than the prediction**: they are refused. I predicted a bad number where the honest
   answer was no number.
4. At least one added row within 3× (35 %) — **right**, decadienal at 1.62×.
5. The within-3× rate falls (65 %) — **wrong**. It rose, 20.5 % → 22.2 %, because five of the eight
   rows became refusals rather than misses.
6. No further bundle has this defect (60 %) — **not resolved by this wave**. The audit covered every
   bundle with a DOI and a dossier; a paper on disk with no dossier, or a bundle with no DOI, could
   still hide one. Recorded as open.

## The clean half of the audit, worth saying

Eleven of the fourteen bundles that declare their source missing are **telling the truth** — no PDF,
no dossier. The provenance notes are mostly right, and the two that were wrong are now corrected with
their false claims retained and labelled.
