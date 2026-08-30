# Kinetic core B4 -- matrix / OAV output layer: HOLD-OUT REPORT

Generated 2026-08-28.

## Flagship: Hong 2020, ten paired soy/water threshold ratios

**Criterion:** >= 7/10 within 5.0x AND correct sign on all 10, including the ethyl-4-methylpentanoate inversion (Amendment 4)

**RESULT: FAIL.** 0/10 within 5x; 3/10 signs correct.

**Not fully blind.** the build brief instructed this wave to read k4b_paired_thresholds_and_browning.md, which prints Hong's aggregate statistics and four per-compound facts outside the paired table. On the 5 rows where nothing leaked: 0/5 within 5x, 2/5 signs correct. the clean-row score is the one that carries evidential weight; the ten-row score is reported because the declaration's criterion is defined over ten rows.

### Scorecard, row by row

| # | compound | class | measured | predicted | fold | sign meas. | sign pred. | sign ok | <=5x | leaked | model state |
|---:|---|---|---:|---:|---:|---|---|---|---|---|---|
| 1 | 2,3-dimethyl pyrazine | alkylpyrazine | 48.88 | 1 | 48.88 | elevated | no_prediction | NO | no | no | no_binding_constant_for_class |
| 2 | ethyl-4-methylpentanoate | ester | 0.1546 | 1.697 | 10.98 | inverted | elevated | NO | no | yes | predicted |
| 3 | 2-pentyl furan | alkylfuran | 1359 | 1 | 1359 | elevated | no_prediction | NO | no | no | no_binding_constant_for_class |
| 4 | 4-vinyl phenol | phenol | 151.4 | 1 | 151.4 | elevated | no_prediction | NO | no | yes | no_binding_constant_for_class |
| 5 | hexanal | n_alkanal | 132.5 | 2.634 | 50.3 | elevated | elevated | yes | no | yes | predicted |
| 6 | 3-methyl butanal | branched_alkanal | 263.2 | 1.582 | 166.4 | elevated | elevated | yes | no | no | predicted_with_chain_length_surrogate |
| 7 | 2-methyl butanal | branched_alkanal | 261.4 | 1.582 | 165.3 | elevated | elevated | yes | no | no | predicted_with_chain_length_surrogate |
| 8 | butyric acid | carboxylic_acid | 2035 | 1 | 2035 | elevated | no_prediction | NO | no | yes | no_binding_constant_for_class |
| 9 | 4-ethyl phenol | phenol | 29.1 | 1 | 29.1 | elevated | no_prediction | NO | no | yes | no_binding_constant_for_class |
| 10 | dimethyl disulfide | disulfide | 1755 | 1 | 1755 | elevated | no_prediction | NO | no | no | no_binding_constant_for_class |

### Residual decomposition, per compound

*measured shift = reversible binding x unsaturation x covalent (inert) x UNEXPLAINED RESIDUAL*

| compound | measured | reversible | unsaturation | covalent | **unexplained residual** | explained share of log |
|---|---:|---:|---:|---:|---:|---:|
| 2,3-dimethyl pyrazine | 48.88x | 1x | 1x | 1x | **48.88x** | 0.0% |
| ethyl-4-methylpentanoate | 0.1546x | 1.697x | 1x | 1x | **0.0911x** | -28.3% |
| 2-pentyl furan | 1359x | 1x | 1x | 1x | **1359x** | 0.0% |
| 4-vinyl phenol | 151.4x | 1x | 1x | 1x | **151.4x** | 0.0% |
| hexanal | 132.5x | 2.634x | 1x | 1x | **50.3x** | 19.8% |
| 3-methyl butanal | 263.2x | 1.582x | 1x | 1x | **166.4x** | 8.2% |
| 2-methyl butanal | 261.4x | 1.582x | 1x | 1x | **165.3x** | 8.2% |
| butyric acid | 2035x | 1x | 1x | 1x | **2035x** | 0.0% |
| 4-ethyl phenol | 29.1x | 1x | 1x | 1x | **29.1x** | 0.0% |
| dimethyl disulfide | 1755x | 1x | 1x | 1x | **1755x** | 0.0% |

### Flags raised by the decomposition

- **2,3-dimethyl pyrazine**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **2,3-dimethyl pyrazine**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.
- **2-pentyl furan**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **2-pentyl furan**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.
- **4-vinyl phenol**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **4-vinyl phenol**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.
- **butyric acid**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **butyric acid**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.
- **4-ethyl phenol**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **4-ethyl phenol**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.
- **dimethyl disulfide**: NO TERM AVAILABLE: the entire measured shift is unexplained residual, because the corpus supplies no constant for this class.
- **dimethyl disulfide**: SIGN NOT PREDICTED: the model emits exactly 1.0, which is the absence of a prediction, not a prediction of no shift.

### Pre-registered expectations, checked against the outcome

| expectation | outcome | detail |
|---|---|---|
| ethyl-4-methylpentanoate: SIGN WRONG (~95 %) | **HELD** | measured inverted, predicted elevated, fold 11 |
| six compounds get NO prediction at all, scored as sign-fail | **HELD** | 6 compounds emitted no prediction: 2_3_dimethylpyrazine, 2_pentylfuran, 4_ethylphenol, 4_vinylphenol, butyric_acid, dimethyl_disulfide. Every one of them is a sign-fail and carries 100 % of its shift as unexplained residual. |
| hexanal and the two butanals: SIGN RIGHT, MAGNITUDE FAR TOO SMALL, fold error of order 20-30x (~85 %) | **HELD ON SIGN, UNDER-STATED ON MAGNITUDE** | all three signs correct; fold errors 50.3x, 166.4x, 165.3x. The pre-registration guessed 20-30x and the truth is 50.3-166x, so the model is WORSE than this wave expected, not better. |
| the two branched butanals get IDENTICAL predictions, because branch position is unmeasured anywhere | **HELD, AND COST ALMOST NOTHING** | predictions identical at 1.582x; and Hong measured the two only 1.01x apart (263.2x vs 261.4x). The layer's structural blindness to branch position turns out to be an accurate blindness on this pair -- an unexpectedly favourable result for a limitation that was declared as a defect. |
| Amendment 6 ruling 2: reversible binding explains ~25 % of a matrix log-shift and no more | **CORROBORATED OUT OF SAMPLE** | on the three rows where the term is active the explained share hexanal 19.8%, 3_methylbutanal 8.2%, 2_methylbutanal 8.2% -- every one BELOW the ~25 % ceiling, on a matrix and a panel the ceiling was not derived from (it came from beef and from dairy protein). This is the strongest positive result in the report. |
| (not pre-registered; recorded on discovery) | **THE UNSATURATION TERM WAS NEVER EXERCISED** | none of Hong's ten compounds is an alpha,beta-unsaturated CARBONYL -- 4-vinyl phenol has a conjugated C=C but no carbonyl -- so the layer's second named term contributed nothing to any row. The flagship hold-out tests one of the three terms, not three. That is a gap in the hold-out's coverage, not in the layer, and it means the ~3.7x penalty remains unvalidated out of sample. |

### What the residual is, and what it is not

The unexplained residual spans **29x to 2 035x** and the three named terms
account for at most **20 %** of any row's log-shift. The residual does not
correlate with anything the layer carries: the largest is a log P 1.0 acid
and one of the smallest is a log P 2.6 phenol, which is Hong's own
`r = -0.05` finding reproduced from the other side. The leading NAMED
candidates for it, none implemented and each with its reason, are in the
fit report's *Named terms NOT implemented* section: lipid-phase partition
(fittable, but the soy fat content is unreported), background-odour masking
(the only mechanism that can produce an inversion, and unimplementable
because Hong measured no background volatiles), delivery kinetics, and the
criterion bias between a 50 % forced-choice and a 75 % uncorrected task.

## Other declared hold-outs

### Zhou 2023 SI Table S2 OAV arithmetic (HOLD-OUT (arithmetic check on the OAV code path))

- potency ratio dimer/monomer reproduced: **True** (15.62x vs expected 15.625x)
- OAV arithmetic exact: **True** (MFT 3.18e+05 vs printed 3.18e+05; dimer 3.21e+05 vs printed 3.21e+05)
- The dimer's OAV matches the monomer's while carrying 6.5-9.6 % of the MFT-equivalents. Mass lost to MFT dimerisation is NOT aroma lost, and any objective that scores the dimerisation channel as a pure loss is wrong by roughly the threshold ratio.

### Brewer 1995 beef (HOLD-OUT, reclassified dose_added_pre_cook (D.6 Module 7))

**not_scored_by_ruling.** Brewer's 'thresholds' are the dose you must add to RAW ground beef BEFORE a 70 C cook for the compound to be detectable AFTER cooking. Every thermally driven loss between dosing and sniffing -- evaporation, Schiff-base condensation with lysine on a denaturing protein, further oxidation -- comes straight off the numerator, and the paper analytically verifies nothing at the moment of sniffing. Scoring the layer against it would be scoring a perception model against a thermal-loss measurement.

Hong's properly-controlled hexanal soy ratio sits an order of magnitude below Brewer's beef ratio, which is exactly what k2 sec. D.2(ii) predicted would happen once the pre-cook thermal loss was removed. That is an out-of-sample confirmation of the reclassification, and it is the one thing this hold-out row genuinely scores.

### Frozen matrix-path bundles

The 4 frozen matrix-path bundles and 17 maillard-path bundles under `data/benchmarks/external_validation/` were NEVER OPENED by this wave -- not to build the layer, not to score it, and not to check a name. They remain `evidence_class: external_validation_only`.

