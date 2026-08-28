# Kinetic core B4 -- matrix / OAV output layer: FIT REPORT

Generated 2026-08-28 by `scripts/generators/generate_kinetic_core_b4_fit.py`.

## The four-line answer

1. **Ratios lead.** The layer's primary output is formulation-vs-formulation
   ratios with a validity bound; absolutes are emitted only as intervals
   carrying the measured reliability band (HS-SPME same-sample dispersion
   10-23x, air/water K +/-0.5 decades).
2. **Three named terms, and nothing beyond them.** Reversible binding
   (capped at ~25 % of a log-shift), the alpha,beta-unsaturation penalty,
   and the covalent ceiling as an INERT bound. Everything else is reported
   as quantified unexplained residual.
3. **The corpus supports a matrix prediction for 4 of the flagship
   hold-out's 10 compounds and no more.** 6 belong to classes for
   which the entire FIT column contains no binding constant.
4. **The gate is expected to fail, by construction.** The declared criterion
   requires the correct sign on all ten including an inversion; nothing in
   the FIT column produces a sign reversal for an ester.

## Fitted constants (FIT rows only)

| class | K_g (L/g) | n FIT rows | members | sign | sources |
|---|---:|---:|---|---|---|
| branched_alkanal | 0.004096 | 0 | (chain-length surrogate) | suppression | Meynier2002_extraction.md |
| diketone | 0.073 | 1 | diacetyl | suppression | leksrisompong2010_extraction.md |
| ester | 0.004908 | 3 | isoamyl_acetate, amyl_acetate, ethyl_pentanoate | suppression | Meynier2002_extraction.md |
| furanone | -0.062 | 1 | furaneol | enhancement_measured | leksrisompong2010_extraction.md |
| lactone | -0.023 | 1 | delta_decalactone | enhancement_measured | leksrisompong2010_extraction.md |
| methyl_ketone | 0.02485 | 3 | 2_heptanone, 2_octanone, 2_nonanone | suppression | andriot2000 via k2 |
| n_alkanal | 0.01151 | 1 | hexanal | suppression | Meynier2002_extraction.md |

**alpha,beta-unsaturation penalty: 3.73x** from 2 FIT rows (unsat_penalty_gelatin 2.81x, unsat_penalty_dairy_headspace 4.95x), spread 1.76x.

> k2 sec. D.3 states the penalty as ~2-3x. This FIT-only estimate is ABOVE that band, and it is above it BECAUSE the observation that would pull it down (beef, 2.01x) is hold-out and excluded. k4b records the same disagreement: 'the term's magnitude is matrix-dependent'. Reported, not smoothed.

## FIT-row reproduction

| compound | matrix | measured | predicted | fold | sign ok | state |
|---|---|---:|---:|---:|---|---|
| hexanal | skim_milk | 1.39 | 1.390189 | 1.000136 | yes | predicted |
| ethyl_pentanoate | skim_milk | 1.33 | 1.166396 | 1.140264 | yes | predicted |
| amyl_acetate | skim_milk | 1.2 | 1.166396 | 1.02881 | yes | predicted |
| isoamyl_acetate | skim_milk | 1.07 | 1.166396 | 1.09009 | yes | predicted |
| diacetyl | caseinate_1pct | 1.73 | 1.73 | 1.0 | yes | predicted |
| delta_decalactone | caseinate_1pct | 0.77 | 0.77 | 1.0 | yes | predicted |
| furaneol | caseinate_1pct | 0.38 | 0.38 | 1.0 | yes | predicted |

## The Vega gelatin ladder as lookup-table entries

| compound | water (ug/L) | gelatin 22 C (ug/L) | measured ratio | predicted | state |
|---|---:|---:|---:|---:|---|
| pentanal | 12.0 | 41.0 | 3.416667 | 1.3453 | predicted |
| hexanal | 4.5 | 58.0 | 12.888889 | 1.3453 | predicted |
| heptanal | 3.0 | 79.0 | 26.333333 | 1.3453 | predicted |
| t_2_hexenal | 3.0 | 109.0 | 36.333333 | 3.729544 | partial_unsaturation_term_only |
| t_2_octenal | 3.0 | 109.0 | 36.333333 | 3.729544 | partial_unsaturation_term_only |
| tt_2_4_decadienal | 0.07 | 64.0 | 914.285714 | 3.729544 | partial_unsaturation_term_only |

## PRE-REGISTERED blind predictions -- Hong 2020 flagship

| # | compound | class | predicted ratio | interval | predicted sign | state |
|---:|---|---|---:|---|---|---|
| 1 | 2,3-dimethyl pyrazine | alkylpyrazine | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |
| 2 | ethyl-4-methylpentanoate | ester | 1.696999 | [1.490845, 1.981689] | elevated | predicted |
| 3 | 2-pentyl furan | alkylfuran | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |
| 4 | 4-vinyl phenol | phenol | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |
| 5 | hexanal | n_alkanal | 2.63442 | [2.151, 3.302] | elevated | predicted |
| 6 | 3-methyl butanal | branched_alkanal | 1.581644 | [1.409609, 1.819217] | elevated | predicted_with_chain_length_surrogate |
| 7 | 2-methyl butanal | branched_alkanal | 1.581644 | [1.409609, 1.819217] | elevated | predicted_with_chain_length_surrogate |
| 8 | butyric acid | carboxylic_acid | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |
| 9 | 4-ethyl phenol | phenol | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |
| 10 | dimethyl disulfide | disulfide | 1.0 | [1.0, 1.0] | flat | no_binding_constant_for_class |

## Pre-registered expected failures

**Criterion:** Declaration Amendment 4: >=7/10 of the paired soy/water ratios within 5x AND correct sign on all 10, including the ethyl-4-methylpentanoate inversion.

**Overall expectation: FAIL, and by construction rather than by accident.**

### ethyl_4_methylpentanoate -- SIGN WRONG (very high (~95%))

The declaration names this compound as carrying an INVERSION. The only term this layer has for an ester is Meynier's measured ester SUPPRESSION (1.07-1.33x at 33.9 g/L dairy protein), which can only push the ratio ABOVE 1. Nothing in the FIT column produces a sign reversal for an ester: the two measured enhancements in the corpus are a lactone and a furanone, not an ester. The one mechanism that could -- background-odour masking -- is unimplementable, because Hong measured no background volatiles in the soy matrix. So this single compound fails the 'correct sign on all 10' clause on its own, and therefore the gate fails whatever the magnitudes do.

### 2_3_dimethylpyrazine, 2_pentylfuran, 4_vinylphenol, butyric_acid, 4_ethylphenol, dimethyl_disulfide -- NO PREDICTION AT ALL -> scored as sign-fail (certain (structural))

Six of the ten compounds belong to classes for which the entire FIT column contains NO per-gram binding constant: alkylpyrazine, alkylfuran, phenol (x2), carboxylic acid, disulfide. The layer reports `no_measured_constant_for_this_class` and emits 1.0, which is the ABSENCE of a prediction, not a prediction of no shift. Imputing a constant for them would be inventing a parameter, and a class-mean over unrelated chemistry would be the general correction factor k2 sec. D.1 refutes. This is the single most informative thing the hold-out will show: the corpus supports a matrix prediction for 4 of Hong's 10 compounds and no more.

### hexanal, 3_methylbutanal, 2_methylbutanal -- SIGN RIGHT, MAGNITUDE FAR TOO SMALL (high (~85%))

These three do get a term. But Amendment 6 ruling 2 already computed the answer: reversible binding explains ~25 % of a hexanal matrix log-shift and covalent ~0.06 %. A model whose only active term is reversible binding must therefore UNDER-PREDICT a real matrix shift by roughly the remaining 0.75 of the log, which for a shift of order 100x is a fold error of order 20-30x -- far outside the 5x window. If any of these three lands inside 5x it will be because the true shift is small, not because the model captured a large one.

### 2_methylbutanal, 3_methylbutanal -- IDENTICAL PREDICTIONS FOR TWO DIFFERENT COMPOUNDS (certain (structural))

Both are C5 branched alkanals reached by the same chain-length surrogate from the C6 n-alkanal. Branch position is unmeasured anywhere in the corpus, so the layer CANNOT distinguish them. Any difference Hong measured between the two is, for this layer, pure unexplained residual by construction.

### dimethyl_disulfide -- carries an unresolved SOURCE CONTRADICTION (certain)

anantharamkrishnan2020b's Table 2 and its own Results text disagree on whether DMDS forms a covalent adduct. Resolved conservatively in favour of the table (no covalent term) and reported on the row rather than decided silently. It changes nothing numerically here, because the covalent term is inert anyway.

### What would falsify this expectation

If the ester's measured ratio turns out to be ABOVE 1 the sign call would pass and this pre-registration would be wrong in the model's favour. If several of the six no-constant compounds happen to sit within 5x of 1.0, the magnitude count could reach 7/10 while the sign clause still fails -- that would be a pass on one half of a criterion that requires both, and it must not be reported as a partial pass.

## Hold-out exposure disclosure

Amendment 2 set the precedent: exposure that occurred during directed reads is disclosed in the report, and its non-use is enforced by a literal-grep firewall test rather than asserted.

**Sealed successfully:**

- results/validation/holdout_frozen/hong2020_paired_thresholds.json was produced by a custodian agent under an explicit no-outcome reporting rule and was NOT read by this wave while the layer was built.
- data/lit/extraction_dossiers/hong2020_extraction.md sec. 4 (the paired BET table) was NOT opened; only the firewalled public manifest was used, and it carries structure, method and matrix description only.
- data/benchmarks/external_validation/ was never opened.

**Exposure incurred:**

- **Hong 2020 AGGREGATE statistics and four per-compound facts** (data/lit/extraction_dossiers/k4b_paired_thresholds_and_browning.md, which the build brief instructed this wave to read). k4b's summary sections print, outside the paired table: the range and geometric mean of the soy shifts; the log10 SD; the fact that one of ten inverts and that it is the ester; the hexanal soy/water value; and the largest and smallest elevated values with the chemical class of each. That is roughly four of the ten rows disclosed in whole or in part.
  - *Consequence:* The Hong hold-out is NOT fully blind for this wave and must not be described as such. It is scored TWO WAYS in the hold-out report: over all ten rows, and over the SIX rows for which nothing leaked. The six-row score is the one that carries evidential weight.
  - *Mitigation:* None of the exposed values enters a parameter, a bound, an initialisation, or a class assignment. Enforced by tests/unit/test_kinetic_core_b4.py::test_holdout_firewall, which greps the executable code of every B4 runtime and fit file for the exposed literals.
- **Brewer 1995's six beef thresholds and their ratios** (k2_matrix_and_thresholds.md sec. A.1 and sec. A.8, also an instructed read). Brewer is D.6 Module 7 HOLD-OUT (reclassified `dose_added_pre_cook`). Its numbers were seen.
  - *Consequence:* The beef alkenal/alkanal ratio (2.01x) is the observation that would pull the unsaturation penalty DOWN toward k2's stated 2-3x band. It is EXCLUDED from the fit, which leaves the fitted penalty above that band. Excluding it is the conservative direction for the fit and the honest one for the hold-out.
  - *Mitigation:* see consequence

## Named terms NOT implemented, and why

- **lipid_phase_partition** -- Leksrisompong's oil arm (FIT: K ratios) measures a 170x suppression for a log P 3.4 lactone and a 1.85-1.94x ENHANCEMENT for two hydrophilic odourants over the same oil. A lipid term is therefore real, sign-switching, and fittable on FIT rows. It is NOT implemented here for two reasons: the build spec fixes this layer at THREE named terms, and the fat content of the one matrix it would be applied to (Hong's soy paste) is unreported, so it would be fitted on FIT rows and then evaluated at an invented loading. Leading candidate for part of the residual on the lipophilic members of the panel.
- **background_odour_masking** -- The only mechanism proposed anywhere in the corpus that is compound-specific in a way that IGNORES hydrophobicity, and the only one that can produce an INVERSION. It is not implementable: Hong measured no background volatile concentrations in the soy matrix, so there is nothing to compute a masking term from. This is the single largest named reason the layer is expected to miss an inversion.
- **delivery_kinetics** -- Baek 1999: perceived intensity tracks the RATE of release (Imax/2(t75-t25), r = 0.968) better than the maximum concentration (r = 0.860). No equilibrium partition model computes a time derivative at all. Would require a release-kinetics state this layer does not have.
- **criterion_bias** -- Vega's 75 % UNCORRECTED criterion yields a systematically HIGHER threshold than a 50 % forced-choice one. Sign known, magnitude not (~2x is a plausible size, k2 sec. D.2(i)). Not implemented because no source measures it; it is a known contributor to every cross-method ratio in the corpus.

## Matrices the layer refuses to tabulate

- **soy_paste_hong** -- Hong 2020's 10 paired soy/water thresholds are the GATING HOLD-OUT of declaration Amendment 4. Ingesting them as table entries would spend the hold-out. The layer PREDICTS this matrix; it does not look it up.
- **caseinate_1pct** -- Leksrisompong 2010's 24 BETs are not fit-eligible: Amendment 4 makes only its K RATIOS FIT, and k4b sec. H.2 proposes the BETs as hold-out. Its K ratios are used as binding constants; its thresholds are not used at all.
- **emulsion_10pct_fat** -- as caseinate_1pct; and fat is confounded with aqueous-phase protein by the paper's own design.
- **emulsion_20pct_fat** -- as emulsion_10pct_fat.
- **soybean_oil** -- as caseinate_1pct.
- **cooked_beef** -- Brewer 1995 is declaration D.6 Module 7 HOLD-OUT and RECLASSIFIED `dose_added_pre_cook`: its numbers are doses added to RAW beef before a 70 C cook, not concentrations present at the moment of perception. They are not thresholds in this repo's sense and are not tabulated here.
- **milk_tian** -- Tian 2020's units cell prints a literal `?`, verified at 900 dpi. A factor-of-1000 basis risk; declared `neither, until the unit is settled`.
- **saliva** -- Starkenmann's thiol binding is STRANDED (no basis, no stoichiometry, mechanism unresolved) and Baek 1999 EXCLUDES mucosal binding for a neutral ester. A single saliva factor applied to all odourants is refuted by the pair.

## Source contradictions carried, not resolved

- **dimethyl_disulfide_adduct** -- anantharamkrishnan2020b Table 2 lists dimethyl disulfide as `no` (unreactive) and is one of the 32. Its own Results text, pp. 12-13, says DMDS DOES form a covalent bond with beta-lactoglobulin and names the adduct (+46 Da, BLG-CysSSMe), while noting the extent is less than dimethyl trisulfide's. THE TABLE AND THE TEXT CONTRADICT EACH OTHER. This layer resolves it CONSERVATIVELY -- DMDS gets no covalent term, matching the tabulated gate -- and reports the contradiction on every DMDS prediction rather than picking a side silently.
- **meynier_partition_absolute_scale** -- Meynier's own measured air/water K sits 6.24x below the Henry's constants printed in its own Table I, systematically across five compounds in two chemical classes with a 1.6x spread. Leksrisompong lands 6-17x low for the same methodological reason (gas-phase sample quantified against a liquid-phase calibration). This is NOT a licence to swap the constant; it is why every absolute here carries +/-0.5 decades and why only RATIOS are fitted.
- **unsaturation_penalty_magnitude** -- k2 sec. D.3 puts the alpha,beta-unsaturation penalty at ~2-3x from two sensory matrices; Meynier's headspace contrast is 4.95x. k4b flags the disagreement explicitly: 'the term's magnitude is matrix-dependent'. The FIT-only estimate here (3.73x) therefore sits ABOVE k2's stated 2-3x band, and it does so because the one observation that would pull it down (beef, 2.01x) is HOLD-OUT and excluded.

