# Sub-model 4: thiol formation and loss (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Do MFT and FFT form from their precursors (norfuraneol / dideoxypentosone and furfural with H2S) and
disappear (decay, dimerisation, binding) at the rates the new time series show, and is H2S supply and loss right?

**Decision.** Whether the unpinned thiol constants (`k_fur_fft_hs`, `k_mft_decay`, `k_ddp_mft`, `k_ddp_mft_hs`,
`k_h2s_loss`; B2.3) can be pinned; how the laboratories' rising-vs-falling thiol curves (CF-013) are explained.

**Rival arms.** Thiol loss as first-order decay vs dimerisation (second order in thiol) vs binding to melanoidins
or protein (Gigl 2021); semi-quantitative series (the Zhai / Kang ladder, CF-004) read by shape only vs the absolute
ones (Feng 2022, isotope-dilution data) alone.

**Expected result.** Formation and loss of FFT and MFT become identified on the time series (~55 %); the H2S budget
stays the weakest link (~60 % unpinned), and the binding question needs the matrix layer (Phase 4).

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 3 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2) -- being filled by the first reader (data/lit/source_verification/thiol_submodel.yml,
data/lit/fit_rows/thiol_submodel.yml, status: draft)
## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
