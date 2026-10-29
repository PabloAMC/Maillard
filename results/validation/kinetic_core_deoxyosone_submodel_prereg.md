# Sub-model 3: the deoxyosones and their fates (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Do the pentose-derived dicarbonyls (1-deoxypentosone, 3-deoxypentosone, 1,4-dideoxypentosone),
norfuraneol and furfural form and disappear at the measured rates, given the route's supply?

**Decision.** Whether the unpinned branch constants (`k_dpo_ddp`, `k_dpo_nf`, `k_nf_decay`; B2.3) can be pinned on
literature data, and whether 1,4-dideoxypentosone, which no laboratory has measured, can stay a modelled species or
must be declared an assumption with a band.

**Rival arms.** Norfuraneol from 1-deoxypentosone (current) vs also from the Cys-Amadori compound directly; furfural
from 3-deoxypentosone only vs also from the pentose's acid-catalysed dehydration.

**Expected result.** Norfuraneol and furfural time courses pin their formation and loss (identified, ~65 %);
1,4-dideoxypentosone stays unpinned (~80 %), and its MFT route is reported as an assumption.

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 2 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2) -- being filled by the first reader (data/lit/source_verification/deoxyosone_submodel.yml,
data/lit/fit_rows/deoxyosone_submodel.yml, status: draft)
## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
