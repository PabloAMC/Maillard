# Sub-model 2: the cysteine-pentose route (pre-registration, in progress)

*Strategy 2026-10-29, Phase 2. Written before any fit, following `docs/guides/OPERATING_PROCEDURE.md`; started in parallel
with sub-model 1 (owner, 2026-10-29: "keep doing what you can"). The unclear steps it must test are in
`results/validation/sulfur_reaction_evidence.md`.*

## 1. Frame (step 1)

**Question.** Does de Roos's measured route (TTCA <-> Cys + pentose; TTCA -> cysteine Amadori compound -> 1- and
3-deoxypentosone + Cys) reproduce the time courses of TTCA, the Cys-Amadori compound, free cysteine and the pentose
across the three laboratories that measured them, and does it replace the shipped step that treats cysteine as an
instant catalyst on the pentose (`r_pent_tdp`, `r_pent_dpo`: 59 % and 45 % of a reference pot's flux, unpinned since
B2.3; `docs/reference/unpinned_constants.yml`)?

**Decision.** Whether the route ships (switch on) and the two catalyst steps retire; which reading of phosphate
catalysis vs pH the data support (stage 1's 1.87x Cao-vs-Zhai tension, CF-011).

**Rival arms.** (a) Phosphate catalyses the Amadori steps (Davidek 2002's 1.3-3.5x) vs (b) the difference is pH
alone vs (c) laboratory. Above pH 6 the decay is extended or flat (CR-14).

**Expected result (before any probe).** The route reproduces Cao 2017 and Zhai 2019's intermediates within the
measurement errors (stage 1 got there with the procedure); the phosphate question stays open on literature data
(~60 %), and the designed experiment's buffer arm is what settles it.

**Done-test** (the operating procedure's mechanism checks, as in sub-model 1): its own species, intermediates included,
across laboratories by leave-one-laboratory-out under the measurement model; no unexplained structural residual flag;
every fitted constant identified and off its bounds, consistent with independent measurements of the same step;
element balances where totals are measured; one out-of-condition prediction frozen before it is scored; the rival
arm below fitted beside the main one. Fits run after sub-model 1 ships its constants (the pH trajectory and the
upstream supply it predicts are inputs here); everything before the fit runs now.

**Budget.** See the strategy's "Next steps" in hours.

## 2. Data (step 2) -- being filled by the first reader (data/lit/source_verification/route_submodel.yml,
data/lit/fit_rows/route_submodel.yml, status: draft)
## 3. Cheap probe (step 3) -- TO DO
## 4. Design (step 4) -- TO DO
## 5. Design review (step 5) -- TO DO
## Amendments
