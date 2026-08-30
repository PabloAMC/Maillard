# Build Wave B6 — PRE-REGISTRATION (the lipid-oxidation module)

**Written 2026-08-29, BEFORE `src/kinetic_core/lipid.py`, `parameters_lipid.py` or
`species_lipid.py` existed, before any fit was run, and before the exam was
re-generated.** Nothing below may be edited after the fit report is written; the
outcome report (`kinetic_core_b6_holdout_report.md`) scores against *this* text.

Module 5 of `docs/reference/FIT_HOLDOUT_DECLARATION.md` (D.6):

| dataset | role |
|---|---|
| Frankel 1989, the 26-column distribution **at zero additive** | **FIT** |
| Frankel 1989, the α-tocopherol arms | ★ HOLD-OUT (two-sided) |
| Frankel 1989, the nonanal ABSENCE | HOLD-OUT (negative test) |
| the 2 frozen Bi 2020 matrix-path bundles | HOLD-OUT (frozen; exam generator only) |

---

## 0. DISCLOSURE, FIRST, BECAUSE IT CHANGES HOW EVERY SCORE BELOW READS

**The builder of this wave has SEEN the α-tocopherol hold-out columns.** This is
not a slip that could have been avoided by better discipline: Frankel 1989's
Table 1 and Table 2 print the zero-additive column (FIT) and the tocopherol
columns (HOLD-OUT) *in the same table rows*, and the paper's ABSTRACT states the
hold-out result in prose ("volatile formation decreased in all tests … the
relative amounts of pentane and methyl octanoate to decrease and hexanal and
methyl 9-oxononanoate to increase"). There is no way to read the FIT column out
of `data/articles/frankel1989.pdf` without reading the hold-out column beside it.

Consequences, adopted here and binding on the outcome report:

1. **The α-tocopherol hold-out is scored `seen_diagnostic`, never `gating`.**
   It follows the precedent of Amendment 9 clause 1 (the Kang 140 °C row) and
   the Amendment 2 disclosure. A pass is *not* evidence of blindness.
2. **The mitigation is structural, not procedural: the donor term has NO FITTED
   PARAMETER.** The prediction registered in §3 is a *monotonicity theorem* that
   holds for every value of the donor-suppression parameter `d ∈ (0, 1)`. There
   is nothing to tune toward the seen numbers, and a unit test asserts the
   monotonicity over the whole range rather than at a fitted point.
3. A **literal-grep firewall test** asserts that none of the hold-out numbers
   appears anywhere in `src/kinetic_core/*lipid*.py` or in the B6 fit report.
4. The two frozen **Bi 2020** bundles were NOT opened. Their directory was never
   read by this wave; they are reached only through
   `scripts/generators/generate_cutover_final_exam.py`.
5. The **nonanal absence** hold-out is unaffected by (1): it is a *structural
   zero*, and the architecture that produces it (nonanal has no edge from any
   linoleate hydroperoxide, because nonanal is the C9 fragment of the OLEATE
   double bond) is fixed by molecular topology, not by any number in the paper.
   It stays **gating**.

---

## 1. WHAT IS BEING BUILT

A `LIPID` lane in `src/kinetic_core`, structured as:

```
lipid pool (declared: lipid fraction x FA profile)
    -> LOOH pools, resolved by POSITION (9- / 13-) and GEOMETRY (cis,trans / trans,trans)
        -> per-(position, geometry) branch simplex over the MEASURED Frankel slate
            13-OOH -> {pentane, hexanal, methyl 13-oxo-9,11-tridecadienoate}
             9-OOH -> {methyl octanoate, 2,4-decadienal, methyl 9-oxononanoate}
        -> consumption: covalent aldehyde-lysine channel (bounded ceiling, INERT by
           Amendment 6 ruling 2) + B4's reversible-binding / threshold layer
```

**The rate architecture, stated before it is built.** The branch DISTRIBUTION is
measured (Frankel 1989, 180 °C). The absolute LOOH decomposition RATE at cooking
temperature is measured NOWHERE — declared gap `research_round3_channels.md`
§F.3, which re-affirms `k3` §C.9. Therefore:

* `k_LOOH_decomp` enters as an **explicitly-labelled bounded input**, anchored
  at Schroën & Berton-Carabin 2022 Table 1 `k4 = 6·10⁻³ h⁻¹` at **25 °C**, pH
  6.7, rapeseed O/W emulsion, **hand-fitted by visual agreement, no standard
  error anywhere in the paper**, and **lumped** (all secondary products).
* Its temperature dependence is an **UNMEASURED assumption**. `Q10` is a
  user-visible parameter with a **declared default** and a **band [2, 3]** taken
  from the authors' own prose licence — applied over a ΔT the authors explicitly
  did not licence (their words cover "adjustment", not 11.5 decades of 10 °C
  steps). **No Q10 number is baked into any rate constant.**
* Every lipid prediction returns `in_envelope_extrapolated`, never
  `in_envelope`, and carries the extrapolation warning.
* **Ratios across formulations are first class** (the rate assumption cancels).
  **Absolute ppb inherits the rate band** into B4's interval, computed by
  RE-INTEGRATING at Q10 = 2 and Q10 = 3 rather than by adding a nominal width.

---

## 2. THE FIT — pre-registered, with its own expected failure

**Fit rows: the three zero-additive columns ONLY** (Table 1 col. "0"; Table 2
`cis,trans`-13 col. "0"; Table 2 `trans,trans` 9+13 col. "0") — 18 numbers, 3
simplex constraints, 15 independent. **Nothing else in Frankel 1989 is read by
the fitter**, and the fitter's input array is asserted to have exactly 18 cells.

Free parameters (12): four 3-way branch simplexes (13-OOH and 9-OOH × ct and tt
geometry) = 8 free; per-system LOOH composition = 4 free (`f13` for each of the
three systems, plus a `ct/tt` mix fraction for the mixed system). Residual df = 3.

**PRE-REGISTERED EXPECTATION F-1 — the fit will NOT be exact, and the reason is
already visible in the fit column itself.** Pentane and methyl
13-oxo-9,11-tridecadienoate are the two halves of ONE β-scission and should be
1 : 1. Their measured ratio across the three zero-additive columns spans
**8.5×**. Either the GC response/recovery factors differ enormously between a C5
alkane and a C14 dienoate (Frankel's own method traps at −65 °C and his
Discussion says the dienals "may be more reactive and less stable"), or the
1 : 1 pairing is wrong. **The model does NOT impose the pairing** — it fits free
shares inside each isomer's simplex — and the wave reports the falsified pairing
as a finding rather than absorbing it into a response factor. Predicted
composition residual: **median |Δshare| ≤ 3 percentage points, worst ≤ 8 points.**

**PRE-REGISTERED EXPECTATION F-2 — the shipped `hexanal 0.37` will be refuted.**
The fitted 13-OOH → hexanal share will land inside Frankel's measured
zero-additive envelope (11–20 % of the six-product sum), i.e. **well below
0.37**, confirming `k3` §A.5's verdict that 0.37 "sits ABOVE the paper's entire
measured range". Confidence ~99 % (arithmetic, not a forecast).

**PRE-REGISTERED EXPECTATION F-3 — `nonanal 0.15` will be refuted structurally,
not numerically.** Nonanal receives NO edge from any linoleate hydroperoxide.

---

## 3. HOLD-OUT 1 — the α-tocopherol two-sided signature (`seen_diagnostic`)

The mechanism assumption, stated so it can be attacked: **a hydrogen donor
quenches the alkoxyl radical before homolytic β-scission, and does not block the
heat-promoted heterolytic (Hock) cleavage.** Basis: Frankel 1989's INTRODUCTION,
which attributes the homolytic/heterolytic product assignment to refs 3–10 (all
pre-1989), plus general radical chemistry. This assumption is NOT taken from the
tocopherol arms and is not fitted to them.

Implementation: one donor-suppression parameter `d ∈ [0, 1)` multiplying the
homolytic channel only. **`d` is never fitted; no value of `d` is stored.**

Registered predictions (all must hold for EVERY `d ∈ (0, 1)`, which is why they
are theorems and not tunings):

| # | prediction | status |
|---|---|---|
| H1-a | **total volatile flux DECREASES** monotonically in `d` | gating on the machinery (unit-tested over the whole range) |
| H1-b | **hexanal SHARE INCREASES** monotonically in `d` | gating on the machinery |
| H1-c | the two move in **OPPOSITE** directions — the two-sided signature a fixed split cannot fake | gating on the machinery |
| H1-d | methyl 9-oxononanoate share **increases**; pentane and methyl octanoate shares **decrease** | `seen_diagnostic` vs Frankel |
| H1-e | methyl 13-oxo-tridecadienoate and 2,4-decadienal shares **decrease in step with pentane** | **registered as EXPECTED WRONG** |

**H1-e is registered as a predicted FAILURE.** The architecture treats those two
as homolytic co-products, so their shares must fall with the donor. Frankel's
Discussion states no trend was observed for them, and his Fig. 4 excludes them
from both the numerator and the denominator. Registering this now means the
outcome report cannot present it as a discovery.

**What this hold-out CANNOT establish:** any magnitude. The wave predicts signs
and monotonicity only. It will NOT report a fold-error against a tocopherol
column, because the donor loading (wt %) has no mapping onto `d` that any source
licenses, and inventing one is exactly the failure mode this module exists to
avoid.

---

## 4. HOLD-OUT 2 — the nonanal ABSENCE (gating, genuinely blind)

**Prediction: exactly 0.0 for nonanal from a pure linoleate-hydroperoxide feed**
— a structural zero, not a small number. Enforced by topology: `NONANAL` has no
incoming edge from `LOOH_13_*` or `LOOH_9_*`; its only edge is from an OLEATE
hydroperoxide pool. Unit-tested as an exact zero, and asserted at import by the
network's structural validator.

**The other half of the same test, registered here so it is not quietly skipped:**
in a REAL matrix (soy, pea) the oleate pool is NOT zero, and **the oleate → nonanal
branch fraction is measured nowhere in the fit corpus.** The engine will therefore
**still refuse absolute nonanal in a real matrix**, with a new and sharper reason.
Confidence that this is the right call: high. A wave that answered nonanal by
carrying `0.15` forward would score better on the exam and be wrong.

---

## 5. HOLD-OUT 3 — the exam (2 frozen Bi 2020 bundles + the wider 8 refusals)

Baseline, from `results/validation/cutover_final_exam.json` at HEAD 7f65cca:
**17 refusals of 40 points**, of which **8 are the matrix/lipid rows**:

| # | bundle | compound | current refusal |
|---|---|---|---|
| 1 | `..._bi_2020_raw_pea_hexanal` | hexanal | unmapped precursor + no lipid path |
| 2 | `..._bi_2020_roasted_pea_hexanal` | hexanal | unmapped precursor + no lipid path |
| 3 | `..._li_2026_spi_wg_hme_control` | hexanal | " |
| 4 | `..._li_2026_spi_wg_hme_control` | nonanal | " |
| 5 | `..._li_2026_spi_wg_hme_control` | 2-pentylfuran | " |
| 6 | `..._li_2026_spi_wg_hme_control` | 1-hexanol | " |
| 7 | `..._liu_2023_ppi_offnote_baseline` | hexanal | " |
| 8 | `..._liu_2023_ppi_offnote_baseline` | nonanal | " |

**PRE-REGISTERED EXAM OUTCOME:**

| row | registered outcome | why |
|---|---|---|
| 1–3, 7 (**hexanal ×4**) | **ANSWERED**, `in_envelope_extrapolated`, with an interval spanning **≥ 2 decades** | the branch fraction is measured; the rate and the LOOH pool are declared assumptions |
| 4, 8 (**nonanal ×2**) | **STILL REFUSED**, new reason: *the oleate → nonanal branch fraction is unmeasured* | §4 |
| 5 (**2-pentylfuran**) | **STILL REFUSED**, new reason: *no measured branch fraction; the alkylfuran route is not in the Frankel slate* | inventing one is prohibited |
| 6 (**1-hexanol**) | **STILL REFUSED**, new reason: *no aldehyde-reduction step is measured anywhere in the corpus* | unchanged in substance, sharper in wording |

So the registered headline is **17 refusals → 13**, not 17 → 9. **A wave that
un-refuses all eight has invented three branch fractions.**

**PRE-REGISTERED ACCURACY EXPECTATION E-1 (the real test).** The four hexanal
point predictions will be **BADLY WRONG at the point value** — the prior for this
lane is the repo's own documented 36× / 1304× hexanal over-prediction
(`matrix-path-accuracy-plan`, Amendment 6 ruling 2). Registered:

* **point fold-error > 3× on at least 3 of the 4** — confidence ~85 %;
* **direction: OVER-prediction on at least 3 of the 4** — confidence ~75 %;
* **the measured value falls INSIDE the reported interval on at least 3 of the
  4** — confidence ~60 %. This, not the point, is the claim the interval makes.

**If the point predictions come out close, that is a warning sign, not a
success**, because none of the three quantities that set the absolute scale
(LOOH pool, Q10, the 25 °C anchor) is measured in these systems.

**PRE-REGISTERED REGRESSION GUARD G-1.** The 23 currently-ANSWERED exam points
must be **byte-identical** after B6. The lipid lane is additive; if any sulfur,
acrylamide or trunk number moves, the wave has broken something and the report
says so instead of re-baselining.

---

## 6. LANE RESOLUTION — the decision, registered before it is exercised

The spec asks for a decision on requests spanning the Maillard and lipid lanes.
**Decision: CO-INTEGRATION, as a direct sum.**

Justification, checkable: the lipid network's species set (LOOH pools, the six
Frankel products, nonanal) is **disjoint** from the trunk / sulfur / acrylamide
species sets — no shared sugar, amine, sulfur or carbonyl pool. The one candidate
coupling is the **aldehyde–lysine covalent channel**, which would put the lipid
aldehydes and the Maillard amine pool in competition. That channel is **INERT BY
RULING** (Amendment 6 ruling 2: it supplies ~0.06 % of the hexanal log-shift, and
its Ea is unmeasured), so it carries zero flux and creates no coupling.

Registered as a **conditional**: a runtime guard asserts the covalent channel's
contribution is zero, and **co-integration must be revisited the moment that Ea
is measured**, because the amine pool then becomes genuinely shared. The guard
raises rather than silently summing.

This is a different ruling from the acrylamide/sulfur `LANE CONFLICT`, and the
difference is the point: those two lanes both consume cysteine, so summing them
spends it twice. These two share nothing.

---

## 7. WHAT THIS WAVE WILL REFUSE TO DO

* No Q10 number is written into any rate constant.
* No branch fraction is invented for 2-pentylfuran, 1-hexanol, nonanal-from-oleate,
  2-nonenal or methyl 12-oxo-10-dodecenoate.
* No activation energy is fitted to Frankel 1989 (one temperature; prohibited by
  `k3` §C.9 and re-affirmed by §F.3).
* The shipped `hexanal 0.37` / `nonanal 0.15` in
  `data/lit/lipid_oxidation_calibration.json` are **not** edited by this wave —
  that file belongs to the FAST lane. B6 declares its own registry and reports
  the contradiction.
* Nothing under `data/benchmarks/external_validation/` is opened.
