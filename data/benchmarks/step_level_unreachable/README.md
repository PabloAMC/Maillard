# Step-level measurements the reaction network cannot execute

Created **2026-08-28 (Wave X)**. Nine verified literature rows live here, and **none of them is
scored**.

## This is not quarantine

`data/benchmarks/quarantined/` holds values for which **no real source could be found**. This
directory holds the opposite case: every number here was read off a rendered page image of
`data/articles/hofmann1998.pdf`, cross-checked against the paper's own `mol %` column by an
independent arithmetic route, and is **true**. What is missing is a **reaction step in `src/`**,
or — for two of the rows — an **observable** and a **seed species**.

Files here are invisible to every panel-wide generator, report and headline count, because
`panel_bundles` in `src/kinetic_core/panel.py` uses a non-recursive `glob("*.json")`.
That is the same mechanism that keeps `external_validation/` and `quarantined/` out.

## Why not just score them and let them fail?

Because they would not fail. **`_mean_abs_log10_error` skips comparisons whose prediction is
non-positive**, so a row the model cannot produce at all scores as *no error*, and a structural
zero would be quietly excluded from the very aggregate that is supposed to notice it. Two of the
rows are worse than that: the mercapto-2-propanone rows produce no entry in `predicted_ppb` at
all, so `_best_prediction_match` returns no match. Scoring these rows would make the panel look
better for being less capable. They are therefore held out of the panel **and named**, which is
the honest version of the same information.

Each file carries a `not_executable` block with three fields: `blocker_class`,
`what_is_missing` (measured, with the function and line of code that stops it) and `step_needed`
(what would have to be written, with the exact balanced stoichiometry where one applies).

## The blockers, cheapest first

| # | blocker | rows blocked | what it needs |
|---|---|---|---|
| 1 | **A fed 3-deoxyosone has no consumer.** `Enolisation_1_2` is emitted only from inside `_enolisation_steps`, which requires an Amadori product. | Table 3 3-deoxyribosulose (1) | A refactor, not new chemistry: move the 1,2-enolisation dehydration into a pool-scanning template keyed on the deoxyosone's canonical SMILES. Both steps, both products and the barrier already exist. |
| 2 | **1-mercapto-2-propanone is not a projected observable.** The step runs and the compound enters the pool; `predict_from_steps` reports only compounds in the desirable/off-flavour/toxic registries. | Table 7 (2) | **One sourced odour threshold.** Deliberately not invented — see the file's own note. |
| 3 | **No C2 + C3 aldol.** The network fragments sugars and never condenses fragments back. Fed hydroxyacetaldehyde + 2-oxopropanal enumerates **zero** steps. | Table 6 (3) | The two branches Hofmann & Schieberle **draw and name** in their Figure 2, both exactly balanced with species the engine already has. Half-blocked by #1. |
| 4 | **No amine-free sugar degradation lane.** The whole sugar trunk is entered through `_amadori_cascade`; sugar + H₂S contains no amine, so the sugar is untouched. | Table 3 ribose, Table 4 ribose (2) | A caramelisation entry step. Not attempted in Wave X: it sits upstream of every sugar row on the panel and no measurement here separates the catalysed from the uncatalysed rate. |
| 5 | **Water is not a tracked species**, so a hydrolytic step cannot be the *first* step of a system. Thiamine alone emits seven steps and produces MFT; the flux propagator returns an empty dict. | Table 10 thiamin (1) | A seed, not a step — but not a local one: seeding water changes the `co_reactant_factor` of every water-consuming step on the panel. Needs its own wave and its own before/after table. |

Blocker #5 is the most general finding in the list and the only one that is not about sulfur
chemistry: **any system whose first step is hydrolytic predicts nothing at all** — not a low
number, nothing.

## Promotion

When a blocker is closed, the affected file moves up one directory into `data/benchmarks/` and
becomes a scored panel row **with no edit to its numbers**. That is the point of storing them
in benchmark format rather than as prose.

## Source

All nine rows: Hofmann, T.; Schieberle, P. *J. Agric. Food Chem.* **1998**, *46* (1), 235–241,
DOI `10.1021/jf9705983`, Tables 3, 4, 6, 7 and 10. Conditions for every row: the stated moles in
50 mL of aqueous buffer, 145 °C, 20 min, laboratory autoclave, at the pH printed in the table.
Precision for every row: Table 1 footnote b, *"Data are mean values of triplicates and differed
by not more than 10%"*, which every one of these tables references.
