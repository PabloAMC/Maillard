# The reaction trees, step by step

*Appendix to the [introduction](INTRODUCTION.md). The three figures are the actual step lists the
model integrates, one per path, drawn from the code by `scripts/generators/build_reaction_tree.py`.
Each box is a molecule; each arrow is one step, coloured by how its rate constant is known. A box is
coloured only when the test panel measures that molecule, by how far the prediction is from the
measurement. Bookkeeping pools (fragment carbon, acid equivalents, oxidant) are hidden, and steps
whose only products are such pools are drawn into one box, "removed into the matrix": those are the
sinks. Steps the shipped model carries at exactly zero (the inert defaults of the two refused sink
variants) are not drawn. Updated 2026-09-09, after the pyrazine step joined the sugar path.*

![How each path's steps are known](../assets/thiol_sink/13_steps_by_status.png)

On the sugar path two thirds of the steps are known with their temperature dependence; the rest
are bands or brackets (the dicarbonyl sinks from a dry glass, the pyrazine condensations declared
fast). On the pentose–cysteine path only four steps are measured with a barrier; thirty-one have a
rate pinned at 145 °C by the fit but no measured temperature dependence, and twenty-five are kept
unchanged from earlier calibrations. That is the whole thiol problem in one bar: away from 145 °C
every constant is an extrapolation.

## The sugar path

![The sugar path, step by step](../assets/thiol_sink/10_tree_sugar.png)

## The pentose–cysteine path

![The pentose-cysteine path, step by step](../assets/thiol_sink/11_tree_sulfur.png)

Read it from the two thiols rightwards: every arrow leaving them, into the disulfides, the
matrix-bound forms and the sink box, is green, a rate fitted at 145 °C with no measured temperature
dependence. Those are the arrows the experiment in the introduction measures.

## The acrylamide path

![The acrylamide path, step by step](../assets/thiol_sink/12_tree_acrylamide.png)

## Other figures kept for reference

The dicarbonyl ordering in water against the model, and the intermediate's decay against the model,
both behind the introduction's section 5: `docs/assets/thiol_sink/06_dicarbonyls_water.png`,
`05_ttca_decay.png`; the 140 and 168 °C shapes: `02_wang2022_shapes.png`, `03_liu2023_168C.png`; the
literature funnel: `07_literature_funnel.png`.
