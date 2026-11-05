"""
src/kinetic_core_jax -- the JAX bake-off port of the sulfur lane (STATUS_AND_PLAN W4 step H2, pulled forward;
2026-11-04). NOT part of the engine: nothing in ``src/kinetic_core`` imports it, it needs the separate
``maillard-jax`` env (environment-jax.yml), and no shipped number comes from it.

What it is: the sulfur lane's 62 species and 132 reactions (stoichiometry and mass-action layout read from
``src.kinetic_core.sulfur``), its rate laws and pH factors, and its inline charge balance (in-situ pH and ionic
strength at the reaction temperature, and the cooled 25 C label pH the measured pH tables read), re-implemented in
JAX with the parameters as a pytree, integrated with diffrax. The pH is solved inside the right-hand side by a
bracketed Newton solve whose derivative is the implicit-function theorem (``chemistry.solve_ph``), so exact
gradients flow through it.

Scope: the configuration of the R3En joint fit (B60 + B58 / B59b / B59a / B61 on, inline coupling, the pentose +
amine entry, the DDP law "none" or "base"); ``lane.build`` refuses any other switch setting rather than computing a
different model silently. Results: results/validation/jax_bakeoff_sulfur_lane.{json,md}.
"""
