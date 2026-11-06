"""
src/kinetic_core_jax -- the engine's lanes in JAX + diffrax (plan step A2 + H3; docs/guides/A2_MODULE_CONTRACT.md
v2). NOT part of the engine: nothing in ``src/kinetic_core`` imports it, it needs the separate ``maillard-jax`` env
(environment-jax.yml), and no shipped number comes from it.

Two layers, side by side:

* the A2 contract v2 runtime, ``network.py`` (Module, compose, integrate, the batch API, the exact Jacobian), and its
  modules: ``lipid.py`` (the lipid lane, the first module; 2026-11-05; results/validation/jax_lipid_module.{json,md}),
  ``trunk.py`` (``network.REACTIONS``, the trunk the acrylamide lane runs, and the trunk lane's declared terms on it)
  and ``acrylamide.py`` (the acrylamide block and its condition terms; composed after the trunk it is the acrylamide
  lane; 2026-11-05; results/validation/jax_trunk_acrylamide_module.{json,md}), ``trunk_lane.py`` (the trunk LANE's
  seven blocks). The network infrastructure (contract v2.2; 2026-11-05; results/validation/jax_network_infra.{json,md}):
  ``forcing.py`` (forced inputs: measured pH / temperature / a_w / species courses), ``fit.py`` (the fit driver:
  least squares on the exact Jacobian) and ``densities.py`` (per-draw predictions for M1.1's predictive densities);
  ``network.compose`` carries the declared-competition and titratable-table hooks. The sulfur lane (2026-11-05;
  results/validation/jax_sulfur_lane_module.{json,md}): ``sulfur_lane.py`` (one module, ``sulfur.FULL_REACTIONS``,
  with its charge balance and both pH couplings in its own driver, ``sulfur_lane.integrate``) and
  ``sulfur_structure.py`` (numpy only: the structure it derives, the hand-ported registry and the check that the
  NumPy lane has not moved under it, run by the core env's unit tier).
* the sulfur-lane BAKE-OFF port (STATUS_AND_PLAN W4 step H2, pulled forward; 2026-11-04), ``lane.py`` / ``model.py``:
  the lane's 62 species and 132 reactions, its rate laws and pH factors, and its inline charge balance (in-situ pH
  and ionic strength at the reaction temperature, and the cooled 25 C label pH the measured pH tables read), with the
  pH solved inside the right-hand side and differentiated by the implicit-function theorem (``model.make_ph_solver``).
  Scope: the configuration of the R3En joint fit only; ``lane.build`` refuses any other switch setting. It is kept as
  the bake-off's evidence (results/validation/jax_bakeoff_sulfur_lane.{json,md}); the contract module
  ``sulfur_lane.py`` reuses its pH solver.
"""
