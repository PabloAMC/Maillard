# A2 module contract (v2, 2026-11-05; draft v1 2026-10-31)

This is the interface through which chemistry modules contribute to ONE composable network, with one frozen configuration and one integrator. It is the target for A2 + H3 (one composable network, compiled with exact gradients) in [STATUS_AND_PLAN.md](STATUS_AND_PLAN.md), and for every standalone module that investigation agents write before the network exists. Sec. 1 is v1's read-only survey of `src/kinetic_core/` (2026-10-31; file:line references are to that date).

**What v2 changes, and why** (the sulfur-lane bake-off, [jax_bakeoff_sulfur_lane.md](../../results/validation/jax_bakeoff_sulfur_lane.md), and the first module port, [jax_lipid_module.md](../../results/validation/jax_lipid_module.md)):

1. **The tool is JAX + diffrax** (float64, CPU first). The contract is now code: [`src/kinetic_core_jax/network.py`](../../src/kinetic_core_jax/network.py) (sec. 2).
2. **Static structure, dynamic numbers.** A module is built from a `Config` into numpy arrays plus one pure rate function; its free numbers are a pytree `theta` that JAX differentiates.
3. **Pot data are arguments**, never constants traced into the compiled code, so one compile serves every pot of a configuration (and every corner and draw).
4. **The pH root goes inside the right-hand side**, differentiated by an implicit-function-theorem custom JVP; no DAE state (sec. 2.4).
5. **Forward mode** for fits of about 20 parameters; the discrete adjoint only above about 50; never the continuous adjoint (sec. 2.7).
6. **Acceptance against today's lanes** is rates to 1e-12 relative and trajectories to 1e-6 dex, not bit-identity, which cannot hold across integrators (sec. 3).
7. **Kvaerno4 vs Kvaerno5: neither dominates.** The solver is a declared choice per composed network (`compose(..., solver=)`), benchmarked on that network's pots, and both are reported for every lane port (sec. 2.6).
8. **Fits run at rtol 1e-6** once a residual check shows the panel moves < 1e-3 dex (sec. 2.6).
9. **A refusal on any unported switch**: a module either implements both values of a switch it reads, or raises `UnsupportedConfiguration` (sec. 2.9).
10. **Measured pH tables with jumps** keep their values exactly; their derivatives miss the crossing term, and pots that cross are flagged (sec. 2.4).
11. **A batch API** (`batched`, `grid`: many draws and pots in lockstep under `jax.vmap`), designed for a later CPU-vs-GPU benchmark (sec. 2.8).
12. **The cases' needs from the network** (M1, M1b, M1c, M3 - M6) are listed with their hook and status (sec. 4).

## 1. Why (what the v1 survey found)

- **Three trunks.**
  - `network.TRUNK_REACTIONS` (`network.py:520`);
  - sulfur `REACTIONS + FED_3DEOXY + SULFUR_REACTIONS` (`sulfur.py:1009`);
  - acrylamide `REACTIONS + ACRYLAMIDE_REACTIONS` (`acrylamide.py:268`).

  The name `TRUNK_REACTIONS` means a smaller set in sulfur / acrylamide than in `network.py`. Three copies of the balance checker, of the stoichiometry builder and of the `r_amadori` derivation exist. Two functions are named `ph_factor`, with different signatures (`trunk_conditions.py:156`, `parameters_sulfur.py:992`).
- **About 30 module-level switches** (sulfur pH coupling, `ph_state.B51`, `TRUNK_SOLUTES_IN_BALANCE`, 13 sulfur dataclass switches, lipid, dicarbonyl, panel, uncertainty). They are set by about 97 assignments in 39 scripts and about 50 in tests; none is set in `src/`. One joint-fit job sets about 15 of them across 3 files. (A3 phases 1 - 2 gave them a `Config`: sec. 6.)
- **Parallel correctness depended on forking** a new pool per call (`parallel.py`). A spawn-started or persistent worker would silently compute with the DEFAULT switches. (Fixed by A3.)
- **Charge closure has a gap.** With `TRUNK_SOLUTES_IN_BALANCE=True`, Gly / AMA / SB / PAMA titrate in `titratable_inventory` (`ph_state.py:653`) but are absent from `TITRATABLE_CENTRES` (`:1069`). The closure check (`validate_charge_closure`, `:1305`) therefore never sees trunk steps that move them. The ledger (validation) and the inventory (physics) are two tables for one job.
- **Amine identity is a convention.**
  - `Gly` stands for "the amine": MET / PRO are folded into Gly on the trunk and popped elsewhere (`engine.py:~2553-2570`).
  - Lysine has two keys (`Lys` in acrylamide, `LYSP` in the trunk).
  - Cys is consumed by two lanes, which `resolve_lanes` therefore refuses to compose.
- **Only the sulfur lane has dynamic pH.** The trunk scales k_ref at the declared pH, acrylamide has no pH argument, and lipid is pH-free.
- **(v2) The lipid lane is not an ODE in the engine.** `lipid.integrate_lipid` solves it analytically per segment. Its arithmetic is the exact solution of a linear mass-action network, so the port declares that network (sec. 3).

## 2. The contract (v2: JAX + diffrax)

### 2.1 Three splits

- **Static vs dynamic.** Species, reactions, stoichiometry, mass-action layout and which switches are on are numpy, frozen when a module is built from a `Config`. Every free number is a leaf of `theta = {module: {coordinate: value}}`, at its shipped value from that `Config` by default.
- **Network vs pot.** A pot's numbers are arguments of the compiled functions: its initial state, its thermal program and, for pH modules, its SID and buffer. One compile serves every pot whose arrays have the same shapes. The bake-off reran four stiff corners with no recompile. The lipid port runs 24 pot records (8 panel pots, each at its point and two corners) through one compile.
- **Time.** Each segment of a thermal program is integrated on a unit interval s in [0, 1], with the right-hand side scaled by the segment's duration. Pots whose segments differ in length therefore share one compiled solve, and run in lockstep under `vmap`. A zero-duration segment is the identity, so padding equalises segment counts.

### 2.2 `Module`

```python
@dataclass(frozen=True)
class Module:                       # src/kinetic_core_jax/network.py
    name: str
    species: tuple[str, ...]                                 # owned species
    elements: Mapping[str, tuple[int, int, int]]             # owned species -> (C, N, S)
    reactions: tuple[str, ...]
    stoichiometry: Mapping[str, Mapping[str, float]]         # reaction -> {species: net coefficient}; keyed by NAME
    reactants: Mapping[str, Mapping[str, int]]               # reaction -> {species: mass-action order (1 or 2)}
    params: Mapping[str, float]                              # free coordinates at their SHIPPED values
    bands: Mapping[str, tuple[float, float] | None]          # from config._BANDS / the engine's own bands; None = none declared
    basis: Mapping[str, str]                                 # coordinate -> dossier / config reference
    rate_constants: Callable[[Mapping[str, Array], RateCtx], Array]   # pure JAX; this module's reactions, in order
    observables: Mapping[str, Callable[[Array, Mapping[str, int]], Array]]
    requires: tuple[str, ...] = ()                           # modules whose species this one may touch
    needs_ph: bool = False
    config_digest: str = ""                                  # Config.digest() the module was built from
    notes: tuple[str, ...] = ()

class RateCtx(NamedTuple):          # what a rate function may read besides theta
    t_c; t_k                        # T_K = T_C + 273.15, as the engine computes it
    ph = None; ph_label = None; ionic_m = None               # pH modules only (sec. 2.4)

def build(cfg: Config) -> Module    # one per module file, e.g. src/kinetic_core_jax/lipid.py
```

Rules for a module:
- **Read the engine, never retype it.** Species, carbon counts, stoichiometry, parameter values, tables and bands are read from `src.kinetic_core` at build time. A literal that has to be re-typed, such as a gas constant inside an engine method, is named and pinned by the rate test.
- **Stoichiometry is static.** A yield that depends on a free coordinate becomes a *yield channel*: a reaction `A -> A + P - C(P) FRAG_C`, whose rate constant carries the yield. The lipid module carries Frankel's branch shares and the donor suppression this way.
- **Rate functions are pure.** They read `theta` and `RateCtx` and nothing else, with no module globals.

### 2.3 `compose`

`compose(modules, layout=None) -> Network`:
- orders the species (module order, or a pinned `layout`);
- builds one stoichiometry matrix, the padded reactant index and power arrays, and the element matrix;
- concatenates the modules' rate functions in module order;
- flattens `theta` in `Network.coordinates` order (`theta_vector` / `theta_from_vector` / `stack_theta`).

**It refuses** (ValueError):
- a duplicate module name, or a `requires` that is not composed;
- a species or reaction key defined twice;
- a reaction that touches a species no module owns, outside its `requires`;
- a reaction that does not balance C / N / S to 1e-9.

**Not yet implemented** (sec. 4 says which case needs each):
- the merged titratable table and its ONE charge-closure check;
- declared competition, so that a species consumed by two modules that do not declare it is refused;
- the forced-input driver.

### 2.4 The pH coupling (from the bake-off)

- **The root is solved inside the right-hand side.** The bake-off solves `F(u; c) = 0`, with u = (pH, I) and c the titratable inventory, at every right-hand-side call: a bisection (the residual is monotone in pH), then Newton on the coupled (pH, I). It solves twice, at the reaction temperature and at the cooled 25 C label. Its derivative is a custom JVP: du = -(dF/du)^-1 (dF/dc) dc, a 2 x 2 solve, and reverse mode transposes it with no extra code ([`model.make_ph_solver`](../../src/kinetic_core_jax/model.py)). There is no DAE state: diffrax has no DAE solver, and a relaxation would change the model.
- **Modes:** `none` (pH-free modules; the lipid module), `inline` (the above), `prescribed` (a measured pH trajectory read as a forced input; M1c: Cao 2017, Kang 2026), and `clamped` (pH-stat).
- **Status.** A network with any `needs_ph` module is refused by `integrate` until the sulfur module is ported, and that port reuses the bake-off's solver.
- **Cost.** The root-find was ~75 % of the bake-off's forward solve (2.1 vs 0.5 ms per step). It gets one focused step: a cheaper residual (a precomputed pK table, fewer exps), and one root per step reused across stages.
- **Measured pH tables with jumps** (Zheng & Ho's mechanism boundary at label pH 5.0, the thiamine share's floor at 6.25, the lipid breakage window at 100 C in temperature) are ported exactly with `jnp.where`. An exact derivative omits the jump-crossing term, which a finite difference across the boundary would see. Rule: a module whose table has a jump declares it, and the fit driver reports every pot whose trajectory crosses one, so a Jacobian there is read as one-sided. No compared column crossed a boundary in the bake-off or in the lipid port: the lipid window switches between segments, never inside one.

### 2.5 The thermal program

- **Segments.** `program_segments(segments, n_segments)` turns the engine's `ThermalProgram.segments` (`((duration_min, T_C), ...)`, B64's staircase for profiles) into an (S, 2) argument, padded with zero-duration segments.
- **Save times.** `refine_segments(segments, save_times)` splits segments so that every save time is a boundary. `integrate` returns the state at the end of every segment.
- **Rate constants** are evaluated once per segment when no module needs pH. A pH module evaluates them in the right-hand side.
- **Not yet ported:** the sulfur lane's PROFILE mode, which integrates T(t) piecewise-linear with its rate constants at T(t). It needs `ClipStepSizeController(jump_ts=breakpoints)`.

### 2.6 Integrator and tolerances

- **The solver is declared per composed network** (`compose(..., solver="Kvaerno5")`, the default; `integrate(..., solver=)` overrides it). Neither solver dominates:
  - On the bake-off's sulfur pot Kvaerno4 matched the engine better at every tolerance (2.7e-9 vs 5.9e-9 dex on the fitted rows at rtol 1e-8). It rejected almost no steps (2 vs 173) and ran ~25 % faster.
  - On the lipid pots Kvaerno5 takes ~5x fewer steps once a hydroperoxide pool is depleted far below atol (160 C, 30 min: 158 vs 794 steps; panel total at rtol 1e-8: 973 vs 4216). Kvaerno4 is slightly closer to the engine there (3.6e-8 vs 1.1e-7 dex at rtol 1e-8, both on a depleted pool).
  
  So each lane port reports both, and the union network's choice is made on the union's own pots.
- **Tolerances.**
  - References and acceptance tests run at rtol 1e-10 / atol 1e-16.
  - Fits run at rtol 1e-6, once a declared check shows the panel's residuals move < 1e-3 dex between rtol 1e-6 and 1e-8. On the sulfur pot the rtol 1e-6 Jacobian differed from the 1e-8 one by <= 4e-6 relative. On the lipid pots the rtol 1e-6 Jacobian differed from the rtol 1e-10 one by <= 1.7e-5 relative, and the rtol 1e-6 trajectories of the scored products by <= 3.6e-6 dex.
  - atol is absolute in mmol/L, so it is set per lane from the smallest species that is scored.
- **Per-pot failure.** `integrate(..., throw=False)` returns `Run.failed` per segment instead of raising, and `max_steps` bounds the work. XLA cannot be interrupted by SIGALRM, so the fit driver reads the flag.

### 2.7 Gradients

- **Forward mode** (`jacobian(net, observe)`: `jax.jacfwd` through `diffrax.ForwardMode()`) for least squares and for HMC with ~20 parameters. One JVP per coordinate.
- **The discrete adjoint** (`RecursiveCheckpointAdjoint`, diffrax's default) only once parameters outnumber observables, at about 50 or more. It is exact for the discretised solve (6e-15 from forward mode on the sulfur pot) but costs about twice the forward mode at 16 parameters.
- **Never the continuous adjoint** (`BacksolveAdjoint`): it did not return in ~23 min on the stiff sulfur pot.
- **The fits' forward-difference Jacobian is replaced.** On the sulfur pot it was up to 7 % off its own central difference. On the lipid pots the exact Jacobian equals the engine's central differences (analytic, so no solver noise) to <= 5e-8 per column.

### 2.8 Batches: many draws and pots in lockstep (the GPU path)

```python
f = network.batched(net, theta_axis=0, pot_axis=None, **solve_kw)   # (B draws) x one pot   -> (B, S, n)
f = network.batched(net, theta_axis=None, pot_axis=0, **solve_kw)   # one theta x (B pots)  -> (B, S, n)
g = network.grid(net, **solve_kw)                                   # (D draws) x (P pots)  -> (P, D, S, n)
thetas = network.stack_theta([theta_1, ..., theta_B])               # a batch of draws as one pytree
```

- **One compile per (network, S, B).** `jax.vmap` of `integrate`: each member keeps its own adaptive steps, and the batched loop runs until the slowest member is done.
- **Not stacked.** The batch is NOT one ODE over a (B, n) state. An implicit solver would then factor one dense (B n) x (B n) Jacobian instead of B small n x n ones, and one stiff member would set every member's step.
- **The benchmark harness** is `stage_batch` of [`build_jax_lipid_module.py`](../../scripts/generators/build_jax_lipid_module.py). It measures cost per member against B (1, 16, 128, 1024), compile time, a pots x draws grid, and that a batch equals its single solves. A GPU run (H2: benchmarked before any Colab run) reuses the same functions with a CUDA jaxlib.
- **What the GPU benchmark must report:**
  - per-member cost against B on CPU and GPU;
  - compile time;
  - float64 throughput. The contract needs float64: rates to 1e-12, and pH roots to 1e-9. Consumer and inference GPUs (T4, L4) run float64 at 1/32 - 1/64 of float32, so my expectation is no gain there (p ~0.8). A100 / H100 run float64 at 1/2 of float32 and are the only candidates.
  - whether the vmapped while loop's idle members (the slowest-member effect) erase the gain on stiff, pH-coupled batches.

On CPU the cost per draw does not fall with B: it is flat from 1 to 1024 draws (~10 - 12 ms on the 160 C lipid pot, Kvaerno5, rtol 1e-8, one core). vmap removes dispatch overhead only; the bake-off saw ~15 % at 8 draws on the sulfur pot. CPU parallelism therefore comes from processes over pots, and a batch's gain is a GPU question.

### 2.9 `Config` and refusals

- **Built from the `Config` it is given.** `build(cfg)` reads every value from a `Config` (`src/kinetic_core/config.py`; default `Config.from_globals()`) and records `cfg.digest()`. It never installs one, so modules built in one process from different configurations coexist. The lipid tests build from `Config.shipped().replace(...)` without touching any global.
- **Two kinds of key.** Build-phase and structural keys become static arrays. Numeric fields that fits or draws sample become `theta` coordinates, with the key's field `bands`.
- **A refusal on any unported switch.** Every switch a module reads is either implemented for every value it can take, or refused with `UnsupportedConfiguration`. Examples:
  - the lipid module implements `lipid.decadienal_breakage.enabled` both ways, and refuses an enabled `lipid.covalent_sink` (the engine raises there too);
  - the bake-off refuses every sulfur switch outside R3En.
  
  A module never computes a different model in silence.

### 2.10 Observables

`Module.observables` maps a name to a function of one state; the lipid module's are `ug_per_L:<species>`. The fit driver takes log10 where the rows are levels. What the cases add (sec. 4):
- phase-tagged observables (oil, water, whole sample: M5);
- observation operators with their own coordinates (an epimer share not integrated: M1c P2);
- the cooled-scale pH as an observable (M1, M1c);
- ratios of isotopologue pools (M4, M6).

## 3. Acceptance against today's lanes (replaces A2a's bit-identity)

A lane port is accepted when, against the engine's own numbers written by a core-env reference generator:
1. **every rate constant** (and, for pH lanes, the charge balance's pH and I) equals the engine's to 1e-12 relative on a grid of conditions and at a perturbed parameter point, before any trajectory is compared;
2. **trajectories** equal the engine's to 1e-6 dex at rtol 1e-10. The comparison runs on the lane's panel pots, their declared-band corners and one stiff or profile corner, for every species above a stated floor, with an absolute bound below it;
3. **the exact Jacobian** equals central differences (the engine's where they are noise-free, the port's own otherwise);
4. **at fit tolerance**, the panel residuals equal the engine's to 1e-4 dex.

The JAX tests skip in the core env (`pytest.importorskip("jax")`). Bit-identity cannot hold across LSODA and ESDIRK.

| lane | status | rates | trajectories (tight) | Jacobian | record |
|---|---|---|---|---|---|
| sulfur | bake-off, R3En configuration only; every other switch refused | 7.1e-15 (252 points x 132 reactions) | 1.8e-11 dex fitted rows (Kvaerno5), stiff corners <= 5e-6 | <= 4e-6 vs own central FD | [jax_bakeoff_sulfur_lane.md](../../results/validation/jax_bakeoff_sulfur_lane.md) |
| lipid | **ported as a contract module** (2026-11-05) | 1.3e-15 (4752 values: (T, Q10, donor) grid x 22 reactions) | 1.1e-9 dex all species, 4.5e-10 on the scored products (8 panel pots x point and 2 corners, Kvaerno5); profile 5.1e-8 | <= 5e-8 vs the engine's central FD (no solver noise), 4 pots | [jax_lipid_module.md](../../results/validation/jax_lipid_module.md) |
| trunk | pending | | | | |
| acrylamide | pending | | | | |

## 4. What the cases need from the network

From each case's "needs from the network" ([cases](../reference/cases/): M4 sec. 6, M5 sec. 6, M6 sec. 6, M3's network part N, M1's and M1b's step tables, and M1c (draft on branch `m1c-case`, sec. 1). "Hook" is where the contract carries the need; "status" is what exists in `network.py` today.

| need | cases | hook | status |
|---|---|---|---|
| one species shared by two modules (1-DG by hexose and sulfur; the Cys-hexose Amadori compound; 1-DPO, DDP, NF) | M4 | one owner per key; the other module `requires` it and touches it by name | done (`requires`, keyed stoichiometry) |
| one pot with every lane (Glc + Cys + Gly + protein) | M4 | `compose` | done for pH-free modules; pH modules wait on 2.4 |
| amine classes as species (free amino acid with identity factors, cysteine, peptide, protein-bound lysine with one key), each declaring which sugars it competes for | M4, M6 | `Module.competes`; `compose` refuses an undeclared shared consumption | design only |
| one Strecker step with two products per amino-acid class (aldehyde and aminoketone from one event) | M3, M6 | one reaction key with both products, per class; never two modules fitting it | design only |
| in-situ pH on the hexose Amadori / enolisation steps; one titratable table with the Cys-hexose Amadori pKa; trunk acids as titratable pools | M4, M1 (D), M1b | `Module.titratables` merged by `compose`, with one charge-closure check (closes v1's Gly / AMA / SB / PAMA gap) | design only (bake-off solver exists) |
| a measured pH trajectory as input, on a declared scale (cooled or in situ) | M1c, M1 | pH mode `prescribed` (a forced input) | design only |
| buffer as a composition, a per-step catalytic coefficient per (step, buffer species) as banded `ConfigKey`s, a hyperbolic form for lumped steps; constants re-referenced to buffer-free values from each calibration pot's buffer | M5, M6, M4 | pot data carry the buffer species; `RateCtx.buffer`; coefficients are `theta` coordinates with bands | design only |
| two phases: an oil volume fraction, fast partition per volatile (declared K_oil/water), reactions confined to their phase, phase-tagged observables, the headspace operator using the same partition once | M5 | species carry a phase tag; fast partition solved like the pH root (algebraic, IFT derivative); observables declare their phase | design only |
| the lipid carbonyls reacting with MFT, Cys and H2S (retiring the direct sum and `lipid.lane_coupling_verdict`) | M5 | a coupling module that `requires` lipid and sulfur; declared competition replaces the guard | lipid now composes as a module (this port); coupling design only |
| water: moisture and a_w as conditions, a whole-phase concentration basis, a Tg gate (warn within 20 K, refuse below), an a_w hook per step | M5, M6 | pot data `aw`, `RateCtx.aw`; the gate is a refusal at pot construction | design only |
| aminoketone, pyrrole and dicarbonyl-amino-acid adduct pools (free and protein-bound), a dihydropyrazine for R3 | M6 | species and reactions of the pyrazine module | design only |
| isotopologue pools (13C sugar, 15N amino acid) so labelling rows are observables | M4, M6 | a module generator that copies species per label and reactions per label combination from an atom map; ratio observables | design only |
| the forced-input driver: measured courses clamped (1-DG / 3-DG, GO / MGO, the ARP-Phe course), with a forced-exchange ledger so element closure still reports | M3, M4, M6 | forced species leave the state; a piecewise-linear interpolant in the right-hand side with `jump_ts` at its breakpoints | design only |
| a fast ring-chain exchange (TTCA_oc with TTCA, f_oc(T) in a band) | M1c | lumped as an algebraic share inside one pool (no stiff pair) | design only |
| observation operators with a band (the minor C-2 epimer share not integrated by an assay) | M1c | observable-level coordinates in `theta` | design only |
| headspace O2 from each pot's vessel volumes (OXR / OXV), a reversible disulfide | M1 | pot data carry the vessel; the sulfur module's switches | with the sulfur port |
| thermal profiles with ramps (B64) | all | profile mode (sec. 2.5) | design only |

## 5. Rules for standalone modules written before the network exists (investigation agents)

1. **Write the module as data plus two functions,** in a standalone script under `scripts/generators/`:
   - `SPECIES`, `REACTIONS` (the existing `Species` / `Reaction` types; keys from the engine's species where they exist, and new keys never colliding);
   - `CONFIG_KEYS` (with bands and basis);
   - `rate_constants(params, ctx)`;
   - `OBSERVABLES`.

   No module-level switches; everything structural is a `ConfigKey`. Rate constants depend on `params` and `ctx` only, so the JAX port is a transliteration.
2. **Forced inputs** are a dict `{species: (t_min, value)}` read from keyed fit rows, interpolated piecewise-linearly.
3. **Carbon / nitrogen / sulfur balance** of every reaction is checked in the module's unit test, with the same conventions (lumped FRAG_C / FRAG_N / FRAG_S, CARRIED carboxyls). Yields on free coordinates are yield channels (sec. 2.2).
4. **Each module's feasibility check carries:**
   - a positive control: the module's own synthetic truth at realistic noise passes;
   - a negative control: a declared known-wrong structure fails;
   - a declared "inconclusive" outcome.

A module that follows 1 - 3 ports onto the network as a `build(cfg) -> Module` in `src/kinetic_core_jax/`, without rewriting its chemistry, and is accepted under sec. 3. It still needs re-validation once the hub replaces its forced inputs.

## 6. Migration order and acceptance

1. **A3 (config): done** (phases 1 - 2, 2026-11-01). [`src/kinetic_core/config.py`](../../src/kinetic_core/config.py):
   - 42 `ConfigKey`s, one per switch, each listing every module binding of it;
   - `Config.shipped()` from declared defaults; `from_globals()`; `install` / `applied()`, which restore on exit;
   - `parallel.map_ordered(..., config=)`, `parallel.Configured`, `parallel.configured`, `MAILLARD_POOL=persistent`;
   - build-phase keys read-only (`BuildPhaseKeyError`); field-level `bands`; `Config.digest()`; `config.update(changes)`;
   - coverage: every UPPERCASE global an engine function reads at call time is a key binding or on the allowlist `tests/engine_constants.py`.
   
   Tests: `tests/unit/test_engine_config.py`. Not migrated: `engine._B2_FIT_REPORT`, scripts' own module globals, the tests' `monkeypatch.setattr` of switches. The JAX modules do not need any of these: they take a `Config` as a value (sec. 2.9).
2. **Lane ports (replaces A2a), each accepted under sec. 3, one writer per lane file.**
   1. lipid: done 2026-11-05 (pH-free; no running chemistry step touches it);
   2. sulfur: the bake-off becomes a contract module, every switch ported or refused, the pH solver moved into `network.integrate` (sec. 2.4), and profile mode added;
   3. trunk;
   4. acrylamide.
   
   A lane whose chemistry a running step still reshapes (M1 / M1b / M2 on sulfur and trunk) is ported after that step's verdict, so the port does not chase it.
3. **A2b (chemistry, judged by evidence, not identity).**
   - one titratable table, which closes the Gly / AMA / SB / PAMA gap;
   - one lysine key;
   - amine classes;
   - declared Cys competition;
   - the union network at a stated tolerance, with the recorded diff.
4. **The fit driver:**
   - least squares on forward-mode Jacobians;
   - per-pot failure flags;
   - one compile cache per configuration;
   - the rtol 1e-6 residual check (sec. 2.6).

**The hardest risks (v2).**
- **Forward-solve cost.**
  - The sulfur port's single solve is ~2.5x slower than LSODA, ~75 % of it the pH root.
  - The lipid port's is 7 - 170x slower than the engine's closed form (its Jacobian 3 - 50x slower than seven engine runs): an ODE solve cannot beat an analytic one. It wins only where Jacobians and batches are needed.
  - The exact Jacobian, not the forward solve, is what makes refits cheaper: 2.5 - 6.6x on the sulfur pot.
- **pH coupling.**
  - One composed balance makes the trunk's acids and amines titrate in pots that never did.
  - Two pH conventions must merge: the label vs the in-situ scale, and declared-pH scaling vs in-situ factors.
  - Each new behaviour gets an explicit config key, off by default.
- **Compile time.**
  - ~12 - 30 s per configuration for the sulfur pot; ~2 s per solve and ~3 s per Jacobian for the lipid module.
  - Pots as arguments keep it to once per configuration and segment count.
  - Every switch a fit varies is a recompile, so a fit's switches are fixed per compile.

## 7. Effort

The bake-off estimated 75 - 115 h for A2 + H3. The lipid port (contract runtime, module, reference generator, builder and tests) took about 1 h of agent work. It is the easiest lane: pH-free, linear, with two switches. So it measures the per-lane fixed cost, about 0.5 h, not the chemistry cost.

Revised estimates, in hours of agent work (finding 6 of [jax_lipid_module.md](../../results/validation/jax_lipid_module.md)):

| work | hours |
|---|---|
| sulfur (every remaining switch, the pH solver in `integrate`, profile mode) | 12 - 20 |
| trunk | 8 - 14 |
| acrylamide | 5 - 9 |
| network infrastructure (merged titratable table, declared competition, forced inputs, the fit driver) | 15 - 25 |
| the pH-root performance step | 8 - 15 |
| **total** | **~50 - 85** (p ~0.6 of landing inside) |
