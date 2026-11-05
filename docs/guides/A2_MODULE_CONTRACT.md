# A2 module contract (v2.2, 2026-11-05; v2 2026-11-05; draft v1 2026-10-31)

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
13. **(v2.2) The network infrastructure** ([jax_network_infra.md](../../results/validation/jax_network_infra.md)): forced inputs (sec. 2.11), the fit driver (sec. 2.12), the bridge to M1.1's predictive densities (sec. 2.13), and the declared-competition and titratable-table hooks of `compose` (sec. 2.3; the charge closure is a typed stub until the sulfur port).

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
    pot_inputs: tuple[str, ...] = ()                         # (v2.1) pot scalars read via RateCtx: "ph_label", "aw"

class RateCtx(NamedTuple):          # what a rate function may read besides theta
    t_c; t_k                        # T_K = T_C + 273.15, as the engine computes it
    ph = None; ph_label = None; ionic_m = None               # pH modules only (sec. 2.4)
    aw = None                                                # (v2.1) the pot's declared a_w; NaN = the engine's None

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

**(v2.2) Declared competition and the titratable table.**
- `Module.competes` names the species OWNED BY ANOTHER MODULE that the module's reactions consume (a mass-action reactant with a negative net coefficient). `compose` records every shared consumption on `Network.competition` (species -> owner, consumers) and the undeclared pairs on `Network.undeclared_competition`. `compose(..., competition="refuse")` refuses an undeclared pair, and a `competes` entry no module owns is always refused. The default is `"report"`, because the ported trunk-lane and acrylamide modules predate the declaration: the trunk lane's 55 steps carry undeclared pairs (Glc, Gly, MGO, AMA, ODG, ... consumed by `trunk_lane` from `trunk`). The sulfur port, which makes Cys a species two lanes consume (sec. 1), is where `"refuse"` becomes the default and every module declares its pairs.
- `Module.titratables` (`Titratable(species, pka_25c, charges, basis)`) are merged by `compose` into `Network.titratables`, refusing a pool declared twice, a pool the module neither owns nor requires, and a charge list that does not have one state more than the pKa list. No ported module declares one yet.
- **The ONE charge-closure check is a typed stub** (`network.charge_closure` raises `UnsupportedConfiguration`). It needs, per reaction, the charge each titratable reactant and product carries at the network's pH: the inline charge balance and the titratable inventory that only the sulfur lane carries (`ph_state.titratable_inventory` / `TITRATABLE_CENTRES`, the bake-off's `model.make_ph_solver`). It lands with the sulfur port.

**(v2.2) The forced-input driver** is sec. 2.11.

### 2.4 The pH coupling (from the bake-off)

- **The root is solved inside the right-hand side.** The bake-off solves `F(u; c) = 0`, with u = (pH, I) and c the titratable inventory, at every right-hand-side call: a bisection (the residual is monotone in pH), then Newton on the coupled (pH, I). It solves twice, at the reaction temperature and at the cooled 25 C label. Its derivative is a custom JVP: du = -(dF/du)^-1 (dF/dc) dc, a 2 x 2 solve, and reverse mode transposes it with no extra code ([`model.make_ph_solver`](../../src/kinetic_core_jax/model.py)). There is no DAE state: diffrax has no DAE solver, and a relaxation would change the model.
- **Modes:** `none` (pH-free modules; the lipid module), `inline` (the above), `prescribed` (a measured pH trajectory read as a forced input; M1c: Cao 2017, Kang 2026), and `clamped` (pH-stat).
- **Status.** A network with any `needs_ph` module is refused by `integrate` until the sulfur module is ported, and that port reuses the bake-off's solver.
- **(v2.1) `prescribed` runs in its constant form** (the trunk and acrylamide ports, 2026-11-05): a module that lists `"ph_label"` in `Module.pot_inputs` reads the pot's DECLARED pH (`pot["ph_label"]`, NaN for None) through `RateCtx.ph_label`, constant over the programme. This is exactly the engine's declared-pH scaling of k_ref (trunk B12 / B40, acrylamide B15), which has no pH trajectory; `Network.ph_mode` reports it. (v2.2) The trajectory form (M1c's measured pH courses) is a forced `ph_label` course (sec. 2.11; `forcing.ph_mode` reports "prescribed (trajectory)"), on a declared scale: "cooled" or "label_25C", both read as the label-scale declared pH; "in_situ" is refused until the sulfur port's inline mode. `"aw"` is read the same way (`RateCtx.aw`), the first piece of sec. 4's water hook. `integrate` refuses a pot that lacks a scalar its modules read.
- **Cost.** The root-find was ~75 % of the bake-off's forward solve (2.1 vs 0.5 ms per step). It gets one focused step: a cheaper residual (a precomputed pK table, fewer exps), and one root per step reused across stages.
- **Measured pH tables with jumps** (Zheng & Ho's mechanism boundary at label pH 5.0, the thiamine share's floor at 6.25, the lipid breakage window at 100 C in temperature) are ported exactly with `jnp.where`. An exact derivative omits the jump-crossing term, which a finite difference across the boundary would see. Rule: a module whose table has a jump declares it, and the fit driver reports every pot whose trajectory crosses one, so a Jacobian there is read as one-sided. (v2.2) The declaration is `Module.breakpoints` (input -> values: the lipid window at 100 C, B18's knot at pH 7, the trunk and acrylamide a_w knots and window edges, acrylamide's pH window edges). For a FORCED course the crossing time is a function of the data alone, so `forcing.forced_pot` puts a knot there and the theta-Jacobian is exact (no missing term); the one-sided caveat stays for an inline pH, whose crossing time depends on theta. No compared column crossed a boundary in the bake-off or in the lipid port: the lipid window switches between segments, never inside one.

### 2.5 The thermal program

- **Segments.** `program_segments(segments, n_segments)` turns the engine's `ThermalProgram.segments` (`((duration_min, T_C), ...)`, B64's staircase for profiles) into an (S, 2) argument, padded with zero-duration segments.
- **Save times.** `refine_segments(segments, save_times)` splits segments so that every save time is a boundary. `integrate` returns the state at the end of every segment.
- **Rate constants** are evaluated once per segment when no module needs pH. A pH module evaluates them in the right-hand side.
- **Not yet ported:** the sulfur lane's PROFILE mode, which integrates T(t) piecewise-linear with its rate constants at T(t). (v2.2) For a pH-free network this is a forced linear `t_c` course (sec. 2.11): each knot is a segment boundary, which is what `jump_ts` would do. The sulfur lane's own profile mode (pH inline) comes with its port.

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

### 2.11 (v2.2) Forced inputs ([`forcing.py`](../../src/kinetic_core_jax/forcing.py))

```python
F = forcing.Forcing(inputs=("t_c", "ph_label"), species=("Gly",),
                    interpolation=(("t_c", "previous"),), ph_scale="cooled")       # static: one compile per F
fp = forcing.forced_pot(net, F, y0, {"t_c": program_course(segments), "ph_label": (t_min, ph), "Gly": (t_min, mM)},
                        scalars={"aw": 0.9}, save_times=[...], n_segments=None)  # data: knots and values
run = forcing.integrate_forced(net, theta, fp.pot, F)    # ys (S, n), exchange (S, n_forced), num_steps, failed
```

- **What can be forced:** a `RateCtx` scalar (`t_c`, `ph_label`, `aw`) or a species. A forcing always carries `t_c`; a thermal program is a 'previous' course (`program_course`).
- **Interpolation** per name: "linear" or "previous" (a staircase). Each interval between knots is one segment of `integrate`'s design (unit interval, duration-scaled right-hand side), with the inputs a + s (b - a) inside it. A knot is therefore an integrator restart. The rate constants are evaluated in the right-hand side when a forced scalar is linear, and once per interval otherwise. Courses with the same number of intervals share one compile; zero-duration intervals pad.
- **A forced species leaves the state.** The interpolant enters every mass-action term, its derivative is zero, and what the network produces into it (or draws from it) is integrated into a forced-exchange ledger (`ForcedRun.exchange`). `forcing.closure` checks E_free . dy_free + E_forced . ledger = 0.
- **Derivatives** are exact through the solve, in theta and in the knot values (both are arguments). `jacobian_forced` and `batched_forced` mirror `jacobian` and `batched`; the fit driver and the density bridge take a `forcing=`.
- **Refused:** an inline-pH network (forcing pH there is "clamped", with the sulfur port); forcing an input no module reads; an in-situ pH scale; a forced species the network does not carry; a non-finite course, non-increasing knots, a negative species level; an acrylamide network whose forced temperature peaks inside CF-106's window (the engine's refusal, kept).
- **Checks** ([jax_network_infra.md](../../results/validation/jax_network_infra.md) sec. 1):
  - a forced constant equals pH mode "prescribed";
  - a 'previous' pH + temperature course equals the engine's own staircase run (the engine's per-segment path with a per-segment declared pH, which is the only pH trajectory the engine can run);
  - linear courses converge on midpoint staircases at second order;
  - a forced constant species equals the network with its row removed, and the ledger closes the elements;
  - the Jacobian through a linear pH course that crosses B18's knot equals central differences.

### 2.12 (v2.2) The fit driver ([`fit.py`](../../src/kinetic_core_jax/fit.py))

- **A `JaxProblem`** has:
  - the free coordinates, a subset of `Network.coordinates`, the others fixed at `theta` (default shipped);
  - their kinds: the engine's (`log10k`, `Ea`, `acid_yield`, `pka`, `ph_slope`, `log10_offset`) plus `linear` (a zero-shipped k_ref, bounded below by 0), `share` and `q10`;
  - their bounds;
  - many pots, plain or forced;
  - `Observation`s: a pot, an output row, a Network observable, a value, a sigma, log10 by default.
- **`fit`** runs scipy's `least_squares` with the engine's `LEAST_SQUARES_SETTINGS` on the JAX residual and its exact forward-mode Jacobian, within a budget (max_nfev), with the engine's `require_converged` semantics (`FitNotConverged`; a pot that failed at the optimum also fails convergence).
- **No central polish:** M1.1's polish corrects the forward-difference Jacobian, and this one is exact.
- **The solve is compiled once per pot shape.** It returns d state / d x by forward mode; the observables are applied after it with their exact gradients, so one compile serves every pot of a shape, whatever is observed on it, and every data set and start.
- **Batching:** `batch="vmap"` runs a group of same-shape pots in lockstep; `batch="loop"` runs one pot at a time on the same compile.
- **Failed pots** read `failed_residual` (25) with a zero Jacobian row, as the engine's `Problem`.
- **Diagnostics:** `laplace` (exact-Jacobian sigmas, eigenbasis, the engine's null-weight rule), `fd_check` (central differences with column-scaled steps, the 2026-11-05 lesson) and `residual_tolerance_check` (sec. 2.6's rtol 1e-6 rule).

### 2.13 (v2.2) The bridge to M1.1's predictive densities ([`densities.py`](../../src/kinetic_core_jax/densities.py))

- **`Predictor(net, pots, observations, keys)`** maps draws (D, p) of a coordinate vector to the unit's predictions (D, n). Every draw of a pot is solved in one compiled batch (`vmap` over draws, chunked, the last chunk padded). Failed draws come back as NaN with a mask.
- **`log_likelihoods(pred, y, cov)`** is `gaussian_log_density` for every draw, with one Cholesky. It refuses failed draws unless the caller declares `on_failed="impossible"` (-inf, which `importance_lppd` counts).
- **`log_integrand`** composes both with a log prior, as `quadrature_lppd`'s callback and `importance_lppd`'s array.
- **Scope:** units whose pots the network covers, which today means the trunk, trunk-lane, acrylamide and lipid lanes. Not sulfur.
- **Cost on CPU:** no gain at equal tolerance. Paired on a trunk-lane pot, the engine ran 9.0 draws/s against 7.5 at rtol 1e-8 and 21 at rtol 1e-6; draws/s are flat from 64 to 2048 draws. The bridge's value is exactness and one code path. Its speed is a GPU question (sec. 2.8).

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
| sulfur (`sulfur.FULL_REACTIONS`, 132 steps, 62 species) | **ported as a contract module** (2026-11-05) at the shipped configuration (pH coupling "fixed_point_nodes"; "inline" ported too); every rival refused | see the record (rates on a (T, pH, label pH, I) grid, shipped and two envelope draws) | every panel pot the lane answers, live against the engine at its shipped rtol 1e-8 (bound 1e-4 dex, sweep counts equal), three envelope draws, profile / staircase / clamped / inline pots; three pots at rtol 1e-10 | `jax.grad` (reverse mode through the fixed point) vs own central FD, 5 coordinates | [jax_sulfur_lane_module.md](../../results/validation/jax_sulfur_lane_module.md) |
| lipid | **ported as a contract module** (2026-11-05) | 1.3e-15 (4752 values: (T, Q10, donor) grid x 22 reactions) | 1.1e-9 dex all species, 4.5e-10 on the scored products (8 panel pots x point and 2 corners, Kvaerno5); profile 5.1e-8 | <= 5e-8 vs the engine's central FD (no solver noise), 4 pots | [jax_lipid_module.md](../../results/validation/jax_lipid_module.md) |
| trunk (`network.REACTIONS`, the 26 steps the acrylamide lane runs) | **ported as a contract module** (2026-11-05), with the trunk LANE's declared terms on these steps | 1.3e-15 (75 816 values: trunk-lane terms, (T, pH, a_w) grid x 26 reactions, shipped / perturbed / exit term off) | accepted inside the acrylamide lane (next row) and inside the trunk lane (row after) | (next rows) | [jax_trunk_acrylamide_module.md](../../results/validation/jax_trunk_acrylamide_module.md) |
| trunk LANE (`network.TRUNK_REACTIONS`, 55 steps: trunk + `trunk_lane` (six trunk-only blocks, 26 steps) + `fed_3deoxy` (3)) | **ported as contract modules** (2026-11-05), every switch ported or refused; pH mode prescribed (the declared pH) | 1.3e-15 (160 380 values: (T, pH, a_w) grid x 55 reactions; shipped / perturbed with B22 / B24 on / exit term off) | 1.2e-9 dex all species (Kvaerno5; Kvaerno4 3.5e-10): 1 panel pot x point and 2 band corners, 4 synthetic isothermal pots (one per idle block), staircase profile shipped and switches-on; panel rows at rtol 1e-6 <= 5.3e-7 dex | <= 3.3e-6 vs own FD (linear zero-k columns one-sided, step scaled to the column); <= 5.2e-5 rel / 3.2e-6 abs vs the engine's (LSODA-noise-limited), 3 pots x 109 coordinates | [jax_trunk_lane_module.md](../../results/validation/jax_trunk_lane_module.md) |
| acrylamide (trunk + `ACRYLAMIDE_REACTIONS`, 42 steps) | **ported as a contract module** (2026-11-05), every `acrylamide_conditions` term; CF-106 refused at the pot | 1.2e-15 (81 648 values: (T, pH, a_w) grid x 42 reactions, shipped and perturbed) | 2.2e-9 dex all species, 1.9e-11 on ACR / HMF (6 panel pots x point and 2 band corners, Kvaerno5); profile with all 42 steps 3.3e-9; panel rows at rtol 1e-6 <= 6e-8 dex | <= 1.3e-6 vs own central FD; <= 8e-5 rel / 3.3e-7 abs vs the engine's (LSODA-noise-limited), 3 pots x 84 coordinates | [jax_trunk_acrylamide_module.md](../../results/validation/jax_trunk_acrylamide_module.md) |

## 4. What the cases need from the network

From each case's "needs from the network" ([cases](../reference/cases/): M4 sec. 6, M5 sec. 6, M6 sec. 6, M3's network part N, M1's and M1b's step tables, and M1c (draft on branch `m1c-case`, sec. 1). "Hook" is where the contract carries the need; "status" is what exists in `network.py` today.

| need | cases | hook | status |
|---|---|---|---|
| one species shared by two modules (1-DG by hexose and sulfur; the Cys-hexose Amadori compound; 1-DPO, DDP, NF) | M4 | one owner per key; the other module `requires` it and touches it by name | done (`requires`, keyed stoichiometry) |
| one pot with every lane (Glc + Cys + Gly + protein) | M4 | `compose` | done for pH-free modules; pH modules wait on 2.4 |
| amine classes as species (free amino acid with identity factors, cysteine, peptide, protein-bound lysine with one key), each declaring which sugars it competes for | M4, M6 | `Module.competes`; `compose` refuses an undeclared shared consumption | hook done (v2.2: `Module.competes`, `compose(competition="report" / "refuse")`, sec. 2.3); the amine-class species are design only |
| one Strecker step with two products per amino-acid class (aldehyde and aminoketone from one event) | M3, M6 | one reaction key with both products, per class; never two modules fitting it | design only |
| in-situ pH on the hexose Amadori / enolisation steps; one titratable table with the Cys-hexose Amadori pKa; trunk acids as titratable pools | M4, M1 (D), M1b | `Module.titratables` merged by `compose`, with one charge-closure check (closes v1's Gly / AMA / SB / PAMA gap) | table merged by `compose` (v2.2, sec. 2.3); the closure check is a typed stub until the sulfur port (bake-off solver exists) |
| a measured pH trajectory as input, on a declared scale (cooled or in situ) | M1c, M1 | pH mode `prescribed` (a forced input) | done on the label scale (v2.2: a forced `ph_label` course, sec. 2.11); in situ waits on the sulfur port |
| buffer as a composition, a per-step catalytic coefficient per (step, buffer species) as banded `ConfigKey`s, a hyperbolic form for lumped steps; constants re-referenced to buffer-free values from each calibration pot's buffer | M5, M6, M4 | pot data carry the buffer species; `RateCtx.buffer`; coefficients are `theta` coordinates with bands | design only |
| two phases: an oil volume fraction, fast partition per volatile (declared K_oil/water), reactions confined to their phase, phase-tagged observables, the headspace operator using the same partition once | M5 | species carry a phase tag; fast partition solved like the pH root (algebraic, IFT derivative); observables declare their phase | design only |
| the lipid carbonyls reacting with MFT, Cys and H2S (retiring the direct sum and `lipid.lane_coupling_verdict`) | M5 | a coupling module that `requires` lipid and sulfur; declared competition replaces the guard | lipid now composes as a module (this port); coupling design only |
| water: moisture and a_w as conditions, a whole-phase concentration basis, a Tg gate (warn within 20 K, refuse below), an a_w hook per step | M5, M6 | pot data `aw`, `RateCtx.aw`; the gate is a refusal at pot construction | design only |
| aminoketone, pyrrole and dicarbonyl-amino-acid adduct pools (free and protein-bound), a dihydropyrazine for R3 | M6 | species and reactions of the pyrazine module | design only |
| isotopologue pools (13C sugar, 15N amino acid) so labelling rows are observables | M4, M6 | a module generator that copies species per label and reactions per label combination from an atom map; ratio observables | design only |
| the forced-input driver: measured courses clamped (1-DG / 3-DG, GO / MGO, the ARP-Phe course), with a forced-exchange ledger so element closure still reports | M3, M4, M6 | forced species leave the state; a piecewise-linear interpolant in the right-hand side with `jump_ts` at its breakpoints | done (v2.2: `forcing.py`, sec. 2.11; every knot is a segment boundary, which does the `jump_ts`' job) |
| a fast ring-chain exchange (TTCA_oc with TTCA, f_oc(T) in a band) | M1c | lumped as an algebraic share inside one pool (no stiff pair) | design only |
| observation operators with a band (the minor C-2 epimer share not integrated by an assay) | M1c | observable-level coordinates in `theta` | design only |
| headspace O2 from each pot's vessel volumes (OXR / OXV), a reversible disulfide | M1 | pot data carry the vessel; the sulfur module's switches | with the sulfur port |
| thermal profiles with ramps (B64) | all | profile mode (sec. 2.5) | done for pH-free networks (v2.2: a forced linear `t_c` course, sec. 2.11); the sulfur profile mode comes with its port |

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
   2. sulfur: the bake-off becomes a contract module, every switch ported or refused, the pH solver moved into `network.integrate` (sec. 2.4), and profile mode added. (2026-11-05) Ported at the shipped configuration as ONE module ([`sulfur_lane.py`](../../src/kinetic_core_jax/sulfur_lane.py), the engine's own copy of the trunk included; structure derived from the engine, hand-ported terms listed in [`sulfur_structure.py`](../../src/kinetic_core_jax/sulfur_structure.py)), with both pH couplings ("fixed_point_nodes", shipped, and "inline"), the clamped pH, staircases and profiles; every off-by-default rival is refused and stays NumPy-only ([jax_sulfur_lane_module.md](../../results/validation/jax_sulfur_lane_module.md)). Still open: the pH solver lives in `sulfur_lane.integrate`, not yet in `network.integrate`, and the charge closure over the merged table. **Any later edit to `src/kinetic_core/sulfur.py` must update the JAX lane (`src/kinetic_core_jax/sulfur_lane.py` and its registry `sulfur_structure.py`) in the same PR: `tests/unit/test_jax_sulfur_structure.py` (core env, every PR) and `tests/scientific/test_jax_sulfur_parity.py` (the `maillard-jax` CI job "JAX sulfur lane parity", live against the engine) enforce it. This replaces the one-writer rule for the sulfur lane while M1b(i) and M1c still edit it;**
   3. trunk: `network.REACTIONS` done 2026-11-05 (the acrylamide lane's trunk); the trunk LANE done 2026-11-05 ([`trunk_lane.py`](../../src/kinetic_core_jax/trunk_lane.py): the six trunk-only blocks as one module, because B18's pyrazine pH slopes are shared by the pyrazine, methionine and proline steps and a rate function reads only its own theta; the fed 3-deoxy triangle as its own module, because the sulfur network carries it too). Steps the engine carries at exactly zero have a linear `k_ref:` coordinate;
   4. acrylamide: done 2026-11-05 (pH mode prescribed: the declared pH; Kvaerno4 takes ~0.6x Kvaerno5's steps on its pots).
   
   A lane whose chemistry a running step still reshapes (M1 / M1b / M2 on sulfur and trunk) is ported after that step's verdict, so the port does not chase it.
3. **A2b (chemistry, judged by evidence, not identity).**
   - one titratable table, which closes the Gly / AMA / SB / PAMA gap;
   - one lysine key;
   - amine classes;
   - declared Cys competition;
   - the union network at a stated tolerance, with the recorded diff.
4. **The fit driver: done 2026-11-05 (v2.2, sec. 2.12)**:
   - least squares on forward-mode Jacobians;
   - per-pot failure flags;
   - one compile per pot shape;
   - the rtol 1e-6 residual check (sec. 2.6).

   On two synthetic toys it reaches `fitting.fit`'s central-polished optimum: cost to 4e-8 relative, theta to 2.3e-4 Laplace sigma. Paired with the engine at the same load it is NOT faster at 5 coordinates: the JAX solve costs ~3x LSODA, the exact Jacobian ~1.6x the forward-difference one and ~0.9x the central one, and a warm fit runs ~1.2x the polished engine fit ([jax_network_infra.md](../../results/validation/jax_network_infra.md) findings 3 - 5). Forward mode pays at many coordinates and on stiff pots where the forward-difference Jacobian is poor (the bake-off's 2.5 - 6.6x at 16 coordinates).

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
| ~~trunk lane's seven blocks~~ (done 2026-11-05: ~1 h of agent work plus ~1.5 h of unattended compute, against 3 - 6 h) | 0 |
| ~~network infrastructure~~ (done 2026-11-05: forced inputs, the fit driver, the density bridge, the competition and titratable hooks; the charge closure is a typed stub that lands with sulfur; ~2.5 h of agent work plus ~0.5 h of unattended compute, against 6 - 12 h) | 0 |
| sulfur (every remaining switch, the pH solver in `integrate`, profile and clamped modes, the in-situ pH scale, the charge closure over the merged table, declared competition switched to "refuse") | 6 - 12 |
| the pH-root performance step | 3 - 6 |
| **total remaining** | **~9 - 18** (p ~0.6 of landing inside); unattended compute on a shared machine adds ~1x |
