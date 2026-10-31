# A2 module contract (draft v1, 2026-10-31)

This is the interface through which chemistry modules contribute to ONE composable network, with one frozen configuration and one integrator. It is the target for A3 (configuration) and A2 (composition) in [STATUS_AND_PLAN.md](STATUS_AND_PLAN.md) sec. 5, and for every standalone module that investigation agents write before A2 exists. It is based on a read-only survey of `src/kinetic_core/` (2026-10-31); file:line references are to that date.

## 1. Why (what the survey found)

- **Three trunks.**
  - `network.TRUNK_REACTIONS` (`network.py:520`);
  - sulfur `REACTIONS + FED_3DEOXY + SULFUR_REACTIONS` (`sulfur.py:1009`);
  - acrylamide `REACTIONS + ACRYLAMIDE_REACTIONS` (`acrylamide.py:268`).
  
  The name `TRUNK_REACTIONS` means a smaller set in sulfur / acrylamide than in `network.py`. Three copies of the balance checker, of the stoichiometry builder and of the `r_amadori` derivation exist. Two functions are named `ph_factor`, with different signatures (`trunk_conditions.py:156`, `parameters_sulfur.py:992`).
- **About 30 module-level switches** (sulfur pH coupling, `ph_state.B51`, `TRUNK_SOLUTES_IN_BALANCE`, 13 sulfur dataclass switches, lipid, dicarbonyl, panel, uncertainty). They are set by about 97 assignments in 39 scripts and about 50 in tests; none is set in `src/`. One joint-fit job sets about 15 of them across 3 files.
- **Parallel correctness depends on forking** a new pool per call (`parallel.py`). A spawn-started or persistent worker would silently compute with the DEFAULT switches.
- **Charge closure has a gap.** With `TRUNK_SOLUTES_IN_BALANCE=True`, Gly / AMA / SB / PAMA titrate in `titratable_inventory` (`ph_state.py:653`) but are absent from `TITRATABLE_CENTRES` (`:1069`). The closure check (`validate_charge_closure`, `:1305`) therefore never sees trunk steps that move them. The ledger (validation) and the inventory (physics) are two tables for one job.
- **Amine identity is a convention.**
  - `Gly` stands for "the amine": MET / PRO are folded into Gly on the trunk and popped elsewhere (`engine.py:~2553-2570`).
  - Lysine has two keys (`Lys` in acrylamide, `LYSP` in the trunk).
  - Cys is consumed by two lanes, which `resolve_lanes` therefore refuses to compose.
- **Only the sulfur lane has dynamic pH.** The trunk scales k_ref at the declared pH, acrylamide has no pH argument, and lipid is pH-free.

## 2. The contract

```python
@dataclass(frozen=True)
class ConfigKey:
    name: str                      # namespaced: "sulfur.tdp_cys_trap.k", "ph.coupling", "ph.b51", ...
    default: Any                   # the SHIPPED value
    band: tuple | None             # declared / measured band (feasibility draws sample inside it)
    phase: Literal["build", "rate", "solver"]
    basis: str                     # dossier / prereg / choice-register reference

@dataclass(frozen=True)
class Config:                      # replaces every module-level switch; hashable; travels with a job
    values: Mapping[str, Any]
    def replace(self, **kw) -> "Config": ...
    @classmethod
    def shipped(cls) -> "Config": ...
    @classmethod
    def from_globals(cls) -> "Config": ...    # A3 migration bridge only; deleted at the end of A2

@dataclass(frozen=True)
class RateCtx:
    T_K: float; ph: float | None; ph_label: float | None; I_M: float
    cfg: Config; k: Mapping[str, float]       # rates already computed (for derived rates)

class Module(Protocol):
    name: str
    requires: tuple[str, ...]                              # e.g. sulfur requires ("trunk",)
    species: tuple[Species, ...]
    reactions: tuple[Reaction, ...]                        # the existing frozen Reaction type
    config_keys: tuple[ConfigKey, ...]
    titratables: Mapping[str, TitratableSpec]              # ONE table: pKa physics + carboxyl / amine centres
    centre_ledger: Mapping[str, Mapping[str, Any]]         # own reactions only; validated against titratables
    observables: Mapping[str, Callable[[Mapping[str, float]], float]]
    def parameters(self, frozen: Mapping, cfg: Config) -> Mapping[str, Any]: ...
    def rate_constants(self, params, ctx: RateCtx) -> Mapping[str, float]: ...   # own reactions only
    def off(self, cfg: Config) -> frozenset[str]: ...

def compose(modules: Sequence[Module], *, layout: Sequence[str] | None = None) -> Network: ...
def integrate(net: Network, params, cfg: Config, initial, program: ThermalProgram, *,
              ph: PhSpec | None = None,
              forced: Mapping[str, tuple[np.ndarray, np.ndarray]] | None = None,
              t_eval=None) -> Run: ...
```

**What `compose` does:**
- orders the species (a `layout` pins the legacy order);
- builds the stoichiometry matrix and reactant index / power arrays;
- merges the ledgers and runs ONE carbon / nitrogen / sulfur balance check plus ONE charge-closure check over the merged titratable table;
- precompiles ln k_ref, Ea / R, the pH-kind index and the caps into arrays, so the only Python per right-hand-side call is the derived-rate closures.

**Duplicate keys refuse to compose.** A species or reaction key defined by two modules, or a species consumed by two modules that do not declare the competition (as with Cys today), is an error. `compose` never merges silently.

**The forced-input driver.** Species named in `forced` leave the ODE state. Their value is the piecewise-linear interpolant of a measured series, substituted into the rates and the pH inventory. The net flux into or out of each clamped species is integrated as a "forced exchange" ledger, so element closure still reports a balance. A Strecker or pyrazine module then runs on measured dicarbonyls, and later on the hub, with no code change.

## 3. Rules for standalone modules written before A2 exists (investigation agents, wave 1+)

1. **Write the module as data plus two functions,** in a standalone script under `scripts/generators/`:
   - `SPECIES`, `REACTIONS` (the existing `Species` / `Reaction` types; keys from the engine's species where they exist, and new keys never colliding);
   - `CONFIG_KEYS` (with bands and basis);
   - `rate_constants(params, ctx)`;
   - `OBSERVABLES`.
   
   No module-level switches; everything structural is a `ConfigKey`.
2. **Forced inputs** are a dict `{species: (t_min, value)}` read from keyed fit rows, interpolated piecewise-linearly.
3. **Carbon / nitrogen / sulfur balance** of every reaction is checked in the module's unit test, with the same conventions (lumped FRAG_C / FRAG_N / FRAG_S, CARRIED carboxyls).
4. **Each module's feasibility check carries:**
   - a positive control: the module's own synthetic truth at realistic noise passes;
   - a negative control: a declared known-wrong structure fails;
   - a declared "inconclusive" outcome.

A module that follows 1-3 ports onto A2 by wrapping it in the `Module` protocol, without rewriting its chemistry. It still needs re-validation once the hub replaces its forced inputs.

## 4. Migration order and acceptance

1. **A3 (config, about 11-13 h).**
   - `Config` and `ConfigKey` are introduced, with `from_globals()` as the bridge.
   - Each job carries a `Config`, which a worker installs at the start of each task: the persistent pool becomes safe.
   - Scripts move from setting globals to building `Config.replace(...)`.
   
   Acceptance:
   - all gates and the full suite pass;
   - the core panel and every signature record are bit-identical;
   - a test that a spawn-started worker reproduces a forked one.
2. **A2a (bit-identical, time-boxed to 3 days).**
   - The trunk, sulfur, acrylamide and lipid lanes become modules.
   - The three legacy lanes become named compositions with their legacy layouts and operation order (derived rates in their current order; the B27 block last).
   
   Acceptance: equality-pinned panel tests per lane.
3. **A2b (chemistry, judged by evidence, not identity).**
   - one titratable table, which closes the Gly / AMA / SB / PAMA gap;
   - one lysine key;
   - amine classes;
   - declared Cys competition;
   - the union network at a stated tolerance, with the recorded diff.

**The three hardest risks.**
- **Bit-identity.** The lanes use different arithmetic paths (dict-loop rates vs vectorised; k_ref pre-scaling vs in-rate pH factors; order-dependent derived rates) and different species orders and shapes. LSODA's step sequence depends on both.
- **pH coupling.** One composed balance makes the trunk's acids and amines titrate in pots that never did, and two pH conventions must merge: label vs in-situ scale, and declared-pH scaling vs in-situ factors. Each new behaviour gets an explicit config key, off by default, and the fixed point's convergence is re-verified.
- **Performance.** `sulfur_rate_constants_at` is a Python loop over about 150 reactions that the inline mode calls on every right-hand-side and Jacobian evaluation. `compose` must precompile, and `Config` is read once per integration.
