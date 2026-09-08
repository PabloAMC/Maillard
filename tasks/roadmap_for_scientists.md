# Roadmap: the tool a plant-based flavour scientist would keep using

*Written 2026-09-08 after the hypothesis layer shipped. Five programmes, in the order of value to a
scientist formulating meaty flavour on pea or soy protein. Each has a design that attaches to a seam
the code already has, a pre-registration, tests, and a success test declared before the work. The
backlog (`data_restructure_plan.md` section 7) points here; this file is the plan.*

## 0. The rules the design keeps

1. **One engine.** Nothing predicts except `src/kinetic_core`. New chemistry enters as species,
   reactions and declared or fitted parameters of that engine, never as a second path.
2. **Every number has a source or a band.** A user's data is a source for the user's calibration
   only; it never enters `data/` and never moves the shipped parameters.
3. **Fit primary evidence, validate levels.** For a user's own pot this becomes: **levels fit the
   laboratory's response factor; contrasts fit the kinetics.** A response factor is a property of
   the measurement, not of the chemistry, so a level may set it. Rate constants move only on
   within-laboratory contrasts (time points, temperatures, pH, recipes).
4. **Overlay, do not fork.** The engine already reads its parameters from frozen fit reports and
   accepts an override mapping (`core_parameters(lane, frozen=...)`). A calibration, a matrix
   declaration and a wave are all overlays on that one mechanism. No new configure/restore dance.
5. **Refuse rather than guess, and say which.** Every new answer either has a rate with a band or is
   a named refusal; the hypothesis layer says whether the route exists.
6. **A pre-registration before a run, an artifact under the freshness gate after it, a test for
   every promise.** Six dependencies stay six.

## 1. Calibrate on the user's own data (`maillard calibrate`)

**Why first.** A laboratory's own GC-MS data on its own matrix is worth more to that laboratory than
any literature pot, and the parity plot says the misses are systematic offsets by path, which is
exactly what a per-laboratory factor absorbs. Today `score` writes a record and refits nothing.

**Design.**
- `src/kinetic_core/calibration.py`: a frozen `Calibration` dataclass: `base_wave` (the shipped
  report it overlays), `lab`, `response_factor_log10` per compound (with its sigma), `overrides`
  (the few kinetic coordinates the user's contrasts identified, each with value, sigma and band),
  `rows` (which user records fitted, which validated), `provenance`. Serialised as
  `results/user/<lab>/calibration_<date>.json` with a markdown card beside it.
- `src/kinetic_core/user_fit.py`: builds residuals from the bundle-shaped records `user_scoring`
  already writes. Two stages, both least squares on log10: (i) response factors from every level
  row, one per compound class the laboratory measured; (ii) kinetic coordinates from contrast rows
  only (ratios between the user's own conditions), restricted to coordinates whose Jacobian column
  is identifiable on those rows (rank test), bands from the shipped Laplace covariance. Everything
  else stays at the shipped value and the card says so. Contrast rows are derived automatically
  from the records: two records that differ in one condition make one contrast.
- Hold-out inside the user's data: records tagged `role: validate` in the spec, or, when untagged,
  every second time point held out, chosen before the fit and written to the card.
- The engine: `ProcessSpec` gains an optional `calibration`; `core_parameters` overlays the
  overrides; the observable layer applies the response factor; the envelope draws the factor's
  sigma. `compare`, `predict` and `score` take `--calibration <file>`; the report names it on every
  page.
- The shipped numbers never move: the panel scorecard is computed without any calibration, and a
  gate asserts no file under `results/user/` is read by the generators.

**Tests.** Synthetic data from the engine with a known threefold response factor and a known shift
in one rate: `calibrate` recovers both within their sigmas and leaves the other coordinates at the
shipped values. Validate-tagged rows never enter a residual (an AST-free runtime guard, like the
hold-out guard). The Yiltirak ladder as a pretend laboratory: factor from two temperatures, the
other two held out, hold-out median must improve. A calibration file from one laboratory cannot be
applied to a spec whose matrix differs from the records' matrix without a printed warning.

**Success.** A user with four measured pots gets a card that says what was pinned, what was not,
and by how much their hold-out improved. About two weeks. Chance it makes the tool useful to a
company: three in four.

## 2. The protein matrix as chemistry, not as a threshold

**Why.** Users cook on pea and soy isolates; the engine's protein disulfide sites and electrophile
pool are zero in every system it has ever run. The adduct dossiers already hold the rates.

**Design.**
- `data/species/protein_matrices.yml`: one entry per matrix (pea isolate, soy isolate, whey,
  casein, BLG): protein fraction, free thiol, disulfide, lysine and arginine sites in mmol per gram,
  each with its dossier anchor and a band. Specs gain `protein_g_per_l` (or `matrix` plus a stated
  fraction); `matrix: water` charges nothing.
- Engine species: `PROT_NH2` beside `PROT_SS`; the electrophile pool `MELE` charged from the
  matrix table and, in a later step, produced by browning.
- Declared reactions with declared bands, no fit: thiol plus protein disulfide (the bracket the
  engine has, now with a charged pool), aldehyde plus lysine (the K6b ladder, barriers 15 to 50
  kJ/mol), HMF plus cysteine (Hamzalıoğlu 2018's rate and barrier), thiol plus electrophile with
  Hofmann 2002's rate and plateau. Evidence class per parameter; the envelope samples the bands.
- Pre-registration B17. Hold-out: the four external matrix pots and the three isolate pots on the
  trust loop. Ship rule: the matrix rows' median fold error improves, and every free-precursor row
  moves less than 0.05 dex (a matrix that is absent must change nothing).

**Tests.** Sites charge to zero at `matrix: water`; sulfur and nitrogen balance close with the new
pools; the declared bands appear in the envelope's prior list; the ship rule as a scientific test.

**Success.** About three weeks. Chance the matrix rows land within tenfold: one in two. Known
limit: the binding data are at thiol loadings thousands of times above food levels, so the bands
will be wide and the card must say so.

## 3. The thiol sink (W7)

**Why.** The meaty character the users want is the thiols, and the model loses them too fast at
every temperature. The candidate table is written (`docs/validation/thiol_sink_candidates.md`).

**Design.** One pre-registration with two variants built and scored against the existing series
(Schieberle 2000 at 100 °C, Hofmann 1998 at 145 °C, Wang 2022 at 140 °C, Yiltirak's ladder):
(a) a saturable covalent sink on a pool browning makes, with Hofmann 2002's rate and plateau;
(b) the disulfide made reversible, or second-order with a depletable oxidant, with Kumazawa's
time-doubling as the check. Ship rule declared before the run, as for every wave. The experiment
in the introduction's section 8 decides between them; the two figure-only Chinese data sets are the
fallback if no laboratory is available.

**Success.** About one week of modelling plus the laboratory. Chance a variant ships: one in three
without the experiment, two in three with it.

## 4. Reach: install, call, click

**Why.** The tool runs only in a container from YAML through a shell script. Right for
reproducibility, wrong for adoption.

**Design.**
- `pip install maillard`: move the argument parser from `scripts/maillard.py` into `src/cli.py`
  and leave the script as a shim; declare `[project.scripts] maillard = "src.cli:main"`; ship
  `data/keys`, `data/species`, `data/lit/*.yml`, `data/benchmarks` and `results/validation` as
  package data (the PDFs are not tracked and not needed). The container stays the reproducible
  route and CI builds both.
- A documented Python API, three functions: `compare(spec_a, spec_b)`, `predict(spec)`,
  `explain(compound)`, returning the same payloads the verbs print.
- `maillard ui`: the standard library's HTTP server serving one page that turns a form into a spec
  and returns the existing HTML report. No new dependency, no accounts, runs on the user's machine.
- One JSON schema for specs (`data/schemas/spec.schema.json`) validated at the front door and used
  by the page to build its form; today validation is hand-rolled in the CLI.

**Tests.** A CLI test per verb against the schema; the API returns the same payload as the verb; the
page round-trips the template spec.

**Success.** About one week. No new science; the effect is on who can use the tool.

## 5. The chemistry users ask about that no lane has

Each starts as a rule in the hypothesis layer and becomes a wave only when a measured rate exists.
- **Pyrazines** (roasted, nutty): the aminoketone condensation after the Strecker step; the
  registry has the species. Rule first; a wave when a rate source is read.
- **The beany note before any heat**: hexanal from lipoxygenase during processing. Different
  chemistry, its own module and its own data programme; the lipid lane must not be stretched to it.
- **Extrusion**: two minutes at 130 to 170 °C at low moisture. The process spec already takes a
  programmed thermal history; what is missing is water-activity terms outside the measured windows
  and any data. Declare the windows, refuse outside them, collect data.
- **The refused routes** the wishlist names: oleate to nonanal, linoleate to 2-pentylfuran,
  aldehyde to alcohol. Rules first, so the refusal reads "no rate".

## 6. Cross-cutting engineering, done once

- **One overlay type** for calibrations, matrix declarations and waves, over the engine's existing
  override mapping, so a fit report, a user calibration and a declared matrix compose in a stated
  order (shipped wave, then matrix, then calibration) and the report prints the stack.
- **One spec schema** for every verb, with the units in the schema, not in comments.
- **Provenance on every user artifact** through `artifact_io`, as for the tracked ones; user
  artifacts are not under the freshness gate because they are not tracked.
- **A CLI test per verb** and a smoke test of the pip package in CI.
- **The molar masses** the structures test recorded, corrected in one change with the guards
  re-pinned, before the matrix layer reports any concentration through them.

## 7. Order, effort, and what would change it

| programme | effort | depends on | declared success |
|---|---|---|---|
| 1 calibrate | 2 weeks | nothing | synthetic recovery; Yiltirak hold-out improves |
| 6 overlay + spec schema | 1 week | do inside 1 | one mechanism used by 1 and 2 |
| 4 reach | 1 week | 6 | pip install, API, page |
| 2 matrix | 3 weeks | 6, the mass fix | matrix median improves; free rows unchanged |
| 3 sink | 1 week + laboratory | none | a variant passes its ship rule |
| 5 chemistry | rules now; waves when data exist | hypothesis layer | each rule cited and controlled |

Calibrate goes first because it makes every other gap something the user's own data can close.
If a laboratory can run the sink experiment soon, programme 3 moves to the front, since its result
changes what programme 2 should bind the thiols to.
