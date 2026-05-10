# Computational Gap Runbook

This guide is for the official selective-QM queue plus the closely related proxy lanes that are being triaged for later promotion.

The current official DFT-ready queue is:

- `hexanal_radical_quench`
- `lysinoalanine_crosslink`
- `aa_ring_open_dicarbonyl`
- `quinone_cys_michael`
- `pe_schiff_base`
- `pe_amadori`

The current proxy or prep-only companion lanes are:

- `asparagine_sugar_explicit_water_cluster`
- `melanoidin_radical_trapping`

Use this when you want exact commands that can be copied and pasted from the repository root.

## What Systems This Runbook Covers

Not every entry in the computational-gap plan is the same kind of thing.

Some entries are:

- direct elementary-step reaction systems that can be run through xTB and then DFT
- bounded surrogate systems where DFT narrows a prior but does not replace the current literature authority
- broader calibration questions that are not a single copy-paste QM job

This runbook is mainly for the first category.

## Systems You Can Replicate Directly Today

These are the targets that are currently prepared as local xTB/DFT jobs with seeded geometries and are part of the official selective-QM queue.

### `hexanal_radical_quench`

- Family: 11
- Chemistry: quenching of hexanal by cysteine-derived sulfur radical chemistry
- Why it matters: this is the cleanest off-note suppression target in the current selective-QM set
- Runtime role today: a successful xTB plus DFT refinement supports ranking-only claims about hexanal suppression, not full benchmark closure
- Practical note: this is often the easiest one to explain to users because it maps directly to hexanal control

### `lysinoalanine_crosslink`

- Family: 12
- Chemistry: dehydroalanine plus lysine epsilon-amine addition to form lysinoalanine
- Why it matters: this is the protein-damage and crosslinking lane; it affects safety and digestibility interpretation under alkaline or high-severity processing
- Runtime role today: the repo still keeps the DHA-crosslinking family rule as the bounded fallback authority, and selective DFT is meant to narrow that prior
- Practical note: this is the most important damage-marker target, but it has also been the most fragile SCF case

### `aa_ring_open_dicarbonyl`

- Family: 14
- Chemistry: dehydroascorbic acid hydration and ring opening toward a reactive dicarbonyl surrogate state
- Why it matters: this controls one upstream browning pressure lane from ascorbic-acid chemistry
- Runtime role today: the live runtime still treats the literature-backed Family 14 source term as the bounded authority until selective DFT tightens it
- Practical note: this is useful when you want to reduce uncertainty around stealth browning and dicarbonyl supply

### `quinone_cys_michael`

- Family: 13
- Chemistry: cysteine capture by polyphenol-derived quinones through Michael addition
- Why it matters: directly affects cysteine depletion in wet, polyphenol-rich systems
- Runtime role today: the lane is active as a bounded upstream precursor sink and the xTB path is now seeded locally
- Practical note: this is the next formal Family 13 candidate after the original three-target batch

### `pe_schiff_base`

- Family: 15
- Chemistry: phosphatidylethanolamine plus sugar condensation at an interface
- Why it matters: gives a route to model sugar loss and interfacial phospholipid-amine chemistry in lecithin-rich systems
- Runtime role today: the lane now has a managed xTB execution artifact and is eligible for bounded-calibration selective DFT
- Practical note: run this before `pe_amadori` if you want the Family 15 pair to stay mechanistically ordered

### `pe_amadori`

- Family: 15
- Chemistry: rearrangement of the PE Schiff base toward a PE-Amadori product
- Why it matters: extends the same phospholipid chemistry lane beyond the initial condensation step
- Runtime role today: the lane now has a managed xTB execution artifact and is eligible for bounded-calibration selective DFT alongside the Schiff-base step
- Practical note: keep the Family 15 pair together in queue planning even if you execute them in separate DFT runs

## Systems In The Same Program But Not Copy-Paste Ready Yet

These belong to the same selective-QM roadmap, but they are not yet at the same "just copy and run" level.

### `asparagine_sugar_explicit_water_cluster`

- Family: 13 safety lane
- Chemistry: explicit-water cluster around a proton-transfer-sensitive asparagine-sugar transition-state motif
- Why it matters: this is the main safety-oriented explicit-solvent target because it can narrow acrylamide uncertainty
- Current blocker: it needs explicit-water cluster seed generation, so it is not a normal xTB path-search target yet

### `melanoidin_radical_trapping`

- Family: 16
- Chemistry: sulfur-radical trapping by high-molecular-weight melanoidin proxies
- Why it matters: this helps explain sulfur loss as browning progresses
- Current blocker: this lane is currently screened as an MLP-first problem and is not in the same ready xTB/DFT copy-paste scope as the official queue

## Systems In The Computational-Gap Plan That Are Not Single QM Jobs

The plan also contains entries such as sulfur selectivity, peptide-bound cysteine reactivity, matrix aldehyde release, and broader retention-calibration questions.

Those are important, but they are not single local xTB-plus-DFT jobs with one reaction directory and one output JSON. They are better understood as model-calibration or triage problems, not as copy-paste replication jobs.

## Before You Start

Requirements:

- You are in the repository root.
- The `maillard` conda environment exists.
- The xTB seeds are already present under `data/geometries/xtb_inputs/`.

Quick verification:

```bash
pwd
conda run -n maillard python --version
```

## Important Rule

Use a separate output JSON per target.

Do not reuse the same DFT output path for multiple targets if you want to preserve the status of each run.

New control rule:

- run the DFT preflight first for every target before starting a heavy calculation
- only launch the heavy DFT job if the execution JSON reports `preflight.quality_gate_passed = true`
- monitor `phase`, `phase_history`, and `elapsed_seconds` in the target-specific execution JSON instead of treating `running` as a black box

## Recommended Docker Workflow

### 1. Preflight First

```bash
./scripts/docker_maillard.sh computational-gap-dft-preflight aa_ring_open_dicarbonyl
```

This writes `results/computational_gap_refinement/aa_ring_open_dicarbonyl_dft_execution.json` without starting the expensive QM run.

What to inspect in that JSON:

- `status`: should be `ready_for_dft`, not `blocked_bad_ts_guess`
- `preflight.blockers`: should be empty
- `preflight.reactant_min_interatomic_distance_angstrom`
- `preflight.ts_min_interatomic_distance_angstrom`
- `preflight.pairwise_distance_rms_delta_angstrom`
- `preflight.pairwise_distance_max_delta_angstrom`

If the status is `blocked_bad_ts_guess`, do not launch DFT yet. Rebuild or replace the TS seed first.

### 2. Launch the Heavy DFT Run

```bash
./scripts/docker_maillard.sh computational-gap-dft aa_ring_open_dicarbonyl
```

The target-specific execution JSON now records live phases such as:

- `preflight_passed`
- `reactant_geometry_optimization`
- `reactant_solvent_scf`
- `reactant_hessian`
- `reactant_single_point`
- `ts_geometry_optimization`
- `ts_solvent_scf`
- `ts_hessian`
- `ts_single_point`
- `barrier_evaluation`
- `completed`

### 3. Decide Whether to Restart

Restart the job when one of these is true:

- the active run predates phase tracking and you need controlled visibility
- the preflight reports `blocked_bad_ts_guess`
- the execution JSON stalls without `updated_at` moving forward
- the phase history gets stuck in an early geometry-optimization phase and the seed is scientifically suspect

For a restart, stop the old terminal/process that owns the run and relaunch through the Docker command above so the target gets a fresh preflight and a clean phase-tracked execution JSON.

## 1. xTB Commands

Run these first if you need to regenerate `xtbpath.xyz` and `xtbpath_ts.xyz` for the official queue.

### Quinone Cys Michael

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target quinone_cys_michael \
  --execute \
  --output results/computational_gap_refinement/quinone_cys_michael_xtb_execution.json
```

### Hexanal Radical Quench

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target hexanal_radical_quench \
  --execute \
  --output results/computational_gap_refinement/hexanal_radical_quench_xtb_execution.json
```

### Lysinoalanine Crosslink

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target lysinoalanine_crosslink \
  --execute \
  --output results/computational_gap_refinement/lysinoalanine_crosslink_xtb_execution.json
```

### AA Ring-Open Dicarbonyl

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target aa_ring_open_dicarbonyl \
  --execute \
  --output results/computational_gap_refinement/aa_ring_open_dicarbonyl_xtb_execution.json
```

### PE Schiff Base

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target pe_schiff_base \
  --execute \
  --output results/computational_gap_refinement/pe_schiff_base_xtb_execution.json
```

### PE Amadori

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target pe_amadori \
  --execute \
  --output results/computational_gap_refinement/pe_amadori_xtb_execution.json
```

## 2. DFT Commands

Run one heavy DFT job at a time unless you intentionally want to saturate the machine.

### Quinone Cys Michael

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target quinone_cys_michael \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/quinone_cys_michael_dft_execution.json
```

### Hexanal Radical Quench

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target hexanal_radical_quench \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/hexanal_radical_quench_dft_execution.json
```

### Lysinoalanine Crosslink

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target lysinoalanine_crosslink \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/lysinoalanine_crosslink_dft_execution.json
```

### AA Ring-Open Dicarbonyl

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target aa_ring_open_dicarbonyl \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/aa_ring_open_dicarbonyl_dft_execution.json
```

### PE Schiff Base

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target pe_schiff_base \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/pe_schiff_base_dft_execution.json
```

### PE Amadori

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target pe_amadori \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/pe_amadori_dft_execution.json
```

## 3. Recommended Order

If you are running the official queue by hand, use this order:

1. `aa_ring_open_dicarbonyl`
2. `hexanal_radical_quench`
3. `quinone_cys_michael`
4. `pe_schiff_base`
5. `pe_amadori`
6. `lysinoalanine_crosslink`

Reason:

- `aa_ring_open_dicarbonyl` is already running in the current queue and is the cleanest bounded-calibration tightening step.
- `hexanal_radical_quench` is the next easiest ranking-only off-note target once the Family 14 lane finishes.
- `quinone_cys_michael` broadens the formal queue into Family 13 without needing the more fragile Family 15 interface proxies first.
- `pe_schiff_base` and `pe_amadori` are now formal candidates because the Family 15 pair gate passed under managed xTB artifacts, but they still sit behind the currently running Family 14 job and the clearer Family 11/13 expansions.
- `lysinoalanine_crosslink` remains highly important, but it is still the most fragile SCF case.

## 4. Check Status

### Check xTB Status

```bash
conda run -n maillard python -c "import json; p=json.load(open('results/computational_gap_refinement/lysinoalanine_crosslink_xtb_execution.json')); print(json.dumps(p['summary'], indent=2)); print(json.dumps(p['jobs'][0], indent=2))"
```

Replace the filename with the target you want to inspect.

### Check DFT Status

```bash
conda run -n maillard python -c "import json; p=json.load(open('results/computational_gap_refinement/lysinoalanine_crosslink_dft_execution.json')); print(json.dumps(p['summary'], indent=2)); print(json.dumps(p['jobs'][0], indent=2))"
```

What the key states mean:

- `running`: the DFT job is live
- `completed`: the DFT job finished and produced a barrier
- `failed`: the DFT calculation crashed or did not converge
- `seed_required`: the required local geometry is missing
- `blocked_missing_xtb_ts`: xTB did not provide `xtbpath_ts.xyz`

## 5. Common Mistake

Wrong:

```bash
conda run -n maillard python run_computational_gap_dft.py --target lysinoalanine_crosslink --execute
```

Correct:

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py --target lysinoalanine_crosslink --execute
```

The script lives under `scripts/`, not in the repository root.

## 6. After DFT Finishes

Once you have completed DFT runs that you want to ingest, use:

```bash
conda run -n maillard python scripts/generators/ingest_computational_gap_dft_results.py \
  --manifest results/computational_gap_refinement/computational_gap_dft_job_manifest.json \
  --execution results/computational_gap_refinement/computational_gap_dft_execution.json \
  --output-dir results/validation
```

Then promote completed non-fast results:

```bash
conda run -n maillard python scripts/generators/promote_computational_gap_dft_results.py \
  --execution results/computational_gap_refinement/computational_gap_dft_execution.json \
  --output-dir results/validation
```

If you ran targets with separate output JSONs, point `--execution` to the file you actually want to ingest or merge first.

## 7. Current Promotion Rule For The Proxy Lanes

- `pe_amadori` can remain as a direct xTB proxy, but it should not enter the formal DFT queue until `pe_schiff_base` also materializes a stable TS guess.
- `pe_schiff_base` is still informative as a Family 15 stress test, but a normal-termination xTB run without a usable `xtbpath_ts.xyz` is not enough for promotion.
- `asparagine_sugar_explicit_water_cluster` stays outside the normal queue until the explicit-water seed exists.
- `melanoidin_radical_trapping` stays MLP-first unless a simpler fragment proxy becomes clearly superior.

Today it is not seed-ready yet, so it is not a copy-paste run like the three targets above.