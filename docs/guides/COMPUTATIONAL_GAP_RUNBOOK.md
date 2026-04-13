# Computational Gap Runbook

This guide is for the three active selective-QM targets that currently add the most value to Maillard:

- `lysinoalanine_crosslink`
- `aa_ring_open_dicarbonyl`
- `hexanal_radical_quench`

Use this when you want exact commands that can be copied and pasted from the repository root.

## What Systems This Runbook Covers

Not every entry in the computational-gap plan is the same kind of thing.

Some entries are:

- direct elementary-step reaction systems that can be run through xTB and then DFT
- bounded surrogate systems where DFT narrows a prior but does not replace the current literature authority
- broader calibration questions that are not a single copy-paste QM job

This runbook is mainly for the first category.

## Systems You Can Replicate Directly Today

These are the targets that are currently prepared as local xTB/DFT jobs with seeded geometries.

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

### `hexanal_radical_quench`

- Family: 11
- Chemistry: quenching of hexanal by cysteine-derived sulfur radical chemistry
- Why it matters: this is the cleanest off-note suppression target in the current selective-QM set
- Runtime role today: a successful xTB plus DFT refinement supports ranking-only claims about hexanal suppression, not full benchmark closure
- Practical note: this is often the easiest one to explain to users because it maps directly to hexanal control

## Systems In The Same Program But Not Copy-Paste Ready Yet

These belong to the same selective-QM roadmap, but they are not yet at the same "just copy and run" level.

### `quinone_cys_michael`

- Family: 13
- Chemistry: cysteine capture by polyphenol-derived quinones through Michael addition
- Why it matters: directly affects cysteine depletion in wet, polyphenol-rich systems, so it is a strong next candidate after the current three
- Current blocker: local balanced xTB seed geometries are not prepared yet

### `pe_schiff_base`

- Family: 15
- Chemistry: phosphatidylethanolamine plus sugar condensation at an interface
- Why it matters: gives a route to model sugar loss and interfacial phospholipid-amine chemistry in lecithin-rich systems
- Current blocker: today it is represented as a literature-transfer prior, not a ready local xTB/DFT run

### `pe_amadori`

- Family: 15
- Chemistry: rearrangement of the PE Schiff base toward a PE-Amadori product
- Why it matters: extends the same phospholipid chemistry lane beyond the initial condensation step
- Current blocker: same as above; it is not yet seeded as a ready copy-paste local job

### `asparagine_sugar_explicit_water_cluster`

- Family: 13 safety lane
- Chemistry: explicit-water cluster around a proton-transfer-sensitive asparagine-sugar transition-state motif
- Why it matters: this is the main safety-oriented explicit-solvent target because it can narrow acrylamide uncertainty
- Current blocker: it needs explicit-water cluster seed generation, so it is not a normal xTB path-search target yet

### `melanoidin_radical_trapping`

- Family: 16
- Chemistry: sulfur-radical trapping by high-molecular-weight melanoidin proxies
- Why it matters: this helps explain sulfur loss as browning progresses
- Current blocker: this lane is currently screened as an MLP-first problem and is not in the same ready xTB/DFT copy-paste scope as the three active targets

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

## 1. xTB Commands

Run these first if you need to regenerate `xtbpath.xyz` and `xtbpath_ts.xyz`.

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

### Hexanal Radical Quench

```bash
conda run -n maillard python scripts/run_computational_gap_xtb.py \
  --target hexanal_radical_quench \
  --execute \
  --output results/computational_gap_refinement/hexanal_radical_quench_xtb_execution.json
```

## 2. DFT Commands

Run one heavy DFT job at a time unless you intentionally want to saturate the machine.

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

### Hexanal Radical Quench

```bash
conda run -n maillard python scripts/run_computational_gap_dft.py \
  --target hexanal_radical_quench \
  --execute \
  --heartbeat-seconds 300 \
  --output results/computational_gap_refinement/hexanal_radical_quench_dft_execution.json
```

## 3. Recommended Order

If you are running all three by hand, use this order:

1. `lysinoalanine_crosslink`
2. `aa_ring_open_dicarbonyl`
3. `hexanal_radical_quench`

Reason:

- `lysinoalanine_crosslink` is the most important damage-marker target and has already shown SCF fragility.
- `aa_ring_open_dicarbonyl` adds bounded upstream browning value.
- `hexanal_radical_quench` is high-value for off-note suppression, but the first two currently do more to tighten the broader scientific surface.

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

## 7. Current Best Candidate After These Three

If you want a fourth target after the current three, the next best candidate is `quinone_cys_michael` because it adds Family 13 coverage and directly affects cysteine depletion in wet, polyphenol-rich systems.

Today it is not seed-ready yet, so it is not a copy-paste run like the three targets above.