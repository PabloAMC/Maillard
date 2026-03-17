# External Matrix Benchmark Unlock Report

## Purpose

This report defines what the repository still needs in order to unlock an honest external scientific assessment for pea and soy matrix benchmarks, especially for meaty-positive systems. The current codebase already separates three states cleanly:

- externally measured matrix off-flavour anchors
- internally reproducible matrix meaty-positive harnesses
- strict-gate-eligible free-precursor benchmarks

What is still missing is not more code first. It is the minimum external evidence package needed to convert matrix meaty-positive assessment from directional support into benchmark-backed scientific assessment.

## What Is Already In Place

- Pea and soy ambient-slurry off-flavour anchors exist through Pratap-Singh 2021.
- Pea and soy meaty-positive matrix harnesses exist as internal reproducibility candidates based on frozen Docker outputs.
- Matrix benchmark evidence is now auditable with the repo tooling.
- Matrix benchmark assertions are now auditable in Docker without opening the strict gate.
- Branch-vs-main matrix benchmark deltas are now auditable before any promotion claim.

## The Scientific Gap

The current blocker is external quantitative evidence for meaty-positive compounds in realistic plant-protein matrices.

Today the repo does not contain an external wet-lab benchmark for pea or soy that quantitatively measures compounds like:

- 2-furfurylthiol
- 2-methyl-3-furanthiol
- bis(2-methyl-3-furyl) disulfide
- 2,5-dimethylpyrazine

under a controlled plant-matrix process state that is comparable to the modelled matrix_precursor_augmented path.

## Minimum Evidence Package Needed

### 1. External Meaty-Positive Matrix Benchmark Pair

At minimum, one pea benchmark and one soy benchmark are needed.

Required structure:

- protein matrix: pea isolate and soy isolate separately
- precursor system: ribose plus cysteine, or another justified meaty-positive system with a clear sulfur route
- controlled process state: e.g. aqueous pre-extrusion slurry, preheated slurry, or post-extrusion hydrated matrix
- fixed pH, water activity, time, and temperature
- replicate quantitative GC-MS or equivalent output

Required outputs:

- absolute or well-normalized concentrations for the meaty-positive compounds above
- at least one adverse or tradeoff marker in the same experiment
- compound-level measurement uncertainty

### 2. Process-State Definition That Maps To The Model

The current code can express process states, but external assessment needs the wet-lab side to define them in a way that maps cleanly into benchmark metadata.

Each external benchmark must specify:

- matrix format
- preheating history
- extrusion history
- moisture state
- whether precursors are exogenous, endogenous, or both
- denaturation or process proxy if directly measured

Without that mapping, any benchmark would be structurally ambiguous even if the compound measurements were good.

### 3. Shared Meaty vs Off-Flavour Tradeoff Measurement

External assessment should not validate only desirable compounds in isolation.

Each candidate external meaty-positive matrix benchmark should also include at least one of:

- hexanal
- 2-pentylfuran
- nonanal
- another justified lipid-coupled off-note marker

This is necessary because P1 is not only about executability or desirable ranking. It is about scientist-useful tradeoff assessment inside a matrix benchmark.

### 4. Reproducibility Thresholds

The external dataset should not be accepted into the benchmark set unless it includes reproducibility information.

Minimum reproducibility package:

- biological or process replicates
- analytical replicate description
- uncertainty or confidence interval per compound
- explicit non-detect policy
- preprocessing and normalization notes

Without this, the benchmark can exist as reference material, but not as a strong external assessment anchor.

### 5. Data Packaging Contract

The benchmark file should include:

- benchmark_id
- source_metadata with lab, date, report or DOI
- process_metadata
- matrix_ranking_contract
- measured_volatiles
- uncertainty fields
- protein_type and denaturation assumptions when applicable

The comparator signal must be measured_volatiles, not reference_volatiles.

## Recommended Wet-Lab Design

### Core Design

- matrices: pea isolate, soy isolate
- precursor load: 1 mM ribose + 1 mM cysteine as the initial narrow unlock case
- conditions: one moderate aqueous-heated state and one more severe process state
- readouts: meaty sulfur targets plus lipid-derived adverse markers
- method: HS-SPME GC-MS or equivalent compound-resolved method with calibration notes

### Minimum Compound Panel

- 2-furfurylthiol
- 2-methyl-3-furanthiol
- bis(2-methyl-3-furyl) disulfide
- 2,5-dimethylpyrazine
- furfural
- hexanal
- 2-pentylfuran

### Why This Panel

- it anchors the current internal meaty harnesses
- it links meaty gain to off-flavour tradeoff in the same experiment
- it is narrow enough to be realistic for a first unlock study

## Acceptance Criteria For External Assessment Unlock

The repo should only consider the external matrix benchmark assessment unlocked when all of the following are true:

- pea has an external quantitative meaty-positive matrix benchmark
- soy has an external quantitative meaty-positive matrix benchmark
- each family still has an off-flavour anchor
- process metadata maps cleanly to a supported model state
- the benchmark passes the matrix ranking assertions in Docker
- branch-vs-main matrix deltas are generated before any promotion claim

Strict gate promotion should still remain disabled until those benchmarks are stable across repeated runs and accepted by the benchmark contract.

## Practical Repo Work Needed Once Data Exists

Once the external data exists, the code work is straightforward:

1. Add new benchmark JSON files with measured_volatiles.
2. Update matrix calibration entries only where the data truly supports compound-specific observability.
3. Regenerate matrix evidence, assertions, readiness, and branch-delta artifacts.
4. Re-run the Docker benchmark lanes.
5. Only then consider changing promotion status or strict-gate policy.

## What Should Not Be Done

- Do not relabel internal reference_volatiles as external data.
- Do not promote matrix_precursor_augmented benchmarks into the strict gate on reproducibility alone.
- Do not treat off-flavour-only anchors as sufficient evidence for meaty-positive external assessment.

## Operational Commands

Current repo commands relevant to this workflow:

```bash
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh matrix-branch-deltas main
./scripts/docker_maillard.sh coverage-gaps
```

These commands do not unlock the science by themselves. They make the current boundary explicit and reproducible.
