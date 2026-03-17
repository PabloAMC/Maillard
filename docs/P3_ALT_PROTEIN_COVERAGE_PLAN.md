# P3 Alt-Protein Coverage Plan

## Objective

P3 should broaden the validated envelope from free-precursor sulfur chemistry plus narrow pea/soy intake modelling into something genuinely useful for alternative-protein R&D.

## First Practical Step

The first step is not to add chemistry blindly. It is to make the coverage gaps explicit and reproducible.

The repo now starts P3 with a generated coverage-gap artifact:

```bash
./scripts/docker_maillard.sh coverage-gaps
```

This artifact is intended to answer three questions before new scientific work is added:

- which protein systems are still completely unbenchmarked
- which matrix process states are still unrepresented
- whether the critical external meaty-positive matrix gap has actually been closed

## Recommended P3 Sequence

1. Close external matrix meaty-positive benchmark gaps for pea and soy.
2. Add at least one broader plant-protein family beyond isolate-only pea/soy.
3. Add a higher-severity process regime that is still benchmark-backed.
4. Expand benchmark compound panels so lipid/Maillard tradeoffs are co-measured.
5. Only then broaden chemistry families further.

## Near-Term Implementation Targets

### Track A: Coverage Expansion

- mycoprotein benchmark candidate definition
- pea_conc or soy_conc benchmark candidate definition
- one heated_matrix external anchor

### Track B: Chemistry Breadth

- broader carbohydrate systems relevant to commercial formulations
- peptide-bound and intact-protein reactivity where it changes predictions materially
- stronger lipid/Maillard coupling validation

### Track C: Safety And Time

- dynamic acrylamide trajectories across realistic heating profiles
- endpoint vs process-path validation for alternative heating schedules

## Success Condition

P3 is meaningfully underway when the repo can show a versioned artifact proving that benchmark coverage is expanding across:

- protein systems
- process states
- matrix tradeoff chemistry
- external benchmark evidence

The coverage-gap report is the first reproducible surface for that work.
