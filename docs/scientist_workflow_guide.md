# Scientist Workflow Guide: Flavour Design with Maillard

This guide provides a minimal starting point for a food scientist to use the Maillard framework for plant-based meat design.

## 1. Rank Your Formulations
To compare different meat-like formulations (e.g., varying ribose/cysteine levels across Pea and Soy):

```bash
./scripts/docker_maillard.sh matrix-family-ranking
```
- **Read:** `results/validation/matrix_family_priority_ranking.md`
- **Goal:** Identify which formulation has the highest "Meaty" aroma signal vs. "Beany" off-notes.

## 2. Inspect Trust & Confidence
Before making a decision based on the ranking, check the evidence quality:

```bash
./scripts/docker_maillard.sh explain-formulation MY_FORM_NAME
```
- **Review:** `results/validation/what_can_i_trust_today.md`
- **Check:** Are the key ingredients (e.g., Ribose, Cysteine) being modeled with `literature_anchor` or only `directional_priors`?

## 3. Assess Process Damage & Safety
Check if your process temperature or shear levels are damaging the protein quality:

```bash
./scripts/docker_maillard.sh literature-learning-loop
```
- **Read:** Section on **Family 12 (Protein Damage Markers)**.
- **Goal:** Ensure Lysine loss and AGE levels remain within the "Green" safety zone.

## 4. Off-Note Suppression
If your formulation uses high-fat plant proteins, check the lipid-oxidation crosstalk:

```bash
./scripts/docker_maillard.sh hexanal-nonanal-resolution
```
- **Read:** `results/validation/hexanal_nonanal_resolution.md`
- **Goal:** Confirm that your Maillard intermediates are effectively quenching lipid-derived off-notes.

---
*For a unified view, run:*
```bash
./scripts/docker_maillard.sh summary
```
