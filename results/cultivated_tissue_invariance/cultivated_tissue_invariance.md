# Does the kinetic layer change what to add to cultivated tissue?

*Generated 2026-09-13 by `src/cultivated_tissue_invariance.py`. Pre-registration: `results/validation/cultivated_tissue_invariance_prereg.md`. The composition box is a sensitivity device; every range is a stub, none is a measurement.*

## Verdict, by programme

| programme | verdict | draws evaluated | engine refused | < 2 restorable | top(E)=top(N) | mean τ | disagreements | dominant reversal | its envelope survival |
|---|---|---|---|---|---|---|---|---|---|
| 100C_20min | **indeterminate** | 191 / 200 | 0 | 9 | 60 % | 0.353 | 40 % | cysteine over glucose (27 % of disagreements) | 100 % |
| 140C_5min | **indeterminate** | 191 / 200 | 0 | 9 | 62 % | 0.330 | 38 % | cysteine over glucose (32 % of disagreements) | 100 % |

Read against section 3 of the pre-registration: T1 needs ≥ 20 % disagreement, one reversal carrying ≥ 50 % of it, and ≥ 80 % envelope survival; T2 needs ≥ 90 % top agreement and mean τ ≥ 0.75; T3 is ≥ 50 % of draws refused. Anything else is indeterminate and is not a licence to build.

## Where each precursor lands

**100C_20min** (mean position, 1 = restore first; over the draws the engine answered)

| precursor | naive N | engine E |
|---|---|---|
| ribose | 1.576 | 1.564 |
| cysteine | 2.597 | 1.903 |
| thiamine | 2.346 | 2.000 |
| glucose | 2.153 | 3.093 |

Disagreements, engine's winner over naive winner: cysteine over glucose ×21; ribose over glucose ×15; thiamine over ribose ×13; cysteine over thiamine ×11; cysteine over ribose ×11; thiamine over glucose ×3; ribose over thiamine ×2; thiamine over cysteine ×1.

**140C_5min** (mean position, 1 = restore first; over the draws the engine answered)

| precursor | naive N | engine E |
|---|---|---|
| ribose | 1.576 | 1.378 |
| cysteine | 2.597 | 1.552 |
| thiamine | 2.346 | 2.774 |
| glucose | 2.153 | 2.933 |

Disagreements, engine's winner over naive winner: cysteine over glucose ×23; cysteine over thiamine ×16; cysteine over ribose ×15; ribose over glucose ×14; ribose over thiamine ×4; glucose over thiamine ×1.

## What the engine refused

Structural, before any draw:

- trunk arm (methional, pyrazines, furaneol on glucose + glycine): methional refused, pyrazines ~1e-13 ug/L, no water threshold on any target; no decision metric (A1)
- cultivated fat: no route from tissue lipid to any lane (A1)
- hexanal: needs a lipid carrier; no tissue-fat carrier exists (A1)
- leucine, IMP, ribose-5-phosphate: not species in any core lane (A2)

## The composition box that was swept

| tissue | precursor | lo mM | hi mM | label | note |
|---|---|---|---|---|---|
| beef | ribose | 0.3 | 5.0 | stub | post-mortem ribose from IMP breakdown in aged beef; order of magnitude only |
| beef | cysteine | 0.05 | 0.5 | stub | free cysteine in raw muscle; order of magnitude only |
| beef | thiamine | 0.001 | 0.01 | stub | beef thiamine ~0.05-0.15 mg/100 g; order of magnitude only |
| beef | glucose | 1.0 | 10.0 | stub | free glucose in post-mortem muscle; order of magnitude only |
| beef | leucine | 0.3 | 1.5 | stub | free leucine; unrankable by the engine, reported for the naive ranking only |
| beef | IMP | 1.0 | 8.0 | stub | inosine monophosphate in aged beef; unrankable by the engine |
| beef | ribose-5-phosphate | 0.01 | 0.1 | stub | unrankable by the engine |
| cultivated_muscle | ribose | 0.01 | 2.0 | stub | no ageing step; IMP-to-ribose route depends on post-harvest handling; deliberately wide |
| cultivated_muscle | cysteine | 0.02 | 0.5 | stub | medium-fed cells; deliberately wide |
| cultivated_muscle | thiamine | 0.0003 | 0.01 | stub | DMEM carries ~12 uM thiamine; intracellular pool unknown; deliberately wide |
| cultivated_muscle | glucose | 0.1 | 10.0 | stub | depends on the harvest wash; deliberately wide |
| cultivated_muscle | leucine | 0.2 | 3.0 | stub | medium is leucine-rich; unrankable by the engine |
| cultivated_muscle | IMP | 0.1 | 3.0 | stub | unrankable by the engine |
| cultivated_muscle | ribose-5-phosphate | 0.01 | 0.1 | stub | unrankable by the engine |

Fixed in every spec: pH 6.0, a_w 0.98, matrix water, phosphate 0.03 mol/L. 200 draws per programme, seed 0, 50 envelope draws per reversal.
