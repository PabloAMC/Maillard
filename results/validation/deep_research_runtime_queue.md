# Deep Research Runtime Queue

Batch: 2026-04-05-runtime-first-batch-03
Source: results/literature/deep_research_backlog.json
Selection policy: Select only Deep Research citations still marked BACKLOG and stage them into non-benchmark runtime lanes: process_state_calibration, safety_reference, and computational_prior. Exclude anything already landed in runtime registries, already runtime-bound, or benchmark-first.

Selected runtime candidates: 0
Process-state calibrations: 0
Safety references: 0
Computational priors: 0
Already landed in runtime registries: 5
Excluded because already encoded or absent: 6

Out of scope for this batch: benchmark payload promotion.

## Process-State Calibration

Landing registry: data/lit/process_state_calibrations.json

| Citation | Families | Score | Theme | Target Modules | Next Action |
| --- | --- | ---: | --- | --- | --- |
| none | n/a | 0 | n/a | n/a | no curated candidates selected for this lane |

## Safety Reference

Landing registry: data/lit/safety_reference_payloads.json

| Citation | Families | Score | Theme | Target Modules | Next Action |
| --- | --- | ---: | --- | --- | --- |
| none | n/a | 0 | n/a | n/a | no curated candidates selected for this lane |

## Computational Prior

Landing registry: data/lit/computational_priors.json

| Citation | Families | Score | Theme | Target Modules | Next Action |
| --- | --- | ---: | --- | --- | --- |
| none | n/a | 0 | n/a | n/a | no curated candidates selected for this lane |

## Excluded Candidates

| Citation | Reason |
| --- | --- |
| Ordoudi et al. (2014 / PMC12484514) | already_landed_in_runtime_registry |
| Hrncirik & Zeelenberg (2014) | already_landed_in_runtime_registry |
| Aliani & Farmer (2005) | already_landed_in_runtime_registry |
| Blank, Devaud & Grosch (2003) | already_runtime_bound |
| Glomb & Monnier (1995) | already_landed_in_runtime_registry |
| Hidalgo & Zamora (2004) | already_landed_in_runtime_registry |
