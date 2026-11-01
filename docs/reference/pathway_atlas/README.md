# Pathway-condition atlas

The atlas answers one question: **which Maillard reactions happen under which conditions, according to the literature.**

- It has 32 reaction families (`families.yml`), each with its own file `F01_*.yml` ... `F32_*.yml`.
- Each file holds **cells**. A cell is one study's observation of one family under one set of conditions. It records:
  - a status;
  - the kind of evidence;
  - the laboratory;
  - the dossiers in `data/lit/extraction_dossiers/` it was read from.
- The engine side lives in `engine_family_map.yml` (which family each engine reaction belongs to).
- The comparison of engine and literature lives in `ENGINE_VS_LITERATURE.md` / `engine_vs_literature.yml`.
- Other files cite a cell as `F25#6`, meaning the 0-based index in that file's `cells:` list. **Cell indices are stable:**
  - add a new cell at the end of its file;
  - delete a cell only when it is the last one, or update every reference to the cells after it.

## Cell rules

**Key and tags.** A cell is keyed on three axes: `ph`, `temperature` and `water`. Five further axes are tags: `sugar`, `amine`, `oxygen`, `buffer` and `metals`. The bins are defined in `families.yml`. A tag may be a list when the cell spans the arms of a contrast.

- **Temperature** is computed from `tags.T_C`:
  - low <= 60 C;
  - moderate > 60 to <= 100 C;
  - high > 100 to <= 140 C;
  - very_high > 140 C.

  A range or list ("80-120", "100 / 130") may sit in any bin it spans. Text in parentheses is detail. So 100 C is `moderate`, 101 C is `high`, 140 C is `high` and 145 C is `very_high`.
- **Water** is binned on all solutes, buffer included. A 0.5 M buffer on its own makes a pot `concentrated`. Record the buffer molarity as `tags.buffer_M`.
- **Oxygen** is `air` or `inert` only where the source states the headspace, the vessel (an open oven or roaster, an open reflux) or the sparge. A sealed vial with no stated gas is `not_stated`, even though it physically holds air.
- **Metals** is `chelated` when a chelator (DTPA, EDTA, DETAPAC) is in the pot. The Glomb lab's deaerated arms carry DTPA, so their air / deaerated contrasts are tagged `[not_controlled, chelated]`: O2 and chelation are confounded.
- **Buffer:** citrate-phosphate (McIlvaine) is `other_buffer` with a `buffer_note`, not `phosphate`.
- **Amine:** amines that are not amino acids are `other_amine`, not `other_free_aa`. Examples: propylamine, t-butylamine, piperidine, methylamine.

**Status is capped by the evidence kind.**

| Status | Meaning | Allowed when |
|---|---|---|
| `dominant` | the main route under these conditions | labelling or fed_intermediate, or trapping with a stated mass balance. A pot whose sugar is 100 % labelled shows only that a product is "sugar-derived", so it does not qualify. Neither does polymer composition. |
| `operating` | shown to occur; share unknown or not dominant | any evidence, provided the family's marker is seen |
| `minor` | shown to occur at a **measured** small share (< ~10 %) | a number or a measured bound. "Negligible, data not shown" is `operating`, not `minor`. |
| `absent` | looked for with a method able to see it, and not found | **requires `detection_basis:`**: a detection limit, a positive control for the same analyte in the same study, or a label pattern that would have shown the route. A positive control does not earn `absent` when the analyte is known to be lost faster than it accumulates (for example FFT in a hot pot, F28). |
| `not_detected` | not seen by a method whose sensitivity is not shown | the default for a non-detect. It never switches a step off: the product may be formed and then lost or captured. |
| `conflicting` | studies disagree | must name a registered CF id |
| `unattributed` | the product forms, but no marker says which family made it | — |

**Quotes.** Every cell carries a `quote` of at most 15 words, copied verbatim from a cited dossier.

## How reliable it is

**Red-team error rate.** On 2026-11-01 a red team graded 98 cells against their dossiers, and every `dominant` and `absent` cell also against the PDFs.

- **Per-cell WRONG rate:** about 9 %, from 6 of 68 random draws (95 % Wilson CI 4 - 18 %). Here WRONG means the status, evidence, family or quantity is materially wrong.
- **Any issue at all**, including tag and bin slips: about 60 %.
- **Worst status:** `absent`, at 5 of 16 sampled (31 %, CI 14 - 56 %).

**Fixes made on 2026-11-01:**

- The 10 WRONG cells were fixed.
- All 27 `absent` cells were re-audited: 14 kept `absent` with a `detection_basis`, 9 became `not_detected`, and 4 became `minor` or `operating` because the product was in fact seen.
- The bins were recomputed.

The roughly 460 cells nobody sampled still carry the pre-fix rate: expect about 40 of them to be wrong (95 % CI about 20 - 80).

**Fit to drive model changes now.** The errors in these cells are mostly tags:

- the `operating` / `dominant` labelling cells of F23 / F24 (CAMOLA 1:1 or positional labels: the intact-C5 vs C2 + C3 split);
- F16 / F17, the fed-dicarbonyl Strecker pots. Their `dominant` means "the only source in a sugar-free pot fed the donor"; do not carry it over to sugar pots;
- F19;
- the F01 - F06 early-chemistry cells;
- the 14 `absent` cells, each of which now states its detection basis.

**Not fit without re-audit:**

- any `not_detected` cell used to switch a step off;
- the F25 `dominant` cells outside Cerny's fed-furfural pot. The three 95 C Cerny cells are one lab and one design, and Cerny 2004's label position is inferred from Cerny 2003;
- F31 `dominant`;
- F30 quantities;
- single-family rates in general. With only 2 - 8 cells sampled per family, one family's rate says little.

**Known limitation: the gate checks the dossier, not the paper.**

- Many quotes are the dossier's own paraphrase, table re-typing or `[Z]` derivation. The gate proves only that the quote is in the dossier.
- Dossier arithmetic flows into cells unchecked. Examples the red team found:
  - Voigt's 0.9 mol/mol, which double-counts the two halves of one cleavage;
  - the zamora2015 dossier's "air multiplies PAC 4-8x". Table 1 gives x1.1 - 9.4.
- Before a cell changes the model, re-read its source table in the PDF.

## Running the gate

```bash
./scripts/docker_maillard.sh run "python scripts/ci/pathway_atlas_gate.py"
```

The gate (`scripts/ci/pathway_atlas_gate.py`) is structural and offline. It checks:

- the file / family ids and the key / tag bins;
- that the temperature key agrees with `T_C`;
- that a buffer of 0.5 M or more is not keyed `dilute`;
- the status / evidence vocabulary and the `dominant` cap;
- that every `absent` cell has a `detection_basis`;
- that each dossier exists, with exact case;
- that each quote appears verbatim in a cited dossier;
- that each CF id is registered.

It does not judge whether a cell reads its paper correctly.
