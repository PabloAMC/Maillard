# Maillard CSV ingest templates

These templates are the canonical starting point for scientists who want to
add a new wet-lab measurement to the calibration / external-validation
surface via `./scripts/docker_maillard.sh ingest`.

## Files

- [hs_spme_gc_ms_template.csv](hs_spme_gc_ms_template.csv) — headspace
  SPME-GC-MS quantitation (one row per detected compound). The most common
  shape for plant-protein flavor work.
- [sida_template.csv](sida_template.csv) — stable-isotope-dilution analysis,
  typically a single odorant with tighter uncertainty.
- (sensory ranking template — coming in a follow-up sprint, see Lane H.)

## How to use

1. Copy the template to a new file (e.g. `data/intake/my_run.csv`).
2. Replace the example rows with your own data. Header names can be loose:
   the ingest CLI fuzzy-matches `Analyte Name`, `Observed ppb`, `RSD %`,
   `Temperature C`, etc. against the canonical schema.
3. Run a dry-run preview (no write):

   ```bash
   ./scripts/docker_maillard.sh ingest \
     --file data/intake/my_run.csv \
     --protein-type pea_iso \
     --process-state heated_matrix
   ```

4. Inspect `*_ingest_preview.md` and `*_support_delta.md` in the output dir.
5. If the impact preview matches what you expected, re-run with `--confirm`
   to persist the canonical YAML payload.

## Mandatory fields

| Field | Aliases the parser accepts | Notes |
| --- | --- | --- |
| `compound` | `compound`, `analyte`, `volatile`, `marker`, `name` | Case-insensitive. Use the IUPAC or common name; the ingest is OK with either. |
| `conc_ppb` | `concentration ppb`, `observed ppb`, `headspace ppb`, `value ppb` | μg/kg of finished matrix. |
| `temp_C` | `temperature c`, `temp`, `temperature deg c` | Cooking / processing temperature, not GC oven temperature. |
| `time_min` | `time min`, `duration min`, `minutes`, `time` | Cooking time at `temp_C`. |
| `ph` | `ph`, `pH` | Maillard pathways are pH-sensitive; do not omit. |
| `water_activity` | `aw`, `water activity` | Required for every row (the ingest treats it as mandatory). Use the post-processing matrix Aw, not raw slurry. |
| `source_reference` | `reference`, `citation`, `source` | Free-text citation; DOI optional. |

## Optional but strongly recommended

| Field | Why it matters |
| --- | --- |
| `uncertainty_pct` (`rsd %`, `cv %`) | Tightens the predicted-vs-measured envelope on the trust loop. Without it, the ingest assumes a conservative default. |
| `target_benchmark_id` | If you are augmenting an existing benchmark, set this to the matching `benchmark_id`. Otherwise the ingest creates a new one. |
