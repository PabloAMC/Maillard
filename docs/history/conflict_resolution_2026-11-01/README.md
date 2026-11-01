# Conflict-resolution programme (2026-11-01; EX-055)

An archived analysis record, not a generator. The records were written in the session scratchpad on 2026-11-01 and
are COPIED here unchanged. Paths inside them (`scratchpad/...`) point at that scratchpad, and the scripts they name
were not archived. Where a red team weakened a claim, the register and the dossiers carry the red team's probability
and wording, not the cluster's.

**How it ran:**
1. A read-only triage of CF-001 to CF-076 (`conflict_triage.md`): register defects first, then the conflicts grouped
   into five clusters by cause class (operator / method, conditions, transcription / duplicates, chemistry, model).
2. Five clusters, each with a pre-registration written before its tests ran (cluster 1 states its decision rule at
   the top of its record). Agents never wrote to the registers; one writer applied the proposals.
3. An independent red team of the clusters' engine claims (`red_team_clusters_2_3_4.md`). Clusters 1 and 5 had none:
   cluster 1 makes no engine claim, and cluster 5's admission is marked provisional until one is run.

## Per cluster

| cluster | record | entries | outcome |
|---|---|---|---|
| 1. same-lab reprints, text / table slips | `cluster1/resolution.md` | CF-018, 019, 020, 029, 040, 057, 067, 070, 074 | A common duplicate rule (>= 3 cells agreeing to the printed precision). Duplicates confirmed: CF-040 (P 0.99), CF-070 (0.9), CF-018 / 029 / 019 in part. CF-057 not an error (0.97). CF-074 explained (one run, 0.8). CF-067 stays open (explanations narrowed). CF-020 split by paper. |
| 2. the deoxyosone hub | `cluster2/prereg.md`, `cluster2/register_edits.yml` (its verdict record) | CF-053, 060, 063, 071 (+ new entries) | CF-053 and CF-063 explained; CF-060 and CF-071 stay open (signature iii is not thiol-diagnostic, P 0.8). Red team (claim 3): the hub's ratios are not separable, a structural fact (P 0.85), but the T8 test is INCONCLUSIVE, not confirmatory, because its negative controls did not fail. Beck's pentose band is weakened (P 0.35), drift-is-required is weakened (0.5), and no trap evidence survives (0.3). Nothing changes in the engine. |
| 3. in-situ pH and speciation | `cluster3/prereg.md`, `cluster3/results.md` | CF-012, 054, 065, 069 (+ CF-084) | The shipped pH model fails a positive (Hou 2026) and a negative (Meynier 1995) control. Red team (claim 1): "fails both" is supported (P 0.85), but "CBX is the cause" is weakened (0.35): the acid is allocated to the wrong source. CF-054 explained (0.75). CF-065 explained in part (conditions ~1.1 dex), residual open. CF-069 explained for the acetic-acid rows (0.6) but stays open. Red team (claim 2): the fed-H2S rows ignore the headspace partition, as the whole sulfur lane does (bias 0.5-0.85 dex, P 0.75); a declared transfer band now, a vessel operator after A3. |
| 4. observation operators, vessel O2, thiol and Strecker levels | `cluster4/prereg.md`, `cluster4/results.md`, `cluster4/results.json` | CF-014, 015, 016, 031, 068, 073; EX-050 | CF-015 explained for the ratio half (sugar : cysteine ratio, P 0.6). CF-068, 073, 016 and 031 stay open with notes: none reaches P 0.6 (CF-073's donor-generic operator is rejected, chi^2 11-22; CF-031's register wording is corrected). CF-014's cause stays unidentified. Grosch's 4-5x is in-pot loss, not extraction recovery (bolton1993, grosch1992, evidence map 2.11). Red team (claim 4): Grosch's standards really did go in after boiling (P 0.85); the Nestle patent bound is not measurement-grade (P 0.15); the proposed MFT loss law is flexible, not identifiable (P 0.1); P(pot-dependent >= 3x) ~0.5, not 0.65. So EX-050 carries an erratum: k_mft_decay is a declared band, Cys-conditional at most, with no law in the engine. |
| 5. cysteine exchange and the Amadori step | `cluster5/prereg.md`, `cluster5/verdict.md`, `cluster5/*.json` | CF-076, CF-011; new CF-091 | The kallen1971 dossier's k_r was the relaxation, 120-270x too fast (corrected with kallen1971b). EX-054 was re-run on the corrected bands with Zhang 2025b's blank as T8, and the re-run supersedes it. The engine exchange is rejected in two laboratories. In-situ pH + Kallen K(pH) / k_r(pH) is feasible: "explained, provisional (P 0.45)", recorded as a note, and CF-076 stays open. CF-011 is explained in part (phosphate 0.5 vs ELSD 0.3) and stays open. Zhai 2021 Fig. 7a's quantitation-by-difference is registered as CF-091. |

## What is still owed

- A red team of cluster 5's admission, and TTCA's ring-N pKa, before any change to the exchange law.
- The co-heated-standard tracer operator (after A3), before the US 10,076,126 shape or level is scored.
- The author requests and designed experiments named in each register entry's remedy.
