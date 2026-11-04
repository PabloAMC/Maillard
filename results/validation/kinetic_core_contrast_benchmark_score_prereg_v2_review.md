# Review: kinetic_core_contrast_benchmark_score_prereg.md (v2)

*Independent reviewer (a fresh agent without the author's context), 2026-11-03. Read-only. Exposure: no `engine.predict`
was run; the two XU07-PH-01 / 02 rows of `results/validation/core_directional_scores.md:153-154` (the Xu MFT / FFT pairs
v2 already tags `dev_exposed`) were seen while grepping claim ids; no Meynier, Grosch or clean-pair engine value was
read. Counts below are from the YAML only, via `docker_maillard.sh`. The v2 document is kept as
`kinetic_core_contrast_benchmark_score_prereg_v2_superseded.md`.*

## v1 items: resolved?

| v1 item | status | where |
|---|---|---|
| B1 exposure | **Resolved.** The `dev_exposed` tag exists, D1 / D2 are decided on clean pairs, and the DEV labs are reported. NF is left clean, unlike the v1 fix. That is correct (the b2 claims are MFT / FFT / FUR from Table 3), but v2 does not say why. | v2:17-21, 51-52, 60-61 |
| B2 names | **Partly resolved.** A map exists, but it contradicts the code for GSH and LYS, omits MPZ / DMP / MTAL / DMDS, and does not name the aliases for MFTD / ACTZ. See new B4. | v2:66-70 |
| B2 buffer / vessel / ph_final / digest / sensitivity / Met | **Resolved.** | v2:65, 71-76 |
| B3 one null | **Resolved.** | v2:82-84 |
| S1 power | **Partly resolved.** The two-of-three-tables rule is in. The power line cannot yet be computed as specified (S1 below). | v2:94-95, 101-102 |
| S2 censored | **Partly resolved.** Censored pairs are in, but the inclusion rule, `tr` cells and the "n.d., no limit" cells are not defined (B6). | v2:53-55 |
| S3 disclosure | **Resolved.** | v2:14-16 |
| S4 pair rule | **Partly resolved.** It works on the three scored tables. It breaks on the other 12 (B5). The swap exclusion is fine. The F26#11 de-duplication is wrong (S3). Engine zeros are "kept" but given no value in S or E (B6). | v2:48-50, 56-58 |
| minor precedence / hash / per-ladder medians | **Resolved.** | v2:32-40, 43-45, 92 |

## Blocking

**B4. The name map contradicts the code, is incomplete, and the per-call target set is not frozen. Each of these changes the item set or refuses whole pots.**

- **GSH is mappable.** `"gsh"` and `"glutathione"` map to `GSH` on the sulfur lane (`engine.py:244-246`), and the Grosch header itself lists GSH among the "Engine keys" (`heldout_grosch1992_tab2.yml:23`). Class 4 ("no engine name", v2:36) is false for it.
  - Leaving it out removes 3 of Grosch's 8 one-factor pot pairs: no. 3-7, 6-8 and 7-8.
  - That drops 3 MFT pairs, 1 detected FFT pair (6-8) and 3 censored pairs (FFT 7-8, MFTD 6-8 and 7-8).
  - The effect is to shrink Grosch, one of the three tables in D2's two-of-three rule, from about 7-11 S items to 3-5.
- **LYS is mappable.** `"lysine"` maps to `Lys`, and the envelope then refuses with "LANE CONFLICT" (`engine.py:553, 975-983`). That is class 6 verbatim, not class 4. Q1 is "a property of the code" (v2:105), so the class must be the code's.
- **MFTD / ACTZ.** "The engine's own target names" (v2:68) are not alias keys. `"mftd"` and `"actz"` are absent from `TARGET_ALIASES`; the aliases are `"mft dimer"` / `"bis(2-methyl-3-furyl) disulfide"` and `"2-acetylthiazole"` (`engine.py:321-323, 351`).
- **Missing species with engine targets.** MPZ, DMP, MTAL and DMDS are rows in Meynier (MPZ / DMP in Tables 1-2, MTAL / DMDS in Table 4), and MPZ is a row in all three Xu pots. All four map to trunk-only targets.
  - Asking for any of them on a sulfur-lane pot (any pot with PENT) adds a whole-request refusal: "run on the trunk lane only" (`engine.py:1303-1311`).
  - Any unmapped name (`V_*`) does the same: "UNREPRESENTED TARGETS" (`engine.py:1189-1193`).
  - So "one engine.predict call per pot" (v2:65) with the pot's full species list refuses every Xu pot, Meynier Gly pH 5.5-6.5 and every Met pot, and FUR / NF / ACTZ / MFTD / MFT go with them.
  - These species are also the E4 consumers the Meynier table was keyed for (README sec. 11). Q1 needs them classified (class 6), and the map as written fails the run on them (v2:69-70).
- **Fix:**
  - Freeze the map by rule: a YAML key K maps to an alias whose value is K (case-sensitive); no such alias means unmapped. This gives GSH -> glutathione, MET -> methionine, THI -> thiamine, PENT -> ribose, MPZ -> 2-methylpyrazine, DMP -> 2,5-dimethylpyrazine, MTAL -> methional, DMDS -> dimethyl disulfide, MFTD -> "mft dimer", ACTZ -> "2-acetylthiazole". It leaves LYS and TPP unmapped. Print the resolved map in the artifact.
  - If GSH is to stay out, put it in an explicitly named exclusion class with its reason, not in class 4.
  - Classify LYS + PENT pots by the envelope (class 6).
  - Make it one `predict` / `declare_envelope` call per (pot, species), so that one target's refusal cannot refuse another's.

**B5. The pair rule (v2:48-50) is used for Q1 over all 15 tables (v2:42), where it mislabels most of the items and loses the contrasts class 5 exists for.**

The rule treats "a factor" as an `initial_mmol_L` key or `label_ph` only. Counted literally on the YAML:
- **811 one-factor pot pairs (2,073 pair x species with rows in both arms).** 709 of the pot pairs also differ in a field the rule does not see:

  | table | flagged pairs | of | the unseen difference |
  |---|---|---|---|
  | `heldout_zamora2015_tab1` | 552 | 576 | `initial_mmol_L` is empty, so every pH 3 / pH 6 pair counts, whatever the lipid (`amounts_umol`) and the N2 / air atmosphere |
  | `heldout_larralde1995` | 129 | 144 | `t_c` |
  | `heldout_thornalley1999` | 28 | 32 | the +BOC-lysine pairs also differ in `phosphate_mM` |

- **Zero items** from the real one-factor contrasts that sit outside `initial_mmol_L` / `label_ph`: Wang 2022 lipid, Zhang 2021 phospholipid and reheating, Zhao 2019b lipid, the Thornalley phosphate series, and the Liu 2023b time and Lee 2024 temperature pairs. The Wang and Zhang pairs are Cys + PENT sulfur-lane pots with MFT / FFT rows: exactly the `factor_not_an_input` items (lipid, reheating, buffer concentration; v2:37-38) that Q1 is meant to count.
- **Unit unstated.** v2 does not say whether a Q1 item is a pot pair or a pot pair x species.
- **Fix:**
  - Define a factor as any differing field among `initial_mmol_L` keys, `amounts_umol` keys, `label_ph`, `t_c`, `minutes`, `buffer`, `phosphate_mM`, `lipid`, `phospholipid_pct`, `aldehyde`, `post_heat`, `atmosphere`, `protein_mg_per_mL`, `matrix`. Keep "exactly one".
  - State that a Q1 item is (pot pair, species with a row or nd / tr in either arm).
  - Keep sec. 2B's narrower rule for Q2 only, where it is exact for the three tables: none of their pairs differs in an unseen field.

**B6. Three scoring choices are left open and can be made after the engine numbers are seen.**

- **Engine zeros and refusals (v2:56-57).** They "stay in the item count", but S (v2:88) and E (v2:90) give them no value.
  - With p = log10(0/x), or no p at all, the implementer can count them as misses or drop them.
  - Under B4 this decides up to 21 MPZ / DMP items in Meynier Gly and Xu, plus any FFT arm refused as unreachable (thiamine-only Grosch no. 1).
  - It also decides the "actual ladder sizes" the power line needs before scoring.
  - *Fix:* refused or zero in either arm counts as a sign miss in S and gets error = +inf in E (median-robust). Both arms zero gives p = 0. Report the counts. Classify refusals by `declare_envelope` before integration, so that the ladder sizes are known before the run.
- **p != 0 (v2:89).** p = 0 counts in the denominator as a miss, which is acceptable, but there is no tolerance. With `rtol=1e-8` (`engine.py:2620`), a factor the engine treats as inert returns |p| ~1e-9 with a random sign, and that sign scores as a hit.
  - *Fix:* treat |p| < log10(1.01) as p = 0.
- **Censored pairs (v2:53-55).** "|m| >= log10(value / LOD)" is a bound, not an inclusion rule. Three points are open:
  - whether a `tr` arm counts: Meynier `trace_below_0p1`, Xu `trace_below_0p3`;
  - what to do with "n.d., no limit printed" (Grosch no. 4 MFTD);
  - whether value / LOD must pass the 1.4x resolution threshold.

  The last choice flips the two clean Grosch FFT censored items (1.9 vs < 1.5, ratio 1.27) in a table with only 3 resolved detected items.
  - *Fix:* the bound is the LOD for nd (Meynier 0.01 ug / 10 mg ribose; Xu 0.07 ng/mL; Grosch `upper_limits_ug`) and the tr ceiling for tr. Include a pair iff value / bound > 1.4. Exclude "n.d., no limit", `present_not_quantified`, and pairs with both arms censored.

## Should-fix

**S1. The power line (v2:94-95) has no generative model.**
- "True per-ladder sign agreement 0.8" can mean iid Bernoulli per pair or a draw per ladder. Within a pH ladder of k pots, the C(k,2) signs are constrained by the orderings, so iid per pair overstates the effective n.
- The line must also apply D2's two-of-three-tables clause and its own resampling counts and seed.
- *Fix:* specify the model, e.g. per ladder an engine ordering that agrees with each measured adjacent step with probability 0.8, with the pair signs derived from it. Use 2,000 outer simulations x the same 2,000-resample bootstrap, seed 20261101, with D2 applied in full. Say that D3 gates D1 on D2's power, since D1 gets no power line of its own.

**S2. D1 (v2:99-100).** "Excludes 0" also passes when the engine is significantly worse. The pooled-median clause usually guards against this, but not always.
- *Fix:* the interval must lie above 0.

**S3. De-duplication (v2:58).**
- **F26#11@meynier1995 is Table 3 (cysteine)** (atlas `F26_thiazoles_thiophenes.yml` cell 11: `where: Table 3`). The keyed file excludes Table 3 (`heldout_meynier1995_tab1_2_4.yml:13`), so there is nothing to de-duplicate; dropping the cell would make the atlas family 64, not 65.
- F27#15 covers Tables II, III and VIII; the keyed file is Table II only.
- F27#16's TPP-vs-thiamine contrast has no keyed pair (it is a swap).
- Other overlaps go unhandled in Q1: Thornalley F05#23 / F06#11 / F08#2, Zamora F30#15, Liu 2023b, Zhao 2019b, Du 2023 and Wang 2022. README sec. 11 says these become "scorable" from the keyed tables.
- *Fix:* keep all 65 atlas cells as one Q1 family and the keyed pairs as a second. Report the two separately and do not pool the shares. Tag the keyed pairs that reproduce an atlas cell.

**S4. Expected counts.** Sec. 7 (v2:110) says "~10 clean ladders". On the YAML, with v2 as written (GSH, TPP and LYS out, Met out, exposed pairs out):

| table | E: detected pairs / ladders | S: resolved detected | S: censored (B6 rule) | S ladders |
|---|---|---|---|---|
| Grosch | 4 MFT / 2 (PENT, Cys) | 3 (no. 1-3 is 4.7 vs 5.2, unresolved) | 0 (2 if the 1.27x FFT pairs are admitted) | 2 (4) |
| Meynier | 10 NF-Gly / 1 | 9 (pH 5.0 vs 5.5 is 5.2 vs 5.5, unresolved) | 0 | 1 |
| Xu | 5 (MFTD 1, ACTZ 3, FUR 1) / 3 | 5 (ACTZ pH 7 vs 8 is 1.417x, 0.005 dex above threshold) | 4 (MFTD, FUR at pH 8 nd) | 3 |
| **total** | **19 / 6** | **17** | **4 (6)** | **6 (8)** |

- Meynier rests on one ladder, with 9 of the 21-23 S items.
- Mapping GSH (B4) adds to Grosch 4 detected pairs (MFT 3-7, 6-8, 7-8; FFT 6-8) and 3 censored pairs (FFT 7-8 at 1.35x, so excluded under the B6 rule; MFTD 6-8 and 7-8 at 8x), plus the MFT/GSH and MFTD ladders.
- *Fix:* print this table, in its final form, in the prereg before the run.

**S5. Plan step 2 departure.** `STATUS_AND_PLAN.md:170` asks for the engine and nulls on DEVELOPMENT and TEST separately. v2:116 scores DEVELOPMENT only for the exposed pairs, and v1's pointer to `contrast_scores.json` was dropped.
- *Fix:* state the departure and its reason in sec. 8.

## Minor

- **Config (v2:65).** The engine reads globals, not a `Config` (`config.py:18`). Record `Config.from_globals().digest()` and assert that it equals `Config.shipped().digest()`.
- **Buffer.** Write the call as `BufferSpec(kind="pyrophosphate", pyrophosphate_mol_l=0.2, declared=True, source=...)`. The positional form in the v1 review sets `phosphate_mol_l` and raises (`ph_state.py:500-529`).
- **Hand classes for atlas cells.** Class 6 for atlas cells (v2:39, 43-44) cannot be "the envelope's reason verbatim" without a spec. Have the hand table give classes 1-5, or a spec that the envelope then classifies.
- **Met pairs.** State their Q1 class: `answerable`, flagged that methionine is recorded and not charged.
- **Hash.** Record the commit sha of `ATLAS_COVERAGE`, not only its sha256, so that "committed before the first scoring run" can be checked.
- **Grosch vessel.** Grosch is an open reflux, and `vessel=None` charges the sealed-pot default reservoir (`engine.py:2488-2490`). State this as declared.

**Verdict: REVISE.**
