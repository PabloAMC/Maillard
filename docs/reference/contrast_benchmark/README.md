# Contrast benchmark: the frozen DEVELOPMENT / TEST split (v1, 2026-11-01)

**Status: the foundation only.** Nothing has been scored. No engine was run to build any of this.

## 1. Purpose

The contrast benchmark will score the engine on **within-study contrasts**: the sign and log-ratio between paired arms of one study, and the route shares from labelling. A ratio inside one study cancels that study's calibration and unreported conditions.

The plan from 2026-11-01 (docs/guides/STATUS_AND_PLAN.md) sets rule (i): **no circular scoring**. The cells that diagnosed an engine error cannot also be its test. This folder fixes, before any score exists:
- which laboratories are independent;
- which atlas cells and contrasts are DEVELOPMENT and which are held-out TEST;
- how much held-out power each open engine error E1 - E20 has.

| file | what it is |
|---|---|
| `lab_ids.yml` | every distinct atlas `lab` string mapped to a canonical laboratory id; the merge rules; collaborations; duplicate data sets |
| `split.yml` | GENERATED. The rule, the DEVELOPMENT / TEST lists (labs, cells, contrasts), the counts, the per-error availability, the provenance, and `sha256` |
| `../../../scripts/generators/build_contrast_split.py` | rebuilds `split.yml` from its inputs and asserts the frozen hash |

## 2. The rule

- **DEVELOPMENT** = every atlas cell cited in [ENGINE_VS_LITERATURE.md](../pathway_atlas/ENGINE_VS_LITERATURE.md) secs. 2 - 5, or in either of the two preregs written from it:
  - [kinetic_core_pyrazine_strecker_step_ph_prereg.md](../../../results/validation/kinetic_core_pyrazine_strecker_step_ph_prereg.md) (E4);
  - [kinetic_core_mft_route_barriers_prereg.md](../../../results/validation/kinetic_core_mft_route_barriers_prereg.md) (E6 + E11);
  
  PLUS every cell of any laboratory that has a cited cell.
- **TEST** = the cells of every laboratory with no cited cell.

A cell is "cited" in three ways:
- by file and index (`F25#9`; lists and ranges such as `F17#0 - 7` and `F03#2 Nestle, #7 Baynes` are expanded);
- by naming its study (author-year, resolved through the atlas dossier stems: "Zhai 2019", "Cerny 2003");
- by naming its laboratory (`text_aliases` in `lab_ids.yml`: "Glomb", "Leahy", "DFA").

How the generator reads the rule:
- **"Secs. 2 - 5"** is every block numbered 2.x - 5.x, wherever it sits in the file. That includes sec. 5.1, which was appended after sec. 6. It also includes Appendix A, the verification probe behind sec. 2's V*p* rows. Sec. 1 (the literature tables) and sec. 6 (the recount) are excluded.
- **A cell carries every lab id of its `lab` string:** the lead lab, plus any collaborating lab and any original-data lab (from `lab_ids.yml` `collaborations` / `data_copies`). The cell is DEVELOPMENT if ANY of its labs is a development lab. For the independent-lab count it belongs to its lead lab only.
- **The preregs are read literally.** Every cell they cite is DEVELOPMENT, including the cells they list as held-out tests. Two consequences (sec. 6):
  - Monforte 2021 (Porto) becomes DEVELOPMENT. The E4 prereg names it, as a non-binding transfer diagnostic.
  - The E6 + E11 prereg's held-out cells H1 - H4 (Liu 2023b and Du 2023 at CAAS; Wang 2012 and Zhao 2019b at BTBU) are NOT held out under this split. Secs. 2 - 5 cite both laboratories (F23#7 BTBU, F23#8 CAAS, for E1), so by the lab rule they are DEVELOPMENT for every error.

**Candidate contrasts** (stable id `F25#9@meynier1995` = family, 0-based cell index, first dossier stem):
- `paired_arms`: evidence `within_study_factor`, with a numeric effect in `contrast` or `quantity` (two or more arm values, a fold or a ratio, after the condition values are stripped), and no comparison against another study.
- `route_share`: evidence `labelling`, with a numeric route share (a % or an n : m split).

There are 16 excluded candidates; each is listed in `split.yml` with its reason.

**Laboratories (lab_ids.yml).** The 200 distinct `lab` strings map to **68 canonical laboratories**. Merges follow five written rules:
1. spelling variants of one group;
2. interlinked co-author groups at one institute;
3. a PI who moved keeps the lab (Cerny Nestle -> Firmenich; Hofmann DFA -> TUM). A junior author who moved does not (Cerny 1994 at DFA; Kerler; Meynier 2004);
4. "X with Y" papers belong to the lead lab and also carry the partner;
5. duplicate data.

Every non-obvious merge carries a `merge_note`. The merges that decide something:
- the Ames papers join the Reading lab, because Elmore co-authors Hincelin 1992;
- Jiangnan (Zhang / Cui) + Rutgers (Ho) are one "Ho network";
- BTBU is one lab, and so is CAAS;
- Wageningen FQD (van Boekel, Fogliano, Hettinga) is one lab;
- Fife's "USC" (Southern California) is NOT Baynes's "USC" (South Carolina).

**Duplicate data sets** come from the register's `duplicate` entries:
- **Cross-lab: one.** CF-018: Cammerer 1994 (Kroh, TU Berlin) reprints Wedzicha 1992 (Leeds) Table 5, so F31#6 also carries Leeds.
- **Same-lab: ten entries** (CF-005 / 006 / 007 / 008 / 019 / 040 / 066 / 070), recorded for the reader. They cannot inflate a laboratory count.

## 3. Frozen

The split's membership (`split.yml` -> `split:`: the labs, the cells and the contrast ids with their partition) is hashed:

    sha256 = 022d0cab5a53ececb27a1fc6837ff275fbab177cf2161ae1f8a4ed0341f273a7
    (sha256 of json.dumps(split, sort_keys=True, separators=(",", ":"), ensure_ascii=True))

`build_contrast_split.py` holds this hash and the atlas extent at the freeze (`FROZEN_EXTENT`). Run it with no flags to check. It **fails** when:
- the membership rebuilt from the inputs no longer hashes to it;
- `split.yml` was hand-edited;
- the cited text changed: ENGINE_VS_LITERATURE secs. 2 - 5 or either prereg, each recorded by sha256. The lab-name aliases were reviewed against that exact text;
- an atlas cell below the frozen extent was deleted, or a lab string has no canonical id.

`--write` refuses to write a split whose membership moved.

**Changing the split needs a NEW, dated, reasoned version, never an edit of this one:**
- a new `VERSION` / `FROZEN_SHA256` in the script;
- the old file kept as `split_v1.yml`;
- an entry in sec. 7 below saying what changed and why.

Atlas cells appended after the freeze (rule (i)'s "newly keyed rows", track P1) are outside v1. `split.yml` -> `annotations.post_freeze_cells` lists them, and they enter only through a new version.

    ./scripts/docker_maillard.sh run "python scripts/generators/build_contrast_split.py"           # check
    ./scripts/docker_maillard.sh run "python scripts/generators/build_contrast_split.py --write"   # regenerate (same membership only)

## 4. Counts

**Laboratories.**
- Before normalisation: 200 lab strings.
- After normalisation: 68 laboratories.
- 35 DEVELOPMENT labs and 33 TEST labs.

**Cells and contrasts.**

| | total | DEVELOPMENT | TEST |
|---|---|---|---|
| cells | 562 | 486 | 76 |
| candidate contrasts | 297 | 265 | 32 |
| of which paired arms | 213 | 187 | 26 |
| of which route shares | 84 | 78 | 6 |

The 32 TEST contrasts come from 18 laboratories.

**By odorant class** (the family's class, `ODORANT_CLASS` in the script):

| odorant class | cells DEV / TEST | contrasts DEV / TEST | TEST labs with a contrast |
|---|---|---|---|
| thiol_h2s_sulfur (F22 - F28) | 150 / 7 | 107 / 5 | 4 |
| strecker_aldehyde_acid (F16 - F18) | 43 / 3 | 22 / 1 | 1 |
| pyrazine_n_heterocycle (F19 - F21) | 45 / 7 | 29 / 2 | 1 |
| furans_furanones (F13) | 13 / 3 | 9 / 2 | 1 |
| acrylamide (F29) | 21 / 5 | 13 / 1 | 1 |
| lipid_crossover (F30) | 17 / 4 | 14 / 1 | 1 |
| precursor_non_odorant (F01 - F12, F14, F15, F31, F32) | 197 / 47 | 71 / 20 | 11 |

Most of the held-out set is precursor chemistry, from glycation powders, acid / base sugar chemistry and melanoidins. The odorant lanes keep **5 + 1 + 2 + 2 + 1 + 1 = 12 TEST contrasts**. The literature the engine's errors were diagnosed on is the same literature that measures odorants.

## 5. Power before score: held-out availability per open engine error

What the columns mean:
- **Scope** of each error (`ERROR_SCOPES` in the script): the families, and where it matters the varied factor, the contrast kind or a text marker, that could test it.
- **Factor tags:** the 32 TEST contrasts' tags were reviewed by hand (`FACTOR_REVIEW`). The DEVELOPMENT tags are a heuristic, so read the DEV columns as context only.
- **Kind:** per rule (ii), the route-share errors E6 and E11 count only labelling (route-share) cells.
- **UNTESTABLE** means fewer than 2 independent TEST labs. Its tests must come from new keying (P1: paired-arm rows from laboratories outside DEVELOPMENT) or from designed data (ENGINE_VS_LITERATURE sec. 3.1).

| error | DEV contrasts / labs | TEST contrasts / labs | TEST labs | TEST labs not named by any other prereg | verdict |
|---|---|---|---|---|---|
| E1 hexose MFT, no intact-C5 route | 8 / 3 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E2 enolisation split frozen in pH | 4 / 2 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E3 no Strecker-aldehyde lane | 9 / 5 | 1 / 1 | TUM Becker | 0 | **UNTESTABLE** |
| E4 pyrazine pH law on the Strecker step | 6 / 4 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E5 no buffer catalysis | 43 / 19 | 1 / 1 | Unilever | 1 | **UNTESTABLE** |
| E6 MFT via free NF (route share) | 5 / 1 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E7 lipid inert | 23 / 7 | 1 / 1 | INRA Nantes | 0 | **UNTESTABLE** |
| E8 FFT supply pH law | 8 / 3 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E9 acrylamide | 13 / 3 | 1 / 1 | Jiangnan (Chen / Zeng) | 0 | **UNTESTABLE** |
| E10 TTCA off / release | 8 / 3 | 1 / 1 | Hokkaido | 0 | **UNTESTABLE** |
| E11 one sulfur barrier (route mix vs T) | 34 / 3 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E12 thiol losses | 15 / 9 | 2 / 2 | Ogawa (Kumazawa 2003), Stack 2018 | 0 | testable (but see below) |
| E13 CF-077 cross pair | 14 / 5 | 2 / 1 | UC Davis (Shibamoto 1977) | 1 | **UNTESTABLE** |
| E14 aminoketone condensation | 29 / 6 | 2 / 1 | UC Davis | 1 | **UNTESTABLE** |
| E15 no Strecker acid / O2 route | 4 / 1 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E16 NF decay base / O2 | 2 / 2 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E17 pentose ARP returns no amine | 3 / 2 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E18 3-MP only from thiamine | 25 / 7 | 1 / 1 | INRA Clermont | 0 | **UNTESTABLE** |
| E19 Cys-Strecker H2S off | 9 / 4 | 2 / 1 | Hokkaido | 0 | **UNTESTABLE** |
| E20 no lipid Strecker / pentylpyridine | 17 / 4 | 0 / 0 | - | 0 | **UNTESTABLE** |

**Reading it.**
- **19 of 20 errors are UNTESTABLE on held-out atlas data.** Ten have no TEST contrast at all: E1, E2, E4, E6, E8, E11, E15 - E17 and E20.
- **E12 is the only error with two TEST labs, and neither is clean:**
  - Kumazawa 2003's FFT pH grid is a fit row of the shipped thiol-loss constant (B9 / B11 preregs);
  - Stack 2018's K(T) is a constant the engine uses as fixed (B17).
  
  The column "not named by any other prereg" counts TEST labs whose papers no other `results/validation/*prereg*.md` names. That column is 0 for E12. Rule (i) does not exclude these labs: it counts cells cited by the two preregs. But a scorer must report TEST results with and without them (`split.yml` -> `annotations`).
- **So with a clean held-out criterion, no error has two independent held-out labs today.** Every mechanism change in the plan's steps 3 - 6 (E4, E6 + E11, E10, E1, E2, E5) can be judged on held-out data only after P1 keying, or after designed data, finds rows.

**Sensitivity to the judgement calls** (all checked; each could only be adopted in a new version):

| call (v1 choice) | the other choice | what moves |
|---|---|---|
| Ames papers merged into Reading (Elmore co-authors Hincelin 1992) | Ames a separate lab | E3 -> 2 TEST labs (Becker + Ames); E18 -> 2 (INRA Clermont + Ames, Hincelin 1992); E20 -> 1 |
| E4 prereg's held-out Monforte 2021 read literally as cited (DEV) | role-aware (held-out citations exempt) | Porto TEST: E3 +1 lab, E4 -> 1 lab (still untestable) |
| "Kallen-shaped" K(pH) (sec. 3, E10) counts as citing Kallen | not a citation | Kallen 1971b F04#2 TEST: E10 -> 2 labs |
| "Namiki C2" (F09-a) counts as citing Nagoya | not a citation | Nagoya TEST; no error scope reaches F09 |
| Meynier 2004 kept at INRA Nantes (junior author) | merged with Reading | E7 -> 0 |

Even under the most generous choices, at most E3, E10, E12 and E18 reach two TEST labs. Each of them has at least one prereg-named lab.

## 6. Caveats

- **The numbers are not yet re-read.** Rule (i) requires the held-out numbers to be re-read in the PDF by two readers before scoring. None of the 32 TEST contrasts had been when v1 was frozen; the 2026-11-01 re-read is in sec. 9. The atlas's per-cell wrong rate is about 9 % (95 % CI 4 - 18 %).
- **Heuristic parts:**
  - the numeric-effect test;
  - the DEVELOPMENT factor tags;
  - the alias detection, reviewed once against the cited text and frozen with its hash.
  
  Each is deterministic, and the generator lists what it excluded and why.
- **Error scopes were declared from secs. 3 and 5.** One was narrowed after the first build: E16's text marker went from "norfuraneol / NF / HDMF / furaneol" to norfuraneol only. E16 is norfuraneol's decay, and furaneol is another compound. The change removed one TEST contrast, Haleva 1999 F13#13, a furaneol pH contrast.
- **Two premortem amendments conflict with the global lab rule:**
  - the E6 + E11 prereg's "held-out tests come from the labelling cells of CAAS / BTBU" (both labs are DEVELOPMENT here);
  - the E4 prereg's own conclusion, which agrees: zero admissible held-out labs.
  
  The coordinator decided: v2 (sec. 8) applies rule (i) per change, beside the frozen v1.

## 7. Versions

- **v1, 2026-11-01:** first freeze, sha256 `022d0cab5a53ececb27a1fc6837ff275fbab177cf2161ae1f8a4ed0341f273a7`.
  - lab_ids_v1;
  - atlas extent as of the 2026-11-01 fix-up (562 cells);
  - ENGINE_VS_LITERATURE including sec. 5.1 (E4 correction);
  - the two preregs as of 2026-11-01 ~03:00.
- **v2 (per change), 2026-11-01:** `split_v2_per_change.yml`, sha256 `dbe32a2b3af83f3bbb603efb391dd8477999263a6b24769ea946fcc0e0fe2e12`, declared before any scoring (sec. 8). It sits beside v1; v1 is unchanged and still checks OK.

## 8. v2: per-change hold-out; why

**Why.** v1 applies rule (i) globally. A lab cited for ANY error is development for EVERY error, so BTBU and CAAS, cited for E1 only, could not test E6. That is stricter than standard hold-out practice. The usual rule is that the data used to diagnose or tune a change cannot test THAT change, and anything else can. The coordinator asked, on 2026-11-01, for a per-change split declared before any score. v1 stays frozen as the conservative global reading. A score should name which split it uses.

**The rule (per error E1 - E20):**
- **DEV(E)** = the cells cited AS EVIDENCE FOR E in ENGINE_VS_LITERATURE secs. 2 - 5, plus E's own prereg cells, plus every cell of those cells' and named laboratories. "Evidence for E" is read from E's own text only:
  - E's sec. 2.1 mismatch rows (`MISMATCH_ERRORS` in the script writes out sec. 3's first column);
  - its rows in secs. 3, 3.1, 5 and 5.1;
  - its recommended-order item in sec. 5 (so the cells its fix was tuned on count, e.g. E4's band on Leahy);
  - its Appendix-A probes;
  - its sec. 2.0 correction (item 1 is E8, item 2 is E10);
  - one unnamed citation, resolved by hand: E17's "hexose ARP 71 - 96 % (3 labs)" = F03#0, #3, #11.
  
  The E4 prereg feeds E4, and the E6 + E11 prereg feeds E6 and E11.
- **TEST(E)** = the candidate contrasts in E's scope from laboratories not in DEV(E), minus any contrast whose paper is a **calibration input** of the shipped engine.

**Choices (all declared before any v2 count was read, except where marked):**
- **Conservative choices kept, as written:** the Ames papers are merged into Reading, and "Kallen-shaped" counts as citing Kallen (so Kallen is DEV(E10)).
- **A prereg's declared held-out tests are not its development.** A target cited in the prereg's sec. 4.1 "Held-out", and not among its sec. 4.2 development checks, does not enter DEV(E). This holds wherever else the prereg mentions it (its summary, its baselines table). Without this, every prereg would turn its own tests into development. The literal reading is reported as an alternative.
- **Calibration inputs** are paper level, which is conservative: one shipped row excludes the whole paper. There are 50 papers:
  1. measurement-ledger entries `in_model_fit` without a not-shipped flag;
  2. the shipped sulfur fit B60's 52 objective rows, mapped by row-id prefix to their papers (the ledger predates part of B60);
  3. measured constants the engine ships fixed (`SRC_FIXED_CONSTANTS`): Chan & Reineccius 1994, Dwivedi & Arnold 1971 (B50) and Stack 2018 (K(T));
  4. a register data copy of any of these (e.g. Schieberle 1998 Table IV reprints Hofmann 1998 T1).
  
  Every other atlas paper the engine's source cites was read and recorded as mention-only, with its reason (`SRC_MENTION_ONLY`). Examples: Du 2023 and Zhao 2019b document the HEXOSE_C5_ENTRY switch, which ships OFF; Quan 2020's constants are refused. The generator fails if a cited paper enters a TEST(E) set unreviewed.
- **Scopes:** v1's, with two widenings that match each error's TEST set to its own prereg's search:
  - E4 also covers F17 and F21;
  - E6 covers every F23 route share, not only those that name NF.
- **Factor tags:** v1's 32 reviewed tags, plus 23 more reviewed for v2 (`FACTOR_REVIEW_V2`). For every factor-scoped error (E2, E4, E5, E7, E8, E16), each in-family candidate whose text carries the factor's vocabulary anywhere had to be reviewed before the build ran.

**Result (power only; nothing scored):**

| error | DEV(E) labs | TEST(E) contrasts / labs | TEST(E) labs | contrasts excluded as calibration | verdict |
|---|---|---|---|---|---|
| E1 hexose MFT, no intact C5 | 4 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E2 enolisation split vs pH | 3 | 0 / 0 | - | 2 (Zhai 2021) | **UNTESTABLE** |
| E3 no Strecker-aldehyde lane | 3 | 5 / 3 | Porto, Reading, TUM Becker | 1 | testable |
| E4 pyrazine pH law on the step | 6 | 1 / 1 | Porto (Monforte 2021) | 0 | **UNTESTABLE** |
| E5 no buffer catalysis | 9 | 5 / 3 | Essex (Thornalley), Kallen, Unilever | 3 | testable |
| E6 MFT via free NF | 5 | 4 / 2 | BTBU, CAAS | 0 | testable |
| E7 lipid inert | 3 | 4 / 3 | Ho network, Nestle (patent), Seville | 1 | testable |
| E8 FFT supply pH law | 3 | 0 / 0 | - | 4 (Cerny 2007) | **UNTESTABLE** |
| E9 acrylamide | 4 | 4 / 2 | Jiangnan (Chen / Zeng), Reading | 0 | testable |
| E10 TTCA | 3 | 1 / 1 | Hokkaido | 0 | **UNTESTABLE** |
| E11 one sulfur barrier | 5 | 9 / 2 | BTBU, CAAS | 0 | testable |
| E12 thiol losses | 3 | 4 / 4 | BTBU, CAAS, Huazhong, Stockholm | 5 (Kumazawa, Stack, ...) | testable |
| E13 CF-077 cross pair | 2 | 8 / 4 | Ghent, Nestle, Reading, UC Davis | 0 | testable |
| E14 aminoketone condensation | 5 | 7 / 3 | Reading, Seville, UC Davis | 5 (Leahy) | testable |
| E15 Strecker acid / O2 | 1 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E16 NF decay | 2 | 0 / 0 | - | 0 | **UNTESTABLE** |
| E17 pentose ARP amine return | 3 | 2 / 1 | Baynes | 0 | **UNTESTABLE** |
| E18 3-MP only from thiamine | 2 | 9 / 5 | BTBU, DFA, Ho network, INRA Clermont, Reading | 2 | testable |
| E19 Cys-Strecker H2S | 1 | 6 / 3 | Copenhagen, Ho network, Hokkaido | 4 (Zheng & Ho 1994, Cerny 2007) | testable |
| E20 lipid Strecker / pentylpyridine | 2 | 8 / 3 | CAAS, Ho network, Reading | 0 | testable |

**12 errors become testable** (>= 2 TEST(E) labs): E3, E5, E6, E7, E9, E11, E12, E13, E14, E18, E19 and E20. **8 stay UNTESTABLE:** E1, E2, E4, E8, E10, E15, E16 and E17. Their tests must come from P1 keying or designed data.

Of the plan's next changes, E6 + E11 and E5 have held-out tests. E4, E1, E2 and E10 do not. E4's single TEST lab, Porto, is the one its own prereg declared non-binding.

**Alternatives** (`split_v2_per_change.yml` -> `alternatives`; TEST labs per error):

| variant | what changes |
|---|---|
| Ames a separate lab | E3 3 -> 4, E18 5 -> 6, E20 3 -> 4; no verdict changes |
| "Kallen-shaped" not a citation of Kallen | E10 1 -> 2: **E10 becomes testable** (Kallen 1971b F04#2 + Hokkaido) |
| prereg read literally (its held-out list is development) | **E6 and E11 drop to 0** (BTBU and CAAS become development); E4 1 -> 0 |
| calibration = ledger + B60 only (no fixed constants) | E3 3 -> 4 (Chan 1994b), E5 3 -> 4 (Arnold 1971), E12 4 -> 5 (Stack 2018); no verdict changes |
| calibration = every paper src/kinetic_core cites | **E9 drops to 1 lab** (Quan 2020); E6 / E11 keep 2 labs on fewer contrasts (Du 2023 and Zhao 2019b go); E20 3 -> 2 |

**Caveats specific to v2:**
- **The per-change scope is not power to discriminate.** The E6 + E11 prereg already expects ~0 power from CAAS's Liu 2023b (no rival touches the fragment route), and it scores Zhao 2019b / Du 2023 under E1, not E6. Only BTBU's Wang 2012 has leverage there. Availability is a necessary condition, not a sufficient one.
- **TEST(E) labs overlap DEV(E') for other errors.** A single regeneration that ships several changes needs the union rule (v1) for the joint claim.
- **Sec. 5's "8+ labs" for E5 names no laboratory** and cannot be resolved. Any of them may sit in TEST(E5) (Essex, Kallen, Unilever).
- **The two scope widenings (E4, E6) were written before the v2 counts were read.** The calibration curation (`SRC_MENTION_ONLY`, the data-copy rule) was written AFTER a first diagnostic pass showed which papers it affects. The `calibration_broad` alternative bounds that choice.

**Two-reader PDF re-read still needed:**
- the 65 contrasts that sit in at least one TEST(E) set;
- v1's 32 TEST contrasts. Nine of them are among the 65, and the other 23 matter only if v1 is scored.

That makes **88 distinct contrasts in all, of which 56 are new to v2** (`split_v2_per_change.yml` -> `reread_list`). Re-read 2026-11-01: sec. 9.

    ./scripts/docker_maillard.sh run "python scripts/generators/build_contrast_split.py --v2"           # check (~2 min)
    ./scripts/docker_maillard.sh run "python scripts/generators/build_contrast_split.py --v2 --write"   # same membership only

## 9. reread_status (2026-11-01)

The 88 `reread_list` contrasts were re-read against the PDFs by two independent readers (A read all 88, B read 58: the priority group, errors E5 / E6 / E7 / E11 / E13, plus 30 others), then adjudicated. Rules: both readers agree -> their grade; they disagree, or a single reader grades MINOR / WRONG -> the adjudicator opens the PDF and decides; single reader MATCH -> accepted, marked single-reader. Corrections were applied to the atlas cells, each recorded in the cell's `note` ("re-read 2026-11-01: ..."); `quote` fields were left alone because they quote the dossier, not the PDF. The per-entry record (grades, readers, basis, action) is in the session's `reread_adjudication.yml` (scratchpad).

```yaml
reread_status:
  date: 2026-11-01
  entries: 88
  final_grades: {MATCH: 70, MINOR: 15, WRONG: 3}
  verification: {two_reader: 58, reader_plus_adjudicator: 6, single_reader_match: 24}
  v2_test_contrasts: {total: 65, two_reader: 51, reader_plus_adjudicator: 3, single_reader_match: 11}   # 54 / 65 verified by two people
  v1_test_contrasts: {total: 32, two_reader: 14, reader_plus_adjudicator: 4, single_reader_match: 14}  # 18 / 32
  wrong: [F05#23@thornalley1999, F13#13@haleva1999, F14#12@haleva1999]
  single_reader_v2_test: [F03#8, F16#11, F16#14, F16#15, F16#16, F16#17, F18#16, F21#4, F21#6, F21#8, F21#9]
```

- **WRONG (3), all corrected:** F05#23 said Thornalley 1999 printed no HSA fructosamine, but Fig. 2d plots it, so the cell is now `operating`, not `unattributed`. F13#13 and F14#12 (Haleva 1999) were tagged `phosphate`; the buffer is 0.1 M citrate-phosphate, now `other_buffer`.
- **MINOR (15), all corrected:** mostly a wrong range or condition. Examples: Kallen's ceiling is x3.1-4.4, not x3.5-5; Stack's kinetics buffer is 0.050 M; Low 2006's ND means not determined; Du 2023's "92 % no glucose carbon" is ~84 % in a 1:1 CAMOLA.
- **Liu 2023b (F23#6, F24#6; E6 / E11):** both readers confirm the printed 41 / 3 / 2 / 1 / 12 / 40. The adjudication (CF-086) finds that the Methods' [M-1]+/[M]+ adjustment did not reach the labelled MFT ion, so the fragment share is ~6 %, not 19 %; P 0.9.
- **Still single-reader:** 11 of the 65 v2 TEST contrasts, including all five of E3's (F16, Strecker aldehydes). A second read of these 11 would finish v2. v1's 14 single-reader contrasts matter only if v1 is scored.
- **Membership is unchanged.** No correction moved a cell into or out of the candidate set; `build_contrast_split.py --v2` still checks OK. The v1 check's `contrasts` block now differs from a rebuild only in the `status` of F05#23 and F24#6. That check was already failing before the re-read (cited-text drift, `development_lab_reasons`), so `--write` waits for that review.

## 10. Addendum: newly keyed TEST rows (2026-11-01)

Three held-out sources were keyed by two independent readers and adjudicated (P1 keying, rule (i)'s "newly keyed rows"). They are **TEST data only, never a fit input**: `fit_rows.load` refuses all three tables (units or species outside its vocabulary; laboratories not in `data/keys/laboratories.yml`).

**They sit outside the frozen splits.** v1 and v2 hash atlas cells and contrast ids, and these are fit-row tables, not atlas cells. `split.yml` and `split_v2_per_change.yml` are unchanged. A score that uses these rows must say so and report them beside the v2 TEST set. To count them as v2 contrasts takes a new, dated version (sec. 3).

| table (`data/lit/fit_rows/`) | error | lab | rows | the contrast it carries |
|---|---|---|---|---|
| `heldout_zhang2021_table3.yml` | E7 lipid inert | China Meat Research Center, Beijing (Zang M.); not in `lab_ids_v1` (proposed id `china_meat_research_center_beijing`); no author shared with CAAS, BTBU or any DEV(E7) lab | 68 (+ 28 not-detected cells listed per pot) | MFT, FFT and ten sulfur volatiles (thiophenes, a thiophenone, 2-methylpentane-3-thiol), with vs without soybean lecithin / cephalin / both, Cys + Xyl 80 + 80 mM, pH 5.5, 120 C / 60 min, with or without reheating |
| `heldout_larralde1995.yml` | E5 buffer catalysis | UCSD (Miller); not in `lab_ids_v1`, and not the University of San Diego (De Haan) | 26 | bicarbonate / phosphate = 0.52 at equal pD 7.4; the phosphate rate law k0 + kP [P] at 100 C; Fig. 2 (Arrhenius) and Fig. 1 (pD profile) rate constants |
| `heldout_lee2024_fig2.yml` | E13 CF-077 cross pair | SayFood, Paris-Saclay (`sayfood_paris_saclay`, already a TEST lab in v1 via Lee 2022; not in DEV(E13)) | 41 | PZ, MPZ, 2,5-DMP and 2,6-DMP in the oven vapour of a glucose + Leu cake at 140 / 170 / 200 C: MPZ^2 / (PZ x DMP) vs time and temperature |

If a new version admits these labs, TEST(E) labs grow by one each:
- E5: 3 -> 4;
- E7: 3 -> 4;
- E13: 4 -> 5.

**Verification.** Each table has a first reading (the adjudicated rows), a second reading under `second/`, two `source_verification` records, and a reconciliation note `results/validation/source_reading_reconciliation_heldout_*.md`. All three reconcile with 0 disagreements.

**Caveats a scorer must carry.**
- **Zhang 2021.** The printed D doses conflict: D1 is 0.24 % + 0.24 %, D2 is 0.07 % + 0.07 %, yet the text calls D2 the reheated D1. So D1 -> D2 is not a clean reheating contrast. The ug/kg values are semi-quantitative (one internal standard, response factor 1).
- **Larralde 1995.**
  - The rate law's printed "s-1" is read as min-1 (P ~0.9). It is consistent with Fig. 2 only in min-1.
  - The bicarbonate temperature (100 C) is inferred.
  - The system is amine-free ribose loss. It tests buffer-species catalysis of sugar degradation, not of the Amadori step.
- **Lee 2024.**
  - These are vapour concentrations in an open, non-isothermal oven. 39 of the 80 points were dropped: at the reading floor, or 140 C points plotted at a time the Methods did not sample.
  - The text says MP is the lowest pyrazine, but the figure shows it the highest (unresolved). So any E13 ratio must also be reported with MPZ / 10.
