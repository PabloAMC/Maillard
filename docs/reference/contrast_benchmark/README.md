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

- **The numbers are not yet re-read.** Rule (i) requires the held-out numbers to be re-read in the PDF by two readers before scoring. None of the 32 TEST contrasts has been. The atlas's per-cell wrong rate is about 9 % (95 % CI 4 - 18 %).
- **Heuristic parts:**
  - the numeric-effect test;
  - the DEVELOPMENT factor tags;
  - the alias detection, reviewed once against the cited text and frozen with its hash.
  
  Each is deterministic, and the generator lists what it excluded and why.
- **Error scopes were declared from secs. 3 and 5.** One was narrowed after the first build: E16's text marker went from "norfuraneol / NF / HDMF / furaneol" to norfuraneol only. E16 is norfuraneol's decay, and furaneol is another compound. The change removed one TEST contrast, Haleva 1999 F13#13, a furaneol pH contrast.
- **Two premortem amendments conflict with the global lab rule:**
  - the E6 + E11 prereg's "held-out tests come from the labelling cells of CAAS / BTBU" (both labs are DEVELOPMENT here);
  - the E4 prereg's own conclusion, which agrees: zero admissible held-out labs.
  
  The coordinator should decide whether rule (i) is applied per change (only the cells that diagnosed THAT change) or globally (v1). A per-change split would return CAAS / BTBU to E6's test set, and would need a v2.

## 7. Versions

- **v1, 2026-11-01:** first freeze, sha256 `022d0cab5a53ececb27a1fc6837ff275fbab177cf2161ae1f8a4ed0341f273a7`.
  - lab_ids_v1;
  - atlas extent as of the 2026-11-01 fix-up (562 cells);
  - ENGINE_VS_LITERATURE including sec. 5.1 (E4 correction);
  - the two preregs as of 2026-11-01 ~03:00.
