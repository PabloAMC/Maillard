# Engine vs literature: which Maillard reactions happen under which conditions

2026-11-01, overnight product. Machine-readable twin: [engine_vs_literature.yml](engine_vs_literature.yml).

**Inputs, read as they stood on 2026-11-01 01:20 - 02:40:** the literature atlas `F01 - F32_*.yml` (561 cells), the vocabulary [families.yml](families.yml), the engine map [engine_family_map.yml](engine_family_map.yml), the agents' mismatch notes in `tasks/overnight_log_2026-11-01.md` (used as LEADS only), and `docs/reference/conflicts_register.yml`. Source code under `src/kinetic_core/` was read, not changed.

**Recount on the fixed atlas (2026-11-01).** Sec. 1 and sec. 4's counts were re-run on the atlas as fixed in commit 6ba6583d (562 cells; [README.md](README.md)), with NO in-analysis overrides. Secs. 2, 3 and 5 are the original analysis; sec. 6 flags which of their mismatches rest on cells the fix-up changed.

**Red-team rules** (coordinator ruling after the 10 % atlas audit: ~9 % of random cells wrong, 95 % CI 4 - 18 %). The original analysis (secs. 2 - 3) applied the first three as overrides; since the fix-up they live INSIDE the files, and the recount applies none of them:
- **Ten known-wrong cells:** F02#12, F07#14, F10#0, F23#15, F24#13, F25#6, F25#7, F29#4, F30#15, F31#3. The original analysis excluded them (551 cells). The fix-up corrected nine and deleted F24#13; F25#6's pH 4 arm became F25#22, and F25#23 (Zhao 2019b, hexose FFT) is new.
- **Temperature bins from `T_C`:** low <= 60, moderate <= 100, high <= 140, very_high > 140 C. The original analysis used range midpoints (48 cells moved); the files now key a range in a bin it spans, checked by the gate. Water is now binned on all solutes, buffer included.
- **`absent` needs a `detection_basis`** (a detection limit, a positive control, or a label pattern), otherwise the cell is `not_detected`. The original analysis approximated this with a text match.
- **`dominant` weighted only when it rests on labelling or fed flux.** F25 dominant cells other than Cerny's fed-furfural pot, F31 dominant cells, and F30 quantities are treated as `operating`. (A reading rule; the classification below pools dominant and operating anyway.)
- **Confidence:** a mismatch is HIGH confidence when it rests on the labelling or fed cells of F23 / F24, F16 / F17, F19 or F01 - F06. Every other mismatch is LOWER confidence.
- **F23 / F24 water bins:** the original caveat (0.5 M buffer omitted) is retired; the fix-up rebinned those pots `concentrated`.

**Verification:**
- "Verified" means three things were checked: I read the cited code lines; I read the cited atlas cells (contrast / quantity / where); and, for the five claims marked *probe*, I ran a small read-only flux-budget probe of the SHIPPED sulfur lane (clamped pH; appendix A).
- "Lead" means an agent's statement I did not re-derive.
- **Two agent leads came out WRONG and are corrected here** (sec. 2.0).

---

## 1. What happens under which conditions (the literature)

### 1.1 Rule and counts

**A regime** is one cell of families.yml's three primary axes: pH (4 bins) x T (4) x water (4) = 64 per family, so 2048 family-regimes. Tags (sugar, amine, O2, buffer, metals) are carried, not binned.

The three statuses:
- **SETTLED:** >= 2 independent laboratories with labelling, fed-intermediate, trapping or within-study-contrast evidence, pointing the same way. The direction classes are {dominant, operating} vs {minor, absent-with-LOD}.
- **UNSETTLED:** a single laboratory; or mixed / conflicting statuses; or kinetic-fit / qualitative / mass-balance evidence only; or only unattributed or not-detected cells.
- **UNKNOWN:** no cell at all.

Independent laboratories were assigned by a fixed lab-name map, kept in the YAML. One lab line each:
- Cerny at Nestle and at Firmenich;
- Hofmann / Schieberle / Grosch / Granvogl at DFA / TUM;
- Mottram / Elmore / Parker at Reading and Bristol;
- Jiangnan and the Ho / Rutgers network;
- Wageningen and RIKILT.

| | count (fixed atlas) | share of 2048 | original analysis |
|---|---|---|---|
| **settled** | **25** (all present) | 1.2 % | 31 (30 present, 1 minor) |
| **unsettled** | **233** (of which 30 are only unattributed / not-detected, 13 mixed or conflicting) | 11.4 % | 217 (33 / 10 on the grid; the earlier "41 / 11" counted off-grid bins too) |
| **unknown** | **1790** | 87.4 % | 1800 |
| also: bins with a primary axis not reported (`not_stated`) | 42 (3 of them settled) | - | 41 (3 settled) |

What moved the counts (31 -> 25 settled):
- **Ten settled bins lost.** Eight are bin splits of claims that still pool (sec. 1.3): the 0.5 M-buffer pots of Nestle / Firmenich and CAAS moved `dilute` -> `concentrated` (F23 x2, F26, F27 x3), and F07#10 / F27#9 moved `high` -> `moderate`. The other two:
  - F16 acid / low: the Porto cell (40 - 80 C) now sits in `moderate`, leaving one lab;
  - F10 neutral / low (was the only `Z`): F10#0 is back as `operating`, because the fixed atlas caps "negligible, data not shown" at `operating`. Its content (MGO negligible) agrees with the minor cells, so the `!` is an artefact of the direction classes, not a disagreement.
- **Four gained**, all from the same rebinning: F07 neutral / moderate, F13 mild-acid / very-high / concentrated, F23 mild-acid and neutral / very-high / concentrated.
- **New `!` bins:**
  - F23 mild-acid / moderate / concentrated mixes sub-claims: intact-C5 MFT `dominant` beside the minor / absent free-NF route (P11 vs P13). Not a disagreement.
  - F16 acid / moderate pairs Quest's fed MGO with Porto's quinone donor: a different donor class with the opposite pH order, so also a sub-claim mix.
  - F27 mild-acid / moderate carries CF-089 (does cysteine stop thiamine degradation?): a real conflict.
  - F24 mild-acid / very-high is the old CF-086 conflict, moved from `dilute` to `concentrated`.

Calibration notes on these counts:
- **The rule is strict**, so the counts under-state how settled the textbook core is. F01 (Schiff), F02 (Amadori) and F05 (glycation) rest on mass balance and kinetic fits in 2 - 8 labs, and none of them reaches SETTLED in any bin.
- **The binning splits single claims.** A claim shown by three labs at 95, 140 and 168 C lands in three bins, one lab each. Sec. 1.3 therefore also lists **pooled claims**: the same claim, same direction, across adjacent bins.
- **Settled bins can mix sub-claims.** A settled bin may combine two labs that showed different sub-claims of the same family, both "operating". I would put P ~0.2 that any given settled bin overstates agreement for that reason.
- **The atlas itself was wrong in ~9 % of cells** (red team). The sampled wrong cells are fixed, but the ~460 unsampled cells still carry that rate (expect ~40 wrong, 95 % CI ~20 - 80), and this propagates.

### 1.2 Compact table (pH x T; water collapsed)

How to read the table:
- **Each pH column** holds four characters, one per temperature bin: L (<= 60 C), M (60 - 100), H (100 - 140), V (> 140).
- **Glyphs:**
  - `S` = settled (present);
  - `Z` = settled (minor / absent with LOD);
  - `!` = mixed / conflicting;
  - `u` = unsettled;
  - `a` = only unattributed / not detected;
  - `·` = unknown.
- **The best glyph across the four water bins** is shown.
- **The water column** lists the bins with any cell: d dilute, c concentrated, i intermediate a_w, y low a_w.
- **The `ns` column** counts bins with an unreported axis.
- **"Tags that matter"** are the tag conditions the atlas shows move the family. They are summarised from each file's `condition_dependence`.

| family | acid `LMHV` | mild_acid `LMHV` | neutral `LMHV` | alkaline `LMHV` | water | ns | tags that matter (evidence) |
|---|---|---|---|---|---|---|---|
| F01 Schiff base | `····` | `··au` | `uu·a` | `····` | dc | 2 | phosphate x6 (not citrate); amine type unknown beyond Gly / Phe / Asn |
| F02 Amadori / Heyns | `····` | `·u··` | `·uu·` | `·uu·` | dc | 4 | phosphate speeds formation AND decay; peptides x5 - 7 (1 lab); Heyns not detected (no LOD) with fructose + casein |
| F03 amine release / reversion | `····` | `·uS·` | `u!u·` | `··u·` | d | 0 | phosphate x1.3 - 12 (4 labs, fades at pH 8); pentose ARP returns less amine (1 lab); O2 x4 at 37 C |
| F04 thiazolidine / TTCA | `uu··` | `uuu·` | `uuS·` | `u··u` | dcy | 2 | amine (Gly / Lys speed drain); phosphate (CF-011); dilute vs concentrated (CF-076); in-situ pH |
| F05 protein glycation | `····` | `··u·` | `uu!·` | `····` | d | 5 | a_w (powders, 3 labs); protein identity; lactose = glucose; phosphate (37 C) |
| F06 1,2-enol -> 3-deoxyosone | `u···` | `·Su·` | `!uu·` | `····` | dc | 2 | O2-INdependent; phosphate / metals at 37 C; amine raises it |
| F07 2,3-enol -> 1-deoxyosone | `a···` | `·uu·` | `uSu·` | `·uu·` | dc | 2 | air lowers 1-DG x3.4 (50 C); amine structure; every flux datum in phosphate |
| F08 oxidative osone | `u···` | `····` | `Suu·` | `····` | dc | 1 | O2 x8; metals (-87 % with chelator); phosphate x14 (GO) |
| F09 Namiki C2 | `····` | `····` | `a·a·` | `·aa·` | dc | 2 | needs amine; aqueous route never marked |
| F10 retro-aldol | `····` | `·ua·` | `!uu·` | `·ua·` | dc | 2 | strong T dependence (bounded <= ~1.7 mol % at 37 - 50 C, material >= 100 C), no single-design measurement |
| F11 beta-cleavage (acetic) | `··u·` | `··ua` | `uu!a` | `·Su·` | dy | 0 | pH (fed 1-DG x4 lower at pH 5); amine not required; total O2-independent; phosphate no yield effect |
| F12 alpha-oxidative cleavage | `····` | `····` | `uu··` | `·uu·` | d | 0 | O2 defining; metal confound |
| F13 furans / furanones | `·u··` | `·uuS` | `·uau` | `··u·` | dcy | 0 | furfural pH (falls >> 10x from 4.5 to 6.5); amine needed for furanones; phosphate x4.3 (HDMF) |
| F14 amine-free degradation | `·S··` | `··uu` | `·Suu` | `·u··` | dcy | 0 | phosphate REQUIRED at pH 7 (1 lab); ketose x3 |
| F15 reductones / redox | `····` | `·au·` | `·uuu` | `·u··` | dcy | 0 | O2 + metals for norfuraneol loss; alpha-amino acid for reduction |
| F16 Strecker aldehyde | `u!··` | `·ua·` | `·uau` | `····` | dciy | 3 | donor identity (MGO > GO > 3-DG); amino-acid identity; low water raises Strecker amine; O2-independent from fed dicarbonyl |
| F17 Strecker acid / ARP O2 route | `·a··` | `····` | `·ua·` | `·u··` | dc | 0 | O2 (77 - 82 % of acid; ARP route x9); Cu x1.3 - single lab |
| F18 lipid Strecker | `uuuu` | `u···` | `u·u·` | `····` | dy | 1 | O2 x4 - 9; donor class - single lab |
| F19 aminoketone pyrazine | `··u·` | `·uuu` | `uu!u` | `·uSu` | dciy | 0 | amine (Lys >> Asn x35; peptides); phosphate x1.8; a_w max ~0.75; supply pH vs step pH differ |
| F20 extra pyrazine routes | `··u·` | `··uS` | `·u·u` | `·au·` | dciy | 0 | free amino acid's aldehyde needed (peptides starve it) |
| F21 pyrroles / pyridines | `uuu·` | `u···` | `·auu` | `u···` | dciy | 0 | NH3 / amide N; oil medium x10; concentration |
| F22 cysteine -> H2S | `·S·u` | `·u··` | `·S·S` | `·u··` | dc | 0 | pH 5 -> 9 x9 - 17; carbonyl partner x4.6 - 47 (1 lab, 1965); reduced thiol required |
| F23 MFT from intact C5 | `···a` | `·!uS` | `··SS` | `····` | dc | 1 | thiamine competes (F27); buffer x2 - 8; pentose ~10x hexose but hexose MFT is ALSO a C5 |
| F24 MFT from C2 + C3 | `···u` | `·u·!` | `··uu` | `····` | dcy | 0 | <= 1 - 6 % of MFT wherever labelled; fed pair rises to pH 7 (opposite of total MFT) |
| F25 FFT (furfural + H2S) | `···u` | `·uuu` | `··Su` | `····` | dcy | 1 | pH (FFT falls 6 - 19x acid -> 6.5 - 7); sugar (pentose 4x hexose); lipid (CF-031); dry / glycerol |
| F26 thiazoles / thiophenes | `····` | `·uuu` | `··!u` | `·uu·` | dcy | 0 | pH flips mercaptoketones down / thiazoles up; PUFA lipid -> alkylthiophenes; dry -> thiazolines |
| F27 thiamine / GSH sources | `··uu` | `·!!u` | `·uuu` | `··u·` | dci | 1 | loading (thiamine share falls as sugar + Cys rise); thiamine H2S needs pH >= 6.5; water needed |
| F28 thiol losses | `····` | `uuSa` | `SSu·` | `····` | dc | 0 | O2 x15 - 45 (25 C only); matrix (melanoidin / protein); cysteine protects; lipid 14 - 143x (MFT-selective) |
| F29 acrylamide | `···u` | `···S` | `u·uu` | `···u` | dcy | 6 | sugar (fructose > glucose, 3 labs); buffer (citrate x0.14 - 2.4 vs phosphate); a_w on elimination; amines compete |
| F30 lipid crossover | `····` | `··S·` | `S··S` | `u···` | dc | 2 | lipid unsaturation / oxidation; oil phase alone ~x10; NH3 / amide N; O2 |
| F31 melanoidins | `····` | `uuu·` | `uuu·` | `·u··` | dc | 4 | water (dry doubles sugar per N); reactant ratio; O2 for the thiol sink (25 - 30 C) |
| F32 DHA cross-links | `····` | `····` | `·u·u` | `····` | d | 1 | nothing beyond whey at pH 7 |

What the table says at a glance:
- The literature is concentrated in **aqueous pots at mild-acid to neutral pH, 60 - 170 C**: dilute, or `concentrated` only through a 0.5 M buffer (the Nestle / Firmenich, CAAS and DFA sulfur pots).
- The **alkaline column is nearly empty** except for F11 / F19.
- **Intermediate a_w (0.5 - 0.9) is almost unstudied.** It has cells in only F05, F16, F19, F20, F21 and F27, and none of them is settled. That is the moisture range of meat analogues and crusts.
- **Low T (<= 60 C) is settled only for physiological / storage chemistry** (F08, F28, F30). F10's minor retro-aldol at 37 - 50 C still pools (P5), but its bin reads `!` (sec. 1.1); F16's low-T bin lost its second lab to `moderate`.

### 1.3 Pooled settled claims

Each claim below has >= 2 independent labs of strong evidence in the same direction, across adjacent bins. Cells are written file#index.

**Recount on the fixed atlas:** all 22 claims survive. P12 and P13 are firmer: their non-detects are now `absent` with a stated label pattern or positive control. P5 is reworded (bounded small, not nil). P16 rests on two lab lines, not three, and no longer has a settled bin under it. P21 regains the Leeds cell.

| # | claim | regime | labs (evidence) | cells |
|---|---|---|---|---|
| P1 | A hexose Amadori compound of a primary amino acid returns most of its amine as it decomposes (71 - 96 %) | pH 5 - 8, 90 - 130 C, dilute | Nestle, Wageningen, Jiangnan (fed) | F03#0,1,3,5,10 - 14 |
| P2 | Phosphate accelerates early condensation and Amadori decomposition (x1.3 - 12), fading at pH 8 | pH 5 - 7.4, 25 - 130 C | Nestle, Baynes, SCUT, Auburn, Jiangnan (within-study, fed) | F01#7,9; F03#2,7,10; F05#22 |
| P3 | The 2,3-enolisation (1-deoxyosone) share RISES with pH; the 3-deoxyosone flux share falls from ~20 - 40 % (pH <= 5) to ~5 % (pH 7). Hexose, phosphate | 100 - 120 C | Ledl school (trapping), Nestle (labelling + fed), Wageningen (fit) | F06#0,1,2,13; F07#0,3,5,6,15 |
| P4 | Hydrolytic beta-cleavage of 1-DG is the main acetic-acid route | pH 6 - 8, 37 - 120 C | Nestle (labelling, 18O2, fed), Glomb (fed), Wageningen | F11#1,5,8,9,12 |
| P5 | Deoxyosone retro-aldol (1-DG, 3-DG -> MGO) is minor at 37 - 50 C: MGO seen at negligible amounts; C3 + C3 cleavage of 1-DG bounded <= ~1.3 - 1.7 mol %; none from fed 3-DG (positive control) | neutral, low T | Glomb, Thornalley (fed) | F10#0 - 4 |
| P6 | Furfural is made from the intact pentose C5 and falls steeply from pH 4.5 to 6.5 | 120 - 145 C | Firmenich, Jiangnan (labelling); Reading, Garching (within-study) | F13#0,1,2,5,6,7 |
| P7 | Fed alpha-dicarbonyls give Strecker aldehydes with ranking MGO > GO > 3-DG, O2-independent; the MAGNITUDE conflicts ~20x (CF-073) | pH 3 - 7, 98 - 170 C | DFA, Quest, Missouri (fed) | F16#0,1,5,7,9 |
| P8 | Pyrazine ring carbons come from dicarbonyl / sugar carbon; ring N from two aminated halves (aminoketone pairing) | 120 - 180 C | Nestle, Reading (13C); Ho network (15N) | F19#0,1,3 |
| P9 | Aldehyde addition to dihydropyrazines makes ethyl / extra-methyl pyrazines in roasts | ~180 C, low a_w | Nestle, Reading (13C), DFA (fed) | F20#0,1,2,4 |
| P10 | Cysteine thermolysis H2S rises with pH 5 -> 9, Ea ~123 - 135 kJ/mol | 80 - 110 C | Rutgers; Hokkaido, Copenhagen confirmations | F22#0 - 5 |
| P11 | **MFT carbon is an intact C5 of one sugar, INCLUDING from glucose** | pH 5 - 6.5, 95 - 168 C; 0.5 M buffer (`concentrated`) except BTBU | Nestle / Firmenich, BTBU, CAAS (labelling) | F23#0,1,6,7,8 |
| P12 | **C2 + C3 recombination carries <= 1 - 6 % of MFT** wherever labelled | 95 - 168 C | Nestle / Firmenich, BTBU, CAAS (labelling; F24#0,1,4,5 `absent` on the label pattern) | F24#0,1,3,4,5 |
| P13 | **Free norfuraneol is a poor MFT source:** 7 % of MFT when fed at the ribose's charge (95 C); 48x below ribose (boiling, concentrated); NF : MFT 3000 in a ribose pot | 95 - 145 C | Nestle (labelling; F23#9 `absent`, positive control in pot D), Givaudan (fed), DFA (within-study) | F23#2,9,11,13 |
| P14 | FFT is made from furfural + H2S in pentose + Cys pots, and FFT falls only modestly from acid to pH 6.5 - 7 (6 - 19x) while furfural collapses | pH 3 - 7, 140 - 145 C | DFA (fed, within-study), Reading (within-study) | F25#4,8,9 |
| P15 | Mercaptoketones collapse and thiazoles / thiophenones rise from pH 4.5 to 7 | 140 - 145 C | Reading, Firmenich (within-study + labelling) | F26 ladders (#6 - 11) |
| P16 | Thiamine is the dominant MFT source at dilute loadings; its share falls as pentose + Cys rise | pH 4 - 7, 120 - 145 C | Firmenich, Minnesota / Nestle (34S), BTBU (13C). The lab map makes Firmenich and Minnesota / Nestle one lab line, so two lines; F27#7 (Nebraska) is a thiamine-sulfur sub-claim | F27#0 - 7 |
| P17 | Furyl disulfides are a REVERSIBLE reservoir (thiol returns in water with a trap, and in beef) | ~100 C | DFA (trap + label), Reading (mass balance) | F28#1,21,22 |
| P18 | Acrylamide: fructose > glucose; formation rises toward pH 8 (in phosphate) | 120 - 200 C | Leuven, Stockholm, Wageningen | F29#1,5,9,21 |
| P19 | Lipid (PUFA, phospholipid, lard) suppresses MFT 3 - 143x in Cys + pentose pots; an oil phase alone ~x10 | pH 5 - 5.7, 100 - 140 C | Reading group, Huazhong (within-study) | F30#0,3,5,6 |
| P20 | 2,4-dienal + free NH3 / amide N -> 2-alkylpyridine | RT - 180 C | Kentucky, Rutgers (13C / 15N labelling); Seville (fed) | F30#11,13,14; F21#5,6 |
| P21 | Melanoidins take sugar C and amine N together, glycine partly decarboxylated (27 - 33 % or ~2/3: Fang 2009 vs Mundt, unregistered) | 70 - 125 C | Iowa State, Wageningen / Leeds (labelling) | F31#0 - 4 |
| P22 | Protein-bound lysine glycates in every regime looked at; flat over a_w 0.3 - 0.7, slower at 0.98 | 37 - 220 C | 8 labs (furosine / available lysine) | F05#0 - 24 |

### 1.4 The main UNSETTLED and UNKNOWN regimes

**Single-lab families** (every quantitative cell from one laboratory):
- F09 Namiki: Nagoya, ethanol only;
- F12 alpha-oxidative cleavage: Nestle for 18O2, Glomb at 37 - 50 C;
- F15 redox;
- F17 Strecker acid / ARP O2 route: DFA;
- F18 lipid Strecker: Seville;
- F21 pyrroles: Ghent;
- F32 DHA: Copenhagen.

**Open conflicts that bear on the engine:**
- CF-073: Strecker magnitude, 20x;
- CF-031: lipid on thiols;
- CF-076 / CF-011: TTCA exchange, phosphate;
- CF-072: glycation lysine return;
- CF-069: acetic-acid loss;
- CF-013: thiol time courses;
- CF-081: melanoidin site density.

**Unknown across nearly all families:**
- intermediate a_w with a measured a_w;
- in-situ pH during heating (almost every pH is a label pH);
- O2 vs inert at >= 100 C (no family has one);
- alkaline in-situ pH at cooking T;
- plant-protein-bound amines;
- pentose formation splits (F06 / F07) by in-situ trap;
- the 100 - 140 C bin with 1:1 pentose labelling (F23 / F25 / F26);
- any labelling of acrylamide.

---

## 2. Engine vs literature, per family

As originally built (overrides as listed at the top). Which entries rest on cells the atlas fix-up changed: sec. 6.

Mismatch kinds:
- **M1** = the literature shows the family operating / dominant, but the engine lacks the route or ships it off;
- **M2** = the engine runs it where the literature shows minor / absent;
- **M3** = opposite sign or grossly different magnitude of a coded condition dependence;
- **M4** = the engine cannot express a condition that matters.

Status: **V** = verified (code read + cells read; *p* = also probed); **L** = lead. Conf. = H (high) / lo (lower), per the red-team rule.

### 2.0 Corrections to the agents' leads

1. **"FFT formation RISES with pH in the engine (a sign error)" is NOT confirmed as a net statement.**
   - The per-furfural HS- branch (`r_fur_fft_hs`, sulfur.py:526, `hs_anion` at :1187) does rise with pH, and so does the literature's step, which is "fast at neutral pH".
   - But net FFT in the shipped engine FALLS with pH, because the furfural supply falls ~1.5 decades per unit: the `acid` factor on `r_pent_tdp` / `r_arp_tdp` / `r_glc_fur`, sulfur.py:1160, :1163, :1175.
   - The real error is the MAGNITUDE (sec. 3, E6). In the clamped probe, FFT made falls x44 from pH 4.5 to 6 and ~x400 to pH 6.5, against x6.1 measured over 4.5 -> 6.5 (Meynier, pyrophosphate held, F25#9) and x19 over 3 -> 7 (Hofmann, F25#8).
   - The MFT / FFT ratio also moves the OPPOSITE way: engine 0.40 -> 3.1 from pH 4.5 to 6.5; measured 1.79 -> < 0.07 (F25#9).
   - At the acid end the engine's FFT is x6 LOWER at pH 3 than at pH 4.5, while DFA measures it x1.9 HIGHER at pH 3 than at pH 5 (F25#8).
2. **"The Cys-ARP route SHIPS dpo_share = 0.30" (F06 / F07 atlas engine notes) is wrong.**
   - `CYS_AMADORI_ROUTE.enabled = False` (parameters_sulfur.py:689; switched off in `inactive_reactions`, sulfur.py:972 - 976).
   - The standing-ratio-used-as-flux-split error (dpo_share, parameters_sulfur.py:704) is LATENT: it bites when the route is enabled, not in shipped predictions.

### 2.1 Mismatch catalogue

| id | fam | kind | literature (cells; labs; evidence) | engine (file:line) | status | conf |
|---|---|---|---|---|---|---|
| F01-a | F01 | M3 | reversible pre-equilibrium, reversion faster than Amadori (F01#0,1,3; Rutgers Ge, Winterswijk; mass balance) | `r_schiff` irreversible, no reverse step (network.py:135); `r_amadori` derived 44.9x (network.py:688 - 689) | V | H (but mass balance only) |
| F01-b | F01 | M4 | phosphate x6, not citrate (F01#7 Auburn; F01#9 SCUT; within-study) | no buffer term anywhere (grep src/kinetic_core: buffer only as titrant, ph_state.py) | V | H |
| F02-a | F02 | M4 | formation ~0.5 dec / pH (F02#2 Wageningen fit) | no pH term on `r_schiff` / `r_amadori`; PH_STEPS covers decays only (trunk_conditions.py:126) | V | lo (fit) |
| F02-b | F02 | M4 | peptides form ARP 5 - 7x faster (F02#5 Jiangnan) | one amine; PRO / MET folded into Gly (engine.py:2560, :2565) | V | lo |
| F03-a | F03 | M3 | GX-ARP returns 25 - 40 % glycine (Hou, Jiangnan; fed; single lab); hexose ARP 71 - 96 % (3 labs) | pentose ARP decays book N to FRAG_N (sulfur.py:337 - 385); disabled PAMA route returns 100 % (parameters_sulfur.py:540): the two pentose routes contradict | V (code) | H (single lab) |
| F03-b | F03 | M4 | phosphate speeds decomposition 1.3 - 12x (F03#2 Nestle, #7 Baynes, #10 SCUT; fed) | no buffer catalysis | V | H |
| F03-c | F03 | M1 | reversion ~20 % at 37 C (F03#7 Baynes), minor 3 - 8 % at 90 - 120 C | no AMA -> Glc + Gly step | V (code) | H (low impact) |
| F04-a | F04 | M1 | Cys + pentose binds 18 - 89 % of Cys as TTCA before / within 0.5 h (F04#11,14,15,19,23; BTBU + Jiangnan; settled bin) | `TTCA_FORMATION.enabled = False` (parameters_sulfur.py:422), off in sulfur.py:972 - 973 | V | H |
| F04-b | F04 | M3 | ring release k_r(100 C, pH 4.5) ~2 - 8e-4 /min (cluster 5, provisional P 0.45; engine exchange rejected in two labs) and fed-TTCA loss x10 over pH 5.5 -> 8 (Zhai 2021) | `k_ttca_cys` 0.0286 /min at 100 C, pH-flat (no REACTION_PH_FACTOR entry) | V (engine value) / L (literature value) | H |
| F05-a | F05 | M4 | glycation x2.25 slower at pH 5.9 than 6.8 (F05#3 Brands); a_w dependence in powders (3 labs) | no pH or a_w term on `k_glyc`: AW_STEPS = (k_schiff,) (trunk_conditions.py:82), PH_STEPS (:126) | V | lo |
| F05-b | F05 | M3 | CF-072 (error_confirmed): constant fitted to a lysine-not-returned network | `r_flp_decay` | L (registered) | lo |
| F06-a | F06 / F07 | M3 | split moves ~0.25 - 0.45 dec / pH toward 1-DG (F06#0,1,2,13; F07#0,3,5,6,15; Ledl trapping, Nestle labelling + fed, Wageningen fit) | ONE exponent 0.69 on k_ama_tdg / k_ama_odg / k_ama_mgo (trunk_conditions.py:121, :126), so the split is pH-INVARIANT; the code's own source table prints per-step slopes 0.37 / 0.92 / 0.72 at 100 C (trunk_conditions.py:116 - 120) and pools them | V | H |
| F06-b | F06 | M3 | absolute 1,2-enolisation (Martins k4) RISES +0.35 dec / pH (F06#13, hexose fit) | pentose / ARP 1,2-enolisation labelled `acid` = [H+], FALLS 1 dec / pH (sulfur.py:1160, :1163) | V (code) / L (transfer hexose -> pentose) | H |
| F06-c | F03 / F06 / F07 / F10 / F11 / F13 | M4 | in-situ pH drifts 1 - 3 units in unbuffered pots | the sulfur lane's 29 trunk steps get no pH: `TRUNK_PH_TERMS = False` (sulfur.py:1278) | V | H |
| F07-a | F07 | M4 | air lowers trapped 1-DG x3.4 (F07#4 Glomb, 50 C) | no O2 term on `r_ama_odg` (network.py:161) | V | H (low T only) |
| F06-d | F06 / F07 | latent | standing 3-DX / 1-DX used as a flux split | dpo_share 0.30 (parameters_sulfur.py:704), route OFF (:689) | V | - |
| F08-a | F08 | M4 | glucosone x8 air vs deaerated; chelator -87 % (F08 neutral / low cells; Glomb, Thornalley, Baynes) | ATMOSPHERE_FACTORS = {"air": 1.0} (trunk_conditions.py:275); argon refused; no metal axis | V | lo (37 - 50 C only) |
| F09-a | F09 | M1 | Namiki C2 only in ethanol; aqueous unattributed | no SB -> C2 step; pentose lane has no glyoxal source | L | lo |
| F10-a | F10 | M2 | 1-DG retro-aldol negligible from fed 1-DG at 37 - 50 C (F10#1 - 4; Glomb, Thornalley) | `r_dpo_c2c3` 1 : 1 HA + MGO from 1-DPO (sulfur.py:425), on the hexose analogue; no pentose fed test at 100 - 140 C | L (regime mismatch: 37 - 50 C vs 100 - 145 C) | H |
| F10-b | F10 | M3 | aldose retro-aldol gives C2 + C4, ketose C3 + C3 | `r_glc_c2c3` Glc -> HA + MGO + C1 (sulfur.py:430) | V (code) | lo |
| F11-a | F11 | M3 | fed 1-DG acetic acid x4 lower at pH 5 than >= 6.5 (F07#5 Nestle) | `r_odg_aa` pH-flat (not in PH_STEPS, trunk_conditions.py:126) | V | lo |
| F11-b | F11 | M3 | <= 1.7 mol % formic from fed 3-DG (Davidek 2006a) | `r_tdg_fa` ~13 % of 3-DG at 120 C (agent's arithmetic) | L | lo |
| F11-c | F11 | M2 | acetic acid plateau (Davidek 2008) | `r_aa_frag` destroys AA (network.py:188; CF-069 open) | L | lo |
| F12-a | F12 | M1 | O2 alpha-cleavage carries ~1/3 of acetic acid under 18O2 (Nestle), and C3 acids x4 - 6 (Glomb) | no O2 acid route; OX consumers only Cys / reductone (engine_family_map) | V (code) | lo |
| F13-a | F13 / F15 | M3 / M1 | NF decay base- and O2-catalysed: 0.032 /min at pH 8, 100 C (Mikami; F13#8,9 Ochanomizu); NF falls 1.7 - 16x pH 4.5 -> 6.5 at 140 C (F13#6 Reading) | `NF_HYDROXIDE_DECAY.enabled = False` (parameters_sulfur.py:447); shipped `k_nf_decay` 4.5e-5 /min at 100 C, pH-free (sulfur.py:537): ~700x slower at pH 8 | V | lo |
| F13-b | F13 | M4 | HDMF x4.3 phosphate over malonate, x43 over water (F13#11 Nestle) | no buffer term | V | lo |
| F13-c | F13 | M3 | furfural Ea 147 is a net-accumulation Ea from hexose + Cys (F13#15) | shipped `r_tdp_fur` shares 64.08 kJ/mol (FURFURAL_BARRIERS.enabled False, parameters_sulfur.py:2860) | V | lo |
| F14-a | F14 | M4 | amine-free browning needs phosphate at pH 7 (F14#4 Rizzi) | amine-free lane runs identically in water and phosphate | V | lo |
| F15-a | F15 | M1 | Strecker-coupled reduction generic to alpha-amino acids (Ledl, qualitative) | `r_dpo_ddp` only with Cys (sulfur.py:393) | V (code) | lo |
| F15-b | F15 | M1 | furfuryl alcohol rises with pH to ~furfural level (Reading, Minnesota) | no FUR -> furfuryl alcohol; lost in `r_fur_decay` (sulfur.py:530) | V (code) | lo |
| F16-a | F16 | M1 | Strecker aldehydes of Leu / Ile / Val / Phe / Met from fed dicarbonyls (F16#0,1,2,5,7; DFA, Quest; fed) | NO Strecker-aldehyde species shipped: Strecker only GO / MGO + Gly -> aminoketone + HCHO (network.py:350, :356); methional off (`METHIONINE_SHIPPED = False`, parameters_methionine.py:59) | V | H |
| F16-b | F16 | M1 | 3-DG and glucosone donors give 73 - 76 % of MGO's aldehyde (F16#0 DFA; F16#9 Missouri) | no 3-DG / glucosone Strecker step (network.py reaction list) | V | H |
| F16-c | F16 / F19 | M3 | fed-dicarbonyl PA peaks at pH 5, falls only 2x to pH 9 (F16#3 DFA); largest at pH 3.2 (F16#5 Quest); fed MGO + Ala 2,5-DMP x2.1 HIGHER at pH 4 than 7 (F19#7 Ghent) | Leahy pyrazine-cascade term on the Strecker steps k_go_ak / k_mgo_ak (parameters_pyrazine.py:153 - 154; applied trunk_conditions.py:357 - 362): x0.091 at pH 5 vs 6.8 | V | H |
| F16-d | F16 | M4 | low water raises the Strecker amine to 19 - 27 % (F16#8 Granvogl) | no a_w term; no Strecker-amine exit | V | H |
| F17-a | F17 | M1 | O2 carries 77 - 82 % of the Strecker acid; fed ARP x9 in air (F17#0 - 7 DFA; single lab) | no Strecker acid species; `r_marp_mtal` not shipped (B22B_SHIPPED False, parameters_methionine.py:103), no O2 term | V | H (single lab) |
| F18-a | F18 | M1 | lipid carbonyls give Phe aldehyde 8 - 17 % at 180 C (Seville only) | lanes disjoint (lipid.py:784 - 830) | V | lo |
| F19-a | F19 / F20 | M3 | statistical cross pair (R = 4 by design; Leahy R 7.5 - 52) | CF-077 error_confirmed: `r_ak_mpz` takes k_cond unscaled, so R = 1 (network.py:370, rate law :695 - 704) | V | H |
| F19-b | F19 | M3 | the only measured aminoketone condensation (ALA analogue) ~1e7x slower, Ea 44, x80 over pH 4.8 -> 7.4 (F19#22 Karolinska; analogue) | k_cond 1e4 L/(mmol min), Ea 0, no pH (parameters_pyrazine.py:52, :121 - 123) | V (code) / L (transfer) | H |
| F19-c | F19 | M4 | phosphate x1.8 over bicarbonate (F19#6 Ghent) | no buffer term | V | H |
| F19-d | F19 | M4 | a_w maximum ~0.75 (F19#11 Leahy) | AW_STEPS = (k_schiff,) only (trunk_conditions.py:82) | V | lo |
| F19-e | F19 | M4 | Lys >> Asn x35 (F19#10); dipeptides 20 - 34x (Ghent) | one amine, alanine data declared as glycine | V | lo |
| F20-a | F20 | M1 | aldehyde addition / extra MPZ route (F20#0,1,4 settled roast) | no dihydropyrazine, no addition step; Strecker HCHO to FRAG_C (network.py:350) | V | H |
| F21-a | F21 | M1 | aminoacetone + MGO pyrrole is the main aminoacetone sink at pH 4, ~0.85 M (Ghent, fed); dienal + NH3 -> 2-pentylpyridine (P20) | no pyrrole / pyridine / imidazole; NH3 not a species | V | H (pyrrole: single lab) |
| F22-a | F22 | M1 | carbonyl partner raises H2S 4.6 - 47x (F22#11 Tokyo, 1965, single) | `STRECKER_SULFUR.enabled = False` (parameters_sulfur.py:657); MGO-Cys Strecker k = 0 | V (code) | lo |
| F22-b | F22 | M3 | 160 C / 160 s H2S 0.05x the extrapolation (F22#8 Copenhagen, one point) | `k_cys_h2s` Ea 133 extrapolated above 110 C | L | lo |
| F23-a | F23 / F24 | M1 + M2 | hexose MFT is an intact glucose C5, < 1 % mixed channels (F23#7 BTBU, F23#8 CAAS, labelling); C2 + C3 <= 1 - 6 % (F24#3 - 5); ribose / glucose MFT 10.4x (F23#21 DFA) | hexose C5 entry OFF: `HEXOSE_C5_ENTRY.k_per_min_100c = 0.0` (parameters_sulfur.py:633; sulfur.py:927, :992, :1473); hexose MFT only via `r_glc_c2c3` (sulfur.py:430, log10 k = -8.59 at 145 C). Probe: glucose MFT ~1e-13 mM, ribose / glucose MFT 4.4e10 (measured 10.4). The envelope flags it only for hexose-only pots (engine.py:1024 - 1043) | V *p* | H |
| F23-b | F23 | M2 (sub-route) | free NF carries 7 % of MFT (F23#2), its marker is absent at 95 C (F23#9), fed NF 48x below ribose (F23#11), NF : MFT 3000 (F23#13) | log10 k(145 C): `k_nf_mft` -2.67 vs `k_ddp_mft` -6.54 (7500x; kinetic_core_b60_fit_report.json). Probe: NF share of MFT 0.999 (95 C, pH 5) | V *p* | H |
| F23-c | F23 | M4 | buffer catalyses the MFT branch 2 - 8x (F23#20 Reading) | no buffer term | V | H |
| F23-d | F22 - F26 | M3 | the route mix shifts with T (NF marker absent at 95 C, on at 145 C pH >= 6: F23#9,10); Cys H2S Ea 123 - 135 (F22#2); furfural 147 (F13#15) | ONE lumped formation barrier 64.08 kJ/mol for every fitted sulfur formation step (b60 report `lumped_formation_Ea_kJ_mol`; no `formation_Ea_by_route`; engine.py:1889 - 1902) | V | H |
| F24-a | F24 | M2 | see F23-a | the engine's hexose MFT is 100 % C2 + C3 (probe C) | V *p* | H |
| F25-a | F25 | M3 (magnitude) | FFT falls 6.1x (pH 4.5 -> 6.5, held, F25#9) and 19x (pH 3 -> 7, F25#8) | probe: x44 (4.5 -> 6), ~x400 (4.5 -> 6.5), x3600 (4.5 -> 7); supply law `acid` on sulfur.py:1160, :1163, :1175 | V *p* | lo |
| F25-b | F25 | M3 (sign) | MFT / FFT 1.79 -> < 0.07 from pH 4.5 to 6.5 (F25#9) | probe: 0.40 -> 3.1 | V *p* | lo (single lab) |
| F25-c | F25 | M3 | pentose FFT 4.3x hexose (F25#10 DFA) | probe: glucose FFT 2.1x ribose (r_glc_fur 10^-4.92 vs r_tdp_fur chain) | V *p* | lo (single lab) |
| F25-d | F25 | M4 | dry 180 C x8, glycerol x56 | no a_w in the sulfur lane | V | lo |
| F26-a | F26 | M1 | 3-mercapto-2-pentanone from the intact ribose chain (F26#0 - 3 Cerny labelling) | MP2P only from thiamine HMP (`r_hmp_mp2p`, sulfur.py:561) | V | lo |
| F26-b | F26 | M1 | thiazoles / thiophenones rise with pH (P15) | `r_nf_shet` inert (k 0 unless a report supplies it); one thiazole (`r_cys_actz`) | V (code) | lo |
| F26-c | F26 / F30 | M1 | PUFA / dienal -> alkylthiophenes / thiapyrans (F26#21,24,25; Reading, CAAS labelling) | no lipid sulfur chemistry | V | lo |
| F27-a | F27 | M4 | thiamine loss 4 - 8 % (pH 3.5) -> 33 - 59 % (pH 8), 121 C (F27#8 - 11 Nebraska) | no pH factor on `r_thi_hmp` | V | lo |
| F27-b | F27 | M1 | thiamine sulfur reaches pentose-skeleton FFT at pH 5, 145 C without Cys (F27#1 Cerny) | `r_thi_h2s` share = 0 below pH 6.25 (sulfur.py:1404) | V (code) | lo |
| F27-c | F27 | M1 | GSH lost by gamma-glutamyl cleavage 7 - 14x faster than its H2S release (F27#25) | GSH has one reaction, `r_gsh_h2s` | V (map) | lo |
| F28-a | F28 | M3 | disulfides revert to thiol (F28#1 DFA, #21, #22 Reading) | `k_dimer_release = 0` (no block in the b60 report; parameters_sulfur.py:1490) | V | lo |
| F28-b | F28 | M3 | in-pot loss cannot be one constant (Grosch 0.044 /min vs patent <= 0.02 - 0.05, F28#12; CF-013) | `ch_mft_decay` first-order, pot-independent: 1.36 /min at 145 C, 0.041 /min at 100 C (sulfur.py:742) | V (code) / L (magnitude) | lo |
| F28-c | F28 | M4 | O2 x15 - 45 on thiol loss (F28#0 Nestle, 25 C) | OX seeded 1.0 mM in every sulfur pot, inert included (engine.py:2533) | V | lo |
| F29-a | F29 | M1 | fructose > glucose 1.1 - 2.9x (F29#1 Wageningen, #9 Leuven, #21 Stockholm) | Asn reacts only with Glc (`a_asn_glc_sb`, acrylamide.py:165); fructose only by isomerising | V | lo |
| F29-b | F29 | M4 / M3 | citrate k_F x0.14 (pH 4), x2.4 (pH 6) vs phosphate (F29#6 Leuven) | phosphate pH law on every buffer (acrylamide_conditions.py:122 - 125) | V | lo |
| F29-c | F29 | M3 | amine adducts ~90 % reversible (F29#17,18 Seville) | irreversible `a_acr_gln` / `a_acr_lys` / `a_acr_ala` (acrylamide.py:251 - 260) | V | lo |
| F30-a | F30 / F28 | M1 / M4 | MFT 3 - 143x down with lipid (F30#0,3,5,6; Reading group, Huazhong); reverse flow, enals x0.1 (F30#1,2) | DIRECT SUM; covalent sink off (lipid.py:784 - 830; parameters_lipid.py:551, :589): lipid changes nothing | V | lo |
| F31-a | F31 | M3 | C/N 6.6 - 21, condition-dependent (F31#5 Leeds, #8 TU Berlin / BAM) | fixed MEL_C 8 : MEL_N 1 (network.py:167) | V (code) | lo |
| F32-a | F32 | M1 | DHA / lanthionine at 80 - 160 C, whey (single lab) | absent | V | lo |

---

## 3. Ranked engine mechanism errors

**Score** = E (evidence: 3 = >= 2 independent labs, labelling / fed / within-study, same way; 2 = 2 labs mixed, or one lab with several marker papers; 1 = one lab or a fit) x O (odorant impact: 3 = moves an MFT / FFT / H2S, Strecker-aldehyde, pyrazine, lipid-aldehyde or acrylamide prediction by >~10x; 2 = 2 - 10x or a supply term; 1 = non-odorant) x R (recipes: 3 = most recipes; 2 = a recipe class; 1 = niche).

**Ranking:** confidence-weighted. Lower-confidence items are multiplied by 2/3 and then sorted.

**Implementability:**
- **band-now:** expressible with EXISTING parameters or switches. A declared band, enabling an existing switch with a declared band, or a fit constraint. Every one is pre-registered, and no new global switch.
- **A3:** needs a new condition axis, term or lane coupling, so it waits for A3's Config.

| rank | id | error (minimum correct structure) | E x O x R | conf | eff. | status | implement |
|---|---|---|---|---|---|---|---|
| 1 | E1 (F23-a / F24-a) | **Hexose MFT has no intact-C5 route.** Glucose + Cys gives ~0 MFT (ribose / glucose 4e10 vs 10.4 measured); what little there is runs 100 % through C2 + C3, which labelling shows carries <= 1 - 6 %. **Fix:** enable the existing hexose -> C5 entry (`r_odg_dpo`, 1-DG -> 1-DPO + C1) with a declared band anchored on Hofmann T1 (ribose / glucose MFT 10.4x) and the labelled C5 pattern; bound C2 + C3 <= 6 % of MFT | 3x3x3 = 27 | H | 27 | V *p* | band-now (HEXOSE_C5_ENTRY exists) |
| 2 | E2 (F06-a, F06-c) | **The hexose 1,2- vs 2,3-enolisation split cannot move with pH,** and the sulfur lane's trunk steps get no pH at all. **Fix:** per-step Martins exponents (k4 ~0.35, k7 ~0.60 dec / pH; the code's own Table 3 slopes 0.37 / 0.92 / 0.72) in place of the pooled 0.69. In the sulfur lane, apply them at in-situ pH (TRUNK_PH_TERMS) | 3x2x3 = 18 | H | 18 | V | band-now (per-step exponents in trunk_conditions; TRUNK_PH_TERMS exists) |
| 3 | E3 (F16-a, -b) | **No Strecker-aldehyde lane.** 3-methylbutanal, 2-methylbutanal, PA and methional are not shipped; donors are GO / MGO only. **Fix:** a per-amino-acid aldehyde product on the existing dicarbonyl + amino-acid steps, with 3-DG / glucosone donors at 0.7 - 0.8x MGO (DFA argon), and a declared x20 magnitude band for CF-073 | 2x3x3 = 18 | H | 18 | V | A3 (new species and amine axis) |
| 4 | E4 (F16-c / F19) | **The pyrazine-cascade pH law sits on the Strecker STEP** (x0.091 at pH 5), against fed-dicarbonyl Strecker peaking near pH 5 and fed MGO + Ala pyrazine x2.1 higher at pH 4. **Fix:** move the below-7 slope to the supply / condensation; the step gets <= 0.07 dec / pH, as the Strecker prereg allows | 2x3x3 = 18 | H | 18 | V | band-now (PYRAZINE_PH_SLOPES are parameters) |
| 5 | E5 (F01-b / F03-b / F13-b / F14-a / F19-c / F23-c / F29-b) | **No buffer-species catalysis anywhere,** while phosphate accelerates Amadori decay, condensation, amine-free degradation, HDMF, the MFT branch, acrylamide and pyrazines 1.3 - 12x (>= 7 labs). The constants were measured in phosphate, so unbuffered and citrate foods are mis-scaled. **Fix:** a buffer-anion factor on the Amadori-decay and enolisation steps (saturating, fading at pH 8) | 3x2x3 = 18 | H | 18 | V | A3 (new axis); interim: a declared transfer band on phosphate-calibrated constants in non-phosphate pots |
| 6 | E6 (F23-b) | **MFT runs 99.9 % through free norfuraneol** (measured 7 %; NF a poor source in 3 labs). **Fix:** re-weight `k_ddp_mft` vs `k_nf_mft` under a fit constraint NF share <= 0.15 at 95 C, pH 5 (Cerny pot D) | 3x3x2 = 18 | H | 18 | V *p* | band-now (fit constraint, no switch) |
| 7 | E7 (F30-a / F28-d / F26-c) | **The lipid lane is a direct sum:** the engine predicts no MFT / FFT change with lipid, vs 3 - 143x suppression (x10 from an oil phase alone) and the reverse removal of enals. **Fix:** a lipid-carbonyl + H2S / thiol capture and an oil-phase partition, jointly integrated | 3x3x3 = 27 | lo | 18 | V | A3 (lane coupling); interim: a declared envelope band on MFT when PUFA lipid is charged |
| 8 | E10 (F04-a, -b) | **TTCA formation is off,** so cysteine is never sequestered (18 - 89 % bound in 2 labs), and the release that IS on is 36 - 140x too fast and pH-flat. **Fix:** enable TTCA_FORMATION with Cao's K band; slow k_r to 2 - 8e-4 /min (cluster 5); K(pH) Kallen-shaped | 3x2x2 = 12 | H | 12 | V | band-now for K and k_r (existing switch / parameter); A3 for the K(pH) term |
| 9 | E13 (F19-a) | **CF-077:** methylpyrazine cross-pair lacks its statistical factor 2 (R = 1) | 2x2x3 = 12 | H | 12 | V | band-now (factor on r_ak_mpz; owner approval: predictions move) |
| 10 | E14 (F19-b, F21-a, F20-a) | **Aminoketone condensation declared instant** (k_cond 1e4, Ea 0); no aminoacetone pyrrole sink; no extra MPZ / aldehyde-addition route | 2x2x3 = 12 | H (pyrrole single lab) | 12 | V (code) / L (transfer of the ALA analogue) | band-now for k_cond band; A3 for new products |
| 11 | E9 (F29-a, -b, -c) | **Acrylamide:** glucose the only carbonyl (fructose 1.1 - 2.9x, 3 labs); a phosphate pH law on every buffer; irreversible amine scavenging. **Fix:** an `a_asn_fru` step at a declared 1.1 - 2.9x band; reversible amine adducts (K from Seville) | 3x3x2 = 18 | lo | 12 | V | band-now for fructose (new step, no switch); A3 for buffer |
| 12 | E11 (F23-d) | **One formation barrier (64.08) for every sulfur formation step,** so the route mix cannot move with T. **Fix:** per-route barriers through the existing B10 `formation_Ea_by_route` map and B62 furfural barriers, bounded by measured classes (Cys H2S 123 - 135, furfural 147, NF formation 81 - 93) | 2x3x3 = 18 | lo (fit) | 12 | V | band-now (B10 map / B62 switch exist) |
| 13 | E12 (F28-a, -b, -c) | **Thiol losses:** disulfide terminal (release 0, vs reversible in 2 labs); one pot-independent first-order MFT decay; O2 a fixed charge. **Fix:** declare `k_dimer_release` from Guth (100 C); a Cys-suppressed loss law | 2x3x3 = 18 | lo | 12 | V | band-now for release; A3 for the loss law and O2 |
| 14 | E15 (F17-a) | **No Strecker acid / Amadori O2 route to Strecker products** (x9 in air, one lab) | 1x3x3 = 9 | H | 9 | V | A3 (O2 axis) |
| 15 | E8 (F25-a, -b, F06-b) | **The FFT / furfural supply pH law is 10 - 70x too steep** (pH 4.5 - 6.5), and the MFT / FFT ratio moves the wrong way. The 1,2-enolisation `acid` factor = [H+] is a full decade per unit, while the hexose k4 RISES with pH. **Fix:** an amine-catalysed 1,2-enolisation term that is pH-flat or weakly rising, beside a weaker acid limb | 2x3x2 = 12 (2 labs within-study; the probe gap is gross but the cells are F25) | lo | 8 | V *p* | A3 (new exponent knob on a factor kind) |
| 16 | E16 (F13-a) | **Norfuraneol decay lacks its base / O2 catalysis** (700x slow at pH 8) | 2x2x2 = 8 | lo | 5 | V | band-now (NF_HYDROXIDE_DECAY exists); A3 for O2 |
| 17 | E17 (F03-a) | **Pentose ARP decays return no amine** (25 - 40 % measured; hexose 71 - 96 %) | 1x2x2 = 4 | H | 4 | V (code) | band-now (amine-return share on the ARP steps as a declared band) |
| 18 | E18 (F26-a, F27-a/b) | **3-mercapto-2-pentanone only from thiamine;** no pH on thiamine degradation; no thiamine sulfur below pH 6.25 | 2x2x2 = 8 | lo | 5 | V | A3 (new route / pH factor) |
| 19 | E19 (F22-a) | **Cysteine-Strecker H2S off** (4.6 - 47x thermolysis, one lab, 1965) | 1x3x2 = 6 | lo | 4 | V (code) | band-now (STRECKER_SULFUR exists) |
| 20 | E20 (F18-a, F21 pyridine) | **No lipid Strecker** (single lab) and **no dienal + NH3 -> 2-pentylpyridine** (3 labs; a plant-protein off-note) | 2x2x2 = 8 | lo | 5 | V | A3 (lane coupling) |

**Lower-ranked leads** (not re-derived):
- formic acid from 3-DG ~8x high (F11-b);
- `r_aa_frag` vs the acetic plateau (CF-069);
- `r_dpo_c2c3` on a refuted analogue: regime mismatch, the refutation is at 37 - 50 C (F10-a);
- `r_g_go` barrier ~0;
- H2S above 110 C ~20x high (F22-b);
- melanoidin C/N fixed (F31-a);
- glycation without pH / a_w (F05-a);
- no Namiki / pentose glyoxal (F09-a).

### 3.1 What would settle each unsettled one

| id | designed experiment or author datum | decides |
|---|---|---|
| E1 | 1:1 [13C6] / [12C6]-glucose + Cys, 120 C, pH 5 and 6.5, time-resolved, with a ribose arm and fed 1-DG + H2S. Or the BTBU / CAAS authors' MFT yields per glucose consumed | the hexose C5 rate and its share vs pentose (sets the band) |
| E2 | o-phenylenediamine in-situ trap on a PENTOSE amino-acid ARP (xylose-Gly) at pH 4.5 / 6 / 7.5 in water vs 0.1 M phosphate, 100 - 120 C | the pentose split law and its buffer dependence (no in-situ pentose trap exists) |
| E3 | fed MGO / 3-DG / glucosone + Leu, Ile, Phe at 100 - 120 C, pH 5 - 7, pH-stat, argon vs air, with SIDA aldehyde AND the aminoketone in one pot. Or DFA's raw Table 4 / Fig. 4 numbers and Quest's SDE recovery (CF-073) | the donor set, the step pH law, the 20x magnitude conflict |
| E4 | same pot as E3, adding 2,5-DMP, pH-stat 4 / 5 / 6 / 7, dilute (Ghent was ~0.85 M) | whether the below-7 slope belongs to the step or to the supply |
| E5 | xylose + Gly and glucose + Gly at fixed in-situ pH 6.0, phosphate 0 / 0.05 / 0.2 M vs citrate, measuring Amadori decay, 1-DG / 3-DG trap, MFT (+Cys) | the buffer factor's form and saturation (CF-011, designed experiment E4) |
| E6 | the Cerny pot D design at 120 and 145 C (fed unlabelled NF vs [13C5]ribose, equimolar), plus 2-mercapto-3-pentanone. Or Cerny's raw 145 C isotopologue data | the free-NF share at process T (unknown above 95 C) |
| E7 | Cys + ribose, pH 5.5, 120 C: buffer vs buffer + fresh oil vs + oxidised oil (PV series) vs + hexanal / 2,4-decadienal dosed, with MFT, FFT, H2S and the aldehydes | chemical capture vs physical partition (Yiltirak is the only separation, one pH) |
| E8 | ribose + Cys, pH-stat 4.5 / 5.5 / 6.5 at 140 C, time-resolved furfural, FFT, furfuryl alcohol, 3-DP (OPD), MFT | the supply law vs the step law (Meynier is one time point) |
| E9 | [13C6]-glucose vs fructose + [15N]-Asn, competition pots, phosphate vs citrate at equal in-situ pH | carbonyl shares and the buffer factor (no labelling on disk) |
| E10 | dilute (10 mM) Cys + xylose formation time course, pH-stat 4.5 / 5.5 / 6.5, phosphate 0 vs 0.2 M, with TTCA + free Cys (CF-076 / CF-011) | K(pH), k_r, catalysis |
| E11 | one pentose + Cys pot at 100 / 120 / 140 / 160 C with the route markers (2-mercapto-3-pentanone, 1:1 label) | per-route barriers |
| E12 | labelled MFT spiked mid-cook, Cys 0 - 150 mM, argon vs air, 110 - 140 C, with the disulfide measured (cluster 4 E3) | the loss law, O2 and reversibility at process T |
| E13 - E14 | fed GO + MGO + Ala dilute, pH-stat, 100 - 130 C, with PZ / MPZ / DMP and the pyrrole. Or Zhou 2024 SI ([M-H]+ correction) | the cross-pair R, the extra MPZ route, the pyrrole sink at dilute concentration |
| E15 | DFA's fed ARP-Phe pot in a second lab, at 120 C, argon / air / chelator | the O2 route's existence above 100 C |

---

## 4. Overall answer

**What we know.**

The early chemistry and its condition dependence are settled in the textbook sense, in dilute phosphate-buffered water at pH 5 - 8 and 90 - 130 C, mostly for glucose + glycine:
- the Amadori compound forms first and, as it decomposes, returns most of its amine;
- phosphate accelerates the early steps;
- the 1-deoxyosone route takes over from the 3-deoxyosone route as pH rises past ~5.5 and carries most of the acetic acid;
- retro-aldol cleavage is negligible at body temperature.

The flavour core is settled by labelling in pentose / glucose + cysteine pots at 95 - 168 C, pH 5 - 6.5:
- MFT is an intact C5 from either sugar, almost never C2 + C3 and rarely free norfuraneol;
- FFT is furfural + H2S;
- thiamine dominates MFT whenever present at low loading;
- pyrazines are aminoketone pairs, plus aldehyde adducts in roasts;
- lipid suppresses meaty thiols several-fold to two orders.

**How sure.**
- Only 25 of 2048 family x pH x T x water regimes (1.2 %) meet the two-independent-labs rule; 233 (11.4 %) have single-lab or conflicting evidence; 1790 (87 %) were never studied (recount on the fixed atlas; the original run gave 31 / 217 / 1800).
- The atlas was ~9 % wrong cell-for-cell (95 % CI 4 - 18 %); the sampled errors are fixed, the unsampled cells still carry that rate.
- So I hold the pooled claims P1 - P22 at roughly 80 - 95 % each. Lower for P7's magnitude (CF-073) and P19's mechanism.
- Any single regime cell outside them deserves at most ~60 %.

**What remains unknown.**
- intermediate water activity (the meat-analogue range);
- in-situ pH during heating;
- any O2 effect above 100 C;
- alkaline cooking;
- plant-protein-bound amines;
- pentose formation splits;
- the 100 - 140 C bin with 1:1 labelling;
- every acrylamide atom map;
- whether lipid's thiol suppression is chemistry or partition.

**Against that literature, the engine's largest errors are structural, not parametric:**
- it has no hexose route to MFT;
- it cannot move the hexose enolisation split with pH;
- it ships no Strecker aldehydes and puts a pyrazine-cascade pH law on the Strecker step;
- it has no buffer catalysis;
- it routes MFT through the wrong C5 intermediate;
- it treats lipid as inert to Maillard odorants.

Seven of the top ten can be started now with existing switches, parameters or declared bands; two of the seven (E10, E14) only in part. The other three (Strecker lane, buffer catalysis, lipid coupling) need A3.

---

## Appendix A. Verification probe (read-only, not a repo artifact)

**Setup:**
- Script: `scratchpad/evl/probe.py` and `probe2.py` (session scratchpad), run through `./scripts/docker_maillard.sh run`.
- Parameters: `engine.core_parameters("sulfur")`, i.e. the shipped b60 frozen report.
- Run: `sulfur.sulfur_flux_budget(..., buffer_spec=BufferSpec(kind="clamped"))`, giving time-integrated fluxes in mM.
- Charges: 30 mM each unless stated.

| probe | conditions | result |
|---|---|---|
| A | ribose + Cys, 95 C, 240 min, pH 5 / 6 | NF share of MFT 0.999 / 0.992; DDP 0.001 / 0.008; C2 + C3 0 |
| B | ribose + Cys, 145 C, 20 min, pH 3 / 4.5 / 5 / 6 / 7 | FUR made 0.105 / 0.096 / 0.055 / 1.2e-3 / 1.5e-5 mM; FFT made 4.0e-3 / 2.5e-2 / 2.3e-2 / 5.7e-4 / 6.9e-6 mM; HS- share of FFT 0.30 -> 1.00 |
| C | glucose + Cys, 145 C, 20 min, pH 6 | MFT made 1.5e-14 mM, 100 % C2 + C3; `r_odg_dpo` flux 0 |
| D | ribose + Cys, 145 C, 20 min, pH 4.5 / 5 / 6 / 7 | MFT made 9.8e-3 / 1.8e-2 / 1.5e-3 / 2.4e-5 mM (NF share 1.00 -> 0.89); MFT / FFT 0.40 / 0.80 / 2.6 / 3.6 |
| E | sugar 100 + Cys 33 mM, pH 5, 145 C, 20 min (Hofmann 1998 T1 design) | MFT ribose 4.2e-3, glucose 9.5e-14 mM (ratio 4.4e10; measured 10.4); FFT ribose 3.5e-3, glucose 7.6e-3 (ratio 0.47; measured 4.3) |

Constants read at 100 / 145 C:
- `k_nf_decay` 4.5e-5 / 6.6e-3 /min;
- `k_mft_decay` 0.041 / 1.36;
- `k_ttca_cys` 0.0286 (100 C);
- `k_dimer_release` 0.

**Caveats:**
- Clamped pH is a pH-stat. The literature pots are buffered with drift <= 0.2 (Meynier) or more (others).
- 30 mM charges are not the papers' charges, except E. The probe tests direction and order of magnitude, not a scored fit.

## 5. Red-team amendment (2026-11-01, night; binding over sec. 3 where they differ)

An independent red team re-ran the probes (they reproduce Appendix A) and tested the top 8 items: the engine claim, the literature strength, the "minimum structure" and the band-now claim. Notes: scratchpad/redteam_engine_errors.md (to be archived with EX-056).

| item | verdict | P (error real and material) | what changes |
|---|---|---|---|
| E1 hexose intact-C5 MFT | engine error SUPPORTED; band-now fix REFUTED | 0.85 | Cysteine is not an amine in the trunk, so Glc + Cys has only amine-free 1-DG supply. HEXOSE_C5_ENTRY on gives ribose / glucose MFT 45 - 78 (measured 10.4), and 0.79 with 33 mM Gly (13x over). Needs a NEW Glc + Cys -> 1-DG step; hexose MFT also inherits E6. |
| E2 enolisation split frozen in pH | WEAKENED | 0.55 | The trapping data are one school; Wageningen's slopes overlap; Martins 2003 at 120 C has equal slopes. No per-step hook exists (a code change). |
| E3 no Strecker aldehydes | a REFUSAL, not a silent error | 0.3 (as a silent error) | Leucine requests are refused. The 3-DG donor band is DFA-only (Quest 0.03x). Stays A3. |
| E4 pyrazine pH law on the Strecker step | SUPPORTED | 0.7 | The supply already carries 0.69 dec / pH, so the step's slope double-counts (shipped methylpyrazine pH 7 / 5 = 170 vs Leahy 14). With the step slope at 0.07, methylpyrazine matches (13.8) but pyrazine goes flat (0.71). |
| E5 no buffer catalysis | SUPPORTED, step-specific | 0.7 | 8+ labs, but FFT and acetic yield are unchanged, and citrate beats phosphate for acrylamide at pH 6. A3. |
| E6 MFT via free norfuraneol | strongly SUPPORTED at <= 100 C | 0.85 | Off on three labs' designs: Cerny pot D NF share 0.98 - 1.00 vs 0.07; Givaudan ribose / NF MFT 0.0066 vs 48; Hofmann 145 C 21 vs 2.6. But F23#13 puts NF at 30 - 100 % at 145 C, and the single 64 kJ/mol barrier makes the share T-invariant. So E6 needs E11 (per-route barriers); k_ddp_mft and k_dpo_nf are frozen keys. |
| E7 lipid inert | SUPPORTED | 0.7 | Reading counts as one lab (+ Huazhong); Du 2025 shows MFT +19 - 30 % with free fatty acids; part is oil partition. A3. |
| E10 TTCA formation off | engine error SUPPORTED; fix REFUTED | 0.75 | B58 already ran TTCA_FORMATION (DO NOT SHIP). The probe moves MFT / FFT <= 4 %. Slowing k_r contradicts Cao's K under detailed balance (k_r >= ~7e-3 /min at 90 C). The real lever is B59b CYS_AMADORI_ROUTE (MFT x0.16 - 3.6, FFT up to x654 at pH 6.5). |

**Recommended order of change.** Each step is a pre-registered declared-band change behind an EXISTING switch or parameter, judged by its evidence test.
1. **E4:** band the below-7 pyrazine Strecker-step slope at [0, 0.15].
   - Test: Leahy's methylpyrazine ratio within 0.3 dex; pyrazine's expected miss pre-registered.
   - Shipped pyrazines move x2.6 - 10 at pH 5 - 6: OWNER APPROVAL.
2. **E6 + E11 as one refit:** unfreeze k_ddp_mft / k_dpo_nf, with per-route barriers via formation_Ea_by_route bounded by measured classes.
   - Tests: pot D NF share <= 0.15; Givaudan >= 10; Hofmann within 3x; no B60 row worse by > 0.3 dex.
3. **E10 via B59b,** after fixing its dpo_share split law. Test: TTCA (Zhai 2019, Zhang 2025b) within 2x.
4. **E1:** a new hexose + Cys -> 1-DG step. Tests: Hofmann T1 MFT 10.4 and FFT 4.3 within 3x.
5. **E2:** per-step exponents, with the 120 C equality declared.
6. **E5, E7, E3:** after A3.

## 6. Recount on the fixed atlas (2026-11-01)

Sec. 5 stays binding. Nothing in secs. 2, 3 or 5 is re-ranked or rewritten; this section only flags entries whose supporting cells the fix-up changed.

**Inputs that changed.** I compared the analysis's overridden view of each cell with the fixed file (scratchpad `evl2/diff.py`):
- **11 cells** back or new: the nine corrected red-team cells, the split F25#22 and the new F25#23. F24#13 was deleted.
- **11 status changes** against the override:
  - eight non-detects now `absent` with a detection basis: F20#3, F23#9, F24#0, F24#1, F24#5, F26#1, F27#10, F27#11;
  - F10#1 and F13#10 became `minor`;
  - F12#5 became `not_detected`.
- **58 water moves**, almost all `dilute` -> `concentrated` (0.5 M buffer).
- **17 temperature moves** against the midpoint rule.

**Mismatches whose supporting cells changed:**

| id (rank) | cells changed | classification | confidence / strength |
|---|---|---|---|
| F23-a / F24-a (E1) | F24#0, #1, #5 `not_detected` -> `absent` (label pattern); F23#8, F24#5 -> `concentrated` | unchanged (M1 + M2, V *p*) | unchanged H; the M2 side is firmer: C2 + C3 is now `absent` on a label basis in all three labelling labs. Sec. 5's E1 verdict is unaffected |
| F23-b (E6), F23-d (E11) | F23#9 -> `absent` (positive control: the same analysis saw the NF marker in fed-NF pot D) | unchanged | unchanged H; firmer at 95 C. Sec. 5's caveat (F23#13: NF 30 - 100 % at 145 C) stands |
| F10-a | F10#0 back as `operating` (MGO "negligible", seen); F10#1 -> `minor` (<= ~1.3 - 1.7 mol %); F10#3 `absent` (positive control, fed 3-DG) | unchanged (M2, L: a 37 - 50 C vs 100 - 145 C regime mismatch) | unchanged; "negligible" now reads "bounded <= ~1.7 mol %" |
| F25-a (E8) | F25#7 back: Cerny 2007's within-study ladder, FFT x2.3 down from pH 4 to 6, and a standing-level zero at pH 7 (`not_detected`) | unchanged (M3 magnitude, V *p*) | stays lo under the F25 rule. A third lab (Firmenich) now agrees with DFA and Reading, so E could read 3, not 2 (E x O x R 18, effective 12 instead of 8). Not re-ranked: these are standing levels, so they bound net accumulation, not formation |
| F25-c | F25#23 new (Zhao 2019b, BTBU CAMOLA): hexose FFT ~94 % an intact glucose C5, `unattributed` | unchanged | unchanged lo. A second hexose intact-C5 datum beside F23#7, so it bears on E1, but from the same lab |
| F20-a (E14) | F20#3 `not_detected` -> `absent` (label): glycine + HCHO addition gives MPZ <= ~1.5 umol/kg at 180 C, while 2,6-DMP takes ~45 % from labelled HCHO | unchanged (M1, V) | **scope narrowed:** the missing addition route is shown for DMP / ethyl pyrazines, not for MPZ |
| F27-b (E18) | F27#10 -> `absent` (methylene-blue assay): thiamine H2S vanishes below pH 6 at ~100 C. With F27#7 (non-thiazole S <= 4 % at pH <= 6, 121 C, 35S) | **narrowed:** M1 only at >= 145 C. At 100 - 121 C the literature SUPPORTS the engine's zero thiamine-H2S share below pH 6.25. What survives is that the cut-off ignores T (F27#1, Cerny 145 C, a single unreplicated control) | lo, and weaker. P(material) lower than the original listing implied |
| F27-a (E18) | F27#11 -> `absent`: 70x cysteine stops thiamine loss at 100 C (CF-089 vs Grosch's x7 thiamine-MFT) | unchanged (M4) | unchanged lo; adds a cysteine dependence the engine does not express (lead; CF-089 open) |
| F26-a | F26#1 -> `absent` (the NF-marker isomer is not made from ribose; positive control); F26#0 - 5 -> `concentrated` | unchanged | unchanged |
| F12-a | F12#5 `absent` -> `not_detected` (erythronic acid, 37 C) | unchanged (the M1 rests on the 18O2 and C3-acid cells) | unchanged |
| F31-a | F31#3 back as `operating` (C/N 7.64) | unchanged (inside 6.6 - 21) | unchanged |
| F04-a, F08-a | F04#11, F08#5 temperature moves | unchanged: the settled bins they cite (F04 neutral / high, F08 neutral / low) are still settled | unchanged |
| F16-a, F17-a, F19-b, F29-a / -c, F30-a, F22-b | water or temperature moves only (F16#2; F17#2, #5 - 7; F19#22; F29#9, #18; F30#2, #6; F22#4 - 6) | unchanged | unchanged |

No mismatch cites the other cells that are back: F02#12, F07#14, F23#15, F25#6 / #22, F29#4 and F30#15.

**Net:** the fix-up changes no mismatch's kind, and no H / lo confidence label.
- Two entries narrow: F27-b to >= 145 C, and F20-a to the non-MPZ adducts.
- One gains evidence that is not used in the ranking: E8, now three labs within-study.
- E1's M2 side and E6 are firmer.

### 5.1 Correction to sec. 5's E4 row (2026-11-01, ~03:40; verified by the coordinator against leahy1989a_extraction.md Table I)

The red team's "shipped methylpyrazine pH 7 / 5 = 170 vs Leahy's 14" misread the table. **14.4 is PYRAZINE's** k(pH 7) / k(pH 5) at 95 C. **Methylpyrazine's is 215** (99 at 85 C). So the shipped engine (170) MATCHES Leahy's methylpyrazine pH dependence, and E4's band [0, 0.15] would make it ~1.2 dex too flat (9.9 - 20).

What remains of E4:
- the Strecker STEP's law fails the fed-dicarbonyl check (Hofmann: PA largest near pH 5; shipped k5 / k7 0.069);
- the same slope cannot serve both Leahy's sugar + lysine cascade and fed-dicarbonyl Strecker.

E4 is therefore NOT a band-now change. Its prereg (kinetic_core_pyrazine_strecker_step_ph_prereg.md) expects REFUSE (P 0.85): a declared band, plus a designed experiment that separates the step law from the supply law (sec. 3.1, E4 row). The recommended order becomes CF-077 first, then A3 phase 2.

Lesson: a red team's numeric claim is a lead too; verify it against the dossier table before it orders work.

## 7. Amendment: batch-2 conflict resolutions (2026-11-01, afternoon; binding over secs. 1.3, 3, 3.1 and 5 where they differ)

Sources: resolve2_pyrbuf (buffer form, Piloty, pyrazine rivals, early steps), resolve2_sulfur (CF-076, van den Ouweland,
Wei, CF-031) and the Wolf 1977a reading, each with an independent red team whose REQUIRED changes are applied here
(records: `docs/history/conflict_resolution_2026-11-01/resolve2_pyrbuf/`, `resolve2_sulfur/`, `wolf1977a/`; register
CF-092 - CF-100, notes on CF-011 / 031 / 072 / 076 / 087; EX-058). Probabilities are the red teams'. Nothing here moves
the engine; every numeric band below is A3-provisional until its own prereg and red team. Earlier sections are not
rewritten (the frozen splits read them at their freeze commits). Sec. 7 is not a cited section of the split builder (secs. 2 - 5 and Appendix A are), so naming a TEST study here does not make it development. The machine-readable twin
`engine_vs_literature.yml` is not regenerated by this amendment.

**Hold-out status, corrected.** Potman, Kallen and Thornalley are E5 TEST in split v2. Larralde 1995 and Lee 2024 are
NOT in any frozen split: they are keyed held-out tables (contrast-benchmark README sec. 10), candidate TEST for a v3.
Every "TEST" label on them in the resolution records reads "keyed held-out". UC Davis (`ucdavis_shibamoto`) is an E13
TEST laboratory; the resolution read Shibamoto 1976, Wong 1988 and Bemis-Young 1993 (all UC Davis), so an E13 score
reports with and without that lab. The resolution also viewed `heldout_zhang2021_table3.yml` (E7 keyed held-out) and
F25#17@us8920862 (E7 TEST) before proposing an E7 structure.

| error | what changed | P (red team) |
|---|---|---|
| **E4** (Strecker-step pH law) | Piloty 1979 is NOT discriminating: the engine's Strecker step is 1e-5 - 1e-3 of Piloty's measured consumption under both rivals; pH-flat consumption allows a Strecker share <= 0.04 (Ala). The piloty1979 dossier's "SUPPORTS E4 (P ~0.6)" is withdrawn. Piloty instead exposes a missing dicarbonyl-amine adduct sink, >= ~1e3x the engine's Strecker consumption, a lower bound (Sedehi's MGO + Gly gap ~2e4); not a k_mgo_ak error (CF-095). Sec. 3.1 E4 row: the designed pot must measure the Strecker aldehyde (and acid) and the aminoketone, not amino-acid loss. | not discriminating 0.85 |
| **E5** (buffer catalysis) | The fix's "saturating" is replaced by: **non-saturating to 0.5 M in dilute pots (order 1.0 - 1.35)**; K <= ~0.3 M rejected on Potman (power >= 0.85 at run-to-run noise; K 0.5 M sits on the threshold). Kallen's ceiling (x3.1 - 4.4 on the re-read) is consistent with a change of rate-determining step (+0.09 / -0.11 / -0.22 dex), but a saturating law fits it equally well. A per-step factor k_s = k_s,u (1 + sum_b kappa_{s,b} c_b) holds only on steps the network RESOLVES; a lumped step whose uncatalysed part can limit (r_schiff) carries a hyperbolic factor; concentrated foods carry a ceiling (Potman Series A, 1.2 M phosphate at 30 - 40 % dry matter: 0.43 dex short of linear). Phosphate kappa [5, 40] /M at pH 5 - 7.4 (DEV: Guan, Davidek, Zyzak, Bell) is NOT temperature-flat (Guan 110 -> 130 C falls 0.4 - 0.75 dex) and is window-dependent in Potman (5-h secant 22 - 29, 3-h 48 - 81 /M). Carboxylates on Amadori steps [0, ~4] /M (Guan's citrate is 3.4): <= 0.1x phosphate at 110 C, <= 0.3x at 130 C; free-sugar steps rho x kappa_P, rho [0.1, 1]; Larralde's bicarbonate consistent with zero. The pH-8 ratio fade is ~60 % k_u rising, ~40 % lost catalysis (Davidek). Browning is not the Amadori step: browning factors (x2 - 13) do not transfer to single steps. Kallen's no-catalysis-below-pH-5 test used carboxylates, so it does not place Cao's catalysis on one step (CF-011). Sec. 3.1 E5 row: add a citrate arm at equal in-situ pH (decider: citrate ~ water on Amadori decay = step-specific; citrate 0.3 - 1x phosphate = general); "saturation" leaves the question for dilute pots, the concentrated-matrix ceiling enters it. CF-093, CF-094. | non-saturating (K <= 0.3 M, dilute) 0.8; "linear, no Michaelis term" as an engine law for lumped steps / foods 0.45; step-specific 0.5; browning is not the Amadori step 0.7 |
| **E6** (MFT via free NF) | Unchanged in substance. van den Ouweland & Peer 1975 propose the NF -> MFT pathway (Fig. 4) but identify no MFT (CF-092); the route's existence rests on four other labs (P 0.97) and E6 (the NF SHARE) never rested on the 1975 paper. Wei 2022's TTCA is not a deoxyosone-trap product: TTCA is C8H15NO6S (M 253); 1-DX cannot give its ring, and either DX can only by rehydration back to the pentose level (redox-neutral); TTCA mainly via ARP reversion / pentose 0.6, direct DX trap 0.15, m/z 254 not (only) TTCA 0.25. | NF route exists without the 1975 paper 0.97; peak 3 is MFT 0.3 |
| **E7** (lipid inert) | CF-031: Zhang 2022's free-MFT rise is a water-phase speciation shift at best: the pre-registered "within x3" rule FAILS on its ENTH arms, and the oil phase holds both MFT and the disulfide, so partition is not excluded. "No partition term" becomes "not needed for the sign data on disk" (cluster 4's <= 0.1 - 0.4 dex). Zhang 2021 (E7 keyed held-out) is removed from the evidence; Xu 2011's lard contrast is E12 TEST F28#20 and must not calibrate an E7 term. Before any E7 engine work: a pre-registered reject / keep rule on the UNSEEN E7 TEST contrasts (F30#11 kim1996b, F30#15 zamora2015, F30#16 hidalgo1993), stating whether each was viewed; F25#17 (viewed, wrong sign) is a carried failure. Sec. 3.1 E7 row: measure free + disulfide (or TCEP-reduced) MFT in both phases. | speciation shift explains the sign 0.4; disulfide partition + more formation with oil ~0.35; other 0.25; an E7 capture term on >= 2 DEV labs, one channel, no TEST cells ~0.3 |
| **E10** (TTCA exchange) | CF-076's in-situ-pH admission is CONDITIONAL on a low release Ea along a (k25, Ea) ridge: Karpinski's printed rate needs ~26 kJ/mol, mid-band ~32, the x6-deflated floor ~41 (E4c) to ~49 (E2c); the x6 deflation is not self-consistent with E2c and exceeds E4c's own phosphate factor (~1.6). "Rejected at Ea >= 50" holds for E4c only (E2c feasible at 48.5). The fitted pKa_T is a T-independent effective value, not comparable with Radomski's room-temperature inference; the pKa is not shown to be non-binding. Sec. 3's "slow k_r to 2 - 8e-4 /min" becomes the run's feasible 0.9 - 3.2e-4 /min, on the ridge only. Decider (sec. 3.1 E10 row, first): one Arrhenius series of cysteine release from TTCA (or AThz) at fixed pH 4.5 and 6, 25 - 90 C, pH recorded at temperature; the pKa titration second. | in-situ pH explanation 0.25 (was 0.45); K-released 0.3; second sink / unlisted 0.3; effective release Ea <= 32: 0.1, <= 45: 0.25, <= 50: 0.35 |
| **E13 / E14 and the pyrazine rivals** | P0 (one constant, R = 4) rejected by **3** independent laboratories (UC Davis, Koehler, Lee), not 5. The P2b rejection (a lower-order extra MPZ route, "pool invariance, P 0.75") is WITHDRAWN: cumulative yields are not pool sizes; with rate-based pool proxies Shibamoto's pools move only 1.1 - 1.4x and P2b gives |z| <= 1.8; after 2 h every species is in net loss. Lee 2024 at face value (R 15 - 26) violates only P0 (under P1x, R = 4 r^2 has no upper bound); only MPZ / 10 (R 0.15 - 0.26) breaks CF-087's floor, so the x10 decides whether any bound beyond P0 fails (CF-097). Low R with glycine is reproduced outside the Ho network (Bemis-Young, UC Davis) and is a condition effect (CF-087: P(low-level data) 0.4 -> 0.2). Same-lab R gap 14.6 vs 357 - 631 (Wong vs Shibamoto 1976 / 1977; CF-098). Survivors: condition-dependent pair constants (P1x / P1y) or a second-order extra route (kinetically identical to P1x), separable only by labels (Zhou 2024 raw ions; EX-052 open). The "~160 kJ/mol 25 - 60 C" constraint rests on 30-day pots at 0.005 - 0.007 % yield (CF-087's low-level caveat). The prereg named four E13 TEST studies (shibamoto1977, low2007, vanlancker2012, amranihemaimi1995; 8 contrasts in split v2); none was scored or excused, against CF-087's no-silent-omission rule: owed before any E13 verdict. For the engine: no structure change; R carries a declared condition band. | P0 rejected (3 labs) upheld; P2b back to ~0.45 |

**Atlas summary P2 (sec. 1.3), restated.** Phosphate accelerates early condensation, Amadori decomposition and
free-sugar degradation (x1.3 - 12). The dose law is non-saturating to 0.5 M in dilute pots (order 1.0 - 1.35); Kallen's
thiazolidine plateau is consistent with a change of rate-determining step; concentrated foods show a ceiling. The
catalysed increment peaks broadly at pH 6 - 7, and the RATIO falls at pH 8, mostly because the uncatalysed rate rises.
The sec. 1.2 F03 row's "fades at pH 8" reads "ratio falls at pH 8 (uncatalysed rate rises)".

**Early steps.** Huber & Ledl 1990: the 1-amino-1,4-dideoxy-2,3-hexodiulose is a C-4 elimination product at the
1-deoxyosone's oxidation level, so the McGill ~5 % share belongs to F07, not F15 (P 0.85; CF-099). The cell is F15#3, a
v1 TEST contrast, so the move waits for a new split version; a Huber cell is appended to F07 after the frozen extent.
a_w: Cui 2026's maxima are 20-min end points of consumed intermediates (no conflict, P 0.85); protein glycation has no
a_w term (F05 (d)), a gap for designed data. Wolf's lysine plateau (CF-096) and the glycation barrier (Wolf 1977a: E
137.7 kJ/mol vs the engine's 96.8, ~5 SE; CF-100) are open; the glucose order is not part of CF-100 (the text's x5.3 /
x25 bracket order 1).
