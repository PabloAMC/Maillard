# Case M1b: pH in the sulfur lane's sugar steps

**Rescoped 2026-11-04.** This case drafted the plan's original M1b claim. Its evidence led the plan to rescope M1b to "the pH balance (CF-084) and the pentose-side steps"; `ph.trunk_terms` alone is not built; see STATUS_AND_PLAN sec. 5, W2. The claim below is kept as drafted.

*Plan step M1b ([STATUS_AND_PLAN sec. 5, W2](../../guides/STATUS_AND_PLAN.md); errors E2 (F06-a, F06-c), E16 (F13-a),
with E8 (F06-b, F25) as context). Items 1 - 5 of sec. 7, written 2026-11-04 before any result of the change exists.
No engine was run for this case. Items 6 - 9 are appended later and never edited into this part.*

**Prompted by a held-out score, declared.** The step entered the plan on 2026-11-04 from the W1.1 contrast score
(`results/validation/kinetic_core_contrast_benchmark_score_v1.md`): on Meynier & Mottram 1995's glycine + ribose pH
ladder (`data/lit/fit_rows/heldout_meynier1995_tab1_2_4.yml`, Table 1; 0.2 M pyrophosphate, 140 C, 1 h) the shipped
engine predicts the same norfuraneol and furfural at every pH from 4.5 to 6.5 (engine log-ratio 0.00 on every pair),
while measured norfuraneol falls with pH (m = log10(NF at the lower pH / NF at the higher pH) from -0.02 to +1.21 dex;
+1.21 for pH 4.5 against 6.5). **This case reads Meynier 1995 only in sec. 7 item 7 (side effects) and never uses it
to choose between rivals.** See sec. 2.8 for what the frozen splits say about that paper.

**Seen by the author.** The W1.1 score file's primary table and its Meynier items (NF and FUR, glycine and
methionine; measured m and engine p), its Grosch 1992 Table 2 items, and the score prereg v6's exposure section; the
header and first two pots of the Meynier keyed table (conditions and nd / tr lists). No other held-out engine result
bears on this case. The E2 and E16 v2 TEST sets are empty, so no v2 TEST result exists to have seen.

## 1. The claim

The sugar-side steps of the sulfur lane move with pH as the measurements of those steps say, so that the pH the pot
drifts to (case M1, part D) reaches norfuraneol and furfural. Today the furfural and norfuraneol a pentose pot makes
are pH-free in glycine pots, and carry a full decade per pH unit through cysteine in cysteine pots.

**What the existing switch touches, read from the code (no run).** `ph.trunk_terms` (`sulfur.TRUNK_PH_TERMS`,
[sulfur.py](../../../src/kinetic_core/sulfur.py) ~1278 - 1295, ~1322) multiplies only the keys in
`trunk_conditions.PH_STEPS` = (`k_ama_tdg`, `k_ama_odg`, `k_ama_mgo`), the three HEXOSE Amadori-decay steps, plus the
shipped 3-DG -> formic exit `k_tdg_fa`, all with the pooled 0.69 dec / pH unless `ph.trunk_exponent_override` gives
per-key slopes. In a pentose-only pot none of these carries flux. Meynier's flat engine answer has two other causes:
`PENTOSE_AMINE_ENTRY` is off in the shipped engine (glycine does not react with ribose), and the amine-free pentose
enolisations `r_pent_caramel` and `r_pent_thermal_dpo` carry no pH factor; `r_nf_decay` is pH-free and
`NF_HYDROXIDE_DECAY` is off.

So the claim has two parts, and the build must say which it ships:

| part | steps | today | M1b |
|---|---|---|---|
| H. hexose trunk (the plan's wording) | `k_ama_tdg`, `k_ama_odg`, `k_ama_mgo` (+ `k_tdg_fa`) | no pH in the sulfur lane | Martins 2003 Part II Table 3 per-step slopes (100 C: 0.37 / 0.92 / 0.72 dec / pH; 120 C: 0.69 / 0.66 / 0.76) through `ph.trunk_exponent_override`, with the 120 C equality declared |
| P. pentose side (where the prompting observation sits) | `r_pent_caramel`, `r_pent_thermal_dpo`; the pentose + amine entry `r_pent_amine` -> `r_pama_dpo` / `r_pama_tdp`; `r_nf_decay` / `r_nf_decay_oh` | pH-free; entry off; hydroxide decay off | a pH law on the steps the chosen rival names (sec. 3), each with a declared source |

## 2. All the evidence, by kind

Split notation as in case M1. For this step the v2 errors are E2 (DEV labs Ledl / Severin school, Nestle / Cerny,
Wageningen; **no TEST contrast**: UNTESTABLE_on_held_out; Zhai 2021's F06#7 and F07#12 excluded as calibration
input), E16 (DEV labs Ochanomizu, Reading; **no TEST contrast**) and E8 (DEV labs DFA, Reading, Wageningen; no TEST
contrast). In v1 every cell below is DEVELOPMENT unless marked.

### 2.1 Isotope labelling

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Davidek 2006a (F07#7) | acetic acid ~70 % from glucose C1 / C2 (the 1-DG cleavage); "pH had only a negligible effect" on the label origin, pH 6 vs 8 | the route of the acid, not the split | DEV |
| Davidek 2008 (F07#8) | xylose + glycine, pH 6 vs 8: acetic 16 vs 28 mol %, >= 12 - 14 % and >= 22 - 25 % of xylose via the 1-deoxypentosone | pentose 2,3-enolisation rises with pH (a lower bound) | DEV; fit rows (EX-013) |
| Cerny 2007b, Zhou 2025b (F13#0, F13#1) | furfural carries the intact C5 (from free sugar or the ARP's own C5) | furfural's parent; no pH | DEV |
| Cerny 2003 / 2004 (F25#1 - #3, F23 pot D) | furfural and FFT intact C5; fed NF gives no FFT | NF and furfural are separate exits; no pH | DEV |

### 2.2 Fed intermediates

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Martins & van Boekel 2003 Part II (martins2003 dossier; `trunk_conditions.PH_SOURCE`) | fed DFG (glucose-glycine Amadori), 100 / 120 C x pH 5.5 / 6.8: Table 3 constants with 95 % HPD; k(6.8) / k(5.5) 3.0 / 15.6 / 8.6 (100 C), 8.0 / 7.3 / 9.8 (120 C) | parts H and rival R4 (hexose) | DEV (Wageningen); cal. input |
| Martins 2003b (F06#5, F07#10) | same fed DFG: standing 3-DG : 1-DG ~12 at pH 5.5, 1.2 (figure) or 3.2 (text) at pH 6.8 | split vs pH, standing levels | DEV |
| Davidek 2006a (F07#5, #6, #9) | fed DFG, 120 C: acetic 4 / 28 / 56 mol % at pH 4 / 6 / 8; fed amine-free 1-DG 0.20 / 0.56 / 0.80 - 0.84 mol / mol at pH 5 / 6 / 6.5 - 8 (60 min, flat above 6.5; corrected 2026-11-05: its grid is not 4 / 6 / 8); 1-DG lifetime < ~2 min at pH >= 6 | the 1-DG cleavage's own pH law (the B7a-2 saturating law) | DEV |
| Mikami 2017 (F15#7, F13#8; B61 prereg) | fed norfuraneol alone, 0.2 M phosphate, 100 C: 0.032 /min at pH 8 (+- 10 %). At pH 7 and 6 only the xylose + lysine pot bounds it (>= 1.8e-3 /min at pH 7; no net loss at pH 6), so those two are not fixed-supply measurements. B61's anchor (`NF_HYDROXIDE_DECAY` is set from this point), so it cannot be decisive for E16 (it bears on D1 below only as the calibration point; corrected 2026-11-05) | rival R3: NF loss at fixed supply | DEV (Ochanomizu); cal. input |
| Nakamura 2020 (F15#8) | fed NF, pH 8, 100 C: decay suppressed by degassing and by DTPA | R3's O2 and metal dependence (couples to case M1 part B) | DEV |
| Deng 2021 (F06#6, F07#11) | fed hexose ARP (Met-glucose), unbuffered, initial pH 5.5 vs 7.5: 3-DG peak ~7 vs <= 1.8 mol %; 1-DG standing lower at 7.5 | Amadori decay and split vs pH, drifting pH | DEV (ho_network) |
| Hou 2026 (F06#8, F07#13; CF-083) | fed GX-ARP (glycine-xylose Amadori), 120 C, unbuffered: standing 3-DX / 1-DX 2.6 -> 10.9 while the measured pH falls 5.05 -> 3.55 | the pentose + glycine analogue closest to the prompting pots | DEV (ho_network) |
| Zhou 2023 ([thiol_zhou2023_ph_series.yml](../../../data/lit/fit_rows/thiol_zhou2023_ph_series.yml)) | fed Ala-xylose ARP +- Cys, initial pH 6 / 7 / 8, 120 C: furfural printed for both arms | furfural vs pH from a pentose ARP | DEV; cal. input |
| Hofmann 1998 fed C5 precursors (F13#2, F25#4) | fed ribose / 3-deoxyribosulose / furfural + H2S, pH 5, 145 C | the furfural exit's size, one pH | DEV |

### 2.3 Trapping

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Beck 1988 (F06#0, #1; F07#0 - #2) | in-situ OPD on a hexose ARP, pH 4.5 vs 7: 3-DG share of trapped C6 osones 38 % -> 5 %; 1-DG : 3-DG x12.5 (~0.44 dec / pH); n = 1 | E2: the flux split moves with pH | DEV (Ledl school) |
| Nedvidek 1992 (F06#2, F07#3, F15#0) | in-situ OPD, glucose + beta-alanine, pH 5 vs 7: 3-DG share 19 -> 4.6 %; temperature never printed; no replicates | E2; the same school as Beck (not independent) | DEV (Ledl school) |
| Nakamura 2020 (F13#9) | in-situ OPD, xylose + Lys, pH 8 only: 1-DX 76 % of xylose, standing 1-DX < 0.2 % of its flux | the pentose 2,3-route dominates at pH 8; no pH contrast | DEV (Ochanomizu) |
| Gobert 2009 (F06#3, F07#4) | in-situ OPD, 50 C: air lowers 1-DG x3.4; 3-DG O2-independent | an O2 term on the 1-deoxyosone, low temperature | DEV (Glomb) |
| none | no in-situ pentose trap across pH (ENGINE_VS_LITERATURE sec. 3.1, E2 row) | gap | -- |

### 2.4 Measured intermediates and lifetimes

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Zhai 2021 (F06#7, F07#12) | TTCA (+- Xyl), unbuffered, label pH 5.5 / 6 / 7 / 8, 100 C: 3-DX / 1-DX at 120 min 19 / 15 / 6.1 / 0.7; ratios robust, levels on a diacetyl surrogate | the pentose standing ratio falls with pH | DEV in v1; excluded from v2 E2 as cal. input |
| Hofmann 1999 (F06#4, F07#20, F13#3) | xylose + Ala, pH 7 reflux: 3-DP peaks at 10 min, furfural 0.06 mol % at 180 min: furfural a small exit of 3-DP | furfural's share of 3-DX, one pH | DEV (DFA) |
| Mikami 2017 (F13#8) | xylose 13.3 + Lys 34 mM, 0.2 M phosphate, 100 C, 60 min: NF 27 / 61 / 100 / 83 / 60 / 21 % of max at pH 5.5 - 8.0; ribose > xylose | NF standing level vs pH: rises to pH 6.5, then falls | DEV |
| Davidek 2006a (F07#9) | 1-DG lifetime < ~2 min at pH >= 6, 120 C | why standing 1-DX is small | DEV |
| Hou 2026 | lifetime bounds 3-DX ~45 - 140 min, 1-DX ~3.5 - 30 min (120 C, dossier-derived) | R5 | DEV |

### 2.5 Mass balance

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Davidek 2006a (F07#5) | >= 67 - 70 % of DFG through 1-DG at pH 8 (dossier-derived) | E2 at pH 8 | DEV |
| Beck 1988 (F07#2) | secondary-amine ARP: 1-DG quinoxaline isolated in 70 % yield even at pH 4.5 | amine structure moves the split | DEV |

### 2.6 One factor varied inside one study

| source (record) | factor -> result | bears on | split |
|---|---|---|---|
| Hofmann 1998 T2 (F25#8; CF-012) | ribose + Cys, 0.5 M phosphate, 145 C, pH 3 / 5 / 7: MFT 55.3 / 19.8 / 2.5, FFT 22.9 / 12.1 / 1.2 ug | NF-side vs furfural-side products vs pH (through H2S) | DEV (DFA); cal. input |
| Cerny 2007 (F13#7, F25#7) | xylose + Cys + thiamin, 145 C, pH 4 - 7: furfural 1.32 / 1.00 / 1.04 / 0 / 0 (rel.), FFT 431 / 368 / 364 / 185 / 0 | furfural standing level vs pH | DEV (Firmenich) |
| Zhou 2023 | initial pH 6 / 7 / 8 (see 2.2) | furfural vs pH | DEV |
| Davidek 2008 (F07#8) | pH 6 vs 8 (see 2.1) | 1-DX route | DEV |
| Mikami 2017 (F13#8) | pH ladder (see 2.4) | NF vs pH | DEV |
| Martins 2005b (F06#13, F07#15) | initial pH 4.8 -> 7.5, glucose + glycine: k1 (glucose + glycine -> DFG) +0.50 +- 0.14, k4 +0.35, k7 +0.60 dec / pH (wide CI on k7). Every Martins slope in this case (secs. 1, 2.2, 2.6) is regressed on the initial label pH of drifting pots; on the time-averaged pH they move (CF-118; corrected 2026-11-05) | E2, hexose | DEV (Wageningen) |
| Davidek 2002 (davidek2002 dossier; sulfur.py comment at `TRUNK_PH_TERMS`) | pH-stat: acid 3 vs 79 mmol/L at pH 5 vs 8 | the acid supply that drives the drift | DEV (Nestle) |
| Blank 1997 (F13#11) | water vs phosphate vs malonate: HDMF x43, x4.3 | buffer (E5), not pH | DEV |
| Rizzi 2010 (F15#2) | amine-free pentoses, phosphate pH 7.1 | the amine-free route's size at one pH | DEV |
| Mottram 2002 (F13#10) | ribose alone vs ribose-5-phosphate: NF trace vs 0.25 mol % | amine-free NF is small | DEV |
| Haleva 1999 (F13#12, #13) | pH 3 vs 5: HMF down, furaneol up (hexoses, 70 C) | opposite pH directions in one study | **v1 TEST** (Hebrew University): not used |

### 2.7 Fitted time courses and fits

| source (record) | what it shows | bears on | split |
|---|---|---|---|
| Martins 2005 (F06#12, F07#16) | Ea k4 97.1, k7 107.3 kJ/mol; initial pH 6.8 drifting to ~5.5 | the split vs T | DEV |
| Desclaux 2006 (F06#14, F07#17) | glycine-xylose network fit, 0.5 M phosphate, pH 6: ARP -> 3-DX 5.50e-4, -> 1-DX 2.79e-2 /min (1 : 51); no CI; furfural and nor-HDF fitted at one pH | the pentose split at pH 6 | DEV (Reading + Quest) |
| EX-013 (pentose + amine entry) | entry identified 10^-4.2 at 100 C / pH 6, slope 0.57; norfuraneol sanity bound FAILED (25 - 36 mol %) | rival R4 already fitted once | engine record |
| EX-001, EX-002, EX-003 (sub-model 1) | the same Martins slopes on the acids: acetic transfers (hold-out 5 / 5), formic does not (0 - 1 / 4); arm C's condensation slope rejected on the frozen hold-out | **against or mixed**: the slopes transfer for one product and not another | engine record |
| B40 / B41 (WAVES.md) | 3-DG exits' pH term: the formic exit's slope shipped (B41); the fragmentation slope rejected (Leitzen MGO 1.3x -> 33x) | **mixed**: a slope from a lumped step does not transfer | engine record |
| EX-033, EX-036, EX-041, EX-043 - EX-046, EX-048, EX-053 (deoxyosone partition, B7a) | the engine's 3-DX / 1-DX rises with pH where Jiangnan's falls (wrong sign); a split pH law (B7S) is identified but leaves the ratio's pH slope ~0 (DPO's sinks pH-dependent the other way; 39 - 45 % of 3-DX from pH-flat caramelisation); a fast base-favoured 3-DX exit is the first structure to pass time and pH jointly (EX-046), but its MGO stand-in overshoots MGO and the cysteine trap version fails (EX-048); the tester's power is limited (EX-053: "not found", not "not feasible") | the strongest DEV evidence on rival R2 vs R5 | engine records |
| ENGINE_VS_LITERATURE sec. 5, E2 row | WEAKENED, P 0.55: the trapping data are one school; Wageningen's slopes overlap; Martins at 120 C has equal slopes | **against** part H's strength | review record |
| B61 (`kinetic_core_b61_prereg.md`) | `NF_HYDROXIDE_DECAY` at Mikami's value: audit "grounded at the measured point (pH 8, 100 C)", barrier a class prior (103.7 kJ/mol); declared tension with Lau 2003 / Pandit 2006 plateaus; fit stopped by the owner, no verdict (`tasks/data_restructure_plan.md` ~990) | rival R3, never scored | engine record |

### 2.8 Meynier 1995: what the splits say

| record | status |
|---|---|
| `heldout_meynier1995_tab1_2_4.yml` (Tables 1, 2, 4) | keyed held-out (not in a frozen split; contrast-benchmark README sec. 11); TEST for E4 / E3 only; scored as TEST by W1.1 under prereg v6, NF "clean", FUR `dev_exposed` |
| F13#6@meynier1995 (NF falls 16.1x Gly, 6.2x Lys, 3.5x Cys, 1.7x Met, pH 4.5 -> 6.5) | **DEVELOPMENT in v1** (Reading is a DEV lab); cited for E16 in ENGINE_VS_LITERATURE sec. 2.1 row F13-a, so Reading is DEV(E16) in v2; also cited in CF-012's cluster-3 test |
| F13#5, F15#5, F25#9, F23#19 @meynier1995 | DEVELOPMENT in v1; F25#9 and the furfural fold already scored since B2 (CF-001) |
| F26#11@meynier1995 | DEVELOPMENT in v1; a **v2 TEST contrast for E18** |

So the prompting observation is, under both frozen splits, development evidence for two of this case's rivals (E16,
and E8 through Reading). The case follows the plan's stricter rule anyway. The cost is real: Meynier is the only
study on disk that prints furfural and norfuraneol across a pH ladder in one pot (sec. 4, D4).

**Sources left out, with reasons.** Low-water and storage cells (F06#15 - #17, F07#18 - #19, v1 TEST; F06#10, F06#11
at 37 C) are outside the regime. The furfural barrier question (B62, E11) is temperature, not pH. Buffer catalysis
(E5) is M5. Dossiers outside the atlas families F06, F07, F13, F15, F25 were not screened one by one: a stated gap.

## 3. The rivals

| id | explanation | furfural vs pH (pentose pot) | norfuraneol vs pH | 3-DX / 1-DX vs pH |
|---|---|---|---|---|
| R0 | **current model** | flat (glycine); falls ~1 dec / pH (cysteine, `acid` factor) | flat (glycine); rises (cysteine, `base` factor) | flat (glycine); rises in Jiangnan's pots (EX-033, scored on R3En, not B60; corrected 2026-11-05) |
| R1 | **M1b as the plan words it** (part H) | as R0 in every pentose-only pot | as R0 | as R0 |
| R2 | the enolisation split moves with pH (E2), carried to the pentose steps | falls | rises (more 1-DX) | falls |
| R3 | base-catalysed norfuraneol decay (E16, `NF_HYDROXIDE_DECAY`) | unchanged | falls above ~pH 6 - 7 at 100 C; at 140 C set by an unmeasured barrier | unchanged |
| R4 | a pH-dependent Amadori decay alone (the pentose + amine entry on, its slope; Martins' k1 - k3 for hexose) | rises with supply | rises with supply | unchanged |
| R5 | a base-favoured 3-DX removal (EX-045, EX-046) | falls (3-DX removed before furfural) | unchanged | falls (standing), with the flux split unchanged |

R2 - R5 can combine; the case asks which the decisive data need.

## 4. The decisive data

All development; Meynier excluded.

| # | measurement (record) | R0 / R1 | R2 | R3 | R4 | R5 | separates |
|---|---|---|---|---|---|---|---|
| D1 | NF loss at fixed supply vs pH: Mikami 2017 fed NF alone at pH 8 (100 C), with the pH 6 / 7 bounds from its xylose + lysine pot; Nakamura 2020 degassed / DTPA arms | pH-free, ~700x slower than measured at pH 8 (F13-a) | as R0 | ~1 dec / pH, O2- and metal-dependent | as R0 | as R0 | R3 from all. Decides R3's existence at 100 C, pH 8; NOT its size at 120 - 145 C or pH <= 7 (no fed NF there) |
| D2 | in-situ trapped split vs pH: Beck 1988, Nedvidek 1992 (hexose) | cysteine pots: steeper than measured (Hofmann pots: 0.4 / 40 / 89 % of pentose to DPO at pH 3 / 5 / 7, SULFUR_FLUX_MAP sec. 5); glycine pots: flat | moves (~0.44 dec / pH in Beck) | no change | no change | **no change** (removal acts after the trap) | R2 from R5 for HEXOSE. One school, n = 1. No pentose in-situ trap across pH exists |
| D3 | standing 3-DX / 1-DX vs pH inside one pentose study: Zhai 2021 (label 5.5 - 8), Hou 2026 (in-pot 5.05 -> 3.55) | rises with pH (wrong sign, EX-033: an R3En result, not B60's; corrected 2026-11-05) | falls | no change | no change | falls | R0 / R1 / R3 / R4 from R2 / R5. R2 vs R5 **not separable** on standing levels (EX-036, EX-046, EX-053) |
| D4 | furfural vs norfuraneol inside one pentose study across pH | -- | furfural / NF falls | furfural / NF rises | flat | furfural / NF falls | **No development study prints both across pH.** Desclaux 2006 has both at pH 6 only. Proxy: FFT / MFT inside Hofmann 1998 T2 (pH 3 / 5 / 7), confounded by H2S speciation and by E6 (MFT runs through NF in the engine) |
| D5 | Amadori decay at fixed charge vs pH: Martins 2003 Part II Table 3 (hexose); Deng 2021 (fed ARP, drifting) | hexose trunk in the sulfur lane pH-free | -- | -- | faster at higher pH, per step | -- | R4's size for hexose. For a pentose ARP no fixed-pH decay ladder exists on disk |
| D6 | acetic acid from the 1-deoxyosone vs pH: Davidek 2008 (xylose + Gly), Davidek 2006a (fed DFG, fed 1-DG) | -- | rises (more 1-DX) | no change | rises (more supply) | no change | R2 / R4 from R5; R2 vs R4 needs the sugar consumed (Davidek 2008 prints xylose) |

**What the data cannot decide.** For pentose pots, R2 against R5 (D2 has no pentose arm; D3 is standing levels). R3's
size where the prompting pots sit (140 C, pH 4.5 - 6.5): D1 measures only 100 C and pH >= 6, and Mikami's own pH 5.5
- 8 ladder shows NF RISING from 5.5 to 6.5 at 100 C with lysine. D4 has no clean development measurement. And part H
alone predicts nothing new for any pentose-only pot, so it cannot be judged on pentose data at all: its decisive data
are hexose (D2, D5, D6 on DFG), which is E2's own evidence and carries the red team's WEAKENED verdict.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1 Mikami 2017 | k measured +- 10 %; Fig. 3 n = 3 with error bars; read error +- 3 % of scale | B61 prereg audit item 3; mikami2017 dossier |
| D1 Nakamura 2020 | "triplicate", no +- value printed; bars readable to +- 20 % of their height | nakamura2020 dossier |
| D2 Beck 1988 | n = 1, no response factors | F06#0; beck1988 dossier |
| D2 Nedvidek 1992 | no replicates, no error, no response factors; temperature unprinted | nedvidek1992 dossier |
| D3 Zhai 2021 | SD read from the error bars; levels on a diacetyl surrogate (ratios robust); CF-064 records an error-bar defect in its TTCA figures | zhai2021 dossier; CF-064 |
| D3 Hou 2026 | triplicate; per-point SD [UNVERIFIED] | hou2026 dossier |
| D4 proxy Hofmann 1998 T2 | triplicates "differed by not more than 10 %" | thiol_submodel.yml header |
| D4 Cerny 2007 | triplicate, "the standard deviation did not exceed 23 %" | cerny2007 dossier |
| D4 Desclaux 2006 | none stated; no CI on the constants | desclaux2006 dossier |
| D5 Martins 2003 Part II | 95 % HPD on each constant; data "at least duplicate", error bars undefined | martins2003 dossier; `trunk_conditions.PH_SOURCE` |
| D5 Deng 2021 | figure read error +- 3 mmol/L | deng2021 dossier |
| D6 Davidek 2006a / 2008 | SIDA acetic acid; replicate SD as keyed in `ph_acid_davidek_acids.yml` [UNVERIFIED: not read for this case] | ph_acid_davidek_acids.yml |
| D3 / D4 across laboratories | between-lab spread from H0 (EX-059, pending) | -- |
| item 7, Meynier 1995 | "triplicate analyses, CV < 20 %" (rel. 0.20 / sqrt 3, as frozen in prereg v6); "except some low-concentration compounds" not applied | score prereg v6 sec. 2D |

## Open questions for the build

1. **Part H or part P?** `ph.trunk_terms` is inert in pentose-only pots. If M1b ships part H only, the case's own
   prompting pots cannot move, and its verdict rests on hexose data (E2). Part P needs a named rival (R2, R3, R4 or R5)
   and a pH law per step; none of `r_pent_caramel`, `r_pent_thermal_dpo` has one today.
2. **Is the pentose + amine entry part of M1b?** R4 is that entry (off in B60, on in R3En with slope 0.57; EX-013's
   norfuraneol bound failed). Turning it on changes every glycine + pentose pot's chemistry, not only its pH.
3. **R3's barrier.** Mikami fixes 100 C, pH 8. The 120 - 145 C value is the class prior; the designed fed-NF arm B
   (`docs/guides/DESIGNED_TIMECOURSE_PROTOCOL_DRAFT.md`) is the decider, and it needs a laboratory.
4. **Cysteine pots.** In cysteine pots the `acid` / `base` factors already give a full decade per pH; adding part H or
   R2 on top double-counts unless one replaces the other. Which factor goes?
5. **The Meynier split question.** The plan's rule (keyed held-out, item 7 only) and the frozen splits (F13#6 DEV in
   v1; Reading DEV for E16 and E8 in v2) disagree. A dated decision should say which governs before the comparison.
6. **The 120 C equality.** Martins' three slopes are equal at 120 C; most sulfur pots run at 120 - 145 C. Declare
   whether per-step slopes are temperature-dependent or the 100 C spread is noise.
