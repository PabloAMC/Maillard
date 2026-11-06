# Case M0a: the trunk's pH law in every lane (CF-133)

**Status: DRAFT, not frozen.** Items 1 - 5 of [STATUS_AND_PLAN sec. 7](../../guides/STATUS_AND_PLAN.md), written
2026-11-06 before any result of the change exists. Items 6 - 9 are appended later and never edited into this part.
Items 1 - 5 freeze only after the sweep's candidates (sec. 2.9) are read or excluded.

*Prompted by M0 integrity batch 1, check CK-06 (misfit MF-28; record
`results/validation/m0_integrity_batch1_strecker_sulfur.md`, draft PR #86, branch `m0-integrity-batch1-ss`) and its
register entry CF-133 (PR #88, branch `m0-data-corrections-batch1`). The coordinator chose CF-133's option (a): the same
pH law in every lane, consistent with the one JAX network (W3), which will carry one pH module. It ships only through
sec. 7.*

**Seen by the author.** CK-06's trace and its glycine stand-in numbers (Schibilsky 2019: DMHF x3.18 and HMF x2.06 from pH
5 to 8 in the trunk lane, against measured x5.1 and x1.76). The two Schibilsky bundles. The printed tables of the
dossiers in sec. 2, including the decisive data's (Yu 2017 Table 2, Cha 2019 Table 1); they were read as evidence, and no
engine result was compared with them. Engine calls for this case were read-only and structural (sec. 1.3, sec. 6): lane
routing, and which of the law's steps can carry flux from each charge (forward reachability). **No rival was integrated
on any pot, and no prediction was compared with any measurement.**

## 1. The claim

The trunk lane's declared pH terms act on the same steps wherever those steps run: the three Amadori-decay steps
(`k_ama_tdg`, `k_ama_odg`, `k_ama_mgo`, 0.69 dec / pH, band 0.37 - 0.92; Martins & van Boekel 2003 Part II Table 3, B12)
and the 3-DG -> formic exit (`k_tdg_fa`, 0.77 dec / pH, band 0.65 - 0.89; B40 / B41). Reference pH 6.8. Today these
terms are applied only `if lane == TRUNK`.

### 1.1 Every lane and routing path that bypasses the trunk pH terms today (read from the code, origin/main e560338e)

| # | path | how it is selected | what the trunk's pH steps get there | code |
|---|---|---|---|---|
| P1 | **acrylamide lane** | any of Asn, Gln, Lys, Ala charged (`engine._PRECURSOR_LANE`), or an ACR target | nothing: `trunk_conditions.apply` is never called. `acrylamide_conditions.apply` scales `k_asn_glc` (0.235 dec / pH) and `k_acr_dp` (0.149) only, De Vleeschouwer 2006 | `engine.py` l. 1568 (declarations, `if lane == TRUNK`), l. 2668 - 2694 (acrylamide branch), l. 2695 - 2710 (trunk branch, the only call of `trunk_conditions.apply`); `acrylamide.py` l. 52 "No pH dependence" |
| P1b | the acrylamide lane's competitor amines | Gln, Lys, Ala react only through `a_gln_glc`, `a_lys_glc`, `a_ala_glc` (amine + Glc -> melanoidin, one step) | **the amine has no Amadori chain at all**: with no glycine charged, `r_schiff`, `r_ama_*` carry no flux, so copying the law does not reach these amines (sec. 1.3) | `acrylamide.py` l. 241 - 261 |
| P2 | **sulfur lane** | any of PENT, Cys, THI, ARP, NF, GSH, H2S, FUR, HA, MP; or a thiol / furfural / NF / H2S / ACTZ target | nothing on the hexose steps: the switch `ph.trunk_terms` (`sulfur.TRUNK_PH_TERMS`) exists and is **off**. When on, it applies the same keys at the lane's IN-SITU pH (sub-model 1, Amendment 1). The lane's own `REACTION_PH_FACTOR` covers pentose, H2S and thiolate steps, not the hexose trunk | `sulfur.py` l. 1271 - 1290, l. 1322, l. 1381 - 1386; `config.py` l. 210 - 213 |
| P3 | refused pairings | a pentose with Asn / Gln / Lys / Ala (e.g. Blank 1997 xylose + Ala, Mikami 2017 xylose + Lys, Meynier 1995's lysine arm): LANE CONFLICT, no answer | not applicable (no number is emitted) | `engine.resolve_lanes` l. 950 - 1006 |
| P4 | lipid lane | a lipid carrier; composes with one Maillard lane | lipid-only: no trunk steps. Lipid + trunk: `lanes[0]` is TRUNK, so the law IS applied. Lipid + acrylamide or sulfur: P1 / P2 | `engine.py` l. 1247, l. 1607 - 1627 |
| P5 | the uncertainty envelope | `CoreDraw.trunk_*` | the Amadori exponent's band is drawn only for TRUNK; acrylamide draws its own two slopes | `uncertainty.py` l. 997, l. 1052 - 1061 |
| P6 | the JAX mirror | `kinetic_core_jax.trunk.build(lane="acrylamide")` (the default) | pH mode "none": mirrors P1 exactly; `lane="trunk"` mirrors B12 / B40 | `src/kinetic_core_jax/trunk.py` l. 18 - 30 |
| P7 | final pH | `ProcessSpec.ph_final` | ignored outside the sulfur lane: where the law runs (trunk lane) it runs at the declared INITIAL pH, held over the programme | `engine.py` l. 1574 - 1578 |

Not bypasses, but scope limits of the same law inside the trunk lane (left out of this claim): the condensation
`k_schiff` and the amine-free sugar steps (`k_glc_tdg`, `k_fru_odg`, `k_fru_int`, ...) carry no pH term in any lane
(`trunk_conditions` docstring: "no source contrasts them across pH"). The trunk lane charges Met and Pro as glycine
(B22, B24), so those amines DO get the law there (`engine.py` l. 2598 - 2606).

### 1.2 What the claim changes, by lane

| lane | steps | pH it reads | parts |
|---|---|---|---|
| acrylamide (P1) | `k_ama_tdg`, `k_ama_odg`, `k_ama_mgo`, `k_tdg_fa` | declared initial pH (as the trunk lane) | part A |
| sulfur (P2) | the same keys, through `ph.trunk_terms` = True | in-situ pH (Amendment 1 of sub-model 1) | part S |
| JAX, envelope (P5, P6) | mirrors | -- | follow parts A and S |

### 1.3 Where the law can act at all (read-only engine call, structural)

Forward reachability (`engine._reachable_species`) on the shipped panel's answerable rows (scratchpad script, 1 core,
niced; no integration). The Amadori-decay keys are **live in no acrylamide-lane or sulfur-lane panel row** (no row
charges glycine there). `k_tdg_fa` is live in every hexose row of both lanes (3-DG is reachable from glucose or fructose
through the amine-free `r_glc_tdg`). So, as written, part A acts on the acrylamide lane's pots only through the 3-DG
exit. In a glucose + alanine pot the Amadori chain the measured x3.2 DMHF response rides on (in the glycine stand-in)
does not exist: path P1b, which is why rival R2 (sec. 3) is needed.

## 2. All the evidence, by kind

Split notation as in case M1b (DEV = development; v1 / v2 = the frozen contrast splits,
`docs/reference/contrast_benchmark/split.yml`, `split_v2_per_change.yml`). "Answerable" is the shipped engine's
declaration for the pot (read-only).

### 2.1 Isotope labelling

| source | what it shows | bears on | status |
|---|---|---|---|
| Davidek 2006a / 2008 (F07#7, #8) | acetic acid's carbon origin; "negligible" pH effect on the label origin, pH 6 vs 8 (glycine, DFG, xylose) | the route, not the pH law of non-glycine Amadori decay | DEV; fit rows (EX-013) |
| none | no labelling study on disk measures a non-glycine hexose Amadori route share across pH | gap | -- |

### 2.2 Fed intermediates

| source | what it shows | bears on | status |
|---|---|---|---|
| Martins & van Boekel 2003 Part II (`martins2003_extraction.md` sec. 5) | fed DFG (glucose-**glycine** Amadori compound), 100 / 120 C x pH 5.5 / 6.8: k1 - k3 x3.0 - 15.6 per 1.3 pH | the law itself; glycine only | calibration input (the law's source): cannot be decisive |
| Deng 2021 (F06#6) | fed Met-glucose ARP, initial pH 5.5 vs 7.5 (unbuffered): 3-DG peak ~7 vs <= 1.8 mol % | a non-glycine Amadori compound's decay vs pH (methionine; trunk lane charges Met as glycine) | DEV |
| Davidek 2006a (F07#5, #6) | fed DFG and amine-free 1-DG at pH 4 / 6 / 8 | the 1-DG exit's own law; glycine | DEV |
| Ge 1997 (`ge1997_extraction.md`) | glucose-**Phe** Amadori compound decay, pH 2 - 12, figure only | a second amine's Amadori decay vs pH: for R2 / R5 | not keyed; evidence |
| none | no fed Lys-, Ala-, Asn- or Gln-hexose Amadori compound across pH on disk | the step R2 / R5 need | gap |

### 2.3 Trapping

| source | what it shows | bears on | status |
|---|---|---|---|
| Beck 1988, Nedvidek 1992 (F06#0 - #2) | in-situ OPD; the 3-DG share falls with pH (pH 4.5 -> 7; 5 -> 7); Nedvidek's amine is **beta-alanine** (unmapped) | the split moves with pH for a non-glycine amine | DEV (Ledl school); one school |

### 2.4 Measured intermediates and lifetimes

| source | what it shows | bears on | status |
|---|---|---|---|
| Yu 2017 Figs. 1, 3 (`yu2017_extraction.md` sec. 4.2 - 4.3) | pH 7.4 only: MGO vs time for Lys / Gly / Arg / Pro; glucose lost 66 - 85 % at 0.5 h while the amine lost 3 - 30 % | amines act largely catalytically: a constant amine factor on condensation is incomplete (for R3 / R5) | unseen by the engine |
| Hayashi & Namiki 1986 (`hayashi1986_extraction.md`) | glucose + **beta-alanine**, unbuffered, pH 3.5 / 6.4 / 9.3: C3 fragment x15, browning x4 from 6.4 to 9.3 | a fragment supply that rises steeply with pH (R3) | unmapped amine; M4 band source |

### 2.5 Mass balance

None on disk at more than one pH for a non-glycine hexose pot.

### 2.6 One factor varied inside one study (the core of this case)

| source | pot (lane) | pH | factor -> result, as printed | status |
|---|---|---|---|---|
| **Yu 2017 Table 2** (10.1002/jsfa.8158) | Glc + Lys, Gly, Pro, Arg 40 + 40 mM, 0.2 M citrate-phosphate (held), 120 C, 0.5 h. **Lys: acrylamide lane; Gly, Pro: trunk lane; Arg: unmapped (refused)** | 5 / 6 / 7 / 7.4 / 8 | MGO rises with pH in every arm: 0.50 - 0.76 dec / pH (5 -> 7), 0.18 - 0.36 (7 -> 8); Lys / Gly 7.4 / 3.7 / 2.2 / 2.3 / 2.1 | **unseen by any engine run**; no split, no fit row; decisive (sec. 4) |
| Cha 2019 Table 1 (10.1038/s41598-019-41824-8) | Glc + Lys 0.1 M, unbuffered, set pH only, 100 - 160 C, ~2 h (acrylamide lane) | 3 / 5 / 7 / 9 | MGO at 160 C 8.6 / 49 / 80 / 705 uM; U-shaped at 140 C; many cells below the lowest calibrator | measurement ledger entry (`trunk.yml`, `cha2019_table1_glc_lys_t_ph_grid`, not simulated); decisive, secondary |
| Schibilsky 2019 (bundles; thesis not on disk) | Glc + Ala 0.25 M, unbuffered, 130 C, 2 h (acrylamide lane) | 5 / 8 | DMHF x5.1, HMF x1.76, furfural x1.47, acetol x32, diacetyl x2.8 from pH 5 to 8 | **SPENT** (CK-06) |
| De Vleeschouwer 2006 (`devleeschouwer2006_extraction.md`) | Asn + Glc 0.1 M, 0.05 M phosphate / citrate, 120 - 200 C (acrylamide lane) | 4 / 6 / 8 | acrylamide formation 0.235, elimination 0.149 dec / pH; **no HMF** | calibration input of the lane's own pH factor |
| Martins 2005b (F06#13) | Glc + **Gly**, initial pH 4.8 -> 7.5 | ladder | k4 +0.35, k7 +0.60 dec / pH | DEV; EX-001 / 002 hold-out |
| Haleva-Toledo 1999 (`haleva1999_extraction.md`) | Glc or rhamnose +- **Arg**, 0.1 M citrate-phosphate, 70 C, 48 h | 3 / 5 / 7 | HMF Glc + Arg 92 -> 6 mg/L (pH 3 -> 5), Glc alone 5 -> 0.2: acid-catalysed sugar dehydration dominates HMF at low pH; furaneol rises (rhamnose) | **v1 TEST** (F13#12 - #13, F14#11 - #12); Arg unmapped; 70 C below the 80 C floor: evidence only |
| Tai & Ho 1998 (`tai1998_extraction.md`) | GSH or glutathione sulfonate + Glc, water, 160 C, 2 h (GSH: sulfur lane) | 8 / 6 / 3 | GSH: HMF only at pH 3; HDMF 10.9 / 14.9 / 9.5. Sulfonate: HMF 159 / 1245 / 896 | semi-quantitative, no replicates: evidence for part S only |
| Yeo & Shibamoto 1991 (`yeo1991_extraction.md`) | Cys + Glc 1 M, microwave, temperature not measured (sulfur lane) | 2 / 5 / 7 / 9 | HDMF U-shaped (area ratios) | not simulable (no temperature) |
| Ajandouz 2008 (`ajandouz2008_extraction.md`) | glucose alone, pH 8.0 phosphate vs 9.7 borate, 60 - 100 C | 2 | glucose loss x3.5 - 6 at 9.7 (buffer changes too) | the amine-free sugar steps (R3); confounded |
| Larralde 1995 (`larralde1995_extraction.md`) | ribose alone, pD 4.4 - 8.6 | ladder | sugar loss ~1 dec / pD to pD 6, then flat | R3 for a pentose; buffer switches |
| Meynier 1995 (F13#5, #6) | ribose + Gly / Lys / Met / Cys, 140 C (Lys arm: LANE CONFLICT) | 4.5 - 6.5 | NF falls 16.1x (Gly), 6.2x (Lys), 1.7x (Met), 3.5x (Cys): **the pH response depends on the amine** | DEV (Reading); W1.1-scored; for R5's existence, pentose |
| Blank 1997 (`blank1997_extraction.md`) | xylose + Ala or Gly, 0.2 M phosphate, 90 C (LANE CONFLICT for Ala) | 5 / 6 / 7 | HDMF (Ala) 0.3 / 0.9 / 2.5; the Gly control flat 2.6 / 2.6 / 3.1 | evidence for R5 (amine-specific pH response), pentose, unanswerable |
| Hofmann 1998 T2 (F25#8) | ribose / glucose + Cys, 0.5 M phosphate, 145 C | 3 / 5 / 7 | thiols vs pH | sulfur-lane panel rows; seen; calibration input |

### 2.7 Fitted time courses and fits (engine records)

| record | what it shows | bears on |
|---|---|---|
| B12 (`trunk_conditions`), B40 / B41 | the law shipped in the trunk lane; the same table's 3-DG -> fragments slope REJECTED by the Leitzen MGO row (1.28x -> 33x): a slope from a lumped step does not transfer | **mixed**: transfer of Martins' slopes is not free |
| EX-001, EX-002 (sub-model 1 arms A, B) | the same slopes inside the sulfur lane (part S): acetic transfers (hold-out 5 / 5), formic does not (0 - 1 / 4); nothing shipped | **mixed** for part S |
| EX-003 (arm C) | a condensation slope on `k_schiff` / the glucose steps: in-sample fine, frozen hold-out rejected (low-pH drifting runs 10 - 200x low in acid) | **against** a pH law on condensation (R3's condensation form) |
| M1b case (rescoped 2026-11-04) | `ph.trunk_terms` alone is inert in pentose-only pots; M1b(i) does not build it | part S has no owner |
| CK-06 | the glycine stand-in moves DMHF x3.18 and HMF x2.06 (pH 5 -> 8); the acrylamide lane makes 74x less DMHF at the same charge | nominating diagnostic (spent) |

### 2.8 Spent data (cannot be decisive here)

- **Schibilsky 2019**, both pH arms, every compound: the glycine stand-in was compared with its measured ratios at
  diagnosis (CK-06). Item 7 only.
- Martins 2003 Part II and Martins 2005b: the law's source and sub-model 1's hold-out.
- De Vleeschouwer 2006: the acrylamide lane's own pH factor (B15).
- The Hofmann 1998 sulfur rows, Leitzen 2021 (B40's rejection), Meynier 1995 (W1.1 scored), sub-model 1's frozen hold-out
  (EX-003), Davidek acids (EX-013 fit rows): seen under the engine.

### 2.9 The citation sweep (sec. 7 item 2), 2026-11-06

**Keyed to the decisive criterion:** does the work hold HMF, DMHF, furfural, furan, MGO or another trunk product at
>= 2 pH inside one study, in a hexose pot whose amino acid routes outside the trunk (Asn, Gln, Lys, Ala), or a hexose +
Cys / GSH pot (part S)? Titles and abstracts only.

**Tool and seeds.** `scripts/tools/openalex.py sweep --direction both`, anonymous, ~1 request per 1.1 s, through the
wrapper: Yu 2017 (10.1002/jsfa.8158), Haleva-Toledo 1999 (10.1021/jf9813788), Cha 2019 (10.1038/s41598-019-41824-8),
De Vleeschouwer 2006 (10.1021/jf0611264). Martins 2003 Part II was not swept (glycine; its citers are the trunk's
literature, outside the criterion).

**Result.** 271 works; 26 already held; **245 unseen**. A keyword screen (pH term; a furanic / dicarbonyl / acrylamide
product; a model pot) and a by-eye read of the 97 title-only works:

| decision | n |
|---|---|
| no pH contrast named in title or abstract (food products, reviews, assays, biomedical; 95 of them title-only, read by eye) | 208 |
| pH named, no furanic / dicarbonyl / acrylamide product | 14 |
| pH and product, no defined model pot | 3 |
| pH, product and pot named, read individually, outside the criterion (pressure, salts, oligosaccharides + whey at 50 C, Boc-Arg + MGO at 37 C, rhamnose, histidine, citrus storage, a seafood condiment, theses and reviews) | 16 |
| **candidates** (below: two from the abstract screen, two title-only) | 4 |

**Candidates (unseen; role fixed now, before reading; P = chance it holds the decisive measurement).** The two Ajandouz
papers are not sweep results: the dossier screen found them cited (k3 inventory) and not on disk; their DOIs were
resolved on OpenAlex.

| candidate (authorship from the OpenAlex record) | DOI | role, fixed now | measurement that would count | P |
|---|---|---|---|---|
| Wang Y. et al. 2019, Food Res. Int., "pH and lipid unsaturation impact the formation of acrylamide and 5-hydroxymethylfurfural in model system at frying temperature" | 10.1016/j.foodres.2019.05.001 | **decisive for part A (HMF) only if** an Asn + hexose pot carries HMF at >= 2 pH with the lipid held fixed | HMF at >= 2 pH, one design | ~0.5 |
| Wang P. et al. 2022, JSFA, galacturonic acid and acrylamide, Asn / Glc at pH 3.8 - 7.8 | 10.1002/jsfa.12149 | decisive for part A only if HMF or a dicarbonyl is printed for the Asn / Glc arm across pH; else evidence | as above | ~0.2 |
| Ajandouz & Puigserver 1999, JAFC, pH on caramelisation and Maillard kinetics of essential amino acids | 10.1021/jf980928z | evidence for R3 (sugar steps vs pH, Lys / Thr / Met); decisive only if a trunk product is quantified at >= 2 pH in a Lys pot | sugar loss / HMF at >= 2 pH, Lys arm | ~0.2 |
| Ajandouz et al. 2001, J. Food Sci., pH on caramelisation and Maillard kinetics in fructose-lysine systems | 10.1111/j.1365-2621.2001.tb08213.x | as Ajandouz 1999 (fructose + Lys, acrylamide lane) | as above | ~0.25 |
| Chen X. et al. 2011, Food Res. Int., alpha-dicarbonyls in different sugar-amino acid model systems (title only) | 10.1016/j.foodres.2011.06.002 | evidence (amine-class factor at fixed pH, Gly vs Lys); decisive only if pH is varied | MGO / 3-DG at >= 2 pH, Lys arm | ~0.1 |
| Hollnagel & Kroh 1998, Z. Lebensm. Unters. Forsch., alpha-dicarbonyl fragments under caramelisation and Maillard conditions (title only; TU Berlin, the Schibilsky laboratory) | 10.1007/s002170050294 | evidence; a same-laboratory source is not independent of the spent Schibilsky pots for the between-lab question | dicarbonyls at >= 2 pH, non-glycine amine | ~0.15 |

Excluded with reasons among the pH-and-product works: Illmann 2009 (10.1021/jf803776k; rhamnose + Lys HDMF, rhamnose
not a model species), Liu 2018 (10.3390/molecules24010056; fructose + histidine, histidine unmapped;
its abstract has HMF falling as pH rises: noted as R3 evidence), De Vleeschouwer 2010 HPHT (pressure, not pH), the 1997
citrus-juice DMHF buffer study (10.1021/jf960724n; rhamnose / Arg in storage). Chance that at least one candidate is decisive for HMF in an acrylamide-lane pot (the two Wang papers): ~0.6
(1 - 0.5 x 0.8). For any part A measurement (adding the other four): 0.82 if independent, ~0.7 allowing for the
correlation between the two Ajandouz papers. A candidate enters only by amendment, declared before its data are
read.

**Not swept, and not on disk** (from the dossier screen): Limacher 2008, Fan 2008, Nie 2013 (furan vs pH;
glycine-dominated designs as cited), Perez Locas & Yaylayan 2008.

## 3. The rivals

| id | explanation | Lys / Ala / Gln pot, MGO or DMHF vs pH | Lys / Gly ratio vs pH (same study) | ACR rows | sugar-alone pot |
|---|---|---|---|---|---|
| **R0** | incumbent: pH-flat outside the trunk | flat, exactly (no pH term on any trunk step; `a_*_glc` pH-free) | falls with pH by the Gly arm's slope (the Gly arm moves, the Lys arm does not) | unchanged | flat |
| **R1** | option (a): the trunk's law as is (part A; part S adds the sulfur switch) | flat or FALLING: only `k_tdg_fa` carries flux, so 3-DG leaves faster to formic acid at high pH and less reaches MGO, DDG / HMF | as R0, steeper | unchanged exactly (3-DG is not upstream of ACR) | as R1's Lys pot |
| **R2** | (a'): the law plus an amine-class factor: the competitor amine enters the trunk's Amadori chain (`r_schiff` -> AMA -> `r_ama_*`) with a constant per-amine factor on condensation; the law on decay | rises, with the Gly arm's slope | **constant** (= the factor) | move (glucose competition), if Asn also gets the chain (open question 2) | flat |
| R2b | CF-133 option (b), the factor fixed at 1: route Lys / Ala / Gln to the trunk lane as glycine (as Met and Pro are) | rises, the Gly slope | constant = 1 | ACR pots keep the acrylamide lane | flat |
| **R3** | the pH law sits on the sugar side: base-catalysed enolisation / fragmentation of the sugar (amine-free and amine-catalysed), not on the amine's Amadori decay | rises, amine-independent slope | constant | small | **rises** (Ajandouz 2008, Larralde 1995) |
| **R5** | amine-specific pH response: R2 with a per-amine SLOPE (amine pKa / ARP-decay chemistry; suggested by Meynier 1995 and Blank 1997, where the pH response differs by amine) | rises, the amine's own slope | **moves with pH** | as R2 | flat |

R1 and R2 are the coordinator's (a) and (a'). R3 and R5 are the rivals the evidence suggests (sec. 2.4, 2.6). R2 - R5
can combine.

## 4. The decisive data

Two unseen, within-study pH series in pots that route outside the trunk exist on disk. **Neither measures HMF or DMHF:
both are MGO.** For HMF / DMHF the only on-disk pair in such a pot is Schibilsky's, which is spent.

| # | measurement | R0 | R1 | R2 / R2b | R3 | R5 | separates |
|---|---|---|---|---|---|---|---|
| **D1** | **Yu 2017 Table 2, the Lys arm**: log10 MGO(pH) / MGO(pH 7), pH 5 / 6 / 7.4 / 8 (acrylamide lane) | 0 at every pH | <= 0 | > 0, equal to the Gly arm's | > 0, amine-independent | > 0, its own | R0 / R1 from R2 / R3 / R5: **the sign of a within-arm contrast** |
| **D2** | **Yu 2017 Table 2, Lys arm minus Gly arm** (Gly: trunk lane, the law as shipped): log10 [MGO_Lys / MGO_Gly] at each pH; the shared trunk-law error and the scale cancel | minus the Gly slope | steeper than R0 | 0 slope | 0 slope | non-zero slope | R5 from R2 / R3. With D1, all four groups |
| D2' | Yu 2017 Pro arm (trunk lane, charged as glycine, B24) minus Gly arm | the trunk lane's own prediction, identical across rivals | same | same | same | same | a **control** on the law itself: every rival inherits its error |
| D3 | Cha 2019 Table 1, Lys, sign of MGO(pH 7) / MGO(pH 5) and MGO(pH 9) / MGO(pH 7), 140 and 160 C cells above the lowest calibrator | 0 | <= 0 | > 0 | > 0 | > 0 | as D1, a second laboratory. Unbuffered and set-pH only: the lane runs at the initial pH, so a pH-trajectory effect is not separable |
| -- | HMF / DMHF vs pH in an Asn / Ala / Lys / Gln hexose pot, other than Schibilsky | -- | -- | -- | -- | -- | **none on disk**. Candidates: sweep (sec. 2.9), first Wang 2019 |
| -- | part S: HMF / DMHF / acids vs pH in a hexose + Cys pot | -- | -- | -- | -- | -- | none answerable and unseen (Tai 1998 semi-quantitative, no replicates; Yeo 1991 no temperature) |

**What the data cannot decide.** R2 against R3: the Lys / Gly ratio is flat under both (D2), and no unseen sugar-alone
pH series in a held buffer is on disk (Ajandouz 2008 changes buffer with pH). Part S: no unseen decisive data. The
HMF / DMHF claim of CF-133: no unseen pair. D1 and D2 decide the law's **transfer to an acrylamide-lane amine** on MGO,
one step from `k_ama_mgo`; whether that transfer carries to HMF / DMHF is the measurement below.

**Smallest measurement that would decide HMF / DMHF** (proposed `decisive_measurements.yml` entry, in the PR body,
not in the file): one laboratory, glucose 0.25 M + L-alanine 0.25 M, the same with L-asparagine, and a glycine control
(Ala and Asn: acrylamide lane; Gly: trunk lane), 0.2 M phosphate held at pH 5.0 / 6.5 / 8.0 (label at 25 C, cooled pH
after heating), 130 C, 2 h, sealed with stated fill; HMF, DMHF, MGO, 3-DG / 1-DG (OPD), acrylamide (Asn arm), glucose;
n = 3 independent heatings per cell (27 pots). A glucose-alone arm at the same three pH (9 pots) separates R3.

**Reservation.** The M4 case (branch `m4-case-refresh`, sec. 2.9) lists Yu 2017 as part A evidence and as a band source
"by amendment" for MGO supply vs pH. Data used as a band source in one change cannot also decide another. This case
asks for D1 / D2 to be reserved as decisive here (they are the only unseen within-study contrast in an acrylamide-lane
amine); the coordinator decides before either case freezes.

## 5. Tolerances from the data

| comparison | noise as the source states it | used |
|---|---|---|
| D1 / D2 Yu 2017 | MGO method RSD +-9 % (instrument; recovery 96 %); Table 2's printed SDs (0.1 - 17 %) are below the method's and are not real replicate noise (dossier sec. 3, 6.1) | 9 % per cell: 0.037 dex per cell, 0.053 dex per within-arm ratio, 0.075 dex per D2 double ratio. Cross-buffer comparison (Table 2 vs Fig. 1) is not used: the two buffers differ x1.4 - 1.8 for Gly / Pro |
| D3 Cha 2019 | n = 3, +- SD printed per cell; concentration basis ambiguous (x2, cancels in ratios); cells below the 0.5 ug/mL calibrator are censored, not scored | printed SD per cell |
| across laboratories | D1 and D3 are separate within-study contrasts; no between-laboratory spread of an MGO pH slope exists (H0 covered thiols) | each scored in its own laboratory; no pooling |
| item 7, Schibilsky 2019 | response-factor RSD: furaneol 0.95 %, HMF 21.3 %, furfural 2.16 %; no replicate SD printed | RF RSD as the floor |
| the law's own band | Martins 2003 95 % HPD: 0.37 - 0.92 dec / pH; `k_tdg_fa` 0.65 - 0.89 | carried as the declared band (envelope), not re-fitted |

## 6. Shipped outputs that would move (read-only, structural)

From the shipped panel's 74 answerable rows (6 trunk, 12 acrylamide, 48 sulfur, 8 lipid); `k_ama_*` live in none of the
non-trunk rows, `k_tdg_fa` live in every hexose row.

| rival | rows that move | rows exactly unchanged |
|---|---|---|
| R1, part A | acrylamide lane, through `k_tdg_fa` only: Schibilsky HMF and DMHF at pH 5 and 8 (4 rows, spent), Lin 2022 HMF (fructose + Asn, pH 6.0; demoted on PR #88), Lin 2021 HMF ("Chang2021" water bundle, pH 6.0) | the 6 acrylamide ACR rows (ACSRef3 extrusion, Lin 2022, Chang2021 x3, Ye 2024): 3-DG is not upstream of acrylamide and glucose consumption does not change |
| R1, part S (`ph.trunk_terms` on) | the 10 hexose sulfur-lane rows: Hofmann 1998 glucose + Cys pH 3 / 5 / 7 and fructose + Cys pH 5 (FFT, MFT), Bolton 1994 MFT, Resconi 2023 furfural; through the 3-DG partition (3-DG -> MGO -> MP / ACTZ) and the formic acid the dynamic pH state titrates. Thiol absolute levels ship relative-only (H0), so their within-laboratory ratios move | the 38 pentose / ARP / NF / H2S rows (no live law step) |
| R2 / R2b | as R1, plus the Amadori chain in every Ala / Lys / Gln row (Schibilsky in the panel); ACR rows move through glucose competition only if Asn gets the chain | trunk and lipid rows |
| all | users' predictions in those lanes at pH != 6.8; the JAX trunk module (lane "acrylamide"); the envelope's trunk band (P5); B12 / B15 pins (`tests/unit/test_kinetic_core_b12.py`) | -- |

The coordinator's note on CF-133 stands: shipped acrylamide-lane HMF / DMHF move, so the owner decides before any
engine change.

## Open questions for the build

1. **Part A alone is nearly inert.** As written it reaches acrylamide-lane pots only through `k_tdg_fa`, and the
   structure says that moves HMF / MGO opposite in sign to R2, R3 and R5 (sec. 4, D1). If the build ships the law, it
   must say whether R2 (the competitor amine's Amadori chain) is part of the change; R1 alone cannot answer CF-133.
2. **Does asparagine get the trunk chain?** Asn has its own Schiff base / INT1 route (De Vleeschouwer). Giving it the
   trunk's Amadori decay as well would count the condensation twice.
3. **Initial or in-situ pH.** The trunk lane reads the declared initial pH; the sulfur lane (part S) the in-situ pH.
   One pH module (W3) needs one choice; Yu 2017 is held (no matter), Cha 2019 and Schibilsky drift.
4. **Part S overlaps sub-model 1** (EX-001 - EX-003: acetic transfers, formic does not; the condensation slope rejected).
   A case for part S needs its own unseen decisive data; none is on disk.
5. **R5's per-amine slope** has no independent hexose source: if it is fitted on Yu 2017, Yu cannot decide it (a
   leave-one-arm-out score at most).
6. **The Yu 2017 reservation** against M4 (sec. 4).
