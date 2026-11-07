# Case B21: the aqueous glucosone step `r_ama_g`, an analogue transfer that B1 and B18 now integrate

**Status: items 1 - 5 written 2026-11-07, before any refit. DRAFT, not frozen.** Freezing waits for three
things: the unread candidates of sec. 2.5 are read (by amendment, declared before their data are read) or excluded
with a reason; an independent premortem; and a second reader for the sec. 2 numbers (S2 (d)). Items 6 - 9 are
appended later and never edited into this part. Under [STATUS_AND_PLAN sec. 7](../../guides/STATUS_AND_PLAN.md), this
is an S3 / S4 application, so items 6 and 9 become: band sources verified; before / after; held-out coverage; review.

*Plan step: the remedy of CF-135 (B1) and CF-136 (B18). The frozen-fit audit (PR #107,
[frozen_fit_reproducibility.md](../../../results/validation/frozen_fit_reproducibility.md)) classes B1 and B18
DRIFTED-MATERIAL. Fourteen shipped coordinates (B1 8, B18 6) rest on objectives that no longer exist. Most of the
drift is one step added after both froze: B21's `r_ama_g` (5c5f2857, 2026-09-09). This case settles what that step
should be **before** B1 and B18 are re-derived. A refit now would bend 14 constants to absorb a step whose constant
is an analogue transfer.*

**Mapping to sec. 7 items 1 - 5.**

| sec. 7 item | where it is in this case |
|---|---|
| claim | sec. 1 |
| evidence by kind, with the citation sweep | sec. 2 |
| rivals | sec. 3 |
| decisive data and tolerances | secs. 3.3 and 4.1 |
| decisive measurement | sec. 5 |

Sec. 4, what ships, is the task's item 4.

## 0. Seen by the author, and what is spent

**Read while writing:**

- STATUS_AND_PLAN sec. 7, items 1 - 9 and S1 - S4.
- `frozen_fit_reproducibility.md`, in full.
- CF-135 and CF-136.
- `incumbent_evidence_tiers.md`: its header, sec. 1, and the `k_ama_g` / `k_g_go` rows.
- The B21 artifacts: prereg, fit report, ship rule, and the generator's rows.
- `parameters_dicarbonyl.py` (B21 block); `network.AQUEOUS_GLYOXAL_REACTIONS`.
- The F08 atlas family (`docs/reference/pathway_atlas/F08_oxidative_osone.yml`, every cell).
- The dossiers named in sec. 2, the relevant sections only.
- OpenAlex titles of the sweep in sec. 2.5.
- Crossref metadata of every DOI cited.

**Held-out results already seen:**

- B21's diagnostics on the B1 browning hold-out (median fold error 1.43 -> 1.31).
- The atlas text of test cells F08#7 (Hamzalioglu 2026, which is B21's own fit rows), F08#8 (Gursul 2020) and
  F08#9 (Lee 2022).

**Spent.** These data cannot be decisive here. They may falsify a constant-free rival, never credit one:

- **Martins 2005** (B1's 354 rows). They nominated the claim through the audit.
- **Zhou 2024 and Leahy 1989** (B18's 10 rows).
- **Hamzalioglu 2026** (B21's six fit rows).
- **Quan 2020, Xia 2022, and Leahy's total pyrazine** (B21's diagnostics).

**Computation.** Run through `./scripts/docker_maillard.sh run`, niced, one process, one thread:

- Sixteen single objective evaluations, through `build_frozen_fit_reproducibility.probe` (the `--probe` path), with
  the aqueous-glyoxal block set by `config.update`. No fit, Jacobian or step.
- Hand arithmetic marked **[Z]**.

## 1. What `r_ama_g` claims

**The step.** `r_ama_g: AMA -> G + Gly`, with constant `k_ama_g`. It is in `network.AQUEOUS_GLYOXAL_REACTIONS`, the
only reaction in that block.

- **Rate law:** first order in the Amadori compound (fructosyl-glycine, DFG). It returns glycine. The step is
  labelled "oxidative cleavage", but it has no O2 term. B29's atmosphere factor is exactly 1 for air or an
  undeclared atmosphere, and its air / argon factors were "fitted on bounds, unidentified, not forwarded"
  (`data/lit/measurement_ledger/trunk.yml`, B29 note).
- **No pH term.** The `rate_transfer = "licensed_at_measurement_ph_only"` flag (pH 6.7) is metadata only;
  `trunk_conditions` applies no pH law to `k_ama_g`. So the same constant runs at Leahy's pH 5 and pH 9.
- **Constant:** log10 k(100 C) = -2.034, i.e. 9.25e-3 /min. Hence **[Z]** 2.3e-3 (80 C), 4.7e-3 (90 C) and
  3.21e-2 (120 C) /min.
- **Barrier:** 75.9 kJ/mol, declared, not fitted (Hamzalioglu 2026 Table 2: 75.9 +/- 21.1).
- **S1 tiers:** both `k_ama_g.log10_k_ref` and `k_ama_g.ea_kj_mol` are `analogue`, with **no declared band**. Both
  are among the 26 "analogue coordinates with no declared band".

**The source.** Hamzalıoğlu, Aktağ & Gökmen 2026, J. Agric. Food Chem. 74, 15776, doi:10.1021/acs.jafc.5c14296;
dossier `hamzalioglu2026_extraction.md` secs. 0, 4, 5.

- **System:** pasteurised whole milk; lactose with casein-bound lysine, 195 mmol/kg dm; pH ~6.7. 1 mL in pre-heated
  glass tubes; headspace not stated; atmosphere not in the dossier. 110 / 120 / 130 / 140 C, 0.5 - 5 min.
- **Lactulosyl-lysine is not measured directly.** It is furosine x 3.1.
- **The constant is the authors' fit.** "LacLys -> glucosone" (step 4) is one of 24 constants of the authors'
  multiresponse model: 1.9e-2 / 2.1e-2 / 2.5e-2 (+/- 7.1e-2) / 1.5e-1 /min.
- **The dossier flags (sec. 5, flag 2) that step 4 is resolved only jointly with its sink.** The glucosone sink
  k19 is 6.5 - 45.9 /min, indeterminate. The 130 C column is "a different solution branch". The dossier's own
  verdict is "Does not transfer: any constant as a number".
- **Its evidence kind is "fitted time courses"**, the lowest rank in item 2.

**The transfer.** B21's prereg sec. 2 (ii) gives the justification for milk lactulosyl-lysine -> fructosyl-glycine
in water: "the oxidative cleavage is on the sugar moiety and the amine leaves, so the products are the same; the rate
is not measured for the pair". It declared a +/- 0.5 dex band on every glyoxal, glucosone and pyrazine answer. That
band has no measured law behind it:

- S1 defines an analogue as "a same-class measurement transferred through a measured law, with its band".
- S2 (b) allows a transfer to add a term "only through a measured law".

So the S1 record carries the coordinate with no band.

The sugar (a 1,4-disaccharide), the amine (protein-bound lysine), the matrix (milk: its metals, ascorbate and
riboflavin are uncontrolled), the O2 supply and the pH all differ from the engine's pots. None of these differences
is bridged by a measurement.

**What B21 fitted.** Two free coordinates, log10 k(100 C) of `k_ama_g` and `k_g_go`, with both barriers declared.
The objective was six of the authors' printed constants, not data:

- the four step-4 values;
- the 120 and 130 C step-10 values;
- sigma 0.15 dex where no interval is printed, 0.5 where the half-width exceeds the value.

No integration was involved. Cost 3.54 on 6 rows; reproduced exactly today (audit: B21 REPRODUCES).

The shipped value is not Hamzalioglu's: **[Z]** k4(120 C) carried back with its own barrier gives 6.0e-3 /min at
100 C; k4(110 C) gives 1.0e-2.

**What B21 did not check.** B21's prereg recorded after its verdict that B1's own rows (Martins 2005) "do not want
that drain". The half-sum of squares rose 2017 -> 2495, and DFG's median error went 0.035 -> 0.093 dex. It
announced a joint fit, "B21b", which never ran. The audit now measures the same thing: `r_ama_g` alone moves B1's S
from 4019.0 to 5025.9 and B18's from 27.2 to 295.1.

**The claim under test.** Is the step's existence, size, temperature law and pH independence, carried from milk,
supported for fructosyl-glycine in water at the engine's rate law and concentrations? If not, what should the
coordinate be, and on what should B1 and B18 be re-derived?

## 2. The literature on Amadori -> glucosone / 1-deoxyglucosone in aqueous model systems

**Search.**

- Corpus first: `data/keys/papers.yml` and every dossier mentioning glucosone, the Amadori compound or 1-DG
  (about 50 dossiers, ranked by mention count).
- Then a web search: Crossref / PubMed / publisher abstracts.
- Then an OpenAlex citation sweep (sec. 2.5).

**Bottom line: no direct measurement of the step exists, on disk or found on the web.**

- No first-order constant for fructosyl-glycine, or any fructosyl-amino acid, -> glucosone in water.
- No barrier for it.
- No fed-Amadori experiment that quantifies glucosone against an authentic standard with the atmosphere stated.
- The only Amadori -> glucosone constant anywhere is Hamzalioglu 2026, in milk.
- **A same-class distribution therefore has n = 1 (one laboratory, one matrix) and cannot be assembled.**

What does exist is a set of **measured upper bounds**: the total loss of a fed Amadori compound, which any single exit
cannot exceed. There are also **O2 contrasts** that bear on the step's existence.

### 2.1 Fed Amadori compound, total loss measured (mass balance / fed intermediate: upper bounds on any one exit)

| source (key; DOI) | system | measured total loss of the Amadori compound (first order, /min) | `r_ama_g` alone **[Z]** | glucosone |
|---|---|---|---|---|
| Martins & van Boekel 2003, Part II (`martins_boekel_2003_dfg_amadori_degradation`; 10.1016/s0008-6215(03)00174-5), with Part I (`martins2003b`; 10.1016/s0008-6215(03)00173-3) | DFG 10 mM, 0.1 M phosphate, initial pH 5.5 / 6.8 (drifting), 100 and 120 C, screw-capped tubes, atmosphere not stated, no chelator | Table 3 k1 + k2 + k3 **[Z]**: 4.7e-3 (100 C, pH 5.5), 3.68e-2 (100 C, 6.8), 2.85e-2 (120 C, 5.5), 0.238 (120 C, 6.8); 95 % HPDs printed (Table 3) | 9.25e-3 / 9.25e-3 / 3.21e-2 / 3.21e-2: **x 2.0** and **x 1.13** the total at pH 5.5; 25 % and 13 % of it at pH 6.8 | no glucosone node; not measured |
| Davidek, Clety, Aubin & Blank 2002 (`davidek2002`; 10.1021/jf025561j) | DFG 0.1 M, 90 C, pH-stat 5 / 6 / 7 / 8 (hot electrode), water or 0.1 M phosphate, 7 h; vessel / atmosphere not described | printed: water pH 5 "about 70 % remained" after 7 h (8.5e-4 **[Z]**). Digitised (dossier sec. 4.3): water 7.8 - 9.2e-4 / 2.4 - 3.0e-3 / 8.7e-3 / 1.7 - 1.8e-2; phosphate 2.2 - 2.9e-3 / 8.7e-3 - 1.0e-2 / 1.7 - 1.9e-2 / 2.2e-2 | 4.7e-3 at 90 C: **x 5 - 6** the total (water pH 5); x 1.6 - 2.0 (water pH 6; phosphate pH 5); about half at pH 6 - 7 phosphate | none (no dicarbonyl data) |
| Deng et al. 2021 (`deng2021`; 10.1021/acs.jafc.1c02023) | methionine-glucose ARP 0.2 M in water, initial pH 5.5 / 7.5, 100 / 120 / 130 C, sealed vials, O2 not stated | 4.8e-3 / 1.6e-2 / 3.2e-2 (pH 5.5); 1.2e-2 / 3.7e-2 / 6.6e-2 (pH 7.5); Ea about 72 - 78 | x 1.9 the total at 100 C, pH 5.5 | **levels only**: 10.6 mM at 60 min, pH 7.5 / 100 C (about 5.3 mol %); <= 3 mM elsewhere; semi-quantified on the 3-DG curve |
| Luo et al. 2023 (`luo2023`; 10.1021/acs.jafc.3c06250) | DFG about 4.2 mM in water, pH 3 - 9, 100 / 130 / 160 C, headspace not stated | not first order: k(0 - 30 min) 6.7 - 8.5e-3, k(0 - 180) 1.5 - 1.8e-3 at 100 C, pH 3 - 7 | above the 0 - 180 min total at every pH 3 - 7 | Table S1, **not on disk**; text: glucosone and 1-DG dominate at neutral / alkaline pH |
| Cui et al. 2022 (`cui2022`; 10.1021/acs.jafc.2c03427) | glycylglycine ARP 5 - 20 mM, pH 7 (25 C) unbuffered, 100 / 120 C | 0.026 - 0.028 /min (100 C, 30 min) | about a third | deoxyosones only (semi-quantitative) |
| Zyzak et al. 1995 (`zyzak1995`; 10.1006/abbi.1995.1073) | N-alpha-formyl-fructoselysine, 37 C, pH 7.4 phosphate | t1/2 4.2 - 10.7 d in air, 16.8 - 57.1 d under N2 + DTPA + phytate (x ~4; O2 and chelators confounded) | — (37 C) | not detected (its own t1/2 <= 12 h) |
| Thornalley et al. 1999 (`thornalley1999`; 10.1042/bj3440109) | fructosyl-lysine 0.5 mM, phosphate pH 7.4, 37 C | 37.5 uM/day initial; DTPA cuts glyoxal from it by about 90 % | — | not measured |

**What this establishes.**

- **The step, as shipped (pH-flat), exceeds the total measured loss of fed DFG at pH <= 6, in two laboratories**:
  - Martins, x 2.0 at 100 C, pH 5.5;
  - Davidek, x 5 - 6 in water at pH 5 and x 1.6 - 2.0 at pH 6.

  Deng's methionine ARP agrees in kind (x 1.9 at pH 5.5).
- At pH 6.8 - 7.5 the step fits inside the measured totals, as a 13 - 50 % share.
- No source measures what share of the total goes to glucosone.

### 2.2 One factor varied inside one study: oxygen and metals (bear on existence, not on size)

- **Hofmann & Schieberle 2000** (`hofmann2000b`; 10.1021/jf000076e). Phe-ARP 0.5 M, 0.5 M phosphate pH 7, 85 C.
  - Under air + Cu2+, glucosone is the main dicarbonyl.
  - Under argon, glucosone is "very low" and 1-DG dominates.
  - Figure only; no numbers on disk.
- **Gobert & Glomb 2009** (`gobert2009`; 10.1021/jf9019085), with Glomb, Gobert & Voigt 2010 (`glomb2010`;
  10.1021/bk-2010-1042.ch004). Glucose + lysine 42 mM, pH 7.4, 50 C, 7 d.
  - Glucosone 44.5 vs 5.5 mmol/mol aerated vs deaerated (Ar + 1 mM DTPA), with in-situ o-phenylenediamine (OPD).
  - 0.30 vs 0.16 mmol/mol when trapped after sampling.
  - The authors warn that continuous trapping makes "false glucosone" from the Amadori product.
  - Fed glucosone: t1/2 8 h at 50 C.
- **Kawakishi, Tsunehiro & Uchida 1991** (10.1016/0008-6215(91)84156-9; not on disk). Glucosone is the main sugar
  product of Cu(II)-catalysed autoxidation of an Amadori compound. Title-level reading only.
- **Baker, Zyzak, Thorpe & Baynes 1994** (10.1093/clinchem/40.10.1950; not on disk).
  - fFL with nitroblue tetrazolium, pH 10.35, 37 C.
  - Glucosone is the oxidation product, about 50 % at 10 min.
  - The oxidant is the dye, not O2: an assay chemistry, not a rate in food.
- **Mossine et al. 1999** (10.1021/tx980209e; not on disk). Superoxide from fructosyl-glycine / lysine in
  Chelex-treated water at pH 7, Cu2+-dependent above 1 uM, 0.9 - 3.6e-10 M/s at 1 mM.
  - **[Z]**, from the web reader: about 3e-4 - 1.3e-3 /h of one-electron oxidation per Amadori molecule (temperature
    not in the abstract).
  - This is an order-of-magnitude upper scale for metal-free oxidation, not a bound on the step.

**What this establishes.** The Amadori -> glucosone route is oxidative and metal-catalysed. It is near-absent under
argon / chelation (Hofmann) and x 2 - 8 lower deaerated (Glomb). So its size in a pot depends on O2 and trace metals,
and neither Martins 2005 (B1) nor Hamzalioglu (milk) states either. The step is real in kind; its rate in the
engine's pots is undetermined.

### 2.3 Levels only (not decisive)

These are reported, not used: Hofmann 1999 (`hofmann1999`; glucose / alanine reflux, 2-hexosulose maximum about
0.03 mol %); Hofmann 2000 (`hofmann2000`); Martins & van Boekel 2005b (`martins2005b`; 10.1016/j.foodchem.2004.08.013;
glucose-glycine at pH 4.8 - 7.5, with no glucosone node; its DFG exits sum to 2.8e-3 to 4.1e-2 /min **[Z]**); and
Kocadagli 2016 (both), Lee 2022, Gursul 2020 and Goncuoglu 2016 / 2017, which are glass / cake / juice / hazelnut
matrices, or are glucosone from the free sugar.

### 2.4 Against and for, in one line

| | direction | sources |
|---|---|---|
| against the shipped size at pH <= 6 | the single exit cannot exceed the measured total | Martins 2003 and Davidek 2002 (strong); Deng 2021 and Luo 2023 (concordant) |
| for the step's existence under air | glucosone is the oxidative product of Amadori compounds | Hofmann 2000b, Gobert 2009, Kawakishi 1991, Baker 1994 |
| against any fixed size | the rate depends on O2 and metals, which are unstated where the constant was measured and where it is applied | Hofmann 2000b, Gobert 2009, Zyzak 1995, Thornalley 1999 |
| for a non-zero size at neutral pH | a glucosone level in a fed hexose ARP | Deng 2021 (methionine ARP, semi-quantified, O2 unstated) |

### 2.5 Citation sweep (W1.4), keyed to the criterion "a measured Amadori -> glucosone / 1-DG rate, or total Amadori loss, in water"

**Method.** OpenAlex forward (citing works) and backward (references) for the decisive sources, titles only. Script
in the scratchpad (`openalex_sweep.py`), not committed. Blinding holds.

| source | citing works | references |
|---|---:|---:|
| Martins 2003 Part I | 95 | 32 |
| Martins 2003 Part II | 57 | 19 |
| Davidek 2002 | 116 | 25 |
| Hamzalioglu 2026 | 1 | 72 |
| Deng 2021 | 22 | 41 |
| Gobert 2009 | 173 | 0 |

Gobert 2009 has no reference list in OpenAlex, so its backward sweep is not done. Keyword-matching titles were
screened by hand.

**Unread candidates that could hold a measurement of the step.** None is on disk. They enter only by amendment,
before their data are read:

- **Greifenhagen, Frolov, Hoffmann et al. 2015**, "Oxidative degradation of N-epsilon-fructosylamine-substituted
  peptides in heated aqueous systems", Amino Acids, 10.1007/s00726-015-1940-2. Fed fructosyl peptides in heated
  water, "oxidative": the closest title to a direct measurement.
- **Staempfli, Blank, Fumeaux, Fay et al. 1994**, "Study on the decomposition of the Amadori compound
  N-(1-deoxy-D-fructos-1-yl)-glycine in model systems: quantification ...", Biol. Mass Spectrom.,
  10.1002/bms.1200231007. Fed DFG.
- Luo 2023's Supporting Information Table S1 (glucosone in fed DFG, pH 3 - 9, 100 - 160 C). The paper is on disk;
  the SI is not.
- Kanzler, Haase & Kroh 2014, "Antioxidant capacity of 1-deoxy-D-erythro-hexo-2,3-diulose and
  D-arabino-hexo-2-ulose", J. Agric. Food Chem., 10.1021/jf404322r. Possibly formation data.
- Zhang, Poojary, Olsen, Ray et al. 2019, "Formation of alpha-dicarbonyls from dairy related carbohydrates with and
  without N-alpha-acetyl-L-lysine during incubation", 10.1021/acs.jafc.9b01532. A lactose / amine bridge to the milk
  constant.
- Nomi, Anazawa, Shinzawa, Tamura et al. 2024, "Identification of lactose-derived alpha-dicarbonyl compounds in dairy
  products and elucidation of their formation mechanism", 10.1021/acs.jafc.4c08966. Relevant to which osone
  Hamzalioglu's "glucosone" is.
- Kaufmann, Meissner, Pelke, Mügge et al. 2016, "Structure-reactivity relationship of Amadori rearrangement products
  compared to related ketoses", Carbohydr. Res., 10.1016/j.carres.2016.04.016. Already queued by case M3b.
- Four more, of lower priority:
  - Yaylayan, Huyghues-Despointes & Feather 1994, a review with a kinetics section, 10.1080/10408399409527667.
  - Henning, Smuda, Girndt, Ulrich et al. 2011, 10.1074/jbc.M111.282442: aerated vs deaerated Amadori at 37 C, with
    amide yields only.
  - Wilker, Heinrich & Kroh 2015, 10.1021/acs.jafc.5b04770: glucose / alanine with and without metals.
  - Glomb & Tschirnich 2001, 10.1021/jf010148h: trapping artefacts.

**Excluded (title):** the analytical-method, food-survey and clinical papers in the sweep, and the pentose / dipeptide
ARP formation papers (no hexose-ARP glucosone).

## 3. Rivals

| rival | what it says | runnable | carried? (S3) |
|---|---|---|---|
| **R0** | the analogue as shipped: k(100 C) 9.25e-3 /min, Ea 75.9, pH-flat, O2-flat | yes | yes: analogue-tier constants only. Contradicted at pH <= 6 by sec. 2.1's bounds, but not formally falsified (sec. 3.3) |
| **R1** | a measured aqueous value | **does not exist** (sec. 2) | — |
| **R2** | an S4 band. With no class distribution (n = 1), the band comes from **measured bounds**: lower edge 0 (no measured lower bound for an amino-acid Amadori compound in water under a stated atmosphere; Hofmann's argon arm is "very low"); upper edge the measured total loss of fed DFG at the pot's pH and T (Davidek 2002's pH law at 90 C; Martins 2003 at 100 / 120 C, pH 5.5 / 6.8) | yes (a grid over the band) | yes |
| **R3** | the step absent: B1's and B18's physics for this step (glucosone only from the dry-glass `k_glc_g`) | yes | yes: structural. Its glyoxal levels are 2.5 - 4.2 dex below Quan 2020 (B21 T3, spent; sec. 3.3) |
| R4 (named, not runnable) | Quan's glyoxal comes from another aqueous route: Namiki C2 from the Schiff base (F09), retro-aldol of the sugar or 3-DG (F10), or Wolff autoxidation of glucose (Thornalley 1999, chelator-sensitive) | **no**: no such step in the engine | not until runnable (S3: carried rivals must run). The reviewer's "missing rival" |

### 3.1 B1's and B18's S under each rival (single `--probe` evaluations at the frozen points; not a refit)

Each row reports one objective evaluation at the frozen B1 / B18 point, with `k_ama_g` replaced and everything else
as today. Validation: R0 and R3 reproduce the audit's numbers to 1e-3 (5025.874 / 295.068 and 4019.001 / 27.191).

| rival / point | log10 k_ama_g (100 C) | B1 S (354 rows; recorded 3979.877) | B18 S (10 rows; recorded 2.644) |
|---|---:|---:|---:|
| R0, shipped | -2.034 | **5025.874** | **295.068** |
| R0', Hamzalioglu k4(110 C) carried with its own barrier, no B21 fit | -1.999 | 5164.118 | 285.801 |
| R2 upper edge, Martins 2003 total at 100 C, pH 6.8 (B1's pH) | -1.434 | 12806.498 | 166.382 |
| R2 interior, Martins 2003 total at 100 C, pH 5.5 | -2.328 | 4374.972 | 379.750 |
| R2 interior, Deng 2021's level read as a minimum formation rate (pH 7.5) | -3.054 | 4057.964 | 606.815 |
| R2 interior, 1 % of R0 | -4.034 | 4022.429 | 693.254 |
| **R3**, step absent (= R2's lower edge) | — (k = 0) | **4019.001** | **27.191** |
| R3b, B21 reverted (k = 0 and the glass `k_g_go`) | — | 4019.001 | 14.471 |
| every post-freeze trunk change reversed (audit) | — | 3979.877 | 5.618 |

**How to read this.**

1. **The frozen points were fitted under R3.** S at a frozen point therefore favours R3 by construction. These numbers
   measure how strongly each frozen coordinate set is conditioned on the step. **They are not a comparison of
   rivals**; that needs each rival's own re-derivation (sec. 4.3).
2. **B1 is monotone in `k_ama_g`** (4019 -> 12806 across the band). The B1 rows are spent and are level rows
   (fitted time courses), so under S4 they cannot contract the band.
3. **B18 is not monotone.** It is worst at small non-zero k (693 at 1 % of R0) and best at k = 0 (27). Its frozen pH
   slopes and aminoketone constants were tuned to an almost glyoxal-free supply, and any aqueous glyoxal supply breaks
   them. B18 therefore cannot be re-derived with `k_ama_g` held at any one value.
4. **`k_g_go`'s aqueous value is irrelevant to B1** (R3 = R3b to 1e-6) **but moves B18** (27.2 -> 14.5). It was
   co-fitted in B21 and is out of this case's scope (sec. 4.4).
5. **B1's residual drift once the step is off** (4019.0 vs 3979.9) is the B7 furanic channels and the fed 3-deoxy
   triangle (the audit's other probes). It is not B21's.

### 3.2 Forms and nuisance coordinates (declared now, M1.1)

- **O2 / metal state of each pot:** a discrete nuisance coordinate (air as stated; air-limited closed tube; anoxic),
  with equal weights unless sourced. Martins 2003 / 2005 and Hamzalioglu state none.
- **Davidek's pH-stat pH (hot electrode) vs Martins' initial pH (25 C, drifting):** a form index on the pH axis of
  R2's upper edge.
- **The R2 upper edge uses the point totals.** Martins' Table 3 95 % HPDs and the Davidek digitisation band (three
  fitting methods) enter when the band is built.

### 3.3 Decisive data and tolerances

**Decisive-type data for this step**, all unspent:

- The **ratio** of the step's flux to the total loss of a fed Amadori compound, inside one pot (a route share, which
  cancels scale).
- The **O2 contrast** of that share inside one study.
- The **sign of the pH dependence** of the share.

**What is on disk.** No source measures the share (sec. 2). The total-loss data bound it from above, but cannot place
it within [0, total]. **On disk the rivals R0, R2 and R3 predict the same thing on every unspent decisive-type
measurement, wherever R0 lies inside the bound (pH >= 6.5). The data cannot decide, and the case says so.**

**R0 at pH <= 6.** R0's pH-flat constant exceeds the measured total in the same pot by:

- 0.66 - 0.84 dex in Davidek's water at pH 5 (robust to the barrier's +/- 21.1 kJ/mol: x 4.6 - 6.7 **[Z]**);
- 0.2 - 0.3 dex in Martins' phosphate at pH 5.5, Davidek's at pH 5 and Davidek's water at pH 6.

This is a **structural, constant-free** test: one exit cannot exceed the sum. The noise rule still applies, though.
Neither source prints a run-to-run SD. The comparison is across laboratories (milk constant vs water pot), and the
trunk lane has no H0, so "the largest noise any listed source supports" is the thiol lane's tau_ABS (0.964 dex).
**Under it, no pH <= 6 contradiction stands as an S3 falsification.** Under the within-lab omega (0.419 dex), only
Davidek's pH 5 water arm would.

R0 is therefore recorded as **contradicted by a structural bound at pH <= 6, not formally falsified**. The formal
test, R0's posterior-predictive p in Davidek's and Martins' fed-DFG pots after Holm, is a forward run with no fitted
constant, and it belongs to item 6.

**R3 against Quan 2020.** R3's glyoxal is 4.2 dex (100 C) and 2.5 dex (130 C) below Quan 2020's level (B21 T3).
That is spent data, and R3 is constant-free in this step, so it may falsify. What it falsifies is "no aqueous glyoxal
supply", not "no Amadori -> glucosone step": R4 explains the same level. **R3 stands as a rival for this step; as
the whole glyoxal supply it does not.** Glyoxal / pyrazine outputs under R3 alone are refused (S3's last clause),
whichever rival this case carries.

## 4. What S2 / S4 say should ship

### 4.1 Rule by rule

- **S2 does not apply.** Its conditions fail on both sides:
  - S2 replaces a contradicted fitted-identified constant with a direct measurement of the step. `k_ama_g` is
    analogue tier, not fitted-identified.
  - No measurement of the step exists: none has every lumped channel at the engine's rate law and concentration
    regime (S1 "measured").
  - The total-loss data measure a superset of channels. Under S1 that is a **one-sided bound**, not a value.
- **S4 applies.** The coordinate has no declared band (S1 record), and with n = 1 there is no class distribution. The
  band comes from measured bounds:
  - **lower edge 0**: no measured lower bound in water under a stated atmosphere; Deng's methionine-ARP level is a
    semi-quantified level with O2 unstated, recorded but not an edge;
  - **upper edge the measured total loss of fed DFG at the pot's pH and T**.

  A transfer across pH or T adds a term only through a measured law. One exists for the total loss: Davidek's
  four-point pH law at 90 C (water and phosphate), and Martins' two temperatures x two pH values. It must be fitted
  as a law with its own interval, not read pointwise.
- **S3: R0, R2 and R3 are carried. The envelope is their unweighted union.** R0 lies inside R2's band at pH >= 6.5
  and outside it at pH <= 6. Since R0 is not formally falsified there, the union is [0, max(k_R0, k_total(pH, T))].
  B18 is non-monotone, so the envelope over the band is computed on a grid, not from the edges.

### 4.2 What ships

1. **`k_ama_g.log10_k_ref`** goes from "analogue, no band" to an **S4 band**.
   - **Edges:** [0, max(R0's point, the fed-DFG total loss at the pot's pH and T)].
   - **Sources:** Davidek 2002 and Martins 2003, read by two readers (S2 (d) applies by analogy) and verified with
     OpenAlex `verify`.
   - **The pH and T law of the upper edge:** fitted once to those two sources, with its interval.
   - **No point value ships on its own.** Glucosone, glyoxal and pyrazine outputs are carried over the band. Trunk
     absolute outputs stay refused as the sec. 7 preamble already requires for a lane without H0.
2. **`k_ama_g.ea_kj_mol`** ships as an analogue band, 75.9 +/- 21.1 kJ/mol (the source's printed interval). The
   upper edge of the k band carries the total-loss temperature law instead **[Z]**: about 110 - 114 kJ/mol between
   100 and 120 C in Martins 2003.
3. **The B21 caveat text** (`AQUEOUS_GLYOXAL_CAVEAT`) is rewritten to the band. Its "+/- 0.5 dex" has no measured law
   and goes.
4. **Citation defect, to correct with the band:** `parameters_dicarbonyl._HAMZALIOGLU` and the B21 generator's
   `HAM_ANCHOR` name "Hamzalioglu, Kocadagli & Gokmen 2026". Crossref and the dossier both give **Hamzalıoğlu,
   Aktağ & Gökmen**.
5. **Nothing moves the 14 B1 / B18 coordinates in this case.** They stay fitted-unidentified, carried under S4 as
   bands (CF-135 / CF-136 remedies), until sec. 4.3 runs.

Before / after, held-out coverage and the review (the S4 replacements of items 6 and 9) are appended after the band
is built. They are not computed here.

### 4.3 The later step: re-deriving B1 and B18, conditioned jointly on this case

This is a separate case-gated step, after this case freezes and its band ships. It is not part of this case.

- **What it is.** A **profile over log10 k_ama_g across the carried union**: a grid from the upper edge down to
  k = 0, dense where B18 is non-monotone. At each grid value, B1's 8 and B18's 6 coordinates are free. CF-107's logic
  (S4) requires every coordinate whose sensitivity on those rows is not negligible to be free in the same fit, and
  `k_ama_g` is such a coordinate for both fits (sec. 3.1).
- **Holding `k_ama_g` fixed is not allowed.** In particular, a refit of B1 or B18 at R0 alone is not a re-derivation
  under these rules.
- **What it reports:** each grid point's converged S, the profile of each shipped coordinate across the band, and the
  union envelope. A coordinate pinned at every grid value is a candidate for fitted-identified; otherwise it stays a
  band.
- **The engine-conditional upper edge.** The total-loss bound minus the engine's other DFG exits is tighter. It
  depends on B1's own coordinates (`k_ama_*`), so it is computed inside this profile, never before it.
- **Contraction.** B1's Martins 2005 rows are level rows and spent: **any contraction of the `k_ama_g` band they
  produce does not count** (S4). B18's Leahy pH ratios are within-study ratios (decisive-type) but spent here. Unspent
  decisive-type rows (sec. 5) are needed to contract the band.
- **B18 also carries CF-077's cross-pair factor**, which moves the same rows (S 295 -> 225 alone). The profile frees
  both, or states CF-077's setting as a form index.
- **If an amendment changes this case** (a sec. 2.5 candidate yields a measured share, or sec. 5's measurement lands),
  the profile is rerun on the new band before any coordinate moves. This dependence is why the step is separate.

### 4.4 Out of scope, recorded

`k_g_go`'s aqueous value (B21's second coordinate; 0.309 /min at 100 C, Ea 4.2 +/- 15.7, "consistent with zero")
rests on the same milk fit, flagged by the dossier as unresolved against its sink (flag 2). It moves B18 (sec. 3.1,
R3 vs R3b) and should get the same treatment. Glucosone's measured lifetime in water (Gobert 2009: t1/2 8 h at 50 C)
is the starting point. The note is for its own case.

## 5. The decisive measurement (the case cannot be decided on disk)

**Before W6, in order:**

1. **Read the sec. 2.5 candidates by amendment.** Luo 2023's SI Table S1 is the first, since it may already hold
   glucosone in fed DFG over pH 3 - 9; then Greifenhagen 2015 and Staempfli 1994.
2. **Then the designed measurement**, under W6 (owner, last).

**The designed measurement.**

- **Design:** fed DFG 10 mM in 0.1 M phosphate at pH 5.5 and 6.8 (Martins 2003's design), at 100 and 120 C, in
  sealed tubes with a stated headspace.
- **Arms:** one factor varied inside the study, **air vs argon + 1 mM DTPA**; a fed-glucosone arm for its loss.
- **Measured:** glucosone and 1-DG quantified against authentic standards by derivatisation **after** sampling (not
  in-situ OPD; Gobert's artefact warning), with DFG, glycine, 3-DG and glucose + mannose measured. Replicate runs on
  separate days, so that a run-to-run SD exists.

**What it decides.**

| outcome | means |
|---|---|
| glucosone formed / DFG lost, corrected for the measured glucosone loss | the step's route share, which cancels scale: decisive for R0 vs R2's interior vs R3 |
| that share, air vs argon | whether the step is oxidative in these pots (R3 = anoxic limit) |
| that share, pH 5.5 vs 6.8 | the sign of its pH law (R0 is flat) |

The result would be the first S1 "measured" value for the step. Entered in `docs/reference/decisive_measurements.yml`
as `amadori_glucosone_route_share`.
