# M1 + M2 pre-registration: MFT route shares with the supply question screened first (v2 DRAFT)

**Status: v2 DRAFT, not frozen.** v1 (9bc55883) went to an independent adversarial review, which returned REDESIGN.
The coordinator decided each finding (sec. 20.1). v2 adds two pre-design checks (sec. 0A, 0B) and applies those
decisions. Freezing needs: a re-review of v2; the synthetic power check (sec. 12.1); the owner's decision on the supply
axis (sec. 16.1 item 1); and a freeze commit. No build, fit or score of any rival here exists. The pre-design screen
(sec. 0B) ran plug-in engine evaluations with no fits, in the session scratchpad.

*Plan steps M1 and M2 ([STATUS_AND_PLAN](../../docs/guides/STATUS_AND_PLAN.md) sec. 5, W2, rows M1 and M2; sec. 7 items
1 - 9; sec. 6, "M1 and M2 in one refit"). Predecessors: M1 (EX-061), its frozen prereg
[kinetic_core_m1_supply_loss_combinations_prereg.md](kinetic_core_m1_supply_loss_combinations_prereg.md) (v3.3), its
[decisive comparison](kinetic_core_m1_decisive_comparison.md), [amendment A1](kinetic_core_m1_amendment_A1.md) and
[results review](kinetic_core_m1_results_review.md); CF-063, CF-076 and CF-107. **This file supersedes the 2026-11-01
E6 + E11 proposal** ([kinetic_core_mft_route_barriers_prereg.md](kinetic_core_mft_route_barriers_prereg.md), never frozen;
its stage 0 is on `origin/mft-route-barriers-2026-11-01`, unmerged). **It reverses that proposal's premortem amendment
(1):** amendment (1) made Cerny pot D, Givaudan and Hofmann development checks only, with CAAS / BTBU as the held-out
tests. Here pot D and Hofmann are decisive (held out of every objective), Givaudan is a reported guard, and the CAAS /
BTBU cells are side effects. Reason: those cells have about zero power on E6 (sec. 2), and BTBU's independence was spent
in M1. There is no M2 case file; secs. 1 - 6 are M2's case items 1 - 5.*

- Written 2026-11-05 on branch `m2-refit-prereg` from main at 46d538c4. The registers and the plan are not edited; sec.
  17 lists the entries for the register session.

## 0. Seen by the author

**The expectations in sec. 15 are not blind.** Everything listed in v1's sec. 0 (kept in the git history at 9bc55883)
was seen: M1's records, the decisive comparison and the review; stage 0; EX-023 / 024 / 030 / 040; the measured cells
of Cerny 2003, de Roos 2005, Hofmann 1998 and Zhai 2023 (the 140 C cells included). Added in v2:

- **The independent review of v1:** every finding (sec. 20.1).
- **Step 1 (sec. 0A):** Zhang 2025b Fig. 2g - h at 400 dpi, and the least-squares refits of its printed lines.
- **Step 2 (sec. 0B): engine predictions on the decisive Zhang 2025b blank pot.** Free cysteine, free xylose, TTCA and
  CAMA at 5 - 50 min, for every screened structure across its bands. These numbers decided only which supply structures
  enter this design. Because they were read, no Zhang 2025b unit can be decisive in any later design (sec. 13 item 2).
- CF-076's battery (EX-054, cluster 5, resolve2_sulfur); `whitfield_mercaptoketone_over_MFT` (16.3); the Cerny 2007
  ledger rows; Mikami 2017.

## 0A. Step 1: is there a Cys / Xyl label swap in Zhang 2025b? Not supported; the cysteine panel's fits are off panel-wide

The review found that the keyed "Cys" points fit the paper's printed xylose line, and the keyed "PENT" points the
printed cysteine slope. Checked against `data/articles/zhang2025b.pdf` p. 7 (rendered at 400 dpi), the dossier
(sec. 3.3) and CF-063:

- **The keying is faithful to the figure.** Panel g is titled "xylose" and panel h "cysteine". The legends match
  (blank black, x1 red, x2 orange, x4 cyan, x6 green, x8 purple, x16 maroon). The black blank points read 8 of 8 in
  each panel as keyed: g -0.78 / -1.25 / -1.57 / -1.85 / -2.2 / -2.45 / -3.15 / -3.5 = PENT; h -0.68 / -1.1 / -1.37 /
  -1.7 / -2.0 / -2.3 / -3.25 / -3.91 = Cys.
- **The printed xylose line reproduces from panel g's own points almost exactly:** ordinary least squares through the
  plotted t = 0 point and 5 - 40 min gives k 0.0729 and b -0.320, against the printed 0.073 and -0.324 (CF-063's
  2026-11-01 finding, re-derived). The review's match (panel h, 5 - 50 min: 0.0715, -0.297) is close but not exact,
  and it needs a chosen window.
- **The printed cysteine line (0.059, -0.760) matches no window of panel h's points** (best: 3 % in k, 0.28 in b).
  The 20 - 50 min window of panel g gives 0.0566 and -0.745, so under a swap the cysteine line would need a window that
  drops the first four points.
- **Panel h's printed fits misfit its lysine arms too.** For example, the x1 line (-0.085x - 0.736, R^2 0.900) misses
  the red points by about +0.4 at 5 - 10 min and -0.6 at 50 min (read at 400 dpi; UNVERIFIED beyond eyeball). The
  misfit is panel-wide, not confined to the blank. That is a property of the authors' cysteine fits, not of the data
  labels.
- **One piece points the other way.** The x1 - x16 lines in panel h share intercepts -0.63 to -0.90. Panel h's blank
  points imply about -0.2 to -0.3, which is panel g's family (-0.30 to -0.34).
- **Verdict: a label swap is not supported. P(swap of the blank series) ~0.15.** The keyed rows stand. CF-063's
  "cysteine fit does not belong to the plotted points" extends to the whole cysteine panel (sec. 17).

## 0B. Step 2: mechanism-signature screen of the supply structures (plug-in, no fits)

**Method.** The engine's own pots (`fit_m1_supply_combinations.run_table_pots`, inline pH, scoring rtol), every
constant at M1's carried C1 optimum except those searched. Each structure was implemented by wrapping
`sulfur.sulfur_rate_constants_at` in the scratch process; nothing in `src/` was edited.

- **Searched coordinates, inside their bands:** `k_ttca_cys`, `ea_form`, `k_pent_caramel`, `k_cys_thermal`,
  `k_ttca_deg`, plus each structure's own.
- **Search:** a bounded least-squares search, 150 evaluations on the three signatures jointly plus 60 on each alone;
  260 - 323 evaluations per structure.
- **Feasibility:** a structure is feasible on a signature if any evaluated point meets its declared screen tolerance.
- **Cost:** 1 778 evaluations in about 0.8 h wall (two processes x 3 workers, nice'd; <= 6 cores). Scripts and records are in the
  session scratchpad (`screen/screen.py`, `screen2.py`, `out_*.json`), not in the repo.

**The signatures and their screen tolerances:**
- **(a) Cao 2017 blank, 90 C:** the mean of K_app = [TTCA] / ([Cys][Xyl]) inside 45 - 75 L/mol, and flat (range <= 0.10
  dex: the measured +-4 % plus two read errors) at 0.5 - 5 h.
- **(b) Zhai 2021 Fig. 4b, TTCA alone at 100 C:** the ln-fraction lost by 80 min at pH 8 over that at pH 5.5 at least
  half the measured 9.6 (0.147 vs 0.819 left).
- **(c) Zhang 2025b blank, 120 C:** TTCA at 30 and 50 min within 20 % of 17.7 / 18.2 mmol/L; CAMA at 50 min within x2
  of 1.61; free Cys and Xyl at 50 min within x2 of 0.80 / 1.22. Zhang's TTCA plateau and free reactants near 1 mM
  together are what imply the ~19 mmol/L 1 : 1 missing pool.

| structure | coordinates and bands (none from Zhang) | (a) | (b) | (c) |
|---|---|---|---|---|
| **S1** (C1) | K 0.060 pH- and T-flat | feasible | **no**: pH 8 / 5.5 loss ratio <= 1.01 at every point (measured 9.6) | **no**: no point puts TTCA within 20 % (best TTCA 11.6 / 8.8); with free Cys low, TTCA is lost too |
| **S2, power law** K ∝ 10^(m pH), slope on the ring opening | m 0 - 1 | feasible in the engine only because the engine's in-situ pH in Cao's pot barely drifts (3.86 -> 3.89, against the measured 4.67 -> 4.15). **On the measured drift (arithmetic): the slope of log10 K_app on pH is +0.013, 90 % CI -0.155 to +0.180** (residual and read errors in quadrature, t4). So m <= 0.18 | **no**: K rising with pH stabilises TTCA; ratio <= 1.0, 0.45 at m = 1 (the wrong sign) | closest of all (chi^2 / point 20) but only at **m = 1**, which (a) excludes; and when TTCA and CAMA pass, free Cys stays >= 3.2 mmol/L |
| **S2, Kallen-shaped** (CF-076's form: K ∝ (1 + 10^(pKa_T - pH)) 10^pH, flat below pKa_T), with S3's drain | pKa_T 5.2 - 6.1 (Radomski 1989, as CF-076); drain as S3 | feasible; on the measured drift K moves <= 0.07 dex for pKa_T >= 5.2 (arithmetic), inside the bound | **no**: ratio <= 1.36 | **no**: when TTCA and CAMA pass, free Cys ~10 mmol/L (chi^2 / point 34) |
| **S3**: K pH-flat in Cao's band + an equal first-order drain on free Cys and free PENT (a 1 : 1 drain at equal free levels) | K 0.045 - 0.075; drain at 90 C <= 3e-3 /min (Cao's pool closes within its 6 - 17 % over-closure over 5 h: arithmetic, about 60 mM / (65 mM x 300 min)), barrier 40 - 135 kJ/mol (`sugar_trunk` class) | feasible | **no**: ratio <= 1.00 | **no**: when TTCA and CAMA pass, free Cys 15 / Xyl 18 mmol/L |
| *diagnostic* S3 + a base-favoured TTCA exit (`r_cama_form`, `r_ttca_deg` x 10^(n (pH - 5.5))) | n 0 - 1 | feasible | **feasible** (ratio up to 12.3) | no |
| *diagnostic* Kallen-shaped + drain + exit | as above | feasible | feasible (13.2) | **no**: each part of (c) is reached alone, never together (TTCA and CAMA pass -> free Cys >= 10; free low -> TTCA lost) |

**What the screen says.**
1. **No screened supply structure is feasible on all three signatures.** Under the coordinator's rule ("only feasible
   structures enter v2") none of S1, S2 or S3 enters as a ship candidate.
2. **(b) needs a TTCA exit that is faster at high pH.** None of S1 - S3 has one. It also cuts against S2's premise:
   any K that rises with pH makes TTCA more stable at pH 7 - 8, the opposite of Zhai's measurement, unless an exit
   outruns it.
3. **(c) is reached by no structure inside bands set from data other than Zhang.** The blank pot's combination
   (TTCA at 18 mmol/L, free reactants at ~1 mmol/L, ~19 mmol/L of each in no measured pool) is outside every banded
   form. This is CF-063 / CF-076's open state; CF-076 gives P ~0.2 to "an unmeasured Cys + Xyl adduct pool or a second
   sink" and ~0.2 to "a much larger, pH-dependent K".
4. **Caveats.** About 40 other constants were held at C1's carried values, fitted under C1. The searches are local and
   budgeted. The tolerances are a screen's, declared here. A real fit could move the held constants. The (c) failure is
   structural, though: with K pH-flat in Cao's band, a TTCA plateau of 18 mmol/L cannot coexist with free reactants
   near 1 mmol/L at any release rate (arithmetic: equilibrium needs K_app >= 18 / (0.8 x 1.22) ≈ 19 L/mmol, 300x Cao's;
   a slow release makes formation equally slow, since k_f = K k_r). P(a full fit inside the same bands reaches (c))
   ~0.1.

**Consequence for this design (the owner decides: sec. 16.1 item 1).** v2 does not judge supply. The route question
(M2) is pre-registered at the incumbent's instant supply (S0), the only supply form on which the shipped constants were
fitted. M1's C1 (S1) is carried as a reported robustness check of the route verdict ("errors stack"), never as a supply
candidate. The supply question returns to the owner with the screen's diagnosis: a candidate needs a base-favoured TTCA
exit and an account of Zhang's missing pool (a measured fourth state, a TTCA quantification bias, or Zhang's pot
declared not representable). That candidate is a new case, not a band change here. The required pair C1F vs C1, S3's
entry, Cao's release floor as a bound and the m-ceiling rule (review MAJOR 2) all wait on it. They are carried, ready to
apply, in sec. 16.2.

## 1. The claim (sec. 7 item 1)

**Route shares (M2, E6 + E11).** At <= 100 C, MFT forms mostly NOT through free norfuraneol (NF). Labelling puts 7 %
through NF. It forms through the 1,4-dideoxyosone (DDP) + H2S route instead. The two routes carry different formation
barriers inside measured classes, so the mix can move with temperature. Changed steps: `r_dpo_nf`, `r_nf_mft`,
`r_ddp_mft`, `r_ddp_mft_hs` (rates and barriers).

**Not changed (review MAJOR 4):** the furfural barriers. B62 is unjudged, so `k_tdp_fur` stays on 64.08 and
`k_fur_decay` on the carbonyl-sink family's 143.54 kJ/mol in every rival. Also unchanged: supply (sec. 0B), the loss
law (C1 vs C2), the disulfide side (C3), the pH balance (M1b), the hexose entry (M4) and the fragment route.

## 2. The evidence, by kind (sec. 7 item 2)

| rank | kind | source (lab) | what it shows | for / against | status |
|---|---|---|---|---|---|
| 1 | labelling | Cerny & Davidek 2003 pot D (Nestle; 95 C, pH 5, 4 h; NF + [13C5]ribose + Cys) | MFT 7 % from fed NF; MP3P 96 % from NF | for E6; the MP3P share also says the NF pool stays fed-dominated | v2 DEV; hold-out, **seen at design**. Decisive: within-lab held-out (sec. 5) |
| 1 | labelling | Cerny & Davidek 2003 pot C | MP3P not detected from ribose | for E6 | blind companion; reported |
| 1 | labelling | Cerny & Briffod 2007 Tables 2 - 4; Cerny 2007b (Firmenich, same first author; 145 C, [13C5]xylose + Cys + thiamine, pH 4 - 7) | xylose vs thiamine share of MFT; C5-intact shares of MFT, MP3P and 3-mercapto-2-pentanone; NF marker (MP3P) absent at pH 4 - 5.5 by TIC (F23#10) | for E6 at 145 C, through the NF marker | **into the objective in v2** (sec. 4.1); the MFT cell at pH 5 duplicates Table 5 and Cerny 2007b's `cerny_ternary_thiamine_share` (CF-007): counted once |
| 1 | labelling | Wang 2012 pot C (BTBU); Liu 2023b, Du 2023 (CAAS) | xylose share of MFT in a thiamine pot; intact-C5 vs fragment shares | about zero power on NF vs DDP (both routes are intact C5); Wang 2012 is confounded by thiamine | v2 TEST; side effects only |
| 2 | fed intermediates | Hofmann & Schieberle 1998 (DFA; 145 C, 20 min, pH 5) | fed NF + Cys gives 5.13x the MFT (ppb) of ribose + Cys | **against at 145 C**: NF is a good precursor there; the T-signature | fed rows FIT (B60); ribose + Cys MFT removed from B60. Decisive (sec. 5) |
| 2 | fed intermediates | Whitfield 1999 (Reading; fed NF, 140 C, pH 4.5) | mercaptoketones : MFT = 16.3 from fed NF | pins NF's branching in-sample | B60 rows |
| 2 | fed intermediates | de Roos 2005 Table I (Givaudan; Laarhoven) | ribose gives 48x the MFT of NF per mol Cys | for E6 | ledger `simulable: false`, `condition_uncertain`; **reported guard** (review MAJOR 7) |
| 2 | fed intermediates | Zhai 2023 Food Chem / JAFC (Jiangnan; fed TTCA and TTCA + Cys, pH 7 water) | the MFT time course at 100 / 120 / 140 C | the route's T-law (E11) | 100 / 120 C in-sample (B60's folds; Kang is a reprint, CF-004); **140 C pots hold-out**; the 140 / 120 MFT ratio seen at design. Decisive as 140 C shapes (sec. 5) |
| 2 | fed intermediates | Mikami 2017 (fed NF, pH 8, 100 C) | NF decays at 0.032 /min | Q2's ceiling | band source |
| 2 | fed intermediates | Zhai 2023 JAFC Fig. 6 (NF + sulfide) | ~34 % MFT | conflicts with Hofmann 0.19 % (CF-065) | not used |
| 4 | measured intermediates | **against the v1 supply design:** Cao 2017's K_app flat while pH drifts 4.67 -> 4.15; Zhai 2021 Fig. 4b's TTCA lost ~10x faster at pH 7 - 8 than at 5.5; Zhang 2025b's mass balance (TTCA plateau 18.2, CAMA 1.6, ~19 mmol/L of each reactant in no measured pool) | bound K's pH slope (m <= 0.18); require a base-favoured TTCA exit; reached by no banded structure | the screen's three signatures (sec. 0B) | Cao, Zhai 2021: route block (in-sample for S1). Zhang: seen in the screen, decisive for nothing in v2 |
| 7 | fitted time courses | Schieberle 2000 Table IV (DFA; 100 C, pH 5) | MFT 360 / 30 rises 1.54x; the FFT shape | R0 fits FFT almost exactly | out of every objective; D1 guard units |
| -- | whole-model fits | B10, R3En, EX-024 / 030 / 035 / 036 / 041, stage 0 | the single barrier cannot move the mix; pot D fails in every joint arm; the NF-step barrier is flat in cost (EX-030: within 4.8 across 64 - 133) | weak (CF-107) | recorded |

**Class evidence for the barriers:**
- **NF formation:** Lau 2003 81.5 +- 3.2 and Pandit 2006 93.0 +- 6.4 kJ/mol. Both are appearance barriers.
- **Sulfur-volatile formation:** Chan & Reineccius 1994 81 - 137; the B55 class (n 15, median 100.9, range
  29.7 - 133.0; UNVERIFIED, from a `source_line` string).
- **Cys -> H2S:** Zheng & Ho 1994, 123 - 135 (unchanged).
- **No step barrier is measured for NF + H2S, DDP + H2S, DPO -> DDP or thiamine -> MFT.** No barrier is computed
  (`assert_no_dft_sulfur`).

**Left out, with the reason:** Zhai 2023 JAFC Fig. 6 (CF-065); van den Ouweland & Peer 1975 (CF-092); Grosch 1992 (E18
TEST); Bornhorst 2017b (pH 8.4 - 9.5); Shu 1988 (area %).

## 3. The rivals (sec. 7 item 3)

### 3.1 Five fitted combinations

| rival | supply | route | role | free coordinates |
|---|---|---|---|---|
| **R0** | S0 instant | Q0 shipped | incumbent, refit | 17 |
| **R0-Q1** | S0 | Q1 (the claim) | **the claim** | 22 |
| **R0-NF** | S0 | Q2 (fed NF lost to its own decay) | the rival explanation | 19 |
| C1 | S1 (M1's reservoir) | Q0 | robustness of the route verdict to supply; M1 comparability. Screened infeasible on (b) and (c), so never a supply candidate | 19 |
| C1-Q1 | S1 | Q1 | the claim at the other supply form | 24 |

Every rival keeps M1's fit-side settings:
- `ph.coupling = "inline"`;
- the dimer encounter ceiling (log10 9.47 at 145 C);
- the `k_mft_decay` weak prior (centre -1.391 at 100 C, sd 1 dex);
- `Ea_decay_thiol_sink` free in Gigl's band (7 - 102);
- M1's switch settings otherwise;
- isothermal fits.

**The base 17 coordinates** are M1 sec. 2.5's: B60's 16 identified, plus `k_mft_decay`. Held at values fitted under
instant supply: `k_dpo_ddp`, `k_dpo_c2c3`, `Ea_lumped_formation`, and in every rival the furfural barriers (sec. 1).
S1 adds `k_ttca_cys` and `ea_form` as in M1. **Cao's floor is not applied to S1** (it stays M1's check), so C1 is M1's
rival.

### 3.2 Q1, the claim (5 coordinates through the B62 hook, in memory; a ship adds a named switch)

Each re-barriered step keeps its fitted log10 k at 145 C.

| coordinate | meaning | band | source | starts 0 / 1 / 2 |
|---|---|---|---|---|
| `k_dpo_nf` | DPO -> NF, log10 k(145 C) (B60 frozen +0.434) | -10 to 0.5 (numerical) | none measured | B60 / B8 perturbation / a seeded uniform draw from the band |
| `d_ddp` | one log10 offset on `k_ddp_mft` and `k_ddp_mft_hs` (B60 -6.541 / -7.693) | -2.31 to +7.04 (each key inside -10 to 0.5; numerical) | none measured | 0 / perturbed / **a seeded uniform draw from the band** (v1's +3 was EX-040's diagnostic pass on pot D, the scored unit: review MAJOR 3) |
| `Ea_dpo_nf` | DPO -> NF barrier | 75.1 - 105.8 kJ/mol | Lau 2003 and Pandit 2006, their +-2 sd envelope (appearance barriers, a declared transfer) | 83.8 / perturbed / seeded draw |
| `Ea_nf_mft` | NF + H2S -> MFT | 40 - 133 | sulfur-volatile class (Chan & Reineccius 1994; B55; the thiol prereg's ceiling 133) | 64.08 / perturbed / seeded draw |
| `Ea_ddp_mft` | DDP + H2S -> MFT (`_hs` tied) | 40 - 133 | as above | 64.08 / perturbed / seeded draw |

- The seeds are printed. Start 0 is the parent rival with `Ea_dpo_nf` moved inside its band.
- The NF - DDP composite gap spans about -82 to +135 kJ/mol across the bands (arithmetic), so the bands do not encode
  E6's answer.
- Held: `k_dpo_ddp` and its barrier, the fragment and thiamine routes, and `k_cys_h2s` (133.0, measured).

### 3.3 Q2, the rival explanation (2 coordinates)

| coordinate | band | source | starts |
|---|---|---|---|
| `k_nf_decay`, log10 k at **100 C** | -10 to log10 0.032 /min | Mikami 2017 Fig. 6B, pH 8 (base-catalysed, so a ceiling for pH 5 - 7) | B60 / perturbed / seeded draw |
| `Ea_nf_decay` (leaves the carbonyl-sink family) | 40 - 135 | `FORMATION_ROUTE_OF` maps `k_nf_decay` to `sugar_trunk`; that class band | 135 / perturbed / seeded draw |

- **What Q2 predicts:** fed NF is spent early, so pot D's MFT and MP3P both come mostly from ribose-made NF. The
  measured MP3P is 96 % from NF while MFT is 7 %. Q1 keeps the NF pool fed-dominated and moves MFT to DDP.
- **Dependency on M1b (sec. 16.1 item 3):** M1b owns the pH law of NF decay (E16). Q2 uses a fixed band with Mikami's
  pH 8 rate as its ceiling. A later M1b verdict that sets the NF-decay law replaces Q2's band. If that law puts NF
  decay at pH 5 below the band's reach of Q2's fitted value, Q2's result here is superseded and the Q1-vs-Q2 contrast is
  re-scored, no refit of Q1, under M1b's law (one amendment). M1b's result never changes Q1's coordinates.

### 3.4 Route identifiability: why the routes may now be identified, and what it means if not

- **Before:** EX-030 found the NF-step barrier flat in cost (within 4.8 across 64 - 133 kJ/mol). R3En freed both route
  offsets and still failed pot D. EX-040 passed pot D only as a diagnostic, with the DDP route raised 3 decades.
  Nothing in B60's objective measures in-situ route shares, so `d_ddp` was a free direction.
- **What v2 adds in-sample (review MAJOR 3 default):** Cerny & Briffod 2007's 145 C labelled shares (Table 4) and the
  NF-marker evidence (Table 2's MP3P vs MFT, pH 4 - 7).
  - Whitfield's in-sample ratio (mercaptoketones : MFT = 16.3 from fed NF) pins NF's branching.
  - So a measured absence of MP3P where xylose-derived MFT is present bounds the NF route's MFT at 145 C, and so sets
    `d_ddp` there.
  - Pot D (95 C) then tests out of sample whether the two barriers carry that mix to 95 C.
- **Why it may still fail:**
  - Cerny 2007's MP3P evidence is TIC areas without an internal standard (semi-quantitative, SD <= 23 %), so the
    absence is a bound only under a declared response-factor ratio (sec. 6.1).
  - 75 - 90 % of that pot's MFT is thiamine-derived, so the xylose-derived MFT the bound acts on is small.
  - The labelled shares mainly inform thiamine vs xylose, already in B60.
- **P(`d_ddp` identified by the Laplace rule at R0-Q1's optimum) ~0.35.**
- **If it is not identified:** `d_ddp` is integrated over its numerical band in the primary density (sec. 7.1). A verdict
  that depends on that integration is **"conditional on a numerical band"**, and it does not ship (sec. 9). The route
  verdict then reads: "E6 is not identifiable on disk at the barrier level; the DDP route's rate is constrained by
  nothing but the scored pot". E6 stays a development finding; the W6 request (Cerny pot D at 120 / 145 C, or Cerny's
  raw 145 C isotopologues) becomes the deciding experiment.
- **Pot D is a within-lab held-out unit.** Cerny 2003 (Nestle) and Cerny 2007 (Firmenich) share a first author and the
  repo's laboratory id `nestle_rc_and_cerny_firmenich`. A shared method bias between the fitted 145 C shares and the
  scored 95 C share is not excluded, and the record says so.

## 4. Objective, fit, numerics

### 4.1 Objective

| block | rows | rivals | measurement model |
|---|---|---|---|
| B60's 52 rows | B60's `row_ids` | all | B60's |
| **Cerny 2007 labelled block (new)** | Table 4: the xylose-derived (M+5) fraction of MFT, MP3P and 3-mercapto-2-pentanone at pH 4 / 5 / 6 / 7, MFT's pH 5 cell counted once with `cerny_ternary_thiamine_share` (CF-007); Table 2: MP3P / MFT TIC area ratio at pH 4 - 7, censored where MP3P is absent | all | logit share; sd from the 23 % area SD (sec. 6.1); the TIC ratio through a declared response-factor band (sec. 6.1) |
| route block | M1's 296 rows minus Zhang and the Zhai 140 C pots | S1 rivals | M1's |
| prior | `k_mft_decay` | all | M1's |

- The Cerny 2007 rows are keyed by the build session (writing under `data/` is the build's, not this file's), checked
  against `cerny2007_extraction.md`, with a schema-gate entry.
- **Spent:** the Cerny 2007 hold-out rows `cerny2007_1x/2x_xylose_share` and `branch_responsiveness` (ledger
  `cerny2007_t5_concentration_pair`) lose their hold-out status. Listed in sec. 10.

Not in any objective: every decisive unit, every side-effect row, Zhang 2025b, Grosch 1992, Zhai 2023 JAFC Fig. 6.

### 4.2 Optimiser, convergence, gradient accuracy

- `fitting.fit`, M1's settings, budget 1 200 evaluations per start. **Budget accounting:** `fitting.fit` passes
  `max_nfev = budget // (n + 1)`, i.e. iterations, each with an (n + 1)-evaluation forward Jacobian.
- **The central-difference polish** (the bake-off's 7 % FD error) is a continuation from the forward endpoint with a
  central Jacobian (2n + 1 evaluations per iteration). It runs with `max_nfev = 300 // (2n + 1)`, i.e. <= 300
  evaluations: 6 iterations for n 22. Status 0: one more call; still 0: the forward endpoint is carried, flagged. A
  central option in `fitting.fit` is added by the build if absent, unit-tested against
  `scipy.optimize.approx_fprime`'s central form.
- **The "forward-difference optimum" form** (a required-agreement form, sec. 7.3) is defined as: the start's forward
  endpoint before the polish, its Laplace recomputed there with central differences, and every decisive density
  recomputed at it exactly as at the polished optimum.
- **Comparison rule, run conditions, A1, series noise, G1:** as v1 sec. 4.2 (M1's, carried).
- **Checkpoints:**
  - every optimiser segment, profile point and IS batch is written to disk on completion;
  - a restarted job resumes from its last checkpoint;
  - long jobs run under `caffeinate -i`, so a laptop sleep cannot enter the objective as a failed pot run.

### 4.3 Verification before any fit

As v1 V1 - V9, adapted:
- **V1:** M1's carried costs at their vectors.
- **V2:** dropped (no S2).
- **V3:** Q1 / Q2 reproduce their parents with bands disabled.
- **V4:** the pot D tracer extended to MP3P, against a two-isotopologue integration.
- **V5:** inline vs fixed point on the decisive pots.
- **V6:** smoothness on `d_ddp`, `Ea_nf_mft`, `Ea_ddp_mft`, `k_nf_decay`.
- **V7:** the IS machinery against M1's R0 Cys quadrature (about -746).
- **V8:** the scorer reproduces stage 0 (pot D 0.981; Hofmann 11.5).
- **V9:** limit numbers per decisive pot.
- **V10 (new):** the Cerny 2007 rows reproduce their dossier values; the censored MP3P likelihood is checked on a toy.

### 4.4 Feasibility probe (review MINOR)

For R0-Q1, R0-NF and C1-Q1: 256 Sobol draws over the new coordinates. Each draw is scored on the rival's claim units,
**its guard units and the Cao / Zhai 2021 blocks** (in-sample signatures), by the plug-in density, then a local search
of <= 200 evaluations.
- **Kill rule:** a rival is killed only if its best value on its claim units is not above R0-shipped (b)'s, with the
  stacking guard (a would-be kill is confirmed by one fit).
- Guard units and the Cao / Zhai blocks are printed, never used to kill.

## 5. The decisive data (sec. 7 item 4)

| block | unit | data | role |
|---|---|---|---|
| **Rt** | Rt-D | Cerny 2003 pot D, MFT fed-NF share 0.07 (logit) | claim, Q1 vs Q0 |
| Rt | Rt-H | Hofmann 1998, 145 C, pH 5: ln[MFT(ribose + Cys) / MFT(NF + Cys)], ppb, panel path; measured ln(1 / 5.13) | claim, Q1 vs Q0; claim, Q1 vs Q2 |
| Rt | Rt-DMP | pot D, MP3P fed-NF share 0.96 | claim, Q1 vs Q2; guard, Q1 vs Q0 |
| **T** (E11) | T-140 | Zhai 2023 JAFC fed TTCA, pH 7, 140 C: the MFT time-course **shape** at 40 / 80 / 120 / 180 min (1.129 / 3.043 / 5.907 / 4.012 ug/L), level free (3 shape degrees of freedom) | claim, temperature law |
| T | T-140c | the same, TTCA + Cys arm (1.576 / 3.126 / 7.189 / 5.037) | claim, temperature law |
| **D1** | D1-MFT, D1-FFT | Schieberle 2000 ratios (M1's) | guards |
| reported | Rt-G | de Roos 2005 ribose / NF = 48 | reported guard (no veto, review MAJOR 7) |
| reported | the reservoir units | Zhang 2025b (M1's S-Cys, S-Xyl; S-TTCA) | reported for M1 comparability only: seen in the screen (sec. 0B) |

- **T scores the 140 C shape, not the 140 / 120 ratio (review MAJOR 5).** The 120 C cells (`zhai2023jafc_ttca_ph7_120C`
  at 120 min, = Food Chem Table 1) are the numerators of B60's `zhai_MFT_fold_120_over_100` and
  `zhai_FFT_fold_120_over_100`: a 140 / 120 ratio would share an in-sample denominator (CF-005's kind of reuse). The
  140 C shape uses only 140 C cells (all hold-out). T-FFT is dropped (MAJOR 4: the FFT side's barriers are held).
- **The T verdict requires the Gigl-ceiling profile to agree (MAJOR 5):** the T units' sides are recomputed at
  `Ea_decay_thiol_sink`'s profile points (sec. 9). If a side changes, T reads "conditional on the Gigl ceiling" and
  E11 is not credited.
- Pot D's concentration is nominal (x0.5 / x2 required-agreement); its thermal reading is isothermal (heating block).
  Hofmann's pots take the six thermal readings (CR-03, CR-04; M1 sec. 6).

### 5.1 Contrasts

| contrast | pairs (a vs b) | claim units (all must favour a) | guard units (sec. 8.1's veto rule) |
|---|---|---|---|
| **route** Q1 vs Q0 | R0-Q1 vs R0 (decides); C1-Q1 vs C1 (robustness) | Rt-D AND Rt-H | Rt-DMP, D1-MFT, D1-FFT, T-140, T-140c |
| **temperature law** (E11) | R0-Q1 vs R0 | T-140 AND T-140c | D1 |
| **route vs rival explanation** Q1 vs Q2 | R0-Q1 vs R0-NF | Rt-DMP AND Rt-H | Rt-D, D1, T |

## 6. Tolerances: the measurement model per unit (sec. 7 item 5), all in ln units

**omega in the primary covariance (review BLOCKER 1).**
- H0's within-laboratory error, omega(all) = 0.390 dex = **0.898 ln** (EX-059, reset set), enters every unit's primary
  covariance: per printed value, independent. H0's ship branch puts omega in the intervals of relative within-lab
  predictions, and these units are such predictions.
- So a ratio of two values carries 2 omega^2, a logit share carries omega^2 once, and a level-free shape of k points
  carries omega^2 on each point.
- **tau** (H0's tau_ABS, 0.964 dex) enters no unit: none is an absolute level compared across laboratories.
- Conversions: 0.041 dex = 0.094 ln; read errors 0.026 / 0.035 dex = 0.060 / 0.080 ln.

| unit | value basis | measurement noise (ln) | primary sd (ln, with omega) | floor |
|---|---|---|---|---|
| Rt-D, Rt-DMP | ln odds of the share, one injection | Cerny 2003 prints none (duplicate); transfer: Cerny 2007's 23 % area SD on each isotopologue, sqrt(2) ln 1.23 = **0.293** | sqrt(0.293^2 + 0.898^2) = 0.945 | none (a share) |
| Rt-H | ln ratio of two pots of one study, SIDA | Hofmann 1998, 10 % per value: 0.095 each | sqrt(2 (0.095^2 + 0.898^2)) = 1.277 | 0.1 ug per pot (the lab's "< 0.1 ug") |
| T-140, T-140c | ln level, 4 points, level integrated out (flat) | per point: the lab's median printed relative SD for these pots' Table 1 cells, 9 % (0.086); declared | per point sqrt(0.086^2 + 0.898^2) = 0.902 | MFT 0.28 ug/L (the printed curve's blank-equivalent) |
| D1-MFT, D1-FFT | ln ratios within one series, bivariate (shared 360 min point) | M1's 0.094 per point | per point sqrt(0.094^2 + 0.898^2) = 0.903; ratios and covariance from it | 1 ug/L ("< 0.1 ug" per 100 mL) |
| Rt-G (reported) | ln ratio of two pots | transfer: 23 % per value | sqrt(2 (0.207^2 + 0.898^2)) | 0.05 mg / mol Cys |
| reported Zhang units | **first-point-normalised** ln(C(t) / C(5 min)) (review MINOR; v1's sensitivity becomes primary), so any multiplicative calibration cancels whatever C0 is | Cao's RSD < 10 % (0.095) with the read error | per point with omega | Zhai 2021's LOD (0.025 / 0.020 mmol/L); 0.223 the alternative |

- **Censored MP3P (the Cerny 2007 objective rows):** where MP3P is "0" (not detected), the likelihood is the
  probability that the predicted MP3P / MFT area ratio lies below a declared detection ratio. No limit is printed; the
  declared value is the smallest MP3P / MFT area ratio the same table reports as detected. That value is taken from
  other cells, not the censored one, and is a required-agreement sensitivity at x0.1.
- The MP3P / MFT response-factor ratio is banded 0.5 - 2 (declared; SPME TIC, no internal standard) and integrated over.
- **Noise sensitivities (required agreement):** each measurement sd at 0.5x and 2x; omega at H0's other set (warm,
  0.419 dex); floors' alternatives; pot D's concentration; Rt-H's thermal readings.

## 7. Densities, ESS, tie band

### 7.1 The predictive (review MAJOR 3)

**Primary:** the integral of the unit's likelihood over the rival's parameter distribution:
- identified coordinates: Laplace (central Hessian);
- bound-active coordinates: their one-sided profile density;
- **unidentified coordinates: integrated over their declared band** (uniform on the band's scale: log k or kJ/mol).
  Not held at the optimum, as v1 had it.

**The plug-in density is printed beside it and is a required-agreement form**, so a sign bought by integrating a wide
band cannot decide. A side that holds under the primary but not the plug-in is a tie. A verdict whose sign depends on
integrating a coordinate over a NUMERICAL band (no measured class) is labelled "conditional on a numerical band" and
does not ship (sec. 9).

**Estimator.**
- Importance sampling with a proposal tuned to the unit: a multivariate Student-t with 5 degrees of freedom, centred
  at the maximum of the unit likelihood x the parameter density, with the inverse central Hessian as its scale;
  unidentified coordinates are drawn from their band.
- Draws are reused across noise sensitivities by re-weighting.
- **Accepted only if both:**
  - the Pareto-smoothed IS diagnostic k-hat < 0.7;
  - the ESS meets max(100, 2 / (0.15 sd_rep(d_u))^2), the second term from v1's rule tying MC error to the tie band's
    own precision.
- S doubles from 256 to a cap of 4 096. If still not accepted: Gauss-Hermite quadrature on the one or two dominant
  coordinates; else "density not resolved" = tie.

### 7.2 The tie band and when a unit favours a rival

- d_u = lppd_u(a) - lppd_u(b).
- sd_rep(d_u) is its spread under replicate redraws of y_u from the unit's primary covariance (R = 4 000).
- b_u = z_0.95 sd_rep(d_u): the 90 % half-width, the plan's interval convention.
- **A unit favours a** iff d_u > b_u under the primary density AND the sign holds under every required-agreement
  form:
  - the plug-in density;
  - the forward-difference optimum (sec. 4.2);
  - each noise sensitivity and floor alternative;
  - each thermal reading;
  - pot D's concentration readings.
- Otherwise a tie.

## 8. Decision rule, outcomes, what ships

### 8.1 Units, blocks, guards

- **Claims:** every claim unit of a contrast must favour a (sec. 5.1); no majority rule.
- **A guard unit vetoes** (it favours b) only if b is adequate on that unit or the rivals differ there by more than
  omega. Adequate means a posterior-predictive p > 0.05 under the primary covariance, two-sided, chi^2 tail for a
  multivariate unit. "Differ by more than omega" means |mean_a - mean_b| > omega on the unit's scale.
- Otherwise the guard reads **"both inadequate"**: no veto, and the record names the unit.
- **A contrast is credited** iff its claim units favour a in every deciding pair, no guard vetoes, and the in-sample
  guard (8.2) holds.
- **Robustness pair** (C1-Q1 vs C1): if it reverses a claim unit's side, the route verdict is "conflict: errors stack"
  and nothing ships. A tie there is recorded as "credited at instant supply; robustness at the reservoir supply not
  shown".

### 8.2 The in-sample guard (review lesson 5, MINOR: Holm)

- Per row group (paper x species) of B60's rows and the Cerny 2007 block: Delta c_g at the carried optima.
- Its replicate sd: sqrt(sum (r_a - r_b)^2), with each row's sigma taken as max(the declared sigma, omega).
- One-sided p_g from z = Delta c_g / sd. Holm's step-down across groups at family-wise 0.10 (the plan's 90 %
  convention).
- A group **lost** under Holm must:
  - be named in sec. 10's table, AND
  - have an evidence kind that does not rank above the best-ranked claim unit won. Labelling ranks first, so a
    labelling win can excuse a fed-intermediate loss; a win on fitted time courses cannot.
- Otherwise: "credited, held: an unexplained in-sample loss" (sec. 7 item 7).

### 8.3 Outcomes and what ships

| route Q1 vs Q0 | Q1 vs Q2 | E11 | meaning | ships |
|---|---|---|---|---|
| credited | credited | credited | a second MFT route with its own barrier carries most MFT at <= 100 C, with a T-law that moves the mix | R0-Q1's structure and constants (the route switch; supply unchanged) |
| credited | credited | not credited | route shares right at 95 - 145 C; T-law not attributed | the same, recorded "T-law not attributed" |
| credited | not separated / Q2 | -- | the route evidence is equally explained by fed-NF loss | nothing; E16 (M1b) decides next |
| not separated / conflict | -- | -- | E6 not decidable on disk at the barrier level | nothing; W6 request |
| Q0 credited | -- | -- | the incumbent explains pot D and Hofmann better | nothing; E6 re-opened as evidence against |

- **Owner sign-off at ship time:** a ship moves every pentose MFT prediction, so it needs the owner's sign-off then
  (merge-authority rule). That is stated here, not decided.
- Constants stay provisional until the joint refit (JR).

## 9. Compensating structure (lesson 6)

- **Bound-active:** profiled at three points inside the bound, others re-optimised (<= 300 + 300 evaluations).
- **Switched off** (a rate at -10, "the step does not happen"): one profile point at the other rival's value. This
  applies to every rival, **R0's `k_arp_dpo` included** (it sat at -10 in M1's R0).
- A side change at any profile point labels the verdict "conditional on ...".
- **Classes:**
  - conditional on a **measured** class bound: reported; ships.
  - conditional on a **numerical** bound or band (the -10 floors, the dimer's ceiling, `d_ddp`'s and `k_dpo_nf`'s bands,
    the declared m ceiling if S2 ever returns): does not ship.
- **The Gigl ceiling** is a measured class bound. Its profile also gates T (sec. 5).
- **The shared-coordinate divergence record:** as v1.

## 10. Side effects (sec. 7 item 7)

As v1 sec. 10: every row before and after; the contrast benchmark; the E6 / E11 TEST cells (CAAS; BTBU tagged lost in
M1); the earlier proposal's H1 / H2; Cerny pot C; the declared hold-outs; M1's disulfide record. Spending added in v2:
Cerny 2007's hold-out rows (sec. 4.1).

**Candidate compensating errors (named before the run):**

| expected loss | why | candidate error |
|---|---|---|
| `zhai_MFT_fold_120_over_100`, Kang MFT rows | the route's T-law moves | one barrier for every route (E11) |
| `fed_nf_h2s_MFT`, `fed_nf_cys_MFT`, Whitfield NF rows | the NF step's barrier and the NF decay move | CF-065 (fed sulfide charged as dissolved H2S) |
| `cerny_ternary_thiamine_share`, `fed_thiamine_MFT` | the non-thiamine MFT changes | the thiamine route fitted against an NF-only network |
| `zhou_pH7_*` | 120 C, off 145 C | E5 (no buffer catalysis); E2 (no pH on the trunk) |
| the dimer rows | MFT supply changes | the declared ambient OX; the irreversible dimer (C3) |
| the 145 C fed + H2S rows | the `k_mft_decay` prior | the frozen fast sink (CF-107); CR-03 |
| Kumazawa FFT retention, every FFT row | none expected (FFT barriers held) | a build error until shown otherwise |

## 11. The adversarial review

- (a) The design review of v1: done (REDESIGN; sec. 20.1).
- (a') **A re-review of v2** by a fresh agent before freezing. Its brief: the screen (sec. 0B) and whether it was fair;
  the identifiability argument (sec. 3.4); the omega-and-guard rule (sec. 8.1); the Cerny 2007 keying; the power
  check's design.
- (b) The results review, as v1.

## 12. Order of operations, blinding, run-once

### 12.1 Freeze prerequisite: a synthetic power check (review BLOCKER 1)

**Question:** under the primary covariance (with omega) and the required-agreement forms, how often does the rule
credit Q1 when Q1 is true, and how often when Q0 is true?

**Design.**
- For the route contrast's claim units, from surrogate predictive distributions built before any fit:
  - **Q0-truth:** the shipped engine's values (stage 0: pot D share 0.981, Hofmann ratio 1 / 11.5), with parameter
    spreads of 0.1 / 0.3 / 0.5 ln.
  - **Q1-truth:** the measured values, with the same spreads.
- 2 000 synthetic datasets per truth, drawn with the primary covariance.
- Each is scored by the rule's actual code, IS included: d_u, b_u, the forms modelled as their declared perturbations
  (noise scales exact; thermal readings and concentration as the stage-0 spread of those forms, about +-0.3 ln on
  Rt-H, recomputed from the engine at B60).
- **Reported:**
  - P(credit Q1 | Q1 true) as a function of the true effect on Rt-D (0.5, 1, 2, 3 dex);
  - P(credit Q1 | Q0 true);
  - P("both inadequate" on each guard).
- The same for Q1 vs Q2 (Rt-DMP) and for T.

**What follows:** reported to the owner before freezing. No threshold: the owner decides whether a low power makes the
run not worth its compute (sec. 16.1 item 2). Arithmetic preview: with omega, Rt-D's shipped miss is 6.5 ln = 6.9 sd,
Rt-H's 0.81 ln = 0.6 sd. So Rt-H alone can rarely favour anything, and "Rt-D AND Rt-H" may be nearly unattainable.
That is the main thing the check measures.

**Budget:** ~3 - 4 h of work (the scorer's code is the build's anyway; the surrogate harness is new), ~0.5 h compute.

### 12.2 Order

0. Re-review; the power check; the owner's supply decision; freeze.
1. Build: the Q1 / Q2 hooks, the MP3P tracer, the Cerny 2007 rows (build session), the scorer (IS, PSIS, quadrature,
   tie bands, forms, the guard adequacy rule, Holm), the central polish, checkpoints; generator
   `scripts/generators/fit_m1m2_supply_routes.py` modelled on M1's. Unit tests. V1 - V10.
2. Probe. 3. Fits: R0, C1 (wave 1); R0-Q1, R0-NF, C1-Q1 (wave 2); polish; records read for `failed_pot_runs`; comparison
   rule; Laplace; profiles; G1. 4. One scorer run after every fit. 5. Side effects. 6. Results review. 7. Decision.

**Run-once and freeze checks:** as v1 sec. 12 (M1's and H0's discipline, numbered amendments).

## 13. Differences from the plan, M1, v1 and the 2026-11-01 proposal

1. **Supply is not judged here** (sec. 0B). The plan's "one refit with M1's combination" becomes: routes refit at
   instant supply, with M1's C1 as a robustness pair. The supply step needs a new case (sec. 16.1 item 1).
2. **Zhang 2025b cannot be decisive in any later supply design** run by this author's session: its engine predictions
   were read in the screen. A later design must take its reservoir block from elsewhere, or have an independent author
   who did not see sec. 0B's numbers.
3. **Cerny 2007 enters the objective;** pot D is within-lab held-out (sec. 3.4).
4. **The furfural barriers are held** (B62 unjudged); T-FFT and v1's `Ea_tdp_fur` / `Ea_fur_decay` are gone.
5. **T scores 140 C shapes**, not a ratio with an in-sample denominator.
6. **Givaudan is a reported guard.**
7. **omega is in every unit's primary covariance;** guards veto only against an adequate rival.
8. **Unidentified coordinates are integrated over their bands** (primary), with the plug-in required to agree.
9. **The first-point-normalised Zhang form** is primary (reported units only now).
10. Superseding the 2026-11-01 proposal and reversing its premortem amendment (1): see the header.
11. v1's corrected numbers:
    - Zhai 2021's 140 C first-order equivalent is 0.0081 /min (dossier), not 0.0087.
    - Cao 2017's first apparent K is 58.6 L/mol (v1 printed 58.5); the series is 58.6 / 64.6 / 60.2 / 59.0 / 58.5 /
      61.1 at 0.5 - 5 h (range 58.5 - 64.6).
    - Cao 2017's free reactants plateau at 15 - 20 %; they do not vanish (unchanged from v1).

## 14. Compute estimate (current engine, <= 6 cores)

**Basis:** the W1.3 profile (one evaluation R0 3.92 s, C1 9.86 s serial; ~3.5 - 4 serial-equivalents machine-wide; the
FD Jacobian ~90 % of a fit) and this screen (3 - 4 s per 5-pot evaluation at 3 workers). The Cerny 2007 block adds one
pot per pH (5 pots, about +0.3 s serial per evaluation). Wall = serial / 3.5 at <= 6 cores.

| step | serial h | wall h |
|---|---|---|
| power check (sec. 12.1) | 0.5 | 0.5 (+ 3 - 4 h work) |
| V1 - V10 | ~4 | ~1.5 (+ 2 h work) |
| probe | ~1 | ~0.3 |
| fits, forward: 3 starts x (R0 18 x 12 x 4.2 s; R0-Q1 23 x 15 x 4.2; R0-NF 20 x 13 x 4.2; C1 20 x 12 x 10.2; C1-Q1 25 x 15 x 10.2) | ~6.5 | ~1.9 |
| central polish (<= 300 evaluations per start) | ~2.3 | ~0.7 |
| continuations (~40 % of starts) | ~1.3 | ~0.4 |
| C1 rivals' series-noise continuation | ~1.7 | ~0.5 |
| Laplace (central) | ~0.4 | ~0.1 |
| **profiles** (review MINOR, recomputed): M1 profiled 2 - 4 coordinates per rival; here 2 - 3 per rival x 5 rivals, 3 - 4 points, <= 600 evaluations per point (typically ~300), at 4.2 / 10.2 s | **~26** | **~7.5** |
| densities: 5 rivals x ~8 units, ESS-driven IS (256 - 1 024 draws, PSIS), x 6 thermal readings on D1 / Rt-H; plug-in, forward-optimum form, band integration | ~10 - 15 | ~3 - 4 |
| side effects | ~4 | ~1.2 |
| **total** | **~58 h serial typical** | **~17 h unattended wall typical; worst ~45 h** (every start at its cap and continuation, 4 profiled coordinates per rival, IS at its cap on half the units) |

Plus **~20 - 26 h of author work** (build, keying with the build session, power check, V-steps, reviews). Long jobs run
under `caffeinate -i` with checkpoints (sec. 4.2), nice'd, memory checked at launch, never paused.

## 15. Expectations (not blind)

| event | P |
|---|---|
| the re-review accepts the screen's conclusion that no supply structure is feasible | 0.7 |
| the power check finds P(credit Q1 \| Q1 true, 1-dex effect on Rt-D) < 0.5 (driven by Rt-H's weakness under omega) | 0.6 |
| `d_ddp` identified at R0-Q1's optimum | 0.35 |
| Rt-D favours Q1 (R0-Q1 vs R0) | 0.55 |
| Rt-H favours Q1 (R0-Q1 vs R0) | 0.15 |
| route contrast credited / not separated or conflict / Q0 credited | 0.12 / 0.83 / 0.05 |
| Q1 vs Q2 credited for Q1 | 0.25 |
| E11 credited | 0.1 |
| at least one verdict "conditional on a numerical band or bound" | 0.6 |
| at least one guard reads "both inadequate" | 0.8 |
| anything ships | 0.1 |

## 16. Owner-level decisions

### 16.1 Decisions

1. **The supply axis.** No screened structure is feasible on Cao's drift, Zhai's pH series and Zhang's mass balance
   together (sec. 0B). The options:
   - **(a), the default:** run M2 at instant supply now, as written, and open a new supply case. The case: a
     base-favoured TTCA exit; an account of Zhang's ~19 mM missing pool (a fourth state, a TTCA quantification bias, or
     the pot declared not representable); a decisive reservoir block that this session has not seen.
   - **(b):** pause M2 until that supply case exists (the plan's "one refit" literally).
   - **(c):** widen bands on Zhang's numbers. Excluded by the coordinator's decision (bands never from Zhang).
2. **Whether to run at all after the power check** (sec. 12.1). If "Rt-D AND Rt-H" is nearly unattainable under omega,
   the run can only return "not separated". The alternative is a claim set without Rt-H, which is a design change for
   the re-review.
3. **The M1b dependency** (sec. 3.3): stated; no decision unless M1b runs first.

### 16.2 Carried for the supply case (review MAJOR 2, ready to apply when a feasible structure exists)

- C1F vs C1 as a required pair, with S-Cys, S-Xyl and S-TTCA all counting.
- S3 added if feasible.
- Cao's floor k_r >= 7e-3 /min at 90 C as a bound (log10 `k_ttca_cys`(145 C) >= -0.94 at 64.08 kJ/mol).
- A verdict conditional on a declared slope ceiling classed as numerical.
- Bands from Cao and the mechanism, never from Zhang.
- Item 2's evidence against: Cao's drift, Zhai's pH series, Zhang's mass balance.

## 17. For the register session (not written here)

**CF candidates:**
- (a) **Cao 2017 vs Zhang 2025b (same laboratory):** no single K holds both; apparent 60 L/mol at 90 C, pH ~4.4,
  against a ratio rising to ~18 600 L/mol at 120 C, pH 6.5. Extends CF-076. The screen (sec. 0B) adds that no banded
  form reaches Zhang's blank at all.
- (b) **Zhai 2021's printed TTCA-degradation Ea** of 80.99 kJ/mol is not reproduced by its own three rate constants
  (~35.5).
- (c) **CF-063 extended:** Zhang 2025b's printed cysteine fits misfit the cysteine panel's points panel-wide, the lysine
  arms included. Not a label swap (sec. 0A, P ~0.15).
- (d) **Zhang 2025b's high-lysine arms** plot free cysteine (~0.03 - 0.06 mmol/L) below the zero-signal point of the
  printed ELSD curve (0.223 mmol/L).
- (e) **The engine's in-situ pH in Cao's pot** barely drifts (3.86 -> 3.89) while the measured pH falls 4.67 -> 4.15.
  This is M1b's balance error (CF-084) showing on a route pot, and it makes the engine blind to K's pH slope there.
- (f) **Zhai 2021 Fig. 4b:** TTCA is lost ~10x faster at pH 7 - 8 than at 5.5, and no engine structure has a
  base-favoured TTCA exit (a new engine-error candidate for the atlas, E10-adjacent).

**Plan wording errors:**
- Row M2's "labelling in 3 labs (Cerny, Givaudan, Hofmann)": only Cerny's is labelling. Givaudan is de Roos 2005 /
  Laarhoven; Cerny 2003 is Nestle.
- Row M1's (M1.0) "in Cao 2017 free cysteine and xylose vanish within 30 min": they plateau at 15 - 20 % of the charge.
- "Zhai TTCA" is ambiguous; this design reads it as the 140 C fed-TTCA shapes.

**CR entries:** the floors; the noise transfers (Cerny 2007's 23 %; Cao's 10 %; the 9 % median for the Zhai cells); the
censored-MP3P detection ratio and response-factor band.

**EX:** this design's ID at freeze. The step-2 screen as its own EX (a pre-design feasibility record, verdict "no
supply structure feasible").

## 18. Open design questions

1. **Is the screen fair to S2 and S3?** Its held constants are C1's. Would a supply case need the screen re-run with
   the held constants freed (a cheap fit on the three signatures alone)?
2. **Cerny 2007's TIC ratio** needs a response-factor band and a detection ratio, both declared without a source. Is
   the MP3P absence then informative enough to identify `d_ddp`, or is the Cerny block only a thiamine-share block?
3. **Rt-H under omega:** 0.6 sd of separation at the shipped values. Keep it as a required claim unit (the review's
   decision), or move it to guards? The power check informs this.
4. **The T shape with a free level and omega per point** has ~3 degrees of freedom at 0.9 ln each. Is it informative,
   or should T be dropped and E11 left to W6?
5. **The robustness pair at S1** carries a supply form screened infeasible. Is a route verdict's agreement at an
   infeasible supply evidence of anything? (Default: it is a check of stacking, not evidence.)
6. **Guard adequacy at p > 0.05** is the coordinator's number. Alternatively, derive it from the plan's 90 % convention
   (p > 0.10).
7. **Holm at family-wise 0.10** across ~20 row groups. Is the in-sample guard then too weak to catch C1-style losses
   (36 nats over 52 rows)?
8. **`k_dpo_ddp` held:** freeing it adds a direction confounded with `d_ddp`. A no-refit sensitivity at its B60 Laplace
   +-1 sd?

## 19. Results

*(appended after the run; nothing above this line is edited after the first number exists)*

## 20. Revision log

- v1 DRAFT, 2026-11-05 (9bc55883): first version.
- **v2 DRAFT, 2026-11-05:**
  - step 1 (no label swap); step 2 (no feasible supply structure);
  - routes at instant supply with an S1 robustness pair;
  - Cerny 2007 into the objective; furfural barriers held; T as 140 C shapes; Givaudan reported;
  - omega primary, with the guard-adequacy rule;
  - band-integrated densities with plug-in agreement; PSIS and the minimum ESS; Holm;
  - the power check as a freeze prerequisite;
  - compute recomputed;
  - the register list. Coordinator's decisions applied as in 20.1.

### 20.1 The review of v1 (REDESIGN) and the decisions applied

| finding | decision (coordinator) | where in v2 |
|---|---|---|
| BLOCKER 1: omega missing from the primary covariance; guards veto against inadequate rivals; no power check | accepted | sec. 6 (omega primary), 8.1 (adequacy rule, "both inadequate"), 12.1 (power check, freeze prerequisite) |
| MAJOR 2: C1F vs C1 not required; Cao's floor only a check; the m ceiling ships; the evidence against not listed; the S2 band's brief error | accepted | the screen (0B) found no feasible supply structure, so the items are carried in 16.2 for the supply case; the evidence against is in sec. 2; bands never from Zhang |
| MAJOR 3: unidentified coordinates held (sharp predictions at arbitrary values); d_ddp start 2 from EX-040's pot D pass; route identifiability unargued | accepted | sec. 7.1 (band integration primary, plug-in required), 3.2 (seeded starts), 3.4 (identifiability, Cerny 2007 into the objective, pot D within-lab held-out, supersession) |
| MAJOR 4: B62 unjudged, yet its barriers fitted | accepted | sec. 1, 3.1 (furfural barriers held in every rival); T-FFT dropped |
| MAJOR 5: T's 120 C denominator shared with B60's folds; single-time ratio; Gigl ceiling | accepted | sec. 5 (140 C shapes; shared denominators declared; Gigl profile gates T) |
| MAJOR 7: Givaudan not simulable; Rt claim too lenient | accepted | sec. 5 (reported guard; Rt-D AND Rt-H) |
| MINORs: Holm; first-point form primary; ln units; PSIS and minimum ESS; probe with guards and the Cao / Zhai blocks; the forward-optimum form defined; polish budget accounting; switched-off rule incl. R0's `k_arp_dpo`; compute (profiles ~26 h serial, power check, checkpoints, caffeinate); Zhai 140 C 0.0081, Cao 58.6 | accepted | secs. 8.2, 6, 7.1, 4.4, 4.2, 9, 14, 13 item 11 |
| author's sec. 16.1 items (v1) | (1) decided: bands from Cao + mechanism only; (2) Givaudan demoted; (3) sign-off stated at ship time; (4) M1b dependency stated | 16.2, 5, 8.3, 3.3 |
