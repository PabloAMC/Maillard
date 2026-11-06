# Case M3b: the anoxic donor supply of ARP-Phe

**Status: PARKED, NOT FROZEN (coordinator decision, 2026-11-06, after the independent premortem of sec. 8).** The power
check is not run, and Figs. 2 - 3 stay blind under the NO-UNBLIND rule of sec. 9.

The case reopens only after four band-narrowing sweep candidates have been read: Kaufmann 2016, Robert 2005, Wu 2022
and Kawakishi 1991 (download queued). It is then reconsidered as **one joint M3 + M3b estimation spend** of the
figures.

Items 1 - 5 of [STATUS_AND_PLAN sec. 7](../../guides/STATUS_AND_PLAN.md) (the text of PR #80, branch
`plan-update-2026-11-05`, with its spent-data sentence) were written before any result of the change existed. The
integrity re-read CK-01 / CK-50 is folded into sec. 0. The premortem's fixes (sec. 8) are recorded but not yet
applied to secs. 3 - 5. Applying them is the first step on reopening. Items 6 - 9 are appended later and never edited
into this part.

*Plan step: an M0 nomination (misfit MF-29 of the misfit inventory, `docs/reference/misfit_inventory.md` on branch
`misfit-inventory`, draft PR #82) on the M3
layer. The layer is `src/m3_strecker_layer.py` on branch `m3-phase1-build` (PR #85, unmerged); its prereg is
`results/validation/kinetic_core_m3_strecker_layer_prereg.md` on that branch, cited below as "M3 prereg". Nominating
records: the argon diagnostic (PR #81, `results/validation/m3_argon_strecker_diagnostic.md` on `m3-argon-diagnostic`)
and the M3 phase-1 close-out (M3 prereg sec. 23 and `kinetic_core_m3_phase1_power`). Parent case:
[M3_strecker_aldehydes.md](M3_strecker_aldehydes.md).*

## 0. Seen by the author, blinding, and the CK-01 result

**Read while writing.** STATUS_AND_PLAN sec. 7 (main and PR #80) and the M0 and M3 rows of PR #80; the M3 case; the M3
prereg's header and secs. 2 - 8, 11 - 13, 22 - 23; the argon diagnostic report and its script header; the misfit-inventory entries
MF-06, -12, -23, -29, -30, -32, -33, -34, batches 0 and G, and its sec. 4 groups; CF-075 and CF-121; the
`strecker_routes` entry of `decisive_measurements.yml` (PR #79); the dossiers named in sec. 2 (relevant sections); the
OpenAlex titles and abstracts of sec. 2.8. The hofmann2000b dossier was read in full, **including its text-level
statements about Figs. 2 - 3** (1-DG the main dicarbonyl under argon, "in particular after longer reaction times";
glucosone "very low"; Phe liberated ~55 % at 320 min; under air + Cu 25 mM glucosone the main dicarbonyl, 1-DG and
3-DG "only low yields", Phe liberation "in the same order of magnitude"). Those statements are SEEN. To list the
paper's references for the backward sweep, the text layer of PDF pp. 4304 - 4305 was extracted. It showed the
reference column, one sentence of the mechanism discussion and the bare label "Figure 2.". No figure content, axis or
plotted value was shown.

**Not opened:** hofmann2000b Figs. 2 - 3 (the PDF figures), the keying directory, and any keyed Fig. 2 - 3 value. No
engine, layer or fit was run. The only computation is the band arithmetic marked **[dry]** (sec. 4.2; scratchpad
scripts `m3b_dry_bands.py` and `m3b_dry_bands_ktr03.py`, run through the wrapper). It reads no data, only the
declared bands, and no band, start or tolerance comes from it.

**Design metadata read from the keyed file (2026-11-06, at the coordinator's request; no value read).**
- **What was looked at.** On the local branch `key-hofmann2000b-fig2-fig3`, `git diff --name-only` against main lists
  the file names. Header and rows share one file, `data/lit/fit_rows/strecker_hofmann2000b_fig2_fig3.yml`, so it was
  not displayed. `git show` wrote it to a scratchpad file that was never printed.
- **The filters**, both through the wrapper:
  - one printed only the key paths and value TYPES;
  - one (`design_only.py`) printed only the pots' design fields (source_figure, atmosphere, metal, trap, t_c,
    value_basis, concentration_basis, the stated and plotted sampling times, vessel, analysis) and the set of
    (pot, species) cells with their plotted times and units.
- **Never printed:** `value`, `read_error`, `share`, the readers' y readings, the notes, the origin-cluster
  (near-zero) points, any flag count. The scratchpad copy was deleted afterwards.
- **What it says:**
  - **Both figures plot five series: ARP-Phe, Phe, 1-DG, 3-DG and glucosone.** **MGO and GO are not plotted.**
  - Arms: Fig. 2 argon (no metal; deoxygenation not stated); Fig. 3 air + Cu(II) 0.05 mmol (25 mmol/L; salt not named;
    headspace not stated). Both at 85 C, OPD 1.2 mmol in situ.
  - Times: stated 20 / 40 / 80 / 160 / 320 min, but plotted (keyed) **0 (ARP-Phe only), 20, 40, 80, 180, 320 min**.
  - y axis "amount [mmol/mmol ARP-Phe]", 0 - 1.00, no error bars; keyed as mmol/L = share x 500 (1.0 mmol in 2 mL;
    ~430 - 450 mmol/L if solids' volume or purity is counted).
- **Completeness of the cells.** Every (arm, species) cell has all five plotted times (20 / 40 / 80 / 180 / 320; ARP-Phe
  also at 0), so no cell is missing a time.
- **Per-cell flags were not read.** The filter never printed `disagreement_flag`, the reader counts or the
  origin-cluster entries. So none appeared, but their absence is not attested.
- The rest of this case uses only these design facts (secs. 4.1, 4.4).

**Spent (development evidence, never decisive here)** under the sec. 7 spent-data rule:
- hofmann2000 Tables 3 - 4 (the fed-donor calibrations), hofmann2000 Fig. 2;
- hofmann2000b Tables 1 - 3 and Fig. 1 (circles and squares), **including Table 2's argon arm, the datum this case
  explains**;
- EX-049's calibration inputs;
- everything PR #81 read: the constants of Ghiron 1988 Fig. 2, Martins 2003 Table 3 col. C, Martins 2005 Table 2 and
  Kanzler 2026 sec. 4.7; Granvogl 2006 (the s_DG1 band);
- the M3 phase-1 checks, screen and power records, and the regression in the M3 prereg sec. 23;
- the text-level Fig. 2 - 3 statements above: the qualitative ordering and the Phe end point at 320 min. CF-121 already
  used the end point. **The Figs. 2 - 3 data that remain blind are:**
  - the plotted yields at 20 - 180 min;
  - the absolute 320-min yields of each osone;
  - Phe at 20 - 180 min;
  - **the whole ARP-Phe series** in both arms; no ARP-Phe value is quoted in any text the author read.

> **CK-01 / CK-50 (integrity re-read of Table 2 and its Methods): result folded in 2026-11-06.** Source: draft PR #86,
> `results/validation/m0_integrity_batch1_strecker_sulfur.md` on branch `m0-integrity-batch1-ss` (read here as the
> coordinator's summary; the record itself was not opened). Table 2 stays SPENT development data.
> - **Integrity holds.** All 12 cells, the umol/mmol basis and the conditions agree with the page; argon PA + PAA is
>   0.09 mM. **Rd loses its basis / keying form.** What is left of Rd is a low single argon vial or an unprinted
>   condition, so Rd is demoted (sec. 3).
> - **Not printed:** vial size, headspace and the deoxygenation method. The scale reads 1 mL in the Methods and 10 mL
>   in the footnotes, and the phosphate 0.5 mmol/L (Methods) vs 0.5 mol/L (footnotes). Purity is ~90 %.
> - **Air, not O2-flushed**, is supported (P ~0.85). sigma_ox stays solved on an air arm; the O2-flushed reading is a
>   form at P ~0.15.
> - **CK-50, O2 available under air:** 85.6 r + 2.6 umol/mmol, with r = headspace / liquid volume. R1 works at
>   r >= 0.06. R2 needs r >= 2.9 - 8.2 (air) and >= 5.8 - 22 (air + Cu), and also k_L a >= ~1 /min against the
>   0.006 - 0.06 /min typical of an unstirred vial: P(feasible) ~0.2. This bears on the host (R1 primary), not on the
>   M3b rivals.
> - **Residual O2 under argon can only RAISE the argon datum.** It cannot explain an over-prediction, so o_Ar > 0
>   works against every rival here. It stays a form only for completeness.
> - **Carried as declared uncertainty:**
>   - the scale (1 vs 10 mL) and r, through the O2 budget of the air and air + Cu arms on which sigma_ox and c_Cu5 are
>     solved, hence k_anox;
>   - the phosphate print, through the fed pots' in-situ pH and k_d;
>   - the purity.
>   Each enters every rival's prediction of the argon datum (development) and the power check's truths (sec. 4.4
>   item 2). Verdicts must hold at both ends of r and of the scale.

## 1. The claim

**Under argon, most of the ARP-Phe that decays does not become a Strecker-active donor at the fed donors'
efficiency.** The M3 layer assumes that every anoxic exit (e1 -> 1-DG, e3 -> 3-DG, eF -> fragments -> MGO at fraction
m) feeds route A at the efficiencies calibrated on 100 mM fed donors. That gives ~0.5 % PA + PAA per ARP decayed,
against the ~0.1 % that Table 2's argon 0.09 mM allows (M3 prereg sec. 23).

The claim adds an explicit **Strecker-inert anoxic exit eI**, with share a_I, to the layer's anoxic exits. The steps
that change, all inside the standalone layer (nothing in `src/kinetic_core/`, nothing shipped):

| part | step (M3 prereg sec. 3.2) | the layer as built | M3b |
|---|---|---|---|
| anoxic exits | e1, e3, eF at a1, a3, aF ~ Dirichlet(1, 1, 1); k_anox = k_d (1 - sigma_ox) | every anoxic exit yields a Strecker-active donor (MGO at fraction m ~ U[0, 1]) | + eI at share a_I, renormalising (a1, a3, aF, a_I). Two physical forms, which differ on the decisive data: **Rb1 bypass** (Phe released, no free C6 osone made) and **Rb2 via free 1-DG** (the measured exit split; 1-DG made, then beta-cleaved to acetic acid + C4 faster than its Strecker step). Rb2 alone fails on development data (sec. 3.1); the measured form is carried as **Rab2**, Rb2 at Ra's slower anoxic decay |
| 1-DG fate | s_DG1 log-uniform [s_G, s_MGO x 10^0.82]; k0_DG1 [0.07, 1] /min (upper end numerical) | 1-DG a Strecker donor up to 6.6x MGO's efficiency | Rb2: k0_DG1 at Martins 2003's k9 (1.909 /min, a fit), with Davidek's measured acetic-acid yields as the check |
| trapped pot (85 C, OPD) | the sub-model of sec. 3.2: eps_X = k_tr / (k_tr + k_loss,X) | -- | eI makes no trappable osone (Rb1); in Rb2, 1-DG is trapped in competition with its beta-cleavage |
| glucose + Phe pot | 1-DG and FOX computed from the forced ARP | -- | eI lowers the 1-DG computed from the ARP (Rb1, Rb2) |

The fed-donor calibrations (hofmann2000 Tables 3 - 4), the route-B structure and the O2 exits are unchanged. The
route question (R1 vs R2 of the M3 case) is not asked here. Each M3b rival is grafted onto the layer's R1 as host
(primary) and onto R4 (form). R2's calibration is feasible in only ~12 % of draws and infeasible under the 20 mL
headspace form. R3 equals R1 on the fed pots.

## 2. All the evidence, by kind

"for / against" refers to the claim (Rb), with the rivals of sec. 3 named where a source bears on one. "Split": every
source below is development; none is in a frozen held-out split for E3, because the M3 case's held-out Strecker cells
carry no ARP and no donor.

### 2.1 Isotope labelling

| source | what it shows | bears on |
|---|---|---|
| Davidek 2006b (glucose + Gly 0.1 M, pH 8, 90 C, 15 h; argon / air / 18O2) | acetic acid **30.5 +- 0.4 (argon), 32.3 +- 3.3 (air), 33.5 +- 0.5 mmol/L (18O2)**, ~32 mol % of glucose, atmosphere-independent; ~1/3 carries O2 oxygen in 18O2 | **for Rb2**: an anoxic Strecker-inert C2 sink of the 1-DG route exists at the same size as in air. Against Ra in its "O2 accelerates the decay" form (weak: a 15 h end point, not a rate) |
| Davidek 2006a (13C glucose; 1-DG without amine, pH 5 - 8, 120 C) | acetic acid from C1 - C2 by hydrolytic beta-dicarbonyl cleavage of 1-DG; 1-DG gives 0.56 - 0.84 mol acetic acid per mol | for Rb2 (1-DG's main fate is Strecker-inert) |
| Zhou 2024 (Ala-xylose ARP + [15N]Ala, 120 C) | about half the aminoketone N comes from free Ala on ARP-made dicarbonyls | **against an all-inert anoxic supply** (pentose ARP, aerobic): route A does run on ARP-made donors |
| none | no label traces the argon PA of a fed hexose ARP to its donor | the gap the claim sits in |

### 2.2 Fed intermediates

| source | what it shows | bears on |
|---|---|---|
| **hofmann2000b Table 2, argon arm (SPENT; the misfit)** | fed ARP-Phe 0.1 M, 0.5 M phosphate, pH 7, 100 C, 120 min, closed vial: PA 0.6, PAA 0.3 umol/mmol | the datum every rival must meet; development only |
| Davidek 2002 (fed DFG 0.1 M, 90 C, pH-stat 5 - 8, water or 0.1 M phosphate; atmosphere not described) | per mol DFG decomposed: **glycine 0.46 - 0.88** (falling with conversion), **acetic acid up to 0.64**, formic up to 0.14. No dicarbonyl measured | **for Rb2** (most of the decay ends in a Strecker-inert acid); glycine < 1 is weak support for an amine-retaining form (Rb-ret) |
| Davidek 2005 (fed DFG, 0.2 M phosphate, 120 C, 240 min) | acetic acid **4 / 28 / 56 mol %** of DFG at pH 4 / 6 / 8; the dossier infers >= 67 - 70 % of DFG through 2,3-enolisation to 1-DG at pH 8 | for Rb2; at Hofmann's pH 7 the interpolation is ~40 mol % (arithmetic, a transfer) |
| Martins 2003 / 2003b (fed DFG 10 mM, 0.1 M phosphate, 100 / 120 C, pH 5.5 / 6.8; headspace not stated) | Model 2 split at 100 C, pH 6.8 (the diagnostic's S1): 1-DG 0.40, 3-DG 0.047, C_n fragments 0.42, reversion to hexoses 0.13; 1-DG -> acetic acid k9 1.909 /min (k8 = 0). At 120 C, pH 6.8 (Fig. 7): Gly ~7.6 mM, **MG ~2.9 mM**, acetic ~5.5, formic ~3.0 mM from 10 mM DFG | for Rb2 (a small 3-DG share; 1-DG goes to acetic acid). **Against Rb1**: the fragments make MGO, ~0.29 per DFG at 120 C, a Strecker-active donor. The enaminol-level reading (M3 prereg sec. 2 item 1: 3-DG route 0.15) is a form |
| Voigt 2009 (fed 1-DG + Ala 42 mM, 37 / 50 C, deaerated with DTPA) | acetic acid **650 mmol/mol 1-DG** at 24 h, 37 C; ~0.9 mol per mol 1-DG accounted for | for Rb2 (1-DG's measured fate in the presence of an amine is beta-cleavage) |
| Voigt 2010 (fed 1-DG + Ala, 37 / 50 C) | acetaldehyde, deaerated, 50 C, 12 h: ~19 mmol/mol 1-DG; beta-Ala halves it, so **~10 mmol/mol is Strecker** (OPD-trap recovery 27 %, correction not stated) | **the only on-disk 1-DG Strecker yield** (~1 %, comparable with DFA's fed 3-DG 0.84 %). The diagnostic recorded 1-DG's efficiency as unsourced; this is a same-class analogue across 50 C and amino acid |
| Smuda 2010 (1-DG + BOC-Lys, 37 C) | carboxylic-acid amides <= 1.5 mmol/mol 1-DG | against an amine-acylation form of Rb-ret |
| Deng 2021 (fed Met-ARP 0.2 M, 100 - 130 C; headspace not stated) | **88 % of the reacted ARP returns free Met** (0.1 M phosphate, 60 min); standing 3-DG <= 1.8 %, 1-DG <= 0.8 %, MGO <= 0.85 % at pH 7.5 | against Rb-ret at 100 C; standing levels are floors on flux, not yields |
| Hou 2026 (fed Gly-xylose ARP 5 mM, 120 C, unbuffered; headspace not stated) | **glycine released 0.40 -> 0.25 per ARP lost** (20 -> 120 min); peak 3-DX 15.5 %, 1-DX 4.8 % of the ARP | for Rb-ret (pentose ARP, 120 C); conflicts with Deng's 0.88 |
| hofmann2000 Table 4 (SPENT) | argon efficiencies (PA + PAA per mmol donor): MGO 15.8, GO 11.4, 3-DG 9.9, glucosone 9.4 umol/mmol | the efficiencies the incumbent carries to ARP-made donors |
| Weenen 2001 Table II | fed Glc-Phe ARP vs glucose + Phe, PA x4.2 at pH 7.0 (air) | no argon arm: not on the criterion |

### 2.3 Trapping

| source | what it shows | bears on |
|---|---|---|
| **hofmann2000b Figs. 2 - 3** (ARP-Phe 0.5 M, OPD 0.6 M, 85 C; argon / air + Cu 25 mM; 20 - 320 min; 1-DG, 3-DG, glucosone, Phe) | **the decisive data** (sec. 4). Text-level ordering SEEN (sec. 0) | every rival |
| Gobert 2009 / Glomb 2010 (glucose + Lys 42 mM, pH 7.4, 50 C, 7 d, continuous OPD) | deaerated vs aerated, mmol/mol glucose: **1-DG 23.5 vs 6.9, 3-DG 5.0 vs 4.6, glucosone 5.5 vs 44.5, MGO 4.0 vs 3.9**. Post-incubation OPD gives far less: 1-DG 0.2 vs 23.5 | a second laboratory for the **anoxic C6-osone pattern** (1-DG dominant, 3-DG small, MGO O2-flat). The continuous vs post-incubation contrast shows that 5 mM OPD outcompetes 1-DG's own loss at 50 C, a possible floor for k_tr (open question 3) |
| Nedvidek 1992 | in-situ OPD, sugars + beta-Ala: 1,4-dideoxyosones higher with alpha-amino acids | 1-DG-derived donors (route A on a 1-DG product) |

### 2.4 Measured intermediates and lifetimes

| source | what it shows | bears on |
|---|---|---|
| hofmann2000b Table 3 (SPENT) | in-pot k_d >= 0.0126 - 0.0144 /min (CF-121) | Ra (MF-32) |
| Davidek 2002 Fig. 2 (digitised) | DFG first order at 90 C, pH 7: 0.52 /h water, 1.01 - 1.13 /h in 0.1 M phosphate (0.009 - 0.019 /min) | the anoxic(?) decay rate of a fed hexose ARP (atmosphere not described) |
| Deng 2021 | k 0.0117 (unbuffered, pH 7.5), 0.020 - 0.022 /min (0.1 M phosphate), 100 C; Ea 75 +- 10 | the layer's k_d band source (spent through the M3 prereg) |
| Martins 2003 Fig. 7 | 1-DG and 3-DG are transients (peaks 0.06 and 0.045 mM from 10 mM DFG, 120 C, gone by 40 min) | in an untrapped pot both osones are short-lived; Figs. 2 - 3 capture what this design cannot |
| Gobert 2009 | 3-DG half-life 40 h and glucosone 8 h with Lys at 50 C | 3-DG is slow to leave: trapping efficiency for 3-DG ~1 (sec. 4.2) |

### 2.5 Mass balance

| source | what it shows | bears on |
|---|---|---|
| Martins 2003 Fig. 3 (molar balance of fed DFG, 100 / 120 C) | 78 - 85 % early, **~100 % within error late**: the identified products (Gly, sugars, deoxyosones, acetic / formic acid, MG, melanoidins) account for the decay | **against Rb1** as an unidentified bypass; consistent with Rb2, whose "inert" carbon is the measured acetic acid |
| Voigt 2009 | ~0.9 mol per mol fed 1-DG accounted for, mostly acetic acid | for Rb2 |
| M3 prereg sec. 2 item 11 | Martins: ~1.08 mol acid per mol DFG (from the pH drop); Davidek 2006a / 2005: 0.28 - 0.84 mol acetic acid per mol | Rb2; y_acid in [0.4, 1.1] (the trapped pot's pH) |

### 2.6 One factor varied inside one study

| source | factor -> result | bears on |
|---|---|---|
| hofmann2000b Table 2 (SPENT) | argon / air / air + Cu: ARP PA 0.6 / 5.5 / 13.8 | the misfit and the air excess that sigma_ox carries |
| **hofmann2000b Table 2, fed ARP-Phe vs glucose + Phe (SPENT)** | ARP-Phe over the glucose pot, PA: **1.5x under argon**, 3.9x in air, 5.3x in air + Cu (PAA 1.5 / 1.7 / 1.9x) | under argon, the fed ARP gives barely more aldehyde than the sugar pot that makes ARP only transiently (<= 4.8 mol %). Its anoxic donor supply per mol is far below what the layer's R0 gives it. Development evidence for the supply class; added after the premortem |
| Davidek 2006b | argon / air / 18O2: acetic acid unchanged | Rb2's sink is anoxic; weakly against Ra |
| Gobert 2009 | deaerated / aerated (DTPA rides with deaeration) | the anoxic osone pattern |
| Voigt 2009 / 2010, Smuda 2010 | aerated / deaerated (DTPA confound) | the 1-DG chemistry |
| Davidek 2002, 2005 | pH 4 - 8, water vs phosphate | the 1-DG route rises with pH: Rb2's share at pH 7 is a transfer |

### 2.7 Fitted time courses and records

| source | what it shows | bears on |
|---|---|---|
| Martins 2003 / 2005 | 1-DG -> AA 1.909 / 1.45 /min (fits; 1-DG unmeasured at that speed) | Rb2's k0_DG1 (above the layer's numerical upper end of 1 /min) |
| EX-049 | matched argon only because a lumped anoxic route (k_a0) was calibrated on these rows (M3 prereg sec. 23) | why the misfit was hidden before |
| argon diagnostic (PR #81) | H1 (first-order donor loss) closes 2 - 81 % without the fed counterpart, 14 - 20 % with it; H2 (donor efficiency) 8 - 37 %; H3 (PA sink) 1 - 10 %; S1 (Martins split) 8 - 11 %; H12_max passes but is inadmissible. Baseline argon shares: 3-DG 42 - 55 %, fragment MGO 22 - 25 %, 1-DG 6 - 14 % | Rc is disfavoured on spent data; S1 alone does not reach the datum because fragment MGO and s_DG1 stay free |
| M3 prereg sec. 23 | regression of the argon total (R^2 0.62 / 0.80) on m, a3, s_DG1 and the amount decayed; the k_d floor (0.005 /min) removes ~2x of the 5x | the nomination: Ra alone is short by ~2.5x on spent arithmetic |
| Piornos 2025 | one shared ARP decay for seven ARPs (malt) | the decay law (M3) |

**Left out, with reasons.** Route B vs route A (the M3 case, unchanged). Sulfur and cysteine ARPs (M1 / M2), except
Hou 2026's amine release. Lipid carbonyls (E20). Pentose ARP flavour studies without an argon arm or a yield per ARP
(Cui 2017 / 2018, Zhou 2022 - 2025, Zhan 2019 - 2025). Food matrices. Kanzler 2026 and Pfeifer 2013 (MGO self-loss)
enter only through Rc's constants, already spent.

**Against the claim, summarised.** Martins' mass balance closes with identified products, and its fragments make MGO
(~0.29 per DFG at 120 C). Zhou 2024 shows route A on ARP-made donors. Deng 2021 returns 88 % of the amine. Rb1 (bypass)
has **no measured share on disk**. Rb2's sink is well measured, but it acts on 1-DG, which carries only 6 - 14 % of the
baseline argon total. The large shares are 3-DG (Rb2 makes it small through the measured split) and fragment MGO
(nothing on disk makes it small at 100 C).

### 2.8 The citation sweep (sec. 7 item 2), 2026-11-06

**Keyed to the decisive criterion, not the topic.** Each work was asked: does it hold (a) a dicarbonyl (osone) yield
per ARP decayed, or per amine released, under argon or deaeration; (b) a fed hexose ARP decayed with and without O2;
(c) a Strecker-inert exit of an ARP or of 1-DG with its yield; (d) a trapping rate of 1-DG / 3-DG by OPD; or (e) the
amine returned per ARP decayed? Titles and abstracts only (blinding holds).

**Tool and seeds.** `scripts/tools/openalex.py sweep` from `origin/openalex-tool` (PR #70), copied to the scratchpad
and run through the wrapper with PYTHONPATH on this worktree, ~1 request per 1.1 s, anonymous.
- hofmann2000b (10.1021/jf000076e), forward and backward. 151 citing works. OpenAlex indexes only 6 of its
  references. The reference column of the PDF's text layer lists ~10: DFA's own papers, Kawakishi 1991 (below), Ledl &
  Schleicher (Angew. Chem., a review) and Strecker 1862.
- Forward only:
  - Martins 2003 Parts I and II (10.1016/s0008-6215(03)00173-3, 10.1016/s0008-6215(03)00174-5), the measured fed-ARP
    split;
  - Gobert 2009 (10.1021/jf9019085), the only other in-situ OPD trap with a deaerated arm;
  - Davidek 2002 (10.1021/jf025561j), fed DFG products per mol decomposed.

**Result.** 508 distinct works, 72 already in the registry by DOI, **436 unseen** (259 with an abstract, 177 title
only). Each unseen work got a decision (scratchpad table `m3b_sweep_screen_final.tsv`, not committed). The first pass
was by keyword, then every title that matched an ARP, osone, O2 or amine term (213) was read by hand:

| decision | n |
|---|---|
| food or matrix product study, no defined fed-Amadori or model pot | 187 |
| read individually: no fed-Amadori pot, no anoxic dicarbonyl yield per ARP decayed | 91 |
| non-supply chemistry (colour, antioxidant, acrylamide, HMF, furan, AGEs, ARP synthesis or application) | 69 |
| biomedical or other field | 43 |
| review, book or chapter without new pot data | 28 |
| computed (DFT, theory): excluded by the no-computed-barrier rule | 7 |
| passed to M3 (route question): Urdaneta Hartmann 2016 (10.1021/bk-2016-1237.ch001), ARP vs free amino acid by isotope enrichment in cocoa | 1 |
| thesis of the held Martins 2003 - 2005 papers (10.18174/121418) | 1 |
| **candidates** (below) | 9 |

Notable exclusions:
- Yuan 2024 (high-O2 Met-glucose pyrazines): a sugar pot, no ARP yield; G9 side.
- Granvogl 2012 (3-oxazolines): the Strecker step, not the supply.
- Zhang 2020 (Amadori derivatives in organic solvents, TGA).
- Akira 2005 (13C NMR of an ARP's fate): aerobic, 37 C, 28 d.
- Kaufmann 2012 / 2013 (1-DG and glucosone isomers by NMR): speciation, no rate; relevant only to k_tr's
  interpretation.
- Bruhns 2018 (2-deoxyglucosone): a Heyns route.

**Candidates (each unseen; role fixed now, before reading; P = chance it holds the named measurement).**

| candidate (authorship from the OpenAlex record) | DOI | role, fixed now | the measurement that would count | P |
|---|---|---|---|---|
| Kawakishi, Tsunehiro & Uchida 1991, Carbohydr. Res. 211:167 (Nagoya), autoxidative degradation of Amadori compounds with Cu | 10.1016/0008-6215(91)84156-9 | band source for Ra (ARP decay with vs without O2 / Cu) and sigma_ox; never decisive (another laboratory, probably 37 C) | ARP loss and glucosone yield under N2 vs O2 | ~0.5 |
| Greifenhagen, Frolov & Hoffmann 2015, Amino Acids (Leipzig), oxidative degradation of Nepsilon-fructosylamine peptides in heated water | 10.1007/s00726-015-1940-2 | evidence for Rb / Rb-ret (products per ARP with and without O2; the amine is retained by construction) | dicarbonyls or acids per ARP decayed, deaerated | ~0.4 |
| Kaufmann, Meissner, Pelke, Muegge & Kroh 2016, Carbohydr. Res. (TU / HU Berlin), structure-reactivity of ARPs vs ketoses | 10.1016/j.carres.2016.04.016 | band source for the anoxic exit split (a1, a3) of a hexose ARP | 1-DG / 3-DG per ARP decayed, atmosphere stated | ~0.4 |
| Wu, Fu, Ma, Dong, Du, Dong & Xu 2022, Foods (Dalian), beta-CD complex of Glc-Phe ARPs: thermal degradation | 10.3390/foods11091309 | band source for k_d of **ARP-Phe itself** (a second laboratory; MF-32) | an ARP-Phe decay constant in water at ~100 C, atmosphere stated | ~0.4 |
| Robert, Arce Vera, Kervella, Davidek & Blank 2005, Ann. N. Y. Acad. Sci. (Nestle), 17O-NMR of Maillard pathways | 10.1196/annals.1333.008 | band source for a1 : a3 from labelling (the top kind); Nestle, as Davidek | the 1-DG vs 3-DG route share from 17O | ~0.2 |
| Deng, Cui, Hussain, Hayat, Liu, Zhang & Ho 2024, Food Chem. (Jiangnan), targeted precursors in thermally processed Met-glucose ARP | 10.1016/j.foodchem.2024.142033 | evidence (Ho network, as Deng 2021; dicarbonyls per ARP; atmosphere probably unstated) | 3-DG / 1-DG / MGO per ARP decayed | ~0.3 |
| Mossine & Mawhinney 2010, Adv. Carbohydr. Chem. Biochem. 64 (Missouri; review) | 10.1016/s0065-2318(10)64006-1 | pointer only: its references may hold a primary anoxic fructosamine degradation not on disk | a cited primary study with argon vs air yields | ~0.3 |
| Glomb & Henning 2016, ACS Symp. Ser. 1237 ch. 9 (Halle; review of reactive fragmentation products) | 10.1021/bk-2016-1237.ch009 | pointer only (fragment yields, Strecker-inert cleavages) | a cited primary fragment yield per ARP | ~0.3 |
| Cui, Yu, Zhai, Feng, Chen, Hayat, Xu, Zhang & Ho 2021, Trends Food Sci. Technol. (Jiangnan; review of ARP formation and fate) | 10.1016/j.tifs.2021.06.055 | pointer only | a cited fed-ARP study under N2 / Ar | ~0.2 |

**None of these can be decisive.** They are other laboratories or other systems, and the decisive data are fixed in
sec. 4. They can only move bands, by an amendment registered before any M3b result exists. The chance that at least
one gives a band that narrows a confound of sec. 4.3 (k_tr, the split, Ra's rate) is **~0.6**. The download list was handed to
the coordinator for the private `data/articles/DOWNLOAD_QUEUE.md` (gitignored, main checkout). It is not committed.
The sweep did not seed Deng 2021, Voigt 2009 / 2010 or Davidek 2005 / 2006a / 2006b. Their forward sweeps overlap the
M4 sweep of 2026-11-05 (Voigt 2009, Davidek 2006a and Deng 2021 were seeds there).

## 3. The rivals

| id | explanation | what it changes in the layer | status on the spent argon datum (development, not a verdict) |
|---|---|---|---|
| **R0** | **the M3 layer as built** (R1 host; R4 form) | -- | **fails**: median 0.39 - 0.48 vs 0.09 mM; 5.7 % of R1 draws within 2x |
| **Ra** | **a slower anoxic ARP-Phe decay under argon** | k_anox from the fed designs' own measured-class rate (CF-121: ~0.006 - 0.010 /min at 100 C), not from the transfer-class band [0.0113, 0.132] x (1 - sigma_ox). Form: an O2 / metal-accelerated decay, k_d = k_anox + k_O2 [O2], so an O2-limited reflux pot decays slower | short: the floor form removes ~2x of 5x (M3 prereg sec. 23) |
| **Rb1** (the claim, bypass form) | an explicit Strecker-inert anoxic exit, inert **before** any free C6 osone; Phe released | + eI, a_I in [0, 0.9] (definitional; **no measured share on disk**; required agreement at both ends) | meets the datum by construction at a_I ~0.75 - 0.8 (full decay), ~0.6 combined with Ra (arithmetic from the 5x) |
| Rb2 (**a declared-failing control**, not a carried rival) | the anoxic decay follows the measured fed-DFG chemistry: a small 3-DG share and a large 1-DG share, with 1-DG beta-cleaved to acetic acid + C4 faster than its Strecker step | (a1, a3, aF, reversion) at Martins 2003 col. C (0.40 / 0.047 / 0.42 / 0.13; enaminol form 3-DG 0.15); k0_DG1 = 1.909 /min (Martins k9); s_DG1 from Voigt 2010 (~1 %, a transfer); fragment MGO at Martins' 100 C yield | **fails on development data** (sec. 3.1): fragment MGO alone gives ~0.10 - 0.14 mM against 0.09 measured, and the total is ~0.16 - 0.36 mM |
| **Rab2** (the claim, measured form, carried) | **Ra and Rb2 together**: the measured anoxic chemistry, at the fed designs' own slower anoxic decay | Rb2's constants plus Ra's k_anox (~0.006 - 0.010 /min at 100 C) | **[dry] reaches the datum only at its low end**: ~0.09 - 0.29 mM (sec. 3.1). Every constant is a measurement or a fit of one; nothing is free |
| Rb-ret (form of Rb1) | the inert exit **keeps Phe** (amine-bound products) | eI releases no Phe | as Rb1 on Table 2; differs on Figs. 2 - 3 (sec. 4) |
| Ra+b (form) | Ra and Rb1 together | both | as Rb1 at the lower a_I |
| **Rc** | **a first-order loss of the ARP-made donors**, carried into the fed calibration | k0 for 3-DG (Ghiron; Martins k5 + k6), MGO (Kanzler) and 1-DG (Martins k9), in the ARP pot AND in the Table 4 inversion | **fails**: 14 - 20 % of the log gap (`H1_max_consistent`) |
| **Rd** (demoted by CK-01) | **the datum**: Table 2's argon vial reads low (a single vial; an unprinted vial or deoxygenation difference). CK-01 found the cells, the basis and the conditions correct, and residual O2 would raise the datum, not lower it | none (the layer is right) | not testable on disk. The diagnostic gave p ~0.10; after CK-01 ~0.05 (my reading) |

Considered and not carried:
- **An ARP-pot-specific donor sink** (donor + ARP-Phe as amine, which the fed-donor pots lack). No constant or
  same-class analogue is on disk. On Figs. 2 - 3 it predicts what Rc predicts (the trap at 0.6 M competes with it).
- **Route A's efficiency failing to transfer to sub-mM donors** (the M3 prereg's second suspect). The diagnostic's
  algebra leaves only the amine level, worth <= 2x. On Figs. 2 - 3 it predicts what R0 predicts.

Both are kept in the record as members of the "downstream efficiency" class with Rc and Rd.

### 3.1 Rb2 on Martins' 100 C methylglyoxal (development data; pre-freeze step 3)

**The datum.** Martins 2003b Fig. 4 (dossier `martins2003b_extraction.md` sec. 8, digitised, +-0.08 mmol/L; DFG from
its Fig. 1, sec. 5, +-0.3 mmol/L) gives, for fed DFG 9.6 mM at 100 C, pH 6.8:

| t (min) | MG (mM) | DFG decayed (mM) | MG per DFG decayed |
|---|---|---|---|
| 30 | ~1.0 | ~5.2 | 0.19 |
| 60 | ~1.5 | ~7.7 | 0.19 |
| 120 | ~2.0 | ~9.4 | 0.21 |

MG accumulates ("no decrease in its amount was observed"). The authors find pH 5.5 -> 6.8 has "almost no influence".
At 120 C, pH 6.8 it is ~0.33 per DFG (3.2 / 9.8). **Fragment MGO is ~0.20 per ARP decayed at 100 C**, against the
incumbent's prior mean aF x m = 1/3 x 1/2 = 0.17.

**[dry] Rb2's argon total, by rescaling the diagnostic's baseline anatomy** (median argon total 0.39 - 0.48 mM; shares
3-DG 42 - 55 %, fragment MGO 22 - 25 %, 1-DG 6 - 14 %, the rest 6 - 30 %):
- fragment MGO x 0.20 / 0.17 = x1.2, i.e. 26 - 30 % of baseline: **0.10 - 0.14 mM, alone at or above the 0.09 mM
  datum** and near the 2 s band's top (0.139 mM);
- 3-DG x 0.047 / 0.33 = x0.14: 6 - 8 %;
- 1-DG beta-cleaved at Martins' k9: at most x0.5, so <= 3 - 7 %;
- the rest unchanged.
- **Total ~0.41 - 0.75 of baseline, ~0.16 - 0.36 mM: 1.8 - 4x high.**

**Verdict on development data: Rb2 alone does not survive.** It is kept as a **declared-failing control**: scored
beside the rivals, reported, never credited. Figs. 2 - 3 can still show its signature (D1 ~0.06), and that would
locate the 3-DG share whatever the argon datum says.

**The composite Rab2 (Ra + Rb2) is carried instead.** Ra's slower anoxic decay (~0.006 - 0.010 /min at 100 C, the
fed designs' own rate, CF-121) leaves 51 - 70 % decayed by 120 min instead of 87 - 99 %, x0.55 - 0.8. That gives
**~0.09 - 0.29 mM**, which reaches the datum only at the slow end of Ra with the fragment-MGO share at its measured
value. This is a rescaling of medians, not a run, and it assumes the anatomy's shares scale linearly. The comparison
(item 6) replaces it.

### 3.2 The OPD trap rate k_tr: a floor from Gobert 2009 / Glomb 2010 (pre-freeze step 2)

**The datum (held dossiers).** Glucose + Lys 42 mM, 0.1 M phosphate pH 7.4, 50 C, 7 d, deaerated. 1-DG trapped with
5 mM OPD present throughout (continuous) is **23.5 mmol/mol glucose** (`gobert2009_extraction.md` sec. 2;
`glomb2010_extraction.md` gives the same as 2.4 mol %). The standing 1-DG at 7 d, OPD added only at sampling, is **0.2
mmol/mol**. Fed 1-DG's half-life there is 0.5 h (Glomb 2010), so k_loss,DG1(50 C) ~0.023 /min.

**Arithmetic.**
- With a constant formation rate F over T = 10 080 min, the 1-DG formed is F T and the standing level is F /
  k_loss. Their ratio is k_loss T ~ 230.
- The trap captures eps1 F T, so trapped / standing = eps1 x 230. Measured: 23.5 / 0.2 = 117, so **eps1(5 mM OPD,
  50 C) ~0.5**, i.e. k_tr(5 mM, 50 C) ~ k_loss ~0.023 /min.
- Three biases push this estimate UP, so the floor is set well below it:
  - formation falls over 7 d (lysine fell 30 %), so the day-7 standing level under-states the mean and the true
    ratio is larger, giving a smaller eps1;
  - the authors warn that continuous OPD makes "false" osone from the Amadori product, inflating the trapped amount;
  - a short post-incubation trap may under-recover the standing 1-DG.
  **Allowing a factor 5 for them, eps1(5 mM, 50 C) >= 0.1**, so k_tr(5 mM, 50 C) >= ~0.0026 /min, and a second-order
  OPD constant >= ~0.5 /(M min).
- Hofmann's trap is 0.6 M OPD at 85 C. Linear in OPD (x120), with an Ea of the condensation >= 0 (x >= 1 from 50 to
  85 C), this gives **k_tr >= ~0.3 /min**. With a typical condensation Ea of 40 - 60 kJ/mol it would be x4 - 9 higher.

**Decision.**
- **k_tr band: [0.3, 10] /min primary**, a declared transfer (50 -> 85 C, 5 mM -> 0.6 M, glucose-Lys -> ARP-Phe).
- [0.1, 0.3) is a sensitivity labelled "below the Gobert floor". A verdict that needs it does not ship.
- The upper end, 10, stays numerical.
- Consequence: eps1 >= 0.45 at the layer's k0_DG1 band (<= 0.36 /min at 85 C) and >= 0.30 at Martins' k9 (0.69 /min
  at 85 C). The 1-DG trap confound shrinks from ~1 dex to <= 0.5 dex.
- The [dry] band arithmetic re-run on [0.3, 10] (`m3b_dry_bands_ktr03.py`) moves the predictions by <= 0.06 (sec.
  4.2).
- **Caveat carried into sec. 4.3:** the same "false osone from the Amadori product" artefact, at 0.6 M OPD against 0.5
  M ARP-Phe, can raise the trapped osones of Figs. 2 - 3 themselves.

## 4. The decisive data

### 4.1 What Figs. 2 - 3 measure, and the units

hofmann2000b Figs. 2 (argon) and 3 (air + Cu 25 mM): ARP-Phe 1.0 mmol + OPD 1.2 mmol in 2 mL phosphate pH 7 (0.5 M
ARP), 85 C. **Design (from the keyed file's metadata, sec. 0):**
- five series in each arm: ARP-Phe, Phe, and the quinoxalines of 1-DG, 3-DG and glucosone;
- plotted at 20, 40, 80, 180 and 320 min, plus ARP-Phe at 0;
- y in mmol/mmol ARP-Phe charged, axis 0 - 1.00;
- **MGO and GO not plotted.**

Because **ARP-Phe itself is plotted**, "per ARP decayed" is measured directly, and Phe liberated per ARP decayed
becomes a unit of its own (D6).

| # | unit (argon unless stated) | basis | why |
|---|---|---|---|
| **D1** | **Y3(t) = trapped 3-DG / ARP-Phe decayed**, t = 20 - 180 min (320 min as a form; per Phe liberated as a second form) | within-run ratio of two analytes (sigma_run cancels; the calibration coordinates do not) | the most robust unit. 3-DG leaves slowly, so its trap efficiency eps3 is 0.70 - 1.0 over the whole k_tr band. It measures a3, the anoxic exit that carries 42 - 55 % of the baseline argon total |
| **D2** | **Yos(t) = (1-DG + 3-DG + glucosone) / ARP-Phe decayed** | as D1 | bounds the free-C6-osone flux, and so the room for a bypass (Rb1). Confounded with 1-DG's trap efficiency eps1 = k_tr / (k_tr + k0_DG1(85)), now >= 0.30 - 0.45 with k_tr in [0.3, 10] (sec. 3.2) |
| **D3** | **Y1 / Y3** | ratio of two quinoxalines (1-DG N(0, 0.20), 3-DG N(0, 0.10)) | the split a1 : a3 x eps1 / eps3. The seen ordering (1-DG main) already fixes its sign, so only its size is blind |
| **D4** | **ARP-Phe(t) / ARP-Phe(0), t = 20 - 320** (blind in full) | within-run course, one analyte (the calibration cancels; the analyser CV does not) | **the anoxic decay itself** (Ra vs R0). Confounded with Ea (85 -> 100 C), order (0.5 vs 0.1 M) and the acid load (0.07 - 0.7 dex, worse under the "0.5 mmol/L" print) |
| **D5** | **ARP-Phe decayed (air + Cu) / ARP-Phe decayed (argon)** at each t | cross-run ratio (each arm keeps its sigma_run) | the O2 / Cu acceleration of the decay (Ra's O2-accelerated form). Confounded with the Cu law (5 -> 25 mM, saturated or linear) |
| **D6** | **Phe liberated / ARP-Phe decayed** (argon; air + Cu as a second arm) | within-run ratio, one analyser (ninhydrin, two channels) | **amine return**: Rb-ret predicts < 1 by its share; Rb1, Rab2 and R0 predict ~1 minus (1 - r_F) aF. Free Phe in the ~10 % impurity biases it up early (imp_Phe). The 320-min Phe point is seen, its ARP partner is not |
| guard | glucosone / Phe in air + Cu (M3's U3b) | -- | M3's route unit. Shared blinding. M3b may not move it outside its band |

### 4.2 What each rival predicts (without reading the data)

The **[dry]** numbers are the declared bands pushed through the trap competition alone, per ARP decayed:
- the Dirichlet shares, r_F and k0_DG1 as the M3 prereg declares;
- k_tr on [0.3, 10] (sec. 3.2);
- eps at 85 C (Ea 75 kJ/mol).

There are no kinetics and no conditioning on Table 2. Values are the median with 5 - 95 % in brackets. The real power
check (sec. 4.4) replaces them. Script: `m3b_dry_bands_perARP.py`.

| rival | D1 Y3 | D2 Y1 + Y3 | D3 Y1 / Y3 | D4 ARP-Phe decay (argon) | D5 | D6 Phe / ARP decayed | Table 2 argon (spent) |
|---|---|---|---|---|---|---|---|
| R0 as built | 0.29 (0.025 - 0.77); 0.14 (0.012 - 0.37) once conditioned on the seen ordering (Y1 > Y3) | 0.67 (0.21 - 0.94) | any (P(Y1 > Y3) 0.48 before conditioning) | first order at k_anox(85) = 0.36 x k_d (1 - sigma_ox), k_d in the transfer band | per sigma_ox and the Cu law | 0.89 (0.48 - 1.0): the eF exit releases Phe only at r_F ~ U[0, 1] | fails |
| Ra | as R0 (shares unchanged) | as R0 | as R0 | **slower than R0** at the same Ea; first order, or an O2-depletion shape | **larger than R0** (O2 / Cu accelerate the decay) | as R0 | short ~2.5x |
| Rb1 | 0.13 (0.01 - 0.53) | **0.31 (0.07 - 0.75)**, low for any k_tr | as R0 | as R0 | as R0 | 0.95 (0.68 - 1.0) | met by a_I |
| Rb-ret | as Rb1 | as Rb1 | as R0 | as R0 | as R0 | **0.44 (0.11 - 0.87)** | met by a_I |
| **Rab2** | **0.047 (0.044 - 0.047)**; enaminol form ~0.15 | 0.33 (0.18 - 0.42) at Martins' k9 | **large** (~6 x eps1 / eps3) | **as Ra** (slower than R0) | as Ra | 0.79 (0.60 - 0.98) | reaches it only at its low end (sec. 3.1) |
| Rb2 (control) | as Rab2 | as Rab2 | as Rab2 | as R0 (Davidek 2002's rate at 90 C is inside R0's band) | as R0 | as Rab2 | fails (sec. 3.1) |
| Rc | as R0 (the trap outcompetes a first-order loss: eps3 >= 0.70) | as R0, slightly lower | as R0 | as R0 | as R0 | as R0 | fails |
| Rd | as R0 | as R0 | as R0 | as R0 | as R0 | as R0 | met only as a low single vial (CK-01 found no integrity error) |

Read as classes, the rivals separate on Figs. 2 - 3 like this:
- **Supply (Rb1, Rab2): low D2.** They split on D1 / D3: Rab2 has a tiny 3-DG share and much 1-DG; Rb1 scales every
  osone down.
- **Decay (Ra, Rab2): D4 and D5.** Now measured on the ARP-Phe series itself.
- **Amine retention (Rb-ret): D6.** R0's free r_F overlaps its lower half.
- **Downstream efficiency (Rc, Rd, the two not carried): indistinguishable from R0** on every unit.

### 4.3 Can Figs. 2 - 3 decide? What they can and cannot

- **They can test the supply class against the downstream class.** D1 is nearly free of the trap confound. If the
  trapped C6 osones per ARP decayed are high (D2 near R0's 0.67), there is no room for an inert bypass: Rb1 is
  refuted, and the misfit's cause is downstream of the osone (Rc, Rd, 1-DG's fate) or in the decay (Ra). If D1 is
  ~0.05, the measured split (Rab2) is favoured over R0's Dirichlet prior.
- **The decay is now measured, not inferred (ARP-Phe is plotted, and its series is blind).** D4 tests Ra / Rab2 against
  R0 on the anoxic decay itself, and D6 tests Rb-ret on the amine return. The seen 320-min Phe point (~55 %) becomes a
  lower bound on ARP decayed under any rival that releases Phe. The Ea, order and acid-load transfers still sit between
  85 C, 0.5 M and 100 C, 0.1 M.
- **They cannot confirm Rb1 against R0 itself.** R0's own freedoms contain an inert corner: aF large with m small, or
  a1 large with k0_DG1 at its end. Conditioned on the spent argon datum, R0 moves there (the M3 prereg sec. 23
  regression: the closest draws have a3 ~0.06, m ~0.07). The fragments are not trapped, and **MGO / GO are not plotted**
  (design, sec. 0). So "ARP decayed, Phe released and nothing trapped" reads the same whether the untrapped carbon is
  Strecker-inert (Rb1) or inert fragments (R0's corner). Their difference on Figs. 2 - 3 is prior mass, not mechanism.
  **P(Rb1 vs R0 ends "not separated") ~0.85.**
- **Fragment MGO is not measured** by Figs. 2 - 3. It is the second-largest piece of the misfit (22 - 25 % of the
  baseline argon total), and Martins puts it at ~0.20 per ARP decayed at 100 C (sec. 3.1). This is now the case's main
  unmeasured piece. It is also why Rb2 fails on development data while Rab2 survives only at its low end.
- **The trap may make osone (new caveat, sec. 3.2).** Glomb's group warns that continuous OPD produces "false"
  glucosone from the Amadori product. At 0.6 M OPD against 0.5 M ARP-Phe, part of the trapped osones (glucosone
  surely; 1-DG / 3-DG possibly) may be OPD-made, not ARP exits. This biases D2 up, i.e. against Rb1 (a spurious
  refutation), and biases M3's U3b. No size is on disk. It must enter the power check as a form: an OPD-induced osone
  flux proportional to [OPD][ARP], from 0 up to the Gobert-floor bias allowance (factor 5).
- **Rc vs R0 and Rd vs R0 cannot be separated** on any unit (identical predictions). Rc is already disfavoured on
  spent data (14 - 20 % closed). CK-01 found no integrity error, so Rd is left as an untestable low vial; Figs. 2 - 3 cannot speak to it.
- **Every transfer is from 85 C, 0.5 M ARP, 0.6 M OPD, a falling pH and 25 mM Cu to 100 C, 0.1 M ARP, pH 7.0 and 5
  mM Cu.** The exit split itself moves with temperature and pH (Martins: k4 / k7 0.37 at 100 C vs 3.9 at 120 C, pH
  6.8; Deng: 3-DG ~7 % at pH 5.5 vs <= 1.8 % at pH 7.5). So **a3 at 85 C under acidification is not a3 at 100 C,
  pH 7.**
- Under M1.1's omega (0.419 dex per printed value; a ratio carries 2 omega^2, 0.59 dex), the largest [dry] separation
  is ~0.8 dex (Rab2 vs unconditioned R0 on D1). The conditioned ones are ~0.4 dex, and D6 (Rb-ret vs the rest) ~0.3
  dex. With five correlated points per arm, **the expected power is low** (the M3 power record found the same for U4
  under omega).

**Overall (revised 2026-11-06 after the design read, the k_tr floor and Martins' 100 C MGO):**
- **P(at least one pair of carried rivals separates at usable power) ~0.4**, up from ~0.25. Reasons: ARP-Phe is
  plotted, so D4 and D6 exist; the trap confound is narrower.
- **P(they decide the claim against R0) ~0.15**, up from ~0.1. Rab2 vs R0 (D1 with D4) ~0.3; Rb1 vs R0 ~0.1.
- P(Ra vs R0 separates on D4) ~0.25: the transfers dominate.
- P(Rb-ret vs the rest on D6) ~0.3.
- The pairs most likely to clear the power check are Rab2 vs R0 on D1 + D4 and Rb1 vs R0 on D2, and only in the form
  that does not condition R0 on the argon datum.

### 4.4 What the power check must simulate (before the freeze)

1. **Truths per rival** (R0, Ra, Rb1, Rb-ret, Rab2, Rc, and Rb2 as the failing control; R1 host, R4 form): draws from each rival's bands **conditioned
   on the development set only**:
   - Table 2 (every arm), Fig. 1, Table 3, hofmann2000 Tables 3 - 4;
   - the seen text constraints on Figs. 2 - 3, since a truth that contradicts seen text is not a possible world: 1-DG
     above 3-DG and glucosone "very low" under argon; glucosone the main osone under air + Cu, with 1-DG and 3-DG low;
     Phe ~55 % at 320 min under argon, the same order under air + Cu;
   - in **two conditioning forms**, with and without Table 2's argon cell. A verdict that flips between them does not
     ship.
2. **The trapped pot as the M3 prereg sec. 3.2 - 3.4 builds it**:
   - 85 C via the exits' Ea [65, 85];
   - 0.5 M ARP with the order n [1, 2] (R4 host);
   - OPD trapping, k_tr on [0.3, 10] /min (sec. 3.2) with both ends as required agreement, and [0.1, 0.3) as a
     sensitivity (a verdict that depends on k_tr's ends does not ship);
   - the OPD artefact form: an OPD-induced osone flux from the ARP, from 0 up to the factor-5 bias allowance of sec. 3.2;
   - eps_X per species (3-DG with Ghiron's and Martins' losses; 1-DG with k0_DG1 [0.07, 1] and Martins' 1.909);
   - the in-situ pH from the acid load, y_acid [0.4, 1.1], under **both phosphate prints** (0.5 M and an unbuffered
     form), with the exit split's pH dependence (Rb2: Martins' pH 5.5 vs 6.8 columns as the two ends);
   - o_Ar {0, 0.1} (it can only raise the argon datum, CK-50);
   - the Table 2 scale (1 vs 10 mL) and headspace ratio r (CK-50: r >= 0.06 for R1; unprinted), through the O2
     budget of the air / air + Cu arms that set sigma_ox, c_Cu5 and so k_anox; both ends required agreement;
   - the Cu law (5 -> 25 mM saturated or linear);
   - arp0 [90, 100] % and imp_Phe [0, 1] (free Phe at t = 0 inflates early Phe liberated).
3. **The measurement model on the keyed design** (sec. 0: five series per arm, ARP-Phe, Phe, 1-DG, 3-DG and glucosone;
   20 / 40 / 80 / 180 / 320 min, ARP-Phe also at 0; y 0 - 1.00 mmol/mmol; no MGO / GO):
   - read error at the Fig. 1 keying level (<= 0.011 umol/mmol plus 0.03 dex), replaced at scoring by the adjudicated
     per-point error;
   - left-censoring at the axis resolution (~0.003 - 0.005 mmol/mmol; the x2 form);
   - quinoxaline calibration coordinates: 3-DG and glucosone N(0, 0.10), 1-DG N(0, 0.20) (form 0.30); Phe N(0, 0.05)
     plus the analyser CV [2, 10] %;
   - sigma_run 0.075 dex per arm, plus the between-points form;
   - **omega 0.419 per printed value in the primary covariance (M1.1), with the measurement-only covariance reported
     beside it**;
   - the ARP-Phe series with the analyser CV and the purity coordinates (arp0, the ~430 - 500 mmol/L charge reading).
4. **The pairs and controls:**
   - Rb1 vs R0 (D2);
   - Rab2 vs R0 (D1, D3, D4);
   - Rb1 vs Rab2 (D1, D3);
   - Ra vs R0 (D4, D5);
   - Rb-ret vs Rb1 and vs R0 (D6);
   - Rc vs R0, expected 0: a negative control for the check itself;
   - Rb2 scored as the declared-failing control.
5. **Scoring:** R = 200 synthetic sets per truth. The densities use M3 v2.1's estimator (reweighting, then the
   validated adaptive-tempering SMC). An unresolved density reads as a tie. Then the binomial demotion rule of M3
   prereg sec. 12.2: power must beat false credit by a one-sided exact test at 0.05, and power below 0.5 is labelled
   low.
6. **The bound profile:** k_tr's ends, k0_DG1's ends, a_I's ends, the Cu law, the OPD artefact form, and the
   conditioning form.

**Compute:** ~1 - 2 h work to adapt the M3 power harness. At 1 core, ~3.5 - 7 h unattended for seven truths x R = 200
under the primary and three forms, by the M3 record's per-density cost. It cannot run before PR #85's layer is
importable.

## 5. Tolerances from the data

| comparison | noise as the source states it | where stated |
|---|---|---|
| D1 - D6 read error | two blind readers, adjudicated by a third: per point mean(A, B), error max(err_A, err_B, abs(A - B) / 2); the axis read resolution on a 0 - 1.00 mmol/mmol axis sets the censoring level | M3 prereg sec. 7 (the rule; values unknown until scoring); axis from the design read (sec. 0) |
| D1 - D6 calibration | quinoxalines without a printed response factor: 3-DG, glucosone N(0, 0.10) dex; 1-DG N(0, 0.20), form 0.30; Phe and ARP-Phe by the same amino-acid analyser N(0, 0.05) each plus CV [2, 10] % (in D6 the shared analyser scale partly cancels: declared form rho 0.5) | M3 prereg sec. 6 |
| replicates | none printed (hofmann2000b: "not stated; no uncertainties printed") | dossier sec. 2 |
| run to run (D5 and every cross-pot use) | sigma_run 0.075 dex per arm: the 7 - 18 % run scatter of this recipe (Table 2 vs Fig. 1 at 120 min) | CF-074, CF-075 |
| model discrepancy | omega 0.419 dex per printed value (H0 EX-059 warm set; 0.390 sensitivity), the M1.1 rule | EX-059; M1.1 |
| purity | arp0 [90, 100] %; imp_Phe [0, 1] | hofmann2000b Methods ("about 90 %") |
| Table 2 (conditioning only) | print +-0.05 umol/mmol (0.036 / 0.072 dex on argon PA / PAA); s_rep 0.025 dex (Granvogl 2006 RSD) | `strecker_hofmann2000b.yml`; M3 prereg sec. 2 item 7 |
| constants from other laboratories (Martins, Davidek, Voigt) | their own HPD or SD where printed (Martins: 95 % HPD; Davidek 2006b: +-0.4 - 3.3 mmol/L; Voigt: SD < 10 mmol/mol); the transfer to 85 / 100 C and to Phe is declared, never a noise term | dossiers |
| between laboratories | not applied: every decisive value is DFA's. H0's tau enters only if a sweep candidate is admitted as a second laboratory's decisive datum (by amendment) | -- |

No round-number threshold is used. The diagnostic's +-2 s band is its own pass rule, not a tolerance here.

## 6. What the M3b rivals say about the G11 and G8 members

The protocol's best evidence for a mechanism is that it predicts beyond the misfit it came from. Checked on
**already-seen data only, by direction, with no run**:

| member | what it is | rival that speaks to it | direction | reading |
|---|---|---|---|---|
| **MF-32** (G11) | fed ARP-Phe decays slower (~0.006 - 0.010) than in-pot (>= 0.0126); the fed early / late shape drifts (early k 0.014, late 0.0078) | **Ra in its O2-accelerated form** | a 0.1 M ARP reflux pot (Fig. 1) and the argon trapped pot have less O2 per ARP than the ~5 mM ARP glucose pot, so they decay slower (right sign); O2 depletion in the fed reflux pot gives early-fast, late-slow (right sign) | **the only cross-misfit prediction found, and not unique.** The acid load (G1, MF-32 cand. (b)) predicts both signs too. A second-order law (R4) predicts the drift but the wrong sign for fed vs in-pot (the dilute glucose pot would decay slower). It does not decide Ra; it is why Ra is carried |
| **MF-06** (G11 / G1) | cooled pH too acidic (buffered Cys pots); wrong shape in unbuffered ARP pots; cand. (a) the amine is returned | Rb-ret vs Rb1 / Rb2 | they differ on the amine returned per ARP decayed, the quantity cand. (a) turns on: Rb1 / Rb2 return it (Deng 0.88, Davidek 0.46 - 0.88); Rb-ret keeps it (Hou 0.25 - 0.40) | speaks to the same quantity, but the sources conflict by ~0.4 dex and MF-06's pH data need the CK-12 machinery (a run). No direction can be checked on seen data. **New:** D6 (Phe liberated per ARP-Phe decayed, blind) measures that quantity in a fed hexose ARP. Whatever it shows also bears on MF-06 cand. (a), but only for the DFA trapped pot, so a cross-misfit claim needs CK-12. Rb2's acetic-acid exit adds acid per ARP, against MF-06's "too acidic" unless the amine return offsets it |
| MF-29 itself | -- | Rb1 by construction; Ra short; Rb2 fails on fragment MGO (sec. 3.1); Rab2 only at its low end; Rc fails | -- | the nomination (spent) |
| **MF-30** (G8) | the B18 Strecker constant ~2.3 dex low; glucose-pot PA lag | Rb1 / Rb2 lower the 1-DG computed from the ARP in the layer's glucose pot | lowers late PA a little; the sign of the effect on the lag needs a run | no seen-data test without a run. A run would spend nothing new (hofmann2000 Fig. 2 is spent) and could be a development diagnostic |
| **MF-34** (G8) | total pyrazines ~2.9 dex low: a supply deficit | every M3b rival lowers the ARP-made donor supply | **opposite sign**: carried into the network (phase 2), it would deepen the deficit | G8's deficits cannot share M3b's piece. They need another (amine identity, the CF-095 adduct sink, the B18 constant) |
| **MF-12** (G8) | glucose + Cys gives almost no MFT; cand. (a) a Glc + Cys -> 1-DG step | Rb2 (1-DG the dominant anoxic osone) is consistent in spirit; Rb1 (bypass) is against | qualitative; another amine (Cys forms thiazolidines) | no test |
| **MF-23** (G8) | glucose alone, autoclave: osones low | none (amine-free pot) | -- | -- |

**Summary.** No M3b rival explains a second member on seen data. Ra predicts MF-32's two signs but shares that with
the acid load. Every supply rival runs against G8's deficits. So the M3b piece is MF-29's alone, or MF-29's plus
MF-32's (G11) if Ra holds. The trunk's MGO excess (B21: Xia 7.8 vs <= 0.5 mM) has the sign Rb predicts, but it is the
trunk's glucose MGO, not an ARP-made donor, and it is not counted.

## 7. Open questions before the freeze

1. **Which series Figs. 2 - 3 plot: done** (sec. 0). ARP-Phe, Phe, 1-DG, 3-DG and glucosone are plotted; MGO and GO
   are not. The fragment-MGO gap stays open. The cheapest design that would close it is a W6 trapped-ARP pot with MGO
   quantified.
2. **CK-01 / CK-50: done** (sec. 0). Table 2 holds; Rd is demoted. The scale, headspace and phosphate prints stay open
   and are carried as declared uncertainty.
3. **k_tr floor: decided** (sec. 3.2). k_tr in [0.3, 10] /min primary, [0.1, 0.3) a sensitivity. It brings in a new
   form, the OPD artefact.
4. **The 1-DG Strecker efficiency.** Voigt 2010's ~10 mmol/mol (50 C, Ala) is an on-disk analogue that the diagnostic
   did not use. The M3 layer's s_DG1 band (up to MGO x 10^0.82) should be re-examined in the M3 record, not here.
5. **Rb2's fragment-MGO share at 100 C: done** (sec. 3.1). Martins 2003b Fig. 4 gives ~0.20 MG per DFG decayed. Rb2
   alone fails on development data and is kept as a declared-failing control. The composite Rab2 is carried.
6. **Merge order.** The layer (PR #85) and the diagnostic (PR #81) are unmerged. This case cites them by branch.
7. **CF entry.** No CF exists for MF-29 (inventory sec. 7). The register session should file the "every M3 rival
   over-predicts argon PA + PAA" conflict, citing this case (one writer per register; not filed here).

**Owner-level decisions: none.** M3b changes only the standalone layer and moves no shipped prediction.

## 8. Premortem (2026-11-06)

An independent premortem, run on this case after pre-freeze steps 1 - 3, gave the ranked failure modes below. The
coordinator parked the case on them. Each fix is to be applied to secs. 3 - 5 on reopening, before any power check.

| # | failure mode | fix |
|---|---|---|
| 1 | **The carried rivals are nested in R0.** Rb1, Rb-ret, Rc and Rd equal corners of R0 on D1 - D6. Conditioning on Table 2's argon cell collapses R0 toward Rab2 | declare those pairs **"cannot decide"**. Keep **one primary pair, Rab2 vs R0 conditioned on the argon cell**, with Holm across the rest. Add prior-sensitivity forms. Prefer an **estimation framing**: R0's posterior over a3, the inert share and k_anox, rather than a rival verdict |
| 2 | **One exit Ea carries the 85 C split to 100 C.** Martins' k4 / k7 (0.37 at 100 C vs 3.9 at 120 C, pH 6.8) implies ~x7 in the split between 100 and 85 C, and Deng's pH dependence (3-DG up at lower pH) runs the other way under the trapped pot's acidification | band the split's temperature and pH dependence for every rival, not only Rb2. Redo the [dry] bands |
| 3 | **OPD side effects.** The OPD artefact must scale like k_tr (x120 from 5 mM to 0.6 M). OPD may catch pre-osone flux, biasing D2 against Rb1. Cu / OPD / phosphate interactions confound D5 | make D1 and D2 one-sided; make D5 report-only; split Rb1 by where its inert branch leaves (before or after a trappable intermediate); add an OPD-perturbed decay form |
| 4 | **Double counting.** Ra's k_anox band comes partly from this trapped pot (CF-121, the seen Phe end point) | condition every scoring density on the seen text as numbers (ordering; Phe ~55 % at 320 min); drop the 320-min points from D4 and D6 |
| 5 | **Early points are poorly read.** At 20 - 40 min the yields sit near the axis floor: ~50 - 200 % relative error. The plotted time is 180 where the Methods say 160, and t = 0 normalises ARP-Phe | forms for 160 vs 180 min and for the t = 0 normalisation (the keyed 0-min ARP vs the nominal charge); read error relative near the floor |
| 6 | **The power-check spec was optimistic** | Holm; misspecified truths; a prior-mass negative control (Rb1 vs R0's inert corner); SMC conditioning wherever ESS < 100; simulate every required-agreement form; never count infeasibility as credit. Realistic cost **10 - 20 h** (not 3.5 - 7 h) |
| 7 | **Opportunity cost** | P(the primary pair passes after these fixes) **~0.15 - 0.2**. Spending blind data on that is poor value alone; hence the joint M3 + M3b spend (status line) |
| 8 | **Rd was demoted too far** | keep **Rd >= 0.10**. A result "R0 favoured" reads as "the downstream class or Rd", never as "R0 confirmed" |
| 9 | **Noise forms** | D6's shared-analyser rho 0.5 becomes required agreement (rho 0 and 0.5); an HPLC CV band for the quinoxalines; power reported under omega and under a DFA-internal discrepancy term |
| 10 | **Evidence and design attestations** | Table 2's glucose + Phe argon contrast added as spent development evidence (sec. 2.6). Every (arm, species) cell is complete in time; per-cell flags were not read (sec. 0) |

**What the premortem changes in sec. 4.3's estimates:** P(primary pair passes) ~0.15 - 0.2 replaces the per-pair
figures, which counted nested pairs as separable.

## 9. The NO-UNBLIND rule (pre-declared, 2026-11-06)

**hofmann2000b Figs. 2 - 3 stay blind** unless the fixed primary pair (Rab2 vs R0 conditioned on Table 2's argon cell,
after the sec. 8 fixes) passes a power check that simulates **every** required-agreement form of secs. 4.4 and 8. That
check must apply the binomial demotion rule, Holm, and never count infeasibility as credit. If it does not pass, the
case closes **"undecidable on disk"** and the figures stay unspent for a later use (the joint M3 + M3b estimation, or a
W6 design). No unit, form or pair may be added after the power check to rescue it.
