# Grosch & Zeiler-Hilgart 1992 — EXTRACTION (MFT, FFT and bis(2-methyl-3-furyl) disulfide from thiamine, ribose, cysteine, GSH, H2S at meat-like concentrations, 0.2 M Na-pyrophosphate pH 5.7, 2 h boiling; SIDA; the same pots refluxed-then-extracted (Method I) vs steam-distilled continuously (SDE, Method II): **SDE gives 4-5x more MFT**)

### The "4-5x SDE effect" the corpus cites second-hand (bolton1993, evidence map 2.11). Read at source it is NOT an extraction-recovery artefact: the deuterated standards were added after boiling, so Method I measures what SURVIVES 2 h in the pot. The ratio is therefore a within-study (kind 6) measure of an in-pot MFT sink at 100 °C, which SDE avoids by stripping MFT as it forms.

**Read 2026-10-31 from `data/articles/grosch1992.pdf`** (10 pp., ACS Symp. Ser. 490, ch. 14, pp. 183-192; scanned
book chapter with a garbled text layer). **Every number below was read from 100 dpi page renders** (scratchpad
`dos7/gro-01..10.png`); the text layer splits words and loses table structure. One figure (Fig. 1, a scheme). No SI.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 1. What this paper answers for the model

| what the model needs | what Grosch & Zeiler-Hilgart give | tag |
|---|---|---|
| SDE vs direct extraction of MFT (evidence map 2.11, "cited only") | **Table VIII, same pots, SIDA:** thiamine 30 µmol: 4.7 ± 1.2 (Method I) vs 25 ± 5 µg (Method II); thiamine + cysteine 325 µmol: 33 ± 3.7 vs 136 ± 0 µg; ribose 670 + cysteine 325 µmol: < 0.1 vs 0.1 µg. **x5.3 and x4.1** | [M] |
| is the effect recovery or chemistry? | d-MFT, d-FFT, d-MF-S2 internal standards; text: their addition "immediately after the boiling period would compensate alterations" (p. 189) [V] → Method I's extraction and storage losses are compensated; **in-pot loss during the 2 h boil is not**. The SDE value is the one that avoids the pot | [V] / [Z] |
| MFT loss law at ~100 °C (EX-050, evidence map 2.4: "lifetime at 95-145 C: NONE") | **an inferred first-order in-pot loss λ ≈ 0.044 /min (thiamine pot) and 0.034 /min (thiamine + Cys)** at 100 °C, pH 5.7, 0.2 M pyrophosphate, under constant formation over 2 h [Z, assumption-laden, sec. 3] | [Z] |
| is the sink the disulfide? | **No.** MF-S2 < 0.2 µg (thiamine, Method I), n.d. with cysteine, 1.6 ± 1.2 µg in the full system (Table II): ≤ 1-5 % of the MFT missing from Method I [M/Z] | [M] |
| MFT sources at meat-like levels | thiamine >> ribose + cysteine; cysteine raises thiamine-MFT x7 (Method I); thiamine pyrophosphate inactive; IMP + Cys only traces | [M] |
| cysteine dose response | Table III: thiamine 30 µmol + ribose 670 µmol, Cys 0 → 2720 µmol: MFT 5.2, 14, 18, 18, 19, 14, 10, 4.8, 2.0 µg — **maximum at 170-325 µmol, falls x9.5 by 2720 µmol** (Method I) | [M] |
| buffer effect | Table IV: pyrophosphate 0.2 → 0.1 M lowers MFT to 35-46 % (Method I) | [M] |
| H2S effect | Table VI: thiamine + H2S (26.8 mmol total, 3 flushes) 15 ± 7 µg vs 4.7 without (x3); ribose, R5P, ribulose, IMP + H2S ≤ 0.4 µg | [M] |
| FFT, MF-S2 | low everywhere (FFT ≤ 2.3 µg, MF-S2 ≤ 1.6 µg) | [M] |

**Answers.**

1. **The 4-5x is a sink-avoidance effect, not a recovery bias** (my reading, **~0.75**). With SIDA standards added
   after boiling, Method I's number is the MFT present in the pot at 2 h; Method II strips MFT into the solvent as it
   forms. If formation F is the same in both and the pot loses MFT at first-order λ, SDE/reflux = λT / (1 - e^-λT)
   [Z]. The residual 0.25: (i) the standard may not equilibrate with reversibly bound MFT (mixed disulfides, adducts)
   at the end of boiling, so part of the "loss" could be release-able binding — still an in-pot state, not an
   extraction loss; (ii) the timing of the standard in Method II is not stated separately ("immediately after the
   boiling period" read for both) [!].
2. **What this means for bolton1993's caveat** (bolton1993_extraction.md, sec. "Arithmetic" point 7: "the absolute
   level, 11.7 ppb, may be a fourfold to fivefold under-estimate"): on this reading, direct extraction after heating
   measures the pot correctly; SDE measures a different quantity (cumulative formation with the sink short-circuited).
   **The bolton1993 caveat is not supported** by the source it cites [Z]. See sec. 6.
3. **λ(100 °C) in a dilute, low-cysteine pot is of the same order as the engine's `k_mft_decay`** [Z]: Grosch gives
   0.034-0.044 /min (range 0.026-0.053 propagating ±1 SD); the engine's 0.095 /min at 110 °C (EX-050 A0) transports
   to 0.031 / 0.045 / 0.068 /min at 100 °C for Ea 133 / 90 / 40 kJ/mol. So this pot does NOT show the engine's sink
   3x too fast. It contrasts with the Nestec patent bound (λ(110 °C) ≲ 0.02-0.05 /min, 0.25 M cysteine, EX-050 Ta1):
   the two are reconcilable only if the loss is pot-dependent (cysteine level, carbonyl load), i.e. it bears on
   A1 (pot-independent) vs A2 / A3 (sec. 5).

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Meatlike Flavor Compounds" |
| Authors | Werner Grosch, Gabriele Zeiler-Hilgart (Deutsche Forschungsanstalt für Lebensmittelchemie, Garching) |
| Venue | Teranishi, Takeoka & Güntert (eds.), *Flavor Precursors: Thermal and Enzymatic Conversions*, ACS Symp. Ser. **490**, ch. 14, pp. 183-192, 1992; received December 18, 1991 |
| DOI | **10.1021/bk-1992-0490.ch014** (printed in the download banner; Crossref 2026-10-31: same title, Grosch / Zeiler-Hilgart, ACS Symposium Series, 183-192, 1992). Matches the brief |
| Related | Zeiler-Hilgart's thesis (ref 26 of Cerny 2007; k3_final_parameter_inventory.md 6c) presumably holds the full data; Gasser & Grosch 1988 / 1990 (refs 1-2, MFT and MF-S2 in boiled beef); Sen & Grosch 1991 (ref 21, the SIDA method) — none on disk |

## 2. Conditions

| item | as printed |
|---|---|
| medium | Na-pyrophosphate buffer **pH 5.7, 0.2 mol/L, 500 mL** unless stated (Table IV: 0.1 M) [M]. pH chosen because meat is pH 5.5-6.0 [Q]. pH not re-measured after boiling [!] |
| reactants (Table II footnote) | thiamine·HCl 10.2 mg (**30 µmol → 60 µM**), thiamine pyrophosphate 13.9 mg (30 µmol), D-ribose 100 mg (**670 µmol → 1.34 mM**), cysteine·HCl 57 mg (**325 µmol → 0.65 mM**), GSH 71 mg (230 µmol → 0.46 mM) [M; concentrations Z] |
| rationale | concentrations "approximated to those occurring in raw beef and pork" (Table I: beef ribose 10-1300 mg/kg, total SH 690 mg/kg (2530 as Cys), non-protein SH 11 (42), GSH 71, thiamine 0.9-2.3; pork thiamine 8-10) [M/C] |
| Method I | refluxed 2 h (~100 °C, open to air through the condenser [Z]), cooled, then extracted [M] |
| Method II | distilled and continuously extracted (Likens-Nickerson SDE) for 2 h [M] |
| quantification | mass chromatography (SIDA) with d-MF-SH, d-FF-SH, d-MF-S2 (refs 20-21) [M]; added after boiling (p. 189) |
| replication | N = 2-4 trials (Table II), "at least two" (Table VIII); mean ± SD |
| Table V no. 2 | Farmer & Mottram's conditions: ribose 30 + cysteine 41 µmol/mL (**30 / 41 mM**), Na-phosphate pH 5.6, **0.5 M, 2 mL, 140 °C, 1 h** [M] |
| Table VI | before boiling and at 30 and 90 min the solution was flushed with H2S (200 mL each); total H2S 26.8 mmol, sized to all beef cysteine (2530 mg/kg) degraded [M] |
| Table VII | 500 g meat, minced, in 500 mL tap water, SDE (Method II), 2 h [M] |
| vessel / O2 | not stated beyond reflux / SDE apparatus; under air [Z] |

## 3. Quantitative results

**Table II (Method I, µg in 500 mL, mean ± SD)** [M]:

| no. | system | N | MF-SH | FF-SH | MF-S2 |
|---|---|---|---|---|---|
| 1 | thiamine | 3 | 4.7 ± 1.2 | 1.9 ± 0.4 | < 0.2 |
| 2 | thiamine pyrophosphate | 2 | < 0.5 | < 0.4 | < 0.2 |
| 3 | thiamine + ribose | 2 | 5.2 ± 0.4 | < 1.5 | n.a. |
| 4 | thiamine + cysteine | 4 | 33 ± 3.7 | < 1.5 | n.d. |
| 5 | thiamine pyrophosphate + cysteine | 2 | < 0.5 | < 0.5 | < 0.2 |
| 6 | thiamine + cysteine + ribose | 3 | 19 ± 2.5 | 1.2 ± 1.1 | < 0.2 |
| 7 | thiamine + GSH + ribose | 3 | 23 ± 1.2 | < 1.7 | < 0.2 |
| 8 | thiamine + GSH + cysteine + ribose | 4 | 9.0 ± 2.8 | 2.3 ± 1.0 | 1.6 ± 1.2 |

Text: no. 1 = 41 nmol MF-SH and 17 nmol FF-SH [M] (checked: 4.7 µg / 114.16 = 41.2 nmol [Z]). Molar yield of MFT
from thiamine: 0.14 % (Method I), 0.73 % (Method II); thiamine + Cys 0.97 % / 3.97 % [Z].

**Table III (Method I; thiamine with ribose 670 µmol, pyrophosphate 0.2 M)** [M]:

| cysteine (µmol) | 0 | 85 | 170 | 255 | 325 | 340 | 680 | 1360 | 2720 |
|---|---|---|---|---|---|---|---|---|---|
| MF-SH (µg), thiamine 30 µmol | 5.2 | 14 | 18 | 18 | 19 | 14 | 10 | 4.8 | 2.0 |
| MF-SH (µg), thiamine 6 µmol | n.a. | n.a. | 2.1 | n.a. | n.a. | 2.8 | 1.6 | 0.5 | n.a. |

**Table IV (Method I; ribose 670 µmol, thiamine 30 µmol)** [M]: Cys 325 / 340 / 680 µmol: 19 / 14 / 10 µg at 0.2 M
pyrophosphate; 8.8 / 7.2 / 3.5 µg at 0.1 M.

**Table V (Method I)** [M]: no. 1 ribose 1.34 + Cys 0.64 µmol/mL, pyrophosphate pH 5.7, 500 mL, reflux 2 h: MF-SH
< 0.1, FF-SH n.d., MF-S2 n.d. (n = 3). No. 2 ribose 30 + Cys 41 µmol/mL, phosphate pH 5.6 0.5 M, 2 mL, 140 °C 1 h:
MF-SH 0.6, FF-SH 0.4, MF-S2 n.d. µg (n = 2) → **MFT 2.6 µM, 8.8 × 10⁻⁵ per mol ribose** [Z].

**Table VI (Method I, + H2S)** [M]: thiamine 15 ± 7 (n = 4); ribose < 0.2; ribose-5-phosphate 0.4; ribulose < 0.2;
5'-inosinic acid 0.2 µg MF-SH. Text: IMP 1.3 g + Cys 2.5 g in 500 mL pyrophosphate: "only traces" (data not shown).

**Table VII (Method II, 500 g meat in 500 mL, 2 h)** [M]: pork shoulder 1.0 / 1.6 / 1.9; beef top round < 0.2 /
0.8; eye of round 2.1; rib wing 2.1 / 2.1 µg MF-SH. FF-SH and MF-S2 below 1.0 and 0.2 µg.

**Table VIII (the method comparison, µg, mean ± SD of ≥ 2)** [M]:

| no. | system | Method I (reflux, then extract) | Method II (SDE) | ratio [Z] |
|---|---|---|---|---|
| 1 | thiamine 30 µmol | 4.7 ± 1.2 | 25 ± 5 | 5.3 |
| 2 | thiamine 30 + cysteine 325 µmol | 33 ± 3.7 | 136 ± 0 | 4.1 |
| 3 | ribose 670 + cysteine 325 µmol | < 0.1 | 0.1 | — |

Method I values of rows 1-2 are the Table II rows 1 and 4 again (same means and SDs), so the Method I side is not
independent of Table II [Z]. Text: MF-SH "increased approximately 5- and 4-fold" with SDE [V]; quantification "was successful only, when
the volatiles formed were continuously removed" [V].

**Inferred in-pot loss (my derivation, not printed)** [Z]. Assume (a) MFT forms at a constant rate F over T = 120 min
in both methods (thiamine at pH 5.7 decays slowly at 100 °C; not measured here); (b) SDE removes MFT as it forms and
recovers it fully; (c) the pot loses MFT first order at λ. Then SDE/reflux = λT / (1 - e^(-λT)):

| system | ratio (±1 SD of SDE) | λ (/min) | t½ (min) |
|---|---|---|---|
| thiamine | 5.32 (4.26-6.38) | **0.044** (0.035-0.053) | 15.7 |
| thiamine + cysteine | 4.12 (3.30-4.95) | **0.034** (0.026-0.041) | 20.5 |

Direction of each assumption's error: incomplete or slow SDE stripping, SDE recovery < 100 % (standard added to the
extract, not the pot), and SDE's removal of H2S (which raises thiamine-MFT x3, Table VI) all make the true ratio
LARGER, so these λ are closer to lower bounds [Z]. A formation rate that decelerates within 2 h would make them upper
bounds. Net: λ(100 °C, this pot) ~0.03-0.06 /min, confidence ~0.6 that the true value lies in that range.

## 4. Mechanism evidence and its kind

- **Kind 6 (within-study factor: work-up method, same pot) + SIDA**, DFA Garching (Grosch group; same institute as
  Schieberle / Hofmann / Guth). Not labelling, not trapping.
- Disulfide is excluded as the main sink by measurement (MF-S2 ≤ 0.2 µg where 20-100 µg MFT is missing) [M/Z].
- The sink itself is not identified. Candidates the paper does not test: covalent coupling with carbonyls / thiamine
  fragments, mixed disulfides with cysteine (only in pots with Cys; the thiamine-only pot loses as much), oxidation by
  air through the reflux condenser.
- Thiamine route to MFT (Fig. 1: thiamine → 5-hydroxy-3-mercaptopentan-2-one → MFT; H2S may add to the oxofuran 3
  to give 4 → MFT) is [P] (van der Linde 1979 [C] plus the authors' extension); cysteine's role "may shift the
  equilibrium of SH-group replacement" [P].
- Pyrophosphate dependence (Table IV) and the cysteine optimum (Table III) are kind-6 within-study factors [M].

## 5. Bearing on open questions

- **Evidence map 2.11, "SDE recovery ... for MFT itself the 4-5x SDE effect is cited only":** now read at source.
  It is a **CHEMISTRY effect of the operator (a sink short-circuited), not a recovery factor.** CONTRADICTS treating
  it as a multiplicative recovery correction for direct-extraction data (sec. 6). SUPPORTS "operators: SDE for MFT"
  moving from NONE to SINGLE (one lab, one within-study comparison, two pots).
- **Evidence map 2.4 / 2.12, "MFT lifetime at 95-145 C: NONE":** gives the first processing-temperature,
  within-study inference: λ(100 °C) ~0.03-0.05 /min in dilute (µM-mM) thiamine pots, pH 5.7, 0.2 M pyrophosphate,
  low cysteine. SINGLE (DFA), and inferred under stated assumptions [Z]. Moves the verdict towards SINGLE, not CAN.
- **EX-050 (MFT loss law, INCONCLUSIVE):**
  - A0 (engine, 0.095 /min at 110 °C): **SUPPORTED for this kind of pot** if Ea ≳ 90 kJ/mol (engine at 100 °C 0.031-0.045
    for Ea 133-90).
    (Caveat added 2026-11-01, cluster 4 and its red team: these are THIAMINE pots. Thiamine cleavage releases the
    pyrimidinylmethyl electrophile that consumes MFT (jhoo2002), so this lambda may be a thiamine-specific capture
    sink, and its agreement with the engine's pentose-pot k_mft_decay (0.041 /min at 100 °C) may be coincidental. The
    red team also notes that the p. 189 sentence on the standards concerns the extract's storage oxidation, that
    Method II's standard goes into the extract (SDE recovery is not compensated), and that lambda assumes equal
    formation under SDE and reflux although SDE strips H2S (Table VI, x3) and the precursor: lambda is model-dependent
    and probably a lower bound. "Standards added after boiling" itself stands, P 0.85.)
  - A1 (pot-independent, slower): CONTRADICTED as a pot-independent constant if the patent's λ(110 °C) ≲ 0.02-0.05
    holds in 0.25 M cysteine; Grosch's 100 °C value already sits at that bound.
  - A3 (oxidative / coupling loss suppressed by free Cys): weakly SUPPORTED in direction (λ 0.034 with 0.65 mM Cys
    vs 0.044 without), but the difference is within 1 SD [Z]. Table III's fall of MFT at high Cys (19 → 2.0 µg,
    0.65 → 5.4 mM Cys) runs AGAINST a pure "Cys protects MFT" reading unless Cys inhibits formation from thiamine;
    the design cannot separate formation from loss [!].
  - A2 (second order in a growing electrophile pool): LEAVES UNCHANGED (no carbonyl series).
  - The red team's P(engine too fast by ≥ 3x as a pot-independent constant up to 145 °C) ~0.45: this paper is evidence
    AGAINST "too fast" in low-Cys dilute pots, so it lowers that probability for pot-independent use; it does not
    touch the cysteine-rich case.
- **CF-013** (sinks too strong; "US 10,076,126 ... the engine's MFT sink is at least 5x too fast"): CONTRADICTS the
  generalisation of that test result beyond cysteine-rich pots [Z]. Does not touch the rising / falling shape.
- **CF-068** (Schieberle 2000 vs Yiltirak 2026, ribose + Cys 100 °C): Table V no. 1 (ribose 1.34 + Cys 0.65 mM,
  pyrophosphate pH 5.7, reflux 2 h: MFT < 0.1 µg / 500 mL, < 1.8 nM) is a third, dilute data point; LEAVES the
  conflict UNCHANGED (concentrations 20-75x lower).
  (Added 2026-11-01, cluster 4: Grosch is the DFA institute, as Schieberle 2000 is, so neither Table V point is a
  third LABORATORY; Table V no. 2, ribose 30 + Cys 41 mM, 140 °C, 1 h, 0.6 µg / 2 mL = 300 µg/L, is a DFA-institute
  point.)
- **Strecker / pyrazine precursors:** not touched.
- **Thiamine as an MFT source** (k3 inventory 6c, Zeiler thesis): this chapter is the printed source of the
  thiamine >> ribose + Cys ranking at meat-like levels; the "45 vs 3.4 µg/L" pair is NOT in this chapter [!] (no
  120 mmol/L thiamine pot here).

## 6. Disagreements found

1. **bolton1993_extraction.md (point 7) and the evidence map 2.11 treat the 4-5x as a recovery bias of direct
   extraction** ("may be a fourfold to fivefold under-estimate"). Source: Method I quantifies with deuterated
   standards added after boiling, so its extraction is compensated; the 4-5x is MFT that did not survive the pot.
   Candidate causes of the misreading: Bolton's citation wording ("recovered"); the SIDA timing is on p. 189, away
   from Table VIII. Consequence: no 4-5x correction should be applied to direct-extraction levels of heated pots; an
   SDE-measured pot should instead be modelled as cumulative formation with the sink short-circuited.
2. **Grosch λ(100 °C) ≈ 0.03-0.05 /min (dilute thiamine pot, ≤ 0.65 mM Cys) vs Nestec patent λ(110 °C) ≲ 0.02-0.05
   /min (0.25 M cysteine + xylose)** — same order at face value, but ×1.5-3 apart once 10 °C is allowed for. Causes:
   cysteine level (A3), carbonyl / sugar load (A2), the patent's single unreplicated series, or assumption (a) above.
   (Added 2026-11-01, cluster 4: also MFT lost through the open reflux condenser, an apparatus loss, small, P ~0.15;
   and a thiamine-specific electrophile sink, jhoo2002. The red team of cluster 4 adds that the two numbers measure
   different channels -- an inferred in-pot loss in a thiamine pot vs a bound derived from 4 single-vial digitised
   points with co-heated standards, not measurement-grade, P 0.15 -- so P(pot-dependent >= 3x at matched T) is ~0.5.)
3. **Same-lab tension with Guth 1995** (MFT-MFT 86 % gone in 2 h at 100 °C, 62 % back as MFT, EX-050 Ta6): Grosch
   finds almost no MF-S2 in these pots, consistent with Guth's disulfide reverting at 100 °C; not a disagreement, a
   consistency.

## 7. Ledger-ready measurements

- **Table VIII rows 1-2** (two readers): the method pair, same pot, SIDA — the one within-study quantity for the
  operator model (SDE vs reflux) and the λ inference. Key both methods with the SIDA timing quoted.
- **Table II rows 1, 4, 6, 7** (Method I): thiamine-MFT levels at meat-like concentrations, 2 h, 100 °C, pH 5.7 —
  data for a thiamine route if one is ever added (the engine has none).
- **Table III** (cysteine dose series, 9 levels): a within-study factor for the thiamine route x cysteine; confounds
  formation and loss.
- **Table V no. 2**: ribose 30 + Cys 41 mM, 0.5 M phosphate pH 5.6, 140 °C, 1 h → MFT 0.6 µg / 2 mL (2.6 µM), Method
  I, n = 2 — a ribose + cysteine pot with SIDA at 140 °C (check against the engine's sulfur lane).
- Table VII (meat, SDE): context only (protein matrix; see mottram1996b for the SDE redox artefact in meat).
