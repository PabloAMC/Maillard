# Thornalley, Langborg & Minhas 1999 — EXTRACTION
DOI: 10.1042/bj3440109 (checked against Crossref, 2026-11-01)
### *"Formation of glyoxal, methylglyoxal and 3-deoxyglucosone in the glycation of proteins by glucose"*
### Glucose 50 mM alone, + Nα-t-BOC-lysine 50 mM, + human serum albumin 50 mg/mL, and fructosyl-lysine 0.5 mM, in sodium phosphate pH 7.4, **37 °C, 21 days**: glyoxal (GO), methylglyoxal (MGO), 3-deoxyglucosone (3-DG) by DDB derivatisation + HPLC. Phosphate **10-100 mM** series (initial rates) and **± 100 µM DETAPAC (= DTPA)** at 100 mM phosphate. Fed 3-DG ± BOC-Lys (2 d). **No argon arm, no other buffer, no 1-DG, no glucosone measured, one pH, one temperature.**

**Read 2026-10-31 from `data/articles/thornalley1999.pdf`** (8 pp., Biochem. J. 344:109-116; text layer clean except
decimal points printed as "n"). Figs. 1-3 read from 220 dpi renders of pp. 111-113 (read error ~ ±0.03 µM/day on the
2-µM/day axes, ~ ±0.15 on the 10-µM/day axes, ~ ±10 µM on Fig. 3a).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Questions: hub H4 / H6 (`docs/history/mechanism_reviews_2026-10-31/hub.md`: "GO / MGO / 3-DG formation ... depends on
phosphate concentration and trace metals (37 C; H4, H6)"), triage items 7 (oxygen) and 9 (phosphate), evidence map
sec. 1.2 ("Deoxyosone split vs phosphate: NO study ... thornalley1999 ... in the queue").

| need | what Thornalley 1999 gives | tag |
|---|---|---|
| O2 / argon contrast | **None.** All pots aerobic (not stated explicitly; no deaeration described). O2 enters only through the proposed autoxidation steps (Scheme 1) | [!] / [P] |
| phosphate vs other buffer at one pH | **None.** Only phosphate 10 / 25 / 50 / 100 mM at pH 7.4 (25 mM absent in Fig. 1b). Ionic strength not held; pH at each phosphate level not reported as measured; no non-phosphate buffer | [M] / [!] |
| phosphate concentration law | supralinear, and in the BOC-Lys and fructosyl-lysine arms nearly flat from 10 to 50 mM then a jump to 100 mM (sec. 3). Glucose alone, 10 -> 100 mM: GO ×14, 3-DG ×3.8, MGO ×2.9 [F] | [F] / [Z] |
| metal / chelator | **DETAPAC 100 µM at 100 mM phosphate cuts every initial rate:** GO −87 / −97 / −90 %, 3-DG −41 / −87 / −63 / −38 % (glucose / +BOC-Lys / +HSA / fructosyl-lysine) | [F] |
| is phosphate confounded with metals? | **Yes, by design and by the authors' own reading:** "trace redox active metal ions in phosphate buffer" promote the oxidative steps [Q]; the DETAPAC test is run at 100 mM only, never crossed with the phosphate series | [Q] / [!] |
| rate constants / Ea | initial rates (µM/day, n = 3) printed; no k, no Ea. My first-order equivalents in sec. 3 | [M] / [Z] |
| isotope labels | none | — |
| transfer 37 °C -> 100-140 °C | **Not shown.** 37 °C only | [!] |

**Answer.** The paper shows that at 37 °C, pH 7.4, the formation of GO, MGO and 3-DG from glucose, a glucose + amine
pot and an Amadori compound **rises with phosphate concentration and falls with a chelator**, and that **the
amine's acceleration of GO and 3-DG formation from glucose is wholly chelator-sensitive**: with DETAPAC the BOC-Lys
pot's GO and 3-DG rates are "not significantly different" from glucose alone [M/Q]. It cannot separate phosphate
catalysis (HPO4²⁻ general base, the authors' proposal [P]) from the trace metals the phosphate salt carries: the
×20 jump between 50 and 100 mM in the BOC-Lys arm (an apparent order ~4 [Z]) is not what general-base catalysis gives
and is what a metal impurity scaling with salt would give (my reading, ~0.6) [Z]. **It does not measure 1-DG, so it
says nothing about the 1-DX / 3-DX split** (H4 stays a declared band; confidence ~0.95 that no split information is
here). Transfer to 100-140 °C: none shown.

## 2. Conditions

| item | as printed |
|---|---|
| pots | (a) glucose 50 mM; (b) glucose 50 mM + Nα-t-BOC-lysine 50 mM; (c) glucose 50 mM + HSA 50 mg/mL; (d) fructosyl-lysine 0.5 mM (synthesised, free α-amine after BOC removal) |
| buffer | 100 mM sodium phosphate pH 7.4 (main runs); 10-100 mM for the initial-rate series; "aseptic conditions"; Chelex treatment NOT stated [!] |
| chelator | DETAPAC 100 µM, at 100 mM phosphate only |
| temperature / time | 37 °C, 21 d (rates over the first 2 d) |
| atmosphere | not stated; vessel and headspace not stated [!] |
| analysis | 1 mL sample + ice-cold HClO4 0.6 M; DDB (1,2-diamino-4,5-dimethoxybenzene) 4 h, RT, dark; SPE; HPLC fluorescence; internal standard DEMQ. Recoveries **52 / 70 / 59 %** (3-DG / GO / MGO); LOD 21 / 11 / 14 pmol; inter-batch CV 15 / 11 / 10 %. Whether results are recovery-corrected: not stated [!]. Derivatisation ex situ (at sampling) |
| artefact checks | no interference by 50 mM glucose or 1 mM fructosamine, except 20 nM 3-DG in glucose at t0 (attributed to contamination) |
| fructosamine | NBT reduction, fructosyl-lysine standard |

## 3. Quantitative results

**Printed [M]** (mean ± SD, n = 3):

| pot (100 mM phosphate) | GO | MGO | 3-DG | other |
|---|---|---|---|---|
| glucose, 21 d (µM) | 59.5 ± 3.7 | 2.77 ± 1.22 | 7.46 ± 1.48 | — |
| glucose, initial rate (µM/day) | 1.72 ± 0.05 | 0.098 ± 0.003 | 0.437 ± 0.032 | — |
| + BOC-Lys, 21 d (µM) | 43.7 ± 4.0 | 10.7 ± 2.0 | 14.8 ± 0.5 | fructosamine 1424 ± 8 |
| + BOC-Lys, initial rate | 8.09 ± 0.78 | 3.13 ± 0.15 | 2.62 ± 0.62 | fructosamine 116 ± 4 |
| fructosyl-lysine 0.5 mM, 21 d (µM) | 34.1 ± 3.5 | 1.61 ± 0.45 | 8.80 ± 2.86 | FL fully degraded |
| fructosyl-lysine, initial rate | 1.66 ± 0.60 | 0.55 ± 0.01 | 1.03 ± 0.15 | FL loss 37.5 ± 10.8 µM/day |
| 3-DG 50 µM, 2 d | — | 0.09 ± 0.03 | 43.1 ± 2.4 | — |
| 3-DG 50 µM + BOC-Lys 50 mM, 2 d | — | 0.55 ± 0.08 | 30.3 ± 2.0 | — |

BOC-Lys raises the initial rates ×5 / ×32 / ×6 (GO / MGO / 3-DG) [M]. HSA: GO and MGO < 0.5 µM throughout [M].

**Phosphate series, initial rates (µM/day)** [F]:

| pot | species | 10 mM | 25 mM | 50 mM | 100 mM |
|---|---|---|---|---|---|
| glucose (Fig. 1b; 10 / 20 / 50 / 100) | GO | 0.12 | 0.27 (20) | 0.60 | 1.73 |
| | 3-DG | 0.13 | 0.18 (20) | 0.22 | 0.49 |
| | MGO | 0.08 | 0.09 (20) | 0.12 | 0.23 |
| + BOC-Lys (Fig. 2e) | GO | ~0.1 | ~0.1 | ~0.4 | 8.1 |
| | 3-DG | ~0.25 | ~0.1 | ~0.25 | 3.2 |
| | MGO | ~0.05 | ~0.1 | ~0.3 | 2.7 |
| + HSA (Fig. 2f) | 3-DG | 0.36 | 0.41 | 0.56 | 1.22 |
| | GO, MGO | ~0.1 | ~0.15 | ~0.1 | ~0.07 |
| fructosyl-lysine (Fig. 3c) | GO | 0.22 | 0.25 | 0.27 | 1.67 |
| | 3-DG | 0.24 | 0.31 | 0.39 | 1.03 |
| | MGO | 0.21 | 0.24 | 0.24 | 0.56 |

Apparent log-log orders 10 -> 100 mM, glucose alone [Z]: GO 1.16, 3-DG 0.58, MGO 0.46. BOC-Lys, 50 -> 100 mM: GO ~4.3,
3-DG ~3.7 [Z] (low-end read error ±0.15 makes the 10-50 mM ratios meaningless).

**DETAPAC 100 µM (100 mM phosphate), initial rates − / + (µM/day)** [F]: glucose GO 1.73 / 0.22, MGO 0.10 / 0.13
(n.s.), 3-DG 0.49 / 0.29; + BOC-Lys GO 8.1 / 0.25, MGO 3.2 / 0.3, 3-DG 2.7 / 0.35; + HSA 3-DG 1.05 / 0.39 (GO, MGO
~0.06 / 0.1); fructosyl-lysine GO 1.67 / 0.16, MGO 0.56 / 0.08, 3-DG 1.21 / 0.75.

**Fructosyl-lysine decay (Fig. 3a)** [F]: 500 / ~425 / ~188 / ~38 / ~0 µM at 0 / 2 / 5 / 10 / 21 d. Not first order:
ln-slope 0.08 d⁻¹ over 0-2 d, ~0.26 d⁻¹ over 5-10 d [Z] (an accelerating loss).

**Derived [Z]:** glucose-alone 3-DG formation ≈ 0.437 µM/day / 50 mM = 8.7e-6 d⁻¹ (6e-9 min⁻¹) at 37 °C; fed 3-DG
loss 0.074 d⁻¹ alone, 0.25 d⁻¹ with 50 mM BOC-Lys (t½ 9.4 and 2.8 d); MGO from fed 3-DG 1.3 % (alone) and 2.8 % (+ BOC-Lys)
of 3-DG lost in 2 d; net 3-DG yield from fructosyl-lysine 1.8 mol % at 21 d, GO 6.8 mol %.

## 4. Mechanism evidence and its kind (Thornalley laboratory, Essex)

- **Within-study factors (kind 6):** phosphate concentration (one buffer, one pH) and chelator (one phosphate level).
  The two are not crossed, so "phosphate catalysis" and "metals carried by phosphate" are not separated [!].
- **Fed intermediates (kind 2):** fructosyl-lysine -> GO / MGO / 3-DG (rates parallel FL loss); 3-DG -> MGO (small,
  amine-accelerated ×6 [Z]).
- **Comparison of pots:** with DETAPAC, glucose + BOC-Lys = glucose alone for GO and 3-DG, so the amine's effect on
  those two depends on metals [M/Q]; MGO stays amine-accelerated with DETAPAC [M].
- **Proposed, not shown [P]:** HPO4²⁻ catalysing the activating C-2 deprotonation; metal-phosphate complexes
  activating 3-DG formation; glyoxal by retro-aldol plus metal-catalysed autoxidation of glycolaldehyde; MGO by
  retro-aldol of 3-DG (or its Schiff base). Glucosone is drawn but not measured.

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| Triage 9 / evidence map 1.2 "Deoxyosone split vs phosphate: NO study" | **LEAVES UNCHANGED** | no 1-DG; 3-DG rate vs phosphate only, metal-confounded, 37 °C |
| Hub H4 (1-DX sinks calibrated in phosphate) | **LEAVES UNCHANGED** | nothing on 1-DX |
| Hub H6 / triage 7 (oxygen) | **weakly SUPPORTS a metal-catalysed oxidative GO source** from glucose and from glucose + amine at 37 °C | chelator, not O2, was varied; O2 dependence is inferred |
| Evidence map "O2: Wolff vs ARP oxidation — NONE" | **partial**: Wolff-type GO from free glucose exists and is ~87 % chelator-sensitive; with amine, the extra GO is entirely chelator-sensitive | no argon arm; no glucosone; 37 °C |
| CF-011 / CF-076 (phosphate catalysis of the Cys exchange / Amadori step) | **LEAVES UNCHANGED** | no thiazolidine, no ring opening; general support that phosphate speeds early-stage chemistry at 37 °C, with a metal confound that matters for the "k_P[P]" form |
| CF-072 / CF-052 (glycation arm) | **LEAVES UNCHANGED** | fructosyl-lysine fate at 37 °C; lysine not measured |
| Amine catalysis (triage list item 1) | relevant: BOC-Lys ×5-6 on GO / 3-DG vanishes with DETAPAC; ×32 on MGO persists | a caution that "amine acceleration" can be metal-mediated in phosphate |

## 6. Disagreements found

1. **With zyzak1995 (Baynes; Nα-formyl-fructoselysine, 37 °C, pH 7.4, phosphate 10-200 mM):** Zyzak's parent decay is
   linear in phosphate (k = 0.0105 + 0.157 [Pi]/M d⁻¹ under N2 + chelators; ×2.3 from 10 to 100 mM). Thornalley's
   product formation from fructosyl-lysine is flat from 10 to 50 mM and jumps ×4-7 from 50 to 100 mM. Quantity differs
   (products vs parent loss) and Zyzak's law is under chelators, Thornalley's without; candidate causes: metal impurity
   scaling with the phosphate salt (Thornalley does not state Chelex), the Nα substituent, product vs parent. Not
   register-grade until put on one quantity; worth a note under hub H4.
2. Internal [!]: the text cites the 1424 µM fructosamine value to "Figure 2b"; it is Fig. 2c (2b is HSA α-OA). The text
   gives the fructosyl-lysine 3-DG initial rate as 1.03 ± 0.15; the Fig. 3d "− DETAPAC" bar reads ~1.21 (Fig. 3c at
   100 mM reads 1.03).

## 7. Ledger-ready measurements (two readers)

1. Printed initial rates and 21-d levels (table in sec. 3) — exact, n = 3, with the 37 °C / metal-confound flags.
2. Figs. 1b / 2e / 2f / 3c phosphate series — the only within-study phosphate concentration law for 3-DG formation
   on disk; keyed with "phosphate co-varies with trace metals".
3. DETAPAC bars (Figs. 1c / 2g / 2h / 3d) — the chelator contrast per species.
4. Fed 3-DG ± BOC-Lys (printed) — 3-DG lifetime and MGO yield at 37 °C.
5. Fig. 3a fructosyl-lysine decay — an Amadori lifetime at 37 °C in 100 mM phosphate.
