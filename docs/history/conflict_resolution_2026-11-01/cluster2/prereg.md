# Cluster 2 pre-registration: the deoxyosone hub conflicts (CF-060, 071, 053, 063) -- 2026-10-31/11-01

Written BEFORE any script in this folder was run. Read-only on tracked files; scratch only.

**Disclosure of what was seen before writing.** I read the dossiers (zhou2022, zhou2023c, zhou2024, zhou2024b,
zhou2025b, cui2022, hou2026, wei2022, beck1988, cao2017, zhang2025b, martins2005 / 2005b, kang2026) and the four
register entries, so the printed tables are known to me. Three numbers were done mentally while reading and are
therefore NOT blind: (a) Zhou 2022's +-Cys double ratio falls from ~1.3 at 30 min to ~0.4 at 90-120 min; (b) the
OPD derivatisation time does not sort the AX-ARP pots into high / low scale (2022 4 h high, 2024b 4 h low, 2023c
by ref = 2022 low, 2024 12 h high, 2025b 12 h low); (c) a rough H-pH prediction for Hou's cystine double ratio
(10^(-0.49 dpH)) rises toward 1 late while the printed one stays ~0.33. Tests built on (a)-(c) are reported as
non-blind.

## 0. Hypotheses (stated neutrally) and rivals

| code | hypothesis | predicts |
|---|---|---|
| H-op | Jiangnan / Xie data are STANDING levels (ex-situ OPD after the quench); Beck 1988 (in-situ OPD) is a cumulative formation split. L3/L1 = (F3/F1)(tau3/tau1) when quasi-steady | flux 1-DX : 3-DX = 20 : 1 and level 3-DX / 1-DX = 2.6 - 22 are reconcilable with tau1/tau3 << 1 |
| H-pH | unbuffered "pH 7 / 8" pots are at pH ~5 -> 3.5 within 20 - 30 min (Hou 2026, Zhou 2024 Fig. 3b, Kang 2026 endpoints); label pH is not reaction pH | label-pH contrasts shrink as in-situ pH converges; +-additive contrasts carry a pH term when the additive changes the drift |
| H-scale | 3-DX / 1-DX are semi-quantitative (no authentic standard; CF-062) | only within-figure, same-analyte contrasts compare across papers; absolute 3-DX levels may differ by a response factor per paper |
| H-trap | a cysteine thiol traps 3-DX | the +-Cys double ratio < 1 for thiols, ~1 for non-thiol additives |
| R-amine (rival to H-trap) | an amine (free alpha-NH2, thiol or not) or the ARP itself removes 3-DX faster than 1-DX (Zhou 2024b ARP + DX condensation; Hou cystine) | double ratio < 1 for any reactive amine; graded by amine reactivity, not by SH |
| R-pHadd (rival to H-trap) | the additive changes the in-situ pH drift, and the double ratio is the split / lifetime pH law acting through dpH | double ratio = 10^(-s dpH) with s from signature (ii); -> 1 when dpH -> 0 |
| R-batch (rival for CF-060) | ARP batch / co-purified catalyst or sample handling differs between papers | calibrated species (ARP, Ala, GO, MGO, furfural) differ along with 3-DX; no pattern by method |
| R-print (rival for CF-060) | a figure or unit slip in one print | that print violates a physical bound (carbon or nitrogen closure) |

## 1. Tests, decision rules (declared now), and expectations

### CF-060 -- the AX-ARP dicarbonyl scale (Zhou 2022 / 2023c / 2024 / 2024b / 2025b)

**T1 carbon closure.** Per pot and time: R_C = (5 [3-DX] + 5 [1-DX] + 2 [GO] + 3 [MGO] + 5 [FUR]) / (5 dARP), all in
mM; dARP from the pot's calibrated ARP curve (Zhou 2024 / 2025b Fig. 3a read from the PDF if not in a dossier).
Also N closure: released Ala / dARP <= 1. **Rule:** a print with R_C > 1.25 at any time (combined read and
calibration error declared +-25 %; 3-DX carries an unknown response factor, so a violation condemns the ABSOLUTE
3-DX print, not the pot) is declared physically wrong on its absolute scale and scored shape-only (it already is);
R_C <= 1 everywhere = closure does not discriminate. **Expect:** no print violates (P 0.6); Zhou 2024 violates at
some time (P 0.3); a low print violates (P 0.03).
**T2 method pattern.** Tabulate per paper: OPD time, sample : reagent volume, IS, ARP lot / synthesis, Ala method,
headspace. **Rule:** OP (protocol) is supported only if a protocol field sorts all five pots into the high (2022,
2024) and low (2023c, 2024b, 2025b) groups. (Non-blind for OPD time: known not to sort.) **Expect:** no field sorts
(P 0.7).
**T3 species-ratio pattern.** high / low ratio per species (3-DX, 1-DX, GO, MGO, Ala, furfural) for the two matched
pairs (2022 vs 2024b at pH 7; 2024 vs 2025b at pH 8). **Rule:** R-print (a dilution / unit slip) predicts a common
factor within 0.15 dex across species; H-scale predicts the factor confined to 3-DX / 1-DX; R-batch predicts all
species move, not necessarily together. **Expect:** neither a common factor nor DX-only (P 0.8).
**T4 H-pH.** Both pH 8 pots carry measured trajectories (2024 Fig. 3b in the fit table; 2025b Fig. 3d "~4.1-4.5" at
120 min). **Rule:** H-pH is excluded for a pair whose in-situ pH trajectories agree within 0.4 pH units at every
time while the dicarbonyl ratio exceeds 0.5 dex. **Expect:** H-pH excluded for 2024 vs 2025b (P 0.85).
**Verdict rule.** RESOLVED if one print is condemned by T1 and the others are mutually consistent; EXPLAINED (P) if
T2 / T3 single out one cause class with P >= 0.6; else OPEN with the deciding measurement named.
**Overall expectation:** OPEN (P 0.6), EXPLAINED as R-batch / handling (P 0.25), RESOLVED (P 0.15).

### CF-071 -- the cysteine double ratio (Zhai 2019 vs Zhai 2023 JAFC)

**T5 matched windows.** Double ratio DR(t) = (TDP / DPO)_{+Cys} / (TDP / DPO)_{-Cys} per time, from the keyed rows;
Zhai 2019 (120 C, 10-60 min) vs JAFC at 120 C <= 60 min, then all JAFC temperatures and times; the copied JAFC
T-140 120 min 3-DX cell dropped; points below the paper's LOQ flagged. Time trend: slope of log10 DR vs t in JAFC.
**Rule (triage's):** EXPLAINED (time window) if the matched-window values agree within 0.15 dex (mean log DR) and
the >1 excursions lie only outside the window. **Expect:** matched windows agree (P 0.4); JAFC DR rises with time
(P 0.55).
**T6 the same-lab additive panel (H-trap vs R-amine vs R-pHadd).** DR for every same-lab additive pot on disk:
cysteine (Zhai 2019, JAFC, Zhou 2022 AX-ARP, Wei 2022 GX-ARP, both panel readings), glutathione (Zhou 2023c), cystine
(Hou 2026), alanine (Zhou 2024b pH 7; Zhou 2024 pH 8). **Rule:** H-trap is supported only if the thiols (Cys, GSH)
give DR <= 0.7 AND the non-thiols (cystine, Ala) give DR within [0.8, 1.25]. R-amine is supported if cystine gives
DR <= 0.7 with Ala ~1 (amine reactivity, not SH). Signature (iii) is "not thiol-diagnostic" if any non-thiol gives
DR <= 0.7. **Expect:** cystine <= 0.7 (P 0.9, it is printed); Ala ~1 (P 0.85); thiol DRs all <= 0.7 at late times
(P 0.5); hence "(iii) not thiol-diagnostic" (P 0.85).
**T7 H-pH on the cystine control (non-blind, see disclosure (c)).** Predicted DR = 10^(-s dpH(t)) with s in [0.25,
0.6] dex/pH (signature ii and Beck's flux slope, bands from data other than Hou) and dpH(t) Hou's measured +-cystine
difference. **Rule:** R-pHadd alone is rejected for cystine if the printed DR lies below the band's lower edge by
> 0.15 dex at >= 3 of 6 times. **Expect:** rejected late (P 0.7).
**Overall expectation for CF-071:** OPEN with the time trend explained in part (P 0.5); EXPLAINED (P 0.35).

### Hub question -- H-op vs Beck, and "the split law near off"

**T8 chain model (standalone module, A2 contract sec. 3 shape).** ARP -> X1 (k_f F1), ARP -> X3 (k_f F3) with ARP
forced from each pot's calibrated curve (piecewise-linear); X_i -> sinks (k1, k3, first order); in-situ pH forced
(measured where printed; else a band between Hou's and Zhou 2024's measured trajectories, mixing weight lambda in
[0, 1] per pot); log10(F1/F3)(pH) = a7 + b (pH - 7); log10 k_i(pH) = log10 k_i,7 + c_i (pH - 7).
Bands from data OTHER than the Jiangnan pots: a7 in [log10 1.45, log10 20] (Martins 2005 fitted 100 C .. Beck
in-situ); b in [0.2, 0.45] (Martins 2005b pD7 - pD4 0.25 .. Beck 0.44); k1,7 in Martins 1-DG -> AA at 120 C
(1.4 /min at 100 C, Ea 76 -> ~4.9 /min) +-1 dex; k3,7 Martins 3-DG -> FA at 120 C (~0.057 /min) +-1 dex; c1, c3 in
[0, 0.6] (Martins 2005b pD bands). Scored observable: log10(3-DX / 1-DX) per point (scale-free; tolerance 0.25 dex,
i.e. 3 x the 0.075 within-pot scatter of prereg sec. 24 plus a transfer allowance). Pots: Hou 2026 (pH measured),
Zhou 2024 (pH measured), Zhou 2022, Zhou 2024b, Zhou 2023c, Zhou 2025b, Wei 2022.
Search: bounded multi-start minimisation of the worst-point residual (sec. 25 decision 2), 20 starts per arm.
**Arms:** A (H-op + H-pH, full bands); B "split off" (b = 0, a7 free in [-1.5, 1.5]); C "no drift" (label pH held);
**negative control N1** "standing = flux" (k1 = k3, i.e. level ratio equals the formation split: must FAIL);
**negative control N2** inverted lifetimes (k1 < k3 forced: must FAIL); **positive control P** synthetic ratios
from one arm-A point inside the bands at sigma 0.075 dex on the same pots and times (must PASS).
**Feasible** = every scored point within 0.25 dex (per pot, then jointly with shared chemistry and per-pot lambda).
**Rules:** (1) H-op reconciles Beck with the Jiangnan levels if arm A is feasible jointly AND N1 fails AND P passes;
(2) "split law near off" is DISFAVOURED if B is infeasible while A is feasible; INCONCLUSIVE if both are feasible
(the levels cannot separate split from lifetime); SUPPORTED only if B is feasible and A's feasible region needs
b < 0.1 (impossible by A's band: so the most it can be is "not excluded"); (3) H-pH is needed if C is infeasible
while A is feasible. Inconclusive is allowed when P fails.
**Expect:** P passes (P 0.85); N1 fails (P 0.9); N2 fails (P 0.8); A feasible per pot (P 0.65), jointly (P 0.4);
B feasible (P 0.5) -> split question inconclusive (P 0.45); C infeasible (P 0.55).

**T9 the label-pH contrast under drift (Zhai 2021 TTCA ladder, 100 C).** H-pH (instantaneous) predicts that the
between-label spread of log10(TDP / DPO) shrinks with time as in-situ pH converges (Kang 2026: 5.5 / 7 / 8 -> 5.1 /
4.9 / 4.7 at 120 C, 120 min). **Rule:** if the pH 5.5 vs 8 spread at 100-120 min is >= 0.8 x that at 20-40 min, an
instantaneous-pH reading of signature (ii) is disfavoured at 100 C (either pH did not converge at 100 C or the
contrast is history-dependent). **Expect:** spread persists (P 0.6).

### CF-053 -- Cao 2017's glycine multiplier

**T10.** (i) Physical: glycine solubility at 90 C (literature 62-67 g / 100 g water at 90-100 C, my values, not on
disk) vs the nominal arms (x8 24, x12 36, x16 48, x24 72 g / 100 g of buffer if "mol/L" = mol per L of buffer
charged). (ii) Kinetic back-solve: excess TTCA first-order drain over the blank (0.5-2 h, Fig. 5a) and peak Gly-ARP
(Fig. 5c) vs the nominal multiplier: proportional (slope ~1 in log-log) = the glycine is in solution at roughly the
stated ratio; flat across x8-x24 = a saturated slurry. **Rule:** the x8-x16 arms re-enter as ordinal (qualitative,
monotone-in-glycine) rows if (i) says dissolved at 90 C AND (ii) is not flat; x24 stays out if (i) puts it above
solubility. EXPLAINED (P) only to that extent; absolute [Gly] remains a declared band. **Expect:** x8-x16 dissolve at
90 C (P 0.7); response not flat (P 0.75); verdict EXPLAINED-partial (P 0.55), OPEN (P 0.45).

### CF-063 -- Zhang 2025b printed fits vs points

**T11.** Refit ln(C/C0) of the blank xylose and cysteine points (route table, 5-50 min) with free slope + intercept,
and with free slope + lag (intercept forced 0 at t = lag). **Rule:** EXPLAINED (fit window / unreported intercept) if
the free-intercept fit over SOME contiguous window of >= 4 points reproduces the printed (k, b) within 15 % in k and
0.1 in b; else the printed constants are a transcription / fitting slip (error_confirmed on the paper side). Mass
balance (Cys + TTCA + CAMA vs 40 mM; Xyl likewise) reported per time; the missing half passes to CF-076 as the sink
term. **Expect:** a window reproduces the printed fits (P 0.55); missing ~45-55 % at 50 min (P 0.8).

## 2. What would change for the hub (declared now)

- If T8 rule (1) holds, Beck's flux split and the Jiangnan standing levels stop being in conflict; the hub's
  dpo_share should be anchored on a FLUX datum, never on a standing ratio (CR-12's basis).
- If T6 says "(iii) not thiol-diagnostic", no cysteine-trap credit survives from signature (iii) or (vii) alone;
  only a thiol-vs-non-thiol contrast at matched in-situ pH could revive it.
- If T1 condemns no print, CF-060 stays OPEN and the absolute Jiangnan DX scale stays unusable (status quo).
</content>
</invoke>

## Amendment A1 (written after reading Zhou 2024 Fig. 3a-b / Fig. 4 and Zhou 2025b Fig. 3a, b, d from the PDFs
## for T1's dARP input, BEFORE any script ran)

Data acquisition for T1 exposed numbers that bear on T3 / T4 (non-blind from here on for CF-060):
- Zhou 2024 (pH 8): ARP 5.08 / 2.84 / 1.65 / 1.25 mM, pH 5.13 / 4.32 / 4.19 / 4.08, released 14N-Ala 5.3 / 6.8 / 7.6
  / 8.4 mM (Fig. 4 green) at 30-120 min. Zhou 2025b (pH 8): ARP 3.57 / 1.79 / 1.25 / 0.95, pH 6.04 / 5.05 / 4.66 /
  4.50, Ala 11.8 / 14.4 / 14.1 / 12.6. My 250-300 dpi reads, +-0.15 mM, +-0.03 pH, +-0.3 mM Ala.
- So CF-060's "pH agree" is loose: 2025b is 0.9 / 0.7 / 0.5 / 0.4 pH units LESS acidic. T4's rule (agree within 0.4)
  fails, so H-pH is NOT excluded by T4 as written. Added magnitude rule **T4b**: H-pH can explain the pair only if
  the required absolute-level slope (dex of dicarbonyl level per pH unit of in-situ difference) lies within the
  on-disk rate slopes (|pD| <= 0.6 dex/pH, Martins 2005b; 0.44 Beck) AND has the sign of the base-catalysed
  formation (the less acidic pot higher). Expect H-pH fails T4b (P 0.85).
- New **T1b nitrogen closure** (every enolisation and the Amadori reversal release the amine 1 : 1): free Ala >=
  standing (3-DX + 1-DX) unless an Ala sink exists. Rule: a pot whose standing DX exceeds free Ala by more than the
  combined read error at >= 2 times needs an Ala sink of that size; if the lab's other pots of the same design show
  free Ala ~0.7-1.0 per ARP lost (no room for such a sink), the high-DX prints and the high-Ala prints cannot share
  one chemistry. Expect: the two high prints (2022, 2024) violate early (P 0.7, non-blind: mentally computed).
- New **T3b**: across the five AX-ARP pots, Spearman rank of log(DX / dARP) vs log(Ala / dARP). A strong negative
  rank (rho <= -0.8) ties the DX scale to a CALIBRATED co-variable, which a DX response factor (H-scale) cannot
  produce: R-batch / pot chemistry then leads. Expect rho <= -0.8 (P 0.6).

## Amendment A2 (before T8 is coded or run; after T1-T7, T9 were run)

- Implementation of T8: exponential-integrator steps of 1 min; ARP forced by PCHIP of log ARP; pH forced by linear
  interpolation; unmeasured pH (Zhou 2022, 2024b, 2023c, Wei 2022) = lambda * (Hou's measured trajectory - 0.3) +
  (1 - lambda) * min(label pH, Zhou 2025b's measured trajectory) (the fastest and slowest measured acidifications on
  disk, the fast edge widened 0.3 for the 4x higher ARP charge). Wei 2022 enters with the axes / text panel reading;
  the caption reading is run as a sensitivity. Search: scipy differential_evolution (bounded, vectorised; popsize 20,
  maxiter 400, seed 20261101) on sum(max(0, |res| - 0.25)^2), then polished; 3 seeds per arm.
- **Added sanity negative control N0** (must fail on arithmetic alone): label pH held AND k1 = k3 -- the Beck flux
  ratio read directly as a level ratio at pH 7 - 8 (F3/F1 <= 0.69 against printed level ratios >= 1.2).
- **Caveat declared now:** with the measured drift to pH ~3.5 - 4.5 and b up to 0.45, F3/F1 can reach ~10 - 26 at
  late times even with equal lifetimes, so N1 (k1 = k3, drift allowed) MAY pass; if it does, the standing levels
  cannot by themselves establish the lifetime asymmetry (H-op) -- the drift alone would reconcile Beck with the
  levels -- and rule (1) is reported as "not established by the levels" rather than "reconciled". P(N1 passes)
  revised from 0.1 to 0.35.
