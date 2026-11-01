# Cluster 4 results (2026-11-01): observation operators and vessel O2, thiol and Strecker levels

Prereg: `prereg.md` (written before the script ran). Script: `cluster4_tests.py` (run in the maillard env, arithmetic
only, no engine call). Numbers: `results.json`. Covariate table: `covariates.md` (plus Zhou rows in sec. 7).

## 1. CF-068 (Schieberle 2000 vs Yiltirak 2026): OPEN, a declared between-lab band (P that O2 is >= half of it: 0.45)

**Gap decomposition.**
- Raw gap: 2.12 dex (Schieberle at 4 h, interpolated as t^1.35, 901 ug/L, against Yiltirak's 6.88).
- Mass action plus pH explain 0.35 dex (first order in Cys), 0.83 dex (first order in pentose) or 0.95 dex (product).
- Residual: Yiltirak sits -1.17 to -1.77 dex below the DFA anchor under every variant (as expected, P 0.8).

**Engine-free residuals of the other isotope-quantified pots** (first order in the limiting precursor; xylose x0.72; pH
0.45 dex per unit above 5; formation Ea 80-130):

| pot | residual | notes |
|---|---|---|
| US 8,920,862 (95 C, 6 h, standards after heating) | -0.55 to -0.65 dex | -1.8 under the product law; the 1.45 M sugar breaks mass action |
| US 10,076,126 (110 C, 4 h, standards co-heated) | -0.07 to -0.64 dex at lambda 0; -0.59 to -1.16 at lambda 0.005 /min; -1.1 to -1.7 at lambda 0.01 | its level is only interpretable with the spike operator |

So both Nestle pots also sit 0.5-1 dex below DFA. Part of CF-068 may therefore be "DFA high" (R68c), not "Reading low".

**Covariate regression (EXPLORATORY).**
- Fit: residual on log10 O2 : (Cys + pentose), n = 4, one point per lab: R^2 0.56, slope -0.73 dex per dex.
- Significance: one-sided permutation p = 0.125 (the minimum possible with n = 4 is 0.042).
- O2 per reductant spans 0.013 (US 8,920,862), 0.065 (DFA), 0.07-0.16 (US 10,076,126), 0.98 (Yiltirak).
- Below the 70 % rule. Lab, concentration and O2 are confounded one-to-one.
- The lowest-O2 pot (US 8,920,862) is not the highest, so O2 alone cannot order the pots.

**Within-study signs.**
- No study varies O2 on thiols at process temperature (evidence map 1.3 still holds).
- Zhang 2024 (BTBU, 115 C), the sulfur-matched cystine vs cysteine arms: the oxidised sulfur raises the disulfide share
  x6.3 (54 % vs 8.6 %) but does NOT lower free MFT (x1.5 higher). The sign agrees for "oxidant -> disulfide", not for
  "oxidant -> less free thiol".

**The IS-timing operator does not separate the two pots.** Both added the standard after heating, so R68b alone cannot
explain the gap. It matters only through the bound share, which nobody measured in Yiltirak's pot.

**Equal-cook series.** Yiltirak's own residual grows from -1.77 (100 C) to -2.17 dex (130 C). That is consistent with a
temperature-activated loss in that vessel (crude: one time shape is used for all four cooks).

**Verdict:** OPEN. The 0.5-0.75 dex (1 SD) between-lab scatter of the absolute MFT levels becomes the declared band
(thiol prereg sec. 6).

| rival | P |
|---|---|
| O2 >= half the residual | 0.45 |
| DFA reads >= 0.5 dex high | 0.35 |
| Yiltirak work-up low | 0.15 |
| other | 0.05 |

**What decides it (designed experiment):** Yiltirak's recipe (25 / 25 mM, pH 5.5, 100 C, 4 h) in 3 mL vs 15 mL fills,
argon vs air, with SIDA for MFT and bis(2-methyl-3-furyl) disulfide. Plus author requests: Elmore for the disulfide or
chromatograms (data "on request"); DFA for the Table IV autoclave fill.

**Model:** an O2-consumer arm (OX charged from the vessel block, k_cys_ox with the Ehrenberg form) is a declared arm,
not fitted. Yiltirak stays a hold-out.

## 2. MFT loss: one constant or a pot-dependent law (EX-050)

**Grosch's in-pot loss, re-derived.** lambda(100 C) is 0.044 (0.035-0.053) /min with thiamine alone and 0.034
(0.030-0.038) with thiamine + Cys. The engine's R3En value at 100 C is 0.041 /min, the same size.

**Transported to 110 C** (Ea 40-133): 0.042-0.163 /min. Against the patent bound (EX-050 red team: 0.02 for slow
exhaustion, 0.05 for fast):
- **L0** (one pot-independent first-order constant): no Ea >= 0 meets the 0.02 bound. The 0.05 bound is met only with
  Ea <= 61 kJ/mol at Grosch's LOW edge, or Ea <= 30 at its centre. The engine's Ea of 101 lies outside both.
  P(L0 serves both pots) ~0.15.
- **L1 / L2** (cysteine suppression, or electrophile capture competed by Cys, lambda / (1 + Cys / K_c)): feasible.
  K_c is >= ~0.9 mM from Grosch (2.1 mM central; the CI includes no effect). The patent caps it at 3-55 mM if the free
  Cys is ~16 mM (TTCA-bound at K 60 L/mol), or at 50-860 mM if the full 250 mM is free.

**Caveat found here.** Grosch's pots are THIAMINE pots. Thiamine cleavage releases a pyrimidinylmethyl electrophile that
consumes MFT (jhoo2002, evidence map 2.4). Grosch's lambda may therefore be an electrophile-capture sink specific to
that pot, which is L2 with E = the thiamine fragment. Its agreement with the engine's pentose-pot constant may be
coincidental.

**Answer: yes, a pot-dependent law.** Proposed form:
`lambda = k_E * E(t) / (1 + Cys(t) / K_c) + k_ox * O2(t)`, with the capture reversible (Esterbauer 1975: 2-alkenal
adducts revert, t1/2 4-6 h at 20 C).

**Data that bound each term:**

| term | bounded by | what is missing |
|---|---|---|
| Cys | Grosch T VIII (weak); the patent (needs free Cys, i.e. TTCA K at 110 C); Grosch T III (confounded with formation) | designed E3 (labelled MFT spiked during heating, Cys 0-150 mM) |
| electrophile | in kind: Elmore 2002 (-1.14 dex MFT with 25 mM linoleate, FFT flat); Xu 2011 (fresh / oxidised lard -1.32 / -2.16); jhoo2002 (thiamine cation); Esterbauer (reversibility, thiolate dependence) | any constant at process T: E3 plus a hexenal or furfural arm |
| O2 | nothing within-study; Ehrenberg 1989 gives the Cys autoxidation form (Cu / Fe, 4 RSH per O2, flat 29-37 C); Yiltirak vs DFA is cross-study only | E3 run under argon vs air |

## 3. EX-050: revise the verdict wording

| proposition | P now | change |
|---|---|---|
| engine k_mft_decay too fast by >= 3x at <= 110 C in a cysteine-rich pot | 0.7 | unchanged in substance (Grosch does not touch it) |
| as a pot-independent constant up to 145 C | 0.3 | down from 0.45 |
| a single constant fits a dilute low-Cys pot and a 0.25 M Cys pot | 0.15 | new |
| the loss is pot-dependent (>= 3x apart at matched T) | 0.65 | new |

**The follow-up must change.** EX-050 planned "k_mft_decay as a declared single-lab band lambda(110 C) <= 0.02-0.03"
pot-independently. That would under-predict the loss 3-6x in Grosch-type pots. Carry the band as Cys-conditional (apply
it at Cys >= ~10 mM), or as L1 with a K_c band of 1-250 mM, until E3.

## 4. Observation operators the engine needs (ranked by the rows they touch)

| # | operator | when | what it does |
|---|---|---|---|
| O1 | standard added after heating (SIDA) | every DFA / Reading / US 8,920,862 row | reports the standing FREE thiol only (no disulfide / adduct) |
| O2 | standard co-heated | US 10,076,126 | reported = native(t) x IS0 / IS(t): run the d3-MFT / d2-FFT spike as a tracer under the same loss law. Mandatory before its shape (CF-013 / 015) or level is scored |
| O3 | SDE during heating (Weenen 2001, Grosch Method II) | rows from such pots | the analyte is stripped at rate >> its in-pot loss and the cumulative amount is reported; volatile precursors (H2S) are stripped too |
| O4 | SDE after heating with protein or thiol present | e.g. Meynier, meat | reducing: y = r (RSH + rho x 2 RSSR), rho ~0.6 in beef (mottram1996b), ~0 in water; r unknown, so scale-free use only |
| O5 | headspace without an isotope standard, lipid present | lipid pots | x 1 / (1 + phi K_ow): <= 0.1-0.4 dex at phi <= 1.5 % |
| O6 | cross-compound SPME / TIC area ratios | any area ratio | not molar: response scales ~K_ow^0.6-1.0, so +-1 dex between furfural and MFT |

No multiplicative 4-5x correction is ever applied to direct-extraction data (the bolton1993 / evidence map 2.11 reading
is wrong; sec. 9).

## 5. CF-015 (MFT/FFT 3.6 vs 0.50 between the Nestle patents): EXPLAINED, P 0.6

**Sugar : cysteine ratio is the main cause.**
- Within study (Xu 2008, 140 C, one lab, RF-free double ratio): tripling ribose at fixed Cys moves MFT/FFT 8.6 -> 0.75,
  i.e. -1.06 dex for 3x sugar. The patents differ by -0.86 dex for 3x xylose. Same sign, similar size.
- Across pots: the 1:1 pots read 3.6-3.8 (US 10,076,126, co-heated standard) and 5.4 (Yiltirak 100 C / 4 h, standard
  after heating). The 3:1 pots read 1.4 (Schieberle 6-12 h), 1.5-1.6 (Hofmann 145 C) and 0.50 (US 8,920,862).
  Three labs.
- The IS-timing rival alone predicts the co-heated pot above ALL post-IS pots. Yiltirak (post-IS) is higher, so that
  prediction fails.

**IS timing is a co-term, not excluded.**
- Taking all 0.86 dex would need lambda_M - lambda_F = 0.0082 /min.
- With FFT's loss from Kumazawa 2003 (0.0115 /min at 121 C, pH 5; 0.004-0.007 at 110 C), that needs lambda_M =
  0.012-0.015 /min, inside the EX-050 bound.
- P(ratio >= half the gap) 0.6; P(IS timing >= half) 0.25.

**Temperature goes the wrong way.** Yiltirak's MFT/FFT falls with T (5.4 -> 1.06), so 95 C should raise US 8,920,862's
ratio, not lower it.

**Shape half:** O2 (the co-heated standard) is required before the US 10,076,126 shape is scored.

## 6. CF-073 (Weenen 2001 vs Hofmann 2000): OPEN, localised. The donor-generic operator is REJECTED

**Per-donor log10(Weenen / Hofmann)** at pH 3 - 3.2 (Hofmann x1.16 for 2 h vs 1 h, from T3):

| donor | log10 ratio |
|---|---|
| MGO | +1.24 |
| GO | +0.89 |
| 3-DG | -0.20 |
| glucose | -0.59 (T III) / -1.29 (T II) |

**The test.** Equal ratios across the three dicarbonyls (R73a: stripping, quantification, a scale) give chi^2 22.5
(sigma_Weenen 0.2) or 11.2 (sigma 0.3), df 2, p = 1e-5 or 0.004. The MGO - 3-DG spread is 1.44 dex. So no donor-generic
operator can explain the gap. Stripping could act only through a donor-specific in-pot PA sink (R73b).

**Concentration** (26 vs 100 mM, a second-order competing donor sink) caps at 0.59 dex. That leaves MGO +0.65 and GO
+0.30 dex unexplained.

**The O2 sign check (Hofmann T4)** gives an aldehyde air / argon ratio of 0.85-0.93. Weenen's boiling pot (O2 ~0) would
gain <= 0.07 dex from it, so O2 is not a cause.

| rival | P |
|---|---|
| R73b, a donor-specific sink avoided by stripping | 0.3 |
| R73c, concentration | 0.3 |
| R73d, MGO / GO speciation or handling | 0.2 |
| R73a, donor-generic operator | 0.05 |
| typo / basis | 0.15 |

**What decides it (designed):** MGO + Phe and 3-DG + Phe at pH 3.2, at 26 and 100 mM, closed vial vs continuous
stripping, SIDA (a 2 x 2 x 2 design). Or a third lab's fed pots.

**Model:** keep CF-073 as a signature, not as rows. Weenen's 3-DG (and C4 / C5 osone) yields are usable cross-lab checks
of the C6 donor (they agree with DFA within 0.2 dex).

## 7. CF-014 (cause half): OPEN, no covariate recorded

Zhou 2022, 2023 and 2023c are all "pressure-resistant bottle, volume / fill / atmosphere not stated", HS-SPME
3 mL + 2 mL NaCl, dichlorobenzene IS, external curves (Table S1 not on disk). No recorded covariate differs. Remedy as
triaged: one lab, a between-paper scale term (explained by remedy). Cause: author request (fill, headspace, storage
time before SPME, the S1 curves).

## 8. CF-031 (lipid lowers vs raises thiols): EXPLAINED in part, P 0.55. Register correction

**Within-study signs:**
- Whitfield 1988, -0.47 dex (lecithin)
- Elmore 2002, -1.14 (linoleate)
- Xu 2011, -1.32 (fresh lard) and -2.16 (oxidised lard)
- Yiltirak buffer + oil, -0.58 to -0.98
- Yiltirak emulsion, -0.07 to -0.27 vs BUFFER (+0.31 to +0.71 vs buffer + oil)
- du2025, +0.10 (free fatty acids)
- US 8,920,862, +1.00 (monoglyceride / soybean mesophase)

**The headspace partition operator** (no isotope standard) can explain at most 0.1-0.4 dex of the 0.5-2.2 dex falls.
So the falls are chemistry: they scale with the electrophile / oxidation load (fresh < oxidised lard; PUFA ester).

**Correction.** The register's "raise them in emulsions (Yiltirak 2026)" is wrong. Yiltirak's emulsion lowers MFT
relative to buffer at every cook and raises it only relative to buffer + oil. Only the patent mesophase (single run,
85 % lipid) and free fatty acids (+0.1) raise MFT.

**Residual rival (P ~0.25): sampling of the non-emulsified oil arm.** A 30 uL aliquot from a creaming oil layer would
under-sample a lipophilic thiol. The MFT fall (-0.98) exceeds FFT's (-0.03 to -0.24) in the direction log K_ow predicts,
but 3M2B (more polar) falls -0.52, so the pattern is mixed. Ask the authors how aliquots were drawn; Table S2 holds the
creaming index.

## 9. CF-016 (furfural x13, same Munich pot): OPEN, probabilities shift

**The outside ratios used on 10-29 are operator-biased.**
- Cerny 2007's FUR/MFT 0.47 is a TIC area ratio from SPME. Corrected for the fibre partition (K_ow^0.6-1.0, log K_ow
  0.41 for furfural vs 1.5-2.0 for MFT), the molar ratio is 2.1-18 (central ~7). That brackets Hofmann's 3.4, not
  Schieberle's 0.26.
- Meynier 1995 (SDE after heating, octadecane in the solvent, RF 1) has no estimable correction.
- So neither outside ratio sides with Schieberle.
- Purge-and-trap at room temperature strips ~1 % of furfural vs most of the MFT. If no furfural isotopologue or matched
  calibration was used, that predicts Schieberle low.

**Revised P:**

| explanation | P now | was |
|---|---|---|
| E1, Schieberle under-recovered | 0.45 | 0.35 |
| E2, a slip | 0.25 | 0.25 |
| Hofmann high | 0.15 | 0.25 |
| E4, run-to-run | 0.15 | 0.15 |

**Decider:** the three-vial Hofmann-pot furfural measurement (as before).

## 10. Pyrazine R (Shibamoto vs Leahy vs Ghent)

**Ghent's 2,5- / 2,6-DMP co-elution.** Leahy (which resolves the two isomers) finds no 2,6-DMP in lysine-glucose /
fructose pots and 2,6 / 2,5 ~1.1 with asparagine. Shibamoto (8 M NH3) has 1.8-4.0. So:
- in Ghent's lysine pots the bias is ~0 dex;
- it is <= 0.33 dex for asparagine-like amines, <= 0.45-0.70 if NH3-like.
- Ghent's R <= 3.5 is not a co-elution artefact for its lysine pots (P 0.75).

**What remains.**
- Lysine-glucose, Leahy 11.7-17.4 vs Ghent <= 2.97: >= 0.60 dex. Candidates: Ghent's theoretical, unmeasured SBSE
  recoveries (R ∝ rec_PZ rec_DMP / rec_MPZ^2; PZ at 0.3 % is the weak link), endpoint vs rate, pH drift to 4-5.
- Shibamoto's R(2,5) of 322-631 (83-138 with 2,6 summed) is a CONDITION (8 M NH3, no amino acid, 2,6 > 2,5), consistent
  with a direct dicarbonyl + NH3 ring route (P2b). It is not an operator (P 0.8).

**2,6 / 2,5 as a fingerprint.** It flags non-aminoketone DMP. But within Leahy it does not order R (Lys 0 vs Asn 1.1 at
similar R), so it is a lead, not a test.

## 11. Data scored / excluded

As in prereg sec. 9. Bolton 1993 is excluded from the level regression: the standard sits in the extract (not the pot),
the matrix is thiamine-dominated, and the work-up is CH2Cl2 extraction concentrated to 0.1 mL (Schieberle 1996 T2:
ether extraction makes MFT "d.n.p."). It is also an O2-rich pot (O2 : Cys ~2), over-predicted by the core.

## 12. Expectations vs outcome

| item | expected | outcome |
|---|---|---|
| CF-068 verdict | OPEN / band P 0.75 | OPEN, met |
| O2 explains >= 70 % | P 0.2 | not met (R^2 0.56, n = 4) |
| Nestle between DFA and Reading | P 0.6 | met (-0.6 to -1.1) |
| L0 confined to a sliver | P 0.7 | met |
| K_c bounds | P 0.6 | met (looser on the upper side) |
| CF-015 Xu sign | P 0.8 | met |
| 1:1 > 3:1 across labs | P 0.7 | met |
| CF-073 donor-generic operator rejected | P 0.8 | met |
| MGO residual >= 0.6 | P 0.75 | met (0.65) |
| CF-031 operator <= 0.3 dex | P 0.8 | met (<= 0.4) |
| register misreading | P 0.9 | confirmed |
| CF-016 Cerny shift | P 0.65 | met |
| Ghent co-elution <= 0.3 dex | P 0.75 | met for lysine pots |

**Unexpected:** the thiamine-electrophile caveat on Grosch, and the "DFA high" component of CF-068 that the Nestle
pots expose.

**Red team:** not yet run on these verdicts. Per lessons 10-31, nothing here moves the engine before one.
