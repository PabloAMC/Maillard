# Pre-registration: wave B25, the thiol sink, variant (c): irreversible addition to the pot's own unsaturated carbonyls (written 2026-09-09, before the fit ran)

## 1. Why

Both pre-registered sink structures of B17 were refused (`kinetic_core_b17_prereg.md`, sections 6
and 6b): the reversible disulfide because the model makes almost no disulfide, the saturable
thioether because its source was dead at the optimum and its measured equilibrium releases the
thiol when hot. Section 6b named the third candidate the lipid papers supply: an IRREVERSIBLE
addition of the thiol to unsaturated carbonyls, the adducts that halve the thiols when a lipid is
present (Farmer 1990: MFT to 0.15 to 0.40 of the lipid-free pot with 6 to 15 g/L phospholipid;
Whitfield 1988: MFT 3.0×, FFT 2.0× lower with lecithin), on a pool sourced from the flux the lane
actually carries. In a pentose and cysteine pot without lipid the unsaturated carbonyls are the
pot's own: the two deoxypentosones, the live intermediates the thiols are made from, whose
1,2-enol forms are Michael acceptors and which the lane carries at millimolar levels while the
sugar lasts. A thiol that adds to them is removed for good, and the removal stops when the sugar is
spent, which is the shape the reference pot shows at 100 °C (thiols still rising at twelve hours)
and the first-order sink cannot.

## 2. The structure

Four steps on the sulfur lane, no new species: the two meaty thiols add to the two deoxypentosones,
irreversibly, and the adduct joins the terminal thiol-oligomer pool:

    ch_add_mft_dpo   MFT + DPO -> OLG + 5 FRAG_C     k_add     second order, L/(mmol min)
    ch_add_mft_tdp   MFT + TDP -> OLG + 5 FRAG_C     k_add
    ch_add_fft_dpo   FFT + DPO -> OLG + 5 FRAG_C     k_add
    ch_add_fft_tdp   FFT + TDP -> OLG + 5 FRAG_C     k_add

One shared constant and one barrier, both FITTED: log10 `k_add` at the lane's 145 °C reference,
band −6 to 0, and its barrier within 10 to 120 kJ/mol (no thiol-Michael barrier on an osone is
measured anywhere on disk; the matrix layer's aldehyde brackets carry 15 to 30 kJ/mol). The
existing first-order thiol sinks (`k_mft_decay`, `k_fft_decay`) stay free in the objective as in
B16: the fit decides how much of the loss the addition takes from them. At the inert default
(`k_add` = 0) the four steps carry no flux and every earlier wave reproduces bit for bit.

## 3. What runs

B16's objective (B9's 54 primary-evidence rows, the seven Schieberle 2000 Table IV within-study
ratios at 100 °C, Zhai 2021's three TTCA rows; 64 rows), every B9 band kept, the thiol-sink ceiling
at 102 kJ/mol; the two new coordinates appended to the vector as B11 and B17 appended theirs (25
free). Two starts (B9's optimum with the addition at the centre of its band; B8's perturbation
protocol), the 600-evaluation budget, quick mode for the search, the Laplace covariance at the
optimum and the cost slice along each new coordinate. Hold-outs never in the objective: the four
Hofmann 1998 pH-5 level bundles, Yiltirak 2026's four pots, Wang 2022's 140 °C series, Bolton 1994,
the returned Hofmann pH-3 and pH-7 rows, Zhou 2023's dimer shares.

## 4. What counts as success, declared before the run

The six tests of B17 (`kinetic_core_b17_prereg.md`, section 4), unchanged, with T6 on the two new
coordinates:

- **T1** the reference pot at 100 °C: MFT and FFT at 360 and 720 min within 0.5 dex of Table IV as
  ratios to the 30-minute point, and both still rising between 6 and 12 hours.
- **T2** every B9 fit row within 0.3 dex of its B9 residual.
- **T3** Zhou 2023's dimer shares within 0.3 dex (Zhang 2024's reported as figure-derived).
- **T4** Yiltirak 2026's median fold below 20.
- **T5** Wang 2022's 140 °C decline from the peak by less than one decade.
- **T6** Laplace sigma below one decade on log10 `k_add` and below 60 kJ/mol on its barrier, neither
  on its bound, the slices not flat.

Ship rule: SHIP if T1, T2 and T6 hold; T3 to T5 are reported. If it ships, the engine reads the B25
report for the sulfur lane (the first change of the shipped thiol parameters since B9) and the
introduction's section 7 is rewritten. If it does not ship, the record joins B17's two and the
laboratory experiment in the introduction's section 8 is the remaining route.

## 6. Outcome (2026-09-09, run the same day) — DO NOT SHIP; the third structure is refused as the first two were

**What was built.** The four addition steps on one constant with its own barrier, inert at zero;
the `thiol_addition` block through the engine and the Laplace vector; generator, ship rule, tests.

**What happened.** Both starts converged to B16's optimum (cost 931.27 against 930.98) with the
addition driven to its floor (log10 k_add -6.00) and its barrier to its ceiling (120 kJ/mol),
which is the optimiser's way of switching the step off at every temperature: sigma unbounded, both
slices bound-limited. The addition helps no row of B16's objective. T1, T2, T3 and T5 fail as for
both B17 variants (the reference pot's MFT still peaks at six hours and falls, the fed-ribose row
still moves 1.9 dex); T4 passes (Yiltirak median fold 13.7). Verdict by the rule: DO NOT SHIP.

**What the three refusals say together.** Three sink structures have now been offered to the same
64 rows and refused for three different reasons: the reversible disulfide because the model makes
almost no disulfide (oxidant-limited); the saturable thioether because its source is dead at the
optimum and its equilibrium releases the thiol when hot; the irreversible addition to the pot's own
osones because the objective does not want it at any temperature. What they share is the objective:
18 rows at 145 °C against 13 at 100 °C, and inside those, **12 fed-intermediate mol % rows from one laboratory against 6 Schieberle within-study ratios**.
*(Count corrected 2026-09-09 by enumerating the objective's own row table. This sentence read "54 fed-intermediate rows at 145 °C … against seven within-study ratios at 100 °C" until then: 54 is the size of the whole B9 block, not of its fed rows, and the seventh Schieberle ratio is attached to the 145 °C system, not the 100 °C one. The imbalance is real and it is 2:1 on the decisive rows, not 8:1. The argument stands; the number that was carrying it did not.)*
Any structure that slows the 100 °C loss also perturbs the 145 °C fed pots, and the fed pots win by
weight of number. The next step is not a fourth structure on the same rows; it is either the
experiment of the introduction's section 8 (the reference pot at 100 and 140 °C with the removal
measured on fed thiols alone) or a pre-registration that re-weights the objective by laboratory
rather than by row, which is a change of rule the owner would have to make, not a wave. The steps
stay in the network at zero and the record is kept.

