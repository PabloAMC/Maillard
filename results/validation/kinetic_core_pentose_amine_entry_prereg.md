# Pentose + amino acid entry and the pentose acetic-acid step (pre-registration, design; 2026-10-30)

*Written before any code or number for it exists. Part of the shared trunk (plan sec. 7, D3, owner 2026-10-28:
"pentoses and hexoses, any amine"). Switch: off until it passes. Owner, 2026-10-30: "you decide on the priorities".*

## 1. Why

Two findings of 2026-10-30 have one cause.

1. **The sulfur lane has no pentose + amino acid condensation.** Pentose is consumed only through cysteine (TTCA,
   the cysteine Amadori compound, the retired catalytic steps), through a pre-made alanine-xylose Amadori compound
   (`ARP`, fed by Zhou's pots), or by caramelisation / thermolysis. Hexoses have a generic glycine route
   (`r_schiff`); pentoses do not. In Davidek 2008's xylose + glycine pots (0.1 M each, 0.2 M phosphate, 120 C) the
   lane leaves 64 % of the xylose after 4 h where 2.5 % (pH 6) and < 1 % (pH 8, by 60 min) remain, and its pH 6 and
   pH 8 runs are identical. Desclaux 2006's xylose + glycine norfuraneol is 1 dex low in the joint fit (sec. 8 of the
   thiol prereg). Hofmann 1999 (xylose + alanine) cannot be simulated at all (CF-067 table, `engine_species_missing`).
2. **The pentose acid supply is missing.** In Zhou 2023's unbuffered ARP pots the lane's final pH is right at initial
   pH 6 and 1 / 2.2 units high at pH 7 / 8 (thiol prereg sec. 7, B1): it makes ~3.5 mM net acid at every pH and has
   no pentose route to acetic or formic acid. Davidek 2008 (isotope dilution) measures acetic acid from xylose +
   glycine rising with pH: 16 mol% at pH 6 and 28 mol% at pH 8 after 4 h at 120 C, 77 - 89 % of it from xylose C-1 /
   C-2 -- the hydrolytic beta-dicarbonyl cleavage of 1-deoxypento-2,4-diulose (the 1-deoxypentosone tautomer).

Both matter for the goal: real alternative-protein matrices are pentose + many amino acids (cysteine is a minor
one), and the in-situ pH of near-neutral pots is upstream of every pH-dependent odorant step. The thiol layer's pH
arm (thiol prereg R2c) waits on (2).

## 2. What changes (switch on)

- **New species `GARP`**: the pentose Amadori compound of a non-sulfur amino acid formed in the pot. `ARP` stays the
  fed alanine-xylose compound of Zhou's pots, unchanged.
- **`r_pent_amine`: PENT + Gly -> GARP**, second order. `Gly` stands for the free non-sulfur alpha-amino acid
  (glycine; alanine only in a declared sensitivity arm, sec. 4). Rate: fitted (`k_pent_amine`, at 100 C / pH 6), with
  a pH slope (`ph_slope`, decades per pH unit on the in-situ pH, as sub-model 1) and Ea held at the hexose
  condensation's measured value (Martins 2003; a declared hexose -> pentose transfer, band +-15 kJ/mol).
- **`GARP` breakdown**: the four `ARP` steps (`r_arp_dpo`, `r_arp_tdp` and their thermal twins) with ARP's shipped
  constants, releasing the amino acid as `Gly` (not CBX + FRAG_N: the amine is regenerated). Declared assumption: an
  Amadori compound's breakdown rate does not depend on which amino acid it carries, within the band between
  Desclaux 2006's measured glycine-xylose constants (ARP -> 1-DP 2.79e-2 /min, -> 3-DP 5.50e-4 /min, pH 6, 100 C) and
  ARP's shipped ones (checked numerically in the probe, sec. 3; if they differ by > x3, that is a finding and the
  measured Desclaux constants are used for `GARP`).
- **`r_dpo_aa`: DPO -> AA + C3 fragment** (hydroxyacetone, FRAG_C 3), hydrolytic beta-cleavage of the 2,4-diulose,
  base-catalysed (HO- attack): fitted rate at 100 C / pH 6 with its own pH slope; Ea held at the carbonyl family's.
  The minor C-3 / C-4 cleavage (9 - 15 % of acetic acid) and the parallel glycolic acid (C-4 attack; "must be
  formed", not measured) are left out and declared.

## 3. Steps (the usual order)

1. **Cheap probe (no fit).** Switch on at a grid of `k_pent_amine` and `k_dpo_aa` (three values each, around the
   values the data imply: xylose half-life 12 min at pH 6 -> k ~ 6e-4 L/(mmol min) at 120 C; acetic 16 mol%):
   residual tables for every pot with PENT + Gly or ARP, the Desclaux constant check, and the Zhou 2023 final pH.
2. **Design review** (independent, before any fit).
3. **Fit** (4 coordinates: `k_pent_amine`, its slope, `k_dpo_aa`, its slope) on: Davidek 2008 Fig. 1 (acetic acid
   and residual xylose, pH 6 / 8: 19 rows, absolute) and Desclaux 2006 Fig. 2 (norfuraneol, 3-DP; the deoxyosone
   table's tau 0.3 scale for 3-DP). Engine state: the joint thiol fit's s042 optimum, route switches on.
4. **Scored, not fitted:** Zhou 2023 final pH (6 rows; ARP-fed, so it tests `r_dpo_aa` and the pH model, not the
   entry), Cao 2017 / Zhao 2017b pentose pH (sub-model 1's table), Davidek 2006a / 2006b glucose + glycine acetic acid
   (the hexose side, `r_odg_aa`: a check that nothing hexose moved). Re-score the joint thiol fit's rows at the new
   optimum (the entry changes no thiol pot without glycine; `r_dpo_aa` drains DPO everywhere, so MFT and MGO move).
5. **Sensitivity arm (not fitted):** alanine mapped to Gly, Hofmann 1999 xylose + alanine methylglyoxal and
   glycolaldehyde scored.

## 4. Done-tests (raw residuals; declared now)

- Davidek 2008 acetic acid and xylose: RMS <= 0.15 dex per pH (isotope dilution; the xylose method unstated, so 0.2).
- Desclaux 2006 norfuraneol mean within 0.3 dex (from -1.0).
- Every new constant off its bound, Laplace sigma < 1 dex; both slopes inside [-0.5, +1.5].
- Zhou 2023 final pH (scored): |mean| <= 0.3 pH and the pH 6 -> 8 ordering right (measured 3.23 / 3.43 / 5.05).
- Joint thiol rows re-scored: no data set worse by > 0.15 dex RMS (else the arm goes back to the joint fit with
  `k_dpo_aa` free there).

## 5. Expected (written now)

- `k_pent_amine` lands near 6e-4 L/(mmol min) at 120 C / pH 6 (within x3: P 0.7); its slope positive, 0.2 - 0.6
  (P 0.6; Davidek's xylose falls ~3x faster at pH 8 than at 6, less than the free-amine fraction's 100x).
- Acetic acid at pH 6 / 8 within 0.15 dex once `r_dpo_aa` competes with DPO -> norfuraneol (P 0.5): the 1-DP pool
  must feed both norfuraneol (Desclaux) and acetic acid (Davidek); if one cannot be met without breaking the other,
  that is the finding (a missing 1-DP source, or the C-2 / C-3 split is not through free 1-DP).
- Zhou 2023 final pH at initial pH 8 moves from 7.24 to within 0.5 of 5.05 (P 0.35): ARP-fed pots get `r_dpo_aa` but
  not the entry, and 20 mM ARP x 28 % = ~6 mM acetic acid may not be enough against the added base.
- MFT in the joint thiol pots falls where DPO is drained (P 0.7 that some data set worsens by > 0.15; the done-test
  then sends `k_dpo_aa` into the joint fit).

## 6. What it does not do

No formic acid from pentoses (unmeasured near pH 5 - 8 in every source found: Novotny 2008 is persulfate / NaOH
only); no amino-acid-specific condensation rates (glycine stands for all); no protein-bound lysine (the matrix layer).
## 7. Cheap probe (step 1) -- done 2026-10-30, nothing fitted (`kinetic_core_pentose_amine_entry_probe.md`)

- **The declared assumption failed its check, and the declared alternative applies.** ARP's derived enolisation
  rates at Desclaux's conditions (100 C, pH 6) are x4.1 (to 1-DP) and x9.1 (to 3-DP) slower than Desclaux 2006's
  measured glycine-xylose constants: PAMA runs at Desclaux's level with ARP's pH and temperature dependence
  (`pama_dpo_factor` 4.14, `pama_tdp_factor` 9.13; computed by the probe from the lane's own rates, not typed in).
  Whether ARP's shipped constants are themselves slow (they rest on Zhou's ARP pots) or the amino acid matters is
  open; it is recorded, not acted on.
- **Grid** (both slopes 0.3; entry 3e-5 - 2e-4 L/(mmol min), acetic 1e-2 - 3e-1 /min at 100 C / pH 6):
  - Davidek 2008 acetic acid: from absent to -0.07 / 0.42 dex (mean / RMS) at acetic 0.3 /min -- only a FAST
    beta-cleavage reaches 16 - 28 mol%.
  - Davidek 2008 xylose: best near entry 7e-5 (+0.25 / 0.55): the level is reachable, the shape is not at slope 0.3
    (the fast first 15 min) -- the entry's slope must be fitted.
  - Desclaux 2006 norfuraneol: from -0.99 to +0.5 - +1.3 -- the entry now OVER-supplies 1-DP to norfuraneol unless
    the acetic step takes its share (+0.53 at acetic 0.3); 3-DP stays +1.2 (the 3-DX excess every table shows).
  - Zhou 2023 final pH: +0.82 -> +0.42 pH (acetic 0.3; RMS 1.36 -> 0.80).
  - Glucose + glycine acetic acid (Davidek 2006a / b): unchanged, as designed (2006b's -1.29 is the hexose side's own
    shortfall, `r_odg_aa`, not this arm's).
- **Consequence for the design (sec. 3 step 3 superseded):** at ~0.3 /min the beta-cleavage takes most of the 1-DP,
  the MFT precursor (via norfuraneol) in every pentose + cysteine pot, so it cannot be fitted apart from the thiol
  layer. The fit is a JOINT-fit arm: the s042 coordinates plus the four entry coordinates (entry rate and slope,
  acetic rate and slope), on the joint fit's rows plus Davidek 2008 (19 rows) and the Desclaux pot (already in);
  scored as sec. 3 step 4. Done-tests unchanged (sec. 4), the thiol re-score now read on the joint fit's own tables.
