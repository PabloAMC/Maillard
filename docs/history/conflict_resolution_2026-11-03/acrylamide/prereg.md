# Pre-registration — acrylamide cluster (CF-039, CF-104), 2026-11-03

Written before any test below was computed. Integrity reads done first (they are inputs, not tests): Gökmen 2012
Figs. 2-3 (p. 170-171, 300 dpi), Bassama 2011 Table 2 (p. 1457, 200 dpi) and its water contents (70 / 41 / 8 % wb,
p. 1453), De Vleeschouwer 2008 Table 1 (water contents per a_w) and Fig. 2 (p. 6464, 110 dpi: the "mM" axis),
Palazoğlu & Gökmen 2008 Eqs. 7-8, Knol 2005 Table 1 Ea (dossier). Scratch: `scratchpad/conflicts/acrylamide/`.

## CF-039 — formation Ea: 84 (Gökmen 2012) vs 154-182 (Leuven, Bassama)

Corpus of formation Ea, sorted by what is estimated (read for this pre-registration):
- (A) single-response lumped "precursors → AA" with co-fitted first-order elimination, measured T(t): Claeys 2005
  168.3 (aqueous 0.01 M, pH 6); De Vleeschouwer 2007 153.6-182.2 (powders a_w 0.34-0.92); Bassama 2011 168.5-177.8
  (plantain, a_w 0.43-0.97). Also Gökmen & Şenyuva 2006 via Palazoğlu 2008 Eq. 7: 93.8 (fructose-Asn, A → C → D).
- (B) multiresponse step constants: Knol 2005 k4 94.4, k3 102 (aqueous 0.2 M); De Vleeschouwer 2008 k_F(Int1 → AA)
  84.6-171.5 over 8 columns, k_INT 105-126 (potato-based); De Vleeschouwer 2009 I k_INTg 117.5, k_F 159.2 (the lane).
- (C) Gökmen 2012 83.94: ln "k (1/min)" vs 1/T, 90-180 °C, silica-supported, k derivation not stated.
- (D) not measurements: Parker 2012b (fixed 100), Chan 2020 (44.5, tuned), Balagiannis 2019 (none printed).

Hypotheses, tests, decision rules:
- **H1 (medium / water state)**: the Ea depends on the reaction medium (aqueous vs low-moisture vs silica).
  Test T1: within estimand class (A), regress Ea_F on a_w / medium over the 15 rows outside Gökmen (Claeys 5 systems
  at a_w ~1, De Vleeschouwer 2007 six a_w, Bassama three a_w). Rule: H1 is supported only if Ea_F changes by
  ≥ 40 kJ/mol across a_w 0.34 → 1.0 or between media within one estimand; falsified if the spread within (A) across
  media is < 30 kJ/mol (half the gap).
- **H2 (estimand: Gökmen's "k" is a yield-per-time under precursor depletion, not a rate constant)**. Test T2a
  (integrity): reconstruct Gökmen's plotted ln k from Fig. 2 yields as (AA / C0) / t_hold; rule: H2's premise holds if
  that reproduces the four plotted ln k within 0.3. Test T2b (outside the conflict): apply the same estimator to
  Bassama Fig. 3 (time courses whose rate-constant Ea_F is 168-178 on the same data) at a 60-min hold and at the
  earliest common time; rule: H2 can explain the gap if this estimator returns ≤ 110 kJ/mol on data whose rate
  constants say ≥ 168. Test T2c: depletion-corrected Ea for Gökmen (Y = Ymax(1 − e^{−kt}), Ymax from the 180 °C
  plateau); rule: if corrected Ea ≥ 140, H2 explains the conflict; if ≤ 115, H2 explains at most part.
- **H3 (heat model: apparent vs intrinsic)**, cf. sprague2011 (75 vs 106 kJ/mol from one data set). Test T3: a
  come-up lag of 1-3 min (sealed tube, 0.5 g silica, oil bath, no T(t) logged) applied to Gökmen's holds; rule: H3
  matters if the Ea moves by ≥ 20 kJ/mol for a 3 min lag.
- **H4 (laboratory / scheme split)**: what remains is a split between single-response lumped fits and the other
  estimands, not a medium effect. Not testable inside the corpus; recorded with its decider.
- Falsifiers: an (A)-class Ea_F < 110 in any medium would break H1's "no medium effect" reading; Gökmen's ln k not
  reproducible from Fig. 2 would void T2.

## CF-104 — Bassama's ×8 k_F at low a_w vs the lane's flat formation multiplier

- **H5 (concentration basis)**: on the basis the lane and De Vleeschouwer use (moles per whole sample, "M" ≈ mol/kg of
  wet sample; checked on DV 2008 Fig. 2: Asn0 ≈ 2650 mM at 29.6 % water, where per-water molarity would be ~7.6 M),
  Bassama's pseudo-first-order k_F = k2 · [reducing sugar] with sugars ~600 × in excess, and [sugar] per wet kg ∝
  (1 − w). Test T4: predicted k_F ratios 1 : (0.59/0.30) : (0.92/0.30) vs Table 2's 1 : 5.1 : 8.1. Rule: H5 explains
  an arm when the residual (observed / predicted) is within the Table 2 95 % CI of the ratio (≈ ×0.55-1.8 at 0.904 /
  0.430 relative to each other). Per-dry-matter and per-water bases are computed too (rule: report which basis, if any,
  makes the residual flat).
- **H6 (moisture range, not a_w)**: the lane's window is in a_w (0.34-0.99), but De Vleeschouwer's powders hold
  ≤ 34.2 % water at a_w 0.99 while plantain holds 70 % at 0.972. Test T5 (outside the conflict, DV 2008 basic and
  potato systems, 10.4-34.2 % water): lumped formation flux k_INT·partition per wet basis across moisture; rule: H6
  holds if DV is flat (within ×2) over 10-34 % water and Bassama's residual is flat between its two arms inside that
  moisture range (41 % and 8 %) but not across 41 → 70 %.
- **H7 (elimination)**: Bassama's k_E rises ×6.1 / ×14.4 dry-side; the lane's elimination multiplier falls (0.70 at
  a_w 0.43). k_E is first order, basis-invariant; on a bimolecular (matrix-amine) reading it scales with (1 − w).
  Test T6: residual on both readings; rule: if the residual is > ×3 on both, the elimination conflict is real (stays
  open).
- Decision for the lane: change `AW_MULTIPLIER` or its band only if T4-T5 leave a residual outside the band
  (0.41, 1.39) at moisture the lane claims to cover. Never tune to Bassama; at most a band edge or a window statement.
