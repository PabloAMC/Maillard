# Pre-registration: wave B20, the glycation arm (written 2026-09-09, before the fit ran)

## 1. Why

The people this tool is for cook protein isolates, and the Maillard chemistry of an isolate is
mostly the glycation of its bound lysine: the sugar joins the protein, the Amadori compound on the
protein breaks down, and two of its products, CML and CEL, are the markers the safety literature
measures and the panel asks for and is refused (`results/validation/data_wishlist.md`, section 3).
The matrix layer (2026-09-08) already charges an isolate's amine pool from measured densities and
uses it only to bind aldehydes. Four dossiers read on 2026-09-08 and 2026-09-09 (`nguyen2016`,
`berk2021`, `troise2015`, `hamzalioglu2026`) hold what a wave needs: one laboratory's rate
constants for every step on casein-bound lysine in water at two temperatures, a second
laboratory's barriers for the two marker-forming steps, and two independent findings that the
glyoxal route to CML fits to zero.

## 2. The arm

Five steps on the trunk lane (glucose and an amine, the lane that carries 3-deoxyglucosone and
methylglyoxal), on four new species: LYSP, the protein-bound lysine residue; FLP, the bound
fructosyl-lysine (Nguyen's AP, furosine's parent); CML; CEL. Appended at the end of the species
table and trunk-only, so the sulfur and acrylamide state vectors keep their shape.

    r_glc_lysp   Glc + LYSP -> FLP                 k_glyc       second order, L/(mmol min)
    r_flp_cml    FLP -> CML + 4 FRAG_C            k_flp_cml    first order
    r_flp_cel    FLP -> CEL + 3 FRAG_C            k_flp_cel    first order (via methylglyoxal, lumped as Nguyen fitted it)
    r_flp_decay  FLP -> TDG + LYSP                k_flp_decay  first order (the Amadori decay that returns the amine)
    r_cml_loss   CML -> 8 MEL_C + 2 MEL_N         k_cml_loss   first order

The glyoxal + lysine route to CML is not written: Nguyen 2016 (casein, water), Berk 2021 (sesame,
dry) and Hamzalioglu 2026 (milk) each fit it to zero or a thousandfold below the Amadori route.

**The pool.** LYSP is charged from the spec's protein loading and the matrix table's amine density
(`data/species/protein_matrices.yml`: soy 0.36, pea 0.47, β-lactoglobulin from its sequence, in
mmol per gram) times the declared available fraction, the centre of the table's availability band
(0.4 to 1.0, so 0.7); the band's corners travel as the interval on every answer. Without a
loading the pool is zero and the five steps carry no flux, so every earlier pot reproduces bit for
bit. A request for CML, CEL or fructosyl-lysine without a protein loading is refused by name:
free lysine resolves to the acrylamide lane as an amine and is not this arm's substrate.

**What fixes the numbers.** Nguyen 2016 Table 1, system M1 (sodium caseinate 30 g/L, about 16
mmol/L lysine residues, glucose 150 mmol/L, 0.1 M phosphate pH 6.8, 120 and 130 °C, sealed, air):
k3 (glucose + lysine residue → AP) 1.5e-4 ± 3.0e-5 and 1.6e-4 ± 2.2e-5 L/(mmol min); k7 (AP → CML)
8.8e-3 ± 6.6e-3 and 6.0e-3 ± 1.6e-3 per minute; k9 (AP → CEL) 2.3e-3 ± 2.7e-3 and 2.0e-3 ± 3.0e-4;
k8 (AP → other products) 5.2e-2 ± 3.3e-2 and 1.5e-1 ± 3.4e-2; k11 (CML loss) 2.9e-1 ± 2.7e-1 and
7.7e-2 ± 3.3e-2. Ten rows, one per constant per temperature, each weighted by its printed interval
(sigma in log10 = log10(1 + half-width / value); a row whose interval spans zero gets 0.5).

**What is fitted.** Five coordinates: the log10 of each constant at the trunk's 100 °C reference.
Bands: two decades either side of the prior centre (Nguyen's values brought to 100 °C with the
declared barriers). Two starts (the prior centre; a seeded perturbation), scipy least squares on
the log10 residuals, then the Laplace covariance at the optimum.

**What is declared.** No barrier is fitted: the two temperatures are ten degrees apart with
overlapping intervals, and the Q10 values run from 0.3 to 2.9, which the authors call apparent.
Each step takes a measured barrier from the nearest measured step: glycation, the trunk's Amadori
formation (Martins 2005, 96.8 kJ/mol); the Amadori decay, the trunk's Amadori → 3-deoxyglucosone
(97.1); CML formation, Berk 2021's fructosyl-lysine → CML (113); CEL formation, Berk 2021's
methylglyoxal + bound lysine → CEL (92); the CML loss, flat (Nguyen's pair falls with temperature,
which no barrier reproduces; flagged, the Kocadagli glyoxal-sink precedent). The available fraction
is declared as above. The rates are one laboratory's on casein; a plant isolate's lysine may
glycate at a different rate, and the wishlist says so.

## 3. What runs

The generator `generate_kinetic_core_b20_fit.py` compares each fitted constant, evaluated at the
row's temperature with its declared barrier, with Nguyen's printed value; no integration is needed
for the fit because the rows are the constants themselves. Beside the fit it integrates Nguyen's
pot (LYSP 16 mmol/L charged directly, glucose 150, pH 6.8, 120 and 130 °C, 30 min) and reports the
CML, CEL and fructosyl-lysine levels and the lysine loss, and the comparators below.

Hold-outs and comparators, never in the objective: Nguyen's printed CML level range (0.025 to
0.135 mmol/L over 0 to 30 min at 120 to 130 °C, the same pot, end-of-cook levels); Berk 2021's
fructosyl-lysine → CML constant at 180 °C (5.54e-3 per minute, a dry seed); Hamzalioglu 2026's
lactulosyl-lysine → CML constants in milk (1.7e-4 to 3.6e-3 per minute, 110 to 140 °C); Troise
2015's lysine loss in expanded soybean (about a quarter in 60 min at 110 °C, a moist solid the
engine cannot charge, so a direction only); the panel's CML and CEL row (a proxy operating point
with no cook; reported, not decisive).

## 4. What counts as success, declared before the run

- **T1, the rows.** Every row whose printed interval does not span zero within 0.3 dex of its
  value (k3 both temperatures, k7 at 130 °C, k8 both, k9 at 130 °C, k11 at 130 °C); the rest
  reported.
- **T2, the panel is untouched.** No currently scored panel row moves by more than 0.05 dex (the
  arm is inert without a protein loading; the two isolate pots on the panel state a matrix but no
  loading).
- **T3, the levels.** In Nguyen's pot at 120 °C and 30 min, the model's CML within the printed
  range 0.025 to 0.135 mmol/L, and CEL below CML (the paper's ordering).
- **T4, another laboratory.** Berk's 180 °C constant and Hamzalioglu's 110 to 140 °C constants
  against the model's, in decades; reported.
- **T5, identification.** Laplace sigma below one decade on every coordinate, none on its bound.

Ship rule: SHIP if T1, T2 and T5 hold; T3 and T4 are reported. If it ships, the engine reads the
frozen literals (the B7 pattern), the answer for CML, CEL and fructosyl-lysine on a loaded pot
carries the availability interval and the caveat, and the roadmap's programme 7 success criterion
("the CML row on the panel becomes evaluable") is checked against what the panel's proxy row can
support. If it does not ship, the steps stay in the network at their frozen literals with a
DO-NOT-SHIP note and the arm is refused by name, as the pyrazine step would have been.

## 6. Outcome (2026-09-09, run the same day) — SHIP

**What was built.** The four species and five steps of section 2, trunk-only and appended last; the
bound-lysine pool charged from the protein loading through the matrix table (the band centre as the
available fraction); the request for a glycation target without a loading refused by name; the
constants as frozen literals in `parameters_glycation.py` asserted equal to the report by
`tests/unit/test_kinetic_core_b20.py`. Generator `generate_kinetic_core_b20_fit.py`, ship rule
`generate_kinetic_core_b20_ship_rule.py`.

**What happened.** Cost 17.2 on ten rows (reduced chi-square
3.4), both starts agreeing; every coordinate identified (Laplace sigma
0.09 to 0.25 decades), none on its bound. T1 passes: the decisive rows within 0.19 dex. The two
reported rows, k7 and k11 at 120 °C, sit 0.46 and 0.45 dex below the model, inside their own printed
intervals: they are the two constants whose pairs fall with temperature, and a declared positive
barrier cannot follow a pair that falls. T2 passes: no scored panel row moves. T3 passes: in
Nguyen's pot at 120 °C and 30 min the model's CML is 0.053 mmol/L against the printed 0.025 to 0.135,
CEL 0.046 below it, 14 % of the lysine consumed. T4, reported: Berk 2021's dry-seed constant at 180 °C
is +1.7 decades from the model's (a sesame seed with sucrose in an open dish is not casein in water,
and the model's 180 °C value is a 100 °C anchor carried up Berk's own barrier); Hamzalioglu 2026's milk
constants at 110 to 140 °C are within 1.2 decades, 0.07 at 120 °C. Verdict by the rule: SHIP.

**What every answer carries.** The availability interval (the matrix band's corners), the caveat
that the rates are one laboratory's on casein and the barriers declared, and the 100 °C extrapolation
warning for any pot below 120 °C. What the arm does not claim: an absolute for a plant isolate (no
plant-protein glycation rate exists on disk; the wishlist names it), anything about free lysine (it
resolves to the acrylamide lane), or the panel's CML row (a proxy operating point with no protein
loading; still refused, now by the arm's own reason). The roadmap's programme 7 criterion, "the CML
row becomes evaluable", is therefore met in the engine and not on the panel: the panel needs a
benchmark with a stated protein loading and a real cook, which the wishlist now asks for.

