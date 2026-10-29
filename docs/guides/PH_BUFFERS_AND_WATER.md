# pH, buffers and water: the physical chemistry behind the model's conditions

*For a reader who knows undergraduate physical chemistry and wants to know exactly how this model
turns "pH 5, 0.5 M phosphate, a_w 0.98, 145 °C" into rate constants, why each piece is built the
way it is, and where it is still an approximation. Written 2026-10-27, during wave B51. Every number
in the worked examples was computed by the code named next to it; the scripts are cited so the
numbers can be regenerated. Companion documents: the [introduction](INTRODUCTION.md) (what the model
is), the [reaction trees](REACTION_TREES.md) (which steps exist), [SOURCES.md](SOURCES.md) (every
paper).*

1. [Why conditions need physics, not multipliers](#1-why-conditions-need-physics-not-multipliers)
2. [pH is a state of the pot](#2-ph-is-a-state-of-the-pot)
3. [Equilibrium constants at reaction temperature](#3-equilibrium-constants-at-reaction-temperature)
4. [Activities and ionic strength](#4-activities-and-ionic-strength)
5. [Rates through speciation](#5-rates-through-speciation)
6. [Buffer catalysis](#6-buffer-catalysis)
7. [Water activity](#7-water-activity)
8. [What is measured and what is fitted](#8-what-is-measured-and-what-is-fitted)
9. [How a change of physics is allowed into the model](#9-how-a-change-of-physics-is-allowed-into-the-model)
10. [The approximations that remain](#10-the-approximations-that-remain)
11. [Where it lives in the code](#11-where-it-lives-in-the-code)

**Status, 2026-10-27.**

| piece | status |
| --- | --- |
| pH as a computed state (charge balance, conserved strong-ion difference, two pH scales) | in the engine since wave B2.2, sulfur lane only; on GKL / Suleimenov constants since B51 |
| heat capacities of ionisation, corrected cysteine and H2S constants, Davies activities, 25 °C scale for measured pH tables, SID carried across heating segments | **wave B51: SHIPPED 2026-10-27** under the correctness rule (§9), with its refit. It first failed its pre-registered score test (the Hofmann pH 5 pots fall ~100x low); the owner then changed the rule so that corrections ship on verification of the change itself, and the rows it made worse are named in the model card with the assumption they expose (the thiolate-acting thiol sink). |
| speciation on the trunk lane's Amadori steps | planned, wave B52 |
| general buffer catalysis | planned, wave B53 |
| water activity: concentration basis and a glass-transition gate | planned, wave B54 |
| trunk and acrylamide pH today | an empirical log-linear slope about pH 6.8, inside a measured window (§5) |
| water activity today | a measured multiplier table on the trunk's condensation step, a declared band on acrylamide; refused on the sulfur and lipid lanes (§7) |

---

## 1. Why conditions need physics, not multipliers

Almost every rate constant in this model was measured or fitted in a **model pot**: a sugar and an
amino acid dissolved in 0.1–0.5 M phosphate buffer, at a pH set with NaOH, heated in a sealed
vessel. The model is used for **foods**: a pea-protein extrudate, a roasted pulse, a meat analogue,
with whatever endogenous phosphate and organic acids the ingredients carry (beef has
26–35 mmol/kg free phosphate, roughly 0.04 M in its water; Koutsidis et al. 2008,
doi 10.1016/j.meatsci.2007.09.010).

A constant fitted in the pot silently contains the pot's conditions. If the step is catalysed by
phosphate, the fitted constant includes 0.5 M phosphate. If the step needs a deprotonated amine, the
fitted constant includes the pot's amine fraction at the pot's pH and temperature. A multiplier
measured in one system ("the rate is 3.7x faster at a_w 0.69 in skim-milk powder") contains that
system's composition. Carrying any of these to a different system is an extrapolation whose error
is invisible.

The alternative is to model each condition as the physical chemistry it is, with the equilibrium
constants taken from measurement, so that what remains to fit is an **intrinsic** rate constant: the
rate of one named molecular species reacting at infinite dilution. That constant transfers. The
rest of this document is how that is done, piece by piece, and where the data stop it.

## 2. pH is a state of the pot

### 2.1 Electroneutrality determines the pH

Proton transfers between oxygen and nitrogen acids happen in microseconds; the Maillard steps take
minutes to hours. So at every instant of a run the solution is at acid–base equilibrium, and the pH
is not a free variable. It is fixed by charge balance:

```
SID + [H⁺] − [OH⁻] + Σᵢ Cᵢ · z̄ᵢ(pH, T, I) = 0                                         (2.1)
```

- **SID**, the strong-ion difference, is the net concentration of ions that never exchange protons
  in this pH range: the Na⁺ of the NaOH used to set the pH, the K⁺ or Na⁺ of a buffer salt, the Cl⁻
  of an HCl adjustment. **No reaction in the network creates or destroys one of these**, so SID is a
  conserved quantity for the whole run. It is the variable that encodes "this pot was taken to pH 8
  with NaOH and that one was not".
- **Cᵢ** is the total concentration of each titratable pool: the buffer, cysteine, the Amadori
  compound, the formic and acetic acid the network makes, the carboxylates that sink reactions carry
  out of consumed molecules.
- **z̄ᵢ** is the pool's mean charge, from its ladder of acid dissociation constants. For a polyprotic
  acid H_nA with charge z_max when fully protonated, the fraction carrying k fewer protons is

```
x_k = Π_{j≤k} 10^(pH − pK_j) / Σ_m Π_{j≤m} 10^(pH − pK_j),     z̄ = z_max − Σ_k k · x_k    (2.2)
```

Equation (2.1) is monotone in pH, so the model solves it by bisection at every point where it needs
the pH. This is the quasi-equilibrium reduction: [H⁺] is an algebraic function of the slow state,
not an extra stiff ODE.

### 2.2 Two pH scales, and why the distinction is not pedantic

A pH meter reads a **cooled sample at 25 °C**. Every pH printed in this corpus is such a reading.
The rate laws need the pH **in the pot at reaction temperature**, which nobody measured. The two
differ because every pK in (2.1) moves with temperature (§3) and water's own ion product moves by
two orders of magnitude between 25 and 145 °C.

The model therefore uses the published pH exactly once: it solves (2.1) **for SID** at 298.15 K with
the declared composition. From then on SID is fixed, and the pH at any temperature and any point of
the trajectory is solved from it. The final state's pH at 25 °C, solved from the same SID, is the only
model quantity comparable with a published endpoint pH.

**Worked example** (`ph_state.strong_ion_difference`, `solve_ph_and_ionic_strength`; 33 mM cysteine,
the declared pH read at 25 °C, solved at 145 °C; script `doc_numbers.py` in the B51 session record):

| pot | label pH | in-situ pH at 145 °C, B50 constants | in-situ pH at 145 °C, B51 constants |
| --- | ---: | ---: | ---: |
| 0.5 M phosphate | 5.0 | 4.94 | 5.27 |
| 0.5 M phosphate | 7.0 | 6.78 | 7.24 |
| unbuffered | 5.0 | 4.49 | 4.45 |
| unbuffered | 7.0 | 5.69 | 5.36 |

An unbuffered pot labelled pH 7 is **one and a half to two units more acidic** at 145 °C than its
label, because cysteine's thiol and ammonium groups release protons as their pK fall with
temperature and nothing absorbs them. A 0.5 M phosphate pot barely moves, and under B51's constants
moves slightly *up* (§3.2). The label is not the reaction pH, and the error is largest exactly in
the unbuffered systems that foods resemble.

### 2.3 Drift

The network makes acid: formic and acetic acid from the deoxyosones (measured by Martins &
van Boekel), and an organic-acid pool from the lumped sinks. Those pools enter (2.1) as titratable
groups, so as they accumulate the solved pH falls, the pH-dependent rate constants change, and the
acid production changes in turn. The integrator solves this as a fixed point: it integrates with the
pH at 41 time nodes, re-solves the pH from the resulting states, and repeats until the node pH moves
by less than 10⁻⁴ (`sulfur.integrate_sulfur`, "dynamic" mode). One defect fixed by B51: on a
multi-segment heating program SID was re-derived from the label at the start of every segment, which
silently re-titrated away the acid made in earlier segments. SID is a property of the pot, and B51
carries it across segments.

Two constants of the drift model are fitted, because no measurement in the corpus fixes them: the
titratable-acid yield of the lumped sink and the pK of the pentose Amadori compound's secondary
amine. Both are calibrated only on three declared endpoint measurements (Zhou 2023), and everything
else in (2.1) is a tabulated equilibrium constant.

## 3. Equilibrium constants at reaction temperature

### 3.1 The integrated van 't Hoff equation with a heat capacity

For an ionisation HA ⇌ H⁺ + A⁻ with standard enthalpy ΔH° and heat capacity ΔCp°, assumed constant,

```
d ln K / dT = ΔH°(T) / RT²,    ΔH°(T) = ΔH°₀ + ΔCp° (T − T₀)
```

integrates to (Clarke & Glew 1966, doi 10.1039/TF9666200539, to first order)

```
pK(T) = pK₀ + ΔH°₀/(R ln10) · (1/T − 1/T₀) − ΔCp°/(R ln10) · [ln(T/T₀) + T₀/T − 1]       (3.1)
```

- The **enthalpy term** is the familiar one. Amine ionisations are strongly endothermic (glycine
  44.2 kJ/mol, cysteine 34–36 kJ/mol): their pK fall by 1.5–2 units between 25 and 145 °C.
  Carboxylic acids and phosphate are nearly thermoneutral (|ΔH°| < 8 kJ/mol): their pK hardly move
  through this term.
- The **heat-capacity term** is the curvature. For every ionisation here ΔCp° is negative (the ions
  formed order the surrounding water; GKL tabulate −57 to −295 J K⁻¹ mol⁻¹), so the bracket is
  positive and the term **raises** pK at high temperature. Over 120 K it is worth 0.3–0.6 units: as
  large as the enthalpy term for phosphate and acetate.

Before B51 the model used (3.1) without the ΔCp° term. The module already knew this was wrong for
water — constant-enthalpy van 't Hoff misses pKw by 0.49 units at 145 °C, which is why pKw is
tabulated (Marshall & Franck 1981, doi 10.1063/1.555643) — but carried the same error for every
other acid. Equation (3.1) with ΔCp° is itself a truncation: applied to water it is 0.04 units high
at 100 °C and 0.11 at 150 °C (pinned in `tests/unit/test_kinetic_core_b51.py`), against the 0.27 and
0.75 that the constant-enthalpy form misses by.

### 3.2 The numbers, before and after B51

pK at 25 / 100 / 145 °C (`ph_state.TitratableGroup.pkas_at`, `parameters_sulfur.pka_at`):

| equilibrium | B50 (constant ΔH) | B51 | source of the B51 value |
| --- | --- | --- | --- |
| H2PO4⁻ ⇌ HPO4²⁻ | 7.20 / 7.07 / 7.02 | 7.20 / 7.35 / 7.63 | GKL Table 7.51: ΔH 3.6, ΔCp −230 |
| acetic acid | 4.76 / 4.77 / 4.78 | 4.76 / 4.94 / 5.16 | GKL Table 7.2: ΔH −0.41, ΔCp −142 |
| H2S ⇌ HS⁻ | 7.05 / 6.28 / 5.94 | 6.98 / 6.47 / 6.49 | Suleimenov & Seward 1997, measured to 350 °C |
| cysteine thiol | 8.33 / 7.41 / 7.02 | 8.36 / 7.17 / 6.72 | GKL Table 7.23: ΔH 36.1, ΔCp −66 |
| cysteine ammonium | 10.28 / 8.77 / 8.12 | 10.75 / 9.80 / 9.58 | GKL Table 7.23: ΔH 34.1, ΔCp −204 |

Three of these are corrections of errors, not refinements:

- **Cysteine.** The B50 pair (26 and 43 kJ/mol) was cited to Goldberg, Kishore & Lennen 2002 (GKL,
  doi 10.1063/1.1416902). GKL print 36.1 and 34.1. The citation had never been checked against the
  table.
- **H2S.** The first ionisation constant of H2S has a shallow **minimum** near 100–150 °C (pK₁ 6.49 at
  both 100 and 150 °C; Suleimenov & Seward 1997, Geochim. Cosmochim. Acta 61:5187, Table 2 and Eq. 22).
  No constant-enthalpy law can produce a minimum; the B50 law put pK₁ at 5.94 and still falling, which
  overstated the HS⁻ fraction about 3.5-fold at pH 5 and 145 °C.
- **Phosphate.** Because its pK2 *rises* at temperature under (3.1), a 0.5 M phosphate pot is
  slightly more alkaline in situ than its label (table in §2.2), not more acidic as the constant-ΔH
  law said.

### 3.3 Why an ionisation enthalpy contaminates an activation energy

Suppose a step proceeds through a minority species, say the free amine of a molecule whose pK lies
above the pH. Then

```
k_obs = k_int(T) · x_reactive ≈ k_int(T) · 10^(pH − pK(T))                              (3.2)

⇒  Ea,apparent = −R d ln k_obs / d(1/T) ≈ Ea,intrinsic + ΔH°_ion   (at fixed in-situ pH)
```

An Arrhenius fit to k_obs therefore returns the intrinsic barrier **plus** the ionisation enthalpy.
If the model's ionisation enthalpy is wrong by 10 kJ/mol, as the cysteine thiol's was, every
activation energy fitted alongside it absorbs 10 kJ/mol of error, and the error reappears the moment
the model is run at a temperature the fit did not see. This is why a correction to an equilibrium
constant is never "just a pH factor": the fitted constants were calibrated on the wrong speciation
and must be refitted (§9).

## 4. Activities and ionic strength

### 4.1 What pH means, and what Davies does

pH is defined through the **activity** of the hydrogen ion, a_H = γ_H [H⁺], and a meter measures
(approximately) −log a_H. Equilibrium constants are ratios of activities. At the ionic strengths of
the corpus's buffers, I = ½ Σ cᵢ zᵢ² ≈ 0.5–1 M, activity coefficients of ions are far from one, so
concentrations and activities must be distinguished.

The model uses the **Davies equation**,

```
log₁₀ γ_z = −A(T) · z² · [ √I/(1 + √I) − 0.3 I ]                                        (4.1)

A(T) = 1.82483 × 10⁶ · √ρ(T) / (ε(T) · T)^{3/2}
```

with the dielectric constant of water ε(T) from Malmberg & Maryott (1956, J. Res. NBS 56:1) and its
density from the steam tables. A is 0.511 at 25 °C, 0.596 at 100 °C and 0.674 at 145 °C: the
electrostatic non-ideality grows as water's dielectric constant falls.

### 4.2 Mixed constants

The charge balance (2.1) is written in concentrations, with the meter's pH. So each acid step is
used as a **mixed constant**,

```
K′ = a_H [A⁻] / [HA] = K · γ_HA / γ_A    ⇒    pK′ = pK + A·f(I) · (z_HA² − z_A²)          (4.2)
```

with f(I) the bracket of (4.1). For H2PO4⁻ → HPO4²⁻ (z −1 → −2) the shift is −3A·f(I), about
−0.3 units at I = 0.5 M and 25 °C; for a neutral acid going to a monoanion (H2S, a thiol) it is
−A·f(I). In (2.1), [H⁺] = a_H/γ₁ and [OH⁻] = K_w/(a_H γ₁). Ionic strength depends on the speciation
and the speciation on ionic strength, so the two are solved together as a fixed point, converged to
10⁻⁹ M (`ph_state.solve_ph_and_ionic_strength`).

### 4.3 Its limits

Davies is semi-empirical: it is good to about ±0.02 in log γ below I ≈ 0.1 M and degrades above
0.5 M. The rigorous treatment at these strengths is Pitzer's ion-interaction model (Pitzer 1973,
doi 10.1021/j100621a026), which needs specific interaction parameters for every ion pair. None exist
for Amadori compounds, sugar–amino-acid mixtures or melanoidin fragments, so a Pitzer model would
be rigorous in form and invented in its numbers. Davies is declared as an approximation, with its
range stated, rather than replaced by an unparameterised rigour.

## 5. Rates through speciation

### 5.1 The rate law

A Maillard step proceeds through specific protonation states of its reactants: the **neutral amine**
is the nucleophile that attacks a carbonyl; the **hydrosulfide anion** HS⁻ is the sulfur nucleophile;
the **thiolate** R–S⁻ is what oxidises to a disulfide; 1,2-enolisation is acid-catalysed. Writing the
rate on the reactive species,

```
rate = [ k₀(T) + k_int(T) · x_reactive(pH, T, I) ] · [A]_total · ...                         (5.1)
```

where x_reactive comes from (2.2) with the pK of §3 and §4, and k₀ is a parallel uncatalysed route
where one exists. The sulfur lane has carried this form since wave B2.1, with five declared kinds
(acid, base, neutral sulfide, hydrosulfide anion, thiolate; `parameters_sulfur.ph_factor`). Each
factor is normalised to 1 at pH 5, the same temperature and infinite dilution, so a fitted constant
always means "k at pH 5, at this temperature, in an ideal solution".

### 5.2 What speciation predicts that a slope cannot

The trunk and acrylamide lanes still use an empirical slope, k × 10^(s·(pH − 6.8)), measured over a
window (pH 5.5–6.8 for the trunk, from Martins & van Boekel 2003). A slope rises without limit; a
speciation factor **plateaus** once the pH passes the reactive species' pK. At cooking temperatures
that pK sits inside the range foods occupy: the amine of the glucose–glycine Amadori compound (pK 8.2
at 25 °C) falls to about 6.9 at 100 °C. Davidek et al. 2002 (doi 10.1021/jf025561j) measured that
compound's degradation at 90 °C with the pH held by a pH-stat at temperature, from pH 5 to 8: the rate
flattens above pH 7, as speciation predicts and a slope does not. Mossine et al. 1994 measured that amine's pK directly (8.18 at 25 °C). Wave B52 will move the trunk's
Amadori steps to (5.1), with Martins & van Boekel 2005's pH 4.8 and 7.5 columns pre-registered as the
test; Davidek's unbuffered series has already been looked at and is reported as a non-blind check, and
its phosphate series is reserved for the buffer term (§6), because one dataset cannot both fit a term
and test it.

### 5.3 What a pH–rate profile cannot tell you

The rate laws k[HA][B] and k′[A⁻][BH⁺] are **kinetically equivalent**: they produce identical pH–rate
profiles, because [HA][B]/([A⁻][BH⁺]) is fixed by equilibrium. A pH profile alone therefore never
identifies which protonation state reacts. That assignment is a chemical claim (the amine is the
nucleophile because a protonated ammonium has no lone pair) and the model records it as such, with
no number attached; the profile can only test whether the claimed species gives the right shape.

## 6. Buffer catalysis

### 6.1 Specific and general catalysis

The complete rate law for a step catalysed by acids and bases is

```
k_obs = k₀ + k_H a_H + k_OH a_OH + Σ_B ( k_HB [HB] + k_B [B⁻] )                              (6.1)
```

- **Specific** acid–base catalysis (the k_H and k_OH terms) depends only on pH. Speciation (§5)
  captures it.
- **General** catalysis (the buffer terms) depends on the buffer's **concentration** at fixed pH: the
  buffer's own acid and base forms transfer protons in the rate-limiting step. Sugar enolisation and
  the Amadori rearrangement are classic cases, and phosphate is the classic catalyst.

The experiment that separates them is buffer dilution at constant pH, constant buffer ratio and
constant ionic strength: specific catalysis gives a flat line, general catalysis a slope. A rate
measured in 0.5 M phosphate and applied to a food with 0.04 M phosphate carries the difference as
an unseen error, which is why this term matters for transfer more than for fitting.

### 6.2 What can be measured, and the Brønsted law

The Brønsted catalysis law, log k_B = log G + β·pK_a(B), relates catalytic constants of different
bases through their pK. It is what would let a constant measured for phosphate be carried to
pyrophosphate or phthalate. Mottram & Nobrega 2002 (doi 10.1021/jf0200826) found 0.3 M phthalate
catalyses MFT formation about as well as 0.2 M pyrophosphate, both 2.2–2.6x over an unbuffered pot
at pH 5.6, with no gain for FFT: catalysis on the 2,3-enolisation branch, not the 1,2 branch.
Reports of a **sub-linear** dose response (rate ∝ √[H2PO4⁻] for glucose–glycine; Reynolds 1959)
mean a constant measured at 0.1 M is not extrapolated linearly to 0.5 M. Wave B53 will add the
buffer term on the Amadori-degradation and 2,3-enolisation steps, with its catalytic coefficient
taken from Davidek 2002 and Mottram & Nobrega 2002 as the hold-out; the numbers are set in that
wave's pre-registration from the dossiers, not here.

## 7. Water activity

### 7.1 Three physical effects

Removing water does three different things to a Maillard step, and a rigorous model treats them
separately:

1. **Concentration.** The reactants are dissolved in less water, so their concentrations rise. The
   moisture content at a given a_w comes from the sorption isotherm, usually the GAB equation:

   ```
   W(a_w) = W_m C K a_w / [ (1 − K a_w)(1 − K a_w + C K a_w) ]                                (7.1)
   ```

2. **Mobility.** As the matrix approaches its glass transition, diffusion slows and the reaction
   becomes diffusion-limited. Tg depends on water content through the Gordon–Taylor mixing rule
   (Gordon & Taylor 1952, doi 10.1002/jctb.5010020901),

   ```
   Tg(w) = (w_s Tg,s + k w_w Tg,w) / (w_s + k w_w)                                             (7.2)
   ```

   and the rate collapse is usually written with the Williams–Landel–Ferry equation (Williams, Landel
   & Ferry 1955, doi 10.1021/ja01619a008), log(k/k_g) = C₁(T − Tg)/(C₂ + T − Tg).
3. **Water as a product.** Schiff-base and Amadori formation release water, so high a_w pushes those
   equilibria back by mass action.

Together these produce the well-known maximum of Maillard rates at intermediate a_w (0.6–0.8).

### 7.2 Why the rigorous model is not built, and what is

Tested against the corpus, the rigorous decomposition fails where it needs constants nobody has
measured for the matrix in question:

- **Concentration per unit water** (molality in the sorbed water) predicts a 4x drop at a_w 0.98
  that Pereyra Gonzales et al. 2010 (skim-milk powder; doi 10.1016/j.idairyj.2009.07.007) do not see,
  and a 4.4x formation-rate change for asparagine–glucose between a_w 0.82 and 0.92 where
  De Vleeschouwer et al. measured 1.09x. **Concentration per total volume** (whole-phase), with no
  fitted constant, reproduces both.
- **WLF with universal constants** misses Pereyra Gonzales's rates by 10^2.7 to 10^4.9; a fitted
  coupling exponent varies from 0.00 to 0.36 between temperatures in the same matrix.
- **Water as a product** has no measured order or equilibrium constant for these systems.

So the defensible model (wave B54) is: declare the concentration basis (moles per kilogram of dry
matter plus moisture, converted with the whole-phase density); re-derive the trunk's measured a_w
table as the residual after that concentration effect, kept as a band; and use Gordon–Taylor Tg as a
**gate**, warning when the process temperature is within 20 K of Tg and refusing below it, instead of
a WLF multiplier with transplanted constants. For high-moisture meat analogues (a_w ≈ 0.97–0.99,
far above Tg at cooking temperatures) all three effects are small; they matter for roasting, dry
extrusion and storage. The sulfur and lipid lanes refuse a_w until a series exists that a scalar
factor can represent (Hofmann 1998's dry arm moves ribose FFT ×8 and hexose FFT ×0.5, which none can).

## 8. What is measured and what is fitted

| quantity | status | source |
| --- | --- | --- |
| pK at 25 °C, ΔH°, ΔCp° of every buffer and solute | measured, tabulated | GKL 2002; CRC Handbook |
| H2S pK₁ from 25 to 350 °C | measured correlation | Suleimenov & Seward 1997 |
| pKw(T) | measured, tabulated | Marshall & Franck 1981 |
| ε(T), ρ(T) of water | measured | Malmberg & Maryott 1956; steam tables |
| which protonation state each step uses | declared chemistry, no number | §5.3 |
| intrinsic rate constants and barriers | measured where a paper measured the step; otherwise fitted on declared rows | each lane's parameter module |
| sink acid yield, pentose Amadori amine pK | fitted, on three declared endpoint measurements | Zhou 2023, `ph_state` §3 |
| class values (ΔCp of formic acid, of TTCA's groups) | declared assumption, flagged in the source string | `ph_state.py` |

## 9. How a change of physics is allowed into the model

Fitted constants absorb the errors of the physics they were fitted under. A model with corrected
physics and old constants is therefore not a valid model, and neither is one where corrected physics
was adopted because it happened to improve a score. Every change described here goes through the
same procedure (`scripts/generators/WAVES.md`):

1. **Pre-registration**, written before any score exists: what changes, which candidates, the
   expected direction and a probability of passing, and the ship rule. For B51:
   `results/validation/kinetic_core_b51_prereg.md`.
2. **Switches**, so the change is attributable. B51 has six (`ph_state.PhStateCorrections`), each
   scored alone against the frozen fit; all off reproduces the previous wave bit for bit.
3. **Refit** on the declared fit rows with the new physics, from the frozen optimum and from a
   perturbed start.
4. **Hold-out tests** the fit never saw: the Maillard-path hold-out panel, named out-of-sample pots,
   and the directional claims panel. A change that makes them worse does not ship, however correct
   it is in principle, and the record says which test failed and why.

B51 is the worked example, and it changed the procedure. Every one of its six corrections is a measured
constant replacing an error, and refitted, the model predicted two out-of-sample pots ~100x too low,
because the corrected speciation exposed a wrong structural assumption (the thiol sink acting through the
thiolate) that the old, wrong constants had been compensating for. Under a pure score gate it would not
have shipped, and the model would have kept constants known to be wrong because they cancel an error
elsewhere. Since 2026-10-27 the rule (`scripts/generators/WAVES.md`, "How a change is accepted") accepts
a change by kind:

- a **correction** (a measured constant replacing a wrong one, a data fix, a scale or unit error, a bug)
  ships on verification of the change itself — against its source, a derivation or limiting case, and the
  code's invariants — and a score that falls is reported and diagnosed, not used as a veto; the rows it
  made worse are listed in the model card with the assumption they expose;
- a **refit** after a correction ships with it;
- a **modelling choice** with real alternatives still needs a pre-registered hold-out test.

The per-switch attribution is what makes the next wave targeted: it names which correction exposed what.

Data corrections are different: a buffer a paper states (the nine trust-loop Hofmann bundles' 0.5 M
phosphate, Mottram's pyrophosphate) is a fact, applied to every arm, under a dated amendment to
`docs/reference/FIT_HOLDOUT_DECLARATION.md` (Amendment 44).

## 10. The approximations that remain

- **Cysteine's constants are macroscopic.** Its thiol and ammonium ionise in overlapping ranges, so
  GKL's second and third constants mix both groups; the speciation factors take them as thiol and
  ammonium respectively, the same approximation the previous pair made.
- **The thiol pK of MFT and FFT is not measured.** Cysteine's thiol stands in for both, although FFT
  (benzylic-type, pK nearer 9.4) and MFT (heteroaryl, nearer 7) differ in opposite directions
  (`parameters_sulfur.THIOL_PKA_CAVEAT`).
- **Constant ΔCp°.** Equation (3.1) truncates the heat capacity's own temperature dependence (0.11
  units on water at 150 °C).
- **Class values.** Formic acid's ΔCp° (not tabulated in GKL; acetic acid's used) and the ΔCp° of
  the Amadori compound's and TTCA's groups (glycine's used) are declared assumptions.
- **Davies above 0.5 M** (§4.3).
- **Pressure** on the pK is ignored: about 0.01–0.03 units at a sealed vessel's few bar
  (Suleimenov & Seward's correction to saturation pressure is below 0.03 to 300 °C).
- **Headspace.** H2S partitions into the gas phase, which lowers its aqueous concentration; the
  partition is only partly modelled.
- **Buffer catalysis is assumed temperature-independent**, carried with a ×/÷1.5 band, until a
  source measures it at two temperatures.

## 11. Where it lives in the code

| piece | file |
| --- | --- |
| charge balance, SID, two pH scales, buffers, ΔCp, Davies, the B51 switches | [`src/kinetic_core/ph_state.py`](../../src/kinetic_core/ph_state.py) |
| sulfur speciation factors, H2S pK₁(T), the measured pH tables | [`src/kinetic_core/parameters_sulfur.py`](../../src/kinetic_core/parameters_sulfur.py) |
| the dynamic pH fixed point, SID carry, the 25 °C scale for measured tables | [`src/kinetic_core/sulfur.py`](../../src/kinetic_core/sulfur.py) (`integrate_sulfur`, `sulfur_rate_constants_at`) |
| trunk pH slope and a_w table | [`src/kinetic_core/trunk_conditions.py`](../../src/kinetic_core/trunk_conditions.py) |
| acrylamide pH and a_w | [`src/kinetic_core/acrylamide_conditions.py`](../../src/kinetic_core/acrylamide_conditions.py) |
| the tests that pin the constants | [`tests/unit/test_kinetic_core_b51.py`](../../tests/unit/test_kinetic_core_b51.py), [`tests/unit/test_kinetic_core_b2_2.py`](../../tests/unit/test_kinetic_core_b2_2.py) |
| the wave records | `results/validation/kinetic_core_b51_prereg.md`, `scripts/generators/WAVES.md`, [`docs/reference/FIT_HOLDOUT_DECLARATION.md`](../reference/FIT_HOLDOUT_DECLARATION.md) |
