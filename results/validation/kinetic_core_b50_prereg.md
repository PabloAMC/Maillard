# Wave B50 pre-registration — thiamine releases H2S (a measured share, no fit)

*Written 2026-10-27 before any B50 score existed. Code: `parameters_sulfur.thiamine_h2s_share`, `sulfur.py`
`r_thi_h2s` (a derived rate); ship rule `scripts/generators/generate_kinetic_core_b50_ship_rule.py`.*

## 1. Why

The sulfur lane makes H2S only from cysteine (and, since 2026-10-27, glutathione), so a pentose + thiamine
pot cannot form 2-furfurylthiol. Dwivedi & Arnold 1971 (J. Agric. Food Chem. 19:923; `arnold1971_extraction.md`)
measured H2S from heated thiamine: about 0.39 (0.34–0.44) mol H2S per mol thiamine degraded at pH 7.0,
~0.04 at pH 6.5, essentially none at pH ≤ 6.0 (1.2e-4 M thiamine HCl, 0.1 M phosphate, ~100 °C, 30 min).
The patent literature the backlog cites (Shi, US 2025/0221436) finds a furfurylthiol disulfide in a
thiamine + xylose flavour.

## 2. The step

`r_thi_h2s`: THI → H2S + C12/N4 fragments, at a DERIVED rate share(pH)/(1 − share) × (k_thi_hmp + k_thi_mesh),
so that the measured share of thiamine's total thermal loss goes to sulfide. share(pH): zero below 6.25,
log-linear through (6.5, 0.04) and (7.0, 0.39), held at 0.39 above 7.0. **Nothing is fitted.** One
temperature: the share is held constant with temperature (declared). Every thiamine pot in the sulfur fit
runs at pH ≤ 5, where the share is zero, so B9 is reproduced exactly and no refit is needed.

## 3. Ship rule

- **T1:** B9's cost on its 54 rows unchanged to 1e-9 relative.
- **T2:** union panel within 3x not lower with the step than without; median |log10 error| on the sulfur
  rows not worse by more than 0.05 dex.
- **T3:** directional panel, independent claims agreeing: not fewer with the step than without.
- **T4 (recorded):** ribose + thiamine (50 / 10 mM, 150 °C, 10 min) FFT at pH 5.0 / 6.5 / 7.0, with and without.

SHIP if T1, T2 and T3 hold.

## 4. Verdict

*(Written after the run, below this line.)*

**Amendment to T1, made after its first evaluation and stated here because it was.** As written, T1
compared B9's cost to 18.74024636, the value recorded in the B49 prereg. Glutathione (a new species,
added the same afternoon before B50 ran) moved B9's cost to 18.74023152 through the solver state, so the
literal T1 failed on a stale reference while B50's branch changed nothing. T1's purpose is "the branch
leaves every fit pot unchanged"; it is now evaluated as B9's cost with the branch off against on, in the
same code: **18.74023152427709 both ways, identical.**

**SHIP** (`kinetic_core_b50_ship_rule.md`). T1 as amended: identical. T2: union 17/50 within 3x both ways,
sulfur median |log10 error| 1.212 both ways (no scored pot runs thiamine above pH 6.25). T3: directional
17 both ways (on the panel as corrected the same day). T4: ribose 50 + thiamine 10 mM, 150 °C, 10 min:
FFT 0 / 1.22 / 10.0 µg/L at pH 5.0 / 6.5 / 7.0 with the branch (0 at every pH without); MFT 52.7 / 44.3 /
37.0. The branch is untested against any scored measurement; its size carries two unmeasured transfers
(temperature, and 0.12 mM to flavour-pot concentrations).
