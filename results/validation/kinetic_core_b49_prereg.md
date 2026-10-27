# Wave B49 pre-registration — norfuraneol's other sulfide products, and the pH they answer to

*Written 2026-10-27 before any B49 number existed. Generator:
`scripts/generators/generate_kinetic_core_b49_fit.py`; ship rule:
`scripts/generators/generate_kinetic_core_b49_ship_rule.py`.*

## 1. Why

Two failures share one cause.

- **B47** (§4a and §6): at B9's optimum the model turns 85–88 % of fed norfuraneol into MFT and then
  destroys ~99 % of it through the thiol sink (`ch_mft_decay` ≈ 1 /min at 140 °C, in a pot with nothing
  for a thiol to react with). The measured yields are 0.05–0.15 mol %. Correcting Whitfield's pots to
  their real 50 mM exposed it, and the fit could only buy the row back by opening the hexose route.
- **B30**: on the only pot measured at two pH values by one laboratory (Whitfield & Mottram 1999 at
  pH 4.5, 2001 at pH 6.5; fed norfuraneol + cysteine or H2S, 140 °C, 60 min), the model makes ~2x MORE
  thiol at pH 6.5 where the measurement drops by ≥150x. B30 found the lane's pH mechanism absent on
  this route and left the cause as "the substrate or the sulfide budget".

Whitfield 2001 Table 2 settles B30's question: with **H2S fed directly** the thiols still collapse at
pH 6.5 (44 → 1 µg per 10 mg norfuraneol) while total volatiles hardly move (616 → 372), shifting to
thiophenes (21 → 128) and thiophenones (39 → 129). So it is neither the sulfide supply nor the
substrate: norfuraneol + sulfide has other products, they dominate even at pH 4.5 (503 of 616 µg; MFT
12), and their share grows with pH. The model has no such branch, so it can only make norfuraneol into
thiols, then needs a fast sink to get the yield down, and has nothing that can fall with pH.

## 2. The structure

One reaction, `r_nf_shet`: NF + H2S → SHET (a lumped C5, one-sulfur pool of dithiolanones,
thiophenones, thiophenes and trithiolanes; terminal). It proceeds through the **hydrosulfide anion**
(the lane's existing `hs_anion` factor), so its rise with pH is the declared mechanism, not a fitted
slope: from pH 4.5 to 6.5 the HS⁻ fraction rises roughly a hundredfold, while the existing
`r_nf_mft` / `r_nf_mp3p` run through neutral H2S. Its barrier is the lumped formation barrier the fit
already carries. **One new coordinate**, `log10 k_nf_shet` (band −8 … 1). The step is inert (zero)
unless a B49 report supplies it. The new species changes the ODE state, which moves B9's cost in the
seventh figure (18.7402487 → 18.7402464); that is numerical, and the live artifacts are regenerated.

## 3. Objective

B47's 52 corrected rows plus **two new FIT rows at pH 4.5** (Whitfield & Mottram 1999 Table 2, via
Whitfield 2001's reprint of it), as molar ratios SHET / MFT in the existing systems:

| row | system | classes (µg / 10 mg NF) | molar ratio | sigma (log10) |
| --- | --- | --- | --- | --- |
| `whitfield_h2s_shet_over_MFT` | `whitfield_nf_h2s` | dithiolanones/dithianones 443, thiophenes 21, thiophenones 39; MFT 12 | 34.9 | 0.4 |
| `whitfield_cys_shet_over_MFT` | `whitfield_nf_cys` | 41, 94, 94; MFT 15 | 14.9 | 0.4 |

Class masses become moles with declared mean molar masses (dithiolanones/dithianones 141, thiophenes
112, thiophenones 114 g/mol); sigma 0.4 covers that, the class assignment and the semi-quantitative
GC. Free set: B47's 23 + 1. Starts: B47's optimum with `k_nf_shet` at its band centre; the perturbation
protocol. Same bands, weighting, optimiser, budget.

**Held out, not fitted: everything at pH 6.5.**

## 4. Ship rule

- **T1, the mechanism (hold-out):** Whitfield 2001's pot, NF + cysteine 50/50 mM, 0.5 M phosphate
  pH 6.5, 140 °C, 60 min. (a) predicted MFT ≤ 0.003 mol % (the measured non-detect is < 0.001; 3x for
  the detection limit's own uncertainty); (b) predicted MFT falls by ≥ 10x from the pH 4.5 twin to the
  pH 6.5 pot (measured ≥ 150x). Under B9 it rises ~2x.
- **T2:** hold-out sulfur rows within 3x not lower than B9's by more than one row; median |log10
  error| not worse by more than 0.1 dex (as B47).
- **T3:** the Hofmann ribose pH-5 rows: within 3x not lower than B9's; median not worse by more than
  0.1 dex.
- **T4:** directional panel, independent claims agreeing: not fewer than B9's 25.
- **T5 (recorded):** Laplace identified count; whether `k_nf_shet` is identified; the MFT sink constant
  at 140 °C against B9's; the Whitfield 50/97 mM H2S row's residual.

**SHIP if T1, T2, T3 and T4 hold.** Probability, stated now: T1 ~0.6 (the anion fraction gives the
right sign and roughly the right size, but the fit may keep the sink fast), the whole rule ~0.35.

## 5. Verdict

*(Written after the run, below this line.)*

**DO NOT SHIP** (`kinetic_core_b49_ship_rule.md`, 2026-10-27). The engine keeps B9; `r_nf_shet` stays in
the network at zero.

| test | B9 | B49 | pass |
| --- | --- | --- | --- |
| T1 Whitfield pot, NF + Cys 50/50 mM, 140 °C: MFT fall pH 4.5 → 6.5 (measured ≥ 150x) | 0.48x (rises) | 1.01x | **no** |
| T2 hold-out sulfur rows within 3x | 1/14 (1.87 dex) | 3/14 (1.68 dex) | yes |
| T3 Hofmann ribose pH 5 within 3x | 0/2 (0.70 dex) | 1/2 (0.36 dex) | yes |
| T4 directional, independent agree | 25/43 | 25/44 | yes |
| T5 Laplace (recorded) | 20/23 | 17/24; `k_nf_shet` identified (σ 0.39 dex) | — |
| union within 3x (recorded) | 12/45 | 15/45 | — |

Both starts at cost 28.99; log10 k_nf_shet −3.45; the excess-H2S row −2.0σ → −1.6σ without the hexose
route opening (`k_glc_ha` +1.0 decade, against B47's +3.6).

**Why T1 failed, probed after the verdict.** The branch does what §2 declared: SHET/MFT rises 85x (Cys pot)
and 44x (H2S pot) from pH 4.5 to 6.5. MFT nevertheless stays flat in the cysteine pot because the model's
pot is H2S-limited (77 % of the norfuraneol unreacted after 60 min) and cysteine releases H2S faster at
pH 6.5 (Zheng & Ho's measured pH dependence), which cancels the branch. In the H2S-fed pot MFT does fall
(2.9x). The fit also moved `k_thiol_decay` +6.65 decades (unidentified, σ 7 dex), so MFT levels still
follow the timing of formation against a fast sink: early-made MFT in the H2S pot is destroyed, the
cysteine pot keeps making it — which is why the two new ratio rows pull opposite ways (+2.0σ, −2.45σ;
model SHET/MFT 225 vs 1.55 at pH 4.5 where the paper has 35 vs 15). **The heterocycle branch is the right
chemistry and the wrong lever:** until the thiol sink is measured (EXPERIMENTS_BRIEF row 1, "where the
thiols go"), the sulfur lane's levels and their pH response are set by a constant no source pins.
