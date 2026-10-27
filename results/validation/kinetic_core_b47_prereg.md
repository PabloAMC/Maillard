# Wave B47 pre-registration — the sulfur objective, corrected against its PDFs

*Written 2026-10-27 before any B47 number existed. Generator:
`scripts/generators/generate_kinetic_core_b47_fit.py`; ship rule:
`scripts/generators/generate_kinetic_core_b47_ship_rule.py`.*

## 1. Why

On 2026-10-27 the 54 rows of B9's objective (the engine's sulfur fit since 2026-09-03) were checked,
row by row, against the source PDFs for the first time (two read-only passes; findings in
`tasks/data_restructure_plan.md`, "Sulfur fit rows checked against the PDFs"). 29 rows reproduce
exactly or by a correct conversion, 15 are figure reads within their sigma, and **six rows and two
systems do not say what their papers say**. B9 was fitted to them. B47 is B9 with those rows
corrected and nothing else changed.

## 2. What changes

| row / system | B9 | B47 | evidence |
| --- | --- | --- | --- |
| `whitfield_mercaptoketone_over_MFT` | target 16.3 | **4.99** | Whitfield 1999 Table 1: 16.3 = (92 + 74.5 + 77.5)/15, three compounds by mass. The row's species is 2-mercapto-3-pentanone alone: (77.5/118.20)/(15/114.16) = 4.99 molar |
| `cerny_isomer_split` | MP3P/MP2P = 1.0 | **removed** | Cerny 2007 Table 2, pH 5.0: 2-mercapto-3-pentanone 0, 3-mercapto-2-pentanone 141 (formed "only at pH 6 and 7"). The 1:1 came from Whitfield's different pot; a log residual cannot hold a zero |
| `yaghmur_fft_share_ceiling` | ceiling on the 145 °C ribose pot | **removed** | the "< 1 %" is furfural 0.3 M + cysteine 0.6 M at 65 °C / 15 h; it does not transfer to a 145 °C / 20 min ribose pot |
| `kang_free_cys_yield_ceiling` | 16.3 % | **12.1 %** | Kang 2026 SI Fig. S3: 16.3 % is the 140 °C curve's peak; the row's own 120 °C system peaks at ~12.1 % (figure read) |
| system `whitfield_nf_cys` | NF 20, Cys 20 mM | **NF 50, Cys 50 mM** | Whitfield 1999 Methods: 1 mL 0.1 M + 1 mL 0.1 M |
| system `whitfield_nf_h2s` | NF 20, H2S 40 mM | **NF 50, H2S 97 mM** | same, with ~0.2 M H2S solution |
| `whitfield_nf_cys_MFT`, `whitfield_nf_h2s_MFT` | mol % basis 20 | **basis 50** | follows the systems; the printed mol % values are unchanged |
| `fed_ribose_h2s_MFT` | 0.01 mol % | **0.0132** | Hofmann 1998: 15.1 µg; the printed % is rounded to one figure |
| `fed_c2c3_MFT_pH3` | 0.01 | **0.0136** | 15.5 µg |
| `fed_thiamine_MFT` | 0.01 | **0.0072** | 8.2 µg |
| `fed_tdp_h2s_FFT` | 0.08 | **0.069** | 78.6 µg; the paper's own % is 16 % higher than its µg value |
| `fed_nf_cys_MFT` | 0.05 | **0.0445** | 50.8 µg; the paper's own % is 12 % higher |
| buffers of the Whitfield, van Seeventer and Cerny 2007b systems | `BUFFER_PHOSPHATE_ASSUMED` | **declared by source**, numerically identical (0.5 M phosphate) | each paper prints 0.5 M phosphate at the matching pH |

Rule for the Hofmann rows, stated once: where a paper prints both a mass and a rounded percentage, the
row uses the mass, converted. Objective: **52 rows** (54 − 2).

## 3. What does not change

Network, species, T-structure, the 23-coordinate free set, the declared bands (`full_bounds`), the
W-HALF pH weighting, every sigma, the optimiser, its 600-evaluation budget and the two-start protocol
(start 0 = B9's frozen optimum, start 1 = the perturbation protocol). The B2.3, B8 and B9 generators
stay frozen; B47 imports B9 and patches.

## 4. Expectations, stated in advance

- Most coordinates move little: 46 of 52 rows are unchanged.
- The mercaptoketone row falls 3.3x and the Hofmann thiamine row falls 1.4x, so the routes feeding
  2-mercapto-3-pentanone and the thiamine → MFT route should come down; the Kang ceiling tightens by 1.3x.
- The two removed rows were both constraints the model could satisfy by being low, so removing them
  should loosen rather than tighten.
- **Probability B47 ships under §5: about 0.75.** The corrections are small against sigmas of 0.3–0.6
  dex, but the B9 optimum sits on four bounds and a bound-limited optimum can move by more than its
  rows suggest.

### 4a. Found while checking the harness, before either start ran

B9's vector scored on B47's 52 rows gives cost 36.46 (B9's own 54 rows: 18.74, reproduced). Almost all
of the rise is one row: `whitfield_nf_h2s_MFT` moves from −2.0σ to −6.3σ. At B9's optimum the model's
MFT yield from norfuraneol + H2S falls ~150x when the pot goes from 20/40 mM to the paper's 50/97 mM
(0.012 → 8e-5 mol %; quick and careful integration agree to five digits), and it is non-monotone in
H2S (50/50 mM: 0.004 mol %). So **excess H2S destroys MFT in the model**, and B9's wrong concentrations
hid it; B9's note "two labs agree within ~1.6x on the H2S channel" rested on them. Expectation revised:
the optimiser will try to buy this row back by moving the H2S-consuming constants, and T2/T3 decide
whether that costs the out-of-sample rows. The ship probability above drops to **about 0.55**.

## 5. Ship rule (evaluated by `generate_kinetic_core_b47_ship_rule.py`, B9 vs B47 by parameter override)

A correction ships unless it makes the out-of-sample evidence worse.

- **T1 (recorded):** B47 cost on its 52 rows, and B9's vector scored on the same 52 rows.
- **T2:** the Maillard-path hold-out panel's sulfur rows: within-3x count under B47 is not lower than
  under B9 by more than one row, and the median |log10 error| does not grow by more than 0.1 dex.
- **T3:** the four Hofmann pH-5 out-of-sample bundles (ribose, xylose, glucose, fructose + cysteine;
  validation rows since B9): within-3x count not lower than B9's, median |log10 error| not worse by
  more than 0.1 dex.
- **T4:** the directional panel on the core: independent claims agreeing not fewer than B9's.
- **T5 (recorded):** Laplace at the B47 optimum, identified coordinates against B9's 20 of 23.

**SHIP if T2, T3 and T4 hold.** Then the engine reads `kinetic_core_b47_fit_report.json` and its
Laplace covariance, and every live artifact is regenerated after the fit is frozen. **DO NOT SHIP** keeps
B9; the record is kept and the report says which test failed, because a corrected objective that
scores worse out of sample means B9's defects were compensating for a structural error somewhere else.

## 6. Verdict

*(Written after the run, below this line.)*
