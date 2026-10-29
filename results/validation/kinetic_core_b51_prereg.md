# Wave B51 pre-registration — the pH state, corrected against its thermodynamic sources

*Written 2026-10-27 before any B51 fit, score or cost existed. The only B51 numbers computed before this
file was written are equilibrium constants (pKa at 145 °C, below). Generators:
`scripts/generators/generate_kinetic_core_b51_fit.py` (fit), `scripts/generators/generate_kinetic_core_b51_ship_rule.py`
(ship rule). Switches: `src/kinetic_core/ph_state.py` `PhStateCorrections` (`ph_state.B51`); all off reproduces B50.*

## 1. Why

The sulfur lane's pH is a computed state (B2.2): a conserved strong-ion difference (SID) from the 25 °C
label, re-solved at reaction temperature, feeding speciation factors on every pH-sensitive step. Checking
its constants against their cited sources on 2026-10-27 found errors in the constants themselves, two
defects in how the state is carried, and three data records that do not say what their papers say.

## 2. What changes — six switches, each attributable

| switch | B50 (now) | B51 | source |
| --- | --- | --- | --- |
| `heat_capacity` | every pK(T) by constant-dH van 't Hoff | + the ionisation heat capacity dCp | Goldberg, Kishore & Lennen 2002, J. Phys. Chem. Ref. Data 31:231 (GKL) Tables 7.2, 7.22, 7.23, 7.31, 7.51, 7.52, 7.56. dCp < 0 everywhere, so pK RISES relative to constant-dH: phosphate pK2 at 145 °C 7.02 → 7.63 |
| `cysteine_gkl` | cysteine thiol 8.33 / 26 kJ/mol, ammonium 10.28 / 43 kJ/mol, cited "[GKL]" | GKL Table 7.23: 8.36 / 36.1 / −66, 10.75 / 34.1 / −204 (pK / dH kJ/mol / dCp J/K/mol) | the B2 pair is NOT in GKL; at 145 °C the thiol moves 7.02 → 6.72 and the ammonium 8.12 → 9.58 |
| `h2s_suleimenov` | H2S pK1 7.05, dH 22 kJ/mol → 5.94 at 145 °C | Suleimenov & Seward 1997 (GCA 61:5187) measured correlation as printed by Ning et al. 2014 (Corrosion 70:375) Table 1 → 6.49 at 145 °C | constant-dH misses the pK1 minimum near 100–130 °C; HS⁻ fraction at pH 5 falls ~3.5x |
| `ionic_strength` | activities = concentrations | Davies, A(T) from Malmberg & Maryott's dielectric constant; mixed constants in the charge balance and the speciation numerators (normaliser at I = 0, so k_ref means infinite dilution) | the buffered pots run at I ≈ 0.5–1.5 M; Davies is a declared approximation above ~0.5 M |
| `label_scale_tables` | Zheng & Ho cysteine and glutathione thermolysis, Dwivedi & Arnold thiamine H2S share read the IN-SITU pH | they read the COOLED (25 °C) pH of the current state — the scale their sources printed | the mismatch disclosed at B50 |
| `sid_across_segments` | SID re-derived from the label at every thermal segment | computed once per pot and carried | on a multi-segment program the drift reset at each boundary (a bug) |

pKa at 145 °C, all off → all on (computed before this file, equilibrium constants only): H2S 5.94 → 6.49,
cysteine ammonium 8.12 → 9.58, thiol 7.02 → 6.72, phosphate pK2 7.02 → 7.63.

**Data corrections, applied to BOTH arms (they are facts, not candidates):**

- the nine trust-loop Hofmann 1998 bundles receive the paper's 0.5 M phosphate block (Amendment 44); they
  were scored unbuffered while every sulfur fit integrated them in 0.5 M phosphate;
- the pH claims PH-01/02 (Meynier & Mottram 1995) and MOT02-PH-01/02 are charged as 0.2 M **pyrophosphate**
  (were orthophosphate), and MOT-01/02/03 (Mottram & Nobrega 2002), which carried no buffer, get theirs;
- `BufferSpec` gains `pyrophosphate` and `phthalate` (GKL Tables 7.56, 7.52); `comparative_cli` reads every
  molarity BufferSpec knows (it silently dropped acetate and citrate).

## 3. Candidates, in a fixed order of preference

- **(A) B51 refit.** B47's 52-row objective (the rows corrected against their PDFs; B47's generator
  installs them) with all six switches on; same network, free set, bands, sigmas, W-HALF weighting,
  optimiser, budget and two-start protocol (start 0 = B9's frozen optimum, start 1 = perturbation).
- **(B) B9 under B51 physics, no refit.** B9's frozen vector with all six switches on.
- **Baseline.** B9's frozen vector, switches off, data corrections applied (what ships today plus facts).

Attribution (recorded only): B9's vector scored on B9's 54 rows with each switch alone.

## 4. Expectations, stated in advance

- The HS⁻ branches lose ~3.5x of their sulfide fraction at pH 5 and the thiolate steps gain ~2x; the fit
  will move `k_*_hs` up and the thiolate sinks down to compensate. Under (B) nothing compensates, so (B)
  will score worse on the fit rows than the baseline. Probability (B) passes §5: **about 0.3**.
- The excess-H2S defect B47 exposed (MFT destroyed when H2S is in excess) runs through the HS⁻ branches;
  a smaller HS⁻ fraction should weaken it, so (A) should need less of B47's hexose-route opening.
  Probability (A) passes §5: **about 0.45**. The main risk is B47's: SUG-09 and the hexose claims.
- The pH 7 hold-out rows (Hofmann ribose pH 7) should move most: phosphate's pK2 rises, so a pH-7-labelled
  0.5 M phosphate pot is more alkaline in situ than B50 thinks.

## 5. Ship rule (`generate_kinetic_core_b51_ship_rule.py`)

Each arm is scored in its own child process (switches and fit report patched), against the baseline.

- **T1 (recorded):** cost on B47's 52 rows for every arm; per-switch attribution on B9's 54 rows.
- **T2:** Maillard-path hold-out sulfur rows: within-3x count not lower than the baseline's by more than
  one row, and median |log10 error| not worse by more than 0.1 dex.
- **T3:** the trust-loop Hofmann ribose-cysteine pH 5 rows plus the xylose pH 5 hold-out: within-3x not
  lower, median |log10 error| not worse by more than 0.1 dex.
- **T4:** directional panel, independent claims agreeing: not fewer than the baseline's.
- **T5 (recorded):** the pH-axis subset of the directional panel; Laplace identified coordinates for (A).

**Verdict:** SHIP (A) if (A) passes T2, T3, T4; else SHIP (B) if (B) passes them; else DO NOT SHIP (switches
stay off; the data corrections stay, since they are facts, and the baseline becomes what ships).
If (A) ships the engine reads `kinetic_core_b51_fit_report.json` and its Laplace covariance; every live
artifact is regenerated afterwards.

## 6. Verdict

*(Written after the run, below this line.)*

**DO NOT SHIP** (`kinetic_core_b51_ship_rule.md`, 2026-10-27). The six switches stay OFF; the engine keeps
reading B9. The data corrections of §2 stay (they are facts), so the baseline column is what ships.

| | baseline (B9, off) | (A) refit | (B) B9, on |
| --- | --- | --- | --- |
| T1 cost, B47's 52 rows (recorded) | 36.46 | **20.41** | 43.29 |
| T2 hold-out sulfur within 3x | 1/14, 1.87 dex | 1/14, 1.61 dex — PASS | 0/14, 1.42 dex — PASS |
| T3 Hofmann ribose + xylose pH 5 | 3/4, 0.37 dex | **0/4, 1.80 dex — FAIL** | 0/4, 2.48 dex — FAIL |
| T4 directional, independent agree | 18/31 | 18/32 — PASS | 17/31 — FAIL |
| T5 pH-axis claims (recorded) | 6/8 | 5/8 | 5/8 |
| Laplace identified (recorded) | 20/23 (B9) | 15/23 | — |

**What happened.** (A) fits its own rows much better (20.41 against B47's 22.90 under the old physics)
and moves every Yiltirak hold-out row toward the measurement (MFT at 120 °C 56x → 5.1x; at 130 °C
100x → 7.4x). But it predicts the Hofmann ribose and xylose pH 5 pots ~100x too LOW (ribose FFT 2.6
against 121 µg/L measured; MFT 1.7 against 198), while pH 3 stays too high (FFT 2530 against 229): the
model's FFT falls ~1000x from pH 3 to pH 5 where Hofmann measured ~2x. Those pH 5 level rows have been
validation rows since B9, so the fit never saw them.

**Why, in the model's own terms.** The per-switch attribution on B9's vector separates the corrections:
`label_scale_tables` alone LOWERS B9's cost (18.74 → 15.47), `sid_across_segments` and `h2s_suleimenov`
are neutral (18.74, 18.88), `ionic_strength` costs little (19.64), and the two that raise the in-situ
THIOLATE fraction in the fit pots — `cysteine_gkl` (thiol pK 7.02 → 6.72 at 145 °C) and `heat_capacity`
(phosphate pK2 up, so the buffered pots' in-situ pH rises, 4.94 → 5.27 for a pH 5 label) — cost ~7 each.
The unidentified thiol sink is declared to act through the thiolate (`ch_thiolate_loss_*`, kind
"thiolate", using cysteine's thiol pK for MFT and FFT because neither is measured). Under correct
speciation that sink destroys the thiols above pH ~5, and the refit cannot hold both Yiltirak and
Hofmann pH 5. The physics corrections are right; the structure they expose is not. This is the same
constant B47 and B49 ended at.

**Follow-ups, each needing its own pre-registration (not adopted here):**
1. The two corrections that do not touch the thiolate structure — the 25 °C scale for the measured pH
   tables (improves B9's own fit) and the SID carried across segments (a bug fix, inert on single-segment
   pots) — as a separate wave.
2. The thiolate assignment of the thiol sink itself: with FFT's thiol pK nearer 9.4 (benzylic) the
   thiolate route would be ~100x weaker for FFT, and MFT's (heteroaryl, nearer 7) stronger. Neither pK is
   measured; the thiol-sink experiment in `docs/guides/EXPERIMENTS.md` would decide it.

## 7. Re-decision under the changed acceptance rule (2026-10-27, owner decision)

The verdict above was correct under the rule this file pre-registered, and it is kept. On the same day the
owner changed how changes are accepted (`scripts/generators/WAVES.md`, "How a change is accepted"):
corrections ship on verification of the change itself, and their scores are reported as diagnosis rather
than used as a veto. Under that rule **B51 SHIPS**, candidate (A):

- **Kind 1, the six switches, each verified independently of the model's scores:** GKL Table 7.23 (cysteine)
  and Tables 7.2 / 7.22 / 7.51 / 7.52 / 7.56 (dCp) checked against the PDF by the `goldberg2002` dossier's
  audit; Suleimenov & Seward's Eq. 22 checked against their Table 2 (within 0.023 everywhere); the dCp law
  checked against the tabulated pKw (0.04 / 0.11 units at 100 / 150 °C); all-off reproducing B50 bit for
  bit and the ionic-strength fixed point converged (`tests/unit/test_kinetic_core_b51.py`); the SID carry a
  bug fix; the 25 °C scale the one printed by the sources of the three tables.
- **Kind 2, the refit (A):** cost 20.41 on B47's 52 corrected rows, Laplace 15/23 identified, profile
  4 quadratic / 10 asymmetric / 1 flat / 8 bound-limited. (B), the old constants under the new physics, is
  not a candidate under the rule: constants calibrated on the old errors do not ship with the new physics.
- **Recorded worse, with the assumption each exposes:** the Hofmann ribose and xylose pH 5 pots (~100x low;
  T3 3/4 → 0/4) and the union panel (18/50 → 16/50) expose the thiol sink declared to act through the
  THIOLATE with cysteine's thiol pK standing in for MFT and FFT; SUG-09 exposes the hexose entry to the
  thiols. Recorded better: every Yiltirak hold-out row (MFT 120 °C 56x → 5.1x), hold-out sulfur median
  1.87 → 1.61 dex. These two assumptions are the next waves.
