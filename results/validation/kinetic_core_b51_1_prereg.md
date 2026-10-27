# Wave B51.1 pre-registration — the two pH-state corrections that do not touch the thiolate

*Written 2026-10-27 after B51's verdict and before any B51.1 score existed. Ship rule:
`scripts/generators/generate_kinetic_core_b51_1_ship_rule.py`. B51's own records:
`kinetic_core_b51_prereg.md` §6 and `kinetic_core_b51_ship_rule.md`.*

## 1. Why

B51 bundled six corrections and did not ship: refitted, the model lost the Hofmann pH 5 pots ~100x
because the corrections that raise the in-situ thiolate fraction (GKL cysteine, ionisation heat
capacities) expose the thiolate-acting thiol sink. B51's per-switch attribution separated the six. Two of
them do not change any speciation factor at all:

- **`label_scale_tables`**: the three MEASURED pH-indexed tables (Zheng & Ho's cysteine and glutathione
  thermolysis, Dwivedi & Arnold's thiamine H2S share) were measured on 25 °C labels and were being read
  at the in-situ pH. Under the switch they read the cooled pH of the current state, the scale their
  sources printed. Speciation factors are untouched.
- **`sid_across_segments`**: on a multi-segment heating program the strong-ion difference was re-derived
  from the label at every segment, re-titrating away the acid made earlier. A bug; inert on every
  single-segment pot, which is every fit system and every panel pot.

## 2. Disclosure

B51's attribution, computed before this file, already scored B9's frozen vector with each switch alone on
B9's own 54 rows: `label_scale_tables` 18.74 → 15.47, `sid_across_segments` 18.74 → 18.74. That is a
FIT-side number and it is the reason this wave exists; it is not evidence about the hold-outs, which no
B51 arm scored with these two switches alone.

## 3. Candidate and rule

- **Candidate.** B9's frozen vector (no refit) with `label_scale_tables` and `sid_across_segments` on, the
  other four B51 switches off. **Baseline.** B9, all off (what ships).
- **T1 (recorded):** B9's cost on its 54 rows, candidate vs baseline.
- **T2:** Maillard-path hold-out sulfur rows: within-3x not lower than the baseline's by more than one row,
  median |log10 error| not worse by more than 0.1 dex.
- **T3:** Hofmann ribose pH 5 (trust loop) + xylose pH 5 (hold-out): within-3x not lower, median
  |log10 error| not worse by more than 0.1 dex.
- **T4:** directional panel, independent claims agreeing: not fewer.
- **SHIP if T2, T3 and T4 hold**: `ph_state.B51` becomes these two switches on. DO NOT SHIP leaves it off.

## 4. Expectation

Every buffered pot reads nearly the same pH on both scales (0.5 M phosphate: 5.00 cooled vs 4.94 in situ
at 145 °C), so the hold-out panel, which is buffered, should barely move; the effect is in unbuffered pots
(an unbuffered pH 7 label is ~5.4–5.7 in situ). Probability of SHIP: **about 0.7**; the risk is a
directional claim on an unbuffered pot flipping.

## 5. Verdict

*(Written after the run, below this line.)*

**DO NOT SHIP** (`kinetic_core_b51_1_ship_rule.md`, 2026-10-27). `ph_state.B51` stays all off.

T1 (recorded) 18.74 → 15.47; T2 PASS (1/14, unchanged; both Hofmann ribose pH 7 rows improve 134x → 110x
and 140x → 115x); T3 PASS (3/4, unchanged); **T4 FAIL, 18/31 → 17/31**: SUG-09 ("FFT from pentoses, NOT
DETECTED from hexoses") flips — pentose 36.9 vs hexose 33.5 µg/L (agree) becomes 39.5 vs 49.9 (disagree).
The claim agreed before only by ordering: the model already makes 33 µg/L of FFT from a hexose that the
paper does not detect at all, the hexose-route defect B47 also ran into. On the unbuffered pots the
cooled pH is higher than the in-situ pH, cysteine thermolysis read at it is faster, and the extra H2S
reaches FFT through the hexose route. The correction is right; the defect it moves is the hexose entry
(`r_glc_c2c3` / `r_glc_fur`, unidentified since B9). Both switches stay in the code, off, until that entry
is fixed; the SID carry (a bug fix, inert on every single-segment pot) is held back with it because the
rule was written for the pair.
