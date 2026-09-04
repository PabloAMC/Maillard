# Hofmann & Schieberle 2002 — Chemical interactions between odor-active thiols and melanoidins involved in the aroma staling of coffee beverages

### Per-paper extraction, 2026-09-04, from the full text supplied by the owner (`data/articles/hofmann2001.pdf`; the article was published on the web 12/07/2001 and printed in the 2002 volume, hence the file name). **Nothing in `src/`, `tests/` or `results/` was changed by this extraction; it records what the code already cites and whether the paper supports it.**

**Provenance codes:** **[M]** measured and printed · **[F]** fitted by the authors · **[D]** derived by this extraction from a printed figure or table · **[NEG]** verified negative.

---

## §0. IDENTITY — CONFIRMED, DOI ON THE PAGE

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/hofmann2001.pdf` (8 pages) | read in full |
| title | "Chemical Interactions between Odor-Active Thiols and Melanoidins Involved in the Aroma Staling of Coffee Beverages" | p. 319 |
| authors | Thomas Hofmann and Peter Schieberle, Deutsche Forschungsanstalt für Lebensmittelchemie, Garching | p. 319 |
| DOI | `10.1021/jf010823n` | p. 319 footer |
| journal | *J. Agric. Food Chem.* 2002, 50, 319–326 | p. 319 header |
| dates | received June 25, 2001; revised October 11, 2001; accepted October 16, 2001; published on web 12/07/2001 | p. 326 |

## §1. WHAT THE CODE CITES, AND WHAT THE PAPER SAYS

| code site | claim in the code | what the paper prints | verdict |
|---|---|---|---|
| `src/kinetic_core/parameters_sulfur.py` (k_thioether, "Hofmann 2002 Table 2 9.8e-4 /s, 30 C, SIDA, CROSSPY model") | pseudo-first-order thiol sink 9.8e-4 /s at 30 °C in the N-acetyl-lysine/glycolaldehyde (CROSSPY) model | Table 2 **[M]**: 2-furfurylthiol stored 30 min at 30 °C with dry-heated N-acetyl-L-lysine/glycolaldehyde (10 mg each): **17 % (15–19) of the free thiol remains**; with albumin/glycolaldehyde 31 % (28–34); albumin/glucose 58 % (55–61); chlorogenic acid 92 % / 86 %. Method: 500 µg thiol in 10 mL 0.1 M phosphate pH 6.0, static headspace, SIDA. | **[D]** ln(100/17)/1800 s = **9.8e-4 /s** — reproduces the code's constant exactly. Note the paper's pH is **6.0** (the code records pH 5.6 for the channel). |
| same site ("Hofmann 2002 Fig. 6 9.4e-4 /s, 30 C, SIDA, 12.5 g/L coffee melanoidin") | 9.4e-4 /s from Figure 6 | Figure 6 **[M]** (thiol 400 µg + melanoidins at 30 °C): ≈240 µg at 15 min, ≈110 at 30, ≈85 at 60, ≈65 at 90 (read off the plot); with 1,4-diethyl-pyrazinium ions ≈190 / 50 / 12 / ≈5. Melanoidins 125 mg in 10 mL (12.5 g/L) per the Static Headspace Analysis section. | **[D]** first-order fits over the first 30 min give 5.7e-4 (0–15 min) to 7.2e-4 /s (0–30 min) for the melanoidin curve and 0.8–1.2e-3 /s for the pyrazinium curve; 9.4e-4 sits inside that spread, i.e. **supported to within the read-off** but not a printed number. **The melanoidin curve plateaus at ≈65–85 µg (≈80 % bound)**: a capacity-limited sink, which is the depletable-electrophile structure the module already encodes. |
| same site (note: "THE DISULFIDE BRANCH IS DEAD AT 30 C … <1.5 % of the thiol flux") | disulfide formation negligible | p. 325 **[M]**: "Although about 400 or 330 µg of 2-furfurylthiol was bound to the coffee melanoidins or the pyrazinium-derived intermediates, respectively, **less than 6 µg of the corresponding bis(2-furfuryl) disulfide was generated**." | **confirmed**: 6/400 = 1.5 %. |
| `scripts/generators/generate_kinetic_core_b2_holdout.py` row `hofmann2002_brew_80C_FFT` ("Hofmann 2002 Fig. 1: FFT loss in a real coffee brew (50 g powder/L, thermos, 80 C) at 0.023 /min, t1/2 ~30 min") | FFT decay 0.023 /min at 80 °C in brew | Figure 1 **[M]** (brew 50 g/L, thermos 80 °C, SIDA): 2-furfurylthiol **16.0 µg/kg at 0 min, ≈10 at 30, ≈4 at 60, ≈1.5 at 90, ≈0.3 at 210**; 3-mercapto-3-methylbutyl formate 8.2 → 6.7 → 5.0 → 3.5 → 1.3. Text p. 322: "After 60 min the 2-furfurylthiol concentration decreased by a factor of more than four". | **[D]** ln(16/4)/60 min = **0.023 /min** — the row's number, from the printed 0 and 60 min points. Brew basis 50 g/L confirmed (Isolation section). |
| B8 prereg H-3 ("hofmann2002_brew_80C_FFT … AT RISK") | the brew is slower than the 30 °C model systems | 80 °C brew: 3.9e-4 /s vs 30 °C melanoidin model 6–9e-4 /s at 12.5 g/L | **confirmed**: the sink in real brew at 80 °C is *not* faster than the 30 °C model, which is what the row exists to test. |

## §2. WHAT THE PAPER DOES NOT GIVE

- **[NEG]** No activation energy and no second temperature for the same matrix: the 30 °C model and the 80 °C brew differ in melanoidin concentration and matrix, so no Ea can be derived from them. The code's "no activation energy for this channel by policy" stands.
- **[NEG]** No MFT (2-methyl-3-furanthiol) kinetics: MFT appears only in the aroma dilution table (rFD 4 → 2 with melanoidins) and in the LC/MS binding experiments; the code's assignment of the FFT channel constant to MFT remains a declared assumption, as it says.
- Table 1 (rFD factors, 30 °C, 30 min, 125 mg melanoidins): FFT 32 → 2, 3-methyl-2-butenthiol 8 → 1, 3-mercapto-3-methylbutyl formate 8 → 2, methanethiol 2 → <1; non-thiol odorants unchanged. Qualitative support for thiol specificity.

## §3. VERDICT

Every number the code attributes to this paper is either printed (Table 2, the <6 µg disulfide) or a one-step derivation from printed points (Figures 1 and 6). The one discrepancy is the pH: the model reactions ran at pH 6.0 in 0.1 M phosphate, while the channel's `ph=5.6` follows Charles-Bernard's coffee conditions. No change to the code is warranted by this paper.
