# Why laboratories disagree: what the method papers say (STATUS_AND_PLAN sec. 0, step 5b, part d)

*Written 2026-11-12. Literature only; no engine output. The numbers behind each statement are in the two reading
notes beside this file ([thiols](disagreement_reading/method_read_thiols.md),
[inter-laboratory](disagreement_reading/method_read_interlab.md)), each written by a separate reader who saw the
seven conflict cells, the 12 method papers and the methods sections of the conflict cells' own sources. The
generated audit is [lever_disagreement_audit.md](../../../results/validation/lever_disagreement_audit.md).*

## The short answer

**Measurement noise cannot flip a sign; how the lever was applied can.** A laboratory repeats itself to 0.01 - 0.03
dex, and laboratories measuring the same samples under one protocol differ by 0.10 - 0.24 dex in absolute level
(two olive-oil ring trials, Aparicio-Ruiz 2022 and Casadei 2021). A within-study ratio cancels a constant bias, so
only a bias that moves WITH the lever survives into a contrast; the documented ones are 0.1 - 0.45 dex (thiol
binding to browned matrix or protein, the matrix's physical state, SPME fibre competition at high volatile load).
The sign flips of the seven conflict cells are larger than that, and the papers point instead to **differences in
what the lever meant** in each laboratory: an initial pH that the reaction did not keep, a "none" arm that held a
different base mixture, an effect too small to separate from zero in either laboratory. A five-laboratory
comparison of identical samples (Eckert 2023) shows the same pattern at its extreme: one treatment changed sign
between laboratories by up to about 2 dex.

## The seven conflict cells

| cell | the two laboratories | best-supported explanation | confidence |
|---|---|---|---|
| MFT x pH (6 - 7.5) | xu2007 (CAU; 0.2 M phosphate) -0.62 against zhou2023 (Ho network; unbuffered) +0.36 | **The pH step is not the same lever.** zhou2023 adjusted the initial pH with NaOH and did not buffer: initial pH 6 and 7 both ended near pH 3.3 - 3.4 (its Fig. 2), so its "6 to 7" step barely moved the reaction pH. Its 7 to 8 step, where the final pH did move (about 3.4 to 5.0), gives MFT -0.48 dex: the same sign as the buffered xu2007. Secondary: xu2007 sampled the browned, undiluted mixture by 30 min headspace at 60 C, and browned matrix binds thiols (Hofmann 2001: 0.3 - 1 dex lost in 20 - 30 min), which would steepen its fall with pH | medium - high |
| bis(2-methyl-3-furyl) disulfide x pH (6 - 7.5) | the same two | The same pH-basis mismatch (zhou2023's 7 to 8 step: -0.31 dex, the sign of xu2007). Inside xu2007, the control arm falls and the less browned supercritical-CO2 arm rises, so browning-linked loss also matters | medium |
| hexanal x sugar identity (xylose vs arabinose, fructose, glucose) | chen2023lycium (Hefei; defatted goji seed meal) + against zhang2024sugars (NUS; chicken carcass hydrolysate) - | **Probably real chemistry, not measurement.** In the fatty chicken hydrolysate hexanal is already present before heating and falls with pentoses; the defatted seed meal starts with little. Measurement adds noise, not the sign: chen2023lycium's own GC-MS and GC-IMS rows of the same pots disagree on the effect by 0.12 - 0.39 dex, and zhang2024sugars quantifies by raw FID area with no internal standard; fibre competition is unlikely (the xylose pot has the lower total area) | low - medium |
| benzaldehyde x sugar identity (xylose vs glucose) | chen2023lycium against theng2024 (A*STAR; nutritional yeast) | **Not a real conflict.** Both laboratories' own tests find no difference between the two sugars; in theng2024 benzaldehyde is present before heating. chen2023lycium prints nearly the same retention index for hexanal and benzaldehyde (872, 873), where standard values on that column differ by about 160 units, so at least one identification is in doubt (the reader's comparison, not the paper's) | medium - high |
| trimethylpyrazine x amino acid identity (none vs arginine) | hwang1995 (Ho network) against park2022mealworm (Wonkwang) | **The "none" arm is not the same pot.** In hwang1995 "none" still holds glucose + glycine; in park2022mealworm it is protein + xylose with no free amino acid, and every added amino acid raises pyrazines over that base. park2022mealworm also added basic arginine with no reported pH correction, and extracted by 2 h of simultaneous distillation-extraction with sugar and amino acid still present, which runs the Maillard reaction during the extraction (Chaintreau 2001) | medium - high |

**Added after this reading (same day).** Keying seven more papers raised the audit to 9 conflict cells. The two new
ones, 2-methylpyrazine and nonanal against cysteine dose (ho_network and nus_liu_shaoquan against ding2023, a fish
protein hydrolysate), were not read against the method papers. Cysteine is known to divert carbonyls from pyrazines
and lipid aldehydes in some matrices and not others; that is a hypothesis, not a reading.

## What this changes

1. **pH cells should not pool buffered and unbuffered initial-pH studies.** In an unbuffered pot the initial pH
   does not describe the reaction. PREREG sec. 4 already separates initial, held and final pH, but an unbuffered
   initial-pH study sits in the same sub-cell as a buffered one. A rule that keys the final pH where the paper
   prints it, or bands initial-pH contrasts by buffered against unbuffered, would remove the two thiol conflicts.
   Because this is seen after the data, it goes to the backlog as a proposed amendment, scored beside the current
   table, not applied silently.
2. **A categorical cell's reference level must be the same pot.** "None vs arginine" compares different base
   mixtures across studies. The categorical network (PREREG sec. 4) assumes the reference arm means the same thing;
   a per-study note of the base mixture would let the generator separate them.
3. **The "what to record" table of the experiments brief is right, and two rows matter most:** the pH before and
   after heating, and the base mixture of the reference arm.
4. **Method, not chemistry, sets the floor on absolute levels:** about 0.1 - 0.25 dex between laboratories even
   under a harmonised protocol, more without one. For ratios the method cancels unless recovery moves with the lever.

## Limits

- Seven cells, each with two laboratories. The explanations are hypotheses that fit the papers; none is tested.
- No inter-laboratory study exists for thiols, Strecker aldehydes or alkylpyrazines; the measurement figures come
  from olive-oil volatiles, cocoa and plant volatiles.
- Reading notes were produced from text layers for most PDFs; the notes say which were read as page images.
