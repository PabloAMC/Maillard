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

**The two cells added the same day** (cysteine dose x a non-thiol; [reading note](disagreement_reading/method_read_cysteine_dose.md)):

| cell | the two laboratories | best-supported explanation | confidence |
|---|---|---|---|
| 2-methylpyrazine x cysteine dose | ding2023 (ZJUT; fish hydrolysate + xylose + glycine) -0.43 against zhou2022 (Ho network; Amadori product, no free amine) +0.20, +0.35 | **Opposite sides of one rise-then-fall curve.** ding2023 rises from 0 to 0.5 % cysteine and falls above it (Table 2); the keyed slope is the falling half, where cysteine is 2.5 - 5 times the xylose (mol/mol) and takes the dicarbonyls that would make pyrazines. zhou2022 stays at or below one cysteine per sugar, in a system with no other nitrogen, where cysteine is the nitrogen donor (Table 3) | medium - high (mechanism: medium) |
| nonanal x cysteine dose | ding2023 +0.22 against li2022heat (NUS; pork-trimming hydrolysate) -0.13 to -0.44 | **Mostly different segments of non-monotone ladders.** ding2023 falls at its top step and li2022heat's hydrolysate rises from 0 to 0.5 %; over a matched dose range what remains is li2022heat's fatty non-enzymatic extract falling at every dose against ding2023's lean fish rising. Fibre competition pushes both the same way, so it cannot flip the sign; formation against thiol-aldehyde trapping is not settled | medium (mechanism: low) |

Both point to the same generator choice: `numeric_observation` keeps only the **longest monotone run** of a ladder
and drops a zero-dose arm (no log of zero). PREREG sec. 4 says a non-monotone ladder's monotone pieces "enter their
own bands"; for dose there are no bands, so the other pieces are silently lost, and which piece survives decides the
sign. The keying itself was checked against the three papers: no transcription error.

## Beyond the sign conflicts: 8 sign changes inside a laboratory, and the 16 widest size spreads

Read 2026-11-12 against their sources ([sign changes inside a laboratory](disagreement_reading/read_sign_changes_inside_a_lab.md);
size spreads [1](disagreement_reading/read_size_spread_1.md) and [2](disagreement_reading/read_size_spread_2.md)).
Every keyed arm that was checked matches its paper; the causes are in the design of the studies and in the table's
own rules.

**The 8 cells where the sign changes inside one laboratory: none needs a laboratory effect.**

| cause | cells | examples |
|---|---|---|
| a segment of a rise-then-fall curve, picked by the longest-monotone-run rule or left by a dropped zero-dose arm | 4 | FFT x time (zhou2023c peaks at 60 min; the table keeps the fall, the rise agrees with feng2022); MFT x cysteine dose (grosch1992 rises then falls, "inhibitory" above 170 - 325 umol); benzaldehyde x xylose dose and phenylacetaldehyde x temperature in li2021ph |
| a second factor the study itself varied | 3 | MFT x pH in hofmann1998: -0.34 per pH unit from cysteine / ribose, +0.33 from fed intermediates, pooled in one cell; pyrazines x sugar identity (the amino acid decides) ; MFT x time |
| scatter of uncorrected peak areas | 1 | phenylacetaldehyde x pH in li2021ph: 15 contrasts, mean -0.04, SD 0.29; its unheated controls already differ by 0.3 dex |

**The 16 cells with the widest size spread among those that agree in direction (5x to 25x):**

| cause | cells | examples |
|---|---|---|
| the reference arm is a different pot | 9 | every amino-acid-identity cell: hwang1995's "none" holds glycine, martin2001's has no amine, park2022mealworm's holds 2 g protein; FFT x glutathione (feng2022's "none" is itself a glutathione Amadori product) |
| a table artefact | 5 | furfural x pH: mottram2002's -2.93 is a detection limit used as a value; FFT x glutathione: grosch1992's -0.05 is a ratio of two "<" limits; furfural x cysteine dose and furfural x temperature: zero or not-detected arms dropped, leaving a plateau or a two-arm slope; phenylacetaldehyde x pH: zamora2015's -0.07 is one piece of a flat ladder |
| where on the reaction curve the window sits | 4 | pyrazines x time (wei2019flax: 0.97 / 0.35 / 0.12 per doubling at 110 / 120 / 130 C over one window); 2,5-dimethylpyrazine x time (a lag step); furfural x temperature (100 - 120 against 150 - 180 C) |
| the lever is not what its label says | 3 | FFT x pH (zhou2023's unbuffered pots; furfural and FFT respond to pH like a switch near pH 5.5 - 6.5, so the final pH decides the slope); a "time" lever that is a pre-heating time (zhou2023c-21); arginine added with no pH correction |
| a second reactant | 1 | phenylacetaldehyde x pH in chan1994b: only that study holds leucine, whose aldehyde condenses with phenylacetaldehyde at higher pH (the authors' explanation) |

(A cell can have more than one cause.) Across all 55 agreed cells, the generated audit (sec. 3b) finds that two
effects from ONE laboratory differ in size by a median factor 1.4 and from two laboratories by 1.6: most of the size
spread is inside laboratories, between conditions. The one condition that separates sizes is the quantification
method (absolute against peak area: +0.48 dex, 90 % +0.25 to +0.72, on three cells only).

**The table's own rules produce part of the disagreement.** Five rules or gaps, each seen in several cells:

1. **A censored arm enters as its limit.** PREREG sec. 3 says a not-detected arm is left-censored AT its limit, a
   bound; `arm_log_and_var` uses the limit as the value. 89 of 1,643 contrasts use such an arm. Leaving them out
   (a check, not a new table) turns the "answered" furfural x pH cell into "too few laboratories" and removes 27
   thin cells.
2. **A non-monotone ladder keeps only its longest monotone run** (see the two cysteine-dose conflicts above).
3. **A zero-dose arm is dropped** (no log of zero), so the 0 to first-dose step, often the largest, is never seen.
4. **Categorical reference arms are treated as the same pot** across studies.
5. **Fed-intermediate and stability experiments share cells with formation from sugar and amino acid.**

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

5. **Non-monotone dose ladders need their pieces kept.** Cysteine both feeds and starves a pathway depending on how
   much there is relative to the sugar; a slope over the longest monotone run discards that. Proposed (backlog):
   keep every monotone piece of a dose ladder as its own contrast, banded by the dose relative to the sugar
   (sulfur donor : sugar below or above about 1 mol/mol), scored beside the current table.

6. **Proposed: a second version of the table's rules** (one amendment bundle, written after the data and scored
   beside the current table, never replacing it silently): censored arms as bounds; every monotone piece of a
   ladder; the zero-dose step as its own "added vs none" contrast; pH by buffered against unbuffered or by final pH;
   reference arms grouped by base mixture; fed-intermediate experiments in their own sub-cells. The test of whether
   it is better: the between-laboratory SD (0.31 dex today) and the leave-one-laboratory-out coverage.

## Limits

- Seven cells, each with two laboratories. The explanations are hypotheses that fit the papers; none is tested.
- No inter-laboratory study exists for thiols, Strecker aldehydes or alkylpyrazines; the measurement figures come
  from olive-oil volatiles, cocoa and plant volatiles.
- Reading notes were produced from text layers for most PDFs; the notes say which were read as page images.
