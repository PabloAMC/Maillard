# Step 7.0: can a perception reading decide anything? (desk check)

*2026-11-13. Plan: [SCOPING_2026-11.md](SCOPING_2026-11.md) sec. 4, row 7.0. Desk arithmetic on numbers already on
disk; nothing is keyed, built or fitted. **Result: 0 of 3 cells decidable for the call that matters ("crosses
threshold"); the kill rule fires and steps 7.2 - 7.3 do not start.***

## 1. The three cells

The answered lever-table cells for the three odorants ([lever_table.md](../../../results/validation/lever_table.md)),
with the interval a new laboratory (a formulator's kitchen) should expect, in log10 units (dex):

| odorant | lever | estimate (90 % CI) | new laboratory, 90 % = L |
|---|---|---|---|
| 2-methyl-3-furanthiol (MFT) | none vs cysteine | -0.78 (-1.16 to -0.39) | -1.77 to +0.22 |
| 2-furfurylthiol (FFT) | pH, +1 unit, pH 6 - 7.5 (initial or label) **[fragile]** | -0.74 (-1.28 to -0.16) | -1.95 to +0.50 |
| methional | xylose vs glucose | +0.51 (+0.45 to +0.56) | +0.41 to +0.60 |

## 2. Thresholds and the assumed spreads

Water thresholds on disk: MFT 0.007 ug/kg (Kerler 2010, Table 3.2); FFT 0.01 ug/kg (Kerler 2010) and 0.005 ppb
(Thong 2024, Table 4); methional 0.2 ug/kg (Kerler 2010) and 0.2 ppb (Thong 2024). Two assumptions, both kinder than
the premortem's 2-dex threshold spread:

- **between-laboratory spread of the threshold: +/-0.5 dex** (the two FFT values differ by 0.3 dex);
- **matrix suppression in a fat- or protein-rich product: 0 to +2 dex** on the effective threshold (scoping sec. 4).

Write x = log10(base level / water threshold). The product's odour-activity value before the lever is then known to
lie in **[x - 2.5, x + 0.5]** (width 3 dex) even when the base level x is known exactly; after the lever it lies in
that interval shifted by L.

## 3. The arithmetic, per call (rules of 7.3)

**"Crosses threshold"** needs the before and the after interval on opposite sides of OAV = 1.

- **MFT, FFT: impossible, whatever the thresholds and base levels.** Their L includes 0 (a new kitchen may see no
  change or the opposite sign), so the after interval always overlaps the before interval. No threshold work and no
  base-level work can change this; only a lever interval that excludes zero can.
- **Methional: impossible with literature inputs.** Upward crossing needs x + 0.5 < 0 (before entirely below) and
  x - 2.5 + 0.41 > 0 (after entirely above): x < -0.5 and x > 2.09 at once. In general a crossing call needs the
  combined width of threshold and base-level uncertainty below L's lower end, 0.41 dex. Even a threshold measured in
  the formulator's own product (about 0.3 dex) with a measured base level (about 0.15 dex each way) gives 0.6 dex.

**"Stays far above (intensity shifts)"** needs the whole after interval more than 1 dex above OAV = 1:

| odorant | condition | base level needed in the product |
|---|---|---|
| MFT | x - 2.5 - 1.77 > 1, so x > 5.27 | above about 1,300 ug/kg (1.3 mg/kg) |
| FFT | x - 2.5 - 1.95 > 1, so x > 5.45 | above about 2,800 ug/kg |
| methional | x - 2.5 > 1 (the lever only raises it), so x > 3.5 | above about 630 ug/kg |

For the two thiols these are levels at which a trace odorant would be a bulk volatile: not a product anyone makes.
For methional the level is not absurd, but the call then says only what the lever row already says (the note rises
by 0.4 - 0.6 dex in concentration); it adds no perception information.

**"Stays below"** is forbidden by the 7.3 rule from a water threshold in a fat- or protein-rich product.

## 4. Verdict and what it implies

**Decidable cells: 0 of 3** for "crosses threshold"; "stays far above" is reachable only at base levels no product
has (MFT, FFT) or restates the lever (methional). **Kill rule (fewer than 3): fires. Steps 7.2 (base levels) and 7.3
(perception column) do not start.** This is the premortem's most likely failure (about 35 %), and the binding
constraint is not the one the premortem expected: before thresholds or base levels, **the lever intervals themselves
include zero for the thiols**. A perception layer on the literature table cannot outrun the table.

- **7.1 (keyed thresholds) is parked too.** With 7.2 - 7.3 stopped nothing would consume it; the uncited thresholds
  in `desirable_targets.yml` already say they are uncited and feed no shipped call.
- **What would reopen the layer:** a thiol cell whose new-laboratory interval excludes zero by more than about
  1 dex (needs more laboratories, the same binding gap as STEP5_RECORD), or a formulator who brings a threshold
  measured in their own product together with a measured base level, and then only for levers stronger than about
  0.6 dex.
- **For a formulator, today:** whether a lever changes what a person smells is decided by a sensory or GC-O test on
  their own product (row 7.4; the factorial's optional GC-O arm), not by any table. The experiments brief's
  experiment 3 (thresholds in pea protein, on a panel) is the measurement that would shrink the threshold side.
