# The perception layer: scoping pass (STATUS_AND_PLAN sec. 0, step 7)

*2026-11-13. Literature and inventory only; no engine constant moves and nothing is built. The plan in sec. 4 was
premortemed by an independent reviewer the same day (ADOPT WITH EDITS) and the edits are folded in.*

## 1. The question

A formulator does not smell log ratios. The lever table says, for example, that leaving out cysteine lowers
2-methyl-3-furanthiol (MFT) about 6-fold. Whether that changes what a person smells depends on (i) MFT's odour
threshold in the product's matrix, (ii) how far above that threshold the product sits, and (iii) how perceived
intensity grows with concentration. This pass asks what the repository has for each, per odorant of the lever table,
and what a perception reading of a lever cell would need.

## 2. Thresholds on disk

| source | odorants | standing |
|---|---|---|
| `data/species/desirable_targets.yml`, `off_flavour_targets.yml` | 27 | **uncited**: the files' declared source (`pathways.md`) exists nowhere; the files say so. |
| `src/kinetic_core/matrix_oav.py`, water (Zhou 2023 SI Table S2) | 15 | **provenance unknown**: the SI prints no reference for any value. |
| `src/kinetic_core/matrix_oav.py`, 3 % gelatin (Vega) | 6 lipid aldehydes x 4 temperatures | **measured, matrix and method stated** (16 panellists, ascending series, 75 % uncorrected). |
| Thong 2024, Table 4 (on disk, `data/articles/Thong2024.pdf`) | 9 GC-O training standards | quoted from the literature, units and medium as printed there. |

**Coverage of the lever table's 42 odorants:** 5 have a measured threshold in a stated matrix (pentanal, hexanal,
heptanal, (E)-2-octenal, (E,E)-2,4-decadienal: all lipid aldehydes, all in gelatin); 25 have only uncited or
provenance-unknown values; 12 have none (benzaldehyde, total pyrazines, 2-ethylpyrazine, 2-furanmethanol,
trimethylpyrazine, 5-methylfurfural, norfuraneol, 2-ethyl-3,6-dimethylpyrazine, octanal, tetramethylpyrazine,
2-methyl-3-(2-furfuryldithio)furan, (E)-2-nonenal). No meaty thiol, Strecker aldehyde
or pyrazine has a threshold the repository can cite.

**The uncited values are not safe to use.** Bis(2-methyl-3-furyl) disulfide reads 0.02 ug/kg in
`desirable_targets.yml` and 0.0007 ppt (7e-7 ug/kg) in Thong 2024's table: a factor of about 30 000. Published
thresholds for one odorant routinely span two or more orders of magnitude between laboratories and media (Czerny
2008 make the same point); a single number without its medium, method and criterion cannot carry an output.

## 3. What a perception reading of a lever cell needs

1. **A threshold with its medium, method and criterion.** The `ThresholdRecord` schema in `matrix_oav.py` already
   holds these fields; what is missing is the data.
2. **A base level in the product.** A lever ratio scales the odour-activity value (OAV) by the same ratio, so a
   ratio alone says how much the OAV moves, not whether it crosses 1. Whether a 6-fold drop in MFT is smelled depends
   on whether the product starts near threshold or a thousand-fold above it. The engine refuses absolute levels
   (sec. 0, step 5a), so the base level must come from measurements: the user's own (`calibrate`), or published
   levels in a comparable product.
3. **The step from OAV to intensity.** Perceived intensity grows sub-linearly with concentration (a power law with
   an exponent typically well below 1 for odorants), so a 6-fold change in OAV is a smaller change in intensity, and
   odorants mask and enhance one another. A perception reading can therefore be honest only as an ordinal call
   ("crosses threshold in a typical product", "stays far above: intensity shifts", "stays below: not smelled"), with
   the interval of the lever effect carried through.

What makes this useful despite the gaps: in commercial plant-based meats the panel's meatiness tracks the odour
activity of MFT, 2-furfurylthiol and dimethyl trisulfide, and not the protein source (Thong 2024, GC-O with a
semi-trained panel). So the cells that matter most for a perception reading are the thiol cells, which is where the
lever table is thinnest.

## 4. Plan (premortem 2026-11-13: ADOPT WITH EDITS; all nine edits folded in)

The premortem's most likely failure (about 35 %) is that every call comes out undecidable: a base-level range of
about 2 dex, a threshold spread of about 2 dex and a lever interval of up to 2.5 dex multiply to an OAV interval of
5 - 6 dex. So the plan starts with a cheap check that this layer can decide anything, and stops if it cannot.

| # | step | done means | hours |
|---|---|---|---|
| 7.0 | **Decidability desk check.** MFT, FFT and methional with rough numbers (thresholds as printed in two or three sources, a plausible base-level range, the lever interval). **Kill rule:** if fewer than 3 cells can be decided, 7.2 - 7.3 do not start and this record says so | the three cells read decidable or undecidable, with the arithmetic | 1 |
| 7.1 | **Thresholds, one panel and one method first.** Verify Czerny et al. 2008's medium and units before keying (air or water). Key into a curated `data/perception/thresholds.yml` with the `ThresholdRecord` fields plus **medium, units, criterion and route (orthonasal or retronasal) as required keys**; one criterion per odorant unless a conversion is cited. Other measured sources only where Czerny has no value. The uncited target files are not used | each of the 42 odorants has a cited threshold with medium, units, criterion and route, or reads "no measured threshold" | 4 |
| 7.2 | **Base levels.** Starts with a 1-hour feasibility search; **kill rule:** fewer than 2 quantified plant-based-meat sources (isotope dilution or stated calibration) per thiol stops the literature route for that thiol. The **primary input is the user's own measured level** (`calibrate`); published ranges are shown as "illustrative, product-dependent" | per odorant: sources and ranges, or "none quantified" | 6 |
| 7.3 | **A perception column in the lever table**, for the cells whose lever a formulator can set (not "none vs cysteine": protein always brings sulfur). Four calls: **crosses threshold**, **stays far above (intensity shifts)**, **stays below**, **undecidable**. A call other than undecidable is made only when the whole propagated interval lies on one side of OAV = 1, or more than 1 dex above it. **Media are matched** (a water threshold against a water-phase level) or the cell refuses, as `select_threshold` does; an air threshold is never compared with a product level without a cited partition coefficient. **"Stays below" is forbidden** from a water or air threshold in a fat- or protein-rich product (matrix suppression of 1 - 2 dex is the wrong way round for that call). Furfural is stated as probably below threshold whatever the lever does | generated with the table; each call carries its medium and route | 6 |
| 7.4 | **Optional, with WUR or a partner:** GC-O or sensory on the factorial corners (P3 already lists it) | -- | -- |

**Out of scope:** intensity functions (power laws), masking and mixture models, release and retronasal transfer
models. **Cap:** 17 h for 7.0 - 7.3.

**Known gap: dimethyl trisulfide.** One of the three odorants whose odour activity tracked meatiness in Thong 2024,
DMTS appears in four keyed source files but enters no cell of the lever table. A perception layer without it cannot
test that finding; the lever table would need its contrasts keyed as a cell first (a PREREG amendment).

## 5. To read (owner's download list)

- Czerny, M. et al. (2008). Re-investigation on odour thresholds of key food aroma compounds and development of an
  aroma language based on odour qualities of defined aqueous odorant solutions. *Eur. Food Res. Technol.* 228,
  265-273. doi:10.1007/s00217-008-0931-x
- The Wageningen sensory and GC-O study of commercial plant-based meat analogues against beef, *Food Research
  International* 177 (2024) 113848 (verify the DOI when downloading).
- For 7.2, papers quantifying the meaty thiols and Strecker aldehydes in plant-based meats by isotope dilution: to
  be searched at the start of 7.2.
