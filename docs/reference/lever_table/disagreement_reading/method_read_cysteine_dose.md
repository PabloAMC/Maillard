# The two cysteine-dose conflict cells: what the sources say

*Reader's note, 2026-11-12. Literature only; no engine output. Effects are log10 change per doubling of cysteine ("dex").*

## What I read

| file | how |
|---|---|
| ding2023.pdf | Text layer, pp. 1-15 in full (methods, Tables 1-2, all results text). References skipped. Figures 1-6 from captions and text only. |
| zhou2022.pdf | Text layer, pp. 1-9 in full. Table 3 (rotated) read as a page image (p. 9, journal p. 15210). Supporting information (calibration Table S1, browning) **not read**. |
| li2022b.pdf | Text layer, pp. 1-9 in full (methods, Tables 1-4B, results). Fig. 1 (pH) and Fig. 2 (xylose) from text only. |
| contrast files | ding2023, zhou2022, li2022heat `__families.yml`: the 2-methylpyrazine and nonanal entries, checked row by row against the tables. |
| generator | `scripts/generators/build_lever_table.py`, `numeric_observation` and `_lever_x` (how a ladder becomes one slope). |
| earlier notes | DISAGREEMENT_READING.md, method_read_thiols.md (in full), method_read_interlab.md (SPME competition section only). |
| method papers | **Not re-read.** Roberts 2000 competition figures are taken from method_read_interlab.md. The zhou2023b final-pH figures are taken from method_read_thiols.md / DISAGREEMENT_READING.md. |

## One finding that bears on both cells: how a ladder becomes a slope

`numeric_observation` takes log2 of the dose. For the 0 % arm, log2(0) is undefined, so **that arm is silently dropped**. The function also keeps only the **longest monotone run** of a non-monotone ladder. All three sources have a 0 arm, and two of them have non-monotone ladders. So each keyed effect is a slope over a *segment* of the ladder, and different laboratories contribute different segments:

| contrast | arms printed | arms used for the keyed slope | segment shape |
|---|---|---|---|
| ding2023 2-MP | 0, .25, .5, .75, 1.0 % | 0.5, 0.75, 1.0 | falling part of a rise-then-fall curve |
| zhou2022 2-MP pH 7 | 0, 6, 14, 20 mM | 14, 20 (0 dropped; 6 is ND) | rising |
| zhou2022 2-MP pH 8 | 0, 6, 14, 20 mM | 6, 14, 20 | rising |
| ding2023 nonanal | 0, .25, .5, .75, 1.0 % | 0.25, 0.5, 0.75 | rising part of a rise-then-fall curve |
| li2022heat nonanal (4 ladders) | 0, 0.5, 1.0 g/100 mL | 0.5, 1.0 only | falling |

This is not a keying error. It is a choice in the generator, and it explains a large part of both conflicts.

---

## Cell 1: 2-methylpyrazine against cysteine dose (ding2023 -0.43; zhou2022 +0.20, +0.35)

### The two designs

| | ding2023 (ZJUT, Zheng/Zhou) | zhou2022 (Jiangnan, Ho network) |
|---|---|---|
| base mixture | Snakehead fish-waste protein hydrolysate (amount unprinted) + D-xylose 0.25 % (16.7 mM). Glycine is named in the title, but no amount is printed and the methods (p. 3) do not list it. Base free amino acids include Gly ~1.3 g/100 g (Table 1; the basis is unclear). | Purified xylose-alanine Amadori product (AX-ARP) 20 mM in water. No free sugar. The only amine is the alanine released from the ARP. |
| cysteine-free arm | **Already makes pyrazines**: total N-compounds 1119 ng/g, 2,5-DMP 852, 2-MP 54.9 ng/g (Table 2) | **Makes none**: pyrazine and 2-MP ND at pH 7 (Tables 1, 3). At pH 8, 2-MP is 3.66 µg/L. |
| cysteine dose | 0 / 0.25 / 0.5 / 0.75 / 1.0 % = 0 / 21 / 41 / 62 / 83 mM. **Cys : xylose = 0 / 1.2 / 2.5 / 3.7 / 5.0 mol/mol** | 0 / 6 / 14 / 20 mM. **Cys : ARP = 0 / 0.3 / 0.7 / 1.0** |
| keyed segment | 0.5 → 1.0 % (Cys : sugar 2.5 → 5) | 14 → 20 mM (pH 7); 6 → 20 mM (pH 8) (Cys : sugar 0.3 → 1) |
| pH | Initial 7.0 (1 M NaOH), unbuffered, final pH not reported | Initial 7 or 8 (3 M NaOH), unbuffered, final pH not reported here. The same lab's identical system (zhou2023b: 20 mM AX-ARP + 20 mM Cys) ended near pH 3.3-3.4 from an initial pH 7, and near 5.0 from an initial pH 8. |
| heat | 110 °C, 60 min, sealed tubes, oil bath | 120 °C, 60 min, pressure bottle, oil bath |
| matrix | Lean protein hydrolysate (centrifuged supernatant) | Water |
| extraction | HS-SPME, 5 mL, 55 °C 10 min + 40 min. **Fibre type not stated.** DB-5MS 60 m. | HS-SPME DVB/CAR/PDMS, 3 mL + 2 mL sat. NaCl, 60 °C 30 min. DB-WAX. |
| quantification | Semi-quantified: peak area against 1,2-dichlorobenzene IS (≈50 ng/mL in the vial). "ng/g". NIST identification only. | Standard curves against the 1,2-dichlorobenzene IS (Table S1, not read), µg/L. Identified with authentic standard + RI. |
| cross-check | GC-IMS of the same pots: 2-MP "firstly increased and then decreased" (p. 13) | — |

### Best-supported explanation: the same curve sampled on opposite sides of its peak, with a different base mixture

- **Inside ding2023, 2-MP rises and then falls.** Per step: 0 → 0.25 % +0.16 dex, 0.25 → 0.5 % +0.05, 0.5 → 0.75 % -0.30, 0.75 → 1.0 % -0.14 (Table 2 row 31: 54.86, 78.98, 88.70, 44.69, 32.75 ng/g, with letters c, b, a, d, e, so every step is significant). The peak sits at Cys : xylose ≈ 2.5. The keyed -0.43 comes only from the falling half. The low-dose half has **zhou2022's sign**, and GC-IMS confirms the shape independently, so it is not a GC-MS peak-area artefact.
- **zhou2022 samples only the low-dose side** (Cys ≤ ARP, mol/mol), and its base has no free amine. The paper's own data and discussion say that cysteine is the missing nitrogen donor. The ARP alone gives no pyrazines because the recovered alanine is out of step with dicarbonyl formation. Added cysteine captures GO and MGO (both "greatly reduced", Fig. 1e,f) and "made up for the deficiency of pyrazines" (pp. 15205-15206). Its degradation products (NH3, H2S) open a second route (p. 15208). So pyrazines rise with cysteine at both pH values. At pH 8 the rise is about 100-fold larger (20 mM Cys: 96.3 against 1.12 µg/L at pH 7). This fits the sister paper's finding that an initial pH of 7 ends near pH 3.3, where pyrazine chemistry is nearly off.
- **ding2023's base already has an excess of amine** (glycine plus peptides), and pyrazine formation is limited by carbonyls. Its own explanation is competition: cysteine competes "with glycine with xylose" (p. 13). Its data agree: A294 and A420 fall with cysteine (Fig. 1a, text p. 6), alanine falls 0.134 → 0.092 (Table 1), total pyrazines fall 1105 → 39 ng/g (p. 12), and the Strecker aldehydes 3-methylbutanal (5.27 → ND) and phenylacetaldehyde (13.5 → 4.9) fall too. Once cysteine is 2.5-5 times the sugar on a molar basis, it takes the dicarbonyls that would have made methylpyrazine.
- **So the lever is not the same.** Cysteine works as an N donor in an amine-free system at Cys ≤ sugar. It works as a carbonyl sink in an amine-rich system at Cys ≫ sugar. Both papers agree with each other where their designs overlap: low Cys : sugar raises 2-MP.

### Measurement: could a lever-linked bias of 0.1-0.45 dex do it?

No. The gap between the laboratories is 0.6-0.8 dex. In ding2023 the bias most likely to move with the lever, SPME fibre competition, points the wrong way. Total volatile load falls with cysteine (about 1.5 µg/g in arm A against about 0.64 µg/g in E, summed from the Table 2 class totals), which would favour 2-MP recovery at high cysteine. The keyed slope is negative. The loads are also near the ~1 ppm level below which Roberts 2000 saw little competition. In zhou2022 the load is of the same order (furfural falls, sulfur volatiles rise), and its calibration is the same in every arm. The GC-IMS agreement inside ding2023 also argues against a GC-MS artefact.

### Keying check

- ding2023-17: values, SDs and arms match Table 2 row 31 (p. 10). Unit ng/g as printed. `peak_area` is correct. Glycine "unstated" is correct (no amount printed anywhere I read). **Add:** the SPME fibre is unstated.
- zhou2022-7 and zhou2022-8: the values match Table 3 (pH 7: ND, ND, 0.88 ± 0.13, 1.12 ± 0.11; pH 8: 3.66, 18.30, 77.87, 96.30). **Flag:** zhou2022-7 rests on the 14 → 20 mM pair alone, and the paper's own Duncan letters do not separate it (both "a"). The real pH 7 effect is ND → detected between 6 and 14 mM, which the slope cannot use.
- No transcription error found.

### What would confirm or refute

- **Confirm:** one design that runs Cys : sugar from about 0.25 to 5 in both an amine-free base (an ARP) and an amine-rich base (sugar + glycine or a hydrolysate), and records the final pH. The prediction: 2-MP peaks, and the peak sits at a lower Cys : sugar when the base is rich in amine. Already within reach: other lever-table cysteine ladders, keyed with their Cys : sugar ratio and the free-amine status of the base, should change sign at that ratio and not by laboratory.
- **Refute:** zhou2022's system still rising at Cys : ARP of 3-5, or an amine-rich system that falls already at Cys : sugar below 1.

**Confidence: medium-high** that this is a dose-regime and base-mixture difference, not measurement. **Medium** on the exact mechanism split: in zhou2022, N donation against the NH3 from cysteine degradation; in ding2023, dicarbonyl capture against browning inhibition.

---

## Cell 2: nonanal against cysteine dose (ding2023 +0.22; li2022heat -0.23, -0.13, -0.44, -0.39)

### The two designs

| | ding2023 | li2022heat (li2022b; NUS, Liu) |
|---|---|---|
| base mixture | Lean fish-waste hydrolysate supernatant + xylose 0.25 % (+ glycine, amount unprinted) | Pork-trimming extract, **10 % fat (dry basis)**. Either NPE (non-enzymatic; free amino acids 222 mg/100 mL) or EPH (Flavourzyme hydrolysate; 3143 mg/100 mL), each with xylose 0 or 3.5 g/100 mL (233 mM) |
| cysteine-free arm | Heated arm A: nonanal 24.3, **half the unheated F (48.1)**. It also has the most browning and the highest pyrazine load. | N0-0 10.3, E0-0 5.4, N35-0 5.7, E35-0 6.7 ng/mL. Lipid oxidation dominates (hexanal 171-345, 2-pentylfuran up to 1802 ng/mL). |
| cysteine dose | 0-1.0 % (Cys : xylose 1.2-5 mol/mol) | 0 / 0.5 / 1.0 g/100 mL = 0 / 41 / 83 mM. Cys : xylose 0.18 / 0.35, or no sugar at all. In NPE, cysteine outnumbers the native free amino acids. |
| keyed segment | 0.25 → 0.75 % (rising run) | 0.5 → 1.0 only (the 0 arm is dropped) |
| pH | Initial 7.0, unbuffered | Initial 5.50 (lactic acid / NaOH), unbuffered. After heating: xylose arms ≈ 5.15 (NPE) / 4.75 (EPH). Cysteine-only arms show no significant change (text p. 4, Fig. 1). |
| heat | 110 °C, 60 min | 100 °C, 60 min (autoclave, 250 mL bottles, 100 mL) |
| extraction | HS-SPME, fibre unstated, 55 °C 50 min total, DB-5MS | HS-SPME CAR/PDMS 85 µm, 60 °C 30 min, agitated; polar column (nonanal LRI 1392) |
| quantification | Peak area against 1,2-dichlorobenzene, ng/g | Peak area against 2-octanol (120 ng/mL), "semi-quantified", ng/mL |
| cross-check | GC-IMS: nonanal "increased firstly and then decreased" (p. 14) | — |

### Best-supported explanation: mostly segment selection on non-monotone curves; a residual matrix difference

- **ding2023 rises and then falls**: 24.33, 37.44, 51.35, 96.64, 75.29 ng/g (Table 2 row 48; letters e, d, c, a, b). Its last step, 0.75 → 1.0 %, is -0.11 dex, i.e. **-0.26 dex per doubling: li's sign.** The keyed +0.22 comes from the rising run only.
- **li2022heat EPH also rises and then falls**: E0 5.39 → 8.69 → 3.12, and E35 6.65 → 15.31 → 6.29 (Table 4B). From 0 to 0.5 % that is **+0.21 and +0.36 dex: ding's sign.** The keyed -0.44 and -0.39 come from the falling 0.5 → 1.0 step only.
- **What remains a real difference.** li's NPE series falls at every dose (N0 10.33 → 7.21 → 4.27; N35 5.69 → 5.40 → 3.96; Table 4A). On the one dose range both studies print, 0.5 → 1.0 %, ding2023 is +0.17 dex (it passes its peak at 0.75 %) against li2022heat's -0.13 to -0.44. So part of the conflict survives a matched segment. It concerns where the peak sits, and the NPE base, which is low in amine and high in fat, where no rise is seen.
- **Formation, trapping or extraction?** The papers cannot settle it. Neither discusses nonanal mechanistically. Ding says only that aldehydes "firstly increased and then decreased… due to the inhibitory effect of excessive cysteine" (p. 13). Li says hexanal fell 30-90 % on heating "depending on the amount of xylose and cysteine" (p. 8). What the data do show:
  - **In ding2023, heating without cysteine *consumes* lipid aldehydes** (F → A: nonanal 48 → 24, hexanal 13.6 → 7.6, heptanal 8.1 → 4.3, decanal 10.0 → 6.3). That happens in the arm that browns most and makes the most pyrazines. Cysteine reduces browning (Fig. 1a), and the lipid aldehydes come back, up to 2× above unheated at 0.75 %. My inference: in the glycine/peptide + xylose base at pH 7, cysteine *spares* the aldehydes from amine-carbonyl consumption, or adds net oxidation. Thiol trapping cannot be the dominant term there, because trapping would lower nonanal as cysteine rises.
  - **In li2022heat, the whole lipid-aldehyde class falls with cysteine**, in NPE at every step (hexanal N0 345 → 61 → 22; pentanal 49 → 12 → 10) and in EPH from 0.5 to 1.0 %. That is consistent with cysteine acting as an antioxidant or forming thiol-aldehyde adducts (thiazolidines) at pH 5.5. The data cannot tell these two apart.
  - **Extraction** is argued against below.
- **Same lever?** Only nominally. In ding2023, cysteine enters a carbonyl-limited, strongly browning base at pH 7 that already uses up aldehydes. In li2022heat, it enters a fatty, oxidising base at pH 5.5, with 0.18-0.35 mol per mol of sugar, or with no sugar at all. Nonanal is a lipid-oxidation product, and what cysteine does to it depends on what else consumes or forms it in that base.

### Measurement

Not the sign. In both laboratories the total volatile load *falls* with cysteine. ding2023: about 1.5 → 0.64 µg/g. li2022heat EPH: 2-pentylfuran 1802 → 748 → 435 and hexanal 171 → 101 → 56 ng/mL. On a Carboxen fibre, less competition at high cysteine raises nonanal's apparent recovery **in both**. That bias has the same direction in both laboratories, so it cannot produce opposite signs. At most it inflates ding2023's rise and damps li2022heat's fall, by an amount the documented 0.1-0.45 dex range bounds from above. li2022heat's EPH load (~2-3 µg/mL) is above Roberts' ~1 ppm onset, so the damping there may be real. Thiol-type matrix binding does not apply to an aldehyde (Hofmann 2001: melanoidin binding of short aldehydes ≤ 0.04 dex, per the earlier notes). The fat content cancels within each study, because cysteine does not change it.

### Keying check

- ding2023-25: matches Table 2 row 48 (p. 10). Correct.
- li2022heat-59, -60, -64, -65: match Tables 4A and 4B exactly (N and E series, xylose 0 and 3.5). Lever values in g/100 mL are correct. The unheated "Controls" are rightly excluded.
- No transcription error found. The flag is the slope rule above: li2022heat contributes only its 0.5 → 1.0 step, and ding2023 only its 0.25 → 0.75 run.

### What would confirm or refute

- **Confirm (from existing data):** rescore the cell on matched dose segments, with the 0 arm kept (for example, slopes against dose or Cys : sugar instead of log2 dose). The prediction: the sign difference shrinks to li2022heat's NPE series against ding2023's 0.5 → 1.0 range.
- **Separate formation, trapping and extraction (new data):**
  1. Spike d-labelled nonanal into each arm *before* heating. Recovery that falls with cysteine means trapping or consumption.
  2. Spike it *after* heating, before SPME. Recovery that changes with cysteine means extraction.
  3. Measure the cysteine-nonanal thiazolidine by LC-MS, and hexanal/nonanal precursors (hydroperoxides or TBARS) for formation.
- **Refute:** a matched-segment rescore that still leaves ding2023's EPH-like rise opposite to li2022heat's EPH at the same Cys : sugar.

**Confidence: medium** that the conflict is mainly the segment rule plus a base-mixture difference, not measurement. **Low** on the chemical mechanism (formation against trapping), which neither paper tests.

---

## For the backlog (proposed, not applied)

1. **Precursor-dose ladders with a 0 arm**: drop-arm-and-longest-run hides non-monotone dose responses, and it chose the conflicting segments in both cells. Options: key Cys : sugar (mol/mol) as the band coordinate; flag non-monotone ladders as their own class ("peaked") instead of a slope; keep the 0 → first-dose step as a separate pair contrast.
2. **Record for each cysteine ladder:** the Cys : reducing-sugar molar ratio, and whether the base holds free amine.
3. **zhou2022-7:** note that the only usable pair is not significant in the source (Duncan letters a, a).
4. **ding2023:** record the SPME fibre as unstated, and note identification by NIST only (the table also prints duplicate rows: (E)-2-octenal and isopinocarveol).
