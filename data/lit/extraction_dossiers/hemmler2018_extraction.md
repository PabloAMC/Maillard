# Hemmler et al. 2018 — EXTRACTION (FT-ICR-MS census of ribose + Gly/Ile/Lys/Cys Maillard products; NO quantitative kinetics)

**Source on disk:** `data/articles/hemmler2018.pdf` (owner's download, 2026-08-28), main text only; the
Supplementary Information (Figs S1–S5, SI Methods) is NOT on disk. Read-only extraction, 2026-09-07, from
`pdftotext -layout`; the paper has no tables — every number below is from the text or a figure caption.
`k1_kinetic_parameters.md` and `k3_final_parameter_inventory.md` already cite this paper (ordinal only) and
`docs/reference/FIT_HOLDOUT_DECLARATION.md` classes it "**neither** — ESI intensity ≠ concentration; pH
uncontrolled and drifting 2–3 units. Ordinal only". This dossier confirms that verdict from the full text.

## 0. Identity

| field | value |
|---|---|
| Title | "Insights into the Chemistry of Non-Enzymatic Browning Reactions in Different Ribose-Amino Acid Model Systems" |
| Authors | Daniel Hemmler, Chloé Roullier-Gall, James W. Marshall, Michael Rychlik, Andrew J. Taylor, Philippe Schmitt-Kopplin (TU Munich / Helmholtz München / Mars Petcare UK) |
| Venue | Scientific Reports (2018) 8:16879 |
| DOI | 10.1038/s41598-018-34335-5 (received 4 May 2018, accepted 12 Sep 2018); CC-BY 4.0 |

## 1. Why it matters — and why it holds no kinetics

The paper is a non-targeted molecular-formula census by direct-infusion ESI(−) FT-ICR-MS. It counts formulae and
reports relative ion intensities; it never converts intensity to concentration, never prints a rate constant,
and its pots are unbuffered and drift 2–3 pH units. **There is nothing in it a rate can be fitted to.** What a
model can use is structural and ordinal: which amino acid makes more products, which sugar reacts faster, which
carbon-backbone families dominate, and that cysteine suppresses sugar degradation and browning.

## 2. Methods as they matter to a model

| item | value | where |
|---|---|---|
| Charge | ribose 0.1 mol/L + amino acid 0.1 mol/L (glycine, isoleucine, lysine or cysteine), equimolar, in Milli-Q water; **unbuffered**; blanks of 0.1 M ribose or amino acid alone | Methods, p. 8 |
| Vessel / atmosphere | 1 mL in 2 mL glass vials, crimp caps, "to exclude additional air/gas exchange" — i.e. ≈ 1 mL air headspace, sealed; no flushing | Methods |
| Temperature / time | 100 °C; sampled at 2, 4, 6, 10 h | p. 2 |
| pH | not stated at start (unbuffered 0.1 M amino acid + ribose); "After ten hours, the pH decreased by approximately 2–3 pH units compared to unheated model systems" | p. 3 |
| Replicates | triplicate; only features present in all three replicates reported | Methods, p. 2 |
| Analysis | direct-infusion FT-ICR-MS, ESI(−), 1:500 in methanol; resolving power 400 000 at m/z 300; > 90 % of formulae within ±0.2 ppm; S/N ≥ 8 for Fig. 5; browning by absorbance at 294 nm (n = 3) | pp. 2–3, Fig. 1 |
| Six-sugar comparison (Fig. S3) | ribose, arabinose, xylose, fructose, galactose, glucose each with glycine, **24 h** (conditions otherwise as above; details in SI, not on disk) | p. 4 |
| Quantification basis | **counts of molecular formulae** and **relative peak intensity** — no concentrations, no response factors | throughout |

## 3. Every number the text holds (there are no tables)

| quantity | value | where |
|---|---|---|
| Distinct molecular formulae after 10 h, four systems together | 1493 | p. 2 |
| MRPs (need both precursors) after 10 h | lysine > 700; glycine, isoleucine, cysteine 300–400 each; order **Lys > Cys > Ile ≈ Gly** | p. 2, Fig. 1a |
| Amino-acid degradation products (no ribose needed) | cysteine 27; lysine "a few"; glycine and isoleucine none | p. 3 |
| Ribose degradation products | 16 ± 2 with Gly, Ile, Lys; **9** with Cys | p. 3 |
| Browning A₂₉₄ after 10 h | order **Lys > Ile > Gly > Cys**; cysteine "only a minor amount of browning over the entire reaction timescale" | p. 3, Fig. 1b (figure only) |
| Glycine-ARP relative intensity, 10 h vs 2 h | −66 % ("degradation rates for the ARPs formed by the other amino acids were significantly lower") | p. 4 |
| N-containing MRPs after 10 h | 89 % (1268/1431) | p. 4 |
| S-free N-containing MRPs in the cysteine pot | 12, < 1 % of intensity | p. 4 |
| Odd-N formulae in the lysine pot | 346, 20 % of intensity (side-chain fragmentation) | p. 4 |
| CHO (N-free) formulae | Gly 34, Ile 39, Lys 130, **Cys 5** | p. 4 |
| Dominant carbon number | ribose + amino acid: Gly C7, Ile/Lys C11, Cys C8; 2 ribose + amino acid: Gly C12, Ile/Lys C16, Cys C13 — the C12 family = 16 % of ribose-glycine intensity | p. 5 |
| "General" MRPs (same core, ≥ 3 of 4 systems) | 73 formulae; explain 45 / 73 / 55 / 46 % of intensity for Gly / Ile / Lys / Cys | p. 6, Fig. 5 |
| Amadori-degradation pathway share, ribose-glycine | 20 % of intensity; the diketosamine pathway explains "up to three times more" in glycine than in the other amino acids | p. 6 |
| ARP as base peak | all except isoleucine-ARP; Cys- and Lys-ARP relative intensity > 30 % | p. 6 |
| Sugar reactivity (formula count, glycine, 24 h) | **ribose > arabinose > fructose ≈ xylose > galactose > glucose** | p. 4 |
| Formula overlap with ribose-glycine, 24 h | arabinose 88 %, xylose 95 %, galactose 75 %, fructose 78 %, **glucose 45 %** | p. 4 |
| Formulae common to all four amino-acid systems | **zero** (apart from ribose degradation products) | p. 4, 8 |
| Lys/Gly MRP count ratio over time | figure only; k3 B4.5 reads ≈ 17× at 2 h, ≈ 2× at 10 h from Fig. 1a | Fig. 1a |

## 4. What the repo could take

No FIT rows. Nothing here is a concentration, a yield or a rate, and the pH is uncontrolled. Directional and
structural items only, all already in `k3_final_parameter_inventory.md` B4.4–B4.6 and confirmed here:

1. **Sugars drive rates, amino acids drive product identity** — 45–95 % formula overlap across six sugars with
   glycine, zero overlap across four amino acids with ribose.
2. Pentoses react faster than hexoses; ribose fastest; glucose slowest (formula count, 24 h, 100 °C, unbuffered).
3. **Cysteine suppresses both browning (A₂₉₄ last of four) and sugar decomposition (9 vs 16 ± 2 ribose degradation
   products), and traps free carbonyls (5 CHO formulae vs 34–130)** — a structural constraint on any sulfur-lane
   step that lets dicarbonyls accumulate in a cysteine pot.
4. The ribose-glycine ARP decays fastest among the four (−66 % relative intensity, 2 → 10 h). The ≈ 0.135 /h k1 once
   derived from that is refused (two points, ESI intensity, drifting pH) — this dossier agrees.
5. Diketosamines (ARP + second sugar) are a major product family (16 % of intensity in ribose-glycine); the trunk
   lane has no such node. Record as a known omission, not a target.
6. Unbuffered 0.1 M ribose + amino acid loses 2–3 pH units in 10 h at 100 °C — the same drift Martins 2005 reports
   for glucose/glycine (6.8 → ≈ 5.5 in 4 h); consistent with the trunk lane's acid production.

## 5. Caveats

- ESI(−) intensity is not concentration; counts of formulae are not amounts. Nothing is scorable as a level.
- pH uncontrolled and unreported at t = 0; "2–3 units" is the only pH datum.
- Sealed vials with ~1 mL air headspace: oxygen neither excluded nor controlled.
- All figures (1a, 1b, 3a, 4, 5 bar charts) are undigitised; the SI with the six-sugar data is not on disk.
- ESI(−) sees polar, oxygen-rich, non-volatile intermediates only; volatiles, pyrazines, furans and melanoidins are
  outside the window, so "reactivity order" here is an order over a proxy the model does not output.
