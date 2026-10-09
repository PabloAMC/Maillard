# Milani & Conti 2024 — EXTRACTION (1.5 % w/w thiamine added to soy protein concentrate before single-screw extrusion, 160 °C final zone; the textured protein used in a meat analogue and a soy burger; hedonic and RATA sensory, no volatiles)

**Source on disk:** `data/articles/Milani2024.pdf` (the publisher's PDF, 10 pages, with a text layer). Read and
checked by eye on 2026-10-09: every number below was read from the page image (pp. 744-750) and
cross-checked against `pdftotext -layout`; Tables 1-4 match the text layer cell for cell. Figs 1 and 2 are PCA
biplots with no numbers taken from them. The supplement (processing flowchart, Supplementary Tables 1 and 2
with water absorption, compression force, yield, shrinkage, TPA) is not on disk; the few supplement numbers
quoted in the text are given as printed. Written for the thiamine row of
`results/validation/cultivated_tissue_invariance_prereg.md` (sections 10-11), not for the core fit.

**Companions on disk.** `conti2025_extraction.md` and `conti2025b_extraction.md` are the same laboratory,
the same SPC (Arcon SM) and the same 1.5 % w/w thiamine, extruded at three moisture/temperature pairs;
Conti 2025b holds the volatile inventory (peak areas scaled to one deuterated standard) that this paper
lacks. The process itself comes from Milani, Menis-Henrique & Conti 2022 (J Food Process Preserv 46:e16731),
which is not on disk.

| field | value |
|---|---|
| Title | "Textured soy protein with meat odor as an ingredient for improving the sensory quality of meat analog and soy burger" |
| Authors | Talita Maira Goss Milani, Ana Carolina Conti (UNESP, IBILCE, São José do Rio Preto, Brazil) |
| Venue | Journal of Food Science and Technology (April 2024) 61(4):743-752; revised 9 Aug 2023, accepted 14 Oct 2023 |
| DOI | 10.1007/s13197-023-05875-0 |

## 1. Methods

**Material (p. 743-744).** Soy protein concentrate Arcon SM (ADM Foods & Wellness); thiamine supplied as
**thiamine hydrochloride** (Sigma-Aldrich, purity > 99 %). SPC moisture adjusted to 28.4 % or 34 % (dry
basis), sealed and refrigerated 48 h, re-measured; then **"1.5% thiamine (w/w) was added to the SPC two hours
before the extrusion"** (planetary mixer, 150 rpm, 1 min), held at room temperature. Whether 1.5 % is on the
moistened or the dry SPC, and whether it counts the hydrochloride or the free base, is not printed.

**Three TSPs (p. 744).** 34 % M without thiamine (control); 28.4 % M with thiamine; 34 % M with thiamine; all
at 216 rpm, selected from the 2022 central composite design. Each processed in two repetitions (six samples,
about 18 kg; 3 kg SPC each); the control extruded first, the others randomised.

**Extrusion (p. 744).** RXPQ Labor 24 single-screw extruder (INBRAMAQ), five heating zones at **30, 60, 80,
145, 160 °C**; helically grooved barrel; screw compression ratio 2.3:1; L/D 15.5:1; pre-die holes 5.5 mm;
die 3.2 mm round; feed rate 170 g/min; screw 216 rpm. **Residence time, melt temperature, pressure and SME
are not printed.** Product ground (30-60 s), sieved to a set granulometry, stored sealed in aroma-barrier
laminate at room temperature in the dark; the storage time before the sensory work is not printed.

**Products (p. 745-746).**
- Meat analogue (%): TSP 18.9, water 75.5, soybean oil 3.5, salt 1.0, MSG 0.4, onion powder 0.3, dried
  parsley 0.2, dried chives 0.2. TSP hydrated 1:4 in water preheated to 100 °C, 15 min; pan-fried to 72 °C
  at the centre, held 2 min; served at ≥ 50 °C within 30 min.
- Soy burger (%): TSP 19.2, water 57.5, soybean oil 9.7, warm water (~50 °C) 9.7, CMC 2.0, salt 1.0, MSG
  0.4, onion powder 0.4. Hydrated 1:3 w/w at 100 °C, 15 min; 56 g patties, 85 mm, frozen at −18 °C; grilled
  from frozen 3 min per side at level 3, then 1 min per side at level 5 (8 min total).

**Sensory (p. 745-746).** Consumers, individual booths, white light, 22 °C. Nine-point structured hedonic
scale (1 disliked extremely to 9 liked extremely) for appearance, odour, texture, flavour, overall; RATA
with attributes from a five-consumer focus group (15 for the analogue, 17 for the burger), intensity
low/medium/high, coded 0-3 with "not applicable" = 0 (Meyners 2016). Monadic, balanced complete blocks
(FIZZ 2.50), Williams Latin square for attribute order; about 10 g (analogue) or one quarter burger (about
14 g) per sample.
- Meat analogue panel: 66 consumers (68 % female, 86 % aged 18-30, 92 % eat meat). Table 1 is n = 132,
  Table 2 n = 66; the text does not explain 132 (two process repetitions × 66 is the obvious reading, not
  stated).
- Soy burger panel: consumer count not printed in the text; Table 4 is n = 70, Table 3 n = 140 (70 % female,
  87 % aged 18-30, 86 % eat meat).

**Statistics (p. 746).** One-way ANOVA for physical properties; two-way ANOVA (sample, consumer) for
acceptance and RATA; Tukey, P ≤ 0.05; PCA on standardised means (Statistica 10).

**No volatile was measured.** The authors say so (p. 751): "we did not evaluate the volatile compounds of
the products".

## 2. Findings that matter

### 2.1 Meat analogue (Table 1, n = 132, and Table 2, n = 66; p. 747)

| | 34.0 % M, no thiamine | 28.4 % M, thiamine | 34.0 % M, thiamine |
|---|---|---|---|
| acceptance, odour (1-9) | 5.6 ± 1.3 b | 6.1 ± 1.6 a | 6.3 ± 1.5 a |
| acceptance, flavour | 6.3 ± 1.4 ns | 6.2 ± 1.8 ns | 6.2 ± 1.8 ns |
| acceptance, overall | 6.2 ± 1.4 ns | 6.2 ± 1.7 ns | 6.2 ± 1.6 ns |
| acceptance, appearance / texture | 6.4 ± 1.4 / 6.3 ± 1.5 ns | 6.2 ± 1.5 / 6.4 ± 1.5 ns | 6.2 ± 1.6 / 6.3 ± 1.5 ns |
| RATA meat odour (0-3) | 0.61 ± 0.63 b | 0.97 ± 0.76 a | 0.94 ± 0.78 a |
| RATA soy odour | 1.89 ± 1.01 a | 1.39 ± 0.91 b | 1.48 ± 0.90 b |
| RATA meat flavour | 0.94 ± 0.78 ns | 1.05 ± 0.77 ns | 1.08 ± 0.79 ns |
| RATA burnt aftertaste | 0.18 ± 0.43 b | 0.70 ± 0.94 a | 0.80 ± 1.03 a |
| RATA aromatic | 1.18 ± 0.72 b | 1.68 ± 0.90 a | 1.56 ± 0.86 a |
| RATA salty taste | 1.56 ± 0.64 b | 1.70 ± 0.76 ab | 1.88 ± 0.87 a |

Other RATA attributes not significant (white and caramel colour, uniform granules, rubber, crumbling/sandy,
wet texture, noodle and seasoning flavour, tasty). Physical (text, p. 747, from Supplementary Table 1): the
28.4 % M TSP had lower water absorption (484 % vs 656 %) and higher compression force (10.6 vs 6.9 N) than
the control; no difference between the two 34 % M samples, so the authors assign the change to moisture,
not thiamine.

### 2.2 Soy burger (Table 3, n = 140, and Table 4, n = 70; p. 749-750)

| | 34.0 % M, no thiamine | 28.4 % M, thiamine | 34.0 % M, thiamine |
|---|---|---|---|
| acceptance, odour (1-9) | 6.7 ± 1.2 b | 7.2 ± 1.2 a | 7.1 ± 1.1 a |
| acceptance, flavour | 6.3 ± 1.5 b | 6.9 ± 1.2 a | 6.8 ± 1.3 a |
| acceptance, overall | 6.4 ± 1.3 b | 6.9 ± 1.2 a | 6.8 ± 1.2 a |
| acceptance, appearance | 7.5 ± 1.2 a | 7.3 ± 1.1 ab | 7.2 ± 1.1 b |
| acceptance, texture | 6.4 ± 1.4 b | 6.8 ± 1.2 a | 6.7 ± 1.3 ab |
| RATA meat odour (0-3) | 0.84 ± 0.83 b | 1.20 ± 0.93 a | 1.06 ± 0.90 ab |
| RATA chicken odour | 0.76 ± 0.79 b | 1.19 ± 1.00 a | 1.04 ± 0.91 a |
| RATA soy/vegetable odour | 1.23 ± 0.90 a | 0.99 ± 0.89 b | 1.00 ± 0.85 b |
| RATA bacon odour | 0.49 ± 0.76 ns | 0.47 ± 0.72 ns | 0.40 ± 0.69 ns |
| RATA chicken flavour | 0.97 ± 0.85 ns | 1.04 ± 0.79 ns | 1.06 ± 0.96 ns |
| RATA soy flavour | 1.31 ± 0.86 ns | 1.10 ± 0.95 ns | 1.11 ± 0.81 ns |
| RATA caramel colour | 2.23 ± 0.75 a | 2.10 ± 0.76 ab | 1.96 ± 0.71 b |
| RATA aromatic | 1.51 ± 0.81 b | 2.10 ± 0.84 a | 1.93 ± 0.87 a |

Yield > 85 %, shrinkage < 1.8 % for all burgers, lower for 28.4 % M with thiamine (0.5 %); no TPA
differences (text, p. 749, from Supplementary Table 2).

### 2.3 What the effect sizes are (derived here)

Thiamine raised RATA meat odour by +0.33 to +0.36 (analogue) and +0.22 to +0.36 (burger) on a 0-3 scale, and
lowered soy odour by −0.41 to −0.50 (analogue) and −0.23 to −0.24 (burger). The control already scored meat
odour 0.61 and 0.84, between "not applicable" and "low"; the authors attribute it to native thiamine in SPC,
citing 0.02 mg/100 g (Deak 2008) to 0.32 mg/100 g (Perkins 1995), not measured here. Meat odour moved; meat
or chicken FLAVOUR did not (ns in both products). Odour acceptance rose by 0.5-0.7 points (analogue) and 0.4-0.5 points (burger); flavour and overall acceptance rose only in the burger (9.7 % oil, grilled 8 min), which the
authors connect to the oil phase carrying the volatiles. Burnt aftertaste is the cost in the analogue
(0.18 → 0.70-0.80).

## 3. What it means for the model

**The dose is not a beef-level dose (derived here).** 1.5 g thiamine HCl (337.27 g/mol) per 100 g SPC is
4.45 mmol. In the water of the moistened SPC that is 131-201 mM, depending on the unprinted basis (on dry SPC:
131 mM at 34 % d.b., 157 mM at 28.4 %; on moist SPC: 175 and 201 mM). The sourced beef range carried by the
cultivated-tissue box is 0.00044-0.0040 mM (`lombardiboccia2005_extraction.md`), so this is 3 × 10⁴ to
5 × 10⁵ times beef; against the native SPC thiamine the authors cite (0.02-0.32 mg/100 g) it is about
5 × 10³ to 8 × 10⁴ times.

**Relation to the thiamine row (prereg sections 10-11).** The row is a claim about restoring thiamine to its
BEEF level in a ribose/cysteine/glucose pot: +24 % on the MFT metric at 100 °C / 20 min, +4 % at 140 °C /
5 min (section 11, corrected box). Milani tests something else: a ten-thousand-fold overdose in a dry-ish
melt (28-34 % d.b.) at 145-160 °C for an unprinted residence time, read by consumers. So it neither confirms
nor refutes the row's magnitude. What it does say:

- **Direction at high temperature, at high dose.** Thiamine produces a consumer-detectable meat (and
  chicken) odour and suppresses perceived soy odour after a short, hot process. That matches section 10's
  reading of Thomas 2014 (a clear MFT rise needed thiamine far above native): high temperature does not
  stop thiamine from making meaty odorants, it is the beef-level dose that makes the 140 °C contribution
  small next to ribose. Milani says nothing about which odorant (MFT, thiazoles or other thiamine
  fragments); Conti 2025b, same material and dose, inventories thiazoles and reports that the HMP
  (5-hydroxy-3-mercapto-2-pentanone) "was not found" (`conti2025b_extraction.md`, quoting that paper's
  section 3.2).
- **Jhoo's sink applies here.** At 131-201 mM thiamine, the pyrimidine fragment that traps MFT as an almost
  odourless thioether (section 10, Jhoo 2002) is at its most relevant; the engine lacks that step, so any
  engine run on this recipe would be biased upward for MFT.
- **Engine arithmetic for scale (live values).** `k_thi_hmp` is frozen at log10 k = −2.588 per minute at
  145 °C (`results/validation/kinetic_core_b9_fit_report.json`, `frozen_parameters.log10_k_ref_at_145C`,
  read by `src/kinetic_core/engine.py`), carried by the lumped formation barrier 64.08 kJ/mol (same file,
  `lumped_formation_Ea_kJ_mol`). Derived here: k = 2.6e-3 min⁻¹ at 145 °C and 4.9e-3 min⁻¹ at 160 °C, so
  about 0.5 % of the thiamine would go to HMP per minute at 160 °C. Even at one minute that is about
  0.7-1 mM HMP from this dose, two to three orders of magnitude above the whole beef thiamine pool. The
  engine has no water-activity term for a 28-34 % d.b. melt; that is outside its validated domain.

**What it gives the model: nothing quantitative.** No concentration, no time course, no residence time,
one process temperature profile, one dose. It is a sensory direction check only, and the thiamine row's
probability in section 10 (about 0.3 of surviving the reference-pot test) is unchanged by it: the row's
open question is magnitude at 100 °C and beef dose, which this paper does not touch.

**For a formulator.** Adding thiamine HCl at 1.5 % of SPC before single-screw extrusion (145/160 °C last
zones) is a cheap, lipid-free way to put a mild meat/chicken odour into TSP and cut soy odour, with no
measurable cost to texture, yield or shrinkage. The gain is small on a 0-3 RATA scale (about +0.3 meat
odour, −0.2 to −0.5 soy odour) and shows up as odour liking; it reaches flavour and overall liking only
where the product carries enough fat (burger, 9.7 % oil). Burnt aftertaste is the side effect.

## What it does not give

- Any volatile measurement (no MFT, no thiazoles, no hexanal): the paper says so.
- Residence time, melt temperature, SME, pressure; the pH of SPC, melt or product.
- The basis of "1.5 % (w/w)" (dry or moistened SPC; salt or free base).
- Thiamine remaining in the TSP after extrusion, or in the SPC before it.
- A dose-response: one dose, one control.
- The storage time between extrusion and sensory testing.
- A separation of thiamine's effect from moisture's at 28.4 % M (no 28.4 % M control was extruded).
- Any beef-level or cultivated-tissue-level condition.
