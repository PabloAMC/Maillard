# Mulley, Stumbo & Hunting 1975 — EXTRACTION (thermal destruction of thiamine hydrochloride and co-carboxylase at 265 °F, pH 4.5 to 6.5)

**Source on disk:** `data/articles/mulley1975.pdf` (a 4-page scan with an OCR text layer, pp. 989-992).
Read and checked by eye on 2026-10-09: every number below was read from the page image (Table 1 also
from a 300 dpi crop) and cross-checked against `pdftotext -layout`. The OCR layer gets one methods number
wrong ("2.5 ml" where the page prints "25 ml", p. 989), and it misses the figures, so the image is the
authority. Written for the thiamine route of the sulfur lane (`src/kinetic_core/sulfur.py`, "THE THIAMINE
ROUTE").

| field | value |
|---|---|
| Title | "Kinetics of thiamine degradation by heat. Effect of pH and form of the vitamin on its rate of destruction" (p. 989) |
| Authors | E. A. Mulley, C. R. Stumbo, W. M. Hunting (Dept. of Food Science & Nutrition, University of Massachusetts, Amherst) |
| Venue | Journal of Food Science, Volume 40 (1975), pp. 989-992 (running heads, pp. 989-992). Issue number not printed. Ms received 5/31/74, revised 3/19/75, accepted 3/22/75 (p. 992). |
| DOI | **Not printed on the PDF.** 10.1111/j.1365-2621.1975.tb02250.x is the DOI as supplied with the task; the page range it should point to (989-992) is confirmed, the DOI string itself is not verifiable from this file. |

The paper says the thermoresistometer modifications were "described earlier" (p. 989); that is the
companion paper (J. Food Sci. 40:985), which is **not on disk**. Anything about temperature dependence
would have been there or in Feliciotti's 1955 thesis (cited p. 992, "temperature range 228° to 300°F");
neither is read here.

## 1. Methods (pp. 989)

- **Forms of the vitamin:** synthetic thiamine hydrochloride ("T-HCl") and co-carboxylase (thiamine
  pyrophosphate, "CoCar"), and two mixtures. Nothing else: **no mononitrate, no thiamine monophosphate,
  no protein-bound or food-matrix thiamine** in this paper.
- **Stocks:** 37.5 mg of each chemical in 25 ml of 25 % ethanol, after ≥ 24 h drying over P₂O₅.
- **Working solutions** (each diluted to 25 ml with buffer of the desired pH):
  100 % T-HCl = 1 ml T-HCl stock; 65 % T-HCl (35 % CoCar) = 0.65 ml T-HCl + 0.35 ml CoCar stock;
  30 % T-HCl (70 % CoCar) = 0.3 ml + 0.7 ml; 0 % T-HCl (100 % CoCar) = 1 ml CoCar stock.
  *Derived here:* 37.5 mg / 25 ml = 1.5 mg/ml stock; ÷ 25 = 60 µg/ml working; × 20 µl = 1.2 µg per
  sample, which matches the ~1 µg starting ordinate of Figs. 1-5. For pure T-HCl (337.26 g/mol) that is
  60 / 337.26 = **0.178 mM**.
- **Buffer:** phosphate, 1/10 M, pH 4.5, 5.0, 5.5, 6.0, 6.5. "The effect at hydrogen-ion concentrations
  of pH 7.0 and above was not investigated" (p. 989).
- **Heating:** Stumbo's thermoresistometer, triplicates of 20 µl, **one temperature only: 265 °F**
  (*derived here:* (265 − 32) × 5/9 = **129.4 °C**). Unheated controls handled alike.
- **Assay:** sample into 2 ml 0.1 N HCl, enzyme hydrolysis ≥ 3 h at 113-122 °F, thiochrome fluorescence
  (AVC 1966). The enzyme step dephosphorylates co-carboxylase, so the assay reads **total remaining
  intact thiamine** (free + phosphorylated). **Loss is total thiamine destruction, not a named product
  channel.**
- **Data treatment:** log(concentration) vs time at 265 °F, lines by linear regression; "At every pH,
  first order rates of reaction were observed" for both forms and the mixtures (p. 989). D values from
  these lines.

## 2. Findings that matter

### Table 1 (p. 991): D values, minutes, phosphate buffer, 265 °F

**How the header must be read.** The table prints a top percentage row (0 %, 30 %, 65 %, 100 %) and a
second one (100 %, 70 %, 35 %, 0 %), with the row labels "T-HCl" and "CoCar" stacked at the left beside
the second and third header lines. Read literally ("T-HCl 100 %" over the first column) the table would
contradict the paper. The top row is the **T-HCl percentage** and the second the **CoCar percentage**,
on four independent pieces of evidence: (i) the methods' mixtures are 65 % T-HCl / 35 % CoCar and 30 %
T-HCl / 70 % CoCar (p. 989), which only the top-row = T-HCl reading reproduces; (ii) the text says
co-carboxylase "is destroyed more rapidly than thiamine hydrochloride" (p. 992), and the first column has
the smallest D at every pH; (iii) Fig. 7, whose legend labels curves by % THCl, puts curve 4 (0 % THCl)
lowest, at ≈ 11 min at pH 6.5, matching the first column's 11.6; (iv) the text says up to 35 % CoCar
does not change the rate (p. 991), and the 65 % T-HCl column tracks the 100 % T-HCl column. Columns
below are labelled accordingly.

| pH | 100 % CoCar (0 % T-HCl) | 70 % CoCar (30 % T-HCl) | 35 % CoCar (65 % T-HCl) | 100 % T-HCl (0 % CoCar) |
|---|---|---|---|---|
| 4.5 | 88.6 | 100.1 | 100.8 | 99.5 |
| 5.0 | 97.6 | 102.2 | 114.3 | 107.3 |
| 5.5 | 64.4 | 76.6 | 97.0 | 107.7 |
| 6.0 | 27.7 | 46.8 | 73.6 | 76.0 |
| 6.5 | 11.6 | 24.6 | 32.8 | 36.0 |

### First-order rate constants at 129.4 °C — derived here, k = ln 10 / D, min⁻¹

| pH | 100 % CoCar | 70 % CoCar | 35 % CoCar | 100 % T-HCl |
|---|---|---|---|---|
| 4.5 | 0.0260 | 0.0230 | 0.0228 | 0.0231 |
| 5.0 | 0.0236 | 0.0225 | 0.0202 | 0.0215 |
| 5.5 | 0.0358 | 0.0301 | 0.0237 | 0.0214 |
| 6.0 | 0.0831 | 0.0492 | 0.0313 | 0.0303 |
| 6.5 | 0.198 | 0.0936 | 0.0702 | 0.0640 |

(e.g. 2.302585 / 76.0 = 0.0303 min⁻¹; half-life = D × log₁₀2 = 0.301 D, so 22.9 min for T-HCl at pH 6.0.)

### Fig. 7 (p. 991) cross-check, D vs pH on a log axis — read from graph, approx.

Curve 1 (100 % THCl): ≈ 106, 112, 98, 73, 35 min at pH 4.5, 5.0, 5.5, 6.0, 6.5. Curve 4 (0 % THCl):
≈ 83, 89, 59, 29, 11 min. Curve 3 at pH 6.0 ≈ 46, at 6.5 ≈ 22; curve 2 at 6.5 ≈ 29. Agreement with
Table 1 is within ≈ 10 % throughout (largest gap: curve 1 at pH 5.5, ≈ 98 against 107.7), which
confirms the column reading above. The authors read the log D vs pH plot as two straight lines with a
break that sits at lower pH for co-carboxylase than for thiamine hydrochloride; the break is drawn,
not tabulated.

### What the text says, without numbers

- D is flat from pH 4.5 to 5.0 and falls sharply above pH 6.0 (p. 991); the authors tie this, after
  Feliciotti 1955, to the hydrochloride being neutralised near pH 6.2.
- Up to 35 % co-carboxylase in a mixture does not change the rate over pH 4.5-6.5; above that, the
  mixture degrades faster (pp. 991-992).
- Mechanism cited, not measured here: cleavage of the C-N bond of the methylene bridge, giving a
  pyrimidine (probably 2-methyl-4-amino-5-hydroxymethylpyrimidine) and 4-methyl-5-(β-hydroxyethyl)
  thiazole (p. 991, after Dwivedi & Arnold 1972). The products were not measured in this paper.

## 3. What it means for the model

1. **There is no Ea, no z value and no second temperature in this paper.** Every run is at 265 °F. The
   paper's own parameters therefore give k at 129.4 °C only; **k at 100 °C and 140 °C cannot be derived
   from this paper.** The borrowed band of the thiol-assembly family (55-145 kJ/mol, centre 100,
   `parameters_sulfur.FORMATION_EA_BOUNDS_BY_ROUTE`) is untouched by it.
2. **What k(100) and k(140) would be if the model's band is borrowed** (derived here, NOT the paper's
   parameters: k(T) = k(402.59 K) · exp[−Ea/R · (1/T − 1/402.59 K)], R = 8.3145 J mol⁻¹ K⁻¹):

   | system, pH | Ea, kJ/mol | k(100 °C), min⁻¹ | k(140 °C), min⁻¹ | lost in 20 min at 100 °C | lost in 5 min at 140 °C |
   |---|---|---|---|---|---|
   | 100 % T-HCl, 6.0 | 55 | 0.0083 | 0.046 | 15 % | 21 % |
   | 100 % T-HCl, 6.0 | 100 | 0.0029 | 0.065 | 5.6 % | 28 % |
   | 100 % T-HCl, 6.0 | 145 | 0.0010 | 0.092 | 2.0 % | 37 % |
   | 100 % T-HCl, 5.5 | 100 | 0.0020 | 0.046 | 4.0 % | 21 % |
   | 100 % T-HCl, 6.5 | 100 | 0.0061 | 0.137 | 11 % | 50 % |
   | 100 % CoCar, 6.0 | 100 | 0.0079 | 0.178 | 15 % | 59 % |

   "Lost" = 1 − exp(−k t), isothermal. **Consequence for the sweep:** for any Ea above 44.4 kJ/mol
   (derived here: Ea = R ln 4 / (1/373.15 − 1/413.15 K) = 44.4 kJ/mol, the barrier at which 20 min at
   100 °C and 5 min at 140 °C destroy the same fraction) — i.e. the model's whole band — 140 °C / 5 min
   destroys MORE thiamine than 100 °C / 20 min, isothermally. So the sweep's "large thiamine boost at
   100 °C/20 min, small at 140 °C/5 min" cannot come from how much thiamine is consumed; it has to come
   from downstream (HMP → MFT vs. MFT sinks) or from the heating ramp. This is an inference, not a
   finding of the paper.
3. **A ceiling, not a value, for `k_thi_hmp`.** The model's thiamine has two fates, `r_thi_hmp` and
   `r_thi_mesh`; its total first-order loss is k_thi_hmp + k_thi_mesh. Mulley's k is total destruction of
   free thiamine in 0.1 M phosphate, with bridge cleavage to pyrimidine + intact thiazole (no HMP) as the
   cited dominant route, so it is an **upper bound** on the sum at 129.4 °C and the matching pH, not a
   measurement of the HMP branch.
4. **pH is a factor of ≈ 3 across food pH** for the free vitamin (k 0.021 → 0.064 min⁻¹ from pH 5.5 to
   6.5, T-HCl) and ≈ 8 for co-carboxylase (0.0236 → 0.198 from pH 5.0 to 6.5). A grep of `src/` found no
   pH factor keyed to `k_thi_hmp` by name; whether the route is pH-dependent through another mechanism
   was not checked.
5. **Form matters only at high pyrophosphate share.** At ≥ 70 % co-carboxylase the rate is 1.4-3.1x the
   T-HCl rate at pH 5.5-6.5 (columns 1-2 vs 4, derived here from the k table). At pH 4.5-5.0 form does
   not matter (all columns within 0.020-0.026 min⁻¹). The phosphorylated share of thiamine in meat is not
   in this paper.

## What it does not give

- No activation energy, no z value, no second temperature, no Arrhenius or TDT plot. k(100 °C) and
  k(140 °C) are not obtainable from this paper alone.
- No data above pH 6.5 (not investigated), none below 4.5 (pH 3.5 is cited from Molitor & Sampson 1936
  only).
- No thiamine mononitrate, monophosphate, protein-bound or in-food thiamine; synthetic vitamins in
  buffer only, at 0.178 mM (derived) — no matrix effect.
- No products measured: nothing on HMP, MFT, 3-mercapto-2-pentanone or the thiazole; loss is total
  thiamine.
- No standard errors or confidence intervals on the D values; triplicates are stated but their spread is
  not printed. No regression statistics.
- No buffer-concentration series (the buffer-salt effect is cited from Farrer, not measured).
