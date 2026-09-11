# Unlu & Faller 2002 — EXTRACTION (measured residence time distribution in a twin-screw extruder at 25 % moisture, with the exact screw speed and feed rate the model's acrylamide row runs at)

**Source on disk:** `data/articles/unlu2002.pdf` (1.3 MB; downloaded 2026-09-11 at this repository's
request). Read 2026-09-11 via `pdftotext -layout`. Wave B45.

| field | value |
|---|---|
| Title | "RTD in twin-screw food extrusion" |
| Venue | Journal of Food Engineering 53 (2002) 115–131 |
| DOI | 10.1016/S0260-8774(01)00148-0 (printed PII S0260-8774(01)00148-0) |
| Extruder | **WP ZSK-30 co-rotating, intermeshing twin-screw**; barrel length 116 cm, diameter 30 mm, **L/D 38.7:1**; free barrel volume 447 cm³ per metre |
| Barrel profile | **five zones at 40 / 80 / 100 / 120 / 140 °C** from feed to die |
| Feed | **degermed yellow cornmeal**, initial moisture 11.24 % (w/w), **total moisture adjusted to 25 % wet basis** by pumping distilled water at the feed inlet |
| Design | 3 × 3 factorial, two replications: feed rate **8.55 / 14.3 / 20.0 kg/h** × screw speed **150 / 250 / 350 rpm**, randomised |
| Method | **KCl tracer** (10 g KCl + 5 g cornmeal dropped at the feed inlet), conductivity measured at flush-mounted electrodes in the die every 5 s; the residence time reported as a **geometric mean** (GMRT) of the log-transformed distribution, with upper and lower limits at one standard deviation |

## 1. Why this paper was fetched

`data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json` — the model's only extrusion acrylamide
row, and a FIT row — declares a hold of **25 s**. Wave B36 read its source (Ma et al. 2024) and found
that **the paper prints no residence time at all**; the 25 s has been carried since as a labelled
assumption. Wave B37 found Yu et al.'s measured table but on a 20 %-protein corn-flour feed. This
paper measures the residence time at **the same screw speed and almost exactly the same feed rate as
Ma's pot, on a machine of almost the same length-to-diameter ratio.**

## 2. Table 5 — the full factorial, verbatim (the cell that matters is the first row)

"LSD analysis results for the RTD parameters and the barrel fill with respect to SFL"

| screw speed (rpm) | feed rate (kg/h) | SFL (kg/(h·rpm)) | **GMRT (s)** | UL (s) | LL (s) | spread (s) | normalised spread | barrel fill (%) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **150** | **8.55** | 0.057 | **123** | 173.5 | 87.2 | 86.3 | 0.70 | 55.5 |
| 250 | 8.55 | 0.034 | 106 | 152.8 | 73.5 | 79.3 | 0.75 | 47.5 |
| 350 | 8.55 | 0.024 | 94.4 | 140.6 | 63.4 | 77.2 | 0.82 | 42.5 |
| 150 | 14.3 | 0.095 | 97.6 | 134.2 | 71.0 | 63.2 | 0.65 | 73.5 |
| 250 | 14.3 | 0.057 | 77.2 | 112.6 | 52.9 | 59.7 | 0.77 | 58.0 |
| 350 | 14.3 | 0.041 | 70.2 | 103.8 | 47.5 | 56.3 | 0.80 | 53.0 |
| 150 | 20.0 | 0.133 | 71.1 | 97.8 | 51.8 | 47.7 | 0.65 | 75.0 |
| 250 | 20.0 | 0.080 | 64.7 | 92.8 | 45.1 | 46.2 | 0.74 | 68.0 |
| 350 | 20.0 | 0.057 | 56.2 | 83.9 | 37.7 | 46.0 | 0.82 | 59.5 |

Marginal means: by feed rate (Table 3) 108 / 81.7 / 64.0 s at 8.55 / 14.3 / 20.0 kg/h; by screw speed
(Table 4) **97.3 / 82.6 / 73.6 s at 150 / 250 / 350 rpm**, with barrel fill 68.0 / 57.8 / 51.7 %.

Verbatim on the ordering: feed rate "had a greater effect on the mean residence time than the screw
speed, with a 41 % reduction in GMRT (108–64.0 s) with a 2.33-fold" increase in feed rate.

## 3. The comparison this repository came for

| | Unlu & Faller 2002 | Ma et al. 2024, the model's row |
|---|---|---|
| screw speed | **150 rpm** | **150 rpm** |
| total feed | **8.55 kg/h** | **8.57 kg/h** (6.0 raw + 2.57 water) |
| L/D | **38.7:1** | **40:1** |
| moisture | 25 % | 30 % |
| feed material | degermed yellow cornmeal | **soy protein isolate : corn starch 9:1** |
| zones | 5, to 140 °C | 10, the 130 °C arm at 80/80/85/90/100/110/120/130/130/130 °C |
| **whole-barrel residence** | **123 s (87.2–173.5)** | **not printed anywhere in the paper** |

**The arithmetic this licenses, and its assumption stated.** Three of Ma's ten zones sit at 130 °C.
If residence time is distributed in proportion to zone length, the material spends about **3/10 of
123 s ≈ 37 s at 130 °C**. The bundle assumes **25 s**. So the unsourced assumption is **about 1.5×
low, not fourfold low** — which is the opposite of what the reading list expected when it asked for
this measurement, and it is the reason the number is recorded here rather than changed.

**Four things that make it an estimate and not a transfer.** (i) The feed is a starch meal, not a
90 %-protein melt, and melt viscosity is the first-order determinant of barrel fill; (ii) 25 % against
30 % moisture, and this paper's own result is that higher moisture shortens residence; (iii) equal
zone lengths is my assumption, not the paper's; (iv) the GMRT spans 87–174 s at one standard
deviation, so even the central figure carries a factor of two.

## 4. What else is printed

Table 7 gives stepwise regressions for barrel fill, die pressure, motor torque, product temperature
and specific mechanical energy against feed rate and screw speed (R² 0.95–0.99), and Table 8 the same
against the specific feeding load. Product temperature rises with screw speed and falls with feed
rate. Those are process relations for a cornmeal extrudate and are not taken.

## 5. Verdict

The best-matched extrusion residence time on disk for the model's acrylamide row: same screw speed,
same feed rate to within 0.2 %, nearly the same machine ratio — and a different feed material. It
does not license replacing the bundle's 25 s, because the melt is not the same melt. It does
establish that the assumption is of the right order, which retires the louder half of a named debt
and replaces it with a smaller one.
