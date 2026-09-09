# Zhang et al. 2020 — EXTRACTION (α-dicarbonyl kinetics in aqueous glucose and glucose-glutamate)
### Glucosone, 1-DG, 3-DG, 3,4-DDG, glyoxal, methylglyoxal and diacetyl at 90-110 C, 0-6 h, in water.

**Source on disk:** `data/articles/zhang2020.pdf` (owner's download, 2026-09-07; open access CC-BY).
Read-only extraction, B13 validation search reopened.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Kinetics of α-dicarbonyl compounds formation in glucose-glutamic acid model of Maillard reaction" |
| Venue | Food Science & Nutrition 2021, 9, 290-302 (received 13 Sep 2020, accepted 25 Oct 2020) |
| DOI | 10.1002/fsn3.1995 |

## 1. Methods

- Glucose-only (0.3 M) and glucose + L-glutamic acid (0.3 M each) **in water** (no buffer; pH NOT
  stated anywhere), in 250 mL four-neck flasks with reflux cooler, silicone oil bath, atmospheric
  pressure, 90 / 95 / 100 / 105 / 110 C for 0-6 h (hourly), triplicates.
- Dicarbonyls by o-phenylenediamine derivatisation (Kocadagli & Gokmen 2016 method: 0.5 mL reaction
  mixture + phosphate buffer pH 7 + o-PDA/DETAPAC) and LC-ESI-MS/MS of the quinoxalines; glucose and
  fructose by HPLC-RID; HMF by HPLC-DAD. Calibration 3.2e-3 to 2 ug/mL for the quinoxalines.
- Kinetics: linear, exponential or logarithmic fits of C(t) per compound and temperature (Tables 1-2);
  NO activation energies reported.

## 2. Tables 1 and 2 — the dicarbonyl rows (C in ug/mL of reaction mixture, t in h)

Glucose-only (Table 1):

| compound | 90 C | 95 C | 100 C | 105 C | 110 C |
|---|---|---|---|---|---|
| glucosone | 0.0003 t + 0.0005 | 0.0004 t + 0.0015 | 0.0002 t + 0.0015 | 0.0005 t + 0.0017 (3-6 h) | 0.0001 t + 0.0014 (3-6 h) |
| glyoxal | 5e-5 t + 0.0003 | 6e-5 t + 0.0002 | 4e-5 t + 0.0002 | 4e-5 t + 0.0001 | 1e-5 t + 0.0002 |
| methylglyoxal | 2e-5 t + 1e-5 | 2e-5 t + 2e-5 | 5e-5 t + 2e-5 | 1e-4 t + 4e-5 | 5e-5 t + 4e-5 |
| diacetyl | 1e-6 t + 1e-6 | 1e-6 t + 4e-6 | 2e-6 t + 4e-6 | 3e-6 t + 8e-6 | (row cut in the text layer) |
| 3-DG | 0.0028 t + 0.0049 | 0.0027 t + 0.0063 | 0.0041 t + 0.0119 | 0.0056 t + 0.0173 | 0.0054 t + 0.017 |
| 3,4-DDG | 0.0005 t + 0.0016 | 0.0006 t + 0.0016 | 0.0010 t + 0.0028 | 0.0018 t + 0.0024 | 0.0015 t + 0.0047 |
| 1-DG | - | 3e-5 ln t - 3e-5 | 7e-5 ln t + 3e-5 | 2e-4 ln t + 1e-4 | 1e-4 ln t + 4e-5 |

Glucose-Glu (Table 2):

| compound | 90 C | 95 C | 100 C | 105 C | 110 C |
|---|---|---|---|---|---|
| glucosone | 0.0004 t + 0.0009 | 0.0004 t + 0.0016 | 0.0003 t + 0.001 | 0.0003 t + 0.0007 | 0.0004 t + 0.0008 |
| glyoxal | 0.0001 t + 6e-5 | ln C = 0.569 t + ln 0.0002 (1-4 h) | ln C = 0.2617 t + ln 0.0004 (1-4 h) | ln C = 0.1225 t + ln 0.0003 (1-4 h) | ln C = 0.1671 t + ln 0.0004 |
| methylglyoxal | 0.0002 t - 1e-4 | 0.0002 t + 1e-4 | 0.0001 t + 1e-4 (1-4 h) | 0.0001 t - 1e-4 | 0.0002 t - 2e-4 |
| diacetyl | ln C = 0.3968 t + ln 8e-6 | ln C = 0.8168 t + ln 6e-6 (1-5 h) | 3e-5 t + 2e-5 | 5e-5 t - 4e-5 | 1e-4 t - 1e-4 |
| 3-DG | 0.0013 t + 0.0037 | 0.0011 t + 0.0061 | 0.0018 t + 0.0094 | 0.0029 t + 0.011 | 0.006 t + 0.008 |
| 3,4-DDG | 0.0002 t + 0.0007 | 0.0002 t + 0.001 | 0.0003 t + 0.0014 | 0.0006 t + 0.0014 | ln C = 0.2785 t + ln 0.0017 |
| glucose | -353 t + 54103 (90 C, glucose-only) ... 4-19 % lost by 6 h | | | | |

## 3. ⚠ The absolute levels do not mass-balance

Glucose falls by 4-19 % of 0.3 M (12-58 mM, 2-10 g/L) while the summed dicarbonyls at 6 h reach
~0.05 ug/mL (3-DG) — nine orders of magnitude below the sugar lost, and three to four below the
mg/L-scale 3-DG that Kocadagli & Gokmen 2016 and Gokmen's aqueous systems report. Either the printed
unit is not ug/mL of reaction mixture (a dilution or an extract-basis unit is likely) or the calibration
is off by a large factor. **The repo therefore uses this paper for ORDERINGS and TIME SHAPES only, never
for a level.** Orderings that survive any common factor: in glucose-only at 105 C / 6 h, 3-DG (0.051)
> 3,4-DDG (0.013) > glucosone (0.0047) > methylglyoxal (0.00064) > glyoxal (0.00034) > diacetyl
(0.000026); in glucose-Glu at 100 C / 6 h, 3-DG (0.020) > glucosone (0.0028) > glyoxal (0.0019) >
methylglyoxal (0.0007) > diacetyl (0.0002).

## 4. What the repo takes

- Directional claims DIC-01 (glucose-only ordering at 105 C / 6 h: 3-DG > glucosone > methylglyoxal >
  glyoxal > diacetyl) and DIC-02 (glucose-Glu ordering at 100 C / 6 h). The glucose-only pot is the
  trunk's amine-free caramelisation set (B7 constants) and is chargeable; the glutamate pot has no
  glutamate species in the core and is recorded as not evaluable.
- pH: unstated. For the glucose-only pot the trunk's pH term touches only Amadori decay, which is
  absent without an amine, so the declared pH changes nothing there and 6.5 (a fresh glucose solution)
  is recorded as the run pH with that note.
