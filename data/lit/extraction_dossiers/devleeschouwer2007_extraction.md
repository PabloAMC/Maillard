# De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2007 — EXTRACTION (the a_w 0.34-0.92 table)
### Acrylamide formation/elimination in low-moisture equimolar asparagine-glucose at six water activities, 120-200 C.

**Source on disk:** `data/articles/devleeschouwer2008b.pdf` (owner's download, 2026-09-07; the file is
named 2008b but the paper is Biotechnol. Prog. 2007). Read-only extraction, wave B15.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Kinetics of Acrylamide Formation/Elimination Reactions as Affected by Water Activity" |
| Authors | Kristel De Vleeschouwer, Iesel Van der Plancken, Ann Van Loey, Marc E. Hendrickx |
| Venue | Biotechnol. Prog. 2007, 23, 722-728 |
| DOI | 10.1021/bp060389f |

## 1. Why the repo needed it

B14 declared the acrylamide lane's a_w term flat inside 0.88-0.99 from the 2008 paper and refused
everything drier. This paper is the same group's earlier series over a_w 0.34 / 0.59 / 0.73 / 0.82 /
0.88 / 0.92 (equilibrated at 4 C over saturated salts; Table 1 of the paper lists the salts), so it
extends the measured window down to 0.34 and adds the elimination and Maillard-competition shapes.

## 2. Methods

Equimolar Asn-Glc freeze-dried powders, dried over P2O5 then equilibrated over saturated salts;
sorption isotherm measured; heated in closed reactor tubes in an oil bath at 120 / 140 / 160 / 180 /
200 C; the model was extended (Scheme 2) with an apparent second-order Maillard competition k_M
because Scheme 1 lacked fit. T_ref 160 C. Values +/- SE.

## 3. Table 2 — Scheme 2 parameters per initial a_w

| a_w | k_Fref (1e-3 /M/min) | k_Eref (1e-3 /min) | k_Mref (1e-3 /M/min) | Ea_F | Ea_E | Ea_M | pseudo-R2 S2 / S1 |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0.34 | 0.647 +/- 0.168 | 371 +/- 333 | 53.6 +/- 52.7 | 153.6 +/- 24.95 | 116.1 +/- 41.2 | 86.50 +/- 64.34 | 0.871 / 0.776 |
| 0.59 | 0.736 +/- 0.110 | 291 +/- 96.4 | 83.3 +/- 44.6 | 161.0 +/- 24.17 | 44.34 +/- 17.21 | 194.2 +/- 22.16 | 0.943 / 0.749 |
| 0.73 | 0.675 +/- 0.138 | 162 +/- 36.5 | 160 +/- 59.9 | 160.8 +/- 11.96 | 60.88 +/- 13.00 | 157.1 +/- 17.53 | 0.967 / 0.801 |
| 0.82 | 0.996 +/- 0.259 | 178 +/- 41.4 | 250 +/- 107 | 170.8 +/- 15.47 | 56.66 +/- 15.14 | 179.8 +/- 20.90 | 0.942 / 0.909 |
| 0.88 | 0.703 +/- 0.058 | 321 +/- 161 | 57.6 +/- 37.3 | 161.9 +/- 7.368 | 97.80 +/- 19.68 | 127.9 +/- 25.95 | 0.966 / 0.773 |
| 0.92 | 0.913 +/- 0.082 | 485 +/- 186 | 58.6 +/- 24.4 | 182.2 +/- 12.05 | 152.2 +/- 18.54 | 104.1 +/- 34.88 | 0.954 / 0.737 |

Authors' reading: k_F "varies only slightly" with a_w (no significant trend); k_M has a MAXIMUM at
a_w 0.82 (the classic Maillard optimum, ~5 % moisture on their isotherm); k_E has a MINIMUM at 0.82.
Relative to the 0.92 column: k_F 0.71 / 0.81 / 0.74 / 1.09 / 0.77 / 1.00; k_E 0.76 / 0.60 / 0.33 /
0.37 / 0.66 / 1.00 (the 0.34 value carries an SE of 90 %).

## 4. What the repo takes (wave B15) and what it does not

- TAKEN: the formation window of B14's flat term extended from 0.88-0.99 down to **0.34**, the band
  widened to the union of both papers' point-estimate spreads (0.41-1.39 already covers 0.71-1.09).
- TAKEN: a DECLARED elimination multiplier on `k_acr_dp` through the dry-side points (0.76 / 0.60 /
  0.33 / 0.37 at a_w 0.34 / 0.59 / 0.73 / 0.82, relative to 0.92), joining 1.0 at a_w 0.88 where the 2008
  series measured the constants flat (this paper's 0.88 point, 0.66 +/- 0.33, is within one SE of 1);
  the envelope scales the deficit by 0-1.2 (from no effect to 1.2x the shape); 1.0 at a_w None.
- NOT taken: k_M (the lane's competition is mass-action, not a multiplier); the Ea-vs-a_w shapes (the
  lane keeps one Ea per step; recorded as a caveat that the elimination barrier spans 44-152 kJ/mol).
- Caveat: these are the 2007 Scheme-2 constants; the lane's shipped constants are the 2008/2009
  two-step scheme. Only within-study ratios transfer.
