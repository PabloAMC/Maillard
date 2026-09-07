# De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2006 — EXTRACTION (the pH tables)
### Acrylamide formation/elimination in buffered equimolar asparagine-glucose at pH 4 / 6 / 8, phosphate vs citrate, 120-200 C.

**Source on disk:** `data/articles/devleeschouwer2006.pdf` (owner's download, 2026-09-07). Read-only extraction, wave B15.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Impact of pH on the Kinetics of Acrylamide Formation/Elimination Reactions in Model Systems" |
| Authors | Kristel De Vleeschouwer, Iesel Van der Plancken, Ann Van Loey, Marc E. Hendrickx (KU Leuven) |
| Venue | J. Agric. Food Chem. 2006, 54, 7847-7855 |
| DOI | 10.1021/jf0611264 |

## 1. Why the repo needed it

`parameters_acrylamide.py` policy 3: "Nothing in the Module 3 fit corpus varies pH ... One pH, therefore
no pH term." This paper varies it: the same equimolar system, the same laboratory and the same
two-step (formation second order, elimination first order) fit at three initial pH values in two buffers.

## 2. Methods, as they matter to the lane

- 0.1 M L-asparagine + 0.1 M D-glucose in 0.05 M citrate or phosphate buffer, initial pH 4, 6 or 8
  measured at room temperature after dissolving (unchanged, except citrate pH 8 -> 6.96).
- Heated in closed reactor tubes in an oil bath at 120 / 140 / 160 / 180 / 200 C for variable times;
  two tubes per point; acrylamide by LC-MS/MS after SPE (recovery-checked). Non-isothermal correction
  as in the group's other papers; T_ref 160 C.
- Potato-based variant: 20 % lyophilised potato matrix, phosphate only, 140-200 C.
- **pH drifts during heating (Table 2 below):** at initial pH 8 the pH falls by up to 2.25 units; at
  pH 4 it rises by up to 1.2. Every constant is therefore indexed to the INITIAL room-temperature pH.

## 3. Table 1 — simplified system, T_ref 160 C (value +/- SE; letters = 95 % groups)

| buffer | pH | k_Fref (1e-3 /M/min) | k_Eref (1e-3 /min) | Ea_F (kJ/mol) | Ea_E (kJ/mol) | pseudo-R2 |
|---|---:|---:|---:|---:|---:|---:|
| phosphate | 8 | 37.5 +/- 4.21 a | 333.6 +/- 41.4 a | 130.3 +/- 6.12 a | 84.19 +/- 6.52 a | 0.951 |
| phosphate | 6 | 8.78 +/- 1.34 b | 175.3 +/- 36.1 b | 190.6 +/- 8.18 b | 128.3 +/- 9.54 b | 0.961 |
| phosphate | 4 | 4.30 +/- 0.636 c | 84.2 +/- 18.8 c | 208.1 +/- 8.53 c | 158.5 +/- 10.0 c | 0.958 |
| citrate | 8 (6.96) | 20.1 +/- 1.73 d | 219.7 +/- 25.0 b,d | 146.1 +/- 4.30 d | 115.3 +/- 5.23 b,d | 0.962 |
| citrate | 6 | 20.8 +/- 2.01 d | 262.8 +/- 31.1 d | 159.1 +/- 5.15 e | 106.3 +/- 5.82 d | 0.975 |
| citrate | 4 | 0.599 +/- 0.134 e | 43.1 +/- 13.4 e | 277.4 +/- 13.1 f | 221.1 +/- 14.2 e | 0.946 |

Authors' "log-linear" fit (phosphate): slope **0.5414 +/- 0.106** per pH unit for formation (r 0.9986)
and **0.3442 +/- 0.013** for elimination (r 0.9627) — in NATURAL-log units: ln(37.5/4.30)/4 = 0.54 and
ln(333.6/84.2)/4 = 0.34 reproduce them, log10 would give 0.235 / 0.149. So **0.235 and 0.149 decades per
pH unit**. Not valid for citrate.
Both barriers RISE as pH falls (Ea_F 130 -> 208, Ea_E 84 -> 159 in phosphate).

## 4. Table 3 — potato-based system, phosphate, 140-200 C

| pH | k_Fref (1e-3 /M/min) | k_Eref (1e-3 /min) | Ea_F | Ea_E | pseudo-R2 |
|---:|---:|---:|---:|---:|---:|
| 8 | 40.9 +/- 4.61 | 600.9 +/- 66.90 | 122.7 +/- 6.81 | 83.9 +/- 6.60 | 0.961 |
| 6 | 9.63 +/- 0.965 | 236.4 +/- 31.7 | 189.1 +/- 5.21 | 139.9 +/- 6.41 | 0.977 |
| 4 | 7.29 +/- 0.607 | 154.4 +/- 19.6 | 178.4 +/- 4.91 | 142.0 +/- 5.88 | 0.975 |

Slopes: formation **0.4312 +/- 0.168**, elimination **0.3397 +/- 0.073** per pH unit (ln units; 0.187 / 0.148 decades).

## 5. Table 2 — pH drift after heating (room temperature, vs initial)

| T (C) | phos pH 8 | phos 6 | phos 4 | cit 8 | cit 6 | cit 4 |
|---:|---:|---:|---:|---:|---:|---:|
| 120 | -1.21 | -0.39 | +0.40 | -1.05 | -0.17 | +0.26 |
| 140 | -1.71 | -0.91 | +0.84 | -1.40 | -0.64 | +0.26 |
| 160 | -2.25 | -1.02 | +1.01 | -1.81 | -1.03 | +0.53 |
| 180 | -2.18 | -1.32 | +1.17 | -1.96 | -1.00 | +0.71 |
| 200 | -1.95 | -1.28 | +1.20 | -1.85 | -0.96 | +1.05 |

## 6. What the repo takes (wave B15) and what it does not

- TAKEN: a DECLARED pH factor on the lane's initiation `k_asn_glc` of 10^(0.235 (pH - 6.8)) and on the
  elimination `k_acr_dp` of 10^(0.149 (pH - 6.8)), phosphate slopes converted to decades, bands spanning the
  potato slopes and the SEs (formation 0.114-0.281, elimination 0.116-0.155 decades per unit), window pH 4-8, reference the
  lane's declared 6.8. The absolute k values are NOT transplanted: the 2006 scheme is one-step
  (Asn + Glc -> AA) while the lane's is two-step (De Vleeschouwer 2009), so only the within-study ratios
  transfer. The barrier shifts with pH are recorded, not modelled (the lane keeps one Ea per step).
- NOT taken: citrate (no log-linear relation; the authors attribute the difference to phosphate
  catalysis); the potato matrix (a validation set for a matrix the lane does not represent).
- Caveat: constants are indexed to INITIAL pH while the pot drifts by up to 2 units; the lane has no
  pH trajectory on this lane, so the declared factor is an initial-pH factor by construction.
