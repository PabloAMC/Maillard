# Pre-registration: wave B40, the 3-deoxyglucosone exits get their pH term, and the fed fit is re-run under it (written 2026-09-11, BEFORE any run)

## 1. Why

B39 fitted the fed 3-deoxyglucosone triangle to 0.09 dex on every maximum and pinned all five
coordinates, and did not ship, because the fed pot's 3,4-DGE still peaked at 12.5 min against a
20–60 min bracket. The timing is set by how fast 3-DG leaves, and its dominant exit — Martins'
formic-acid step, `k_tdg_fa` — was measured at pH 6.8 and is applied at pH 5 with no correction.
The same paper prints the correction.

## 2. The source, already on disk, and the declaration

Martins & van Boekel 2003, Part II, Table 3 (`martins2003_extraction.md` §5.2): rate constants of
the glucose–glycine network at 100 and 120 °C, pH 5.5 and 6.8, with 95 % HPD.

| step | 100 °C pH 5.5 | 100 °C pH 6.8 | 120 °C pH 5.5 | 120 °C pH 6.8 |
|---|---:|---:|---:|---:|
| k6, 3-DG → formic acid (the trunk's `k_tdg_fa`) | 1.9 × 10⁻³ | 2.74 × 10⁻² | 4.30 × 10⁻² | 3.04 × 10⁻¹ |
| k5, 3-DG → fragments (the trunk's `k_tdg_mgo`, Kocadagli's amine-free route, is the nearest step) | 1.38 × 10⁻² | 9.07 × 10⁻² | 2.23 × 10⁻² | 5.07 (flagged by the dossier as a fitting artefact) |

Slopes in decades per pH unit over the 1.3-unit interval: **k6: 0.89 at 100 °C, 0.65 at 120 °C;
k5: 0.63 at 100 °C** (its 120 °C pair is not used, per the dossier's flag). B12 declared the Amadori
steps' term from the same table at 0.69 with band (0.37, 0.92) and reference pH 6.8.

**Declared here, not fitted:** `k_tdg_fa` scales by 10^(0.77·(pH − 6.8)), band (0.65, 0.89), the
mean and range of k6's two temperatures; `k_tdg_mgo` scales by 10^(0.63·(pH − 6.8)), band (0.37,
0.92) — k5's one usable slope inside B12's declared band. Measured window pH 5.5–6.8, as B12; outside
it the engine says so, as it already does for the Amadori term. The term applies on the trunk lane
wherever pH is declared, so it touches every pot below 6.8: that is the point, and it is what the
ship rule must judge.

## 3. The fit

B39's generator, rows, sigmas, bounds and starts, unchanged, run with the term on. The Leitzen 2021
hold-out is never read.

## 4. Predictions

- **P1.** The fed rows fit as in B39 and the fed-3-DG peak now lands **inside 20–60 min**; χ²_red
  below 3.
- **P2.** `k_tdg_ddg` moves **less** than B39's +0.57 dex, because the slower exits do part of the
  work the forward rate was doing.
- **P3 — the risk, named.** Leitzen's autoclaved glucose is at pH 4.36, where the term slows the 3-DG
  exits by roughly 70×. Its 3-DG row, 1.11× today, **rises**. I predict it ends **between 2× and 6×**;
  the ship rule needs it at or under 3×, so this prediction is compatible with either verdict and is
  written so that the number cannot be re-read after the fact.
- **P4.** Leitzen's 3,4-DGE improves from 32× to under 10× and HMF does not worsen beyond 3 % of
  its current 11.9× — the tolerance B39's rule did not have and the reason B39's P4 failed by 2 %;
  declared here before the run, not after.
- **P5.** Martins' own pH 6.8 rows (B1's fit corpus) are untouched, since 6.8 is the reference.
- **P6.** No row now within 3× leaves the band; the two rows B39's candidate brought in stay in.
- **P7.** Every other trunk-lane pot below pH 6.8 (Schibilsky pH 5, Hofmann pH 3 and 5, Yiltirak
  5.5, the acrylamide pots at 6.0) is reported row by row; none is a test, because none measures a
  3-DG exit.

## 5. Ship rule

SHIP if P1, P4 and P6 hold **and** Leitzen's 3-DG is at or under 3× (the decisive half of P3). If
the 3-DG row alone fails, the record says: a pH term measured on a glycine pot at pH 5.5–6.8 does
not transfer to an amine-free pot at pH 4.4, and the term and the triangle both stay uninstalled.
Frozen pair under `_b40_baseline/`.

If it ships: the term's flag defaults on, B39's frozen literals become the B40 optimum
(`SHIPPED_B39` True), the envelope gains the five `b39.` rows from the fit's Laplace σ and drops
the ENV-B34 printed band on `k_tdg_ddg` (superseded by data), and the panel headline is re-pinned
with the reason.

## 6. Outcome (written 2026-09-11, after the run) — DO NOT SHIP, and the hold-out said exactly which declaration was wrong

Artifacts: `kinetic_core_b40_fit_report.*`, `kinetic_core_b40_ship_rule.*`, frozen pair under `_b40_baseline/`.

| prediction | result |
|---|---|
| P1 fed rows fit, peak inside 20–60 min | **HELD**: peak at **25.5 min**, cost 1.04 on twelve rows, χ²_red 0.15 |
| P2 `k_tdg_ddg` moves less than B39's +0.57 | HELD: +0.54 dex |
| P3 Leitzen 3-DG ends between 2× and 6× | **REFUTED, in the model's favour**: 1.11× → **1.36×**. The 70× slowing of the exits at pH 4.36 did not pile 3-DG up, because the reversible triangle drains it |
| P4 Leitzen 3,4-DGE < 10×, HMF within 3 % | HELD: 32.4× → **6.0×**; HMF 11.9× → **7.2×** |
| P5 Martins' pH 6.8 rows untouched | HELD (reference pH) |
| P6 no row leaves the band | **REFUTED**: Leitzen's **methylglyoxal 1.28× → 33.2×** left it; within-3× net 10 → 11 (Schibilsky pH-8 HMF 3.6× → 1.6× and Chang-water HMF 3.2× → 2.6× entered) |
| P7 other pots reported | Schibilsky pH-5 HMF 2.05× → 1.07×, Lin 2022 HMF 6.3× → 6.2× |

**The finding.** Every gain came from the formic-acid exit's term — Martins' k6, measured at two pH
values on the very step the trunk carries. The one loss came from the term on `k_tdg_mgo`, which
was declared from Martins' k5, a lumped "3-DG → fragments" step, and *transferred* to Kocadagli's
amine-free methylglyoxal route. At pH 4.36 that transfer starved the methylglyoxal route 35-fold,
and Leitzen's methylglyoxal row — right to 28 % before — went to 33×. The hold-out rejected one of
two declarations and kept the other, by name. **B41 keeps what the data kept.** Nothing installed.
