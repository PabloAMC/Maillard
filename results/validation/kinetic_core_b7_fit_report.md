# Build Wave B7 — the furanic channel — FIT REPORT

Generated 2026-08-29 at `49df685` on `audit-remediation`.
Pre-registered in `results/validation/kinetic_core_b7_prereg.md`, written before this script existed. Roles: docs/reference/FIT_HOLDOUT_DECLARATION.md Amendments 8 and 12.

## The whole fit, in one line

**One free parameter.** `k_dpo_af` = **4.028960e-06** L/(mmol*min), fitted to **three** cells of Blank 1997 Table 1.

Everything else in the channel is ingested (Kocadagli's seven glucose constants, Hamzalioglu's sink), derived (the hexose transfer), a digitised prior (Edge B's level) or **exactly zero** (Edge C).

### The two seeded starts

| start `log10 k` | converged `log10 k` | objective |
|---|---:|---:|
| -8.0 | -5.394807040 | 24.9936 |
| -4.0 | -5.394807040 | 24.9936 |

Agreement: **0.00e+00** decades (IDENTIFIED; the pre-registration's threshold is 1e-6).

### The three fit rows

| row | measured µg/mmol | predicted µg/mmol | residual (decades) |
|---|---:|---:|---:|
| arabinose/glycine HDMF | 5.10 | 3.628 | -0.1480 |
| xylose/glycine HDMF | 2.60 | 3.628 | +0.1446 |
| ribose/glycine HDMF | 3.60 | 3.628 | +0.0033 |

RMS residual **0.1195 decades = 1.317×**. The pre-registration predicted ≈0.12 decades ≈1.33× and said a materially SMALLER value would mean something had been fitted beyond the one declared parameter.

> The model carries ONE generic aldopentose, so it emits one number for three sugars and the residual IS Blank's 1.96x sugar spread. The pre-registration predicted ~0.12 decades on exactly this ground; a materially smaller RMS would mean something beyond the one declared parameter had been fitted.

### The nine declared-FIT cells the core cannot represent

| table | system | compound | µg/mmol | why not represented |
|---|---|---|---:|---|
| T1 | arabinose/glycine | EHMF | 1.3 | HEMF needs a C2 Strecker donor. The core has no route to it: the sulfur lane (which owns the pentose) carries no alanine, and the acrylamide lane (which owns alanine) carries no pentose. Reported as an unrepresentable FIT row rather than dropped. |
| T1 | xylose/glycine | EHMF | 0.3 | as above. NOTE that it is 0.3 and NOT ZERO -- this single cell is why Amendment 8's 'HEMF requires alanine' truth table had to be demoted to a preference by Amendment 12. |
| T1 | ribose/glycine | EHMF | 0.7 | as above |
| T1 | arabinose/alanine | HDMF | 1.2 | no lane carries pentose + alanine; see the lane-conflict note |
| T1 | xylose/alanine | HDMF | 0.9 | as above |
| T1 | ribose/alanine | HDMF | 1.6 | as above |
| T1 | arabinose/alanine | EHMF | 6.8 | as above |
| T1 | xylose/alanine | EHMF | 7.5 | as above |
| T1 | ribose/alanine | EHMF | 10.0 | as above |

## The derived hexose transfer

`k_odg_af` = **4.028960e-03** 1/min. Rule: k_odg_af [1/min] = k_dpo_af [L/(mmol*min)] x 1000 mmol/L, the pentose edge's pseudo-first-order constant at Blank's 1 M amine loading.

> DECLARED TRANSFER, and the weakest link in the DMHF node. There is NO absolute hexose DMHF yield in any of the five papers of the cluster; only the intact-C6 STRUCTURE is measured, and it is measured twice.

## Ingestion check 1 — Hamzalioglu's HMF + cysteine sink

the DOSSIER'S REFIT prefactor (24115.1 day^-1 pseudo-first-order, /0.020 M) and Ea (29.675 kJ/mol), evaluated against the three second-order constants derived from Hamzalioglu's own Table 1.

| T (°C) | measured (M⁻¹ day⁻¹) | measured (core units) | refit line (core units) | fold |
|---:|---:|---:|---:|---:|
| 5 | 3.95 | 2.7431e-06 | 2.2400e-06 | 1.225× |
| 25 | 5.15 | 3.5764e-06 | 5.2974e-06 | 1.481× |
| 50 | 23.30 | 1.6181e-05 | 1.3374e-05 | 1.210× |

Max fold **1.481×**. NOT 1.00. The three-point refit has R^2 = 0.874 and its slope rests on the single 50 C point at which the authors declare the pseudo-first-order assumption compromised.

**Why the refit prefactor and not the published one:** All six of that paper's activation energies reproduce from its own Table 1 to four decimal places and only TWO of its six pre-exponentials do -- a sign flip on every negative intercept and a SWAP of the Coffee-Cys/Coffee-Lys pair. It is the third audited case in the corpus of a correct Ea bolted to a defective prefactor. HMF-Cys happens to be one of the two that pass (published 23980.59 vs refit 24115.1, 0.56 % apart); the refit is used anyway, because using the published value for the rows that pass and the refit for the rows that fail would be a selection nobody could audit. Amendment 12 makes it binding.

## Ingestion check 2 — Kocadagli's seven glucose constants

each ingested Kocadagli glucose-system constant, re-referenced from the paper's T_b = 180 C to the core's 100 C at import and evaluated back at 180 C here. A pure arithmetic round trip, which is exactly why it is worth doing: the exponent's sign is the one place this ingestion could silently invert.

| parameter | published k_b (min⁻¹×10³) | recovered | rel. error | Ea |
|---|---:|---:|---:|---:|
| `k_glc_tdg` | 4.19 | 4.1900 | 0.00e+00 | 107.2 |
| `k_tdg_ddg` | 30.50 | 30.5000 | 0.00e+00 | 36.9 |
| `k_ddg_hmf` | 119.00 | 119.0000 | 1.19e-16 | 0.0 |
| `k_fru_int` | 330.00 | 330.0000 | 0.00e+00 | 100.4 |
| `k_int_hmf` | 1.84 | 1.8400 | 1.21e-16 | 151.4 |
| `k_fru_odg` | 2.11 | 2.1100 | 0.00e+00 | 99.3 |
| `k_tdg_mgo` | 304.00 | 304.0000 | 1.87e-16 | 84.8 |

Max relative error **1.87e-16**.

## The sourceless constant, swept

``k_af_dmhf`` (acetylformoin -> DMHF) has NO SOURCE OF ANY KIND. It encodes the assumption 'acetylformoin does not accumulate'. The pre-registration requires that sweeping it over three decades move the DMHF prediction by less than 1.2x; if it moves more, the constant is rate-limiting after all and the assumption is false.

| ×`k_af_dmhf` | DMHF (mmol/L) |
|---:|---:|
| 0.1 | 2.767759e-02 |
| 1 | 2.831161e-02 |
| 10 | 2.837498e-02 |
| 100 | 2.838132e-02 |

Span over three decades: **1.0254×** against a pre-registered threshold of 1.2× — **PASS**.

## No fixed branch fraction — demonstrated, not asserted

K5a MUST-NOT #1: the HMF node may not carry a fixed branch fraction. It does not -- the share below is computed from the Fru and 3,4-DG pools at the end of each run, and it MOVES with sugar identity, amine presence and temperature.

| system | HMF (µg/L) | fructose-limb share | 3-DG-limb share |
|---|---:|---:|---:|
| glucose only, 160 C | 2.449e+04 | 0.8441 | 0.1559 |
| fructose only, 160 C | 1.863e+04 | 0.9556 | 0.0444 |
| glucose + glycine, 160 C | 1.196e+04 | 0.9317 | 0.0683 |
| glucose only, 120 C | 2063 | 0.0211 | 0.9789 |

Share range **0.0211 – 0.9556**; it MOVES: **True**.

> Six papers, six matrices, verdicts spanning 'fructose limb dominant' to '3-DG limb dominant by infinity', and EVERY paper that names why its limb wins names a SUPPLY reason -- pool size, a starved 3-DG source, a drained cation pool -- never a terminal rate constant. In four of six published comparisons the terminal constant points the other way.

## Declared gaps carried out of this wave

* NO usable activation energy exists for either HMF-forming edge in any REAL FOOD MATRIX. The only reproducible one (Int -> HMF, 145-152 kJ/mol) is from an amine-free freeze-dried glucose melt over an UNMEASURED intermediate pool, and all three real-matrix systems in the corpus destroy it (wheat flour Ea = -97.6 with R^2 = 0.322; hazelnut non-monotonic; five nuts/seeds with k = 0 in six of ten cells). K5a G1.
* NO HMF-sink constant exists above 50 C that survives audit. Hamzalioglu covers 5-50 C; the hazelnut first-order sink covers 150-160 C but is a lumped decay with no amine dependence. THE 50-150 C WINDOW IS EMPTY, so this module has no effective HMF sink at cooking temperature and must be expected to OVER-PREDICT HMF. K5a G2.
* The fructose-limb intermediate has NEVER been measured in any of the five published networks. Until [Int] or [FFC] is quantified, no absolute rate constant on that limb is recoverable from any of them. K5a G3.
* NO rate constant, NO activation energy, NO time course and NO temperature series exists anywhere in the DMHF cluster. All five papers are single-temperature. K5b sec. 10.
* There is NO absolute HEXOSE DMHF yield in any of the five K5b papers. The intact-C6 STRUCTURE is settled twice over and the MAGNITUDE is measured nowhere; this module transfers the pentose level by a declared rule and flags every hexose answer.
* There is NO measurement of DMHF CONSUMPTION. Edge C ships at exactly zero. Haleva-Toledo et al. 1999 (JAFC 47:4140-4145) is the identified source that would close it, and it would serve the HMF node at the same time.
* NO pH ladder exists anywhere in the HMF cluster. Six distinct pH values appear across the seven papers and NO SINGLE PAPER VARIES pH, so any pH term on the HMF node would be fitted across labs and matrices at once -- which k3 sec. B.2 already forbids at family level. This module therefore carries NO pH term on any HMF edge. K5a G8.
* 3-deoxyglucosone and 3-deoxyGALACTOSONE were not chromatographically resolved in ANY of the five Gokmen multiresponse papers (author-declared once, same method throughout). Every 3-DG number ingested here carries that ambiguity. K5a G6.
* 3,4-dideoxyglucosone is SEMI-QUANTITATED against the 3-DG response factor, so both edges that touch it inherit an unknown multiplicative scale. K5a C22.

