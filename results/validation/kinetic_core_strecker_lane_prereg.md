# Strecker aldehydes lane (plan step C9) -- pre-registration, step 1: inventory and first-arm sketch (DRAFT, 2026-10-31)

*A read-only inventory of what the engine and the corpus hold, written while the sulfur lane's round 3 ran. No code,
no fit. The design (step 4) and its independent review (step 5) come before any number. Target species: 3- and
2-methylbutanal (Leu / Ile), 2-methylpropanal (Val), phenylacetaldehyde (Phe), methional (Met), acetaldehyde (Ala).*

## 1. What the engine has (nothing Strecker-aldehyde is shipped)

- **B18 pyrazine step** (`network.py` PYRAZINE_REACTIONS; `parameters_pyrazine.py` FROZEN_B18): GO + Gly -> AKG and
  MGO + Gly -> AKM, second order, log10 k -6.54 / -7.53 L/(mmol min) at 100 C, Ea 103 / 115, fitted on Zhou 2024's
  PYRAZINE rates under "pyrazine rate = Strecker rate / 2"; pH term from Leahy's pyrazines (+0.58 dex / pH below 7,
  +0.20 above). The trunk carries one amine, `Gly`.
- **B22 methionine chain** (`METHIONINE_SHIPPED = False`): GO / MGO + MET -> methional -> MeSH -> DMDS; failed (identity
  ratio at its x100 ceiling, Pan 2025 still 3.6 - 5.6 dex low, Deng 2022 1.5 - 2.9 dex high). B22b (Amadori ->
  methional) failed its time-course shape. B24 (proline, `r_mgo_pro`) unshipped.
- The sulfur lane's cysteine Strecker (`r_dpo_ddp`) and the MGO + Cys Strecker (rejected, experiment register EX-011).
- Dicarbonyl donors on the trunk: Amadori -> MGO (Martins, 125 kJ/mol), -> 3-DG (97), GO via glucosone (B21 fitted),
  1-DG and diacetyl (Kocadagli 160 - 200 C, extrapolated). **No Strecker step uses 3-DG, glucosone, 1-DG or diacetyl.**

## 2. Data on disk (no keyed fit rows yet)

- **Hofmann 2000** (anchor; aqueous pH 7, 98 C, SIDA): Table 3 fed MGO + Phe time course (7 points to 500 min,
  phenylacetaldehyde and Phe consumed); Table 4 phenylacetaldehyde from MGO / GO / 3-DG / glucosone, argon vs air;
  Table 2 glucose + Phe / Leu / Ala at 30 min; Fig. 3 the dicarbonyls in the glucose + Phe pot (figure only); Figs.
  4 - 5 a pH 3 - 9 series.
- **Hofmann 2000b** (Phe Amadori compound vs glucose + Phe; argon / air / air + Cu; 100 C; SIDA).
- **Weenen 2001** (Phe + 15 carbonyls, pH 3.2, steam-distillation extraction; semi-quantitative endpoints).
- **Chan 1994b** (glucose + Met / Phe / Pro / Leu, pH 6 - 8, 75 - 115 C; averaged Ea only: 3-methylbutanal 80,
  phenylacetaldehyde 90 kJ/mol); **Chan 1994** (methional zero-order rate ladders, digitised, semi-quantitative).
- **Cremer 2000** (glucose + Ala / Val / Ile / Leu on cellulose, aw 0.52, 80 - 110 C; Ea 115 - 124, text only; the
  aldehyde lag tied to the Fru-Leu steady state).
- **Balagiannis 2009** (liver extract, 120 - 140 C, calibrated): glucose -> intermediate pool -> aldehyde constants;
  Ile : Leu yield 1.6 per mole; measured aldehyde sinks (3-methylbutanal 2.7e-3 /min at 130 C, Ea 78 +- 39;
  2-methylbutanal ~0).
- **Huang 2017** (glucose / maltose + Leu / Ile, pH 5.2, 90 - 130 C; constants and Ea, ratios only usable).
- Methional: Pan 2025 (MGO and 3-DG measured in the same pot, figure only), Deng 2022, Meynier 1995 (pH 4.5 - 6.5),
  Zhang 2023b, Yu 1995. Lipid-carbonyl Strecker of Phe: Zamora / Hidalgo 2007 - 2016 (mostly pH 3).
- **The donor measured in the same pot:** only Hofmann 2000 Fig. 3 and Pan 2025 Fig. 5 (both figures) -- digitising
  Hofmann Fig. 3 is the highest-value single task.

## 3. First findings (the inventory's own arithmetic, to be re-derived in the design)

- **B18's constant is not the Strecker step.** Hofmann Table 3's fed MGO + Phe pot makes phenylacetaldehyde at
  ~7e-6 L/(mmol min) (98 C, pH 7); B18 gives ~3e-8 there: ~2.3 dex low -- a pyrazine-yield constant. Phe is
  consumed ~15x faster still: a large non-Strecker amino-acid sink (Strecker products 18 % -> 6 % of the Phe used).
  0.5 M phosphate is a confounder (Rizzi 2004: phosphate catalysis).
- **The missing donors.** At pH 7, MGO : GO : 3-DG : glucosone = 1 : 0.83 : 0.76 : 0.73 in phenylacetaldehyde
  (Hofmann Table 4): the absence of 3-DG and glucosone donors is a likely second cause of B22's sugar-pot shortfall
  (Hofmann's "3-DG era" in the first hour).

## 4. Mechanism signatures the data offer

Donor identity at pH 7 (Hofmann: MGO ~ GO ~ 3-DG ~ glucosone) vs pH 3.2 (Weenen: MGO : GO : 3-DG = 35 : 15 : 1);
pH direction -- fed pots peak at pH 5 (pH 5 / 9 = 2) while sugar pots give pH 7 / 3.2 = 25 (the step and the supply
have opposite signs, and both oppose B18's pyrazine term); the Strecker share of Phe consumed; oxygen independence of
the aldehyde from dicarbonyls (0.85 - 0.93; the acid branch parallel, not an oxidation); amino-acid identity (Leu /
Phe 2.1, Hofmann; Ile / Leu 1.6, Balagiannis, Huang); Chan's Ea order Pro < Leu < Phe; aldehyde losses (plateau and a
10 % fall by 500 min; 3-methylbutanal sinks, 2-methylbutanal does not -- aldol); Cremer's lag following Fru-Leu.

## 5. Gaps for a targeted search (DOIs for the owner, to be verified before any download)

G1 the Strecker step's own barrier (a fed-dicarbonyl + Leu / Phe aldehyde ladder at >= 3 temperatures in water,
dicarbonyl followed) -- nothing on disk; G2 diacetyl as a donor (Pripis-Nicolau 2000, JAFC 48:3761,
10.1021/jf991024w, unverified at doi.org); G3 3-DG / 1-DG donors at pH 5 - 7 (Ghiron 1988; Fujimaki 1968 -- DOIs not
found); G4 Strecker-aldehyde loss kinetics vs pH; G5 Val -> 2-methylpropanal in water (Kocadagli 2021's supplement).

## 6. First-arm sketch (C9-A, trunk, switch off) -- for the design step, not yet a design

Species Leu / Ile / Val / Phe / Ala and their aldehydes (Met / methional reused). Steps: donor + amino acid ->
aldehyde + aminoketone + CO2 for donors MGO, GO, 3-DG, glucosone (rate k_MGO,Phe x donor ratio x identity ratio x
f(pH)); amino acid + donor -> melanoidin (the non-Strecker sink, Hofmann Table 3); first-order aldehyde loss.
Coordinates: k_MGO,Phe at 100 C, three donor ratios, a two-slope pH term (knot near pH 5), Leu identity; Ile = Leu x
1.6 (+-0.15 dex, Balagiannis), Val / Ala declared (+-0.5); losses at Balagiannis's; Ea declared 80 - 115 (G1 open).
Fit rows: Hofmann 2000 Tables 2 - 4 (argon). Hold-outs: Weenen 2001, Hofmann 2000b, Balagiannis 2009, Huang 2017
ratios, Chan 1994b, Cremer 2000's lag, Deng 2022. Expected: >= 7 of the 11 signatures (P 0.55); Weenen's pH x donor
contrast (P 0.25); the step / supply pH sign flip (P 0.4); sugar-pot absolute levels within x3 across laboratories
(P 0.2: the dicarbonyl supply limits it).

## 7. Step 2 done (2026-10-31): Hofmann 2000 keyed by two readers and adjudicated

`data/lit/fit_rows/strecker_hofmann2000.yml`, status ready, 111 rows (Tables 2 - 4, Fig. 3's donor time courses, Figs.
4 - 5's pH series); reconciliation `results/validation/source_reading_reconciliation_strecker_hofmann2000.md`;
experiment register EX-026. It loads as literature only: the engine has no species for Phe, phenylacetaldehyde /
phenylacetic acid, Ala, Leu, 3-methylbutanal / 3-methylbutanoic acid (provisional key MEBUAC3), which arrive with the
lane on the composable network (STATUS_AND_PLAN A2, then C9). For the design: Table 3's phenylalanine rows are
semi-quantitative and probably duplicate hofmann2000b Table 3 (CF-070) -- at most one enters a fit; Table 4's buffer
is 0.1 mol/L as printed (alternative 0.5); the Figs. 4 - 5 pH points carry declared conditions (60 min, air), not
printed ones, so they are signatures (pH shape), not level rows.
