# Heng 2005 (PhD thesis) — EXTRACTION (pea vicilin and legumin binding of aldehydes and ketones, heat and pH; volatiles carried by pea flour, isolate, legumin and vicilin)

**Source on disk:** `data/articles/Heng2005.pdf` (the Wageningen thesis, 176 pp.; PDF page = printed page + 10).
Table of contents read first (`pdftotext -l 12`); then **Chapter 5** (pp. 95-108), **Chapter 6** (pp. 109-130)
and the volatile half of **Chapter 7, General discussion** (pp. 131-136) read by eye on 2026-10-09. Every
number below was read on the page image (Read tool) and cross-checked against `pdftotext -layout`; the two
agree on every value quoted. Chapters 1-4 (introduction; saponins) were not read. Written for the
plant-matrix binding gap and the carried-volatile question for pea ingredients, not for the core fit.

Already on file before this read: `data/lit/binding_constants.yml` source `heng_2005_wur_thesis` with four
`percent_bound_at_conditions` records from Ch. 6 Table 1 (2-octanone, 2-heptanone, 2-hexanone at 10 g/L;
octanal marked `usable_for_model: false` for an ambiguous loading) and one denaturation record
(`heng2005_denaturation_vicilin`); `data/keys/papers.yml` `doi_10_18174_121674` with `dossier: null`. No
published paper from this thesis is in the dossier folder. Three dossiers cite "Heng 2004" second-hand
(`guo2019`, `fischer2021`, `wang2014`); Ch. 6 says it "is based on" that 2004 review (Trends Food Sci.
Technol. 15:217-224), which is unread here. The second-hand figures ("75-88 % of aldehydes", "pentanal 75 %")
do not match the thesis's own Table 1 (pentanal 52 %, octanal 96 %); the thesis is the primary.

| field | value |
|---|---|
| Title | "Flavour aspects of pea and its protein preparations in relation to novel protein foods" |
| Author | Lynn Heng |
| Degree | PhD thesis, Wageningen University, defended 2 June 2005; ISBN 90-8504-198-8 |
| DOI | 10.18174/121674 (as in `papers.yml`, crossref-verified there 2026-08-27) |
| Promotors | A.G.J. Voragen, M.A.J.S. van Boekel; co-promotor J.-P. Vincken |
| Chapter 6 basis | Heng, van Koningsveld, Gruppen, van Boekel, Vincken, Roozen, Voragen (2004) TFST 15:217-224 (p. 109) |
| Chapter 5 | no publication note |

## 1. Methods

**Proteins (Ch. 5, pp. 98-99; Ch. 6, p. 113).** Dried de-hulled green split peas (*Pisum sativum*, cv.
Solara), milled with dry ice 1:1 "to avoid possible heat denaturation". Ch. 5: extraction pH 8 (100 mM
Tris-HCl), acid precipitation at pH 4.8 ("protein isolate"), re-solubilisation, re-precipitation, DEAE anion
exchange to legumin and vicilin; freeze-dried. So the "isolate" is a laboratory, never-heated,
freeze-dried acid-precipitate, not a commercial spray-dried isolate. Ch. 6: vicilin by batch DEAE then
DEAE column; protein by Dumas, N × 5.7.

**Volatiles carried by the preparations (Ch. 5, pp. 99-100).** Dynamic headspace: 100 mg dry sample in
10 mL McIlvaine buffer, **pH 4 (pea flour only) or pH 8**, hydrated 12 h at 4 °C, purged with N2 at
45 mL/min for **1 h at 37 °C** onto Tenax TA; thermal desorption, GC-FID (Supelcowax 10), identity by GC-MS
library and retention time. **Quantities are arbitrary units (AU) per g dry matter**, printed only as
+ to +++++ bands; no internal standard and no response factors ("the response factors of all aldehydes are
similar", p. 101). No lipid content of any preparation is measured in Ch. 5.

**Binding (Ch. 6, pp. 116-117).** Static headspace GC-FID. 0.5 mL ligand solution + 0.5 mL protein solution
in an **11.5 mL** vial, so 1 mL liquid and a gas-to-liquid volume ratio of **10.5** (derived here:
(11.5 − 1)/1). 10 mM potassium phosphate **pH 7.6, I = 0.024 M**; protein solutions filtered 0.2 µm.
Incubation **37 °C, 10 min, 500 rpm**; 1 mL headspace injected. Final concentrations: aldehydes
0.006-0.04 mM (0.8-3.3 ppm), ketones 0.03-1.2 mM (4-103 ppm), protein 0.1-1 % w/v; triplicates; blank =
ligand in buffer. "% bound" is "calculated as a proportion of the total volatiles added to the system"
(Table 1 footnote, p. 119); the headspace-to-bound arithmetic is not printed. **Heat pre-treatment: capped
vials 90 °C for 30 min**, then read at 37 °C (p. 116). Ligands: pentanal, octanal, 2-pentanone,
2-hexanone, 2-heptanone, 2-octanone (vicilin); octanal, heptanal, hexanal, pentanal (legumin, Ch. 7).
**No binding constant (K, n, Klotz) is fitted anywhere**; only bound amounts (nmol/mg) and % bound.

**pH and non-protein components (Ch. 6, pp. 114-116).** 0.2 % w/v vicilin taken to pH 4.5, 5.5, 6.5,
centrifuged 1500 × g; pellet and supernatant re-suspended at pH 7.6 and tested with octanal 0.026 mM.
Defatting of freeze-dried vicilin by 6 h Soxhlet with hexane or chloroform:methanol 1:1. Lipid = sum of
fatty acids by FAME GC-FID (margaric methyl ester internal standard); carbohydrate by methanolysis/TFA,
HPAEC-PAD.

## 2. Findings that matter

### 2a. Aldehyde and ketone binding to vicilin, pH 7.6, 37 °C (Ch. 6 Table 1, p. 119; Figs 1-2, pp. 118, 120)

| ligand at 0.025 mM | vicilin loading | "octanal affinity" (nmol/mg) | % bound, native | % bound, 90 °C/30 min | K_eff native (L/g, derived here) | K_eff heated (L/g, derived here) |
|---|---|---|---|---|---|---|
| octanal | 0.1 % w/v (Fig. 1) | ~23 | 96 | 32 | 24 | 0.47 |
| pentanal | 0.1 % w/v (Fig. 1) | ~15 | 52 | — | 1.08 | — |
| 2-octanone | 1 % w/v (Fig. 2) | <1.5 | 33 | 16 | 0.049 | 0.019 |
| 2-heptanone | 1 % w/v | <1 | 19 | — | 0.023 | — |
| 2-hexanone | 1 % w/v | <1 | 13 | — | 0.015 | — |
| 2-pentanone | 1 % w/v | <0.5 | 10 | — | 0.011 | — |

The column header "Octanal affinity" is printed for every row. K_eff = [b/(1 − b)] / c_protein, with b the
printed % bound and c_protein in g/L (derived here; it ignores the headspace share, which the thesis's own
"% of total added" also does not separate).

**The loading ambiguity recorded in `binding_constants.yml` is resolved by the thesis's own arithmetic**
(derived here). At 0.1 % w/v (1 mg protein in the 1 mL vial) and 0.025 mM (25 nmol added), 23 nmol/mg × 1 mg
= 23 nmol = 92 % bound, against the printed 96 %; at 1 % w/v it would be 230 nmol bound out of 25 nmol added,
which is impossible. For the ketones the reverse holds: 33 % of 25 nmol over 10 mg is 0.83 nmol/mg, inside
the printed "<1.5", whereas 0.1 % w/v would need 8.3 nmol/mg. So the aldehyde rows are at **1 g/L** and the
ketone rows at **10 g/L**, matching the figure captions. Pentanal is less clean: 15/25 = 60 % against the
printed 52 % (the "~15" is approximate).

**Heat (pp. 120-121).** The lower-chain aldehydes and ketones "showed no binding to vicilin after heating
(results not shown)". Gel filtration puts aggregates at ≥ 4200 kDa (≈ 83 monomers of 50 kDa) after
90 °C/30 min, and a sphere model gives a 77 % smaller surface. The thesis says binding "decreased by 17 % and
64 %" for 2-octanone and octanal; those are **percentage-point** drops (33 − 16, 96 − 32, derived here). The
relative drops are 52 % and 67 %, and in K_eff the drops are **2.6×** (2-octanone) and **51×** (octanal),
derived here. So the thesis's comparison with the 77 % surface loss, and its conclusion that affinity is
unchanged, rest on a points-versus-percent mix-up. Read from Fig. 3 (approx.): heated octanal binding is
~2.7, ~8.1 and ~8.5 nmol/mg at ~0.006, ~0.019 and ~0.026 mM, flattening, while native binding keeps rising to
~23 nmol/mg.

### 2b. pH and co-purified lipid (Ch. 6 Tables 2-3, pp. 123, 125; Fig. 5a, p. 123)

| pH | protein in pellet (% of total, text p. 122) | octanal bound, supernatant (nmol/mg, ± SD) | pellet | supernatant after CHCl3:MeOH | pellet after CHCl3:MeOH | pellet lipid, µg (± SD) | pellet carbohydrate, µg |
|---|---|---|---|---|---|---|---|
| 4.5 | 17 | 5.8 (0.5) | 23.5 (2.3) | 1.4 (0.2) | 18.6 (1.1) | 151 (10) | 334 (15) |
| 5.5 | 61 | 4.8 (0.5) | 6 (0.9) | 5.8 (0.5) | 8.7 (0.7) | 248 (31) | 592 (15) |
| 6.5 | 6 | 4.4 (0.4) | 41.4 (1.2) | 4.9 (0.2) | 19.3 (1) | 96 (6) | 126 (12) |
| 7.6 | — (all soluble) | 22.1 (0.9) | — | 9.3 (0.9) | — | — | — |

Octanal 25 µM; affinity per mg of protein. Hexane-extracted values (Table 2) are 2.3 / 6.5 / 5 / 19.7
(supernatants) and 27.4 / 8.3 / 40.1 (pellets): hexane removed ~25 % of pellet lipid and most carbohydrate
without lowering binding; chloroform:methanol removed ~75 % of the lipid and lowered pellet binding at pH 4.5
and 6.5 (pp. 125). The thesis attributes this to polar lipids, probably DGDG-type galactolipids, from a
removed lipid-to-carbohydrate weight ratio of about 2 (p. 126); that assignment is an inference, not an
identification.

Three derived readings (derived here). (i) Soluble vicilin binds octanal **3.8-5.0× less** at pH 4.5-6.5
than at pH 7.6 (22.1 / 5.8 to 22.1 / 4.4). (ii) Even the fully soluble pH 7.6 fraction loses **58 %** of its
octanal binding after chloroform:methanol ((22.1 − 9.3)/22.1), so more than half of the "protein" binding of
a chromatographically purified globulin sits on co-purified polar lipid or on what the solvent did to the
protein; the thesis cannot separate the two. (iii) The pellets alone carry 96-248 µg fatty-acid lipid from
20 mg vicilin (10 mL × 0.2 %), i.e. **≥ 0.48-1.24 %** of the vicilin mass (supernatant lipid was not
measured, so this is a floor).

### 2c. Legumin, aldehydes including hexanal (Ch. 7 Fig. 1, p. 133; text p. 132)

0.1 % w/v legumin, pH 7.6. Printed: "At 0.13 mM of octanal, about 70 nmoles ... and at 0.18 mM of pentanal,
about 50 nmoles ... constituting ~51% and ~28%". Hexanal and heptanal are figure-only:

| ligand | added (mM), read from graph, approx. | bound (nmol/mg), read from graph, approx. | fraction bound (derived here) | K_eff, L/g (derived here) |
|---|---|---|---|---|
| hexanal | 0.040 / 0.081 / 0.121 / 0.162 | 13.5 / 22 / 35.5 / 61 | 0.34 / 0.27 / 0.29 / 0.38 | 0.51 / 0.37 / 0.42 / 0.60 |
| heptanal | 0.036 / 0.072 / 0.108 / 0.144 | 13 / 21.5 / 56.5 / 63.5 | 0.36 / 0.30 / 0.52 / 0.44 | 0.57 / 0.43 / 1.10 / 0.79 |
| octanal | 0.032 / 0.064 / 0.096 / 0.128 | 16 / 27 / 41.5 / 69.5 | 0.50 / 0.42 / 0.43 / 0.54 | 1.00 / 0.73 / 0.76 / 1.19 |
| pentanal | 0.047 / 0.094 / 0.141 / 0.188 | 21 / 29.5 / 43 / 51.5 | 0.45 / 0.31 / 0.30 / 0.27 | 0.81 / 0.46 / 0.44 / 0.38 |

Fraction bound = nmol/mg × 1 mg ÷ (mM × 1000 nmol), K_eff as in 2a. **Pea legumin binds hexanal at
K_eff ≈ 0.4-0.6 L/g** (read from graph, approx.; derived). The text's "~51 %" for octanal checks at 70/130 =
54 %, and "~28 %" for pentanal at 50/180 = 28 % (derived here). Legumin "was found to have affinity for
aldehydes only" (no ketone binding). The text adds "At 0.03 mM octanal, vicilin retains ~3500 nmoles/mg,
whereas legumin retains ~16 nmoles/mg" (p. 132); **3500 is a printed error**: 0.03 mM in 1 mL is 30 nmol in
total over 1 mg protein, so the ceiling is 30 nmol/mg, and Ch. 6 gives ~23 at 0.025 mM. Per mg, vicilin binds
octanal ~1.5× legumin at ~0.03 mM (23 against ~16), not 200×.

### 2d. Volatiles carried by pea preparations (Ch. 5 Tables 1-2, pp. 102, 104; text pp. 101, 103, 105-106)

| AU per g dry matter | pea flour pH 4 | pea flour pH 8 | lab isolate pH 8 | legumin pH 8 | vicilin pH 8 |
|---|---|---|---|---|---|
| hexanal (band) | +++++ | +++++ | +++ | + | + |
| total aldehydes | 199110 | 355510 | 18142 | 26669 | 11926 |
| total ketones | 25900 | 19690 | 1598 | 22279 | 26413 |
| total alcohols | 28200 | 98690 | 3855 | 256421 | 140052 |
| total volatiles | 255800 | 477390 | 25987 | 327224 | 202634 |
| number of volatiles | 18 | 21 | 10 | 23 | 18 |

Bands: + 500-5000, ++ 5000-10000, +++ 10000-50000, ++++ 50000-100000, +++++ 100000-500000 AU. Hexanal is
">35 %" of total release from flour at both pHs and ">34 %" from the isolate; it fell ">90 %" from flour to
isolate and further to legumin and vicilin. Pentanal, 2-hexenal, 2-heptenal and 2-octenal were released from
flour only at pH 8; decanal and undecanal only at pH 4. Flour at pH 8 released 1.87× the total volatiles and
1.79× the aldehydes of pH 4 (derived here). High-log P compounds (2-nonenal, undecanal, 2-decanone, decane,
2-ethyl-1-hexanol, which rose > 30× and > 15× into legumin and vicilin) were enriched by purification; the
thesis explains release versus retention by log P (pp. 105-106). 2-Pentylfuran appears only in legumin (+).
The 3-alkyl-2-methoxypyrazines reported for green peas were **not found** in these dried peas (p. 133).

### 2e. Sensory and thresholds

None for any volatile in a pea matrix. The only sensory data in the thesis are for saponin **bitterness**
(Ch. 3, not read; summarised p. 135-136): DDMP saponin perceived below 2 mg/L, saponin B at ~8 mg/L in
water; Ch. 7 then estimates 4.4 g saponin per kg of a 55 %-protein meat analogue, its own back-of-envelope.
The odour-threshold content is limited to citing Buttery 1973 (oil lowers aldehyde thresholds) on p. 135.

## 3. What it means for the model and for formulation

**Binding layer.** The live pea constant is `kg_hexanal_pea` = **2.537e-1 L/g**
(`src/kinetic_core/parameters_matrix.py`, anchor "bi2022_extraction.md sec. 3.8 / Table S3", pea isolate
10 g/L, 0.01 M potassium phosphate pH 7.6, 37 °C, FIT). Heng's legumin hexanal at the same buffer, pH and
temperature gives **K_eff ≈ 0.37-0.60 L/g** (2c, graph reading, derived), **1.5-2.4×** the live value, a
second laboratory, a purified globulin, 1 g/L instead of 10 g/L and a different data reduction. That is
agreement inside the factor-of-2.5 method spread the registry already flags for Bi 2022
(`within_paper_method_spread_x: 2.5`). It **tests** `kg_hexanal_pea` and does not move it; it cannot pin it
because it is figure-only, has no replicate errors on most points, and its "% bound" ignores the 10.5 phase
ratio.

The vicilin aldehyde rows are a different story: octanal at K_eff 24 L/g native and 1.08 L/g for pentanal
(derived) are 4-95× the pea hexanal constant, and 2c shows legumin at ~1 L/g for octanal. Two things in the
thesis say this is not protein alone: 58 % of the soluble fraction's octanal binding goes with
chloroform:methanol (2b), and the insoluble, lipid-rich pellets bind up to 41.4 nmol/mg. For a formulator
the useful number is the **ratio**, not the absolute: an aldehyde-binding constant measured on a purified pea
globulin can be dominated by co-purified polar lipid, and an isolate's binding will move with its residual
polar-lipid content.

**Heat (the gap `BINDING_AT_PROCESS_TEMPERATURE` names).** The registry says "No paper in the corpus
measures an AQUEOUS binding constant at a temperature above 60 C" and lists cooked-then-read-cold
experiments that split by ligand. Heng is another cooked-then-read-cold one, purified vicilin, 90 °C for
30 min, read at 37 °C: **both carbonyl classes fall**, octanal by 51× and 2-octanone by 2.6× in K_eff
(derived), and C5-C7 ligands lose all measurable binding. Unlike Xu 2022 (2-methylpyrazine rises on pea
heating) there is no ligand here that gains. It does not license a binding constant at 90 °C either.

**pH.** Every pea constant in the registry is at pH 7.6. Soluble vicilin binds octanal 3.8-5.0× less at
pH 4.5-6.5 (2b). A meat analogue at pH 5.5-6.5 would then sit well below the pH 7.6 value for the soluble
protein, but the insoluble fraction at those pHs (lipid-rich) binds more, so the sign of the net change for a
whole isolate is not fixed by this thesis.

**Lipid lane and carried volatiles.** The live carrier `lipid.pea_protein_isolate.lipid_mass_fraction` is a
`declared_band` [0.01, 0.06] kg/kg, centre 0.025 (`results/validation/core_prediction_uncertainty.json`,
source `parameters_lipid.LIPID_CARRIERS['pea_protein_isolate'] (declared_assumption)`). Heng's ≥ 0.48-1.24 %
fatty-acid lipid in the insoluble part of a **purified** vicilin (2b) is a floor that sits at or below the
band's lower corner; it is consistent with the band and does not narrow it (wrong material, FA-sum basis,
supernatant not measured). Ch. 5's hexanal cannot fill `conditions.carried_volatiles`, which the engine asks
for when it refuses an uncooked pea hexanal row: the levels are AU bands, not concentrations. What Ch. 5 does
give the formulation side is a ranking: acid precipitation strips >90 % of flour hexanal and all the
short-chain alkenals, while purification **concentrates** hydrophobic ketones, long aldehydes and
2-ethyl-1-hexanol into legumin and vicilin. For hexanal and the alkenals the more purified fraction is the
cleaner base for added meat aroma; for those hydrophobic off-notes it is the dirtier one.

## What it does not give

- No Klotz, Scatchard or partition constant; no binding isotherm fitted; no binding at any temperature
  above 37 °C. Heat is a 90 °C/30 min pre-treatment, read cold.
- No hexanal binding to vicilin and no ketone binding to legumin (legumin binds "aldehydes only"); no
  commercial isolate in any binding experiment; no sulfur compound, pyrazine or furan ligand.
- No absolute volatile concentration in any pea preparation (AU bands only), no peroxide value, no lipid
  content of flour, isolate, legumin or vicilin in Ch. 5, no heat-treated preparation in Ch. 5.
- No odour threshold, GC-O or sensory test of any volatile in a pea matrix.
- Ch. 6's affinity drop with heat is mis-stated in percentage points (2a), and Ch. 7's "~3500 nmoles/mg" is
  impossible (2c); both are flagged, not corrected in the source.
