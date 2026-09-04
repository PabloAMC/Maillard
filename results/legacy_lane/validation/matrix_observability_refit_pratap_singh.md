# Matrix observability refit against the content-corrected Pratap-Singh anchors (Wave O)

Generator: `scripts/generators/refit_matrix_observability_pratap_singh.py` — 2026-08-27, owner-approved.

Fit targets (both already `fitted_to_benchmark` in the registry, both excluded from
the honest literature-coverage count before AND after this refit):

- `pea_isolate_40C_PratapSingh2021` — hexanal **1138.00 ± 297.30 ppb**
- `soy_isolate_40C_PratapSingh2021` — hexanal **1621.71 ± 159.69 ppb**

Molecules 2021, 26, 4104, Table 1 (Europe PMC PMC8271896). The repo's previous
260 / 380 ppb appear nowhere in the paper; the shipped constants reproduced them to
four significant figures, i.e. the model reproduced the transcription error.

## What was fitted

ONE parameter: a common multiplicative scale on the ambient-slurry **hexanal**
observability factors of the pea and soy lanes. No barriers, no projection constants,
no marker yields, no 2-pentylfuran factor (those anchors are verified verbatim).

- objective at incumbent: **0.807024 dex²**
- objective at the fit:  **0.000048 dex²**
- adopted shared scale:  **4.31725x**
- search range 0.001 – 30 (4.5 decades); optimum sits 3.64 dex above the floor and 0.84 dex below the ceiling — **no bound was hit, nothing saturated**

## Mutual consistency of the two corrected anchors

| lane | required scale |
| --- | --- |
| `pea_iso/ambient_slurry/hexanal` | 4.36606x |
| `soy_iso/ambient_slurry/hexanal` | 4.26899x |

A single shared scale of 4.31725x satisfies both to within **1.0113x**. The correction is therefore
almost purely an absolute-scale error on the ambient hexanal lane; the pea-vs-soy
release ratio the registry encodes survives it. That residual is the informative
part of this fit, and it only exists because the second degree of freedom was
declined — see `alternative_per_lane_fit` in the JSON.

## Adopted constants

| lane | role | previous | refit | move |
| --- | --- | --- | --- | --- |
| `pea_iso/ambient_slurry/hexanal` | fitted | 1 | 4.31725 | 4.3172x |
| `soy_iso/ambient_slurry/hexanal` | fitted | 2.20976 | 9.54007 | 4.3172x |
| `soy_iso/heated_matrix/hexanal` | propagated | 0.649668 | 2.80478 | 4.3172x |

## Rows

| benchmark | compound | measured ppb | predicted (incumbent) | predicted (fit) | fold before | fold after | in objective |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `pea_isolate_40C_PratapSingh2021` | hexanal | 1138 | 260.647 | 1125.28 | 4.3661x | 1.0113x | yes |
| `pea_isolate_40C_PratapSingh2021` | 2-pentylfuran | 638 | 638.267 | 638.267 | 1.0004x | 1.0004x | no |
| `pea_isolate_uht_140C_Trikusuma2019` | hexanal | 782 | 782.006 | 782.006 | 1.0000x | 1.0000x | no |
| `pea_isolate_uht_140C_Trikusuma2019` | 2-pentylfuran | 163 | 163 | 163 | 1.0000x | 1.0000x | no |
| `pea_isolate_uht_140C_Trikusuma2019` | nonanal | 24 | 24 | 24 | 1.0000x | 1.0000x | no |
| `soy_isolate_40C_PratapSingh2021` | hexanal | 1621.71 | 379.882 | 1640.04 | 4.2690x | 1.0113x | yes |
| `soy_isolate_40C_PratapSingh2021` | 2-pentylfuran | 2492 | 2492.29 | 2492.29 | 1.0001x | 1.0001x | no |

## Reported, not fitted

- `pea_iso/ambient_slurry/2-pentylfuran` — implied scale 0.999581x. 638 ppb verified VERBATIM against Molecules 26:4104 Table 1 (Wave K). The shipped factor already recovers it; moving it would be fitting rounding.
- `soy_iso/ambient_slurry/2-pentylfuran` — implied scale 0.999882x. 2492 ppb verified VERBATIM against the same table. Same reading.
- `pea_iso/heated_matrix/hexanal` — implied scale 0.999993x. Trikusuma 2019 was NOT content-corrected — it is also the last unverified pillar of the matrix lane (Wave K: full text not retrieved). Refitting it here would move a constant against an anchor nobody has read.
- `pea_iso/heated_matrix/2-pentylfuran` — implied scale 0.999997x. Same as the heated-pea hexanal entry.
- `pea_iso/heated_matrix/nonanal` — implied scale 1x. Same as the heated-pea hexanal entry.

## Unanchored after the content correction

- `pea_iso/ambient_slurry/1-hexanol` — The paper reports n.d. for pea 1-hexanol and its text states pea proteins 'contained no alcohol compounds'. The 80 ppb this factor was built from is not a measurement. There is nothing to fit against; the factor is left at 1.0.
- `soy_iso/ambient_slurry/1-hexanol` — Soy 1-hexanol is n.d. as well (the paper's soy total alcohols are 40 +/- 9 ppb, 1-octen-3-ol only). The shipped factor 0.143/0.063 is a RATIO OF TWO FABRICATED NUMBERS (120 and 80 ppb), neither of which appears in the paper. It is left untouched because there is no anchor to fit it to -- but it is not a calibration, and it is the factor behind the hold-out's 1117x 1-hexanol miss on Li 2026. Retiring it is a science decision outside this refit's approved scope.
- `pea_iso/ambient_slurry/nonanal` — Nonanal is not reported for pea in this benchmark. Factor stays 1.0.
- `soy_iso/ambient_slurry/nonanal` — Nonanal is not reported for soy either; 0.160/0.150 is lane-internal. Untouched.

## Hold-out

This script never reads the hold-out. Scored AFTER the constants landed, by the normal generator, the hold-out median |log10| error moves from 1.1849 dex (15.31x) to a worse value, because the pea ambient lane contains two mutually contradictory external measurements (Bi 2020: 1260 ppb; Liu 2023 band midpoint: 51.96 ppb, a 24x spread at nominally the same conditions) and the erroneous 260 ppb constant happened to sit almost exactly at their geometric mean (sqrt(1260 x 51.96) = 255.9). Moving to the verified anchor moves the lane onto the Bi side of that disagreement. See tasks/audit_remediation.md, Wave O, for the 8-point table.

## Shipped-vs-refit check

tolerance 0.0005 dex — **all shipped constants agree with this refit**

- `pea_iso/ambient_slurry/hexanal`: shipped 4.31725, refit 4.31725, delta 9.27e-08 dex — ok
- `soy_iso/ambient_slurry/hexanal`: shipped 9.54007, refit 9.54007, delta 1.15e-07 dex — ok
- `soy_iso/heated_matrix/hexanal`: shipped 2.80478, refit 2.80478, delta 1.15e-07 dex — ok
