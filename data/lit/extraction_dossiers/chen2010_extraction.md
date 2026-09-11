# Chen, Wei, Zhang & Ojokoh 2010 — EXTRACTION (soy protein isolate extruded across a wide moisture range: system parameters including mean residence time)

**Source on disk:** `data/articles/chen2010.pdf` (0.52 MB; downloaded 2026-09-11 at this
repository's request). Read 2026-09-11 via `pdftotext -layout`. Wave B37.

| field | value |
|---|---|
| Title | "System parameters and product properties response of soybean protein extruded at wide moisture range" |
| Venue | Journal of Food Engineering 96 (2010) 208–213 |
| DOI | 10.1016/j.jfoodeng.2009.07.014 — read from the printed footer |
| Group | Chinese Academy of Agricultural Sciences (Wei Yimin) |
| System | **soy protein isolate alone** (Yuwang Group Ltd.), twin-screw extruder |
| Design (Table 1) | 3 × 5 factorial: **moisture 28, 36, 44, 52, 60 %** × **cooking temperature at the middle zone 140, 150, 160 °C**; screw speed fixed at **160 rpm**, feed rate fixed at **20 g/min** |
| What is measured | in-line viscosity at the die, **mean residence time**, specific mechanical energy (SME), and product texture |

## 1. Why this paper was fetched, and the limitation that caps its use

It was fetched as the second anchor for the extrusion residence time that
`acrylamide_spi_extrusion_130C_ACSRef3` assumes (see `yu2012_extraction.md`), and it is the closer
match on composition: **pure soy protein isolate**, where Yu et al. run 20 % SPI in corn flour, and
it reaches the model's 28–30 % moisture and 140–160 °C.

**But the mean residence time is FIGURE ONLY.** It appears as Fig. 3b — "In-line viscosity (a), mean
residence time (b) and SME (c) versus moisture content and cooking temperature with screw speed
160 rpm and feed speed 20 g/min" — and **no residence-time value is printed anywhere in the text or
in either table**. Table 1 is the factorial design; Table 2 is a correlation matrix between system
parameters and product properties. This dossier therefore transcribes **no residence time from this
paper**, and the repository takes none.

## 2. Table 2 — correlation analysis between system parameters and product properties (verbatim)

| property | in-line viscosity at die | mean residence time | SME |
|---|---:|---:|---:|
| Tensile strength | 0.86** | 0.82** | 0.96** |
| Hardness | 0.72** | 0.76** | 0.76** |
| Chewiness | 0.70** | 0.75** | 0.74** |
| Degree of texturization | 0.36 | 0.45 | 0.49 |

`**` significant at p < 0.01, as printed.

## 3. What the repository takes

A **direction only**: verbatim, *"increasing moisture content could result in accelerating the flow
speed of extrudate coming out from extruder"*, confirmed by the effect of moisture on mean residence
time — so residence time **falls** as moisture rises, over 28–60 % in pure SPI at 160 rpm. That
agrees with Yu et al.'s measured table and with Ma et al. 2024's own reasoning about moisture, feed
rate and residence time. It supports no number.

## 4. Verdict

A supporting directional source, not a quantitative one, because its residence times were never
printed. Recorded so that a future reader does not fetch it twice expecting a table.
