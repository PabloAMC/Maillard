# The papers the model rests on, and the overlap with an expert reading list

*Written 2026-10-27 after an outside expert shared a reading list of 24 papers on volatile
formation from amino acids, sugars, thiamine and glutathione. Refreshed 2026-11-24: counts recomputed on main,
every paper on the list now read, and a section on the papers behind the thiol forecasts. Counts are the
repository's own: rate constants taken from a paper, fit rows a calibration was tuned on, benchmark pots the model
is scored against, and directional claims it is scored on (`build_thiol_sink_figures.paper_usage_counts`; every
paper and what is taken from it: [SOURCES.md](SOURCES.md)).*

## 1. The fifteen papers the model rests on most

Across the whole model (sugar breakdown, sulfur, acrylamide and lipid chemistry). Counts recomputed 2026-11-24.

| # | paper | records | what the model takes | on the expert list |
| ---: | --- | ---: | --- | :---: |
| 1 | Hofmann & Schieberle 1998 (cysteine + sugars, isotope dilution) | 77 | 41 fit rows, 23 benchmark pots, 13 claims | |
| 2 | Zhou et al. 2023 (Amadori compound + cysteine, pH drift) | 16 | 10 fit rows, 4 benchmark pots, 2 claims | |
| 3 | Kang et al. 2026 (cysteine–xylose, TTCA) | 14 | fit rows | |
| 4 | Whitfield & Mottram 1999 (norfuraneol + cysteine / H2S) | 14 | fit rows | |
| 5 | **Cerny 2015** (sulfur chemistry in thermal aroma generation, book chapter) | 11 | directional claims | **yes** |
| 6 | Martins & van Boekel 2003 (Amadori degradation kinetics) | 10 | rate constants | |
| 7 | Martins & van Boekel 2005 (glucose–glycine kinetics) | 10 | rate constants | |
| 8 | Kocadağlı & Gökmen 2016 (sugar degradation) | 9 | rate constants | |
| 9 | Trikusuma et al. 2018 (pea protein, lipid volatiles) | 9 | claims | |
| 10 | Mottram & Nobrega 2002 (buffered vs unbuffered cysteine–ribose) | 8 | claims | |
| 11 | De Vleeschouwer et al. 2009 (acrylamide kinetics) | 8 | constants, fit rows | |
| 12 | Yiltirak et al. 2026 (ribose–cysteine, temperature ladder) | 8 | hold-out pots, fit rows, claims | |
| 13 | Ma et al. 2024 | 7 | claims | |
| 14 | Zhang et al. 2024 | 6 | fit rows, claims | |
| 15 | Ni et al. 2021 | 6 | claims | |

Further down the same ranking, also on the expert list: **Cerny 2007 (LWT, origin of carbons)** — a fit row (the
54:46 thiamine:xylose carbon split of MFT) and two hold-outs since wave B2; **Cerny & Briffod 2007 (JAFC,
pH)** — a claim and fit-system conditions; **Meynier & Mottram 1995** — the conditions of the two pH claims
PH-01/02.

## 2. The papers behind the thiol forecasts

The ranking above counts records across the whole model. For the meaty thiols (2-methyl-3-furanthiol,
2-furfurylthiol and their disulfides) from cysteine, thiamine, glutathione and a pentose, these are the papers the
sulfur lane's constants come from. Several are stored as module constants in
`src/kinetic_core/parameters_sulfur.py`, which the counter in section 1 does not see, so this list is curated
from that file.

| paper | what the sulfur lane takes | DOI |
| --- | --- | --- |
| Hofmann & Schieberle 1998 | thiol yields from each precursor system (isotope dilution): the main fit rows | [10.1021/jf9705983](https://doi.org/10.1021/jf9705983) |
| Zheng & Ho 1994 | H2S release from cysteine and glutathione: measured rate constants and activation energies | [10.1021/bk-1994-0564.ch011](https://doi.org/10.1021/bk-1994-0564.ch011) |
| Dwivedi & Arnold 1971 | share of degraded thiamine released as H2S, by pH | [10.1021/jf60177a033](https://doi.org/10.1021/jf60177a033) |
| Kang et al. 2026 | cysteine–xylose intermediates: fit rows | [10.1039/d5fb00932d](https://doi.org/10.1039/d5fb00932d) |
| Whitfield & Mottram 1999 | norfuraneol + cysteine or H2S: fit rows | [10.1021/jf980980v](https://doi.org/10.1021/jf980980v) |
| Kumazawa & Masuda 2003 | thiol loss against pH (fit) | [10.1021/jf021025f](https://doi.org/10.1021/jf021025f) |
| Zhang et al. 2024 | thiamine + methionine + xylose: disulfide and methylthio-furan formation (fit) | [10.1016/j.foodres.2024.114149](https://doi.org/10.1016/j.foodres.2024.114149) |

The thiamine H2S share is the least settled of these: Dwivedi & Arnold 1971 and their 35S sulfur balance of 1972
([10.1111/j.1365-2621.1972.tb03694.x](https://doi.org/10.1111/j.1365-2621.1972.tb03694.x)) used different
methods and concentrations, and near pH 6 and 120 °C they imply shares from about 1 % to over 30 %. Outputs that
depend on it there are refused ([STEP5_RECORD.md](../reference/lever_table/STEP5_RECORD.md), amendment 2026-11-17).

## 3. The expert list against the corpus

24 unique papers. **8 were already known to the corpus (33 %); 5 had been read.** By 2026-11-12 every paper on the
list had been read into an extraction record. None changed a rate constant: the measured values agreed with those
already in the model, or concerned something the model does not compute.

| status | papers |
| --- | --- |
| already used by the model | Cerny 2015 (claims; the chapter itself was not on disk, its statements were quoted through the claims panel); Cerny 2007 LWT (fit row, hold-outs); Cerny & Briffod 2007; Meynier & Mottram 1995 |
| data the model is now fitted on or tested against | Zhai et al. 2019; Cao et al. 2017; Zhang et al. 2025 (lysine); Zhao et al. 2019; Du et al. 2025 |
| directional claims the model is scored on | Modi, Linforth & Taylor 2008; Martins et al. 2010 |
| read; little or nothing taken | Wang, Yang & Song 2012 (glutathione–xylose; directional tests only); Hou et al. 2017; Liu et al. 2023 (LWT); Wang et al. 2026 (FRI); Zheng et al. 1997; Cerny & Guntz-Dubini 2013; Farmer, Mottram & Whitfield 1989; Cha, Debnath & Lee 2019; Zhai et al. 2023; Xu et al. 2013; Garbusov et al. 1976; Sun et al. 2022 (review); Sohail et al. 2022 (review) |

Cao et al. 2017 gave one new step (the cysteine–xylose thiazolidine equilibrium); it failed its acceptance test and
is switched off.

## 4. The three gaps Cerny 2007 (LWT) showed

Cerny 2007 (LWT) labels which carbons of each sulfur volatile come from xylose and which from thiamine.
Against it, the model's network has three gaps:

1. 3-mercapto-2-pentanone is made in the model only from thiamine; the paper shows >95 % of it from
   xylose when thiamine is absent. A sugar + H2S route is missing. **Open**: no shipped output depends on it yet,
   and its constant is confounded with parked routes (STEP5_RECORD, amendment 2026-11-17).
2. The methyl of 2-methyl-3-(methylthio)furan comes in the model only from thiamine; the paper shows
   64–70 % of it from xylose. **Open.**
3. At pH 5 without cysteine, thiamine supplies the sulfur of 2-furfurylthiol; the model's thiamine
   H2S release is zero below pH 6.25. **Partly addressed 2026-11-17**: the zero is measured only near 100 °C, so
   above 105 °C and below pH 6.25 thiamine-H2S-dependent outputs are now refused instead of set to zero. A fitted
   share waits for data that measure it directly.
