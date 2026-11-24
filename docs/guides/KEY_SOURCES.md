# The papers the model rests on

*2026-11-24. The repository has read about 350 papers, but most give nothing quantitative (the full list, and
what is taken from each: [SOURCES.md](SOURCES.md)). This page lists the few that carry the model. Counts are the
repository's own: rate constants taken from a paper, fit rows a calibration was tuned on, benchmark pots the model
is scored against, and directional claims it is scored on (`build_thiol_sink_figures.paper_usage_counts`).*

## 1. The fifteen papers the model rests on most

Across the whole model (sugar breakdown, sulfur, acrylamide and lipid chemistry).

| # | paper | records | what the model takes | DOI |
| ---: | --- | ---: | --- | --- |
| 1 | Hofmann & Schieberle 1998 (cysteine + sugars, isotope dilution) | 77 | 41 fit rows, 23 benchmark pots, 13 claims | [10.1021/jf9705983](https://doi.org/10.1021/jf9705983) |
| 2 | Zhou et al. 2023 (Amadori compound + cysteine or glutathione, pH drift) | 16 | 10 fit rows, 4 benchmark pots, 2 claims | [10.1021/acs.jafc.3c04424](https://doi.org/10.1021/acs.jafc.3c04424) |
| 3 | Kang et al. 2026 (cysteine–xylose intermediates) | 14 | fit rows | [10.1039/d5fb00932d](https://doi.org/10.1039/d5fb00932d) |
| 4 | Whitfield & Mottram 1999 (norfuraneol + cysteine or H2S) | 14 | fit rows | [10.1021/jf980980v](https://doi.org/10.1021/jf980980v) |
| 5 | Cerny 2015 (sulfur chemistry in thermal aroma generation, book chapter) | 11 | directional claims | [10.1016/b978-1-78242-103-0.00009-6](https://doi.org/10.1016/b978-1-78242-103-0.00009-6) |
| 6 | Martins & van Boekel 2003 (Amadori degradation kinetics) | 10 | rate constants | [10.1016/s0008-6215(03)00174-5](https://doi.org/10.1016/s0008-6215(03)00174-5) |
| 7 | Martins & van Boekel 2005 (glucose–glycine kinetics) | 10 | rate constants | [10.1016/j.foodchem.2004.08.013](https://doi.org/10.1016/j.foodchem.2004.08.013) |
| 8 | Kocadağlı & Gökmen 2016 (sugar degradation) | 9 | rate constants | [10.1016/j.foodchem.2016.05.150](https://doi.org/10.1016/j.foodchem.2016.05.150) |
| 9 | Trikusuma, Paravisini & Peterson 2020 (pea protein beverages, aroma; the claims cite the 2018 thesis behind it) | 9 | claims | [10.1016/j.foodchem.2019.126082](https://doi.org/10.1016/j.foodchem.2019.126082) |
| 10 | Mottram & Nobrega 2002 (cysteine with three forms of ribose) | 8 | claims | [10.1021/jf0200826](https://doi.org/10.1021/jf0200826) |
| 11 | De Vleeschouwer et al. 2009, Part II (acrylamide kinetics) | 8 | constants, fit rows | [10.1016/j.foodchem.2008.09.084](https://doi.org/10.1016/j.foodchem.2008.09.084) |
| 12 | Yiltirak et al. 2026 (ribose–cysteine, temperature ladder) | 8 | hold-out pots, fit rows, claims | [10.1016/j.foodres.2026.118600](https://doi.org/10.1016/j.foodres.2026.118600) |
| 13 | Ma et al. 2024 (extrusion, CML, CEL and acrylamide) | 7 | claims | [10.3390/ijms25168668](https://doi.org/10.3390/ijms25168668) |
| 14 | Zhang et al. 2024 (thiamine + methionine + xylose, sulfide radicals) | 6 | fit rows, claims | [10.1016/j.foodres.2024.114149](https://doi.org/10.1016/j.foodres.2024.114149) |
| 15 | Ni et al. 2021 (sugar type in hydrolysate Maillard products) | 6 | claims | [10.1016/j.fochx.2021.100127](https://doi.org/10.1016/j.fochx.2021.100127) |

The counter keys papers by first author and year, so two papers sharing both are counted together; each row above was
resolved to the paper the records actually cite.

## 2. The papers behind the thiol forecasts

The ranking above counts records across the whole model. For the meaty thiols (2-methyl-3-furanthiol,
2-furfurylthiol and their disulfides) from cysteine, thiamine, glutathione and a pentose, these are the papers the
sulfur lane's constants come from. Several are stored as module constants in
`src/kinetic_core/parameters_sulfur.py`, which the counter in section 1 does not see, so this list is curated
from that file. Routes that are off by default are not listed.

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

## 3. Known gaps in the thiol network

Cerny 2007 (LWT, [10.1016/j.lwt.2006.09.008](https://doi.org/10.1016/j.lwt.2006.09.008)) labels which carbons of
each sulfur volatile come from xylose and which from thiamine. Against it, the model's network has three gaps:

1. 3-mercapto-2-pentanone is made in the model only from thiamine; the paper shows >95 % of it from
   xylose when thiamine is absent. A sugar + H2S route is missing. **Open**: no shipped output depends on it yet,
   and its constant is confounded with parked routes (STEP5_RECORD, amendment 2026-11-17).
2. The methyl of 2-methyl-3-(methylthio)furan comes in the model only from thiamine; the paper shows
   64–70 % of it from xylose. **Open.**
3. At pH 5 without cysteine, thiamine supplies the sulfur of 2-furfurylthiol; the model's thiamine
   H2S release is zero below pH 6.25. **Partly addressed 2026-11-17**: the zero is measured only near 100 °C, so
   above 105 °C and below pH 6.25 thiamine-H2S-dependent outputs are now refused instead of set to zero. A fitted
   share waits for data that measure it directly.
