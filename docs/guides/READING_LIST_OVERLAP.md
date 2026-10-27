# The papers the model rests on, and the overlap with the WUR reading list

*Written 2026-10-27 after Bei Wang (Wageningen University) shared a reading list of 24 papers on volatile
formation from amino acids, sugars, thiamine and glutathione. Counts are the repository's own: rate
constants taken from a paper, fit rows a calibration was tuned on, benchmark pots the model is scored
against, and directional claims it is scored on (`build_thiol_sink_figures.paper_usage_counts`; every
paper and what is taken from it: [SOURCES.md](SOURCES.md)).*

## 1. The fifteen papers the model rests on most

| # | paper | records | what the model takes | on the WUR list |
| ---: | --- | ---: | --- | :---: |
| 1 | Hofmann & Schieberle 1998 (cysteine + sugars, isotope dilution) | 77 | 41 fit rows, 23 benchmark pots, 13 claims | |
| 2 | Kang et al. 2026 (cysteine–xylose, TTCA) | 14 | fit rows | |
| 3 | Whitfield & Mottram 1999 (norfuraneol + cysteine / H2S) | 14 | fit rows | |
| 4 | **Cerny 2015** (sulfur chemistry in thermal aroma generation, book chapter) | 11 | directional claims | **yes** |
| 5 | Zhou et al. 2023 (Amadori compound + cysteine, pH drift) | 10 | fit rows | |
| 6 | Martins & van Boekel 2003 (Amadori degradation kinetics) | 10 | rate constants | |
| 7 | Martins & van Boekel 2005 (glucose–glycine kinetics) | 10 | rate constants | |
| 8 | Trikusuma et al. 2018 (pea protein, lipid volatiles) | 9 | claims | |
| 9 | Kocadağlı & Gökmen 2016 (sugar degradation) | 9 | rate constants | |
| 10 | De Vleeschouwer et al. 2009 (acrylamide kinetics) | 8 | constants, fit rows | |
| 11 | Yiltirak et al. 2026 (ribose–cysteine, temperature ladder) | 8 | hold-out pots, claims | |
| 12 | Mottram & Nobrega 2002 (buffered vs unbuffered cysteine–ribose) | 8 | claims | |
| 13 | Ma et al. 2024 | 7 | claims | |
| 14 | Zhang et al. 2024 | 6 | fit rows, claims | |
| 15 | Ni et al. 2021 | 6 | claims | |

Further down the same ranking, also on the WUR list: **Cerny 2007 (LWT, origin of carbons)** — a fit row (the
54:46 thiamine:xylose carbon split of MFT) and two hold-outs since wave B2; **Cerny & Briffod 2007 (JAFC,
pH)** — a claim and fit-system conditions; **Meynier & Mottram 1995** — the conditions of the two pH claims
PH-01/02.

## 2. The WUR list against the corpus

24 unique papers. **8 were already known to the corpus (33 %); 5 had been read.** Two of the eight were on
disk without an extraction record and were read on 2026-10-27 (Wang, Yang & Song 2012; Cerny 2007 LWT).

| status | papers |
| --- | --- |
| already used by the model | Cerny 2015 (claims; the chapter itself was not on disk, its statements were quoted through the claims panel); Cerny 2007 LWT (fit row, hold-outs); Cerny & Briffod 2007; Meynier & Mottram 1995 |
| read, little or nothing taken | Hou et al. 2017; Liu et al. 2023 (LWT); Wang et al. 2026 (FRI) |
| read on 2026-10-27 | Wang, Yang & Song 2012 (glutathione–xylose): directional tests only, no rates; Cerny 2007 LWT: confirms the rows in use and shows three gaps in the network (below) |
| new, being added | Zheng et al. 1997; Cerny & Guntz-Dubini 2013; Modi, Linforth & Taylor 2008; Farmer, Mottram & Whitfield 1989; Martins et al. 2010; Cha, Debnath & Lee 2019; Cao et al. 2017; Zhao et al. 2019; Zhai et al. 2019; Zhai et al. 2023; Zhang et al. 2025 (lysine); Du et al. 2025; Xu et al. 2013; Sun et al. 2022 (review); Sohail et al. 2022 (review); Garbusov et al. 1976 |

## 3. What the first two reads already changed

Cerny 2007 (LWT) labels which carbons of each sulfur volatile come from xylose and which from thiamine.
Against it, the model's network has three gaps, now in the backlog:

1. 3-mercapto-2-pentanone is made in the model only from thiamine; the paper shows >95 % of it from
   xylose when thiamine is absent. A sugar + H2S route is missing.
2. The methyl of 2-methyl-3-(methylthio)furan comes in the model only from thiamine; the paper shows
   64–70 % of it from xylose.
3. At pH 5 without cysteine, thiamine supplies the sulfur of 2-furfurylthiol; the model's thiamine
   H2S release is zero below pH 6.25.

## 4. Next

Every paper in section 2's last row is read into an extraction record, in the order the backlog gives
(quantitative papers that could anchor or test a step first; reviews last), and what each can anchor or
test enters the model through the usual pre-registered waves.
