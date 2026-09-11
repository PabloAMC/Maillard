# Yu, Meng, Ramaswamy & Boye 2014 — EXTRACTION (measured residence time distribution of a soy-protein-isolate feed in a twin-screw extruder)

**Source on disk:** `data/articles/yu2012.pdf` (1.0 MB; downloaded 2026-09-11 at this repository's
request). Read 2026-09-11 via `pdftotext -layout`. Wave B37.

| field | value |
|---|---|
| Title | "Residence time distribution of soy protein isolate and corn flour feed mix in a twin-screw extruder" |
| Venue | Journal of Food Processing and Preservation (Wiley); the repository's reading list carried it as 2014 |
| DOI | 10.1111/jfpp.12005 |
| Group | McGill University (Ramaswamy) with Agriculture and Agri-Food Canada (Boye) |
| System | **corn flour : soy protein isolate 4:1 (20 % SPI)** in a co-rotating twin-screw extruder |
| Design | full factorial: screw speed **75, 100, 125 rpm** × feed moisture **25, 30, 35 % (w/w)** × die diameter **3 and 5 mm**; tracer method, samples collected at 10-s intervals |
| What is measured | first passage residence time (FPRT), exit completion time (ECT), **mean residence time tm (s)** and its variance tv |

## 1. Why this paper matters to this model

`data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json` — the model's only extrusion acrylamide
row, and a FIT row — declares a hold of **25 s**. Wave B36 read its source (Ma et al. 2024) in full
and found that **the paper prints no residence time at all**: the 25 s is an assumption with no
source, now labelled one in the bundle's vessel note. This paper is the nearest measured residence
time distribution for a soy-protein-isolate feed in a twin-screw extruder.

## 2. Table 1 — the full design and result, verbatim ("EXPERIMENT DESIGN AND RESULTS")

| run | feed moisture (w/w %) | screw speed (rpm) | die (mm) | FPRT (s) | ECT (s) | **tm (s)** | tv |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 30 | 75 | 5 | 30 | 100 | **67.4** | 33.4 |
| 2 | 30 | 100 | 5 | 20 | 80 | **48.2** | 20.9 |
| 3 | 30 | 125 | 5 | 10 | 70 | **39.9** | 13.4 |
| 4 | 30 | 75 | 3 | 30 | 120 | **76.0** | 42.6 |
| 5 | 30 | 100 | 3 | 20 | 100 | **55.3** | 26.7 |
| 6 | 30 | 125 | 3 | 10 | 80 | **44.0** | 16.6 |
| 7 | 35 | 75 | 3 | 30 | 90 | 61.5 | 26.6 |
| 8 | 35 | 100 | 3 | 20 | 80 | 45.7 | 17.6 |
| 9 | 35 | 125 | 3 | 10 | 60 | 35.0 | 9.79 |
| 10 | 25 | 75 | 3 | 40 | 100 | 86.6 | 52.3 |
| 11 | 25 | 100 | 3 | 20 | 110 | 64.05 | 42.00 |
| 12 | 25 | 125 | 3 | 10 | 80 | 54.15 | 25.28 |
| 13 | 35 | 75 | 5 | 30 | 80 | 61.14 | 16.46 |
| 14 | 35 | 100 | 5 | 20 | 60 | 41.64 | 13.50 |
| 15 | 35 | 125 | 5 | 15 | 70 | 34.81 | 12.40 |
| 16 | 25 | 75 | 5 | 30 | 100 | 77.40 | 34.82 |
| 17 | 25 | 100 | 5 | 20 | 100 | 59.94 | 31.35 |
| 18 | 25 | 125 | 5 | 15 | 90 | 44.75 | 21.42 |

Verbatim on the trends: *"higher screw speed, higher initial moisture content and larger die
diameter resulted in a shorter mean residence time"*, and *"increasing screw speed from 75 to
125 rpm resulted in decreasing tm from 87 to 54 s at a feed moisture level of 25%."*

## 3. The rows that bear on the model's extrusion pot

Ma et al. 2024's acrylamide series ran at **30 % moisture and 150 rpm**. The two rows at 30 %
moisture and the highest speed measured here are **runs 3 and 6: tm = 39.9 s (5 mm die) and 44.0 s
(3 mm die) at 125 rpm**. The trend with speed is monotone and steep (67.4 → 48.2 → 39.9 s for the
5 mm die), so 150 rpm would sit below 39.9 s — but that is off the end of the measured range and this
dossier does not extrapolate it.

**The comparison is nevertheless not like for like, and the difference is large:**

| | Yu et al. | Ma et al. 2024 (the model's row) |
|---|---|---|
| composition | corn flour : SPI **4:1** (20 % SPI) | SPI : corn starch **9:1** (90 % SPI) |
| moisture | 25–35 % | 30 % |
| screw speed | 75–125 rpm | 150 rpm |
| what tm covers | the whole extruder, feed to exit | — |

A 20 %-protein starch feed and a 90 %-protein feed are different melts with different viscosities,
so the measured tm does not transfer as a number. What it establishes robustly is the **order of
magnitude and the direction**: a twin-screw extruder at this scale holds material for **tens of
seconds, and at 30 % moisture and the fastest screw measured, still 40 s** — against the bundle's
unsourced 25 s, and covering the whole barrel rather than one zone.

## 4. Verdict

The measured anchor for a condition the model currently assumes. It does not license replacing the
bundle's 25 s with a number from a different feed — but it does mean the assumption can now be
priced rather than merely flagged. That pricing is a pre-registered question for a later wave,
because the row is a FIT row whose report is frozen; nothing was changed in B37.
