# Zamora, Navarro, Aguilar & Hidalgo 2015 — EXTRACTION (thermal breakage of 2-alkenals and 2,4-alkadienals: 2,4-decadienal makes hexanal, with yields and barriers)

**Source on disk:** `data/articles/Zamora2015b.pdf` (0.74 MB; downloaded 2026-09-11 at this repository's
request). Read 2026-09-11 via `pdftotext -layout`. Wave B42. Not the `zamora2015.pdf` already on disk
(Strecker aldehydes versus amines), a different paper by the same group.

| field | value |
|---|---|
| Title | "Lipid-derived aldehyde degradation under thermal conditions" |
| Venue | Food Chemistry 174 (2015) 89–96 |
| DOI | 10.1016/j.foodchem.2014.11.034 — read from the printed footer |
| Group | Instituto de la Grasa, CSIC, Seville |
| Systems | 2-pentenal, 2-octenal, 2,4-heptadienal and 2,4-decadienal (0–80 µmol; 4 µmol standard) in 80 µL tetrahydrofuran + 420 µL of 0.2 M buffer (pH 2.15–11), closed tube under air, **120, 160 or 200 °C**, up to 60 min; products by GC-MS and LC-MS/MS after derivatisation, against seven-level standard curves |

## 1. What it is, and what it is not

It was fetched as a candidate law for the **loss** of hexanal that four of the sixteen B37 papers show
and the lipid lane cannot represent. **It is not that.** It is the thermal breakage of the
carbon–carbon double bonds of *unsaturated* aldehydes: 2-alkenals give the alkanal two carbons
shorter, 2,4-alkadienals give the alkanal and the 2-alkenal. The saturated aldehydes are the
**products**, not the substrates. For this model the paper is a hexanal **source** the lipid lane
lacks — 2,4-decadienal, which the lane already makes and a B35 bundle already scores, breaks to
hexanal at cooking temperature.

## 2. Numbers printed in the prose (all at pH 8, 1 h, 200 °C unless stated)

| reaction | yield | where |
|---|---:|---|
| 2-pentenal → propanal | **12.5 %** (slope 0.125, r = 0.993) | §3.2 |
| 2-octenal → hexanal | **18.0 %** (slope 0.180, r = 0.999), constant over 0–80 µmol | §3.3 |
| 2,4-heptadienal → propanal / 2-pentenal | 9.8 % / 1.0 % | §3.4 |
| **2,4-decadienal → hexanal / 2-octenal** | **11.5 % / 0.8 %** (slopes 0.1154 / 0.00821, r > 0.998) | §3.5 |

Barriers, from Arrhenius plots of the initial linear formation rates at 120/160/200 °C: propanal
from 2-pentenal **25.2 kJ/mol**; hexanal from 2,4-decadienal **21.3 kJ/mol**, 2-octenal from
2,4-decadienal **29.6 kJ/mol**; "always very similar and about 25 kJ/mol" for all four aldehydes.

Time scales, verbatim: 2-pentenal "less than 10 % ... after 25 min at 200 °C and after 45 min at
160 °C. When 2-pentenal was heated at 120 °C, 17 % of the initial aldehyde was still present after
60 min"; 2-octenal "less than 10 % ... after 10 min heating at 200 °C, 50 min heating at 160 °C, and
about 60 min when heating at 120 °C". pH: hexanal from 2-octenal is maximal around pH 10; hexanal
and 2-octenal from 2,4-decadienal peak "at about pH 8", the alkenal slightly more acidic.

The rate constants themselves and the disappearance curves are **figure-only** (Figs. 2–5); the
yields, barriers and half-life brackets above are the printed numbers.

## 3. What the repository can take, and the limit

- A **route the lipid lane lacks**: 2,4-decadienal → hexanal (11.5 %) + 2-octenal (0.8 %), with a low
  barrier (21 kJ/mol) and a half-life of tens of minutes at 120–200 °C. On Trikusuma's UHT pot (140 °C,
  6 s) it is negligible; on a 30-minute roast it is not. A candidate step, not built here.
- A **directional fact** the lane's decadienal rows should carry: the compound is unstable on the
  cook's own time scale, so a measured 2,4-decadienal level is a net of formation and breakage.
- **Not a hexanal loss law.** The loss the sixteen papers show is a different question — enzymatic
  formation switched off, stripping, and binding — and this paper does not speak to it.

Limit: tetrahydrofuran–buffer at 4 µmol in 0.5 mL, closed tube; the medium is not a food.
