# Gobert & Glomb 2009 — EXTRACTION (the α-dicarbonyls of glucose degradation with lysine at 50 °C, and two fed half-lives: 3-deoxyglucosone 40 h, glucosone 8 h)

**Source on disk:** `data/articles/Gobert2009.pdf` (0.95 MB; downloaded 2026-09-11 at this
repository's request). Read 2026-09-11 via `pdftotext -layout`. Wave B42.

| field | value |
|---|---|
| Title | "Degradation of Glucose: Reinvestigation of Reactive α-Dicarbonyl Compounds" |
| Venue | Journal of Agricultural and Food Chemistry 2009, 57, 8591–8597 |
| DOI | 10.1021/jf9019085 — read from the printed header |
| Group | Martin Luther University Halle-Wittenberg (Glomb) |
| System | **glucose 42 mM + lysine 42 mM in 0.1 M phosphate, pH 7.4, 50 °C, 7 days**, aerated and deaerated; dicarbonyls trapped with o-phenylenediamine (5 mM) either during incubation or added at sampling; quinoxalines by HPLC-UV against isolated, NMR-verified standards; ¹³C-labelled glucose for the carbon origin |
| Fed experiments | **glucosone (42 mM) and 3-deoxyglucosone (42 mM) each incubated with lysine (42 mM) at 50 °C for 0–48 h**, aerated and deaerated, OPD added at sampling |

## 1. Why it matters here

The B37–B41 finding was that the model destroys the 3-deoxy pool too fast, and that its aqueous
glucosone → glyoxal barrier (B21) is "consistent with zero". This paper feeds both compounds and
prints their **half-lives at 50 °C** — a temperature 70 °C below the model's window, which is
exactly where a wrong barrier shows.

## 2. Table 1 — quinoxalines after 7 days from glucose + lysine (+ OPD), mmol per mol glucose

| dicarbonyl | aerated | deaerated |
|---|---:|---:|
| glucosone | 44.5 | 5.5 |
| 1-deoxyglucosone | 6.9 | 23.5 |
| 3-deoxyglucosone | 4.6 | 5.0 |
| Lederer's glucosone | 4.1 | 8.2 |
| 1-deoxypentosone | 1.1 | 0.6 |
| 3-deoxypentosone | 0.5 | 0.2 |
| 1-deoxythreosone | 4.4 | 14.8 |
| 3-deoxythreosone | 1.5 | 4.6 |
| threosone | 8.4 | 2.0 |
| methylglyoxal | 3.9 | 4.0 |
| glyoxal | 3.3 | 1.6 |

With OPD added only at sampling (no trapping during the week) the same products appear "at much
lower concentrations": 3-deoxyglucosone 2.2 mmol/mol at 7 days, "the major product ... and
increased independent from oxygen"; glucosone preferentially under air; 1-deoxyglucosone 0.05 vs
0.2 (deaerated); Lederer's glucosone 0.04 vs 0.2; methylglyoxal 0.08 vs 0.08; glyoxal 0.12 vs 0.06;
threosone 0.25 vs 0.08 (Figs. 3–4, figure-only beyond these prose values). Lysine fell by 30 % in
7 days.

## 3. The two fed half-lives, verbatim

> "In this experiment **3-deoxyglucosone had a half-life of 40 h**. However no additional
> quinoxalines were formed. In contrast from glucosone, **with a half-life of 8 h**, the formation of
> 1-deoxypentosone-Q, 3-deoxypentosone-Q, pentosone-Q, threosone-Q, methylglyoxal-Q and glyoxal-Q
> was established (Table 3)."

Table 3 — degradation of glucosone (42 mM, with lysine, 50 °C, 8 h), mmol per mol glucosone:

| product | aerated | deaerated |
|---|---:|---:|
| threosone | 7.5 | 3.5 |
| 3-deoxypentosone | 2.9 | 6.1 |
| 1-deoxypentosone | 0.3 | 0.9 |
| methylglyoxal | 0.3 | 0.3 |
| glyoxal | 0.7 | 0.7 |

So in 8 h at 50 °C half the glucosone is gone and **about 1 % of it is found as any smaller
dicarbonyl**; glyoxal is 0.07 %. The rest went to the amine (lysine) or to products the assay does
not see. And fed 3-DG loses half in 40 h **without making any smaller dicarbonyl at all**: at 50 °C
with an amine, the 3-DG pool leaves by the amine route, not by fragmentation.

## 4. What the repository can take

Two fed rows at 50 °C, pH 7.4, with an amine present — under the owner's rule, fed yields are FIT
evidence, and these are also the lowest-temperature aqueous measurements of these two pools on disk:

- **3-DG half-life 40 h** with 42 mM lysine. The model's 3-DG exits at 50 °C are Martins' formic-acid
  step (Ea 30 kJ/mol, which the parameter table already flags as conflicting with Knol 2010's 84)
  and Kocadagli's methylglyoxal step; the pH term B41 installed raises the formic-acid exit at 7.4.
- **glucosone half-life 8 h**, and glyoxal 0.07 % of it. The model's aqueous glucosone → glyoxal
  constant (B21) carries a barrier of 4.2 kJ/mol, "consistent with zero", fitted at 110–140 °C in
  milk; a near-zero barrier extrapolated to 50 °C keeps that step fast.

Both are probed in `results/validation/kinetic_core_b42_prereg.md` with the prediction written
before the run. Limit: 50 °C with lysine is a glycation pot, not a cook; what transfers is the
temperature dependence, not the level.
