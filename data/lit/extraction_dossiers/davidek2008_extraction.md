# Davidek, Gouezec, Devaud & Blank 2008 — EXTRACTION
### *"Origin and Yields of Acetic Acid in Pentose-based Maillard Reaction Systems"*
### Xylose + glycine (0.1 M each), 0.2 M phosphate, pH 6.0 / 8.0, 120 C 4 h (and 90 C 15 h, labelling only); acetic acid by isotope dilution, residual xylose

**Read 2026-10-30 from `data/articles/davidek2008.pdf`** (3 pp.) by two independent readers for the keyed table
`data/lit/fit_rows/ph_acid_davidek_acids.yml` (ledger `davidek2008_fig1_xylose_gly_acetic_residual_xylose`; both
readings agree, reconciliation `results/validation/source_reading_reconciliation_ph_acid_davidek_acids.md`). This
dossier records what the paper is; the numbers live in the keyed table.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Tomas Davidek, Elisabeth Gouezec (Nestle Product Technology Centre Orbe), Stephanie Devaud (Nestle Research Center), Imre Blank |
| Venue | *Ann. N.Y. Acad. Sci.* **2008**, 1126, 241-243 |
| DOI | 10.1196/annals.1433.053 (printed) |

## 1. What it measures

- **Model reactions:** glycine 0.35 mmol + D-xylose 0.35 mmol in 3.5 mL phosphate buffer (0.2 mol/L), 0.4 mL per
  1.5 mL screw-cap vial, silicone bath, 120 C for 4 h or 90 C for 15 h, pH 6.0 and 8.0 (initial). Unlabelled and
  [1-13C] / [2-13C] / [5-13C]-xylose.
- **Acetic acid:** SPME-GC/MS with [2H3]acetic acid as internal standard (isotope dilution) [M].
- **Fig. 1:** acetic acid and residual xylose vs time at 120 C, pH 6 (squares) and pH 8 (triangles), in mol% of the
  initial xylose [F; digitised by both readers]. The abstract's 16 / 28 mol% are its 240-min points.
- **Labelling:** 77-89 % of acetic acid from C-1 / C-2 (abstract 77-87 %, text 85-89 %), 8-15 % from C-4 / C-5; the
  authors assign it to the hydrolytic beta-dicarbonyl cleavage of 1-deoxypento-2,4-diulose [Q].
- "Formation almost stopped when no more xylose was available" [Q]: acetic acid plateaus once xylose is gone.

## 2. Flags

- [!] How residual xylose was measured is never stated.
- [!] Methods say 90 C for 15 h; the Fig. 4 legend says 16 h (labelling data only, no amounts).
- [!] At pH 8, xylose at 60 - 240 min sits on the x-axis (censored, not keyed).
- [!] CF-069: the plateau caps acetic-acid loss far below the trunk's fitted `r_aa_frag`.
