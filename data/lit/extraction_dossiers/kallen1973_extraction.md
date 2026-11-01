# Kallen 1973 (Additions and Corrections to Kallen 1971b) — EXTRACTION
DOI: 10.1021/ja00800a602 (not printed; Crossref 2026-11-01: "Additions and Corrections - The Mechanism of Reactions Involving Schiff Base Intermediates. Thiazolidine Formation from L-Cysteine and Fromaldehyde [sic]", Kallen, J. Am. Chem. Soc. 95(19), 6510, 1973-09)
### A five-line erratum: it replaces the printed steady-state rate equations 9 and 10 of Kallen 1971b (JACS 93, 6236). No table, figure, constant or data point is changed. **The k_r and K used in `kallen1971b_extraction.md` and CF-076 are unaffected.**

**Read 2026-11-01 from `data/articles/kallen1973.pdf`** (1 p., JACS 95, p. 6510; Acrobat Capture OCR). The page
also carries the end of an unrelated nootkatone communication (Dastur) and two other corrections (Stewart &
MacPhee; Reynolds et al.). Both equations were read from a 200 dpi render (scratchpad `ing1101/ka-1.png`) and
compared with the original eqs 9-10 on JACS 93 p. 6244 (render of `kallen1971b.pdf` p. 9, `kb-09.png`).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** printed · **[Z]** derived by me · **[!]** flag.

**Laboratory.** `kallen_pennsylvania` (lab_ids.yml), University of Pennsylvania (Kallen).

## 1. What it corrects, exactly

The erratum's only sentence: "Equations 9 and 10 should read as shown below." [V]

| | as printed in 1971b (p. 6244) | as corrected (1973) |
|---|---|---|
| eq 9 | k_obsd = (k1a α3 + k1b α13)[F] / { (k1a α3/K1 + k1b α13/K1) / (k2 + k2′a_H+) + 1 } | k_obsd = (k1a α3 + k1b α13)[F] **α_RNH2T / (α3 + α13)** / { (k1a α3/K1 + k1b α13/K1) / **((α3 + α13)**(k2 + k2′a_H+)**)** + 1 } |
| eq 10 | the same with k2 + k2′a_H+ + k2″[HA] | the same correction, with k2 + k2′a_H+ + k2″[HA] |

Two changes, both normalisations by the free-base amine fraction [M]:
1. **Numerator** multiplied by α_RNH2T / (α3 + α13) (total amine-free-base fraction over the fraction in the two
   reacting species 3 and 13).
2. **Partition term** in the denominator divided by (α3 + α13): the reverse attack k−1 (= k1/K1) is now averaged
   over the free-base species, not weighted by their fraction of TOTAL cysteine.

## 2. Does it change the k_r / K in `kallen1971b_extraction.md` and CF-076?

**No** (confidence ~95 %).
- The dossier's k_r(pH) = k_f(pH) / K_app(pH) uses **k_f from Table I (printed k_obsd/[F]) and Fig. 6 data points**,
  and **K_app from Kallen 1971 Table IV**. Neither is computed from eq 9 or 10 [Z, checked against the dossier's
  sec. 4 derivation].
- Tables I-V, Figs. 3, 6, 7 and the constants of Table II are untouched by the erratum; it prints no new value
  [M]. The fitted lines and constants were presumably computed with the correct form (the erratum corrects the
  print, not the analysis); the paper's text agrees with the corrected form (eq 8a carries α_RNH2T) [Z].
- CF-076 uses Kallen 1971b's Fig. 7 saturation through an empirical ceiling (`F_cat = 1 + rho x / (rho + x)`,
  cluster-5 prereg sec. 2), not eq 10 itself.

**What the correction does mean [Z]:** only the corrected eq 10 produces the phosphate saturation Fig. 7 shows. At
pH 6.3, with Table II's constants (k−1 ≈ 42-105 s⁻¹, k2′ 1.4 × 10⁸ M⁻¹ s⁻¹, phosphate k2″ 1.2 × 10⁴ M⁻¹ s⁻¹), the
corrected partition term is ~0.6-1.5 without buffer and ~0.02-0.06 at 0.15 M phosphate monoanion: a ceiling of
~x1.6-2.4 at that pH, the order of Fig. 7's x3.5. The uncorrected print multiplies that term by α3 + α13
(~2.5 × 10⁻³ at pH 6.3), i.e. makes it ~400x smaller: no attack limitation, so no saturation. Anyone who implemented eq 9 / 10 as printed in 1971 would have lost the change of
rate-determining step; nothing in the repo does (grep of docs/history/cys_exchange_2026-10-31 and
conflict_resolution_2026-11-01/cluster5).

## 3. Bearing on open questions

| item | verdict |
|---|---|
| CF-076 (k_r band 2-8 × 10⁻⁴ min⁻¹ at 100 °C, Kallen-shaped K(pH), k_r(pH)) | LEAVES UNCHANGED |
| `kallen1971b_extraction.md` flag 7.6 ("the 1973 Additions and Corrections may revise any of the above") | **resolved: it revises only eqs 9-10; every number in that dossier stands** |
| saturation of buffer catalysis (E1c rejection) | consistent: the corrected eq 10 is what makes it saturate |

## 4. Disagreements found

None with other dossiers. The erratum's Crossref title misspells "Formaldehyde" ("Fromaldehyde"); the printed
page spells it correctly [!].

## 5. Ledger-ready measurements

None (no data). Note beside the kallen1971b ledger rows: eqs 9-10 as corrected in JACS 95, 6510.
