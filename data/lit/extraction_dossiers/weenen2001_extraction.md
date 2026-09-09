# Weenen 2001 — EXTRACTION (phenylalanine + fifteen carbonyl partners, ~26 mmol/L each in 50 mmol/L citrate pH 3.2, reflux 2 h under Likens-Nickerson extraction; Glc-Phe ARP at pH 7.0 and 3.2; Glc-Pro and glyceraldehyde-Pro ARPs at pH 7.2)
### A single-temperature yield ranking of Strecker co-reactants (pyruvaldehyde 21 %, glyoxal 8.8 %, 3-deoxyglucosone 0.6 %, glucose 0.03 %) for phenylacetaldehyde; no rate, no barrier, no time series.

**Source on disk:** `data/articles/weenen2001.pdf` (13 pp., ACS Symposium Series 794, chapter 15, pp.
183-195; owner's download, 2026-09-08). Read from the text layer (`scratchpad/articles/weenen2001.txt`)
and re-extracted with `pypdf` for pages 8-9 (Tables II and III), which came through cleanly. Figures 1-4
are reaction schemes only. Table I is an odour-threshold list (secondary). No value was read off any
figure.

## 0. Identity

| field | value |
|---|---|
| Title | "The Formation of Strecker Aldehydes" |
| Authors | H. Weenen, J. G. M. van der Ven (Bio-Organic Chemistry Section, Quest International, Bussum, NL) |
| Venue | ACS Symposium Series 794, "Aroma Active Compounds in Foods" (Takeoka, Güntert, Engel, eds.), chapter 15, pp. 183-195; published 14 Aug 2001 |
| DOI | 10.1021/bk-2001-0794.ch015 |
| Naming | ARP = Amadori rearrangement product; Glc-Phe ARP = N-(1-deoxy-D-fructos-1-yl)-phenylalanine; Glc-Pro ARP, Glc-Val ARP likewise; pyruvaldehyde = methylglyoxal; 3-deoxyerythrosone / 3-deoxyxylosone / 3-deoxyglucosone = the C4 / C5 / C6 3-deoxyosones; ACP = 2-acetyl-1-pyrroline; ACTP = 2-acetyl-3,4,5,6-tetrahydropyridine / 2-acetyl-1,4,5,6-tetrahydropyridine; "yield" = per cent, basis not defined (read as mol % on the amino acid or the ARP) |
| Quirks | Experimental says citrate "pH 3.5"; Tables II-III and the discussion say "pH 3.2". Reference numbering: the text cites "Fujimaki et al. (27)" but ref 27 is Kerler 1997; the Fujimaki 1968 reactivity paper is ref 26 |

## 1. Why it matters

Programme 6 needs the Strecker partner named and ranked. This chapter ranks fifteen carbonyl
partners for one amino acid (phenylalanine) under one condition, with pyruvaldehyde (methylglyoxal)
at 21 % and glyoxal at 8.8 % against 3-deoxyglucosone at 0.6 % and glucose at 0.03 %. Those are
within-study ratios the trunk can be checked against: the model's relative use of methylglyoxal,
glyoxal and the deoxyosone as Strecker donors at low pH. It also gives the ARP-vs-free-precursor
comparison at two pH values (the Glc-Phe ARP is 4 times glucose + Phe at pH 7.0 and 270 times at pH
3.2), and the proline branch (ACP, ACTP from the Glc-Pro ARP). What it cannot give is a rate: one
temperature (reflux), one time (2 h), and the aldehyde is continuously stripped by the distillation.

## 2. Methods as they matter to a model

- **Phenylacetaldehyde runs (verbatim).** "The carbonyl containing compound (~1.3 mmol) and an
  equimolar amount of phenylalanine, were dissolved in citrate buffer (50 mL, 50 mM, pH 3.5), and
  continuously extracted for 2 h with dichloromethane to which undecane was added as the internal
  standard. After 2 h, heating of the flask containing the aqueous phase was stopped and extraction
  with dichloromethane continued for another 20 min." So **~1.3 mmol / 50 mL = 26 mmol/L** carbonyl
  and 26 mmol/L phenylalanine, 50 mmol/L citrate; Likens-Nickerson simultaneous distillation-extraction
  (the aqueous flask refluxes at ~100 C; the aldehyde is removed into dichloromethane as it forms).
  The pH 7.0 entries of Table II used 100 mM phosphate (volume not stated).
- **Proline / ARP runs.** ARP 150 µmol-2.7 mmol, or the parent sugar + amino acid at the same molar
  amounts, in 35 mL 50 mM phosphate pH 7.2 (or 50 mL 50 mM citrate pH 3.5), 4 h extraction + 20 min;
  IS 2-acetylpyrazine; "Quantitation was based on comparison with external standard consisting of
  solution of ACTP, ACP and acetylpyrazine. Concentration-response correlation was linear in the
  concentration range in which measurements took place."
- **Ethanol run.** ~1 mmol Glc-Phe ARP + KH2PO4 170 mg (1.25 mmol) + K2HPO4·3H2O 285 mg (1.25 mmol) in
  ethanol, reflux 2 h; IS 2-acetylpyrazine.
- **GC.** Carlo Erba HRGC 5300, FID, HP-5 50 m × 0.32 mm × 1.05 µm, 75 → 150 C at 3 C/min, → 300 C at
  40 C/min, 5 min hold; injector 125 C, detector 210 C. For phenylacetaldehyde only the internal
  standard (undecane) is named; a response factor or external standard for PAC is not stated.
- **ARP syntheses.** Glc-Pro after Vernin 1992; glyceraldehyde-Pro after Huyghues-Despointes &
  Yaylayan 1996; Glc-Phe after Sosnovsky 1993; Glc-Val after Xenakis 1983. Purities not stated.
- **Replicates.** Not stated; no error bars or SD anywhere.
- **Yield definition.** "Yield of phenylacetaldehyde" in per cent; the basis (mol PAC per mol Phe, or
  per mol carbonyl; the two are equimolar except in the ARP-alone entries) is not written.

## 3. Tables re-typed

### Table II. "Formation of phenylacetaldehyde from the glucose-phenylalanine ARP, during a Likens-Nickerson steam distillation-extraction procedure."

| entry | starting materials | conditions | pH | yield of phenylacetaldehyde |
|---|---|---|---|---|
| 1 | glucose + phenylalanine | phosphate buffer (100 mM), reflux, 2 h | 7.0 | 0.15 % |
| 2 | ARP from glucose and phenylalanine | phosphate buffer (100 mM), reflux, 2 h | 7.0 | 0.63 % |
| 3 | glucose + phenylalanine | citrate buffer (50 mM), reflux, 2 h | 3.2 | 0.006 % |
| 4 | phenylalanine | citrate buffer (50 mM), reflux, 2 h | 3.2 | 0.006 % |
| 5 | ARP from glucose and phenylalanine | citrate buffer (50 mM), reflux, 2 h | 3.2 | 1.6 % |
| 6 | ARP from glucose and phenylalanine | EtOH, phosphate salts (100 mM), reflux, 2 h | 7.0 | "very small amount, also furfural formed" |

Ratios: ARP/(Glc + Phe) = 4.2 at pH 7.0 and 267 at pH 3.2; ARP pH 3.2 / pH 7.0 = 2.5; (Glc + Phe) pH
7.0 / pH 3.2 = 25; Phe alone = Glc + Phe at pH 3.2 (0.006 %), i.e. at pH 3.2 free glucose contributes
nothing measurable in 2 h.

### Table III. "Formation of phenylacetaldehyde at pH 3.2."

Conditions: ~26 mmol/L carbonyl + 26 mmol/L phenylalanine, 50 mM citrate, reflux 2 h under
Likens-Nickerson. Footnote: "* Glc-Phe ARP was reacted without added phenylalanine."

| carbonyl compound | yield of phenylacetaldehyde |
|---|---:|
| pyruvaldehyde (methylglyoxal) | 21.0 % |
| 3-deoxyerythrosone | 17.2 % |
| dihydroxyacetone | 9.0 % |
| glyoxal | 8.8 % |
| erythrose | 5.7 % |
| glyceraldehyde | 3.7 % |
| 3-deoxyxylosone | 2.6 % |
| hydroxyacetone | 1.2 % |
| Glc-Pro ARP | 0.9 % |
| glycolaldehyde | 0.8 % |
| 3-deoxyglucosone | 0.6 % |
| xylose | 0.6 % |
| Glc-Phe ARP* | 0.5 % |
| fructose | 0.13 % |
| glucose | 0.03 % |

Within-study ratios: pyruvaldehyde/glyoxal = 2.39; pyruvaldehyde/3-deoxyglucosone = 35;
3-deoxyerythrosone/3-deoxyxylosone/3-deoxyglucosone = 17.2 : 2.6 : 0.6 = 28.7 : 4.3 : 1;
dihydroxyacetone/glyceraldehyde = 2.4; fructose/glucose = 4.3; Glc-Pro ARP / Glc-Phe ARP = 1.8 (the
text says "almost twice"); glucose pH 3.2 here (0.03 %) vs Table II entry 3 (0.006 %) differ 5-fold for
what should be the same experiment (see Flags).

**Unit reconciliation and a derived average rate (mine, assumption-laden).** If yield is mol PAC per
mol Phe and the 26 mmol/L is right, pyruvaldehyde gave 0.21 × 26 = 5.5 mmol/L of PAC in 120 min, an
average formation rate of **45 µmol L^-1 min^-1** at ~100 C, pH 3.2; glyoxal 0.088 × 26 / 120 = **19
µmol L^-1 min^-1**; 3-deoxyglucosone 1.3 µmol L^-1 min^-1. Expressed as a second-order constant at the
initial concentrations, k2 = rate / (26 × 26 mmol²/L²) = 6.7e-5 (MGO) and 2.8e-5 (GO) L mmol^-1 min^-1,
assuming the reactants are not depleted (they are, by up to 21 %). These are averages over an open,
product-stripping distillation, not measured rates.

### Table IV. "Formation of ACP and ACTP from Amadori rearrangement products."

Conditions: 50 mM phosphate pH 7.2, 35 mL, 4 h Likens-Nickerson; ARP or parent sugar + proline at the
same molar amounts (150 µmol-2.7 mmol).

| starting compound | yield ACP | yield ACTP |
|---|---:|---:|
| Glc-Pro ARP | 0.016 % | 0.04 % |
| Glc + Pro | 0.017 % | 0.18 % |
| glyceraldehyde-Pro ARP | 0.033 % | 0.24 % |
| glyceraldehyde + Pro | 0.034 % | 0.29 % |

Ratios: ACTP (Glc + Pro)/(Glc-Pro ARP) = 4.5 ("more than 4 x as low" from the ARP); ACP unchanged
(1.06); glyceraldehyde systems give 2 × the ACP of the glucose systems.

### Secondary numbers quoted in the Discussion (pointers, not this chapter's data)

- Fujimaki et al. (1968, ref 26 by content; cited as 27): 3-deoxyglucosone + L-leucine in distilled
  water at 80 C, 30-60 min, gave 560 times less isovaleraldehyde (3-methylbutanal) than pyruvaldehyde.
- Ghiron et al. 1988 (ref 30): 3-deoxyglucosone + phenylalanine in water at 100 C, 30 min: Strecker
  aldehyde formed.
- Chuyen, Kurata & Fujimaki 1972 (ref 19): the decarboxylation rate of alanine + glyoxal increases
  with pH.
- Hofmann & Schieberle (Weurman 1999): oxygen catalyses carboxylic-acid formation in the Strecker
  degradation (cf. `hofmann2000b_extraction.md`).

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`): phenylacetaldehyde → `phenylacetaldehyde`; 2-acetyl-1-pyrroline,
ACTP, methylglyoxal, glyoxal, the deoxyosones, the sugars and the ARPs → not in registry
(reaction_rules.yml uses MGO, GO, 3-DG, Glc).

| quantity | value | unit | conditions | reaction order (authors) | source location | evidence class |
|---|---|---|---|---|---|---|
| PAC yield from Phe + pyruvaldehyde | 21.0 | % (basis undefined; read as mol/mol Phe) | 26 + 26 mmol/L, 50 mM citrate pH 3.2, reflux (~100 C) 2 h, aldehyde stripped continuously | — | Table III | level_only (yield on a fed intermediate) |
| PAC yield from Phe + glyoxal | 8.8 | % | same | — | Table III | level_only |
| PAC yield from Phe + 3-deoxyglucosone / 3-deoxyxylosone / 3-deoxyerythrosone | 0.6 / 2.6 / 17.2 | % | same | — | Table III | level_only |
| PAC yield from Phe + dihydroxyacetone / glyceraldehyde / hydroxyacetone / glycolaldehyde / erythrose / xylose / fructose / glucose | 9.0 / 3.7 / 1.2 / 0.8 / 5.7 / 0.6 / 0.13 / 0.03 | % | same | — | Table III | level_only |
| PAC yield from Glc-Pro ARP (+ Phe) / Glc-Phe ARP (alone) | 0.9 / 0.5 | % | same | — | Table III | level_only |
| k(MGO)/k(GO) as Strecker donors for Phe | 2.39 | — | pH 3.2, ~100 C | — | derived from Table III | within_study_ratio |
| MGO : GO : 3-DG donor ratio | 35 : 14.7 : 1 | — | same | — | derived | within_study_ratio |
| C4 : C5 : C6 3-deoxyosone ratio | 28.7 : 4.3 : 1 | — | same | — | derived | within_study_ratio |
| Glc-Phe ARP vs Glc + Phe, pH 7.0 / pH 3.2 | 4.2 / 267 | — | 100 mM phosphate / 50 mM citrate, 2 h | — | derived from Table II | within_study_ratio |
| Glc + Phe, pH 7.0 vs pH 3.2 | 25 | — | different buffers | — | derived from Table II | within_study_ratio (buffer-confounded) |
| Glc-Phe ARP, pH 3.2 vs pH 7.0 | 2.5 | — | same | — | derived | within_study_ratio (buffer-confounded) |
| PAC yields, Table II | 0.15 / 0.63 / 0.006 / 0.006 / 1.6 | % | as tabulated | — | Table II | level_only |
| ACP and ACTP yields from Pro systems | Table IV (0.016-0.034 % ACP; 0.04-0.29 % ACTP) | % | pH 7.2 phosphate, 4 h | — | Table IV | level_only |
| average PAC formation rate, Phe + MGO / + GO (mine) | 45 / 19 | µmol L^-1 min^-1 | 26 + 26 mmol/L, pH 3.2, ~100 C, 0-120 min average, open distillation | assumed constant rate | derived from Table III | derived_assumption |
| second-order re-expression (mine), MGO / GO | 6.7e-5 / 2.8e-5 | L mmol^-1 min^-1 | same, initial concentrations, no depletion | assumed second order | derived | derived_assumption |
| 3-DG vs MGO with leucine, 80 C (Fujimaki) | 1/560 | — | water, 30-60 min | — | Discussion, secondary | pointer (not admissible as this paper's data) |

Cross-reference inside the repo: `hofmann2000_extraction.md` is the source of rule R07 (the dicarbonyl
Strecker) and `hofmann2000b_extraction.md` the ARP-Phe oxidative route; Zhou 2024
(`zhou2024_extraction.md`) prints alanine + MGO / GO pyrazine rates at pH 8 whose GO/MGO ordering
(pyrazine 8 times faster than 2,5-dimethylpyrazine) is the reverse of this chapter's MGO > GO for the
aldehyde at pH 3.2; the two measure different products (ring vs aldehyde) at different pH. The B18
Strecker constants (`parameters_pyrazine.py`: 10^-6.54 and 10^-7.53 L mmol^-1 min^-1 at 100 C, pH 6.8)
are 2-3 decades below the derived averages here (6.7e-5 and 2.8e-5), which is the expected gap
between an aldehyde yield under product stripping at pH 3.2 and a pyrazine rate that needs two
Strecker events and a condensation; it is recorded, not reconciled.

## 5. Flags

1. **No rate, no barrier, no replicate.** One temperature (reflux), one time (2 h), no SD, no n. Every
   number is a yield; only the within-study ratios travel.
2. **Yield basis undefined.** "Yield of phenylacetaldehyde" in per cent without saying per mole of
   what; equimolar runs make Phe and carbonyl bases coincide, but the ARP-alone entries (Table II 2, 5,
   6; Table III Glc-Phe ARP) are per mole of ARP.
3. **PAC quantification underdocumented.** Undecane as IS is named; no response factor, external
   standard or calibration for PAC is stated (the external standards named are for ACP/ACTP).
4. **Open system.** Likens-Nickerson strips the aldehyde into dichloromethane as it forms, so the
   yields are cumulative formation with the aldehyde protected from further reaction (aldol,
   oxidation to phenylacetic acid): higher than a closed-pot level, and not comparable to Zhou 2024's
   or Zamora 2015's closed tubes without a declared correction.
5. **pH 3.2 vs 3.5.** The Experimental says citrate pH 3.5; the tables and discussion say 3.2. Carry
   as "pH 3.2-3.5".
6. **Internal inconsistency on glucose.** Glucose + Phe at pH 3.2 is 0.006 % in Table II (entry 3) and
   0.03 % in Table III; the same experiment in the two tables differs 5-fold, which is the
   reproducibility band of this chapter's smallest yields.
7. **Buffer confounds pH.** The pH 7.0 runs are 100 mM phosphate, the pH 3.2 runs 50 mM citrate; the
   authors argue the difference in catalysis "is not expected to be significant", which is an
   assertion.
8. **Commercial carbonyls, nominal amounts.** "~1.3 mmol" of pyruvaldehyde and glyoxal (aqueous
   commercial solutions) — hydrate/oligomer content not assayed; the deoxyosones were the authors'
   own preparations (refs 28, 31), purity not stated.
9. **Reference numbering slip** (Fujimaki cited as 27, ref 27 is Kerler); the 560-fold quote is
   secondary and must not be entered as a number.
10. Registry gaps: 2-acetyl-1-pyrroline, ACTP, MGO, GO, the deoxyosones and ARPs have no
    `compounds.yml` rows.
