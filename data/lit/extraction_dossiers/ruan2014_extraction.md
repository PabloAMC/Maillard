# Ruan, Chen, Kong & Hua 2014 — EXTRACTION (laboratory soy protein, DTT-reduced to 7.5-75.9 µmol SH/g, heated 100 C / 30 min)
### Soy protein free-SH, total-cysteine and S-S bookkeeping before and after heating, by 4,4'-dithiodipyridine in 2 % SDS.

**Source on disk:** `data/articles/ruan2014.pdf` (owner's download, 2026-09-08). Read from the pdftotext
layer in the scratchpad (`articles/ruan2014.txt`); the text layer is clean (ligatures and degree signs
come out as `/C176`, `lmol` = µmol). The paper has NO numeric tables: every SH / S-S number lives in the
running text of §3.1.1 and §3.2.3 or in Figures 1a and 3a. Only the text numbers are re-typed below;
the per-DTT-level intermediate points are FIGURE-ONLY. Repo status before this dossier: no soy site
density on file (`data/species/protein_matrices.yml` carries only BLG).

## 0. Identity

| field | value |
|---|---|
| Title | "Heat-induced aggregation and sulphydryl/disulphide reaction products of soy protein with different sulphydryl contents" |
| Authors | Qijun Ruan, Yeming Chen, Xiangzhen Kong, Yufei Hua (State Key Lab of Food Science and Technology, Jiangnan University, Wuxi) |
| Venue | Food Chemistry 156 (2014) 14-22; received 1 Sep 2013, available online 6 Feb 2014 |
| DOI | 10.1016/j.foodchem.2014.01.083 (printed; PII S0308814614001150 per the brief) |
| Naming | "SH" = sulphydryl, "SS" = disulphide bond; "total free SH" = SH measured in 2 % SDS (unfolded, NOT reduced); "total cysteine" = cysteic acid after performic oxidation (= half-cystine). A = acidic, B = basic glycinin polypeptide; a', a, b = beta-conglycinin subunits; KTI = Kunitz trypsin inhibitor. |
| Companion method paper | Ruan, Chen, Kong & Hua 2013, JAFC 61, 2661-2668 (the DPS-vs-DTNB comparison the SH method cites); not on disk |

## 1. Why it matters

The matrix layer (`src/kinetic_core/matrix_sites.py`) charges free thiol and disulfide sites per gram of
protein, and the pre-registration says no dossier states them for soy. This paper gives a laboratory
soy protein isolate's native free SH (7.5 µmol/g protein), its heated (100 C, 30 min) free SH and S-S,
and, by difference, its total half-cystine (~114 µmol/g protein, derived here). It also supplies a
controlled series of the SAME isolate with the SH pool artificially opened by DTT reduction (up to
75.9 µmol/g), which is a direct lever on the question "how much of the thiol pool is available" that
the sulfur lane's thiol-to-protein exchange channel depends on. Its central finding, that opening more
SH does NOT drive S-S polymerisation on heating (S-S content falls with degree of reduction), is a
directional hold-out shape for any model that couples soy aggregation to thiol density.

## 2. Methods as they matter to a model

- **Protein (§2.2):** laboratory isolate, not commercial. Dehulled, milled soybean flour, hexane-defatted
  5x, then hexane/ethanol 1:2 (v/v) 4 C 1 h, then 95 % ethanol 4 C 1 h, dried; dispersed 1:10 (w/v) in
  water, pH 7.0 (2 M NaOH), 1 h 20 C, 15 800 g 30 min; supernatant to pH 4.5, 6000 g; pellet washed twice,
  resuspended 5-fold (w/w) water, neutralised to pH 7.0, 15 000 g; supernatant freeze-dried. "The protein
  content in the prepared soy protein was 90% (w/w), as determined by the micro-Kjeldahl method."
  Nitrogen factor NOT stated. The alcohol wash is an important quirk: it removes soluble sugars but is
  known to partially denature soy protein.
- **Working solution (§2.3):** powder at 30 g/L in 10 mM sodium phosphate pH 7.0, centrifuged 40 000 g
  30 min to discard insolubles. So all SH numbers refer to the SOLUBLE fraction of the isolate.
- **Reduction:** "incubated with freshly prepared DTT between 0.1 and 10 mM under nitrogen overnight at
  25 C" (levels used: 0, 0.1, 0.5, 2.5, 5, 10 mM); excess DTT removed on a HiTrap 5 mL desalting column
  "equilibrated and eluted with distilled water". ⚠ After desalting the protein is in WATER, not the
  phosphate buffer; §3.2.1 nevertheless states the heated solutions were "22.5 mg/ml, pH 7.0". Ionic
  strength during heating is therefore essentially zero (unbuffered).
- **Heating:** screw-capped tubes, water bath 100 C, times 0.5, 1, 2, 3, 30 min, then ice-water. The SH /
  S-S numbers are for 30 min only.
- **Free SH method (§2.4), verbatim:** "The total free sulphydryl (SH) content was determined according
  to the method from a previous study (Ruan, Chen, Kong, & Hua, 2013). Soy protein solutions were
  diluted with the SDS-buffer (pH 7.0) to give a final SDS concentration of 2% (w/v). The SH detecting
  reagent, 4,4'-dithiodipyridine (DPS) was added prior to SDS addition to avoid SH oxidation. The samples
  were vortexed and detected immediately at 324 nm against the SDS-buffer blank in a UV-2450 UV-VIS
  spectrophotometer (Shimadzu, Kyoto, Japan) until the absorbance reached the maximum value (this time
  period was recorded). [...] The SH content was expressed as µmol SH/g protein." Reagent is DPS
  (4,4'-dithiodipyridine, 4-thiopyridone read at 324 nm), NOT Ellman's DTNB; denaturant 2 % SDS, no
  urea/GuHCl, no reducing agent. Extinction coefficient not printed (in the 2013 paper). "Total free
  SH" here means SDS-exposed free SH, i.e. buried + surface, but NOT S-S-derived.
- **Total cysteine (§2.5), verbatim:** "Performic acid was added and incubated at 0 C for 20 h. Both the
  cysteine and cystine could be converted into cysteic acid by performic acid oxidation [...]. Cysteic
  acid was separated from the other amino acids using a Hitachi L-8900 amino acid analyser after
  digestion in 6 M HCl at 110 C for 22-24 h. The total cysteine content was determined in the form of
  cysteic acid. The SS content = (total cysteine residue content - total free SH content)/2."
  ⚠ The total-cysteine value itself is NEVER printed; it is recoverable from the printed SH and SS pairs
  (§4).
- **Protein basis:** "µmol SH/g protein". How protein was quantified in the desalted solutions is not
  stated (the 90 % Kjeldahl figure is for the powder). Treat "per g protein" at face value.
- **Replicates:** "Three separate soy samples were used, and each sample was run in triplicate. [...]
  Data were expressed as the mean ± SD (n = 3)."
- **Units in the repo:** 1 µmol/g = 0.001 mmol/g. No powder-to-protein conversion is needed (values are
  already per g protein); if a spec states soy powder loading, multiply by 0.90 g protein / g powder.

## 3. Tables re-typed

The paper has no numeric tables. Text-stated numbers, by section:

### §3.1.1 Unheated, total free SH (Fig. 1a; text values only)

| DTT (mM) | total free SH (µmol/g protein) | source |
|---:|---:|---|
| 0 (unreduced) | **7.5 ± 0.26** | text |
| 0.1, 0.5 | "slightly increased" | FIGURE-ONLY |
| 2.5, 5 | "greatly increased" | FIGURE-ONLY |
| 10 | **75.9 ± 0.5** | text |

Unheated S-S is not printed for any level (Fig. 1b is a non-reducing SDS-PAGE, not an S-S plot).
Unreduced particle size 23.0 ± 1.6 nm (DLS), unchanged by reduction.

### §3.2.3 Heated 100 C / 30 min, 22.5 mg/mL (Fig. 3a; text values only)

| DTT (mM) | total free SH (µmol/g protein) | S-S (µmol/g protein) | source |
|---:|---:|---:|---|
| 0 (unreduced) | **0.75 ± 0.1** | **56.6 ± 0.56** | text ("increased from 0.75 ± 0.1 to 50.2 ± 0.3"; "decreased from 56.6 ± 0.56 to 31.9 ± 0.5") |
| 0.1 - 5 | monotone between the ends | FIGURE-ONLY | |
| 10 | **50.2 ± 0.3** | **31.9 ± 0.5** | text |

### Other heated-state numbers (text)
- Aggregate size after 30 min: 40 ± 2 nm (unreduced) to 70 ± 2 nm (10 mM DTT).
- Power-law consistency index k: 0.20 ± 0.02 to 0.40 ± 0.01 Pa s^n (n < 1 throughout).
- Ultracentrifugation (270 000 g) supernatant protein: 14 ± 0.7 to 9 ± 0.6 mg/mL (of 22.5).
- Subunit partition: supernatant a'+a+A 53 -> 60 %, b+B 33.7 -> 25.8 %; precipitate a'+a+A 25 -> 12 %,
  b+B 51 -> 60 % (with increasing reduction).
- SS-linked products (Fig. 6, text): polymer + (dimer of B + monomer of A) 45 -> 20 %; polymer alone
  40 -> 8 %; correlation with heated S-S content R² = 0.9. The b subunit "had no SH" and does not take
  part.
- Non-reducing SDS-PAGE: glycinin AB and A5B3 bands vanish after 3 min (0-0.5 mM DTT), 2 min (2.5-5 mM),
  within 1 min (10 mM).

## 4. Site densities the repository can use

Protein basis: all printed values are per g PROTEIN (paper's unit); no conversion assumption needed
beyond taking the paper at its word. Derived rows use the paper's own formula
SS = (total Cys - free SH)/2, inverted.

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN | conditions | source | evidence |
|---|---|---:|---|---:|---|---|---|
| soy protein isolate, lab-made, alcohol-washed, soluble fraction | free SH (SDS-exposed, unreduced) | 7.5 ± 0.26 | µmol SH/g protein | **0.0075** | native, 30 g/L, 10 mM phosphate pH 7.0, 25 C | §3.1.1 | measured |
| same | total half-cystine (cysteic acid) | ~114 (113.95 from 0 mM pair; 114.0 from 10 mM pair) | not printed | **~0.114** | independent of heating | back-calculated: 2 x 56.6 + 0.75; 2 x 31.9 + 50.2 | inferred (arithmetic on printed numbers; sd not propagable) |
| same | S-S, native | (114 - 7.5)/2 = **~53** | not printed | **~0.053** | native | derived from the two rows above | inferred |
| same | free SH after heating | 0.75 ± 0.1 | µmol SH/g protein | **0.00075** | 100 C, 30 min, 22.5 mg/mL, water (nominal pH 7.0), unreduced | §3.2.3 | measured |
| same | S-S after heating | 56.6 ± 0.56 | µmol SS/g protein | **0.0566** | same | §3.2.3 | measured (via the SS formula, using an unprinted cysteic-acid total) |
| same, 10 mM DTT-reduced then desalted | free SH, unheated | 75.9 ± 0.5 | µmol SH/g protein | 0.0759 | native after reduction | §3.1.1 | measured |
| same, 10 mM DTT | S-S, unheated | (114 - 75.9)/2 = ~19 | not printed | ~0.019 | | derived | inferred |
| same, 10 mM DTT | free SH / S-S after heating | 50.2 ± 0.3 / 31.9 ± 0.5 | µmol/g protein | 0.0502 / 0.0319 | 100 C, 30 min | §3.2.3 | measured |
| same, 0.1-5 mM DTT | all of the above | — | — | — | — | Figs 1a, 3a | FIGURE-ONLY |

Fraction of half-cystine present as free SH in the native isolate: 7.5 / 114 = **6.6 %** (cf. Shimada &
Cheftel 1988: 8 %). Heating the unreduced isolate at 100 C for 30 min consumed **90 %** of the free SH
(7.5 -> 0.75), i.e. 6.75 µmol/g SH -> 3.4 µmol/g new S-S (from ~53 to 56.6).

## 5. Flags

1. **Every S-S number rests on an unprinted total-cysteine value.** The paper defines SS from cysteic
   acid but never prints the cysteic-acid content. The two printed (SH, SS) pairs invert to 113.95 and
   114.0 µmol/g, so a single value ~114 µmol/g protein was used; that agreement is the only evidence
   for it. The repo should carry ~114 as *inferred*, with Shimada & Cheftel's 104.0 ± 6.6 (amino-acid
   analysis, commercial SPI) as the independent cross-check.
2. **0.75 vs 7.5 after heating.** The heated unreduced free SH (0.75 ± 0.1) is exactly one decade below
   the native value (7.5 ± 0.26), which invites suspicion of a misplaced decimal. It is, however,
   arithmetically consistent with the printed SS of 56.6 under the ~114 total (with 7.5 the SS would be
   53.25), so it is taken as printed. A 90 % loss is much larger than Shimada & Cheftel's 40 % maximum
   at 80 C; the higher temperature (100 C) and the unbuffered water medium may explain it. Treat the
   heated free-SH value as measured-but-unconfirmed.
3. **Method is DPS in 2 % SDS, not Ellman's DTNB.** Ruan et al. 2013 (not on disk) found reagent-
   dependent SH values for soy protein; direct comparison with DTNB numbers (Shimada; Chihi) carries a
   method offset of unknown sign.
4. **Soluble fraction only.** The 40 000 g pre-centrifugation discards insoluble protein; the site
   densities are for what dissolves at 30 g/L, pH 7.0, I ~ 0.01, of an alcohol-washed isolate, not for
   the whole powder.
5. **Medium during heating is distilled water** (desalting eluent), despite the "pH 7.0" statement;
   ionic strength during the 100 C step is unstated and near zero.
6. **Protein basis of the solution assays is not documented** (how "g protein" in "µmol SH/g protein"
   was measured for the desalted solutions). The powder is 90 % protein by micro-Kjeldahl, factor
   unstated (soy convention 6.25 or 5.71).
7. **Intermediate DTT levels are FIGURE-ONLY** (Figs 1a, 3a); the dossier records only the two ends.
8. **Directional hold-out:** more free SH -> larger aggregates, higher viscosity, but LESS S-S and fewer
   SS-linked polymers after heating (R² = 0.9 between S-S content and polymer yield). Any coupling in
   the repo that makes soy aggregation scale with thiol density runs against this paper.
