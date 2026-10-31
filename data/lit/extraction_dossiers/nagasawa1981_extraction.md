# Nagasawa, Goon & Shirota 1981 — EXTRACTION (C-2 epimerisation of 2,5,5-trimethylthiazolidine-4-carboxylic acid from PENICILLAMINE + acetaldehyde, followed by mutarotation and ¹H NMR at room temperature: ~24 h to equilibrium in neutral D₂O / methanol, < 1 h in 1 N HCl, no change in dilute Na₂CO₃; acid catalysis inferred; NO rate constant, NO equilibrium constant, NO cysteine or sugar thiazolidine, NO buffer or temperature series)

### A STEREOCHEMISTRY PAPER WITH ONE KINETIC OBSERVATION. The ring-chain step (C-S cleavage to an iminium and re-closure) is shown to run on an hours scale at room temperature in neutral protic solvent, to be acid-catalysed (minutes in 1 N HCl) and to stop in weak base, WITHOUT release of free amino thiol (TLC). It gives a lower bound on a ring-opening rate for an analogue — a gem-dimethyl penicillamine ring with a methyl at C-2 — not for TTCA, and nothing on the equilibrium with free cysteine + aldehyde.

**Source on disk:** `data/articles/nagasawa1981.pdf` (5 pp., J. Heterocycl. Chem. pp. 1047-1051; Acrobat Paper Capture
OCR of a scan). Prose read from the `pdftotext -layout` layer; **Figures 1-2 and Schemes 1-4 read from 110 dpi page
renders** (scratchpad `pdf2/nag-2.png`, `nag-3.png`), the kinetic passage on p. 1047 checked against a 200 dpi crop
(`pdf2/nag1crop-1.png`). The OCR layer garbles stereodescriptors ("2(R,5')" = 2(R,S)) and loses Fig. 1's rotations.

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[F]** read from a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed / hypothesised, not
measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question put to it (plan step B7b, `results/validation/kinetic_core_dpo_partition_prereg.md` secs. 7 and 12,
blocking item B2): the engine's TTCA ring opening `k_ttca_cys` (~0.028 min⁻¹, fitted, no pH factor; `r_ttca_form`
by detailed balance with Cao's K ~60 L/mol) cannot reconcile Cao 2017 (`cao2017_extraction.md`: 0.4 M, 0.2 M
phosphate pH 5.5, equilibrium within 30 min at 90 °C) with Zhai 2021 (`zhai2021_extraction.md`: 10 mM TTCA,
unbuffered water, natural pH 5.5, loss ≤ 0.0033 min⁻¹ at 100 °C).

| what the model needs | what Nagasawa 1981 gives | tag |
|---|---|---|
| rate of ring opening | **bounds only, for an analogue.** 1a (2S-rich) in water / methanol "reaching equilibrium essentially within 24 hours" (rotation +108° → +152° in D₂O by 24 h, unchanged at 48 h); in 1 N HCl "Equilibration occurred rapidly (within 1 hour)"; in D₂O / 2 % Na₂CO₃ rotations and NMR constant over 24 h. Only 0 / 24 / 48 h are shown | [M] / [F] / [Q] |
| equilibrium constant with free thiol + aldehyde | **none.** "tlc does not show the presence of the fully dissociated penicillamine" during epimerisation | [Q] |
| pH dependence | **three points, qualitative:** fast in 1 N HCl, ~24 h in neutral protic solvent, none in weak base → "the epimerization is acid catalyzed" | [Q] / [P] |
| temperature dependence | **none** (room temperature; ice-bath synthesis) | — |
| buffer / general-acid catalysis | **not tested** (HCl, water, methanol, Na₂CO₃, pyridine only) | — |
| mechanism | Scheme 4: zwitterion 1d (2R,4S) ⇌ transient open iminium-thiolate 1e ⇌ 1f (2S,4S); no ²H exchange at the C-2 methyl in D₂O, so no enamine; authors cite strongly basic solutions as favouring ring-to-chain reversal (refs 10a, 11) | [P] / [Q] / [C] |
| epimer ratios | 1a cis:trans 3.3:1; equilibrium product 1b 1.8:1; 1c (from ethyl acetate) 69:31 % (pyridine-d₅ NMR) | [M] |

**Answer to B2: no rate constant and no K.** Confidence ~99 % (all five pages read). What it gives, weighed:

1. **A lower bound on the ring-opening rate of an analogue** [Z]. "Equilibrium essentially within 24 h" ≥ ~3
   relaxation times: k_relax ≳ 3 / 1440 min ≈ **2 × 10⁻³ min⁻¹** in neutral D₂O at room temperature (~22 °C
   assumed; not printed), 1a at c ~1.0 g/100 mL ≈ 57 mM [Z: M 175.2]. Epimerisation relaxation = k_open × (sum of the
   closure branching), so k_open is of the same order. In 1 N HCl: ≳ 3 / 60 ≈ **0.05 min⁻¹**. Both are bounds (no
   intermediate time point), not measurements; D₂O may carry a solvent isotope effect.
2. **Transfer to the engine's conditions** [Z, assumption-laden]: with an assumed Ea 50-90 kJ/mol (none printed;
   nothing in the corpus pins it for a thiazolidine ring opening), 22 → 100 °C multiplies by ~70-2,000, giving
   k_open(100 °C, neutral) ≳ 0.15-4 min⁻¹ — 50-1,000× above Zhai 2021's 0.0033 min⁻¹ cap and consistent with Cao's
   ≥ 0.013 min⁻¹. **Structural caveats** cut both ways: the 5,5-gem-dimethyl ring (penicillamine) should favour
   closure and slow opening relative to cysteine (Thorpe-Ingold); a polyhydroxyalkyl at C-2 (TTCA) withdraws
   electrons and should destabilise the iminium and slow opening relative to a C-2 methyl. Net direction unknown;
   order of magnitude uncertain by ≥ 10×.
3. **What it does not bound.** The step observed is C-S ring opening + re-closure (no free penicillamine on TLC). The
   engine's `r_ttca_cys` is full dissociation (TTCA → Cys + PENT), which also needs hydrolysis of the iminium. That
   hydrolysis — and any general-acid (phosphate) catalysis of it — is untouched here.
4. **pH sign.** Acid accelerates the ring-chain step (specific acid; general acid untested). Cao's pot (pH 4.7-4.15
   after heating) is more acidic than Zhai's (initial pH 5.5): a specific-acid law could give ≤ ~3-6× [Z: 0.5-0.8 pH
   units, first order in H⁺], not the ≥ 4-30× the B2 contradiction needs if it were a rate problem.

**Bearing on the Cao-vs-Zhai contradiction (my reading, confidence ~0.65):** together with Schubert 1939 (sugar
thiazolidines yield their thiol to iodine as fast as free cysteine, `schubert1939_extraction.md`) and Roberts 1987
(minutes-scale mutarotation, `roberts1987_extraction.md`), this makes "the ring opening is intrinsically slow in
water" unlikely (P ~0.2 that a slow C-S opening is the cause). The contradiction then sits in (a) the **equilibrium
position** (Zhai's 10 mM pot would need K ≳ 10⁴ L/mol at pH 5.5 / 100 °C to stay ≥ 90 % closed if equilibrated,
against 60 L/mol in Cao, where K = 60 would leave only ~23 % closed [Z]; Kallen 1971's pH law, ~20× over 1.9 pH
units, gives only ~3-4× over Cao-to-Zhai's 0.8 pH units, so pH alone is unlikely to close it), (b) a
**slow, possibly buffer-catalysed hydrolysis** of the open iminium (rival R-P relocated from the ring opening to the
hydrolysis), or (c) Zhai's analytics re-forming TTCA after quench. This paper cannot rank them.

**Can anything anchor the engine's constant?** No: a bound on an analogue's C-S opening at one temperature, with no
K. It could at most back a declared lower bound on a split-out ring-opening step if the engine ever separates
opening from hydrolysis. **Keying: not worth it.** Fig. 1's three time points per series are rotations, not
concentrations, and the endpoints are printed.

**Leads (none on disk; DOIs checked against Crossref 2026-10-31):**
- **Pesek & Frost 1975**, *Tetrahedron* 31:907-913, "Decomposition of thiazolidines in acidic and basic solution" —
  DOI 10.1016/0040-4020(75)80099-8 (ref 10a). **The most direct lead for ring-opening / hydrolysis rates and their pH
  profile.**
- Kallen 1971b, *J. Am. Chem. Soc.* 93:6236-6248 — DOI 10.1021/ja00752a040 (ref 10b; already in
  `docs/guides/LITERATURE_SEARCH_2026-10-29_phosphate_thiazolidine.md`).
- Schubert 1936, *J. Biol. Chem.* 114:341-350 — DOI 10.1016/s0021-9258(20)65226-7 (ref 1a).
- Luhowy, Cieciuch & Meneghini 1973, *Tetrahedron Lett.* p. 1285; Stacey & Strong 1967, *J. Org. Chem.* 32:1487
  (ref 11: base-favoured ring-chain reversal). DOIs not verified.
- Carroll, Jung & Sklavounos 1977, *J. Heterocycl. Chem.* 14:503 (ref 5: penicic acid C-2 epimerisation at elevated pH
  or temperature). DOI not verified.

## 0. Identity

| field | value |
|---|---|
| Title | "Epimerization at C-2 of 2-Substituted Thiazolidine-4-carboxylic Acids" |
| Authors | Herbert T. Nagasawa, David J. W. Goon, Frances N. Shirota — VA Medical Center and Dept. of Medicinal Chemistry, University of Minnesota, Minneapolis |
| Venue | *J. Heterocycl. Chem.* 1981, **18**(5), 1047-1051 (Aug. 1981); received October 17, 1980. Printed: "J. Heterocyclic Chem., 18, 1047 (1981)"; code 0022-152X/81/051047-05$02.25 |
| DOI | 10.1002/jhet.5570180538 (not printed; Crossref returns this title, these three authors, J. Heterocycl. Chem. 18(5) 1047-1051, 1981-08; checked 2026-10-31). **Matches the brief.** |
| Funding | Veterans Administration |
| Paper type | Physical-organic stereochemistry: polarimetry, ¹H NMR (60 MHz, Varian T-60), chemical correlation via münchnone cycloadditions |
| Same lab | Roberts 1987 (`roberts1987_extraction.md`); Nagasawa 1978, J. Med. Chem. 21:1274 (ref 3b); 1984, J. Med. Chem. 27:591 — not on disk |

## 2. Methods as they matter to a model

- **Compounds** [M]: 2(R,S),5,5-trimethylthiazolidine-4(S)-carboxylic acid (**1**: 1a precipitate, 1b mother-liquor
  equilibrium product, 1c ethyl-acetate recrystallisate) from **D-penicillamine + acetaldehyde**; 2(R,S)-phenyl-5,5-
  dimethyl analogue (**5**) from benzaldehyde; N-acetyl derivatives 2a / 2b / 6a-c. **No cysteine, no sugar.**
- **Synthesis of 1** [M]: D-penicillamine 0.081 mol + acetaldehyde 0.085 mol, aqueous, ice bath; 1a precipitates in
  15-20 min (32 %), total 30 min; mother liquor 2 h at room temperature gives 1b (65 %). Left at room temperature, the
  precipitate "redissolved after 1 hour" and work-up gave the 1c-like mixture [Q].
- **Kinetic observation** [M/F]: [α]D in water or methanol (c ~1.0) and ¹H NMR of the C-4 methine (δ 4.12 / 4.02) at
  0, 24, 48 h (D₂O, CD₃OD) and 0, 24, 72 h (pyridine-d₅). Temperature: room (not printed). No pH measured.
- **Analytics** [M]: Perkin-Elmer 141 polarimeter; TLC on silica GF-254 (n-butanol : acetic acid : water); chiral
  shift reagent Eu(hfbc)₃ for the pyrrolothiazoles.

## 3. Numbers printed (none is a kinetic or thermodynamic constant)

Fig. 1 / Fig. 2 rotations [F] (read from the render; time in h; rotations measured in the protic solvent named):

| series | 0 h | 24 h | 48 h | 72 h |
|---|---:|---:|---:|---:|
| 1a, D₂O | +108° | +152° | +152° | — |
| 1b, D₂O | +147° | +149° | +149° | — |
| 1a, CH₃OD | +103° | +132° | +132° | — |
| 1a, pyridine-d₅ (NMR) | +80.7° | +89.7° | — | +93.7° |
| 1b, pyridine-d₅ (NMR) | +130° | +108° | — | +93.5° |

Text [M]: in D₂O / 2 % Na₂CO₃, constant 1a +146°, 1b +133° (over 24 h). 1c in methanol +120.9° → +137.1° (24 h) →
same at 47 h; recrystallised 1a +121.7° → +137.8° (24 h) → +138.5° (47 h). Epimer ratios: 1a 3.3:1 cis, 1b 1.8:1,
1c 69:31 %; 5 (phenyl) 3:1 cis as isolated, 5:1 after acetonitrile.

## 4. Figures

Figs. 1-2 (NMR stacks with rotations at three times): the rotations are tabulated above; no digitising beyond
reading the printed labels. **No two-reader keying needed.**

## 5. What the repository takes

Nothing quantitative: no FIT row, no constant, no pot. Context only, never a prior: (i) the C-2 ring-chain step of a
2-alkyl thiazolidine-4-carboxylic acid equilibrates within ~24 h at room temperature in neutral protic solvent and
within 1 h in 1 N HCl, and is frozen in weak base [Q/M]; (ii) it proceeds without detectable release of the free
amino thiol [Q] — ring opening ≠ dissociation.

## 6. Flags

1. **Penicillamine, not cysteine; acetaldehyde / benzaldehyde, not an aldose.** Every statement is an analogue.
2. **Only bounds.** Three time points per series; the 24 h and 48 h values coincide, so the relaxation time is
   anywhere below ~8 h. Do not fit a rate to Fig. 1.
3. **Room temperature and pH not printed.** "Neutral" means the zwitterion dissolved in water or methanol.
4. **Epimerisation ≠ the engine's `r_ttca_cys`** (full dissociation to Cys + PENT).
5. **OCR** garbles "2(R,S)" as "2(R,5')" and drops the Fig. 1 rotations; use the renders.
6. **Registry:** penicillamine and the 2,5,5-trimethyl thiazolidine are not compound keys (not needed).
