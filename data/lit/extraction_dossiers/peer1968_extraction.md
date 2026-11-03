# Peer & van den Ouweland 1968 — EXTRACTION (file stem `peer1968`)
DOI: 10.1002/recl.19680870906 (**not printed** on the scan; found by Crossref bibliographic query 2026-11-03):
"Synthesis of 4-hydroxy-5-methyl-2,3-dihydrofuran-3-one from D-ribose 5-phosphate", authors Peer, van den Ouweland,
*Recueil des Travaux Chimiques des Pays-Bas* **87**:1017-1020, 1968. **Matches the PDF and the brief.**
Not to be confused with the companion paper immediately before it, Peer, van den Ouweland & de Groot 1968, Recl. Trav.
Chim. 87:1011-1016 (10.1002/recl.19680870905; aldopentoses + secondary amine salts → the same furanone; its ref. 1,
**not on disk**). The compound "I" (4-hydroxy-5-methyl-2,3-dihydrofuran-3-one) is norfuraneol (NF,
4-hydroxy-5-methyl-3(2H)-furanone); the text layer OCRs "5-methyl" as "S-methyl".

### D-ribose 5-phosphate (R5P, 1.77 mmol; 11-30 mM) heated 2 h at 100 °C in four concentrated buffers at pH 5.5, no amine: NF yields 43-47 %; preparative run (0.18 mol R5P in sodium acetate / acetic acid, 2 h reflux) 43 % (UV) / 45 % (GC) / 33.8 % isolated. ¹⁴C tracer: ribose C1 does not end in the NF methyl. **One time point, one temperature; no rate, no time course, no pH series in numbers, no H2S.**

**Read 2026-11-03 from `data/articles/peer1968.pdf`** (4 pp., journal pp. 1017-1020, Acrobat Paper Capture OCR of a
scan). Prose from `pdftotext -layout`; **the Table and the Experimental (p. 1019) read from a 250 dpi render**
(`ingest/peer1968_work/pg-3.png`; the text layer drops the Table's yield column), the mechanism scheme (p. 1018) from
a 120 dpi render (`pg-2.png`). No figures with data.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** **Unilever Research Laboratory Duiven**, The Netherlands (H.G. Peer, G.A.M. van den Ouweland; no
corresponding author marked). **Not in `data/keys/laboratories.yml`.** Same proposal as `vandenouweland1975_extraction.md`:
**NEW `unilever_duiven_van_den_ouweland`**, strings "Unilever (van den Ouweland)" (already used in sulfur.yml for
van den Ouweland & Peer 1975 and US 3,904,655) plus a new "Unilever Duiven (Peer / van den Ouweland)". The ledger's
"Unilever Duiven (Tonsbeek)" (`tonsbeek1969_extraction.md`) is the same site and the same NF-in-meat-broth programme
(this paper's "recent work in this laboratory" on NF in meat broth is Tonsbeek 1969; the acknowledgements thank
Losekoot, a Tonsbeek co-author): merging it into the same id for leave-one-lab-out is the conservative choice
(~75 % that they are one group in the fitting sense). Same-lab dossiers: `vandenouweland1975_extraction.md`,
`us3904655_extraction.md`, `tonsbeek1969_extraction.md`.

## 1. What this paper answers for the model

| model need | what Peer & van den Ouweland 1968 give | tag |
|---|---|---|
| NF yield from R5P vs buffer | Table: 43 / 43 / 45 / 47 % (acetate / citrate / phthalate / phosphate), pH 5.5, 2 h, 100 °C | [M] |
| NF yield vs time | **none** (2 h only) | [!] |
| NF yield vs pH | "highly dependent on the pH", best at 5.5 [Q]; **no numbers** | [Q] |
| NF yield vs temperature | **none** (100 °C / reflux only) | [!] |
| rate constant | **none**; lower bound on R5P consumption k ≥ 0.005 min⁻¹ [Z] | [Z] |
| amine catalysis | not required: buffer anions alone give 43-47 %; pyridine ~7 % [M]; many bases and anions "active" [Q] | [M] / [Q] |
| label position | [1-¹⁴C]R5P → NF keeps 98 % of the specific activity; Kuhn-Roth acetic acid (methyl + C5) carries 4.5 % | [M] |
| mechanism | base-catalysed β-elimination of phosphate (C5), then acid-catalysed dehydration of the 5-deoxy-pentosulose | [P] supported by tracer |
| NF UV constants | λmax (MeOH) 289 nm, ε = 9,200; m.p. 126-127 °C | [M] |
| NF + H2S → MFT | **nothing** (no sulfur in this paper) | — |

**Answers.**
1. **Conditions** (confidence ~95 %). Table: 1.77 mmol R5P in 60-155 mL of very concentrated buffer, pH 5.5, 2 h at
   100 °C; [R5P] 11-30 mM, buffer 0.6-5 M, buffer : R5P ~100-280 : 1 [Z] (sec. 2). Preparative: 0.18 mol R5P (made by
   acid-resin hydrolysis of IMP/GMP) in ~5.3 L containing 20.3 mol NaOAc + 4.32 mol AcOH, refluxed 2 h. No amine in the
   Table runs. pH: "~5.5 at the end of the reaction" is the stated criterion for the nucleophile screen [Q]; when and
   at what temperature the Table's "pH 5.5" was read is not stated.
2. **Yield vs time / conditions: one time point** (confidence ~99 %). 43-47 % across four chemically different buffers
   at the same pH and time, so the buffer identity barely matters (spread within what a single UV/GC determination
   resolves). Preparative scale reproduces the acetate row (43 % UV, 45 % GC). Pyridine "approximately 7 %" with no
   time, concentration or pH value. No yield-vs-time curve exists, so it is unknown whether 45 % is a plateau (R5P
   exhausted, NF selectivity ~45 %) or a still-rising point.
3. **Rate: none printed.** Only a bound [Z]: if R5P is lost first-order and at most all of it becomes NF, yield Y at
   t = 120 min requires k_loss ≥ −ln(1 − 0.47)/120 = **0.0053 min⁻¹** (t½ ≤ 131 min) at ~100 °C in 0.6-5 M buffer, pH
   5.5. NF's own loss under these conditions is unmeasured, which would raise the bound. The repo's Larralde 1995
   reading (R5P decomposes ~12 × faster than ribose at pD 7.4, 0.05 M phosphate, 100 °C [C]) is consistent with it.
   Confidence ~90 % that no rate can be extracted beyond this bound.
4. **Relevance to E6 (MFT via NF + H2S): source term and labelling, not the H2S step** (confidence ~85 %). (a) It
   establishes an amine-free, phosphate-leaving-group route R5P → NF at ~45 % molar yield in 2 h at 100 °C in
   concentrated buffer, i.e. R5P is a far richer NF source than free ribose (which gives NF in the companion paper
   only with secondary amine salts [C]). In beef, where R5P is a measured pool (aroma_other.yml composition entries),
   this route can feed the NF → MFT step; in plant-based recipes R5P is absent unless yeast extract nucleotides are
   cleaved (cf. `kavitha2007_extraction.md`), so P(the route matters for a plant-based recipe) ~15 %. (b) **Labelling**:
   the R5P route puts ribose C1 in the ring (C2 of NF, "we may assume" [Q]) and the methyl from C5, the opposite of the
   amine (Amadori / 1-deoxyosone) route, where the companion paper found the methyl from C1 [Q, C]. Any E6 / E11
   labelling cell with [1-¹³C]ribose would read the two routes differently. The corpus's labelled pots use free
   ribose, not R5P, so P(this changes an E6 cell today) ~5 %. (c) It gives no H2S, no MFT, no barrier.
5. **Ribose-5-phosphate as a pentose source** (confidence ~80 %): a fast, amine-independent NF route at pH 5.5 that
   the engine (no R5P species; `src/cultivated_tissue_invariance.py` marks R5P "UNMAPPED PRECURSOR") cannot represent.
   If R5P is ever added, this paper sets the NF selectivity ceiling (~0.45-0.47 per R5P at 2 h, 100 °C, pH 5.5) and
   the label topology, not the rate.

## 2. Conditions

- **R5P source** [M]: Khym 1954 method [C]: 150 g (0.3 mol) commercial 50/50 disodium 5′-inosinate / 5′-guanylate
  (20 % H2O) added to 900 mL Bio-Rad AG 50W-X8 resin in 450 mL water at 100 °C, 4 min, ice-cooled, filtered:
  1000 mL containing 0.18 mol R5P (60 % of the nucleotide [Z]). Purity, counter-ion and free phosphate not stated.
- **Table runs** [M] (1.77 mmol R5P each, pH 5.5, 2 h, 100 °C); concentrations [Z]:

| buffer | total volume | printed amounts | [R5P] | buffer species [Z] | nominal pH [Z, infinite dilution] |
|---|---|---|---|---|---|
| NaOAc / AcOH | 87 mL | 254 mmol NaOAc, 54 mmol AcOH | 20.3 mM | 2.92 M OAc⁻, 0.62 M AcOH | 5.43 |
| NaOH / citric acid | 60 mL | 300 mmol NaOH, 117 mmol CA | 29.5 mM | 1.95 M citrate, 5.0 M Na⁺ | ~6.4 (high ionic strength lowers it) |
| NaOH / KH-phthalate | 155 mL | 175 mmol NaOH, 330 mmol KH-pht | 11.4 mM | 2.13 M phthalate, 1.13 M NaOH | 5.46 |
| Na2HPO4 / KH2PO4 | 100 mL | 103 mmol Na2HPO4, 382 mmol KH2PO4 | 17.7 mM | 4.85 M phosphate | 6.63 (pKa2 7.2; ~5.5-6 at I ≈ 5-7 M) |

- **Preparative run** [M]: R5P solution 1000 mL + 1646 g (20.3 mol) NaOAc + 259 g (4.32 mol) AcOH + 4300 mL water,
  refluxed 2 h; [Z] ~5.5-6.5 L, [R5P] ~28-33 mM, [NaOAc] ~3.1-3.7 M; NaOAc / AcOH = 4.70, the Table's acetate ratio
  (254 / 54 = 4.70) scaled ~80 ×. Reflux of ~3.5 M NaOAc boils a few K above 100 °C (ebullioscopic estimate ~103 °C
  [Z]).
- **Vessel / headspace / O2**: open reflux in air (preparative); Table vessel not stated. **Replicates**: none stated.

## 3. Analytical method

- NF in the preparative residue [M]: UV spectrophotometry (85 % of the 10.3 g residue) and GC (90 %), giving 43 % and
  45 %; crystallised 6.93 g = 33.8 %. **The method behind the Table yields is not stated** [!] (the acetate row equals
  the preparative UV figure). No calibration curve, no error, no replicate.
- Tracer [M]: Packard Tri-Carb 3375 liquid scintillation, external-standard efficiency; Kuhn-Roth oxidation [C].
- NF identity: m.p. 126-127 °C, λmax (MeOH) 289 nm, ε 9,200 [M]. Useful for UV quantification of NF elsewhere
  (e.g. Tonsbeek 1969's UV screen).

## 4. Data

### 4a. Table (p. 1019) [M]: NF yield from 1.77 mmol R5P, pH 5.5, 2 h, 100 °C

| buffer | yield (mol % on R5P) |
|---|---|
| NaOAc / AcOH | 43 % |
| NaOH / citric acid | 43 % |
| NaOH / KH-phthalate | 45 % |
| Na2HPO4 / KH2PO4 | 47 % |

### 4b. Preparative yield (p. 1019) [M], checked [Z]

10.3 g residue × 0.85 = 8.76 g = 76.7 mmol NF (114.10 g/mol) = **42.6 %** of 180 mmol (printed 43 %); × 0.90 = 81.2 mmol
= **45.1 %** (printed 45 %); 6.93 g = 60.7 mmol = **33.7 %** (printed 33.8 %). All reproduce.

### 4c. Other text numbers [M]

- Pyridine, "weakly acidic aqueous medium at 100°", NF ~7 % (p. 1017); time and concentrations not given.
- Nucleophiles screened, all "active" when the end pH is ~5.5 [Q]: diisopropylamine, dicyclohexylamine, pyridine,
  collidine; acetate, citrate, phthalate; carbonate (Pb), hydroxide (Na, Ca, Ba), azide, sulfite, thiosulfate (Na).
  Yields for these are not printed.
- Tracer: [1-¹⁴C]ribose 1.51 × 10⁶ dpm/mmol → [1-¹⁴C]R5P → NF 1.48 × 10⁶ dpm/mmol (98 % retained [Z]); Kuhn-Roth
  acetic acid 6.88 × 10⁴ dpm/mmol = "4.5 %" of NF (my ratio 4.65 %; 4.56 % against the ribose [Z]).

## 5. Flags and internal consistency

1. [!] Single time point and temperature; yields are net of any NF loss; no plateau test.
2. [!] Table yield method not stated; no replicates or errors anywhere.
3. [!] "pH 5.5": the time / temperature of the reading is not stated; the four buffers are 0.6-5 M at extreme ionic
   strength, where nominal Henderson-Hasselbalch pH differs by up to ~1 unit (phosphate, citrate rows) [Z]. KH2PO4 at
   3.8 M and KH-phthalate at 2.1 M exceed their room-temperature solubility (~1.6 M and ~0.4 M [Z, handbook values
   from memory, ±20 %]): the mixtures were presumably prepared or kept hot.
4. [!] The C2 assignment of the label is inferred ("we may assume" [V]); only "not in the methyl / C5" is measured.
5. [!] No DOI printed; OCR text layer garbles "5-methyl", "¹⁴C", "D-".

## 6. Existing records in the repo that cite this paper

grep 2026-11-03 (`peer`, `ouweland`, `recl.1968*`, `Recueil`, `Rec. Trav`, `ribose-5-phosphate` across `data/ docs/
src/ scripts/ tasks/`). Many records cite van den Ouweland & Peer **1975** (JAFC 23:501) and US 3,904,655; those are
out of scope. **Only `data/lit/extraction_dossiers/tonsbeek1969_extraction.md` cites the 1968 Recueil paper** (and the
`pdf_dossier_aliases.yml` pending line, not checked as a record). No ledger, registry, prior, code or plan record cites
it; `data/keys/papers.yml` has no entry for it.

### 6a. `tonsbeek1969_extraction.md` (l. 73, 87, 95, 105)

| claim | paper | verdict |
|---|---|---|
| Peer 1968 = "Rec. Trav. Chim. 87:1017" | 87 (1968) 1017-1020 | **matches** |
| NF "in good yield from ribose-5-phosphate in concentrated buffer systems at pH 5.5" (Tonsbeek's words) | 43-47 % in 0.6-5 M buffers, pH 5.5 | **matches** |
| mechanism "base-catalyzed β-elimination of phosphoric acid", then acid-catalysed dehydration | p. 1018, proposed and tracer-supported | **matches** |
| "that paper, not this one, carries yields" | yes, at one time point; no rates | **matches** (add: no rates) |
| Peer "report the opposite in concentrated buffer, i.e. buffer catalysis is the operative variable" | Peer differs from Tonsbeek in buffer (0.6-5 M vs none) **and** [R5P] (11-30 mM vs ~0.4 mM, ~50 ×) and in detection (isolation vs UV screen); no dilute-buffer arm | **partly**: "buffer catalysis is the operative variable" is an inference [P], not shown by either paper |
| "norfuraneol from R5P needs a nitrogen partner" (point ii, from Tonsbeek's pots) | acetate, citrate, phthalate, phosphate (no nitrogen) give 43-47 % | **contradicted by Peer** under concentrated buffer; true only of Tonsbeek's dilute unbuffered pots |
| "Engine relevance ... a ribose-5-phosphate → 1-deoxyosone shortcut" | tracer: methyl not from C1, β-elimination at C5 → a **5-deoxy** route, unlike the amine / 1-deoxyosone route | **mismatch**: the R5P route is not a 1-deoxyosone shortcut |
| "Action: if an R5P route is ever wanted, the quantitative source is Peer & van den Ouweland 1968" | yields yes; no rate, no time course | **matches**, with the caveat that it bounds selectivity and topology, not kinetics |

### 6b. `data/lit/measurement_ledger/trunk.yml` → `tonsbeek1969_text_r5p_model`

Does not cite Peer; its note ("R5P alone gives no detectable NF") is Tonsbeek's dilute-pot result and is not
contradicted as a statement about that pot. A cross-reference to `peer1968_*` (sec. 7) would prevent its being read as
"R5P → NF needs an amine".

## 7. Proposed ledger entries (trunk.yml; not written)

```yaml
  - id: peer1968_table_norfuraneol_yield_r5p_four_buffers
    dossier: peer1968_extraction.md
    paper: "Peer & van den Ouweland 1968 Recl Trav Chim Pays-Bas 87:1017"
    lab: "Unilever Duiven (Peer / van den Ouweland)"
    system: "D-ribose 5-phosphate 1.77 mmol (11-30 mM), no amine"
    medium: "NaOAc/AcOH 3.5 M, NaOH/citric acid 1.95 M citrate, NaOH/KH-phthalate 2.1 M, Na2HPO4/KH2PO4 4.9 M (60-155 mL)"
    temperature_c: [100]
    ph: {initial: [5.5], measured_during: false, final: null}
    time_min: [120]
    species: [norfuraneol]
    kind: endpoint
    location: "Table p. 1019 (yield %, read from a 250 dpi render); dossier sec. 4a"
    status: digitised_in_dossier
    quant_class: unknown
    simulable: false
    condition_quote: "heating 1.77 mmole of D-ribose 5-phosphate in various buffer solutions at pH 5.5 for 2 hours"
    flags: [single_time_point, yield_method_not_stated, no_replicates, ph_reading_time_not_stated, concentrated_buffer_extreme_ionic_strength]
    note: "NF 43-47 % per R5P, amine-free; buffer identity barely matters. Bounds R5P loss k >= 0.005 min-1 at 100 C. Engine has no R5P species. Contrast with tonsbeek1969_text_r5p_model (dilute, unbuffered: no NF)."
  - id: peer1968_text_preparative_acetate_reflux_yield
    dossier: peer1968_extraction.md
    paper: "Peer & van den Ouweland 1968 Recl Trav Chim Pays-Bas 87:1017"
    lab: "Unilever Duiven (Peer / van den Ouweland)"
    system: "D-ribose 5-phosphate 0.18 mol (~30 mM) from acid-resin hydrolysis of IMP/GMP, no amine"
    medium: "~3.5 M sodium acetate / 0.7 M acetic acid, ~5.5-6.5 L, open reflux in air"
    temperature_c: [100]
    ph: {initial: [5.4], measured_during: false, final: null}
    time_min: [120]
    species: [norfuraneol]
    kind: endpoint
    location: "Experimental p. 1019 (UV 43 %, GC 45 %, isolated 33.8 %); dossier sec. 4b"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "After refluxing for 2 hours the reaction mixture was cooled and extracted"
    flags: [single_time_point, ph_nominal_from_recipe, reflux_above_100c, no_replicates]
    note: "Reproduces the Table acetate row at 100x scale. pH 5.4 is Henderson-Hasselbalch from the recipe, not measured. NF lambda max 289 nm, epsilon 9200 (MeOH)."
  - id: peer1968_text_1_14c_ribose5p_norfuraneol_label_position
    dossier: peer1968_extraction.md
    paper: "Peer & van den Ouweland 1968 Recl Trav Chim Pays-Bas 87:1017"
    lab: "Unilever Duiven (Peer / van den Ouweland)"
    system: "D-[1-14C]ribose 5-phosphate (1.51e6 dpm/mmol), conditions as the acetate preparation (not restated)"
    medium: "NaOAc/AcOH buffer (presumed; not restated for the tracer run)"
    temperature_c: [100]
    ph: {initial: [5.5], measured_during: false, final: null}
    time_min: []
    species: [norfuraneol]
    kind: other
    location: "text pp. 1018 and 1020 (specific activities, Kuhn-Roth acetic acid 4.5 %); dossier sec. 4c"
    status: text_only
    quant_class: other
    simulable: false
    condition_quote: "the methyl group does not originate from the aldehyde group of the ribose"
    flags: [labelling, condition_not_stated, c2_assignment_inferred]
    note: "R5P route: ribose C1 kept (98 %), not in methyl/C5 (4.5 %); methyl from ribose C5 via beta-elimination of phosphate. Opposite topology to the amine / 1-deoxyosone route. Relevant to E6/E11 labelling cells only if R5P is present."
```

## 8. Supplementary material

None referenced. Nothing model-relevant is missing; the companion paper (Recl. Trav. Chim. 87:1011, aldopentose +
secondary amine salts, the methyl-from-C1 result this paper contrasts with) is a separate article, not SI, and is
not on disk.
