# Pesek & Frost 1975 — EXTRACTION (2-substituted thiazolidine-4-carboxylic acids from L-cysteine, D,L-penicillamine and 2-aminoethanethiol with aliphatic and aromatic aldehydes, 0.5 M in D₂O, followed by 60 MHz ¹H NMR in 8 / 4 M HCl, neutral D₂O and 0.25-3 M NaOH, over hours to days; α-H/D exchange, line broadening and a directly observed ArCH=N imine; temperature NOT stated)

### A QUALITATIVE STABILITY MAP AT pH EXTREMES. Aliphatic 2-substituted TCAs (the class TTCA belongs to) are unchanged for days in 8 M HCl and in neutral D₂O, and decompose to cysteine + aldehyde only in strong base (> 3 M NaOH). Three rate constants, all at pD ≥ 11.8 or in 0.5-3 M NaOH, all ring-chain or Schiff-base steps of ARYL or C1-C2 alkyl substrates. Nothing at pH 4-8, nothing at any stated temperature.

**Read 2026-11-01 from `data/articles/pesek1975.pdf`** (7 pp., Tetrahedron 31, 907-913; Acrobat Capture OCR of a
scan). Prose from the `pdftotext -layout` layer, **checked against 130 dpi renders of all pages** (scratchpad
`dos6/pes-*.png`); Fig. 5 re-read from a 300 dpi crop (`pesf5-5.png`). The OCR garbles every exponent; exponents
below are from the renders. Figure values read by me: ±0.02 on Fig. 5, ±0.3 s⁻¹ on Fig. 7, ±0.005 on Fig. 8.

**Legend.** **[V]** verbatim · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative claim ·
**[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| model need | what Pesek 1975 gives | tag |
|---|---|---|
| TTCA-class (aliphatic C-2) stability vs pH | unchanged by NMR "even after several days" in neutral D₂O; no change in up to 8 M HCl; decomposition to cysteine + aldehyde in > 3 M NaOH; with isobutyraldehyde + cysteine in 1-3 M NaOH "decreasing amounts of ring compound with increasing base concentration" | [Q] |
| forward rate, acid | aliphatic aldehyde + cysteine (0.5 M) in 4 M HCl/D₂O reaches the thiazolidine spectrum in **12-24 h** (temperature not stated) | [M] (time only) |
| ring-chain (C-S / ring ⇌ imine) rate | aryl TCAs (2-thienyl, 2-phenyl), 0.5-1.0 M NaOH: interconversion k₁ = **25 M⁻¹ s⁻¹** (OH⁻-catalysed, NMR line broadening, Fig. 7) | [M] |
| α-H exchange via the open imine (aliphatic) | 2-methyl TCA pD 11.8: **6.0 × 10⁻⁵ s⁻¹**; 2-ethyl TCA pD 12.2: **7.3 × 10⁻⁵ s⁻¹** (Fig. 5) | [M] |
| Schiff-base hydrolysis | 2-thienyl, 2 and 3 M NaOH: k_h = **2.2 × 10⁻³ s⁻¹** (Fig. 8; initial-rate plot) | [M] |
| K, Ea, buffer catalysis | **none** | — |
| temperature | **not stated** (Varian A-60A "equipped with a variable temp probe"; "a lowering of the temperature produces sharper lines") | [!] |

**Answers.**

1. **`k_ttca_cys`: no bound** (confidence ~97 %). Nothing is measured at pH 4-8 or at a stated temperature; the three
   constants are for aryl rings or for alkaline D₂O (pD ≥ 11.8, 0.5-3 M NaOH).
2. **Direction of the pH law of RELEASE, qualitatively** [Q → Z]: the aliphatic TCAs do not release cysteine in 8 M HCl
   (days), consistent with Kallen 1971b's TC k_r falling with acidity (a ≈ −0.5 between pH 2.8 and 5.6,
   `kallen1971b_extraction.md`) and against an acid-accelerated release (E2 / E4 acid limb, `a` > 0). Weak: NMR at
   0.5 M, where the equilibrium lies on the ring side whatever the rate, so "no change" also reflects K.
3. **Ring opening to the imine is base-catalysed** [M/Q] at pH ≥ 13 (aryl, k₁ = 25 M⁻¹ s⁻¹ in [OH⁻]); Pesek's
   summary item (3). This is the C-S ring-chain step, not release.
4. **Order of thiazolidine stability: formaldehyde ≫ aliphatic aldehydes > aromatic aldehydes** [Q, with C refs 26-27]:
   the direction of the ~10⁴ gap between Kallen's TC K and Cao's TTCA K, no magnitude.
5. **CF-076 E1 / E4 (catalysis), E3 (adducts), MFT loss law**: nothing. Cysteine in 3 M NaOH is unchanged for two weeks
   by NMR (no disulfide) unless H₂O₂ is added [Q]: oxygen alone did not oxidise it on that timescale at that pH.

## 2. Identity

| field | value |
|---|---|
| Title | "Decomposition of Thiazolidines in Acidic and Basic Solution. Spectroscopic Evidence for Schiff Base Intermediates" |
| Authors | J. J. Pesek and J. H. Frost, Department of Chemistry, Northern Illinois University, DeKalb |
| Venue | *Tetrahedron* **31**, 907-913 (1975); received USA 20 February 1974, UK 24 October 1974 |
| DOI | not printed; PDF metadata PII 0040-4020(75)80099-8. Crossref (2026-11-01): **10.1016/0040-4020(75)80099-8**, same title, Pesek & Frost, Tetrahedron 31(8) 907-913. Matches the brief |
| Funding | Petroleum Research Fund (ACS) |
| Paper type | NMR (60 MHz) and preparative; no UV kinetics, no buffered pH |

## 3. Conditions

- **Substrates** [M]: 2-substituted (Ph, cyclohexyl, Me, Et, i-Pr, Pr) TCAs from cysteine; 2-(2-thienyl) TCA and TC
  (Aldrich); 2-(Ph, thienyl, Pr) 5,5-dimethyl TCAs from penicillamine; mercaptals (cysteine or 3-mercaptopropionic acid
  + benzaldehyde / thiophene carboxaldehyde, 4 M HCl).
- **Media** [M]: all 0.5 M in D₂O; 8 M HCl, 4 M HCl, neutral, 0.25 / 0.5 / 0.75 / 1 / 3 M NaOH; formation studies
  0.5 M in 4 M HCl; NaOH 1-3 M for isobutyraldehyde + cysteine. **No buffer; pH / pD measured only for Fig. 5
  (11.8, 12.2).** Isobutyraldehyde ~30 % as aldol under the formation conditions [M].
- **Temperature, vessel, O₂**: not stated (NMR tube) [!].
- **Analysis**: ¹H NMR, chemical shifts vs external TMS; MS and elemental analysis of the mercaptals.

## 4. Quantitative results

- **Fig. 5** [M/F]: α-H peak area vs time, 2-methyl TCA (pD 11.8) and 2-ethyl TCA (pD 12.2), 0-6 h; printed first-order
  constants 6.0 × 10⁻⁵ and 7.3 × 10⁻⁵ s⁻¹. My read of the plotted slopes: 2.5 × 10⁻⁵ and 2.8 × 10⁻⁵ s⁻¹ if the
  ordinate is ln as labelled; **5.7 × 10⁻⁵ and 6.5 × 10⁻⁵ s⁻¹ if it is log₁₀** [Z][!]. Comparison values cited:
  0.2-8.9 × 10⁻⁵ s⁻¹ for isobutyraldehyde-amino-acid imines [C, ref 17].
- **Fig. 7** [M/F]: 1/τ of the C-2 proton rises linearly from ~8 to ~13.4 s⁻¹ over [OH⁻] 0.70-0.93 M (2-thienyl, open
  circles) and ~7.5 to ~12.5 s⁻¹ (2-phenyl); slope printed 25 M⁻¹ s⁻¹ for both (my read ~22-23).
- **Fig. 8** [M/F]: d[9]/dt · 1/[7] ≈ 0.13 min⁻¹ flat for the first ~2-3 min in 2 and 3 M NaOH, then falling; printed
  k_h 2.2 × 10⁻³ s⁻¹ (= 0.13 min⁻¹, consistent). Cited ranges 1.7 × 10⁻²-1.7 × 10⁻¹ s⁻¹ and 5 × 10⁻⁵-5 × 10⁻¹ s⁻¹
  for other aryl Schiff bases [C].
- **Fig. 1** [Q]: 0.5 M 2-ethyl TCA in 0.75 M NaOH: after 24 h the C-2 multiplet collapses to a singlet and the side-
  chain α-H signal disappears (exchange via the open imine); no loss of the C-4 proton.
- Eq (1)-(2) [M]: d[7]/dt = k₁[6] − k_h[7] − k₋₁[7]; d[9]/dt = k_h[7] − k₄[9]²[OH⁻]² (Cannizzaro removal of aldehyde).
  k₁ and k₋₁ separately "no way to estimate" [V-short].

## 5. Mechanism evidence and its kind

- Northern Illinois group (Pesek, Frost). Kind 4 (directly observed intermediate: the ArCH=N imine at ~8.3 ppm in
  1-1.5 M NaOH; the mercaptal isolated and characterised) and kind 3 (trapping: phenylhydrazine precipitates the
  aldehyde from neutral solution, showing ring-chain tautomerism) [M/Q]. Kind 6 only across extreme media.
- **Proposed, not shown** [P]: a sulfonium-ion intermediate (-S⁺=CH-) for mercaptal formation in strong acid with aryl
  rings (Fig. 4); Kallen 1971b had excluded it for formaldehyde at ordinary acidity.

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| CF-076 E2 / E4 acid limb | **weakly CONTRADICTS** (qualitative) | no acid-driven release of aliphatic TCAs, days in 8 M HCl |
| CF-076 E1 catalysis | LEAVES UNCHANGED | no buffer used |
| CF-076 E3 / O₂ sink | LEAVES UNCHANGED | cysteine stable in 3 M NaOH (no oxidant) is not a pot condition |
| evidence map 2.12 cysteine fate / TTCA | LEAVES UNCHANGED | |
| engine `k_ttca_cys`, MFT loss law | LEAVES UNCHANGED | no constant at relevant pH / T |

## 7. Disagreements found

1. **Fig. 5 axis vs printed constants** [Z][!]: slopes reproduce 6.0 / 7.3 × 10⁻⁵ s⁻¹ only if the "ln peak area" axis
   is log₁₀ (x2.30). Most likely a labelling error; either way the constants are as printed.
2. **With Kallen's alkaline equations** [Q]: the authors state Kallen's equations "which imply alkaline stability do not
   apply to 2-substituted" TCAs (TC itself is stable in strong base here, consistent with Kallen). Chemistry (C-2
   substitution), not a contradiction of Kallen's data.
3. With Fife 1991 (`fife1991_extraction.md`): Fife writes ring opening is favourable only below pH ~10 (thiol
   un-ionised); Pesek sees OH⁻-catalysed ring ⇌ imine exchange at pH ~14 for aryl TCAs. Different quantity (rate of
   exchange vs fraction open) and substrate; not a contradiction.

## 8. Ledger-ready measurements

None worth keying for the engine. If a thiazolidine-stability row is ever wanted: the qualitative "aliphatic 2-substituted
TCA unchanged in 8 M HCl and neutral D₂O for days at 0.5 M" (sec. 1), with its missing temperature flagged.
