# Cämmerer & Kroh 1995 (received 1994) — EXTRACTION (elementary composition of model melanoidins vs reaction conditions)

**Source on disk:** `data/articles/Cammerer1994.pdf` (5 pp., scanned journal pages with an OCR text
layer; the OCR is noisy — "SO.28" for 50.28, "500" for 5.00, "4290" for 42.90, "17O"C" — so Tables 1-3
were verified against 150/200-dpi rasters of pp. 3-4, from which the values below are taken). Read-only
extraction, 2026-09-07. Filed under the owner's name `cammerer1994` although the printed year is 1995.

## 0. Identity

| field | value |
|---|---|
| Title | "Investigation of the influence of reaction conditions on the elementary composition of melanoidins" |
| Authors | B. Cämmerer & L. W. Kroh (Technical University Berlin, Institute of Food Chemistry) |
| Venue | Food Chemistry 53 (1995) 55-59; received 8 April 1994; revised/accepted 4 July 1994 |
| DOI as printed | **none printed** (pre-DOI Elsevier issue; only the code 0308-8146/95/$9.50). The registered DOI is 10.1016/0308-8146(95)95786-6 — not verified on the page. |

## 1. Why it matters to the model

It is **not a kinetics paper**: it contains no rate constant, no activation energy and no time course.
What it gives is the **stoichiometry of the melanoidin sink** — how many moles of sugar (or C3 fragments)
are incorporated into the non-dialysable polymer per mole of amino acid, and how many water molecules are
lost per sugar — as a function of sugar type, sugar:glycine ratio, temperature and water (solvent-free
170-180 C vs aqueous 60-100 C, pH 5 / 7). The repo's melanoidin step (Martins 3-DG + Gly -> Mel, and the
extinction-coefficient bookkeeping from Brands 2002b / Martins 2003c) needs an assumed sugar:amine
stoichiometry to convert "sugar units incorporated" to mass; this paper is the primary source for that
number in glucose/glycine, and it shows that in water at 100 C the ratio is ~1.2-1.3 sugar per glycine
while dry at 170-180 C it is ~2.2. Also the structural claim: the polymer is built from dicarbonyl
(3-deoxyhexosulose) units, consistent with the trunk routing melanoidin through 3-DG.

## 2. Methods as they matter to a model

- **Solvent-free:** carbohydrate + amino acid mixed 1:1 molar, heated "in a flat sheet" for 10-40 min at
  170 or 180 C; solid ground and dialysed. (Table 1/2 footnote: "solvent free, 10 min, 180 C"; Table 3
  first row: "170 C/20 min" — see caveats.)
- **Aqueous (a):** 0.1 M sugar + 0.1 M amino acid (glycine or DL-phenylalanine) refluxed 10 h (i.e. 100 C);
  pH 5 or 7.3 held constant with 0.1 M NaOH via a sterilisable electrode ("10-20 ml during the total
  heating time"). Volume of the charge not stated.
- **Aqueous (b):** 0.1 M solution 160 h at 60 C, pH 5 held constant. Also a 90 C / 22 h / pH 5 row
  credited to Wedzicha & Kaputo (1992).
- Lyophilised residues dialysed (cellulose, MWCO 12-14 kDa, 5 g in 1 L water, water changed every
  8-10 h, 136 h total), then freeze-dried. **Only the non-dialysable (> 12 kDa) fraction is analysed.**
- CHNS microanalysis (Leco CHNS 932, ~1.5 mg); HPLC (Nucleogel aqua-OH 40, water eluent, DAD 190-450 nm) to
  follow dialysis; IR (KBr); CP-MAS NMR cited from Engelke et al. 1994.
- Calculation (Wedzicha & Kaputo 1992): melanoidin formula C(la+pb-x) H(ma+qb-2y) O(na+rb-2x-y) H(b) with a
  carbonyl molecules (C_l H_m O_n) and b amino acids (C_p H_q O_r N); x, y = mol CO2, H2O lost; solved per
  atom N to give a (mol sugar per mol amino acid) and y/a (mol water per mol sugar).
- Replicates not stated; no uncertainties printed anywhere.

## 3. Tables (verbatim, %, from the raster)

**Table 1.** "Microanalysis data (%) of non-dialysable melanoidins from D-carbohydrates/glycine (1:1) model
systems^a" — a: "Reaction conditions: solvent free, 10 min, 180 C." b: "Mol sugar (C6- or 2C3-fragments)
which is incorporated into the polymer per mol amino acid." c: "Mol water which is liberated per mol sugar."

| Sugar | C | H | N | O | a^b | y/a^c |
|---|---:|---:|---:|---:|---:|---:|
| Ribose | 50.28 | 6.28 | 2.41 | 41.03 | 4.74 | 1.69 |
| Glucose | 53.42 | 5.38 | 4.26 | 37.61 | 2.19 | 3.05 |
| Fructose | 47.39 | 6.42 | 5.00 | 41.19 | 1.53 | 3.54 |
| Maltose | 42.34 | 6.09 | 1.61 | 49.96 | 2.34 | 0.76 |
| Lactose | 42.03 | 5.41 | 4.37 | 48.19 | 0.71 | 2.31 |
| Sucrose | 42.90 | 5.53 | 5.80 | 45.77 | 0.52 | 2.98 |

**Table 2.** "Microanalysis data (%) of non-dialysable melanoidins from glucose/glycine model systems^a"
(same footnotes; molar ratio = glucose:glycine).

| Molar ratio | C | H | N | O | a | y/a |
|---|---:|---:|---:|---:|---:|---:|
| 8:2 | 50.32 | 6.19 | 4.15 | 39.34 | 2.16 | 2.31 |
| 1:1 | 53.42 | 5.39 | 4.26 | 37.61 | 2.19 | 3.05 |
| 2:8 | 47.81 | 6.21 | 3.50 | 42.48 | 2.43 | 1.92 |

**Table 3.** "Microanalysis data (%) of non-dialysable glucose/glycine (1:1) melanoidins produced under
various reaction conditions" — a: mol sugar per mol amino acid; b: mol water per mol sugar; c: "Wedzicha &
Kaputo (1992)".

| Conditions | C | H | N | O | a | y/a |
|---|---:|---:|---:|---:|---:|---:|
| 170 C / 20 min | 53.42 | 4.38 | 4.26 | 37.61 | 2.19 | 3.05 |
| 100 C / 10 h / pH 5 | 55.58 | 5.38 | 6.97 | 32.07 | 1.28 | 3.75 |
| 100 C / 10 h / pH 7 | 49.05 | 5.25 | 6.12 | 39.57 | 1.20 | 3.10 |
| 90 C / 22 h / pH 5^c | — | — | — | — | 1.06 | 2.96 |
| 60 C / 160 h / pH 5 | 43.02 | 4.78 | 6.94 | 45.25 | 0.75 | 2.89 |

**Table 4** (IR, cm^-1; aqueous 100 C/10 h vs solvent-free 170 C/20 min): 3340/3420 (O-H, N-H, s);
2920/"1930" (aliphatic C-H, m — the 1930 is almost certainly a misprint for 2930); 1700/1710 (C=O, m sh);
1615/1630 (C=N, C=C); 1420/1420; 1375/1380; 1215/1200; 1050/1020. No kinetics.

**Fig. 1:** HPLC of dialysate (8 h) and retentate (112 h) for glucose/glycine 1:1, 10 h / 100 C, pH 5 —
qualitative. **Fig. 2:** proposed polymer backbone from 3-deoxyhexosuloses + amino acids (enamine/imine
chain, R = H or saccharide). No scheme with rate constants exists in the paper.

Text statements with numbers: "In all experiments carried out, three mol water were eliminated per mol
carbonyl compound"; "in boiling reaction mixtures (100 C) the ratio of amino to carbonyl compound falls to
around 1:1"; "in solvent-free milieu ... at least two mol glucose ... per mol amino acid"; ribose "more
than 4 mol"; "Even a significant excess of glycine (2:8) ... can hardly affect the molecular composition".

## 4. What the repo could take

No FIT rows (no measured rate). Stoichiometric constants and directional claims for the melanoidin sink:

| item | value | use |
|---|---|---|
| sugar units per glycine in glucose/glycine melanoidin, water, 100 C, 10 h | 1.28 (pH 5), 1.20 (pH 7) | the mass/N bookkeeping of the Mel species: ~1.2 glucose-C6 per amine at cook temperatures in water |
| same, 60 C / 160 h / pH 5 | 0.75 | ratio falls at lower T |
| same, solvent-free 170-180 C | 2.19 | dry/high-T doubles sugar incorporation |
| water lost per incorporated sugar | 2.9-3.75 (all conditions; "three") | melanoidin C6 unit ~ C6H6O3 + amine, i.e. a 3-DG-derived (triple-dehydrated) unit — consistent with routing Mel through 3-DG |
| ribose vs glucose incorporation (dry, 180 C) | 4.74 vs 2.19 (2.2x) | pentose melanoidins are sugar-richer: relevant to the sulfur/pentose lane's Mel stoichiometry |
| fructose vs glucose (dry) | 1.53 vs 2.19 | |
| glucose:glycine 8:2 vs 1:1 vs 2:8 | a = 2.16 / 2.19 / 2.43 | polymer composition insensitive to reactant ratio (directional: the Mel stoichiometry can be a constant, not a function of charge) |
| pH 5 vs pH 7 at 100 C | a = 1.28 vs 1.20; N 6.97 vs 6.12 % | weak pH dependence of composition |
| N content of glucose/glycine melanoidin, water 100 C | 6.1-7.0 % w/w | converts melanoidin mass to bound-glycine equivalents (~0.45-0.5 mmol N per 100 mg) |

Directional claims: sugar incorporation per amine rises with temperature (0.75 at 60 C -> 1.06 at 90 C ->
1.2-1.3 at 100 C -> 2.2 at 170-180 C) and with dryness; C content rises and N falls with harsher
conditions (43 % C at 60 C vs 53-56 % at 100-180 C); the 100 C polymer is nitrogen-richer (7 % N) than the
dry one (4.3 %).

## 5. Caveats

- Scanned PDF; OCR text layer is unreliable for numbers — all table values above were read from the
  raster. Table 3 row 1 repeats Table 1's glucose row (C, N, O, a, y/a identical) but with H = 4.38 vs
  5.38 and conditions "170 C/20 min" vs "180 C, 10 min": one of the two captions or the H value is a
  misprint; the paper does not resolve it.
- Only the > 12 kDa retentate is characterised; the authors say part of "the higher and also the highest
  molecular weight fraction were washed out" during dialysis, so a and y/a describe a dialysis-defined
  fraction, not all colour.
- No replicates, no uncertainties, no time resolution; aqueous runs are 10 h at reflux with pH held by
  NaOH titration (10-20 mL of 0.1 M NaOH added to an unstated volume), i.e. sodium accumulates.
- The 90 C row is another group's data (Wedzicha & Kaputo 1992).
- a and y/a are model-derived from CHN via the Wedzicha-Kaputo formula, which assumes losses only as CO2
  and H2O; the "2C3-fragments" alternative means a counts C6 equivalents, not intact glucose.
- Nothing about rates, so nothing here can move a trunk constant; use only for the Mel species' formula
  weight and nitrogen bookkeeping.
