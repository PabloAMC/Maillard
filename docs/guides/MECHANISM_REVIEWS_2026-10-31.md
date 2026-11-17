# Mechanism reviews, 2026-10-31: what the literature says about the six core problems

*Six independent read-only reviews (one per core problem), each reading the dossiers and PDFs on disk and searching
for what is missing. The raw reports are in `docs/history/mechanism_reviews_2026-10-31/`. Numbers marked "derived"
are the reviewers' own arithmetic on keyed data; nothing here is a fit. This synthesis ranks what changes a decision
in `STATUS_AND_PLAN.md` sec. 5. Probabilities are the reviewers', kept where they were argued.*

## 1. The findings that change the plan

| # | problem | finding | P | consequence |
|---|---|---|---|---|
| 1 | deoxyosone hub | both deoxyosones are short-lived and quasi-steady next to their parent (Hofmann 1999, Zhai 2023 JAFC, Zhai 2021 pH 5.5; Martins M4: 1-DG ~1 min, 3-DG ~30 min); the engine's 3-DX lives ~35 h, so the ratio climbs | 0.65 | the hub step starts from Martins & van Boekel 2005 M4 as its null |
| 2 | hub / cysteine | free cysteine traps the aldehydic 3-DX, faster at high pH -- the pre-declared signature (vii) (+xylose / -xylose double ratio of TDP / DPO, measured 0.85 / 0.92 / 1.22 / 2.71 at pH 5.5 - 8) is reproduced at no pH 7 / 8 point by any arm (all ~0.73) | 0.45 | a cysteine-3-DX adduct with a free-amine factor enters the hub step |
| 3 | cysteine fate | the dilute fed-TTCA pots are kinetically far from equilibrium: adding xylose RAISES free cysteine (Zhai 2021 Fig. 6), which no equilibrium model allows; the same lab's 83 mM pot forms TTCA 8 - 30x faster per mol -- ring hydrolysis / recombination catalysed by solutes and buffer | 0.55 | B7b becomes catalysed rate laws k = k0 + k_aa [AA] + k_P [phosphate], dilute and concentrated pots fitted together, the xylose contrast must pass |
| 4 | thiols | the Cerny / Schieberle conflict most likely comes from the MFT LOSS law (the engine destroys 93 % of the MFT it makes at 100 C, so only an accelerating precursor -- norfuraneol -- can give Schieberle's rise); with an hours-long MFT lifetime (CF-013) a steady intact-C5 supply already gives the lag | 0.40 | a cheap dry check first: k_mft_decay x0.1 - 0.2 at the Cerny-passing corner |
| 5 | thiols | Liu 2023's "19 % fragment recombination" at 168 C is an EI-MS artefact: its MFT channels match the NIST spectrum's [M-H]+ / M+ (0.30 vs 0.28), as do furfural's and 2-methylfuran's | 0.8 | Liu's signature (thiol prereg sec. 11) is withdrawn; the "route mix changes with temperature" finding loses its high-T anchor |
| 6 | pyrazines | the lane uses a pyrazine-YIELD constant as the Strecker FLUX: fed-dicarbonyl pots turn ~1 % of the aminoketone into pyrazine (the rest to pyrroles, Adams 2008); Hofmann 2000 puts the true Strecker flux ~1.9 dex above B18's constant | 0.55 | pyrazines follow straight after the forced-input Strecker layer (roadmap step 2), with an aminoketone -> pyrrole sink |
| 7 | Strecker | a fed-donor yield is a branching share, not a rate (MGO is gone in ~15 min); the rate law is mass action on TOTAL amino acid with a shallow pH factor (<= 0.07 dex / pH); O2 opens a SEPARATE acid route (PA -> PAA < 1 %); lipid enals need O2 for the aldehyde itself | 0.5 | adopted in the forced-input prereg; a slope-1 free-amine law is refuted by > 3 dex |
| 8 | matrix | protein-bound lysine is the dominant amine in a real alt-protein product (~60 mmol/kg at 16 % pea isolate vs ~1 mmol/kg free amino acids; ~1.9x glycine's reactivity per amine): ~90 % of the sugar-amine flux | 0.40 (largest gap) | the composable network's amine classes (roadmap step 6); product numbers need designed data (P 0.85) |
| 9 | matrix | **code bug:** the B20 glycation arm returns lysine on fructoselysine decay; Nguyen 2016's model does not | -- | CF-072; fixed in its own PR |
| 10 | process | sec. 13 of the B7 prereg wrongly said Zhai 2021's pots were phosphate-buffered (they are water): phosphate catalysis is a live rival for the Davidek / Zhai tension | -- | corrected (dpo prereg, CORRECTION note) |

## 2. The cheapest discriminating checks (no new data)

1. The MFT-loss dry check (finding 4): ~30 min, one core.
2. Signature (vii) and the +-xylose free-cysteine contrast as must-pass tests of the hub step (findings 2, 3): declared.
3. The methylpyrazine ratio in Leahy Table I: methylpyrazine^2 / (pyrazine x 2,5-DMP) = 11.7 / 14.4 / 17.3 at 95 / 85 /
   75 C against 4 for one shared condensation constant -- methylpyrazine has an extra route (finding 6).
4. 2-mercapto-3-pentanone / MFT from ribose as a presence / absence signature across temperature (absent at 95 C,
   present at 145 C).
5. The Weenen 2001 vs Hofmann 2000 Strecker yield from MGO at pH ~3 (21 % vs ~1 %): a cross-laboratory signature
   (CF-073). Corrected by the forced-input design review: an aldehyde-loss term CANNOT explain it (k_loss ~4e-4 /min
   removes <= 5 % in 2 h); donor concentration or a reversible aldehyde-amine pool drained by steam stripping can.

## 3. Papers to download

32 candidates; the list is kept privately.

## 4. Gaps no search filled (designed data or authors)

No time-resolved MFT isotope labelling at 95 - 100 C; no fed-DDP experiment; no H2S time course in a sugar + cysteine
pot; no DPO decay in water vs phosphate; no glycation kinetics on pea or soy protein; no thiol capture by a plant
protein at process temperature; no measured aminoketone-dimerisation or dihydropyrazine-oxidation constant; no
Strecker step measured in water (G1). Zhai 2021's 13C5-xylose label fraction in TTCA would give the formation constant
at 10 mM directly (only MS fragments are printed -- ask the authors).
