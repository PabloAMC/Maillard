# Voigt, Smuda, Pfahler & Glomb 2010 — EXTRACTION
### *"Oxygen-Dependent Fragmentation Reactions during the Degradation of 1-Deoxy-D-erythro-hexo-2,3-diulose"*
### Synthetic 1-deoxyglucosone (1-DG) + L-alanine, 42 mM each, 0.1 M phosphate pH 7.4, **37 and 50 °C**, 0-48 h, aerated vs deaerated: lactic acid (Fig. 1), γ-pyranone (Fig. 2), acetaldehyde (Fig. 3), glyceric acid (text only). [13C]glucose isotopomers + alanine (50 °C, 7 d, air): lactic acid carbon origin (Table 1) and glyceric acid origin. Fed γ-pyranone gives the same acids, also O2-dependent. **No acetic acid, no 100-140 °C, and the deaerated arm carries 1 mM DTPA that the aerated arm lacks** (protocol of ref 8).

**Read 2026-10-31 from `data/articles/voigt2010.pdf`** (7 pp., JAFC 58:5685-5691; born-digital, clean text layer).
Figs. 1-3 and Table 1 checked against 200 dpi renders of pp. 5687-5688 (read error ~ ±3 mmol/mol on Figs. 1-2,
~ ±0.5 mmol/mol on Fig. 3). Figs. 5-7 are schemes. The deaeration protocol was checked in ref 8 (Voigt & Glomb 2009,
`data/articles/voigt2009.pdf`, no dossier yet).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Questions: triage item 7 (oxygen as a reagent) and the evidence-map row "O2: oxidative alpha-cleavage acids —
SINGLE (Nestlé); second lab voigt2010" (`docs/guides/EVIDENCE_MAP_2026-10-31.md` sec. 1.3, 2.12).

| need | what Voigt 2010 gives | tag |
|---|---|---|
| O2 / argon contrast | **Lactic acid at 48 h: 166 vs 26 mmol/mol 1-DG (37 °C), 118 vs 28 (50 °C)**, aerated vs deaerated (×6.4 and ×4.2 [Z]). Glyceric acid: 266 mmol/mol at 37 °C aerated; "significant increase" under aeration, deaerated value not printed. γ-pyranone runs the other way: 146 (deaerated) vs 67 (aerated) at 12 h, 37 °C; 114 vs 103 at 50 °C. Acetaldehyde ×1.5-1.6 with air [F]. Fed γ-pyranone: acids "preferred under aerated conditions" | [M] / [F] / [Q] |
| what "deaerated" means | Not described here ("according to ref 8"). Ref 8 (Voigt & Glomb 2009): phosphate + **1 mM DTPA**, helium-degassed buffer, samples argon-flushed. Aerated pots: no DTPA stated. **O2 is confounded with chelation**, as in Smuda 2010 | [C] / [!] |
| phosphate vs other buffer | **None.** All pots 0.1 M phosphate pH 7.4 | [M] |
| metals / chelators | None varied in this paper, except the DTPA that rides with deaeration [C]. Authors cite metal-catalysed oxidation (ref 17) | [C] / [!] |
| rate constants / Ea | **None.** 1-DG half-life "about 0.5 h" at 37 °C is cited from ref 7. Two temperatures, but yields fall with temperature: lactic acid 166 -> 118 (aerated) | [C] / [M] |
| isotope-label shares | [13C]glucose (C-1, 2, 3, 5, 6; no C-4) + Ala, 50 °C, 7 d, air. Lactic acid: **48 % from C1-C3 (COOH at C-3), 36 % from C3-C5 (COOH at C-5), 8 % + 8 % from C4-C6 (COOH at C-4 or C-6)** (Table 1). Glyceric acid: label from C-5 and C-6, never from C-1 to C-3, so "exclusively" from the lower three carbons | [M] |
| transfer to 100-140 °C | **Not shown.** 37 and 50 °C only. The O2-dependent lactic excess (aerated minus deaerated) falls from 140 to 90 mmol/mol between 37 and 50 °C [Z]; the γ-pyranone O2 contrast nearly vanishes at 50 °C. Any statement at 100-140 °C is an extrapolation | [Z] / [!] |

**Answer.** Voigt 2010 is a **second laboratory (Glomb, Halle) for the DIRECTION "O2 raises the C3 acids from
1-DG"**, at 37-50 °C (confidence ~0.85 in the direction, given the DTPA confound below). It is **not** a second
laboratory for Davidek 2006b's quantity: Davidek measures acetic acid (argon 30.5, air 32.3 mmol/L, O2-independent) and
18O incorporation; Voigt measures lactic and glyceric acid and no acetic acid. The two do not conflict:

1. Davidek 2006b's ~20 % (lactic) and 1/3 (glyceric) 18O labelling from 18O2 is not a measure of how much of the acid
   is O2-dependent if the route is Voigt's [P] Fig. 5: O2 oxidises the reductone to a tricarbonyl (the O2 oxygen leaves
   as H2O2), then hydrolytic cleavage supplies the carboxyl oxygen from water. Under that route an O2-dependent acid can
   carry no 18O [Z]. So "only a minor share is oxidative" (Davidek) and "×6 with air" (Voigt) can both hold.
2. The DTPA confound: the deaerated arm also removes catalytic metals. The ×6 is an upper bound on the pure O2 effect if
   metals matter, and the right number for "air + adventitious metals vs argon + chelator" (~0.7 that metals carry a
   material part, by analogy with thornalley1999's DETAPAC effects [Z]).
3. Temperature runs against transfer: from 37 to 50 °C the O2-dependent excess falls by ~35 % and γ-pyranone loses its
   O2 dependence. Nothing here supports carrying ×6 to 100-140 °C; ~0.6 that the oxidative share is smaller at 120 °C
   in an air-headspace vial (lower O2 solubility, faster non-oxidative 1-DG sinks) [Z].

## 2. Conditions

| item | as printed |
|---|---|
| 1-DG | synthesised per ref 7 (Glomb & Pfahler 2000); deprotection details not repeated; 42 mM nominal [!] |
| amine | L-alanine 42 mM (β-alanine in one acetaldehyde control) |
| buffer / pH | 0.1 M phosphate (K2HPO4 / KH2PO4), pH 7.4; pH not re-measured [!] |
| temperature / time | 37 and 50 °C; 0-48 h (Figs. 1-2), 0-12 h (Fig. 3) |
| atmosphere | aerated vs deaerated "according to ref 8" = 1 mM DTPA + helium-degassed buffer + argon (ref 8) [C] |
| label pots | [1-, 2-, 3-, 5-, 6-13C]glucose 12.6 µmol + Ala 12.6 µmol in 0.3 mL 0.1 M phosphate pH 7.4 (42 mM each [Z]), 50 °C, 7 d, aerobic |
| γ-pyranone pots | synthetic γ-pyranone + Ala, 42 mM, same buffer, 37 / 50 °C, aerated / deaerated |
| analytics | acids, γ-pyranone: GC-MS after silylation (SD < 10 mmol/mol, CV < 5 %); glyceric label: GC-LCI-MS (m/z 323 / 324); lactic label: LC-MS/MS MRM (m/z 89/90 -> 43/44 and 45/46 locate the label in CHOH-CH3 vs COOH); aldehydes: post-incubation OPD trapping to benzimidazoles, HPLC-UV 272 nm; ≥ 3 determinations |
| aldehyde recovery | acetaldehyde 27 %, glycolaldehyde 38 %, glyceraldehyde 76 % by the OPD trap; whether yields are recovery-corrected is not stated [!] |

## 3. Quantitative results

**Fig. 1, lactic acid (mmol/mol 1-DG)** — printed [M]: 37 °C aerated 87 at 2 h, 166 at 48 h; deaerated 26 at 48 h;
50 °C aerated 118, deaerated 28 at 48 h. Read [F] (±3):

| h | 37 °C air | 50 °C air | 37 °C deaer. | 50 °C deaer. |
|---|---|---|---|---|
| 0.5 | ~55 | ~13 | ~8 | ~9 |
| 1 | ~76 | ~31 | ~9 | ~10 |
| 2 | ~87 | ~43 | ~11 | ~10 |
| 3 | ~95 | ~59 | ~13 | ~13 |
| 5 | ~101 | ~64 | ~12 | ~13 |
| 12 | ~115 | ~79 | ~19 | ~19 |
| 24 | ~135 | ~92 | ~22 | ~23 |
| 48 | 166 | 118 | 26 | 28 |

**Glyceric acid** [M]: 266 mmol/mol at 37 °C aerated (48 h implied; "266 vs 166" lactic); "significant increase"
under aeration, "significant decrease" at 50 °C; no figure, no deaerated number. Ref 8 [C]: 171 (air) vs 40
(deaerated) at 12 h, 37 °C.

**Fig. 2, γ-pyranone (mmol/mol 1-DG)** — printed [M]: max 173 at 2-3 h (37 °C, deaerated); 12 h: 146 vs 67 (37 °C),
114 vs 103 (50 °C) deaerated vs aerated. Read [F]: 37 °C air 46 / 60 / 72 / 79 / 76 / 68 / 53 / 33 at 0.5 / 1 / 2 /
3 / 5 / 12 / 24 / 48 h; 37 °C deaer. ~76 / 111 / 142 / 173 / 162 / 146 / 136 / 116; 50 °C air ~68 / 98 / 122 / 123 /
121 / 103 / 76 / 43; 50 °C deaer. ~76 / 115 / 146 / 140 / 129 / 115 / 92 / 65.

**Fig. 3, acetaldehyde (mmol/mol 1-DG)** [F] (the caption gives no symbol key; read with Figs. 1-2's convention, full =
aerated, circle = 37 °C [!]): 50 °C air ~14.5 / 20 / 25.5 / 30.3 / 30.7 at 1 / 2 / 3 / 5 / 12 h (text: 30 at 12 h
[M]); 37 °C air ~9.5 / 14.7 / 18.6 / 23 / 26.5; 50 °C deaer. ~8.8 / 12.5 / 15.2 / 18.4 / 19.2; 37 °C deaer. ~5 / 10 /
12 / 14.2 / 17.8. Glycolaldehyde 29 (50 °C, air, 12 h), glyceraldehyde ≤ 13 [M]. β-alanine instead of α-alanine
halves acetaldehyde (so ~half is alanine's Strecker aldehyde) [M].

**Fed γ-pyranone, 12 h, aerated** [M]: glyceric 30 vs lactic 23 mmol/mol γ-pyranone (37 °C); 66 vs 61 (50 °C).

**Table 1, lactic acid carbon origin** [M]: 48 % C1-C3 (COOH = C-3), 0 % reverse; 36 % C3-C5 (COOH = C-5), 0 %
reverse; 8 % C4-C6 (COOH = C-6) and 8 % C4-C6 (COOH = C-4). The natural-abundance correction is not described [!].

## 4. Mechanism evidence and its kind (Glomb laboratory, Halle)

- **Within-study factor (O2 × DTPA, 2 temperatures):** lactic, glyceric, aldehydes up with air; γ-pyranone down. Kind 6,
  confounded by the chelator [!].
- **Isotope labelling (kind 1), but from glucose, not 1-DG:** lactic and glyceric acid origins. Labelled 1-DG was not
  available, so the label routes include every glucose -> C3 route, not only 1-DG's [!].
- **Fed intermediate (kind 2):** γ-pyranone gives the same acids and furan-3-one, O2-dependent; the authors read this as
  rehydration of γ-pyranone back to open-chain 1-DG [Q] (furan-3-one "can in fact only occur via" open-chain 1-DG [V
  paraphrase]).
- **Proposed, not shown [P]:** the tricarbonyl (from 1-DG or its α-oxo-enediol) and its β-cleavage via a ketene
  (Figs. 5-6); 36 % of lactic acid via oxidation of 1-deoxythreosone; retro-aldol of 1-DG "ruled out" as an important
  route (MGO only negligible, quinoxaline data not shown).
- **Self-contradiction [!]:** the Results say all three aldehydes "responded significantly to aeration"; the
  Discussion dismisses a retro-aldol route to acetaldehyde partly for its "independence from oxygen". Fig. 3 shows
  ×1.5-1.6 with air.

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| Triage 7, oxygen as a reagent | **SUPPORTS** (direction, 37-50 °C) | ×4-6 lactic, glyceric up, aldehydes ×1.5; O2 also SUPPRESSES γ-pyranone, i.e. it re-partitions 1-DG, not only adds acids |
| Evidence map "O2: oxidative alpha-cleavage acids — SINGLE (Nestlé)" | **upgrades to two labs for the direction on C3 acids; LEAVES UNCHANGED for acetic acid and for magnitude** | Voigt has no acetic acid and is DTPA-confounded; magnitude at 100-140 °C is NONE |
| Evidence map "O2 -> glucosone CAN" | unchanged | not measured |
| Hub H6 (O2 diverts ARP away from 1-DX, P 0.15) | **LEAVES UNCHANGED** | starts from 1-DG; says what O2 does to 1-DG's fate, not to its formation |
| B7a-2 amine acceleration of 1-DG β-cleavage | **LEAVES UNCHANGED** | every pot has alanine; no amine-free arm |
| CF-069 (acetic acid loss) | unchanged | no acetic acid here (ref 8 has it: 550 air / 650 deaerated at 24 h [C]) |
| Phosphate | unchanged | one buffer |

## 6. Disagreements found

1. **NEW: same laboratory, same year, same nominal pot, amine swapped: yields differ ~2-4×.** Smuda 2010
   (smuda2010_extraction.md: 1-DG + Nα-BOC-lysine 42 mM, 0.1 M phosphate pH 7.4, 37 °C, aerated) gives lactic **~30**
   and glyceric **~79** mmol/mol at 12 h, 58 / 88 at 7 d. Voigt 2010 (Ala) gives lactic **~115** at 12 h and **166** at
   48 h, glyceric **266** at 48 h; ref 8 (Ala) glyceric **171** at 12 h. Acetic acid: Smuda ~250 at 24 h vs ref 8 550
   (air, 24 h). Candidate causes: the amine (α-amino acid vs BOC-protected ε-amine: Strecker-capable, different pKa);
   the nominal 1-DG concentration (protected-precursor moles, deprotection yield unprinted in both); time-zero / batch.
   Same GC-MS silylation method, so the method is an unlikely cause.
2. **With Davidek 2006b (davidek2006b_extraction.md):** Voigt states Davidek found the oxidative route to glyceric /
   lactic acid "negligible", against Voigt's "significant yield". Quantities differ (18O share vs aerated / deaerated
   yield; 90 °C pH 8 glucose + Gly vs 37-50 °C 1-DG + Ala); see sec. 1 point 1 for a reading under which both hold.
   Not a register-grade conflict until both are put on one quantity.
3. Internal: aldehydes "responded significantly to aeration" vs "independence from oxygen" (sec. 4).

## 7. Ledger-ready measurements (two readers)

1. Fig. 1 lactic acid, 4 series × 8 times (text anchors at 2 and 48 h) — the cleanest air / deaerated time course for
   1-DG; key with the DTPA flag.
2. Glyceric acid 266 (text) — a single cell; better paired with ref 8 Fig. 3 (voigt2009) when that is read.
3. Table 1 lactic-acid label distribution (6 cells) — a scale-free route signature (C1-C3 vs C3-C5 vs C4-C6).
4. Fig. 2 γ-pyranone, 4 series — O2 shifts 1-DG from cyclisation to acids.
5. Fig. 3 acetaldehyde, 4 series — only after the symbol assignment is confirmed by a second reader.
