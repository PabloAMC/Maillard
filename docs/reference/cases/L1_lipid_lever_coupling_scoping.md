# L1 scoping: do the Maillard levers move the lipid aldehydes? (2026-11-14, literature only)

**Status: scoping note, not a sec. 7 case.** Nothing here changes the engine; the engine is paused (STATUS_AND_PLAN
row 6, STEP5_RECORD). It records what the corpus on disk says, so the case can be written when the pause ends.

## Why

In the P3 engine dry run (`results/validation/p3_engine_dry_run.md`) the lipid lane predicts hexanal, nonanal and
2-pentylfuran at exactly 0 dex for every lever: sugar identity, sulfur source and pH. The zero is structural. The lane
(`src/kinetic_core/lipid.py`) decomposes a pre-existing hydroperoxide pool and has no initiation term. Its only
coupling to the Maillard lanes, the aldehyde-lysine `COVALENT_SINK`, is inert by ruling
(FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2). `lane_coupling_verdict` considers lysine only, and never looks at
cysteine, H2S or pH.

## What the corpus on disk says (ratios: with the lever / without; dex = log10)

| source | system | lever | hexanal | others |
|---|---|---|---|---|
| farmer1992 | 4 lipids +/- cysteine + ribose, pH 5.6, 140 C 1 h | Cys + ribose | -0.08 to -0.49 | nonanal -0.31 to -1.11; 2,4-decadienal -0.23 to -1.30; 2-pentylfuran +0.34 to +0.51 (3 of 4) |
| elmore2002 | methyl linoleate +/- cysteine + ribose, pH 5.5, 140 C 30 min | Cys + ribose | -0.23 (not significant) | 2-enals and dienals 5x to >250x lower; 2-pentylfuran +0.74 |
| mandin1999 | linoleic acid +/- glucose or methionine, pH 6, ~100 C 2 h | glucose; methionine | +0.24; ~0 | decadienal +0.18; -0.6 |
| xin2026 | extruded plant-protein analogue, 8 carbohydrates + control | carbohydrate identity | spread <= 0.38 over all arms | nonanal spread 0.22 |
| madruga1995 | beef, pH 4.0 -> 5.6, 140 C 30 min | pH | ~ -0.14 per pH unit | 2-pentylfuran -0.52 overall |
| desousafontes2024 | goat hydrolysate + xylose + Cys + thiamine, pH 4 vs 6, 121 C 1 h | pH | +0.15 (opposite sign to Madruga) | nonanal +0.25; 2-pentylfuran +0.98 |
| zamora2015b | 2,4-decadienal -> hexanal | pH | yield peaks near pH 8 | |
| zhang1989 | 2,4-decadienal + cysteine, 180 C | cysteine | cysteine also FORMS hexanal (278.6 mg/mol vs 45.4 with glutathione) | |
| zhang2021, du2023, du2024 | aldehyde + cysteine or H2S | thiol sink | thiazolidines, 2-alkylthiophenes (figures only) | |
| fischer2021 | pea isolate, unheated | pH (release) | the isolate carries ~5 ug/g hexanal; release +59 % at pH 4.5 vs 6.5 | |
| elmore1999 | beef, linseed / fish-oil feeding, 140 C 30 min | fatty-acid profile | +0.37 to +0.40 (not tracking linoleic share) | nonanal, decadienal, 2-pentylfuran n.s. |
| morales2001 | 0.1 M sugar + amino acid, 100 C 3 - 24 h | sugar, amino acid (radical scavenging) | not measured | glucose vs lactose <= 0.15 dex; lysine fastest |
| hidalgo2006 | bulk olive oil, Rancimat 110 C | PE, lysine (antioxidant) | not measured | induction period +0.07 to +0.17 dex |
| whitfield1992, zamora2005 | reviews (secondary) | mechanisms | | cysteine vs glycine 2-pentylfuran -0.31 dex (whitfield1988, cited); glutathione diverts decadienal from hexanal |
| alaiz1999 | ribose-browned BSA, 100 ppm in soybean oil, 60 C 264 h | browning pH (buffer-confounded) | not measured | TBARS vs control protein -0.07 / -0.20 / -0.08 dex at pH 4 / 7 / 10 (derived from figures); antioxidant, tracks lysine loss |
| hidalgo1999 | BSA + oxidised-linoleate products, 80 C | pH | not measured | lysine consumed +0.33 dex at pH 7 vs 4: the amine sink peaks near neutral |
| ohnishi1984 | dry beef fat +/- glycine, 200 C 4 h (area %, single runs) | amino acid (NH3 sink) | -0.38 | nonanal -1.14; heptanal, decanal <= -1.9; 2-alkylpyridines form |
| mottram1989 | cysteine + ribose +/- lipid, pH 5.7, 140 C 1 h (review, relative values) | lipid on the sulfur lane | not measured | lipid lowers H2S-derived thiols -0.14 to -1.52 dex; cysteine >= +1.0 dex 2-pentylpyridine vs glycine |
| farag1982 | MRPs of 16 amino acids x 5 sugars on linoleate emulsion, 50 C | amino acid; sugar (induction period) | not measured | cysteine-HCl products the weakest antioxidant (-0.05 to -0.14 dex vs median; pH 2-confounded); glucose vs pentose median +0.03 dex |

**Added 2026-11-15** (the last five rows; dossiers in `data/lit/extraction_dossiers/`). None supports a declared band:
- Maillard products act both ways on lipid oxidation: some studies find them pro-oxidant (glycated phospholipid,
  Amadori compounds), others antioxidant (browned protein). No heated-pot measurement on disk settles which.
- pH is almost unstudied (Whitfield names it a gap).
- No paper on disk tests thiamine against cysteine in a lipid pot.
- elmore1999 and whitfield1992's numbers come from the Reading laboratory and pot design of the two hold-outs, so a
  leave-one-laboratory-out split puts them in the same fold.

**Added 2026-11-15, second batch** (alaiz1999, hidalgo1999, ohnishi1984, mottram1989, farag1982). No probability in the Reading changes; none measures hexanal with and without a P3 lever:
- Sulfur: ohnishi1984 and mottram1989 support an amine / H2S carbonyl sink, strongest on dienals and long aldehydes. farag1982 adds a small opposite channel (cysteine products the weakest antioxidant, <= 0.14 dex and pH-confounded), so the sign is a little less certain; the magnitude is unchanged.
- pH: alaiz1999 and hidalgo1999 lean weakly to hexanal lower at pH 7 than pH 4 (more antioxidant browning, a stronger amine sink); no pH 6 arm, and pH is confounded with buffer.
- Sugar: farag1982's glucose-vs-pentose median is about 0 with +/- 0.1 dex scatter.
- Folds: alaiz1999 and hidalgo1999 join zamora2005 and hidalgo2006 in the Seville laboratory; mottram1989 joins the Reading fold of the hold-outs.

## Reading (calibrated; the probability that the lever moves hexanal by more than 0.3 dex in the P3 pot)

- **Sulfur source (cysteine vs thiamine): P ~ 0.30, cysteine lower.** Thiazolidines and H2S trapping are specific to
  the thiol. Nonanal P ~ 0.45, cysteine lower; 2-pentylfuran P ~ 0.3, sign unclear (up with cysteine + ribose in farmer1992 and elmore2002, down for cysteine vs glycine in whitfield1988). A zero is not
  defensible here.
- **Sugar (glucose vs ribose): P ~ 0.15, sign unknown.** Xin 2026 bounds the spread across carbohydrates at 0.38 dex.
- **pH (6 vs 7): P ~ 0.12, sign unknown.** The two on-disk pH series disagree in sign.

Two things in this pot lower all three probabilities. The isolate's carried-over hexanal dominates, and the protein's
amines already make up most of the carbonyl sink, so swapping one additive changes less.

## What the case would be, when the pause ends

- A **thiol / H2S carbonyl sink**, carbonyl-class specific: saturated aldehydes ~2x, enals and dienals 5x - 80x,
  ketones and 1-octen-3-ol untouched (Farmer 1992's pattern). Its structure comes from elmore1997, zhang1989 and
  du2023/2024, with a rival in which cysteine also forms hexanal (zhang1989).
- Until it is built, the shipped zeros should carry declared bands, not point zeros:
  - sulfur source on hexanal [-0.5, 0] dex;
  - sugar and pH +/- 0.3 dex.

  This is a presentation change, made only when the pause ends.
- `lane_coupling_verdict` should name cysteine as a shared species alongside lysine.
- **The hold-out conflict comes first.** farmer1992 and elmore2002 are designated hold-outs in their dossiers, so a
  band set from them would be fitted to data read for scoring. Either the band comes from the mechanism papers alone,
  with those two kept as the test, or they are formally re-designated before use.
- P3 itself measures hexanal, nonanal and 2-pentylfuran on all three levers, so it is the designed test of this.
