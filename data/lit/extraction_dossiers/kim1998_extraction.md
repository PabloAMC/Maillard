# Kim & Ho 1998 — EXTRACTION (five amino acids heated in medium-chain triglyceride oil, 180 C / 1 h: free ammonia by electrode; 2- and 3-pentylpyridine from each amino acid + 2,4-decadienal by acid-trap GC-FID)
### The oil-medium companion to Kim, Hartman & Ho 1996 (aqueous): the pentylpyridine data are FIGURE-ONLY, but the paper is the source for "ammonia, not the alpha-amino group, is the nitrogen of 2-pentylpyridine" and for "oil gives > 10x more 2-pentylpyridine than water".

**Source on disk:** `data/articles/kim1998.pdf` (owner's download, 2026-09-08). Read from the `pypdf` text
layer (4 pages, words run together in places but every sentence is recoverable). There is **no table** in
the paper. Figure 1 (ammonia per amino acid, bar chart) and Figure 4 (alkylpyridines per amino acid, bar
chart) carry the results: **FIGURE-ONLY**; no value was read off them. The ammonia numbers quoted below
are the ones printed in the Results text. Figures 2 and 3 are drawn deamidation schemes, described from
the text. A layout-mode re-extraction of page 2 confirmed that the amount of 2,4-decadienal charged is
not printed anywhere.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Pentylpyridines in an Oil Medium" |
| Authors | Young-Suk Kim, Chi-Tang Ho (Rutgers, Dept. of Food Science) |
| Venue | J. Agric. Food Chem. 1998, 46(2), 644-647 |
| DOI | 10.1021/jf970719z |
| Companions | Kim, Hartman & Ho 1996, JAFC 44, 3906-3908 (same reaction in water; 3-pentylpyridine "considerable" there); Sohn & Ho 1995, JAFC 43, 3001 (ammonia release from amino acids in water); Zhang & Ho 1989, JAFC 37, 1016 (2,4-decadienal + cysteine / glutathione); Schieberle 1993 (ammonium sulfate raises 2-pentylpyridine in roasted sesame). **Not** the Kim & Ho 1998 J. Food Lipids 5, 173 paper that Zamora 2020 cites for the dienal-imine mechanism. |

## 1. Why it matters

For Programme 7 the paper settles two qualitative points the rule layer needs: (i) the ring nitrogen of
2-pentylpyridine comes from free ammonia, not from the alpha-amino group: glycine, aspartic and
glutamic acid, which release little ammonia, gave "much less" 2-pentylpyridine than asparagine and
glutamine, whose amide side chains deamidate; (ii) the medium matters more than the ammonia supply:
glutamine in oil releases 13x less ammonia than in water yet makes > 10x more 2-pentylpyridine, which
the authors attribute to NH3 (nucleophile) vs NH4+ (not) and to the dienal surviving longer in oil
(no retro-aldol). It also ranks the amide amino acids in oil: glutamine >> asparagine for ammonia
(9x), the reverse of water. What it does not give is any number for 2-pentylpyridine.

## 2. Methods as they matter to a model

- **Ammonia release runs:** amino acid "0.5 mol" (glycine, L-aspartic acid, L-glutamic acid,
  L-asparagine, L-glutamine; Sigma) in 100 mL medium-chain triglycerides (saturated C8 and C10 acids;
  Stepan), 250 mL round-bottom flask, 180 C, 1 h; N2 purge 20 mL/min; ammonia trapped in 100 mL 2 N
  HCl on ice behind a dry-ice/acetone condenser; stored at 4 C. (0.5 mol glutamine is 73 g in 100 mL
  oil, a suspension, not a solution; the charge is as printed, see flag 3.)
- **Ammonia assay:** 5 mL of trap diluted to 100 mL; Orion ISA added to alkaline range; Orion 95-12
  gas-sensing ammonia electrode; NH4Cl standards 10^-1 to 10^-6 M, semilog calibration each run.
- **Pentylpyridine runs:** amino acid 0.005 mol in 100 mL MCT in a 0.3 L Hoke stainless-steel
  cylinder (closed), 180 C oil bath, 1 h. **The 2,4-decadienal charge is not stated.** After cooling,
  2.5 mL of 1000 ppm 2-ethoxy-3-ethylpyrazine (2.5 mg) as internal standard; purge N2 400 mL/min, 6 h,
  sample at 70 C with agitation; trap in 70 mL 11.7 % (w/v) HCl on ice (basic volatiles only); wash the
  acid 3 x 100 mL CH2Cl2; make pH 12.5 with 30 % NaOH; extract 3 x 100 mL CH2Cl2; Kuderna-Danish to
  5 mL, N2 to 0.1 mL.
- **GC-FID:** Varian 3400, DB-1 60 m x 0.25 mm x 0.25 µm; injector 270 C, detector 300 C; He 1.2
  mL/min at 40 C; 40 -> 280 C at 2 C/min, hold 30 min; 1 µL split 50:1; RI vs C6-C19 n-paraffins.
- **GC-MS:** Varian 3400 / Finnigan MAT 8230 magnetic sector, same column; EI 70 eV, source 250 C;
  NIST and Wiley 138.1 libraries plus literature.
- **Quantification:** "peak area relative to that of the internal standard" (FID, no response factor,
  no calibration curve). Results shown as bar charts (Figure 4) only; units on the axis were not read.
- **Replication:** none stated.
- **Medium:** MCT oil, no water added, no buffer, no pH.

## 3. Tables re-typed

There are no tables. The numbers printed in the text:

### Ammonia released at 180 C / 1 h (Results, paragraph 1; Figure 1 holds the bar chart)

| amino acid | medium | ammonia, as printed | note |
|---|---|---|---|
| glutamine | MCT oil | 7.36 x 10^-3 "mol of NH3 / 0.1 M of glutamine" | this study |
| asparagine | MCT oil | 8.01 x 10^-4 "mol of NH3 / 0.1 M of asparagine" | this study |
| glutamine | water | 9.48 x 10^-2 "mol of NH3 / 0.1 M" | "the same reaction in aqueous system" (Sohn & Ho 1995 conditions) |
| asparagine | water | 1.25 x 10^-1 "mol of NH3 / 0.1 M" | idem |
| glycine, aspartic acid, glutamic acid | MCT oil | "a small amount", from deamination of the alpha-amino group | Figure 1 only |

The printed unit "mol of NH3 / 0.1 M of amino acid" is not dimensionally clean; read it as a fraction
of the amino acid charge on the authors' basis. Ratios are safe: oil/water = 0.078 (Gln), 0.0064 (Asn);
Gln/Asn in oil = 9.2; Asn/Gln in water = 1.3.

### Pentylpyridines at 180 C / 1 h in MCT oil (Results, paragraphs 5-8; Figure 4 holds the bar chart)

| statement, as printed | status |
|---|---|
| "all five amino acids generated 2-pentylpyridine but not 3-pentylpyridine (Figure 4)" | FIGURE-ONLY; contradicted two paragraphs later (below) |
| "The relative amount of 2-pentylpyridine produced from asparagine and glutamine ... was proportional to the amount of free ammonia available" | FIGURE-ONLY |
| "much less 2-pentylpyridine was formed from ... glycine, aspartic acid, and glutamic acid" | FIGURE-ONLY |
| "Glutamine produced > 10 times the amount of 2-pentylpyridine in oil systems, with less free ammonia, than in aqueous systems" (aqueous = Kim et al. 1996) | cross-study ratio, no numbers |
| "Only a small amount of 3-pentylpyridine was produced from the asparagine and glutamine, which generated a large amount of 2-pentylpyridine" | FIGURE-ONLY; the "small amount" contradicts the "not 3-pentylpyridine" sentence |
| 3-pentylpyridine "considerable" in water (Kim 1996): "other intermediates degraded from 2,4-decadienal in the presence of water were essential to the formation of 3-pentylpyridine" | interpretation |

### Figures 2 and 3 (drawn deamidation schemes, described)

- Fig. 2, glutamine: the alpha-amino N attacks the side-chain amide C, closing a five-membered ring
  (pyrrolidone carboxylic acid = pyroglutamic acid) with loss of NH3; alternatively the carboxyl O
  attacks the amide giving a six-membered cyclic anhydride + NH3. No water needed.
- Fig. 3, asparagine: the carboxyl O attacks the side-chain amide (five-membered cyclic anhydride) +
  NH3; the amino-N route is disfavoured (would be a four-membered ring).
- No scheme is drawn for 2-pentylpyridine formation. The text's mechanism is one sentence: ammonia
  acts "as a nucleophile to attack the carbonyl center of 2,4-decadienal"; NH4+ cannot.

## 4. Routes and numbers the repository can use

| route | reactant -> product | mechanism as drawn / stated | measured numbers (units, conditions) | evidence class |
|---|---|---|---|---|
| KH-GLN-NH3 | glutamine -> NH3 + pyroglutamic acid (or cyclic anhydride) | Fig. 2 (drawn, described above); intramolecular, anhydrous | 7.36 x 10^-3 (oil) vs 9.48 x 10^-2 (water), authors' unit, 180 C / 1 h | mechanism_drawn; level_only (text) |
| KH-ASN-NH3 | asparagine -> NH3 + cyclic anhydride | Fig. 3 | 8.01 x 10^-4 (oil) vs 1.25 x 10^-1 (water) | mechanism_drawn; level_only (text) |
| KH-DEAMIN | Gly, Asp, Glu -> NH3 (C-N cleavage, deamination) | stated, citing Lien & Nawar 1974 | "small amount", Figure 1 | figure_only |
| KH-2PP | **2,4-decadienal + NH3 (from Asn, Gln) -> 2-pentylpyridine, in oil** | one sentence: NH3 nucleophile on the dienal carbonyl; no scheme | Figure 4 bar chart; ordering Gln, Asn >> Gly, Asp, Glu; oil > 10x water for Gln (vs Kim 1996) | figure_only |
| KH-3PP | 2,4-decadienal + N source -> 3-pentylpyridine | "exact mechanism ... not clear"; needs water-borne dienal fragments | "only a small amount" in oil (Asn, Gln); "considerable" in water (Kim 1996) | figure_only |
| KH-MEDIUM | same reactants, oil vs water | NH3/NH4+ speciation; dienal stability (no retro-aldol in oil, Josephson & Lindsay 1987) | qualitative | within_study statement, no numbers |

Nothing in this paper can be a rate or a yield. Its use is (a) the N-source ranking (amide side chain,
not alpha-amine) as a rule constraint, and (b) the qualitative medium effect as a flag on any rate fitted
in water (Zhou 2000) or on silica (Zamora 2020) when applied to an oil phase.

## 5. Rule sketches (repository suggestions, not the paper's)

Registry: no key for 2-pentylpyridine or 3-pentylpyridine in `data/keys/compounds.yml`; `Gln`, `Asn`,
`Glu`, `Gly`, `Lys` and `DECADIENAL` exist in `data/species/structures.yml`; ammonia does not.

**S1. Deamidation of the amide amino acids to ammonia (net, anhydrous; Fig. 2/3).**
- positive: `NC(CCC(N)=O)C(=O)O` (Gln) -> `O=C1CCC(N1)C(=O)O` (pyroglutamic acid) + `N`; `NC(CC(N)=O)C(=O)O` (Asn) -> `NC1CC(=O)OC1=O` (aminosuccinic anhydride, as Fig. 3 draws) + `N`
- negative: `NC(CCC(=O)O)C(=O)O` (Glu) -> no fire; `NCC(=O)O` (Gly) -> no fire (no side-chain amide). The paper's Figure 1 says these do release a little NH3 by deamination; a separate, slower deamination rule would be `proposed`, not this one.
- Order-of-magnitude constraint from the text: in oil Gln : Asn ≈ 9 : 1; in water Asn ≥ Gln.

**S2. 2,4-decadienal + NH3 -> 2-pentylpyridine.** Same rule as `zamora2020_extraction.md` S1;
this paper adds only the medium flag and the N-source constraint. Positive: `CCCCC/C=C/C=C/C=O` + `N`
-> `CCCCCc1ccccn1`. Negative: `CCCCC/C=C/C=C/C=O` + `NCC(=O)O` (glycine's alpha-amine) must NOT give
2-pentylpyridine directly (the paper: alpha-amino groups "seem to be less involved"); the amine gives a
Schiff base (R02) at most.

**S3. 3-pentylpyridine `CCCCCc1cccnc1`:** no mechanism here; leave unrouted or as `proposed` via the
Zamora acrolein + alkanal route.

## 6. Flags

1. **All pentylpyridine results are in a bar chart (Figure 4) with no numbers in the text**; there is
   no table; quantification is FID area relative to one internal standard without calibration. Evidence
   class figure_only throughout. Do not cite this paper for any 2-pentylpyridine amount.
2. **The 2,4-decadienal charge is not printed** (checked in layout mode on page 2). Even a reader of
   Figure 4 could not compute a yield.
3. **"0.5 mol" amino acid in 100 mL oil** for the ammonia runs (73 g glutamine) against 0.005 mol for
   the pyridine runs: a 100-fold different loading, and the ammonia unit "mol NH3 / 0.1 M" is not
   dimensionally interpretable. Use ratios only.
4. **Internal contradiction on 3-pentylpyridine** ("not 3-pentylpyridine" vs "only a small amount").
5. **No replication, no statistics, no blank** (amino acid without dienal, dienal without amino acid).
6. **The "> 10x more in oil" claim compares two papers** (this one vs Kim 1996, aqueous) with different
   isolation trains; it is a cross-study statement.
7. **Ammonia at 180 C in MCT with N2 purge** measures ammonia that left the oil, not ammonia
   available in it; the paper itself argues that the speciation (NH3 vs NH4+) rather than the total
   controls the pyridine yield.
8. **Zamora 2020's mechanism attribution** ("Kim & Ho, 1998") points to the J. Food Lipids paper, not
   this one; this paper draws no pyridine mechanism.
9. Medium-chain triglyceride (C8/C10 saturated) is inert to oxidation: no lipid-derived carbonyls other
   than the added dienal; the roadmap's pea/soy lipid is polyunsaturated and would add its own.
