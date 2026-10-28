# Literature search 2026-10: mechanism-first search for the sulfur lane's missing routes

Run 2026-10-28. Earlier rounds searched by rate-data keyword ([LITERATURE_SEARCH_2026-10.md](LITERATURE_SEARCH_2026-10.md),
[ROUND3](LITERATURE_SEARCH_2026-10_ROUND3.md), [VARIABLES](LITERATURE_SEARCH_2026-10_VARIABLES.md),
[PATHS](LITERATURE_SEARCH_2026-10_PATHS.md)) and by forward citation ([FORWARD](LITERATURE_SEARCH_2026-10_FORWARD.md)). This
round starts from the **routes** the model lacks. It reads the engine's own hypothesis layer
(`results/validation/network_hypotheses.json`), then searches eight mechanism items (M1-M8) against the errors in
[SULFUR_FLUX_MAP.md](SULFUR_FLUX_MAP.md) sections 4-6. It records the search only. Nothing was added to `data/`, and no paper
was downloaded.

## Bottom line

- **The modern literature still has no aqueous rate law for any of the eight routes at 90-150 °C with cysteine or
  glutathione.** Most of what the searches return is identification work, sensory work or applications. The numbers that do
  exist sit in three neighbouring literatures the food-flavour rounds had not touched:
  1. **Microwave and aseptic sterilisation "chemical markers" (Kim & Taub; Tang).** Norfuraneol, called "M-2" there, forms
     from ribose plus protein amines. These papers report **first-order formation rate constants with Arrhenius parameters
     at about 100-130 °C** (Lau 2003; Pandit 2006; Kim 1996). They give the only measured activation energy for the
     DPO → norfuraneol branch. That bears on finding 2 (the trunk's shared 64 kJ/mol is too low) and on the FFT/MFT
     temperature mis-weighting (finding 3).
  2. **Physical-organic thiol chemistry (Whitesides; Jencks).** These are rate and equilibrium constants for thiolate +
     disulfide exchange (with a Brønsted law in thiol pKa) and for thiol + aldehyde hemithioacetals. They would put numbers
     on two of the hypothesis layer's unrun steps: the R12 mixed disulfides, which are the suspected MFT sink in Zhou's
     glutathione pot, and the R14 hemithioacetals. The numbers are at 25 °C.
  3. **Environmental and dialysis-fluid chemistry.** Chen & Morris 1972 gives aqueous sulfide autoxidation kinetics.
     Kjellstrand 2001 gives the pH profile of glucose degradation at 121 °C without amine.
- **One modern-looking aqueous kinetic paper on cysteine + sugar:** Schirle-Keller & Reineccius 1992 (glucose + cysteine in
  water, 80/100/120/150 °C). Furfural, **furfuryl alcohol**, HMF and other furans form at zero order with **Ea 28-33 kcal/mol
  (117-138 kJ/mol)** [abstract]. That is a direct, multi-temperature check on the trunk Ea in a cysteine pot, but with a
  hexose. **Correction after reading the paper (2026-10-28, `schirle1992_extraction.md`):** the abstract's range
  contradicts the paper's own Table 2, which prints furfural at 35.18 kcal/mol = **147 kJ/mol** (a refit of Table 1
  gives 156-189). The charge is 0.5 M glucose + 0.1 M cysteine (0.2 and 0.04 mol in 400 mL), not "0.5 + 0.1 mol".
  No furfuryl-alcohol Ea is printed. B62 uses the Table 2 value.
- **Furfuryl alcohol → FFT (M1)** has route evidence beyond Cerny 2021 and Liu 2021, which the corpus holds: Rizzi 1995 and
  2000 (furfuryl and allylic alcohols plus sulfur nucleophiles at acid pH), and Parliment & Stahl 1995 (FFT in coffee model
  systems). **No rate constant was found** for FA + H2S or FA + Cys → FFT in water. I am about 80 % confident none is
  published.
- **Plain negatives** (searched, nothing beyond the corpus or earlier documents):
  - M2: MFT oxidative dimerisation kinetics in water at any temperature.
  - M4: the mercaptoacetaldehyde → H2S + acetaldehyde split, and H2S yield per dicarbonyl.
  - M5: cysteamine formation from cysteine.
  - M6: glutathione degradation kinetics with sugars beyond the Xia/Zhang group papers already listed.
  - M8: H2S + furfural or dicarbonyl rate or equilibrium constants in water.
- **Calibration.** I put about 70 % on Lau 2003 and Pandit 2006 each reporting a usable Ea for norfuraneol formation (both
  abstracts promise Arrhenius parameters but print no number). I put about 35 % on that Ea transferring to a free-ribose +
  cysteine pot within ±20 kJ/mol, because a protein gel at pH ~6-7 is a different amine. I put about 15 % on having missed
  an aqueous FA → FFT or MFT → disulfide rate law that exists (Chinese-language journals are the main blind spot).

## Method

- **Engines.**
  - **Europe PMC REST**, restricted to `TITLE_ABS:` with boolean phrasing. This was the main engine. It is unbudgeted and
    returns abstracts.
  - **OpenAlex** `title_and_abstract.search`, about 30 queries, plus `referenced_works` and `cites:` for one-hop chaining.
    OpenAlex answered this round. It reports a shared $0.10/day budget per IP at about $0.001 per search, so it was used
    sparingly.
  - **Semantic Scholar** search mostly returned HTTP 429, and its keyword ranking was noisy (a "furfuryl alcohol H2S" query
    returned catalysis papers). It contributed nothing unique.
  - **WebSearch**, three queries, for full-text snippets of chapters without abstracts.
  - **Crossref** `/works/<doi>` verified every DOI below.
  - Everything was anonymous: no key, no e-mail. The responses are cached under `scratch/mech_search/cache/` (gitignored).
- **Exclusion set.** 921 DOIs: every DOI in `data/keys/papers.yml`, in the five `LITERATURE_SEARCH_2026-10*.md` documents and
  in `CITATION_CONTEXT_2026-10.md`, plus the `data/articles/` file names.
  - Excluded hits are not re-reported. Among them were Cerny 2021, Cerny 2013 (S-conjugates), Hofmann 1996 (thiol
    stability), Mottram 1995 (disulfide stability), Grosch & Zeiler-Hilgart 1992, Wang 2020 (2-acetylthiazole pathways),
    Liu 2023b (carbon-module labelling, `liu2023b`) and the Tang/Ho 2012 MFT review.
- **Queries.** About 45 Europe PMC queries, 30 OpenAlex queries and 12 Semantic Scholar queries across M1-M8, 4-8 phrasings
  per item.
- **One-hop chains** (references plus up to 100 citers): Schirle-Keller 1992, Albouchi 2018, Cerny 2021, Cerny 2013,
  Hofmann 1996, Lau 2003, Pandit 2006, Kobayasi 1965, Huang 1998, Larralde 1995, Kaufmann 2018, Rizzi 2010 and the GSH-xylose
  Amadori paper (10.1021/acs.jafc.9b04694).
- **Triage.** About 600 titles were screened and about 110 abstracts read.
  - **Classes.** **Q** = measured numbers usable as a parameter or hold-out. **M** = mechanism or labelling evidence for a
    route. **B** = identification or application; dropped.
  - **Tags.** [abstract] = read in the abstract. [snippet] = full-text search excerpt. [inferred] = from the title or venue.

---

## Step 1: the hypothesis layer against the sulfur lane's errors

`network_hypotheses.json` has 38 rules and 270 steps: 26 modelled, 172 "mechanism_known" and 72 "proposed". The table lists
the sulfur-lane steps that bear on the current misses (flux map sections 4-5b).

| step (rule) | placement | error it bears on | does literature measure it? |
|---|---|---|---|
| Cys + PENT → Cys-Amadori; Cys-Amadori → DPO + Cys (R04, 2,3-enolisation) and → TDP + Cys (R03) | mechanism_known | Finding 2 (trunk 10³-10⁴x too fast); the pH 3 arm (MFT 22x low) | Only indirectly: Cao 2017 (TTCA K ≈ 60 L/mol), Zhai 2019-23 and Kang 2026, all in the corpus. **No Cys-ARP enolisation rate constant anywhere.** New this round: the M-2 (norfuraneol) Arrhenius papers measure the lumped amine-catalysed pentose → NF rate (Lau 2003, Pandit 2006, Kim 1996) |
| Cys + PENT → TTCA (R15) | mechanism_known in this artifact (B58 tested it) | Finding 2 (the reservoir) | Cao 2017 equilibrium. Rates: de Roos 1992 (round 3). Nothing new. Dworkin & Miller 2000 (open-chain aldehyde fraction of aldoses) sets the free-carbonyl term of any Schiff-base/TTCA rate law |
| Cys + MGO / DA → aminoketone + mercaptoacetaldehyde (R07, Strecker) | mechanism_known; the engine's `r_dpo_ddp` sends the S to an inert pool | pH 7 arm (H2S starvation, 4 900x low) | **No.** The split mercaptoacetaldehyde → H2S + acetaldehyde vs ring products is unmeasured. Kobayasi 1965 (corpus) gives total H2S with glucosamine as the reductant (4.2-8 % of cysteine), not per dicarbonyl |
| MFT + Cys, FFT + Cys, MFT + FFT → mixed disulfides (R12) | mechanism_known | Zhou GSH pot: MFT collapse, ~81 % of the cysteine moiety unaccounted for; Zhang dimer 3e-5x | Identification only (Mottram 2002, corpus). **New:** Szajewski & Whitesides 1980 and Lees & Whitesides 1993 give k and K for RS⁻ + RSSR exchange as a Brønsted law in thiol pKa, so MFT/FFT + GSSG/cystine can be estimated, but **no measured pKa for MFT or FFT was found** |
| Cys + FUR → furyl-thiazolidine; Cys + MGO / GO → thiazolidines (R15) | mechanism_known | FFT/MFT ratio 6-60x high; Zhou furfural 25-45x low | Mizutani 1960 and Yaghmur 2002/05 (corpus). Huang, Su & Ho 1998 (new, M): phosphate strongly promotes thiazolidine formation and a phosphate-mediated thiazolidine → thiazoline redox |
| MFT / FFT / Cys + FUR, MGO, GO, PENT → hemithioacetals (R14) | proposed (no source) | a reversible MFT/FFT sequestration term in carbonyl-rich pots | **New:** Lienhard & Jencks 1966 measures equilibria and rates of thiol addition to aldehydes. These are class constants, not MFT-specific |
| NF + H2S → 4-mercapto-5-methyl-3(2H)-furanone (R36) | mechanism_known | competes with MFT for NF and H2S (MFT low at pH 3 and 7) | No rate anywhere. Whitfield 1999/2001 identify mercaptoketones and thiols from NF + H2S/Cys (corpus) |
| NOT IN THE RULE SET: furfuryl alcohol and FA → FFT; cysteamine → 2-acetyl-2-thiazoline from cysteine; glutathione as a species (GSSG, Cys-Gly); phosphate catalysis; H2S oxidation | — | M1, M5, M6, M7, M8 | see the per-item sections |

**The steps that matter most**, in order of how much of the scorecard they could move:

1. The R07 Strecker sulfur split (the pH 7 arm).
2. Cys-Amadori 2,3-enolisation (R04), together with the pentose enolisation pH law (finding 2 and the pH 3 arm).
3. The R12 MFT / GSH mixed-disulfide exchange (the Zhou series).
4. R15 Cys + furfural thiazolidine as a furfural sink (FFT/MFT ratio).

The first two have **no direct measurement in any literature found** across six rounds. The third can now be estimated from
physical-organic constants. The fourth has only Yaghmur's 40-70 °C data.

---

## Per-item findings

### M1: furfuryl alcohol formation, and FA → FFT

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Schirle-Keller, J.-P. & Reineccius, G. A. (1992).** Reaction kinetics for the formation of oxygen-containing heterocyclic compounds in model systems. ACS Symp. Ser. 490:244-258 | [10.1021/bk-1992-0490.ch019](https://doi.org/10.1021/bk-1992-0490.ch019), paywalled chapter | [abstract] glucose 0.5 mol + cysteine 0.1 mol in 400 mL water, **80, 100, 120 and 150 °C**, 2-10 h, GC-AED. Furfural, furfuryl alcohol, 2-acetylfuran, 5-methylfurfural, DDMP and HMF form at **zero order**, with **Ea 28-33 kcal/mol** | **Q.** (i) Formation rates of FA and furfural in a cysteine + sugar pot at four temperatures. (ii) An Ea of 117-138 kJ/mol against the trunk's 64 kJ/mol (finding 2). Companion to `chan1994` (same lab, sulfur volatiles). Hexose, unbuffered |
| **Parliment, T. H. & Stahl, H. D. (1995).** Formation of furfuryl mercaptan in coffee model systems. Dev. Food Sci. 37:805-813 | [10.1016/s0167-4501(06)80197-4](https://doi.org/10.1016/s0167-4501(06)80197-4), paywalled | [inferred] no abstract online. Cited by Cerny 2021 for FFT precursors. Likely compares furfural, FA, cysteine and H2S as FFT sources in aqueous models | M, possibly Q (0.4 for yields per precursor) |
| **Rizzi, G. P. (1995).** Formation of sulfur-containing flavor compounds from allylic alcohol precursors. Dev. Food Sci. 37:289-302 | [10.1016/s0167-4501(06)80162-7](https://doi.org/10.1016/s0167-4501(06)80162-7), paywalled | [inferred] allylic and furfurylic alcohols + H2S or thiols → thiols by SN1 substitution (an acid-catalysed carbocation route) | M (route; acid dependence) |
| **Rizzi, G. P. (2000).** Formation of sulfur-containing volatiles under coffee roasting conditions. ACS Symp. Ser. 754:210-215 | [10.1021/bk-2000-0754.ch022](https://doi.org/10.1021/bk-2000-0754.ch022), paywalled | [abstract] "furfuryl mercaptan and related compounds are readily formed by reactions of Maillard-derived furfuryl alcohols and amides of S-containing amino acids at acid pH" | M. The FA → FFT route via a cysteinyl intermediate at acid pH (cf. Cerny 2021's S-furfuryl-L-cysteine) |
| Albouchi, A. & Murkovic, M. (2018). Formation kinetics of furfuryl alcohol in a coffee model system. Food Chem. 243:91-95 | [10.1016/j.foodchem.2017.09.112](https://doi.org/10.1016/j.foodchem.2017.09.112), paywalled | [abstract] FA formation during roasting of green coffee, temperature and time effects, up to 57 % lost by evaporation | Q but dry roasting. Low transfer. SKIP unless FA enters the model |
| Swasti & Murkovic 2012, Food Funct. ([10.1039/c2fo30020f](https://doi.org/10.1039/c2fo30020f)) | — | [abstract] FA "is formed via the degradation of quinic acid or 1,2-enediols", and polymerises in acid | M. Note: FA need not come from furfural reduction |

**Verdict.**

- **FA formation.** Schirle-Keller is the only aqueous, multi-temperature FA formation rate found in a cysteine pot. The route
  furfural → FA by reductones or Cannizzaro in Maillard water was **not found measured anywhere**. The biomass literature is
  all heterogeneous-catalysis hydrogenation.
- **FA + S → FFT.** Route evidence exists: Rizzi 1995/2000, Parliment 1995, and Cerny 2021 and Liu 2021 in the corpus. No
  rate exists.
- **For the model.** An FA branch would enter as a structural hypothesis, with Hofmann's FFT/MFT pH series and Yiltirak's
  temperature ladder as pre-registered hold-outs. Its strength would be declared, not fitted from these papers.
- **Ties to the error.** If FFT comes partly via FA by an acid-catalysed SN1 route, FFT formation would carry an acid
  dependence and a different Ea from furfural + H2S. That is one candidate for finding 3's FFT/MFT mis-weighting. Poisson
  2009 already argues against free furfural as the FFT route in coffee.

### M2: MFT/FFT oxidative coupling, mixed disulfides, thiol sinks

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Szajewski, R. P. & Whitesides, G. M. (1980).** Rate constants and equilibrium constants for thiol-disulfide interchange reactions involving oxidized glutathione. JACS 102:2011-2026 | [10.1021/ja00526a042](https://doi.org/10.1021/ja00526a042), paywalled | [known content, not re-read] second-order k for RS⁻ + GSSG across thiols of pKa about 7-10, at 25 °C in water. Brønsted β_nuc ≈ 0.5 | **Q.** Parameterises R12 exchange (MFT + GSSG ⇌ MFT-SG + GSH) for the Zhou pot from the MFT thiolate fraction. Needs an MFT pKa (none measured; heteroaryl thiols are about 6-7, an assumption) and a T-extrapolation assumption |
| Lees, W. J. & Whitesides, G. M. (1993). Equilibrium constants for thiol-disulfide interchange reactions: a coherent, corrected set. J. Org. Chem. 58:642-647 | [10.1021/jo00055a016](https://doi.org/10.1021/jo00055a016), paywalled | [inferred] corrected K values for GSH/cysteine/disulfide pairs | Q (equilibria for the GSH ⇌ GSSG ⇌ Cys mixed-disulfide pool) |
| Sun, Z. et al. (2020). Enhancement of coffee brew aroma through control of the aroma staling pathway of 2-furfurylthiol. Food Chem. 322:126754 | [10.1016/j.foodchem.2020.126754](https://doi.org/10.1016/j.foodchem.2020.126754), OA status per Semantic Scholar | [abstract] hydroxyhydroquinone-derived quinones bind FFT. Cysteine releases it. Storage temperature | M (thiol-quinone sink; cysteine as a competitive releaser). B for numbers |
| Franco-Luesma, E. & Ferreira, V. (2014). Quantitative analysis of free and bonded forms of volatile sulfur compounds in wine … reversible cation-complexed forms. J. Chromatogr. A 1359:8-15 | [10.1016/j.chroma.2014.07.011](https://doi.org/10.1016/j.chroma.2014.07.011), paywalled | [inferred] H2S and thiols reversibly bound by Cu/Zn | M (a metal-complex sink; relevant only if pot metals are recorded) |
| Zhang, Kuang, Wang & Cao 2024, Food Res. Int. ([10.1016/j.foodres.2024.114149](https://doi.org/10.1016/j.foodres.2024.114149)) | already in an earlier document | Cys/GSH 3:7 "controls the free-radical reaction of sulfides". Metabolomics | B |

**Verdict.**

- **Measured MFT oxidation kinetics in water: none found.** The corpus already holds every measured MFT/FFT stability datum
  in the open literature: Hofmann 1996 in ether, Blank 2002 under Fenton conditions, Mottram 1995 for disulfide stability,
  Hofmann 2001/02 and Müller 2005-07 in coffee, and Charles-Bernard 2005. I am about 85 % confident of this.
- **Oxidant identity.** No paper found says which oxidant (O2, dicarbonyl or metal) couples MFT in a sealed 120-145 °C
  buffer.
- **The Zhou sink.** The physically grounded way to parameterise it is thiol-disulfide exchange with GSSG, using the
  Whitesides constants. That makes MFT-SG the sink: reversible, and fast at pH 7.5. A test of that structure: the model
  should reproduce the MFT plateau-then-collapse once GSSG accumulates.

### M3: norfuraneol, furanones and thiols

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Lau, M. H., Tang, J., Taub, I. A., Yang, T. C. S., Edwards, C. G. & Mao, R. (2003).** Kinetics of chemical marker formation in whey protein gels for studying microwave sterilization. J. Food Eng. 60:397-405 | [10.1016/s0260-8774(03)00063-3](https://doi.org/10.1016/s0260-8774(03)00063-3), paywalled (a WSU repository record exists) | [abstract] M-2 (= norfuraneol) "formed from D-ribose and amines through … enolization under low acid conditions (pH > 5)". **First-order formation; Arrhenius k(T); Ea by three regression methods**. 20 % whey protein gel | **Q.** The only measured Ea for the amine-catalysed pentose → NF branch at sterilisation temperatures. Sets a band on the NF branch's Ea (`k_dpo_nf` currently rides the shared trunk Ea). Tests finding 3's "the furfural branch's Ea is wrong relative to the NF branch" |
| **Pandit, R. B., Tang, J., Mikhaylenko, G. & Liu, F. (2006).** Kinetics of chemical marker M-2 formation in mashed potato. J. Food Eng. 76:353-361 | [10.1016/j.jfoodeng.2005.05.032](https://doi.org/10.1016/j.jfoodeng.2005.05.032), paywalled | [abstract] 1.5 % D-ribose in mashed potato (83 % water), four temperatures, first-order, Arrhenius; "amino acid was the limiting element"; salt has no effect | **Q.** An independent second Ea in a different matrix |
| **Kim, H.-J., Taub, I. A., Choi, Y.-M. & Prakash, A. (1996).** Principles and applications of chemical markers of sterility in HTST processing of particulate foods. ACS Symp. Ser. 631:54-69 | [10.1021/bk-1996-0631.ch006](https://doi.org/10.1021/bk-1996-0631.ch006), paywalled | [abstract] M-2 forms from D-ribose **or ribose 5-phosphate** with amines via 2,3-enolisation under weakly acidic or neutral conditions. M-1 (DDMP) from glucose. The 1998 sequel (10.1007/978-1-4899-1925-0_9) uses "the rate constants and the activation energies" of the markers | Q/M. The precursor chemistry and probably the original k and Ea, with the pH window. M-1 also bears on the trunk lane |
| Gupta, R. et al. (2011). Combined pressure-temperature effects on M-2 formation in whey protein gels. LWT 44:2141-2146 | [10.1016/j.lwt.2011.05.007](https://doi.org/10.1016/j.lwt.2011.05.007) | [abstract] 105 °C, pH 6.1 vs 8.3: first-order; **M-2 yield decreases with increasing pH** | Q (a shape point: NF yield pH 6.1 > 8.3 at 105 °C; compare Mikami 2017) |
| Hauck, T. et al. (2003). Alternative pathway for the formation of 4,5-dihydroxy-2,3-pentanedione … BBA 1623:109-119 | [10.1016/j.bbagen.2003.08.002](https://doi.org/10.1016/j.bbagen.2003.08.002), paywalled | [abstract] ribulose 5-phosphate → 1-deoxypentosone (DPD) → norfuraneol **non-enzymatically**, by phosphate elimination, 1.3 % in 15 h at 35 °C and pH 7.5. 13C labelling | M/Q. A phosphate-leaving-group NF route (IMP and ribose-5-P matrices). Not for free-ribose pots |

**Verdict.**

- **NF + cysteine/H2S.** No new quantitative paper beyond the corpus (Whitfield 1999/2001, Hofmann 1998, Zheng 1997,
  Shu 1985/88, Nakamura 2020/21). Thiol addition to 3(2H)-furanones has no measured rate.
- **The new handle is the M-2 literature.** It measures how fast norfuraneol forms from ribose + amines vs temperature,
  which is exactly the quantity that splits the NF/MFT branch's temperature dependence from the furfural/FFT branch's.
- **Caveats.** Protein amines rather than cysteine; gels at about pH 6-7; and the exact temperature ranges are not in the
  abstracts.

### M4: Strecker degradation of cysteine; the fate of mercaptoacetaldehyde

- **Searched:** mercaptoacetaldehyde (Europe PMC and OpenAlex return almost only synthetic [3+3] annulation chemistry);
  1,4-dithiane-2,5-diol; dithiane and trithiolane with cysteine; "Strecker degradation" AND cysteine AND H2S; cysteine +
  MGO/GO/diacetyl + H2S yield; β-elimination; and a chain from Kobayasi 1965.
- **Found nothing quantitative beyond the corpus and earlier documents.** That covers Kobayasi 1965, Edwards & Wedzicha
  1992, Wedzicha 1991, Pripis-Nicolau 2000, Du 2023 and Liu 2023b's carbon-module labelling (3-thiophenethiol carbon mostly
  from cysteine via Strecker).
- **Correction to PATHS gap 1.** The 0.42/0.80 mmol H2S numbers in Kobayasi 1965 are from cysteine + **glucosamine**, the
  enaminol reductant test, not from an α-dicarbonyl. The corpus dossier `kobayasi1965_extraction.md` already records this.
- **New M-class items:**
  - Li et al. 2022, Food Chem. Mol. Sci. 5:100120 ([10.1016/j.fochms.2022.100120](https://doi.org/10.1016/j.fochms.2022.100120),
    **open**). In heated whey protein, H2S comes from cysteine residues by **β-elimination**, with no free thiol amino acid
    required. That is an H2S source independent of Strecker, for protein matrices.
  - Rizzi 1992, "Formation of dialkylthiophenes in Maillard reactions involving cysteine" (Dev. Food Sci., no DOI found):
    a ring-product sink for the Strecker fragments [inferred].
- **I am about 85 % confident that no measurement of the H2S : ring-product split of mercaptoacetaldehyde exists.** The
  split stays a declared branching fraction, with Hofmann's pH sweep as the shape hold-out.

### M5: cysteamine, 2-acetyl-2-thiazoline, 2-acetylthiazole

| paper | DOI / access | what it has | class |
|---|---|---|---|
| Huang, T.-C., Su, Y.-M. & Ho, C.-T. (1998). Mechanistic studies on thiazolidine formation in aldehyde/cysteamine model systems. JAFC 46:224-227 | [10.1021/jf9705633](https://doi.org/10.1021/jf9705633), paywalled | [abstract] "Buffer dramatically promotes thiazolidine formation". **Phosphate**-catalysed thiazolidine → thiazoline redox (hydride transfer) | M (phosphate catalysis of the thiazoline step; overlaps M7) |
| Huang, Su & Ho 1998b, cysteamine/2,3-butanedione, JAFC 46:664-667 | [10.1021/jf970602a](https://doi.org/10.1021/jf970602a) | [abstract] phosphate strongly enhances 2-methyl-2-acetylthiazolidine | M |
| Novotny, O., Dufossé, T. & Davidek, T. (2019). Unexpected potential of iso-oligosaccharides in the generation of important food odorants. JAFC 67:2954-2962 | [10.1021/acs.jafc.8b07100](https://doi.org/10.1021/acs.jafc.8b07100), paywalled | [abstract] yields of 2-acetylthiazole, 2-acetyl-2-thiazoline and HDMF from sugar + Cys/Gly/Pro binary mixtures (probably SIDA; Nestlé) | Q? (0.4 that it is aqueous with per-mol yields) |

**Verdict.**

- **Cysteamine formation from cysteine under Maillard conditions has no measurement in any source found.** I am about 85 %
  confident of this. The corpus's Hofmann & Schieberle 1995b remains the only kinetic data on the cysteamine route.
- **2-Acetylthiazole.** Its pathway papers (Wang 2020 Food Chem; the 2022-23 JAFC papers on cysteine vs deoxyosone
  fragments) are already in earlier documents.
- **For the model.** Zhang 2025b's 2-acetylthiazole (1 750x low) can only enter as a declared route:
  Cys + MGO → (Strecker) aminoketone/H2S/NH3 → thiazole. It has no measured constant. Novotny 2019 is the best chance of a
  per-mol yield.

### M6: glutathione Maillard chemistry

- **Everything relevant is already in the corpus or earlier documents.** That includes:
  - the Xia/Zhang Jiangnan group's GSH-Amadori work (N-(1-deoxy-D-xylulos-1-yl)-GSH, 2019-24);
  - Zhang & Ho 1988/89/91 and Zheng & Ho 1994;
  - Wang 2012, Zhou 2023c and Du 2025 (GSH-Amadori + fatty acids);
  - the GSH preformulation degradation-kinetics paper (10.1080/03639045.2020.1752708).
- **The chain from the GSH-xylose Amadori paper adds only applications.**
- **What GSH-as-a-species needs beyond those:**
  - the exchange and equilibrium constants of M2 (Szajewski 1980; Lees 1993);
  - a GSSG-formation rate under the pot's oxidant. No aqueous GSH autoxidation constant at above 90 °C was found. The
    biochemical metal-catalysed GSH autoxidation literature is at 25-37 °C and was not pursued.
- **Class: M6 is a negative for new sources.**

### M7: buffer (phosphate) catalysis of sugar enolisation and Amadori degradation; amine catalysis

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Kjellstrand, P. et al. (2001).** Degradation in peritoneal dialysis fluids may be avoided by using low pH and high glucose concentration. Perit. Dial. Int. 21:338-344 | [10.1177/089686080102100402](https://doi.org/10.1177/089686080102100402), paywalled | [abstract] glucose 1-60 %, **pH 1-8, 121 °C**, no amine: **degradation minimum at initial pH ≈ 3.5**; the HMF route dominates below pH 3.5 and 3-DG-type products above | Q (shape). Independent support, in a hexose, for PATHS's inference that pH 3 sits at the sugar pH-rate minimum, where the specific-acid term should not dominate |
| **Rizzi, G. P., Amba, E. E. & Heineman, W. R. (2010).** Quantification of chemically reducing species in the phosphate ion catalyzed degradation of reducing sugars. JAFC 58:9739-9743 | [10.1021/jf102433u](https://doi.org/10.1021/jf102433u), paywalled | [abstract] reductones titrated (DCPIP) and redox potential in **0.1 M phosphate at 100 °C**, time-resolved | Q/M. Phosphate catalysis makes reductones without amine. These are the reductant for furfural → FA (M1) and a redox buffer against thiol oxidation (M2) |
| Kaufmann, M., Krüger, S., Mügge, C. & Kroh, L. W. (2018). General acid/base catalysis of sugar anomerization. Food Chem. 265:216-221 | [10.1016/j.foodchem.2018.05.101](https://doi.org/10.1016/j.foodchem.2018.05.101), paywalled | [abstract] a rate law for general acid/base catalysis of ring opening (fructose). "The most effective catalysts … are weak acids and their conjugated bases". Sugar reactivity in amino-acid systems is "primarily due to carboxylic acid catalysis" | M/Q (a general-catalysis rate-law form; fructose only) |
| Zhang, W. & Serianni, A. S. (2012). Phosphate-catalyzed degradation of D-glucosone … JACS 134:11511-11524 | [10.1021/ja3020296](https://doi.org/10.1021/ja3020296), paywalled | [abstract] Pi (and arsenate) catalyse osone degradation at pH 7.5, 37 °C; mechanism through phosphate adducts; also shown for D-xylosone | M (phosphate is a nucleophilic catalyst, not just general base) |
| Dworkin, J. P. & Miller, S. L. (2000). A kinetic estimate of the free aldehyde content of aldoses. Carbohydr. Res. 329:359-365 | [10.1016/s0008-6215(00)00204-4](https://doi.org/10.1016/s0008-6215(00)00204-4), paywalled | [inferred] open-chain fractions of the pentoses and hexoses (Larralde's lab) | Q. The free-carbonyl term of any Schiff/TTCA rate law (finding 2); ribose vs xylose |
| Borenfreund, E. & Dische, Z. (1957). Instability of ketopentose-5-phosphates in buffer solutions of varying pH. BBA 25:215-216 | [10.1016/0006-3002(57)90454-7](https://doi.org/10.1016/0006-3002(57)90454-7) | [inferred] decomposition vs pH of ribulose/xylulose 5-phosphate (the Hauck 2003 NF route) | M (IMP and ribose-5-P matrices only) |

**Verdict.**

- **Corpus coverage.** The corpus already holds the food-side buffer papers: Bell 1997, Potman 1989, Rizzi 2004 and
  Davidek 2002, plus Larralde 1995 from PATHS. A one-hop chain from Larralde returned only prebiotic chemistry.
- **What is still missing:**
  - **A measurement of the 1-deoxy- vs 3-deoxypentosone partition below pH 5 with or without phosphate.** None was found;
    Nedvidek 1992 (pH ≥ 5) remains the only partition data.
  - **Amine (or ammonia) catalysis of pentose enolisation at low pH.** None was found beyond Kim 1996's "2,3-enolisation
    under weakly acidic conditions" statement and Nomi 2021, which observed 1-DP, 3-DP and pentosone by LC-MS/MS in
    carnosine + xylose at pH 5.0 (10.1093/bbb/zbab119; M, no rates in the abstract).
- **For the model.** Kjellstrand's hexose pH profile and Lu & Mosier's xylose profile (PATHS) are shape constraints on the
  acid term. Phosphate's catalytic role (Serianni, Rizzi, Huang) argues for a buffer term in the enolisation law, carried as
  a declared band, rather than an [H+]-only law.

### M8: H2S chemistry in heated sugar solutions

| paper | DOI / access | what it has | class |
|---|---|---|---|
| **Chen, K. Y. & Morris, J. C. (1972).** Kinetics of oxidation of aqueous sulfide by oxygen. Environ. Sci. Technol. 6:529-537 | [10.1021/es60065a008](https://doi.org/10.1021/es60065a008), paywalled | [known content] rate law in [S(-II)] and [O2] vs pH (maximum near pH 8), temperature dependence, trace-metal catalysis, products (S⁰, S2O3²⁻, SO4²⁻) | **Q.** An H2S autoxidation sink competing with thiol oxidation for the B11 oxygen pool. Extrapolating from 25-35 °C to 145 °C is an assumption. Hydroquinone also catalyses it (10.31788/rjc.2020.1315612), which links to melanoidin pots |
| **Lienhard, G. E. & Jencks, W. P. (1966).** Thiol addition to the carbonyl group. Equilibria and kinetics. JACS 88:3982-3995 | [10.1021/ja00969a017](https://doi.org/10.1021/ja00969a017), paywalled | [known content] K and rate constants for hemithioacetal formation of thiols, including H2S-type nucleophiles, with aldehydes, and their acid/base catalysis | **Q.** The class constants for R14 (proposed). Says whether a hemithioacetal or gem-thiol pool of FUR/MGO with H2S, MFT or Cys is large enough to matter |
| Mouret, J.-R. et al. (2015). Gas-liquid partitioning of H2S in model solutions simulating winemaking. JAFC 63:3271-3278 | [10.1021/jf505733a](https://doi.org/10.1021/jf505733a) | [abstract] measured H2S K_i vs temperature and composition | Q, but Sander's compilation already covers H2S Henry constants (VARIABLES). SKIP |

**Verdict.**

- **H2S + furfural and H2S + dicarbonyl rate or equilibrium constants in water at 90-150 °C:** none found. The only
  quantitative H2S-aldehyde scavenging literature is oilfield glyoxal/triazine scavenger work (industrial, mostly
  proprietary).
- **Sealed vs open vessel losses** are a Henry-law and headspace question that VARIABLES already covers.

---

## Ranked download list (new papers only; at most 12)

File names follow the earlier rounds (`<firstauthor><year>.pdf`). None clashes with `data/articles/`.

| # | paper | DOI | item | access | would give | likely numbers? |
|---|---|---|---|---|---|---|
| 1 | Lau et al. 2003, J. Food Eng. 60:397 | [10.1016/s0260-8774(03)00063-3](https://doi.org/10.1016/s0260-8774(03)00063-3) | M3 (finding 2/3) | paywalled; WSU repository record | Ea and k(T) for ribose + amine → norfuraneol; a band on the NF branch's Ea, separate from the furfural branch | yes (0.7) |
| 2 | Schirle-Keller & Reineccius 1992, ACS Symp. Ser. 490:244 | [10.1021/bk-1992-0490.ch019](https://doi.org/10.1021/bk-1992-0490.ch019) | M1, finding 2 | paywalled | zero-order rates at 80-150 °C for furfural, FA and HMF in glucose + Cys water; Ea 117-138 kJ/mol as a hold-out against the 64 kJ/mol trunk | yes (0.9 for Ea, in abstract) |
| 3 | Szajewski & Whitesides 1980, JACS 102:2011 | [10.1021/ja00526a042](https://doi.org/10.1021/ja00526a042) | M2/M6 | paywalled | Brønsted law for thiolate + GSSG exchange, to parameterise MFT/FFT + GSSG (the Zhou sink) | yes at 25 °C (0.9); transfer to 120 °C is an assumption |
| 4 | Pandit et al. 2006, J. Food Eng. 76:353 | [10.1016/j.jfoodeng.2005.05.032](https://doi.org/10.1016/j.jfoodeng.2005.05.032) | M3 | paywalled | a second, independent NF-formation Ea (ribose in potato, amino acid limiting) | yes (0.7) |
| 5 | Kim, Taub, Choi & Prakash 1996, ACS Symp. Ser. 631:54 | [10.1021/bk-1996-0631.ch006](https://doi.org/10.1021/bk-1996-0631.ch006) | M3/M7 | paywalled | M-2 from ribose / ribose-5-P via 2,3-enolisation, pH window; the M-1 (glucose) constants for the trunk | probably (0.55) |
| 6 | Lienhard & Jencks 1966, JACS 88:3982 | [10.1021/ja00969a017](https://doi.org/10.1021/ja00969a017) | M8 / R14 | paywalled | hemithioacetal K and rates: is R14 sequestration of MFT/FFT/H2S by FUR/MGO material? | yes at 25 °C (0.85) |
| 7 | Parliment & Stahl 1995, Dev. Food Sci. 37:805 | [10.1016/s0167-4501(06)80197-4](https://doi.org/10.1016/s0167-4501(06)80197-4) | M1 | paywalled chapter | FFT yields by precursor (furfural vs FA, Cys vs H2S) in aqueous coffee models | maybe (0.4) |
| 8 | Rizzi 1995, Dev. Food Sci. 37:289 (with Rizzi 2000, [10.1021/bk-2000-0754.ch022](https://doi.org/10.1021/bk-2000-0754.ch022)) | [10.1016/s0167-4501(06)80162-7](https://doi.org/10.1016/s0167-4501(06)80162-7) | M1 | paywalled chapters | the acid-catalysed SN1 route from FA/allylic alcohols + S nucleophiles to thiols: a route with an acid dependence furfural + H2S lacks | route mainly (0.3 for yields) |
| 9 | Chen & Morris 1972, ES&T 6:529 | [10.1021/es60065a008](https://doi.org/10.1021/es60065a008) | M8 | paywalled | H2S autoxidation rate law (pH, T, metals): an H2S sink and a claimant on the oxygen pool | yes (0.9); at 25-35 °C |
| 10 | Rizzi, Amba & Heineman 2010, JAFC 58:9739 | [10.1021/jf102433u](https://doi.org/10.1021/jf102433u) | M7/M1 | paywalled | time course of reductones from phosphate-catalysed sugar degradation, 0.1 M phosphate, 100 °C | partly (0.5) |
| 11 | Kjellstrand et al. 2001, Perit. Dial. Int. 21:338 | [10.1177/089686080102100402](https://doi.org/10.1177/089686080102100402) | M7 | paywalled | pH 1-8 profile of glucose degradation at 121 °C, minimum at pH ≈ 3.5; a shape check on the acid term (hexose) | shape (0.6) |
| 12 | Dworkin & Miller 2000, Carbohydr. Res. 329:359 | [10.1016/s0008-6215(00)00204-4](https://doi.org/10.1016/s0008-6215(00)00204-4) | finding 2 (TTCA/Schiff) | paywalled | open-chain aldehyde fraction of ribose, xylose and others, the free-carbonyl factor in the TTCA/Schiff rate law | yes (0.7) |

**Also-rans (M-class; fetch only if the item is opened):**

- Huang, Su & Ho 1998 ([10.1021/jf9705633](https://doi.org/10.1021/jf9705633)): phosphate and the cysteamine thiazolidine.
- Hauck 2003 ([10.1016/j.bbagen.2003.08.002](https://doi.org/10.1016/j.bbagen.2003.08.002)): the ribulose-5-P → NF route.
- Zhang & Serianni 2012 ([10.1021/ja3020296](https://doi.org/10.1021/ja3020296)): phosphate catalysis of osones.
- Kaufmann 2018 ([10.1016/j.foodchem.2018.05.101](https://doi.org/10.1016/j.foodchem.2018.05.101)): general-catalysis rate
  law.
- Novotny 2019 ([10.1021/acs.jafc.8b07100](https://doi.org/10.1021/acs.jafc.8b07100)): 2-acetylthiazole yields.
- Li 2022 ([10.1016/j.fochms.2022.100120](https://doi.org/10.1016/j.fochms.2022.100120), open): β-elimination H2S from protein.
- Lees & Whitesides 1993 ([10.1021/jo00055a016](https://doi.org/10.1021/jo00055a016)): exchange equilibria.
- Gupta 2011 ([10.1016/j.lwt.2011.05.007](https://doi.org/10.1016/j.lwt.2011.05.007)): M-2 vs pH at 105 °C.

**Where nothing exists (after six rounds):**

- the mercaptoacetaldehyde H2S/ring split;
- the Cys-Amadori enolisation rate;
- cysteamine formation from cysteine;
- MFT oxidation kinetics in water;
- FA → FFT and furfural + H2S → FFT rate constants;
- the deoxypentosone partition below pH 5.

These are experiment-shaped gaps, not search-shaped ones. The cheapest decisive measurement is still PATHS's three-temperature
Hofmann-T3 protocol, extended with an FA arm (FA + Cys and FA + Na2S in 0.5 M phosphate pH 3/5/7, sealed).
