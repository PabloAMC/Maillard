# Reaction-network topology vs. isotope-labelling literature

**Wave L1, 2026-08-27.** Independent validation of the *topology* of the model's reaction
network — which atoms come from which precursor, through which intermediate — against
isotope-labelling and pathway-elucidation studies. Written for external review by flavour
chemists.

Rationale: the two worst defects found by this audit campaign so far were topological, not
parametric (a fabricated one-step MFT shortcut; a spurious H2 gate whose removal moved
Hofmann1998 MFT from 5.58x under to 1.45x under with **zero** parameter changes). Barrier
values can be refitted; a wrong carbon skeleton cannot.

## How to read this document

**Verdicts**

| Verdict | Meaning |
|---|---|
| CONFIRMED | Isotope/labelling evidence supports the topology the model encodes. |
| CONTRADICTED | Labelling evidence shows a different atom origin or a different intermediate from the one the model encodes. The discrepancy and a candidate fix are spelled out. |
| PARTIALLY SUPPORTED | The endpoint or one limb is supported; a specific limb, intermediate or co-product is not evidenced. |
| UNTESTED | No isotope or pathway-elucidation study was found. Not a criticism — a statement about the evidence base. |

**Access levels.** Every claim below is tagged with how the source was actually obtained:

* `[FULL]` — full text retrieved and read.
* `[ABS]` — publisher/PubMed abstract retrieved verbatim (this is most of the corpus; ACS
  abstracts carry the labelling results in the abstract itself).
* `[META]` — bibliographic record verified via Crossref/OpenAlex only; no abstract exists or
  none was retrievable (typical for pre-1980 JAFC and ACS Symposium Series chapters).
* `[SECONDARY]` — the finding reached me through a search-engine summary of the abstract, not
  by direct retrieval. **Never used to set a verdict**, only as corroboration.

Nothing below is asserted from memory. Where no source was found, the row says UNTESTED.

## Measured route inventory (what the network actually emits)

Enumerated live rather than read off the source, so the verdicts below attach to routes that
really fire. `SmirksEngine.enumerate(..., max_generations=3)` at pH 5.5 / 150 °C, three
canonical systems; the cell gives the reaction family (or families) that produce the compound:

| Product | ribose + cysteine | glucose + cysteine | glucose + glycine |
|---|---|---|---|
| 2-methyl-3-furanthiol | `Thiol_Addition_Norfuraneol` (**sole producer**) | `Thiol_Addition_Hexose_Legacy_Shortcut` | absent |
| 2-furfurylthiol | `Thiol_Addition_H2`, `Thiol_Dehydration` | `Thiol_Addition_H2`, `Thiol_Dehydration` | absent |
| DMHF (furaneol) | absent | `Furanone_Cyclisation` | `Furanone_Cyclisation` |
| norfuraneol | `Furanone_Cyclisation` | absent | absent |
| furfural | `Enolisation_1_2` | `Enolisation_1_2` | `Enolisation_1_2` |
| HMF | absent | `Enolisation_1_2` | `Enolisation_1_2` |
| 2,5-dimethylpyrazine | `Aminoketone_Condensation` | `Aminoketone_Condensation` | `Aminoketone_Condensation` |

Two facts to carry into the reading below. First, in the model's flagship system
(ribose + cysteine) the route this document finds **contradicted** — norfuraneol → MFT — is
the *only* producer of MFT. Second, DMHF from a hexose is reachable in both hexose systems,
but see topology risk 3 for what it is gated on.

---

## Summary table

| # | Product | Model route (file:line) | Isotope / pathway evidence | Verdict | Citation |
|---|---|---|---|---|---|
| 1 | **2-Methyl-3-furanthiol** — carbon origin | intact pentose C5 throughout `_norfuraneol_mft_route` (`src/reaction_templates.py:275-373`) | [13C5]ribose 1:1 mix: ribose fragmentation "did not play a significant role"; C skeleton intact. [1-13C]ribose gives 2-[13C]**methyl**-3-furanthiol — model's atom mapping (sugar C-1 → 1-deoxyosone CH3 → norfuraneol 5-CH3 → MFT 2-CH3) agrees | **CONFIRMED** | Cerny & Davidek 2003 `[ABS]`; Cerny & Davidek 2004 `[ABS]` |
| 2 | **MFT** — *norfuraneol* as the intermediate | `Furanone_Cyclisation` → `Thiol_Addition_Norfuraneol` (`reaction_templates.py:333-371`); the primary MFT lane since Wave G1 | Spiking experiment: norfuraneol + [13C5]ribose + cysteine gave MFT that was *mainly 13C5-labelled* ⇒ "4-hydroxy-5-methyl-3(2H)-furanone is **unimportant as an intermediate**". Independently: "Thiamin and **norfuraneol/cysteine were less effective precursors of MFT**" | **CONTRADICTED** | Cerny & Davidek 2003 `[ABS]`; Hofmann & Schieberle 1998 `[ABS]` |
| 3 | **MFT** — hexose legacy shortcut | one-step hexose 3-deoxyosone + H2S + 2[H] → MFT + CH2O + 3 H2O (`reaction_templates.py:870-969`) | Hexoses *do* yield MFT with cysteine ("glucose and rhamnose also gave significant yields"), so a hexose channel must exist; but no labelling study proposes a 3-deoxyosone route, and the evidenced hexose intermediates are the 1-deoxyosone/acetylformoin family | **PARTIALLY SUPPORTED** (endpoint) / mechanism UNTESTED | Hofmann & Schieberle 1998 `[ABS]`; Wang & Ho 2008 `[ABS]` |
| 4 | **MFT** — thiamine cascade | thiamine → thiazole → 5-hydroxy-3-mercapto-2-pentanone → MFT (`reaction_templates.py:1433-1516`) | [13C5]xylose/cysteine/thiamin: xylose-unlabelled sulfur volatiles "apparently stemmed from thiamin"; MFT gets carbon from **both** xylose and thiamin, "therefore different formation pathways must exist". The cascade's key intermediate is a literature-identified compound of exactly this system | **CONFIRMED** | Cerny & Briffod 2007 `[ABS]`; Cerny & Guntz-Dubini 2008 `[META]`; Bolton et al. 1994 34S `[SECONDARY]` |
| 5 | **2-Furfurylthiol** | furfural + H2S (+2[H]) → FFT (`reaction_templates.py:782-867`) | "confirm furan-2-carbaldehyde as an intermediate of 2-furfurylthiol"; [1-13C]ribose → 2-[13CH2]furfurylthiol, i.e. sugar C-1 = furfural's aldehyde C = FFT's exocyclic CH2SH; xylose "the exclusive carbon source" for both furfural and FFT; furfural/H2S "10 times higher efficiency in generating FFT" than the alternatives tested | **CONFIRMED** | Cerny & Davidek 2004 `[ABS]`; Cerny & Briffod 2007 `[ABS]`; Hofmann & Schieberle 1998 `[ABS]` |
| 6 | **FFT** — thiohemiacetal sub-route | `Thiohemiacetal_Formation` → `Thiol_Dehydration` (`reaction_templates.py:849-865`) | No labelling study addresses the intermediate | **UNTESTED** | — |
| 7 | **Furfural** from pentose | 3-deoxypentosone → furfural + 2 H2O (`reaction_templates.py:186-199`) | [13C5]xylose: "for 2-furaldehyde ... xylose was the exclusive carbon source"; deoxyosone intermediacy supported by phosphorylated-glucose labelling work | **CONFIRMED** (origin) / PARTIALLY (intermediate) | Cerny & Briffod 2007 `[ABS]`; Yaylayan et al. 2003 `[ABS]` |
| 8 | **Furfural** from hexose (secondary branch) | hexose 3-deoxyosone → furfural + **formaldehyde** + 2 H2O (`reaction_templates.py:201-207`) | Labelled glucose/alanine: glucose "can lose the C-6 atom to produce a pentose moiety responsible for the formation of furfural" — one terminal carbon is lost, as the model encodes; the specific 3-deoxyosone precursor and the formaldehyde identity of the lost C1 are not established | **PARTIALLY SUPPORTED** | Yaylayan & Keyhani 2000 `[ABS]` |
| 9 | **HMF** from hexose | hexose 3-deoxyosone → HMF + 2 H2O (`reaction_templates.py:190`), reached by **both** glucose (Amadori) and fructose (Heyns) | 3-Deoxyglucosone converts to HMF far faster than glucose does (supports the glucose route), **but** "both fructose and sucrose showed much higher conversion rates than 3-deoxyglucosone thus **precluding it as a major precursor of HMF in fructose** and sucrose solutions" | **PARTIALLY SUPPORTED** (glucose) / **CONTRADICTED** (fructose) | Perez Locas & Yaylayan 2008 `[ABS]` |
| 10 | **Furaneol / HDMF (DMHF)** from hexose | hexose 1-deoxyosone + 2[H] → DMHF + 2 H2O (`reaction_templates.py:343-350`) | CAMOLA-style 1:1 [13C6]/[12C6]glucose gave a 1:1 mixture of [13C6]- and [12C6]-DMHF ⇒ **intact C6 skeleton**, no C3+C3 recombination; acetylformoin (the cyclised 1-deoxyosone) detected as the DMHF precursor in glucose/amino-acid systems | **CONFIRMED** | Wang & Ho 2008 `[ABS]`; corroborated Spreng et al. 2021 `[ABS]`, Schieberle 2005 `[ABS]` |
| 11 | **DMHF / HEMF** from pentose + amino acid | pentose + glycine → DMHF; pentose + alanine → HEMF (`reaction_templates.py:1519-1575`) | 13C-labelled glycine and alanine: "incorporation of the Strecker degradation products **formaldehyde and acetaldehyde** into the pentose moiety, forming the furanones HDMF and HEMF"; mechanism = Amadori → 2,3-enolisation → chain elongation by the Strecker aldehyde → **reduction** of acetylformoin-type intermediates | **CONFIRMED** (C5+C1 / C5+C2 stoichiometry exactly as modelled) | Blank & Fay 1996 `[ABS]` |
| 12 | **Norfuraneol** (as a species / cyclisation step) | pentose 1-deoxyosone → norfuraneol + H2O (`reaction_templates.py:333-342`) | Norfuraneol is a genuine Maillard product and reacts with H2S (van den Ouweland & Peer's original synthesis); in situ it is the demonstrated precursor of **2-mercapto-3-pentanone**, not of MFT | **PARTIALLY SUPPORTED** (keep the species and the step; retire its MFT role) | Van den Ouweland & Peer 1975 `[META]`; Cerny & Davidek 2003 `[ABS]` |
| 13 | **Retro-aldol fragmentation** of the deoxyosones | 3-deoxyosone → pyruvaldehyde + glyceraldehyde / glycolaldehyde; second channel → glyoxal + ketol (`reaction_templates.py:561-606`) | Labelled glucose: pyruvaldehyde and acetol "incorporated intact C1-C2-C3 and C4-C5-C6 carbon chains"; glycolaldehyde "incorporated intact C5-C6 and C1-C2 chains" (70/30 split) — i.e. both C3+C3 and C2+C4 cleavages are real | **CONFIRMED** | Yaylayan & Keyhani 2000 `[ABS]`; Yaylayan & Keyhani 2001 `[ABS]` |
| 14 | **Strecker aldehydes** (methional, 3-methylbutanal, phenylacetaldehyde) | α-dicarbonyl + amino acid → aldehyde + aminoketone + CO2 (`reaction_templates.py:376-485`) | [15N]- and [methyl-13C]methionine confirm the amino-acid origin of the Strecker aldehyde. **But** the α-dicarbonyl is not the only initiator: lipid-derived 4-oxo-2-alkenals, 2,4-alkadienals and 4,5-epoxy-2-alkenals raise phenylacetaldehyde yields by 300-900 % and have measured Ea 28-67 kJ/mol | **PARTIALLY SUPPORTED** (route real, initiator set incomplete) | Yaylayan & Keyhani 2001 `[ABS]`; Zamora et al. 2013 `[ABS]`; Zamora et al. 2015 `[ABS]` |
| 15 | **Methanethiol, DMDS, DMTS** | methional → MeSH + acrolein; 2 MeSH → DMDS; DMDS + MeSH → DMTS (`reaction_templates.py:1016-1060`) | Precursor-addition studies: "methionine was the precursor of dimethyl disulfide and dimethyl trisulfide"; "methionine was the major source of methanethiol". The acrolein co-product is not isotope-traced anywhere I found | **PARTIALLY SUPPORTED** | Pan et al. 2021 `[ABS]`; Cheng et al. 2020 `[ABS]` |
| 16 | **Dimethyl sulfide (DMS)** | **not implemented** | DMS's precursor is S-methylmethionine, not methional/methionine: SMM→DMS yields 24-27 mol %; "S-methyl methionine was the precursor of dimethyl sulfide, and methionine was the precursor of dimethyl disulfide and dimethyl trisulfide" | **CORRECT OMISSION** — SMM is not a registered precursor, so DMS *should* be absent | Scherb et al. 2009 `[ABS]`; Pan et al. 2021 `[ABS]` |
| 17 | **Acrylamide** | Asn + reducing sugar → acrylamide + CO2 + NH3 + 3-deoxyosone (`reaction_templates.py:1316-1397`); curated layer emits the substituted-imine branch (`src/curated_pathways.py:153-156`) | "The mechanism involves formation of a Schiff base followed by decarboxylation and elimination of **either ammonia or a substituted imine**... Isotope substitution studies and mass spectrometric analysis of heated model systems confirm the presence of key reaction intermediates." Both literature branches are represented across the two layers | **CONFIRMED** | Zyzak et al. 2003 `[ABS]` |
| 18 | **Acrylamide** — carbonyl-partner identity | no identity term in either layer (single lumped `reducing_sugar_mM`) | The N-glycosyl of Asn gives ~2.4 mmol/mol vs 0.1-0.2 for the Amadori compound, and α-hydroxycarbonyls (hydroxyacetone, >4 mmol/mol) beat α-dicarbonyls by an order of magnitude | **CONTRADICTED as a modelling assumption** (sugar/carbonyl identity is not neutral) | Stadler et al. 2004 `[ABS]` |
| 19 | **2,5-Dimethylpyrazine** | 2 aminoacetone → 2,5-DMP + 2 H2O + 2[H] (`reaction_templates.py:530-558`) | Label-incorporation patterns *of the methyl/dimethylpyrazines themselves* were used to establish the C2/C3 sugar-fragment origins, i.e. dimethylpyrazines are built from two small sugar fragments plus amino nitrogen. The aminoacetone → 3,6-dihydro-2,5-DMP → (oxidation) → 2,5-DMP sequence is independently established in isotope-traced microbial systems. CAMOLA on coffee: sugar is the predominant backbone source but only 33-55 % of pyrazine nitrogen comes from free amino acids | **PARTIALLY SUPPORTED** | Yaylayan & Keyhani 2000 `[ABS]`; Zhang et al. 2020 `[ABS]`; Molla et al. 2015 `[ABS]`; Spreng et al. 2021 `[ABS]` |
| 20 | **2,3-Dimethylpyrazine, 2-ethyl-3,5-dimethylpyrazine** | **not implemented** (both are scored targets in `data/species/desirable_targets.yml`) | Isotope-traced route for EDMP exists: 2,3-pentanedione + aminoacetone | **COMPLETENESS GAP** | Zhang et al. 2021 `[ABS]` |
| 21 | **2-Alkylthiazoles** | Strecker aldehyde + glyoxal + NH3 + H2S → thiazole + 3 H2O (`reaction_templates.py:705-779`); lipid analogue at `:1166-1214` | No isotope study found for this condensation. Note the one thiazole that *is* isotope-resolved in a comparable system — 5-(2-hydroxyethyl)-4-methylthiazole — comes from **thiamin**, and the model does get that one from thiamin | **UNTESTED** | (Cerny & Briffod 2007 `[ABS]` for the thiamin thiazole only) |
| 22 | **Hexanal** (empirical lipid lane) | hydroperoxide pool × 0.37 (`src/lipid_oxidation.py:398-453`); SMIRKS lane = alkoxy β-scission (`src/smirks_engine.py:157-172`) | Full-text isomer-resolved study: **13-HpODE → hexanal, vinyl hexanoate, 2-pentylfuran, 4,5-epoxy-2-decenal**; 9-HpODE → 2,4-decadienal, octanoic acid, 2-octenal, hexanal. Hexanal + 2-pentylfuran co-occurring from 13-HpODE matches the model's branching set. The same paper argues hexanal there is *not* formed by simple β-scission: "3-nonenal was not detected. Therefore ... the unfavorable β-scission would not be occurred" | **PARTIALLY SUPPORTED** (products right; the model's β-scission mechanism is questioned by the one full-text study found) | Miyazaki et al. 2023 `[FULL]` |
| 23 | **Nonanal** | 12-15 % of a hydroperoxide pool computed **from the linoleic fraction only** (`lipid_oxidation.py:437` uses `linoleic_acid_pct`; `:451` and `:395` emit nonanal) | Nonanal appears in the decomposition products of **neither** linoleate hydroperoxide isomer. Nonanal is the C9 fragment of the **oleate** Δ9 double bond (recovered at 30 ± 3 % carbon yield from oleic acid cleavage). `LipidProfile.oleic_acid_pct` is declared at `lipid_oxidation.py:145` and **never read anywhere in the repo** | **CONTRADICTED** | Miyazaki et al. 2023 `[FULL]`; Zahardis & Petrucci-type oleate cleavage product studies `[ABS]` (Reynolds et al. 2006; Hearn & Smith 2007) |
| 24 | **CML / CEL** | lysine + glyoxal → CML; lysine + methylglyoxal → CEL, single condensation, no oxidant (`reaction_templates.py:1400-1428`) | In glucose/lysine, "approximately 50 % of the CML ... originates from oxidation of Amadori product, and 40-50 % originates from a pre-Amadori stage"; glyoxal trapping blocked 50 % of CML. Crucially for the model's *balance*: "Initial CML formation rate from glyoxal was not dependent on oxidation, suggesting an intramolecular Cannizzaro reaction" — so a non-oxidative one-step lump is right **for the glyoxal half** | **PARTIALLY SUPPORTED** (~half the literature flux is unmodelled) | Glomb & Monnier 1995 `[ABS]` |
| 25 | **2-Acetyl-1-pyrroline** | **not implemented**; proline and ornithine are not in `data/species/precursors.yml` | CAMOLA **in an extrusion cooker** (135 °C, 20 % moisture, rice model): 2-AP forms "(i) by acylation of 1-pyrroline via **C2 sugar fragments** (major pathway) and (ii) via ring-opening of 1-pyrroline incorporating **C3 sugar fragments** (minor)" | **COMPLETENESS GAP** (and the study is in this repo's own process domain) | Davidek et al. 2013 `[ABS]` |
| 26 | **Cysteine-S-conjugates** (bound-thiol sink) | **not implemented** | S-furfuryl-L-cysteine and S-(2-methyl-3-furyl)-L-cysteine identified as Maillard products of xylose + cysteine (100 °C, 2 h); odourless, with furfuryl alcohol proposed as the FFT-S-Cys intermediate | **COMPLETENESS GAP** (affects free-thiol mass balance) | Cerny & Guntz-Dubini 2013 `[ABS]` |
| 27 | **Reducing-equivalent token** `[HH]` and the disulfides | 2 cysteine → cystine + 2[H] (`reaction_templates.py:645-702`); 2 MFT → bis(2-methyl-3-furyl) disulfide + H2 (`:972-1013`) | No isotope study bears on the model's bookkeeping token; thiol→disulfide oxidation is textbook but untraced here | **UNTESTED** | — |
| 28 | **Hexanal trapping by amines** (off-note masking) | hexanal + Gly/Lys → Schiff base (`curated_pathways.py:143-148`; SMIRKS `Lipid_Schiff_Base`) | Carbonyl-amine reactions between aldehydes and amines are real and consume the carbonyl, but the imine is an *intermediate* that goes on to further products, not a stable sink | **PARTIALLY SUPPORTED** | Zamora & Hidalgo 2012 `[ABS]` |

---

## Per-product detail

### 1-2. 2-Methyl-3-furanthiol — the flagship, and the campaign's largest surviving topology risk

**What the model does.** Since Wave G1 the primary MFT lane is the three-step route
Amadori → (2,3-enolisation) → 1-deoxy-2,3-pentodiulose → (cyclodehydration) →
**norfuraneol** → (+ H2S + 2[H]) → MFT, cited to van den Ouweland & Peer 1975. Wave H fitted
`thiol_addition_norfuraneol` (28.60 → 26.85 kcal/mol) against `cys_ribose_140C_Hofmann1998`,
which Wave I then recorded as *the only surviving literature anchor on the entire sulfur
branch*. [SUPERSEDED 2026-08-27, Waves S2b/S2c: that benchmark's 342/200 ppb targets were
later proven to be a repo-internal derivation from an abstract-reconstructed brief
(`maillard_validation_benchmarks.md` §1.3), not literature values. The benchmark is retired
to REFERENCE tier and the sulfur branch has ZERO absolute literature anchors. The isotope
findings in this document are unaffected — they concern route topology, not yields.]

**What the isotope evidence says.** Two independent results, from two laboratories, both
say the norfuraneol → MFT step is not the in-situ route.

Cerny & Davidek 2003 (JAFC 51:2714-2721) ran the decisive competition experiment — spiking
authentic norfuraneol into a [13C5]ribose/cysteine system, so that the two candidate
precursors carry different masses:

> "In another trial cysteine, 4-hydroxy-5-methyl-3(2H)-furanone and [13C5]ribose were
> reacted under the same conditions. The resulting 2-methyl-3-furanthiol was mainly
> 13C5-labeled, suggesting that it stems from ribose and that **4-hydroxy-5-methyl-3(2H)-furanone
> is unimportant as an intermediate**. ... A new reaction pathway from ribose via its
> **1,4-dideoxyosone** is proposed, which explains both the formation of 2-methyl-3-furanthiol
> without 4-hydroxy-5-methyl-3(2H)-furanone as an intermediate and a new way to form
> 3-mercapto-2-pentanone." `[ABS]`

Cerny & Davidek 2004 (JAFC 52:958-961) then confirmed the proposed intermediate positionally,
with [1-13C]ribose:

> "The results confirm furan-2-carbaldehyde as an intermediate of 2-furfurylthiol, as well as
> **1,4-dideoxypento-2,3-diulose as an intermediate of 2-methyl-3-furanthiol** and
> 3-mercaptopentan-2-one." `[ABS]`

And Hofmann & Schieberle 1998 (JAFC 46:235-241) — *the repo's own fit target* — says in its
abstract:

> "Studies on several intermediates indicated the highest yields for MFT (1.4 mol %) when
> **hydroxyacetaldehyde and mercapto-2-propanone** were reacted for 6 min at 180 °C in the
> absence of water. ... **Thiamin and norfuraneol/cysteine were less effective precursors of
> MFT.** The results imply that different formation pathways may run in parallel during food
> processing." `[ABS]`

**Discrepancy, stated precisely.**

* The model's *carbon-skeleton* claim is right. Cerny & Davidek 2003 show ribose fragmentation
  "did not play a significant role" and the C5 skeleton stays intact for MFT, FFT and
  3-mercapto-2-pentanone. The model's route is intact-C5 end to end. The 2004 positional
  result (`2-[13C]methyl-3-furanthiol` from `[1-13C]ribose`) also matches the model's implicit
  atom mapping: sugar C-1 → the 1-deoxyosone methyl → norfuraneol's 5-methyl → MFT's 2-methyl.
* The model's *intermediate* is contradicted. The evidenced intermediate is
  1,4-dideoxypento-2,3-diulose (CH3-CO-CO-CH2-CH2OH, C5H8O3), which is **absent from the
  network** — grep `_DEOXYOSONE_1_PENTOSE` and there is only the 1-deoxy species,
  `CC(=O)C(=O)C(O)CO` (C5H8O4).

**What a fix looks like.** Adding the 1,4-dideoxyosone is unusually cheap, and it removes a
piece of fiction rather than adding one:

```
  1-deoxy-2,3-pentodiulose  C5H8O4 + 2[H] -> 1,4-dideoxy-2,3-pentodiulose C5H8O3 + H2O
  1,4-dideoxy-2,3-pentodiulose C5H8O3 + H2S -> MFT C5H6OS + 2 H2O          [EXACT, no 2[H]]
```

The second step balances with **no reducing-equivalent token at all** — the `[HH]` fiction
that caused red-team finding H4 is not needed on the sulfur-incorporation step under the
literature topology; it moves upstream onto a dehydration/reduction of the 1-deoxyosone,
where an enediol/reductone donor is chemically ordinary. Norfuraneol should be *kept* (it is
a real product, and Cerny & Davidek show it is the demonstrated precursor of
2-mercapto-3-pentanone) but demoted out of the MFT lane.

**Exposure.** Measured, not inferred: in ribose + cysteine at pH 5.5 / 150 °C,
`Thiol_Addition_Norfuraneol` is the **sole** producer of MFT in the enumerated network (see
the route inventory above). There is no second channel to fall back on, so whatever fraction
of the real MFT flux the literature assigns to the 1,4-dideoxyosone route is currently being
carried by a step the literature says is unimportant.

**Why this matters beyond chemistry.** `thiol_addition_norfuraneol` was fitted, by this
campaign's own generator, to Hofmann1998 — a paper whose abstract explicitly ranks
norfuraneol/cysteine as a *less effective* MFT precursor. The fitted barrier is therefore
absorbing a route-selection error. Any re-run of the sulfur refit (open item 1 in the Wave I
carry-forward list) should be done **after** the route change, not before, or it will fit
the wrong route again more precisely.

### 3. MFT from hexoses (the demoted legacy shortcut)

Restricting the one-step lump to hexoses (Wave I fix 12) was the right call directionally,
but it is not the literature's hexose route. Hofmann & Schieberle 1998 confirm hexoses reach
MFT ("glucose and rhamnose also gave significant yields"), and their *highest-yielding* MFT
system of all was a **C2 + C3 recombination** (hydroxyacetaldehyde + mercapto-2-propanone,
1.4 mol %) in the absence of water — a topology the network cannot express at all, since
mercapto-2-propanone is not a species anywhere. Under low-moisture conditions (extrusion,
roasting), that is likely the *dominant* channel and the model has no representation of it.

### 4. MFT and FFT from thiamine

Cerny & Briffod 2007 (JAFC 55:1552-1556), with [13C5]xylose, cysteine and thiamin at pH 4-7:

> "For 2-furaldehyde and 2-furfurylthiol, which were favored at low pH, the labeling
> experiments clearly indicated that xylose was the exclusive carbon source. On the other
> hand, xylose was virtually not involved in the formation of 3-mercapto-2-butanone,
> 4,5-dihydro-2-methyl-3-furanthiol, and 5-(2-hydroxyethyl)-4-methylthiazole, which apparently
> stemmed from thiamin degradation. **Both xylose and thiamin seemed to significantly
> contribute to the formation of 2-methyl-3-furanthiol**, 3-mercapto-2-pentanone, and
> 2-mercapto-3-pentanone, and therefore **different formation pathways must exist for each of
> these molecules**." `[ABS]`

The model's thiamine cascade (bridge hydrolysis → thiazolium ring opening →
5-hydroxy-3-mercapto-2-pentanone → MFT) is therefore the right *lane*, and its key
intermediate is a compound identified in exactly this system (Cerny & Guntz-Dubini 2008,
"Identification of 5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine,
Cysteine, and Xylose", JAFC 56:10679-10682 `[META]` — title and metadata verified via
Crossref; the abstract was not retrievable).

**Consequence for the panel.** The repo's `thiamine_cys_glucose_120C_Bolton1994` benchmark
(current status: 724x under-prediction, pinned as an honest failure) is sourced from
Bolton, Reineccius, Liardon & Ba 1994, *"Role of Cysteine in the Formation of
2-Methyl-3-furanthiol in a Thiamine—Cysteine Model System"* (ACS Symp. Ser. 543:270-278,
10.1021/bk-1994-0543.ch022 `[META]` — metadata verified via Crossref). That chapter is itself
a **34S-labelling study**: reported secondary-source figures are that only ~7.5-8.0 % of the
MFT carried 34S from labelled cysteine, i.e. thiamine supplied the great majority of the
sulfur `[SECONDARY]` — I could not retrieve the chapter text and this number must be
verified against the original before being used for anything. If it holds, it says the
model's 724x miss on that benchmark is a *magnitude* problem in a lane whose topology is
correct — which is a materially different diagnosis from the sulfur-branch story the repo
currently tells, and it is checkable with one library visit.

### 5. 2-Furfurylthiol — the best-evidenced route in the model

Three independent labelling results converge on exactly what the model encodes:

* Cerny & Davidek 2004: "confirm **furan-2-carbaldehyde as an intermediate of
  2-furfurylthiol**", and the product from [1-13C]ribose is `2-[13CH2]furfurylthiol` — the
  sugar's C-1 ends up as the exocyclic CH2-SH carbon. In the model, furfural's aldehyde carbon
  is the one that becomes CH2SH. Same atom. `[ABS]`
* Cerny & Briffod 2007: xylose is the *exclusive* carbon source for both furfural and FFT. `[ABS]`
* Hofmann & Schieberle 1998: "the system **furan-2-aldehyde/H2S showed a 10 times higher
  efficiency in generating FFT**" than the C2+C3 intermediates tested. `[ABS]`

The thiohemiacetal intermediate the model interposes (`Thiohemiacetal_Formation` →
`Thiol_Dehydration`) is chemically reasonable and **untested** by any labelling study I found.
Not a defect; just unevidenced, and it is a second parallel emission of the same product,
which is worth remembering when reading FFT flux attributions.

### 8-9. Furfural and HMF: where the model's hexose chemistry diverges from the labels

Yaylayan & Keyhani 2000 (JAFC 48:2415-2419), labelled glucose + alanine:

> "glucose in the presence of L-alanine can lose either **C-1** atom to produce a pentitol
> moiety responsible for the formation of **furanmethanol** or it can lose the **C-6** atom to
> produce a pentose moiety responsible for the formation of **furfural**." `[ABS]`

The model implements the second (furfural + one lost carbon, emitted as formaldehyde) and not
the first — **furfuryl alcohol is absent from the network**, which is a double loss because
Cerny & Guntz-Dubini 2013 propose furfuryl alcohol as the intermediate of the odourless
S-furfuryl-cysteine conjugate.

For HMF, Perez Locas & Yaylayan 2008 (JAFC, 10.1021/jf8010245):

> "Glucose exhibited a much lower conversion rate than 3-deoxyglucosone, however, both
> **fructose and sucrose showed much higher conversion rates than 3-deoxyglucosone thus
> precluding it as a major precursor of HMF in fructose and sucrose solutions**." `[ABS]`

The model sends fructose through Heyns → 3-deoxyosone → HMF, i.e. through the precursor this
labelling study excludes for fructose. Fructose is a registered precursor
(`data/species/precursors.yml:152`), so this is live. The literature alternative is a direct
fructofuranosyl-cation route, which is a different topology, not a different rate constant.

### 10-11. Furaneol / HDMF: the model gets this one right, twice

Wang & Ho 2008 (JAFC 56:7405-7409) ran the CAMOLA design directly:

> "Experiments using a 1:1 mixture of [13C6]glucose and [12C6]glucose indicate that in the
> presence of glycine or cysteine, **glucose skeleton kept intact during DMHF formation** since
> a 1:1 mixture of [13C6]DMHF and [12C6]DMHF was formed. **Acetylformoin was detected in the
> glucose with amino acid reaction system as a precursor of DMHF**, while in the MG reaction
> systems, acetylformoin could not be identified." `[ABS]`

A 1:1 labelled/unlabelled product distribution with **no intermediate isotopologues** is the
signature of an intact skeleton; had DMHF been assembled from two C3 fragments, the CAMOLA
experiment would have produced the 1:2:1 mixed pattern. The model's hexose route
(1-deoxyosone → DMHF, where the cyclised 1-deoxyosone *is* acetylformoin) is exactly this
topology, including the reduction step.

The pentose route is equally well evidenced. Blank & Fay 1996 (JAFC 44:531-536):

> "Experiments using 13C-labeled glycine and alanine suggest the incorporation of the
> Strecker degradation products **formaldehyde and acetaldehyde** into the pentose moiety,
> forming the furanones HDMF and HEMF, respectively. The presence of 12C-HDMF, which was
> approximately **30 % of the total HDMF** amount found in xylose/glycine, indicates that HDMF
> is partly formed by sugar fragmentation. The proposed mechanism ... is based on decomposition
> of the Amadori compound via 2,3-enolization, chain elongation by the Strecker aldehydes, and
> **reduction** of the resulting acetylformoin-type intermediates." `[ABS]`

The model's `_furanone_generation` stoichiometry (pentose + glycine − CO2 → C6 DMHF; pentose +
alanine − CO2 → C7 HEMF) matches the C5+C1 / C5+C2 accounting. Two caveats worth recording:
the ~30 % sugar-fragmentation channel is not modelled, and Blank & Fay report HDMF from
pentose/**alanine** as well as pentose/glycine, whereas the model's alanine branch emits only
HEMF.

### 14. Strecker aldehydes: right route, incomplete initiator set

The amino-acid origin of the aldehyde is not in doubt (Yaylayan & Keyhani 2001 used
[15N]methionine and [methyl-13C]methionine and detected the Strecker aldehyde). What the model
misses is that in any system containing oxidising lipid — which is the *entire point* of this
repo's off-note lane — the α-dicarbonyl is not the only Strecker initiator:

> "all other lipid oxidation products assayed increased the amount of phenylacetaldehyde
> produced by **300-900 %**, with [4-oxo-2-nonenal] being the most reactive compound"
> (Zamora & Hidalgo 2012, `[ABS]`)

and the kinetics have been measured: Ea 28-38 kJ/mol for alkadienals, 55-67 kJ/mol for
oxoalkenals and epoxyalkenals (Zamora et al. 2013 `[ABS]`). The model's only lipid×Maillard
crosstalk family (`Lipid_Strecker_Synergy`) makes **thiazoles**, not Strecker aldehydes, so
this coupling is absent in the direction the literature says is strongest.

### 22-23. The lipid lane: hexanal is defensible, nonanal is not

Miyazaki et al. 2023 (*Biosci. Biotechnol. Biochem.* 87:179-190, 10.1093/bbb/zbac189) is the
one source here read in **full text**. Isomer-resolved decomposition:

* 13-HpODE → **hexanal**, vinyl hexanoate, **2-pentylfuran**, 4,5-epoxy-2-decenal
* 9-HpODE → 2,4-decadienal, octanoic acid, 2-octenal, hexanal

Two consequences for the model.

*Hexanal and 2-pentylfuran*: co-production from 13-HpODE supports the model's fixed branching
set {hexanal 0.37, 2-pentylfuran 0.08} at the level of *which* products appear. But the paper
disputes the mechanism the SMIRKS lane uses: "3-nonenal was not detected. Therefore ... the
unfavorable β-scission would not be occurred", with hexanal instead attributed to
decomposition of `12-((1-hydroperoxyhexyl)oxy)dodeca-9,11-dienoic acid` formed from an epoxy
allyl radical. The model's `beta_scission_alkoxy` rule (widened in Wave G1 specifically so
hexanal became structurally reachable) is therefore *classical* but not confirmed by the one
isomer-resolved study found. Marked PARTIALLY SUPPORTED, not CONTRADICTED — the classical
β-scission literature (Frankel) was not retrieved and one paper does not overturn it.

*Nonanal*: **contradicted**. Nonanal is in neither isomer's product list. Nonanal is the C9
aldehyde produced by cleaving the **oleate** Δ9 double bond — recovered at 30 ± 3 % carbon
yield among the four canonical C9 products of oleic acid cleavage (nonanal, nonanoic acid,
9-oxononanoic acid, azelaic acid; Reynolds et al. 2006 `[ABS]`, Hearn & Smith 2007 `[ABS]`
— both atmospheric-chemistry ozonolysis studies of the same double bond, cited here only for
the structural provenance of the C9 fragment, not for food-relevant kinetics). The model
nonetheless computes its entire hydroperoxide pool from `linoleic_acid_pct` and then assigns
12-15 % of it to nonanal, while `LipidProfile.oleic_acid_pct` — declared, populated for pea
and soy, and carried through the calibration registry — is **read by no code path in the
repo** (verified by grep across `src/`, `scripts/` and `tests/`). Same for
`alpha_linolenic_pct`.

This is a topology error with live consequences: nonanal is one of four
`LIPID_OXIDATION_MARKER_NAMES`, it is a registered precursor, and the off-note lane's
composition claims flow from it. Minimal honest fix: compute a second hydroperoxide pool from
`oleic_acid_pct` and source nonanal from it, or (if that is too large a change) mark the
nonanal branching ratio as unanchored in `data/lit/lipid_oxidation_calibration.json` and say
in the docs that nonanal is currently a proxy, not a mechanistic prediction.

### 24. CML/CEL: the right half, modelled correctly

Glomb & Monnier 1995 (JBC 270:10017-10026):

> "**Initial CML formation rate from glyoxal was not dependent on oxidation, suggesting an
> intramolecular Cannizzaro reaction.** CML formation from glucose/lysine or Amadori product of
> both was strongly dependent on oxidation. ... it can be estimated that approximately **50 %
> of the CML forming in a glucose/lysine system originates from oxidation of Amadori product,
> and 40-50 % originates from a pre-Amadori stage** largely independent from glucose
> autoxidation." `[ABS]`

The model's `Safety_Risk_AGE` step (lysine + glyoxal → CML, atom-balanced, no oxidant, no
water) is a faithful lump of the glyoxal half — including, correctly, the absence of an
oxidant, which the Cannizzaro finding justifies. The Amadori-oxidation half (~50 % of the
flux) has no representation. Given that Wave I measured model CML at ~1204x below its
literature band, "half the routes are missing" is a partial but real explanation that should
be tried before any parameter is touched.

### 25. 2-Acetyl-1-pyrroline: the most domain-relevant completeness gap

Davidek et al. 2013 (JAFC 61:10215-10219) ran CAMOLA **inside a twin-screw extruder** on a
rice model recipe with [U-13C6]-D-glucose, proline and ornithine:

> "2-Acetyl-1-pyrroline is formed (i) by acylation of 1-pyrroline via **C2 sugar fragments**
> (major pathway) and (ii) via ring-opening of 1-pyrroline incorporating **C3 sugar fragments**
> (minor pathway) ... In addition, it has been shown that the formation of 2-acetyl-1-pyrroline
> in low-moisture systems depends on the pH value of the reaction mixture." `[ABS]`

The model has neither proline nor ornithine in its precursor registry and no 2-AP route. 2-AP
is a defining roast odorant of extruded cereal/protein systems — this repo's stated
application domain — and the labelling work needed to implement it correctly (which sugar
fragment, which pH dependence, at extrusion conditions) already exists. The C2 acyl donor
(glycolaldehyde/hydroxyacetaldehyde) and the C3 donor (pyruvaldehyde) are both already
species in the network, so the missing pieces are the amino-acid precursors and one
condensation family.

---

## Summary

**Counts over 28 assessed routes / route-absences** (29 verdicts: row 9 splits, HMF being
supported for glucose and contradicted for fructose):

| Verdict | N | Rows |
|---|---|---|
| CONFIRMED | 8 | MFT intact-C5 origin (1); MFT thiamine lane (4); FFT (5); furfural-from-pentose (7); HDMF hexose intact-C6 (10); DMHF/HEMF pentose+AA (11); retro-aldol fragmentation (13); acrylamide Asn route (17) |
| CONTRADICTED | 4 | norfuraneol as the MFT intermediate (2); HMF-from-3-deoxyosone **for fructose** (9); acrylamide carbonyl-identity neutrality (18); nonanal sourced from the linoleate pool (23) |
| PARTIALLY SUPPORTED | 10 | hexose MFT shortcut (3); furfural-from-hexose (8); HMF for glucose (9); norfuraneol species/cyclisation (12); Strecker aldehydes (14); MeSH/DMDS/DMTS (15); 2,5-DMP (19); hexanal (22); CML/CEL (24); hexanal-amine trapping (28) |
| UNTESTED | 3 | FFT thiohemiacetal sub-route (6); 2-alkylthiazole condensation (21); the `[HH]` token and disulfide dimerisation (27) |
| CORRECT OMISSION | 1 | DMS absent because SMM is not a precursor (16) |
| COMPLETENESS GAP | 3 | other dimethylpyrazines (20); 2-acetyl-1-pyrroline (25); cysteine-S-conjugate sink (26) |

**Top 3 topology risks**

1. **The flagship MFT route is contradicted by two independent laboratories, and the model's
   only surviving sulfur-branch anchor is one of them.** Cerny & Davidek's spiking experiment
   says norfuraneol is "unimportant as an intermediate"; Hofmann & Schieberle 1998 — the paper
   `refit_sulfur_barriers_hofmann.py` fits against — calls norfuraneol/cysteine "less
   effective". The evidenced intermediate, 1,4-dideoxypento-2,3-diulose, is absent from the
   network, and adopting it would *remove* a reducing-equivalent token rather than add one.
   Any re-run of the sulfur refit before the route change will fit the wrong route more
   precisely.
2. **Nonanal is generated from the wrong fatty acid.** The whole hydroperoxide pool is scaled
   by `linoleic_acid_pct`, 12-15 % of it is emitted as nonanal, and nonanal is not a product of
   either linoleate hydroperoxide isomer; `oleic_acid_pct` is dead code. This is a
   substrate-to-product mis-assignment in one of the four shipped lipid markers, in the lane
   the hold-out anchors sit on.
3. **The `[HH]` reducing-equivalent gate that red-team finding H4 exposed is only half fixed —
   measured.** Wave I gave MFT a cysteine-derived donor, but the *hexose furanone* step
   (`reaction_templates.py:343`, `elif can == hexose_do1 and h2 is not None`) is still
   pool-gated on `[HH]`, and in a cysteine-free system the only `[HH]` producer that runs
   before it is `Aminoketone_Condensation` (engine ordering: 3b-1 needs cysteine, 3c is the
   pyrazine step, 3e-1 is the furanone route). Measured on glucose + glycine, pH 5.5, 150 °C,
   `max_generations=3`, by monkey-patching `_aminoketone_condensation` to return no steps
   (in-process; no repo file touched):

   ```
   [baseline]                    16 steps | DMHF steps 1 | H2 producers: ['Aminoketone_Condensation']
   [aminoketone cond. DISABLED]  14 steps | DMHF steps 0 | H2 producers: []
   ```

   So **predicted furaneol from glucose is still contingent on pyrazine chemistry** — the exact
   failure mode Wave I fixed for MFT, surviving one lane over. Wang & Ho's CAMOLA result makes
   this worse rather than better: it establishes the hexose → DMHF route as real and
   intact-skeleton, so the model is gating a *confirmed* route behind an unphysical dependency.
   The same donor fix already applied to MFT (a real reductone couple, or the cysteine/cystine
   token where cysteine exists) is not available in a cysteine-free system, so this one needs
   the sugar-derived enediol/dehydro-reductone couple the Wave I note declined to invent — or
   the reduction has to be lumped into the cyclisation step's barrier and the token dropped.

**Top 3 completeness gaps**

1. **2-Acetyl-1-pyrroline and the proline/ornithine lane** — CAMOLA-elucidated *in extrusion*,
   at this repo's own process conditions, and entirely absent (no precursors, no route).
2. **Lipid-derived Strecker initiation** — 4-oxo-2-alkenals, 2,4-alkadienals and
   4,5-epoxy-2-alkenals raise Strecker-aldehyde yields 300-900 % with measured Ea, and none of
   them are species in the network; 2,4-decadienal, the principal 9-HpODE product, is missing
   too. The model's only lipid×Maillard family produces thiazoles instead.
3. **The low-moisture C2+C3 recombination route to MFT** (hydroxyacetaldehyde +
   mercapto-2-propanone, Hofmann & Schieberle's highest-yielding system at 1.4 mol %) — likely
   the dominant channel under roasting/extrusion conditions, and inexpressible in the current
   network. Runner-up: ~50 % of CML flux (the Amadori-oxidation route) is unmodelled while the
   CML predictor sits 1204x low.

**One thing this exercise did *not* find:** no evidence contradicts the core sugar trunk
(Schiff base → Amadori/Heyns → 1-/3-deoxyosone → retro-aldol fragments). Yaylayan's labelled
glucose work independently reproduces every fragment species the model emits from that trunk,
with the intact carbon chains the model assigns them. The Wave-G1 chemistry review's verdict
that the trunk is sound survives isotope scrutiny.

---

## Sources

Listed with the access level actually achieved. All DOIs verified against Crossref during this
review; none is asserted from memory.

| Source | DOI | Access |
|---|---|---|
| Cerny C., Davidek T. (2003) *Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction.* JAFC 51:2714-2721 | 10.1021/jf026123f | `[ABS]` PubMed 12696962 |
| Cerny C., Davidek T. (2004) *α-Mercaptoketone Formation during the Maillard Reaction of Cysteine and [1-13C]Ribose.* JAFC 52:958-961 | 10.1021/jf035265m | `[ABS]` PubMed 14969557 |
| Cerny C., Briffod M. (2007) *Effect of pH on the Maillard Reaction of [13C5]Xylose, Cysteine, and Thiamin.* JAFC 55:1552-1556 | 10.1021/jf062874w | `[ABS]` PubMed 17243706 |
| Cerny C. (2007) *Origin of carbons in sulfur-containing aroma compounds from the Maillard reaction of xylose, cysteine and thiamine.* LWT 40:1309-1315 | 10.1016/j.lwt.2006.09.008 | `[META]` (publisher abstract not retrievable; not used for any verdict) |
| Cerny C., Guntz-Dubini R. (2008) *Identification of 5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine, Cysteine, and Xylose.* JAFC 56:10679-10682 | 10.1021/jf801762c | `[META]` |
| Cerny C., Guntz-Dubini R. (2013) *Formation of cysteine-S-conjugates in the Maillard reaction of cysteine and xylose.* Food Chem 141:1078-1086 | 10.1016/j.foodchem.2013.04.043 | `[ABS]` PubMed 23790889 |
| Hofmann T., Schieberle P. (1998) *Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol.* JAFC 46:235-241 | 10.1021/jf9705983 | `[ABS]` PubMed 10554225 |
| Van den Ouweland G.A.M., Peer H.G. (1975) *Components contributing to beef flavor. Volatile compounds produced by the reaction of 4-hydroxy-5-methyl-3(2H)-furanone and its thio analog with hydrogen sulfide.* JAFC 23:501-505 | 10.1021/jf60199a045 | `[META]` (title/authors/pages verified; the Wave I DOI repair is **correct**) |
| Bolton T.A., Reineccius G.A., Liardon R., Ba T. (1994) *Role of Cysteine in the Formation of 2-Methyl-3-furanthiol in a Thiamine—Cysteine Model System.* ACS Symp. Ser. 543:270-278 | 10.1021/bk-1994-0543.ch022 | `[META]`; its 34S percentages `[SECONDARY]` — must be verified against the chapter |
| Blank I., Fay L.B. (1996) *Formation of 4-Hydroxy-2,5-dimethyl-3(2H)-furanone and 4-Hydroxy-2(or 5)-ethyl-5(or 2)-methyl-3(2H)-furanone through Maillard Reaction Based on Pentose Sugars.* JAFC 44:531-536 | 10.1021/jf950439o | `[ABS]` via OpenAlex |
| Wang Y., Ho C.-T. (2008) *Formation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone through Methylglyoxal: A Maillard Reaction Intermediate.* JAFC 56:7405-7409 | 10.1021/jf8012025 | `[ABS]` PubMed 18593173 |
| Yaylayan V.A., Keyhani A. (2000) *Origin of Carbohydrate Degradation Products in L-Alanine/D-[13C]Glucose Model Systems.* JAFC 48:2415-2419 | 10.1021/jf000004n | `[ABS]` PubMed 10888560 |
| Yaylayan V.A., Keyhani A. (2001) *Carbohydrate and Amino Acid Degradation Pathways in L-Methionine/D-[13C]Glucose Model Systems.* JAFC 49:800-803 | 10.1021/jf000986w | `[ABS]` PubMed 11262032 |
| Yaylayan V.A. et al. (2003) *Thermal Decomposition of Specifically Phosphorylated D-Glucoses and Their Role in the Control of the Maillard Reaction.* JAFC 51:3358-3366 | 10.1021/jf034037p | `[ABS]` PubMed 12744667 |
| Perez Locas C., Yaylayan V.A. (2008) *Isotope labeling studies on the formation of 5-(hydroxymethyl)-2-furaldehyde (HMF) from sucrose by pyrolysis-GC/MS.* JAFC | 10.1021/jf8010245 | `[ABS]` PubMed 18611024 |
| Zyzak D.V. et al. (2003) *Acrylamide Formation Mechanism in Heated Foods.* JAFC 51:4782-4787 | 10.1021/jf034180i | `[ABS]` PubMed 14705913 |
| Stadler R.H. et al. (2004) *In-Depth Mechanistic Study on the Formation of Acrylamide and Other Vinylogous Compounds by the Maillard Reaction.* JAFC 52:5550-5558 | 10.1021/jf0495486 | `[ABS]` PubMed 15315399 |
| Davidek T. et al. (2013) *Study To Elucidate Formation Pathways of Selected Roast-Smelling Odorants upon Extrusion Cooking.* JAFC 61:10215-10219 | 10.1021/jf4004237 | `[ABS]` PubMed 23621440 |
| Schieberle P. (2005) *The Carbon Module Labeling (CAMOLA) Technique.* Ann. N.Y. Acad. Sci. 1043:236-248 | 10.1196/annals.1333.029 | `[ABS]` PubMed 16037244 (method review; no compound-specific data in the abstract) |
| Spreng S. et al. (2021) *Discovery of Polyhydroxyalkyl Pyrazine Generation upon Coffee Roasting by In-Bean Labeling Experiments.* JAFC 69:6636-6649 | 10.1021/acs.jafc.1c01894 | `[ABS]` PubMed 34097401 |
| Zhang L. et al. (2020) *An Alkylpyrazine Synthesis Mechanism Involving L-Threonine-3-Dehydrogenase...* Appl. Environ. Microbiol. | 10.1128/AEM.01807-19 | `[ABS]` PubMed 31585995 |
| Zhang L. et al. (2021) *The Biosynthesis Mechanism Involving 2,3-Pentanedione and Aminoacetone Describes the Production of 2-Ethyl-3,5-dimethylpyrazine...* JAFC | 10.1021/acs.jafc.9b07809 | `[ABS]` PubMed 32065523 |
| Molla G. et al. (2015) *Aminoacetone oxidase from Streptococcus oligofermentans...* Biochem. J. | 10.1042/BJ20140972 | `[ABS]` PubMed 25269103 |
| Zamora R., Hidalgo F.J. (2012) *Chemical conversion of phenylethylamine into phenylacetaldehyde by carbonyl-amine reactions in model systems.* JAFC | 10.1021/jf301258s | `[ABS]` PubMed 22578256 |
| Zamora R., Alcón E., Hidalgo F.J. (2013) *Strecker-Type Degradation of Phenylalanine Initiated by 4-Oxo-2-alkenals...* JAFC 61:10231-10237 | 10.1021/jf305007y | `[ABS]` PubMed 23360317 |
| Zamora R., León M.M., Hidalgo F.J. (2015) *Oxidative versus Non-oxidative Decarboxylation of Amino Acids...* JAFC 63:8037-8043 | 10.1021/acs.jafc.5b02619 | `[ABS]` PubMed 26189462 |
| Pan X. et al. (2021) *Isolation and identification of putative precursors of the volatile sulfur compounds ... in heat-sterilized melon juices.* Food Chem 343:128459 | 10.1016/j.foodchem.2020.128459 | `[ABS]` PubMed 33158672 |
| Cheng Y. et al. (2020) *Methanethiol, an Off-Flavor Produced from the Thermal Treatment of Mandarin Juices.* JAFC 68:1030-1037 | 10.1021/acs.jafc.9b06647 | `[ABS]` PubMed 31903752 |
| Scherb J. et al. (2009) *Quantitation of S-Methylmethionine in Raw Vegetables and Green Malt...* JAFC 57:9091-9096 | 10.1021/jf901765q | `[ABS]` PubMed 19754146 |
| Miyazaki R. et al. (2023) *Elucidation of decomposition pathways of linoleic acid hydroperoxide isomers by GC-MS and LC-MS/MS.* Biosci. Biotechnol. Biochem. 87:179-190 | 10.1093/bbb/zbac189 | `[FULL]` |
| Glomb M.A., Monnier V.M. (1995) *Mechanism of Protein Modification by Glyoxal and Glycolaldehyde...* J. Biol. Chem. 270:10017-10026 | 10.1074/jbc.270.17.10017 | `[ABS]` PubMed 7730303 |
| Reynolds J.C. et al. (2007) *Structural analysis of oligomeric molecules formed from the reaction products of oleic acid ozonolysis.* Environ. Sci. Technol. | 10.1021/es060942p | `[ABS]` PubMed 17144295 (cited only for the C9-fragment identity) |
| Hearn J.D., Smith G.D. (2007) *Products and mechanisms of the reaction of oleic acid with ozone and nitrate radical.* J. Phys. Chem. A | 10.1021/jp0500900 | `[ABS]` PubMed 16833788 (same; 1-nonanal at 30 ± 3 % carbon yield) |

### Reproducing the two measurements in this document

Both were run in-process against the working tree at Wave L1 time, with no repo file modified.

* Route inventory: `SmirksEngine(ReactionConditions(pH=5.5, temperature_celsius=150))
  .enumerate([sugar, amino_acid], max_generations=3)`, then group each step's products by
  `reaction_family`.
* Furaneol/pyrazine coupling: same call for D-glucose + glycine, run twice, the second time
  with `src.smirks_engine._aminoketone_condensation` monkey-patched to return `[]`.

### Methodological caveats for the reviewer

* **Coverage is uneven by design of the literature, not by choice.** JAFC before 1996 and the
  ACS Symposium Series are not indexed in PubMed, so the two oldest and most-cited sources here
  (van den Ouweland & Peer 1975; Bolton et al. 1994) could be verified bibliographically but not
  read. Any claim resting on them is labelled `[META]` or `[SECONDARY]` and no verdict depends
  on one.
* **Abstract-level evidence is unusually strong in this field.** ACS abstracts in flavour
  chemistry state the labelling outcome explicitly ("mainly 13C5-labeled", "1:1 mixture of
  [13C6]DMHF and [12C6]DMHF"). Where a verdict turns on a number rather than a direction, the
  full text should still be obtained before the finding is acted on.
* **Two verdicts would change with better access**: the hexanal β-scission mechanism (one
  full-text study against a large classical literature not retrieved here), and the Bolton 34S
  percentages (which, if confirmed, re-diagnose the panel's worst quantitative failure).
