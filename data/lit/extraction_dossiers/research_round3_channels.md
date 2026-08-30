# Round-3 literature dossier — the four DECLARED MISSING CHANNELS

**Date:** 2026-08-29. **Agent:** literature-research (round 3).
**Scope:** HMF · DMHF/furaneol · furfural without cysteine · lipid-oxidation → hexanal.
**Nothing in the repository was written or modified.** All downloads landed in this scratchpad.

## HOW TO READ THIS FILE — verification codes

Every claim carries one of these. **Nothing here is a guess.** Where I could not read a thing,
the row says so and is not counted as a finding.

| code | meaning |
|---|---|
| **[FT-local]** | I read the **full text of a PDF held on this machine** (repo `data/articles/` or this scratchpad) and extracted the number from the text layer myself |
| **[FT-web]** | I read the **open-access full text on the web** (PMC / HAL / publisher) and extracted the number from the article body or an HTML table |
| **[ABS]** | I read only the **abstract** (Europe PMC `resultType=core` or Crossref). Numbers, if any, are abstract-level only |
| **[META]** | Crossref-verified **bibliographic metadata only** (DOI, exact title, journal, volume, pages, authors). No content read |
| **[NEG]** | A **verified negative** — I read the source and it does NOT contain what was hoped for |

**Every DOI and every exact title in the DOWNLOAD LIST at the bottom was resolved through the
Crossref API (`api.crossref.org/works/<doi>`) in this session.** None is reconstructed from memory.

---

# §0. THE FIVE THINGS THE BUILD WAVE MUST NOT GET WRONG

**0a. ★ THE BRIEF'S PREMISE FOR CHANNEL 1 IS CONTRADICTED BY THE STRONGEST KINETIC LITERATURE.**
The brief asks for "HMF formation from 3-DG (dehydration)". **Four independent multiresponse
kinetic studies from the Gökmen school, plus the AgroParisTech cake series, place the dominant
HMF flux on the FRUCTOSE / fructofuranosyl-cation dehydration limb, not on 3-DG** — and place
**FURFURAL, not HMF, downstream of 3-deoxyglucosone.** See §A.2. Verbatim, from an open-access
full text I read:

> "5–HMF and F can form from fructose and 3,4–deoxyosone or 3–dexyosone, respectively. The
> pathways are independent of Leu, but the amounts of these precursors vary according to the
> presence of the amino acid" — Lee et al. 2024, Food Res. Int. 183:114183, §3.2 **[FT-web]**

and

> "5-hydroxymethylfurfural formation mainly proceeded via fructofuranosyl cation dehydration
> **rather than 3-deoxglucosone**" — Göncüoğlu Taş & Gökmen 2017 abstract **[ABS]**

This is the same finding the repo's own audit already recorded from Perez Locas & Yaylayan 2008
(`tasks/audit_remediation.md` §(g) ITEM 5). **The 3-DG → HMF edge is a minority channel at best;
building the HMF node on it will reproduce the SUG-12 loss.** The one dissent is Şen & Gökmen
2022 (nuts/seeds), where "5-hydroxymethylfurfural formation **from 3-deoxyglucosone**" was
"quantitatively important" — a sucrose-rich, low-moisture, restricted-reducing-sugar matrix
**[ABS]**. Treat the branch fraction as matrix-dependent, exactly as §0 of `k3` already says
about branch fractions generally.

**0b. ★ THE VAN BOEKEL / MARTINS / BRANDS / KNOL CORPUS CONTAINS NO HMF CONSTANTS, AND THIS IS NOW
VERIFIED BY READING THE PDFs WE ALREADY HOLD.** The brief asked me to check before proposing
downloads. I checked all of them. **They deliberately EXCLUDED HMF because it was negligible at
their pH.** Full receipts in §A.1. Do not order any of these papers for HMF.

**0c. ⚠️ NAMING TRAP — "HMF" MEANS NORFURANEOL IN TWO PAPERS THE REPO ALREADY HOLDS.**
`whitfield1999.pdf` (48 hits) and `whitfield2001.pdf` (65 hits) use **"HMF" = 4-hydroxy-5-methyl-
3(2H)-furanone = NORFURANEOL**, not hydroxymethylfurfural. Verified from their titles
**[FT-local]**: *"Reaction of 4-hydroxy-5-methyl-3(2H)-furanone (HMF) with cysteine or hydrogen
sulfide at pH 4.5"*. Any grep-driven ingestion of "HMF" across the corpus will silently merge the
furanone lane into the furanic lane. `k3` §D already flags the Whitfield 2001 "Methods omits HMF"
line as "likely a typo" — **it is not a typo; it is this abbreviation collision.**

**0d. DMHF/FURANEOL HAS NO RATE CONSTANTS IN THE ACCESSIBLE LITERATURE. AT ALL.** I ran three
independent Europe PMC queries and four web searches. There is no Arrhenius fit, no rate constant
and no activation energy for HDMF formation or degradation in any Maillard model system I could
reach. What exists is **quantified yields at single temperatures** plus **two independent
isotope-labelling structural constraints** that are strong enough to build the node's topology.
See §B. **10.1021/jf950439o is already in the repo (`data/articles/blank1996.pdf`) and quantifies
nothing** — receipts in §B.1.

**0e. LIPID → HEXANAL AT COOKING TEMPERATURE DOES NOT EXIST AS MEASURED KINETICS.** The single
best source has an explicit `k_hexanal`, is open access, is in a **protein-stabilised aqueous
emulsion at pH 6.7** — and is at **25 °C**, with parameters fitted "based on visual agreement of
fit", not regression. It also hands you an author-stated temperature-extrapolation factor. That
combination is either a gift or a trap depending on how it is labelled. See §D.

---

# §A. CHANNEL 1 — HMF

## A.1 — What the corpus already holds: five verified negatives

I read every candidate PDF already on this machine. **Every one is a NEGATIVE for HMF constants.**

### A.1.1 Martins & Van Boekel 2005, Food Chem. 90:257–269 (`data/articles/martins2005.pdf`) — **[FT-local] [NEG]**

This is the flagship glucose/glycine multiresponse network. It **names HMF in the abstract** and
then **deletes it from the model**. Verbatim, §3.3 "Building a reaction network", p. 262:

> "Under the pH conditions used in the present study formation of HMF was consequently very low
> (0–20 µmol L⁻¹). **This reaction pathway was therefore not considered in Scheme 1.**"

Scheme 1 has 15 numbered steps (I read the ASCII rendering of the scheme in the text layer);
**there is no HMF node and no k for one.** The paper also states the 3-deoxy-2-hexosulose "can
cyclise to form HMF" and that "This 3-deoxy-2-hexosulose route is favoured under slightly acidic
conditions" — i.e. the authors themselves say the 3-DG→HMF route needs acid, and their pH 6.8
system does not supply it. **This is a directly usable directional constraint** and it is the
mechanistic reason 0a is true.

### A.1.2 Martins & Van Boekel 2005b, Food Chem. 92:437–448 (the pH ladder) — **[FT-local] [NEG]**

Already dossiered as `martins2005b_extraction.md`. Re-confirmed by grep of the PDF:

> "HMF | '**one order of magnitude lower … µmolar range**'; increases with **decreasing** pH |
> µmol/l | all pHs | §3.1, p. 440 | [M] — **no numeric value published**"
> "**HMF** | **not usable — no numbers published** … Measured, judged negligible … and discarded
> without publication. Only the *direction* (increases with decreasing pH) is quotable."

### A.1.3 Martins, Marcelis & van Boekel 2003 Part I, Carbohydr. Res. 338:1651–1663 — **[FT-local] [NEG]**

Existing dossier §13, quoting the paper: *"HMF was not formed at pH 6.8 and only µmole amounts at
pH 5.5"*; *"quantification was not possible since no reference material was available"*.
**HHMF / DDMP / HMF → "not usable — no numbers published".**

### A.1.4 Brands & van Boekel 2001, JAFC 49:4667–4675 (`data/articles/brands2001.pdf`) — **[FT-local] [NEG]**

The sugar-casein network paper, 120 °C, pH 6.7. Verbatim, p. 4670:

> "Other identified compounds were HMF, furfuryl alcohol, HHMF, and DDMP. **However, HMF and
> furfuryl alcohol were formed in very low amounts (0-40 µM).** HHMF and DDMP could not be
> quantified because no reference material was available … These compounds were presumably also
> formed in low amounts (as judged using the response factor of HMF)."

HMF is not in the reaction network. **The 0–40 µM ceiling at 120 °C / pH 6.7 is itself a
quotable [M] bound.**

### A.1.5 Knol 2005 / 2009 / 2010 (`data/articles/knol*.pdf`) — **[FT-local] [NEG]**

Case-insensitive grep for `HMF|hydroxymethyl` across all three full text layers returns
**zero hits in all three**. The glucose–asparagine / fructose–asparagine acrylamide networks do
not contain an HMF node. Do not re-open them for this.

**Corpus-wide sweep [FT-local]:** I grepped all 101 PDFs in `data/articles/` for
`hydroxymethylfurfural|HMF` with a threshold of >3 hits. Only six files pass:
`Wang2025`(4), `brands2001`(4), `martins2003b`(7), `martins2005`(4) — all covered above — and
`whitfield1999`(48) / `whitfield2001`(65), which are the **norfuraneol naming collision of §0c**.
**The corpus's HMF gap is total.**

## A.2 — The real HMF kinetics literature: the Gökmen multiresponse school

This is one lab (Hacettepe FOQUS), one method (multiresponse Bayesian/least-squares fitting of a
full elementary-step network to 8–12 simultaneously measured responses), applied across six
matrices. **It is the only body of work that fits a genuine HMF node inside a Maillard network.**
None of the six is open access; all six abstracts read **[ABS]**, all six DOIs Crossref-verified
**[META]**.

| # | paper | system | T | what the abstract states about the HMF node |
|---|---|---|---|---|
| 1 | **Kocadağlı & Gökmen 2016**, JAFC 64:6333–6342, `10.1021/acs.jafc.6b01862` | **glucose ALONE, caramelization, no amine** ± NaCl | 180, 200 °C | *"HMF formation was revealed to be mainly **via isomerization to fructose and dehydration over cyclic intermediates**, and the rate constants increase 4-fold in the presence of NaCl."* Also: *"A decrease in rate constants of 3-deoxyglucosone and 1-deoxyglucosone formations by the presence of NaCl"*; *"Interconversion between glucose and fructose became 2.5 times faster in the presence of NaCl"* |
| 2 | **Kocadağlı & Gökmen 2016**, Food Chem. 211:892–902, `10.1016/j.foodchem.2016.05.150` | glucose/wheat flour, low moisture | not in abstract | *"**Formation of 5-hydroxymethyl-2-furfural from fructose was found to be a key step.**"* Also: *"Formation of 3-deoxyglucosone proceeded directly from glucose and also Amadori product degradation."* Ten responses incl. glucosone, 1-DG, 3-DG, 3,4-DGE, HMF, GO, MGO, diacetyl |
| 3 | **Göncüoğlu Taş & Gökmen 2017**, Food Chem. 221:1911–1922, `10.1016/j.foodchem.2016.11.159` | hazelnut roasting | 150/160/170 °C | *"5-hydroxymethylfurfural formation mainly proceeded **via fructofuranosyl cation dehydration rather than 3-deoxglucosone**"*. ⚠️ **and**: *"The temperature dependence of the reactions was complicated and **could not be explained by the Arrhenius equation**."* — a pre-registered warning that this matrix yields no Ea |
| 4 | **Hamzalıoğlu & Gökmen 2020**, Food Chem. 318:126467, `10.1016/j.foodchem.2020.126467` | coffee roasting | 200/220/240 °C | *"fructofuranosyl cation contributed mostly to the formation of 5-hydroxymethylfurfural which was found to be the most important intermediate precursor of acrylamide"* |
| 5 | **Gürsul Aktağ & Gökmen 2020**, Food Chem. 320:126620, `10.1016/j.foodchem.2020.126620` | apple juice, orange juice, peach nectar, **storage** | low T, acidic | *"The contribution of fructose dehydration through fructofuranosyl cation on the formation of 5-hydroxymethylfurfural was **significantly higher (p < 0.05) than 3-deoxyglucosone dehydration**."* **The only acidic, food-pH, low-temperature member of the set** |
| 6 | **Şen & Gökmen 2022**, Food Chem. 395:133583, `10.1016/j.foodchem.2022.133583` | 5 nuts/seeds, sucrose-rich low moisture | 160/180 °C | ⚠️ **the dissent**: *"3-deoxyglucosone formation via sugar degradation; **5-hydroxymethylfurfural formation from 3-deoxyglucosone** … were found to be quantitatively important"* |

Plus **Berk, Hamzalıoğlu & Gökmen 2020**, Eur. Food Res. Technol. 246:2399–2410,
`10.1007/s00217-020-03583-z`, sesame at 180/200/220 °C: *"5-hydroxymethylfurfural formation was
mostly originated from fructofuranosyl cation"* **[ABS]**.

**Verdict.** For a *3-DG → HMF* edge specifically, **#1 is the single highest-value download**: it
is the only member with **no amine at all**, so the fitted constants are pure sugar chemistry and
cannot be confounded by the amino-acid lane; and its abstract is the only one that states the
HMF rate constants' *sensitivity* (4× on NaCl). **#5 is the best pH match to the repo** (acidic
juice, storage temperatures) and the only one where the fructose-vs-3-DG partition was tested for
statistical significance. **#3's abstract pre-declares that no Arrhenius Ea is recoverable from
it** — that is a reason to rank it lower, and to expect the same problem elsewhere in the set.

## A.3 — HMF constants I could actually READ, and what they are worth

### A.3.1 ★ Han et al. 2025, *Foods* 14:2136, `10.3390/foods14122136` — **[FT-web], OA, numbers extracted** — **VERDICT: STRUCTURE = USE, k VALUES = REFUSE**

Membrane-clarified sugarcane juice + a synthetic model system (sucrose, glucose, fructose, lysine,
histidine, proline in 0.07 M Na-acetate / 0.03 M acetic acid, 15 °Bx), **pH 6.50 (model) / 6.42
(juice), 60 / 70 / 80 / 90 °C, vacuum 0.08 MPa.** This is the closest temperature and pH match to
the repo's window of anything I found with a published HMF ODE. I pulled Table 5 and the appendix
equations out of the PMC HTML myself.

**★ The ODE is directly transplantable — this is the paper's real value:**

> `d[5-HMF]/dt = k7·[Fru] + k8·[3-DG] − k15·[5-HMF]`   (their Eq. A6)
> `d[melanoidin]/dt = k12·[CML] + k14·[CEL] + k15·[5-HMF]`   (their Eq. A11)

i.e. **HMF has two parallel sources (fructose AND 3-DG) and exactly one sink (melanoidins)**, and
the sink is first order in HMF. Fifteen elementary steps; k7 = Fru→5-HMF, k8 = 3-DG→5-HMF,
k15 = 5-HMF→melanoidin.

**⚠️ THE FITTED k8 AND k15 ARE UNIDENTIFIABLE AND MUST BE REFUSED.** Verbatim from their Table 5
(k in h⁻¹, ± is the 95 % HPD half-width):

| step | 60 °C | 70 °C | 80 °C | 90 °C |
|---|---|---|---|---|
| k7 Fru → 5-HMF | 0.0100 ± 0.0001 | 0.0100 ± 0.0001 | 0.0149 ± 0.0149 | 0.0100 ± 0.0001 |
| **k8 3-DG → 5-HMF** | **38,991 ± 37,170** | **20,453 ± 10,810** | **27,521 ± 21,480** | **28,918 ± 22,910** |
| **k15 5-HMF → melanoidin** | **31,945 ± 26,290** | **23,069 ± 12,100** | **31,232 ± 24,990** | **28,874 ± 19,900** |
| k6 1,2-enediol → 3-DG | 0.0445 ± 0.0375 | 0.0343 ± 0.0196 | 0.5366 ± 0.3574 | 0.0506 ± 0.0252 |
| k10 3-DG → MGO | 46.4663 ± 29.7000 | 46.6424 ± 21.4401 | 51.9015 ± 29.9000 | 56.3724 ± 35.2000 |

Three disqualifying defects, all visible in the table itself: **(i)** the HPD half-width is
50–95 % of the estimate on k8 and k15 — the same SD ≥ estimate pathology the repo already refused
in Knol (`k3` §C.6); **(ii)** k8 is **non-monotonic and actually DECREASES** from 60 → 70 °C
(38,991 → 20,453), so no Arrhenius fit is possible and the temperature signal is noise;
**(iii)** k8 and k15 are ~10⁶× every other constant in the same table, which is the signature of
a **stiff, unidentified pair** — HMF is at quasi-steady state and only the *ratio* k8/k15 is
constrained, not either value. Note k7 is pinned at its bound (0.0100 ± 0.0001) in 3 of 4 columns.
**k12 and k14 show the same bound-pinning.** This is a fit that did not converge on these nodes.

**What IS ingestable [M/F]:** their separate "simple kinetic model" (Table 4, zero/first/second
order fits to each marker independently), which gives
**Ea(5-HMF) = 18.836 ± 5.755 kJ/mol** and **Ea(3-DG) = 25.691 ± 3.081 kJ/mol**, 60–90 °C, pH 6.5,
4 points. R² for the 5-HMF first-order rows: 0.991 / 0.980 / 0.987 / 0.990. **These are LUMPED
apparent accumulation Ea, exactly the same species as the Bornhorst norfuraneol Ea the repo
already carries** — and they must carry the same label discipline: *"apparent lumped
approach-to-accumulation Ea for 5-HMF in a vacuum-evaporating acetate-buffered sugarcane-juice
model, pH 6.5, 60–90 °C, 4 points"*, never "the HMF barrier". ⚠️ Their same Table 4 prints
**negative** k for every second-order row (e.g. −350.2164 for 5-HMF), which is nonsense and shows
the table is a mechanical three-order sweep, not a curated fit. Ingest the Ea with that caveat
attached or not at all.

### A.3.2 Ma et al. 2022, *Front. Nutr.* 9:940202, `10.3389/fnut.2022.940202` — **[FT-web], OA, numbers extracted** — **VERDICT: REFUSE the Ea, keep the 3-DG row as a weak prior**

Glucose (60 mg) + asparagine (50 mg) + linoleic acid (16 mg) in 2 mL 0.1 M phosphate buffer
**pH 7.4**, sealed tubes, **160 / 180 / 200 °C**, 1–15 min. I read Tables 1–3 including the R²
column off the Frontiers HTML.

| Table 3 | k(160) | k(180) | k(200) | R² | **Ea (kJ/mol)** |
|---|---|---|---|---|---|
| acrylamide | 0.268 | 0.282 | 0.364 | 0.980 / 0.944 / 0.926 | **12.87** |
| **5-HMF** | 0.285 | 0.305 | 0.405 | 0.996 / 0.974 / 0.999 | **14.85** |

| Table 2 | k(160) | k(180) | k(200) | R² | **Ea (kJ/mol)** |
|---|---|---|---|---|---|
| **3-DG** | 0.103 | 0.235 | 0.758 | 0.915 / 0.917 / 0.917 | **84.55** |
| MGO | 0.147 | 0.149 | 0.153 | 0.991 / 0.819 / 0.819 | **1.84** |
| GO | 0.067 | 0.155 | 1.128 | 0.961 / 0.922 / 0.922 | **119.75** |

All k are **first-order FORMATION** constants in min⁻¹ (their Origin "Exp2PMod1"/"Exponential"
fits to each compound's own accumulation curve, three temperatures). **They are NOT network
constants and there is no reaction scheme behind them.**

**⚠️ WHY THE 5-HMF Ea MUST BE REFUSED: the same table produces Ea(MGO) = 1.84 kJ/mol.** That is
a chemically impossible barrier — it says methylglyoxal formation is temperature-independent over
a 40 °C span — and it comes from three k values (0.147 / 0.149 / 0.153) that are flat because the
compound has **plateaued**, not because the barrier is small. An approach-to-plateau curve fitted
with a pure exponential returns an apparent k that saturates, and its Arrhenius slope is then an
artefact of the plateau, not a barrier. **Ea(5-HMF)=14.85 and Ea(acrylamide)=12.87 sit in exactly
that regime** (k moves only 1.4× over 40 °C). Ea(3-DG)=84.55 and Ea(GO)=119.75 come from k values
that move 7× and 17× and are the only two rows in the paper not obviously plateau-limited.
**Verdict: Ea(3-DG) = 84.55 kJ/mol is a defensible weak prior [F]; Ea(5-HMF) = 14.85 is REFUSE,
and the reason must be recorded so a later wave does not re-ingest it.** Note independently that
Ma's Ea(5-HMF)=14.85 and Han's Ea(5-HMF)=18.836 agree — but **both are plateau-limited lumps from
different directions, so their agreement is not corroboration.**

### A.3.3 Lee et al. 2022, Food Chem. 376:131917 — **[FT-web via HAL, OA, downloaded]** — **VERDICT: ★ BENCHMARK / HOLD-OUT, not a rate source**

Downloaded from `agroparistech.hal.science/hal-03819621` → `scratchpad/lee2022.pdf`. Solid sponge-
cake model containing **only known precursors: glucose (G) or glucose+leucine (G+L)**, baked at
**140 / 170 / 200 °C**, two fan frequencies, 170 °C in triplicate. **12 markers quantified
absolutely in µmol·gDM⁻¹**: glucose, leucine, fructose, glucosone, 1-DG, 3-DG, 3,4-DGE, GO, MGO,
diacetyl, **furfural, 5-HMF**, plus browning.

**[FT-web] absolute numbers, verbatim from p. 14:**
> "the kinetics of 5-HMF and furfural formation were very similar, but 5-HMF was found in much
> larger quantities than furfural (**20 and 50 µmol.gDM-1, in G and G+L, respectively, for 5-HMF,
> against 2 and 8 µmol.gDM-1, in G and G+L, respectively, for furfural, after 80 min at 200°C**)
> … A short lag phase of approximatively 20 min was seen, similar to that seen for fructose"

⇒ **HMF : furfural = 10 : 1 (G, no amine) and 6.25 : 1 (G+L)**, and **adding leucine multiplies
HMF by 2.5× and furfural by 4×**. Four ratios, all unit-independent, all from one lab / one
method / one matrix — **immune to the response-factor caveat that disqualifies most of §B of
`k3`.** This is exactly the kind of same-method pair `k3` §C.2 complains does not exist.

**⚠️ NO FITTED RATE CONSTANTS EXIST IN THIS PAPER.** I grepped the full text: there is no k table,
no Arrhenius fit, no regression section. The abstract's phrase "compared to experimental kinetic
data" means concentration-versus-time curves, not constants. **Its value is as a hold-out
benchmark, and it is a very good one**: a two-formula (amine on/off) × three-temperature ×
12-marker absolute dataset in a solid matrix, with **the underlying data published under their own
open-data DOIs** (`doi.org/10.15454/IICJV3` for the furanics, `10.15454/HV3UTK` for the
α-dicarbonyls, `10.15454/GUXKYQ` for precursors, `10.15454/QYCBIW` for browning) — so the raw
curves are recoverable without the paper.

### A.3.4 Lee et al. 2024, Food Res. Int. 183:114183 — **[FT-web via HAL, OA, downloaded]** — **VERDICT: ⚠️ TRAP — its "rate constants" are TRANSPORT, not chemistry**

Downloaded from `hal.inrae.fr/hal-05318480` → `scratchpad/lee2024.pdf`. The volatile companion to
A.3.3: same cakes, same three temperatures, **furfural and 5-HMF quantified in BOTH the matrix and
the oven vapour**, with a fitted **reparametrised Arrhenius** `k(T) = k_ref·exp(−Ea/R·(1/T − 1/T_ref))`
for two constants `ka` and `kd`.

**★ THIS IS THE MOST DANGEROUS ITEM IN THIS DOSSIER AND I AM FLAGGING IT BEFORE THE VALUES ARE
EVEN READ.** A wave scanning for "Ea for HMF" will find this paper, see an Arrhenius fit on 5-HMF,
and ingest it. **It is not a chemical barrier.** Verbatim from their §2, defining `ka`:

> "It should be noted that ka finally expresses several contributions: – one linked to the
> **partition between the matrix and the vapor** that could resemble a partition coefficient …;
> – a second linked to the **water evaporation flow rate**; – a third giving the kinetic part"

and on `kd`, verbatim from §3:

> "No perfect explanation could be found regarding the exact mechanisms behind this disappearing
> rate but some hypotheses can be put forward, including **the deposition of molecules on the oven
> walls, some further reactions occurring in the vapor or even gas leaks** that were not considered
> in the air renewal."

⇒ **`ka` is a lumped partition + evaporation + generation term, and `kd` is an oven-wall
deposition / leak term.** Neither is an HMF formation or degradation barrier. The authors also
state `ka` "displayed considerable dependence on the formulae" and attribute that to melanoidins
physically **trapping** the volatile — a matrix-retention effect, not chemistry.
**RECORD AS A PROHIBITED DERIVATION**, in the same register as `k3` §C.1's thiol-Ea prohibition.

⚠️ Second defect: **the numeric k_ref and Ea values live only inside Fig. 4 panel (C), a raster.**
They are not in the PDF text layer. I could not read them and **I am not reporting numbers I did
not read.** If a wave wants them, they require reading the figure image from
`scratchpad/lee2024.pdf` — and per the above, they should not want them.

**What IS usable [FT-web]:** the reaction-scheme sentence quoted in §0a, and the fact that the
furanic pathways are **"independent of Leu"** — i.e. the amino acid changes precursor pools but
opens no new HMF/furfural route. That is a clean structural constraint and it is the amine-on/off
control the corpus otherwise lacks.

### A.3.5 The HMF SINK — the one channel with a dedicated paper — **[ABS] only**

*"Investigation and kinetic evaluation of the reactions of hydroxymethylfurfural with amino and
thiol groups of amino acids"*, Food Chem., ScienceDirect PII `S0308814617312852`.
⚠️ **I did not verify its authors — do not attribute it.** I could **not** verify a DOI for it through Crossref in this session and I am
therefore **NOT** putting it in the download list with a DOI I have not resolved — it is listed at
the bottom as a *lead requiring DOI resolution first*. What I verified is only that the title
exists and that the topic is HMF's reaction with amino and thiol groups **[ABS-adjacent]**.
⚠️ **The repo needs an HMF sink and currently has none** — Han 2025's k15 (§A.3.1) is the only
published one and it is unidentifiable. If this paper contains second-order HMF + amine/thiol
constants it is the highest-value HMF-sink source in existence; **but its identity must be
Crossref-verified before it is ordered.**

---

# §B. CHANNEL 2 — DMHF / FURANEOL

## B.1 — 10.1021/jf950439o, as the brief asked: what it actually quantifies — **[FT-local] [NEG]**

**★ FIRST: THIS PAPER IS ALREADY IN THE REPO.** `data/articles/blank1996.pdf`, 137,885 bytes,
added 2026-08-29. It does **not** yet have an extraction dossier. Verified identity from its own
first page: **Blank, I. & Fay, L. B., JAFC 1996, 44(2), 531–536**, *"Formation of 4-Hydroxy-2,5-
dimethyl-3(2H)-furanone and 4-Hydroxy-2(or 5)-ethyl-5(or 2)-methyl-3(2H)-furanone through Maillard
Reaction Based on Pentose Sugars"*, Nestlé Research Center Lausanne.

**IT QUANTIFIES NOTHING.** The paper has exactly two tables and I read both:

- **Table 1** — "Sensory Contribution of HDMF and HEMF Estimated by GC-Olfactometry". The entire
  data content is **odour intensities on a `+` / `++` / `+++` scale** at the sniffing port, plus
  linear retention indices on OV-1701 and FFAP. Verbatim body: HDMF `++` (pentose/Gly), `+`
  (pentose/Ala), `+` (pentose alone); HEMF 2a `−` / `+++` / `−`; HEMF 2b `−` / `+`(d) / `−`.
  Footnote a: *"GC-O data are presented in terms of odor intensities perceived at the sniffing port
  (+: weak; +++: intense)."*
- **Table 2** — "Mass Spectral Data of Nonlabeled and ¹³C-Labeled HDMF": m/z fragmentation
  patterns only.

**There is no concentration, no yield, no molar conversion and no rate constant anywhere in the
paper.** Conditions, verbatim from Experimental: *"5 mmol of pentose (xylose, ribose, or
arabinose) and 5 mmol of amino acid (glycine or alanine) … in 5 mL of phosphate buffer (0.2 mol/L
Na₂HPO₄, pH 6.0) … heated at 90 °C for 1 h"*; *"During the reaction, the pH dropped to 5.0
(xylose/glycine) and 5.3 (xylose/alanine)."* ⚠️ **the pH labels are initial, not held** — the same
defect `k3` §C.11 records for Zhou 2023.

**The ONE number in the paper**, verbatim from p. 534:
> "preferentially, but not exclusively, produced by incorporation of the labeled carbon of glycine.
> **The remaining 30% might be formed by sugar fragmentation.**"

⇒ **In xylose/glycine, ~70 % of HDMF carries the glycine-derived C1 and ~30 % comes from sugar
fragmentation.** That is a **branch fraction [M]**, and it is the only quantitative thing here.

**★ THE STRUCTURAL YIELD IS LARGE, THOUGH.** Three facts, all [M], all falsifiable:
1. **HEMF is ALANINE-SPECIFIC.** Detected only in pentose/alanine; absent (`−`) from pentose/glycine
   and from pentose alone. This is an **on/off switch**, the strongest kind of constraint.
2. **HDMF forms from pentose ALONE** (`+`), with no amino acid at all — so the DMHF node needs a
   sugar-only route, not just an amine route.
3. The mechanism is **Strecker-aldehyde chain elongation**: formaldehyde (from glycine) → HDMF,
   acetaldehyde (from alanine) → HEMF, onto a C5 skeleton via acetylformoin-type intermediates.
   For a **pentose**, DMHF is C6 and therefore *requires* a one-carbon addition.

**Verdict: USE for topology and the on/off switch. There is no later Blank/Schieberle paper with
rate constants — I searched for one (§B.3) and it does not exist.**

## B.2 — The two papers that DO carry DMHF numbers — **[ABS]**

### B.2.1 ★ Wang & Ho 2008, JAFC 56:7405–7409, `10.1021/jf8012025` — **the best DMHF download**

*"Formation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone through Methylglyoxal: A Maillard Reaction
Intermediate"*. Verbatim abstract **[ABS]**:

> "DMHF was **quantified and verified by HPLC and GC-MS** … The reaction was performed in the
> **0.5 M phosphate buffer** by heating MG with or without either glycine or cysteine at
> **120 °C for 1 h**. MG alone or MG with cysteine could produce **increased level of DMHF with pH
> increased**, whereas MG with glycine had contrary trend. Experiments using a 1:1 mixture of
> [¹³C₆]glucose and [¹²C₆]glucose indicate that in the presence of glycine or cysteine, **glucose
> skeleton kept intact during DMHF formation** since a 1:1 mixture of [¹³C₆]DMHF and [¹²C₆]DMHF
> was formed. **Acetylformoin was detected in the glucose with amino acid reaction system as a
> precursor of DMHF, while in the MG reaction systems, acetylformoin could not be identified.**
> It is suggested different pathways of DMHF formation via MG and glucose."

Why this is the top DMHF pick: **(a)** actual quantification (HPLC + GC-MS, not GC-O intensities);
**(b)** a **pH ladder**, and a *sign-reversing* one — DMHF rises with pH under MG+cysteine and
falls with pH under MG+glycine, which is a two-sided directional test a model can fail in either
direction; **(c)** it carries **cysteine**, so it plugs straight into the repo's sulfur lane;
**(d)** the ¹³C₆/¹²C₆ scrambling test is a hard **structural** result: **DMHF from hexose is NOT
formed by C3+C3 fragment recombination.**

### B.2.2 Poisson et al. 2019, JAFC 67:13829–13839, `10.1021/acs.jafc.9b00770` — **[ABS]**, corroborates B.2.1

*"Generation of α-Diketones and 4-Hydroxy-2,5-dimethyl-3(2H)-furanone upon Coffee Roasting—Impact
of Roast Degree on Reaction Pathways"*. CAMOLA labelled/unlabelled sucrose, biomimetic in-bean.
Verbatim: *"**HDMF was generated from sucrose almost exclusively via cyclization of an intact
skeleton, irrespective of the roast time.**"* — and, as an internal control in the same experiment,
2,3-butanedione and 2,3-pentanedione **do** shift to fragment recombination with roast degree.
**So the intact-skeleton result for HDMF is not an artefact of the method: the method demonstrably
detects recombination when it happens, and finds none for HDMF.**

**★ RECONCILE B.1 WITH B.2.1/B.2.2 — this is the rule the build wave should encode:**
- **hexose → HDMF: intact C6 skeleton, no fragment recombination** (Wang & Ho ¹³C₆; Poisson CAMOLA)
- **pentose → HDMF: REQUIRES a C1 addition**, ~70 % from the Strecker aldehyde of the amino acid,
  ~30 % from sugar fragmentation (Blank & Fay)

These do not conflict — a C5 skeleton *cannot* give C6 intact — and together they fully specify
the node's carbon topology on both sugar classes. **That is the whole channel-2 answer that
literature can currently give.**

## B.3 — The verified negative: no DMHF rate constants exist — **[NEG]**

I ran, and read the returns of:
- Europe PMC: `("dimethyl-3(2H)-furanone" OR HDMF OR furaneol) AND ("rate constant" OR "activation energy" OR "first order")` → **23 hits, zero relevant** (returns are porphyrin electrochemistry, biofilm aerogels, azo dyes, cheese sensory)
- Europe PMC: `furaneol AND kinetic` → 58 hits; the only food-chemistry hit is Poisson 2019 (B.2.2)
- Europe PMC: `acetylformoin` → 10 hits total in the entire database; the relevant ones are Blank & Fay 1996, Wang & Ho 2008, and Hofmann 1998 `10.1021/jf980512l`
- Europe PMC: `"dimethyl-3(2H)-furanone" AND (quantit* OR yield) AND (Maillard OR "model system")` → 98 hits, all aroma-profiling/sensomics surveys, no kinetics
- four web searches on furaneol/HDMF formation and degradation kinetics

**There is no Arrhenius fit, no rate constant, and no activation energy for HDMF anywhere I could
reach. The corpus's own DMHF holdings are consistent with this**: I grepped all 101 repo PDFs for
`dimethyl-3(2H)-furanone|furaneol|HDMF|DMHF` (>2 hits) and got `leksrisompong2010`(53, sensory
thresholds), `shu1988`(22), `Wang2025`(10), `bi2020`(9), `blank1996`(69), `blank2002`(5),
`hofmann1998`(5), `anantharamkrishnan2020b`(4), `hofmann1995`(4), `kumazawa2003`(3) — **[FT-local]**.
Of these, `shu1988.pdf` is **Shu & Ho 1988, JAFC 36:801–803, "Effect of pH on the Volatile
Formation from the Reaction between Cysteine and 2,5-Dimethyl-4-hydroxy-3(2H)-furanone"** — a
**fed-DMHF + cysteine pH series at 2.2 / 5.1 / 7.1 in a Parr bomb**, read from its first page
**[FT-local]**. It is qualitative GC identification only (no response factors stated on the page I
read), **but it is a fed-precursor sulfur experiment on DMHF that the repo already owns and has
not dossiered** — the exact analogue of the Whitfield fed-norfuraneol experiments the repo does
use. Worth a dossier before any new DMHF download.

**Honest bottom line for channel 2: the DMHF node can be built with correct topology and a pH
sign, but it cannot be given a barrier from literature. Any Ea the repo assigns to DMHF is the
repo's own assumption and must be labelled as such**, in the same register as `k3` §C.1.

---

# §C. CHANNEL 3 — FURFURAL WITHOUT CYSTEINE, AND THE ALANINE SYSTEMS

## C.1 — ★ The best single result: furfural comes from 3-DG, and it is the 3-DG product, not HMF

Both open-access cake papers state it and one of them is the only place I found it asserted from a
network fitted to a full marker set:

> "**Furfural was mainly formed via 3-deoxyglucosone.**" — Lee et al. 2022 abstract **[FT-web]**
> "5–HMF and F can form from fructose and 3,4–deoxyosone or **3–dexyosone**, respectively."
> — Lee et al. 2024 §3.2 **[FT-web]**

**This is the structural swap the repo needs**: the edge it wanted to build (3-DG → HMF) is, in the
best-measured systems, **3-DG → FURFURAL**, with HMF hanging off fructose / 3,4-DGE instead.

**Quantitative, no cysteine anywhere in the system** (§A.3.3, **[FT-web]**): furfural =
**2 µmol·gDM⁻¹ (glucose only) and 8 µmol·gDM⁻¹ (glucose + leucine)** after 80 min at 200 °C, with
5-HMF at 20 and 50 in the same runs. **The amine-on/off contrast is a 4× furfural gain from
leucine.** ⚠️ Note this is a **hexose** system — furfural from glucose is a C5 fragment route, not
pentose dehydration. If the repo's furfural feed is pentose-based, these numbers are a
*different* production route and are not commensurable; use the G-vs-G+L **ratio**, not the level.

## C.2 — Pentose + amine, no cysteine: Apriyantono & Ames 1993 — **[ABS]** — **VERDICT: RATIO-ONLY**

**JSFA 61(4):477–484, `10.1002/jsfa.2740610416`**, *"Xylose-lysine model systems: the effect of pH
on the volatile reaction products"* — Crossref-verified **[META]**. Verbatim abstract:

> "Aqueous molal solutions of xylose and lysine (initial pH 4.9) were refluxed either with control
> of the pH at 5.0 or without pH control (final pH 2.6) … **2-furfural alone comprised 522 and 999
> g kg⁻¹ of the volatiles**, respectively, from the systems with final pH values of 5.0 and 2.6.
> Maintaining the pH at 5.0 resulted in a **higher yield and greater numbers of nitrogen-containing
> compounds**"

⚠️ **522 and 999 g/kg are shares OF THE VOLATILE FRACTION, not absolute yields.** They mean
"furfural was 52 % and 99.9 % of the volatiles". **RATIO-ONLY**, and the 999 figure means the
uncontrolled-pH system produced essentially nothing but furfural, which is a statement about the
denominator collapsing as much as about furfural rising. **Its real value is the pH-control
design**: it is a **held-pH vs drifting-pH pair in one lab**, which is exactly the control
`k3` §C.11 says Zhou 2023 lacks. Pentose + lysine, no cysteine, no alanine.

## C.3 — Alanine specifically — **[FT-local]**

The corpus's alanine content for this channel is **Blank & Fay 1996 (§B.1), already in the repo**:
xylose/ribose/arabinose **+ L-alanine**, 90 °C, 1 h, pH 6.0 → 5.3. It gives **HEMF as an
alanine-only product** and HDMF in the alanine system too, but **no furfural data and no numbers**.
**I found no paper reporting furfural yields from a pentose/alanine system.** That is a genuine
gap and I am reporting it as one rather than filling it with a guess.

## C.4 — The acidic-juice route — **[META] only, abstract not retrievable**

**Agcam, E. 2022**, *"A Kinetic Approach to Explain Hydroxymethylfurfural and Furfural Formations
Induced by Maillard, Caramelization, and Ascorbic Acid Degradation Reactions in Fruit Juice-Based
Mediums"*, **Food Analytical Methods 15:1286–1299, `10.1007/s12161-021-02214-x`** — Crossref-
verified title/journal/volume/pages/author **[META]**. ⚠️ **I could NOT read the abstract**: no
Europe PMC record exists, Crossref carries no abstract, and the Springer page 303-redirects to an
IdP. **I therefore know only its title.** The title promises HMF **and** furfural kinetics in an
acidic food medium with three named source reactions separated — which would be the single best
pH match in this whole dossier — **but I have verified nothing about its contents and it must be
ordered on that basis, not on mine.**

## C.5 — Verified negative: the biorefinery furfural literature is not transferable — **[NEG]**

Two web searches for furfural/xylose dehydration kinetics return an entirely different field:
xylose → furfural in **acetic acid, formic acid, sulfuric acid, pTSA, or beta zeolite** at
**pH 0.9–1.7 and 150–200 °C**, with Ea ≈ 44.7–81.8 kJ/mol. **These are specific-acid-catalysis
constants at proton concentrations 4–6 orders of magnitude above the repo's pH 5–7 systems.** The
rate law is explicitly first-order in [H⁺]. **Do not transfer any of these.** I am recording this
so a later wave does not re-derive the same dead end.

---

# §D. CHANNEL 4 — LIPID OXIDATION → HEXANAL

## D.1 — ★ Schroën & Berton-Carabin 2022 — **[FT-local, downloaded and read]** — the only source with an explicit k_hexanal

**Food Research International 160:111621, `10.1016/j.foodres.2022.111621`.** Crossref-verified
**[META]**. **Open access (CC BY-NC-ND)**, downloaded from `edepot.wur.nl/573932`; text at
`scratchpad/schroen2022.txt`. I extracted everything below from the PDF text layer myself.

**System:** rapeseed-oil O/W emulsions, **five emulsifiers — Tween 20, Tween 80, and three
proteins: β-lactoglobulin, bovine serum albumin, β-casein** — in buffer at **pH 6.7**, droplets
1.4–1.8 µm, incubated under rotative agitation at **25 °C**, oxygen-to-lipid ratio strictly
controlled. Oil composition given as C18:1n-9 / C18:2n-6 / C18:3n-3 at **618.5 / 191.7 / 92.2
mg·g⁻¹ oil**. Measured responses: headspace O₂, conjugated dienes, total hydroperoxides,
**propanal and hexanal**.

**The scheme (their Fig. 1 / Eqs. III–VI):** LH →(k1) L• →(k2, +O₂) LOO• →(k3, +LH′) LOOH;
LOOH →(k4) LO• secondary products; LOOH →(k5) LOO• (hydroperoxide-driven re-initiation).

**Table 1, verbatim [FT-local]** (conjugated-diene basis; all "other hydroperoxide" parameters are
0.95× these, per their footnote):

| | β-lactoglobulin | BSA | β-casein | Tween 20 | Tween 80 |
|---|---|---|---|---|---|
| k1,CD (h⁻¹) | 6.5·10⁻⁵ | 2.7·10⁻⁵ | 1.4·10⁻⁵ | – | – |
| k2,CD ((mol/kg oil)·h)⁻¹ | 10 | 10 | 10 | 10 | 10 |
| k3,CD ((mol/kg oil)·h)⁻¹ | 1 | 1 | 1 | 1 | 1 |
| **k4,CD (h⁻¹)** | **6·10⁻³** | **6·10⁻³** | **6·10⁻³** | **6·10⁻³** | **6·10⁻³** |
| k5,CD (h⁻¹) | – | – | – | 3.4·10⁻³ | 2.7·10⁻³ |
| **k_propanal (h⁻¹)** | **3.5·10⁻⁴** | 3.5·10⁻⁴ | 3.5·10⁻⁴ | 3.5·10⁻⁴ | 3.5·10⁻⁴ |
| **k_hexanal (h⁻¹)** | **6·10⁻⁵** | **6·10⁻⁵** | **6·10⁻⁵** | **6·10⁻⁵** | **6·10⁻⁵** |
| initial L• (µmol/kg oil) | 500 | 10 | 100 | 1 | 1 |
| droplet size (µm) | 1.5 | 1.8 | 1.7 | 1.4 | 1.7 |

**★ The four things that make this usable:**
1. **`k_hexanal = 6·10⁻⁵ h⁻¹`, first order in [LOOH]**, and it is **identical across all five
   emulsifiers** — three of them proteins. Verbatim: *"For propanal and hexanal formation, we use
   the same reaction as for the overall reaction for secondary oxidation products (equation VI),
   but use lower albeit constant k-values for all emulsions."*
2. **The branch fraction is stated explicitly**, verbatim: *"The corresponding constants in Table 1
   can be interpreted as a percentage of the total aldehydes formed related to k4 (**7 % and 1.2 %,
   for propanal and hexanal, respectively**)."* ⇒ **hexanal is 1.2 % of the total secondary-product
   flux, and propanal is 5.8× hexanal.** Two unit-independent ratios.
3. **`k4 = 6·10⁻³ h⁻¹` is the LOOH → secondary-products lump**, constant across all five systems —
   the paper's central claim is that only the *initiation* terms (k1, k5) differ between
   emulsifiers, and k1 is 4.6× higher for β-lactoglobulin than for β-casein.
4. **An author-stated temperature-extrapolation rule**, verbatim (×2, §3.2.2 and §4.2):
   *"In general the reaction rate constant would need to be decreased by a factor of 2–3 for every
   10 °C temperature difference"* and *"When using other temperatures, parameters will need to be
   adjusted according to their temperature-sensitivity, which typically is in the order of a factor
   of 2–3 per 10 °C temperature difference."* ⇒ **Q10 = 2–3**, which is Ea ≈ 53–75 kJ/mol over
   298–308 K. **This is the authors' own statement, [M/C], not my derivation.**

**⚠️ THREE CAVEATS, ALL AUTHOR-STATED, ALL DISQUALIFYING IF DROPPED:**
- **25 °C.** This is a storage experiment. The repo's window is 90–180 °C. **Applying the Q10 = 2–3
  rule from 25 °C to 140 °C is 11.5 decades of extrapolation and a factor of 2¹¹·⁵–3¹¹·⁵ ≈
  3·10³–8·10⁵.** The authors licensed adjustment, not extrapolation of that span. If used above
  ~60 °C the constant must be labelled **"extrapolated far beyond the measured range"**, in the
  same register as `k3`'s "Extrapolation of these Ea above 100 °C | NOT LICENSED".
- **The parameters were not regressed.** Verbatim: *"The parameter values (k1-k5) were determined
  based on **visual agreement of fit**."* There are **no standard errors anywhere in the paper**.
  These are hand-tuned values, honestly labelled as such by the authors. ⇒ **[F], but with no
  uncertainty attached, ever.**
- **`k4` is explicitly a lump.** Verbatim: *"k4 is a lumped reaction rate, ultimately leading to
  formation of secondary oxidation products"* and *"k4 corresponds to the formation of all
  secondary oxidation products (Table 1), not just to hexanal and propanal"*. Same species as the
  Yaghmur furfural-sink Ea the repo already carries with a lump label.

**Verdict: ★ USE, as the ONLY measured hexanal-formation constant in a protein-containing aqueous
system at the repo's pH — under the step name "first-order hexanal formation from lipid
hydroperoxide in a protein-stabilised O/W emulsion at 25 °C and pH 6.7, hand-fitted, no SD".
NEVER as "the hexanal barrier". The 1.2 % branch fraction and the propanal:hexanal = 5.8:1 ratio
are the more robust exports.**

## D.2 — Koelsch, Downes & Labuza 1991 — **[FT-local, downloaded and read]** — **VERDICT: DO NOT ORDER**

**J. Food Sci. 56:816–820, `10.1111/j.1365-2621.1991.tb05389.x`**, *"Hexanal Formation via Lipid
Oxidation as a Function of Oxygen Concentration: Measurement and Kinetics"*. Crossref-verified
**[META]**; full text read from an open PDF at `packagingtechnologyandresearch.com`.

Named in searches as the obvious hexanal-kinetics paper. **It is not usable here.** From its own
methods **[FT-local]**: soybean oil in a **cracker matrix**, at **23 ± 2 °C** and **23 ± 1 % RH**,
in the absence of light, at four constant oxygen levels (**1.2, 4.5, 10.0, 15.4 %**). It derives a
**cubic** model for the monomolecular period and an extended model for the accelerated stage, and
its Tables 1 and 2 give rate constants **as a function of oxygen concentration only**.

**Disqualifying: ONE temperature (23 °C) ⇒ no activation energy exists in it; a dry cracker matrix
(23 % RH) ⇒ not aqueous, not protein; and the reported k's are cubic/extended-model constants in
concentration/(time)³ units that are not commensurable with a first-order sink.** I am recording
this as a **[NEG]** specifically to stop a later wave ordering it on the strength of its title.

## D.3 — Verified negative: no 13-HPODE → hexanal rate constant at cooking temperature — **[NEG]**

I searched for measured decomposition kinetics of the linoleate 13-hydroperoxide to hexanal at
60–100 °C in aqueous or protein systems. **It does not exist in a form the repo can use.** What
the searches return instead:
- **α-alkoxyalkyl-hydroperoxides** in aqueous solution, k(288 K) = 5.3·10⁻⁴ s⁻¹, k(298 K) =
  1.2·10⁻³ s⁻¹, k(308 K) = 2.1·10⁻³ s⁻¹, Ea = 12.3 ± 0.6 kcal/mol — **a different compound class
  (atmospheric chemistry), 15–35 °C** ⇒ not transferable
- **13S-HPOD → 13S-HOD in 3.75–5 M KOH**, Ea 15.3–15.6 kcal/mol — **a strongly alkaline hydroxide
  reduction, not a cleavage to hexanal** ⇒ not transferable
- **peroxide formation** in linseed oil, 71 ± 1 kJ/mol — **formation, not decomposition; bulk oil**
- **Frankel-type volatile distributions** at 180 °C — which is exactly the source the repo already
  holds and `k3` §C.9 already declared *"NOT a yield source … no absolute yield and no Ea exist in
  it (one temperature, 180 °C)"*

**The honest statement for the build wave: `k3` §C.9's declaration still stands after this round.
The lipid channel's temperature dependence is the repo's own assumption, and the only measured
alternative is a 25 °C hand-fitted constant (§D.1) plus its authors' Q10 = 2–3 rule.**

---

# §E. RANKED DOWNLOAD LIST

**Every DOI and exact title below was resolved through `api.crossref.org/works/<doi>` in this
session unless the row says otherwise.** Titles are as Crossref returns them.

## TIER 1 — order these

| # | DOI | exact title | why, in one line |
|---|---|---|---|
| **1** | `10.1021/acs.jafc.6b01862` | **Effect of Sodium Chloride on α-Dicarbonyl Compound and 5-Hydroxymethyl-2-furfural Formations from Glucose under Caramelization Conditions: A Multiresponse Kinetic Modeling Approach** — Kocadağlı T, Gökmen V, *J. Agric. Food Chem.* **64**(32):6333–6342, 2016 | **THE HMF PICK. Glucose alone, no amine** ⇒ the fitted HMF constants are pure sugar chemistry and cannot be confounded by the amino-acid lane. Its abstract already states the fructose-route dominance and a 4× sensitivity — the constants themselves are in the tables |
| **2** | `10.1021/jf8012025` | **Formation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone through Methylglyoxal: A Maillard Reaction Intermediate** — Wang Y, Ho C-T, *J. Agric. Food Chem.* **56**:7405–7409, 2008 | **THE DMHF PICK, and the only quantitative one.** HPLC+GC-MS quantification, a **sign-reversing pH ladder** (MG+Cys rises with pH, MG+Gly falls), **cysteine in the system**, and the ¹³C₆ intact-skeleton test |
| **3** | `10.1016/j.foodchem.2020.126620` | **Multiresponse kinetic modelling of α-dicarbonyl compounds formation in fruit juices during storage** — Gürsul Aktağ I, Gökmen V, *Food Chemistry* **320**:126620, 2020 | **The only Gökmen network at FOOD pH and LOW temperature.** Three matrices; the fructose-vs-3-DG partition for HMF is the one place it was tested for statistical significance (p < 0.05) |
| **4** | `10.1016/j.foodchem.2016.05.150` | **Multiresponse kinetic modelling of Maillard reaction and caramelisation in a heated glucose/wheat flour system** — Kocadağlı T, Gökmen V, *Food Chemistry* **211**:892–902, 2016 | Ten simultaneous responses incl. glucosone/1-DG/3-DG/3,4-DGE/HMF/GO/MGO/diacetyl. **The most complete single α-dicarbonyl + furanic network published**, and the paper Lee 2022/2024 benchmark against |

## TIER 2 — order if the HMF node needs more, or the sulfur/furfural lane is being extended

| # | DOI | exact title | why |
|---|---|---|---|
| 5 | `10.1002/jsfa.2740610416` | **Xylose-lysine model systems: The effect of pH on the volatile reaction products** — Apriyantono A, Ames J, *J. Sci. Food Agric.* **61**:477–484, 1993 | Pentose + amine, **no cysteine**; a **held-pH-5.0 vs drifting-to-2.6 pair in one lab**. RATIO-ONLY — the 522/999 g·kg⁻¹ are shares of the volatile fraction |
| 6 | `10.1016/j.foodchem.2016.11.159` | **Maillard reaction and caramelization during hazelnut roasting: A multiresponse kinetic study** — Göncüoğlu Taş N, Gökmen V, *Food Chemistry* **221**:1911–1922, 2017 | The clearest statement of HMF-via-fructofuranosyl-cation. ⚠️ **its own abstract warns the temperature dependence "could not be explained by the Arrhenius equation"** — order for structure, not for Ea |
| 7 | `10.1016/j.foodchem.2022.133583` | **Kinetic modeling of Maillard and caramelization reactions in sucrose-rich and low moisture foods applied for roasted nuts and seeds** — Şen D, Gökmen V, *Food Chemistry* **395**:133583, 2022 | **The dissent** — the one member of the set where HMF-from-3-DG is quantitatively important. Order it to bound the matrix-dependence of the branch fraction, not to settle it |
| 8 | `10.1021/acs.jafc.9b00770` | **Generation of α-Diketones and 4-Hydroxy-2,5-dimethyl-3(2H)-furanone upon Coffee Roasting—Impact of Roast Degree on Reaction Pathways** — Poisson L, Schaerer A, Spreng S, Mestdagh F, *J. Agric. Food Chem.* **67**:13829–13839, 2019 | CAMOLA corroboration of B.2.1's intact-skeleton result, with a **built-in positive control** (the diketones *do* show recombination) |
| 9 | `10.1016/j.foodchem.2020.126467` | **5-Hydroxymethylfurfural accumulation plays a critical role on acrylamide formation in coffee during roasting as confirmed by multiresponse kinetic modelling** — Hamzalıoğlu A, Gökmen V, *Food Chemistry* **318**:126467, 2020 | Only if the **acrylamide** lane wants the HMF→acrylamide edge; 200–240 °C is out of the repo's window |
| 10 | `10.1007/s00217-020-03583-z` | **Multiresponse kinetic modelling of 5-hydroxymethylfurfural and acrylamide formation in sesame (Sesamum indicum L.) seeds during roasting** — Berk E, Hamzalıoğlu A, Gökmen V, *Eur. Food Res. Technol.* **246**:2399–2410, 2020 | Same architecture as #9 in a simpler matrix; lowest marginal value of the Gökmen set |

## TIER 3 — leads, NOT verified enough to order blind

| # | identifier | what I know | what must happen first |
|---|---|---|---|
| L1 | *"Investigation and kinetic evaluation of the reactions of hydroxymethylfurfural with amino and thiol groups of amino acids"*, *Food Chemistry*, ScienceDirect PII **S0308814617312852** | **The only dedicated HMF-SINK kinetics paper I found**, and the repo has no HMF sink at all. Kocadağlı/Gökmen group | ⚠️ **I could not resolve a DOI for it.** Resolve the DOI and read the abstract before ordering. Do not cite it until then |
| L2 | `10.1007/s12161-021-02214-x` — **A Kinetic Approach to Explain Hydroxymethylfurfural and Furfural Formations Induced by Maillard, Caramelization, and Ascorbic Acid Degradation Reactions in Fruit Juice-Based Mediums** — Agcam E, *Food Analytical Methods* **15**:1286–1299, 2022 | Metadata Crossref-verified. Title promises **HMF and furfural kinetics in an acidic food medium with three source reactions separated** — potentially the best pH match in this dossier | ⚠️ **I read NOTHING of its content** — no PMC record, no Crossref abstract, Springer 303-redirects. Order on the title alone, knowing that |

## ALREADY ON DISK — dossier these instead of downloading anything

| file | what it is | status |
|---|---|---|
| `data/articles/blank1996.pdf` | **Blank & Fay 1996, JAFC 44(2):531–536** = `10.1021/jf950439o` | **★ Already in the repo, NO dossier exists.** §B.1 is the extraction. Quantifies nothing; supplies the HEMF/alanine on-off switch and the 70/30 branch fraction |
| `data/articles/shu1988.pdf` | **Shu & Ho 1988, JAFC 36:801–803** — fed-DMHF + cysteine at pH 2.2 / 5.1 / 7.1 | **Already in the repo, no dossier.** The DMHF analogue of the Whitfield fed-norfuraneol runs |
| `scratchpad/lee2022.pdf` | **Lee et al. 2022, Food Chem. 376:131917** (HAL OA, downloaded this session) | **★ HOLD-OUT candidate.** 12 markers, absolute µmol·gDM⁻¹, 140/170/200 °C, amine on/off, raw data at `doi.org/10.15454/{IICJV3,HV3UTK,GUXKYQ,QYCBIW}` |
| `scratchpad/lee2024.pdf` | **Lee et al. 2024, Food Res. Int. 183:114183** (HAL OA, downloaded this session) | ⚠️ **PROHIBITED-DERIVATION source.** Its Arrhenius ka/kd are matrix→vapour transport and oven-wall deposition, not chemistry (§A.3.4). Keep only the scheme sentence |
| `scratchpad/schroen2022.txt` | **Schroën & Berton-Carabin 2022, Food Res. Int. 160:111621** (WUR OA, downloaded this session) | **★ The channel-4 answer.** `k_hexanal = 6·10⁻⁵ h⁻¹`, 25 °C, pH 6.7, three protein emulsifiers, hand-fitted, no SD |
| `data/articles/{martins2005,martins2005b,martins2003b,brands2001,knol2005,knol2009,knol2010}.pdf` | the van Boekel school | **[NEG] — verified this session to contain NO HMF constants.** Do not re-open for HMF |

## DO NOT ORDER — verified negatives

| item | why |
|---|---|
| `10.1111/j.1365-2621.1991.tb05389.x` — Koelsch, Downes & Labuza 1991 | One temperature (23 °C) ⇒ no Ea; dry cracker matrix at 23 % RH ⇒ not aqueous/protein; cubic-model constants in concentration/(time)³ (§D.2) |
| `10.3389/fnut.2022.940202` — Ma et al. 2022 | Already fully read (§A.3.2). Its Ea(5-HMF)=14.85 is plateau-artefact, proven by Ea(MGO)=1.84 in the same table. Only Ea(3-DG)=84.55 survives, as a weak prior |
| `10.3390/foods14122136` — Han et al. 2025 | Already fully read (§A.3.1). k8/k15 unidentifiable (SD ≈ estimate, non-monotonic in T, 10⁶× stiff). Keep only the ODE form and the lumped Ea(5-HMF)=18.836 ± 5.755 |
| any xylose→furfural biorefinery kinetics (acetic/formic/sulfuric/pTSA/zeolite) | pH 0.9–1.7, explicit first-order-in-[H⁺] rate law, 150–200 °C. 4–6 orders of magnitude off the repo's proton concentration (§C.5) |

---

# §F. DECLARED GAPS FROM THIS ROUND — carried in the `k3` §C register

**F.1 — There is no rate constant, no activation energy and no Arrhenius fit for DMHF/furaneol
formation or degradation anywhere in the accessible literature.** Five database queries and four
web searches, all returned in §B.3. The node can be given correct carbon topology (hexose = intact
C6; pentose = C1 addition, 70 % Strecker / 30 % fragmentation) and a pH sign, **but any barrier the
repo assigns to DMHF is the repo's own assumption and must be labelled as such.**

**F.2 — There is no measured furfural yield from any pentose/alanine system.** Blank & Fay 1996 is
the only pentose+alanine paper in reach and it reports no furfural at all. The brief's request to
"extend the sulfur lane's furfural feed" with alanine-system yields **cannot be met from
literature**.

**F.3 — `k3` §C.9's Frankel declaration survives this round.** No measured rate constant exists for
linoleate hydroperoxide → hexanal at cooking temperatures in aqueous or protein systems. The only
measured `k_hexanal` is at **25 °C** and was **hand-fitted without standard errors**. Extrapolating
it to the repo's window using its authors' own Q10 = 2–3 spans 11.5 decades of 10 °C steps
(≈3·10³–8·10⁵×) and is **the repo's assumption, not the authors'**.

**F.4 — The 3-DG → HMF edge the brief asked me to source is, in the best-measured systems, the
wrong edge.** Five independent studies place HMF on the fructose/fructofuranosyl-cation limb and
**furfural** on 3-DG. One study (Şen & Gökmen 2022, sucrose-rich low-moisture) dissents. **The
branch fraction is matrix-dependent and must not be hard-coded**, per `k3` §0.

**F.5 — I could not read the contents of two of the twelve items in the download list**
(Tier 3 L1 and L2). They are listed as leads with that stated on their face. **Nothing else in this
dossier is unread.**
