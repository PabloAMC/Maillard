# Round-4 literature dossier — ADVERSARIAL TEST OF FOUR DECLARED NULLS

**Date:** 2026-08-29. **Agent:** literature-research (round 4).
**Mission:** try to FALSIFY four declared literature-nulls. Owner's question: *"are you sure these
don't exist?"*
**Nothing in the repository was written, staged or modified.** Every download landed in this
scratchpad.

---

## HOW TO READ THIS FILE — verification codes (same register as round 3)

| code | meaning |
|---|---|
| **[FT-local]** | I read the **full text of a PDF held on this machine** and extracted the number from the text layer myself |
| **[FT-web]** | I read the **open-access full text on the web** (PMC / repository HTML) and extracted the number from the article body or a table |
| **[ABS]** | I read only the **abstract** (Europe PMC `resultType=core`). Numbers, if any, are abstract-level only |
| **[META]** | **Crossref-verified bibliographic metadata only** (`api.crossref.org/works/<doi>`). No content read |
| **[NEG]** | A **verified negative** — I read the source and it does NOT contain what was hoped for |
| **[REPO]** | Read from an existing repo dossier in `data/lit/extraction_dossiers/` |

**Every DOI and every author list in this file was resolved through the Crossref API in this
session** unless the row says otherwise. None is reconstructed from memory.

### ⚠️ A prerequisite that could not be met

The brief told me to read `research_round2_thiols.md` first. **That file does not exist on this
machine.** I ran `find /Users/pabloantoniomorenocasares -iname "*research_round*"` — the only
matches are `research_round3_channels.md` (repo) and its copy inside the OrbStack container.
I substituted the two places round-2's conclusions are actually carried: `k3_final_parameter_
inventory.md` §A.4 + §B.7 + §C.1/C.3/C.17 (the thiol-Ea refusals) and §C.2 + §A.7.4/A.7.5 +
`k4b_paired_thresholds_and_browning.md` (the threshold refusals), all read in full **[REPO]**.
**If round-2 searched something I re-searched, this dossier duplicates it and I could not know.**

---

# §0. THE FOUR VERDICTS, ON ONE PAGE

| # | null | verdict | the single item that decides it |
|---|---|---|---|
| **1** | **ADDUCT-FORMATION Ea** — no temperature-resolved kinetics of aldehyde→amino-nucleophile covalent adduct formation | **★ NULL FALSIFIED** (as literally stated) — **but the narrower protein-specific null SURVIVES** | **Hamzalıoğlu & Gökmen 2018** (HMF + Cys/Arg/Lys, 5/25/50 °C, Michael adduct + Schiff base, Arrhenius obeyed) **[ABS]**, plus an entire **Zamora–Hidalgo Ea ladder, 22–81 kJ/mol** **[ABS]** and one **[FT-web]** |
| **2** | **THIOL-LOSS Ea** (Chinese journals / theses / patents) | **NULL SURVIVES in those three channels — but ONE LIVE FALSIFIER exists in the mainstream literature and must be read** | **Zhou et al. 2025, LWT 218:117469** — first-order + **Arrhenius** fitted to **2-furfurylthiol** across storage temperatures. **Gold OA.** I could not get past Elsevier's 403 |
| **3** | **SECOND SULFUR TEMPERATURE LADDER** (MFT/FFT/methional, absolute, 3+ T, one system) | **NULL SURVIVES for MFT and FFT. FALSIFIED for METHIONAL** | **Chan & Reineccius 1994**, ACS Symp. Ser. **564**:127–137 — methional + DMDS + 2-acetylthiophene, aqueous, **75→115 °C**, **pH 6/7/8**, 5 min–7.5 h, *"Kinetic data are presented"* **[ABS]** |
| **4** | **MORE PAIRED SAME-PANEL MATRIX/WATER THRESHOLDS** | **★ NULL FALSIFIED — four new sets, one read in full text** | **Perry & Hayes 2016**, *Foods* **5**(2):35, **open access, full text read [FT-web]** — water vs model wine vs wine, ASTM E-679, **water and matrix measured inside the same experiment**, n = 36–44 |

**The two things a build wave should act on first:**

**0a. ★ THE ZAMORA–HIDALGO Ea LADDER RESOLVES `anantharamkrishnan2020b`'s CONDITIONAL — IN THE
NEGATIVE.** That dossier's §6 says the covalent lysine sink can matter *"only if its activation
energy is at the top of the plausible range (**Eₐ ≳ 70 kJ mol⁻¹**). At Eₐ ≲ 60 kJ mol⁻¹ the channel
is negligible on any food-process timescale"* **[REPO]**. The measured Ea for carbonyl–amine
reactions of lipid-derived aldehydes with amino acids, from one lab across seven papers, is
**22–81 kJ/mol and mostly 25–67** **[ABS]**. Every alkadienal value (**22–38**) is *below* the
threshold; the closest analogue to a lysine adduct, **4,5-epoxy-2-heptenal + lysine, is 66.5 kJ/mol
for browning and 50 kJ/mol for fluorescence** **[ABS]**. **The 50–80 kJ/mol band the repo assumed
is roughly right in span but sits at its own upper edge, and the measurement does not support the
≳70 case the sink needs.** Caveats in §A.4 — these are overall-conversion Ea, not the adduction
elementary step, and the repo must not paper over that.

**0b. ★ ALL FOUR NEW PAIRED THRESHOLD SETS GIVE SMALL RATIOS (1.4× – 12×).** Not one approaches
the 900× beef geometric mean or the 6 714× decadienal outlier in `k3` §A.7.5. This is the third
independent confirmation of `k4b`'s finding #2 — *"most of K2's dispersion was method noise while
its tails were real"* **[REPO]** — and it pushes further: **when water and matrix are measured in
the same experiment, the shift is one order of magnitude at most.** Perry & Hayes measured the same
water threshold three times in three independent experiments and got **7.51 / 8.10 / 7.57 µg/L —
a 1.08× spread** **[FT-web]**. Cross-study threshold dispersion in the repo's tables is
overwhelmingly method noise.

---

# §A. NULL 1 — ADDUCT-FORMATION Ea

## A.0 What the null says, verbatim

`anantharamkrishnan2020b_extraction.md` **[REPO]**:

> "**This paper cannot distinguish those cases — it reports one temperature and never states what
> it is (§2b [!]).** Its companion (`anantharamkrishnan2020_extraction.md` §C3) establishes only
> that adduct formation increases from 4 °C to 20 °C, **directionally, with no number and no upper
> temperature.**"

and the derivation it forced:

> "| assumed Eₐ | rate acceleration 20 → 140 °C | t½ at 140 °C | … | **70 kJ mol⁻¹** | ≈ 4,300× |
> **≈ 65 minutes** | **80 kJ mol⁻¹** | ≈ 13,900× | **≈ 20 minutes** |"

⇒ the repo's aldehyde–protein sink currently runs on an **assumed** Ea, and its size is a pure
function of that assumption.

## A.1 ★ THE FALSIFIER — a food-chemistry Ea corpus the repo has never touched

**Rosario Zamora & Francisco J. Hidalgo, Instituto de la Grasa (CSIC), Sevilla.** One lab, one
method (Arrhenius fit to formation rate constants of a measured product across a temperature
series), applied to carbonyl–amine reactions for thirty years. I found it by an author-scoped
Europe PMC query (`(AUTH:"Zamora R" OR AUTH:"Hidalgo FJ") AND ABSTRACT:"activation energy"`,
**18 hits**), then Crossref-verified every DOI and author list.

| paper | system | **Ea (kJ/mol)** | code |
|---|---|---|---|
| **Hidalgo & Zamora 1993**, *J. Food Sci.* **58**(3):667–670, `10.1111/j.1365-2621.1993.tb04352.x` | **(E)-4,5-epoxy-(E)-2-heptenal + LYSINE** — an oxidised-lipid/protein browning model. Zero-order in colour and fluorescence, followed *"as a function of time, pH, temperature and epoxy-aldehyde/lysine ratio"* | **66.5** (colour difference) · **50** (fluorescence intensity) | **[ABS]** |
| **Zamora, Alcón & Hidalgo 2013**, *JAFC* **61**(43):10231–10237, `10.1021/jf305007y` | Strecker degradation of **phenylalanine** initiated by four lipid-oxidation carbonyl classes; Ea of phenylacetaldehyde formation | **55–64** (4-oxo-2-alkenals) · **58–67** (4,5-epoxy-2-alkenals + 4-HNE) · **28–38** (2,4-alkadienals) | **[ABS]** |
| **Hidalgo, Navarro, Delgado & Zamora 2013**, *Food Res. Int.* **52**(1):206–213, `10.1016/j.foodres.2013.03.031` | **histidine → histamine** by 2,4-alkadienals, 2-alkenals, 4-oxo-2-alkenals, 4,5-epoxy-2-alkenals, 4-HNE, hydroperoxides | **22–29** (2,4-alkadienals, the most reactive class) | **[ABS]** |
| **Zamora, Delgado & Hidalgo 2012**, *Food Res. Int.* **46**(1):321–325, `10.1016/j.foodres.2011.12.029` | **phenylalanine** decarboxylation by 2,4-decadienal → β-phenylethylamine | **54** | **[ABS]** |
| **Hidalgo, Delgado, Navarro & Zamora 2010**, *JAFC* **58**(19):10512–10517, `10.1021/jf102026c` | **asparagine** decarboxylation by 2,4-decadienal | **81.0** | **[ABS]** |
| **Hidalgo, Delgado & Zamora 2010**, *Food Chem.* **122**(3):596–601, `10.1016/j.foodchem.2010.03.016` | **acrylamide + benzyl mercaptan / N-acetylcysteine** — the anaerobic Michael addition adduct | **28–30** | **[ABS]** |
| **Zamora, Delgado & Hidalgo 2010**, *JAFC* **58**(3):1708–1713, `10.1021/jf903378x` | **acrylamide + amines / amino acids / polypeptides** — Michael adduct formation vs its reverse elimination | **no numbers in the abstract**; the abstract states only: *"acrylamide seems to be converted into its Michael adduct with a **lower** activation energy than the elimination reaction of the Michael adduct"* | **[ABS]** |
| **Zamora, Navarro, Aguilar & Hidalgo 2015**, *Food Chem.* **174**:89–96, `10.1016/j.foodchem.2014.11.034` | thermal degradation of 2-alkenals / 2,4-alkadienals — **not an adduct reaction**, see §A.3 | **21.3–29.6** | **[FT-web]** — downloaded and read |

**And the second, independent falsifier — and it also closes round-3's Tier-3 lead L1:**

**Hamzalıoğlu, A. & Gökmen, V. 2018**, *"Investigation and kinetic evaluation of the reactions of
hydroxymethylfurfural with amino and thiol groups of amino acids"*, **Food Chemistry 240:354–360**,
**`10.1016/j.foodchem.2017.07.131`**, PII **S0308814617312852** — Crossref-verified **[META]**.
Verbatim abstract **[ABS]**:

> "reactions of hydroxymethylfurfural (HMF) with selected amino acids (**arginine, cysteine and
> lysine**) were investigated in HMF-amino acid (high moisture) and Coffee-amino acid (low
> moisture) model systems at **5, 25 and 50 °C** … High-resolution mass spectrometry (HRMS)
> analyses of the reaction mixtures **confirmed the formations of Michael adduct and Schiff base of
> HMF with amino acids**. Calculated **pseudo-first order reaction rate constants** were in the
> following order; k_Cysteine > k_Arginine > k_Lysine for high moisture model systems … **The
> temperature dependence of the rate constants was found to obey the Arrhenius law in a temperature
> range of 5–50 °C under both low and high moisture conditions.**"

⇒ **This is a temperature-resolved, three-point kinetic study of an ALDEHYDE forming BOTH a Michael
adduct and a Schiff base with three amino-acid nucleophiles including a thiol, with an explicit
Arrhenius claim, in both a high-moisture and a low-moisture (coffee) matrix.** It is the exact
shape of measurement the null said did not exist. **⚠️ The Ea VALUES are not in the abstract and I
have not read them.** The paper is `closed` on OpenAlex (checked this session); it is **download
item #1**.

> ★ **BONUS FOR ROUND 3.** `research_round3_channels.md` Tier-3 lead **L1** is this paper, listed
> as *"⚠️ I could not resolve a DOI for it. Resolve the DOI and read the abstract before ordering.
> Do not cite it until then."* **The DOI is now resolved, the authors are now known (Hamzalıoğlu &
> Gökmen — the same FOQUS lab as the round-3 HMF picks), and the abstract is now read.** It is not
> only the HMF-sink paper round 3 hoped for; it is also the adduct-Ea paper round 4 was sent to
> find. **The repo has no HMF sink and this paper is the only published kinetics for one.**

**A third, adjacent falsifier — same measurement class, wrong electrophile:**
**Marianne N. Lund's lab (Copenhagen)** does quinone–amine kinetics with explicit Arrhenius fits.
**Liu, J.; Poojary, M. M.; Thygesen, M. B.; Andersen, M. L.; Lund, M. N. (2023)**, *"Temperature
affects the kinetics but not the products of the reaction between 4-methylbenzoquinone and lysine"*,
*Food Res. Int.* **163**:112187, `10.1016/j.foodres.2022.112187` — Crossref-verified **[META]**:
stopped-flow, pseudo-first-order, **15–45 °C**, Lys / Nᵅ-acetyl-Lys / Nᵋ-acetyl-Lys,
*"the temperature dependence of reaction rates followed the Arrhenius law, with activation energies
in the order: Lys < Nᵋ-acetyl Lys < Nᵅ-acetyl Lys"* — **ordering only, no values in the abstract.**
Companion: `10.1016/j.foodchem.2023.137473`, 4MBQ + **β-lactoglobulin** and amino acids, apparent
k₂ at **25 °C only**, order β-LG > His > Trp > Arg **[ABS]**.
⚠️ **A quinone is not an aldehyde.** Transferring a quinone Ea to a Schiff-base sink would be
exactly the class of substitution `k3` §C.14 refuses for the Yaghmur lump. **Ingest as an analogue
bound, never as the aldehyde value.**

## A.2 ⚠️ WHERE THE NULL SURVIVES — and it is the part the brief bet on

**THE CROSS-DOMAIN BET FAILED. The biomedical 4-HNE / alkenal / acrolein / MDA protein-adduct
literature measures second-order rate constants at ONE temperature and publishes no activation
energies.** I searched it hard and this is a **[NEG]**:

- Europe PMC `ABSTRACT:"activation energy" AND ABSTRACT:"4-hydroxynonenal"` → **0 hits**
- Europe PMC `ABSTRACT:"activation energy" AND ABSTRACT:acrolein` → 16 hits, **zero biological**
  (Criegee intermediates, Diels–Alder DFT, glycerol dehydration over zeolites, PU foam autoxidation)
- Europe PMC `ABSTRACT:"rate constant*" AND ABSTRACT:(aldehyde OR alkenal OR acrolein) AND
  ABSTRACT:(cysteine OR glutathione) AND ABSTRACT:temperature*` → **0 hits**
- Europe PMC `TITLE:(kinetic* OR "rate constant*") AND TITLE:(acrolein OR hydroxynonenal OR
  oxononenal OR crotonaldehyde OR "unsaturated aldehyde*") AND TITLE:(protein* OR albumin OR lysine
  OR cysteine OR glutathione OR "amino acid*")` → **2 hits**, both single-temperature
- three WebSearches on HNE-protein adduct kinetics / Michael-addition Arrhenius

The three canonical papers, abstracts read **[ABS]**, all single-temperature:

| paper | what it gives | what it does NOT give |
|---|---|---|
| **Doorn, J. A. & Petersen, D. R. (2002)**, *Chem. Res. Toxicol.* **15**(11):1445–1450, `10.1021/tx025590o` — "Covalent Modification of Amino Acid Nucleophiles by the Lipid Peroxidation Products 4-Hydroxy-2-nonenal and 4-Oxo-2-nonenal" **[META]** | the reference k₂ set for HNE/ONE + Cys/His/Lys | **no temperature series, no Ea** |
| **Szapacs, M. E.; Riggins, J. N.; Zimmerman, L. J.; Liebler, D. C. (2006)**, *Biochemistry* **45**(35):10521–10528, `10.1021/bi060535q` **[META]** | **site-resolved** k for six HSA histidines, *"Rate constants ranged over 4 orders of magnitude, with the order of reactivity being H242 > H510 > H67 > H367 > H247 ≈ K233"* — H242 is a fatty-acid-binding-cavity hot spot | **one temperature, no Ea** |
| **Sauerland, M.; Mertes, R.; Morozzi, C.; Eggler, A. L.; Gamon, L. F.; Davies, M. J. (2021)**, *Free Radic. Biol. Med.* **169**:1–11, `10.1016/j.freeradbiomed.2021.03.040` **[META]** | k₂ for acrolein/crotonaldehyde/DMF/cyclenones with N-Ac-Cys, GSH, **BSA, creatine kinase, papain, GAPDH, Keap-1**; *"k₂ values for N-Ac-Cys and GSH vary by > 250-fold"*; *"A linear inverse correlation for acrolein with the thiol pKa indicates that the **thiolate anion is the reactive species**"* | **one temperature, no Ea** |

**⇒ THE NARROW NULL, RESTATED AND STILL STANDING: there is no measured activation energy for the
covalent adduction of a saturated n-alkanal (hexanal above all) to a PROTEIN.** Everything found in
§A.1 is an amino acid or a small amine, not a protein, and — except the 1993 lysine paper — the
measured product is a *downstream* Strecker/decarboxylation product, not the adduct.

**The three specific gaps the round-4 search could not fill, reported as gaps rather than filled:**
1. **No hexanal + lysine/BSA/whey kinetics at more than one temperature exists.** Searched via
   `TITLE:(hexanal OR "2-hexenal" OR aldehyde*) AND TITLE:(binding OR adduct* OR covalent) AND
   TITLE:(protein OR "beta-lactoglobulin" OR whey OR soy) AND ABSTRACT:temperature*` (2 hits, neither
   relevant) and four WebSearches. **[NEG]**
2. The one modern near-miss is **`10.1038/s41538-026-00874-9`** (*npj Sci. Food* 2026, *"Novel
   insights into the binding mechanisms of selected aldehydes during heat-induced protein
   unfolding"*, six aldehydes + myofibrillar proteins, stage-heating) **[ABS]**. Its
   *"thermodynamic parameters"* are **binding** ΔH/ΔS, not an **Ea** for adduct formation, and it
   is the non-covalent lane. Its proteomics finding is usable as topology: *"saturated aldehydes
   predominantly formed Schiff bases with the lysine ε-amino group. Unsaturated aldehydes,
   especially (E,E)-2,4-decadienal, undergo **both** Schiff base reactions **and** Michael addition
   with cysteine, histidine, and tryptophan residues."*
3. The nearest thing to a Schiff-base Ea in the crosslinking literature is a **REVERSAL** Ea, not a
   formation Ea: **Kennedy-Darling & Smith 2014**, *Anal. Chem.* **86**(12):5678–5681,
   `10.1021/ac501354y`, PMC4063333 **[FT-web]** — formaldehyde protein–DNA cross-link **reversal**,
   4 / 23 / 37 / 47 °C, half-lives **179 / 45.9 / 22.7 / 11.3 h**, *"the activation energy (Eₐ) was
   determined to be **47 kJ mol⁻¹**, while from the Eyring plot, the enthalpy of activation (ΔH‡)
   and entropy of activation (ΔS‡) were determined to be **44 kJ mol⁻¹** and **−0.1 kJ mol⁻¹ K⁻¹**"*.
   The paper states no formation Ea. **⚠️ Do not file this as an adduct-formation barrier.**

## A.3 ★ AN UNSOUGHT FINDING THE LIPID LANE SHOULD HAVE — full text, numbers extracted

**Zamora, Navarro, Aguilar & Hidalgo 2015**, *Food Chemistry* **174**:89–96,
`10.1016/j.foodchem.2014.11.034`. **Green OA at DIGITAL.CSIC** (`hdl.handle.net/10261/116078`);
postprint downloaded to `scratchpad/zamora2015.pdf`, text at `scratchpad/zamora2015.txt`. Every
number below I pulled from the text layer myself **[FT-web]**.

**System, verbatim (p. 5, lines 85–88):** *"a solution of the aldehyde (0–10 µmol) in tetrahydrofuran
(80 µL) was treated with 420 µL of **0.2 M buffer (pH 2.15–11)** and, then, heated for the indicated
time and temperature in **closed test tubes under air**."* Quantification by **LC-MS/MS after
dansylhydrazine derivatisation with a formaldehyde-d₂ internal standard**. Temperatures, from the
figure captions: **200, 160 and 120 °C** — three points.

| aldehyde heated | product | **yield** | **Ea (kJ/mol)** | anchor |
|---|---|---:|---:|---|
| 2-pentenal | propanal | ~10 % | **25.2** | line 221 |
| **2-octenal** | **hexanal** | **18.0 %** | **25.3** | lines 228, 238 |
| 2,4-heptadienal | propanal | 9.8 % | **25.2** | lines 254, 261 |
| 2,4-heptadienal | 2-pentenal | 1.0 % | **22.5** | lines 254, 261 |
| **2,4-decadienal** | **hexanal** | **11.5 %** | **21.3** | lines 278, 285 |
| 2,4-decadienal | 2-octenal | 0.8 % | **29.6** | lines 278, 285 |

Verbatim summary, p. 14 line 309: *"The Ea for the breakage of the different carbon-carbon double
bonds was always very similar and was about **25 kJ/mol**. However, alkanals were produced to a much
higher extent than 2-alkenals. Thus, **10–18 % of the initial either 2-alkenal or 2,4-alkadienal was
converted into alkanal after 1 h heating at 200 °C** and only about 1 % of the initial 2,4-alkadienal
was converted into 2-alkenal."* pH optimum for hexanal from 2-octenal: *"maximum around pH 10"*
(line 225); from 2,4-decadienal: *"maximum at about pH 8"* (line 270).

> ★ **WHY THIS MATTERS TO THE MATRIX-PATH ACCURACY PLAN.** The repo currently routes hexanal from
> **13-hydroperoxide cleavage only**. Here are **two additional hexanal sources measured
> absolutely, in aqueous buffer, at the repo's own processing temperatures**: 2-octenal → hexanal at
> **18 %** and 2,4-decadienal → hexanal at **11.5 %**, each with an Ea near **25 kJ/mol** — i.e. a
> **weakly temperature-dependent** route that will *not* switch off when the model cools. If the
> repo's network already carries 2-octenal and 2,4-decadienal as lipid products, it is currently
> **under**-counting hexanal formation while the calibration is fighting an **over**-prediction.
> ⚠️ **This cuts against the over-prediction fix, not for it — record it before it is discovered by
> a fit.**
> ⚠️ Second caveat: the maxima are at **pH 8–10**, above the repo's pH 5–7 window, and the pH
> profiles are figure-bound (Figs. 3A / 5A) — I did not digitise them.

## A.4 ⚠️ THE LABEL DISCIPLINE THE ZAMORA Ea MUST CARRY — four caveats, all real

1. **They are OVERALL-CONVERSION Ea, not adduction-step Ea.** Every Zamora Ea in §A.1 is fitted to
   the appearance of a *downstream* product (phenylacetaldehyde, histamine, 3-aminopropionamide,
   β-phenylethylamine). The carbonyl–amine condensation is the *first* step of that sequence; the
   *rate-limiting* step may be the decarboxylation. **Same species as the Bornhorst norfuraneol Ea
   and Han's lumped Ea(5-HMF) the repo already carries with a lump label** (`k3` §C.13, round-3
   §A.3.1). Label: *"apparent overall Ea for carbonyl-amine-initiated conversion of amino acid X by
   lipid carbonyl Y"*, never *"the Schiff-base barrier"*.
2. **Amino acids, not proteins.** Free histidine, phenylalanine, asparagine, lysine — never a
   protein-bound ε-NH₂. The same objection the repo raises against Ajandouz/Bornhorst in reverse
   (`k3` §C.13 complains those use proteins where the repo wants free amino acids; here it is the
   other way round).
3. **The 1993 lysine paper is the ONLY true aldehyde+lysine member, and its measured quantity is
   COLOUR and FLUORESCENCE, not adduct concentration.** Ea = 66.5 / 50 kJ/mol are browning-pigment
   and fluorophore development rates. The aldehyde is **(E)-4,5-epoxy-(E)-2-heptenal**, an ω-3
   epoxyalkenal — far more electrophilic than hexanal.
4. **Temperature ranges are unread for every [ABS] row.** I read the temperature ladder only for
   the 2015 paper (120/160/200 °C). For the other seven the range is unknown to me and **must be
   checked before any of these Ea is extrapolated** — `k3` §C.13's *"Extrapolation of these Ea above
   100 °C | NOT LICENSED"* rule applies until the range is read.

---

# §B. NULL 2 — THIOL-LOSS Ea, THE THREE COVERAGE HOLES

## B.0 What the null says, verbatim (`k3` §C.1) **[REPO]**

> "**No activation energy for thiol consumption exists anywhere in this basket.** … **The literature
> reviewed supports NO activation energy for thiol consumption.**"

## B.1 ★ THE ONE LIVE FALSIFIER — and it is not in any of the three holes I was sent to search

**Zhou X, Diao Y, Zhang J, Yu X, Wang G, Gu Y, Ren D, Li S, Dong W, Yi L. (2025)**
*"Sulfury/roasty fading indicators in roasted coffees: Their contribution and applicability in coffee
freshness perception and prediction"*, **LWT 218:117469**, `10.1016/j.lwt.2025.117469` —
Crossref-verified **[META]**. Verbatim abstract **[ABS]**:

> "the variation of volatile profiles and aroma perceptions of coffee brews of roasted coffee beans
> treated with **accelerated storage experiments** were analyzed. Besides, the **first-order kinetic
> reaction model and Arrhenius formula were adopted for the establishment of prediction models**.
> The results show that **2-furfurylthiol**, methanethiol and 2-ethyl-3-methyl-pyrazine are reliable
> compounds that not only contribute sulfury/roasty perception to coffee brews, but can also act as
> reliable indicators for the shelf life prediction … **during different storage temperature**."

⇒ **An explicit Arrhenius fit to first-order 2-FFT loss across ≥2 storage temperatures.** If it
carries an Ea, `k3` §C.1's *"NO activation energy for thiol consumption"* is false as written.

**⚠️ I DID NOT READ IT AND I AM NOT REPORTING NUMBERS I DID NOT READ.** OpenAlex reports
`oa_status: gold`, but ScienceDirect returned **HTTP 403** to both WebFetch and curl (with and
without a browser UA), and the SSRN preprint (`10.2139/ssrn.4836239`, Crossref-verified) also
403s. **It is gold OA — the owner can open it in a browser in ten seconds.** Download item #2.

**⚠️ AND THE CAVEAT THAT WILL PROBABLY DEFUSE IT.** This is an **accelerated-storage** Ea for
whole/ground roasted coffee, i.e. the same species as van Seeventer's 50 °C storage rates. `k3`
§C.17's finding stands over it before it is even read:

> "**the processing sink and the storage sink are DIFFERENT CHANNELS.** That qualitative conclusion
> is the robust result" — `vanseeventer2001_z3_addendum.md` **[REPO]**

**So the honest prediction is: this paper will falsify the null's LETTER and leave its SPIRIT
intact.** It gives a storage-regime Ea for FFT; the repo's open question is a
115–180 °C processing-regime Ea. **Read it, amend `k3` §C.1's wording, and do NOT bolt the Ea onto
the processing sink** — the exact error `k3` §C.17 already names.

## B.2 HOLE (a) — CHINESE-LANGUAGE FOOD JOURNALS: **NULL SURVIVES** **[NEG]**

Searched: two Chinese-language WebSearches
(`2-甲基-3-呋喃硫醇 降解 动力学 活化能 热反应肉味香精`;
`糠硫醇 2-呋喃甲硫醇 稳定性 动力学 半衰期 温度 咖啡 香气 研究`), plus Europe PMC on the
English-indexed abstracts of Chinese work (`TITLE:("2-furfurylthiol" OR "2-methyl-3-furanthiol" OR
"furfuryl mercaptan")` → 31 hits, read all titles).

**What came back:** the Chinese FFT/MFT literature is enormous but it is **Baijiu biosynthesis and
analysis**, not degradation kinetics — `10.1021/acs.jafc.7b06125` (C–S lyase in Baijiu yeast),
`10.1021/acs.jafc.7b01359` (STR3/CYS3), `10.1016/j.foodres.2020.109757` (raising FFT by inoculating
a cysteine producer), `10.1021/acs.jafc.0c04170` (MFT/FFT as Baijiu aroma-type markers). **Not one
degradation-rate paper.** No Chinese source with a thiol Ea was found.

**The strongest negative control I can offer**, and it is a good one: **Zhang, G.; Xiao, P.; Yuan,
M.; Li, Y.; Xu, Y.; Li, H.; Sun, J.; Sun, B. (2023)** — Crossref-verified author list **[META]** —
*"Roles of sulfur-containing compounds in fermented beverages with 2-furfurylthiol as a case
example"*, *Front. Nutr.* **10**:1196816, `10.3389/fnut.2023.1196816`, PMC10348841 — a dedicated,
open-access, Chinese-authored (⚠️ affiliations **not** checked by me)
review of FFT with a whole section (§4.3.1) on the **staling pathway**. I read it **[FT-web]**.
**It reports no rate constant, no half-life, no activation energy and no temperature-dependent
kinetics anywhere.** Its only quantitative-adjacent sentence is *"The degradation rate of FFT was
positively correlated with the activity of hydroxyl radical"* — qualitative, citing Blank.

## B.3 HOLE (b) — THESES: **NULL SURVIVES, but the thesis trick was not properly testable**
**[NEG, with a caveat about method]**

Searched ProQuest-adjacent and repository-indexed routes via WebSearch (English and Chinese). **No
dissertation with thiol-stability kinetics surfaced.** Two named leads worth a manual library check
that I could NOT verify and am therefore listing as leads, not findings:

- **Zhenchun Sun** is first author of both FFT-stability papers in §B.5, whose Crossref-verified
  co-author lists include **Xiaoming Zhang** (Jiangnan University) and **Ian D. Fisk** (Nottingham).
  ⚠️ **His affiliation and supervisor are my INFERENCE from the co-author list, not verified.** If
  the inference is right, a Chinese-language thesis on FFT stability likely exists on CNKI. **I
  could not reach CNKI from this session and have verified nothing about any such thesis.** Do not
  cite it until the record is in hand.
- **F. Chan**, co-author of the methional kinetics chapter (§C.1), worked in Reineccius's group at
  the **University of Minnesota**, whose **conservancy.umn.edu** hosts theses. A Chan thesis would
  carry the full kinetic tables the 11-page chapter compresses. **Unverified — I did not find the
  record.**

## B.4 HOLE (c) — PATENTS: **NULL SURVIVES** **[NEG]**

One dedicated WebSearch on flavour-house patent disclosures of thiol stability kinetics. Returns
were product listings, a β-cyclodextrin inclusion patent, and a coffee-aroma stabilisation patent
(sodium thiosulfate/sulfite) reporting **retention percentages after one year at room temperature**,
not rates. **No patent disclosing a rate constant or an Ea for thiol loss was found.**

## B.5 THE MAINSTREAM RE-CHECK — three verified negatives worth recording

| source | why it looked promising | verdict |
|---|---|---|
| **Weerawatanakorn, M.; Wu, J.-C.; Pan, M.-H.; Ho, C.-T. (2015)**, *"Reactivity and stability of selected flavor compounds"*, *J. Food Drug Anal.* **23**(2):176–190, `10.1016/j.jfda.2015.02.001`, PMC9351765 — Crossref-verified **[META]** | **the standard review of exactly this question**, by Ho's group | **[FT-web] [NEG]** — **no rate constant, no half-life, no Ea, no Arrhenius parameter for any thiol.** Its only quantitative sentence is *"The degradation is temperature dependent, which caused **20 % loss of FFT at room temperature compared to 90 % loss at 37 °C after 1 hour**"* — the Blank 2002 two-point contrast, restated |
| **Blank et al. 2002**, `10.1021/jf011329m` | a **same-lab, same-system 22 °C vs 37 °C pair** — better than K1's two-lab pair | **ALREADY IN THE REPO AND ALREADY REFUSED. [REPO]** `blank2002_extraction.md`: *"There are **NO rate constants, NO half-lives, NO Arrhenius parameters and NO activation energies** anywhere in the paper"*, and *"the 37 °C and 22 °C arms of Table 3 were run at **different initial FFT concentrations** (3.8 mM vs 3.2 mM), which confounds the temperature contrast."* **Round 4 adds nothing here** |
| **Sun, Yang, Liu, Linforth, Zhang & Fisk 2019**, *"Aroma binding and stability in brewed coffee: A case study of 2-furfurylthiol"*, *Food Chem.* **295**:449–455, `10.1016/j.foodchem.2019.05.175` | abstract says *"The reduction in available 2-FFT was investigated at **different pH and temperatures**"* | **UNRESOLVED — [ABS] only.** Green OA at Repository@Nottingham (`worktribe.com/output/2226865`, Unpaywall-confirmed) but **Cloudflare-blocked from this session** (curl and WebFetch both). **Download item #3.** Whether it contains ≥3 temperatures with fitted rates is **unknown to me** |

---

# §C. NULL 3 — A SECOND SULFUR TEMPERATURE LADDER

## C.1 ★ FALSIFIED FOR METHIONAL

**Chan, F. & Reineccius, G. A. (1994)**, *"Kinetics of the Formation of Methional, Dimethyl
Disulfide, and 2-Acetylthiophene via the Maillard Reaction"*, in **Mussinan & Keelan (eds.),
*Sulfur Compounds in Foods*, ACS Symposium Series 564, ch. 10, pp. 127–137**,
**`10.1021/bk-1994-0564.ch010`** — Crossref-verified title, authors, book, series volume, pages and
date **[META]**. Verbatim abstract **[ABS]**:

> "There is **a lack of information on the reaction kinetics of these Strecker aldehydes** as well as
> other flavor compounds. Thus a kinetic study on the formation of **methional and two secondary
> products (dimethyl disulfide and 2-acetylthiophene)** from the reaction of **amino acids
> (0.075 mole) and glucose (0.5 mole) in aqueous model systems** was conducted. Systems were heated
> at **temperatures from 75 to 115 °C** at times from **5 min to 7.5 h** and **pH's of 6, 7, and 8**.
> **Kinetic data are presented and discussed.**"

**Why this is the round-4 find for the sulfur lane, if it holds up:**
- **A genuine temperature ladder** (75→115 °C) crossed with **a pH ladder at exactly the repo's
  pH 6/7/8**, in **one lab, one aqueous system** — the design `k3` §B.2 says nothing in the corpus
  has.
- **Its upper end, 115 °C, is exactly Zhang 2024's single temperature** (`k3` §A.4). If Chan &
  Reineccius quantified absolutely, the repo gains its first sulfur-lane Arrhenius axis.
- It couples a **Strecker aldehyde to two secondary sulfur products** (DMDS, 2-acetylthiophene) in
  the same experiment — a formation/consumption pair the repo has for no sulfur compound.

**⚠️ THREE THINGS I DID NOT VERIFY AND WILL NOT ASSERT:**
1. **Whether quantification is absolute with an internal standard.** The abstract does not say. A
   1994 ACS-symposium chapter could easily be peak areas — the exact defect `k3` §C.15 records for
   Amrani-Hemaimi and Cerny.
2. **How many temperatures.** *"from 75 to 115 °C"* could be two points or five.
3. **Whether an Ea is printed.** *"Kinetic data are presented"* is not *"activation energies are
   reported"*. ACS blocks the chapter (403 to WebFetch this session).
**Do not record this as a ladder until someone reads pp. 127–137.** Download item #4.

## C.2 NULL SURVIVES FOR MFT AND FFT — and one hopeful route was already closed by the repo

**No absolute, internal-standard quantification of MFT or FFT at 3+ temperatures in one system was
found.** Queries run, all **[NEG]**:
- `ABSTRACT:("2-methyl-3-furanthiol" OR "2-furfurylthiol") AND ABSTRACT:temperature* AND
  ABSTRACT:(quantit* OR "internal standard" OR "isotope dilution")` → **6 hits**, listed below
- `ABSTRACT:("volatile sulfur"…) AND ABSTRACT:(kinetic* OR Arrhenius OR "activation energy") AND
  ABSTRACT:(Maillard OR "model system") AND ABSTRACT:temperature*` → 1 hit, an ion-trap physics paper
- `ABSTRACT:("meaty flavor" OR "process flavor" OR "reaction flavor") AND ABSTRACT:temperature AND
  ABSTRACT:(kinetic* OR "activation energy")` → 1 hit, a 1997 process-engineering chapter (§C.3)
- one WebSearch on cysteine/xylose model systems at 100–130 °C

**The one hit that looked like the answer was already checked and refused by the repo:**
**Hofmann & Schieberle 1998** (`10.1021/jf9705983`) — whose abstract promises *"model systems
**varying in temperature**, pH value, or water content … determined by using **stable isotope
dilution assays**"* **[ABS]**. `hofmann1998_reconciliation.md` §1 **[REPO]**:

> "## 1. ★ ANSWER TO CROSS-CHECK #1b — THE TEMPERATURE / WATER-CONTENT SERIES **DOES NOT EXIST**"
> … *"varying in temperature, pH value, or water content" (p. 235) **promises more than the
> Experimental delivers**"* — the entire dataset is **145 °C vs 180 °C**, with a_w, temperature,
> time and buffer molarity all varying together.

**Round 4 confirms the repo was right, and adds that the abstract is the trap.** Anyone re-finding
this paper by keyword will believe it has a temperature series. It does not. **Record the abstract's
wording as a known snare.**

**The three remaining leads, none verified as ladders:**

| lead | why | status |
|---|---|---|
| **`Kang2026.pdf`** — already on disk | `k4c_wave_consolidated.md` §1.2 **[REPO]**: *"⚠⚠ Kang 2026 — the 100/120/140 °C per-compound data are **NOT IN THE PDF**. The ladder **exists as an experiment** but its numeric table is [absent]"* | **★ THE HIGHEST-VALUE RETRIEVAL FOR THIS NULL IS AN SI, NOT A NEW PAPER.** A three-temperature Cys–Xyl sulfur ladder already exists and the repo owns the paper. Chase the SI |
| **Zhang 2024 SI** | `k3` §C.18 already ranks it *"the single highest-value retrieval in the corpus"* — Fig. S2 is *"the temperature sweep"* | unchanged by round 4 |
| **`10.1021/jf4008086`** — Schieberle-school hazelnut roasting, SIDA + GC×GC-TOF-MS, 2-furfurylthiol among six marker odorants **[ABS]** | *"roasted under **defined conditions**"* — cultivar comparison, not obviously a temperature ladder | weak; do not order on this basis |

## C.3 One lead I could not evaluate

**"Implementation of process kinetics to scale up automated thermally reacted flavor processes"**,
in *Flavor Technology: Physical Chemistry, Modification, and Process* (ACS Symp. Ser., 1995/1997)
**[ABS]** — batch reaction-flavour processing, *"process kinetics was sequentially developed from
laboratory, pilot plant, and scale up in production"*, *"a **semi-empirical process model**"*,
*"the non-linear correlation between the **virtual reaction temperature** with the apparent probe
temperature"*. **This is process-control engineering, and "virtual reaction temperature" is a
control variable, not a chemical one.** I found no DOI for it and read no content beyond the
abstract. **Listed for completeness; I would not order it.**

---

# §D. NULL 4 — PAIRED SAME-PANEL MATRIX/WATER THRESHOLDS

## D.0 What the null says (`k3` §C.2, from K2) **[REPO]**

> "**NO SAME-METHOD MATRIX-VS-WATER PAIR EXISTS ANYWHERE IN THESE TEN PAPERS.** Every ratio in §A.8
> crosses labs, decades and criteria."

`k4b` **[REPO]** already softened this to *"the corpus now has TWO same-panel, same-method paired
sets"* (Hong 2020; Leksrisompong 2010). **Round 4 adds four more.**

## D.1 ★ THE BEST NEW SET — read in full text, numbers extracted myself

**Perry, D. M. & Hayes, J. E. (2016)**, *"Effects of Matrix Composition on Detection Threshold
Estimates for Methyl Anthranilate and 2-Aminoacetophenone"*, ***Foods* 5(2):35**,
**`10.3390/foods5020035`** — Crossref-verified **[META]**; **open access**, full text read from
PMC5302346, tables transcribed by me **[FT-web]**.

**Design (verbatim, §2.1):** *"A total of **five experiments** … See Table 1 for a summary of the
combinations of delivery matrices (**water, model wine, and wine**) and odorants (MA or 2AAP)."*
**Experiments 1, 4 and 5 measure water AND a matrix inside the same experiment** — Experiment 5 is
explicitly *"a **within-subjects** design"*. Method: *"the forced-choice methodology outlined in
**ASTM method E-679**"*, six ascending triads, two blanks + one spike, orthonasal.
Matrices: RO water; **model wine** (9.65 % v/v ethanol, TA 8.15 g/L tartaric, **pH 3.10**, matched
to the wine by measurement); **wine** (bulk Franz Reinhart Riesling).

**Table 2 — Detection Thresholds (µg/L) for Methyl Anthranilate [FT-web]:**

| | water DT | water 95 % CI | n | wine DT | wine CI | n | model wine DT | CI | n |
|---|---:|---|---:|---:|---|---:|---:|---|---:|
| **Experiment 1** | **7.51** | 3.1–17.9 | 38 | **45.0** | 21.2–95.5 | 36 | – | | |
| Experiment 2 | 8.10 | 3.82–17.2 | 43 | – | | | – | | |
| **Experiment 4** | **7.57** | 2.30–25.0 | 42 | – | | | **89.4** | 28.2–283 | 40 |

**Table 3 — Detection Thresholds (µg/L) for 2-aminoacetophenone [FT-web]:**

| | water DT | CI | n | wine DT | CI | n | model wine DT | CI | n |
|---|---:|---|---:|---:|---|---:|---:|---|---:|
| Experiment 2 | 1.00 | 0.685–1.47 | 43 | – | | | – | | |
| Experiment 3 | – | | | **10.5** | 3.79–29.3 | 38 | – | | |
| **Experiment 5** | **1.17** | 0.614–2.24 | 44 | – | | | **5.56** | 2.94–10.5 | 43 |

**The ratios the repo can actually use — `cross_study_cross_method: FALSE`:**

| pair | design | ratio |
|---|---|---:|
| MA, wine / water (**Exp 1, same experiment**) | between-subjects within one experiment | **6.0×** |
| MA, model wine / water (**Exp 4, same experiment**) | between-subjects within one experiment | **11.8×** |
| **2AAP, model wine / water (Exp 5, WITHIN-SUBJECTS)** | **the strongest design in this dossier** | **4.75×** |
| 2AAP, wine / water (Exp 2 vs Exp 3) | between-subjects, **across** experiments | 10.5× |

**★ FOUR STRUCTURAL RESULTS, all quotable [FT-web]:**
1. **The water threshold reproduces to 1.08× across three independent experiments** (7.51 / 8.10 /
   7.57 µg/L, n = 38/43/42). Authors: *"the detection thresholds for MA in water are **consistent
   and reproducible across experiments**, suggesting the methods used here to determine the
   detection threshold were **robust**."* ⇒ **the ~30× dispersion in `k3` §A.7.5's cross-study
   ratios is method noise, not chemistry.**
2. **The matrix ORDERING INVERTS between two structurally similar compounds.** MA: model wine
   (89.4) > wine (45.0). 2AAP: wine (10.5) > model wine (5.56). Two anthranilate-type aromatics,
   one panel protocol, opposite matrix rankings. **No scalar matrix factor and no partition-based
   correction can produce that** — the same conclusion `k4b` §B reached three ways.
3. **A model wine matched on ethanol, pH and TA does NOT reproduce the real wine's threshold shift.**
   Authors' own limitation: the model wine *"was included here as an attempt to differentiate
   perceptual masking effects from physical partitioning effects."* It failed to separate them.
4. **Author caution, verbatim:** *"we caution that the apparently higher value for model wine versus
   the wine should not be over interpreted, as the confidence intervals for these two matrices show
   substantial overlap."* **Carry that caution with result 2.**
⚠️ **Table 4's "mean across all matrices" water values (MA 7.73, 2AAP 1.09) are means of individual
sigmoidal fits pooled across experiments and carry NO confidence interval.** Prefer Tables 2/3.

## D.2 ★ THE CLEANEST SINGLE-PANEL PAIR IN THE LITERATURE

**Eisele, T. A. & Semon, M. J. (2005)**, *"Best Estimated Aroma and Taste Detection Threshold for
Guaiacol in Water and Apple Juice"*, ***J. Food Sci.* 70(4)**,
**`10.1111/j.1365-2621.2005.tb07201.x`** — Crossref-verified **[META]**. Verbatim abstract **[ABS]**:

> "determine the **best estimate threshold (BET)** for detection of guaiacol in **water and
> commercial pasteurized apple juice from concentrate** using the **forced-choice ascending
> concentration method of limits** with an **experienced 17-member sensory panel**. The mean BET for
> **aroma** detection of guaiacol in water and apple juice was **0.48 ppb and 0.91 ppb**,
> respectively. The mean BET for **taste** detection … was **0.17 ppb and 0.24 ppb** … Individual
> aroma BET values ranged from 0.06 to 4.71 ppb in water and 0.17 to 4.71 ppb in apple juice."

⇒ **One panel, 17 people, one method, both matrices, both modalities. Aroma juice/water = 1.9×.
Taste juice/water = 1.4×.** This is precisely the design `k3` §C.2 says does not exist, and the
matrix penalty is **under 2×**. ⚠️ Guaiacol is a phenol, not an aldehyde or a thiol, and apple juice
is a clear, low-protein, low-fat matrix — the mildest possible test. **Its value is as the LOWER
BOUND of the matrix-shift distribution**, which the repo currently has no anchor for.

## D.3 A THIRD SET — a matrix LADDER for one compound, one trained panel

**Ziegler, M., Gök, R., Bechtloff, P., Winterhalter, P., Schmarr, H.-G. & Fischer, U. (2019)**,
*"Impact of matrix variables and expertise of panelists on sensory thresholds of
1,1,6-trimethyl-1,2-dihydronaphthalene known as petrol off-flavor compound in Riesling wines"*,
***Food Quality and Preference* 78:103735**, `10.1016/j.foodqual.2019.103735` — Crossref-verified
**[META]**. Verbatim abstract **[ABS]**:

> "the impact of matrix variables such as **ethanol and carbonation** on the odor detection threshold
> of TDN. **Increasing carbonation nearly doubled the detection threshold in water**, however, this
> effect could not be observed in alcoholic matrices. **Ethanol enhanced detection thresholds
> presumably due to diminished volatilization, which was only partially supported by measurement of
> partition coefficients using GC–MS.** Overall, the impact of matrix changes was **small ranging
> from 1.1 µg/L in still water to 4.0 µg/L in Riesling sparkling wine**. **Consumer detection
> threshold of 14.7 µg/L … exceeded the threshold of the trained panel by a factor of five**."

**Three things the repo should take from this [ABS]:**
1. **Water → real sparkling wine = 3.6×.** A fourth small ratio.
2. **A DIRECT TEST OF THE `f_free` LOGIC, AND IT PARTLY FAILS.** The authors attributed the ethanol
   effect to reduced volatilisation and then **measured the partition coefficients** — which *"only
   partially supported"* the perceptual shift. **This is a fourth independent demonstration of
   `k2` §D.4.5 / `k4b` §B: you cannot get matrix odour activity from a partition coefficient.**
3. **★ A PANEL-TYPE FACTOR OF 5, MEASURED ON THE SAME COMPOUND IN THE SAME STUDY.** Consumers
   (n = 156) sit **5× above** the trained panel. `k3` §A.7.4's mandatory qualification already warns
   that *"Criteria differ: Brewer 50 % forced-choice; Vega 75 % detection UNCORRECTED for chance;
   Tian 50 % 3-AFC"*. **This paper puts a number on panel expertise alone: 5×.** That is the same
   order as most of the new matrix ratios in this section — i.e. **panel type is as large an effect
   as matrix for these compounds.** Any threshold row without a `panel_type` field is under-specified.

## D.4 A FOURTH SET — and it points straight at an OPEN repo gap

**Tian, H., Yu, B., Yu, H. & Chen, C. (2020)**, *"Evaluation of the synergistic olfactory effects of
diacetyl, acetaldehyde, and acetoin in a yogurt matrix using odor threshold, aroma intensity, and
electronic nose analyses"*, ***J. Dairy Sci.* 103(9):7957–7967**, `10.3168/jds.2019-17495` —
Crossref-verified **[META]**. Verbatim abstract **[ABS]**:

> "the odor thresholds of diacetyl, acetaldehyde, and acetoin **in the yogurt matrix** were
> **5.43, 15.4, and 29.0 mg/L**, respectively, which were **significantly higher than the
> corresponding values in water**."

**⚠️ I COULD NOT VERIFY WHETHER THE WATER VALUES WERE MEASURED IN THIS STUDY OR CITED.** The
abstract says *"the corresponding values in water"* without saying whose. **Until that is checked
this is a LEAD, not a pair.** I am flagging it rather than counting it.

> ★ **BUT NOTE THE AUTHOR AND THE JOURNAL.** `k3` §C.18 lists as OPEN: *"**Tian et al. 2019**
> (`10.3168/jds.2019-16796`) — the only way to settle the `?/kg` unit; without it **seven milk
> thresholds carry a factor-of-1000 basis risk**"* **[REPO]**. This is **the same corresponding
> author (Tian Huaixiang), the same journal, and an adjacent submission number** — 2019-17495 vs
> 2019-16796, same submission year. **Its Methods section will almost certainly describe the same
> threshold protocol and the same units.** ⇒ **Ordering `10.3168/jds.2019-17495` may settle the
> factor-of-1000 milk-threshold basis risk as a side effect**, and it is a *J. Dairy Sci.* paper
> (that journal's back catalogue is broadly open). **This is the cheapest route to a gap `k3` has
> carried as OPEN.**

## D.5 What round 4 did NOT find, reported as a gap

- **No Korean or Japanese same-panel water/matrix pair beyond Hong 2020.** Searched
  `ABSTRACT:("odor threshold" …) AND ABSTRACT:("in water") AND ABSTRACT:(soy OR sake OR tea OR rice
  OR miso OR soybean OR kimchi OR seaweed)` → 14 hits, none a paired-matrix threshold study. **[NEG]**
- **No beer or cheese paired set.** The ASTM E679 beer tradition (Meilgaard-style hop thresholds in
  unhopped pale ale) measures **in beer only**, with no water arm in the same study. **[NEG]**
- **One unverified adjacent lead:** *"Best estimated taste detection threshold for cardamom
  (Elettaria cardamomum Maton) aroma **in different media**"*, Senthil & Bhat, *J. Sensory Studies*
  **26**(1):48–53, `10.1111/j.1745-459x.2010.00320.x` — Crossref-verified metadata only **[META]**,
  **content unread**. A spice oleoresin, so probably of low value to the repo.

---

# §E. RANKED DOWNLOAD LIST

**Every DOI, title and author list below was resolved through `api.crossref.org/works/<doi>` in this
session.** Nothing here is reconstructed from memory.

## TIER 1 — order these; each one decides a null

| # | DOI | exact citation | why, in one line |
|---|---|---|---|
| **1** | `10.1016/j.foodchem.2017.07.131` | **Hamzalıoğlu, A. & Gökmen, V.**, *"Investigation and kinetic evaluation of the reactions of hydroxymethylfurfural with amino and thiol groups of amino acids"*, **Food Chemistry 240:354–360 (2018)** | **★ FALSIFIES NULL 1 AND CLOSES ROUND-3 LEAD L1 IN ONE PAPER.** Aldehyde + Cys/Arg/Lys, **5/25/50 °C**, Michael adduct **and** Schiff base confirmed by HRMS, pseudo-first-order k, Arrhenius obeyed, high- **and** low-moisture arms. **Also the only published HMF-sink kinetics, a node the repo has none of.** Closed access |
| **2** | `10.1016/j.lwt.2025.117469` | **Zhou, X.; Diao, Y.; Zhang, J.; Yu, X.; Wang, G.; Gu, Y.; Ren, D.; Li, S.; Dong, W.; Yi, L.**, *"Sulfury/roasty fading indicators in roasted coffees…"*, **LWT 218:117469 (2025)** | **★ THE ONLY LIVE FALSIFIER OF NULL 2.** First-order + **Arrhenius** on **2-furfurylthiol** across storage temperatures. **GOLD OPEN ACCESS — open it in a browser.** ⚠️ Expect a *storage* Ea; `k3` §C.17 says that is a different channel from the processing sink |
| **3** | `10.1021/bk-1994-0564.ch010` | **Chan, F. & Reineccius, G. A.**, *"Kinetics of the Formation of Methional, Dimethyl Disulfide, and 2-Acetylthiophene via the Maillard Reaction"*, **ACS Symp. Ser. 564 (*Sulfur Compounds in Foods*), pp. 127–137 (1994)** | **★ THE ONLY CANDIDATE SULFUR TEMPERATURE LADDER.** Aqueous, **75–115 °C**, **pH 6/7/8**, 5 min–7.5 h. ⚠️ **Verify absolute quantification and the number of temperatures BEFORE recording it as a ladder** |
| **4** | `10.3390/foods5020035` | **Perry, D. M. & Hayes, J. E.**, *"Effects of Matrix Composition on Detection Threshold Estimates for Methyl Anthranilate and 2-Aminoacetophenone"*, ***Foods* 5(2):35 (2016)** | **★ FALSIFIES NULL 4. ALREADY READ IN FULL — §D.1 IS THE EXTRACTION.** Order only to hold the PDF; no new information will come from it |
| **5** | `10.3168/jds.2019-17495` | **Tian, H.; Yu, B.; Yu, H.; Chen, C.**, *"Evaluation of the synergistic olfactory effects of diacetyl, acetaldehyde, and acetoin in a yogurt matrix…"*, ***J. Dairy Sci.* 103(9):7957–7967 (2020)** | Yogurt thresholds **+** possibly a same-study water arm, **and** — the real prize — **the same lab, journal and submission year as the Tian 2019 paper `k3` §C.18 needs to settle the milk `?/kg` factor-of-1000 basis risk** |

## TIER 2 — order to put numbers on the Ea ladder of §A.1

| # | DOI | exact citation | why |
|---|---|---|---|
| 6 | `10.1021/jf305007y` | **Zamora, R.; Alcón, E.; Hidalgo, F. J.**, *JAFC* **61**(43):10231–10237 (2013) | **The Ea LADDER across four lipid-carbonyl classes: 55–64 / 58–67 / 28–38 kJ/mol.** The single most informative row for whether the covalent aldehyde sink can matter at process temperature. Closed |
| 7 | `10.1111/j.1365-2621.1993.tb04352.x` | **Hidalgo, F. J. & Zamora, R.**, *J. Food Sci.* **58**(3):667–670 (1993) | **The only true ALDEHYDE + LYSINE Arrhenius study found.** Ea **66.5** (colour) / **50** (fluorescence). ⚠️ measures browning and fluorescence, not adduct concentration. Closed |
| 8 | `10.1016/j.foodchem.2010.03.016` | **Hidalgo, F. J.; Delgado, R. M.; Zamora, R.**, *Food Chem.* **122**(3):596–601 (2010) | **Ea = 28–30 kJ/mol for a THIOL + Michael-acceptor addition** (benzyl mercaptan / N-acetylcysteine + acrylamide). The nearest measured analogue in existence to `k3` §A.4's thiol sink. Closed |
| 9 | `10.1021/jf903378x` | **Zamora, R.; Delgado, R. M.; Hidalgo, F. J.**, *JAFC* **58**(3):1708–1713 (2010) | The **forward-vs-reverse** Ea comparison for a Michael adduct — the reversibility datum `meynier2004_extraction.md` X-5 says the corpus has **ZERO** of. Closed |
| 10 | `10.1016/j.foodchem.2019.05.175` | **Sun, Z.; Yang, N.; Liu, C.; Linforth, R. S. T.; Zhang, X.; Fisk, I. D.**, *Food Chem.* **295**:449–455 (2019) | *"investigated at different pH and temperatures"* for **2-FFT in real coffee brew**. **Green OA at Repository@Nottingham but Cloudflare-blocked from this session.** Unknown whether ≥3 temperatures |

## TIER 3 — leads, NOT verified enough to order blind

| # | identifier | what I know | what must happen first |
|---|---|---|---|
| L1 | `10.1016/j.foodres.2022.112187` — **Liu, J.; Poojary, M. M.; Thygesen, M. B.; Andersen, M. L.; Lund, M. N.**, *Food Res. Int.* **163**:112187 (2023) | Stopped-flow Arrhenius, **15–45 °C**, quinone + Lys/Nᵅ-Ac-Lys/Nᵋ-Ac-Lys. Ea **ordering** published, **values not in the abstract** | ⚠️ **A quinone is not an aldehyde.** Order only as an analogue bound, and label it as one |
| L2 | `10.1038/s41538-026-00874-9` — *npj Sci. Food* (2026), aldehydes + myofibrillar protein under stage heating | Proteomics shows saturated aldehydes → **lysine ε-NH₂ Schiff bases**; unsaturated → **both** Schiff base and Michael addition on Cys/His/Trp | ⚠️ its *"thermodynamic parameters"* are **binding** ΔH/ΔS, **not an Ea**. Order for topology only |
| L3 | Chinese-language thesis by **Zhenchun Sun** (Jiangnan Univ., supervisor Xiaoming Zhang) on FFT stability | Inferred from his two FFT-stability papers. **CNKI unreachable this session** | ⚠️ **I verified NOTHING about it.** Do not cite until the record is in hand |
| L4 | A **University of Minnesota** thesis by **F. Chan** underlying download #3 | conservancy.umn.edu hosts UMN theses; a thesis would carry the tables an 11-page chapter compresses | ⚠️ **I did not find the record.** Search the Conservancy before assuming it exists |
| L5 | `10.1111/j.1745-459x.2010.00320.x` — Senthil & Bhat, *J. Sensory Studies* **26**(1):48–53 | *"Best estimated taste detection threshold for cardamom aroma **in different media**"*. Metadata verified | ⚠️ **content unread**; a spice oleoresin, probably low value |

## ALREADY ON DISK OR ALREADY SETTLED — do not re-open

| item | status |
|---|---|
| `scratchpad/zamora2015.pdf` — **Zamora, Navarro, Aguilar & Hidalgo 2015**, *Food Chem.* **174**:89–96 | **★ DOWNLOADED AND READ THIS SESSION. §A.3 is the extraction.** Two new hexanal routes (2-octenal 18 %, 2,4-decadienal 11.5 %) with Ea ≈ 21–30 kJ/mol at **120/160/200 °C**, aqueous, internal standard |
| `data/articles/blank2002.pdf` | **[REPO] [NEG]** — the 22/37 °C FFT pair is already dossiered **and already refused** for a concentration confound. Round 4 adds nothing |
| `data/articles/hofmann1998.pdf` | **[REPO] [NEG]** — the "temperature series" its abstract promises **does not exist**; already verified by `hofmann1998_reconciliation.md`. ⚠️ **The abstract is a snare for keyword search** |
| `Kang2026.pdf` | **[REPO]** — a 100/120/140 °C sulfur ladder **exists as an experiment** but its numeric table is not in the PDF. **The SI is the retrieval, not a new paper** |

## DO NOT ORDER — verified negatives from this round

| item | why |
|---|---|
| `10.1016/j.freeradbiomed.2021.03.040` (Davies group), `10.1021/tx025590o` (Doorn & Petersen), `10.1021/bi060535q` (HNE + HSA) | **The entire biomedical alkenal–protein kinetics literature is single-temperature.** Excellent k₂ sets, **zero** activation energies. §A.2 |
| `10.1021/ac501354y` (Kennedy-Darling & Smith 2014) | Ea = 47 kJ/mol is for formaldehyde crosslink **REVERSAL**, on protein–**DNA**, 4–47 °C. The authors state **no formation constant**. Filing it as an adduct-formation barrier would invert the paper |
| `10.1016/j.jfda.2015.02.001` (Weerawatanakorn/Ho 2015 review) | **[FT-web] [NEG]** — the standard review of flavour-compound stability contains **no rate constant, no half-life and no Ea for any thiol.** Its only number is Blank 2002's 20 %/90 % contrast |
| `10.3389/fnut.2023.1196816` (Ye et al. 2023 FFT review) | **[FT-web] [NEG]** — a dedicated OA review of FFT with a staling-pathway section and **zero kinetic parameters** |
| `10.1016/j.foodchem.2023.137473` (4MBQ + β-lactoglobulin) | k₂ at **25 °C only** — no temperature axis. Order L1 instead if the quinone analogue is wanted |

---

# §F. DECLARED GAPS FROM THIS ROUND — for the `k3` §C register

**F.1 — `k3` §C.1's sentence "*No activation energy for thiol consumption exists anywhere*" is
UNSAFE AS WRITTEN and should be scoped before it is quoted again.** One paper
(`10.1016/j.lwt.2025.117469`) explicitly fits Arrhenius to 2-FFT loss, and one
(`10.1016/j.foodchem.2010.03.016`) reports **Ea = 28–30 kJ/mol** for a thiol adding to a Michael
acceptor. **Neither is a processing-temperature thiol sink.** Proposed rewording: *"No activation
energy for thiol consumption **at Maillard processing temperatures (≥100 °C)** exists in the
literature reviewed. One accelerated-storage Arrhenius fit for 2-FFT exists (Zhou 2025) and, per
`k3` §C.17, storage and processing are different channels."*

**F.2 — The adduct-formation Ea null is falsified in the AMINO-ACID lane and survives in the PROTEIN
lane, and the repo needs both statements, not one.** Measured carbonyl–amine Ea span **22–81 kJ/mol**
across seven Zamora–Hidalgo papers **[ABS]**; **no** Ea exists for any aldehyde adducting a protein.
**Any repo Ea for the covalent aldehyde–protein sink remains the repo's own assumption** — but it is
now an assumption with a *measured neighbourhood*, and that neighbourhood does **not** support the
≳70 kJ/mol case `anantharamkrishnan2020b` says the sink needs to matter.

**F.3 — Round 4 found NO absolute, multi-temperature MFT or FFT dataset.** `k3` §C.3's *"No
temperature dependence for MFT"* **survives intact**. The only ladder that exists for the sulfur lane
is Chan & Reineccius 1994 and it is for **methional**, not a furanthiol — and its quantification
basis is unverified.

**F.4 — `k3` §C.2's "no same-method matrix-vs-water pair exists" is now falsified four times over
and should be retired.** Six such pairs now exist (Hong 2020, Leksrisompong 2010, Perry & Hayes 2016,
Eisele & Semon 2005, Ziegler 2019, and provisionally Tian 2020). **Every new ratio falls in
1.4×–12×.** ⇒ **The 900×–6 714× matrix ratios in `k3` §A.7.5 are, on the weight of all same-panel
evidence now available, artefacts of cross-study method mismatch and of Brewer's
`dose_added_pre_cook` reclassification — not of matrix chemistry.** The `no_factor_available` lookup
table should stay, but its *tail* rows now look much weaker than its *body*.

**F.5 — A NEW measured quantity for the threshold table that the repo does not model: panel type is
worth ~5×.** Ziegler 2019 measured a **consumer** threshold of **14.7 µg/L** (n = 156) and states it
*"exceeded the threshold of the trained panel by a **factor of five**"* **[ABS]**. (⚠️ the trained
value of ≈2.9 µg/L is **my arithmetic from the stated 5×**, not printed in the abstract.) That is the
same magnitude as most matrix shifts.
**A `panel_type` field belongs beside `criterion` in the threshold schema.**

**F.6 — Two items in this dossier are unread and are labelled as such on their face:** download #2
(Elsevier 403 despite gold OA) and Tier-2 #10 (Cloudflare 403 despite green OA). **Nothing else in
this file is unread beyond its stated verification code.** No number appears here that I did not
either read myself or copy from an abstract I read, and every such number carries **[ABS]**.
