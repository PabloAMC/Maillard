# Literature search 2026-10: measured kinetics for the three sulfur-lane path gaps

Run 2026-10-28. This search targets the three path-level errors in [SULFUR_FLUX_MAP.md](SULFUR_FLUX_MAP.md) section 5:
(1) where cysteine's Strecker sulfur goes, (2) the pentose enolisation pH law at low pH, and (3) the temperature dependence
of furfural + H2S → FFT. It records the search only. Nothing was added to `data/`, and no paper was downloaded.
It follows [LITERATURE_SEARCH_2026-10.md](LITERATURE_SEARCH_2026-10.md) (rounds 1-2),
[LITERATURE_SEARCH_2026-10_ROUND3.md](LITERATURE_SEARCH_2026-10_ROUND3.md) (G1-G10) and
[LITERATURE_SEARCH_2026-10_VARIABLES.md](LITERATURE_SEARCH_2026-10_VARIABLES.md) (V6, buffer catalysis). Papers those
documents or the corpus already hold appear here only as cross-references.

## Bottom line

- **Gap 1, the Strecker sulfur split.** I found no measurement of the mercaptoacetaldehyde → H2S fraction, its rate, or its
  temperature or pH dependence. **The one new quantitative source is Kobayasi & Fujimaki 1965**, a two-page note titled for
  exactly this question: mercaptoacetaldehyde, H2S and acetaldehyde from cysteine boiled with carbonyl compounds. A
  search-engine snippet of its full text gives 0.42 and 0.80 mmol H2S from 10 mmol cysteine with 10 and 20 mmol of carbonyl,
  i.e. **about 4 % and 8 % of the cysteine**. I have not verified this, and the snippet does not say which carbonyl was used.
  The rest of what bears on gap 1 is already on the round-3 list: Edwards & Wedzicha 1992, Mizutani 1960 and de Roos 1992.
- **Gap 2, the low-pH enolisation law.** No paper gives the 1- vs 3-deoxypentosone partition below pH 5. Three new sources
  constrain the pieces around it:
  - **Larralde, Robertson & Miller 1995 (PNAS, open access)**: ribose loss at pH 4-8 and 40-120 °C, with a measured
    phosphate-catalysis law.
  - **Lu & Mosier 2008** (conference abstract): xylose loss is *slowest at pH 2.2*, is about 9x faster at pH 4.2, and
    gives mainly furfural only at pH ≤ 3.
  - **Mikami et al. 2017 (open access)**: norfuraneol yields from pentose + lysine in 0.2 M phosphate at pH 5.5-8 and
    100 °C, about 10 % w/w at pH 7.

  Together these argue that at pH 3 the specific-acid term is close to its floor and the uncatalysed and buffer terms set the
  rate. That is inference, explained below.
- **Gap 3, the temperature dependence of FFT.** I found no measured Ea or multi-temperature rate constants for furfural
  + H2S → FFT in water (about 75 % confident that none exists in the open literature at 90-150 °C). The new item is
  path-level: **Poisson et al. 2009** used labelled precursors inside coffee beans and *ruled out the furfural route to FFT*
  that model studies propose. The temperature data that do exist are already in the corpus: Yaghmur 2005 (40-70 °C), the
  Kang 2026 SI (100/120/140 °C) and Yiltirak 2026 (100-130 °C). **The quickest fix is an experiment, not a paper** (end of
  the gap 3 section).

## Method

- **Tools.** WebSearch (standard and extended) and WebFetch on abstract or open-access pages: J-STAGE, PMC, NTRS, AGRIS
  and a conference-abstract server. Crossref `query.bibliographic` (anonymous, no e-mail) resolved titles and DOIs. OpenAlex
  and Semantic Scholar returned HTTP 429 for this IP (the OpenAlex budget was spent in earlier rounds), so this round has no
  citation chaining. That is the main caveat: pre-1990 papers without abstracts are under-found.
- **Corpus check first.** Every candidate was checked against `data/keys/papers.yml`, the `data/articles/` file names, the
  extraction dossiers in `data/lit/extraction_dossiers/`, and the three earlier search documents.
- **Confidence** below is my probability that the paper contains the number the gap needs, in a form the engine can use. It
  is not the probability that the paper is good. **[abstract]** marks what I read in an abstract or full-text HTML.
  **[snippet]** marks what came only from a search-engine excerpt of the full text. **[inferred]** marks what I infer from the
  title, the venue or other papers.

---

## Gap 1: cysteine's Strecker sulfur (mercaptoacetaldehyde vs H2S vs ring products)

**What the corpus and earlier documents already cover.** Use these first:

- Shu 1985 (`shu1985`): cysteine thermolysis at pH 2.2/5.1/7.1, with the mercaptoacetaldehyde self-condensation products.
- Zheng & Ho 1994 (`zheng1994`): H2S release kinetics from cysteine and GSH vs pH, by sulfide electrode.
- Du 2023 (`du2023`): 84 % of the thiophene C2 unit comes from cysteine-derived acetaldehyde.
- Elmore 2002 (`elmore2002`): a side experiment with mercaptoacetaldehyde (from 1,4-dithiane-2,5-diol) + enals at
  140 °C, amounts not stated.
- Mottram & Whitfield 1995b (`mottram1995b`) and Zheng 1997 (`zheng1997`): thiirane and thiophenes from mercaptoacetaldehyde.
- Round 3 G6/G3/G2: Edwards & Wedzicha 1992 (DOWNLOAD #2), Mizutani 1960 (#6), de Roos 1992 (#7). Round 2:
  Pripis-Nicolau 2000 (10.1021/jf991024w). Round 2 SKIP: Parker, Bishara, Baines & Mottram 2014.

| paper | DOI / access | what it measures | serves | confidence |
|---|---|---|---|---|
| **Kobayasi, N. & Fujimaki, M. (1965).** "On the formation of mercaptoacetaldehyde, hydrogen sulfide and acetaldehyde on boiling cysteine with carbonyl compounds." *Agric. Biol. Chem.* 29(7):698-699 | [10.1271/bbb1961.29.698](https://doi.org/10.1271/bbb1961.29.698) (also 10.1080/00021369.1965.10858450). **Free on J-STAGE** (128 KB PDF) | [abstract] no abstract; the title names all three products. [snippet] 10 mmol cysteine boiled with 10 or 20 mmol carbonyl gave 14.3 mg (0.42 mmol) and 27.2 mg (0.80 mmol) H2S. The proposed mechanism is cysteine → mercaptoacetaldehyde by reaction with an α-diketone. Boiling (about 100 °C), one temperature; pH not seen | gap 1: an **H2S yield per cysteine vs dicarbonyl dose**, and possibly the mercaptoacetaldehyde : acetaldehyde split | H2S yield vs carbonyl dose: **0.75**. A mercaptoacetaldehyde/H2S partition: 0.35. Any temperature or pH dependence: 0.05 |
| Kato, Kurata & Fujimaki (1973). "Volatile compounds produced by the reaction of L-cysteine or L-cystine with carbonyl compounds." *Agric. Biol. Chem.* 37(3):539-544 | [10.1271/bbb1961.37.539](https://doi.org/10.1271/bbb1961.37.539), free on J-STAGE | [abstract] cysteine/cystine + glucose or pyruvaldehyde, "various temperatures" for sensory, then GC-MS identification at 160 °C. Identification only | gap 1 (products only) | 0.05. **SKIP** |
| Parker, Bishara, Baines & Mottram (2014). "Comparison of ribose and ascorbic acid in model process reactions containing cysteine." *Flavour Science* (Elsevier) 211-214 | [10.1016/b978-0-12-398549-1.00040-4](https://doi.org/10.1016/b978-0-12-398549-1.00040-4), paywalled chapter | [abstract] ascorbic acid + Cys gives more MFT than ribose + Cys. The authors propose that a fairly stable thiazolidine holds the cysteine in the ribose system, "limiting the formation of H2S" and NF → MFT | gap 1 (mechanism: thiazolidine vs Strecker as the cysteine sink) | 0.15. Round 2 already SKIPped it; the abstract gives a mechanism claim and no rates |

**What this means for B59.**

- The Strecker split is still unmeasured. Kobayasi & Fujimaki is the only primary number found: **about 4-8 % of cysteine as
  H2S at boiling with an α-dicarbonyl**, if the snippet is right.
  - That number is an *upper bound on nothing and a lower bound on little*. It is total H2S from all routes (thermolysis,
    thiazolidine release and Strecker together), with the H2S lost to the open boiling vessel.
  - Carry it as a declared-assumption band, not a fitted value.
- **[inferred]** The textbook step "mercaptoacetaldehyde → acetaldehyde + H2S" is formally a reduction. In the papers I read,
  its mechanism is drawn but not established. The literature schemes put the H2S/acetaldehyde exit either at the
  decarboxylated Schiff base (elimination before hydrolysis) or after the free aldehyde forms.
  - Which one holds changes the model structure: a branch at `r_dpo_ddp` vs a downstream first-order decay of a
    mercaptoacetaldehyde pool.
  - Kobayasi's note proposes the latter. Nothing found measures the rate of either.
- Hofmann's pH sweep remains the best available constraint, as a pre-registered shape hold-out, on a split carried as a
  free branching fraction with a declared band.

---

## Gap 2: pentose 1,2- vs 2,3-enolisation at low pH; phosphate catalysis; norfuraneol

**What the corpus and earlier documents already cover.**

- **Deoxyosone partition vs pH**: Nedvidek 1992 (`nedvidek1992`): α-dicarbonyls trapped as quinoxalines from xylose at pH
  5, 6.5 and 7. This is the only partition measurement, and none of it is below pH 5.
- **Thiol pH ladders**:
  - Meynier & Mottram 1995: pH 4.5-6.5 in pyrophosphate.
  - Cerny & Briffod 2007: pH 4.0-7.0, labelled xylose.
  - Madruga & Mottram 1995.
- **Fed norfuraneol**: Whitfield & Mottram 1999 (pH 4.5) and 2001 (pH 6.5).
- **Low-pH MFT from the C2+C3 route**: Hofmann & Schieberle 1998 Table 8 (hydroxyacetaldehyde + mercapto-2-propanone, pH
  series).
- **Norfuraneol and furaneol in phosphate at pH 5/6/7**: Blank 1997.
- **Norfuraneol stability and browning**: Nakamura 2020/21.
- **Buffer catalysis**: Rizzi 2004, Bell 1997, Potman 1989 and Davidek 2002 are in the corpus. VARIABLES V6 lists Saunders
  1966, Huang 2011, Pepa 2023 and Reynolds 1959.

| paper | DOI / access | what it measures | serves | confidence |
|---|---|---|---|---|
| **Larralde, R., Robertson, M. P. & Miller, S. L. (1995).** "Rates of decomposition of ribose and other sugars: implications for chemical evolution." *PNAS* 92(18):8158-8160 | [10.1073/pnas.92.18.8158](https://doi.org/10.1073/pnas.92.18.8158). **Open access** (PMC41115; NASA NTRS 19980033941, a scanned PDF) | [abstract] ribose loss by 1H NMR at **pH 4-8, 40-120 °C**; t½ = 73 min at pH 7.0 / 100 °C. Other aldopentoses, aldohexoses and ribose 5-phosphate are within 10x of ribose. "Considerable general base catalysis." [snippet] at pD 7.4, 100 °C, an empirical law k = 2.86e-3 + 9.36e-2·([DPO4²⁻] + [D2PO4⁻]), with the printed unit garbled by OCR. Only **min⁻¹** reproduces the 73-min half-life at their 0.05 M buffer (k = 7.5e-3 min⁻¹ gives t½ = 92 min), so read it as min⁻¹ | gap 2: **a measured phosphate-catalysis law for pentose loss**, and a pH 4-8 × T profile with Arrhenius data | phosphate law and pH/T profile: **0.85**. Deoxyosone partition: 0.05 (it measures sugar loss, not products) |
| **Lu, Y. & Mosier, N. S. (2008).** "pH and buffer effects on xylose degradation rates and products." Poster, SIM 30th Symposium on Biotechnology for Fuels and Chemicals | no DOI. The abstract is online (sim.confex.com P5869); no journal version found | [abstract] **pH 0.5-7.0**, McIlvaine (phosphate/citrate) buffer, HCl/NaCl below pH 1.5; temperature not stated. **Xylose loss is slowest at pH 2.2, about 9x faster at pH 4.2, and complete within 5 min at pH 7.0. Furfural is the main product only at pH 0.5-3.0**; at pH 4.5-7.0 the products are "complex" and polymeric. Added phosphate does not change the rate at pH 1.0 | gaps 2 and 3: **the shape of the pentose pH-rate profile and the furfural share vs pH**, without amine | that this is the profile: 0.8. That usable numbers can be found: 0.3 (look for Y. Lu's Purdue thesis) |
| **Mikami, Y., Nakamura, M., Yamada, S. & Murata, M. (2017).** "4-Hydroxy-5-methyl-3(2H)-furanone (HMFO) contributes to browning in the xylose-lysine Maillard reaction system." *Food Sci. Technol. Res.* 23(2):283-289 | [10.3136/fstr.23.283](https://doi.org/10.3136/fstr.23.283). **Open access** (J-STAGE, CC BY-NC-SA) | [full-text HTML] 13.3 mM pentose + 34 mM Lys in **0.2 M phosphate, pH 5.5-8.0, 100 °C**. HMFO (= norfuraneol) peaks at **pH 6.5-7** and is at its maximum at 1 h at pH 7, then falls; it forms and decays more slowly at pH 6. Ribose is the best precursor: about 30 mg/100 mL HMFO from 100 mg/100 mL ribose. About 10 % w/w from xylose. HMFO decay ± Lys at pH 6/7/8, 100 °C, 0-6 h. The values are in figures only. Acetate pH 3.5-6 is listed in the methods, but no HMFO data are given for it | gap 2: **a norfuraneol yield and decay vs pH at 100 °C in phosphate**, the upper half of the pH law. Not pH 3 | norfuraneol vs pH 5.5-8 (figure digitisation): **0.9**. Anything at pH < 5.5: 0.05 |
| Garrett, E. R. & Dvorchik, B. H. (1969). "Kinetics and mechanisms of the acid degradation of the aldopentoses to furfural." *J. Pharm. Sci.* 58(7):813-820 | [10.1002/jps.2600580703](https://doi.org/10.1002/jps.2600580703), paywalled (Wiley) | [inferred from title and venue] Garrett's pharmaceutical-kinetics papers usually report pH-rate profiles and Arrhenius parameters. Expect pentose loss and furfural appearance vs [H+] and temperature; likely mostly below pH 2 | gaps 2 and 3: **the specific-acid term of pentose → furfural and its Ea** (the model's `r_pent_tdp` acid factor and `r_tdp_fur`) | acid rate law + Ea: 0.6. Data reaching pH 3-5: 0.35 |
| Li, M.-C., Ko, W.-C. & Hsieh, C.-W. Kinetics of ribose/xylose/glucose + cysteine at 95 °C (Da-Yeh University, Taiwan; venue and year not resolved) | Open PDF on people.dyu.edu.tw/paper/354509_e.pdf, **unreachable at search time** (connection refused) | [snippet] pentoses or glucose + Cys, **95 °C, 120 min, pH 5/7/8 (or 9)**: sugar loss ribose > xylose > glucose, with kinetic parameters for reactant loss. Amadori content highest at pH 7; MFT forms more readily at pH ≤ 7 | gap 2 (sugar loss vs pH with cysteine, at a temperature below Hofmann's) | sugar-loss rate constants vs pH: 0.5. Product partition: 0.1 |
| Severin, T. & Seilmeier, W. (1967, 1968). "Studien zur Maillard-Reaktion" (parts; subtitles not resolved). *Z. Lebensm. Unters. Forsch.* 134:230-232 and 137:4-6 | [10.1007/bf01787086](https://doi.org/10.1007/bf01787086), [10.1007/bf01842036](https://doi.org/10.1007/bf01842036); Springer, paywalled, German | [inferred] Severin's group made norfuraneol from pentoses with amine salts in mildly acidic media; these are the likely primary reports. Product yields vs amine and pH are possible | gap 2: **norfuraneol from pentoses below pH 5**, amine-catalysed | 0.25 |
| Hicks, K. B. & Feather, M. S. (1975). "Mechanism of formation of 4-hydroxy-5-methyl-3(2H)-furanone, a component of beef flavor, from Amadori products." *J. Agric. Food Chem.* 23:957-960 | [10.1021/jf60201a005](https://doi.org/10.1021/jf60201a005), ACS | [inferred] mechanism, probably by deuterium/13C labelling, of norfuraneol from pentose-type Amadori compounds; 1-deoxyosone vs other routes | gap 2 (the 2,3-enolisation route of the Amadori compound to norfuraneol, the B59 de Roos path) | rates: 0.1. Mechanism/route shares: 0.5 |
| Tonsbeek, C. H. T. et al. (1969; co-authors not verified). "Components contributing to beef flavor: natural precursors of 4-hydroxy-5-methyl-3(2H)-furanone in beef broth." *J. Agric. Food Chem.* 17:397-400 | [10.1021/jf60162a022](https://doi.org/10.1021/jf60162a022), ACS | [inferred] ribose 5-phosphate (from IMP) + amine → norfuraneol. The phosphate ester's leaving group is a phosphate-assisted 1-deoxy route that free ribose lacks | gap 2 (an alternative NF source; matters for meat and IMP-containing matrices, not for Hofmann's pots) | 0.2 |

**Inference for the model (my reasoning, not a measurement; please check):**

1. **The pH 3 arm is near the pentose pH-rate minimum.**
   - Lu & Mosier put the xylose minimum at pH 2.2, with the rate rising about 9x by pH 4.2.
   - Larralde's ribose rates rise from pH 4 to 8 by general-base catalysis.
   - So between about pH 2 and 4, *neither* specific-acid nor specific-base catalysis dominates pentose loss. A law in which
     the 1,2-enolisation scales with [H+] and the 2,3-enolisation with [OH-] (or the free-amine fraction) gives a ratio that
     moves by about 10x per pH unit right where the data say the acid term is at its floor.
   - That is the 0.4 % DPO share at pH 3. The B60 law may be extrapolating the acid branch through a region where the
     uncatalysed and buffer terms set both branches.
   - **Confidence that this explains a material part of the pH 3 MFT miss: about 0.5.** The alternative is a low-pH MFT
     source the model lacks: the C2+C3 route (Hofmann 1998 T8) or Cys-Amadori 2,3-enolisation.
2. **Phosphate at 0.5 M is not a spectator.**
   - Read in min⁻¹, Larralde's law at 100 °C makes 0.5 M phosphate speed ribose loss about 17x over buffer-free, and about
     7x over their own 0.05 M.
   - Every Hofmann arm is at 0.5 M phosphate. At pH 3 the buffer is about 87 % H2PO4⁻ and 13 % H3PO4 (pKa1 2.15).
   - Larralde's law sums the two phosphate species at pD 7.4. So it cannot say whether H2PO4⁻, the species at pH 3, is as
     active as HPO4²⁻. That needs Rizzi 2004 or Bell 1997 (corpus), or the Larralde full text if it has a pH series.
   - Phosphate as a bifunctional enolisation catalyst would raise the 2,3-branch at low pH relative to a pure [H+] law.
3. **Norfuraneol's own pH window.** Mikami 2017 has NF peaking at pH 6.5-7 and decaying faster at higher pH in 0.2 M
   phosphate at 100 °C. That is an upper-pH constraint on `k_dpo_nf` against NF loss, usable as a shape check.

---

## Gap 3: furfural + H2S → FFT, temperature dependence; FFT time courses; FFT degradation

**What the corpus and earlier documents already cover.**

- Hofmann & Schieberle 1998 (`hofmann1998_reconciliation`): Table 3 (FFT from C5 precursors + H2S) and Table 5 (the FA/FFT
  pair), both at 145 °C only. The dossier shows that the promised temperature series does not exist.
- Yaghmur 2002/2005: furfural loss in the presence of cysteine at 40-70 °C, Ea 46.5 kJ/mol in water, with FFT under 1 %
  of the flux.
- Kang 2026 SI Table S4: MFT, FFT and furfural at 100/120/140 °C from TTCA.
- Yiltirak 2026: MFT/FFT at 100 °C/4 h to 130 °C/0.5 h.
- Liu 2021 (`liu2021`, 10.1016/j.foodres.2021.110295): furfuryl alcohol as an FFT intermediate.
- Cerny & Davidek 2003: FFT keeps the intact ribose skeleton at 95 °C.
- FFT stability: Blank 2002 (Fenton, 22/37 °C) and Guth 1995 (disulfide ⇌ thiol at 100 °C).
- Round 3 G1: no buffer stability data for FFT/MFT at two or more temperatures.
- Rounds 1-2: Shibamoto 1977 (furfural + H2S + NH3). Round 3: Mizutani 1960 (the furfural-cysteine thiazolidine releases
  the most H2S at 120 °C).

| paper | DOI / access | what it measures | serves | confidence |
|---|---|---|---|---|
| **Poisson, L., Schmalzried, F., Davidek, T., Blank, I. & Kerler, J. (2009).** "Study on the role of precursors in coffee flavor formation using in-bean experiments." *J. Agric. Food Chem.* 57(21):9923-9931 | [10.1021/jf901683v](https://doi.org/10.1021/jf901683v), ACS, paywalled | [abstract, via a vendor reprint of it] labelled and unlabelled precursors infused into green beans, then roasted. **"Ruled out the 2-furaldehyde route to FFT as suggested by model studies"**; no evidence that arabinose's C5 skeleton enters FFT | gap 3, **path level**: if FFT does not come from free furfural + H2S in a real matrix, then the model's FFT branch (95-99 % furfural + H2S) may have the wrong route as well as the wrong Ea. That favours a route through a cysteine-carbonyl intermediate (the Mizutani thiazolidine, or Liu 2021's furfuryl alcohol) with a different temperature law | the route claim: 0.9. Kinetics: 0.05. Transfer from coffee bean to aqueous pots: low (in-bean ≠ buffer; Cerny 2003's intact skeleton argues for a pentose-derived C5, which free furfural is) |
| Garrett & Dvorchik 1969 (table above) | see above | pentose → furfural, acid kinetics, probably with Arrhenius data | gap 3 via `r_tdp_fur` / `r_pent_tdp`: if furfural formation's Ea is underestimated, FFT falls too slowly as temperature drops | 0.35 (for data relevant to pH 4-6) |
| Lu & Mosier 2008 (table above) | see above | furfural is the dominant xylose product only at pH ≤ 3 | gap 3: the **furfural share** vs pH; above pH 4.5 most pentose goes elsewhere, which caps FFT | 0.3 for usable numbers |

**Verdict and the cheap way out.**

- I found no aqueous measurement of the temperature dependence of furfural + H2S (or HS⁻) → FFT. I am about 75 % confident
  none exists at 90-150 °C in the open literature. The wine literature (Tominaga 2000, Blanchard 2001) is yeast-mediated and
  at cellar temperature. The biorefinery literature measures furfural *formation* at pH < 2 and above 150 °C.
- Two cheap measurements would close the gap, in a SIDA lab such as the WUR contact:
  - Hofmann 1998 Table 3's own protocol (furfural + Na2S or cysteine, 0.5 M phosphate pH 5, sealed) at 100, 120 and 145 °C.
  - FFT alone, heated at the same temperatures, for its loss.
  - Together these give k(T) for formation and for the FFT-specific sink, a few runs each.
- Until then, **Yaghmur 2005's 46.5 kJ/mol is the only Ea on the furfural side, and it is for furfural loss, not FFT
  formation**. Its 40-70 °C window needs extrapolation, and the flux it measures is dominated by a sink the paper never
  identifies (see its dossier).
- **[inferred]** A low Ea on the FFT branch relative to the norfuraneol branch is the simplest reading of the "excess grows as
  temperature falls" signature. Poisson 2009 also allows a route error. Kang 2026 SI furfural at 100/120/140 °C, against FFT,
  can tell the two apart within the corpus: if FFT/furfural is flat in temperature, the route is fine and the Ea is the issue.

---

## Ranked download list (new papers only; at most 10)

File names follow round 3: lowercase `<firstauthor><year>.pdf`. I checked them against `data/articles/`; there are no
clashes (`poisson2019.pdf` is a different paper).

| # | paper | DOI | gap | access | why |
|---|---|---|---|---|---|
| 1 | Larralde, Robertson & Miller 1995, PNAS | [10.1073/pnas.92.18.8158](https://doi.org/10.1073/pnas.92.18.8158) | 2 | open (PMC41115) | the only measured phosphate-catalysis law for pentose loss, with a pH 4-8 × 40-120 °C profile; decides whether 0.5 M phosphate moves Hofmann's pots (`larralde1995.pdf`) |
| 2 | Kobayasi & Fujimaki 1965, Agric. Biol. Chem. 29:698 | [10.1271/bbb1961.29.698](https://doi.org/10.1271/bbb1961.29.698) | 1 | free (J-STAGE) | the only primary H2S + mercaptoacetaldehyde + acetaldehyde yield from cysteine + carbonyl; two pages (`kobayasi1965.pdf`) |
| 3 | Mikami et al. 2017, Food Sci. Technol. Res. 23:283 | [10.3136/fstr.23.283](https://doi.org/10.3136/fstr.23.283) | 2 | open | norfuraneol formation and decay vs pH 5.5-8 in 0.2 M phosphate at 100 °C; a shape check on `k_dpo_nf` against NF loss (`mikami2017.pdf`) |
| 4 | Garrett & Dvorchik 1969, J. Pharm. Sci. 58:813 | [10.1002/jps.2600580703](https://doi.org/10.1002/jps.2600580703) | 2, 3 | Wiley | the acid term and Ea of pentose → furfural; tests the acid factor and `r_tdp_fur` (`garrett1969.pdf`) |
| 5 | Poisson et al. 2009, JAFC 57:9923 | [10.1021/jf901683v](https://doi.org/10.1021/jf901683v) | 3 | ACS | labelling evidence against furfural + H2S as the FFT route in a real matrix; a path-level check on the FFT branch (`poisson2009.pdf`) |
| 6 | Li, Ko & Hsieh, Da-Yeh Univ., pentoses/glucose + Cys kinetics, 95 °C, pH 5/7/8-9 | none found | 2 | open PDF (server down at search time) | sugar-loss rate constants vs pH with cysteine at a sub-Hofmann temperature (`li_dyu.pdf`; year unresolved) |
| 7 | Severin & Seilmeier 1968, Z. Lebensm. Unters. Forsch. 137:4 (with the 1967 part, 134:230) | [10.1007/bf01842036](https://doi.org/10.1007/bf01842036) | 2 | Springer | the likely primary reports of norfuraneol from pentoses with amine salts in acidic media, i.e. a low-pH NF route (`severin1968.pdf`) |
| 8 | Hicks & Feather 1975, JAFC 23:957 | [10.1021/jf60201a005](https://doi.org/10.1021/jf60201a005) | 2 | ACS | the route from Amadori compound to norfuraneol (for B59's Cys-Amadori 2,3-enolisation) (`hicks1975.pdf`) |
| 9 | Tonsbeek et al. 1969, JAFC 17:397 | [10.1021/jf60162a022](https://doi.org/10.1021/jf60162a022) | 2 | ACS | ribose 5-phosphate as the NF precursor: phosphate as a leaving group, relevant to IMP matrices (`tonsbeek1969.pdf`) |

**Not a downloadable paper, but worth a lookup:** Y. Lu (Purdue, Mosier lab) PhD thesis, about 2009-2011, for the full xylose
pH 0.5-7 rate and product data behind the 2008 abstract.

**Already on the round-3 list and directly relevant here; fetch them in the same batch:**

- Edwards & Wedzicha 1992 and Wedzicha 1991: cysteine + 3-deoxyhexosulose, fast capture then slow Strecker release
  (gap 1).
- Mizutani 1960: H2S from cysteine-aldehyde thiazolidines at 120 °C; the furfural adduct releases the most (gaps 1 and 3).
- de Roos 1992: reactivity scales with the neutral thiazolidine, i.e. the pH law of the TTCA step (gap 2).

**The most decisive item per gap:**

- Gap 1: Kobayasi & Fujimaki 1965 (new). Edwards & Wedzicha 1992 is the best kinetic partner.
- Gap 2: Larralde 1995.
- Gap 3: no paper settles it. The best new item is Poisson 2009, at path level only. The decisive move is the
  three-temperature Hofmann-T3 experiment above, plus the within-corpus FFT/furfural check on the Kang 2026 SI.
