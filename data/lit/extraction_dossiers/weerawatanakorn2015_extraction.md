# Weerawatanakorn, Wu, Pan & Ho 2015 — EXTRACTION (review of flavor-compound stability; the MFT/FFT, 2-acetyl-2-thiazoline, methional and furaneol sections, with their primary sources)
### A secondary source. Nothing here is primary evidence; its value is the map of which primary papers carry the thiol-stability numbers and which of those are not yet on disk.

**Source on disk:** `data/articles/weerawatanakorn2015.pdf` (owner's download, 2026-09-08; 15-page Elsevier
"article in press" PDF, born-digital). Read from the pdftotext extraction in the scratchpad; the text layer
is clean except that the degree sign is rendered "/C14" and, in several threshold statements, the unit
prefix reads "mg" where the original almost certainly prints "µg" (the Symbol-font µ was lost; see §5.1).
Figures 1-9 are structures and mechanism schemes with no numbers (FIGURE-ONLY, nothing lost).

## 0. Identity

| field | value |
|---|---|
| Title | "Reactivity and stability of selected flavor compounds" |
| Authors | Monthana Weerawatanakorn (Naresuan Univ.), Jia-Ching Wu (NCKU), Min-Hsiung Pan (NTU), Chi-Tang Ho* (Rutgers) |
| Venue | Journal of Food and Drug Analysis 2015 (in press at the time of the PDF; vol. 23, pp. 176-190 in print), Review Article |
| DOI | 10.1016/j.jfda.2015.02.001 |
| Scope | citral, allyl isothiocyanate, furaneol/homofuraneol/norfuraneol, vanillin, MFT and FFT, methional, 2-acetyl-1-pyrroline, 6-acetyl-tetrahydropyridine, 2-acetyl-2-thiazoline, 5-acetyl-2,3-dihydro-4H-1,4-thiazine; 110 references |
| Type | narrative review; **no tables**; no numbers of its own |

## 1. Why it matters

The sink re-fit (`results/validation/kinetic_core_b17_prereg.md`) needs to know every way MFT and FFT leave
a pot and which measurement pins each. `docs/validation/thiol_sink_candidates.md` §2 already lists the
primary sources on disk (Kumazawa 2003, Mottram 2002, Zhou 2023, Zhang 2024, Hofmann 2002, Gigl 2021,
Blank 2002, van Seeventer 2001, Anantharamkrishnan 2020b, Starkenmann 2008). This review was read to check
whether any thiol-stability measurement exists in the literature that the corpus has not seen. The answer
is: its MFT/FFT section (§2.5) rests on ten primaries, seven of which are on disk; the three not on disk
are Mottram, Szauman-Szumski & Dodson 1996 (thiol/disulfide loss to egg albumin at 100 °C), Kerscher &
Grosch 1998 (SIDA levels of MFT, FFT and the two mercaptopentanones in heated meat) and Kumazawa et al.
1998 (Japanese, coffee sterilisation). The review adds one number the corpus lacks on 2-acetyl-2-thiazoline
(96 % loss in 60 min at 100 °C in water; 99 % at 145 °C in pH-5 buffer; 2.5 % in oil), from Hofmann &
Schieberle 1995 JAFC 43, 2946, also not on disk. It carries nothing on sotolon and nothing quantitative on
Strecker aldehydes.

## 2. Methods as they matter to a model

Not applicable: a narrative review without a stated search protocol. Every number below is the review's
paraphrase of a primary paper; where the primary is on disk, the dossier for that paper is the authority
and the row here only records what the review says and whether it agrees. The review's reference numbers
are kept in brackets so the fetch list (§6) can be checked against the printed reference list (pp. 13-15).

## 3. Tables re-typed

The review contains no tables. §4 below is the tabulation of its quantitative statements.

## 4. Numbers the repository can use (all evidence class **secondary_review**; use only to locate the primary)

### 4.1 MFT and FFT stability statements

| compound | condition (as stated in the review) | statement / number | primary as cited in the review | on disk? |
|---|---|---|---|---|
| MFT, FFT | diethyl ether, 6 °C | "Oxidation of MFT and FFT in diethyl ether occurs even within 1 day at 6 °C"; products the corresponding disulfides or mixed disulfides; "the oxidation rate of MFT was higher than that of FFT"; "higher temperatures ... increased the MFT oxidation rate" | [84] Hofmann, Schieberle, Grosch, JAFC 1996, 44, 251-255 | yes (`hofmann1996.pdf`, dossier `hofmann1996_extraction.md`; declared "neither", organic solvent) |
| thiol and disulfide flavor compounds | aqueous solution, 100 °C, with egg albumin | "causes a decrease in the concentration of flavor", attributed to thiol/disulfide interchange with protein; "depends on the structure of protein and the number and position of sulfhydryl groups"; no number given | [83] Mottram, Szauman-Szumski, Dodson, JAFC 1996, 44, 2349-2351 | **no** |
| FFT (coffee drink) | canned, sterilised 121 °C / 10 min; pH 3-7 | roasty flavor "significantly decreased"; "Higher pH, particularly in the range of 5.0-7.0, and higher temperatures significantly reduced the FFT concentration"; major volatile product difurfuryl disulfide, then furfural and furfuryl alcohol; non-volatile products "presumably produced through a Fenton-type reaction"; loss attributed to melanoidin binding, Fenton reaction and pH-dependent degradation | [91] Kumazawa & Masuda, JAFC 2003, 51, 2674-2678 | yes (`kumazawa2003_extraction.md` has the residual-FFT ladder 99.5 → 0.1 % from pH 3 to 7) |
| coffee flavor | sterilisation 121 °C / 15 min or 134 °C / 3 min | flavor "changed"; no number | [90] Kumazawa, Masuda, Nishimura, Hiraishi, Nippon Shokuhin Kagaku Kogaku Kaishi 1998, 45, 108-113 | **no** (Japanese) |
| MFT, FFT | cysteine + ribose process flavoring, 50 °C accelerated storage, with and without air | both decrease with time; "MFT was found to be less stable than FFT"; MFT loss "not because of the oxidative pathway ... to form disulfide, but ... oligomerization/polymerization" via protonation at C-2 and attack of the thiol and the furan ring at C-5; "MFT can react with other thiols such as cysteine"; no rate quoted | [87] van Seeventer, Weenen, Winkel, Kerler, JAFC 2001, 49, 4292-4295 | yes (`vanseeventer2001.pdf`; ingested in `k3_final_parameter_inventory.md` and FIT_HOLDOUT_DECLARATION as zero-order 50 °C rows, HOLD-OUT; no standalone dossier) |
| FFT | Fenton system (Fe(II) + H2O2), 37 °C, 1 h | "90 % degradation of FFT within 1 hour at 37 °C"; "20 % loss of FFT at room temperature" after 1 h; products difurfuryl disulfide > bifurfuryl > difurfuryl monosulfide, MW 124-262 Da; Fe/Fe(II) more effective than Mn or Cu; C-centred and S-centred radical intermediates | [86] Blank, Pascual, Devaud, Fay, Stadler, Yeretzian, Goodman, JAFC 2002, 50, 2356-2364 | yes (`blank2002_extraction.md`) |
| FFT (and MFT) | coffee beverage kept warm with coffee melanoidins | FFT "rapidly reduced (50 % after 20 minutes and almost 100 % after 30 minutes)"; FFT "most affected among other coffee thiols"; covalent binding through CROSSPY-derived pyrazinium ions to the thioether 2-(2-furyl)methylthiol-1,4-dihydropyrazine; "Without oxygen, degradation of FFT was very slow" | [94] Hofmann & Schieberle, JAFC "2001;50:319-26" (= JAFC 2002, 50, 319-326); [95] Hofmann, Czerny, Calligaris, Schieberle, JAFC 2001, 49, 2382-2386 | [94] yes (`hofmann2001.pdf`, dossier `hofmann2002_extraction.md`, which records 80 % bound in 30-90 min at 30 °C, 12.5 g/L melanoidin, a plateau); [95] **no** |
| volatile thiols in coffee | general ROS statement only | no number | [93] Charles-Bernard, Roberts, Kraehenbuehl, JAFC 2005, 53, 4426-4433 | yes (`charles-bernard2005.pdf`; cited in `thiol_sink_candidates.md` §1 for the thioether channel; no standalone dossier) |
| MFT, FFT, 3-mercapto-2-pentanone, 2-mercapto-3-pentanone | heated meat | cited only for occurrence and thresholds; the primary is a SIDA quantification in meat | [85] Kerscher & Grosch, JAFC 1998, 46, 1954-1958 | **no** |
| MFT, FFT formation | cysteine + pentose/hexose, NF + cysteine or H2S, thiamine | cited for routes only | [82] Hofmann & Schieberle, JAFC 1998, 46, 235-241 | yes (`hofmann1998.pdf`) |
| MFT, methional | stored orange juice | off-flavors; no rate | [88] Bezman, Rouseff, Naim, JAFC 2001, 49, 5425-5432 | no (low priority) |
| thresholds | water / air | MFT "0.007 mg/kg in water and 0.0025 ng/L in air"; FFT "0.01 mg/kg in water of 0.01 ng/kg in air" (as extracted; the water values are µg/kg in the primaries, see §5.1) | [82, 85, 86] | — |

Statements on **disulfide formation**: [84] (ether, 6 °C), [86] (Fenton, difurfuryl disulfide the major
product), [91] (canned coffee, difurfuryl disulfide the major volatile product). On **melanoidin binding**:
[94], [95] (covalent, oxygen-dependent). On **protein binding**: [83] (thiol-disulfide interchange, egg
albumin). On **metals**: [86] (Fe(II) > Mn, Cu; Fenton), [91] ("Fenton-type"). On **oxygen**: [87] (loss with
and without air; MFT loss not oxidative), [94] ("without oxygen ... very slow"). On **pH**: [91] only. On
**thiol + aldehyde adducts (hemithioacetal, thiazolidine)**: none for MFT or FFT anywhere in the review.
On **methanethiol**: only as a methional breakdown product, [96], [97], no rate. **2-mercapto-3-butanone
and 3-mercapto-2-pentanone**: not discussed beyond the citation of [85].

### 4.2 2-Acetyl-2-thiazoline and the thiazine

| compound | condition | statement / number | primary | on disk? |
|---|---|---|---|---|
| 2-acetyl-2-thiazoline (2-AT) | reflux in tap water, 100 °C, 60 min | "the degradation was 96 % after 60 minutes" | [109] Hofmann & Schieberle, JAFC 1995, 43, 2946-2950 | **no** (`hofmann1995.pdf` is the ribose/cysteine AEDA paper, JAFC 43, 2187) |
| 2-AT | sunflower oil, 100 °C | "only 2.5 % 2-AT degradation" (time not restated; presumably the same 60 min) | [109] | no |
| 2-AT | phosphate buffer pH 5, autoclave, 145 °C | "99 %" degradation (time not restated) | [109] | no |
| 2-AT | general | "unstable during heat treatment in the presence of water. Fat-containing food systems can stabilize 2-AT"; formed from cysteine + methylglyoxal | [108] Meynier & Mottram, Food Chem 1995, 52, 361-366 (on disk `meynier1995.pdf`, no dossier); [110] Pripis-Nicolau et al., JAFC 2000, 48, 3761 (no) | — |
| 5-acetyl-2,3-dihydro-4H-1,4-thiazine | — | "first identified from a model ribose-cysteine reaction system and has not been identified in food systems"; threshold 0.05 ng/L air; no stability number | [103] Adams & De Kimpe, Chem Rev 2006, 106, 2299-2319 (review) | no |
| 2-AT threshold | water | "1 mg/L" as extracted (µg/L in the primaries) | [103] | — |

### 4.3 Furanones and Strecker aldehydes (brief)

| compound | condition | statement / number | primary | on disk? |
|---|---|---|---|---|
| furaneol (HDMF) | aqueous buffer | degradation pH-dependent, "optimum stability being at pH 4", first-order; sucrose and ethanol 0-20 % no effect | [53] Hirvi, Honkanen, Pyysalo, Lebensm.-Wiss. u.-Technol. 1980, 13, 324-325 | **no** |
| furaneol | pH 2.0-8.0, 23 °C | greatest stability at pH 3.5 | [54] Roscher, Schwab, Schreier, Z. Lebensm. Unters. Forsch. 1997, 204, 438-441 | **no** |
| furaneol | closed system, 160 °C, pH 2.2 / 5.1 / 7.1 | degradation "preferred at a lower pH"; ring opening then retro-aldol to acetaldehyde, hydroxyacetone, 1-hydroxy-2-butanone, acetoin, 2,3-butanedione | [50] Shu, Mookherjee, Ho, JAFC 1985, 33, 446-448 | no (`shu1988.pdf` is [56], DMHF + cysteine pH study) |
| furaneol + cysteine | 160 °C, pH 2.2 / 5.1 / 7.1 | more thiophenes at pH 2.2 than 5.1; pyrazines only above the pI of cysteine (pH 7.1) | [56] Shu & Ho, JAFC 1988, 36, 801-803 | yes (`shu1988_extraction.md`) |
| furaneol + cysteine / glutathione / Na2S / alanine | 130 °C | H2S availability "might be the limiting factor" | [49] Zheng, Brown, Ledig, Mussinan, Ho, JAFC 1997, 45, 894-897 | no (`zheng1994.pdf` is the H2S-release kinetics chapter) |
| furaneol, homofuraneol, norfuraneol thresholds | water | "0.04 mg/kg", "20 mg/kg", "23,000 mg/kg" as extracted (µg/kg in the primaries) | [34-38] | — |
| methional | heat, light | "heat labile and readily decomposes to methanethiol, which oxidizes to dimethyl disulfide"; light gives methanethiol and dimethyl sulfide; reacts with hydroxyl radical to give ethylene; "Data on the mechanism of methional degradation are still lacking"; no rate | [96] Di et al. JAFC 2003, 51, 5695; [97] Jung et al. J Food Sci 1998, 63, 408; [100-102] | no |
| 2-acetyl-1-pyrroline, 6-acetyl-tetrahydropyridine | popcorn, polyethylene bag, 1 week | decreased "75 and 69 %" | [107] Schieberle, JAFC 1995, 43, 2442-2448 | no |
| vanillin + amino acids | 55 / 65 / 75 °C | first-order loss, rate rises with temperature; Schiff base "not a major reaction" | [74] Chobpattana, Jeon, Smith, JAFC 2000, 48, 3885 | no (outside scope) |
| sotolon | — | **not mentioned anywhere in the review** | — | — |

## 5. Flags

1. **Unit-glyph loss in the text layer.** Every water threshold reads "mg/kg" or "mg/L" in the extraction
   (MFT 0.007, FFT 0.01, furaneol 0.04, 2-AT 1, 2-AP 0.1) where the primaries (and Schieberle & Hofmann
   1998 Table IV on disk) give µg/L. Treat every "mg" prefix in this text layer as suspect; none of these
   thresholds should be transcribed from the review in any case.
2. **Citation year/volume slips**: [94] is printed "2001;50:319-26" — JAFC volume 50 is 2002 (web release
   Dec 2001), the paper on disk as `hofmann2001.pdf` and dossiered as `hofmann2002`. [55] Van den Ouweland
   & Peer is printed "1975;2:501" (volume 23).
3. **Attribution merging.** The melanoidin numbers ("50 % after 20 minutes and almost 100 % after 30
   minutes") are attributed jointly to [94] and [95]; `hofmann2002_extraction.md` records a different
   shape from [94] (80 % bound within 30-90 min at 30 °C then a plateau). The 20/30-min figures are most
   likely from [95] (Hofmann, Czerny, Calligaris, Schieberle 2001, coffee beverage kept warm), which is
   not on disk. Do not attach them to hofmann2002.
4. **The review's MFT mechanism paragraph** (oligomerisation via C-2 protonation, C-5 attack) is van
   Seeventer's proposal; the review states it as fact. The repo already holds van Seeventer's 50 °C rows
   as HOLD-OUT with the oligomer channel at zero (`thiol_sink_candidates.md` §1).
5. **Nothing quantitative** on 2-mercapto-3-butanone, 3-mercapto-2-pentanone, methanethiol, disulfide
   yields, thiol-aldehyde adducts or sotolon. The review does not shorten the fetch list for those.
6. **Registry keys**: MFT `2_methyl_3_furanthiol`, FFT `2_furfurylthiol`, bis(2-methyl-3-furyl) disulfide
   `bis_2_methyl_3_furyl_disulfide`, methanethiol `methanethiol`, methional `methional`, dimethyl disulfide
   `dimethyl_disulfide`, HDMF `hdmf`, HEMF (homofuraneol) `hemf`, norfuraneol `norfuraneol`. No key for
   3-mercapto-2-pentanone, 2-mercapto-3-butanone, 2-acetyl-2-thiazoline, 5-acetyl-2,3-dihydro-1,4-
   thiazine, 2-furfuryl methyl disulfide, difurfuryl disulfide, or sotolon.

## 6. Fetch list (primary papers not on disk, citation as printed in the review's reference list)

Ranked by what they would give the sink re-fit.

1. **[85] Kerscher R, Grosch W. Quantification of 2-methyl-3-furanthiol, 2-furfurylthiol, 3-mercapto-2-
   pentanone, and 2-mercapto-3-pentanone in heated meat. J Agric Food Chem 1998;46:1954-8.** SIDA levels
   of the four sulfur odorants in real cooked meat (boiled, roasted, several temperatures): end-of-cook
   validation levels in a food matrix for the same four species the fit carries, from the Garching group
   whose model-pot numbers anchor the fit.
2. **[83] Mottram DS, Szauman-Szumski C, Dodson A. Interaction of thiol and disulfide flavor compounds
   with food components. J Agric Food Chem 1996;44:2349-51.** The only primary on thiol/disulfide loss to
   a protein at cooking temperature (100 °C, egg albumin) in the corpus's target range; would give a
   protein-sink magnitude at 100 °C to set against Anantharamkrishnan 2020b's ambient bracket.
3. **[109] Hofmann T, Schieberle P. Studies on the formation and stability of the roast-flavor compound
   2-acetyl-2-thiazoline. J Agric Food Chem 1995;43:2946-50.** 96 % loss in 60 min at 100 °C in water,
   99 % at 145 °C / pH 5, 2.5 % in oil: the time course and products of 2-AT hydrolysis, needed if the
   thiazoline/thiazine pair (both quantified in Schieberle & Hofmann 1998 dry vs aqueous) enters the
   sulfur lane.
4. **[95] Hofmann T, Czerny M, Calligaris S, Schieberle P. Model studies on the influence of coffee
   melanoidins on flavor volatiles of coffee beverages. J Agric Food Chem 2001;49:2382-6.** The
   time-resolved thiol loss in a warm coffee beverage (the "50 % at 20 min, ~100 % at 30 min" statement),
   the companion of hofmann2002; gives the melanoidin sink at drinking temperature with several thiols.
5. **[53] Hirvi T, Honkanen E, Pyysalo T. Stability of 2,5-dimethyl-4-hydroxy-3(2H)-furanone and
   2,5-dimethyl-4-methoxy-3(2H)-furanone in aqueous buffer solutions. Lebensm Wiss u Technol
   1980;13:324-5** and **[54] Roscher R, Schwab W, Schreier P. Stability of naturally occurring
   2,5-dimethyl-4-hydroxy-3[2H]-furanone derivatives. Z Lebensm Unters Forsch 1997;204:438-41.**
   First-order HDMF loss rates versus pH (optimum pH 3.5-4) at ambient: the pH shape of the DMHF sink,
   which the repo's k5b synthesis lacks below cooking temperature.
6. [50] Shu CK, Mookherjee BD, Ho CT. Volatile components of the thermal degradation of
   2,5-dimethyl-4-hydroxy-3(2H)-furanone. J Agric Food Chem 1985;33:446-8. Product spectrum of HDMF at
   160 °C by pH; qualitative.
7. [49] Zheng Y, Brown S, Ledig WO, Mussinan C, Ho CT. Formation of sulfur-containing flavor compounds
   from reactions of furaneol and cysteine, glutathione, hydrogen sulfide, and alanine/hydrogen sulfide.
   J Agric Food Chem 1997;45:894-7. HDMF + sulfur donors at 130 °C; the H2S-limitation claim.
8. [103] Adams A, De Kimpe N. Chemistry of 2-acetyl-1-pyrroline, 6-acetyl-1,2,3,4-tetrahydropyridine,
   2-acetyl-2-thiazoline, and 5-acetyl-2,3-dihydro-4H-thiazine: extraordinary Maillard flavor compounds.
   Chem Rev 2006;106:2299-319. A review; only as a pointer to 2-AT/thiazine primaries.
9. [90] Kumazawa K, Masuda H, Nishimura O, Hiraishi S. Change in flavor of coffee drink during heating.
   Nippon Shokuhin Kagaku Kogaku Kaishi 1998;45:108-13. Japanese; superseded by kumazawa2003 on disk.

**On disk but without a standalone dossier** (no download needed): `vanseeventer2001.pdf` (ingested via the
inventory), `charles-bernard2005.pdf`, `meynier1995.pdf`.
