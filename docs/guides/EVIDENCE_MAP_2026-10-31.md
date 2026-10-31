# Evidence map: within-study factor matrix and mechanism-discriminating experiments (2026-10-31)

READ-ONLY audit of `data/lit/extraction_dossiers/*.md` (475 files) and `data/lit/fit_rows/*.yml`, for
STATUS_AND_PLAN roadmap row 3 (identifiability audit) and sec. 5 "How a mechanism is identified" / "Mechanism triage".
Nothing was run or edited in the repo.

## How this was made, and how far to trust it

- **Method.** Every dossier was grepped. The candidate dossiers' Identity, Methods and contrast lines were read by four
  read-only search passes, one per factor group or mechanism group, then merged here. `file:line` means
  `data/lit/extraction_dossiers/<file>_extraction.md:<line>`.
- **Quote checks.** I re-checked ten load-bearing quotes against the files with `sed -n`; all ten match:
  hofmann2000:112-113, cerny2003:64, amranihemaimi1995:46, davidek2006b:22, gobert2009:24-26, nedvidek1992:241,
  vanseeventer2001:52, zhao2019:33, zyzak1995:22, cao2017:83. The other rows' quotes come from the search passes and were
  not re-opened; treat them as accurate to the dossier, not as checked against the PDF.
- **Keyed status** comes from the pot names in `data/lit/fit_rows/*.yml`:
  - **K** = the contrast's arms are keyed fit rows.
  - **k** = the paper is keyed, but the arms that make this contrast are not, or only one arm is.
  - **–** = not keyed.
- **Counting labs.** Labs are counted strictly:
  - The Ho network (Rutgers before 2000; Zhang/Cui/Ho at Jiangnan; Zhai/Fu at Zhejiang Gongshang) counts as ONE lab.
  - Hofmann, Schieberle, Granvogl and Grosch (DFA Garching / TUM) count as one.
  - Mottram, Parker, Elmore, Ames and Whitfield (Reading / IFR) count as one.
  - Davidek, Blank, Cerny-at-Nestlé, Poisson and Kerler count as one (Nestlé).
  - Cerny's Firmenich papers are listed separately but flagged as overlapping with Nestlé.
  - van Boekel, Martins, Brands and Knol count as one (Wageningen).
- **"Weak" (W)** means the contrast is confounded, with the reason given in a few words.
- **Response tags:**
  - DO = deoxyosones (1-DG / 3-DG / 1-DX / 3-DX / glucosone)
  - DC = other dicarbonyls (MGO, GO, diacetyl)
  - SA = Strecker aldehydes; SAc = Strecker acids
  - PZ = pyrazines
  - TH = thiols (MFT, FFT, ...); H2S
  - OA = organic acids (formic, acetic)
  - Cys = free cysteine; TTCA / ARP = thiazolidine / Amadori compounds
  - SL = sugar or amine loss; BR = browning; OT = other
- **Quantification:** SIDA = isotope dilution; cal = calibrated; IS = one internal standard; area = peak area.

---

# TABLE 1 — Within-study factor matrix (roadmap row 3)

The decision rule: a factor gets a FITTED coordinate for a response class only if at least 2 studies from at least 2
independent labs vary it inside one study, with the other conditions fixed, AND measure that response class.
Otherwise it gets a declared band. Factor-level totals hide this, because a factor varied by 40 labs for browning can
still have one lab for deoxyosones. The verdict is therefore given per response class (sec. 1.9).

## 1.1 Free amine identity / concentration (incl. peptides, protein-bound lysine)

| dossier | lab | levels | responses | reps / quant | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| martins2005b | Wageningen | Glc:Gly 1:1 / 2:1 / 1:2 | DO, OA, ARP | >= 2 / cal | k (pH arms keyed, ratio arms not) | "molar ratios of 2:1 and 1:2" :111 |
| xia2022 | Ho network (Jiangnan) | Gly / diGly / triGly + Glc | ARP, DO (1-DG, 3-DG), DC | ns / cal | – | "Glycine, Diglycine, and Triglycine Exhibit Different Reactivities" :16 |
| deng2021 | Ho network | MG-ARP +/- Met 0.1 M | ARP decay, DO, DC | 3 / cal | K (phosphate arms) | "+- DL-Met 0.1 mol/L" :1 |
| zhou2024 | Ho network | AX-ARP +/- 15N-Ala | DO, DC, PZ | 3 / cal | K (DO / DC) | "Ala-xylose Amadori compound +- 15N-alanine at 120 C" :1 |
| edwards1992 | Wedzicha, Leeds | Gly 0-1.0 M; Cys vs homocysteine vs cysteamine | 3-DG (DH) consumption, thiol | n=6 / cal | – | "initial rate = k1(1 + k2[glycine])[DH]" :61 |
| wedzicha1991 | Wedzicha, Leeds | Gly 0-1.0 M | 3-DG (DH), thiols | 5-15 / cal | – | "+- glycine 0-1.0 M, 50 mM acetate pH 5.5" :1 |
| weenen1998 | Weenen, Quest | Ala / cyclohexylamine / none | DO, DC | ns / uncalibrated OPD | – | "heated with alanine (0.6 M), cyclohexylamine, or no amine" :3 |
| potman1989 | Unilever | glycine-order series | SL, ARP | none / cal | – | "Series C -- order in glycine (Figs. 6-8)" :61 |
| nedvidek1992 | Ledl, Stuttgart | beta-Ala vs hydrolysed wheat protein | DO (quinoxalines) | none / area | – | W: temperature never stated; :1 |
| hofmann2000 | DFA | Phe / Ala / Leu + Glc | SA, SAc | ns / SIDA | K | "Generated from Three Different alpha-Amino Acids" :62 |
| granvogl2006 | DFA | Leu / Ile / Phe / Met + MGO | SA, Strecker amines | 3 / SIDA | – | "Leu / Ile / Phe / Met (not Val) + methylglyoxal" :3 |
| huang2017 | Becker, TUM Freising | Leu vs Ile (x Glc / maltose) | SA, SL | 3 / unclear | – | "glucose or maltose + leucine or isoleucine, 10:1 molar" :1 |
| martin2001 | Reading (Ames) | potato, one amino acid changed at a time | SA, PZ | 3 / IS | – | "with one amino acid changed at a time" :3 |
| kocadagli2021 | Reading (Parker) | 7 amino acids + whey hydrolysate dosed | SA, PZ, acrylamide | ns / SI | – | :1 |
| balagiannis2009 | Reading | 3 liver batches (free amino acids 3-4x) | SA | ns / cal | – | W: natural batches; :34 |
| solina2007 | Whitfield, FSA | +/- acid-HVP | SA, PZ | 3 / IS | – | "Strecker aldehydes rising 25-fold when free amino acids are added" :3 |
| lee2024 | INRAE | Glc vs Glc+Leu cake | SA, PZ | ns / on-line MS | – | "found ONLY in the glucose + leucine cake" :3 |
| adams2008 | De Kimpe, Ghent | 20 amino acids + DHA; Ala:MGO 1:1 / 2:1 / 1:3 | PZ, pyrrole | ns / cal (T2-3) | – | "Ala : MGO 1:1, 2:1, 1:3" :45 |
| vanlancker2012 | De Kimpe, Ghent | 9 di- + 3 tripeptides vs free amino acids x Glc / MGO / GO | PZ | none / area | – | :1 |
| amranihemaimi1995 | Nestlé (Cerny / Fay) | Ala vs Gly | PZ (labelled) | none / area | – | :1 |
| hwang1995 / 1995b | Ho network (Rutgers) | Glc + 15N-Gly + 9 second amino acids | PZ / pyridines, pyrroles | ns / IS | – | :1 |
| leahy1989 | Reineccius, Minnesota | Lys vs Asn x 3 sugars | PZ | 2-3 / cal | – | :1 |
| meynier1995 | Reading | Gly / Lys / Cys / Met + Rib (x 5 pH) | TH, PZ | 3 / IS | – | :1 |
| whitfield1988 | Reading / IFR | Gly / Lys / Cys + Rib +/- lecithin | TH, PZ | 3 / cal ng | – | :1 |
| cao2017 | Xie, BTBU | Gly 0 / 1 / 8 / 12 / 16 / 24 x Cys | Cys, TTCA, ARP, volatiles | 3 / mixed | k (blank and 1x keyed) | "at 0, 1, 8, 12, 16, or 24 times of the cysteine" :83 |
| zhang2025b | Xie, BTBU | Lys 0-16 x Cys (7 levels) | TH, TTCA, ARP, Cys | 3 / IS | k (blank only) | :35 |
| feng2022 | Ho network | RG-ARP; +GSH / +Cys / +Glu / +Gly | TH, PZ, ARP | ns / IS | k (+/-GSH keyed) | :176 |
| kang2026 | Ho network (Zhai) | TTCA + 13 amino acids + Cys | TH, PZ | 3 / IS | – | :130 |
| zhang2023b | Cao, BTBU | Met +/- 16 amino acids | DMDS, DMTS, MeSH | 3 / IS | – | :53 |
| zheng1994 | Ho network | Cys vs GSH, 16 matched conditions | H2S | none / cal | – | :167 |
| brands2001 | Wageningen | sugars +/- casein (protein-bound Lys) | OA, SL | >= 2 / cal | k (casein arms) | "with and without an amine" :27 |
| ajandouz2008 | Marseille | Glc alone vs + BSA / + casein | SL, BR | >= 2 / cal | – | :18 |
| tang2024 | Hettinga, Wageningen | soy vs pea protein + Glc | Lys, CML, CEL | 2 / LC-MS | – | W: whole-protein swap |
| goncuoglu2026 | Reading / Hacettepe | 0-70 % pea isolate in cornmeal | PZ, acrylamide | ns | – | W: starch and moisture co-vary |

Also checked (amine or peptide rows outside the hub and odorant responses):
- **Acrylamide:** claeys2005, devleeschouwer2009 (W).
- **Other:** quan2020 (GO / MGO / CML); hamzalioglu2018; zamora2010 / 2020; hidalgo1993 / 2016; kim1998; zhang1988 / 1989 / 1991 / 1991b; tai1997 / 1998; yu1995; chen2000 / 2002; wang2008; wang2012 (W); xu2013; hou2017; zhao2019 (Amadori-bound, W); zhou2023c (k); haleva1999; friedman1990; rizzi2004; cammerer1994; wedzicha1992 (W); hayashi1980 (W, ethanol); zhang2020; wang2025 (W); bornhorst2017 / 2017b (W, pH rises with Lys); hemmler2018; shu1999; poisson2009; lin2021 / 2022; coukos2021; whitfield1999 (W); us3904655 (W); tonsbeek1969 (W).

**Count: about 79 rows; about 29 strict independent labs.** The count that matters, by response:

- **Deoxyosones vs amine concentration.**
  - Wageningen (martins2005b, ratio), Leeds (edwards1992, wedzicha1991: glycine order of 3-DG consumption), Quest (weenen1998, uncalibrated), Unilever (potman1989, sugar loss only).
  - That is at least 2 labs, so a FITTED coordinate is allowed for the glycine order of the early steps (hexose).
  - Pentose: Ho network only (deng2021, zhou2024). Declared band.
- **Peptide vs free amine on the hub:** xia2022 only (one lab). Declared band.
  - cui2022 (glycylglycine, 10.1021/acs.jafc.2c03427) is in the private download queue, tier 2.
- **Peptide vs free amine on pyrazines:**
  - De Kimpe (vanlancker2012, peak areas) only on disk; scalone2015 and vanlancker2010 are in the queue.
  - That leaves one lab, so a declared band.
- **Protein-bound lysine vs free amine** (same sugar, same pot design):
  - Wageningen (brands2001) and Marseille (ajandouz2008) give sugar loss / acids only. Nothing on disk measures DO or odorants.
  - The sugar-loss coordinate can be fitted. DO and odorants need a declared band (as sec. 4 item 8 already says).
- **Amino-acid identity in Strecker:**
  - DFA (hofmann2000, granvogl2006), Becker (huang2017), Reading (martin2001, kocadagli2021, balagiannis2009 ratio), Reineccius (chan1994b, all in one pot).
  - That is at least 3 labs, so a FITTED identity factor is allowed for Leu / Ile / Phe / Met / Ala.
  - Caveat: only hofmann2000 and granvogl2006 fix the dicarbonyl donor, and they are the same lab.
- **Amino-acid identity on pyrazines:** Ghent, Nestlé (labelled), Ho network, Reineccius, Reading. Fitted for the
  identity ranking. Peak areas dominate, so it is directional only.

## 1.2 Phosphate / buffer

| dossier | lab | levels | responses | reps / quant | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| davidek2002 | Nestlé | water vs 0.1 M phosphate, pH-stat 5 / 6 / 7 / 8 | ARP (DFG) decay, OA, SL | ns / cal | **K** (8 pots) | "in water or in 0.1 M phosphate" :3 |
| zyzak1995 | Baynes, USC | phosphate 10 / 50 / 100 / 200 mM; phosphate vs HEPES | Amadori (fFL) decay + products | 3 / cal | – | "at 10 / 50 / 100 / 200 mM phosphate (triplicates)" :22 |
| deng2021 | Ho network | water vs 0.1 M phosphate, 100 C | ARP decay, DO | 3 / cal | k (phosphate arm only) | W: unbuffered drift; "~1.8x faster in phosphate" :132 |
| potman1989 | Unilever | NaH2PO4 0-0.5 M (6 levels) | SL, ARP | none / cal | – | :59 |
| bell1997 | Bell, Auburn | phosphate or citrate at 0.02 / 0.05 / 0.2 / 0.5 M, pH held | Gly loss, BR | 95 % CL / cal | – | :35 |
| rizzi2004 | Rizzi, P&G | phosphate vs bis-tris / citrate / succinate / acetate | BR | range / A420 | – | :3 |
| larralde1995 | Miller, UCSD | phosphate vs bicarbonate / acetate (no amine) | ribose loss | none / NMR | – | :67 |
| blank1997 | Nestlé | 0.2 M phosphate / malonate / water, pH 6 | HDMF, EHMF (furanones) | 2 / SIDA | – | water arm W; :606 |
| mottram2002 | Reading | pyrophosphate vs phthalate vs unbuffered, pH 5.6 | TH, norfuraneol, furfural | 3 / IS | – | unbuffered arm W (5.6 -> 3.9); :41 |
| parker2014 | Reading | pyrophosphate 5.5 vs water | TH (MFT), furfural | 2 / area | – | W: unbuffered pH 3.0-4.7 |
| adams2008 | Ghent | phosphate vs bicarbonate, pH 7 | PZ, pyrrole | ns / cal | – | :149 |
| devleeschouwer2006 | Leuven | citrate vs phosphate, pH 4 / 6 / 8 | acrylamide, SL | 2 / cal | – | :23 |
| hidalgo2010 | Seville | citrate vs phosphate, pH 6 | acrylamide | >= 2 / SIDA | – | :227 |
| fife1991 | Fife, USC | acetate / phosphate / cacodylate buffer series | thiazolidine ring-opening / hydrolysis k | ns / UV | – | :29 |
| lievonen2002 | Roos | water vs 0.05 vs 0.5 M citrate (glassy solids) | BR | 3 | – | :146 |
| martins2005b / brands2002c | Wageningen | pH-stat vs free drift | DO, OA | >= 2 / cal | K (pH-stat in hold-out) | W: pH control, not buffer |
| us10076126 / yaghmur2005 | Nestec / Garti | phosphate vs glycerol / propylene glycol / microemulsion | TH; furfural loss | none / 3 | k / – | W: solvent swap |

Also W, because pH co-varies: hidalgo1993, ge1997, ajandouz2008, weenen2001.

**Count: about 23 rows; about 12 strict labs with a clean contrast.** By response:

- **Amadori / ARP decay:** Nestlé (davidek2002, K), Baynes (zyzak1995, 37 C), Ho network (deng2021, W).
  - That is 2 clean labs, so a FITTED phosphate term on ARP decay is allowed. The temperature gap (37 vs 90-120 C) should be declared.
- **Deoxyosone split vs phosphate: NO study.** Declared band.
  - This is exactly the Cao-vs-Zhai / B7b question.
  - thornalley1999 and guan2012 are in the queue tier 2.
- **Thiols vs phosphate:** Reading only (mottram2002 clean, parker2014 W). Declared band.
- **Cysteine / thiazolidine hydrolysis vs buffer:** fife1991 only (non-cysteine rings). Declared band in form and magnitude;
  pesek1975 is in the queue.
- **Sugar loss, browning, acrylamide:** many labs. Fitted coordinate allowed, but this is not a hub response.

## 1.3 O2 (air vs inert)

| dossier | lab | levels | responses | reps / quant | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| hofmann2000 | DFA | argon / air / air + Cu, Phe + MGO / GO / 3-DG / glucosone | SA, SAc | ns / SIDA | **K** | "PAA air / argon \| 4.34" :113 (checked) |
| hofmann2000b | DFA | argon / air / air + Cu, ARP-Phe vs Glc+Phe | SA, ARP decay, DO (OPD) | ns / SIDA | **K** | "9x more phenylacetaldehyde under air than under argon" :25 |
| davidek2006b | Nestlé | argon / air / 18O2, Glc+Gly pH 8, 90 C | OA, C4-aldonic acids | +- SD / SIDA | **K** | "argon 30.5 +- 0.4, air 32.3 +- 3.3, 18O2 33.5 +- 0.5" :22 (checked) |
| gobert2009 | Glomb, Halle | aerated vs deaerated, Glc+Lys, 50 C, OPD in situ | DO (glucosone, 1-DG, 3-DG), DC | ns / cal | – | "glucosone \| 44.5 \| 5.5" :26 (checked) |
| smuda2010 | Glomb | aerated vs deaerated, 1-DG + Boc-Lys | OA, amides | ns / cal | – | W: DTPA only in deaerated pots |
| zyzak1995 | Baynes | air vs N2 + chelators, 37 C | fFL decay | 3 / cal | – | W: chelators co-varied |
| zamora2015, hidalgo2013, hidalgo2016, zamora2013, zamora2008 | Seville | N2 vs air; O2 0-100 % | SA (from lipid carbonyl or keto-acid donors), amine | >= 3 / cal or IS | – | "Air multiplies PAC by 4-8" zamora2015:124 |
| charlesbernard2005 | Nestlé | air vs N2, coffee 25 C | TH (10 thiols, no MFT) | 2 / area | – | :216 |
| chin1994 | Lindsay, Wisconsin | aerobic vs anaerobic x Cu / Fe | MeSH -> DMDS / DMTS | 2 / cal | – | :27 |
| ehrenberg1989 | Stockholm | air vs O2; anaerobic arm | Cys oxidation | parallel / Ellman | – | :25 |
| vanseeventer2001 | Quest | air vs argon, stored flavouring 50 C | TH | – | – | W: "data not presented" :272 |

Also W or out of scope: wondrak2002, arnold1971, friedman1990, kreitman2016, xu2010 (pressure, not O2),
koelsch1991 (lipid).

**Count: about 22 rows; 7 strict labs with clean contrasts.** By response:

- **O2 on the hub dicarbonyls:**
  - Glomb (gobert2009, glucosone up 8x in air, 50 C) and DFA (hofmann2000b, ARP + O2 + Cu -> glucosone).
  - That is 2 labs, but both have metal / chelator confounds and neither is at 100-140 C in a sugar + amino-acid pot.
    PROVISIONAL fitted coordinate at best; a declared band is the safer call.
- **O2 on acids:** Nestlé only (davidek2006b, K). Declared band (single lab).
  - voigt2010 (10.1021/jf100140h) is queue tier 2 and would be the second lab.
- **O2 on Strecker aldehyde from sugar-derived dicarbonyls:** DFA only. Single-lab, so PROVISIONAL.
  - Seville covers lipid / keto-acid donors, a different donor class.
- **O2 on Strecker ACID:** DFA only. Declared band.
- **O2 on MFT / FFT loss at process temperature: NO clean study.**
  - vanseeventer2001 shows no data; charlesbernard2005 has no MFT, at 25 C. Declared band.

## 1.4 pH

Core rows only (sugar-amine, sulfur, Strecker). "pH-stat" and "buffered" are clean; "initial pH, unbuffered" is W.

| dossier | lab | levels | responses | TR | reps / quant | keyed | quote (file:line) |
|---|---|---|---|---|---|---|---|
| martins2005b | Wageningen | initial 4.8 / 5.5 / 6.0 / 6.8 / 7.5 + pH-stat | DO, DC, OA, SL, BR | Y | >= 2 / cal | **K** | :62 |
| martins2003b (+2003, 2003c) | Wageningen | 5.5 / 6.8 x 100 / 120 C, fed DFG | DO, SL, OA | Y | >= 2 / cal | **K** | :71 |
| davidek2002 | Nestlé | pH-stat 5 / 6 / 7 / 8 | ARP, OA | Y | ns / cal | **K** | :3 |
| davidek2005 | Nestlé | 4 / 6 / 8 (0.2 M phosphate) | OA, SL | Y | ns / cal | – | :3 |
| davidek2006a | Nestlé | 4 / 6 / 8; fed 1-DG pH 5-8 | OA, DO (fed) | Y | ns / cal | **K** (Glc / 3-DG arms) | :3 |
| davidek2008 | Nestlé | 6.0 / 8.0 | OA, xylose | EP | ns / SIDA | **K** | :3 |
| nedvidek1992 | Ledl | 5 / 6.5 / 7 | DO, DC (OPD in situ) | EP | none / area | – | :43 |
| zhai2021 | Ho network | 5.5 / 6 / 7 / 8 x 3 T, fed TTCA +/- Xyl | TTCA, Cys, DO, BR | Y | 3 / cal | **K** | :2 |
| tang2019 | Ho network | initial 5.4 / 7.4 / 9.4 (GSH + Xyl) | DO, ARP, SL | Y | 3 / cal | – | :3 |
| deng2021 | Ho network | initial 5.5 / 7.5 (water) | DO / DC, ARP | Y | 3 / cal | k | W: unbuffered |
| zhou2023 | Ho network | 6 / 7 / 8 x +/-Cys, ARP | TH, PZ | Y | 3 / IS | **K** | :97 |
| wedzicha1991 | Leeds | 3.85-6.00 scan | 3-DG + thiol | Y | / cal | – | :1 |
| hofmann2000 | DFA | 3 / 5 / 7 / 9 | SA, SAc | EP | ns / SIDA (figure) | **K** | :46 |
| weenen2001 | Quest | 7.0 / 3.2 (fed ARP) | SA | EP | ns / cal | – | W: buffer changes with pH; :71 |
| chan1994 / chan1994b | Reineccius | 6 / 7 / 8 x 5 T | SA, PZ, methional | Y | 2 / IS | – | W: drift not measured |
| leahy1989a | Reineccius | 5 / 7 / 9 x 3 T (Lys + Glc) | PZ | Y | 2-3 / cal | – | :1 |
| hofmann1998 (reconciliation T4 / T8 / T9) | DFA | 3 / 5 / 7 | TH (MFT, FFT), mercapto-propanone | EP | ns / SIDA | **K** | thiol_submodel pots `_ph3/_ph5/_ph7` |
| schieberle1996 | DFA | 3 / 5 / 7 | TH (OAVs) | EP | ns / SIDA | – | :35 |
| hofmann1998b / schieberle2005 | DFA | 3 / 5 / 7 / 9 | 2-AP, ATHP | EP | ns / SIDA | – | :37 |
| blank2003 | Nestlé | 6 / 7 / 8 | 2-AP, ATHP | Y | >= 2 / SIDA | – | :1 |
| cerny2007 | Firmenich | 5 levels 4.0-7.0 (13C5-xylose ternary) | TH, furfural, mercaptoketones, label shares | EP | ns / area | – | :2 |
| meynier1995 | Reading | 4.5-6.5 (5 levels, pyrophosphate) | TH, PZ, Lys / His loss | EP | 3 / IS | – | :1 |
| madruga1995 / mottram2002 / ames2001 | Reading | 4.0-5.6 / 4.2 vs 5.6 / 5.5-7.5 | TH, PZ | EP | 3 / IS | – | :3 / :209 / :15 |
| whitfield1999 vs whitfield2001 | Reading | 4.5 vs 6.5 (two papers, same design) | TH from fed norfuraneol | EP | 2 / cal | – | whitfield2001:27 |
| kumazawa2003 | Ogawa | 3.0-6.4 (6 levels, FFT survival, 121 C) | TH loss | EP | none / IS | **K** (thiol_submodel) | :31 |
| us9005689 | Nestec | 3 / 5 / 7 / 9 | H2S, cysteamine, DC | Y | none | – | :32 |
| zheng1994 | Ho network | 3 / 5 / 7 / 9 x 4 T | H2S | Y | none / cal | – | :253 |
| mikami2017 / nakamura2020 | Murata, Ochanomizu | 5.5-8.0 / 6.5 vs 8.0 | norfuraneol formation and decay | Y | 3 / cal | k (pH 8 fed-NF only) | :32 / :141 |
| wang2026, modi2008, xu2007, yu2010, kang2026 | Ningxia / CFTRI / CAU / Hubei / Ho network | 3-5 levels | TH, PZ, S-volatiles | EP | IS or area | – | wang2026:290 |
| tai1998, chen2000, chen2002, shu1985, shu1988 | Ho network | 2-3 initial-pH levels | TH, S-volatiles | EP | IS | – | W: unbuffered |

Also: haleva1999 (furaneol), monforte2021 (SA, o-quinone route), adams2008 (W), Seville x 6 (lipid-donor SA),
fife1991 / kallen1971 / lienhard1966 / kim1986 (thiol-carbonyl physical chemistry).
Sugar-only: larralde1995, garrett1969. Browning: potman1989, cammerer1994, ajandouz2008, friedman1990 (W), cha2019 (W).

**Count: about 64 core rows; about 30 strict labs.** By response:

- **Hexose deoxyosones vs pH:** Wageningen (calibrated, pH-stat), Nestlé (fed 1-DG; acids by label), Ledl (trap,
  semi-quantitative). FITTED.
- **Pentose deoxyosones vs pH:** Ho network only (zhai2021 fed TTCA; tang2019 GSH-Xyl; deng2021 is Met-Glc).
  - That is ONE lab, so a declared band, or PROVISIONAL with one cross-lab check owed.
  - Nestlé's davidek2008 gives pentose ACIDS at pH 6 / 8, not deoxyosones.
- **Acids vs pH:** Wageningen + Nestlé. FITTED.
- **MFT / FFT vs pH:** DFA (SIDA), Reading, Firmenich, Ho network, Ogawa (FFT loss). FITTED in direction.
  - The calibrated subset is DFA only (hofmann1998 / schieberle1996).
- **Strecker aldehyde vs pH, sugar-derived donor:** DFA (SIDA) and Quest (two points, buffer co-varies), plus Reineccius (W).
  - Effectively ONE clean lab. PROVISIONAL.
- **Pyrazines vs pH:** Reineccius (leahy1989a) and the Ho network (zhou2023); Reading at end point only.
  - That is 2 labs, so FITTED in direction.
- **H2S vs pH:** Ho network (zheng1994) and Nestec (us9005689, no replicates). PROVISIONAL.

## 1.5 Temperature (isothermal series)

| dossier | lab | levels (C) | responses | reps / quant | keyed |
|---|---|---|---|---|---|
| martins2005 | Wageningen | 80 / 90 / 100 / 110 / 120 | DO, DC, OA, SL, BR | ns / cal | **K** (acids at 80 / 100 / 120) |
| martins2003b | Wageningen | 100 / 120 (fed DFG) | DO, SL, OA | >= 2 / cal | **K** |
| brands2002c, knol2005, knol2010, nguyen2016 | Wageningen | 90-130 / 120-200 / 120-200 / 120-130 | SL, OA, acrylamide, CML | 2-3 / cal or SIDA | nguyen2016 **K**; others – |
| desclaux2006 | Reading (Mottram) | 80 / 100 / 120 | DO, DC (uncalibrated OPD) | ns | k (100 C only) |
| balagiannis2009 | Reading | 120 / 130 / 140 | SA (2- / 3-MB), SL | >= 2 / cal | – (W: batch per temperature) |
| yiltirak2026 | Reading (Elmore) | 100 / 110 / 120 / 130 (time 240 -> 30 min) | TH (MFT, FFT) | / SIDA | **K** (W: time co-varies with temperature) |
| hamzalioglu2026, berk2021, goncuoglu2016, kocadagli2016 x 2, gursul2020 | Gökmen, Hacettepe | 3-4 levels, 4-220 C | DO (3-DG, 1-DG, glucosone), DC, SL | / cal | – |
| xia2022, wei2022, zhang2025, zhai2021, tang2019, deng2021 | Ho network | 2-6 levels, 50-140 C | DO (pentose and hexose), ARP | 3 / cal | zhai2021 **K**; zhang2025, deng2021 k |
| zhai2019, zhai2023foodchem, zhai2023jafc, kang2026 | Ho network | 70-130 / 100-140 | TTCA, Cys-ARP, TH ladder | 3 | first three **K** |
| zhou2024 / zhou2025b / deng2022 | Ho network | 100-120 / 70-90 / 100-130 | PZ, SA (methional) | 3 / IS | k (pyrazine rows not keyed) |
| li2011 | Chen, Jiangnan | 100 / 110 / 120 / 130 | Cys, Glc loss | ns / cal | – |
| yu2020 | Yao, Jiangnan | 100 / 110 / 120 | 3-DG, MGO, SL | 3 / uncalibrated | – |
| zhang2020 | lab not printed | 90 / 95 / 100 / 105 / 110 | DO (glucosone, 1-DG, 3-DG, 3,4-DDG), DC | 3 / cal | – |
| zhang2026 | BTBU | 130 / 140 / 150 (GSH + Glc) | SL, TH | 3 / cal | **K** |
| lee2022 / lee2024 | INRAE | oven 140 / 170 / 200 | DO, DC, SA, PZ | 1-3 / cal or IS | – (W: non-isothermal) |
| leahy1989 / 1989a, chan1994 / 1994b, schirle1992 | Reineccius | 75-150 | PZ, SA, furfural, S-volatiles | 1-2 / cal or IS | – |
| cremer2000 | Eichner, Münster | 80 / 90 / 100 / 110 (a_w 0.52) | SA (4), ARP | ns / standard addition | – |
| huang2017 | Becker, TUM | 90-130 (5 levels) | SA | 3 / IS | – |
| pan2025, Cai2024 | Wu, China Agricultural University | 100 / 120 / 140 | methional, MeSH; 2-AP | 3 / cal | – |
| zheng1994 | Ho network | 80 / 90 / 100 / 110 | H2S | none / cal | – |
| hidalgo2013, zamora2013, zamora2015 | Seville | 80-190 | SA (lipid or keto-acid donors) | / IS or cal | – |
| lau2003, pandit2006, bornhorst2017b | Tang, WSU | 80-131 | norfuraneol | 2-3 / cal | – |
| wang2022, wang2026, meng2017, sun2019, Zhou2025 | BTBU / Ningxia / Iwate / Jiangnan + Fisk / Kunming | 3-5 levels | TH | IS | – (wang2026, meng2017 W) |

**Count: about 75 core rows; about 45 labs.** By response:

- **Deoxyosones vs temperature:** Wageningen, Gökmen, Ho network, Reading (uncalibrated), INRAE, zhang2020's lab. FITTED.
  - Pentose deoxyosones specifically: Ho network (calibrated) plus Reading desclaux2006 (uncalibrated OPD, a 50:1 fitted
    split that contradicts the standing levels). PROVISIONAL.
- **Strecker aldehydes vs temperature:** Reineccius, Eichner, Becker, Reading, CAU, INRAE. FITTED.
  - All are pseudo-zero-order or yield rates in mixed pots.
- **Pyrazines vs temperature:** Reineccius, Ho network, NUS (yu2018), INRAE. FITTED.
- **Thiols vs temperature:** Ho network, Reading, BTBU. FITTED in direction.
  - Calibrated MFT ladders exist only for yiltirak2026 (time co-varies) and the zhai2023foodchem / kang2026 ladder
    (one lab, IS). PROVISIONAL for barriers.
- **H2S vs temperature:** zheng1994 only (Ho network). Declared band.

## 1.6 Water activity / moisture

| dossier | lab | levels | responses | reps / quant | keyed |
|---|---|---|---|---|---|
| devleeschouwer2007 / 2008 | Leuven | 6 a_w (0.34-0.92) x 5 T; 4 initial a_w | SL, acrylamide, BR | ns / cal | – |
| warmbier1976 | Labuza, Minnesota | a_w 0.14-0.87 (glycerol) x 3 T | Lys, Glc, BR | ns / cal | – |
| pereyragonzales2010 | Malec, Buenos Aires | 6 a_w 0.33-0.98 x 3 T | available Lys | 2 / UV | – |
| bell1995 | Bell | 6 a_w 0.33-0.85 (PVP) | BR | ns / UV | – |
| lievonen2002, miao2004 | Roos | a_w / water content x 6-7 T | BR | 2-3 / UV | – |
| leahy1989a | Reineccius | a_w 0.32 / 0.58 / 0.75 / 0.85 (nonfat dry milk, 95 C) | PZ | 2-3 / cal | – |
| modi2008 | Modi / Taylor | a_w 0.60 / 0.80 / 0.90 | TH, S-odorants | 4 batches / area or IS | – |
| hartman1984b | Ho network | water / propylene-glycol mixes, a_w 0-0.75 | TH (MFT) | 1 / IS | – (W: thiamine, PG solvent) |
| kang2026_SI | Ho network (Zhai) | storage a_w 0.113 / 0.432 / 0.843 | TTCA / ARP retention | 3 / cal | – |
| moisio2015, wang2026c, goncuoglu2026, solina2007 | Helsinki / BTBU / Reading / FSA | extrusion moisture 13-70 % | volatiles, SA, PZ | 3 / area or IS | – (last two W) |
| schieberle1998, hofmann1998c, hofmann1998b, engel2004 | DFA | dry 180 C / 6 min vs aqueous 145 C / 20 min | TH, 2-AP (SIDA) | ns / SIDA | – (W: temperature and time differ) |
| tang2019, zhai2023c | Ho network | atmospheric vs vacuum dehydration | DO, ARP, TTCA | 3 / cal | – (W: not an a_w ladder) |
| us10076126, yaghmur2005 | Nestec / Garti | water vs glycerol / propylene-glycol media | TH; furfural loss | none / 3 | k / – |

Correction to the brief's hints: cremer2000 runs at a single a_w (0.52; cremer2000:34), so it is not an a_w series.

**Count: 26 rows; 14 clean same-recipe series; about 12 labs.** By response:

- **Sugar / lysine loss and browning vs a_w:** Leuven, Labuza, Malec, Bell, Roos. FITTED.
- **Deoxyosones on an a_w ladder: NO study.** Declared band.
- **Pyrazines:** Reineccius only. Declared band.
- **Thiols:** Modi / Taylor (peak area) plus Ho-network W. Declared band.
- The DFA dry-vs-aqueous contrasts (MFT x ~6 from fed C2 + C3, hofmann1998_reconciliation:286; 2-AP x 53) are confounded
  with temperature and time. They are directional only.

## 1.7 Thermal history (ramp vs isothermal)

Only five studies vary the history at a fixed recipe:
- **bornhorst2017b** (WSU): microwave vs preheated water vs non-preheated, at equal F90; browning and norfuraneol. The
  only clean one. Quote: "F90 at the cold spot: MAPS / hot water preheated / not preheated \| 10.9 / 11.0 / 10.9" (:437).
- **lee2022 / lee2024** (INRAE): oven fan 25 vs 50 Hz with logged core temperature; DO, SA, PZ.
- **Zhang2024** (BTBU Cao): one-stage vs two-stage. W: an arm also adds Cys + GSH.
- **zhang2023b** (BTBU Cao): sequential addition. W.
- **us9005689** (Nestec): multi-step vs single-step. W: the steps differ in pH, not temperature.

About 18 more labs MEASURE or correct a come-up time without varying it: Leuven's thermocouple-integrated fits,
Wageningen (heat-up included in the fit), Reading (desclaux2006, balagiannis2009, martin2001's frying profile),
Reineccius, Gökmen, Roos, Unilever and others. The DFA autoclave ramp (20 -> 145 C in 20 min: hofmann1995:123,
hofmann1997, engel2004) is the single history behind most keyed thiol rows (hofmann1998). It is not logged and not varied.

**Verdict:** a declared band for every odorant response. No keyed contrast exists. The one clean study (WSU) measures
browning and norfuraneol.

## 1.8 Cysteine concentration

| dossier | lab | levels | responses | reps / quant | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| zhai2019 | Ho network | TTCA 5 mM +/- Cys 5 mM, 120 C | DO (1-DX, 3-DX), GO | 3 / cal | **K** | :54 |
| zhai2023jafc | Ho network | TTCA + Cys 0-10 mol/mol (browning); 0 vs 1:1 (other responses) x 100 / 120 / 140 C | TH, PZ, DO, TTCA / ARP | 3 / cal | **K** (0 vs 1:1) | :108 |
| zhou2022 | Ho network | AX-ARP + Cys 0 / 6 / 14 / 20 mM | PZ, diacetyl, Cys, DO | 3 / cal | **K** (0 vs 20) | :253 |
| liu2022 | Ho network | GR-ARP + Cys 0 / 10 / 20 / 40 mM | 2-acetylthiazole, MGO, GO, PZ, Cys | 3 / cal | – | :3 |
| wei2022 | Ho network | GX-ARP +/- Cys | TTCA, DO, Cys | 3 / cal | – | :26 |
| zhou2023 | Ho network | ARP +/- Cys x pH 6 / 7 / 8 | PZ, TH | 3 / IS | **K** | :97 |
| feng2022, kang2026 | Ho network | ARP or TTCA +/- Cys | TH, PZ | ns or 3 / IS | k / – | |
| wedzicha1991 | Leeds | Cys 10-50 mM with 3-DG | Cys binding, 3-DG | / cal | – | "independent of [cysteine] over 10-50 mM" :29 |
| Zhang2024 (= zhang2024b) | Cao, BTBU | Cys 0.5-75 mg/mL (10 levels); Cys / GSH / cystine | TH (MFT, MMFT, MFT-MFT), MeSH, DMDS, DMTS | 3 / external curves (figure only) | – | :62 |
| zhang2023 / zhang2023b | Cao, BTBU | Met / Xyl +/- Cys | DMDS, DMTS | 3 / IS | – | |
| hou2017 | Xie, BTBU | Gly-Amadori +/- Cys | TH | ns / area | – | :13 |
| vanseeventer2001 | Quest | flavouring +/- Cys; Cys series (MFT minimum near 50 mM) | TH (MFT, FFT), furanones | 2 / IS | – | :176 |
| thomas2014 | Berdagué, INRA Theix | Cys 0 / 100 / 1000 / 3000 / 10000 mg/kg with Fru or Xyl | TH (MFT + disulfides) | 3 / area | – | :30 |
| xu2008 | Gao, China Agricultural University | Cys 16.7 vs 50 mM at Rib 50 mM | TH, norfuraneol | 2 / IS | – | :35 |
| madruga1997 | Madruga, UFPB | beef +/- Cys (vs IMP, thiamine) | S-furans | 3 / area | – | :1 |
| poisson2009 | Nestlé | recombinate +/- Cys spike | FFT, MeSH, furfural | 2 / SIDA | – | :95 |
| cerny2007b | Firmenich | ternary vs no-Cys control | TH label shares | ns | – | :43 |
| haleva1999 | Naim, Hebrew University | Cys 0 / 1 / 5 / 10 mM | furaneol, HMF, BR | 3 / cal | – | :36 |
| brehm2020 | DFA (Hofmann, TUM) | thiamine vs thiamine / Cys | thioethers | 2 / SIDA | – | :1 |
| arnold1971 | Nebraska | thiamine +/- Cys | H2S | 2 | – | :57 |

Also W (multifactor or response-surface designs, or co-varying sugar): raza2020, gao2020b, ma2026, martins2010, meng2017,
jhoo2002. Acrylamide-only: devleeschouwer2009, claeys2005.

**Count: about 31 rows; about 15 strict labs; the Ho network supplies 8 and BTBU 4.** By response:

- **Cys on deoxyosones (the 3-DX trap):** Ho network (zhai2019, zhou2022, wei2022, zhai2023jafc; all K except wei2022)
  and Leeds (wedzicha1991 / edwards1992, hexose 3-DG, 55 C).
  - That is 2 labs, so a FITTED trap coordinate is allowed for the hexose / 3-DG reference.
  - The pH law of the trap is single-lab (zhai2021).
- **Cys dose on MFT, calibrated:** Zhang2024 (figure-only), vanseeventer2001 (IS). Others are peak area.
  - PROVISIONAL: 2 labs, semi-quantitative.
- **Cys dose on pyrazines and thiazoles:** Ho network only (zhou2022, liu2022). Declared band.
- **Cys on free-Cys or TTCA kinetics:** a single lab family. Declared band.

## 1.9 Summary: which hub and odorant coordinates can be fitted

| factor | studies (strict labs) | deoxyosone hub | Strecker / pyrazine | thiols / H2S | decision for step 5 (hub) |
|---|---|---|---|---|---|
| amine conc. (glycine order) | ~79 (~29) | Wageningen + Leeds (+ Quest) | Strecker identity: DFA + Becker + Reading | Reading, BTBU (Gly / Lys dose on thiols) | **fit** the hexose order; pentose = band |
| amine identity (peptide / protein-bound) | subset | xia2022 only | Ghent only (peak area) | – | **band** |
| phosphate / buffer | ~23 (~12) | ARP decay: Nestlé (K) + Baynes; split: none | Ghent (pyrazine vs bicarbonate) | Reading only | **fit** the ARP-decay term; split = **band** |
| O2 | ~22 (7) | Glomb + DFA (confounded, 50-100 C) | DFA only (sugar donors) | none clean for MFT | **band** (provisional coordinate on glucosone only) |
| pH | ~64 (~30) | hexose: Wageningen + Nestlé + Ledl; pentose: Ho network only | DFA (+ Quest W); pyrazines Reineccius + Ho network | DFA + Reading + Firmenich + Ho network | **fit** hexose; pentose provisional with one cross-lab check |
| temperature | ~75 (~45) | 6 labs; pentose: Ho network + Reading (uncalibrated) | 6 labs | Ho network + Reading + BTBU | **fit** |
| water activity | 26 (~12 clean) | none | Reineccius only | Modi only | **band** |
| thermal history | 5 vary it (1 clean) | INRAE (oven, W) | INRAE (W) | BTBU (W) | **band**; log the come-up |
| cysteine conc. | ~31 (~15) | Ho network + Leeds | Ho network only | BTBU + Quest (semi-quantitative) | **fit** the trap existence; pH law = band |

---

# TABLE 2 — Mechanism-discriminating experiments on disk (kinds 1-5)

K1 = isotope labelling of a reactant (not SIDA standards); K2 = fed intermediate; K3 = trapping; K4 = directly measured
intermediate or lifetime; K5 = mass balance.

**Verdict labels:**
- **CAN** = literature can discriminate: at least one K1-K4 experiment from at least 2 labs on the SAME rival question.
- **SINGLE** = single-lab.
- **NONE** = no discriminating experiment; needs designed data.

## 2.1 Deoxyosone hub: split and lifetime

| kind | dossier | lab | conditions | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| K1 | davidek2006a / davidek2005 | Nestlé | Glc + Gly, [1- ... 6-13C]glucose, pH 4-8, 90 / 120 C | Acetic acid ~70 % from C1 / C2, so it comes from 1-DG beta-cleavage (2,3-enolisation) | K (acids) | "pairs C1/C2 ~70 %, C3/C4 ~8-10 %, C5/C6 ~20 %" davidek2006a:30 |
| K1 | davidek2008 | Nestlé | Xyl + Gly, 13C xylose, pH 6 / 8, 120 C | Pentose acetic acid 77-89 % from C1 / C2, via 1-deoxypentosone | K | davidek2008:26 |
| K2 / K4 | davidek2006a | Nestlé | synthetic 1-DG, no amine, pH 5-8, 120 C | 1-DG gone in < 5 min at pH >= 6; 0.8 mol/mol acid with no amine (amine not REQUIRED; acceleration untested) | K (no-amine arm not keyed) | :117 |
| K2 | smuda2010 | Glomb | 1-DG + Boc-Lys, 37 C | Amine-induced cleavage < 1 % of the acetyl flux; no amine-free arm | – | :34 |
| K2 / K4 | martins2003b | Wageningen | fed DFG, pH 5.5 / 6.8, 100 / 120 C | Time-resolved 1-DG and 3-DG; higher pH favours 1-DG; 1-DG decays faster | K | "An increase in pH seems to favour the 1-DG formation." :445 |
| K5 | martins2003 | Wageningen | DFG, pH 6.8 | Molar balance dips to 72-81 % early: a transient unmeasured pool | k | :644 |
| K4 | martins2005 | Wageningen | Glc + Gly, 80-120 C | Multiresponse fit: 1-DG lives ~1 min, 3-DG slower; both quasi-steady | K | :35 |
| K4 | hofmann1999 | DFA | Xyl or Glc + Ala, phosphate pH 7, reflux | Pentosones peak at 10 min, then 3-DX : 1-DX ~3 | keyed table, engine_species_missing | :37-38 |
| K3 | nedvidek1992 | Ledl | OPD in situ, pH 5-7 | 3-DG falls 39x and 1-DG 8x as pH rises (semi-quantitative; temperature unstated) | – | :3 |
| K3 / K2 | gobert2009 | Glomb | Glc + Lys, 50 C, OPD in situ; fed 3-DG | 1-DG 23.5 vs 3-DG 5.0 mmol/mol (deaerated); fed 3-DG half-life 40 h at 50 C | – | :47 |
| K1 | hicks1975 | Feather, Missouri | Amadori compounds, 14C and D2O, 100 C | D exchange during the reaction, so a 1-deoxyosone enol intermediate | – | :36-37 |
| K2 / K3 | hofmann2000b | DFA | ARP-Phe + OPD in situ, argon, 85 C | Under argon the ARP goes mainly to 1-deoxyosone | K (Strecker rows) | :107 |
| K2 / K4 | zhai2019, zhai2021, wei2022, zhang2025, deng2021, xia2022, tang2019 | Ho network | fed TTCA / ARP, water, pH 5.5-8, 100-130 C | Pentose: 3-DX > 1-DX standing levels (2-2.6x); TTCA's 3-DX share falls steeply with pH; peptide amines raise both 2-8x; free Met does not catalyse ARP decay; phosphate 1.8x | zhai2019 / zhai2021 / zhang2025 / deng2021 K | "3-DX > GO > 1-DX at every time, with 3-DX / 1-DX = 2.1-2.6" zhai2019:39 |
| K4 (fitted) | desclaux2006 | Reading | Xyl + Gly, 0.5 M phosphate pH 6, 80-120 C | Fitted ARP -> 1-DX : 3-DX flux split ~50 : 1, OPPOSITE to the standing levels (uncalibrated) | k | :4 |
| K2 / K4 | zyzak1995 | Baynes | fFL, 37 C, phosphate 10-200 mM vs HEPES | Phosphate catalyses Amadori decomposition on every channel | – | :31 |

**Verdict.**
- **Hexose hub: CAN.** 1-DG is short-lived and the split shifts to 1-DG with pH (Wageningen, Nestlé, DFA, Ledl, Glomb).
- **Pentose flux split in water vs pH: SINGLE** (the Ho network's standing levels; DFA at one pH). Reading's fitted split
  contradicts both.
- **Amine ACCELERATION of 2,3-enolisation / beta-cleavage: NONE.** Davidek and Smuda only exclude an amine REQUIREMENT.
  Voigt & Glomb 2009 (now verified, below) is the lead.
- **Phosphate effect on the split: NONE.**

## 2.2 The 3-DX sink

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K2 | zhai2019 | Ho network | +Cys cuts 3-DX to 20 % (vs 1-DX less): Cys traps the aldehydic 3-DX | K | :45 |
| K2 / K5 | zhou2022 | Ho network | Cys cuts 3-DP ~5-8x and 1-DP ~3x; ~11 mM of ARP carbon missing (an unmeasured Cys adduct) | K | :44 |
| K2 | wei2022 | Ho network | 1-DX responds first: CONFLICTS with 3-DX being the preferred target | – | :32 |
| K2 | zhai2021 | Ho network | -Xyl / +Xyl double ratio of 3-DX / 1-DX = 0.85 / 0.92 / 1.22 / 2.7 at pH 5.5 / 6 / 7 / 8 (reviewer arithmetic) | K | docs/history/mechanism_reviews_2026-10-31/hub.md:16 |
| K2 | wedzicha1991, edwards1992 | Leeds | Cys binds 3-DG 1 : 1 within minutes, then slowly releases thiol; glycine catalyses 3-DG activation | – | wedzicha1991:3; edwards1992:27 |
| K2 | wondrak2002 | Arizona | Cys + phenylglyoxal k2 = 0.63 M-1 s-1 at 37 C (a model alpha-oxoaldehyde) | – | :3 |
| K2 | lu2019 | Ho network | MS evidence of a GSH + deoxypentosone adduct | – | :43 |
| K2 | gobert2009 | Glomb | Fed 3-DG + Lys leaves without making smaller dicarbonyls (an amine route, not fragmentation) | – | :64 |
| K2 / K3 | weenen1998 | Quest | 3-DG -> MGO is amine-independent (56 vs 57) | – | :98 |
| K2 | davidek2005 | Nestlé | 3-DG gives hardly any acetic acid | – | :101-102 |
| K2 | mittelmaier2010 | Erlangen | Fed 3-DG -> 3,4-DGE, reversible (the dehydration route) | – | :45 |
| K1 | cerny2004, zhou2025b | Nestlé; Ho network | FFT carbon via 3-deoxypentosone -> furfural; furfural mostly from the ARP's own 3-DX | – / K | cerny2004:29; zhou2025b:46 |

**Verdict.**
- **That thiols trap 3-deoxyosones: CAN** (Ho network for pentose, Leeds for hexose).
- **That 3-DG leaves by adducts and dehydration rather than cleavage: CAN** (Glomb, Nestlé, Quest, Erlangen).
- **The pH law of the Cys trap and the pentose 3-DX lifetime in water: SINGLE** (zhai2021). wei2022 conflicts.
  No fed 3-DX + Cys x pH experiment exists.

## 2.3 MFT routes (and FFT)

| kind | dossier | lab | conditions | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|---|
| K1 | cerny2003 | Nestlé | [13C5]ribose + Cys, pH 5, 95 C, 4 h | MFT 95 % intact C5: C2 + C3 recombination rejected at 95 C | – | "pathways via ribose fragmentation were not relevant." :64 (checked) |
| K1 + K2 | cerny2003 pot D / E | Nestlé | NF vs [13C5]ribose; furfural vs ribose | Only 7 % of MFT via norfuraneol; FFT 92 % from fed furfural | – | :69, :195 |
| K1 | cerny2004 | Nestlé | [1-13C]ribose | Ribose C-1 -> MFT 2-methyl; C2 + C3 "not relevant" | – | :30 |
| K1 | zhao2019 | Xie, BTBU | [13C5]-TTCA + Gly + Cys, 140 C | 97 % of MFT intact TTCA skeleton | – | :33 (checked) |
| K1 | cerny2007b / cerny2007 | Firmenich | [13C5]xylose + Cys + thiamine, 145 C, pH 4-7 | MFT 46 % xylose-derived, the rest thiamine | – | cerny2007b:3 |
| K1 (34S) | bolton1993 | Minnesota + Nestlé | bouillon, 34S-Cys, 120 C | MFT sulfur 92 : 8 thiamine : Cys | – | :3 |
| K1 | xu2013 | Song, BTBU | [13C5]xylose + Cys + thiamine, 140 C | No xylose label in MFT / FFT (contradicts Cerny; detection limit suspected) | – | :3 |
| K1 | liu2023b | lab unrecorded | 168 C CAMOLA | "19 % fragments", re-read as an EI-MS [M-H]+ artefact | – | :27 |
| K2 / K4 | hofmann1998_reconciliation | DFA | fed precursors + H2S, pH 3 / 5 / 7, 145 C; in-situ NF / MFT | C2 + C3 is the strongest FED route at 145 C; NF / MFT ratio 2754x (ribose) vs 7x (glucose): NF not the sole precursor; in-situ mercapto-2-propanone | **K** (thiol_submodel) | :328, :232, :408 |
| K2 | whitfield1999 / 2001 | Reading | fed NF + Cys or H2S, pH 4.5 / 6.5, 140 C | NF -> MFT needs only H2S; none at pH 6.5 | – | :95; :27 |
| K2 | deroos2005 | Givaudan | Cys + ribose / xylose / NF, ~100 C | NF gives ~48x LESS MFT than ribose: second lab against the NF route near 100 C | – | :144 |
| K4 | schieberle2000 | DFA | ribose + Cys 100 C, 0.5-12 h, SIDA | MFT rises 35x from 30 to 360 min, still rising at 12 h | **K** | :31 |
| K1 (2H) + K2 | poisson2009, cerny2021 | Nestlé; Firmenich | coffee + labelled furfural / furfuryl alcohol | In roast, FFT does not come from furfural (0 % label); it comes from furfuryl alcohol (d2-FFA -> d2-FFT) | – | poisson2009:3; cerny2021:3 |

**Verdict.**
- **Intact C5 vs C2 + C3 for pentose at 95-140 C: CAN** (Nestlé K1, BTBU K1; DFA fed pots show feasibility only).
- **Norfuraneol vs non-NF intact C5 at 95-100 C: CAN** (Nestlé K1 + K2, Givaudan K2).
- **The same at 140-170 C: SINGLE.** Only DFA's in-situ ratios and fed pots exist; there is no labelled NF-vs-sugar
  competition.
- **Thiamine route: CAN** (Minnesota / Nestlé 34S, Firmenich 13C). Caveat: Bolton includes Nestlé co-authors, and
  Firmenich is Cerny.
- **FFT: CAN, by regime** (furfural + H2S in water; furfuryl alcohol in roast).
- **The time-resolved route mix: NONE.**

## 2.4 MFT loss

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K4 + K5 | vanseeventer2001 | Quest | MFT lost 59 %/day at 50 C, the same in air and argon; not disulfide; electrophilic coupling | – | :52 (checked) |
| K5 | whitfield1999 | Reading | ~35 % of MFT in disulfides at 140 C; the oxidant is the dicarbonyl redox couple, not air | – | :61 |
| K2 | guth1995 | DFA | MFT-MFT splits back to thiols at 100 C: the disulfide is a reservoir | – | :3 |
| K2 + K4 | hofmann2001 (= hofmann2002) | DFA | FFT bound covalently to coffee melanoidin (~80 %); FFT only | – | :3 |
| K2 | mottram1996 | Reading | Protein reduces furyl disulfides back to thiols during SDE | – | :74 |
| K2 / K3 | jhoo2002 | Ho network | Thiamine's pyrimidinyl cation consumes MFT (MAMP) | – | :116 |
| K3 (identification) | cerny2013 | Firmenich | MFT-S-Cys / FFT-S-Cys conjugates isolated; no yield | – | :3 |
| K4 | Zhang2024, zhou2023 | BTBU Cao; Ho network | MFT-MFT 54 % of the MFT pool at 115 C; disulfide share 6.5-9.6 % (derived) | – / K (zhou2023 thiols) | Zhang2024:128; zhou2023:206 |
| K4 | liu2023b, us10076126 | ?; Nestec | MFT halves 20 -> 60 min at 168 C; rises with upward curvature to 4 h at 110 C (labelled standard added BEFORE heating) | – / k | liu2023b:36; us10076126:3 |
| matrix | elmore2002 | Reading | PUFA esters cut MFT 14 to > 80x | – | :36 |

**Verdict.**
- **Disulfide as a reversible reservoir vs a terminal sink: CAN** (DFA guth1995, Quest mass balance, Reading).
- **MFT lifetime at 95-145 C: NONE.** Only bounds exist (schieberle2000's rise; 50 C storage). The loss law behind
  finding 4 of the 10-31 reviews cannot be fitted from the literature.
- **Adduct sinks (melanoidin, protein, thiamine, Cys conjugates): SINGLE per adduct type.**

## 2.5 Cysteine fate / TTCA

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K4 | cao2017 / zhao2017b | Xie, BTBU | K_app = [TTCA] / ([Cys][Xyl]) constant at 58.5-64.6 L/mol over 6 times: a fast pre-equilibrium; TTCA drains to Cys-ARP | K / K | cao2017:41-42 |
| K4 | deroos2005 | Givaudan | Thiazolidine drains into the Cys-Amadori compound (35 % maximum, 90 C, pH 5) | – | :40 |
| K4 | zhai2019 | Ho network | TTCA relaxes in 20 min at 70 C (83 mM pot) | K | review cysteine.md:5 |
| K2 (+K1) | zhai2021 | Ho network | Added xylose RAISES free Cys in dilute pots: far from equilibrium | K | review cysteine.md:3 |
| K2 | zhai2023jafc | Ho network | Added Cys delays ARP and osones: TTCA <-> Cys + Xyl is reversible | K | :1150 |
| K1 | zhai2023foodchem, zhao2019 | Ho network; BTBU | Exogenous 13C-xylose enters only above 100 C; TTCA's C5 reaches MFT intact | K / – | zhai2023foodchem:380 |
| K2 | kang2026_SI | Ho network (Zhai) | Ea of free-Cys depletion from TTCA = 55.1 kJ/mol | – | SI:603 |
| K2 | wei2022 | Ho network | Most Cys loss is NOT TTCA (a fast loss, flat in temperature) | – | :47 |
| K4 | li2011 | Chen, Jiangnan | Free-Cys first-order loss constants at 100-130 C | – | :3 |
| K4 | zheng1994 | Ho network (Rutgers) | H2S from Cys slow, Ea ~130 kJ/mol | – | :209 |
| K2 | mizutani1960 / obata1960 | Hokkaido | Thiazolidine-bound Cys releases H2S only 1.2-3.5x slower than free Cys | – | :3 |
| K4 | fife1991, kallen1971 | USC; Penn | Ring opening fast, hydrolysis rate-limiting and general-base (phosphate) catalysed; K(pH) rises ~20x from pH 5.6 to 7.5 (formaldehyde) | – | fife1991:3; kallen1971:3 |

**Verdict.**
- **TTCA <-> Cys-ARP isomerisation in concentrated pots: CAN** (BTBU + Givaudan).
- **H2S from Cys and thiazolidines is minor and slow: CAN** (Rutgers + Hokkaido).
- **Why dilute pots stay off equilibrium (catalysed hydrolysis vs K(pH) vs adducts vs O2): SINGLE** (Zhai).
  No H2S time course in a sugar + Cys pot; no cystine / total-thiol balance.

## 2.6 Strecker (donor identity; the O2 acid route)

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K2 | hofmann2000 T4 | DFA | Aldehyde O2-independent (air / argon 0.85-0.93); acid 1.6-5.5x higher in air; donor order MGO > GO > 3-DG > glucosone | K | :112 (checked) |
| K1 + K2 | hofmann2000 T5 | DFA | Fed [13C2]phenylacetaldehyde: < 1 % becomes acid, so the acid comes by a SEPARATE route | – (T5 not keyed) | :132 |
| K2 | hofmann2000b | DFA | ARP oxidation gives the aldehyde without free dicarbonyl (9x in air) | K | :28 |
| K2 | weenen2001 | Quest | Same donor order, but 21 % yield from MGO (~20x Hofmann's); ARP 270x glucose at pH 3.2 (SDE during reaction) | – | :24 |
| K2 | granvogl2006 | DFA | Amino-acid identity with MGO; the Strecker amine is negligible in water | – | :37 |
| K2 | hidalgo2013 | Seville | alpha-Keto acids act as Strecker donors | – | :2 |
| competition | martin2001 | Reading | Mixed amino acids are not additive (2.4-5.2x falls) | – | :65 |

**Verdict.**
- **Donor ranking: CAN** (DFA + Quest). The magnitude is unresolved: a ~20x cross-lab gap (CF-073), plausibly the
  observation operator (SDE stripping).
- **Amino-acid identity at a fixed donor: SINGLE** (DFA x 2 papers).
- **O2 acid route: SINGLE** (DFA).
- **Direct oxidative ARP route on the O2 question: SINGLE.**

## 2.7 Pyrazines (aminoketone, pyrrole sink, methylpyrazine extra route)

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K2 + K5 | adams2008 | Ghent (+ Wageningen) | Ala + MGO: 2,5-diacetyl-3-methylpyrrole up to 3.7 mol % vs 2,5(6)-DMP 0.02 %, so an aminoketone -> pyrrole sink | – | :32 |
| K1 | amranihemaimi1995 | Nestlé | Methylpyrazine and 2,5-DMP carry 0 % amino-acid carbon; 2,6-DMP 25-30 %; trialkyl 20-100 % (aldehyde addition) | – | :46 (checked) |
| K1 | low2007 | Reading | [2-13C]glycine: "x + x + y" formaldehyde addition up to ~45-50 % of 2,6- / 2,3-DMP and TMP | – | :30 |
| K1 (15N) + K2 | zhou2024 | Ho network | Exogenous Ala gives ~55 % of pyrazine N but ~37 % of methylpyrazine N | k (pyrazines not keyed) | :112 |
| K1 (13C) | zhou2025b | Ho network | One GO pool, 6-8 % xylose-derived; a late M+3 methylpyrazine channel | K (DO / DC); pyrazines k | :194 |
| K1 (15N) | hwang1995 | Ho network (Rutgers) | Ring N per amino acid; shares in figures only | – | :76 |
| K2 | cerny1994 | DFA | Fed MGO + Ala -> ethyldimethylpyrazines via acetaldehyde addition | – | :27 |
| K2 | vanlancker2012 | Ghent | Lys dipeptides 4-34x the DMP of free amino acids (areas) | – | :31 |

**Verdict.**
- **Carbon origin and aldehyde addition: CAN** (Nestlé + Reading K1).
- **Nitrogen source: CAN** (Rutgers + Jiangnan, but both are the Ho network). Strictly SINGLE-network; flag it.
- **Aminoketone -> pyrrole sink: SINGLE** (Ghent).
- **Methylpyrazine extra route: NONE** (the Leahy ratio test is a fitted-time-course signature, not K1-K4).
- **Aminoketone time course: NONE.**

## 2.8 Oxygen routes (Wolff, oxidative alpha-cleavage, Namiki C2)

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K1 (18O2) | davidek2006b | Nestlé | Total acid unchanged without O2; under O2 ~1/3 of acetic acid comes by oxidative alpha-cleavage of 1-DG | K | :26-27 |
| K3 | gobert2009 | Glomb | O2 diverts flux from 1-DG to glucosone (glucosone 44.5 vs 5.5 mmol/mol) | – | :26 (checked) |
| K2 / K3 | hofmann2000b | DFA | ARP + air + Cu -> glucosone (Hodge type); argon -> 1-DG | K | :149 |
| K2 | gobert2009 | Glomb | Fed glucosone is not a major glyoxal source at 50 C | – | :63 |
| K2 | zyzak1995, smuda2010 | Baynes; Glomb | Oxidative ARP channel ~4x (chelator-confounded); small O2 effect on amides | – | zyzak1995:25 |
| K4 | hayashi1980 | Namiki, Nagoya | Amine-bound C2 appears before 3-DG (refluxing ethanol) | – | :25 |
| K2 / K3 | weenen1998 | Quest | Glyoxal "cannot be explained" from deoxyosones or ARPs, so a pre-Amadori source | – | :33-34 |

**Verdict.**
- **O2 diverts the ARP / 1-DG flux to glucosone: CAN** (Glomb + DFA; both carry metal or chelator confounds).
- **Oxidative alpha-cleavage to acids: SINGLE** (Nestlé).
- **Namiki C2:** two labs, both indirect (Quest K2 aqueous, Namiki K4 in ethanol). Treat it as SINGLE for an aqueous
  quantitative rate.
- **Wolff free-sugar autoxidation vs ARP oxidation as the glucosone source: NONE.**

## 2.9 Redox (1,4-dideoxyosones, reductones)

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K3 | nedvidek1992 | Ledl | 1,4-Dideoxyosones form in higher yield with alpha-amino acids than with beta-alanine (Strecker-coupled reduction); no numbers, temperature unstated | – | :241 (checked) |
| K3 | hollnagel2000 | Kroh, TU Berlin | 1,4-Dideoxyhexosulose < 1 mol % even without amine (dry, 100 C) | – | :20 |
| K1 | cerny2004 | Nestlé | The label pattern is consistent with a 1,4-dideoxypentosone route to MFT (proposed only) | – | :96 |
| K3 / K4 | rizzi2010 | Rizzi | A mM reductone pool, identified as 1-deoxyosone, in amine-free aldose (phosphate, 100 C) | – | :39 |
| K2 (weak) | kobayasi1965 | Tokyo | Reductive H2S mechanism proposed, not shown | – | :31 |

**Verdict: SINGLE** (Ledl) for the Strecker-coupled reduction. The amine-free trace (Kroh) shows a second route but
does not test the coupling. **Needs designed data.**

## 2.10 Lipid carbonyls as Strecker donors

| kind | dossier | lab | what it shows | keyed | quote (file:line) |
|---|---|---|---|---|---|
| K2 | zamora2007, zamora2008, zamora2013, zamora2015, hidalgo2004, hidalgo2016 | Seville (all) | 2,4-Decadienal -> phenylacetaldehyde ~8 %/h at 180 C; hydroperoxides up to 12 % at pH 2.15; O2 lowers Ea 2.8x; aldehyde vs amine split depends on O2 | – | zamora2013:152; zamora2015:100 |
| negative control | mandin1999 | Reading | Linoleic ACID + Met at 100 C: no methional (the acid, not its carbonyls) | – | :96 |
| K1 / K2 (thiophenes, pyridines) | du2023, elmore1997, zhang1989, zhou2000, kim1998, zamora2020 | CAAS; Reading; Ho network; Kentucky; Seville | Lipid dienal + H2S -> 2-pentylthiophene (92 % unlabelled); dienal + NH3 -> 2-pentylpyridine | – | du2023:224; zhou2000:25 |

**Verdict.**
- **Strecker-donor role: SINGLE** (Seville only). hidalgo2005 and zamora2006, in the queue, are also Seville.
- **Lipid aldehyde + H2S -> thiophenes: CAN.**
- **Dienal + NH3 -> pyridines: CAN.**

## 2.11 Observation-operator artefacts

| kind | dossier | lab | what it shows | quote (file:line) |
|---|---|---|---|---|
| method comparison | schieberle1996 | DFA | Same pot: ether extraction destroys MFT; 2-vinylpyridine derivatisation and purge-and-trap give 17-21 ug/L; FFT within +-7 % | :37 |
| spike recovery | schieberle1998 | DFA | Purge-and-trap SIDA recovery 95-104 %; MFT and its standard are completely degraded by distillation work-up | :66 |
| method comparison | mottram1996 | Reading | Thiol lost on mixing or extraction; protein converts disulfide to thiol during SDE | :72 |
| OPD timing | gobert2009 | Glomb | OPD in situ measures cumulative flux; added at sampling, dicarbonyls are "much lower" | :38 |
| OPD artefact | nakamura2020 | Murata | Norfuraneol + OPD forms quinoxalines at room temperature: NF fakes a dicarbonyl signal | :317 |
| spike control | nedvidek1992 | Ledl | Quinoxaline 11 is not an artefact of 32 | :353 |
| spike recovery | hamzalioglu2026 / buttery1986 / fan2018 | Gökmen / USDA / BTBU | OPD-LC-MS recovery 99.5-103.6 %; 2-AP through SDE only 28 %; SAFE recovery of MFT 0.94 / 0.89 | :74 / :3 / :69 |
| fibre competition | cerny2003, fischer2021 | Nestlé; Dijon | Excess furfural saturates the PDMS-DVB fibre; non-linear calibration in mixtures | cerny2003:196 |
| SIDA flaw | us10076126 | Nestec | Labelled standards added BEFORE 110 C heating, so they degrade | :3 |

**Verdict.**
- **Thiol work-up artefacts: CAN** (DFA + Reading).
- **OPD in situ vs ex situ: CAN** (Glomb + Murata). The NF -> quinoxaline artefact is SINGLE (Murata).
- **SDE recovery: CAN** for 2-AP and SAFE. For MFT itself the 4-5x SDE effect is cited only (Grosch & Zeiler-Hilgart,
  below).
- **SPME competition: NONE** (two passing observations).
- **The EI-MS label artefact (liu2023b): NONE.** It rests on reviewer arithmetic.

## 2.12 Verdict table

| mechanism | verdict | labs with K1-K4 | open rival question |
|---|---|---|---|
| hub split / lifetime, hexose | **CAN** | Wageningen, Nestlé, DFA, Ledl, Glomb, Feather | -- |
| hub split, pentose in water x pH; amine acceleration; phosphate | **SINGLE** / **NONE** | Ho network (+ DFA one pH; Reading contradicting) | flux split vs standing level; k0 + k_aa[AA] + k_P[P] |
| 3-DX sink (Cys trap) | **CAN** (existence) / **SINGLE** (pH law) | Ho network, Leeds | trap grows with pH? pentose 3-DX lifetime |
| MFT routes 95-140 C | **CAN** | Nestlé, BTBU, DFA, Reading, Givaudan, Minnesota | NF vs intact C5 at 140-170 C (**SINGLE**); time-resolved mix (**NONE**) |
| MFT loss | **NONE** (lifetime at process temperature); **CAN** (disulfide reversible) | DFA, Quest, Reading | first-order loss law and its O2 / carbonyl dependence at 95-145 C |
| cysteine fate / TTCA | **CAN** (concentrated); **SINGLE** (dilute off-equilibrium) | BTBU, Givaudan, Ho network, Rutgers, Hokkaido | catalysed hydrolysis vs K(pH) vs adducts; H2S time course |
| Strecker donor identity | **CAN** (ranking); magnitude conflict | DFA, Quest | 20x yield gap (operator?) |
| Strecker O2 acid route | **SINGLE** | DFA | cross-lab replicate |
| pyrazines: carbon origin | **CAN** | Nestlé, Reading | -- |
| pyrazines: pyrrole sink | **SINGLE** | Ghent | aminoketone partition rates |
| pyrazines: methylpyrazine extra route | **NONE** | -- | GO + aldehyde vs aminoketone dimer |
| O2 -> glucosone | **CAN** | Glomb, DFA | -- |
| O2: oxidative alpha-cleavage acids | **SINGLE** | Nestlé | second lab (voigt2010 in queue) |
| O2: Namiki C2 | **SINGLE** (aqueous) | Quest (indirect), Namiki (ethanol) | rate from the Schiff base in water |
| O2: Wolff vs ARP oxidation | **NONE** | -- | free sugar +/- amine under air with a trap |
| redox, 1,4-dideoxyosones | **SINGLE** | Ledl (+ Kroh amine-free) | quantified alpha- vs beta-amino acid x O2 |
| lipid carbonyls as Strecker donors | **SINGLE** | Seville | second lab; lipid vs sugar donor in one pot |
| operators: thiol work-up, OPD timing | **CAN** | DFA, Reading, Glomb, Murata | -- |
| operators: SPME competition, SDE for MFT, EI-MS label | **NONE** | -- | method-comparison pots |

---

# 3. Missing experiments with the highest value (each discriminates at least two mechanisms)

Ranked by (mechanisms discriminated) x (P that the result changes a structural choice) / cost. Probabilities are my
judgement. The owner said "no designed experiment for now", so these go to the designed-experiment brief (`rank` verb),
not to the bench.

**E1. Pentose hub ladder with matched in-situ and ex-situ OPD.**
- **Discriminates:** hub split / lifetime (pentose); Cys 3-DX trap and its pH law; phosphate catalysis; OPD operator.
- **Pots:** xylose 50 mM + glycine 50 mM in water under pH-stat.
- **Factors:**
  - pH 5.5 / 6.5 / 7.5;
  - phosphate 0 vs 0.2 M at the held pH (pH-stat in both arms, so buffer is not confounded with pH);
  - Cys 0 / 5 / 25 mM;
  - 100 C, plus a 120 C corner. That is 18 pots + 2 corners, in duplicate.
- **Sampling:** 0, 2, 5, 10, 20, 40, 80, 160, 240 min (quench on ice).
- **Analytics:**
  - 1-DX, 3-DX, GO, MGO by ex-situ OPD-LC-MS/MS with quinoxaline standards;
  - a parallel in-situ-OPD pot per condition for cumulative flux, so the split is a FLUX, not a standing level;
  - furfural and norfuraneol (HPLC-DAD);
  - acetic / formic acid (SIDA);
  - free Cys, TTCA, Cys-ARP, Gly-ARP (LC-MS);
  - an NF + OPD blank (the nakamura2020 artefact).
- P(changes the hub structure) ~0.6.

**E2. Fed 1-DX / 3-DX (or 1-DG / 3-DG) decay with and without amine and phosphate.**
- **Discriminates:** amine acceleration of beta-cleavage (B7a); deoxyosone lifetimes; 3-DX sink identity (Cys vs amine
  vs dehydration); phosphate.
- **Pots:** synthetic 1-DG and 3-DG (and 3-DX if it can be made) at 2 mM.
- **Factors:**
  - glycine 0 / 20 / 100 mM;
  - phosphate 0 vs 0.1 M (pH-stat 7.0);
  - +/- Cys 5 mM;
  - 100 C.
- **Sampling:** 0-60 min in 8 points.
- **Analytics:** substrate (OPD ex situ), acetic / formic / glyceric acid (SIDA), GO / MGO, HMF or furfural, adducts by
  LC-HRMS.
- Cheap: the substrates are published syntheses (Glomb, Davidek). P ~0.55.

**E3. Spiked labelled MFT in a running ribose + Cys pot, with CAMOLA.**
- **Discriminates:** MFT loss law (lifetime, O2 and carbonyl dependence); MFT route mix over time (NF vs intact C5 vs
  C2 + C3); thiol observation operator (standard added before vs after heating).
- **Pots:** ribose (1 : 1 12C / [13C5]) 100 mM + Cys 33 mM in 0.5 M phosphate pH 5.0 (the hofmann1998 / schieberle2000
  recipe).
- **Factors:**
  - 100 vs 140 C;
  - headspace air vs argon;
  - [2H3]-MFT spike at t = 0 vs at sampling only.
- **Sampling:** 15, 30, 60, 120, 240, 480, 720 min.
- **Analytics:**
  - MFT / FFT / NF / mercaptoketones by SIDA HS-SPME-GC-MS with [M-H]+ correction;
  - total thiol after TCEP reduction (disulfide pool);
  - MFT-S-Cys by LC-MS/MS.
- Decides the 10-31 finding 4 (the MFT-loss lead). P ~0.65.

**E4. Cysteine and sulfur mass balance with an H2S time course.**
- **Discriminates:** cysteine fate (dilute off-equilibrium: catalysed hydrolysis vs K(pH) vs adducts vs O2); H2S supply
  to the MFT / FFT routes; cysteine's own Strecker (mercaptoacetaldehyde).
- **Pots:** xylose + Cys at 10 / 10 mM and 83 / 165 mM (dilute and concentrated, as Zhai).
- **Factors:**
  - water vs 0.2 M phosphate (pH-stat 5.5 and 7.4);
  - air vs argon.
- **Sampling:** 0-120 min in 8 points at 100 and 120 C.
- **Analytics:**
  - headspace and dissolved H2S (3-buten-2-one derivatisation SIDA, the li2016b method);
  - free Cys, cystine, TTCA, Cys-ARP (LC-MS/MS);
  - 1-DX, 3-DX;
  - MFT, FFT;
  - a total-S closure.
- Fills the reviews' gap "no H2S time course in a sugar + cysteine pot". P ~0.55.

**E5. Forced-input Strecker / pyrazine pot, cross-lab replicate with an operator arm.**
- **Discriminates:** Strecker donor flux vs yield; the O2 acid route (second lab for hofmann2000 T4 / T5); aminoketone
  -> pyrazine vs -> pyrrole partition; methylpyrazine extra route; the SDE-vs-closed-vial operator (the CF-073 20x gap).
- **Pots:** MGO 10 mM or 2,3-butanedione 10 mM + an amino-acid mix (Phe, Leu, Ala, Gly, 5 mM each).
- **Factors:**
  - pH 3.2 and 7.0 (citrate-phosphate at matched ionic strength);
  - argon vs air;
  - closed vial vs open reflux under Likens-Nickerson;
  - one arm with [13C2]glyoxal 2 mM;
  - one arm with 15N-Ala replacing Ala.
- **Sampling:** 5, 15, 30, 60, 120 min at 100 C.
- **Analytics:**
  - Strecker aldehydes and acids (SIDA);
  - aminoacetone / aminoketones (derivatised LC-MS/MS);
  - pyrazines with isotopomers (GC-MS, corrected);
  - 2,5-diacetyl-3-methylpyrrole;
  - residual MGO (OPD).
- P ~0.6; feeds roadmap step 2 directly.

**E6. Free sugar vs ARP under controlled O2, with alpha- vs beta-amino acid.**
- **Discriminates:** Wolff vs ARP oxidation as the glucosone / glyoxal source; Namiki C2 timing (glyoxal before 3-DG);
  redox (1,4-dideoxyosones, Strecker-coupled reduction); oxidative alpha-cleavage acids (second lab after Davidek
  2006b).
- **Pots:** glucose 100 mM alone, + glycine, + beta-alanine, and the fed Gly-ARP.
- **Factors:**
  - argon / air / O2 headspace;
  - DTPA 0 vs 0.1 mM, crossed so metals are not confounded with O2;
  - pH-stat 7.0;
  - 100 C.
- **Sampling:** 2-240 min in 9 points.
- **Analytics:**
  - in-situ OPD for glucosone, 1-DG, 3-DG, 1,4-dideoxyosones, GO, MGO (LC-MS/MS);
  - acids (SIDA);
  - reductone titre (DCIP, as rizzi2010);
  - dissolved O2 (probe).
- P ~0.45.

**E7. Lipid vs sugar dicarbonyl as Strecker donor in one pot, with a thiol trap arm (second lab for Seville).**
- **Discriminates:** lipid carbonyls as Strecker donors (cross-lab); lipid carbonyls as thiol traps (the elmore2002
  14-80x MFT drop); the O2 dependence of both.
- **Pots:** Phe 10 mM + (E,E)-2,4-decadienal 2 mM, + MGO 2 mM, or both, in 0.1 M phosphate pH 6 with Tween-stabilised
  dispersion.
- **Factors:**
  - argon vs air;
  - +/- MFT 0.05 mM spike;
  - 100 and 140 C.
- **Sampling:** 15-120 min in 5 points.
- **Analytics:** phenylacetaldehyde and phenylethylamine (SIDA); residual decadienal and MGO; MFT (SIDA); MFT-dienal
  adducts (LC-HRMS); hexanal and 2-pentylfuran.
- P ~0.4 for the odorant set of plant-based meat.

**E8. Water-activity x thermal-history block on one recipe.**
- **Discriminates:** water activity (no deoxyosone ladder exists); thermal history (no odorant contrast exists); the
  C2 + C3 MFT route's low-water enhancement (DFA 5.8x).
- **Pots:** ribose 50 + Cys 25 + Gly 25 mmol/kg on microcrystalline cellulose.
- **Factors:**
  - a_w 0.3 / 0.6 / 0.9 vs aqueous;
  - isothermal 120 C vs a logged ramp 20 -> 145 C in 20 min, matched on equivalent time at 120 C with the measured
    barrier;
  - one [13C5]ribose CAMOLA arm.
- **Sampling:** 5 points to 60 min, core temperature logged.
- **Analytics:** 1-DX / 3-DX (ex-situ OPD), MFT / FFT / pyrazines (SIDA), label shares.
- P ~0.35. Lower, because product-level matrix work is out of scope until designed data exist anyway.

---

# 4. Papers the dossiers cite that would fill a gap, NOT on disk and NOT in the private download queue

Each candidate's absence from the main repo's `data/articles` (492 files) was checked by name. DOIs were verified
against Crossref today (title, authors, journal, pages matched); "lead" means not verified. Already queued and so left
out: Smuda 2010 group tier 1 (on disk), and tier 2 / 3 (pesek1975, schubert1936, zhou2024b, cui2022, voigt2010, hou2026,
thornalley1999, guan2012, pripisnicolau2000, hidalgo2005, zamora2006, scalone2015, vanlancker2010, lu2005, adams2001,
lo1994, shibamoto1977, vandenouweland1975, huber1990, and others).

| paper | DOI | verified | gap it fills | dossier line relied on |
|---|---|---|---|---|
| Voigt & Glomb 2009, J. Agric. Food Chem. 57:4765-4770, "Reactivity of 1-deoxy-D-erythro-hexo-2,3-diulose: a key intermediate..." | 10.1021/jf900459x | Crossref | Amine-free vs amine 1-DG cleavage: the amine-acceleration (B7a) question, 2.1 | "The lead that could answer acceleration is Voigt & Glomb 2009" smuda2010:55-56 (the queue lists it as "DOI to verify"; now verified) |
| Beck, Ledl & Severin 1988, Carbohydr. Res. 177:240-243, "Formation of 1-deoxy-D-erythro-2,3-hexodiulose from Amadori compounds" | 10.1016/0008-6215(88)85058-4 | Crossref | 1-DG : 3-DG split at pH 7 (hub, 2.1) | martins2003b:462-464; "It has not been read." blank1997:744 |
| Huyghues-Despointes & Yaylayan 1996, J. Agric. Food Chem. 44:672-681, "Retro-aldol and redox reactions of Amadori compounds: mechanistic studies with variously labeled D-[13C]glucose" | 10.1021/jf9502921 | Crossref | K1 labelling for the REDOX mechanism (2.9) and ARP retro-aldol (2.1) | docs/history/mechanism_reviews_2026-10-31/hub.md:81 ("Unverified content"); ARP rates cited at jousse2002:271 |
| Kallen 1971, J. Am. Chem. Soc. 93:6236-6248, "Mechanism of reactions involving Schiff base intermediates. Thiazolidine formation from L-cysteine and formaldehyde" | 10.1021/ja00752a040 | Crossref | Forward / reverse thiazolidine RATE constants (2.5). The on-disk kallen1971 is the equilibria paper, 93:6227, 10.1021/ja00752a039 | "Not on disk; this is the paper that would carry forward/reverse rate constants." kallen1971:19 |
| Charles-Bernard, Kraehenbuehl, Rytz & Roberts 2005, J. Agric. Food Chem. 53:4417-4425, "Interactions between volatile and nonvolatile coffee components. 1. Screening..." | 10.1021/jf048021q | Crossref | Which matrix fraction binds thiols (MFT loss, 2.4). The on-disk charles-bernard2005.pdf is Part 2 | "Part 1 of the pair ... -- not on disk" charlesbernard2005:28 |
| Ghiron, Quack, Mawhinney & Feather 1988, J. Agric. Food Chem. 36:677-680, "Studies on the role of 3-deoxy-D-erythro-glucosulose (3-deoxyglucosone) in nonenzymic browning..." | 10.1021/jf00082a002 | Crossref (title) | 3-DG as a Strecker donor (2.6, 3-DG sink 2.2) | docs/history/mechanism_reviews_2026-10-31/strecker.md:113 |
| Grosch & Zeiler-Hilgart 1992, ACS Symp. Ser. 490, ch. 14, 183-192, "Formation of meatlike flavor compounds" | 10.1021/bk-1992-0490.ch014 | Crossref | The 4-5x SDE vs direct-extraction MFT effect (operator, 2.11); thiamine / Cys as MFT source | "a 4-5 fold increase in recovered MFT by Likens-Nickerson distillation-extraction over direct extraction" bolton1993:203 |
| Weenen & Tjan 1992, ACS Symp. Ser. 490, ch. 17, 217-231, "Analysis, structure, and reactivity of 3-deoxyglucosone" | 10.1021/bk-1992-0490.ch017 | Crossref | Cysteamine trapping of 3-DG (K3, 2.2). The dossier cites it as "Weenen & Tjan 1994" | "(Weenen & Tjan 1994 with cysteamine trapping; NMR)" weenen1998:114 |
| Kim, Hartman & Ho 1996, J. Agric. Food Chem. 44:3906-3908, "Formation of 2-pentylpyridine from the thermal interaction of amino acids and 2,4-decadienal" | 10.1021/jf950623+ | Crossref (title) | Companion to kim1998 (lipid carbonyl + amine, 2.10). The on-disk kim1996 is a different paper | "The oil-medium companion to Kim, Hartman & Ho 1996 (aqueous)" kim1998:2 |
| Low, Koutsidis, Parker, Elmore, Dodson & Mottram 2006, J. Agric. Food Chem. 54:5976-5983 | 10.1021/jf060328x | Crossref | Concentrations behind low2007's label shares (2.7); queued as "optional, DOI to verify", now verified | "Low 2006 (not on disk)" low2007:162 |
| Mottram, Nobrega & Dodson 1996, Flavour Science, 413-418, "Changes in thiol and disulfide flavour compounds resulting from the interaction with proteins" | 10.1533/9781845698232.7.413 | Crossref (title) | Thiol / disulfide exchange with protein (2.4). Lead that it is the "Mottram et al., 1998" whitfield1999 cites: the years differ | "thiols are not converted to disulfides during the analysis procedures ... (Mottram et al., 1998)" whitfield1999:60 |
| Weenen & Apeldoorn 1996, Flavour Science (RSC) 211-216 | -- | lead | Primary data behind weenen1998 (hub / 3-DG / glyoxal) | weenen1998:21 |
| de Roos 1992, ACS Symp. Ser. 490:203-216 | -- | lead (Crossref search failed) | Thiazolidine -> Cys-Amadori kinetics (2.5) | deroos2005:20 |
| Mottram & Leseigneur 1990 (chapter) | -- | lead | Norfuraneol consumption vs pH in ribose-Cys (2.3 / 2.4) | "(Mottram & Leseigneur 1990, norfuraneol decreasing with pH in ribose-cysteine)" apriyantono1993:461 |
| Desclaux 2006 PhD thesis | -- | lead | The full pentose network parameters, to adjudicate its 50 : 1 split (2.1) | desclaux2006:208 |
| Cerny & Davidek follow-up CAMOLA with norfuraneol + cysteine | -- | lead (may never have been published) | The NF share of MFT at 95 C, directly (2.3) | "'CAMOLA experiments with 4-hydroxy-5-methyl-3(2H)-furanone/cysteine ... Further experiments are currently under way'" cerny2003:357 |

**Data to ask authors or SI for (not papers):**
- Zhai 2021: the 13C5-xylose label fraction in TTCA (MECHANISM_REVIEWS sec. 4).
- Xia 2022 SI Fig. S5 (120 C 1-DG / 3-DG; xia2022:231).
- Zhou 2024 SI Table S1 (glycine isotopomers; zhou2024:228).
- Zamora 2015 SI Figs. S5-S6 (zamora2015:69).
- Wei 2022 SI Fig. S3 is already asked.

**Housekeeping found on the way (no action taken):**
- hofmann2001 and hofmann2002 dossiers are the same paper; so are Zhang2024 and zhang2024b, and huang2016 and huang2017.
- liu2023b's laboratory is not recorded in its dossier.
- The ledger-prefix grep shows zhang2025b's blank pot is keyed (route_submodel), but its lysine series is not.
