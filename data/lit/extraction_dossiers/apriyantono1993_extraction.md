# Apriyantono & Ames 1993 (`10.1002/jsfa.2740610416`) — single-paper extraction, Wave K5b, 2026-08-29

> Apriyantono, A.*; Ames, J. M.‡
> "**Xylose−Lysine Model Systems: The Effect of pH on the Volatile Reaction Products**"
> *J. Sci. Food Agric.* **1993**, 61 (4), 477−484.
> Department of Food Science and Technology, University of Reading, Whiteknights, PO Box 226,
> Reading RG6 2AP, UK.
> *Present address (AA): Department of Food Technology and Human Nutrition, Bogor Agricultural
> University, Kampus IPB Darmaga, PO Box 122, Bogor 16001, Indonesia. ‡ To whom correspondence
> should be addressed.
> Received 14 September 1992; revised 25 January 1993; accepted 26 February 1993.
> Funding: Indonesian government, Interuniversity Project World Bank XVII Scholarship (AA).

**Source PDF:** `data/articles/apriyantono1993.pdf` (8 pp., 748,351 bytes, on disk 2026-08-29 09:42).
**Read method:** OCR'd scan. `pdftotext -layout` recovered the body text but **mangled Table 1's
superscript exponents and several decimal points** (it rendered `12.9 × 10³` as `129 x 103`,
`3.53 × 10⁶` as `3.53 x 108` in one row, `66.2` as `662`, and `6.8` as `68`). **Every numeric cell
of Table 1 reproduced below was therefore re-read visually from a 260 dpi render of pages 3 and 4**,
and the totals were checked arithmetically against the component rows (§3.4). **Do not trust a
naive text-layer harvest of this PDF.**

---

## 0. VERDICT UP FRONT

### 0a. ★★ THE ROUND-3 "RATIO-ONLY" VERDICT IS WRONG. THIS PAPER CARRIES ABSOLUTE MOLAR YIELDS.

`research_round3_channels.md` §C.2 records, from the **abstract only** (`[ABS]`):

> ⚠️ **"522 and 999 g/kg are shares OF THE VOLATILE FRACTION, not absolute yields.** … **RATIO-ONLY**,
> and the 999 figure means the uncontrolled-pH system produced essentially nothing but furfural"

**That is true of the abstract's `g kg⁻¹` figures and FALSE of the paper.** Table 1 — 58 compounds ×
2 systems — reports every yield in **`nmol mol⁻¹ xylose`**, i.e. **an absolute molar yield on a
declared sugar basis**, quantified against an internal standard (*n*-tetradecane), with a stated
precision. The `g kg⁻¹` numbers in the abstract and Conclusion are a *second*, derived presentation
of the same data.

| | round-3 `[ABS]` read | **what the paper actually has** |
|---|---|---|
| basis | share of the volatile fraction | **`nmol mol⁻¹ xylose`** — absolute, molar, on the limiting sugar |
| quantification | not stated | **internal standard (*n*-tetradecane)**, HP 3390A integrator |
| precision | not stated | **±16 %**, stated in Table 1's footnote *d* |
| replication | not stated | **means of two runs**; isolates prepared **in at least triplicate** |
| 2-furfural, pH 2.6 | "999 g/kg of volatiles" | **3.53 × 10⁶ nmol mol⁻¹ xylose = 0.353 mol %** |
| 2-furfural, pH 5.0 | "522 g/kg of volatiles" | **12.9 × 10³ nmol mol⁻¹ xylose = 0.00129 mol %** |
| grand total volatiles | — | **3.53 × 10⁶ (pH 2.6) vs 2.47 × 10⁴ (pH 5.0) nmol mol⁻¹ xylose** |

**★ THE PAPER IS THEREFORE A LEVEL SOURCE, NOT A RATIO SOURCE, AND THE FURFURAL RATIO IT SUPPORTS IS
274× — NOT THE ≈2× A READER WOULD INFER FROM `999 / 522`.** The `g kg⁻¹` shares moved by 1.9× while
the underlying absolute yields moved by 274×, because **the denominator (total volatiles) collapsed
by 143× at the same time.** This is exactly the trap `research_round3_channels.md` §C.2 warned about
("a statement about the denominator collapsing as much as about furfural rising") — and the paper
itself resolves it. **The round-3 register entry should be corrected.**

### 0b. THE OTHER HALF OF THE BRIEF: ⚠️ **THERE IS NO DMHF IN THIS PAPER AT ALL**

The brief filed this paper under the DMHF/furanone cluster. **DMHF/HDMF is not among the 58
identified compounds, is not among the 14 unidentified ones, and is never claimed.** The only
furanone in the entire volatile profile is **compound 13, 4-hydroxy-5-methyl-3(2H)-furanone
(NORFURANEOL)**, at **`tr` (< 2 nmol/mol xylose) at pH 5.0** and **not detected at pH 2.6**.

The only mention of DMHF anywhere in the paper is in the Introduction, **citing Shu & Ho 1988** —
i.e. citing the other paper in this same wave. See §6.

**What the paper IS for this repo: a pH-trajectory validation pair, a furfural level source, and a
NORFURANEOL near-null that has to be reconciled against the `norfuraneol ≫ DMHF` ordering.** See
§0c and §8.

### 0c. ★★ THE PH-TRAJECTORY PAIR — the reason to keep this paper

> **"Aqueous molal solutions of xylose and lysine (initial pH 4.9) were refluxed either with control
> of the pH at 5.0 or without pH control (final pH 2.6)."** — abstract, verbatim

**One lab, one apparatus, one hour, one sugar, one amine, one initial pH — and the ONLY difference
between the two arms is whether the pH was held.** That is the exact control
`k3` §C.11 records as missing from Zhou 2023, and the exact control that
`blank1996_extraction.md` §1.2, `blank1997_extraction.md` §1.1 and `wang2008_extraction.md` §1.3
each had to flag as *absent* from their own papers. **In this cluster of five DMHF papers, four
report pH as an unheld set-point and one — this one — runs the held-vs-drifting comparison
explicitly.**

**Its use for the pH-state model, stated explicitly as the brief asked:**

> **This paper is the corpus's only single-lab, single-protocol HELD-pH vs DRIFTING-pH pair.** It
> should be given a **named validation role of its own** — not folded into a general "pH
> sensitivity" bucket — because it is the only dataset that can distinguish **"the model gets pH 5
> chemistry right"** from **"the model gets the pH TRAJECTORY right."** A model that carries an
> internal pH state must reproduce **(i)** the ≈143× rise in total volatiles when the pH is allowed
> to fall 4.9 → 2.6, **(ii)** the 274× rise in furfural specifically, **(iii)** the near-total
> extinction of the nitrogen-containing products (10.1 → 0.13 µmol/mol, §4.2), and **(iv)** the
> disappearance of norfuraneol, which forms only in the held arm. A model that predicts the pH 5.0
> arm correctly and the pH 2.6 arm wrongly has a **pH-trajectory** defect, and this is the only
> experiment in the corpus that can say so. **Proposed role: `HOLDOUT — pH-trajectory pair`,
> scored as a PAIR of arms, never as two independent rows.** See §11.

---

## 1. SYSTEM DEFINITION — verbatim (Experimental, p. 478)

### 1.1 Materials
> "D-(+)-Xylose (**99 %, Gold Label**) was obtained from the Aldrich Chemical Company Ltd
> (Gillingham, Dorset, UK) and **L-lysine monohydrochloride (chromatographically homogeneous)**,
> diethyl ether (AnalaR) and ***n*-tetradecane** from BDH Chemicals Ltd (Poole, Dorset, UK)."

### 1.2 ★ The two arms — verbatim, and this is the whole point of the paper
> "**Solutions of xylose (0·5 mol) and lysine monohydrochloride (0·5 mol) in degassed distilled
> water (500 ml) were refluxed for 1 h in a modified Likens and Nickerson apparatus** (Maarse and
> Kepner 1970). **Diethyl ether (15 ml) was the extraction solvent. The pH of the systems were
> initially 4·9 and it was either uncontrolled during heating (pH 2·6 system) or was adjusted to and
> maintained at pH 5·0 by the addition of 3 M NaOH solution before and during heating (pH 5·0
> system), the pH being monitored with a Russell ACWL autoclavable electrode** (Russell pH Ltd,
> Auchtermuchty, Fife, UK), connected to a Kent EIL 7045/46 meter."

*(⚠️ The `0·5 mol` values were re-read from the page image — the text layer rendered the xylose one
as `0 3 mol`. The image is unambiguous: **`xylose (0·5 mol)` and `lysine monohydrochloride
(0·5 mol)`.**)*

**Derived, arithmetic only: 0.5 mol in 500 mL = 1.0 mol/L xylose and 1.0 mol/L lysine·HCl, 1:1.**
The abstract's phrase *"aqueous molal solutions"* is consistent with this. **Reflux of an aqueous
solution → ≈100 °C at 1 atm. The paper never states a temperature; "refluxed" is all there is.**

> "Each solvent extract was concentrated to a volume of **0·5−1·0 ml using a Vigreux column under a
> reduced pressure of 8·0 kPa (60 mmHg)**. Concentration, to a final volume of **0·05 ml**, was
> achieved using a gentle stream of nitrogen. **Blank isolates were prepared using only distilled
> water in the sample flask. All isolates were prepared in at least triplicate and were stored at
> −20 °C prior to analysis.**"

### 1.3 ⚠️ THREE THINGS THAT LIMIT WHAT TRANSFERS

1. **The isolation is SIMULTANEOUS DISTILLATION−EXTRACTION (Likens−Nickerson), not a static
   headspace or a cold solvent extraction.** The sample is held at reflux *while* being continuously
   stripped into ether for an hour. **SDE is thermally aggressive and strongly biased toward
   volatile, thermally robust, ether-soluble species.** Polar, thermally labile, water-miscible
   compounds — **which is precisely what the 3(2H)-furanones are** (Blank 1997 p. 2642, verbatim:
   *"The quantitative analysis of dihydrofuranones is a challenging task because of their rather
   fragile nature, water miscibility, and delicate gas chromatographic behavior"*) — are recovered
   poorly and are **partly generated or destroyed inside the apparatus.** This is the primary reason
   the norfuraneol result here cannot be compared directly with Blank's (§8.2).
2. **No recovery correction.** One internal standard (*n*-tetradecane, a C14 alkane) is used as a
   response reference for 58 chemically diverse compounds. **No compound-specific response factors
   are reported.** The paper does not claim otherwise.
3. **⚠️ THE HELD ARM IS NOT ISOTONIC WITH THE DRIFTING ARM.** Holding pH 5.0 required **adding
   3 M NaOH throughout the hour**, so the held arm accumulates sodium and dilution that the drifting
   arm does not. Neither the total NaOH added nor the final volume is reported. **This is a real
   confounder in an otherwise excellent design, and the paper does not address it.**

### 1.4 Analytics — verbatim
**GC:** *"fused silica capillary column (**25 m × 0·32 mm id**) coated with **SE 52/54 (1 µm film
thickness)** … Perkin Elmer Sigma 3B … **helium, 1·5 ml min⁻¹; temperature programme, 60 °C for
5 min followed by a gradient of 3 °C min⁻¹ to a final temperature of 200 °C, held for 20 min;
injector temperature, 225 °C; detector temperature, 250 °C; injection volume, 1 µl; injection
technique, Grob splitless.*** *"Retention times and quantitative data were obtained from a
**Hewlett Packard 3390A integrator** … **Individual components were quantified by the use of an
internal standard (tetradecane).**"*

**★ THE FURFURAL EXCEPTION, verbatim — this matters for the biggest number in the paper:**
> "**Quantitative data were obtained for all compounds using the concentrated isolates, apart from
> 2-furfural, for which data from the isolate obtained by heating without pH control were obtained
> BEFORE concentration.**"

**→ The pH 2.6 furfural value (3.53 × 10⁶) was measured on an UNCONCENTRATED isolate, on a different
aliquot from every other number in that column.** Sensible (it would have saturated the detector
otherwise) and honestly disclosed — but it means **the single largest number in the paper is not
strictly commensurable with the rest of its own column**, and the grand total for the pH 2.6 arm is
dominated by it.

**GC-MS:** *"Kratos MS 80 RFA … Carlo Erba MFC 500 GC … Kratos DS 90 data system. **electron impact;
electron energy, 70 eV; ion source temperature, 200 °C; separator temperature, 250 °C; ionisation
current, 100 µA; accelerating voltage, 4 kV; scan speed, 1 s per decade; mass range, 30−500 amu;
interscan time, 0·2 s; resolution, 1000.*** For the isolates prepared with control of the pH at 5·0,
GC-MS was also performed using **chemical ionisation (CI), with ammonia as the reagent gas (60 kPa,
0·6 bar); mass range, 50−500 amu; emission current, 500 µA** … runs were also performed using a
**fused silica column (50 m × 0·32 mm id) coated with CP-WAX 52 CB (0·2 µm film thickness)**."

**Precision, verbatim (Table 1 footnote *d*):**
> "**The yields quoted are the means obtained from two runs. The precision is ±16 % and was obtained
> from the means of the standard deviations of the values obtained for each run.**"

**★ A stated ±16 % on every cell.** That is a real, usable σ — better than Poisson 2019 (which
states none) and looser than Blank 1997 (≤10 %).

**Trace convention, verbatim (footnote):** *"**tr = trace, < 2 nmol mol⁻¹ xylose.**"*
**`—` = "compound not detected in this isolate". `NQ` = "component not quantified in this sample; it
eluted under the solvent peak."**

---

## 2. ★★ TABLE 1 — verbatim, both arms, all 58 compounds

> **TABLE 1. Volatile compounds identified from a xylose−lysine model system heated either without
> pH control (pH 2·6 system) or with control of the pH at 5**
>
> Columns: No. · Identity and MS reference^a · **LRI** Experimental^b (SE 52/54 · CP WAX 52 CB) ·
> **LRI** Reference (CWX 20M)^c · **Yield (nmol mol⁻¹ xylose)^d** for the **pH 2·6 system** and the
> **pH 5·0 system**.

*(Transcribed in full. `*` before a name = **tentatively identified** — footnote a: *"Compounds
whose mass spectra and LRI agree with those in the literature are classed as positively identified.
Other compounds are classed as tentatively identified (*)."* Superscript letters after a name are
the paper's own footnotes e/f/g and are reproduced as `^e` etc. where present. **Compound number 46
is absent from the printed table** — the numbering jumps 45 → 47.)*

| No. | identity | LRI SE 52/54 | LRI CP WAX 52 CB | LRI ref (CWX 20M) | **pH 2·6** | **pH 5·0** |
|---|---|---|---|---|---|---|
| | ***Aliphatic compounds*** | | | | | |
| 1 | Hydroxyacetone | < 800 | 1295 | 1272¹ | — | **NQ** |
| 2 | Butanedione | < 800 | 989 | 963 | **tr** | **NQ** |
| 3 | 2,3-Pentanedione | < 800 | 1056 | 1044 | — | **NQ** |
| 4 | Acetic acid | 820 | 1430 | 1425¹ | **tr** | **NQ** |
| 5 | *2-Oxopropylacetate | 887 | — | — | — | **27·6** |
| | **Total aliphatic compounds** | | | | **tr** | **27·6** |
| | ***Alicyclic compounds*** | | | | | |
| 6 | *2-Cyclohexen-1-ol | 943 | 1790 | — | **42·8** | tr |
| 7 | *Cyclopropylethylketone | < 800 | — | — | tr | — |
| | **Total alicyclic compounds** | | | | **42·8** | **tr** |
| | ***Benzene derivatives*** | | | | | |
| 8 | *2,5-Dihydroxyacetophenone^e | 1287 | 1702 | — | — | **98·6** |
| | ***Monocyclic furans — Type I furans*** | | | | | |
| 9 | 2(5H)-Furanone^e | 926 | 1737 | 1742² | **121** | **50·0** |
| | ***Type II furans*** | | | | | |
| 10 | 2-Methylfuran^f | < 800 | — | 838² / 866 | tr | NQ |
| 11 | 2-Furanmethanol | 866 | 1653 | 1650¹ / 1665³ | — | **492** |
| **12** | **2-Furfural^e** | **843** | **1459** | **1449²** | **★ 3·53 × 10⁶** | **★ 12·9 × 10³** |
| **13** | **4-Hydroxy-5-methyl-3(2H)furanone** *(= norfuraneol)* | **1076** | **2073** | — | **★ —** | **★ tr** |
| | **Total Type II furans** | | | | **3·53 × 10⁶** | **13·4 × 10³** |
| | ***Type III furans*** | | | | | |
| 14 | 3(or 4)-Methyl-2-furanmethanol | 892 | 1587 | 1587¹ | — | **9·0** |
| 15 | 5-Methyl-2-furanmethanol | 962 | 1710 | 1714¹ | tr | **7·2** |
| 16 | 5-Methyl-2-furfural^e | 964 | 1561 | 1562¹ / 1563 | **1410** | **8·2** |
| 17 | A dimethyl-2-furfural^e | 1058 | 1597 | 1594¹ / 1607³ | **80·6** | **286** |
| 18 | 2-Acetylfuran^e | 913 | 1495 | 1491 / 1494¹ | **865** | **111** |
| 19 | 2-Propionylfuran^e | 1020 | 1563 | 1563 / 1564¹ | **19·4** | **53·2** |
| 20 | 1-(2-Furyl)-2-propanone^e | 954 | 1511 | 1513³ / 1525¹ | **111** | **14·6** |
| 21 | 1-(2-Furyl)-3-butanone^e | 1087 | 1634 | 1636² / 1639³ | **350** | **350** |
| 22 | *cis-4-(2-Furyl)-3-buten-2-one^e | 1128 | 1722 | — | **372** | **8·0** |
| 23 | trans-4-(2-Furyl)-3-buten-2-one^e | 1197 | 1885 | 1895¹ | **194** | **440** |
| 24 | 1-(2-Furyl)-1,2-propanedione^e | 1075 | 1762 | 1765¹ / 1772³ | **143** | **28·0** |
| 25 | *1-(2-Furyl)-1,2-butanedione^e | 1182 | — | — | **52** | **14·0** |
| 26 | 2-Furfurylacetate | 1006 | 1528 | 1518 | — | **80·0** |
| 27 | A furan — 109(100), 124(43), 43(10), 53(8) | 974 | 1470 | — | — | tr |
| 28 | A furan — 109(100), 124(45), 53(28), 43(15), 110(15) | 1041 | 1479 | — | — | tr |
| 29 | A furan — 81(100), 109(88), 124(66), 53(17) | 1120 | 1534 | — | — | tr |
| 30 | A furan — 109(100), 110(15), 81(13) | 1084 | 1610 | — | tr | **200** |
| 31 | A furan — 81(100), 43(60), 109(50), 53(28), 109(8) | 1256 | — | — | — | tr |
| 32 | A furan — 81(100), 109(23), 123(22), 43(14), 53(8), 67(5) | 1175 | — | — | — | **4·0** |
| 33 | A furan — 123(100), 109(10), 67(8), 43(7) | 1323 | — | — | — | tr |
| | **Total Type III furans** | | | | **3550** | **1220** |
| | **Total monocyclic furans** | | | | **3·53 × 10⁶** | **14·7 × 10³** |
| | ***Dimeric furans*** | | | | | |
| 34 | 2,2′-Bifuran | 1059 | 1619 | — | **38·8** | — |
| 35 | 2,2′-Difurylmethane^e | 1107 | 1626 | 1622⁴ | **9·4** | **13·6** |
| 36 | 2,2′-Difurylethane^e | 1186 | 1669 | 1673⁴ | **13·6** | **37·0** |
| 37 | cis-2,2′-Difurylethylene^e | 1341 | 1990 | 2007⁴ | **226** | **111** |
| 38 | trans-2,2′-Difurylethylene^e | 1320 | 1968 | 1974⁴ | — | **12·6** |
| 39 | *5-(2-Furfuryl)-2-furfural | 1381 | 2148 | — | **2·2** | — |
| | **Total dimeric furans** | | | | **290** | **174** |
| | ***Benzofurans*** | | | | | |
| 40 | Benzofuran | 998 | 1493 | 1486³ | tr | **42·4** |
| 41 | *Benzodifuran | 1364 | #^(eluted too late) | — | — | **31·2** |
| | **Total benzofurans** | | | | **tr** | **73·6** |
| | ***Monocyclic pyrroles*** | | | | | |
| 42 | 2-Pyrrolealdehyde | 1035 | — | 1005 | — | **14·8** |
| 43 | 2-Acetylpyrrole | — | 1949 | 1935 | — | tr |
| 44 | A pyrrole^e — 137(100), 94(58), 109(32), 66(12) | 1098 | 1621 | 1629³ | — | **10·0** |
| 45 | A pyrrole — 151(100), 136(63), 108(42), 80(18), 43(18), 53(12) | 1274 | — | — | — | tr |
| | **Total monocyclic pyrroles** | | | | **—** | **24·8** |
| | ***1-(2-Furfuryl)pyrroles*** | | | | | |
| 47 | *1-(2-Furfuryl)pyrrole^e | 1189 | — | 1822⁵ | **13·6** | **19·0** |
| 48 | 1-(2-Furfuryl)-2-pyrrolealdehyde^e | 1425 | 2174 | 2234² | **91·4** | **286** |
| | **Total 1-(2-furfuryl)pyrroles** | | | | **105** | **305** |
| | ***Pyridines*** | | | | | |
| 49 | *An ethylmethylpyridine | 1140 | 1579 | — | — | **66·2** |
| 50 | *An ethyldimethylpyridine | 1264 | 1715 | — | — | tr |
| | **Total pyridines** | | | | **—** | **66·2** |
| | ***Pyrazines*** | | | | | |
| 51 | Pyrazine^e | 781 | 1208 | 1194⁶ | — | **★ 7910** |
| 52 | 2-Methylpyrazine | 824 | 1252 | 1251 / 1257⁶ | — | **1020** |
| 53 | 2,5(or 2,6)-Dimethylpyrazine^e | 910 | 1309 | 1306 / 1315⁶ | — | tr |
| 54 | A propenylpyrazine^e | 1067 | 1534 | 1546³ | **21·6** | **10·0** |
| 55 | An ethyl-2-vinylpyrazine | — | 1521 | 1528⁶ | — | tr |
| | **Total pyrazines** | | | | **21·6** | **★ 8940** |
| | ***2,3-Dihydro-1H-pyrrolizines*** | | | | | |
| 56 | *5-Formyl-2,3-dihydro-1H-pyrrolizine | 1402 | — | 1993⁷ | — | **3·0** |
| 57 | *5-Formyl-6-methyl-2,3-dihydro-1H-pyrrolizine^e | 1463 | — | 2157⁷ | — | **6·8** |
| 58 | *5-Acetyl-7-methyl-2,3-dihydro-1H-pyrrolizine^e | 1522 | 2144 | 2126⁷ | — | **112** |
| | **Total 2,3-dihydro-1H-pyrrolizines** | | | | **—** | **122** |
| | **Unknowns** | | | | **52·2** | **52·2** |
| | **★ GRAND TOTAL** | | | | **★ 3·53 × 10⁶** | **★ 2·47 × 10⁴** |

Table 1 footnote *b*, verbatim: *"LRI quoted were normally calculated using the retention time data
for components of the pH 5·0 system, using a standard solution of C₆−C₂₄ *n*-alkanes. **Values in
excess of 2000 were obtained by extrapolation** of the retention data."*
Footnote *e*: *"The identity is supported by CI MS data."* Footnote *f*: *"The identity was confirmed
by comparing the retention time of this compound with that of the authentic compound run under the
same conditions on the SE 52/54 column."* Literature LRI sources ¹−⁷: ¹ Jennings & Shibamoto 1980;
² Vernin et al. 1988; ³ Baltes & Bochmann 1987b; ⁴ Salter et al. 1988; ⁵ Farmer et al. 1989;
⁶ Baltes & Bochmann 1987c; ⁷ Tressl et al. 1985.

---

## 3. ARITHMETIC VERIFICATION — I checked the totals, and they close

Because the text layer mangled the exponents (§ preamble), I re-derived every group total from its
component rows against the visually-read values:

| group total | components summed | **computed** | **printed** | ✔ |
|---|---|---|---|---|
| Type II furans, pH 5·0 | 492 (#11) + 12 900 (#12) + tr (#13) | **13 392** | **13·4 × 10³** | ✔ |
| Monocyclic furans, pH 5·0 | 50·0 + 13 400 + 1220 | **14 670** | **14·7 × 10³** | ✔ |
| 2,3-Dihydro-1H-pyrrolizines, pH 5·0 | 3·0 + 6·8 + 112 | **121·8** | **122** | ✔ |
| Pyrazines, pH 5·0 | 7910 + 1020 + tr + 10·0 + tr | **8940** | **8940** | ✔ |
| Type III furans, pH 5·0 | 9·0+7·2+8·2+286+111+53·2+14·6+350+8·0+440+28·0+14·0+80·0+200+4·0 | **1613** | **1220** | ⚠️ *see below* |
| **Grand total, pH 5·0** | 27·6 + tr + 98·6 + 14 670 + 174 + 73·6 + 24·8 + 305 + 66·2 + 8940 + 122 + 52·2 | **24 554** | **2·47 × 10⁴** | ✔ |
| **Grand total, pH 2·6** | dominated by #12 (3·53 × 10⁶); all other groups sum to ≈ 4060 | **3·53 × 10⁶** | **3·53 × 10⁶** | ✔ |

**★ Every total closes to within rounding except the Type III furan subtotal at pH 5·0, where my
component sum (1613) exceeds the printed subtotal (1220) by 393.** Since the *monocyclic furans*
total and the *grand total* both close using the **printed 1220**, the discrepancy is internal to
the Type III subtotal, not to the rest of the table. **Most likely explanation: one or two of the
"A furan" rows (#30 at 200, #17 at 286) were excluded from the subtotal, or my reading of one of
them is off by a factor.** Recorded as an unresolved minor inconsistency; **it does not affect the
furfural numbers, the grand totals, or anything used below.**

**★ The g/kg cross-check confirms the whole scale.** MW(2-furfural) = 96.08.
pH 5·0: 12 900 nmol × 96.08 g/mol = **1.239 mg furfural per mol xylose**. For that to be
**522 g kg⁻¹** of the volatiles (abstract), the total volatile mass must be **2.374 mg per mol
xylose**, i.e. an average MW of 2.374 mg / 24 554 nmol = **96.7 g/mol** across the whole isolate.
**That is almost exactly furfural's own MW and is entirely plausible for a furan-dominated
profile.** **The abstract's mass shares and Table 1's molar yields are mutually consistent, which
independently validates both readings.**

---

## 4. THE PH-TRAJECTORY RESULT, IN NUMBERS

### 4.1 ★★ The headline contrast

| quantity | **pH held at 5·0** | **pH uncontrolled (4·9 → 2·6)** | **fold** |
|---|---|---|---|
| **Total volatiles** | **2·47 × 10⁴ nmol/mol xylose** | **3·53 × 10⁶** | **★ 143× UP when pH drifts** |
| **2-Furfural** | **1·29 × 10⁴** | **3·53 × 10⁶** | **★ 274× UP** |
| **2-Furfural, as mol % on xylose** *(derived)* | **0·00129 mol %** | **0·353 mol %** | 274× |
| **Furfural share of the volatiles** *(abstract)* | **522 g kg⁻¹** | **999 g kg⁻¹** | 1·9× |
| **Compounds identified** | **58** | **28** | 2·1× more at pH 5 |
| **Total N-containing volatiles** | **≈ 9·46 µmol/mol** *(derived, §4.2)* | **0·127 µmol/mol** | **★ 75× DOWN when pH drifts** |
| **Total pyrazines** | **8940 nmol/mol** | **21·6** | **★ 414× DOWN** |
| **Pyrazine (parent) alone** | **7910 nmol/mol** | **not detected** | **★ ∞ (on/off)** |
| **Monocyclic pyrroles** | 24·8 | **not detected** | on/off |
| **Pyridines** | 66·2 | **not detected** | on/off |
| **2,3-Dihydro-1H-pyrrolizines** | 122 | **not detected** | on/off |
| **★ Norfuraneol (4-hydroxy-5-methyl-3(2H)furanone)** | **tr (< 2 nmol/mol)** | **not detected** | **on/off, at the floor** |
| **2-Furanmethanol** | 492 | **not detected** | on/off |
| Type III furans (total) | 1220 | 3550 | 2·9× up when pH drifts |
| Dimeric furans | 174 | 290 | 1·7× up |
| Benzofurans | 73·6 | tr | ≥37× down |

### 4.2 ⚠️ A UNIT ERROR IN THE PAPER'S OWN TEXT — and the arithmetic that exposes it

The Results (p. 478) state, verbatim:
> "**The total yield of the 16 nitrogen-containing compounds from the system held at pH 5·0 was
> 10·1 mmol mol⁻¹ xylose compared with a total yield of 0·127 mmol mol⁻¹ xylose for the three
> compounds from the pH 2·6 system.**"

**Summing Table 1's own N-containing rows** (compounds 42−45, 47−50, 51−55, 56−58 = **16 compounds**,
exactly as the text says):
`14·8 + 10·0 + 19·0 + 286 + 66·2 + 7910 + 1020 + 10·0 + 3·0 + 6·8 + 112` = **9458 nmol/mol =
9·46 µmol/mol = 0·00946 mmol/mol.**
pH 2·6: `13·6 + 91·4 + 21·6` = **126·6 nmol/mol = 0·127 µmol/mol.**

**★ The pH 2·6 figure matches the text EXACTLY at 0·127 — but in µmol, not mmol.** The text's unit
is wrong by a factor of 1000 in both figures; it should read **µmol mol⁻¹ xylose**. (The pH 5·0
figure, 10·1 vs my 9·46, differs by 6·8 % — **inside the paper's own ±16 % precision**, and probably
an inclusion/exclusion of `tr` cells.)

**→ Use `9·46 µmol/mol` and `0·127 µmol/mol`, or better, sum Table 1 yourself. Do NOT ingest
"10·1 mmol mol⁻¹" — it is 1000× too large and would put the N-compounds at 1 mol % of the xylose.**

### 4.3 The authors' own explanation, verbatim (p. 478 and p. 483)
> "**The main impact of heating without pH control was an increase in the total yield of reaction
> products by a factor of ~100, due mainly to an increased yield of 2-furfural. Degradation of the
> ARP via 1,2-enolisation is favoured in the pH 2·6 system, and 2-furfural is the main compound
> formed from a pentose ARP by this route.** The other principal effect of pH was that in the pH 2·6
> system both total yields and numbers of some nitrogen-containing reaction products were lower
> while the monocyclic pyrroles, pyridines and 2,3-dihydro-1H-pyrrolizines were not detectable."

> "Heating without pH control increased the total yield of volatile reaction products by a factor of
> about 100, due mainly to an increased yield of 2-furfural. **Both total yields and numbers of
> nitrogen-containing reaction products were greater when the pH was held at 5·0 and pyrazines were
> the second most abundant class of volatile compound in that system, accounting for 360 g kg⁻¹ of
> the total isolate.**"

**★ THE 1,2- vs 2,3-ENOLISATION FORK IS THE MECHANISM THE AUTHORS NAME, AND IT IS THE SAME FORK
BLANK 1997 §5 INVOKES (via Beck 1988's 1-DG:3-DG = 20:1 at pH 7 / 8:5 at pH 4.5).** Two papers on
two different products, ten years apart, attributing their pH response to the same partition, in the
same direction: **acid favours 1,2-enolisation → 3-deoxyosone → furfural; neutral/alkaline favours
2,3-enolisation → 1-deoxyosone → furanones and N-heterocycles.** ★ This is a genuine cross-paper
mechanistic convergence and it is the most transferable thing in this paper. See
`k5b_dmhf_synthesis.md` §4.

### 4.4 The Introduction's framing, verbatim (p. 477) — includes the paper's own DMHF citation
> "pH has an important effect on the Maillard reaction. The rate of browning … **increases with pH
> over the pH range 4−8** (Kawashima et al. 1980; Westphal et al. 1988). **Degradation of the Amadori
> rearrangement product (ARP) via 2,3-enolisation is favoured by increased pH values, and this route
> possesses a greater browning potential than 1,2-enolisation** (Nursten 1986). As far as the
> volatile reaction products are concerned, some are produced only over a certain pH range … For
> example, **a wider range of heterocyclic reaction products was formed in a model system of
> 2,5-dimethyl-4-hydroxy-3(2H)-furanone and cysteine at pH 7·1 rather than at pH 2·2 or 5·1 and, in
> particular, pyrazines were only detected at pH 7·1 (Shu and Ho 1988).** In a study of a
> glucose−proline model system **Tressl et al. (1989) showed that large amounts of
> 5-hydroxymethyl-2-furfural were formed at pH 3, but that this compound could not be detected at
> pH 5 or above.** Yields of selected nitrogen-containing components formed in a range of
> sugar−amino acid model systems increased with pH and some could not be detected when the pH was
> below 7 (Tressl et al. 1989)."

**★ Note two things.** (i) **The paper's only DMHF mention is a citation of Shu & Ho 1988** — the
other paper dossiered in this wave (`shu1988_extraction.md`), and it is cited for exactly the result
that dossier transcribes. (ii) **The Tressl 1989 5-HMF claim — "large amounts at pH 3, undetectable
at pH 5 or above" — is directly relevant to the HMF channel** (`research_round3_channels.md` §A) and
is **[C], not read.**

---

## 5. ★★ THE NORFURANEOL RESULT — and why it does NOT contradict `norfuraneol ≫ DMHF`

### 5.1 What the paper says, verbatim (p. 481)
> "**4-Hydroxy-5-methyl-3(2H)-furanone (HMFone) was also detected only at the higher final pH**, and
> is of interest since it is a precursor of several coloured Maillard reaction products (Ames 1992).
> **It is formed by the degradation of a pentose ARP via 2,3-enolisation, a pathway favoured by
> increased pH values. The identification of HMFone in only trace amounts may be due to its
> participation in subsequent reactions, including those resulting in the formation of coloured
> compounds.** This suggestion is consistent with the results of **Mottram and Leseigneur (1990) who
> identified HMFone in a ribose−cysteine system in decreasing amounts when the pH was increased over
> the range 4·5−6·5**, probably due to the formation of increased amounts of coloured compounds."

### 5.2 The apparent conflict, and its resolution

`blank1996_extraction.md` item #9 and `blank1997_extraction.md` item #32 both record that
**norfuraneol is the MAJOR volatile of a pentose + amino-acid system.** Here, in a pentose + amino
acid system, **norfuraneol is at the detection floor (`tr`, < 2 nmol/mol) in one arm and absent in
the other**, while **furfural is 6450× larger** in the same isolate.

**This is NOT a contradiction of the ordering, for four independent reasons, and all four should be
carried:**

1. **DMHF was not detected either.** The ordering `norfuraneol ≫ DMHF` is untestable here because
   **both** terms are at or below the floor. `tr > not-detected` is, if anything, weakly consistent
   with it. **The paper is SILENT on the ordering, not contradictory.**
2. **The amine is LYSINE, not glycine or alanine.** Lysine has a side-chain ε-amino group, forms
   pyrrolizines and pyrroles (compounds 56−58, 42−45 — the paper's headline novelty), and provides
   **no Strecker C1 or C2 unit of the kind Blank's furanone route needs.** A different amine gives a
   different product spectrum; that is the whole point of the amino-acid switch in Blank 1997 §2.1.
3. **★ SDE destroys and diverts furanones** (§1.3 item 1). An hour of continuous steam distillation
   at reflux into ether is close to the worst possible isolation for a water-miscible, thermally
   labile 3(2H)-furanone. Blank's ether perforation was at **pH 4 and 4 °C with an isotope-labelled
   internal standard**; this is at reflux with an alkane standard. **The two methods are not
   comparable for this compound class.**
4. **The authors themselves attribute the trace level to CONSUMPTION**, not to non-formation, and
   support it with an independent citation (Mottram & Leseigneur 1990, norfuraneol *decreasing* with
   pH in ribose−cysteine). **Norfuraneol is being made and then eaten by the browning chemistry.**
   ★ That is a **sink** statement, and it is the same class of statement `shu1988_extraction.md`
   makes for DMHF.

**→ FOR THE REGISTER: this paper must NOT be scored against the `norfuraneol ≫ DMHF` ordering. Its
norfuraneol cell is a method-limited near-null in a different amine system, and the authors read it
as consumption.**

---

## 6. THE NAMING TRAP — what "HMF"/"HDMF"/"DMHF" denotes IN THIS PAPER

**⚠️ THIS PAPER CARRIES THE TRAP, IN ITS MOST DANGEROUS FORM.** There is **no abbreviation table**
anywhere; abbreviations are introduced inline, and the key one is introduced **on page 481, four
pages after the compound first appears in Table 1**.

| token | occurrences | denotes, in this paper | defined where |
|---|---|---|---|
| **★★ `HMFone`** | 3 (all on p. 481) | **4-HYDROXY-5-METHYL-3(2H)-FURANONE = NORFURANEOL.** Verbatim: *"4-Hydroxy-5-methyl-3(2H)-furanone (**HMFone**) was also detected only at the higher final pH."* | p. 481, Results, "Monocyclic furans" |
| **`HMF` standing alone** | **ZERO** | — | — |
| **`HDMF`** | **ZERO** | — | — |
| **`DMHF`** | **ZERO** | — | — |
| **2,5-dimethyl-4-hydroxy-3(2H)-furanone** (furaneol, spelled out) | **2** | Furaneol — but **only in the Introduction (p. 477) citing Shu & Ho 1988** and in the reference list. **It is NOT a compound in this study.** | Introduction |
| `ARP` | many | Amadori Rearrangement Product | p. 477 |
| **5-hydroxymethyl-2-furfural** (the *other* HMF) | **1** | **5-HMF** — in the Introduction, citing Tressl et al. 1989 (pH 3 vs pH 5 in glucose−proline). **Also not a compound in this study.** | p. 477 |

**★★ THE TRAP: `HMFone` looks like an HMF abbreviation and is not one.** A reader or a harvester
scanning for HMF-family strings in this paper will find `HMFone` (= **norfuraneol**, a C5 furanone),
`2,5-dimethyl-4-hydroxy-3(2H)-furanone` (= **furaneol**, a C6 furanone, cited not measured) and
`5-hydroxymethyl-2-furfural` (= **5-HMF**, cited not measured) — **three different molecules, in one
paper, none of which is what the others are.** Combined with Blank & Fay 1996's `HMF (3)` =
norfuraneol (`blank1996_extraction.md` §7), **this makes TWO papers in the corpus where an
"HMF"-shaped token means norfuraneol.**

**Cluster-wide summary of what the token means, per paper, as the brief asked:**

| paper | `HMF` / `HMF`-shaped token | means |
|---|---|---|
| `blank1996` | **`HMF (3)`** | **NORFURANEOL** |
| `blank1997` | *token never used* | — (norfuraneol always spelled out) |
| **`apriyantono1993`** | **`HMFone`** | **NORFURANEOL** |
| `wang2008` | *body: never used*; **reference title: `5-(hydroxymethyl)furfural`** | **5-HMF** (a cited paper's title only) |
| `poisson2019` | *never used* | — |
| `shu1988` | *never used* | — |

---

## 7. KINETICS — **NONE, and the temperature is never even stated**

No rate constant, no activation energy, no time course, no half-life. **One heating time (1 h) and
one heating mode ("refluxed") for both arms.** The paper **never states a temperature anywhere** —
"refluxed" in an aqueous Likens−Nickerson apparatus implies ≈100 °C at 1 atm, but that is an
inference, not a datum, and the boiling point of a 1 M xylose + 1 M lysine·HCl solution is elevated
above 100 °C by an unstated amount. **Nothing here can be made into a kinetic parameter, and any
"100 °C" attached to this dataset must be flagged as inferred.**

---

## 8. CROSS-CHECKS AGAINST THE CLUSTER

| constraint | this paper's evidence | verdict |
|---|---|---|
| **`blank1996` #9 / `blank1997` #32: norfuraneol ≫ DMHF** | norfuraneol `tr`; **DMHF not detected at all** | **★ SILENT — must not be scored.** §5.2. Both terms at/below floor; different amine; SDE isolation; authors read the trace as consumption. |
| **`blank1997` §5.2 (via Beck 1988): 2,3-enolisation is favoured at higher pH; 1,2-enolisation at lower** | **§4.3 — the authors' explicit explanation for their entire result**, and the 274× furfural rise on acidification is exactly what a 1,2-enolisation shift predicts | **★★ CORROBORATED, on a different product, by a different lab, ten years apart.** The strongest mechanistic convergence in this wave. |
| **`blank1997` Table 4: HDMF/EHMF FAVOURED at pH 7 over pH 5** | This paper's furanone (norfuraneol) forms **only in the held-pH-5·0 arm** and not in the acid arm — same directional sign | **★ CORROBORATED directionally**, though at the detection floor. |
| **`wang2008`: MG + glycine DMHF RISES as pH FALLS** | This paper: everything *except* furfural falls as pH falls | **NO CONTACT.** Different precursor (MG vs pentose), different amine (glycine vs lysine), and this paper has no DMHF. |
| **`shu1988`: pyrazines from DMHF + cysteine only at pH 7·1** | This paper: **pyrazines 8940 nmol/mol at pH 5·0 vs 21·6 at pH 2·6, and pyrazine itself entirely absent in the acid arm** | **★ CORROBORATED as a direction** — pyrazine formation is strongly base-favoured in both systems. **This paper cites Shu & Ho 1988 for precisely this** (§4.4). Two systems, one direction. |
| **`research_round3_channels.md` §C.2: "RATIO-ONLY"** | Table 1 is in `nmol mol⁻¹ xylose` | **★★ CONTRADICTED. See §0a. The register entry needs correcting.** |
| **`research_round3_channels.md` §C.1: furfural comes from 3-DG, not HMF** | The authors state exactly that: *"2-furfural is the main compound formed from a pentose ARP by [the 1,2-enolisation] route"* | **★ CORROBORATED**, and here for a **pentose** — which is the case round-3 §C.1 flagged as *missing*, since the Lee et al. cake papers are hexose systems. **This is the pentose furfural source round-3 said did not exist.** |

**★★ THE FURFURAL FINDING IS THE UNDER-VALUED RESULT HERE.**
`research_round3_channels.md` §C.1 records, verbatim: *"⚠️ Note this is a **hexose** system — furfural
from glucose is a C5 fragment route, not pentose dehydration. If the repo's furfural feed is
pentose-based, these numbers are a *different* production route and are not commensurable."*
**This paper supplies the missing pentose case, with absolute molar yields, at two pH states:**
**0·353 mol % on xylose (drifting to pH 2·6) and 0·00129 mol % (held at pH 5·0).** That is a
**274× pH lever on furfural from a pentose**, and it is a level, not a ratio.

---

## 9. CITED, NOT MEASURED

| quantity / claim | source, as this paper cites it | page |
|---|---|---|
| **A wider range of heterocycles from DMHF + cysteine at pH 7·1 than at 2·2 or 5·1; pyrazines only at pH 7·1** | **Shu & Ho 1988**, JAFC **36**, 801−803 — *`data/articles/shu1988.pdf`, dossiered in this wave* | 477 |
| **Large amounts of 5-hydroxymethyl-2-furfural at pH 3; undetectable at pH 5 or above** (glucose−proline) | **Tressl, Helak, Martin & Kersten 1989**, in *Thermal Generation of Aromas*, ACS, pp. 156−171 | 477 |
| Browning rate increases with pH over pH 4−8 | **Kawashima, Itoh & Chibata 1980**, *Agric. Biol. Chem.* **44**, 1595−1599; **Westphal, Kroh & Follmer 1988**, *Nahrung* **32**, 117−120 | 477 |
| **2,3-enolisation has a greater browning potential than 1,2-enolisation** | **Nursten 1986**, in *Concentration and Drying of Foods*, pp. 53−68 | 477 |
| **Norfuraneol DECREASES as pH rises over 4·5−6·5 in a ribose−cysteine system** | **Mottram & Leseigneur 1990**, in *Flavour Science and Technology*, pp. 121−124 | 481 |
| Norfuraneol is a precursor of coloured Maillard products | **Ames 1992**, in *Progress in Food Proteins − Biochemistry*, pp. 99−153 | 481 |
| Sugar-fragmentation reactions are favoured by increased pH | **Hodge 1953**, *JAFC* **1**, 928−943 | 481 |
| cis-2,2′-difurylethylene forms by **self-condensation of 2-furfural**; 2-acetylfuran from furfural + H₂S + NH₃ | **Shibamoto 1977**, JAFC **25**, 206−208 | 481−482 |
| 2-acetylfuran and 5-methyl-2-furfural from a 2-furfural−cysteine−methionine mixture | **Silwar & Tressl 1989**, *Z. Lebensm. Unters. Forsch.* **189**, 205−211 | 481 |
| Pyrazine 632 g kg⁻¹ and 2-methylpyrazine 368 g kg⁻¹ of total pyrazines in a **glucose−lysine system at 95 °C, pH 5, 2 h** | **Leahy & Reineccius 1989**, in *Thermal Generation of Aromas*, pp. 196−208 | 483 |
| Hydroxyacetone is a precursor of 2,5- and 2,6-dimethylpyrazine | **Rizzi 1988**, JAFC **36**, 349−352 | 483 |
| 2,3-dihydro-1H-pyrrolizines are proline-specific Maillard products; the mechanism adapted for lysine | **Tressl, Rewicki, Helak, Kamperschroer & Martin 1985**, JAFC **33**, 919−923 | 482 |
| Pyrrolidine from a lysine ARP by intramolecular cyclisation | **Yaylayan & Sporns 1989**, *Food Chem.* **33**, 81−91; **Yaylayan, Paré, Laing & Sporns 1990** | 482 |
| 2-acetylfuran and 1-(2-furyl)-2-propanone from **1 M aqueous xylose refluxed at pH 7 for 1 h** | **Hincelin, Ames, Apriyantono & Elmore 1992**, *Food Chem.* **44**, 381−389 | 481 |

**★ The Leahy & Reineccius 1989 comparison is worth flagging**: the authors compare their pyrazine
**885 / 114 g kg⁻¹** (pyrazine / 2-methylpyrazine, as a share of total pyrazines at pH 5·0) against
**632 / 368 g kg⁻¹** for a **glucose**−lysine system at 95 °C, pH 5, 2 h. **A pentose-vs-hexose
pyrazine-distribution contrast at nearly matched conditions**, and it is the only near-matched
hexose comparator in the paper. **[C], not read.**

---

## 10. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured here · **[C]** cited · **[D]** derived by me. Conditions for all **[M]**
rows: **1·0 M xylose + 1·0 M lysine·HCl in 500 mL water, refluxed 1 h (temperature never stated;
≈100 °C inferred), modified Likens−Nickerson SDE into 15 mL diethyl ether, internal standard
*n*-tetradecane, means of two runs, precision ±16 %, isolates prepared in at least triplicate.**

| # | quantity | value | class | verdict |
|---|---|---|---|---|
| 1 | **★★ Total volatiles, pH held 5·0** | **2·47 × 10⁴ nmol mol⁻¹ xylose** | **[M]** | ★★ USE. |
| 2 | **★★ Total volatiles, pH drifting 4·9 → 2·6** | **3·53 × 10⁶ nmol mol⁻¹ xylose** | **[M]** | ★★ USE. |
| 3 | **★★ Total-volatile pH-trajectory ratio** | **143× (drifting / held)** | **[D]** | ★★ USE. |
| 4 | **★★ 2-Furfural, pH held 5·0** | **12·9 × 10³ nmol/mol = 0·00129 mol % on xylose** | **[M]** | ★★ USE. **A pentose furfural LEVEL, which round-3 §C.1 says did not exist.** |
| 5 | **★★ 2-Furfural, pH drifting to 2·6** | **3·53 × 10⁶ nmol/mol = 0·353 mol % on xylose** | **[M]** | ★★ USE. ⚠️ measured on an **unconcentrated** isolate (§1.4). |
| 6 | **★★ Furfural pH-trajectory ratio** | **274×** | **[D]** | ★★ USE. **NOT the ~2× the abstract's g/kg suggests.** |
| 7 | **Furfural as a share of the volatiles** | **522 g kg⁻¹ (pH 5·0) / 999 g kg⁻¹ (pH 2·6)** | **[M]** | USE only as a share; §0a. |
| 8 | **★ Norfuraneol** | **`tr` (< 2 nmol/mol) at pH 5·0; NOT DETECTED at pH 2·6** | **[M]**, at floor | ★ USE **only as an on/off**. **Do NOT score against `norfuraneol ≫ DMHF`.** §5.2. |
| 9 | **★ DMHF / HDMF** | **NOT DETECTED, NOT SOUGHT, NOT MENTIONED as a product** | — | **ABSENT. §0b.** |
| 10 | **Total pyrazines** | **8940 (pH 5·0) vs 21·6 (pH 2·6) nmol/mol** — **414×** | **[M]** | ★★ USE. Strongest base-favoured class. |
| 11 | **Pyrazine (parent) alone** | **7910 (pH 5·0) vs not detected (pH 2·6)** | **[M]** | ★★ USE as an on/off. |
| 12 | 2-Methylpyrazine | **1020 (pH 5·0) vs not detected** | **[M]** | ★ USE. |
| 13 | Pyrazine : 2-methylpyrazine share of total pyrazines, pH 5·0 | **885 : 114 g kg⁻¹** | **[M]** | USE as a ratio. Compare **632 : 368** for glucose−lysine (**[C]** Leahy 1989). |
| 14 | **Total N-containing volatiles (16 compounds)** | **9·46 µmol/mol (pH 5·0) vs 0·127 µmol/mol (pH 2·6)** — **75×** | **[D]** from Table 1 | ★★ USE. ⚠️ **The paper prints these as `mmol`; that is a 1000× unit error. §4.2.** |
| 15 | **Monocyclic pyrroles / pyridines / pyrrolizines** | **24·8 / 66·2 / 122 (pH 5·0); ALL NOT DETECTED (pH 2·6)** | **[M]** | ★★ USE — **three clean on/off classes.** |
| 16 | 2-Furanmethanol | **492 (pH 5·0) vs not detected** | **[M]** | ★ USE as on/off. |
| 17 | 5-Methyl-2-furfural | **8·2 (pH 5·0) vs 1410 (pH 2·6)** — **172× up on acidification** | **[M]** | ★ USE. Tracks furfural, as expected for a 1,2-enolisation product. |
| 18 | 2-Acetylfuran | **111 (pH 5·0) vs 865 (pH 2·6)** | **[M]** | USE. |
| 19 | Type III furans (total) | **1220 (pH 5·0) vs 3550 (pH 2·6)** | **[M]** | USE. ⚠️ subtotal does not close on my component sum (§3). |
| 20 | Dimeric furans (total) | **174 (pH 5·0) vs 290 (pH 2·6)** | **[M]** | USE. |
| 21 | Benzofurans | **73·6 (pH 5·0) vs tr (pH 2·6)** | **[M]** | USE. **"the current study is the first reporting benzodifuran in a sugar−amino acid model system."** |
| 22 | Compound counts | **58 identified (pH 5·0) vs 28 (pH 2·6)**; plus 14 unidentified, **2 g kg⁻¹ (pH 5·0) and a trace (pH 2·6)** | **[M]** | USE. |
| 23 | **★★ Precision** | **±16 %**, means of two runs, isolates in ≥ triplicate | **[M]** | ★★ USE as σ on every cell. |
| 24 | Trace threshold | **`tr` = < 2 nmol mol⁻¹ xylose** | **[M]** | ★★ **CRITICAL** — defines the floor for item #8. |
| 25 | **★★ The pH-trajectory design** | **held at 5·0 by 3 M NaOH addition throughout vs uncontrolled 4·9 → 2·6**, otherwise identical | **[M]** | ★★★ USE — **the reason this paper matters.** ⚠️ §1.3 item 3: the held arm accumulates unreported NaOH. |
| 26 | Loading | **1·0 M xylose : 1·0 M lysine·HCl**, 500 mL | **[M]**/**[D]** | Context. No dose−response. |
| 27 | **Temperature** | **NEVER STATED.** "Refluxed", 1 h. ≈100 °C **inferred**. | — | **⚠️ ABSENT. §7.** |
| 28 | Isolation | **Likens−Nickerson SDE, 1 h, into 15 mL Et₂O**, concentrated to 0·05 mL | **[M]** | ★ **CRITICAL CAVEAT for furanones. §1.3.** |
| 29 | Quantification basis | **one internal standard (*n*-tetradecane) for all 58 compounds; no response factors** | **[M]** | **⚠️ Weight accordingly.** |
| 30 | rate constant, Ea, time course, temperature series | **NONE** | — | **ABSENT. §7.** |
| 31 | **1,2- vs 2,3-enolisation as the pH mechanism** | qualitative, authors' own explanation | **[M]** interpretation / **[C]** Nursten 1986 | ★★ USE as topology. **Converges with Blank 1997 §5.2 / Beck 1988.** |
| 32 | Norfuraneol is **consumed** into coloured products | qualitative | **[M]** interpretation, supported by **[C]** Mottram & Leseigneur 1990 | ★ USE as a **sink** statement. |
| 33 | 5-HMF: large at pH 3, undetectable at pH ≥ 5 (glucose−proline) | qualitative | **[C]** Tressl 1989 | **Secondary. Not read. Relevant to the HMF channel.** |
| 34 | Pyrazine/2-methylpyrazine 632/368 g kg⁻¹, glucose−lysine, 95 °C, pH 5, 2 h | ratio | **[C]** Leahy & Reineccius 1989 | **Secondary. A near-matched HEXOSE comparator.** |

---

## 11. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION.**

### 11.1 ★★★ THE NAMED ROLE THIS PAPER DESERVES

**Proposed name: `apriyantono1993_xylose_lysine_pH_trajectory_pair`.**
**Proposed role: `HOLD-OUT — pH-TRAJECTORY VALIDATION PAIR`. Scored as ONE paired test, never as
two independent rows.**

**Why it needs its own named role rather than a slot in a general pH bucket:**

1. **It is the only held-vs-drifting pair in the corpus.** `k3` §C.11 records the absence of exactly
   this control for Zhou 2023. Four of the five papers in this DMHF wave report pH as an **unheld
   set-point** with no final pH; **this one runs both states, in one lab, with one protocol, from
   one initial pH.**
2. **It tests a different thing from a pH ladder.** Blank 1997 Table 4 (pH 5/6/7, all buffered)
   tests *"does the model get the pH-5 rate right, and the pH-7 rate right?"*. This pair tests
   *"does the model's internal pH STATE evolve correctly, and does the chemistry follow it?"* **A
   model can pass every point of a pH ladder and still fail this**, if it holds pH as a fixed input
   rather than tracking acid generation.
3. **It has four independent scoreable channels in one experiment**, so it is hard to pass by luck:
   **(i)** total volatiles **143×** up on drift; **(ii)** furfural specifically **274×** up;
   **(iii)** N-heterocycles **75×** down, with **four whole compound classes going from present to
   not-detected**; **(iv)** the only furanone present goes to zero.
4. **Its σ is stated (±16 %)**, so a likelihood can be written honestly — unlike Poisson 2019.

**Suggested scoring shape (proposal only):** score the **log-ratio between the arms**, not the two
levels, for each of the four channels. That form is immune to the internal-standard/response-factor
weakness (§ item 29) and to the SDE recovery bias, because both arms went through the identical
isolation. **Passing requires getting the DIRECTION and the ORDER OF MAGNITUDE of all four right
simultaneously.**

**⚠️ The one caveat that must ride with the role:** the held arm received **unreported 3 M NaOH
throughout the hour** (§1.3 item 3), so it is not sodium- or volume-matched to the drifting arm.
**Any model that represents ionic strength must be told the held arm's Na⁺ is unknown.**

### 11.2 Other proposed **HOLD-OUT** candidates

| candidate | proposed role | why |
|---|---|---|
| **★ Item #4/#5/#6 — pentose furfural at two pH states, in mol %** | **HOLD-OUT, LEVEL.** | Round-3 §C.1 declares a gap: no pentose furfural yields exist. **This closes it, with an absolute basis and a stated σ.** ⚠️ carry the "unconcentrated isolate" caveat on the pH 2·6 cell. |
| **★ Item #15 — four compound classes present at pH 5·0 and NOT DETECTED at pH 2·6** | **HOLD-OUT, structural on/off, ×4.** | Free of any calibration; a model that emits pyrazines in the acid arm fails. |
| **Item #10/#11 — pyrazines 414× down on acidification; pyrazine itself on/off** | **HOLD-OUT, directional.** | Corroborates Shu & Ho 1988's pyrazines-only-at-pH-7·1 in a completely different system. **Score the two together.** |

### 11.3 Proposed **FIT**

| candidate | recommendation |
|---|---|
| **Nothing from this paper as a magnitude fit target** | **RECOMMEND AGAINST.** Three transmission defects stack: **(a)** one alkane internal standard for 58 diverse compounds with no response factors; **(b)** Likens−Nickerson SDE, which for the furanone class is actively destructive (§1.3); **(c)** **the temperature is never stated** (§7), so the condition cannot be reproduced exactly. The **ratios** survive all three because both arms share them; the **levels** do not. |
| **The four log-ratios of §11.1** | **Fit-eligible in principle** if a pH-state term ever needs calibrating — but **strongly preferred as a hold-out**, because it is the only test of that term the corpus has and burning it in a fit set would leave the pH-state model unvalidated. |
| **Any barrier or Ea** | **PROHIBITED.** One time, one (unstated) temperature. |

### 11.4 Required register corrections

1. **`research_round3_channels.md` §C.2 and its §E download-list row must be corrected from
   "RATIO-ONLY" to "ABSOLUTE MOLAR YIELDS, `nmol mol⁻¹ xylose`, ±16 %, plus derived mass shares."**
   The `[ABS]` verdict was correct for the abstract and wrong for the paper.
2. **`research_round3_channels.md` §C.1's declared gap — "no pentose furfural yields" — is closed by
   items #4/#5.**
3. **The paper's own "10·1 mmol mol⁻¹" must never be ingested** (§4.2); it is 1000× too large.

---

## 12. DECLARED GAPS — verbatim, for `k3` §C

> **"There is no DMHF, HDMF or furaneol IN this paper.** 2,5-Dimethyl-4-hydroxy-3(2H)-furanone
> appears twice, both times in the Introduction citing Shu & Ho 1988, and is not among the 58
> identified or 14 unidentified volatiles. The only furanone measured is 4-hydroxy-5-methyl-3(2H)-
> furanone (norfuraneol, which the paper abbreviates **HMFone**), at `tr` (< 2 nmol mol⁻¹ xylose) in
> the held-pH arm and not detected in the drifting arm."

> **"The temperature is never stated anywhere in the paper.** Both arms were 'refluxed for 1 h' in an
> aqueous Likens−Nickerson apparatus. Any temperature attached to this dataset — including 100 °C —
> is inferred, and the boiling-point elevation of a 1 M xylose + 1 M lysine·HCl solution is not
> accounted for."

> **"⚠️ THE ROUND-3 'RATIO-ONLY' VERDICT ON THIS PAPER IS WRONG.** `research_round3_channels.md`
> §C.2 was written from the abstract alone and concluded that the 522/999 g kg⁻¹ figures made this a
> ratio source. **Table 1 reports every yield in `nmol mol⁻¹ xylose`**, an absolute molar yield on
> the limiting sugar, quantified against an internal standard with a stated ±16 % precision. The
> furfural pH effect is **274×**, not the ≈2× the g kg⁻¹ shares imply, because the denominator
> collapsed 143× at the same time."

> **"The paper's own text carries a 1000× unit error.** p. 478: 'The total yield of the 16
> nitrogen-containing compounds from the system held at pH 5·0 was 10·1 **mmol** mol⁻¹ xylose
> compared with a total yield of 0·127 **mmol** mol⁻¹ xylose'. Summing Table 1's own rows gives
> **9·46 µmol** and **0·127 µmol** respectively — the second matching the text's digits exactly.
> The unit should be µmol. Do not ingest the printed mmol values."

> **"The pH 2·6 furfural value was measured on a DIFFERENT ALIQUOT from every other number in its
> column.** Experimental, verbatim: 'Quantitative data were obtained for all compounds using the
> concentrated isolates, apart from 2-furfural, for which data from the isolate obtained by heating
> without pH control were obtained before concentration.' The largest number in the paper —
> 3·53 × 10⁶ — and therefore the pH 2·6 grand total, is not strictly commensurable with the rest of
> its own column."

> **"The held-pH arm is not sodium- or volume-matched to the drifting arm.** Holding pH 5·0 required
> adding 3 M NaOH before and during heating; neither the total NaOH added nor the final volume is
> reported. This is an unaddressed confounder in an otherwise unique design."

> **"Isolation was by Likens−Nickerson simultaneous distillation−extraction, an hour of continuous
> steam stripping at reflux into diethyl ether.** For water-miscible, thermally labile 3(2H)-
> furanones this is close to a worst-case isolation, and it is the primary reason this paper's
> norfuraneol near-null cannot be compared with the cold, pH-4, isotope-corrected ether perforation
> of Blank et al. 1997. **This paper must NOT be scored against the `norfuraneol ≫ DMHF` ordering.**"

> **"Quantification used a single internal standard, *n*-tetradecane, for all 58 compounds, with no
> compound-specific response factors reported.** Ratios between the two arms survive this; absolute
> levels carry an uncorrected, compound-dependent response bias."

> **"Table 1's Type III furan subtotal at pH 5·0 (1220) does not equal my sum of its own component
> rows (1613).** The monocyclic-furan total and the grand total both close on the printed 1220, so
> the discrepancy is internal to that one subtotal. Unresolved; it affects nothing used elsewhere.
> **Compound number 46 is also absent from the printed table** — the numbering runs 45 → 47."

> **"NAMING TRAP: `HMFone` in this paper means NORFURANEOL, not 5-hydroxymethylfurfural and not
> furaneol.** It is defined on p. 481, four pages after the compound first appears in Table 1. This
> is the second paper in the corpus where an 'HMF'-shaped token means norfuraneol (the other is
> Blank & Fay 1996's `HMF (3)`). The paper separately mentions **5-hydroxymethyl-2-furfural** (the
> real 5-HMF, cited to Tressl 1989) and **2,5-dimethyl-4-hydroxy-3(2H)-furanone** (furaneol, cited
> to Shu & Ho 1988) — three different molecules with confusable short names, in one paper, none of
> them measured except norfuraneol."

---

## 13. WHAT TO FETCH NEXT — ranked

| # | paper | why | confidence |
|---|---|---|---|
| 1 | **Mottram, D. S.; Leseigneur, A. 1990**, *"The effect of pH on the formation of aroma volatiles in meat-like model systems"*, in *Flavour Science and Technology* (Bessière & Thomas, eds.), Wiley, pp. 121−124 | **Norfuraneol quantified across pH 4·5−6·5 in a RIBOSE−CYSTEINE system, DECREASING with pH.** That is (a) a pentose, (b) with cysteine, (c) with a pH axis, (d) on the exact compound the repo's `Furanone_Cyclisation` produces — **four matches to the repo's sulfur lane at once**, and it is the only source in this paper that quantifies a furanone against pH. | Full citation transcribed from the reference list. High on identity; conference-volume availability uncertain. |
| 2 | **Tressl, R.; Helak, B.; Martin, N.; Kersten, E. 1989**, *"Formation of amino acid specific Maillard products and their contribution to thermally generated aromas"*, in *Thermal Generation of Aromas*, ACS, pp. 156−171 | The **5-HMF pH claim** ("large amounts at pH 3, undetectable at pH 5 or above" in glucose−proline) and amino-acid-specific yield data. **Directly serves the HMF channel** (`research_round3_channels.md` §A), which round-3 left with no usable constants. | High. |
| 3 | **Leahy, M. M.; Reineccius, G. A. 1989**, *"Kinetics of the formation of alkylpyrazines"*, in *Thermal Generation of Aromas*, ACS, pp. 196−208 | The title says **kinetics**, on **glucose−lysine at 95 °C and pH 5** — a near-matched hexose comparator to this paper's pentose system, and one of very few pyrazine-kinetics sources. **Same ACS volume as item 2 — one acquisition serves both.** | High. |
| 4 | **Hincelin, O.; Ames, J. M.; Apriyantono, A.; Elmore, J. S. 1992**, *"The influence of xylose on the volatile thermal degradation products of thiamine"*, *Food Chem.* **44**, 381−389 | Same lab, and it reports **1 M aqueous xylose refluxed at pH 7 for 1 h** — i.e. a **third pH point on the same sugar with the same apparatus**, which would extend this paper's two-arm pair into a three-point ladder. | High. |
| 5 | **Nursten, H. E. 1986**, in *Concentration and Drying of Foods*, pp. 53−68 | The source for "2,3-enolisation has a greater browning potential than 1,2-enolisation" — the mechanism both this paper and Blank 1997 lean on. ⚠️ `src/barrier_constants.py` already cites "Nursten 2005" for the `dehydration` class value 28.0; **check whether these are the same author's book and whether the repo's citation is to the 1986 chapter or the 2005 monograph.** | High on identity. |

---

## 14. SOURCES USED IN THIS EXTRACTION BEYOND THE PDF

**None.** Every number and quotation is from `data/articles/apriyantono1993.pdf` — the text layer for
prose, and **260 dpi renders of pages 2, 3 and 4 read visually for every numeric cell of Table 1 and
for the reagent quantities**, because the OCR text layer corrupted exponents and decimal points
(see the preamble). Comparisons are to sibling dossiers in this directory. No web lookup was
performed. Nothing in `src/`, `tests/`, `results/` or `data/benchmarks/` was modified.

*End of dossier. Nothing outside this file was created or modified.*
