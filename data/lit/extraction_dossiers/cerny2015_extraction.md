# Cerny 2015 — EXTRACTION (review chapter) + AUDIT OF THE 11 `ref: cerny2015` CLAIMS
### *"The role of sulfur chemistry in thermal generation of aroma"*
### A 24-page review chapter (Firmenich) with no new measurements. It summarises cysteine Maillard chemistry, lipid–Maillard interaction, thiol stability, cysteine-S-conjugates, Strecker degradation of methionine/SMM, thiamine degradation, Allium and sesame. Every number in it is cited from a primary. **All 11 `docs/validation/directional_claims_panel.yml` quotes attributed to it are verbatim. Two claim statements the chapter does not contain (PH-05, CYS-02). Two misrepresent it: SUG-06, already retired, and PH-04, built from three spliced sources. The PH-01/PH-02 audit opened the wrong "Farmer & Mottram 1990".** Read 2026-10-27.

**Full read of `data/articles/cerny2015.pdf`** (24 pages, Elsevier chapter PDF with a clean text layer
for prose and tables; figures are structure drawings with no text layer). Quotes were checked by
script against two extractions (`pdftotext -layout` and plain `pdftotext`), after normalising
whitespace and dashes and removing running headers. Each "…"-separated fragment was matched
separately. The two quotes that cross a page header or Fig. 9.6 were checked by eye on the page.
Before today the chapter was **not on disk**. The panel's `sources.cerny2015.retrieved_from` names a
scratchpad text ("scratchpad/cerny_chapter.txt") that is not in the repo, and the 2026-10-27 audits
of PH-01, CYS-01 and CYS-02 say "chapter not on disk". This read closes those.

**Provenance:** **[M]** a number printed in the chapter (always **cited** from a primary; the chapter
measures nothing) · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Author | C. Cerny, Firmenich SA, Geneva |
| Title | "The role of sulfur chemistry in thermal generation of aroma" (Chapter **9**) |
| Book | *Flavour Development, Analysis and Perception in Food and Beverages* (Elsevier/Woodhead), pp. **187–210**. **Editors are not printed on the chapter PDF** |
| DOI | **10.1016/B978-1-78242-103-0.00009-6** (printed p. 187 footer: "http://dx.doi.org/10.1016/B978-1-78242-103-0.00009-6"; confirmed on Crossref 2026-10-27: title, container, pp. 187–210, 2015) |
| Year | "Copyright © 2015 Elsevier Ltd." The PDF metadata says "(2014)" (online Nov 2014). Cite as 2015 |
| Keywords (PDF metadata) | Maillard reaction; Cysteine; Thiamine; Allium chemistry; Sesame |
| PDF character | Tables 9.1–9.3; Figures 9.1–9.19 (schemes); ~70 references. **No experimental section, no new data** |

---

## 1. WHAT IT IS AND WHAT IT CONTAINS

A narrative review. It measures nothing, so nothing in it can anchor a constant directly. Its value
to the repo is (a) as the cited source of 11 directional claims (audited in §4), and (b) as a
pointer to primaries. Section map:

| § | topic | primaries leaned on |
|---|---|---|
| 9.1 | thiol vs alcohol odour thresholds; thialdine as an SDE artefact | van Gemert 2003; Siegmund 1997; Mottram & Mottram 2002 |
| 9.2.1 | cysteine Maillard: 2-acetylthiazole/-thiazoline, acetyl-dihydrothiazine, ribose/Cys key odorants, 1,4-dideoxypentosone → MFT and 3-mercapto-2-pentanone, the dione + H2S route to α-mercaptoketones, IMP vs ribose, furfural → FFT, glucose/Cys, norfuraneol and DMHF + H2S/Cys, pH effects | Mulders 1973; Hofmann & Schieberle 1995, 1997, 1998a,b; Cerny & Davidek 2003, 2004; Mottram et al. 1995; Mottram & Nobrega 2002; Shibamoto 1977; Whitfield & Mottram 1999, 2001; Shu & Ho 1988; Shu et al. 1985; Madruga & Mottram 1995; Farmer & Mottram 1990 (book chapter, see §4); Cerny & Briffod 2007 |
| 9.2.2 | lipid–Maillard interaction | Farmer et al. 1989; Farmer & Mottram 1990; Whitfield 1988; Elmore et al. 1997, 2002; Elmore & Mottram 2000; Mottram & Whitfield 1995 |
| 9.2.3 | thiol and disulfide stability | Hofmann et al. 1996; van Seeventer 2001; Hofmann & Schieberle 2002; Guth 1995; Mottram et al. 1996 |
| 9.2.4 | cysteine-S-conjugates | Cerny & Guntz-Dubini 2013 (`cerny2013_extraction.md`) |
| 9.3 | Strecker: methionine, 3-oxazolines, SMM → DMS | Buhr 2010; Granvogl 2012; Scherb 2009; Cerny et al. 2011 |
| 9.4 | thiamine → 5-hydroxy-3-mercapto-2-pentanone → MFT; alkaline route to sulfurol/thiazoles | Grosch et al. 1993; Van der Linde 1978; Guentert 1990, 1993, 1994; Cerny & Guntz-Dubini 2008 |
| 9.5–9.6 | Allium; roasted sesame | Granvogl 2004; Krammer 2006; Roessner 2000; Schieberle 1996; Tamura 2010, 2011 |

---

## 2. TABLES AS PRINTED

### Table 9.1 "Odour thresholds of thiols and their corresponding alcohols (van Gemert, 2003)" [M, cited]

| thiol | odour threshold (µg/kg) | alcohol | odour threshold (µg/kg) |
|---|---|---|---|
| Ethanethiol | 2 | Ethanol | 3500–190,000,000 |
| 2-Furanmethanethiol (furfural mercaptan) | 0.005–0.12 | 2-Furanmethanol | 1900–2000 |
| 3-Mercapto-2-butanone | 3 | 3-Hydroxy-2-butanone | 14–10,000 |
| Methanethiol | 0.02–4 | Methanol | 14–10,000 |
| 2-Phenylethanethiol | 0.05^a | 2-Phenylethanol | 0.015–75,000 |
| 3-Methyl-2-buten-1-thiol (prenyl thiol) | 0.0002–0.005 | 3-Methyl-2-buten-1-ol (prenyl alcohol) | 250–7800 |
| 4-Mercapto-4-methyl-2-pentanone | 0.0001–0.005 | 4-Hydroxy-4-methyl-2-pentanone | 44,000–100,000 |
| 1-p-Menthen-8-thiol | 0.0001^b | α-Terpineol | 4.6–15,000 |

^a Schieberle (1996). ^b Buettner and Schieberle (2001). **[!]** The methanol range (14–10,000) is
identical to the 3-hydroxy-2-butanone range printed directly above it, which looks like a
copy-down error.

### Table 9.2 "Aroma-active sulfur compounds in thermally processed foods (examples)" [Q]

Qualitative (compound / odour / foods). The model-relevant rows: **2-furanmethanethiol**
"coffee-like" in beef, chicken, coffee, peanut, pork, salmon, sesame. **2-Methyl-3-furanthiol**
"meaty" in beef, chicken, coffee, pork, sesame. **2-Methyl-3-(methyldithio)furan** "meaty" in beef,
cocoa. 2-Acetylthiazole "roasted" in asparagus, beef, chicken, pork, salmon, sesame. Methional,
methanethiol, DMS and dimethyl trisulfide also appear, with food-source footnotes a–n.

### Table 9.3 "Aroma-active sulfur compounds from Maillard reaction and their precursors in food" [Q]

| compound | precursor(s) |
|---|---|
| 2-Acetylthiazole | Cysteine/monosaccharides |
| 2-Acetyl-2-thiazoline | Cysteine/monosaccharides |
| Dimethyl sulfide | S-Methylmethionine |
| 2-Furanmethanethiol | Cysteine/ribose |
| Methanethiol | Methionine/monosaccharides |
| 3-Mercapto-2-butanone | Cysteine/monosaccharides |
| 2-Mercapto-3-pentanone | Cysteine/ribose, thiamine |
| 3-Mercapto-2-pentanone | Cysteine/ribose |
| Methional | Methionine/monosaccharides |
| 2-Methyl-3-furanthiol | Thiamine, cysteine/ribose |

**[!]** Table 9.3 gives thiamine as a precursor of **2-mercapto-3-pentanone** and not of
**3-mercapto-2-pentanone**. Its own §9.2.1 text (p. 196, on Cerny & Briffod 2007) says
3-mercapto-2-pentanone comes "primarily from thiamine". The model labels `MP2P`
(3-mercapto-2-pentanone) "THIAMINE-diagnostic" and `MP3P` (2-mercapto-3-pentanone)
"NORFURANEOL-diagnostic", which is the **opposite** assignment to Table 9.3. The primaries side with the model and the text. Cerny & Briffod 2007 Table 4 finds
2-mercapto-3-pentanone absent at pH 4–5.5 and 94–95 % xylose-derived at pH 6–7, never
thiamine-derived. Cerny 2007b finds 3-mercapto-2-pentanone >99 % thiamine-derived without cysteine.
So Table 9.3 puts "thiamine" on the wrong isomer. Do not use it.

---

## 3. NUMBERS IN THE TEXT [M, all cited]

| statement | value | cited primary | p. |
|---|---|---|---|
| FFT odour threshold vs 2-furanmethanol | 0.005–0.12 µg/kg; "10,000–100,000-times lower" | van Gemert 2003 | 187 |
| 2-acetyl-2-thiazoline in pressure-cooked chicken / fried beef | 14 / 28 µg/kg | Kerler & Grosch 1997; Cerny & Grosch 1993 | 190 |
| ribose + Cys, labelled ribose | carbon skeleton of ribose intact in MFT, FFT and 3-mercapto-2-pentanone at pH 5 | Cerny & Davidek 2003 | 191 |
| 2-mercapto-3-pentanone absent | ribose + Cys, **95 °C, pH 5, 4 h** | Cerny & Davidek 2003 | 192–193 |
| only 3-mercapto-2-pentanone formed | Cys : ribose **1 : 3**, 95 °C 4 h (Cerny & Davidek 2003) or **145 °C 20 min** (Hofmann & Schieberle 1995) | as stated | 193 |
| both isomers formed, the 2,3-pentanedione route dominant | **equimolar** Cys + ribose, **140 °C, 30 min** | Mottram & Nobrega 2002 | 193 |
| raw meat pH | "generally around 5.6" | — | 195 |
| minced meat pre-adjusted to pH 4.0–5.6 | more volatiles at lower pH; MFT and FFT "in general, higher" at low pH | Madruga & Mottram 1995 | 195 |
| MFT in diethyl ether | −30 % after 1 day at 6 °C, with bis(2-methyl-3-furyl) disulfide forming | Hofmann et al. 1996 | 198 |
| process flavour, 50 °C, 1 day | MFT −59 %, FFT −28 %; "L-Cysteine had a stabilising effect on the thiols" | van Seeventer et al. 2001 | 198 |
| FFT in a coffee brew | −75 % after 60 min in a thermos flask | Hofmann & Schieberle 2002 | 198 |
| **Cys 0.65 mM + ribose 1.34 mM, boiled** | **MFT 0.2 µg/L**; ribose replaced by **thiamine 0.06 mM**: "the yield increased **1500 times** to **300 µg/L**" | Grosch et al. 1993 | 201 |
| S-methylmethionine | tomato 2.8, malt 6, cabbage 81, celery 176, asparagus 53–252 mg/kg | Scherb 2009 | 201 |
| alliin / isoalliin | garlic up to 2.4 g/kg; onion (~)2 g/kg | Whitaker 1976 | 203 |
| 3-mercapto-2-methylpentan-1-ol | raw onion 8–32, cooked 34–246 µg/kg; threshold 0.15 µg/L | Granvogl 2004 | 203 |
| 2-methylthioacetaldehyde / 2-allylthioacetaldehyde | broccoli 0.5 mg/kg / garlic 1 mg/kg | Roessner 2000 | 205 |
| sesame alkenethiols | 20–800 µg/kg, OAV 27–2400 | Tamura 2011 | 206 |

**Mechanistic statements the model uses** [Q]: DDP + H2S → 5-hydroxy-3-mercapto-3-penten-2-one →
MFT, and its Strecker-type branch → 3-mercapto-2-pentanone (Fig. 9.5, Cerny & Davidek 2003).
Thiamine → 5-hydroxy-3-mercapto-2-pentanone (isolated by Cerny & Guntz-Dubini 2008) →
4,5-dihydro-2-methyl-3-furanthiol → MFT (Fig. 9.14). Furfural → FFT with H2S (Shibamoto 1977).
Hydroxyacetaldehyde + mercaptopropanone → MFT and FFT in glucose/Cys (Fig. 9.6, Hofmann & Schieberle
1998). Dione + H2S → α-mercaptoketones (Mottram et al. 1995). FFT irreversibly bound by pyrazinium
radical cations in coffee melanoidins (Fig. 9.9). Disulfide → thiol by exchange with protein cysteine
(Fig. 9.10, Mottram 1996). FFA + Cys → S-conjugates (Fig. 9.11).

---

## 4. AUDIT OF EVERY CLAIM CITING `ref: cerny2015` (`docs/validation/directional_claims_panel.yml`)

**Two separate questions per claim:** (1) **Quote fidelity**: is the `quote` in the chapter?
(2) **Claim support**: does the chapter say what the claim's `statement` and scored arms assert?
Verdicts: **VERBATIM** (quote exact), **PARAPHRASE OK** (statement is a fair reading),
**NOT FOUND** (the statement's content is not in the chapter), **MISREPRESENTS** (the quote is
real, but the statement says something different or is assembled out of context). The panel was
not edited.

| claim_id | fit_status (panel) | quote fidelity | claim support | primary the chapter cites for the quoted text | notes |
|---|---|---|---|---|---|
| **SUG-03** | independent | **VERBATIM** (p. 193; one contiguous sentence; the `verbatim_page_elided` tag is unnecessary because the quote stops before "(Hofmann and", i.e. before the page break) | **PARAPHRASE OK, stretched** | Hofmann & Schieberle **1997** (glucose/cysteine and rhamnose/cysteine, AEDA) | "played a more important role" is **aroma importance (AEDA/FD)**, not concentration. The statement turns it into "FFT > MFT" in amount. The quantitative primary Hofmann & Schieberle 1998b Table 1 agrees (glucose FFT 2.8 > MFT 1.9 µg; `hofmann1998_reconciliation.md` l. 153). The claim's arms (10 : 10 mM, 140 °C, 30 min) are not that primary's pot |
| **SUG-04** | fit_adjacent | **VERBATIM** (p. 193–194; the ellipsis start is fine; "(Hofmann and \| [p. 194 header + Fig. 9.6] \| Schieberle, 1997)" is verified by eye) | **PARAPHRASE OK, stretched** | in-text: Hofmann & Schieberle **1997**. **[!]** That paper is the glucose/rhamnose study; the pentose statement belongs to Hofmann & Schieberle **1995** (ribose/Cys AEDA, cited on p. 191) | "key" ≠ "more than FFT". The arms' primary (Hofmann 1998b T1: ribose MFT 19.8 > FFT 12.1 µg) supports A > B at that pot. The audit's caveat stands: the ordering flips dry-heated and at pH 5.5–6.0 |
| **SUG-05** | fit_system_overlap | **VERBATIM** (two fragments, p. 193, ~4 lines apart) | **PARAPHRASE OK** | Shibamoto 1977 (furfural + H2S + NH3, no cysteine); the glucose-for-furfural substitution is **uncited in the text** (the reference list carries Shibamoto & Russell 1977, glucose–H2S–NH3, which is presumably it); "less effective … than pentoses (ribose, xylose)": **Hofmann & Schieberle 1998b** | grammatically "less effective as a precursor" covers "furfural and 2-furanmethanethiol", so furfural is included. The audited arms use Hofmann 1998b Table 5 (ribose 67.5 vs glucose 6.9 µg furfural): fine |
| **SUG-06** | retired | **VERBATIM** (p. 193) | **MISREPRESENTS** | Hofmann & Schieberle 1998b | the chapter **groups** ribose and xylose ("pentoses (ribose, xylose)"); it never says they are not interchangeable. The 2026-10-27 retirement is **confirmed** on the full text |
| **PH-01** | fit_system_overlap | **VERBATIM** (p. 195, two consecutive sentences) | **PARAPHRASE OK** | sentence 1: the **Madruga & Mottram 1995** meat-pH context (minced meat pH 4.0–5.6); sentence 2: **Farmer & Mottram 1990 = the BOOK CHAPTER** "Recent studies on the formation of meat-like aroma compounds", in Bessière & Thomas (eds), *Flavour Science and Technology*, Wiley, Chichester, **pp. 113–116** | **[!] The 2026-10-27 audit opened the wrong paper.** It says "the cited Farmer & Mottram 1990 (data/articles/farmer1990.pdf) is a lipid study at fixed pH 5.6". That PDF is *J. Sci. Food Agric.* **53**, 505–525 (`farmer1990_extraction.md`), a **different 1990 publication**. The chapter's reference list (p. 207) cites the Wiley proceedings chapter, which is **not on disk** (Zheng 1997 cites the same one). The re-sourcing of the arms to Meynier & Mottram 1995 Table 3 is still sound, but the stated reason ("the cited paper has no pH series") is not established: the cited paper was never read |
| **PH-02** | fit_system_overlap | **VERBATIM** (p. 195) | **PARAPHRASE OK** | Farmer & Mottram 1990 **book chapter** (as PH-01) | "their corresponding disulfides" covers bis(2-methyl-3-furyl) disulfide. Same wrong-paper issue as PH-01. Meynier 1995 gives only trace-level disulfide evidence (as the audit says) |
| **PH-03** | not_independent | **VERBATIM** (p. 195–196; a page header sits between "3-mercapto-2-pentanone" and "and 4,5-dihydro…", verified by eye) | **PARAPHRASE OK** | **Cerny & Briffod 2007** (JAFC 55, 1552) | faithful. The chapter's next sentence ("3-mercapto-2-pentanone primarily from thiamine") is right for C&B's dilute pot (77–90 % thiamine) and wrong for Cerny 2007b's concentrated pot (60 %; `cerny2007b_extraction.md`), but it is not part of the claim |
| **PH-04** | independent | **VERBATIM fragments, but stitched**: three sentences from **three different paragraphs and model systems**, quoted **out of page order** (all on p. 195, where fragment 2 comes first, fragment 3 second and fragment 1 last, in a later paragraph) | **MISREPRESENTS (splice)** | fragment 1 "Higher pH produced more thiazoles": **Farmer & Mottram 1990** (Cys + ribose); fragment 2 "At pH 4.5, sulfur compounds dominated … pH 6.5 nitrogen compounds prevailed": **Whitfield & Mottram 1999/2001, NORFURANEOL** + H2S or Cys; fragment 3 "at pH 7.1, heterocyclic compounds (thiazoles, thiazolines, pyrazines …) dominate": **Shu & Ho 1988, DMHF + Cys** | the claim scores **2,5-dimethylpyrazine in Cys + ribose** (pH 4.5/5.6/6.5). The chapter **never says pyrazines rise with pH in Cys + ribose**: the only pyrazine mention is the DMHF + Cys system, and the "nitrogen compounds" sentence is about norfuraneol. The ellipses hide that the three subjects differ. **The direction is nonetheless supported by a primary on disk:** Meynier & Mottram 1995 **Table 3 (Cys + ribose)**: dimethylpyrazine nd / nd / nd / 0.3 / 1.5 µg per 10 mg ribose at pH 4.5 / 5.0 / 5.5 / 6.0 / 6.5 (`meynier1995_extraction.md` §3). **Recommend re-sourcing PH-04 to Meynier 1995 T3** and moving its arms to that pot (50 : 35 mM, 0.2 M pyrophosphate, 140 °C, 1 h), as PH-01 was |
| **PH-05** | independent | **VERBATIM** (p. 195) | **NOT FOUND** | Whitfield & Mottram 1999/2001 (norfuraneol) | the quote is about sulfur vs nitrogen volatiles from **norfuraneol**. **The chapter says nothing anywhere about furfural (or HMF) versus pH.** The furfural direction is the claim-writer's "acid-catalysis corollary", as the panel's own `fit_note` admits. It should not be counted as sourced to this chapter (or at all) in the independent headline |
| **CYS-01** | not_independent | **VERBATIM** (three fragments: p. 187 ×1 with an internal elision, p. 193 ×2) | **PARAPHRASE OK (weak)** | Mottram & Mottram 2002 (review); Shibamoto 1977 | the chapter calls cysteine "the most important precursor for Maillard reaction-derived sulfur volatiles" (p. 188) but **contains no ± cysteine comparison**. The audit's "(chapter not on disk)" inference is now **confirmed**. The tautology verdict stands |
| **CYS-02** | retired | **VERBATIM** (p. 193; "α" written as "alpha-") | **NOT FOUND** | Mottram, Madruga & Whitfield 1995 (alkanediones + H2S) | the chapter has **no mention of HMF**; "browning" appears only in the title of the Mulders 1973 system. Nothing about cysteine suppressing dehydration. The retirement's "(inferred; the chapter is not on disk)" is now **confirmed** on the full text |

**Tally.** Quote fidelity: **11/11 VERBATIM**. No quote is fabricated or altered beyond declared
elisions; SUG-05, PH-04 and CYS-01 join non-contiguous fragments with "…". Claim support:
**PARAPHRASE OK 7** (SUG-03, SUG-04, SUG-05, PH-01, PH-02, PH-03, CYS-01; SUG-03/04 stretched from
"importance" to amount), **NOT FOUND 2** (PH-05, CYS-02), **MISREPRESENTS 2** (SUG-06, PH-04).
**Of the four NOT FOUND / MISREPRESENTS rows, only PH-04 and PH-05 still sit in the `independent`
headline** (SUG-06 and CYS-02 are retired; SUG-03 is independent and PARAPHRASE OK). PH-04's direction is
rescuable from Meynier 1995 T3. PH-05 has no source.

**Citation-level finding to carry to the panel's source notes:** the "Farmer and Mottram, 1990" that
the chapter cites for the Cys + ribose pH 4.5–6.5 thiol series (PH-01, PH-02, PH-04 fragment 1) is
the Wiley *Flavour Science and Technology* proceedings chapter (pp. 113–116). It is not
`data/articles/farmer1990.pdf` (JSFA 53:505). That primary has never been read in this repo.

---

## 5. CHAPTER ERRORS FOUND ON THE WAY [!]

1. **Table 9.3 swaps the mercaptopentanone precursors** relative to the chapter's own text and the
   primaries (§2).
2. **Fig. 9.10 caption** says "3-methyl-3-furanthiol (16)"; compound 16 is 2-methyl-3-furanthiol.
3. **§9.2.4** says rearrangement gives "S-(2-methyl-3-furyl)-L-cysteine (36)". 36 is the
   intermediate; the product is **38**.
4. **p. 195** gives Shu & Ho's main product as "2,5-dimethyl-4-hydroxy-3(2H)-thiophene"; it is the
   **thiophenone** (the model's `DMHFS`).
5. The acetyl-thiazine is "5-acetyl-3,6-dihydro-2H-1,4-thiazine" on pp. 190–191 and
   "5-acetyl-3,4-dihydro-2H-1,4-thiazine" on pp. 191 and 193.
6. **p. 191** calls 1,4-dideoxypentosone "the alleged key intermediate for 2-methyl-3-furanthiol
   **and 2-furanmethanethiol**". On p. 193 the chapter makes furfural the key FFT intermediate,
   and Fig. 9.5 shows only MFT and 3-mercapto-2-pentanone from DDP.
7. References: Cerny & Davidek 2003 is printed without "cysteine" in the title and with pp.
   2714–2721; Cerny & Guntz-Dubini 2013 gives pp. 2714–2718. Starkenmann et al. 2008 is printed as
   vol. "96"; Cerny 2013 has 56.
8. Table 9.1's methanol threshold range duplicates the row above (§2).
9. SUG-04's source sentence attributes a **pentose** statement to the **glucose/rhamnose** paper
   (Hofmann & Schieberle 1997), as noted in §4.

---

## 6. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Cannot anchor anything numerically.** It is a review, and every number is a citation. Use it only
to route to the primary.

**Useful pointers (not anchors):**
1. **Grosch et al. 1993** (Progress in Flavor Precursor Studies, Allured, pp. 329–342): boiled
   Cys 0.65 mM + ribose 1.34 mM gives MFT 0.2 µg/L; Cys 0.65 mM + thiamine 0.06 mM gives 300 µg/L
   (×1500). This is a **low-concentration thiamine-vs-ribose MFT contrast** at near-meat precursor
   levels, the regime where the model's first-order thiamine route should dominate the second-order
   sugar route. Time, pH and buffer are not given in the chapter. The primary is not on disk. If
   retrieved, it is a candidate HOLD-OUT for the route mix (cf. the dilution trend in
   `cerny2007b_extraction.md` §5).
2. **van Seeventer 2001** (MFT −59 %, FFT −28 %, 1 day, 50 °C, process flavour) and **Hofmann &
   Schieberle 2002** (FFT −75 % in 60 min, coffee thermos) are already the model's
   `ch_oligomer_*` / coffee-binding evidence. The chapter adds that "L-Cysteine had a stabilising
   effect on the thiols" (van Seeventer). That is a cysteine-dependent **suppression** of thiol loss,
   the opposite sign to what a cysteine-S-conjugation sink would give (`cerny2013_extraction.md` §4).
3. **Mottram & Nobrega 2002 vs Cerny & Davidek 2003 / Hofmann & Schieberle 1995**: whether
   2-mercapto-3-pentanone forms depends on the Cys : ribose ratio (1 : 1 → both isomers; 1 : 3 → only
   3-mercapto-2-pentanone). That is a ratio-dependent isomer switch the `MP3P`/`MP2P` split could be
   tested on, from the primaries.
4. The **Farmer & Mottram 1990 Wiley chapter** is the actual cited primary for the Cys + ribose
   pH 4.5–6.5 thiol series; worth fetching.

---

## 7. CAVEATS

1. Review chapter: no conditions, units or uncertainties for most cited numbers. Never use a number
   from here without its primary.
2. Several mechanism drawings (Figs. 9.4–9.6, 9.11, 9.14) are "hypothesised" schemes, as the chapter
   says.
3. Quote matching used two text extractions plus page inspection for the two page-spanning quotes.
   Typographic normalisation (en-dash ↔ hyphen, α ↔ "alpha-") was applied, consistent with the
   panel's declared `quote_form` values.
