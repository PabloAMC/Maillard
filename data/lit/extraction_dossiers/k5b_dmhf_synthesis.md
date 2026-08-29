# Wave K5b — DMHF / furanone cluster synthesis

**Date:** 2026-08-29. **Wave:** K5b (extraction). **Scope:** the five-paper DMHF/furanone cluster.
**Nothing outside `data/lit/extraction_dossiers/` was created or modified. No commit.**

**Per-paper dossiers written by this wave:**
`blank1997_extraction.md` · `wang2008_extraction.md` · `poisson2019_extraction.md` ·
`apriyantono1993_extraction.md` · `shu1988_extraction.md`
**Read first:** `blank1996_extraction.md`, `research_round3_channels.md` §B.

---

# §0. FILE IDENTITIES — **ALL FIVE FILES ARE THE CORRECT PAPERS. NO WRONG-FILE PROBLEM.**

The brief asked for wrong-file identities to be reported prominently. **There are none.** Every file
was verified from its own header/footer, not from a lookup:

| file | expected | verified from the document itself | ✔ |
|---|---|---|---|
| `blank1997.pdf` | Blank, Fay, Lakner & Schlosser 1997, JAFC 45:2642−2648, IDA companion to blank1996 | Title, four authors, both affiliations, journal, volume 45, issue 7, pages 2642−2648. **Article ID `JF960997I` and CCC line `S0021-8561(96)00997-1` printed on the paper.** | ✔ |
| `wang2008.pdf` | `10.1021/jf8012025` Wang & Ho, DMHF via methylglyoxal | Title, both authors, JAFC 56(16):7405−7409, **article ID `JF8012025`** and the DOI printed in the first-page footer. | ✔ |
| `poisson2019.pdf` | `10.1021/acs.jafc.9b00770` coffee CAMOLA | Title, six authors, three affiliations, **DOI printed on every page footer**. ⚠️ It is an ASAP proof: volume/issue/pages read `XXXX, XXX, XXX−XXX`. | ✔ |
| `apriyantono1993.pdf` | `10.1002/jsfa.2740610416` xylose−lysine pH | Title, both authors, *J. Sci. Food Agric.* 1993, 61, 477−484, Reading. **No DOI printed (1993).** | ✔ |
| `shu1988.pdf` | Shu & Ho 1988, JAFC 36:801−803, fed-DMHF + cysteine | Title, both authors, JAFC 36(4):801−803, Rutgers. **No DOI printed (1988); none looked up.** | ✔ |

### ★ ONE POSITIVE IDENTITY UPGRADE

`blank1996_extraction.md` §0a/§12 flagged the companion's DOI as *"~85 % confidence on the
DOI/pagination — taken from a web bibliographic lookup … NOT CrossRef-verified and NOT read."*
**`10.1021/jf960997i`, JAFC 45:2642−2648, is now VERIFIED FROM THE PRIMARY DOCUMENT** (article ID
`JF960997I`). That flag can be cleared.

### ⚠️ THREE PDFs REQUIRE VISUAL READING — a transmission warning for any future wave

| file | what the text layer gets wrong | what I did |
|---|---|---|
| **`wang2008.pdf`** | **Figure 1 — the paper's ONLY quantitative content — has no text layer at all.** `pdftotext -bbox` shows the plot region is empty of glyphs. | Digitised by pixel measurement from a 300 dpi render; axis calibration verified against the tick marks to **0.4 %**. |
| **`apriyantono1993.pdf`** | OCR corrupts Table 1's **superscript exponents and decimal points**: `12.9 × 10³` → `129 x 103`; `3.53 × 10⁶` → `3.53 x 108`; `66.2` → `662`; `6.8` → `68`. | Every numeric cell re-read visually from 260 dpi renders of pp. 3−4; all group totals re-derived arithmetically. |
| **`shu1988.pdf`** | OCR **collapses Table I's three pH columns** wherever a row has empty cells — e.g. it renders the trithiolane row as `3.1 / 3.2` with `24.7` orphaned onto the previous line. **Figure 1 has no numeric text layer.** | Table I re-read visually from a 300 dpi render; Figure 1 digitised, axis calibration verified to **0.2 %**. |

**A naive text-layer harvest of any of these three would produce wrong numbers.** Recorded so it is
not repeated.

---

# §1. ★★ THE HEADLINE: WHAT CHANGED

**Before this wave**, `research_round3_channels.md` §B concluded, honestly and correctly on the
evidence then available:

> *"the DMHF node can be built with correct topology and a pH sign, but **it cannot be given a
> barrier from literature**. Any Ea the repo assigns to DMHF is the repo's own assumption."*
> — §B.3

**After this wave, that verdict splits in two:**

| | round-3 position | **K5b position** |
|---|---|---|
| **Topology** | correct, from Blank & Fay 1996 + two abstracts | **★★ Much stronger.** Two formation edges structurally separated by a measured intermediate test; one sink edge structurally established. |
| **Magnitude (formation, pentose)** | **none** | **★★ 39 absolute SIDA cells in µg/mmol sugar, σ ≤ 10 %, across five parameter axes** (Blank 1997). |
| **Magnitude (formation, MGO)** | **none** | **⚠️ Nine cells, digitised from a bar chart, external-standard HPLC** (Wang & Ho). |
| **Magnitude (sink)** | **none** | **❌ STILL NONE.** Shu & Ho give the product spectrum and the pH shape, not the consumption. |
| **Rate / Ea** | **none** | **❌ STILL NONE, on every edge.** All five papers are **single-temperature**. |

**★ THE ONE-LINE ANSWER: the DMHF channel is now CALIBRATABLE IN SELECTIVITY AND MAGNITUDE ON ITS
PENTOSE FORMATION EDGE, PARTIALLY ON ITS MGO EDGE, AND NOT AT ALL ON ITS SINK OR ON ANY RATE.
Round-3's verdict on barriers stands unchanged and unchallenged: no Arrhenius parameter for any
furanone family exists in this literature.**

---

# §2. ⚠️ TWO PRIOR HOLD-OUT PROPOSALS ARE CONTRADICTED — READ BEFORE DECLARING ANYTHING

`blank1996_extraction.md` §10.2 proposed six hold-outs from GC-O data. **The 1997 isotope-dilution
companion measures through two of them and reverses their sign.** Both failures have the same
cause: **a `−` in a 1996 GC-olfactometry table is a non-detection at a sniffing port, not a zero.**

| prior proposal | prior role | K5b evidence | **status** |
|---|---|---|---|
| **#7 — "HEMF requires alanine; ABSENT with glycine"** | HOLD-OUT, structural truth table: *"a model that emits HEMF from pentose+glycine fails"* | **Blank 1997 Table 1: EHMF = 0.3 (xylose/Gly), 0.7 (ribose/Gly), 1.3 (arabinose/Gly) µg/mmol.** Non-zero in all three glycine systems. | **★ MUST NOT BE DECLARED AS WRITTEN.** A model emitting HEMF from pentose+glycine is **correct**. Demote to a **≈10−25× preference**. |
| **#8 — "DMHF forms from pentose ALONE"** | HOLD-OUT, structural: *"a zero-amino-acid POSITIVE control"* | **Blank 1997 Table 1 fn c and Table 3 rows 1−2: `< 0.01 µg/mmol` for both furanones, in phosphate and in water** — **≥260× below** the xylose/glycine cell, reported only as an upper bound. | **★ SIGN REVERSED.** Re-specify as a **CEILING / negative control**. |
| #9 — norfuraneol ≫ DMHF | HOLD-OUT, ordinal ceiling | Blank 1997 p. 2646: *"in all samples analyzed, 4-hydroxy-5-methyl-3(2H)-furanone was the main reaction product (**data not shown**)"* | **✔ CORROBORATED across two papers and six systems. Still unquantified in both.** Keep as an **ordinal-only** ceiling. |
| #4/#5/#6 — 70/30 Strecker/fragmentation split | FIT, as a branch prior | **Blank 1997's 4-aminobutyric-acid control gives 73/27 by a NON-ISOTOPIC method** (§3.3 below) | **★★ CORROBORATED INDEPENDENTLY. Strengthen, don't weaken.** |
| #3 — C5+C1/C5+C2, backbone intact | HOLD-OUT, atom-mapping | Blank 1997 supports it quantitatively; adds no new isotope evidence | ✔ unchanged |
| #14 — pH drift 6.0 → 5.0/5.3 | HOLD-OUT, weak | Not repeated; **Blank 1997 reports no final pH at all** | unchanged, still a single observation |

### ⚠️ AND ONE ROUND-3 REGISTER ENTRY IS WRONG

`research_round3_channels.md` §C.2 (and its §E download-list row) records Apriyantono & Ames 1993 as
**"RATIO-ONLY"**, on the basis of the abstract's `522 / 999 g kg⁻¹` volatile shares. **The paper's
Table 1 reports all 58 compounds in `nmol mol⁻¹ xylose` — absolute molar yields on the limiting
sugar, against an internal standard, with a stated ±16 % precision.** The furfural pH effect is
**274×**, not the ≈1.9× the mass shares imply, because the denominator collapsed 143× at the same
time. **`research_round3_channels.md` §C.1's declared gap — "I found no paper reporting furfural
yields from a pentose … That is a genuine gap" — is closed by this paper.** See
`apriyantono1993_extraction.md` §0a.

---

# §3. ★★ THE CHANNEL DESIGN THE EVIDENCE SUPPORTS

Three edges. **They are structurally separated by measurements, not by assumption.**

```
                                       ┌──── EDGE A ────────────────────────────────────┐
   pentose ─(Amadori, 2,3-enol)→ 1-deoxypentosone ─+Strecker aldehyde(aldol)→ acetylformoin ─(reduction)→ DMHF
   hexose  ─(Amadori, 2,3-enol)→ 1-deoxyhexosone  ────────────────────────────→ acetylformoin ─(reduction)→ DMHF
                                       └───── intact sugar skeleton; acetylformoin IS on this route ─────┘

                                       ┌──── EDGE B ────────────────────────────────────┐
   2 × MGO ─(Cannizzaro, base-cat.)→ MGO + acetol ─(aldol)→ 2,5-dioxo-3,4-dihydroxyhexane ─(cyclise)→ DMHF
                                       └───── C3+C3 recombination; acetylformoin is NOT on this route ────┘

                                       ┌──── EDGE C (SINK) ─────────────────────────────┐
   DMHF + H2S/cysteine ──→ 2,5-dimethyl-4-hydroxy-3(2H)-THIOPHENONE + H2O   (+ 2,4-hexanedione, thiophenes)
                                       └───── acid/neutral-favoured; markers vanish by pH 7.1 ───────────┘
```

### 3.1 ★★ THE MEASUREMENT THAT SEPARATES EDGE A FROM EDGE B

**Wang & Ho 2008, experiment B** (`wang2008_extraction.md` §4.2): feeding `[¹³C₆]glucose +
[¹²C₃]MG (1:1)` gives a **4:1 mixture of `[¹³C₆]`DMHF and `[¹²C₆]`DMHF** — so 20 % of the DMHF came
from two MG units recombining — **and NO `[¹²C₆]`acetylformoin was observed.**

> "**no [¹²C₆]acetylformoin was observed suggesting that acetylformoin was not a precursor during
> the DMHF formation from MG** … Glucose kept carbon skeleton intact during DMHF formation even its
> fragment MG was present, which indicated that **MG and cyclization of intact glucose pathways were
> parallel since the precursors of these two pathways were different.**"

**★ THE TWO EDGES DO NOT SHARE AN INTERMEDIATE.** That is a hard structural result, and it also
**resolves an ambiguity Poisson 2019 explicitly leaves open**: Poisson's Figure 5 proposes *both*
a C3/C3 route through acetylformoin (pathway b, DHA + 2-oxopropanal) *and* one bypassing it
(pathway c, acetol + 2-oxopropanal), and its CAMOLA cannot tell them apart because both produce the
same [M+3] peak. **Wang & Ho's acetylformoin null selects pathway (c).**

**A third confirmation that the edges are non-mixing:** in both of Wang & Ho's CAMOLA experiments,
**only `[¹³C₆]` and `[¹²C₆]` isotopomers were ever observed — never `[¹³C₃]`.** No DMHF is assembled
from one hexose-derived C3 plus one MG-derived C3.

### 3.2 ★★ EDGE B IS BOUNDED FROM BOTH SIDES BY THREE INDEPENDENT MEASUREMENTS

| system | conditions | **C3/C3 share of DMHF** | source |
|---|---|---|---|
| **glucose only, aqueous** | 1.4 M glucose + 1 M Gly or Cys, 0.5 M phosphate, **120 °C**, pH 5, 60 min | **below detection** — `[¹³C₁]`−`[¹³C₅]` all absent | Wang & Ho, exp. A |
| **coffee bean, real matrix** | biomimetic recombinate, pH 5.5 initial, 10 % moisture, **171−220 °C dry** | **8.2 / 12.7 / 11.4 %** of the sucrose-derived HDMF, **at 210/260/300 s ONLY**; below the 0.5 % floor at all other roast times | Poisson, Table 10 |
| **glucose + spiked MG, aqueous** | 1.4 M glucose **+ 1.4 M MG**, 120 °C, pH 5, 60 min | **20 %** — and this is a **LOWER bound** (authors attribute the shortfall to Strecker consumption of MG) | Wang & Ho, exp. B |

**★★ THREE INDEPENDENT POINTS, TWO LABS, TWO MATRICES, A 100 °C TEMPERATURE SPREAD — and all three
say the same thing: Edge B is REAL, it is CONCENTRATION-CONTROLLED, and at in-situ MGO levels it is
MINOR.** A model that gives Edge B a large share in a plain sugar+amine system contradicts Wang &
Ho's experiment A; a model that omits Edge B entirely contradicts both his experiment B and
Poisson's mid-roast window. **The edge is bracketed, which is a better position than a point
estimate.**

### 3.3 ★★ EDGE A'S INTERNAL BRANCH FRACTION IS NOW A TWO-METHOD RESULT

`blank1996_extraction.md` items #4/#6 established, by ¹³C labelling, that in xylose/glycine at 90 °C
**~70 % of HDMF carries the Strecker-derived C1 and ~30 % comes from sugar fragmentation.**

**Blank 1997 Table 2 reproduces that non-isotopically.** 4-Aminobutyric acid is a γ-amino acid and
*"does not decompose by Strecker deamination"*, so `xylose/GABA` is the fragmentation-only baseline
and the increment on adding a Strecker-active amino acid at the same molarity is the Strecker
channel:

| product | GABA only | + Strecker-active AA | **Strecker share** | **fragmentation share** |
|---|---|---|---|---|
| **HDMF** (+ glycine) | 0.4 | 1.5 | **73 %** | **27 %** |
| **EHMF** (+ alanine) | 0.1 | 3.2 | **97 %** | **3 %** |

**★★ 73/27 by amino-acid substitution vs 70/30 by isotopomer distribution — three percentage points
apart, two completely different methods, on the same system.** This is the single strongest
corroboration in the cluster and it makes the 70/30 split safe to use as a branch prior.

Two weaker estimates from the same data bracket it: **15 %** (GABA/xylose-Gly ratio, 0.4/2.6) and
**24−44 %** (alanine-baseline). **Consolidated band for the fragmentation share of HDMF in a
pentose/glycine system at 90 °C: 15−44 %, centred near 30 %.**

### 3.4 EDGE C — established structurally, unmeasured in magnitude

Shu & Ho 1988 feed **0.1 M DMHF + 0.1 M cysteine·HCl** and read out 57 products at three pH values.
**`2,5-dimethyl-4-hydroxy-3(2H)-thiophenone` — DMHF with S in place of the ring O — is an
unambiguous, structurally specific tracer of the sink**, and it behaves as a clean gate:

| marker | pH 2.2 | pH 5.1 | pH 7.1 |
|---|---|---|---|
| **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** (GC area %) | **6.0** | **5.8** | **NOT DETECTED** |
| **2,4-hexanedione** | 4.7 | 6.7 | **NOT DETECTED** |
| the two minor thiophenones | T | T | **not detected** |
| **total volatiles from the reaction** (mg, digitised) | ≈178 | **≈321 (peak)** | ≈247 |

**⚠️ The pH 7.1 zeros are ambiguous by the authors' own reading** — they argue *"secondary reactions
occurred readily at this pH"* and that the primary products were consumed rather than never formed.
**The paper cannot separate "the sink is off" from "the sink ran and its products were eaten."**

**The balanced edge this licenses:**
```
DMHF + H2S -> 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone + H2O
C6H8O3 + H2S -> C6H8O2S + H2O                                   [balances exactly]
```
**This is the exact C6 counterpart of the repo's existing C5 edge** `furanone_reductive_sulfhydrylation`
(`src/barrier_constants.py:324`, norfuraneol + H₂S + 2[H] → 2-methyl-3-furanthiol + 2 H₂O). **The
repo has the C5 thiol chemistry and no C6 counterpart at all: `DMHF` appears in
`src/reaction_templates.py` exactly twice, both times as `products=[dmhf, ...]`. DMHF is never a
reactant.**

### 3.5 ★★ THE FORMATION/SINK CROSS-VALIDATION — the most satisfying result in the wave

| paper | side | measurement | pH direction |
|---|---|---|---|
| **Shu & Ho 1988** | **SINK** | DMHF+cysteine coupling markers at 6.0 % (pH 2.2), 5.8 % (pH 5.1), **zero (pH 7.1)** | **sink is acid/neutral-favoured, OFF by pH 7** |
| **Wang & Ho 2008** | **FORMATION (net)** | DMHF from MG+cysteine: **20 → 34 → 70 mg/mol MG** at pH 3 → 5 → 8 | **net DMHF RISES with pH** |

**These are the same phenomenon seen from opposite ends.** If the cysteine sink is acid-favoured and
extinguishes by pH 7−8, then the **net** DMHF surviving in a cysteine system must rise with pH —
which is precisely Wang & Ho's measured curve. **Two experiments, two decades apart, same Rutgers
group (Wang & Ho cite Shu & Ho 1988 as their ref 10), one from each side of the node, giving the
same pH dependence with the same sign.** ★ A model that reproduces one and not the other has the
sink or the formation edge mis-signed, and the pair will say which.

---

# §4. ★★ THE pH PROBLEM — FOUR LADDERS THAT DO NOT AGREE, AND WHY THAT IS INFORMATION

The brief asked whether the quantitative papers corroborate or contradict. **On pH they appear to
contradict, and the resolution is that pH acts at four different places.**

| # | ladder | system | direction | magnitude | source |
|---|---|---|---|---|---|
| **L1** | Edge A, pentose + **glycine** | xylose/Gly, 0.2 M PO₄, 90 °C | **flat, then weakly UP** | 2.6 → 2.6 → 3.1 µg/mmol at pH 5/6/7 = **1.19×** | Blank 1997 T4 |
| **L2** | Edge A, pentose + **alanine** | xylose/Ala, same | **strongly UP** | HDMF 0.3 → 0.9 → 2.5 = **8.3×**; EHMF 2.0 → 7.5 → 13.5 = **6.75×** | Blank 1997 T4 |
| **L3** | Edge B, **MG alone** and **MG + cysteine** | 1.4 M MG, 0.5 M PO₄, 120 °C | **UP** | MG alone **49×**; MG+Cys **3.5×** over pH 3→8 | Wang & Ho F1 |
| **L4** | Edge B, **MG + glycine** | same | **★ DOWN** | 72 → 18 → 8 mg/mol = **9× down** over pH 3→8 | Wang & Ho F1 |

### 4.1 ★★ THE RESOLUTION — pH acts at four distinct places, and they need four terms

1. **The enolisation fork (base-favoured).** Higher pH shifts Amadori decomposition from
   1,2-enolisation (→ 3-deoxyosone → **furfural/HMF**) toward 2,3-enolisation (→ 1-deoxyosone →
   **furanones**). **Two papers name this as their own explanation, ten years apart, on different
   products, in the same direction:** Blank 1997 p. 2647 (citing Beck 1988's 1-DG:3-DG = **20:1 at
   pH 7 vs 8:5 at pH 4.5**), and Apriyantono & Ames 1993 p. 478 (*"Degradation of the ARP via
   1,2-enolisation is favoured in the pH 2.6 system, and 2-furfural is the main compound formed from
   a pentose ARP by this route"*). ★ **Apriyantono then MEASURES it: furfural 274× UP and the
   furanone + all four N-heterocycle classes gone, when the pH is allowed to fall from 4.9 to 2.6.**
   This is the cluster's clean, quantified enolisation-fork evidence. **It explains L1, L2 and L3 —
   everything base-favoured.**
2. **The Cannizzaro step (base-catalysed).** Edge B's reductant, acetol, is made by a base-catalysed
   disproportionation of MG. Wang & Ho: *"DMHF formation from MG was pH-dependent because the
   Cannizzaro reaction is a base preferential reaction."* **Explains the steepest ladder, L3's 49×.**
3. **Strecker competition for the precursor (base-favoured, but SUBTRACTIVE).** At high pH, Strecker
   degradation consumes MG before it can dimerise. **This is Wang & Ho's own explanation for L4** —
   and note that **the authors decline to explain L4 fully** (*"Further studies are required to
   explain this phenomenon"*). ★ **L4 is a NET curve on Edge B in the presence of a competitor, not
   a curve on the edge itself.**
4. **The cysteine sink (acid/neutral-favoured).** §3.4/§3.5.

**★★ THE MODELLING CONSEQUENCE, STATED PLAINLY: L1 (1.19×) AND L2 (8.3×) COME FROM THE SAME
REACTION FAMILY, THE SAME SUGAR, THE SAME BUFFER AND THE SAME LAB, DIFFERING ONLY IN THE AMINO
ACID.** A single pH factor attached to a `Furanone_Formation` family **cannot** produce a 1.2× and an
8.3× response. **That is the sharpest falsification target the wave produced, and it does not need
an absolute scale — it is a contrast of two slopes.**

### 4.2 ⚠️ AND EVERY ONE OF THESE LADDERS IS AN INITIAL-pH LADDER

| paper | buffer | pH labels are | final pH reported? |
|---|---|---|---|
| Blank & Fay 1996 | 0.2 M Na₂HPO₄ | initial | **YES — 6.0 → 5.0 (Gly) / 5.3 (Ala)** ← the only one |
| **Blank 1997** | 0.2 M phosphate / malonate / **water** | **set-point; hold NOT stated** | **NO** |
| **Wang & Ho 2008** | 0.5 M phosphate, adjusted with 1 N NaOH | **set-point; hold NOT stated** | **NO** |
| **Poisson 2019** | none (a dry roast) | **initial BRE pH 5.5** | **N/A** |
| **Shu & Ho 1988** | **★ NONE — unbuffered, titrated with 10 % Na₂CO₃** | **initial, room temperature** | **NO** |
| **Apriyantono & Ames 1993** | none — **actively titrated with 3 M NaOH throughout** | **★ HELD at 5.0 in one arm, DRIFTING 4.9 → 2.6 in the other** | **★ YES, both arms** |

**★ Shu & Ho 1988 is the worst case in the corpus: an unbuffered 0.1 M cysteine·HCl solution
titrated to pH 7.1 with sodium carbonate and then sealed at 160 °C for 30 minutes while generating
H₂S, CO₂ and NH₃.** Its "pH 7.1" cannot have survived the heat-up. **Every use of a Shu & Ho pH
label must carry `initial, unbuffered, final unknown`.**

**★★ Apriyantono & Ames 1993 is the ONLY paper in the corpus that runs the held-vs-drifting
comparison.** That is why it gets a named role of its own — §7.1.

---

# §5. ★★ THE NAMING TRAP — WHAT "HMF"/"HDMF"/"DMHF" DENOTES, PER PAPER

The brief asked for this per paper, from each paper's own abbreviation table. **I checked each
one explicitly. Only ONE of the six papers has a formal abbreviations section, and it does not list
any furanone.**

| paper | has an abbreviation table? | **`HMF`-shaped token** | **what it denotes** | furaneol is written | norfuraneol is written |
|---|---|---|---|---|---|
| **`blank1996`** | no | **`HMF (3)`** | **★ NORFURANEOL** ← trap | `HDMF` | **`HMF`** |
| **`blank1997`** | no | **never used** | — | `HDMF` | spelled out in full |
| **`wang2008`** | no | **body: never used**; **reference title: `5-(hydroxymethyl)furfural`** | **5-HMF** (in a cited title only) | `DMHF` | **not mentioned at all** |
| **`poisson2019`** | **★ YES** — `EB, BRE, RB, CAMOLA` **and nothing else** | **never used** | — | `HDMF` (defined inline only) | spelled out in full |
| **`apriyantono1993`** | no | **★★ `HMFone`** | **★★ NORFURANEOL** ← trap, defined on p. 481, **four pages after the compound first appears in Table 1** | spelled out, **cited only** | **`HMFone`** |
| **`shu1988`** | no | **never used** | — | **`DMHF`** | in one reference title only |

**★★ TWO PAPERS IN THE CORPUS USE AN "HMF"-SHAPED TOKEN TO MEAN NORFURANEOL** (`blank1996`'s
`HMF (3)` and `apriyantono1993`'s `HMFone`). **Neither means 5-hydroxymethylfurfural.**
**Apriyantono is the more dangerous of the two**, because the same paper separately mentions
**5-hydroxymethyl-2-furfural** (the real 5-HMF, cited to Tressl 1989) *and*
**2,5-dimethyl-4-hydroxy-3(2H)-furanone** (furaneol, cited to Shu & Ho 1988) — **three confusably
named molecules in one paper, only one of which (norfuraneol) it actually measures.**

**Letter-order warning:** furaneol is written **`HDMF`** by Blank 1996, Blank 1997 and Poisson 2019,
and **`DMHF`** by Wang & Ho 2008 and Shu & Ho 1988. **The repo uses `DMHF`.** And within the Blank
pair, homofuraneol is **`HEMF`** in 1996 and **`EHMF`** in 1997 — **same lab, same first author, one
year apart, same molecule.**

---

# §6. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured · **[M-fig]/[D-fig]** published only as a figure / my digitisation of it ·
**[C]** cited from elsewhere · **[D]** derived arithmetic · **[P]** proposed by authors.
Full transcriptions, conditions and caveats live in the per-paper dossiers; this is the index.

## 6a. EDGE A — formation from an intact sugar skeleton via acetylformoin

| # | quantity | value | conditions | class | source |
|---|---|---|---|---|---|
| A1 | **HDMF, arabinose / xylose / ribose + glycine** | **5.1 / 2.6 / 3.6 µg per mmol sugar** (0.0040 / 0.0020 / 0.0028 mol %) | 90 °C, 1 h, 0.2 M PO₄, pH 6, 1:1 | **[M]** SIDA, **SD ≤ 10 %** | Blank 1997 T1 |
| A2 | **HDMF, same three + L-alanine** | **1.2 / 0.9 / 1.6 µg/mmol** | same | **[M]** | Blank 1997 T1 |
| A3 | **EHMF (=HEMF), three sugars + glycine** | **1.3 / 0.3 / 0.7 µg/mmol** | same | **[M]** | Blank 1997 T1 |
| A4 | **EHMF, three sugars + L-alanine** | **6.8 / 7.5 / 10.0 µg/mmol** (0.0048/0.0053/0.0070 mol %) | same | **[M]** | Blank 1997 T1 |
| A5 | **Both furanones, pentose ALONE** | **< 0.01 µg/mmol**, phosphate **and** water | same | **[M]** upper bound | Blank 1997 T1 fn c, T3 |
| A6 | **Sugar ordering, HDMF** | **arabinose > ribose > xylose**, 1.96× | Gly, pH 6 | **[M]** | Blank 1997 |
| A7 | **Sugar ordering, EHMF** | **ribose > xylose > arabinose**, 1.47× — **DIFFERENT from A6** | Ala, pH 6 | **[M]** | Blank 1997 |
| A8 | **Selectivity swing (xylose)** | HDMF/EHMF **8.7 (Gly) → 0.12 (Ala) = 72×** | pH 6 | **[D]** | Blank 1997 |
| A9 | **★ Strecker / fragmentation split, HDMF** | **73 / 27 %** (GABA control) | 90 °C, pH 6 | **[D]** from **[M]** | Blank 1997 T2 |
| A10 | **Strecker / fragmentation split, EHMF** | **97 / 3 %** | same | **[D]** from **[M]** | Blank 1997 T2 |
| A11 | **Isotopomer branch fraction, xylose+Gly** | **~70 % labelled / ~30 % fragmentation**; double-label **4/12/84 %** | 90 °C, 1 h | **[M]** | Blank & Fay 1996 |
| A12 | **Isotopomer branch fraction, xylose+Ala** | **15 / 47 / 38 %** | same | **[M]** | Blank & Fay 1996 |
| A13 | **★ pH ladder, xylose/glycine** | HDMF **2.6 / 2.6 / 3.1**; EHMF **0.3 / 0.3 / 0.7** at pH **5/6/7** | 0.2 M PO₄ | **[M]** | Blank 1997 T4 |
| A14 | **★ pH ladder, xylose/L-alanine** | HDMF **0.3 / 0.9 / 2.5**; EHMF **2.0 / 7.5 / 13.5** | same | **[M]** | Blank 1997 T4 |
| A15 | **★★ Amino-acid-specific pH sensitivity** | pH5→7: **Gly 1.19× vs Ala 8.3×** (HDMF) | same | **[D]** | Blank 1997 |
| A16 | **★★ Apparent reaction order in amino acid** | **n = 0.35 (HDMF/Gly), 0.37, 0.42, 0.71 (EHMF/Ala)** — **all sub-linear** | 1:1/1:2/1:4, 90 °C | **[D]** from **[M]** | Blank 1997 T5 |
| A17 | **Phosphate catalysis vs 0.2 M malonate** (clean) | **3−9.4×** | pH 6, both buffered | **[D]** | Blank 1997 T3 |
| A18 | Phosphate catalysis vs water (confounded) | **43−150×** ⚠️ water arm pH-uncontrolled | pH 6 nominal | **[D]** | Blank 1997 T3 |
| A19 | **★★ Hexose DMHF is 100 % intact-C6** | `[¹³C₆]:[¹²C₆] = 1:1`, **`[¹³C₁]`−`[¹³C₅]` ABSENT** | 120 °C, pH 5, 60 min, +Gly or Cys | **[M]** | Wang & Ho §4.1 |
| A20 | **★★ Acetylformoin is ALSO intact-C6** | m/z **144 : 150 = 1:1**, equal intensity | same | **[M]** | Wang & Ho §4.1 |
| A21 | **★★ Sucrose → HDMF in a real bean is intact-skeleton** | intact = **87.3−100 %** of the sucrose share at **all nine** roast times | 171−220 °C, 50−400 s, pH 5.5 init. | **[M]** | Poisson T9/T10 |
| A22 | **★★ POSITIVE CONTROL for A21** | 2,3-butanedione's intact share **collapses 25.4 → 0.4 %** and 2,3-pentanedione's C2/C3 **rises 4.8 → 18.8 %** in the SAME runs | same | **[M]** | Poisson T5/T8 |
| A23 | No glucose fragmentation even at **165 °C** | qualitative, **"data not shown"** | — | **[M]** | Wang & Ho |
| A24 | **norfuraneol ≫ DMHF** | ordinal, **"main reaction product … (data not shown)"** | 90 °C, all six systems | **[M]** ordinal | Blank 1997 p. 2646 |
| A25 | 1-DG : 3-DG = **20:1 (pH 7)**, **8:5 (pH 4.5)** | ratio, **hexose** | — | **[C]** Beck 1988 | via Blank 1997 |
| A26 | **rate constant / Ea / time course** | **NONE anywhere** | — | — | all five papers |

## 6b. EDGE B — formation from methylglyoxal (C3 + C3), acetylformoin-free

| # | quantity | value | conditions | class | source |
|---|---|---|---|---|---|
| B1 | **DMHF from MG alone** | **4 ± 2 / 13 ± 1 / 199 ± 19 mg per mol MG** at pH 3/5/8 | 1.4 M MG, 0.5 M PO₄, 120 °C, 60 min, n=3 | **[D-fig]** | Wang & Ho F1 |
| B2 | **DMHF from MG + cysteine** | **20 ± 2 / 34 / 70 ± 16 mg/mol MG** | same +1.0 M Cys | **[D-fig]** | Wang & Ho F1 |
| B3 | **DMHF from MG + glycine** | **72 ± 21 / 18 ± 3 / 8 mg/mol MG** | same +1.0 M Gly | **[D-fig]** | Wang & Ho F1 |
| B4 | **pH trends** | MG alone **↑49×**; MG+Cys **↑3.5×**; **MG+Gly ↓9×** | pH 3→8 | **[M]** sign / **[D-fig]** magnitude | Wang & Ho |
| B5 | **★★ Amino-acid effect swing across pH** | glycine **×18 at pH 3, ÷26 at pH 8** = a **460× swing** | — | **[D]** | Wang & Ho |
| B6 | **★★ MG-spike branch fraction** | **`[¹³C₆]`:`[¹²C₆]` = 4:1 → 20 % from MG**, 80 % from glucose. **A LOWER BOUND.** | 1.4 M glucose + 1.4 M MG, 120 °C, pH 5 | **[M]** | Wang & Ho §4.2 |
| B7 | **★★ Edge B produces NO acetylformoin** | `[¹²C₆]`acetylformoin **ABSENT** in the MG-spiked run | same | **[M]** | Wang & Ho §4.2 |
| B8 | **★ In-situ Edge B share, aqueous hexose** | **below detection** (no `[¹³C₃]`DMHF) | 120 °C, pH 5, no MG spike | **[M]** ceiling | Wang & Ho §4.1 |
| B9 | **★★ Edge B share in a real matrix** | **8.2 / 12.7 / 11.4 %** of sucrose-derived HDMF at **210/260/300 s only**; **< 0.5 % floor elsewhere** | coffee bean, 199−204 °C | **[M]** | Poisson T10 |
| B10 | **The edges do not mix** | only `[¹³C₆]` and `[¹²C₆]` isotopomers ever seen — **no `[¹³C₃]`DMHF** | both CAMOLA runs | **[M]** | Wang & Ho |
| B11 | **Route** | 2 MG →(Cannizzaro, base-cat.)→ MG + acetol →(aldol)→ 2,5-dioxo-3,4-dihydroxyhexane → DMHF | — | **[P]**/**[C]** Schieberle 2005 | Wang & Ho F2 |
| B12 | Blank's lab saw the same edge in 1997 | *"some of the well-known C3 fragments do generate HDMF, e.g. **methylglyoxal, acetol, and dihydroxyacetone** (data not shown)"* | 90 °C, pentose systems | **[M]**, unpublished | Blank 1997 p. 2646 |

## 6c. EDGE C — consumption by cysteine / H₂S

| # | quantity | value | conditions | class | source |
|---|---|---|---|---|---|
| C1 | **★★ 2,5-dimethyl-4-hydroxy-3(2H)-THIOPHENONE** (the DMHF O→S swap) | **6.0 % / 5.8 % / NOT DETECTED** GC area, pH 2.2/5.1/7.1 | 0.1 M DMHF + 0.1 M Cys·HCl, **UNBUFFERED**, 160 °C, 0.5 h, Parr bomb | **[M]**, area % | Shu & Ho T1 |
| C2 | **2,4-hexanedione** — second marker, same gate | **4.7 / 6.7 / NOT DETECTED** | same | **[M]** | Shu & Ho T1 |
| C3 | Two minor thiophenones | **T / T / not detected** | same | **[M]** | Shu & Ho T1 |
| C4 | **3,5-dimethyl-1,2,4-trithiolane** | **3.1 / 24.7 / 3.2 %** — **peaked 8× at cysteine's pI** | same | **[M]** | Shu & Ho T1 |
| C5 | **★ Pyrazines (9 species)** | **ALL absent at pH 2.2 and 5.1; ALL present at 7.1** | same | **[M]** on/off | Shu & Ho T1 |
| C6 | **Total volatiles, DMHF + cysteine** | **≈178 / 321 / 247 mg** — **PEAKED at pH 5.1** | same | **[D-fig]** | Shu & Ho F1 |
| C7 | Total volatiles, **DMHF alone** | **≈160 / 122 / 93 mg** — monotone, acid-favoured | — | **[C]** Shu 1985b, **[D-fig]** | Shu & Ho F1 |
| C8 | Total volatiles, **cysteine alone** | **≈57 / 202 / 25 mg** — peaked 8× at pI | — | **[C]** Shu 1985a, **[D-fig]** | Shu & Ho F1 |
| C9 | **★ Cysteine, not DMHF, controls the sink's pH shape** | qualitative — the reaction curve resembles C8, not C7 | — | **[M]** interpretation | Shu & Ho p. 802 |
| C10 | **★ Cysteine is BOTH a reductant feeding Edge B AND a sink consuming DMHF; which wins is pH-set, crossover between pH 3 and 5** | qualitative | 120 °C | **[P]** authors, **[C]** Haleva-Toledo 1999 | Wang & Ho §5.3 |
| C11 | **DMHF consumed / conversion / mass balance / rate** | **NONE** | — | — | **ABSENT** |
| C12 | ⚠️ pH 7.1 zeros are ambiguous | authors argue **secondary consumption**, not absence of coupling | — | **[P]** | Shu & Ho p. 801 |

## 6d. THE pH-TRAJECTORY PAIR and the surrounding context

| # | quantity | value | conditions | class | source |
|---|---|---|---|---|---|
| P1 | **★★ Total volatiles, held pH 5.0 vs drifting 4.9→2.6** | **2.47 × 10⁴ vs 3.53 × 10⁶ nmol/mol xylose** — **143× UP on drift** | 1 M xylose + 1 M lysine·HCl, reflux 1 h, SDE | **[M]**, **±16 %** | Apriyantono T1 |
| P2 | **★★ 2-Furfural, same pair** | **12.9 × 10³ vs 3.53 × 10⁶ nmol/mol = 0.00129 vs 0.353 mol %** — **274×** | same | **[M]** | Apriyantono T1 |
| P3 | **★ N-containing volatiles (16 cpds)** | **9.46 vs 0.127 µmol/mol** — **75× DOWN on drift** | same | **[D]** from **[M]** | Apriyantono T1 |
| P4 | **★ Pyrazines** | **8940 vs 21.6 nmol/mol** — **414× down**; pyrazine itself **on/off** | same | **[M]** | Apriyantono T1 |
| P5 | **★ Four classes go present → NOT DETECTED** | monocyclic pyrroles (24.8), pyridines (66.2), pyrrolizines (122), 2-furanmethanol (492) | same | **[M]** | Apriyantono T1 |
| P6 | **Norfuraneol** | **`tr` (< 2 nmol/mol) held; NOT DETECTED drifting** | same | **[M]**, at floor | Apriyantono T1 |
| P7 | **⚠️ NO DMHF anywhere in that paper** | not among 58 identified or 14 unidentified compounds | — | — | Apriyantono |
| P8 | **★ The enolisation fork, MEASURED** | acid → 1,2-enolisation → furfural; neutral → 2,3-enolisation → furanones + N-heterocycles | — | **[M]** interpretation, **[C]** Nursten 1986 | Apriyantono p. 478 |
| P9 | **★ Matrix realism: sugar : glycine ≈ 620 : 1 in a coffee bean** | sucrose 4593 mg vs glycine 1.62 mg, alanine 19.4 mg, per 50 g bean | — | **[M]** | Poisson T1 |
| P10 | Chlorogenic acid = the largest BRE component | **3229 mg / 50 g bean**, 1.3× the total sucrose | — | **[M]** | Poisson T1 |
| P11 | Bean moisture at roast start | **10 ± 0.5 %** | — | **[M]** | Poisson |

---

# §7. ★★ IS CALIBRATION NOW POSSIBLE? — **yes / no / partially, per edge**

| edge | **magnitude** | **selectivity / branching** | **pH response** | **rate (k, Ea)** | **verdict** |
|---|---|---|---|---|---|
| **A — pentose, via Strecker elongation + acetylformoin** | **★★ YES.** 39 SIDA cells, µg/mmol sugar, σ ≤ 10 %, internal standard added pre-workup (A1−A5) | **★★ YES.** Two-method 73/27 and 70/30 agreement (A9, A11); sugar and amino-acid resolved (A6−A8) | **★★ YES**, and amino-acid-resolved (A13−A15) | **❌ NO** — one T, one t | **★★ CALIBRATABLE.** |
| **A — hexose, intact C6** | **⚠️ PARTIALLY.** No absolute hexose DMHF yield exists in any of the five papers. Structure is settled (A19−A22); level is not | **★★ YES**, structurally: 100 % intact-skeleton, with an internal positive control (A21/A22) | ⚠️ only via Wang & Ho's MG-containing systems, which confound Edge B | **❌ NO** | **PARTIALLY — structure yes, magnitude no.** |
| **B — from MGO, C3+C3, acetylformoin-free** | **⚠️ PARTIALLY.** Nine cells, but **digitised from a bar chart**, external-standard HPLC, no recovery (B1−B3) | **★★ YES**, and bracketed from both sides by three independent measurements (B6, B8, B9) | **★★ YES**, and sign-reversing (B4, B5) | **❌ NO** — one T, one t | **PARTIALLY — the branch is well bounded, the level is the weakest anchor in the cluster.** |
| **C — sink, by cysteine / H₂S** | **❌ NO.** No residual DMHF, no conversion, no mass balance, no molar yield of any product. Table I is GC area % with no internal standard | **✔ YES** structurally: an unambiguous O→S tracer plus a second marker (C1−C3) | **✔ YES** as a shape (C1, C6) — ⚠️ but unbuffered, initial-pH only, and the pH 7.1 zeros are ambiguous (C12) | **❌ NO** | **❌ NOT CALIBRATABLE. Structure only.** |
| **Any barrier / activation energy, any edge** | — | — | — | **❌ NO** | **❌ IMPOSSIBLE FROM THIS LITERATURE.** **All five papers are single-temperature.** Round-3 §B.3's verdict stands. |

### 7.1 ★ THE HONEST SUMMARY

**What this wave DOES enable:** calibrating the pentose formation edge's magnitude and its five-way
response surface against real SIDA numbers with a stated σ; bounding the MGO edge's branch share
from both sides; pinning three structural results (intact-C6, acetylformoin-on-A-not-on-B,
sink-exists) as hold-outs; and validating the model's **pH-trajectory** behaviour against the only
held-vs-drifting pair in the corpus.

**What it does NOT enable, and cannot be made to:** any rate constant, any activation energy, any
barrier fit, any hexose DMHF magnitude, and any sink magnitude.

**⚠️ THE PROHIBITION THAT MATTERS MOST.** Every one of the five papers is **single-temperature**
(90 / 120 / 160 °C / reflux / a confounded roast ramp). **A barrier fitted to a single-temperature
dataset is a rate-constant fit wearing an Arrhenius costume: it will absorb every unmodelled factor
in the system and will not transfer.** This is precisely the failure mode
`src/barrier_constants.py:307` documents at length for `thiol_addition_pentodiulose` — a constant
fitted to a target that turned out not to be a measurement, then reverted in Wave S2c. **The
difference here is that these ARE measurements; the defect would be the single-temperature design,
not fabrication. That still does not make them barrier-grade.** Record as a **prohibited
derivation**, by name, before anyone builds these edges.

---

# §8. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION. NOTHING WAS WRITTEN TO ANY REGISTRY, BENCHMARK FILE, CONFIG OR THE DECLARATION.**

## 8.1 ★★★ THE NAMED ROLE THE BRIEF ASKED FOR

**Name:** `apriyantono1993_xylose_lysine_pH_trajectory_pair`
**Role:** **HOLD-OUT — pH-TRAJECTORY VALIDATION PAIR.** Scored as **ONE paired test**, never as two
independent rows.

**Why it needs a named role of its own rather than a slot in a general pH bucket:**

1. **It is the ONLY held-vs-drifting pair in the corpus.** §4.2. Four of the five papers in this
   wave report pH as an unheld set-point with no final pH; `k3` §C.11 records the same absence for
   Zhou 2023. **This paper runs both states, in one lab, one apparatus, one hour, one sugar, one
   amine, from one initial pH (4.9).**
2. **It tests a different thing from a pH ladder.** Blank 1997 Table 4 (pH 5/6/7, all buffered) asks
   *"does the model get the pH-5 rate right, and the pH-7 rate right?"* This pair asks *"does the
   model's internal pH STATE evolve correctly, and does the chemistry follow it?"* **A model can
   pass every point of a buffered ladder and still fail this**, if it treats pH as a fixed input
   rather than tracking acid generation.
3. **Four independent scoreable channels in one experiment**, so it is hard to pass by luck:
   **(i)** total volatiles **143×** up on drift (P1); **(ii)** furfural specifically **274×** up
   (P2); **(iii)** N-heterocycles **75×** down with **four whole classes going to not-detected**
   (P3−P5); **(iv)** the only furanone present goes to zero (P6).
4. **Its σ is stated (±16 %)** — unlike Poisson 2019, which states none.
5. **It is mechanistically anchored** to the same 1,2-vs-2,3-enolisation fork that Blank 1997 and
   Beck 1988 invoke from the other direction (§4.1 item 1).

**Suggested scoring shape (proposal only):** score the **log-ratio between the two arms**, per
channel — not the two levels. The ratio form is immune to the single-alkane-internal-standard
weakness and to the SDE recovery bias, because **both arms went through the identical isolation.**
Passing requires getting the direction and order of magnitude of all four right simultaneously.

**⚠️ The caveat that must ride with the role:** the held arm received **unreported 3 M NaOH
throughout the hour**, so it is not sodium- or volume-matched to the drifting arm.

## 8.2 Proposed **FIT** set

| # | candidate | edge | why it is fit-grade | σ |
|---|---|---|---|---|
| **F1** | **Blank 1997 Table 1 — the six sugar × amino-acid HDMF/EHMF cells (A1−A4)** | A (pentose) | The only absolute furanone anchors in the corpus. SIDA, internal standard added pre-workup, calibration r² = 0.999 over 3−50 µg/mL. | **10 %** (stated) |
| **F2** | **The 73/27 and 97/3 Strecker/fragmentation splits (A9, A10)**, jointly with `blank1996` A11/A12 | A | Two independent methods agreeing to 3 pp (§3.3). **Ratio-form, so immune to absolute-scale error.** The safest fit target in the cluster. | ratio |
| **F3** | **Phosphate/malonate catalysis = 3−9.4× (A17)** — **ONLY if a buffer term is ever added** | A | Both arms buffered at the same pH. **Do NOT use A18 (the water comparison): its pH was uncontrolled.** | 10 % on each arm |
| **—** | **NO barrier or Ea, on any edge, from any of these five papers** | all | **§7.1. Single-temperature. Record as a prohibited derivation.** | — |

## 8.3 Proposed **HOLD-OUT** set

Ordered by how much information each carries per unit of implementation risk.

| # | candidate | edge | proposed role | why |
|---|---|---|---|---|
| **H1** | **★★ Poisson's PAIRED structural test: HDMF intact-skeleton ≥ 87 % at all nine roast times WHILE 2,3-butanedione's intact share collapses 25.4 → 0.4 % in the same runs (A21 + A22)** | A | **HOLD-OUT, structural, PAIRED. Score both halves or neither.** | **The best-designed test available anywhere in this cluster.** The pairing is the point: a model with no fragmentation chemistry at all would pass the HDMF half and fail the diketone half. Scoring both makes it discriminating instead of vacuous. |
| **H2** | **★★ Blank 1997's amino-acid-specific pH sensitivity: glycine 1.19× vs alanine 8.3× over pH 5→7 (A15)** | A | **HOLD-OUT, contrast of two slopes.** | Survives any absolute-scale error. **No single pH factor on a shared reaction family can produce both.** §4.1. |
| **H3** | **★★ `apriyantono1993_xylose_lysine_pH_trajectory_pair` (P1−P6)** | — | **HOLD-OUT, pH-trajectory. See §8.1.** | The corpus's only held-vs-drifting control. |
| **H4** | **★★ Wang & Ho's intact-C6 CAMOLA at BOTH the product and the intermediate (A19 + A20)** | A | **HOLD-OUT, structural.** **Re-declare against the PRIMARY source** — `src/reaction_templates.py:440−442` currently cites it from an abstract. | Isotopomer ratios; response-factor-immune; already implemented, so it is a free regression. |
| **H5** | **★★ Edge B is acetylformoin-free (B7)** | B | **HOLD-OUT, structural NEGATIVE.** If an MGO edge is added it must **not** route through acetylformoin. | A negative structural constraint is cheap to test and expensive to discover late. It also resolves an ambiguity Poisson explicitly leaves open (§3.1). |
| **H6** | **★★ The Edge B bracket: below detection in-situ (B8), 8−13 % in a real matrix (B9), 20 % at a 1.4 M MG spike (B6)** | B | **HOLD-OUT, two-sided bracket.** | Bounds the edge from both directions without needing an absolute yield. Three independent measurements, two labs, two matrices, a 100 °C spread. |
| **H7** | **★★ The formation/sink pH cross-validation: Shu's sink markers OFF by pH 7.1 (C1, C2) AND Wang & Ho's net DMHF RISING with pH in a cysteine system (B2)** | B + C | **HOLD-OUT, PAIRED, cross-paper.** | §3.5. Two experiments from opposite sides of the node giving the same pH dependence. **A model that reproduces one and not the other has the sink or the formation edge mis-signed, and the pair will say which.** ⚠️ Only meaningful if Edge C is implemented. |
| **H8** | **Blank 1997's sub-linear reaction order in amino acid, n = 0.35−0.71 (A16)** | A | **HOLD-OUT — and it is EXPECTED TO FAIL TODAY.** | `_furanone_generation` (`src/reaction_templates.py:2218−2274`) emits a **bimolecular** step, order 1 in amino acid **by construction**. Declaring this records a **known, diagnosed** structural defect instead of discovering it later. |
| **H9** | **Pentose alone < 0.01 µg/mmol, both furanones, two media (A5)** | A | **HOLD-OUT, structural CEILING.** ⚠️ **REPLACES `blank1996` §10.2 #8, which proposed the opposite sign.** §2. | Cleanest negative control in the cluster. |
| **H10** | **The 72× HDMF/EHMF selectivity swing on xylose (A8)** and **the 460× glycine-effect swing across pH (B5)** | A, B | **HOLD-OUT, ratios.** | Scale-free; two-sided directional tests. |
| **H11** | **The two DIFFERENT sugar orderings (A6, A7)** | A | **HOLD-OUT, ordinal, hard.** | arabinose > ribose > xylose for HDMF **and** ribose > xylose > arabinose for EHMF, simultaneously. A single "pentose reactivity" scalar cannot satisfy both. |
| **H12** | **norfuraneol ≫ DMHF and ≫ HEMF (A24 + `blank1996` #9)** | A | **HOLD-OUT, ORDINAL ONLY.** | Now supported by two papers and six systems. **Still unquantified in both — the hold-out can only ever be `norfuraneol > DMHF`, never a ratio.** ⚠️ **Do NOT score Apriyantono against it** (P6 is a method-limited near-null in a lysine system; §8.5). |
| **H13** | **Pyrazines on/off at pH 7 — Shu & Ho (C5) scored JOINTLY with Apriyantono (P4)** | — | **HOLD-OUT, structural on/off, cross-paper.** | Two systems, two labs, two precursor sets, same direction and threshold. **A joint hold-out across two papers is far harder to pass by accident than either alone.** |
| **H14** | **Shu & Ho's sink markers: thiophenone + 2,4-hexanedione present at pH 2.2 AND 5.1, absent at 7.1 (C1, C2)** | C | **HOLD-OUT, structural on/off, PAIRED.** Only if Edge C is implemented. | Two independent structural markers of one gate. ⚠️ Must carry C12: the pH 7.1 absence is, on the authors' own reading, **secondary consumption**, so this tests *net survival*, not *coupling rate*. |
| **H15** | **Poisson's BRE stoichiometry, sugar : glycine ≈ 620 : 1 (P9)** | — | **HOLD-OUT of a different kind: a CONDITION-VALIDITY GATE.** | Not a result to score, but a **precondition test**: if a coffee-lane prediction uses branch fractions derived at Blank's 1:1 stoichiometry, that is a category error and this number detects it. |

## 8.4 ★★ THE INDEPENDENCE AND LEAKAGE WARNINGS — read before splitting

1. **Blank 1997's five tables are NOT 39 independent measurements.** The conditions
   `xylose/glycine, pH 6, phosphate, 1:1` (HDMF 2.6 / EHMF 0.3) and
   `xylose/L-alanine, pH 6, phosphate, 1:1` (HDMF 0.9 / EHMF 7.5) **each appear in FOUR of the five
   tables.** **Any split putting one instance in the fit set and another in the hold-out set has
   leaked.** The safe unit of splitting is the **experimental axis**, or the **table**.
2. **F1 and H2 come from the same paper and share the reference cell.** If F1 (Table 1) is used as a
   fit set, **H2 must be scored on the pH SLOPE CONTRAST, not on the pH 6 levels** — the pH 6 column
   of Table 4 *is* Table 1.
3. **F2 and `blank1996`'s A11/A12 are two methods on one system, not two systems.** Fitting a branch
   prior to both is legitimate; treating them as independent evidence of magnitude is not.
4. **H7 and H13 are cross-paper pairs by design.** Do not decompose them; the pairing is what makes
   them discriminating.

## 8.5 ★ EXPLICITLY **NEITHER** — excluded from both roles, with reasons

| item | why excluded |
|---|---|
| **`blank1996` items #12/#13 — the "low mg/kg" rough estimate and my derived 10⁻³−10⁻² mol %** | **SUPERSEDED** by Blank 1997's measured µg/mmol. `blank1996_extraction.md` §3.3 already excluded them; keep them excluded and retire them. *(For the record: the measured range, 7 × 10⁻⁴ − 1.4 × 10⁻² mol %, VALIDATES that dossier's derivation. That is a reason to trust its method, not to ingest its numbers.)* |
| **`blank1996` §10.2 items #7 and #8 AS WRITTEN** | **CONTRADICTED. §2.** #7 would fail a correct model; #8 has the wrong sign. |
| **Apriyantono's norfuraneol cell (P6) against the `norfuraneol ≫ DMHF` ordering** | Both terms are at/below the detection floor; the amine is **lysine**, not glycine or alanine; the isolation is **Likens−Nickerson SDE**, close to a worst case for a water-miscible labile furanone; and **the authors themselves read the trace as consumption into coloured products**, supported by an independent citation. **The paper is SILENT on the ordering, not contradictory.** |
| **Poisson 2019, anything, as a FIT target** | No absolute quantity; **no replicate count and no SD anywhere in the paper**; T and t confounded. Structural hold-out only. |
| **Shu & Ho 1988, anything, as a FIT target** | GC area % with no internal standard; **no replication stated**; **UNBUFFERED with initial-only pH labels**; one T, one t; **no DMHF conversion**. §8.6. |
| **Wang & Ho's nine bar values (B1−B3) as a FIT target** | Digitised from a chart by an agent, external-standard HPLC with no recovery, pH hold unstated — **three transmission defects deep.** Usable as a **loose (≥3×) hold-out band**, not as a contract. |
| **Beck 1988's 20:1 / 8:5 (A25)** | **[C], hexose, unread.** Attribute to Beck, never to Blank 1997. |
| **The Huber 1992 odour thresholds** (`blank1996` §7) | Trade-magazine secondary source. Unchanged. |

## 8.6 ★★ A PROHIBITED-DERIVATION WARNING, BY NAME, BEFORE ANYONE BUILDS EDGE C

Shu & Ho 1988 is a **near-perfect setup for the failure mode `src/barrier_constants.py:307`
records.** The sulfur branch has, on that entry's own words, **"ZERO ABSOLUTE LITERATURE ANCHORS"**,
and this wave has just delivered a fed-precursor sulfur paper that *looks* like it supplies one.
**It does not.** `thiol_addition_pentodiulose` was fitted to `cys_ribose_140C_Hofmann1998`'s
342/200 ppb, which Wave S2b proved was a repo-internal derivation rather than a measurement, and
Wave S2c had to revert it.

**A `thiol_addition_dmhf` (or similarly named) constant fitted to Shu & Ho's `6.0 % GC area` would
repeat that exactly — with the added defect that 6.0 % is not even a concentration.** If Edge C is
built:

- its barrier must be an **ESTIMATE, explicitly labelled UNCONSTRAINED**, in the same register as
  `thiol_addition_pentodiulose` and `deoxyosone_reduction`;
- the natural neutral starting value is the un-fitted **`thiol_addition` class value 28.60** — the
  same value Wave N gave its siblings — **not** a value tuned to make C1 come out at 6 %;
- C1/C2/C5 become **hold-outs on the edge, never fit targets**;
- **and it should be built only if the orchestrator wants it.** Adding a sink with a free rate to a
  channel whose formation edges are only now becoming calibratable puts **two unconstrained
  parameters in series.** `blank1997` calibrates formation; **nothing calibrates this.**
  **Sequencing matters more than completeness here.**

---

# §9. WHAT THE REPO'S CURRENT CODE LOOKS LIKE AGAINST THIS — read-only observations

Recorded for the orchestrator. **Nothing in `src/` was modified.**

| # | observation | severity |
|---|---|---|
| **R1** | **★ `blank1996_extraction.md` defect D2 is now PARTLY OBSOLETE.** That dossier states *"There is no pentose → DMHF step and no pentose → HEMF step anywhere in the network."* **That is no longer true.** `src/reaction_templates.py::_furanone_generation` (lines 2218−2274) emits both `pentose + glycine -> DMHF + NH₃ + CO₂ + 2 H₂O` and `pentose + alanine -> HEMF + NH₃ + CO₂ + 2 H₂O` as `Furanone_Formation`, gated at T ≥ 90 °C, `source_quality="literature"`, `barrier_uncertainty_kcal=6.0`. **I verified both balance exactly** (Gly: C₅H₁₀O₅ + C₂H₅NO₂ → C₆H₈O₃ + NH₃ + CO₂ + 2 H₂O, C/H/O all close). **The channel exists and is scorable against Blank 1997 Table 1 today.** The 1996 dossier appears to have inspected `_furanone_and_mft_route` only. | **Correct the record before it is read as evidence.** |
| **R2** | **★ That same template is FIRST ORDER in amino acid by construction**, and Blank 1997 Table 5 measures **n = 0.35−0.71** (A16, H8). A known structural mis-specification, testable today. | Medium — a shape error, not a wrong route. |
| **R3** | **DMHF is only ever a PRODUCT.** `grep -n "dmhf" src/reaction_templates.py` → two hits, both `products=[dmhf, ...]` (lines ~490 and ~2264). **No DMHF consumption edge, no thiophenone species, no MGO → DMHF edge.** Edge B and Edge C are both absent. | Medium — two missing edges, both now evidenced. |
| **R4** | **`MGO` already exists as a first-class species** (`src/kinetic_core/species.py:77`) and is produced in **all three lanes** — `r_dpo_c2c3`, `r_glc_c2c3` (`src/kinetic_core/sulfur.py:386,392`) and `r_ama_mgo` (`src/kinetic_core/network.py:139`). **Edge B's reactant is already tracked; only the edge is missing.** | ★ An opportunity, not a defect. |
| **R5** | **`src/literature_runtime.py:3386` already claims** *"the lane now carries an explicit methylglyoxal-to-HDMF C3 anchor so fragmentation-heavy browning is not treated as amino-acid-dependent only."* **I could not find a corresponding edge in `reaction_templates.py` or `kinetic_core/`.** Either the anchor lives somewhere I did not look, or the summary string overstates what ships. **Worth checking before Wang & Ho is cited for it.** | **Flag — verify.** |
| **R6** | **`blank1996_extraction.md` defect D1 stands and this wave points at its fix.** `src/barrier_constants.py:325` attributes *"the accepted mechanism names the reductant and it is the amino acid"* partly to Blank & Fay 1996, which does not say it. **Both Wang & Ho 2008 (ref 17) and Poisson 2019 (ref 19) cite the same primary source for the reductant claim: Hofmann & Schieberle 2001, "Acetylformoin — an important progenitor…", Flavour 2000 Proceedings, pp. 311−322.** Fetching that could resolve D1 properly. | Unchanged; **a route to resolution identified.** |
| **R7** | **`src/barrier_constants.py:322` `furanone_reductive_opening` (norfuraneol → 2,3-pentanedione, cited to Whitfield & Mottram 1999) is INDEPENDENTLY CORROBORATED** by Poisson 2019's Figure 2 legend, which names the same route, cites the same Whitfield & Mottram 1999 paper, and adds a second proposed mechanism (Cerny & Davidek 2003). **The repo's route choice is supported by a source it does not currently cite.** | ★ Good news; a citation the repo could add. |
| **R8** | **`k3` §B2.2 already records Shu 1988 as a shape source** (*"three shapes in ONE experiment"*) and `k3` §1243 already classes it *"GC area %, no internal standard; pH-shape source only."* **Both are correct.** `shu1988_extraction.md` adds the numbers (C1−C5) and **Figure 1's three yield curves, which `k3` §B2.2 does not mention.** | ★ Extends, does not contradict. |

---

# §10. DECLARED GAPS FROM THIS WAVE — for the `k3` §C register

> **"No rate constant, no activation energy, no time course and no temperature series exists
> anywhere in the DMHF cluster.** All five papers are single-temperature: Blank 1997 at 90 °C/1 h
> (39 cells), Wang & Ho at 120 °C/60 min (9 cells, plus one unquantified 'data not shown' excursion
> to 165 °C), Shu & Ho at 160 °C/0.5 h, Apriyantono at reflux (**temperature never stated**)/1 h,
> and Poisson on a roast ramp where **temperature and time are confounded by design**. Round-3
> §B.3's verdict stands unchanged: **any Ea the repo assigns to a furanone family is the repo's own
> assumption and must be labelled as such.**"

> **"There is NO absolute hexose DMHF yield in any of the five papers.** The hexose intact-C6
> structure is settled twice over (Wang & Ho CAMOLA; Poisson coffee CAMOLA with an internal positive
> control) and its magnitude is measured **nowhere**. Blank 1997's 39 cells are all **pentose**.
> Wang & Ho's nine cells are all **per mole of MG**. Poisson reports no concentration at all."

> **"There is NO measurement of DMHF CONSUMPTION.** Shu & Ho 1988 establishes the sink structurally
> (an unambiguous O→S ring-swap tracer plus a second marker), gives its complete 57-compound product
> spectrum at three pH values and its total-volatile pH shape — and reports **no residual DMHF, no
> conversion, no mass balance and no molar yield of any product.** A sink edge built on it would
> have a correct structure and a completely free rate constant."

> **"FOUR OF THE FIVE PAPERS REPORT pH AS AN UNHELD SET-POINT WITH NO FINAL pH.** Blank 1997
> (0.2 M phosphate, hold not stated), Wang & Ho (0.5 M phosphate adjusted with NaOH, hold not
> stated), Poisson (BRE pH 5.5 initial, a dry roast) and **Shu & Ho (UNBUFFERED, titrated with 10 %
> Na₂CO₃, then sealed at 160 °C for 30 min while generating H₂S, CO₂ and NH₃ — the worst case in the
> corpus)**. Only Apriyantono & Ames 1993 reports both a held and a drifting trajectory, and only
> Blank & Fay 1996 reports a measured drift (6.0 → 5.0/5.3)."

> **"REPLICATION IS ABSENT FROM TWO OF THE FIVE PAPERS.** Poisson 2019 and Shu & Ho 1988 state **no
> replicate count, no SD and no error bars anywhere** — every value must be treated as a single
> determination of unknown precision, and **no σ may be assigned to them.** The other three do state
> one: Blank 1997 'maximum SD ≤ 10 %, n ≥ 2 assays × 2 injections'; Wang & Ho 'SD, n = 3
> independent analyses'; Apriyantono '±16 %, means of two runs, isolates in at least triplicate'."

> **"The norfuraneol : DMHF ratio does not exist.** Both Blank papers state that norfuraneol is the
> major product ('data not shown' in 1997) and neither quantifies it. Wang & Ho, Poisson and Shu &
> Ho do not measure norfuraneol at all. **The ordinal ceiling `norfuraneol ≫ DMHF` is supported
> twice and quantified zero times.**"

> **"⚠️ TWO PRIOR HOLD-OUT PROPOSALS ARE CONTRADICTED BY THE 1997 ISOTOPE-DILUTION DATA.**
> `blank1996_extraction.md` §10.2 item **#7** ('HEMF requires alanine; a model that emits HEMF from
> pentose+glycine fails') is wrong — Blank 1997 measures EHMF at 0.3−1.3 µg/mmol in all three
> glycine systems. Item **#8** ('DMHF forms from pentose alone — a zero-amino-acid POSITIVE
> control') has the wrong sign — the 1997 control gives < 0.01 µg/mmol for both furanones, in
> phosphate and in water. **A `−` in a 1996 GC-olfactometry table is a non-detection at a sniffing
> port, not a zero.**"

> **"⚠️ `research_round3_channels.md` §C.2's 'RATIO-ONLY' verdict on Apriyantono & Ames 1993 is
> WRONG.** It was written from the abstract. Table 1 reports all 58 compounds in `nmol mol⁻¹
> xylose` — absolute molar yields on the limiting sugar, against an internal standard, ±16 %. The
> furfural pH effect is **274×**, not the ≈1.9× the g kg⁻¹ shares imply. **§C.1's declared gap ('I
> found no paper reporting furfural yields from a pentose … That is a genuine gap') is closed.**"

> **"THREE OF THE FIVE PDFs CANNOT BE HARVESTED FROM THEIR TEXT LAYER.** `wang2008.pdf`'s only
> quantitative content is a bar chart with no glyphs in the plot area; `apriyantono1993.pdf`'s OCR
> corrupts Table 1's exponents and decimal points; `shu1988.pdf`'s OCR collapses Table I's three pH
> columns. All three were read visually from renders. **A naive text-layer harvest of any of them
> would produce wrong numbers.**"

> **"NAMING: two papers in the corpus use an 'HMF'-shaped token to mean NORFURANEOL** — Blank & Fay
> 1996's `HMF (3)` and Apriyantono & Ames 1993's `HMFone`. Neither means 5-hydroxymethylfurfural.
> Apriyantono is the more dangerous, because the same paper separately mentions the real 5-HMF
> (cited to Tressl 1989) and furaneol (cited to Shu & Ho 1988), neither of which it measures.
> Furaneol itself is written `HDMF` by Blank 1996/1997 and Poisson 2019 and `DMHF` by Wang & Ho 2008
> and Shu & Ho 1988; homofuraneol is `HEMF` in Blank 1996 and `EHMF` in Blank 1997."

> **"Three papers carry internal errors that must not propagate:** Wang & Ho 2008 prints its
> statistical criterion inverted ('A *p* value of >0.05 was considered statistically significant',
> in both the Methods and the Figure 1 caption) and mis-cites CAMOLA to Namiki & Hayashi 1982;
> Apriyantono & Ames 1993 reports its nitrogen-compound totals in **mmol** where the arithmetic
> gives **µmol** (a 1000× error); Shu & Ho 1988 calls DMHF 'an important α-dicarbonyl', which it is
> not (it is an enolic 3(2H)-furanone/reductone)."

---

# §11. RANKED FETCH LIST FROM THIS WAVE

Consolidated across the five dossiers. Items 1−3 are named independently by two or more papers.

| # | paper | why | named by |
|---|---|---|---|
| **1** | **★★ Haleva-Toledo, E.; Naim, M.; Zehavi, U.; Rouseff, R. L. 1999**, *"Effects of L-cysteine and N-acetyl-L-cysteine on 4-hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol), 5-(hydroxymethyl)furfural, and 5-methylfurfural formation and browning in buffer solutions containing either rhamnose or glucose and arginine"*, **JAFC 47, 4140−4145** | **THE MISSING SINK MAGNITUDE.** The only identified source that quantifies DMHF inhibition by cysteine against pH **in a buffer** — i.e. the number Shu & Ho structurally establishes and cannot measure. **It also carries 5-HMF and 5-methylfurfural in the same buffer, so it serves the HMF channel (`research_round3_channels.md` §A) at the same time.** | **Wang & Ho 2008 (ref 20)**; raised further by `shu1988_extraction.md` |
| **2** | **★★ Poisson, L.; Auzanneau, N.; Mestdagh, F.; Blank, I.; Davidek, T. 2018**, *"New Insight into the Role of Sucrose in the Generation of α-Diketones upon Coffee Roasting"*, **JAFC 66, 2422−2431** | **THE ABSOLUTE QUANTIFICATION `poisson2019.pdf` LACKS**, same lab, same in-bean system, one year earlier. Everything the brief hoped to find in the 2019 paper and did not. | **Poisson 2019 (ref 3)** |
| **3** | **★ Hofmann, T.; Schieberle, P. 2001**, *"Acetylformoin — an important progenitor of 4-hydroxy-2,5-dimethyl-3(2H)-furanone and 2-acetyltetrahydropyridine during thermal food processing"*, Flavour 2000 Proceedings, 6th Wartburg Aroma Symposium, pp. 311−322 | The **primary source for the acetylformoin → DMHF reduction and for the identity of the reductant** — i.e. for the claim `src/barrier_constants.py:325` currently mis-attributes to Blank & Fay 1996 (**defect D1**). **Fetching it could resolve D1.** ⚠️ Conference proceedings, no DOI; availability doubtful. | **Wang & Ho (ref 17) AND Poisson (ref 19)** |
| **4** | **★ Schieberle, P. 2005**, *"The carbon module labeling (CAMOLA) technique"*, **Ann. N.Y. Acad. Sci. 1043, 236−248** | Primary source for *"HDMF is formed solely via the intact C6-glucose skeleton"* **and** for the acetol + 2-oxopropanal C3/C3 route (Edge B's pathway c) **and** for the claim that MG + 1-hydroxy-2-propanone recombination is *"a major pathway when glucose is reacted with proline"* — a quantified-sounding statement about Edge B's share that both papers only paraphrase. | **Wang & Ho (ref 18) AND Poisson (ref 8)** |
| **5** | **★ Beck, J.; Ledl, F.; Severin, T. 1988**, *"Formation of 1-deoxy-D-erythro-2,3-hexodiulose from Amadori compounds"*, **Carbohydr. Res. 177, 240−243** | The **pH-dependent 1-DG : 3-DG partition (20:1 at pH 7, 8:5 at pH 4.5)** — the enolisation fork that explains three of the four pH ladders in §4 and that both Blank 1997 and Apriyantono & Ames invoke. **Currently a secondary, unread citation doing load-bearing work.** | Blank 1997 |
| **6** | **Shu, C.-K.; Mookherjee, B. D.; Ho, C.-T. 1985b**, *"Volatile Components of the Thermal Degradation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone"*, **JAFC 33, 446−448**, and **Shu, Hagedorn, Mookherjee & Ho 1985c**, *"Two Novel 2-Hydroxy-3(2H)-thiophenones from the Reaction between Cystine and DMHF"*, **JAFC 33, 638−641** | 1985b is **the DMHF self-degradation paper with a pH axis** (Figure 1's middle curve); 1985c is **the source of every "formed by the interaction of cysteine and DMHF" attribution** in Shu & Ho 1988. Without them the sink's route assignments are unsupported. | Shu & Ho 1988 |
| **7** | **Mottram, D. S.; Leseigneur, A. 1990**, *"The effect of pH on the formation of aroma volatiles in meat-like model systems"*, in *Flavour Science and Technology*, Wiley, pp. 121−124 | **Norfuraneol quantified across pH 4.5−6.5 in a RIBOSE−CYSTEINE system, DECREASING with pH.** A pentose, with cysteine, with a pH axis, on the exact compound `Furanone_Cyclisation` produces — **four matches to the repo's sulfur lane at once.** | Apriyantono & Ames 1993 |
| **8** | **Hirvi, T.; Honkanen, E.; Pyysalo, T. 1980**, *"Stability of 2,5-dimethyl-4-hydroxy-3(2H)-furanone and 2,5-dimethyl-4-methoxy-3(2H)-furanone in aqueous buffer solutions"*, **LWT 13, 324−325** | A DMHF stability-vs-pH study **in buffer** — the buffered counterpart of Shu & Ho's Figure 1 middle curve, and the thing that would settle the pH-4-optimum tension in `shu1988_extraction.md` §4.5. Short paper. | Blank 1997 |
| **9** | **Tressl, R.; Helak, B.; Martin, N.; Kersten, E. 1989**, in *Thermal Generation of Aromas*, ACS, pp. 156−171, **and Leahy, M. M.; Reineccius, G. A. 1989**, *"Kinetics of the formation of alkylpyrazines"*, same volume, pp. 196−208 | Tressl carries the **5-HMF pH claim** ('large amounts at pH 3, undetectable at pH 5 or above'), directly serving the HMF channel; Leahy is titled **kinetics** on glucose−lysine at 95 °C/pH 5 — one of very few pyrazine-kinetics sources, and a near-matched hexose comparator to Apriyantono. **Same ACS volume — one acquisition serves both.** | Apriyantono & Ames 1993 |
| **10** | **Hofmann, T.; Schieberle, P. 1997**, **JAFC 45, 898−906** | The rhamnose ≫ hexose ≫ pentose DMHF reactivity ordering, and the claim that the 2,3-dioxo-4,5-dihydroxyhexane route **cannot** operate from pentoses. Would give the 6-deoxysugar arm, on which the repo has nothing. | Wang & Ho (ref 13) |

---

# §12. SOURCES USED BEYOND THE FIVE PDFs

The **volume/issue/pages of `poisson2019.pdf` (67:13829−13839)** are not printed in that PDF (it is
an ASAP proof reading `XXXX, XXX, XXX−XXX`); they were taken from `research_round3_channels.md`
§B.2.2. **That is the only bibliographic fact in this wave not read from a primary document.**
Everything else — including all five DOIs or their explicit absence — is from the PDFs themselves,
from the sibling dossiers in this directory, or from **read-only** inspection of
`src/reaction_templates.py`, `src/barrier_constants.py`, `src/kinetic_core/species.py`,
`src/kinetic_core/sulfur.py`, `src/kinetic_core/network.py`, `src/literature_runtime.py` and
`src/curated_pathways.py` in the working tree. **No web lookup or CrossRef query was performed in
this wave.**

**Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was modified, and
nothing was committed.**

*End of synthesis.*
