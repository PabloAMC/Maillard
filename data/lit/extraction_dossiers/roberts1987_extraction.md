# Roberts, Nagasawa, Zera, Fricke & Goon 1987 — EXTRACTION (eight L-cysteine + aldose thiazolidine-4(R)-carboxylic acids incl. D-ribose-cysteine "RibCys", synthesised as cysteine prodrugs; rat-hepatocyte GSH and mouse acetaminophen-protection data; NO rate constant, NO equilibrium constant, NO pH / temperature / buffer series for ring opening)

### A PHARMACOLOGY PAPER, NOT A KINETICS PAPER. It asserts that the sugar-cysteine thiazolidines release cysteine "by nonenzymatic ring opening and hydrolysis" but never measures how fast, at what pH, at what temperature, or in which buffer. The only physical-chemical content is qualitative synthesis-bench observation (C-2 mutarotation on a minutes scale at room temperature, free cysteine on TLC and in NMR, a pyridine/pH effect on yield). Read it for leads (its refs 3-5), not for numbers.

**Source on disk:** `data/articles/roberts1987.pdf` (6 pp., J. Med. Chem. pp. 1891-1896; p. 1891 opens with the tail of
an unrelated dopamine paper). **Scanned page images (300 dpi CCITT) with an OCR text layer.** Prose and Tables I-III
read from the `pdftotext -layout` layer and **checked against 130 dpi page renders** (scratchpad `img/roberts-1.png` ...
`roberts-6.png`): the OCR layer misreads the bold compound labels (`1d` -> `la`/`Id`, `1g` -> `lg`), so every
compound label below was taken from the renders. Supplementary material (Figs. 1 and 3, 13C NMR of GlcCys) is **not
on disk** and carries no kinetics.

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[F]** read from a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed / hypothesised, not
measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question put to it (plan step B7, `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 7, blocking
item B2): the engine's TTCA ring-opening `k_ttca_cys` (~0.028 min⁻¹, fitted, no pH factor) cannot reconcile Cao 2017
(`cao2017_extraction.md`: 0.4 M Cys + Xyl, 0.2 M phosphate, equilibrium within 30 min at 90 °C) with Zhai 2021
(`zhai2021_extraction.md`: water, TTCA loss <= 0.0033 min⁻¹ at 100 °C). Is there a measured ring-opening rate, its
pH / T dependence, an equilibrium constant, and is buffer (general-acid) catalysis reported?

| what the model needs | what Roberts 1987 gives | tag |
|---|---|---|
| rate of ring opening / cysteine release from RibCys (or any aldose-Cys TCA) | **none.** Release is asserted (Scheme I: TCA ⇌ Schiff base ⇌ carbinolamine ⇌ Cys + aldehyde) but never timed | [P] |
| pH dependence | **none measured.** One directional remark: added pyridine (1 equiv) "increased the overall yields somewhat", attributed to a higher solution pH that shifts the equilibrium toward the thiazolidine (citing Schubert 1936 and Bognar 1970) | [Q] / [C] |
| temperature dependence | **none.** Syntheses at room temperature (aqueous, RibCys) or in refluxing methanol; no rate at any T | — |
| equilibrium constant | **none.** No K, no equilibrium free-Cys fraction | — |
| buffer / general-acid catalysis | **not reported, not tested.** No buffer is used in any chemistry experiment; the only buffers are cell-culture media and a PBS rinse in the hepatocyte assay | — |
| ring-chain dynamics | **qualitative only:** C-2 mutarotation observed, so all specific rotations were taken after a fixed 3-min mixing interval (RibCys: [α]²⁵D −103.1°, c 0.52, H₂O) | [Q] / [M] |
| dissociation visible at room temperature | TLC "invariably" showed free L-cysteine (authors: dissociation on the plate); 1H NMR of all eight TCAs "usually" showed trace free cysteine; recrystallisation gave "dramatic loss of product" | [Q] |
| epimer ratio | two C-2 epimers, not separated; cis:trans "around 1:2" for most compounds (1H NMR, C-2 doublets ~5 ppm, cis 0.1-0.2 ppm downfield) | [M] |

**Answer to B2: nothing quantitative.** Confidence that the paper holds no rate, K, pH series or buffer comparison:
~99 % (all six pages read, including the Experimental). What it does say, weighed:

1. **The ring-chain step is facile in unbuffered water at ~25 °C** [Z, from Q]. Mutarotation that has to be
   pinned to a 3-min convention means C-2 epimerisation (which needs ring opening to the iminium / Schiff base and
   re-closure) proceeds on a minutes scale at room temperature, with no buffer present. At 90-100 °C the same step
   should be far faster. **Caveat:** epimerisation needs only ring opening, not hydrolysis to free cysteine + sugar;
   the engine's `r_ttca_cys` is the full dissociation. The solution pH of the polarimetry samples (~20 mM RibCys
   [Z: 5.2 g/L ÷ 253.3 g/mol], self-buffered zwitterion, pH not stated) is unknown.
2. **Full dissociation is also not slow at room temperature** [Q]: free cysteine on TLC and in NMR, losses on
   recrystallisation. Again no time, pH or concentration.

**Bearing on the Cao-vs-Zhai contradiction (my reading, low confidence ~0.6):** Roberts mildly disfavours the reading
that the uncatalysed ring opening is intrinsically slow in water, i.e. it gives no support to "phosphate makes the
ring-chain step fast, water leaves it slow" as the *sole* explanation of B2. It is equally consistent with (a) a
pH-dependent K (Kallen 1971: apparent K rises ~20× from pH 5.6 to 7.5) controlling how much cysteine Zhai's pot
releases, and (b) buffer catalysis of the slower *hydrolysis* step (carbinolamine ⇌ Cys + aldose) rather than of the
ring opening. It cannot rank these. The rival R-P (general-acid buffer catalysis) remains untested by any paper in
the corpus; this one does not move it.

**Leads it points to (for the owner's download list; none on disk):**
- Nagasawa, Goon & Shirota 1981, *J. Heterocycl. Chem.* 18:1047-1051, "Epimerization at C-2 of 2-substituted
  thiazolidine-4-carboxylic acids" — DOI 10.1002/jhet.5570180538 (Crossref match on title / authors / year,
  2026-10-31). The mechanism paper behind the mutarotation remark; the most likely place for a ring-chain rate.
- Nagasawa, Goon, Muldoon & Zera 1984, *J. Med. Chem.* 27:591 — DOI 10.1021/jm00371a006 (as recorded in
  `docs/guides/LITERATURE_SEARCH_2026-10-29_phosphate_thiazolidine.md`); 2-alkyl TCA prodrugs, possibly
  hydrolysis kinetics.
- Schubert 1939, *J. Biol. Chem.* 130:601-603, "The combination of cysteine with sugars" — DOI
  10.1016/s0021-9258(18)73532-1 (Crossref match, 2026-10-31); Schubert 1936 *J. Biol. Chem.* 114:341 (thiol acids +
  aldehydes; the pH-shifts-equilibrium claim is cited to it). Possible source of sugar-Cys dissociation constants.
- Bognar, Somogyi & Györgydeák 1970, *Justus Liebigs Ann. Chem.* 738:68 (the RibCys synthesis followed here; sugar-
  cysteine ring-chain tautomerism in solution). DOI not verified.
- Pesek 1978, *Anal. Chem.* 50:787 (NMR of thiazolidines). DOI not verified.
- Kallen 1971, *J. Am. Chem. Soc.* 93:6227 — already read (`kallen1971_extraction.md`).

## 0. Identity

| field | value |
|---|---|
| Title | "Prodrugs of L-Cysteine as Protective Agents against Acetaminophen-Induced Hepatotoxicity. 2-(Polyhydroxyalkyl)- and 2-(Polyacetoxyalkyl)thiazolidine-4(R)-carboxylic Acids" |
| Authors | Jeanette C. Roberts, **Herbert T. Nagasawa** (corresponding), Richard T. Zera, Robert F. Fricke, David J. W. Goon |
| Affiliations | Dept. of Medicinal Chemistry, University of Minnesota; Medical Research Laboratories, VA Medical Center, Minneapolis; Hennepin County Medical Center; USAMRIID, Frederick MD |
| Venue | *J. Med. Chem.* 1987, **30**(10), 1891-1896. Received May 4, 1987. Code printed: 0022-2623/87/1830-1891$01.50/0 |
| DOI | 10.1021/jm00393a034 (not printed in the PDF; Crossref returns this title, these five authors, J. Med. Chem. 30(10) 1891-1896, 1987, checked 2026-10-31). **Matches the brief.** |
| Funding | USPHS 5-T32-GM07994 (University of Minnesota), Veterans Administration, USAMRIID |
| Paper type | Medicinal chemistry: synthesis + in vitro (rat hepatocyte GSH) + in vivo (mouse acetaminophen survival) |
| Same lab | Nagasawa 1981 (J. Heterocycl. Chem., epimerisation), 1982 / 1984 (J. Med. Chem., 2-alkyl TCA prodrugs) — none on disk |

## 2. Methods as they matter to a model

- **Compounds** [M]: 2(RS)-polyhydroxyalkyl thiazolidine-4(R)-carboxylic acids from L-cysteine and D-glyceraldehyde
  (GlycCys, **1a**), D-arabinose (AraCys, 1b), D-lyxose (LyxCys, 1c), **D-ribose (RibCys, 1d)**, D-xylose (XylCys, 1e),
  D-galactose (GalCys, 1f), D-glucose (GlcCys, 1g), D-mannose (ManCys, 1h); peracetylated Rib / Gal / Glc analogues
  2a-c. (The OCR layer prints RibCys as "1a"/"Id"; the render reads **1d** throughout.) **XylCys (1e) is the TTCA of
  Cao 2017 and the Zhai papers** (2-(1,2,3,4-tetrahydroxybutyl)thiazolidine-4-carboxylic acid from xylose); RibCys is
  its ribo diastereomer (Parker 2014's "semi-stable thiazolidine", `parker2014_extraction.md`). The engine's TTCA is
  one lumped pentose species.
- **RibCys synthesis** [M]: "as described" in Bognar 1970 (ref 4d), aqueous medium under N₂, absolute ethanol added to
  precipitate, 1 equiv pyridine; 4.71 g, **92.2 % yield**, pale yellow, mp 149-151 °C dec, [α]²⁵D −103.1° (c 0.52,
  H₂O), IR 3220 (br, OH, COO⁻), 1610 cm⁻¹ (COO⁻); anal. C₈H₁₅NO₆S. XylCys: Weitzel procedure (methanol), 74.8 %,
  mp 119-123 °C dec, [α]²⁵D −97.8° (c 1.01, H₂O). Yields across 1a-h 75-95 % except GlycCys 40.5 %.
- **NMR / IR** [M]: 300 MHz; 13C in DMSO-d₆; 1H with TMS or DSS. IR: no NH stretch, carboxylate band ~1610 cm⁻¹
  (zwitterion in the solid); some preparations show a 1700-1750 cm⁻¹ carbonyl plus a sharp NH at ~3250 cm⁻¹
  (the un-ionised form). 13C of GlcCys: C-2 at 70.18 / 71.76 ppm (two epimers), upfield of glucose's anomeric C
  (92.34 ppm), confirming the thiazolidine (not glycosylamine) form.
- **Acetylated analogues** [Q]: "appeared to be more stable" — TLC showed no free L-cysteine.
- **Hepatocyte assay** [M]: primary rat hepatocytes (Seglen), cultured 24 h, then 4 h with prodrug 1.0 mM (NAC 2.5
  mM); total GSH (GSH + GSSG) by DTNB-glutathione reductase recycling; BSO 0.20 mM pre-exposure. Temperature of
  incubation not printed.
- **Mouse assay** [M]: male Swiss-Webster, fasted, acetaminophen 750 mg/kg ip then prodrug 2.45 mmol/kg ip 30 min
  later; survival at 48 h; necrosis score.

## 3. Numbers printed (none is a kinetic or thermodynamic constant)

| quantity | value | units | conditions | where | tag |
|---|---|---|---|---|---|
| RibCys synthesis yield | 92.2 | % | aqueous, N₂, 1 equiv pyridine, EtOH precipitation | p. 1895 | [M] |
| yields, series | 75-95 (GlycCys 40.5) | % | per compound | pp. 1892, 1895 | [M] |
| C-2 cis:trans | ~1:2 | — | 1H NMR, most compounds, solvent not stated per compound | p. 1892 | [M] (approximate) |
| polarimetry convention | 3 | min after mixing | because of mutarotation (C-2 epimerisation) | p. 1892 | [M] |
| RibCys [α]²⁵D | −103.1 | deg | c 0.52 g/100 mL, H₂O (~20 mM [Z]) | p. 1895 | [M] |
| GSH in hepatocytes, RibCys | 61.2 ± 1.52 (1.7× control 35.4 ± 0.78) | nmol/10⁶ cells | 1.0 mM, 4 h | Table I | [M] |
| GSH, XylCys | 58.3 ± 0.99 (1.6×) | nmol/10⁶ cells | 1.0 mM, 4 h | Table I | [M] |
| GSH with BSO, RibCys | 16.2 ± 3.60 (0.5×; BSO alone 18.4 ± 2.08) | nmol/10⁶ cells | + 0.20 mM BSO | Table II | [M] |
| mouse survival, RibCys | 12/12 (100 %) vs 2/12 (17 %) acetaminophen alone | — | 48 h | Table III | [M] |

The hepatocyte numbers show cysteine becomes available to GSH synthesis within 4 h at 1 mM in culture, but release
may be cellular, enzymatic or chemical; there is no cell-free control. **Not usable as a release rate.**

## 4. Figures

Fig. 2 (DEPT spectra of GlcCys) only; Figs. 1 and 3 are supplementary 13C spectra. **Nothing digitisable for the
model; no two-reader keying needed.**

## 5. What the repository takes

Nothing quantitative: no FIT row, no constant, no pot. Two qualitative statements may be cited as context for the
TTCA ring-chain step: (i) sugar-cysteine thiazolidines epimerise at C-2 (open-ring intermediate) on a minutes scale at
room temperature in unbuffered water [Q]; (ii) higher pH (added pyridine) shifts the Cys + aldose ⇌ TCA equilibrium
toward the thiazolidine [Q, C] — the same sign as Kallen 1971's measured K(pH). Neither is a measured rate or
constant; neither may be promoted to a prior.

## 6. Flags

1. **OCR labels are wrong in the text layer** (RibCys printed as "1a" in the abstract layer and "Id" in tables; GlcCys as
   "lg"). Use the renders.
2. **No temperature, pH or buffer is stated for any observation of dissociation** (TLC, NMR, mutarotation,
   recrystallisation). Do not back-compute a rate from "3-min mixing interval".
3. **"Mutarotation" ≠ cysteine release.** C-2 epimerisation proceeds through the open iminium; it does not require
   hydrolysis to free cysteine. A fast ring-chain step and a slow hydrolysis are compatible.
4. **TLC artefact.** The authors themselves attribute free cysteine on TLC to dissociation on the silica surface, not
   in solution.
5. **Registry:** compounds `cysteine`, `ribose` present in `data/keys/compounds.yml`; RibCys / XylCys (TTCA) are not
   registered as compound keys.
