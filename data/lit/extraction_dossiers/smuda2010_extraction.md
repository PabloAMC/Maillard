# Smuda, Voigt & Glomb 2010 — EXTRACTION
### *"Degradation of 1-Deoxy-D-erythro-hexo-2,3-diulose in the Presence of Lysine Leads to Formation of Carboxylic Acid Amides"*
### Synthetic 1-deoxyglucosone (1-DG) + Nα-t-BOC-lysine, 42 mM each, 0.1 M phosphate pH 7.4, **37 °C**, 0-7 days, aerated and deaerated. Four acids (acetic, glyceric, lactic, formic) and their four Nε-lysine amides quantified. **No amine-free arm**: the paper cannot say whether an amine accelerates the beta-cleavage; it does bound the amine-as-nucleophile channel at < 1 % of the acetyl flux.

**Read 2026-10-31 from `data/articles/smuda2010.pdf`** (7 pp., JAFC 58:6458-6464; born-digital, clean text layer;
Fig. 2 read from a 300 dpi render of p. 6461, read error ~ +-5 mmol/mol on panel A, ~ +-0.02 on panel B).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Question (B7a-2, `results/validation/kinetic_core_dpo_partition_prereg.md` secs. 11 and 13; finding 2 / the Davidek
tension in `docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 1): does an AMINE accelerate the beta-dicarbonyl
cleavage of 1-DG (rate or yield with vs without amine; amine-concentration, pH, temperature dependence), and what are
the products' shares?

| need | what Smuda 2010 gives | tag |
|---|---|---|
| (a) with vs without amine | **None.** Every degradation pot contains 42 mM Nα-t-BOC-lysine; no 1-DG-only control is reported. The only amine-free incubations are acid + BOC-lysine controls (the reverse question) | [M] |
| (b) amine concentration | **None** (one ratio, 1-DG : BOC-Lys = 1 : 1, 42 mM each) | [M] |
| (c) pH / temperature | **None** (pH 7.4, 37 °C only; pH not re-measured) | [M] |
| (d) acid yields, 7 d (aerated) | acetic **310**, glyceric **88**, lactic **58**, formic **27** mmol / mol 1-DG | [M] |
| (e) amide yields, 7 d (aerated) | Nε-formyl **1.48**, Nε-acetyl **1.01**, Nε-lactoyl **0.28**, Nε-glycerinyl **0.16** mmol / mol 1-DG | [M] |
| (f) amide : acid ratio (amine vs water / hydroxide as the attacking nucleophile) | acetyl 1.01 / 310 = **0.33 %**; glycerinyl 0.16 / 88 = 0.18 %; lactoyl 0.28 / 58 = 0.48 %; formyl 1.48 / 27 = **5.5 %** | [Z] |
| (g) time scale of the acid release | acetic acid ~50 / ~122 / ~175 / ~220 / ~250 / ~290 / 310 mmol/mol at ~1 / 3 / 5 / 12 / 24 / 72 / 168 h; half of the 7-d value by ~4 h at 37 °C | [F] / [M] |
| (h) O2 | amides only (Table 3): formyl 1.47 vs 1.08 and glycerinyl 0.15 vs 0.09 mmol/mol at 72 h (aerated vs deaerated); acetyl 0.78 vs 0.75 and lactoyl 0.25 vs 0.26 (no O2 effect). Acids not run deaerated here ("data will be published elsewhere" = Voigt, Smuda, Pfahler & Glomb 2010, JAFC 58:5685, 10.1021/jf100140h, in the download queue as voigt2010) | [M] / [Q] |

**Answer.** Smuda 2010 **does not test amine acceleration**: there is no 1-DG run without lysine, no lysine
concentration series, one pH, one temperature. Confidence that the paper holds no with / without-amine comparison:
**~98 %** (all seven pages read; the acid controls are acid + amine, not 1-DG alone). What it does settle:

1. **The amine-INDUCED cleavage (amine as the nucleophile that splits the 2,4-diulose) is a minor channel.** At 42 mM
   lysine side chain, pH 7.4, 37 °C, Nε-acetyl lysine is 0.33 % of the acetic acid formed [Z from M]; even the most
   amine-favoured pair (formyl) is 5.5 %. The control (13 mM acids + 42 mM BOC-Lys, 7 d) gave **no** amides, so the
   amides are not condensation artefacts [M]. So a law in which the amine itself cleaves DPO (rate proportional to free
   amine) is not what produces the acetic acid; **~90 %** confidence that this transfers as "amine-as-nucleophile < 1 %
   of the beta-cleavage flux" to glycine pots (glycine's alpha-amine pKa ~9.6 vs lysine epsilon ~10.5, so more free
   base at pH 7: maybe x5-10, still < 5 %) [Z].
2. **It does not exclude amine CATALYSIS** (e.g. an amine speeding the 2,3- to 2,4-diulose isomerisation that precedes
   the hydroxide attack, Fig. 3 [P]). Nothing here measures it. With Davidek 2006a (amine-free 1-DG gives
   0.80-0.84 mol/mol acetic acid at pH >= 6.5, 120 °C; davidek2006a_extraction.md) the joint reading is unchanged:
   an amine is not required, its acceleration is untested.
3. **Yield context, not comparable across papers.** 0.31 mol acetic acid / mol 1-DG at 37 °C, 7 d, with lysine,
   vs Davidek's 0.80-0.84 at 120 °C, 1 h, without amine: different temperature, buffer (0.1 vs ~0.16 M phosphate),
   time and 1-DG assay (both nominal). One cannot read an amine effect off this pair. If anything the lower 37 °C yield
   says other 1-DG sinks (amine condensation, cyclisation) compete more at 37 °C with amine present; it is not
   evidence about acceleration. Confidence in "no inference possible": ~85 %.
4. **Product shares for the hexose beta-cleavage at pH 7.4 [M]:** acetic : glyceric : lactic : formic = 310 : 88 : 58 : 27
   (64 / 18 / 12 / 6 % of the four acids; sum 0.48 mol/mol). Acetic acid dominates as in Davidek; glyceric acid is
   ~0.28 x acetic (attack at C-4 of the 2,4-diulose). Lactic and formic acid are assigned (Fig. 3) to beta-cleavage of
   the 3,5- and 4,6-diuloses, and partly (Fig. 4) to an O2 route via tricarbonyls [P].
5. **For B7a-2:** this paper neither supports nor refutes an amine-multiplied term; it removes one mechanism for it
   (direct aminolysis). The lead that could answer acceleration is Voigt & Glomb 2009 (JAFC 57:4765-4770,
   10.1021/jf900459x, Crossref 2026-10-31), the paper cited as ref 20 for 1-DG reactivity; not on disk.

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Mareen Smuda, Michael Voigt, Marcus A. Glomb* (Martin-Luther-University Halle-Wittenberg) |
| Title | as in the header |
| Venue | *J. Agric. Food Chem.* **2010**, 58 (10), 6458-6464 |
| DOI | **10.1021/jf100334r** (printed p. 6458; Crossref 2026-10-31: same title, authors, vol. 58, issue 10, pp. 6458-6464) |
| Dates | received 26 Jan 2010; revised 14 Apr 2010; accepted 17 Apr 2010; published on web 29 Apr 2010 |
| Content | Tables 1 (MRM parameters), 2 (CID fragments), 3 (amides aerated vs deaerated, 24 / 72 h); Figs. 1 (syntheses of lactoyl / glycerinyl lysine), 2 (acids A, amides B, 0-168 h), 3 (beta-cleavage mechanism), 4 (tricarbonyl route) |
| Brief's description | "amine-induced beta-dicarbonyl cleavage of 1-DG": correct for the amide channel; the paper's own title is about amide formation |

## 3. Methods [M]

| item | as printed |
|---|---|
| 1-DG | protected 1-deoxy-4,5-O-isopropylidene-D-erythro-hexo-2,3-diulose 0.35 mmol in 3 mL water + Dowex 50W x8 (H+) 4 mL, 4 h, argon; resin filtered and washed with MeOH; solvents evaporated; residue in 4.2 mL 0.1 M phosphate pH 7.4 |
| Incubation | 150 µL 1-DG solution + 150 µL BOC-Lys (0.35 mmol in 4.2 mL same buffer) in screw-cap vials -> **42 mM each** (nominal; deprotection yield and evaporation losses not printed [!]); shaken at **37 °C**; sampled over 7 d; each sample "at least three times" |
| Deaerated | phosphate buffer + 1 mM DTPA, degassed with helium; samples deaerated with argon. Aerated pots have **no DTPA** [!]: aerated vs deaerated confounds O2 with chelation |
| Control | acetic, lactic, glyceric, formic acid 13 mM each + BOC-Lys 42 mM, same buffer, 37 °C, 7 d: corresponding amides "not detected" |
| Analytics | amides: HPLC-MS/MS (MRM) after 6 M HCl deprotection, standard addition, SD < 0.07 mmol/mol (CV < 5 %); lactic, glyceric: GC-MS (TMS), SD < 5 mmol/mol; acetic: GC-FID (decyl chloroformate), SD < 10 mmol/mol (CV < 3.5 %); formic: enzymatic kit, CV < 2.5 % |
| Not measured | 1-DG decay, pH drift, C4 counterparts (erythrulose, etc.) |

## 4. Results

**Fig. 2A, acids (mmol / mol 1-DG), aerated, 37 °C** (7-d values printed [M]; earlier points [F]; time points
~1, 3, 5, 12, 24, 72, 168 h):

| h | acetic | glyceric | lactic | formic |
|---|---|---|---|---|
| ~1 | ~50 | ~29 | ~10 | ~3 |
| ~3 | ~122 | ~50 | ~13 | ~6 |
| ~5 | ~175 | ~68 | ~21 | ~8 |
| 12 | ~220 | ~79 | ~30 | ~12 |
| 24 | ~250 | ~82 | ~33 | ~17 |
| 72 | ~290 | ~84 | ~47 | ~24 |
| 168 | **310** | **88** | **58** | **27** |

"Within the first twelve hours a rapid accumulation of all carboxylic acids was observed" [V]; then a marginal rise.

**Fig. 2B / Table 3, amides (mmol / mol 1-DG):** 7-d values formyl 1.48, acetyl 1.01, lactoyl 0.28, glycerinyl 0.16
[M]. Formyl ~0.36 / 0.53 / 0.65 / 0.90 / 1.12 / 1.48 at ~1 / 3 / 5 / 12 / 24 / 72 h; acetyl ~0.17 / 0.21 / 0.27 /
0.38 / 0.49 / 0.78 [F].

| Table 3 | 24 h aerated | 24 h deaerated | 72 h aerated | 72 h deaerated |
|---|---|---|---|---|
| Nε-acetyl lysine | 0.48 | 0.45 | 0.78 | 0.75 |
| Nε-formyl lysine | 1.12 | 0.71 | 1.47 | 1.08 |
| Nε-lactoyl lysine | 0.19 | 0.16 | 0.25 | 0.26 |
| Nε-glycerinyl lysine | 0.11 | 0.05 | 0.15 | 0.09 |

**Mechanism.** Proposed [P] (Fig. 3): 1-DG isomerises to 1-deoxyhexo-2,4-diulose (and 3,5- / 4,6-diuloses); the
lysine epsilon-amine attacks C-2 (-> Nε-acetyl lysine + C4 enediol) or C-4 (-> Nε-glycerinyl lysine), "analogous" to
hydroxide attack giving the acids. Fig. 4 [P]: O2 converts 1-DG or 1-deoxy-2,3-tetrodiulose to tricarbonyls whose
cleavage gives glyceric, lactic, formic acid / amides. Shown: the amides exist, track the acids in time, are not made
from free acids, and formyl / glycerinyl amides rise with O2 (with the DTPA confound). The 2,4-diulose intermediate
and the amine attack are not observed.

## 5. Supplementary material

**None.** The paper mentions no Supporting Information (text layer searched for "supporting" / "supplement": no hit;
no SI paragraph before the references). Every number the model could use is in Fig. 2 and Table 3.

## 6. Flags

- [!] No amine-free 1-DG pot; no concentration, pH or temperature series. The paper cannot answer the B7a-2
  acceleration question.
- [!] 1-DG concentration nominal (protected-precursor moles; deprotection / evaporation losses unprinted); yields
  are per nominal 1-DG.
- [!] Aerated and deaerated pots differ in DTPA as well as O2.
- [!] BOC-lysine, not free lysine or glycine: a single, alpha-protected epsilon-amine (pKa ~10.5 [C, general]).
- The authors' statement that formic acid is lowest among acids but formyl amide highest is read by them as amine
  attack competing "more strongly" for the C-1 formyl fragment [Q]; the absolute amide yields remain <= 1.5 mmol/mol.
