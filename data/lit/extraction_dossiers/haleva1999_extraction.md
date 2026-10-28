# Haleva-Toledo, Naim, Zehavi & Rouseff 1999 — EXTRACTION (glucose or rhamnose 55 mM ± L-arginine 55 mM, ± L-cysteine or N-acetyl-L-cysteine 1 / 5 / 10 mM, 0.1 M citrate-phosphate pH 3 / 5 / 7, 70 °C, 48 h, sealed brown bottles: furaneol, HMF, 5-methylfurfural (HPLC-UV 290 nm), A420; plus furaneol / HMF / MF 20 mg/L alone at pH 3.5 + Cys or NAC 0-2.5 mM, 70 °C, 48 h)

### ONE TEMPERATURE (70 °C), ONE TIME POINT (48 h), NO RATE CONSTANTS. What it gives: (1) cysteine and N-acetylcysteine consume pre-formed furaneol DIRECTLY in sugar-free buffer at pH 3.5 — a thiol + furanone sink, with a dose-response from which an apparent second-order constant of **≈ 2 M⁻¹ min⁻¹ (range 1.2-3.9; 70 °C, pH 3.5)** can be backed out (my arithmetic); HMF and MF are NOT attacked; (2) furaneol "self-degradation" of **20 % in 48 h at 70 °C, pH 3.5 (≈ 7.5 × 10⁻⁵ min⁻¹)**; (3) in sugar + arginine pots the same thiols do NOT suppress furaneol at pH 3, cysteine RAISES it 2.5-3× at pH 5 and 4-5× at pH 7 (it is also an amine), and only NAC suppresses it (to ~0.4× at 10 mM). Nothing on the temperature dependence of furanone loss, nothing on its pH dependence, nothing on norfuraneol.

**Source on disk:** `data/articles/haleva1999.pdf` (6 pp., publisher PDF from pubs.acs.org, clean text
layer). Prose, Table 1 and Table 2 from `pdftotext -layout`. **Figs. 1, 2, 4 and 5 are figure-read** from
200 dpi page renders and 400 dpi crops (scratchpad `lit5/hal-3.png`, `hal-4.png`, `hal_fig1-3.png`,
`hal_fig2-3.png`) and tagged **[D]** with read error. No SI.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of L-Cysteine and N-Acetyl-L-cysteine on 4-Hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol), 5-(Hydroxymethyl)furfural, and 5-Methylfurfural Formation and Browning in Buffer Solutions Containing either Rhamnose or Glucose and Arginine" |
| Authors | Einat Haleva-Toledo, **Michael Naim** (corresponding, naim@agri.huji.ac.il), Uri Zehavi (Hebrew University of Jerusalem, Rehovot); Russell L. Rouseff (Citrus REC, Univ. of Florida, Lake Alfred) |
| Venue | J. Agric. Food Chem. **47** (10), 4140-4145, 1999. Received 23 Dec 1998, revised/accepted 20 Jul 1999, web 18 Sep 1999 |
| DOI as printed | **10.1021/jf9813788** ✔ matches the brief. **Identity verified.** |
| Funding | BARD I-1967-91, US-2914-97 |
| Companion | Haleva-Toledo et al. 1997, JAFC 45:1314 (furaneol formation in buffers and citrus models; NOT on disk) — the source of the claims "furaneol from glucose only at pH > 6" and "from rhamnose + Arg at all pH" |

**Compound.** Furaneol = 4-hydroxy-2,5-**dimethyl**-3(2H)-furanone (DMHF, HDMF; C6H8O3, M 128.13). The
engine's `NF` is norfuraneol = 4-hydroxy-5-**methyl**-3(2H)-furanone (C5H6O3). **Norfuraneol is not
studied in this paper.**

## 1. Why it matters

Three questions from wave B61 (`results/validation/kinetic_core_b61_prereg.md`): (a) any furanone loss
kinetics, especially at two or more temperatures (the declared 103.7 kJ/mol barrier of `r_nf_decay_oh`);
(b) does cysteine consume the furanone directly or act upstream, and by how much; (c) pH dependence of
furanone stability. The paper answers (b) qualitatively and semi-quantitatively, gives one point for (a)
at one temperature and one pH, and says nothing usable on (c).

## 2. Methods as they matter to a model

- **Sugar pots.** 0.1 M citrate-phosphate buffer, **pH 3, 3.5, 5 and 7** (pH 3.5 used only for the furan
  stability test), with **55 mM rhamnose or glucose**, ± **55 mM L-arginine**, ± **1, 5 or 10 mM** Cys or
  AcCys. "Immediately stored in sealed 20-mL brown bottles for **48 h at 70 °C**." Fill volume and
  headspace atmosphere **not stated** (air presumed; an argon-sparged repeat of the browning arm "was
  not modified", data not shown). **Controls stored at 4 °C** (the "no thiols" bars are 70 °C samples
  without thiol; the 4 °C controls are not shown).
- **Browning at pH 3:** 138 mM sugar + 138 mM Arg, Cys or AcCys 1 / 5 / 10 / 25 / 50 mM, same storage.
- **Direct thiol-furan test:** buffer **pH 3.5**, **20 mg/L furaneol** (= 0.156 mM), HMF or MF, with
  **0, 0.01, 0.05, 0.15, 0.5 or 2.5 mM** Cys or AcCys, "incubated under the same conditions" (70 °C,
  48 h). No sugar, no amino acid.
- **pH after heating:** **not measured / not reported.** Citrate-phosphate 0.1 M is a good buffer; Cys up
  to 10 mM is unlikely to shift it, but 55 mM Arg (free base vs HCl salt not stated) could.
- **Analysis.** C18 Sep-Pak SPE (2 mL sample → 1.5 mL MeOH; Walsh 1997 method), RP-18 HPLC, UV
  290 nm, external standards. **SPE recovery not reported.** n = 3 samples × 2 injections, mean ± SEM.
  A420 for browning (dilution not stated; Fig. 4A pH 7 ODs of ~12 imply back-calculation from dilution).
- **Sampling:** **a single endpoint (48 h).** No time courses anywhere in the paper.

## 3. Tables and figures

### Table 1 (verbatim, mg/L, mean ± SEM, n = 3) — HMF and MF after 48 h at 70 °C

| pH | pot | no thiol | Cys 1 mM | Cys 5 mM | Cys 10 mM | AcCys 1 mM | AcCys 5 mM | AcCys 10 mM |
|---|---|---|---|---|---|---|---|---|
| **HMF** | | | | | | | | |
| 3 | Glc | 5 ± 1 | 1.4 ± 0.1* | 0.7 ± 0.2* | 0.2 ± 0.1* | 2.7 ± 0.4* | 2 ± 0.5* | 0.6 ± 0.2* |
| 3 | Glc + Arg | 92 ± 4 | 90 ± 5 | 25 ± 6* | 6 ± 1* | 118 ± 5 | 77 ± 7* | 13 ± 1* |
| 5 | Glc | 0.2 ± 0.1 | 0.1 ± 0.1* | ND | ND | 0.1 ± 0.04* | ND | ND |
| 5 | Glc + Arg | 6 ± 1 | 7 ± 1 | 3.4 ± 0.6* | ND | 6 ± 1 | 3 ± 1* | 2.5 ± 0.3* |
| **MF** | | | | | | | | |
| 3 | Rhm | 3 ± 0.1 | 1 ± 0.1* | 0.4 ± 0.1* | 0.3 ± 0.1* | 1.4 ± 0.5* | 1 ± 0.1* | 0.7 ± 0.1* |
| 3 | Rhm + Arg | 27 ± 3 | 7 ± 0.4* | 1.4 ± 0.2* | 1 ± 0.1* | 12 ± 1* | 7 ± 1* | 4.5 ± 0.3* |
| 5 | Rhm | 2.4 ± 0.6 | 0.5 ± 0.1* | 0.5 ± 0.1* | 0.3 ± 0.1* | 1.1 ± 0.1* | 0.8 ± 0.1* | 0.1 ± 0.01* |
| 5 | Rhm + Arg | 8 ± 1 | 2.9 ± 0.5* | ND | ND | 4.8 ± 0.5* | 4.2 ± 0.8* | 2.6 ± 0.6* |

\* p ≤ 0.01 below control. pH 7 not tabulated (HMF/MF "formed mainly under acidic conditions").
Glc + Arg at pH 3 gives 92 mg/L HMF = 0.73 mM = **1.3 mol % of the 55 mM glucose** in 48 h (mine).

### Fig. 2 — furaneol 20 mg/L (0.156 mM) alone, pH 3.5, 70 °C, 48 h, ± thiol [D]

Linear axis 0-20 mg/L; read error **± 0.3 mg/L** (bars are flat-topped; SEM bars shown).

| thiol, mM | 0 | 0.01 | 0.05 | 0.15 | 0.5 | 2.5 |
|---|---|---|---|---|---|---|
| Cys: furaneol, mg/L | **16.1** | 16.1 | 12.4* | 9.4* | 0.4* | ND |
| AcCys: furaneol, mg/L | 16.1 | 15.3 | 13.6* | 3.0* | 1.6* | ND |
| Cys: fraction of thiol-free control | 1 | 1.00 | 0.77 | 0.58 | 0.025 | < LOD |
| AcCys: fraction of control | 1 | 0.95 | 0.84 | 0.19 | 0.10 | < LOD |

Text (p. 4142): "Regardless of thiol fortification, **self-degradation of Furaneol, HMF, and MF amounted to
20 % during 48 h of incubation at 70 °C**"; HMF and MF "contents were not affected" by Cys or AcCys (data
not shown). The 16.1 mg/L control (80.5 % of nominal) agrees with the 20 %.

### Fig. 1 — furaneol formed in rhamnose pots, 70 °C, 48 h [D]

Panels A and C on **log axes** (read error ± 0.05 dex); panel B linear 0-0.30 mg/L (read error
± 0.008 mg/L). Open bar = no thiol; then 1 / 5 / 10 mM thiol.

| panel | pot | no thiol | Cys 1 | Cys 5 | Cys 10 | AcCys 1 | AcCys 5 | AcCys 10 |
|---|---|---|---|---|---|---|---|---|
| A, pH 5 | Rhm | ND | ND | 0.35 | 0.38 | ND | ND | ND |
| A, pH 5 | Rhm + Arg | 9.4 (Cys panel) / 10.3 (AcCys panel) | 24 + | 30 + | 24 + | 7.2 | 6.6 | 4.3 * |
| B, pH 3 | Rhm + Arg | 0.168 | 0.119 | 0.119 | 0.178 | 0.158 | 0.071 * | 0.056 * |
| C, pH 7 | Rhm | 1.37 | 2.3 | 6.5 + | 5.7 + | 1.16 | 1.31 | 0.62 |

Text, not plotted: pH 7 Rhm + Arg — furaneol "formed but its concentration was not affected by the presence
of either Cys or AcCys (data not shown)"; pH 7 glucose ± Arg — "small amounts (∼2 mg/L)", unaffected by
either thiol. Conversions (mine): 30 mg/L = 0.23 mM = **0.43 mol % of 55 mM rhamnose**; the pH 5 Rhm + Arg
control 9.4-10.3 mg/L ≈ 0.14 mol %; pH 3 0.168 mg/L ≈ 0.002 mol %.

### Figs. 4, 5 and Table 2 — browning (A420), qualitative use only [D, ± 5 % of each panel's scale]

- pH 7, Glc + Arg (Fig. 4A): control ~12; Cys 1 / 5 / 10 mM ~13 / 10.5* / 7*; AcCys ~7* / 7.5* / 6.8*.
  Rhm + Arg (Fig. 4B): control ~0.77; Cys 10 mM ~0.51*; AcCys 10 mM ~0.48*.
- pH 5: Cys raises A420 in Glc + Arg (~0.27 → 0.55-0.62 at 1-5 mM); AcCys lowers it.
- pH 3 (Fig. 5, 138 mM sugar + Arg): any thiol raises A420 from ~0.015 to 0.06-0.15, **not
  dose-monotonic** (maximal at 1-5 mM). Table 2: Cys adds a 320-330 nm band (glucose) and a 380-390 nm band
  (rhamnose), unidentified sugar-Cys products.

## 4. Arithmetic (mine)

1. **Furaneol thermal loss, pH 3.5, 70 °C:** k₀ = ln(20/16.1)/2880 min = **7.5 × 10⁻⁵ min⁻¹** (t½ ≈ 150 h).
   This is an **upper bound** on the true thermal loss: the comparison is against nominal 20 mg/L, and SPE
   recovery is not reported (a 4 °C control for this arm is not shown). One temperature, one pH.
2. **Thiol + furaneol apparent second-order constant (70 °C, pH 3.5).** ODE dA/dt = −k₀A − k₂AB with k₀
   from (1), fitted pointwise to each Fig. 2 bar (scratchpad `lit5/haleva_k2.py`, run in the maillard
   env). Two limits: B constant (thiol not consumed) or B consumed 1 : 1.

   | thiol, mM | Cys: k₂ (B const / B 1:1) | AcCys: k₂ (B const / B 1:1) |
   |---|---|---|
   | 0.01 | — (no loss read) | 1.8 / 3.2 |
   | 0.05 | 1.8 / 3.1 | 1.2 / 1.6 |
   | 0.15 | 1.3 / 1.6 | 3.9 / 10 |
   | 0.5 | 2.6 / 3.3 | 1.6 / 2.0 |

   **k₂ ≈ 2 M⁻¹ min⁻¹ (pointwise range 1.2-3.9 with B constant; 1.6-10 with 1:1 consumption)**, the same
   for Cys and AcCys within the scatter. A single endpoint cannot distinguish first order in thiol from a
   saturating or catalytic law; the 0.01 mM Cys bar (no detectable loss where k₂ = 2 predicts 5 %) is
   within read error. The thiol pool itself is not measured (Cys oxidation to cystine over 48 h at 70 °C in
   a sealed air-headspace bottle is plausible and would make these k₂ lower bounds).
3. **What that sink would do in a pot, and why it does not:** at 10 mM free thiol, k₂B ≈ 0.02 min⁻¹ at
   70 °C — t½ ≈ 35 min. If it acted on freely in the Rhm + Arg pots, furaneol accumulating over 48 h would
   be held at ~1/(k₂B·t) ≈ 1-2 % of the thiol-free level. **Measured: 33 % with 10 mM AcCys at pH 3, 100 %
   with 10 mM Cys at pH 3, 250-300 % with Cys at pH 5.** So in the presence of 55 mM sugar the direct sink
   is suppressed by 15-50× or more — most plausibly because the thiol is tied up (Cys as the thiazolidine /
   hemithioacetal with the 55 mM aldose; AcCys as a hemithioacetal) or consumed by the HMF/MF precursors
   (Table 1 shows both thiols suppress those upstream). **Do not apply k₂ to total charged cysteine.**
4. **Cysteine as an amine outweighs cysteine as a thiol sink at pH 5-7.** At pH 5 (Rhm + Arg) Cys 1-10 mM
   raises furaneol 2.5-3.2×; at pH 7 (Rhm alone) 5-10 mM Cys raises it 4-5× — above what AcCys does
   (0.45-1.0×). The authors attribute it to the free amino group (Cys as Maillard amine, Hofmann &
   Schieberle 1997). **Net effect of Cys on furanone accumulation is positive at pH ≥ 5, null at pH 3**;
   NAC (thiol without a reactive amine) gives the pure-thiol effect: 0.42× at 10 mM and pH 5, 0.33× at
   pH 3, 0.45× at pH 7 (n.s.).
5. **Cross-compound comparison with the engine (illustrative only, NOT evidence):** the engine's NF loss at
   pH 3.5, 70 °C is the neutral `r_nf_decay` (log10 k(145 °C) = −2.18, 143.5 kJ/mol) → 8 × 10⁻⁷ min⁻¹,
   plus B61's hydroxide term (−4.495 at 100 °C, pH 5; 103.7 kJ/mol; ×10^(−1.5)) → 5 × 10⁻⁸ min⁻¹: total
   ≈ 9 × 10⁻⁷ min⁻¹. Haleva's furaneol bound (1) is ~90× faster. Were furaneol a faithful proxy for NF and
   the 20 % pure chemistry, an Arrhenius line through this point and the engine's 145 °C value would
   need ≈ 71 kJ/mol. Both "were"s are unsupported; recorded so nobody re-derives it as a finding.

## 5. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| furaneol thermal loss | first-order k (bound) | ≤ 7.5 × 10⁻⁵ | min⁻¹ | 0.156 mM, 0.1 M citrate-phosphate pH 3.5, 70 °C, 48 h, one endpoint | text p. 4142 + Fig. 2 [D]; arithmetic 1 | derived upper bound (recovery confounded) |
| HMF, MF thermal loss | fraction lost | ~20 % | — | 20 mg/L, same | text p. 4142 | stated, rounded |
| thiol + furaneol (direct sink) | apparent k₂ | ≈ 2 (1.2-3.9; ≤ 10) | M⁻¹ min⁻¹ | Cys or NAC 0.01-0.5 mM, no sugar, pH 3.5, 70 °C | Fig. 2 [D]; arithmetic 2 | derived_from_digitised, single endpoint, rate law assumed |
| thiol + HMF / MF | — | no reaction | — | same | text p. 4142 | stated, data not shown |
| Cys effect on furaneol (Rhm + Arg) | ratio to control | 2.5-3.2 (pH 5); ~1 (pH 3, n.s.) | — | 55 mM Rhm + 55 mM Arg, 1-10 mM Cys, 70 °C, 48 h | Fig. 1A, 1B [D] | digitised |
| NAC effect on furaneol (Rhm + Arg) | ratio to control | 0.42 (pH 5, 10 mM); 0.33 (pH 3, 10 mM) | — | same | Fig. 1A, 1B [D] | digitised |
| Cys / NAC on HMF, MF | dose-response | Table 1 | mg/L | pH 3, 5, 70 °C, 48 h | Table 1 | tabulated |
| furaneol formation / loss Ea | — | **NOT REPORTED** (70 °C only) | — | — | — | — |
| furaneol stability vs pH | — | **NOT REPORTED** (pH 3.5 only) | — | — | — | — |
| sink product | — | **NOT IDENTIFIED** (proposed 3(2H)-thiophenone, citing Zheng 1997) | — | — | Fig. 3, text | hypothesis |

## 6. What is NOT reported

- No time course: every number is a 48 h endpoint. No formation or loss rate constant.
- No second temperature (70 °C throughout): **no activation energy for anything**.
- Furanone stability only at pH 3.5; no pH series of furaneol stability. (The usual citation for
  furaneol's pH-dependent stability is Hirvi, Honkanen & Pyysalo 1980, Z. Lebensm. Unters. Forsch. —
  **not on disk, not verified here.**)
- No norfuraneol, no pentose.
- No H2S, no sulfur volatiles, no product of the thiol + furaneol reaction identified.
- No free-thiol or cystine measurement; post-incubation pH not measured; SPE recovery not reported; 4 °C
  controls not shown.

## 7. Flags

1. **Fig. 1A has two different "no thiol" bars for the same Rhm + Arg pH 5 control** (9.4 in the Cys
   panel, 10.3 in the AcCys panel). Same samples presumably; ± 5 % is the replicate-level noise.
2. **Fig. 1B AcCys group has no 0 mM bar**; read against the Cys-panel control (0.168).
3. The authors' mechanism (thiol attacks the carbonyl of the ring-opened form, O → S exchange to a
   3(2H)-thiophenone) is proposed by analogy to Zheng et al. 1997 and the open/closed equilibrium of Fig. 3;
   **no product was detected**. It implies the sink is a sulfur-volatile *source*, not an inert loss.
4. "Under acidic conditions ... readily exchanged with the sulfur atom of the sulfide" — the citation
   (Zheng 1997) is about furaneol + H2S / cysteine at roasting temperatures, not 70 °C.

## 8. Use in the model

**For B61 (norfuraneol hydroxide-catalysed decay, declared barrier 103.7 kJ/mol): NO barrier evidence.**
Every number is at one temperature (70 °C); the only furanone-loss number is at pH 3.5, where B61's
hydroxide term is ~10⁻³·⁵ of its pH 7 value and irrelevant. It neither supports nor contradicts the
declared 103.7 kJ/mol, and it says nothing about the one-decade-per-pH-unit law (one pH). Its only
bearing is indirect: furaneol at pH 3.5 / 70 °C loses ≤ 20 % in 48 h, which bounds an *acid-side* or
pH-independent furanone loss, and the engine's NF decay at that condition is ~90× slower than that bound
(arithmetic 5). That is a separate question (the frozen neutral `r_nf_decay` and its 143.5 kJ/mol), weak
because of the recovery confound and the dimethyl → monomethyl transfer, and **not** a reason to change
B61's pre-registration.

**Cysteine + furanone sink: YES, qualitatively; magnitude ≈ 2 M⁻¹ min⁻¹ at 70 °C, pH 3.5 (factor ~2-5
uncertainty), and operative only on FREE thiol.** Cys and NAC both consume pre-formed furaneol in
sugar-free buffer with nearly the same k₂ (so it is the thiol, not the amine). In sugar pots the same
thiols at 10 mM do not suppress furaneol at pH 3 and Cys *raises* it at pH 5 and 7, so upstream amine
catalysis by Cys dominates and the direct sink is quenched ≥ 15-50× when 55 mM sugar is present. The engine
has no Cys + NF reaction (NF is consumed only with H2S: `r_nf_mft`, `r_nf_mp3p`, `r_nf_shet`, and by
`ch_red_ox_nf`); a candidate `r_nf_cys` with k₂ ≈ 2 M⁻¹ min⁻¹ at 70 °C would be an **unmeasured-T, unmeasured-pH
constant on a different furanone**, acting on free Cys only — and since the B60 flux map shows free cysteine
is gone within seconds in the pentose pots, its effect there would be small. Its relevance is to fed-NF + Cys
pots and to glutathione pots (Zhou 2023c), where free thiol persists; it is a mechanistic candidate for
the missing MFT/NF sink in Zhou, not a calibrated one.

**Transfer furaneol → norfuraneol: moderate for mechanism, weak for numbers.** Both are
4-hydroxy-3(2H)-furanone enols with the same ring-opening equilibrium (Fig. 3) that the proposed thiol
attack and the base-catalysed ring opening both go through; norfuraneol lacks the 2-methyl, so its
C2 is a CH2 and the ring-opened form differs (a C5 vs C6 tricarbonyl), which plausibly changes rates by
more than the factor-2 scatter here. Use as order-of-magnitude only (± 1 dex).

**Also usable:** Table 1 as a qualitative check that thiols suppress the 1,2-enolisation exit (HMF/MF)
dose-dependently at pH 3-5 while leaving the 2,3-enolisation furanone (with Arg) intact or enhanced — the
same split the engine encodes as `r_pent_tdp` vs `r_pent_dpo`, here for hexose / 6-deoxyhexose.
