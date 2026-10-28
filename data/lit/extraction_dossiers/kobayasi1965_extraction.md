# Kobayasi & Fujimaki 1965 — EXTRACTION (cysteine + ninhydrin 10 + 10 mmol in 50 mL water, steam-distilled 100 °C, 1 h → mercaptoacetaldehyde isolated as its bis-2,4-DNPH; cysteine + glucosamine 10 + 10 / 10 + 20 mmol, 100 °C, 1 h → 0.42 / 0.80 mmol H2S; a two-page short communication, two numbers, no rates)

### THE PRIMARY SOURCE FOR "CYSTEINE + α-DICARBONYL → MERCAPTOACETALDEHYDE + H2S + ACETALDEHYDE", AND IT SAYS SOMETHING THE MODEL'S COMMENTARY DOES NOT: the H2S is proposed to come from a REDUCTIVE cleavage of the C–S bond by the dicarbonyl's enaminol, not from a spontaneous breakdown of mercaptoacetaldehyde. Mercaptoacetaldehyde → H2S + acetaldehyde is a two-electron reduction (C2H4OS + 2 H → C2H4O + H2S); it needs a reductant.

**Source on disk:** `data/articles/kobayasi1965.pdf` (3 pp.: a Taylor & Francis cover page + journal pp. 698-699).
Read from 150-dpi page renders through the Read tool; the passage carrying the numbers (p. 698, right column)
was re-checked on a 250-dpi crop (scratchpad `img/kob_crop.png`). Text layer (`pdftotext -layout`) agrees.
There are no tables; Fig. 1 is an IR spectrum (identity only); Charts 1-2 are schemes.

## 0. Identity

| field | value |
|---|---|
| Title | "On the Formation of Mercaptoacetaldehyde, Hydrogen Sulfide and Acetaldehyde on Boiling Cysteine with Carbonyl Compounds" ✔ |
| Authors | Noboru Kobayasi, Masao Fujimaki — Department of Agricultural Chemistry, Faculty of Agriculture, The University of Tokyo ✔ |
| Venue | Agricultural and Biological Chemistry **29**(7), 698-699 (1965), "Short Communication" (letter, opens "Sir:"); received 17 May 1965 ✔ |
| DOI as printed | **10.1080/00021369.1965.10858450** (Taylor & Francis cover page). The brief's 10.1271/bbb1961.29.698 is the J-STAGE DOI for the same article (same journal, volume, issue, pages). **Same paper; two registrar DOIs.** |
| Paper type | Letter: qualitative identification + two H2S masses + one isolated-derivative mass. No kinetics, no pH, no time course. |
| Cited by | Wedzicha & Edwards 1991 (`wedzicha1991_extraction.md`), who note the 3-deoxyglucosone / dehydroascorbic acid claim is made "without supporting evidence" |

## 1. Why it matters

Gap (a) in `docs/guides/SULFUR_FLUX_MAP.md` §5-6: `r_dpo_ddp` (Strecker of the 1-deoxyosone with cysteine)
routes the cysteine sulfur, as mercaptoacetaldehyde (MAA), to the inert `FRAG_S`, and the flux map calls
"mercaptoacetaldehyde → H2S + acetaldehyde" the textbook route that would rescue the Hofmann pH 7 arm.
This is the paper the textbook route traces to. What it actually shows:

1. **MAA is formed** from cysteine + ninhydrin at 100 °C and can be isolated (as the disulfide's bis-DNPH):
   ≥ 3.6 mol % of the cysteine recovered in the steam distillate in 1 h (arithmetic 1).
2. **H2S and acetaldehyde are formed** in the same system (GC of the headspace; lead sulfide; DNPH) — qualitative only.
3. **The proposed H2S mechanism is reductive** (Chart 1, "Second step"): the Strecker by-product of the
   α-diketone is an **enaminol** "which possesses strong reducing power"; "this enaminol fulfills the reductive
   cleavage of the C-S bond of mercaptoacetaldehyde or cysteine" (p. 698). Chart 1 writes R–SH + 2H → R–H + H2S,
   the enaminol going to the ketimine (–2H).
4. **A test of the reductant hypothesis**: glucosamine (an amino sugar that tautomerises partly to an enaminol, Chart 2)
   boiled with cysteine gives H2S **in proportion to the glucosamine charge** — 0.42 mmol with 10 mmol and
   0.80 mmol with 20 mmol, from 10 mmol cysteine, 100 °C, 1 h. **Ascorbic acid**, "a typical aci-endiol
   reductone", gives H2S "in a very small amount" (not quantified).

## 2. Methods as they matter to a model

- **MAA isolation.** 1.76 g (10 mM) cysteine-HCl + 1.78 g (10 mM) ninhydrin in 50 mL water,
  "steam-distilled at 100 °C for an hour" (p. 698). In this paper **"mM" means millimoles**
  (1.76 g cysteine·HCl = 10.0 mmol; 1.78 g ninhydrin = 10.0 mmol), so the charge is **0.20 mol/L each in 50 mL**.
  pH not stated (cysteine hydrochloride, unbuffered: strongly acidic, likely pH ≤ 2 initially — my inference).
  Distillate treated with iodine "under acetic acid-acidity" (oxidising MAA to the disulfide,
  "dithioglycollic aldehyde"), then 2,4-DNPH: **92 mg** hydrazone precipitate, recrystallised from nitrobenzene
  ×3, m.p. 205.5-206.0 °C, C/H/N/S within 0.6 % of C16H14N8O8S2.
- **Headspace.** Equimolar cysteine + ninhydrin boiled at 100 °C: "The hydrogen sulfide-like aroma evolved
  immediately"; GC of the surface vapour showed H2S, acetaldehyde and MAA. No quantities.
- **Glucosamine series.** 1.21 g (10 mM) cysteine gave 14.3 mg (0.42 mM) and 27.2 mg (0.80 mM) H2S with
  1.79 g (10 mM) and 3.58 g (20 mM) glucosamine, "at 100 °C for an hour" (p. 698-699). Free cysteine (1.21 g = 10.0 mmol) and glucosamine as free base (1.79 g = 10.0 mmol;
  MW 179.2). **Volume NOT stated** (if the 50 mL of the ninhydrin run: 0.20 / 0.40 mol/L). **pH NOT stated**
  (cysteine + free-base glucosamine, unbuffered: plausibly near-neutral to mildly basic; my inference, not reported).
  Open or closed vessel NOT stated ("boiled").
- **H2S quantitation.** "determined with gas chromatography, gravimetry as lead sulfide and iodine titration"
  (p. 699) — which method gave the printed masses is NOT stated; no replicates, no error.
- **Not reported anywhere:** a cysteine-only blank, H2S from the ninhydrin run, acetaldehyde amounts, MAA yield as
  such, any time course, any temperature other than 100 °C, any pH, the 3-deoxyglucosone or dehydroascorbic-acid
  data (asserted in the introduction, never shown).

## 3. Numbers, re-typed (all from p. 698-699 text)

| system | charge | conditions | observable | value as printed | mol % of cysteine (mine) |
|---|---|---|---|---|---|
| Cys·HCl + ninhydrin | 10 + 10 mmol in 50 mL water | steam distillation, 100 °C, 1 h | bis-2,4-DNPH of oxidised MAA, crude precipitate | **92 mg** | **≥ 3.6 % as MAA** (arith. 1) |
| Cys + glucosamine | 10 + 10 mmol, volume n.r. | boiled, 100 °C, 1 h | H2S | **14.3 mg (0.42 mmol)** | **4.2 %** |
| Cys + glucosamine | 10 + 20 mmol, volume n.r. | boiled, 100 °C, 1 h | H2S | **27.2 mg (0.80 mmol)** | **8.0 %** |
| Cys + ascorbic acid | n.r. | "under the same condition" | H2S | "a very small amount" | n.r. |
| cystine / methionine / protein + glucosamine | n.r. | boiling | H2S | "also produced" | n.r. |

### Arithmetic (mine)

**1. MAA recovered.** MW of C16H14N8O8S2 = 510.5 g/mol; 92 mg = 0.180 mmol of bis-hydrazone = **0.360 mmol MAA**
(two MAA units per disulfide) = **3.6 mol % of the 10 mmol cysteine**. A **floor**: only what steam-distilled in
1 h, crude precipitate before three recrystallisations, iodine oxidation and hydrazone precipitation not
quantitative, and MAA dimerises to 2,5-dihydroxy-1,4-dithiane in solution (Wedzicha 1991, p. 85).

**2. H2S per glucosamine.** 0.42 → 0.80 mmol on doubling glucosamine: 1.9× — **close to first order in the
reductant**, as the reductive-cleavage proposal predicts. One pair of points, no error.

**3. Against neat thermolysis.** The engine's measured cysteine thermolysis (`cysteine_thermolysis_k`,
Zheng & Ho 1994) at 100 °C gives, over 1 h, **0.17 % H2S at pH 5 and 0.91 % at pH 7**. The glucosamine pots give
**4.2 % and 8.0 %: 4.6-47× more** depending on the unreported pH. Even at the most generous pH (7), the
carbonyl/enaminol partner raises H2S ≥ 4.6×. (Caveat: pH 9 not checked; Zheng's pH-9 rate at 100 °C is higher still.
The comparison is order-of-magnitude.)

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | evidence class |
|---|---|---|---|---|---|
| Cys + α-dicarbonyl → MAA | MAA recovered in distillate | **≥ 3.6** | mol % of Cys | ninhydrin, 0.2 M each, 100 °C, 1 h, pH n.r. (acidic) | **derived floor** (isolated derivative) |
| Cys + enaminol (glucosamine) → H2S | H2S yield | **4.2 / 8.0** | mol % of Cys | 1 / 2 eq glucosamine, 100 °C, 1 h, pH and volume n.r. | **level_only**, single determinations |
| H2S vs reductant | order in glucosamine | **~1** (1.9× on 2×) | — | same | within-study ratio, n = 2 |
| MAA → H2S + acetaldehyde | rate, yield | **NOT REPORTED** | — | — | qualitative only (products co-detected by GC) |

No rate constant, no Ea, no pH dependence.

## 5. What this can and cannot anchor

**Can:**
- A **qualitative licence** for a Strecker sulfur split: MAA, H2S and acetaldehyde co-occur when cysteine meets
  an α-dicarbonyl at 100 °C.
- A **mechanistic constraint**: the MAA → H2S + acetaldehyde step is a reduction, so its rate should scale with a
  reductant pool (reductones/enaminols, which the Strecker step itself makes 1:1), not be a unimolecular decay.
  A model that writes `MAA → H2S + CH3CHO` as first order is not mass-balanced in hydrogen; a 2 e⁻ donor is
  implied.
- A **one-sided magnitude**: H2S from cysteine with an enaminol partner at 100 °C reaches 4-8 mol % of the
  cysteine in 1 h.

**Cannot:**
- No split between MAA → H2S and MAA → other fates (thiophenes, 1,4-dithiane, thiazolines).
- No rate, no temperature dependence, no pH.
- Glucosamine is not a deoxyosone; the H2S may also come from reductive cleavage of **cysteine itself**
  (Chart 1 allows either substrate), so the glucosamine numbers are not an MAA → H2S yield.
- The ninhydrin system (a tricarbonyl, a very fast Strecker reagent) is not a sugar dicarbonyl.

## 6. Flags

1. **Two DOIs for one paper** (T&F 10.1080/00021369.1965.10858450 printed; J-STAGE 10.1271/bbb1961.29.698 in the brief).
2. **"mM" means mmol** throughout; reading it as mmol/L would misstate every charge by 50×.
3. **No blanks.** Without a cysteine-alone control the glucosamine H2S cannot be split into thermolysis + partner.
4. **pH and volume of the glucosamine runs are not reported**; "boiled" (open/closed unknown) — H2S losses to air possible,
   so the yields are floors if the vessel was open, and the trapping method is not described.
5. **The 3-DG / dehydroascorbic acid claim is unsupported in this paper** (no data shown), as Wedzicha & Edwards 1991 note.
6. Registry: paper absent from `data/keys/papers.yml` (not checked for compounds; MAA is not an engine species).

## Use in the model

- **Gap (a), Strecker sulfur.** Supports making the Strecker sulfur split explicit rather than sending MAA to inert
  `FRAG_S` — but it does **not** give a split fraction. The quantitative content is a floor (≥ 3.6 mol % MAA from
  ninhydrin in 1 h at 100 °C) and a magnitude (4-8 mol % H2S with an enaminol partner at 100 °C, 1 h). The
  mechanistic content matters more: if B59 adds `MAA → H2S + acetaldehyde`, write it as **reductant-dependent**
  (second order with the Strecker co-product, or with the reductone pool) and declare its constant an assumption;
  this paper licenses the route's existence and its reductant dependence, not its rate.
- **Contradiction with the flux map's wording:** §5/§6 describe MAA → H2S + acetaldehyde as a plain decomposition.
  Kobayasi's own scheme makes it a reduction. That changes the pH-7 rescue argument: the H2S supply would scale
  with the reductone produced by the same Strecker event (1:1), so it cannot exceed the Strecker flux and is
  bounded by reductant survival.
- **Caveats:** ninhydrin and glucosamine stand in for sugar dicarbonyls; unbuffered, pH unknown; 100 °C only;
  single determinations.
