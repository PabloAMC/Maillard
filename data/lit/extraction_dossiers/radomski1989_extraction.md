# Radomski & Temeriusz 1989 — EXTRACTION
DOI: 10.1016/0008-6215(89)80005-9 (printed as the PII; Crossref 2026-11-01: "Thiazolidine-4(R)-carboxylic acids derived from sugars: Part I, C-2-epimerisation in aqueous solutions", Radomski & Temeriusz, Carbohydrate Research 187(2), 223-237, 1989-04)
### Fifteen 2-(polyhydroxyalkyl)thiazolidine-4(R)-carboxylic acids from L-cysteine + aldoses / glycolaldehyde / glyceraldehydes, including the D-xylose one (TTCA, compound 10): synthesis yields, 1H / 13C NMR, and C-2 EPIMERISATION in water by polarimetry and NMR. Epimer equilibria at 37 °C (TTCA 2R : 2S = 4.06), relaxation in minutes, solution pH 3.7-4.0. **No formation K, no hydrolysis (release) rate, no temperature series.**

**Read 2026-11-01 from `data/articles/radomski1989.pdf`** (15 pp., Carbohydr. Res. pp. 223-237; Acrobat Capture
OCR, text layer partly garbled). **Tables I-V and Figs. 1, 2, 4, 5 read from 120 dpi renders; Fig. 5 digitised from
a 300 dpi crop** (scratchpad `ing1101/ra-*.png`, `f5.png`; read error +-0.003 in molar fraction, +-0.3 min).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a
figure · **[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived
by me, never printed · **[!]** flag.

**Laboratory.** Department of Chemistry, Warsaw University (Temeriusz); high-field NMR at Southampton. Not in
`lab_ids.yml`; proposed id `warsaw_temeriusz`. Independent of every Maillard group in the corpus.

**Compound key (as printed).** 6 glycolaldehyde; 7 D-, 8 L-glyceraldehyde; 9 L-xylose; **10 D-xylose (= TTCA,
2-(D-xylo-1,2,3,4-tetrahydroxybutyl), the repo's TTCA)**; 11 D-ribose; 12 D-lyxose; 13 D-arabinose; 14
L-arabinose; 15 D-glucose; 16 D-mannose; 17 D-galactose; 18 L-rhamnose; 19 L-fucose; 20 D-fucose.

## 1. What this paper answers for the model

CF-076 asks for K, the epimerisation / ring-opening rate, pH and T of the Cys + pentose ⇌ TTCA exchange.

| model need | what Radomski 1989 gives | tag |
|---|---|---|
| formation K (Cys + Xyl ⇌ TTCA) | **none** | — |
| hydrolysis / release rate k_r | **none** measured. Alkaline decomposition only: "only ~50% of the starting compounds remained after 10 min" at pH 12; "considerable amounts" of hydrolysis products after 5 min at pH 9 (T not stated, compounds not specified) | [M] / [Q] |
| ring opening (ring-chain) rate | **indirect, fast.** C-2 epimerisation needs ring opening: no H-2 deuterium exchange (excludes a C-2 carbanion), acyclic form not seen by NMR ("short lived on the n.m.r. time-scale"). Mutarotation complete in < 1 h (Fig. 1); for 8 (L-glyceraldehyde) the relaxation is **k ≈ 0.20 min⁻¹ (t½ ≈ 3.5 min) at 37 °C** [Z from F, Fig. 5] | [Q] / [F] / [Z] |
| epimer equilibrium, TTCA | **2R : 2S = 4.06** (equilibrium after 40 min), [α]eq -92.0°, saturated D2O / water at 37 °C (Table III) | [M] |
| pH of the solution | "Aqueous solutions (c 0.5) of 6-20 have a pH in the range 3.7-4.0"; Fig. 2 equilibrium at pH 3.80 (14) | [M] |
| pH effect on the ring equilibrium | equilibrium rotation of 14 falls with acidity (pH 3.8, 2.1, 0, 2 M, 4 M HCl; Fig. 2), "due mainly to changes in composition rather than conformation" | [F] / [Q] |
| formation under strong acid | in the synthesis mixture (Cys·HCl 20 mmol + aldehyde 20 mmol in 2-5 mL, [H+] up to 4 M) 13C NMR showed "no traces of product after 20 h"; weak signals after 2 weeks; with base the products crystallise "shortly" | [Q] / [M] |
| temperature | none (37 °C for Table III / Fig. 5; Fig. 1 / 2 temperature not printed) | [!] |

**Answers.**

1. **K: not here** (confidence ~98 %). The paper never measures free cysteine or free sugar at equilibrium.
2. **Ring opening is fast; release is not measured.** At 37 °C and pH ~3.8 the ring opens and recloses (C-2
   epimerisation) with a relaxation of ~0.2 min⁻¹ for 8 [Z]; the 15 compounds of Fig. 1 all equilibrate within
   20-35 min (TTCA's curve included) [F]. Any Ea of 40-80 kJ/mol carries this to ~3-40 min⁻¹ at 100 °C [Z]: ~10⁴-10⁵ x
   faster than CF-076's declared k_r(100 °C, pH 4.5) band 2-8 × 10⁻⁴ min⁻¹. So **ring opening cannot be the slow
   step of TTCA release**, as Fife 1991 (aryl thiazolidines) and Kallen 1971b (formaldehyde TC) already show on
   analogues; this is the first sugar-thiazolidine measurement of it on disk. It **supports the structure** CF-076's
   leading explanation assumes (fast ring-chain tautomerism, slow C-N release), not a number.
3. **The TTCA ring-N pKa, bounded from the solution pH [Z]:** a pure zwitterion at c 0.5 g/100 mL (~20 mM for
   TTCA, M 253) sitting at pH 3.7-4.0 implies, with a carboxyl pKa1 of 1.5-2.2 (TC 1.63, Kallen 1971), a ring-N pKa2
   of **~5.2-6.1** (charge-balance solution of the ampholyte, computed in the env). This is inside CF-076's admitted
   region (feasible up to pKa 5.9) and below TC's 6.24. Caveats: the 3.7-4.0 range is for all of 6-20 together, not
   TTCA alone; temperature unstated (presumably room); no ionic-strength control. P(pKa2 of TTCA in 5.2-6.1 at
   ~25 °C) ~0.6.
4. **Acid limb [Q]:** no product in 20 h with the amine fully protonated (4 M H+), fast crystallisation once base is
   added: formation is slow when the free-base amine fraction is tiny, the attack-limited acid limb of Kallen 1971b.
   Qualitative only.

## 2. Conditions

| item | as printed |
|---|---|
| synthesis | Cys·HCl·H2O 3.61 g (20 mmol) + aldehyde 20 mmol in water 2-5 mL, room temperature, 20 h; then 20 mmol base (method A pyridine + 200 mL ethanol after ~20 h more; method B potassium acetate, no ethanol); yields Table V |
| epimerisation | solids dissolved in water (polarimetry, c 0.25-0.5) or D2O (NMR, saturated); Table III "Saturated solutions at 37°" [V]; Fig. 5 [α]D at 37 °C |
| pH | unbuffered; solutions pH 3.7-4.0; Fig. 2 also 0.5 M phthalate pH 2.1, 1 M, 2 M, 4 M HCl; alkaline runs pH ≥ 8.5 (base not stated for these) |
| O2 | water "thoroughly washed with nitrogen" for syntheses and polarimetry |
| analysis | Perkin-Elmer 241 / Polamat A polarimeters; 1H NMR 90 / 100 / 360 MHz (H-2 doublet ratio for epimers), 13C NMR 22.53 MHz; t.l.c. and g.l.c. for hydrolysis products / residual substrates |

## 3. Quantitative results

**Table III** (90 MHz 1H NMR + polarimetry, saturated solutions, 37 °C; ratio 2R : 2S at the time in parentheses,
min; [α] in degrees) [M]:

| cpd (sugar) | 2R : 2S (min) | [α]eq | [α]2S | [α]2R |
|---|---|---:|---:|---:|
| 6 glycolaldehyde | 0.72 (30) | -150.3 | -178.6 | -108.8 |
| 7 D-glyceraldehyde | 0.96 (no change) | -67.5 | | |
| 8 L-glyceraldehyde | 0.54 (16) | 39.0 [sic; -39 in Fig. 5] | -127.0 | 123.0 |
| 9 L-xylose | 0.58 (40) | -97.0 | -130.2 | -39.8 |
| **10 D-xylose (TTCA)** | **4.06 (40)** | **-92.0** | **-213.7** | **-62.0** |
| 11 D-ribose | 0.77 (41) | -120.5 | -155.4 | -75.5 |
| 12 D-lyxose | 1.40 (42) | -80.5 | -92.8 | -71.6 |
| 13 D-arabinose | 0.59 (36) | -107.5 | -133.2 | -64.2 |
| 14 L-arabinose | 1.92 (35) | -82.3 | -125.3 | -63.5 |
| 15 D-glucose | 0.99 (28) | -95.5 | -117.2 | -73.7 |
| 16 D-mannose | 3.20 (45) | -83.2 | -200.5 | -46.5 |
| 17 D-galactose | 2.83 (30) | -71.0 | -113.8 | -55.9 |
| 18 L-rhamnose | 0.25 (40) | -75.5 | | |
| 19 L-fucose | 0.41 (41) | -97.0 | | |
| 20 D-fucose | 6.10 (30) | -87.0 | | |

H-2 shifts for 10: δ 2(S) 5.12, 2(R) 4.97 ppm; J 3.72 / 5.86 Hz. 13C (Table IV) for 10: ring C-2 67.83 / 66.97,
C-4 65.56 / 65.13, C-5 33.65, COOH 172.35 / 172.42 (major / minor) [M].

**Fig. 5** (8, 37 °C, right axis labelled molar fraction of the 2R epimer, probably the 2R : 2S ratio, sec. 6;
polarimetry and NMR) [F]: ~0.709 (3 min), 0.627 (4), 0.609
(5), 0.589 (7), 0.579 (8), 0.565 (10), 0.554 (13), 0.549 (15), 0.543 (18), 0.541 (20), plateau 0.540. Departures
from 0.540 fall e-fold every ~5 min (5 -> 10 -> 15 min: 0.069 -> 0.025 -> 0.009): **k_relax ≈ 0.20 min⁻¹** [Z].

**Table V** (yields, %, methods A / B) [M]: 10 (TTCA) 82 / syrup; 11 57 / syrup; 15 43 / syrup; 9 72 / -; others
53-93 (A), 20-54 (B). Melting points: 10 152 °C (lit. 158-159 / 152).

**Alkaline** [M/Q]: pH ≥ 8.5 mutarotation "seems to be slower"; pH 9, 5 min: considerable hydrolysis products
(t.l.c.); pH 12: ~50 % remaining at 10 min, i.e. ~0.07 min⁻¹ loss [Z], T unstated.

## 4. Mechanism evidence and its kind

- **Measured equilibrium and relaxation (kind 4-like)**, one laboratory: epimer ratios by two independent methods
  (polarimetry and H-2 NMR, Fig. 5 agree) [M].
- **Ring opening shown by exclusion** [Q]: no deuterium exchange at H-2 rules out C-2 carbanion epimerisation; the
  open chain (Schiff base / iminium) must be the route and is too short-lived for NMR.
- **Configuration depends on work-up** (ethanol gives 2R, its absence 2S, for several sugars), i.e. the isolated
  solid is a kinetic / crystallisation product, not the solution equilibrium [M/Q].

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| CF-076 explanation 1 (in-situ pH, Kallen-shaped K(pH), k_r(pH); ring-N pKa ≤ ~5.9) | **SUPPORTS, weakly** | solution pH 3.7-4.0 implies TTCA ring-N pKa ~5.2-6.1 [Z]; ring opening fast, so the slow step is release, as the explanation assumes |
| CF-076 remedy (1): measure TTCA's ring-N pKa | partly answered (bound, not a titration) | see sec. 1 answer 3 |
| CF-076 explanation 2 (slow exchange with larger K) | LEAVES UNCHANGED | no K |
| CF-076 explanation 3 (catalysed exchange) | LEAVES UNCHANGED | no buffer series |
| CF-091 (13C5-TTCA by difference) | LEAVES UNCHANGED | — |
| analytical note for every TTCA row | **new fact**: TTCA is two C-2 epimers, 2R : 2S ≈ 4 at 37 °C, interconverting in minutes; HPLC "TTCA" is their sum unless separated | [M] |

## 6. Disagreements found

1. **TTCA's "natural pH" (Zhai 2021, `zhai2021_extraction.md` sec. 2: 10 mM TTCA "the natural pH of 5.5") vs
   Radomski's 3.7-4.0 for the same class of compounds at ~20 mM.** For a pure zwitterion the solution pH barely
   depends on concentration (computed: 3.7-4.5 at 10 mM for any pKa1 1.5-2.2, pKa2 5.2-6.5), so pH 5.5 cannot be
   pure TTCA [Z]. Candidate causes: Zhai's TTCA (made at pH 7.4, purified "98 %") is partly a sodium salt or carries
   base; electrode / ionic strength; or a different compound composition. This bears directly on CF-076 (the
   dilute pots' in-situ pH is the explanation's key variable) and on the lesson about Zhai's pH labels.
2. **Internal, Table III vs Fig. 5 for compound 8**: Table III prints the 2R : 2S RATIO 0.54; Fig. 5's right axis
   ("molar fraction of 2(R),4(R) epimer") plateaus at 0.540 (a ratio of 1.17). Table III is self-consistent if its
   [α]eq "39.0" lost a minus sign (Fig. 5's left axis plateaus at -39°): the lever rule between [α]2S -127.0 and
   [α]2R +123.0 then gives a 2R fraction of 0.35, i.e. a ratio of 0.54 [Z]. So Fig. 5's right axis most likely
   plots the 2R : 2S ratio under a "molar fraction" label [!]. The relaxation constant (sec. 3) barely depends on
   which (ratio and fraction are near-linear over 0.54-0.71).

## 7. Ledger-ready measurements

1. **Table III, compound 10 (TTCA): 2R : 2S = 4.06 at equilibrium, 37 °C**, and the other pentoses (9, 11-14): the
   epimer composition of the analyte (two readers).
2. **Fig. 5 relaxation, k ≈ 0.20 min⁻¹ at 37 °C** (compound 8, [F] -> [Z]): a lower bound on ring-chain exchange,
   transfer to TTCA assumed.
3. **Solution pH 3.7-4.0 at c 0.5 (6-20)**: an input for a TTCA pKa bound, not itself a pKa.
4. The alkaline decomposition sentence (pH 12, ~50 % in 10 min) only as a qualitative row: T and compound unstated.
