# US 10,076,126 B2, Nestec, 2018 — EXTRACTION (xylose + L-cysteine ~0.25 mol/kg each in aqueous phosphate pH 5.5, glycerol, or choline chloride/glycerol 1:1; crimped glass vials, oil bath 110 °C, 0.5 / 1 / 2 / 4 h; MFT and FFT by HS-SPME-GC-MS with deuterated internal standards; Figs. 1-2 digitised)

### A SECOND HOUR-SCALE MFT / FFT TIME COURSE IN WATER FROM CYSTEINE + PENTOSE — four points, 110 °C, 0.5-4 h. MFT in the buffer arm rises with upward curvature through 4 h (≈15 → 124 → 575 → 1855 ng/g), no plateau. FFT in buffer stays at the axis until 1 h and reaches only ≈520 ng/g at 4 h (MFT/FFT ≈ 3.6 at 4 h). BUT the labelled standards were added BEFORE heating, so the printed amounts are ratios to a standard that itself degraded at 110 °C for up to 4 h (see §6 flag 1): the upward curvature is partly what that artefact produces.

**Source on disk:** `data/articles/US10076126.pdf` (11 pp.: front page, four drawing sheets, specification cols 1-12).
Text read with `pdftotext -layout -q` (prose reliable; figure pages garbage). Example 1 and the figure legends re-checked on
300-dpi renders (`[2H3]`/`[2H2]` labels and the ◇ / □ / ▲ legend do not survive the text layer). Figs. 1 and 2 (sheets 1-2,
PDF pp. 2-3) digitised by eye from 200-dpi page renders and 400-dpi crops around each buffer-arm marker
(scratchpad `img/us126-02.png`, `us126-03.png`, `f1_*.png`, `f2_*.png`); Fig. 4 (p. 5) at 150 dpi.

**Status of the evidence.** A patent, not peer reviewed; examples are written to support claims (here: the deep eutectic
solvent beats water and glycerol). No replicate count, no error bars, no standard amounts, no recovery data are printed.

## 0. Identity

| field | value |
|---|---|
| Title | "Process for preparing a flavour composition" ✔ |
| Patent number | **US 10,076,126 B2** (kind B2; 13 claims, 4 drawing sheets) |
| Applicant / assignee | Nestec S.A., Vevey (CH) |
| Inventors | Karin Kraehenbuehl (Fully); Tomas Davidek (Correvon); Florian Viton (Lausanne); Heiko Oertling (Lausanne); Rachid Belrhlid (Savigny) — all CH |
| Appl. No. | 15/023,836; PCT/EP2014/070247, PCT filed **Sep. 23, 2014**; §371 date Mar. 22, 2016 |
| PCT publication | WO2015/044136, Apr. 2, 2015 |
| Priority (as printed) | **EP 13185655(.1), Sep. 24, 2013** |
| Prior publication | US 2016/0242442 A1, Aug. 25, 2016 |
| Date of patent (grant) | **Sep. 18, 2018** |
| Related literature by the same group | Cerny et al. 2006 (glycerol in Maillard, cited by the examiner); Belrhlid is also an author of `belrhlid2002_extraction.md` |

## 1. Why it matters

The model has one hour-scale MFT/FFT time course in water from cysteine + pentose: Schieberle 2000
(`schieberle2000_extraction.md`, ribose + cysteine, 100 °C, 30-720 min), in which MFT keeps rising for 12 h while the engine
stops forming it within an hour. This patent gives a second, independent laboratory's pot (Nestlé Research, xylose not
ribose, 110 °C, 0.5-4 h) where **MFT in buffer is still accelerating at 4 h**. It is the same qualitative message —
MFT formation in water is slow and sustained over hours, not a burst — from a different sugar, temperature and lab.

## 2. Methods as they matter to a model (Example 1, col. 8 l. 50 – col. 9 l. 30)

- **Charge.** "For each solvent, 75 mg of xylose (Fluka) and 60 mg of L-cysteine (Fluka) were weighed into 2.0 g of the
  solvent". Arithmetic (mine; MW xylose 150.13, L-cysteine 121.16):
  - xylose 75 mg = **0.4996 mmol**; cysteine 60 mg = **0.4952 mmol** → **molar ratio 1.01 : 1**.
  - per kg solvent: xylose **0.250 mol/kg**, cysteine **0.248 mol/kg**.
  - per kg reaction mixture (2.135 g, the basis of the figures' "ng/g reaction mixture"): xylose **0.234 mol/kg**,
    cysteine **0.232 mol/kg**. ≈ 0.25 mol/L in water (density ≈ 1.0; buffer salts and solute volume ignored).
- **Solvents.** A) choline chloride : glycerol 1:1 molar; B) glycerol; C) **"an aqueous phosphate buffer at pH 5.5"**.
  ⚠ **The buffer concentration is NOT printed in this patent** (the 0.2 mol/L sodium phosphate pH 5.5 figure belongs to the
  sibling patent US 8,920,862, Example 5; same group, plausibly the same stock, not stated). Counter-ion not stated.
- **pH.** 5.5 is the buffer's nominal pH before adding precursors. **pH after adding 0.25 mol/kg cysteine not reported;
  pH after heating not measured / not reported.**
- **Internal standards.** [2H3]-2-methyl-3-furanthiol (aromaLAB AG, no. 2119) and [2H2]-2-furfurylthiol (aromaLAB AG,
  no. 2032) added **into the vials together with the precursors, before heating**. Amounts NOT printed.
- **Vessel.** "a series of glass vials", crimped (sealed). Vial volume NOT printed (the SPME was done "in the glass vials",
  so headspace vials, plausibly 10-20 mL — my inference). One vial per time point implied; replicates not stated.
- **Mixing.** Vortexed at room temperature until homogeneous.
- **Heating.** "heated to 110° C. in an oil bath for periods between 30 minutes and 4 hours", then "cooled in a crushed ice
  bath and stored at -20° C. until analyzed". **Come-up time not stated**; for ~2 g in a glass vial dropped in a 110 °C
  oil bath it is of order 2-5 min (my estimate, not measured), material only for the 0.5 h point. Times are bath times.
- **Analysis.** HS-SPME (PDMS/DVB fibre, Supelco 57345-U, "75 um" as printed), 15 min equilibration at 40 °C with agitation
  (5 min pre-incubation), desorption 5 min at 250 °C splitless; GC 6890A / MSD 5973, HP-5MS 30 m × 0.25 mm × 0.25 µm;
  oven 25 °C (5 min) → 6 °C/min → 180 °C → 10 °C/min → 240 °C (5 min); EI 70 eV, m/z 29-350. Quantification by the ratio
  to the deuterated standards (stable isotope dilution in principle; ions not printed here — the sibling patent uses
  m/z 114 vs 117 for MFT and 114 vs 116 for FFT).
- **Not reported:** standard amounts, calibration response factors, replicate number, error, xylose/cysteine consumption,
  H2S, any other compound in Example 1, pH drift, browning.

## 3. Numbers, re-typed

### 3a. Fig. 1 — MFT (ng/g reaction mixture) vs time (h); legend (col. 2): aqueous buffer ◇, glycerol □, choline chloride/glycerol ▲

**Axis calibration.** y: gridlines every 500 ng/g, 0-3500, linear; at 400 dpi the 0 and 500 ng/g gridlines are 250 px apart
(≈ 2.0 ng/g per px). x: linear, 0-4 h; markers sit at 0.5, 1, 2 and 4 h (no other times plotted).
**Read-off uncertainty:** marker-centre judgement ±5 px at 400 dpi → **±10 ng/g** at the low points, **±15-20 ng/g**
at 2-4 h (from the gridline interpolation). No error bars are drawn; these are read-off errors only, not experimental error.

| time (h) | **MFT, aqueous buffer ◇ (ng/g mixture)** | ± read-off | µmol/kg mixture | µmol/kg solvent | µmol/mol xylose (= per mol Cys within 1 %) |
|---|---|---|---|---|---|
| 0.5 | **15** | 10 | 0.13 | 0.14 | 0.56 |
| 1 | **124** | 10 | 1.09 | 1.16 | 4.6 |
| 2 | **575** | 15 | 5.04 | 5.38 | 21.5 |
| 4 | **1855** | 20 | 16.25 | 17.35 | 69.4 |

Conversion: ng/g ÷ 114.17 g/mol (C5H6OS) = µmol/kg mixture; × 2.135/2.000 = per kg solvent; ÷ 0.234 mol/kg = per mol xylose.

Other arms (coarser read, ±50 ng/g, for context only):

| time (h) | glycerol □ | ChCl/glycerol ▲ |
|---|---|---|
| 0.5 | ≈130 | ≈280 |
| 1 | ≈430 | ≈680 |
| 2 | ≈1480 | ≈1810 |
| 4 | ≈2950 | ≈3210 |

### 3b. Fig. 2 — FFT (ng/g reaction mixture) vs time (h); same legend

**Axis calibration.** y: gridlines every 5000 ng/g, 0-30000; at 400 dpi the 0 and 5000 gridlines are ≈290 px apart
(≈ 17 ng/g per px). The buffer arm sits on or just above the zero line, so the axis is 10× coarser than Fig. 1.
**Read-off uncertainty ±90 ng/g per point** (±5 px), i.e. the 0.5 and 1 h points are only an upper bound.

| time (h) | **FFT, aqueous buffer ◇ (ng/g mixture)** | ± read-off | µmol/kg mixture | µmol/kg solvent | µmol/mol xylose |
|---|---|---|---|---|---|
| 0.5 | **≈0 (< 90)** | — | < 0.8 | < 0.85 | < 3.4 |
| 1 | **≈0 (< 90)** | — | < 0.8 | < 0.85 | < 3.4 |
| 2 | **150** | 90 | 1.3 | 1.4 | 5.6 |
| 4 | **520** | 90 | 4.6 | 4.9 | 19.5 |

Other arms (±400 ng/g): glycerol □ ≈100 / 2600 / 7000 / 29200; ChCl/glycerol ▲ ≈400 / 5800 / 7900 / 16300 at 0.5 / 1 / 2 / 4 h.

### 3c. Fig. 4 (Example 4, for completeness) — residual rhamnose after 5 min at 120 °C, closed vessel

Precursor mix (0.57 g rhamnose, 1.69 g fructose, 1.18 g lysine·HCl, 0.47 g glycine, 0.106 g Na2HPO4, 2 g water; 0.547 g of it)
+ 0.8 g solvent + 3.2 g sunflower oil. Read (±2 %): **water 88 %**, betaine/glycerol (A1) 82 %, ChCl/glycerol (A2) 60 %.
Rhamnose by HPAEC (CarboPac PA1, 20 mM OH−). One bar each, no error.

### Arithmetic (mine)

1. **Shape.** MFT buffer increments: 0.5→1 h +109 ng/g (218 ng/g/h); 1→2 h +451 (451 ng/g/h); 2→4 h +1280 (640 ng/g/h).
   The apparent rate rises ~3× over the run: a lag/induction followed by sustained formation, no sign of a maximum by 4 h.
   Log-log slope 1→4 h: ln(1855/124)/ln 4 = **1.95** (t² -like).
2. **Yield.** 69 µmol MFT per mol cysteine at 4 h = **0.007 mol %** of the cysteine; sulfur balance is not constrained by MFT.
3. **MFT/FFT in buffer at 4 h:** 1855/520 = **3.6** (range 2.9-4.4 on read-off). In glycerol it is ≈0.10, in ChCl/glycerol ≈0.20.
4. **Against the sibling patent** (US 8,920,862, buffer arm, 95 °C, 6 h, xylose:Cys 3:1, 0.48 mol/kg Cys): 28 µmol MFT/mol Cys
   and MFT/FFT = 0.50. Here (110 °C, 4 h, 1:1, 0.23 mol/kg): 70 µmol/mol Cys and MFT/FFT = 3.6. **The MFT/FFT ratio in buffer
   disagrees by 7× between the two Nestlé patents** — different temperature, ratio and concentration, and different
   standard-addition timing (§6 flag 1), so no single cause can be assigned.

## 4. Kinetic numbers the repository can use

Grades: **M** = printed number; **F** = figure-digitised (read-off error stated); **Z** = not reported / zero information;
**!** = usable only with the flagged caveat.

| # | observable | value | unit | conditions | grade | source |
|---|---|---|---|---|---|---|
| 1 | charge, xylose / L-cysteine | 0.250 / 0.248 (0.234 / 0.232 per kg mixture) | mol/kg solvent | 75 / 60 mg in 2.0 g | M (masses) | col. 8 |
| 2 | MFT, buffer, 0.5 h | 15 ± 10 | ng/g mixture | 110 °C, pH 5.5 nominal, crimped vial, oil bath | F ! | Fig. 1 |
| 3 | MFT, buffer, 1 h | 124 ± 10 | ng/g mixture | " | F ! | Fig. 1 |
| 4 | MFT, buffer, 2 h | 575 ± 15 | ng/g mixture | " | F ! | Fig. 1 |
| 5 | MFT, buffer, 4 h | 1855 ± 20 | ng/g mixture | " | F ! | Fig. 1 |
| 6 | FFT, buffer, 0.5 h / 1 h | < 90 / < 90 | ng/g mixture | " | F ! (upper bounds) | Fig. 2 |
| 7 | FFT, buffer, 2 h | 150 ± 90 | ng/g mixture | " | F ! | Fig. 2 |
| 8 | FFT, buffer, 4 h | 520 ± 90 | ng/g mixture | " | F ! | Fig. 2 |
| 9 | within-study ratios MFT(4 h)/MFT(1 h), MFT(2 h)/MFT(1 h) | 15.0 / 4.6 | — | " | F ! | Fig. 1 |
| 10 | MFT/FFT, buffer, 4 h | 3.6 (2.9-4.4) | — | " | F ! | Figs. 1-2 |
| 11 | glycerol and ChCl/glycerol arms | §3a-3b | ng/g | 110 °C, non-aqueous | F (not an aqueous pot) | Figs. 1-2 |
| 12 | residual rhamnose, water arm | 88 ± 2 | % | 120 °C, 5 min, oil + water, closed | F | Fig. 4 |
| 13 | buffer concentration, final pH, standard amounts, replicates, xylose/Cys consumption, H2S | — | — | — | Z | — |

## 5. What this can and cannot anchor

**Can:**
- A **second hour-scale aqueous MFT time course from cysteine + pentose**, from a different lab and sugar than Schieberle 2000:
  formation continues and accelerates over 0.5-4 h at 110 °C. A model that stops forming MFT within the first hour
  fails this pot on shape, whatever the absolute level.
- **Within-study ratios** (rows 9-10) are the safest use: they are independent of the unprinted standard amounts and of
  matrix response, but not of the standard-degradation artefact (flag 1).
- A **lag** (MFT at 0.5 h ≈ 1 % of the 4 h value; FFT below detection on this axis until 1 h).

**Cannot:**
- No replicate or experimental error; the stated ± are read-off only. Treat experimental CV as unknown (SIDA-SPME
  typically 5-15 %; my prior, not the source's).
- No pH trajectory (cysteine at 0.25 mol/kg in a phosphate buffer of unstated strength can shift pH; not measured).
- Absolute levels depend on the standard-addition-before-heating design (flag 1) and are not directly comparable with
  Schieberle 2000 or with US 8,920,862 (standards added after heating).

## 6. Flags

1. **! Internal standards added BEFORE heating.** [2H3]-MFT and [2H2]-FFT were heated with the reaction for up to 4 h at 110 °C.
   Thiols are not stable there (oxidation to disulfides, addition to carbonyls). The reported amount is
   native(t) × IS0 / IS(t). If both the native and the labelled thiol are lost first-order with the same k and native is
   formed at a steady rate F, reported(t) = (F/k)(e^{kt} − 1): this exceeds both the true amount present, (F/k)(1 − e^{−kt}),
   and the cumulative formation F·t, and it **curves upward**. So the t²-like shape (arithmetic 1) can be partly an
   artefact; the data read most safely as "loss-corrected cumulative formation, possibly over-corrected". Isotope effects
   and label exchange (methyl-d3 and CH2-d2 labels are non-exchangeable in principle) would not remove this.
2. **Buffer concentration not printed** (the brief's 0.2 M is from US 8,920,862). Nominal pH 5.5 before precursors;
   final pH not reported.
3. **Vial size, headspace, come-up time, replicates not printed.** Times are bath residence times.
4. **The FFT axis is 10× coarser than MFT's;** the buffer FFT 0.5 and 1 h points are upper bounds, and the 2 h point is
   1.7σ of read-off from zero.
5. **Patent selection.** The example exists to show the DES advantage; the aqueous arm is the comparator.
6. **SPME fibre "75 um"** in Example 1 but the same Supelco catalogue number (57345-U) is "65 um" in Example 2 and in the
   sibling patent; almost certainly 65 µm PDMS/DVB (typo); immaterial with isotope standards.
7. Registry: patent not checked against `data/keys/papers.yml` (no edit made).

## Use in the model

- **Hold-out time course: yes, as a shape test, graded F !.** Score the model's MFT(t) for xylose 0.25 + cysteine 0.25 mol/kg,
  pH 5.5 phosphate, 110 °C, sealed, at 0.5 / 1 / 2 / 4 h on the **within-study ratios** MFT(t)/MFT(4 h)
  = 0.008 / 0.067 / 0.31 / 1 (read-off ±0.006 / 0.006 / 0.01 / —), and on FFT(4 h)/FFT(2 h) ≈ 3.5 (wide: 1.8-10 on read-off).
  Directional claim that is robust to flag 1: **MFT(4 h) > MFT(2 h) > MFT(1 h) in water at 110 °C** — no plateau within 4 h.
  Even if flag 1 inflates late points, the artefact factor e^{kt} cannot create a 15× rise from 1 to 4 h unless the
  standard lost > 90 % (k > 0.9 h⁻¹); a monotone rise survives any plausible k.
- **Absolute levels:** use as an order-of-magnitude level (≈ 2 µg/g MFT, ≈ 0.5 µg/g FFT at 4 h) with a declared ×2 band
  for flag 1 and the unknown buffer strength; do not fit to them.
- **MFT/FFT ratio:** report, do not fit; it contradicts the sibling patent's buffer arm by 7× (arithmetic 4).
- Pair with Schieberle 2000 (ribose, 100 °C, 0.5-12 h): two pots, two sugars, two temperatures, same qualitative behaviour.
