# Bork, Stobernack, Rohn & Kanzler 2024 (Food Chem. 449:139189) — EXTRACTION (file stem `bork2024`)
### *"Browning reactions of hydroxycinnamic acids and heterocyclic Maillard reaction intermediates – Formation of phenol-containing colorants"*
### Caffeic acid (CA) or ferulic acid (FA), 0.05 mmol, heated **neat (no water, no buffer)** alone or 1:1 with furfural (FF), 5-hydroxymethylfurfural (HMF) or pyrrole-2-carbaldehyde (PA) in **sealed glass vials at 220 °C, 0-30 min**; reactants by HPLC-DAD **relative to t = 0**; three isolated products (DHMF, FAHMF, DVG) by external calibration; colour index, pH, TEAC, APCI-HRMS. **No norfuraneol (NF) anywhere in the experiments**: NF appears once, as a citation of `bork2022` in the introduction. A dry-roasting system outside the engine's aqueous scope.

**Read 2026-11-05 (first reader) from `data/articles/bork2024.pdf`** (15 pp., Elsevier open-access PDF, clean text layer).
Prose from `pdftotext -raw`; figures from 130 dpi page renders and **300 dpi crops of Fig. 3B-D (p. 6)** (scratchpad
`sweep_img/b24f3b-06.png`, `b24f3cd-06.png`). Bars read by eye against the axis ticks: **read error ±2 % of scale**
(conversion), **±0.05 pH unit** (pH). Supplementary material (Figs. S-1 to S-8, Table S-1; `bork2024_supplementary.pdf`,
7 pp.) read the same day into **sec. 8**.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Leon V. Bork** (corresponding, l.bork@tu-berlin.de), Tobias Stobernack (BfR), Sascha Rohn, Clemens Kanzler |
| Affiliation | TU Berlin, Institute of Food Technology and Food Chemistry, Dept. of Food Chemistry and Analysis; BfR Berlin, Dept. of Chemical and Product Safety |
| Venue | *Food Chemistry* **449** (2024) 139189 |
| DOI | **10.1016/j.foodchem.2024.139189** (printed p. 1 and in the SI note) — matches the brief and the PDF's title metadata |
| Dates | received 2 Feb 2024; revised 26 Mar 2024; accepted 28 Mar 2024; online 2 Apr 2024 |

**Laboratory.** Kanzler group, TU Berlin — **the same laboratory as `bork2022`, `bork2022b`, `bork2023`** (proposed
`tu_berlin_kanzler`, not yet in `data/keys/laboratories.yml`). The BfR co-author ran HRMS [Z, as in `bork2023`].

## 1. What this paper answers for the model

| model need | what Bork 2024 gives | tag |
|---|---|---|
| NF decay (E16) | **nothing**: no NF pot, no NF measurement | — |
| FF / HMF / PA stability | dry, 220 °C: **HMF alone falls to ~2 % in 30 min; FF and PA alone do not decay in 10 min** (Fig. 3B, S-4B) | [F] |
| phenol + heterocycle cross reactions | HMF consumed to ~60 % (FA) / ~45 % (CA); PA to ~36 % with CA; FF hardly | [F] / [M] |
| any aqueous, buffered or pH-controlled pot | **none**: neat solids; the "pH value" is of an extract whose solvent is not stated | [!] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | 0.05 mmol each: CA 9.0 mg, FA 9.7 mg, FF 4.8 mg, HMF 6.3 mg, PA 4.8 mg; heated individually and as equimolar binary phenol + carbonyl mixtures. **No solvent** |
| Medium / buffer | none (dry); no water added |
| pH | **not set**; a "pH value" was "monitored" at 0, 2.5 and 10 min (Fig. 3C) — **method and solvent not stated** [!]; samples are taken up in 1.0 mL methanol for colour / HPLC / TEAC, so this is probably an apparent pH of the methanolic (or a diluted) extract, not of a reaction medium |
| Temperature / vessel | **220 °C**, "sealed glass vials" (volume, heating medium not stated); times 0, 2.5, 5.0, 7.5, 10, 15, 20, 30 min; quench in a −20 °C freezer |
| Headspace / O2 / metals | **not stated**: air presumed; no inert gas; no chelator |
| Replicates | triplicate pots, mean ± SD; two-way ANOVA + Tukey |

## 3. Assays [M]

- **Reactants (CA, FA, FF, HMF, PA):** HPLC-DAD (Prodigy ODS-3 C18, 0.05 % acetic acid / MeOH gradient), 270 nm (FF,
  HMF, PA) or 324 nm (CA, FA); **quantified relative to the unheated sample (t = 0)**, dilution 1:90 (0 min) or 1:30.
  No calibration, **no LOD / LOQ**.
- **DHMF, FAHMF, DVG:** external calibration with the isolated products (identity by 1H / 13C NMR and HRMS). No LOD / LOQ.
- **Colour index** = A420 × dilution / mmol of reactants (methanol, quartz cuvette).
- **Replicate noise** (SD, from error bars): HMF alone ±1.5-8 %, FF / PA ±1-12 %, pH ±0.05-0.3 unit.
- [!] **FF reads above 100 %** of its unheated reference (FF alone 125 ± 1 % at 2.5 min, FA/FF 115-128 %): the t = 0
  reference is biased low by ~20-25 % (volatility of FF during the unheated work-up, or a matrix effect), so the
  relative assay carries a systematic error of that size for FF at least.

## 4. Data

### 4.1 Fig. 3B / 3D and Fig. S-4B — heterocycle remaining (% of t = 0), 220 °C, dry [F ±2 %; SD in brackets]

| system | 2.5 min | 5 min | 7.5 min | 10 min | 20 min | 30 min |
|---|---|---|---|---|---|---|
| **HMF alone** (Fig. 3B) | 89.7 (1.5) | 60 (8) | — | 14.5 (6) | 4.2 (4) | 2.2 (2) |
| HMF in FA/HMF | 85 (1.5) | 80 (4) | — | 75 (2) | 65.5 (3) | 59.5 (4) |
| FF in FA/FF | 127.5 (1) | 123 (5) | — | 125.5 (5) | 120.5 (5) | 115.5 (6) |
| PA in FA/PA | 96.5 (6) | 91 (8) | — | 88.5 (4) | 86 (6) | 88 (3.5) |
| HMF in CA/HMF (Fig. 3D) | 64 (9) | 55.5 (10) | 45.5 (5) | 44.5 (8) | — | — |
| FF in CA/FF | 90 (2) | 83.5 (5.5) | 87.5 (4.5) | 85 (6.5) | — | — |
| PA in CA/PA | 91 (6) | 84 (9) | 63 (6) | 36 (3.5) | — | — |
| **FF alone** (Fig. S-4B) | 125.5 (1) | 117 (12) | 111 (10) | 114 (9) | — | — |
| **PA alone** (Fig. S-4B) | 90.5 (6) | 93.5 (3) | 95 (6) | 88 (11) | — | — |

Text agrees: FA/HMF HMF "around 60%" from 20 min; CA/HMF ~45 % after 7.5 min; FA/PA PA ~88 % at 2.5 min; CA/PA PA ~91 %
at 2.5 min and ~36 % at 10 min [M].

### 4.2 Fig. 3A — FA remaining (% of t = 0) [F ±2 %]

FA alone 87 / 77 / 56 / 25 / 12 % at 2.5 / 5 / 10 / 20 / 30 min (text: "around 12%" at 30 min [M]); FA/HMF 75 / 55 / 31
/ 8 / 3; FA/FF 70 / 66 / 61 / 44 / 25; FA/PA 47 / ~1 / traces. CA is gone by 2.5 min in every pot (Fig. S-4D) [M].

### 4.3 Fig. 3C and Fig. S-4C — "pH value" at 0 / 2.5 / 10 min [F ±0.05]

| pot | FA | FA/HMF | FA/FF | FA/PA | CA | CA/HMF | CA/FF | CA/PA | HMF | FF | PA |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 min | 3.83 | 3.65 | 3.25 | 3.75 | 4.45 | 4.50 | 4.32 | 4.30 | 4.7 | 3.75 | 6.1 |
| 2.5 min | 4.33 | 4.10 | 3.78 | 4.80 | 4.88 | 5.05 | 4.53 | 6.40 | 4.85 | 4.6 | 6.8 |
| 10 min | 4.45 | 4.43 | 4.03 | 6.12 | 5.43 | 5.10 | 4.78 | 7.17 | 4.95 | 4.55 | 6.95 |

The authors read the rise as decarboxylation of the phenolic acid (loss of COOH) [Q]. The heterocycles alone come from
Fig. S-4C (0 / 2.5 / 5 / 7.5 / 10 min; only 0 / 2.5 / 10 shown here).

### 4.4 Fig. 5 — products, mol % of FA [M from text]

FA alone: VG max ~10 mol % (5 min), DVG ~18 mol % (10 min). FA/HMF: FAHMF ~5 mol % at 2.5 min, ~7 mol % from 7.5 min;
DHMF ~2 then ~3 mol %; VG ~10 → ~3 mol % (20-30 min); DVG ~21 mol %. FA/FF: VG rising to ~11, DVG ~19 mol % at 30 min.
FA/PA: VG ~63 mol % at 5 min, then ~7 mol %; DVG ~14 mol %.

### 4.5 Colour and TEAC [M from text]

Colour index: FA/PA 117 mmol⁻¹ (30 min, linear), CA/PA 388 mmol⁻¹ (max at 7.5 min), FA/HMF 90, CA/HMF 109 (7.5 min),
FA/FF ~53, HMF alone 29 (7.5 min). TEAC: FA systems ~150 mmol TE/L initially, falling to 83 (FA), 40 (FA/HMF), 56
(FA/PA), 76 (FA/FF); CA systems 95 initially. No engine observable.

### 4.6 Arithmetic [Z]

1. **HMF alone is not first order**: pointwise k = ln(100/C)/t is 0.044 / 0.10 / 0.19 / 0.16 / 0.13 min⁻¹ at 2.5 / 5 /
   10 / 20 / 30 min. The slow start is the heat-up of a dry vial (no come-up stated); no single k describes it.
2. FF and PA alone: no loss within the read and replicate error over 10 min at 220 °C (FF's >100 % is the assay bias of
   sec. 3).

## 5. What it can and cannot test in the model

**Can:** nothing quantitative in the engine's scope. At most a **qualitative** statement that HMF, not FF or PA, is
thermally labile when dry at 220 °C, and that phenolic acids act as HMF / PA sinks in roasting.

**Cannot:**
- **No NF, no water, no buffer, no controlled pH, no O2 statement**: none of E16's blockers is touched.
- Relative assay (t = 0 reference) with a demonstrated ~20-25 % bias for FF; no LOD / LOQ.
- 220 °C neat solids are outside every lane's domain (the engine is a dilute-solution mass-action model).

## 6. Conflicts and agreements with the repo's records

1. **No conflict with any record**: no ledger, fit table or engine rate covers dry HMF / FF / PA at 220 °C.
2. [!] **FF > 100 % of t = 0** (sec. 3): an internal assay defect worth recording before anyone uses the FF rows.
3. [!] **pH method not stated** (sec. 2): the values cannot be read as an aqueous pH.
4. Consistent with `bork2023` Fig. S4 (FF and PA do not decay alone, there at 130 °C in water): neither heterocycle is
   lost on its own in either medium.

## 7. Proposed ledger entries — none

Nothing here is simulable or within a sub-model's scope; no entry is proposed. **Provenance tiers (first reader only):**
the heterocycle and FA fractions of sec. 4.1-4.2 would be `direct_measurement` (digitised from a figure, relative
assay), confidence low (t = 0 bias, dry system); the pH values of sec. 4.3 `direct_measurement` of an undefined
quantity, confidence very low. Second reader pending.

## 8. Supplementary material (read 2026-11-05, first reader)

**Source.** `data/articles/bork2024_supplementary.pdf` (7 pp.). Identified by its first page: the paper's title, the same
four authors and affiliations; figures and table refer to this paper's pots. Text by `pdftotext`; figures from 130 dpi
renders (Fig. S-4 read by eye, ±2 % / ±0.05 pH).

| item | content | tag |
|---|---|---|
| Fig. S-1 / S-2 / S-3 | NMR assignments of DHMF, FAHMF, DVG | [M] |
| **Fig. S-4** | heterocycles heated alone, 220 °C, 0-10 min: (A) colour index FF 3.9 → 0.7, PA 0.4 → 1.4 mmol⁻¹; (B) remaining FF / PA (sec. 4.1); (C) pH HMF / FF / PA (sec. 4.3); (D) CA remaining: 0 in every pot at 2.5-10 min | [F] |
| Fig. S-5 | APCI(+)-HRMS scan spectra of the six binary pots (5 min) | [Q] |
| Fig. S-6 / Table S-1 | MS² of m/z 259, 151, 246, 381, 351, 368 and their fragment assignments | [M] |
| Fig. S-7 | TEAC of the heterocycles heated alone (negligible) | [Q] |
| Fig. S-8 | TEAC vs colour index, FA systems (R² 0.68-0.95) | [M] |

[!] **Caption / panel mismatch in Fig. S-4**: the caption calls panel (D) "antioxidant activity"; the panel plots CA
conversion. **Nothing in the SI on vessel volume, atmosphere, the pH method, water or NF.**

## 9. Open questions it bears on

**E16 (is NF decay hydroxide-only, partial order, or pH-flat?).** Blockers: (1) O2 status, (2) pH drift in unbuffered
pots, (3) no temperature series for NF alone, (4) NF self-reaction with its own fragments at 20 mM.

| item looked for | result |
|---|---|
| NF-alone control at a new pH, temperature or concentration | **NOT FOUND** — no NF pot in the paper or the SI |
| stated headspace or atmosphere | **NOT FOUND** — "sealed glass vials" only |
| pH after heating | **NOT FOUND for any aqueous pot.** A "pH value" after 2.5 and 10 min is printed for the dry pots (sec. 4.3), method and solvent unstated; it is not an NF pot and not an aqueous pH |
| blocker 1 / 2 / 3 / 4 | **none touched** |
| new NF-alone k | **none** (nothing to compare with `bork2022`'s 4.91 ± 0.23 × 10⁻³ min⁻¹) |

**Net:** E16 stays "undecidable on disk"; this paper is neutral for it.

## Scope note (coordinator, 2026-11-05)

"Not simulable / outside every sub-model" above is a judgement, not the engine's rule, and is narrowed here:
- **Phenolic arms (caffeic / ferulic acid + partner):** outside the engine's declared envelope BY RULE: the phenolic acids
  are not species in any core lane (`declare_envelope` refuses an unnamed reactant).
- **Single-compound controls (HMF, furfural, pyrrole-2-carbaldehyde, glyoxal, MGO, diacetyl alone, heated neat or near-dry
  at 220 C):** these species exist in the engine, and no envelope rule refuses them. Only the acrylamide lane has an a_w
  window (0.34-0.99); the trunk's a_w term has no declared window. They are outside the engine's VALIDATED domain (dilute
  aqueous kinetics; no solvent, a_w undefined), not outside its rules. Tagged as **candidate data for plan step W5.2
  (crust, low-moisture high-temperature)**, to be judged there. Register item: give the trunk lane an a_w window, so a
  dry-state question is refused rather than extrapolated.
