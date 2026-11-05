# Bork, Haase, Rohn & Kanzler 2022b (Food Chem. 395:133592) — EXTRACTION (paper key `bork2022b`; the PDF on disk is `data/articles/bork2022.pdf`, see sec. 0)
### *"Structural characterization of polar melanoidins deriving from Maillard reaction intermediates – A model approach"*
### Norfuraneol (NF) 20 mM alone or 1:1 with glyoxal (GO), glycolaldehyde (GA), acetaldehyde (AA) or furfural (FF) in **unbuffered water adjusted to pH 5** (1 M NaOH / HCl), sealed ampoules, **130 °C, 30-300 min**; low-water arm (40 wt% solids); colour index, SEC, HRMS / Kendrick analysis of the oligomers; NF's proposed degradation scheme (Fig. 3). **Every NF time course is in Fig. S-2 (Supplementary, NOT on disk)**; the main text prints only the NF left at 30 min in each binary pot. **No NF-alone number is printed: this paper adds nothing quantitative to `bork2022`'s NF-alone decay.**

**Read 2026-11-05 (first reader) from `data/articles/bork2022.pdf`** (pp. 1-10, Elsevier publisher PDF, clean text layer).
Prose from `pdftotext -raw`. No figure in the main text carries an NF concentration, so nothing was digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Leon V. Bork, Paul T. Haase, Sascha Rohn, **Clemens Kanzler** (corresponding) — TU Berlin, Food Chemistry and Analysis |
| Venue | *Food Chemistry* **395** (2022) 133592 |
| DOI | **10.1016/j.foodchem.2022.133592** (printed p. 1) — matches the brief's second Bork 2022 |
| Dates | received 30 Mar 2022; revised 20 Jun 2022; accepted 26 Jun 2022; online 28 Jun 2022 |

[!] **The file `bork2022.pdf` is THIS paper, not Food Chem. 380:131852** (that one is `bork2022b.pdf`, dossier
`bork2022_extraction.md`). Keys follow the brief and publication order; the crossing is recorded in
`data/keys/pdf_dossier_aliases.yml`.

**Laboratory.** Same group as `bork2022` and `bork2023` (proposed `tu_berlin_kanzler`); one laboratory.

## 1. What this paper answers for the model

| model need | what Bork 2022b gives | tag |
|---|---|---|
| NF decay alone vs time, pH 5, 130 °C | **Not printed** (Fig. S-2, not on disk). The time grid is 30 / 60 / 120 / 180 / 300 min | [!] |
| NF left at 30 min with a partner (1:1, 20 mM each) | **GO 0.7 %**, **MGO 19 %**, **GA 60 %**, **FF 80 %**, **AA 83 %** of the initial NF | [M] |
| reactivity ranking in water (colour index + NF conversion) | GO > MGO > AA > GA > FF | [Q] |
| NF degradation products | proposed: hydration / ring opening to a 2,4-diketone → retro-aldol to **MGO + GA**; oxidative hydrolysis to a **triketone** → retro-aldol to **formaldehyde + acetol**; HRMS of dry-heated NF (125 °C, 10 min) shows NF oligomers with MGO, acetaldehyde and formaldehyde units and oxidised NF | [P] / [Q] |
| pH, O2 | as `bork2022`: unbuffered, pH set to 5 before heating; no degassing, no chelator stated | [!] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | NF 57 mg (0.50 mmol) with GO 40 wt% 72.5 µL, GA 30 mg (0.50 mmol), AA 22 mg (0.50 mmol) or FF 48.0 mg (0.50 mmol) in **25 mL water** → **20 mM each** [Z]; alone or binary |
| pH | "at pH 5" [V]; 1 M NaOH / 1 M HCl; no buffer; **initial, room temperature; not re-measured** |
| Heating | sealed ampoules, 130 °C, 30 / 60 / 120 / 180 / 300 min; freezer quench −20 °C; centrifuged, diluted |
| Low-water arm (sec. 2.6) | same masses at 60 wt% water, 130 °C, 10 / 30 / 60 / 120 min; colour of water and ethyl acetate extracts |
| Replicates | triplicate, mean ± SD |
| NF assay (sec. 2.4) | HPLC-DAD 285 nm, ion-pair (phosphate pH 6 / EDTA / TBAHS) on Nucleosil C18; **external calibration and relative to t = 0** |
| LOD / LOQ | **not given** |

## 3. Data [M]

- Text, sec. 3.1: NF left after **30 min at 130 °C, pH 5**: "0.7 % with regard to the initial concentration" with GO;
  "19 % (methylglyoxal), 60 % (glycolaldehyde), 80 % (furfural), and 83 % (acetaldehyde)" [V, abridged].
- Fig. 1 (colour index after 120 min, aqueous pH 5 and 40 wt% water): bars, not digitised (colour, not a model observable
  here).
- SEC (Fig. 2, 300 min): NF/GO colorants span > 800 kDa to < 0.34 kDa; NF with AA / GA / MGO / DA limited to ~1-10 kDa
  and ~0.34 kDa fractions [Q].

**Arithmetic [Z].** NF alone at 30 min in the same nominal pot is **85.9 %** (`bork2022` Fig. 1D). So the partners add an
apparent first-order NF loss over 0-30 min of ln(0.859/x)/30: **GO 0.16, MGO 0.050, GA 0.012, FF 2.4 × 10⁻³, AA 1.2 ×
10⁻³ min⁻¹** (FF and AA at the read-error limit of the NF-alone value). These are pot-specific (20 mM partner) and not
rate constants of a step.

## 4. What it can and cannot test

**Can:** rank C2 / C3 carbonyls as NF sinks at pH 5 / 130 °C (GO ≫ MGO ≫ GA > FF ≈ AA); supply the authors' product
scheme for NF degradation (MGO + GA by retro-aldol; formaldehyde + acetol via an oxidised triketone) — mechanism [P],
consistent with Mikami 2017's MGO / diacetyl from NF at pH 8 only in the MGO part (diacetyl is not in this scheme).

**Cannot:** give NF-alone kinetics (supplementary only); give fragment yields (none measured); separate oxidative from
hydrolytic NF loss (no O2 control).

## 5. Conflicts and agreements with the repo's records

1. [!] **duplicate_suspected**: "19 % (methylglyoxal)" at 30 min equals `bork2022` Fig. 1C (19.2 %). Probably the same
   incubation reported twice; never enter both.
2. **Engine:** the NF sinks are `r_nf_decay` (→ inert FRAG_C), B61's `r_nf_decay_oh` (+ optional MGO share) and the H2S
   routes; there is no NF + GO / GA / MGO aldol sink. The fast NF + GO loss (t½ ≈ 4 min at 20 mM each [Z]) is a candidate
   missing sink wherever GO and NF co-exist (pentose pots), not tested here.
3. The proposed retro-aldol NF → MGO + GA agrees in direction with B61's `mgo_share` structure (`parameters_sulfur.py`
   sub-model 3); no yield is given to set it.

## 6. Proposed ledger entry (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: bork2022b_text_norfuraneol_left_30min_with_carbonyls_ph5_130c
    dossier: bork2022b_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 395:133592"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol 20 mM + glyoxal, methylglyoxal, glycolaldehyde, furfural or acetaldehyde 20 mM (1:1)"
    medium: "water, no buffer, pH 5 set with 1 M NaOH/HCl before heating; sealed ampoules"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [30]
    species: [norfuraneol]
    kind: endpoint
    location: "text sec. 3.1 p. 3 (Fig. S-2 not on disk)"
    status: text_only
    quant_class: calibrated
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "only traces of norfuraneol (0.7 % with regard to the initial concentration)"
    flags: [ph_start_only, unbuffered, o2_not_stated, single_time_point, si_not_on_disk, duplicate_suspected_bork2022_fig1c]
    note: "NF left at 30 min: GO 0.7, MGO 19, GA 60, FF 80, AA 83 %. NF alone in the same pot: bork2022 Fig. 1D (85.9 % at 30 min)."
```

## 7. Supplementary material

Not on disk. **Fig. S-2 (A-E): NF vs time in each binary pot** — the only place this paper's NF kinetics live; fetch it
if a second NF-alone series (or the GO / GA / AA / FF sinks as time courses) is wanted. Tables S-1 to S-4 (HRMS
assignments), Fig. S-11 (Kendrick plots of NF alone).
