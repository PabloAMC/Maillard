# Yaghmur, Aserin, Abbas & Garti 2005 — EXTRACTION (furfural 3 mmol + L-cysteine 6 mmol in 10 g phosphate 0.5 M pH 5.0, 40-85 C; the same reaction run in water, in a water/propylene-glycol binary, and in five-component food-grade O/W microemulsions; pseudo-first-order furfural disappearance at four temperatures in three media, plus three apparent Arrhenius barriers)

### THE ONLY MEASUREMENT IN THE CORPUS OF HOW FAST FURFURAL DISAPPEARS IN THE PRESENCE OF CYSTEINE, AND THE SOURCE OF THE 1.2 % CEILING THE SULFUR LANE FITS AGAINST: the water arm prints four first-order constants and a 46.50 kJ/mol barrier, but the FFT branch is under a hundredth of that flux, so what is being measured is overwhelmingly a sink the paper never identifies — which is exactly what `r_fur_decay` says it is.

**Source on disk:** `data/articles/yaghmur2005.pdf` (12 pp., Colloids and Surfaces A:
Physicochemical and Engineering Aspects 253 (2005) 223-234, doi 10.1016/j.colsurfa.2004.10.114).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/yaghmur2005.txt`, 654 lines), which came through clean, **Table 1 included**;
journal page 233 was additionally rendered with `pdftoppm` to confirm that Table 1's last column
genuinely has **no r² sub-column** (it does not — see section 3). Figures 1-15 are all images and
every concentration-, conversion- and velocity-versus-time datum in this paper lives in them: **there
is exactly one table, and it is the paper's only tabulated data.** There is no supplementary
material. Repo status before this dossier: the paper is already load-bearing —
`src/kinetic_core/parameters_sulfur.py` carries `YAGHMUR_FFT_SHARE_CEILING = 0.012` with its anchor
(line 1855), a `MEASURED_SULFUR`-adjacent **audit_flag** record of the 46.50 kJ/mol barrier marked
`operative: False` (line ~1457), and the `k_fur_decay` note (line 1976); `src/kinetic_core/sulfur.py`
line 500 names it in `r_fur_decay`'s comment; and `k3_final_parameter_inventory.md` quotes a
`yaghmur2005_extraction.md` at sections B10.15-B10.17, C.14 and F row 29 — **but that dossier is not
on disk.** This file restores it from the primary paper and checks every number the inventory
carries.

## 0. Identity

| field | value |
|---|---|
| Title | "Reactivity of furfural-cysteine model reaction in food-grade five-component nonionic O/W microemulsions" |
| Authors | Anan Yaghmur (present address given as Institute of Chemistry, Physical Chemistry, University of Graz, Austria), Abraham Aserin, Atallah Abbas, **Nissim Garti** (corresponding; tel. +972 2 6586574/5, fax +972 2 6520262, garti@vms.huji.ac.il) — Casali Institute of Applied Chemistry, Givat Ram Campus, The Hebrew University of Jerusalem, 91904 Jerusalem, Israel |
| Venue | Colloids and Surfaces A: Physicochem. Eng. Aspects 253 (2005) 223-234. Received 13 July 2004, accepted 22 October 2004, available online 22 December 2004 |
| DOI | 10.1016/j.colsurfa.2004.10.114 |
| Paper type | **Physical chemistry of a reaction medium.** The chemistry is a vehicle: the subject is microemulsion catalysis. The kinetics is real (section 3.7, Table 1, Figs. 13-15) but the reaction is followed **only through the disappearance of furfural**, and the products are quantified only as figure traces. |
| Compound letters | **A** = 2-furfurylthiol (FFT); **B** = 2-(2-furanyl)-thiazolidine; **C** = 2-(2-furanyl)-thiazoline; **D** = N-(2-mercaptovinyl)-2-(2-furanyl)-thiazolidine |
| The three media compared | (i) **water** (phosphate buffer, pH 5.0); (ii) **"aqueous solution"** = water/propylene glycol at a constant 1/1 weight ratio; (iii) **O/W microemulsion** = R(+)-limonene / ethanol / water / propylene glycol / Tween 60 |
| Lineage | the same group's earlier O/W and W/O studies, refs [19] Yaghmur, Aserin & Garti, J. Agric. Food Chem. 50 (2002) 2878 and [20] Fanun et al., Colloids Surf. A 194 (2001) 175 (the sugar-ester W/O system, source of the 38.8 kJ/mol comparison); the microreactor idea from [18] Vauthey et al. (Nestlé), J. Agric. Food Chem. 48 (2000) 4808 |
| Companions on disk | `meynier2002_extraction.md` / `meynier2004_extraction.md` (ref [34] Meynier & Mottram 1995 is the pH-and-cysteine-degradation source this paper leans on), `whitfield2001_extraction.md` and `whitfield1988_extraction.md` (refs [16], [17]), `hofmann1996_extraction.md` / `hofmann1998b_extraction.md` (refs [11]-[13]), `schieberle1998_extraction.md` / `schieberle2000_extraction.md` |

## 1. Why it matters

The engine's furfural node has three outward edges and this paper bears on the ratio between them:

| row in `sulfur.py` | rate key | order | pH tag | what this paper says |
|---|---|---|---|---|
| `r_fur_fft`: FUR + H2S -> FFT | `k_fur_fft` | 2 | `neutral_h2s` | the branch exists; **<1 % of the furfural charged reaches it in water** |
| `r_fur_fft_hs`: FUR + HS- -> FFT | `k_fur_fft_hs` | 2 | `hs_anion` | same branch, other nucleophile; not separated here |
| `r_fur_decay`: FUR -> FRAG_C 5 | `k_fur_decay` | 1 | — | **this is what the paper actually measures**: ~85 % of the charged furfural leaves by a route the water arm never identifies |

**The ceiling.** Section 3.1 (p. 226) prints, for the **water** reference reaction at 65 °C and 15 h,
that the conversion to FFT "was less than 1 %". Table 1's water arm gives k_obs = 12.6 × 10⁻² h⁻¹ at
65 °C, so over the same 15 h the furfural conversion is 1 − exp(−0.126 × 15) = **84.9 % (mine)**. The
FFT share of the furfural that disappeared is therefore **< 1 / 84.9 = < 1.18 %**, which is the
`YAGHMUR_FFT_SHARE_CEILING = 0.012` the fit already carries and which this dossier now confirms from
the primary text and table rather than from the inventory's summary. It is a **one-sided constraint**
and the corpus has no other on this branch.

**The barrier, and why it must stay non-operative.** Section 3.7 (p. 232) prints three apparent
activation energies: **32.3 ± 1.0** (O/W microemulsion), **46.50 ± 1.0** (water), **56.72 ± 1.0**
kJ/mol (water/PG 1/1). The registry's audit_flag record of 46.50 ± 1.0 and its four objections stand
up against the primary source: it is (i) a lump over ≥ 98.8 % non-FFT flux, (ii) measured at
**40-70 °C** against a module window of 95-145 °C, (iii) an Arrhenius Ea, not the Eyring ΔG‡ the
module's barriers are, and (iv) — a fourth the registry does not state and this dossier adds — it is
a **pseudo**-first-order constant, k_obs = k₂[Cys]ₜ, so its temperature dependence contains
whatever the cysteine pool is doing as well, and cysteine is thermally unstable in water by the
paper's own account (refs [34-37]).

**Something the engine's comment slightly overstates.** `r_fur_decay`'s note says "nearly all the
furfural that disappears goes somewhere the corpus never identifies." That is right **for the water
arm**, where section 3.1 says the reaction "leads to the formation of one main sulfur compound
(compound A)" — FFT, at under 1 %. In the **microemulsion** arm, three further products are named
and are the majority: 2-(2-furanyl)-thiazolidine (B, "the main product"), the thiazoline (C) and the
mercaptovinyl-thiazolidine (D). So the corpus *does* have candidate identities for the furfural sink;
what it does not have is any evidence that they are the sink **in water**, and none of B, C or D is
quantified anywhere except in figures. The engine carries no species for any of the three (its only
thiazolidine is `TTCA`, the Kang 2026 tetrahydroxybutyl compound from a pentose, a different
chemistry). Recorded as Flags 6, not as a change.

**Three further findings, in order of usefulness to the repository.**

1. **The furfural sink is first order in cysteine** — k_obs rises linearly with total cysteine over
   0-12 mmol at r² = 0.98, three replicates, relative error ±1.0 % (Fig. 13). This licences the
   second-order form the engine's furfural steps already use. **But it was measured in the
   microemulsion only**; the paper does not repeat it in water.
2. **The medium moves the rate by 3.7-6.3× and the barrier by 24 kJ/mol** (Table 1, section 3.7).
   An interfacial matrix accelerates furfural loss and lowers its apparent barrier; a 50 % w/w polyol
   *slows* it and raises the barrier. The engine has one furfural sink constant with one barrier and
   no matrix term. This is a bound on how much of the corpus's cross-laboratory furfural scatter is
   matrix rather than chemistry — a quarter of the barrier, on one axis, in one study.
3. **Lower pH accelerates furfural consumption** (Figs. 5, 6; pH 3.5, 4.0, 6.0, 8.0 at 65 °C), which
   the authors attribute to faster H2S release from cysteine at low pH. The engine's furfural-to-FFT
   pair already splits `neutral_h2s` from `hs_anion` for this reason. **Every number in that pH study
   is figure-only** and none of it can be typed.

What this paper does NOT give: any concentration of any product in mass or molar units; any absolute
yield of B, C or D; any H2S measurement; any temperature above 85 °C; any measurement of MFT; any
rate for the FFT branch itself; any identification of the ~85 % of furfural that vanishes in water;
and any replicate count outside Fig. 13.

## 2. Methods as they matter to a model

- **The reaction, as printed (section 2.3).** "Reactions were carried out at 65 °C in Tween 60-based
  O/W microemulsions and were performed as follows: **cysteine (6 mmol)** was dissolved in the
  **phosphate buffer (10 g, 0.5 M, pH 5.0)** while **furfural (3 mmol)** was dissolved in the oil
  phase (R(+)-limonene plus ethanol). The nonionic surfactant was dissolved in propylene glycol (PG)
  at 37 °C, and added to the aqueous phase. The oil phase was added to the aqueous phase, transferred
  to a 100 mL vial, and mixed to obtain a single-phase O/W microemulsion with the desired
  compositions. The microemulsion was **heated at 65 °C in a water bath while stirring with a
  magnetic stirrer**. The reaction was **quenched by rapid cooling on ice**."
- **The water reference arm (the arm the repository uses).** "For the reference reactions (reactions
  carried out in aqueous phase), the same heating procedure was applied to **10 g of the buffer
  solution containing the reactants**. After heating, the mixture was cooled and **other components
  of the microemulsion were added** in order to have similar workup conditions for the isolation of
  volatiles." So the water arm is **6 mmol cysteine + 3 mmol furfural in 10 g of 0.5 M phosphate at
  pH 5.0** — a **2 : 1 cysteine : furfural** charge — and the microemulsion components are added only
  *after* the reaction, purely to equalise the extraction.
- **Concentrations (mine, and approximate).** Treating 10 g of 0.5 M phosphate as ~10 mL:
  **cysteine ~600 mmol/L, furfural ~300 mmol/L**. The buffer's density is above 1.0, so these are
  slight over-estimates of volume and therefore slight under-estimates of concentration. The paper
  never prints a molarity for either reactant (Flags 3).
- **pH.** 5.0, set by 0.5 M phosphate, in every experiment except the pH study of section 3.2, which
  used **3.5, 4.0, 6.0 and 8.0** (Figs. 5, 6). pH is **never re-measured** after heating and no drift
  is reported. In the microemulsion the "pH" is the pH of the aqueous phase (water plus PG) as made
  up, which is not the pH at the interface where the authors argue the reaction happens.
- **Temperatures.** 65 °C throughout the product work; **85 °C** for section 3.1's temperature test;
  **40, 50, 65, 70 °C** for the kinetics of Table 1 and Fig. 15. **The maximum temperature in this
  paper is 85 °C, and the kinetic series tops out at 70 °C.**
- **Work-up.** "After cooling, the reaction mixture was extracted with diethyl ether. The extract
  (organic phase) was dehydrated by anhydrous MgSO4 at 4 °C. **Naphthalene was used as an internal
  standard.** The same extraction method was applied for the isolation of the volatiles from
  reference reaction." One internal standard, added at extraction; no recovery correction is
  described and no response factors are printed.
- **Quantification — GC-FID.** Hewlett-Packard 5880 with FID; **Rtx-1701** fused silica capillary
  30 m × 0.32 mm × 0.25 µm (Thames Restek); nitrogen carrier, column head pressure 107 kPa; injector
  170 °C, detector 250 °C; oven 50 °C (2 min) then 6 °C/min to 240 °C (10 min). **This is the
  quantitative instrument** ("was used for quantification").
- **Identification — HRGC/MS.** HP 5890 GC-MS, HP5 column 30 m × 0.32 mm × 0.25 µm, helium at
  103 kPa (15 psi); EI at 70 eV, scan **40-250 amu**; 1 µL of concentrated sample; oven 50 °C then
  4 °C/min to 240 °C (10 min); injector 190 °C. **"The compounds were tentatively identified by
  comparing their mass spectra with those contained in the mass spectrometer data system library and
  in previously published literature [18]."** — **tentative library identification, no authentic
  standards, no retention-index confirmation** (Flags 5).
- **PGSE NMR (self-diffusion).** Varian Inova 400 MHz, 5 mm indirect-detection PFG probe, 20 A
  Highland gradient amplifier, **25 ± 0.5 °C**. Equation (1) as printed:
  `E/Eo = exp[−γ² g² δ² (Δ − δ/3) D]`. This is structural characterisation of the medium, not
  chemistry, and all its results (Figs. 11, 12) are figure-only.
- **Phase diagrams.** Pseudo-ternary, constructed at 25 °C; "the accuracy in the location of the
  phase boundaries is within 4 wt.%". Dilution lines **T64** (surfactant/oil 6/2), **T73** (7/1.5),
  **T82** (8/1); microemulsion droplet sizes quoted from ref [25] as **10-12 nm**.
- **The kinetic model, as printed (section 3.7).** Equations (2)-(6):
  `rate = k2 [Fur_t][Cys_t]` (2); with cysteine in excess, `rate = k_obs [Fur_t]` where
  `k_obs = k2 [Cys_t]` (3); `rate = d[Fur_t]/dt = −k_obs [Fur_t]` (4);
  `ln([Fur_t]/[Fur_t]0) = ln(1 − X_Fur) = −k_obs t` (5); `ln k_obs = ln A − Ea/RT` (6), "where the
  overall rate constant (k_obs) is in h⁻¹". **k_obs values are the slopes of −ln(1 − X_Fur) against
  t** (Fig. 14). Note that equation (4) as printed omits the minus sign on the left-hand side that
  the derivation requires; the sense is unambiguous from equation (5).
- **Replication.** Stated **only** for the cysteine-dependence runs of Fig. 13: "Three replicates were
  prepared for these runs and a relative experimental error of ±1.0 % was found." No replicate count
  is given for Table 1, for the conversion figures, or for anything else. The ±1.0 on the three
  activation energies is not explained and no n is attached to it.

## 3. Tables re-typed

There is exactly **one** table. It is re-typed in full below, with the column structure verified
against a 100 dpi render of journal page 233 — **the "Aqueous solution" column has no r² of its
own**, which is how the paper prints it, not a transcription loss.

### Table 1. "Pseudo first-order rate constants, k_obs, and linear regression (r²) for water and O/W microemulsions obtained from data in Fig. 14"

| Temperature (°C) | Microemulsion (O/W) k_obs (h⁻¹) × 10² | r² | Water k_obs (h⁻¹) × 10² | r² | Aqueous solution ᵃ k_obs (h⁻¹) × 10² |
|---|---|---|---|---|---|
| 40 | 20.1 | 0.99 | 3.2 | 0.97 | 1.7 |
| 50 | 29.1 | 0.97 | 6.7 | 0.96 | 4.7 |
| 65 | 48.6 | 0.98 | 12.6 | 0.94 | 9.3 |
| 70 | 59.8 | 0.97 | 16.0 | 0.95 | 12.2 |

Footnote ᵃ as printed: "The aqueous solution composed of water/PG at a constant weight ratio of 1/1."
Caption cross-reference: the data come from **Fig. 14**, whose caption states the microemulsion arm
"contains 60 wt.% aqueous phase (buffered water/PG with a weight ratio of 2/1)" and that the water
arm is "water (pH 5.0)", with the microemulsion reaction "carried out at different temperatures
varying from 40 to 70 °C".

### Numbers printed in the running text (everything else in this paper is figure-only)

| quantity | value | where |
|---|---|---|
| **conversion to FFT, O/W microemulsion, 15 h, 65 °C**, as aqueous phase rises 65 -> 90 wt.% | **4.9-7.3 %** | sec. 3.1 p. 226 |
| **conversion to FFT, WATER, 15 h, 65 °C, same conditions** | **"less than 1 %"** | sec. 3.1 p. 226 — **the ceiling's numerator** |
| in water the reaction "leads to the formation of **one main sulfur compound** (compound A)" | — | sec. 3.1 p. 226 |
| in the O/W microemulsion "the main product ... is 2-(2-furanyl)-thiazolidine (compound B)" | — | sec. 3.1 |
| compound B's time course | "fast (intermediate formation) and the very slow (intermediate disappearance)"; conversion rises to a plateau "followed by a very slight decrease" | sec. 3.1 |
| compounds A, C, D formation rate | "increased slightly with time (Fig. 2a, c and d)" | sec. 3.1 |
| effect of replacing ethanol by butanol on the initial velocity | **"an increase of 28.5 % in the value of V₀"** | sec. 3.3 |
| alcohol ordering | highest V₀ with **hexanol**; partition into the aqueous phase in the order hexanol < butanol < propanol < ethanol | sec. 3.3 |
| effect of the surfactant/oil ratio at 90 wt.% aqueous phase, 3 h | changing T64 -> T73 and T64 -> T82 **decreases** the conversion of furfural by **7 %** and **37 %** | sec. 3.4 |
| the three structural regions along the T64 dilution line | Region I 10-30 wt.% aqueous (V₀ unchanged, W/O); Region II 30-60 wt.% (V₀ rises slightly, bicontinuous); Region III 60-90 wt.% (**sharp linear increase** in V₀, O/W) | sec. 3.5 |
| **Ea, O/W microemulsion** | **32.3 ± 1.0 kJ/mol** | sec. 3.7 p. 232 |
| **Ea, water** | **46.50 ± 1.0 kJ/mol** | sec. 3.7 p. 232 |
| **Ea, aqueous solution (water/PG 1/1)** | **56.72 ± 1.0 kJ/mol** | sec. 3.7 p. 232 |
| Ea, W/O microemulsion on sugar esters (**cited, ref [20], not measured here**) | 38.8 kJ/mol | sec. 3.7 |
| k_obs is **linear in total cysteine over 0-12 mmol** | r² = 0.98; **three replicates**; relative experimental error **±1.0 %**; furfural held at 3 mmol; 90 wt.% aqueous phase; 65 °C | sec. 3.7, Fig. 13 |
| microemulsion droplet size (**cited, ref [25]**) | 10-12 nm | Conclusions |
| phase-boundary accuracy | within 4 wt.% | sec. 2.2 |

**Everything else is figure-only.** Fig. 1 phase diagram; **Figs. 2a-d** conversion (%) to A, B, C, D
against time at 65 °C for 60/65/75/90 wt.% aqueous phase; **Fig. 3** the proposed pathways in
microemulsion (a) and water (b); **Figs. 4a-d** the same four conversions at **85 °C** for 75 and
90 wt.%; **Fig. 5** furfural conversion against time at pH 3.5/4.0/6.0/8.0; **Fig. 6** V₀ against pH
for four media; **Figs. 7a,b** alcohol chain length; **Figs. 8, 9** surfactant/oil ratio; **Fig. 10**
V₀ against aqueous-phase content 10-90 wt.%; **Figs. 11, 12** self-diffusion coefficients;
**Fig. 13** k_obs against cysteine; **Fig. 14** −ln(1 − X_Fur) against t; **Fig. 15** the semi-log
Arrhenius plot. Per house rule none of these is typed as a number.

### Arithmetic on the printed constants (all mine)

**1. Furfural conversion in water at 65 °C over 15 h — the ceiling's denominator.**
From Table 1, water k_obs(65 °C) = 12.6 × 10⁻² h⁻¹ = 0.126 h⁻¹. Equation (5) gives
X_Fur = 1 − exp(−0.126 × 15) = 1 − exp(−1.89) = 1 − 0.1511 = **0.849 = 84.9 %**. Combined with the
printed "less than 1 %" FFT conversion of the furfural charged, the FFT share of the furfural **that
actually reacted** is < 1/84.9 = **< 1.18 %**, i.e. the registry's **≤ 1.2 %** ceiling.
**This is a composite of one printed number and one number I computed from a different figure's
data**, and its two halves come from two different experiment sets in the paper (section 3.1's 15 h
product run and Fig. 14's four-temperature kinetic run). Both are nominally water at pH 5.0 with the
standard charge, but the paper never states the cysteine charge of the Fig. 14 water run
(Flags 2).

**2. The same in the microemulsion, for contrast.** Microemulsion k_obs(65 °C) = 0.486 h⁻¹ over 15 h
gives X_Fur = 1 − exp(−7.29) = **99.93 %**, against a printed FFT conversion of 4.9-7.3 %, so the FFT
share of the furfural flux there is **4.9-7.3 %** — **four to six times the water share**. Caveat:
Table 1's microemulsion is 60 wt.% aqueous phase while the 4.9-7.3 % range is quoted for 65-90 wt.%,
so this is a cross-panel comparison and should be treated as indicative.

**3. Two-point Arrhenius check on all three media (mine).** Using the 40 °C (313.15 K) and 70 °C
(343.15 K) columns, Ea = R ln(k₇₀/k₄₀) / (1/313.15 − 1/343.15) with
1/313.15 − 1/343.15 = 2.7918 × 10⁻⁴ K⁻¹:

| medium | k₄₀ | k₇₀ | Ea, two-point (mine) | Ea, printed |
|---|---|---|---|---|
| O/W microemulsion | 20.1 | 59.8 | **32.5 kJ/mol** | 32.3 ± 1.0 |
| water | 3.2 | 16.0 | **47.9 kJ/mol** | 46.50 ± 1.0 |
| aqueous solution (water/PG 1/1) | 1.7 | 12.2 | **58.7 kJ/mol** | 56.72 ± 1.0 |

All three reproduce to within 2 kJ/mol, i.e. within the scatter a four-point regression on data with
r² = 0.94-0.99 would give. **The table and the three barriers are mutually consistent** — unlike some
papers in this corpus, the temperature columns here are genuine independent measurements, not
back-calculations from a fitted (A, Ea) pair, and the two-point values differing from the printed
ones by 1-2 kJ/mol is the signature of a real four-point fit.

**4. Medium ratios at each temperature (mine).** Microemulsion / water: **6.28** (40 °C), **4.34**
(50 °C), **3.86** (65 °C), **3.74** (70 °C) — this is the inventory's "3.7-6.3×", confirmed, and note
that it **shrinks as temperature rises**, which is the direct consequence of the microemulsion's
lower barrier. Water / aqueous-solution: 1.88, 1.43, 1.35, 1.31 — the polyol slows the reaction and
the penalty also shrinks with temperature.

**5. Furfural half-lives from the water arm (mine, first order).** t₁/₂ = ln2/k_obs: **21.7 h**
(40 °C), 10.3 h (50 °C), **5.50 h** (65 °C), 4.33 h (70 °C). For scale, extrapolating the printed
46.50 kJ/mol to 145 °C would give a half-life of about 12 min — **and that extrapolation is exactly
what the registry forbids** (Flags 1); it is quoted here only to show how far outside the measured
window the module operates.

**6. An apparent second-order constant, and why it is not one (mine, and conditional).** If the
cysteine linearity of Fig. 13 held in water, then k₂ = k_obs/[Cys]ₜ = 0.126 h⁻¹ / 0.600 mol/L =
**0.21 L mol⁻¹ h⁻¹ = 3.5 × 10⁻³ M⁻¹ min⁻¹ at 65 °C**. Three separate assumptions stand between the
printed table and that number: (i) the linearity was demonstrated **only in the microemulsion**;
(ii) the cysteine concentration is my ~600 mmol/L from a mass-based recipe; (iii) the Fig. 14 water
run's cysteine charge is not restated in the paper. **Do not enter this number anywhere.** It is
recorded so that nobody re-derives it and thinks it printed.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Keyed: `furfural`, `2_furfurylthiol`,
`2_methyl_3_furanthiol`, `hydrogen_sulfide`. **Not keyed:** cysteine, 2-(2-furanyl)-thiazolidine,
2-(2-furanyl)-thiazoline, N-(2-mercaptovinyl)-2-(2-furanyl)-thiazolidine, propylene glycol,
R(+)-limonene, Tween 60, naphthalene.

Every row below shares: **furfural 3 mmol + L-cysteine 6 mmol** (2 : 1) in **10 g of 0.5 M phosphate
at pH 5.0**, stirred in a water bath, quenched on ice, diethyl-ether extraction with **naphthalene**
internal standard, quantified by GC-FID, products identified only tentatively by MS library.
Temperatures **40-85 °C**; the module's own window is 95-145 °C.

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| **furfural disappearance in the presence of cysteine, WATER** | **k_obs at 40 / 50 / 65 / 70 °C** | **3.2 / 6.7 / 12.6 / 16.0** (r² 0.97 / 0.96 / 0.94 / 0.95) | **10⁻² h⁻¹** | water, pH 5.0, phosphate 0.5 M | **pseudo-first order in furfural**, with k_obs = k₂[Cys]ₜ | Table 1 p. 233 | **measured_rate** — but a **lump over ≥ 98.8 % non-FFT flux**; never file under a step named "furfural -> FFT" |
| " | **Ea (water)** | **46.50 ± 1.0** | kJ/mol | 40-70 °C, water pH 5.0 | — | sec. 3.7 p. 232 | **measured_barrier** — already carried in `parameters_sulfur.py` as an **audit_flag, `operative: False`**, and this dossier confirms every digit and adds a fourth objection (it is a pseudo-order constant) |
| furfural disappearance, **O/W microemulsion** | k_obs at 40 / 50 / 65 / 70 °C | 20.1 / 29.1 / 48.6 / 59.8 (r² 0.99 / 0.97 / 0.98 / 0.97) | 10⁻² h⁻¹ | Tween 60 / R(+)-limonene / EtOH 3/1/1, 60 wt.% aqueous (water/PG 2/1), pH 5.0 | pseudo-first order in furfural | Table 1 | **measured_rate** — a **different matrix**, not transportable to an aqueous pot |
| " | Ea (microemulsion) | 32.3 ± 1.0 | kJ/mol | 40-70 °C | — | sec. 3.7 | **measured_barrier** (matrix-specific) |
| furfural disappearance, **water/PG 1/1** | k_obs at 40 / 50 / 65 / 70 °C | 1.7 / 4.7 / 9.3 / 12.2 (**no r² printed**) | 10⁻² h⁻¹ | 50 % w/w propylene glycol, pH 5.0 | pseudo-first order in furfural | Table 1 | **measured_rate** (matrix-specific) |
| " | Ea (water/PG) | 56.72 ± 1.0 | kJ/mol | 40-70 °C | — | sec. 3.7 | **measured_barrier** (matrix-specific) |
| **FUR -> FFT branch share, WATER** | conversion of charged furfural to FFT at 65 °C / 15 h | **< 1** | % of furfural charged | water, pH 5.0, 6 mmol Cys : 3 mmol Fur | — | sec. 3.1 p. 226 | **fed_intermediate_yield** (an upper bound, printed as an inequality) |
| **the ceiling the fit uses** | FFT share of the furfural **flux** | **≤ 1.2** | % | as above, over 84.9 % furfural conversion | — | printed "< 1 %" ÷ my 84.9 % (mine) | **derived_assumption** built on a `fed_intermediate_yield` and a `measured_rate` — this is `YAGHMUR_FFT_SHARE_CEILING = 0.012`, and it is the **one FIT row** this paper supplies |
| FUR -> FFT branch share, microemulsion | conversion to FFT at 65 °C / 15 h, aqueous phase 65 -> 90 wt.% | 4.9-7.3 | % of furfural charged | O/W, pH 5.0 | — | sec. 3.1 | **fed_intermediate_yield** (matrix-specific; **not** a licence to raise the water ceiling) |
| order in cysteine | k_obs linear in [Cys]ₜ over 0-12 mmol | r² = 0.98, 3 replicates, ±1.0 % relative error | — | **microemulsion only**, 90 wt.% aqueous, 65 °C, furfural fixed at 3 mmol | **first order in cysteine** ⇒ second order overall | sec. 3.7, Fig. 13 | **within_study_ratio** (a structural/order result, not a constant) |
| matrix acceleration | k_obs(microemulsion)/k_obs(water) | 6.28 / 4.34 / 3.86 / 3.74 at 40 / 50 / 65 / 70 °C | — | as above | — | derived from Table 1 (mine) | **within_study_ratio** |
| polyol retardation | k_obs(water)/k_obs(water-PG) | 1.88 / 1.43 / 1.35 / 1.31 at 40 / 50 / 65 / 70 °C | — | as above | — | derived from Table 1 (mine) | **within_study_ratio** |
| furfural half-life, water | t₁/₂ | 21.7 / 10.3 / 5.50 / 4.33 | h at 40 / 50 / 65 / 70 °C | water, pH 5.0 | — | derived from Table 1 (mine) | **derived_assumption** (arithmetic only) |
| butanol vs ethanol on V₀ | +28.5 | % | 90 wt.% aqueous O/W, 65 °C | — | sec. 3.3 | **level_only** (matrix effect, no engine analogue) |
| surfactant/oil ratio T64 -> T73 / T82 on furfural conversion at 3 h | −7 / −37 | % | 90 wt.% aqueous, 65 °C | — | sec. 3.4 | **level_only** |
| Ea in a sugar-ester W/O microemulsion | 38.8 | kJ/mol | — | — | sec. 3.7, **citing ref [20]** | **level_only** — **not measured here** |
| products B, C, D (thiazolidine, thiazoline, mercaptovinyl-thiazolidine) exist and B is the main microemulsion product | — (presence) | — | O/W, 65 and 85 °C | — | sec. 3.1, Figs. 2b, 4b | **threshold** — and only **tentatively** identified (Flags 5) |
| all conversion-, velocity-, diffusion- and Arrhenius-vs-time/pH/composition traces | — | — | — | — | Figs. 1-15 | **figure_only** |

### Can these be put on the same basis as the sulfur lane's constants? Step by step.

**(a) The ceiling — yes, and it is already in.** A branch *share* is dimensionless and matrix-robust
in a way a rate is not, and it is the one quantity here that survives the transport objections. It
constrains `k_fur_fft · [H2S] / k_fur_decay` at 65 °C and pH 5.0. **One caveat the registry does not
state:** no H2S was charged — the sulfide comes from cysteine degradation in situ and is never
measured — so the ceiling constrains a **flux ratio at an unknown sulfide level**, not a rate ratio.
If the engine's H2S pool at 65 °C and pH 5.0 is wrong, the ceiling will be met for the wrong reason.

**(b) The water rate — record it, do not score against it.** Four temperatures, a real four-point
Arrhenius fit that my two-point check reproduces, an internal standard, and a plain pseudo-first-order
form: as a *measurement* this is sound. As a *constant for this engine* it fails on four counts, of
which the registry names three. It is a lump over ≥ 98.8 % of the flux; the module has no species for
the ~85 % that leaves; 40-70 °C is 25-75 °C below the module's floor; and k_obs is not an elementary
constant. Its most defensible use is a **residual check**: an engine run at 65 °C, pH 5.0, with this
charge should lose furfural at roughly 0.126 h⁻¹ and should put under 1.2 % of it into FFT. That is a
prediction the module can be held to without ingesting a number.

**(c) The matrix result — no transport, but a warning worth keeping.** 24 kJ/mol of barrier and a
factor 3.7-6.3 of rate move with the medium, at one pH, one charge, one reaction. Whenever two
laboratories' furfural numbers disagree by that much, this paper is the reason not to conclude the
chemistry differs.

**(d) What cannot be transported at all.** The alcohol-chain-length series, the surfactant/oil
series, the aqueous-phase-content series, the self-diffusion coefficients and the whole pH study —
the first four because the engine has no matrix axis and no species for any microemulsion component,
the last because every number in it is inside a figure.

## 5. Flags

1. **The 46.50 kJ/mol barrier must stay `operative: False`, and there is a fourth reason.** The
   registry's three objections (a lump over ≥ 98.8 % non-FFT flux; 40-70 °C against a 95-145 °C
   window; an Arrhenius Ea is not an Eyring ΔG‡) are confirmed against the primary text. **Add a
   fourth: k_obs is a pseudo-first-order constant, k_obs = k₂[Cys]ₜ**, so its Arrhenius slope carries
   the temperature dependence of the cysteine pool as well as of the reaction — and the paper itself
   says cysteine is thermally unstable in water and degrades to cysteamine, H2S and
   mercaptoacetaldehyde (refs [34-37]). The measured barrier is the barrier of a composite.
2. **The 1.2 % ceiling is a composite of a printed inequality and my arithmetic, and its two halves
   come from different runs.** "< 1 %" is printed for the water reference at 65 °C / 15 h in
   section 3.1. The 84.9 % denominator is **mine**, computed from Table 1's water k_obs, whose run
   (Fig. 14) is described only as "water (pH 5.0)" — **the paper never restates its cysteine
   charge**. If the Fig. 14 water run used a different cysteine loading than section 3.1's reference
   reaction, k_obs would differ and so would the denominator. The ceiling is sound as an
   order-of-magnitude bound and should keep the "≤" that it has; it should not be tightened.
3. **All molar concentrations are mine.** The recipe gives **masses and millimoles into 10 g of
   buffer**, never a molarity: my ~600 mmol/L cysteine and ~300 mmol/L furfural treat 10 g as 10 mL
   and are therefore slight under-estimates. The **2 : 1 cysteine : furfural ratio is printed** and is
   the robust part.
4. **Nothing here is at a cooking temperature.** 40-70 °C for every rate, 85 °C for the highest
   product run. The module's fit rows sit at 95-145 °C. Any Arrhenius extrapolation from this paper
   into the module's window crosses 25-75 °C of unmeasured ground on a lumped barrier.
5. **Every product identification is tentative.** Section 2.6: "The compounds were tentatively
   identified by comparing their mass spectra with those contained in the mass spectrometer data
   system library and in previously published literature [18]." No authentic standards, no retention
   indices, no NMR. The **quantification** is separately by GC-FID against naphthalene with no
   response factors printed, so the conversion percentages in Figs. 2, 4, 5, 7-10 rest on an
   uncorrected FID response. **The "< 1 %" and "4.9-7.3 %" FFT conversions inherit that.**
6. **The furfural sink is partly named here, but only in the microemulsion.** Compounds B, C and D
   (2-(2-furanyl)-thiazolidine, its thiazoline, and the mercaptovinyl-thiazolidine) are the majority
   products **in the O/W medium**; the water arm is said to give "one main sulfur compound", FFT,
   at under 1 %. So `r_fur_decay`'s comment ("goes somewhere the corpus never identifies") remains
   correct for water, and the three thiazolidine-family compounds are **candidates**, not the
   answer. The engine has no species for any of them and no rate exists for any of them. Do not
   read Figs. 2b/4b as evidence about the water sink.
7. **The cysteine order was measured in the microemulsion only.** Fig. 13's linearity (r² = 0.98,
   three replicates, ±1.0 %) is the paper's only replicated result and its only order determination,
   and it is in the 90 wt.% aqueous O/W system. The water arm's order in cysteine is **assumed**, not
   measured. The engine's second-order furfural steps are consistent with it; they are not confirmed
   by it.
8. **The pH study is entirely figure-only.** Figs. 5 and 6 give furfural conversion against time at
   pH 3.5, 4.0, 6.0 and 8.0, and V₀ against pH for four media, at 65 °C. The qualitative result —
   lower pH, faster furfural consumption, attributed to faster H2S release — supports the engine's
   `neutral_h2s` / `hs_anion` split in direction. **No number from it may be typed**, and it would in
   any case be a formation-side pH response measured through a disappearance, not a pH dependence of
   any single step.
9. **The "aqueous solution" arm has no r² column at all** (verified on the page render), and the
   ±1.0 attached to all three activation energies has no stated n, no confidence level and no
   derivation. Three barriers all carrying exactly ±1.0 is a suspiciously uniform uncertainty.
10. **No H2S is charged or measured anywhere.** The sulfide that makes FFT is generated in situ from
    cysteine, and its level at 65 °C and pH 5.0 is never quantified. The ceiling is therefore a
    constraint on a flux ratio at an unknown sulfide concentration (see section 4a).
11. **The dossier the inventory cites was missing from disk.** `k3_final_parameter_inventory.md`
    quotes `yaghmur2005_extraction.md` at sections B10.15, B10.16, B10.17, C.14 and F row 29, but no
    such file existed in `data/lit/extraction_dossiers/` before this one. **All five of the
    inventory's claims check out against the primary paper**: the ≲1.2 % share, the first-order
    cysteine dependence with its microemulsion-only caveat, the 3.7-6.3× acceleration, the
    46.50 -> 32.3 -> 56.72 barrier ordering, and the four k_obs values per medium.
12. **What this paper does not contain**: any absolute concentration of any product; any yield of B,
    C or D; any H2S measurement; any MFT; any rate for the FFT-forming step itself; any temperature
    above 85 °C; any identification of the ~85 % of furfural that disappears in water; any replicate
    count outside Fig. 13; any error bar on any k_obs; any supplementary material.
13. **What to request from the authors**: (i) the cysteine and furfural charges for the Fig. 14 water
    runs specifically, which is what would firm up the ceiling's denominator; (ii) the numeric data
    behind Figs. 2 and 4 — the four conversion-vs-time series in both media at 65 and 85 °C — which
    would turn the paper's only product kinetics from figure-only into rows; (iii) the FID response
    factors used for FFT and for compounds B, C, D against naphthalene; (iv) whether the water arm
    was ever assayed for compounds B, C and D at all, and at what detection limit — this is the
    single question that decides whether the water sink is the thiazolidine family or something else;
    (v) a mass balance on furfural in water, i.e. what the other ~85 % is; (vi) the n behind the
    ±1.0 kJ/mol on the three barriers.
14. **Registry gaps against `data/keys/compounds.yml`**: `furfural`, `2_furfurylthiol` and
    `hydrogen_sulfide` are present. **Absent: cysteine** (a reactant in every row of this paper),
    **2-(2-furanyl)-thiazolidine**, **2-(2-furanyl)-thiazoline** and
    **N-(2-mercaptovinyl)-2-(2-furanyl)-thiazolidine** — the three named candidates for the furfural
    sink, none of which the engine carries as a species either.
