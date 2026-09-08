# Where the meaty thiols go: the candidate sinks, from the dossiers on disk

*Hand-curated on 2026-09-08 from the extraction dossiers named in each row, for the sink-structure
re-calibration (plan section 7, W7). It answers three questions for every way MFT and FFT can leave a
pot: what the literature measured, what the model has, and what one measurement would decide. The
hypothesis layer (`results/validation/network_hypotheses.md`) lists the same steps from the reaction
rules; this table adds the numbers. Nothing here is a rate the model may use without a
pre-registered wave.*

## 1. What the model has today

Read from `src/kinetic_core/sulfur.py` and the shipped fit report. At trace thiol levels the model's
removal is one first-order decay per thiol with one shared barrier on the ceiling of its band; the
only reversible sink is inert because its partner pool is never charged.

| step in the model | order in thiol | rate at 145 °C | barrier | reversible | partner pool | source |
|---|---|---|---|---|---|---|
| first-order decay (`ch_mft_decay`, `ch_fft_decay`) | 1 | MFT 1.36 /min, FFT 0.22 /min | 102 kJ/mol, on the ceiling of the Gigl band (7 to 102) | no | none | fitted; barrier band from Gigl 2021 |
| thiolate loss (`ch_thiolate_loss_*`), pH-shaped | 1 (thiolate fraction) | 0.015 /min at pH 5 | 102 kJ/mol (shared) | no | none | fitted on Kumazawa 2003's pH grid |
| oxidative dimerisation (`ch_dimer_*`), 2 thiol + oxidant | 2 | 3.2 L²/(mmol²·min), on its upper bound | 122 kJ/mol from Zhang 2026's peptide dimer | no; dimers decay at 2e-10 /min (never) | oxidant equivalents, 1 unit, depletable | fitted on Zhang 2024's dimer shares |
| methanethiol coupling (`ch_mmft`) | 1 (× MeSH) | 0.010 L/(mmol·min) | none (held) | no | methanethiol, from thiamine only | fitted on Zhang 2024 |
| thioether with a matrix electrophile (`ch_thioether_*`) | 1 (× pool) | 1.6e-3 L/(mmol·min) | 10.8 kJ/mol | **yes** (release from Stack 2018's K) | MELE, **never charged** in any fit or deployment | Hofmann 2002, Charles-Bernard 2005 |
| protein disulfide exchange (`ch_protein_ss_*`) | 1 (× sites) | 6e-6 L/(mmol·min) | none | no | PROT_SS, zero in every panel system | Anantharamkrishnan 2020b's 6 to 24 h bracket |
| oligomer (`ch_oligomer_*`) | 1 | 0 | none | no | none | van Seeventer 2001, declared hold-out |

What the model lacks as a step: thiol plus aldehyde or dicarbonyl (hemithioacetal, thiazolidine),
thiol addition to an unsaturated aldehyde, disulfide reduction back to thiol, mixed MFT-FFT
disulfides, FFT plus methanethiol, a reversible non-covalent reservoir, and any sink that scales
with a pool the pot itself makes.

## 2. What the literature measured, sink by sink

| sink | thiol, partner, product | conditions | what was measured | reversible? saturates? | in the model | dossier |
|---|---|---|---|---|---|---|
| **oxidation to the symmetric disulfide** | FFT + O2 → difurfuryl disulfide (major product), plus furfural, furfuryl alcohol, unaccounted non-volatiles | citrate-phosphate buffer, 121 °C 10 min, air in the can | residual FFT 99.5 / 96 / 89 / 80 / 45 / 11 / 0.1 % at pH 3 / 4 / 5 / 5.4 / 6 / 6.4 / 7; disulfide/FFT area ratio rises from <0.02 to >17 across the same ladder; the volatile products do not close the mass balance | not tested; the apparent rate at pH 6 halves when the cook doubles (0.073 vs 0.042 /min), read by the dossier as saturation of the oxidant in the can | yes, 2nd order with an oxidant pool; the pH slope is fitted | kumazawa2003 |
| the same, in a Maillard pot | MFT, FFT, mercaptoketones → symmetric and mixed disulfides | cysteine + ribose, 140 °C 30 min, pH 4.2 to 5.6, sealed air headspace | disulfides are 3 to 10 % of their parent thiols everywhere (MFT dimer / MFT 0.016 to 0.04); at pH 4.2 the dimer tracks MFT | single time point | yes | mottram2002 |
| the same, at 120 °C | MFT, FFT → dimers; FFT → di-2-furfuryl sulfide | Cys-Amadori + cysteine, water, 120 °C 60 min, pH 6 to 8 | dimer holds 6.5 to 9.6 % of MFT and 0 to 4.4 % of FFT as thiol equivalents, nearly flat in pH while MFT swings threefold | single time point | yes | zhou2023 |
| the same, with an oxidised additive | MFT → MFT dimer | thiamine + xylose, pH 4.9, 115 °C 60 min, cystine vs cysteine vs glutathione | cystine sends 54 % of MFT to the dimer, cysteine 9 %, glutathione 7 %: the dimer share is set by what oxidises the thiol | single time point; six of seven dose panels turn over above 50 mg/mL | yes | zhang2024 |
| **mixed disulfide with methanethiol** | MFT + MeSH → MMFT | as above, 115 °C, methionine present | MMFT/dimer ratio falls from about 1.9 to 0.2 as cysteine rises; two-regime kinetics with a break at 90 min; units of the printed rate constants never stated | not tested | yes for MFT; no FFT analogue | zhang2024; kumazawa2003 (FD of FFT-methyl disulfide rises 4 → 20 in canned coffee) |
| **binding to melanoidins, covalent** | FFT + coffee melanoidin → bound thiol (thioether-type); very little disulfide (<6 of 400 µg) | 0.1 M phosphate pH 6, 30 °C, 12.5 g/L melanoidin | 80 % bound within 30 to 90 min and then a **plateau**: a capacity-limited sink; first-order 6 to 7e-4 /s over the first 30 min; in real brew at 80 °C the loss is not faster (0.023 /min) | plateau = saturation of sites; reversibility not tested | the thioether channel, with its pool never charged | hofmann2002 |
| **binding to melanoidins, non-covalent** | FFT ⇌ melanoidin π-complex, then slow covalent incorporation and the disulfide | 0.1 M phosphate pH 5.5, 5 mmol/L FFT, 6 / 27 / 60 °C, up to 96 h, NMR | free FFT after 5 min is 24 / 32 / 55 %: the reversible pool **shrinks with temperature** (ΔH about −19.5 kJ/mol); the slow covalent channel then takes everything at 27 and 60 °C but only 2 % at 6 °C; the dossier's covalent barrier 60 kJ/mol with a defensible range 7 to 102 | the fast pool is reversible and feeds the slow one; capacity per gram not measured | declared not implemented; its barrier band is the model's sink band | gigl2021 |
| **binding to protein disulfides** | FFT (and propanethiol) + BLG disulfides → 1:1 and 2:1 mixed disulfides | 1 % BLG, water, ambient, 12 g/L thiol | mass shifts only; the second exchange takes 6 to 24 h at 158 mmol/L thiol; ordinal "partial" at 22, 63, 72 and 130 °C, flat with temperature | asserted stable, never tested; stoichiometric in disulfides (two per protein) | yes, as a bracket, with the site pool zero everywhere | anantharamkrishnan2020b, yuan2023 |
| **thiol-Michael addition to an unsaturated carbonyl** | cysteine (and thiols generally) + acrylamide, HMF, enals → β-thioether | acrylamide + cysteine 80 to 180 °C; HMF + cysteine 5 to 50 °C pH 3.5 | barriers 28 to 30 kJ/mol (acrylamide + benzyl mercaptan / N-acetylcysteine, Hidalgo 2010 via the K6b ladder); HMF + cysteine second-order 4 / 5 / 23 M⁻¹ day⁻¹ at 5 / 25 / 50 °C, Ea 29.6 kJ/mol, HMF 97 % gone in 7 days at 50 °C | Hidalgo 2010 shows the thioether does **not** release on heating; HMF + Cys is 2:1 as well as 1:1 | acrylamide + cysteine only; nothing for MFT or FFT with any enal or with HMF | hamzalioglu2018, k6b_adduct_kinetics_synthesis, hidalgo2010 |
| **hemithioacetal, thiazolidine** | thiol + aldehyde; cysteine + aldehyde → thiazolidine (TTCA with a pentose) | 25 to 140 °C | TTCA is 94 % of the "Cys-Amadori" pool; its ring opening is fitted, its formation is charged, not modelled | equilibria | TTCA charged only; no thiol + aldehyde step for MFT/FFT | zhai2020, zhai2021 |
| **saliva and protein quench** | thiol + saliva proteins → undetectable | 22 °C, 45 min | detection limit rises 60 000-fold in crude saliva; mechanism unknown | not tested | no | starkenmann2008 |
| **survival** | MFT, FFT in the reference pot | ribose + cysteine, pH 5, 100 °C, autoclave | both thiols still rise between 6 and 12 h (MFT ×1.15, FFT ×1.20): removal is slow against formation at 100 °C | — | the model peaks at 1 h | schieberle2000 |

## 3. What this says about the sink the model needs

1. **Every measured thiol sink is weaker than the model's decay at cooking temperature, and every one has a partner.** Kumazawa's loss at 121 °C is set by pH and the oxidant in the can; Hofmann's and Gigl's binding plateaus when the melanoidin sites fill; the protein exchange is stoichiometric in disulfides; the Michael additions need an electrophile the pot must first make. The model's dominant sink has no partner and never runs out.
2. **The barrier band is not a barrier.** The only measured temperature dependence of a thiol sink is Gigl's covalent channel at 6 to 60 °C, and even there a single first-order fit fails inside the window. The model's 102 kJ/mol on the ceiling is the fit asking for a steeper slope than any sink shows, to reconcile 100 °C with 145 °C. A sink whose partner is made by browning rises with temperature through the partner, not through its own barrier.
3. **The reversible pool has the wrong sign for an Arrhenius sink.** Gigl's free thiol rises with temperature over the first minutes. No first-order sink with a positive barrier can do that. A reversible reservoir with a negative enthalpy can.
4. **The candidates, in order of evidence.** (a) A saturable covalent sink on a pool the pot makes (melanoidin-type electrophiles from the sugar branch), with Hofmann 2002's rate and plateau as the anchor and the pool charged from browning rather than as an input. (b) The disulfide made reversible or made second-order with a depletable oxidant, which the model already has in form, with Kumazawa's time-doubling test as the check. (c) Thiol-Michael addition to the enals and to HMF, with the 28 to 30 kJ/mol barrier and the HMF + cysteine rate as analogues, for the pots with fat or a hexose. (d) Gigl's reversible pool as a declared term with its measured enthalpy.
5. **The measurement that decides it** is the one in the introduction's section 8: the reference pot and the fed thiols on the same grid at 100 and 140 °C, with the disulfides quantified in the same run. If the fed thiol levels off with its disulfide, (b) is right; if it keeps falling while a browning marker rises, (a) is.

### 2b. The dry regime, from the same laboratory (read 2026-09-08)

Schieberle & Hofmann 1998 (`schieberle1998_extraction.md`) heat the same cysteine + ribose charge
dry on silica at 180 °C for 6 min and in 0.5 M phosphate at 145 °C for 20 min, both by stable
isotope dilution. The aqueous column is Hofmann 1998's Table 2 reprinted at tenfold scale (the
same experiment, not a new run; the fit must not enter it twice). What is new is where the sulfur
goes when the water leaves: the mercaptoketone falls from 59.9 to 10.1 µg per pot and the thiazine
from 42.4 to 1.0, while FFT rises eightfold, MFT by a third and 2-acetyl-2-thiazoline sevenfold; the
MFT to FFT ratio flips from 1.6 to 0.26. Four things change at once (water activity, temperature,
time, buffer), so this is a hold-out shape for any sink structure, not a coefficient: a sink that
scales with a browning-made pool should lose MORE thiol in the dry pot, where furfural rises from
5 to 7900 µg, and the pot loses less. The same chapter is the only 145 °C anchor for
3-mercapto-2-pentanone and the thiazine, the two compounds Schieberle 2000's 100 °C series follows.

The review read the same evening (`weerawatanakorn2015_extraction.md`) adds no number of its own;
seven of its ten thiol sources are on disk. The three it names that are not, and that bear on the
sink: Mottram, Szauman-Szumski & Dodson 1996 (thiol and disulfide loss to egg albumin at 100 °C,
the only protein sink at cooking temperature the corpus knows of), Hofmann, Czerny, Calligaris &
Schieberle 2001 (time-resolved thiol loss with coffee melanoidins), and Hofmann & Schieberle 1995
(2-acetyl-2-thiazoline loss in water at 100 and 145 °C and in oil).

## 4. What no dossier measures

An activation energy for any thiol sink printed by its authors; a time-resolved thiol loss in a
Maillard pot at 100 to 180 °C; any kinetics of an MFT sink (every rate is FFT's); thiol release from
a bound state other than Gigl's fast pool; oxygen dependence in water; the order in thiol; a
melanoidin binding capacity per gram; thiol + methanethiol kinetics in stated units.
