# Li et al. 2020 — EXTRACTION (hexanal in a soy-protein emulsion; the temperature axis is on the protein, not on the reaction — REFUSED as a rate source)

**Source on disk:** `data/articles/li2020.pdf` (Elsevier accepted manuscript, 49 pp.; downloaded
2026-09-11 at this repository's request). Read 2026-09-11 via `pdftotext -layout` plus a 400 dpi
render of the figure page. Wave B45.

| field | value |
|---|---|
| Title | "Impact of heating treatments on physical stability and lipid-protein co-oxidation in oil-in-water emulsion prepared with soy protein isolates" |
| Authors | Qingyun Li, Jiabao Zheng, Ge Ge, Mouming Zhao, Weizheng Sun (South China University of Technology) |
| Venue | Food Hydrocolloids — **accepted manuscript; volume, pages and article number are NOT printed on this version and are not inferred here** |
| DOI | 10.1016/j.foodhyd.2019.06.012 (printed on the cover sheet; PII S0268-005X(19)30971-3) |
| System | 10 % (w/w) corn oil in 90 % (w/w) aqueous phase; soy protein isolate 20 g/L, 10 mM phosphate buffer, **pH 7**, 0.02 % sodium azide |
| Reaction temperature | **37 °C, one value, for every sample** |

## 1. Why this paper was fetched, and why it does not answer

The reading list asked for a **hexanal temperature dependence in a wet protein system** — the missing
rung between the 114–122 kJ/mol bulk-oil barrier and the 61–65 kJ/mol moist-nut-paste barrier the
lipid lane already carries. A Google Scholar Labs snippet suggested this paper heated a soy-protein
emulsion at several temperatures and measured hexanal.

**It does not.** The 70 / 90 / 120 °C axis is applied to the **dry-weight SPI dispersion for 15 min
before emulsification**; the protein is then centrifuged, diluted, and used as an emulsifier. Every
emulsion then oxidises at **37 °C**. Verbatim, §2.2:

> "Three aliquots of native SPI dispersions were heated for 15 min using a water bath (at 70°C and
> 90°C) or autoclaved at 120°C, respectively. Native and heated SPI solutions were centrifuged
> (10,000 g for 20 min) and diluted to 20 g/L."

and §2.3: "The resulting fine emulsion containing 0.02% sodium azide were stored at 37°C for the
following analysis." Emulsification itself is mechanical only (8,000 rpm shear, then 30 MPa
high-pressure homogenisation) — no thermal step touches the emulsion.

So the four numbers are a **protein-denaturation axis at fixed reaction temperature**, not a
temperature axis. Reading them as an Arrhenius series would be a category error, and the data say so
themselves: the largest hexanal belongs to the **unheated** protein and the series is non-monotonic.

## 2. Fig. 5B — the four printed values, transcribed

Caption verbatim: "Lipid oxidation of SPI-stabilized emulsions during storage. (A) Hydroperoxide
content and (B) hexanal level (using 2-methyl-3-heptanone as internal standard). Abbreviations
native SPI, 70-SPI, 90-SPI, and 120-SPI represent the emulsions prepared with native, 70 °C heated,
90 °C heated, and 120 °C heated SPI, respectively."

| bar | SPI pre-treatment | reaction T | time | hexanal, "Internal stardard equivalent" (*sic*, as printed) |
|---|---|---:|---:|---:|
| Native SPI | none | 37 °C | 21 d | **16.3** |
| 70-SPI | 70 °C, 15 min | 37 °C | 21 d | **2.8** |
| 90-SPI | 90 °C, 15 min | 37 °C | 21 d | **3.2** |
| 120-SPI | 120 °C, autoclave | 37 °C | 21 d | **5.1** |

The labels are printed above the bars and were read at 400 dpi; they are transcribed, not estimated.
**No error bars are drawn on panel B and no ± values are printed**, although §2.8 claims triplicates.

**The units are not concentrations.** HS-SPME/GC-MS with a single internal standard
(2-methyl-3-heptanone) and **no external calibration and no standard addition** is described; the
y-axis is a dimensionless internal-standard ratio. Under this repository's rule these are one step
above peak areas and **must not be scored as a concentration**.

**One time point only** — 21 days. A single time cannot give a rate even at 37 °C.

## 3. The hydroperoxide pool — figure only

Fig. 5A is "Hydroperoxide content (µmol/kg)" against storage day (ticks 1–17, all at 37 °C), four
series with error bars, calibrated against cumene hydroperoxide. **No value is printed in the text or
in any table**, so nothing is transcribable. No peroxide value, no conjugated dienes, no TBARS.

## 4. What the repository can bank

One qualitative fact, and it is worth carrying: in this hot-processed protein-stabilised emulsion the
hydroperoxide pool and the hexanal output move in **opposite directions** — the emulsion with the
**lowest** hydroperoxide content produced the **highest** hexanal. The authors' own reading, verbatim:

> "This demonstrated that the lowest hydroperoxide content of native SPI-stabilized emulsion resulted
> from degradation of hydroperoxides."

For a lipid lane that routes a hydroperoxide pool into hexanal, this is a warning that the two are
not in fixed ratio and that the emulsifier's interfacial state — not temperature alone — gates the
decomposition step. It is a **structure claim, not a number**, and is recorded as such.

## 5. Kinetics

**None.** No rate constant, no activation energy, no half-life, no Arrhenius treatment, no kinetic
model (grep-verified for arrhenius / activation energ / rate constant / half-life / kinetic / Q10 /
first-order / zero-order — zero hits). The only statistics are Duncan's multiple range test at
p < 0.05.

## 6. Verdict

**REFUSED as a rate or barrier source**, on three independent grounds, any one of which is
sufficient: a single reaction temperature; a single time point; relative internal-standard units
rather than concentrations. Kept as a structure-only corroboration of the pool-versus-aldehyde
decoupling described in §4.
