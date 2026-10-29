# Hicks & Feather 1975 — EXTRACTION (three Amadori compounds — 1-deoxy-1-dibenzylamino-D-fructuronic acid, 1-benzylamino-1-deoxy-D-fructuronic acid, 1-benzylamino-1-deoxy-D-xylulose — heated at 100 °C at pH 7, in 2 N acetic acid at pH 2.5, and in 2 N H2SO4: norfuraneol, furfural and reductic acid by TLC/UV; ¹⁴C and D2O tracer on the norfuraneol route)

### QUALITATIVE MECHANISM PAPER WITH THREE NUMBERS: furfural yields of **16 / 23 / 41 %** from the three Amadori compounds in boiling 2 N H2SO4 (3 h, distilled), an upper bound of **< 6 %** norfuraneol from the secondary-amine Amadori compounds at pH 2.5 (1 h, 100 °C), and a ¹⁴C result placing C-1 of the Amadori sugar in norfuraneol's **methyl group (~89-96 % of the label)**. No rate constant, no yield at pH 7, no time course.

**Source on disk:** `data/articles/hicks1975.pdf` (4 pp.; p. 1 is the tail of the preceding paper,
JAFC-75-950, and p. 4 the head of the following one, Milić et al.). The paper itself is J. Agric. Food
Chem. 23 (5), 957-959. Read from the `pdftotext -layout` text layer; the results page (p. 959) was
checked against a 130 dpi render (scratchpad `fork_nf/hicks-3.png`). No tables; Figs. 1-3 are
structures.

## 0. Identity

| field | value |
|---|---|
| Title | "Studies on the Mechanism of Formation of 4-Hydroxy-5-methyl-3(2H)-furanone, a Component of Beef Flavor, from Amadori Products" |
| Authors | Kevin B. Hicks, **Milton S. Feather** (corresponding) — Department of Biochemistry, University of Missouri, Columbia |
| Venue | J. Agric. Food Chem. 1975, **23** (5), 957-959. Received 15 Nov 1974, accepted 12 May 1975. Missouri Agric. Exp. Station journal paper 7146 |
| DOI | **not printed**. Brief's 10.1021/jf60201a005 is consistent with JAFC 23(5) p. 957. **Identity verified** on title, authors, journal, volume, issue, pages. |
| Compound numbers | **1** = norfuraneol (4-hydroxy-5-methyl-3(2H)-furanone, engine `NF`); **2** = 1-deoxy-1-dibenzylamino-D-fructuronic acid (tertiary amine, glucuronic-acid Amadori); **3** = 1-benzylamino-1-deoxy-D-fructuronic acid (secondary amine); **4** = 1-benzylamino-1-deoxy-D-xylulose (secondary amine; **the xylose Amadori compound**, used as its oxalate salt) |

## 1. Why it matters

Gaps (c) and (d). The engine's pentose trunk splits by pH into `r_pent_tdp` (acid-catalysed
1,2-enolisation → 3-deoxypentosone → furfural) and `r_pent_dpo` / `r_arp_dpo` (2,3-enolisation →
1-deoxypentosone → `r_dpo_nf` → NF). `SULFUR_FLUX_MAP.md` §5 asks whether that pH law is too steep
at pH 3 (Hofmann's MFT is highest at pH 3, the model sends 0.4 % to DPO). Hicks & Feather is the
classic experimental statement of the **Anet rule**: protonation of the Amadori nitrogen promotes
1,2-enolisation (furfural, reductic acid), an unprotonated nitrogen allows 2,3-enolisation (NF);
hence **amine basicity** decides the split at intermediate acidity. Their data:

- **2 N H2SO4:** all three give furfural; NF "either undetectable or found in trace amounts".
- **pH 2.5 (2 N acetic acid):** NF "in major amounts only from 2" (the tertiary, weaker-base amine);
  from **3** and **4** (secondary amines, incl. the xylose Amadori compound) NF **< 6 %** (TLC
  detection limit).
- **pH 7:** "all three compounds gave 1 as the major identifiable product". No yield.
- **Tracer:** C-1 of the Amadori sugar → NF's 5-methyl; methyl H fully exchanged in D2O during
  reaction (not after) → a 1-deoxy-2,3-dicarbonyl (1-deoxyosone) intermediate in extensive
  equilibrium with its enol. This is the DPO → NF step the engine uses.

## 2. Methods as they matter to a model

| arm | substrate charge | medium | T, t | analysis |
|---|---|---|---|---|
| NF preparation | 10 g of **2** in 180 mL | 2 N acetic acid, "which had a pH of 6.8" | reflux, 1 h | isolation, sublimation (mp 129 °C); **no yield given** |
| pH 7 | ~4 g of **2**, **3**, **4** in 50 mL (~80 g/L) | "pH 7 buffer" (**buffer identity not stated**) | 100 °C, 1 h | Dowex 50, CHCl3 extraction, TLC, MS — identification only |
| pH 2.5 | ~0.2 g in 5 mL (~40 g/L) | 2 N acetic acid, adjusted to pH 2.5 with glacial acetic acid | 100 °C, 1 h | 2 µL spotted on TLC directly; detection-limit → "< 6 %" |
| strong acid, products | ~2 g in 25 mL (~80 g/L) | 2 N H2SO4 | 100 °C, 45 min | TLC, UV (reductic acid, NF, 5-formyl-2-furoic acid) |
| strong acid, furfural | ~0.5 g in 25 mL (~20 g/L) | 2 N H2SO4, distilled until ~80 % over | boiling, ~3 h | UV (λmax 278 nm) + phenylhydrazone (mp 94 °C) |
| ¹⁴C | glucuronic acid-1-¹⁴C → **2** (0.027 µCi/mmol) → NF | as NF preparation | reflux 1 h | Kuhn-Roth oxidation to acetic acid, scintillation |
| D2O | ~10 g of **2** in 160 mL 98 % D2O, pH 6.8 with 20 mL acetic acid | — | 100 °C, 1 h | ¹H NMR integrals |

## 3. Numbers extracted (p. 959)

| quantity | value | conditions | locator | class |
|---|---|---|---|---|
| furfural yield from **2** | **16 %** | 2 N H2SO4, boiling, distilled over ~3 h | p. 959 col. 1 | measured (UV, spectrophotometric; basis — molar — not stated) |
| furfural yield from **3** | **23 %** | same | same | measured |
| furfural yield from **4** (xylose Amadori) | **41 %** | same | same | measured |
| NF from **3**, **4** at pH 2.5 | **< 6 %** | 2 N acetic acid pH 2.5, 100 °C, 1 h | p. 959 col. 1 | detection-limit bound |
| NF from **2** at pH 2.5 | "major amounts" | same | p. 959 | qualitative |
| NF at pH 7 | "major identifiable product" from 2, 3, 4 | pH 7 buffer, 100 °C, 1 h | p. 959 | qualitative |
| reductic acid in 2 N H2SO4 | from 2 and 3, not from 4 | 45 min, 100 °C | p. 959 | qualitative |
| ¹⁴C specific activity | 0.027 (Amadori **2**) → 0.026 (NF) → 0.024 µCi/mmol (acetic acid from C-5 methyl + C-4) | — | p. 959 | measured → **≥ 89 %** (0.024/0.027) of C-1 label in the methyl-bearing C2 fragment; 96 % retention into NF |
| D exchange | CH3 : H-2 = 3 : 2 (H2O product) → ≈ 0 : 1 (D2O product) | 100 °C, 1 h, pH 6.8 | p. 959 | measured; product itself does not exchange CH3 under the same conditions |

**Not reported:** any rate constant; any NF yield (at any pH); conversion of the Amadori compounds;
time courses; temperature other than 100 °C; buffer composition at pH 7; furfural yields at pH 2.5 or 7.

## 4. What this can and cannot anchor

**Can:**
1. **Structure of the NF route**: Amadori → 2,3-enol → 1-deoxyosone ⇌ enol → NF, with C-1 → CH3.
   This supports the engine's `DPO → NF` topology (`r_dpo_nf`), and the carbon mapping a ¹³C hold-out
   would use.
2. **The acid limit of the 1,2-branch from a pentose Amadori**: 41 % furfural from the xylose
   Amadori in boiling 2 N H2SO4 — a ceiling-type reference for `r_arp_tdp` → `r_tdp_fur` at
   pH ≈ 0 (not a rate).
3. **Direction of the amine-basicity effect**: weaker-base (tertiary) Amadori nitrogen → more NF at
   pH 2.5; secondary-amine Amadoris (the class cysteine's Amadori compound belongs to) → NF < 6 %.

**Cannot:**
- Give any rate or Ea; the pH law is three qualitative points (pH ~0, 2.5, 7).
- Be read as a branch fraction: "< 6 %" is a product yield after 1 h with unknown conversion.
- Transfer to cysteine without care: benzylamine/dibenzylamine Amadori compounds, no thiol, no
  thiazolidine; the glucuronic-acid Amadoris decarboxylate to reach C5.

## 5. Flags

1. Furfural yields are spectrophotometric from a **distillate** (~80 % of volume taken over), so they
   lump formation and the distillation recovery; basis (mol % of Amadori) is implied, not stated.
2. The pH 2.5 bound rests on an unspecified TLC detection limit ("from previous studies").
3. In 2 N H2SO4 "Only from compound 3 could furanone be detected" (45-min run) — the abstract's
   statement that NF forms at pH 2.5 "as a major product from 2" refers to the acetic-acid arm. Do
   not conflate the two acid arms.
4. "pH 7 buffer" composition not given — could be phosphate (catalytic) or not.
5. NF preparation pot is 10 g of **2** in 2 N acetic acid at a measured **pH 6.8** (the amine buffers
   it), not pH ~2.4 as 2 N acetic acid alone would be.

## 6. Use in the model

- **Gap (c), the pH 3 MFT problem (Hofmann pH 3: model MFT 22× low because DPO gets 0.4 %).** Hicks
  & Feather do **not** contradict a steep 2,3-enolisation pH law for a secondary-amine Amadori: at pH
  2.5 the xylose Amadori compound gives < 6 % NF. That bound is loose (a 0.4 % model branch is
  consistent with it), so it neither confirms nor refutes the model; it does say that a pH-3 MFT
  source via NF from a secondary-amine Amadori is not expected to be large, which shifts weight toward
  the other two candidates in §5 (phosphate catalysis; a non-NF MFT route) or toward a tertiary-like
  (less basic) nitrogen — e.g. the thiazolidine (TTCA) ring nitrogen, pKa 6.24, which is far less
  basic than an α-amino group. **That last point is an inference, not in the paper.**
- **Gap (d).** Supports the DPO → NF mechanism (C-1 → CH3, enol equilibration). No numbers.
- **Furfural acid ceiling.** 41 % furfural from the xylose Amadori in 2 N H2SO4 could serve as a
  qualitative bound on `r_arp_tdp`'s furfural throughput at the acid extreme; outside every fit pot's
  pH range, so not a hold-out.
- **Caveats:** non-food amines, 100 °C only, TLC-level quantitation, distillation-dependent furfural yields.
