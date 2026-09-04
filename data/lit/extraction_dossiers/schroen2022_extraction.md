# Schroën & Berton-Carabin 2022 — A unifying approach to lipid oxidation in emulsions: modelling and experimental validation

### Per-paper extraction, 2026-09-04, from the full text supplied by the owner (`data/articles/Schroen2022.pdf`, 15 pages, open access CC BY-NC-ND). **Nothing in `src/`, `tests/` or `results/` was changed by this extraction; it records what the lipid lane already cites and whether the paper supports it.**

**Provenance codes:** **[M]** measured and printed · **[F]** fitted by the authors · **[D]** derived by this extraction · **[NEG]** verified negative.

---

## §0. IDENTITY — CONFIRMED, DOI ON THE PAGE

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/Schroen2022.pdf` (15 pages); a text dump `data/articles/schroen2022_fulltext.txt` predates it | read in full (pp. 1–8 in detail) |
| title | "A unifying approach to lipid oxidation in emulsions: Modelling and experimental validation" | p. 1 |
| authors | Karin Schroën (Wageningen; Twente) and Claire C. Berton-Carabin (Wageningen; INRAE UR BIA Nantes) | p. 1 |
| DOI | `10.1016/j.foodres.2022.111621` | p. 1 footer |
| journal | *Food Research International* 160 (2022) 111621 | p. 1 header |
| dates | received 22 March 2022; revised 11 June 2022; accepted 1 July 2022; available online 15 July 2022 | p. 1 |
| data source | the kinetics re-analyse emulsions published earlier (Berton, Genot, Viau et al. 2011): rapeseed oil O/W emulsions (30 g oil / 100 g), five emulsifiers, FeSO4/EDTA initiator 200 µM, incubated at **25 °C**, PIPES 10 mM / NaCl 80 mM buffer at **pH 6.7** | pp. 3–4 |

## §1. WHAT THE LIPID LANE CITES, AND WHAT THE PAPER PRINTS

| code site | claim in the code | what the paper prints | verdict |
|---|---|---|---|
| `src/kinetic_core/parameters_lipid.py` `K_LOOH_DECOMP_ANCHOR` (k_LOOH_decomp = 6.0e-3 /h, "Table 1, k4,CD — identical across all five emulsifiers") | first-order decomposition of the hydroperoxide pool to all secondary products | Table 1 **[F]**: **k4,CD = 6·10⁻³ h⁻¹ for β-lactoglobulin, BSA, β-casein, Tween 20 and Tween 80** (only k1,CD and k5,CD vary between emulsifiers, printed in italics). Reaction scheme p. 2: "(4) Hydroperoxide decomposition to form alkoxyl radicals, and in turn secondary oxidation products". Text p. 4: "k4 is a lumped reaction rate, ultimately leading to formation of secondary oxidation products, which includes formation of LO• radicals". | **confirmed**, including the "identical across all five" statement. |
| same file, `K_HEXANAL_SCHROEN` (k_hexanal = 6.0e-5 /h, "Table 1, k_hexanal — identical across all five emulsifiers") | hexanal-specific first-order constant, carried for the ratio 6e-5 / 6e-3 = 1.2 % | Table 1 **[F]**: **k_hexanal = 6·10⁻⁵ h⁻¹ for all five emulsifiers**; k_propanal = 3.5·10⁻⁴ h⁻¹ for all five. Text p. 4: "For propanal and hexanal formation, we use the same reaction as for the overall reaction for secondary oxidation products (equation VI), but use lower albeit constant k-values for all emulsions to do justice to the relative contribution of components". | **confirmed**; the 1.2 % ratio is the paper's own two numbers. |
| flags `hand_fitted_by_visual_agreement`, `no_standard_error_anywhere_in_the_source` | | p. 5, procedure step 5: **"The parameter values (k1–k5) were determined based on visual agreement of fit"**; no standard errors or confidence intervals are printed anywhere (Table 1 and the supporting Table A1 carry point values only). | **confirmed verbatim.** |
| flags `storage_temperature_25C_only`, `temperature_dependence_UNMEASURED` | | p. 3, 2.2.2: "Incubation was held at 25 °C"; the paper's only temperature discussion (p. 7) is a literature remark that constants would be "decreased by a factor of 2–3 for every 10 °C" relative to 37 °C sources. | **confirmed**: no measured temperature dependence. |
| conditions "rapeseed-oil O/W emulsion, droplets 1.4–1.8 µm, pH 6.7 buffer, rotative agitation, controlled oxygen-to-lipid ratio" | | pp. 3–4: droplet d3,2 1.4–1.8 µm (Table 1 last row: 1.5 / 1.8 / 1.7 / 1.4 / 1.7 µm); PIPES/NaCl buffer pH 6.7; 3 mL emulsion in 20.5 mL sealed headspace vials, 154 mmol O2 per kg oil. | **confirmed.** |
| `evidence_class="hand_fitted_no_uncertainty"`, `rate_transfer="licensed_as_an_anchor_only"` | | Consistent with §3.2.2 of the paper, which compares the fitted constants with literature ranges rather than deriving them from a measurement of the step. | appropriate. |

## §2. NUMBERS THE PAPER PRINTS THAT THE LANE DOES NOT USE

- Table 1 **[F]**: k1,CD 6.5e-5 / 2.7e-5 / 1.4e-5 h⁻¹ (β-lg / BSA / β-casein; not used for the surfactant emulsions); k2,CD 10 (mol/kg oil)⁻¹ h⁻¹; k3,CD 1 (mol/kg oil)⁻¹ h⁻¹; k5,CD 3.4e-3 / 2.5e-3 h⁻¹ (Tween 20 / Tween 80; zero for the proteins). Initial L• 500 / 10 / 100 / 1 / 1 µmol per kg oil.
- Other-hydroperoxide constants are 0.95 × the conjugated-diene ones (Table 1 footnote, reactivity ratio).
- Oxidizable PUFA 4.5 mol double bonds per kg oil; oxygen 154 mmol per kg oil.

## §3. VERDICT

Both constants the lipid lane takes from this paper are printed in Table 1 exactly as the code carries them, with the authors' own statement that they were fitted by visual agreement and with no uncertainty anywhere in the source. The lane's treatment (a floor anchor, no rate transfer, no temperature dependence) is the strongest claim the paper supports. No change to the code is warranted.
