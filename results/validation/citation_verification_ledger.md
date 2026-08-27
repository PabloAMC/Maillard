# Citation verification ledger

**Generated:** 2026-08-26 · **Branch:** `audit-remediation` · **Source of truth:** CrossRef REST API (`api.crossref.org`) plus CrossRef bibliographic title search.

Machine-readable companion: `citation_verification_ledger.json` (same directory).

## What this is

Every DOI used as an *anchor* in the repo's literature and benchmark registries was extracted together with the
context that surrounds it — the repo id, the citation string, the chemistry/matrix family, and the note describing
what the anchor is supposed to supply. CrossRef metadata was then fetched for each DOI and compared against that
context. This catches the class of failure that no link-checker finds: **a DOI that resolves perfectly well, to
the wrong paper.**

| Class | Meaning |
|---|---|
| `MATCH` | The DOI resolves to a paper that plausibly supports the claim made at the anchor site. |
| `METADATA-MISMATCH` | Right paper, wrong bibliographic claim — wrong author, year, venue, or an internally conflicting `paper_id`. |
| `TOPIC-MISMATCH` | The DOI is live and resolves to a real paper **on a different subject** than the anchor claims. |
| `DEAD` | The DOI does not resolve (CrossRef HTTP 404). |

## Headline

**225 unique DOIs** across **421 anchor sites** in 8 file globs.

| Class | Unique DOIs | Share | Anchor sites |
|---|---:|---:|---:|
| `MATCH` | 88 | 39% | 168 |
| `METADATA-MISMATCH` | 47 | 21% | 93 |
| `TOPIC-MISMATCH` | 45 | 20% | 80 |
| `DEAD` | 45 | 20% | 80 |

> **137 of 225 anchors (61%) are defective in some way.**
> The prior audit sampled 28 anchors and estimated ~20% dead DOIs; that dead rate is confirmed exactly (45/225 = 20%),
> but the full sweep shows an equally large *live-DOI-wrong-paper* population (45 TOPIC-MISMATCH) that sampling
> under-weighted, plus 46 entries where the paper is right but the repo's bibliographic claim about it is not.

**No numeric value anywhere in the repo was changed by this pass.**

### Where the damage is concentrated

| File | MATCH | METADATA | TOPIC | DEAD |
|---|---:|---:|---:|---:|
| `data/lit/benchmark_intake_registry.json` | 72 | 40 | 39 | 41 |
| `data/lit/slr_incorporation_matrix.json` | 66 | 38 | 30 | 29 |
| `data/lit/computational_priors.json` | 3 | 4 | 1 | 1 |
| `data/lit/arrhenius_params.yml` | 4 | 1 | 3 | 1 |
| `data/benchmarks/quarantined/cys_ribose_150C_Mottram1994.json` | 1 | 0 | 1 | 0 |
| `data/benchmarks/spi_hvp_xylose_120C_PMC9905368.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/wheat_gluten_hvp_xylose_120C_PMC9905368.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/quarantined/cys_glucose_150C_Farmer1999.json` | 0 | 0 | 0 | 1 |
| `data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/external_validation/external_validation_bi_2020_raw_pea_hexanal.json` | 0 | 0 | 1 | 0 |
| `data/benchmarks/external_validation/external_validation_bi_2020_roasted_pea_hexanal.json` | 0 | 0 | 1 | 0 |
| `data/benchmarks/acrylamide_asparagine_glucose_Parker2012.json` | 0 | 0 | 0 | 1 |
| `data/benchmarks/thiamine_cys_xylose_145C_Cerny2008.json` | 0 | 0 | 1 | 0 |
| `data/benchmarks/quarantined/thiamine_cys_ribose_100C_Hofmann1996.json` | 0 | 0 | 0 | 1 |
| `data/benchmarks/cys_ribose_140C_Hofmann1998.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/furosine_extrusion_crossover_140C_RamirezJimenez2000.json` | 0 | 1 | 0 | 0 |
| `data/benchmarks/cml_cel_commercial_pbma_Foods2023.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/external_validation/external_validation_li_2026_spi_wg_hme_control.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/pea_isolate_40C_PratapSingh2021.json` | 1 | 0 | 0 | 0 |
| `data/benchmarks/soy_isolate_40C_PratapSingh2021.json` | 1 | 0 | 0 | 0 |

Two rows above need a caveat: `10.1021/acs.jafc.9b07711` (Bi et al. 2020) is the **correct** anchor for both
`external_validation_bi_2020_*` benchmarks and for `acs_2020_raw_pea_hexanal_baseline`; it is counted as
TOPIC-MISMATCH only because the same DOI is *also* reused in the intake registry for an unrelated EGCG claim.
Likewise the `cys_ribose_150C_Mottram1994` quarantine file scores a MATCH because its own audit note correctly
cites the candidate replacement paper.

---

## 1. TOPIC-MISMATCH — live DOI, wrong paper

The most serious class: the citation resolves, so every automated check passes, but the paper it resolves to
cannot support the claim. Ordered worst-first (wrong discipline entirely, then wrong subject, then partial).

### Top 10 worst

**1. `10.1016/j.mineng.2010.05.003`** — used as `li_2010_phytate_chelation_kinetics`  
*Actually is:* Smith, Hadler et al. (2010), "Flotation bank air addition and distribution for optimal performance", *Minerals Engineering*  
DOI is Smith, Hadler & Cilliers 2010, 'Flotation bank air addition and distribution for optimal performance' in MINERALS ENGINEERING. Repo uses it as 'li_2010_phytate_chelation_kinetics'. Wrong discipline entirely.

**2. `10.3390/plants13020274`** — used as `mdpi_plants_2024_hemp_volatiles`  
*Actually is:* Chen, Zhang et al. (2024), "RNA-Seq-Based WGCNA and Association Analysis Reveal the Key Regulatory Module and Genes Responding to Salt Stress in Wheat Roots", *Plants*  
DOI is Chen et al. 2024, 'RNA-Seq-Based WGCNA ... Salt Stress in WHEAT ROOTS' (Plants). Repo uses it as 'mdpi_plants_2024_hemp_volatiles' (hemp protein off-note flavour reference). Wrong field entirely.

**3. `10.1021/acs.jafc.9b06882`** — used as `ref41_ppi_sulfur_binding`  
*Actually is:* Zhao, Li et al. (2019), "Lycopene Prevents DEHP-Induced Leydig Cell Damage with the Nrf2 Antioxidant Signaling Pathway in Mice", *Journal of Agricultural and Food Chemistry*  
DOI is Zhao et al. 2019, 'Lycopene Prevents DEHP-Induced Leydig Cell Damage with the Nrf2 Antioxidant Signaling Pathway in Mice' - a rodent reproductive-toxicology paper. Repo uses it as 'ref41_ppi_sulfur_binding' (pea-protein sulfur-volatile binding order prior). Wrong field entirely.

**4. `10.1016/j.carbpol.2022.120468`** — used as `tran_2023_reducing_sugar_hme`  
*Actually is:* Liu, Xu et al. (2023), "Self-healing, antibacterial, and conductive double network hydrogel for strain sensors", *Carbohydrate Polymers*  
DOI is Liu et al. 2023, 'Self-healing, antibacterial, and conductive double network hydrogel for strain sensors' (Carbohydr Polym). Repo uses it as 'tran_2023_reducing_sugar_hme' / 'starch reducing-sugar liberation during HME'. Completely unrelated field.

**5. `10.3390/foods7080126`** — used as `pmc6104182_soybean_fermentation`  
*Actually is:* Chambers, Godwin et al. (2018), "Recipes for Determining Doneness in Poultry Do Not Provide Appropriate Information Based on US Government Guidelines", *Foods*  
DOI is Chambers, Godwin & Terry 2018, 'Recipes for Determining Doneness in POULTRY Do Not Provide Appropriate Information Based on US Government Guidelines'. Repo uses it as 'pmc6104182_soybean_fermentation'. Unrelated.

**6. `10.1016/j.foodchem.2011.07.037`** — used as `wang_2012_gsh_xylose_sulfur_uplift`  
*Actually is:* Koblar, Tavčar et al. (2012), "Fluoride in teas of different types and forms and the exposure of humans to fluoride with tea and diet", *Food Chemistry*  
DOI is Koblar et al. 2012, 'Fluoride in teas of different types and forms and the exposure of humans to fluoride with tea and diet'. Repo uses it as 'wang_2012_gsh_xylose_sulfur_uplift' (glutathione/xylose sulfur-uplift prior). Unrelated.

**7. `10.3390/foods12061349`** — used as `pmc10056349_rubisco_amadori`  
*Actually is:* Simkova, Veberic et al. (2023), "Variability in ‘Capri’ Everbearing Strawberry Quality during a Harvest Season", *Foods*  
DOI is Simkova et al. 2023, 'Variability in Capri Everbearing STRAWBERRY Quality during a Harvest Season'. Repo uses it as 'pmc10056349_rubisco_amadori' (Rubisco amino acid composition and Amadori kinetics). Unrelated.

**8. `10.1016/j.lwt.2021.111352`** — used as `huang_2021_sulfur_oav_support`  
*Actually is:* Chraibi, Fadil et al. (2021), "Antimicrobial combined action of Mentha pulegium, Ormenis mixta and Mentha piperita essential oils against S. aureus, E. coli and C. tropicalis: Application of mixture design methodology", *LWT*  
DOI is Chraibi et al. 2021, antimicrobial mixture-design study of Mentha/Ormenis essential oils vs S. aureus, E. coli, C. tropicalis. Repo uses it as 'huang_2021_sulfur_oav_support' (MFT/FFT sulfur prominence prior). Unrelated.

**9. `10.1021/jf00080a032`** — used as `arabshahi_1988_aw_thiamine_kinetics`  
*Actually is:* Schaefer, Sandermann (1988), "Metabolism of pentachlorophenol in cell suspension cultures of wheat (Triticum aestivum L.). Tetrachlorocatechol as a primary metabolite", *Journal of Agricultural and Food Chemistry*  
DOI is Schaefer & Sandermann 1988, metabolism of PENTACHLOROPHENOL in wheat cell suspension cultures. Repo uses it as 'arabshahi_1988_aw_thiamine_kinetics' (water-activity-dependent thiamine kinetics). Unrelated.

**10. `10.1021/jf00058a011`** — used as `glomb_1995_3dg_fragmentation_stoichiometry`  
*Actually is:* Cabras, Garau et al. (1995), "Fate of Some Insecticides from Vine to Wine", *Journal of Agricultural and Food Chemistry*  
DOI is Cabras et al. 1995, 'Fate of Some Insecticides from Vine to Wine'. Repo uses it as 'glomb_1995_3dg_fragmentation_stoichiometry' (3-deoxyglucosone fragmentation prior). Unrelated.

### Remaining topic-mismatches

| DOI | Used in repo as | The DOI actually is |
|---|---|---|
| `10.1002/recl.19871060202` | `de_bruijn_1987_monosaccharide_alkaline_degradation` | Chittenden, Regeling (1987), "An improved procedure for the synthesis of sugar chloroacetates and some substitution reactions thereof" |
| `10.1006/fstl.1996.0092` | `ilo_1996_maize_sme_lysine_damage` | Ilo, Tomschik et al. (1996), "The Effect of Extrusion Operating Conditions on the Apparent Viscosity and the Properties of Extrudates in Twi" |
| `10.1016/S0260-8774(00)00047-9` | `quintas_2000_sucrose_caramelisation` | Ramaswamy, Zareifard (2000), "Evaluation of factors influencing tube-flow fluid-to-particle heat transfer coefficient using a calorimetric t" |
| `10.1016/j.foodchem.2013.11.082` | `hrncirik_2014_coconut_oil_thermal_profile` | Mashhadizadeh, Amoli-Diva et al. (2014), "Solid phase extraction of trace amounts of silver, cadmium, copper, mercury, and lead in various food samples " |
| `10.1016/j.foodchem.2014.05.097` | `ordoudi_2014_hmf_peak_window` | Agcam, Akyıldız et al. (2014), "Effects of PEF and heat pasteurization on PME activity in orange juice with regard to a new inactivation kinet" |
| `10.1016/j.foodchem.2020.128670` | `comert_gokmen_2021_epicatechin_cysteine_mgo_synergy` | Azevedo, Cerqueira et al. (2021), "Rhamnolipids-based nanostructured lipid carriers: Effect of lipid phase on physicochemical properties and stab" |
| `10.1016/j.foodchem.2023.136009` | `yeo_mottram_2023_lecithin_crosstalk_anchor` | Bae, Choi et al. (2023), "Ethanol organosolv lignin as a substitute for commercial antioxidants, focusing on the structural properties a" |
| `10.1016/j.foodchem.2023.136200` | `poulsen_2023_pbma_cml_cel` | Weligama Thuppahige, Moghaddam et al. (2023), "Investigation of critical properties of Cassava (Manihot esculenta) peel and bagasse as starch-rich fibrous ag" |
| `10.1016/j.foodchem.2025.145815` | `siripitakpong_2026_fft_retention` | Sharma, Keast et al. (2025), "Plant-protein isolates and flavour perception: Understanding mechanisms and strategies to balance flavour rete" |
| `10.1016/j.foodhyd.2024.110509` | `malia_2025_pea_free_sh_crosscheck` | Gong, Xu et al. (2025), "Effects of esterification and enzymatic modification on the properties of wheat starch and dough" |
| `10.1016/j.jfoodeng.2014.06.032` | `luna_aguilera_2014_molten_sugar_color_kinetics` | Gelaw, Güell et al. (2014), "Use of attenuated total reflectance infrared microspectroscopy combined with multivariate analysis to study me" |
| `10.1016/j.lwt.2006.07.014` | `mundt_wedzicha_2007_biscuit_browning` | Kavitha, Modi (2007), "Effect of water activity and temperature on degradation of 5′-inosine monophosphate in a meat model system" |
| `10.1016/j.lwt.2021.111802` | `yu_2021_corn_hydrolysate_kinetics` | Patel, Andrade et al. (2021), "Fat crystal-stabilized water-in-oil emulsion breakdown and marker release during in vitro digestion" |
| `10.1016/j.lwt.2024.117316` | `pruteanu_2023_glucose_phe_arrhenius_browning` | Bolchini, Angeli et al. (2025), "Free radical scavenging kinetics of Maillard reaction products: A glucose-glycine model system" |
| `10.1016/j.meatsci.2020.108151` | `sun_2020_solid_matrix_cml_cel_accumulation` | Jo, Lee et al. (2020), "Utility of winter mushroom treated by atmospheric non-thermal plasma as an alternative for synthetic nitrite a" |
| `10.1021/acs.jafc.9b07711` | `acs_2020_raw_pea_hexanal_baseline`, `jafc_2020_egcg_deoxyosone_trapping`, `external_validation_bi_2020_raw_pea_hexanal`, `external_validation_bi_2020_roasted_pea_hexanal` | Bi, Xu et al. (2020), "Characterization of Key Aroma Compounds in Raw and Roasted Peas (<i>Pisum sativum</i> L.) by Application of In" |
| `10.1021/acsomega.7b00321` | `nashalian_yaylayan_2014_cu_catalyzed_strecker` | Jefferson, Hu et al. (2017), "New Insight into and Characterization of the Aqueous
Metal-Enol(ate) Complexes of (Acetonedicarboxylato)copper" |
| `10.1021/bi00270a010` | `beta_elimination_dha` | Helmerhorst, Stokes (1983), "Generation of an acid-stable and protein-bound persulfide-like residue in alkali or sulfhydryl-treated insulin" |
| `10.1021/jf00002a024` | `yeo_shibamoto_1991_wof_hexanal` | Ramarathnam, Rubin et al. (1991), "Studies on meat flavor. 1. Qualitative and quantitative differences in uncured and cured pork" |
| `10.1021/jf00009a012` | `blank_grosch_1991_hdmf_anchor` | Mattheis, Buchanan et al. (1991), "Change in apple fruit volatiles after storage in atmospheres inducing anaerobic metabolism" |
| `10.1021/jf00049a015` | `cys_ribose_150C_Mottram1994`, `$.source_metadata.quarantine_reason` | Zarkadas, Yu et al. (1995), "Assessment of the Protein Quality of Beefstock Bone Isolates for Use as an Ingredient in Meat and Poultry Prod" |
| `10.1021/jf00065a009` | `tannenbaum_1985_thiamine_ea` | Al-Wandawi, Abdul-Rahman et al. (1985), "Tomato processing wastes as essential raw materials source" |
| `10.1021/jf00083a026` | `nakamura_1988_imp_ribose_release` | Epstein, Randecker et al. (1988), "Influence of heat and cure preservatives on residues of sulfamethazine, chloramphenicol, and cyromazine in mus" |
| `10.1021/jf00111a008` | `grosch_1982_hexanal_odt` | Purcell, Walter (1982), "Stability of amino acids during cooking and processing of sweet potatoes" |
| `10.1021/jf026123f` | `blank_devaud_grosch_2003_g6p_hdmf_uplift_v1`, `blank_devaud_grosch_2003_g6p_hdmf_prior` | Cerny, Davidek (2003), "Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction" |
| `10.1021/jf047903m` | `martins_van_boekel_2005_ascorbic_amino_browning` | Adams, Borrelli et al. (2005), "Thermal Degradation Studies of Food Melanoidins" |
| `10.1021/jf800268t` | `thiamine_cys_xylose_145C_Cerny2008` | Limacher, Kerler et al. (2008), "Formation of Furan and Methylfuran by Maillard-Type Reactions in Model Systems and Food" |
| `10.1080/10408398.2017.1378865` | `retro_aldol` | ALjahdali, Carbonero (2017), "Impact of Maillard reaction products on nutrition and health: Current knowledge and need to understand their f" |
| `10.1111/j.1365-2621.1987.tb14251.x` | `buera_1987_maillard_caramelisation_ea` | KOCH, SEIB et al. (1987), "Incorporating L‐Ascorbyl 6‐Palmitate in Bread and Its Shortening‐Sparing and Anti‐Firming Effects" |
| `10.3390/foods10123140` | `karolkowski_2021_ppi_ph_release` | Karolkowski, Guichard et al. (2021), "Volatile Compounds in Pulses: A Review" |
| `10.3390/foods13091393` | `xu_2024_soybean_pbma_hexanal` | Spizzirri, Restuccia (2024), "Advances in Food Waste Biomass Transformation into High-Value Products" |
| `10.3390/foods13162590` | `pmc11353891_lentil_deflavoring` | Streule, André et al. (2024), "Influences of Depulping, Pod Storage and Fermentation Time on Fermentation Dynamics and Quality of Ghanaian Co" |
| `10.3390/foods14111881` | `pmc12154226_crosspy_thiol_adduction` | El Hosry, Elias et al. (2025), "Maillard Reaction: Mechanism, Influencing Parameters, Advantages, Disadvantages, and Food Industrial Applicati" |
| `10.3390/foods14193453` | `wang_2024_glucosamine_synergy` | Wu, Wang et al. (2025), "Physicochemical and Flavor Characteristics of Maillard Reaction Products from Nile Tilapia Fish Skin Collagen " |
| `10.5194/egusphere-2026-321` | `cysteine_thermolysis` | Wu, Jaccard et al. (2026), "Technical note: Kinetically resolved volatile and redox fingerprints of geologic materials by TGA/DSC-MicroGC" |

---

## 2. METADATA-MISMATCH — right paper, wrong bibliographic claim

These do not invalidate the science, but they make the provenance untrustworthy: a reader chasing the cited
author-year will not find the paper, and several are cases where one DOI has been split across multiple
supposedly-distinct repo ids.

| DOI | Repo claims | CrossRef says | Issue |
|---|---|---|---|
| `10.1002/jsfa.10528` | Zhu et al. (2021) | Zhu, Fang et al. (2020) | Right paper (Zhu et al., CML/CEL kinetics in chicken braising) but published 2020; repo cites 'Zhu et al. (2021)'. |
| `10.1007/s10068-016-0185-5` | Chen et al. (2016) | Yu, Gao et al. (2016) | Right topic (CML/CEL formation in meat products) but authors are Yu, Gao, Zeng et al. 2016; repo cites 'Chen et al. (2016)'. |
| `10.1007/s10068-024-01633-w` | Choi et al. (2024), roasted garden peas | Jung, Baek et al. (2024) | Right topic (acrylamide in air-fryer roasted legumes vs asparagine/free sugars) but authors are Jung, Baek, Ma et al. 2024; repo cites 'Choi et al. (2024)'. |
| `10.1007/s13197-021-05084-7` | Arsa & Theerakulkait (2022), J Food Science and Technology | Arsa, Puechkamutr (2021) | Right topic (pyrazine yield vs pH, rice bran protein hydrolysate) but authors/year are Arsa & Puechkamutr 2021; repo cites 'Arsa & Theerakulkait (2022)'. |
| `10.1016/0163-7827(83)90002-4` | Frankel et al. (1982) | Frankel (1983) | Right topic (Frankel, volatile lipid oxidation products) but it is Frankel 1983 sole-author Prog Lipid Res; repo cites 'Frankel et al. (1982)'. |
| `10.1016/j.foodchem.2009.11.049` | Knol et al. (2009), Food Chemistry 113:103 | Knol, Linssen et al. (2010) | Right paper (Knol et al., fructose/asparagine acrylamide multiresponse kinetics) but published 2010 in Food Chem 120; repo cites 'Knol et al. (2009), Food Chemistry 113:103'. |
| `10.1016/j.foodchem.2022.134998` | Liu (2023) | Liu, Cadwallader et al. (2023) | Correct paper for the PPI off-note baseline (Liu, Cadwallader & Drake 2023, pea protein aroma), but this ONE DOI is reused for FOUR distinct repo ids (liu_2023_ppi_offnote_baseline, liu_cadwallader_2023_pea_aeda, liu_2022_ppi_oav_ |
| `10.1016/j.foodchem.2024.138795` | Shu et al. (2024), Food Chemistry | Kong, Wu et al. (2024) | Right topic (SPI volatile flavour under ultrasonic-thermal treatment) but authors are Kong, Wu, Li et al. 2024; repo cites 'Shu et al. (2024)'. |
| `10.1016/j.foodhyd.2024.110543` | Ince et al. (2024), Int J Biological Macromolecules | Ince, Condict et al. (2025) | Right topic (11S glycinin - hexanal binding) but the paper is Ince et al. 2025 in Food Hydrocolloids; repo cites 'Ince et al. (2024), Int J Biological Macromolecules' - wrong year and wrong journal. |
| `10.1016/j.foodhyd.2026.112497` | Sun et al. (2026) | Sun, Qiu et al. (2026) | Correct paper (Sun et al. 2026, VSC-soy protein binding) but slr_incorporation_matrix stores it under paper_id 'zhang_2026_spi_vsc_retention' while its own citation field says Sun et al. Internal id/citation conflict. |
| `10.1016/j.foodres.2016.12.009` | Fratianni et al. (2016) | Fratianni, Niro et al. (2017) | Right paper (Fratianni et al., carotenoid degradation and furosine in dried apricots) but published 2017; repo cites 'Fratianni et al. (2016)'. |
| `10.1016/j.foodres.2022.112187` | Liu et al. (2022) | Liu, Poojary et al. (2023) | Right paper (Liu et al., 4-methylbenzoquinone/lysine kinetics) but published 2023; repo cites 'Liu et al. (2022)'. |
| `10.1016/j.lwt.2010.02.016` | Biondi et al. (2010) | Chiavaro, Rodriguez-Estrada et al. (2010) | Topic matches (microwave heating of vegetable oils, thermal/chemical parameters) but authors are Chiavaro, Rodriguez-Estrada, Vittadini et al. 2010; repo cites 'Biondi et al. (2010)'. |
| `10.1016/j.lwt.2022.113651` | Cui et al. (2022) | Liu, Bau et al. (2022) | Topic broadly matches (shiitake drying: volatiles + taste properties, i.e. nucleotide/EUC context) but authors are Liu, Bau, Jin et al. 2022; repo cites 'Cui et al. (2022)' in both the registry and computational_priors. |
| `10.1016/j.tifs.2020.01.021` | Yu et al. (2017) | Zhu, Huang et al. (2020) | Topic matches (CML/CEL in thermally processed meat products review) but authors are Zhu, Huang, Cheng et al. 2020; repo cites 'Yu et al. (2017)'. |
| `10.1016/s0924-2244(01)00022-x` | Martins, Jongen & van Boekel (2001) | Martins, Jongen et al. (2000) | Correct paper (Martins, Jongen & van Boekel, TIFS Maillard kinetic-modelling review) but CrossRef dates it 2000; computational_priors.json and benchmark_intake_registry.json both cite it as (2001). The repo also stores this DOI in |
| `10.1021/acs.est.5b02902` | Shirai et al. (2015) | Kampf, Liu et al. (2015) | Topic partly matches (dityrosine cross-linking) but the paper is Kampf, Liu, Reinmuth-Selzle et al. 2015 on protein cross-linking upon OZONE exposure (Environ Sci Technol); repo cites 'Shirai et al. (2015)' and claims a BSA diffus |
| `10.1021/acs.jafc.0c04738` | Ding et al. (2020) | Troise, Fogliano et al. (2020) | Topic matches (Maillard/lipid-oxidation interplay in emulsions) but authors are Troise, Fogliano, Vitaglione et al. 2020; repo cites 'Ding et al. (2020)'. |
| `10.1021/acs.jafc.2c03427` | Bao et al. (2022) | Cui, Ma et al. (2022) | Topic matches (glycylglycine-catalysed Amadori degradation to deoxyosone) but authors are Cui, Ma, Wang et al. 2022; repo cites 'Bao et al. (2022)'. |
| `10.1021/acs.jafc.2c04919` | Ohsu et al. (2025) | Yang, Liao et al. (2022) | Topic matches (kokumi gamma-glutamyl peptides, CaSR/T1R1-T1R3) but the paper is Yang, Liao, Dong et al. 2022 JAFC; repo cites 'Ohsu et al. (2025)'. The stored EC50 constants are attributed to a paper that is not this DOI. |
| `10.1021/acs.jafc.3c08432` | Rizzello et al. (2024) | Flores, Comes et al. (2024) | Correct paper (Flores et al. 2024, fermented texturised pea protein aroma) for id acs_jafc_3c08432_crosstalk_cleanup_link, but the SAME DOI is also anchored to id 'rizzello_2024_fermentation_cleanup' citing 'Rizzello et al. (2024) |
| `10.1021/acs.jafc.4c05736` | Huang et al. (2024) | Zhang, Cui et al. (2024) | Topic matches (DiXyl-ARP degradation with extra xylose) but authors are Zhang, Cui, Xia et al. 2024; repo cites 'Huang et al. (2024)'. |
| `10.1021/acs.jafc.6b05811` | Li et al. (2017) | Liu, Xia et al. (2017) | Topic matches (quercetin/MGO adducts and alpha-dicarbonyls in lysine/glucose) but authors are Liu, Xia, Lu et al. 2017; repo cites 'Li et al. (2017)'. |
| `10.1021/acs.jafc.7b00264` | Monforte et al. (2018) | Monforte, Martins et al. (2017) | Right paper (Monforte et al., gallic acid/glucose/metals in phenylacetaldehyde formation) but published 2017; repo cites 'Monforte et al. (2018)'. |
| `10.1021/acs.jafc.9b07752` | J. Agric. Food Chem. 2019 (Ref. 24) | Arsad, Zainudin et al. (2020) | Topic matches (protein cysteine-phenol adducts in minced beef, i.e. polyphenol-thiol capping) but the paper is Arsad et al. 2020; repo cites it as 'J. Agric. Food Chem. 2019 (Ref. 24)'. |
| `10.1021/bi051902l` | Smagghe et al. (2006) | Smagghe, Sarath et al. (2005) | Right paper (Smagghe et al., ferrous hexacoordinate haemoglobin ligand-binding kinetics incl. leghaemoglobin) but published 2005; repo cites 'Smagghe et al. (2006)'. |
| `10.1021/ef101678s` | Spier et al. (2010) | West, Adams et al. (2011) | Topic matches (homogeneous catalysis of hydroperoxide decomposition) but authors/year are West, Adams & Zabarnick 2011 (Energy & Fuels); repo cites 'Spier et al. (2010)'. |
| `10.1021/jf0100797` | Adams et al. (2001) | Adams, Mottram et al. (2001) | Correct paper (Adams, Mottram, Parker et al. 2001, disulfide interchange between ovalbumin and volatile disulfides) and the citation field says Adams, but the repo paper_id is 'mottram_2001_bmfd_retention'. Also note the matrix is |
| `10.1021/jf025561j` | amadori | Davidek, Clety et al. (2002) | Davidek, Clety, Aubin et al. 2002 Amadori degradation. ALREADY FLAGGED in-file (comment claims Martins 2003). |
| `10.1021/jf051197n` | De Vleeschouwer et al. (2006) | Claeys, De Vleeschouwer et al. (2005) | DOI is Claeys, De Vleeschouwer & Hendrickx 2005 - correct for slr id claeys_2005_acrylamide_temperature_time, but the SAME DOI also anchors 'de_vleeschouwer_2006_acrylamide_aqueous' citing 'De Vleeschouwer et al. (2006)'. Duplicat |
| `10.1021/jf102026c` | Zamora et al. (2010) | Hidalgo, Delgado et al. (2010) | Right paper (asparagine decarboxylation by lipid oxidation products) but first author is Hidalgo (Hidalgo, Delgado, Navarro, Zamora 2010); repo cites 'Zamora et al. (2010)'. Same group - minor. |
| `10.1021/jf3032342` | Jian et al. (2012) | Hsu, Tsai et al. (2012) | Topic matches (degradation of ascorbic acid in ethanolic solutions) but authors are Hsu, Tsai, Fu et al. 2012; repo cites 'Jian et al. (2012)'. |
| `10.1021/jf9013394` | Richards et al. (2009) | Carvajal, Rustad et al. (2009) | Topic matches (haemoglobin-induced lipid oxidation in a liposome model) but authors are Carvajal, Rustad, Mozuraityte et al. 2009; repo cites 'Richards et al. (2009)'. |
| `10.1021/jf9708251` | Hendrickx et al. (1998) | Van den Broeck, Ludikhuyze et al. (1998) | Right topic (isobaric-isothermal L-ascorbic acid degradation kinetics) but first authors are Van den Broeck, Ludikhuyze, Weemaes et al. 1998; repo cites 'Hendrickx et al. (1998)' (Hendrickx is a later co-author) and claims a 'sque |
| `10.1021/jf970892v` | Hofmann & Schieberle (1997/1998), JAFC | Kerscher, Grosch (1998) | Topic matches exactly (quantification of MFT, FFT, mercaptopentanones in heated meat) but the paper is Kerscher & Grosch 1998; repo cites 'Hofmann & Schieberle (1997/1998)' under paper_id hofmann_schieberle_1997_meat_anchor. |
| `10.1021/jf990120u` | Hidalgo & Pompei (2000) | Hidalgo, Pompei (1999) | Right paper (Hidalgo & Pompei, HMF and furosine kinetics in tomato products) but published 1999; repo cites 'Hidalgo & Pompei (2000)'. |
| `10.1021/jf9907687` | Ramírez-Jiménez et al. (2000) | Ramírez-Jiménez, Guerra-Hernández et al. (2000) | Correct paper (Ramirez-Jimenez et al. 2000, 'Browning Indicators in Bread'); slr row correctly says toasted_bread, but the intake row records matrix_family='mild_legume_extrudate' for a bread study, and the benchmark file is named |
| `10.1093/chemse/bjaf043` | Tanaka et al. (2025) | Tanaka, Itoh et al. (2025) | Correct paper for the umami-threshold anchor (Tanaka, Itoh & Kondoh 2025) and the citation field says Tanaka, but the repo paper_id is 'mouritsen_2024_umami_thresholds'. |
| `10.1590/1678-457X.08717` | Yu et al. (2018) | YU, LI et al. (2017) | Right paper (Yu et al., browning kinetics in L-ascorbic acid/basic amino acid systems) but published 2017; repo cites 'Yu et al. (2018)'. Also stored twice in two DOI case variants. |
| `10.1590/1678-457x.08717` | Yu et al. (2018) | YU, LI et al. (2017) | Lowercase duplicate of the same Yu et al. 2017 paper; repo cites 2018. |
| `10.3389/fnut.2022.940202` | Zhu et al. (2022) | Ma, Long et al. (2022) | Topic matches (acrylamide/HMF in glucose-asparagine-linoleic acid kinetic model) but authors are Ma, Long, Li et al. 2022; repo cites 'Zhu et al. (2022)'. |
| `10.3390/antiox9080756` | Pereira et al. (2020) | García-Díez, Mora-Diez (2020) | Topic adjacent (theoretical study of iron-aminoguanidine complexes, secondary antioxidant activity) but authors are Garcia-Diez & Mora-Diez 2020; repo cites 'Pereira et al. (2020)' for a 'metal PM Haber-Weiss chelation' prior. |
| `10.3390/foods9030251` | Soladoye et al. (2020) | Hwang, Ismail et al. (2020) | Topic matches (sous-vide beef EUC / umami) but authors are Hwang, Ismail & Joo 2020; repo cites 'Soladoye et al. (2020)' in both the registry and computational_priors. |
| `10.3390/molecules25020373` | Zamora et al. (2020) | Goritschnig, Tadus et al. (2020) | Topic matches (carbonyl-amine adducts in PE-fortified soybean oil, radical scavenging) but authors are Goritschnig, Tadus, Koenig et al. 2020; repo cites 'Zamora et al. (2020)'. |
| `10.3390/molecules27196182` | Zhao et al. (2022) | Jo, An et al. (2022) | Topic matches (moromi fermentation metabolomics) but authors are Jo, An, Kim et al. 2022; repo cites 'Zhao et al. (2022)'. |
| `10.3390/molecules28073151` | Resconi et al. (2023) | Hernandez, Woerner et al. (2023) | Correct paper (Hernandez, Woerner, Brooks et al. 2023, sensory/volatiles of PBMAs vs ground beef) - the slr rows cite Hernandez correctly, but the intake row and the benchmark file both call it 'Resconi et al. (2023)'. |
| `10.4014/jmb.2207.07057` | Ahlberg & Mohammadi (2021) | Tao, Yuan et al. (2022) | Topic matches (yeast extract characteristics/production review) but authors are Tao, Yuan, Liu et al. 2022; repo cites 'Ahlberg & Mohammadi (2021)'. |

### Duplicate-anchor collapse

Cases where a single DOI is doing duty for several nominally-distinct references — the registry over-counts its
own literature coverage as a result:

| DOI | Distinct repo ids | Paper |
|---|---:|---|
| `10.1007/s10068-022-01194-w` | 4 | Cho, Park et al. (2022), "Study on volatile compounds formed from the thermal interaction of hydrolyzed vegetable protei |
| `10.1016/j.foodchem.2022.134998` | 4 | Liu, Cadwallader et al. (2023), "Identification of predominant aroma components of dried pea protein concentrates and is |
| `10.1021/acs.jafc.9b07711` | 4 | Bi, Xu et al. (2020), "Characterization of Key Aroma Compounds in Raw and Roasted Peas (<i>Pisum sativum</i> L.) by Appl |
| `10.3390/foods12101967` | 4 | Fu, Ma et al. (2023), "Contents and Correlations of Nε-(carboxymethyl)lysine, Nε-(carboxyethyl)lysine, Acrylamide and Nu |
| `10.1002/sfp2.1044` | 3 | Lincoln, Girard (2025), "Exogenous Polyphenols Inhibited Flavor Degradation in Legume Protein Matrices", *Sustainable Fo |
| `10.1016/j.foodchem.2019.126082` | 3 | Trikusuma, Paravisini et al. (2020), "Identification of aroma compounds in pea protein UHT beverages", *Food Chemistry* |
| `10.3390/ijms25168668` | 3 | Ma, Fu et al. (2024), "Impact of Extrusion Parameters on the Formation of Nε-(Carboxymethyl)lysine, Nε-(Carboxyethyl)lys |
| `10.3390/molecules26134104` | 3 | Singh, Shi et al. (2021), "A Rapid Gas-Chromatography/Mass-Spectrometry Technique for Determining Odour Activity Values  |
| `10.3390/molecules28073151` | 3 | Hernandez, Woerner et al. (2023), "Descriptive Sensory Attributes and Volatile Flavor Compounds of Plant-Based Meat Alte |

---

## 3. DEAD DOIs

45 DOIs return CrossRef HTTP 404. Repairs were applied **only** in
`benchmark_intake_registry.json` and `slr_incorporation_matrix.json`; each repaired site carries a `doi_repair`
record next to the corrected DOI. Benchmark JSONs and quarantined files were left untouched — see §4.

### 3a. Repaired (17 DOIs, 30 anchor sites)

| Old (dead) | New | Basis |
|---|---|---|
| `10.1007/s00217-003-0685-6` | `10.1007/s00217-002-0649-0` | CrossRef title search: exact match on authors, year, journal and subject - Krause, Knoll & Henle (2003), 'Studies on the formation of furosine and pyridosine during acid hydrolysis of different Amadori products of lysine', Eur Food Res Technol. |
| `10.1007/s00217-004-1065-x` | `10.1111/j.1365-2621.2005.tb11434.x` | CrossRef title search: Hidalgo, Nogales & Zamora (2005), 'Nonenzymatic Browning, Fluorescence Development, and Formation of Pyrrole Derivatives in Phosphatidylethanolamine/Ribose/Lysine Model Systems' - exact author/year/system match; journal is J Food Sci, not Eur Food Res Technol as the dead DOI implied. |
| `10.1007/s00726-005-0187-z` | `10.1007/s00726-005-0200-2` | CrossRef title search: Henle (2005), Amino Acids - exact author, year and journal match. Topic is protein-bound AGEs in foods, broader than the entry's 'soy glycinin CML' claim, so the payload content still needs human verification. |
| `10.1007/s11745-003-1082-x` | `10.1002/ejlt.200390032` | CrossRef title search: Kamal-Eldin, Velasco & Dobarganes (2003), 'Oxidation of mixtures of triolein and trilinolein at elevated temperatures' - exact author/year/substrate match; journal is Eur J Lipid Sci Technol, not Lipids. |
| `10.1016/0308-8146(91)90013-I` | `10.1016/0308-8146(91)90103-u` | CrossRef title search: Farmer & Patterson (1991), 'Compounds contributing to meat flavour', Food Chemistry - exact author/year/journal match (dead DOI differed only in the article suffix). |
| `10.1016/0891-5849(91)90081-A` | `10.1016/0891-5849(91)90192-6` | CrossRef title search: Esterbauer, Schaur & Zollner (1991), 'Chemistry and biochemistry of 4-hydroxynonenal, malonaldehyde and related aldehydes', Free Radic Biol Med - exact author/year/journal/topic match (dead DOI differed only in the article suffix). |
| `10.1016/S0260-8774(98)00085-4` | `10.1021/jf970605n` | CrossRef title search: Bhandari, D'Arcy & Bich (1998), 'Lemon Oil to beta-Cyclodextrin Ratio Effect on the Inclusion Efficiency of beta-Cyclodextrin and the Retention of Oil Volatiles in the Complex' - matches the cited author, year and the beta-cyclodextrin volatile-inclusion claim; journal is JAFC, not J Food Eng. |
| `10.1021/jf000900i` | `10.1007/s11745-001-0781-x` | CrossRef title search + author verification: Lin, Fay, Welti & BLANK (2001), 'Quantification of key odorants formed by autoxidation of arachidonic acid using isotope dilution assay', Lipids. Blank confirmed as co-author; trans-4,5-epoxy-2-decenal is a principal odorant quantified there, matching the 'epoxydecenal guardrail' claim. |
| `10.1021/jf015902m` | `10.1021/jf011164h` | CrossRef title search: Brands & van Boekel (2002), 'Kinetic Modeling of Reactions in Heated Monosaccharide-Casein Systems', JAFC - exact author/year/journal match and the source of the methylglyoxal/HDMF routing this entry describes. (Distinct from the Brands/Wedzicha/van Boekel melanoidin paper already in the registry.) |
| `10.1021/jf020082z` | `10.1016/s0308-8146(02)00155-3` | CrossRef title search: Rawel, Rohn & Kruse (2002), 'Structural changes induced in bovine serum albumin by covalent attachment of chlorogenic acid', Food Chemistry - matches the cited author/year and the CGA lysine/cysteine blocking claim; protein is BSA, not soy, so the transfer needs human sign-off. |
| `10.1021/jf020582d` | `10.1021/jf035265m` | CrossRef title search: 'alpha-Mercaptoketone Formation during the Maillard Reaction of Cysteine and [1-13C]Ribose' - an exact match to the described content (ribose-to-mercaptoketone branching from 13C labelling), and this DOI is already cited elsewhere in the repo for exactly that result. NOTE the true authors are Cerny & Davidek (2004), NOT 'Blank & Mottram (2002)': the citation string still needs correcting. |
| `10.1021/jf025830v` | `10.1016/s0008-6215(03)00174-5` | CrossRef title search: Martins & Van Boekel (2003), 'Kinetic modelling of Amadori N-(1-deoxy-d-fructos-1-yl)-glycine degradation pathways. Part II - Kinetic analysis', Carbohydrate Research - exact author/year/topic match for the 'Amadori degradation prior'. |
| `10.1021/jf0401831` | `10.1021/jf050085t` | CrossRef title search: Aliani & Farmer (2005), 'Precursors of Chicken Flavor. I. Determination of Some Flavor Precursors in Chicken Muscle', JAFC - exact author/year/journal match and the source of the ribose/nucleotide donor-potency context this entry describes. |
| `10.1021/jf0494452` | `10.1021/jf040334y` | CrossRef title search: 'Release of Nucleotides and Nucleosides during Yeast Autolysis: Kinetics and Potential Impact on Flavor', JAFC - an exact match to the described content. NOTE the true authors/year are Charpentier et al. (2005), NOT 'Perdiguero et al. (2004)': the citation string still needs correcting. |
| `10.1074/jbc.266.18.11644` | `10.1016/s0021-9258(18)99006-x` | CrossRef title search: Grandhee & Monnier (1991), 'Mechanism of formation of the Maillard protein cross-link pentosidine. Glucose, fructose, and ascorbate as pentosidine precursors', JBC - exact author/year/journal/topic match (the dead form used the legacy volume-page DOI pattern). |
| `10.1080/17415993.2012.715206` | `10.1080/17415993.2012.708933` | CrossRef title search: Tang, Jiang, Yuan & Ho, 'Flavor chemistry of 2-methyl-3-furanthiol, an intense meaty aroma compound', Journal of Sulfur Chemistry - same journal and same DOI prefix as the dead anchor, exact author and topic match (online 2012, issue 2013 as the repo cites). |
| `10.1111/j.1745-4530.1993.tb00179.x` | `10.1111/j.1745-4549.1993.tb00730.x` | CrossRef title search: Zhang, Lee & Ho (1993), 'Kinetics and Mechanism of Nonenzymatic Deamidation of Soy Protein' - exact author/year/topic match; the correct journal is J Food Processing and Preservation (ISSN 1745-4549), not J Food Process Engineering (1745-4530) as the dead DOI implied. |

> Two repairs fix the DOI but leave a **wrong citation string** in place, flagged in their `doi_repair.basis`:
> `10.1021/jf020582d` → the real authors are Cerny & Davidek (2004), not "Blank & Mottram (2002)"; and
> `10.1021/jf0494452` → the real authors are Charpentier et al. (2005), not "Perdiguero et al. (2004)".

### 3b. Unrepaired (28 DOIs)

No confidently identifiable replacement was found. Several are not DOIs at all but fabricated-looking strings
(`10.1021/acs.jafc.de_leyn_2019`, `10.1021/acs.jafc.liardon_1991`, `10.1016/j.foodchem.2015.00000`,
`10.1016/j.foodres.2025.001279`) — a confabulation signature rather than a transcription error.

| Dead DOI | Used in repo as | Claimed source |
|---|---|---|
| `10.1007/s11746-001-0402-y` | `tan_2001_oil_dsc_oxidation` | Tan et al. (2001) |
| `10.1016/0022-1416(84)90226-9` | `thiamine_degradation` | — |
| `10.1016/S0308-8146(98)00174-8` | `cys_glucose_150C_Farmer1999`, `$.source_metadata.quarantine_reason` | — |
| `10.1016/j.foodchem.2015.00000` | `fadel_2015_mft_retention` | Fadel et al. (2015) |
| `10.1016/j.foodchem.2022.132009` | `zhang_2022_unsaturated_aldehyde_potency_anchor` | Zhang et al. (2022) |
| `10.1016/j.foodchem.2025.142222` | `pmc11889959_spi_tvp_volatiles` | PMC11889959 (2025) |
| `10.1016/j.foodres.2021.110345` | `comunian_2021_thiamine_encapsulation` | Comunian et al. (2021) |
| `10.1016/j.foodres.2025.001279` | `hao_2025_sph_pyrazines` | Hao et al. (2025), Food Research International |
| `10.1016/j.lwt.2022.113009` | `wang_2022_lab_hexanal_cleanup_anchor` | Wang et al. (2022) |
| `10.1021/acs.jafc.2c08283` | `pmid36878579_pe_stoichiometry` | Fujimoto et al. (2023) |
| `10.1021/acs.jafc.5b02009` | `pmc9351765_crosspy_trapping_anchor` | Weerawatanakorn et al. (2015) |
| `10.1021/acs.jafc.9b05898` | `liu_2020_egcg_arp_kinetics` | Liu et al. (2023/2024) |
| `10.1021/acs.jafc.de_leyn_2019` | `de_leyn_2019` | De Leyn & Muylle (2019) |
| `10.1021/acs.jafc.liardon_1991` | `liardon_1991_r5p_donor_potency` | Liardon et al. (1991) |
| `10.1021/bk-1999-0740.ch012` | `hauck_tressl_1999_hdmf_non_amino` | Hauck & Tressl (1999) |
| `10.1021/jf001302l` | `hofmann_2001_melanoidin_thioether` | Hofmann et al. (2001) |
| `10.1021/jf010391+` | `mottram_2001_mft_quench_buffering_anchor` | Mottram et al. (2001) |
| `10.1021/jf021111m` | `martins_2003_lys_glucose_kinetics` | Martins & van Boekel (2003) |
| `10.1021/jf035223w` | `hidalgo_zamora_2004_4hne_pentylpyrrole` | Hidalgo & Zamora (2004) |
| `10.1021/jf053246m` | `davidek_2006_thr_glucose_furan` | Davidek et al. (2006) |
| `10.1021/jf3032779` | `acrylamide_asparagine_glucose_Parker2012` | — |
| `10.1021/jf502636m` | `marquez_ruiz_2014_oleic_oav_anchor` | Marquez-Ruiz et al. (2014) |
| `10.1021/jf960062o` | `hofmann_schieberle_grosch_1996`, `$.eligible_references[12].quarantine.quarantine_reason`, `$.entries[167].notes_on_limits`, `$.entries[167].quarantine.quarantine_reason` … | Hofmann et al. (1996) |
| `10.1021/jf960777q` | `blank_1997_rhamnose_proline_hdmf_uplift_v1`, `blank_1997_rhamnose_proline_hdmf_anchor` | Blank et al. (1997) |
| `10.1021/jf990111u` | `grosch_1999_c183_propanal_scission` | Grosch & Wieser (1999) |
| `10.1039/c9fo00878e` | `finnigan_2019_mycoprotein_rna` | Finnigan et al. (2019) |
| `10.1186/s12934-025-02688-w` | `xu_2025_peptide_hierarchy` | Xu, H. et al. (2025), PMC11743841 |
| `10.3989/gya.1051211` | `bayram_2023_stripped_soybean_oil_hexanal` | Bayram et al. (2023) |

Low-confidence candidates found but deliberately **not** applied:

- `10.1021/acs.jafc.5b02009` — Possible: Weerawatanakorn, Wu & Pan (2015), 'Reactivity and stability of selected flavor compounds', J Food Drug Anal, 10.1016/j.jfda.2015.02.001 - author/year match but the topic is only adjacent to the claimed melanoidin carbonyl-trapping prior. NOT applied.
- `10.1021/jf053246m` — Nearest: Davidek et al. (2006) 'Sugar Fragmentation in the Maillard Reaction Cascade' Parts I/II (10.1021/jf060667q, 10.1021/jf060668i) - author/year/journal match but neither concerns furan from threonine. NOT applied.
- `10.3989/gya.1051211` — Possible: Kavran, Yucel & Bakkalbasi (2024), hexanal/peroxide formation kinetics in ascorbyl-palmitate sunflower oil, 10.3989/gya.0320231 - same journal and endpoint (hexanal Arrhenius) but different authors and oil. NOT applied.

---

## 4. Proposed repairs requiring human sign-off

Benchmark JSONs and quarantined files were **not** modified. These are the changes recommended there:

### `data/benchmarks/thiamine_cys_xylose_145C_Cerny2008.json`

- **Field:** source_doi
- **Current:** `10.1021/jf800268t`
- **Proposed:** `10.1021/jf801762c`
- **Basis:** Current DOI is Limacher/Kerler/Davidek 2008 'Formation of Furan and Methylfuran'. The benchmark is named Cerny2008 and concerns thiamine+cysteine+xylose; the matching paper is Cerny & Guntz-Dubini (2008), 'Identification of 5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine, Cysteine, and Xylose', which is ALREADY registered in benchmark_intake_registry.json under cerny_guntz_dubini_2008. Values unverified against either paper.
- **Status:** NOT APPLIED - benchmark JSONs require human sign-off

### `data/benchmarks/resconi_2023_pbma_beef_identity_benchmark.json`

- **Field:** citation/author label (DOI is correct)
- **Current:** `Resconi et al. (2023)`
- **Proposed:** `Hernandez, Woerner, Brooks et al. (2023)`
- **Basis:** DOI 10.3390/molecules28073151 resolves correctly to the intended study but its authors are Hernandez et al., not Resconi. The slr_incorporation_matrix rows already cite Hernandez; the benchmark filename and the intake row do not.
- **Status:** NOT APPLIED - benchmark JSONs require human sign-off

### `data/benchmarks/acrylamide_asparagine_glucose_Parker2012.json`

- **Field:** source_doi
- **Current:** `10.1021/jf3032779`
- **Proposed:** _none identified_
- **Basis:** DOI is dead. No Parker 2012 acrylamide asparagine/glucose kinetics paper could be located via CrossRef bibliographic search. Nearest same-year JAFC hit is Halford/Muttucumaru/Powers (2012) 10.1021/jf3037566 on free amino acids and sugars in potato varieties, which is NOT a kinetic acrylamide model - not proposed as a replacement.
- **Status:** NOT APPLIED - unidentified; recommend quarantine review

### `data/benchmarks/quarantined/cys_ribose_150C_Mottram1994.json`

- **Field:** source_doi
- **Current:** `10.1021/jf00049a015`
- **Proposed:** `10.1021/jf00052a027 (candidate only)`
- **Basis:** Confirms the existing in-file note: the DOI resolves to Zarkadas et al. (1995) beefstock bone protein quality. Mottram & Whitfield (1995) 10.1021/jf00052a027 verified live and on-topic, but the reported values remain unverified.
- **Status:** NOT APPLIED - already quarantined; needs human sign-off

### `data/benchmarks/quarantined/cys_glucose_150C_Farmer1999.json`

- **Field:** source_doi
- **Current:** `10.1016/S0308-8146(98)00174-8`
- **Proposed:** _none identified_
- **Basis:** Dead DOI confirmed. No Farmer 1999 cysteine/glucose volatile source located. Nearest real Farmer works are on cysteine/RIBOSE/phospholipid (10.1002/jsfa.2740490311, 10.1002/jsfa.2740530409), a different system.
- **Status:** NOT APPLIED - already quarantined

### `data/benchmarks/quarantined/thiamine_cys_ribose_100C_Hofmann1996.json`

- **Field:** source_doi
- **Current:** `10.1021/jf960062o`
- **Proposed:** _none identified_
- **Basis:** Dead DOI confirmed. Nearest real Hofmann & Schieberle 1996 item found is the Flavour Science chapter 10.1533/9781845698232.4.182 on MFT/2-acetyl-2-thiazoline/sotolon intermediates - not a quantitative thiamine/cys/ribose benchmark. The pi-valued MFT finding stands.
- **Status:** NOT APPLIED - already quarantined

---

## 5. Audit flags added by this pass

Values UNCHANGED everywhere. One-line `audit_flag` fields (dated 2026-08-26, matching the existing in-file style)
were added to entries that did not already carry one:

**`data/lit/arrhenius_params.yml`** (3 new; 5 already flagged by the prior pass)

| Entry | Class | Why |
|---|---|---|
| `beta_elimination_dha` | TOPIC-MISMATCH | Helmerhorst & Stokes 1983 is a mechanistic insulin-persulfide characterisation, not a kinetics study; the 19 kcal/mol midpoint is untraceable to it. |
| `thiamine_degradation` | DEAD | `10.1016/0022-1416(84)90226-9` 404s and the source string itself hedges "or similar lit". Real thiamine kinetics anchors (Voelker 2018/2021) already exist in the intake registry. |
| `strecker` | non-citable | Source is the bare string "T.P. Labuza 1998" — no DOI, title or venue; the `A_value` comment separately credits an uncited "Hofmann 2000". |

**`data/lit/computational_priors.json`** (6 new; none previously flagged)

| Prior id | Class |
|---|---|
| `blank_devaud_grosch_2003_g6p_hdmf_uplift_v1` | `TOPIC-MISMATCH` — DOI is Cerny & Davidek 2003, 'Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction'. Repo uses it as 'blank_devaud_grosch_2003_g6p_hdmf' carrying a g6 |
| `ohsu_2025_kokumi_casr_support_v1` | `METADATA-MISMATCH` — Topic matches (kokumi gamma-glutamyl peptides, CaSR/T1R1-T1R3) but the paper is Yang, Liao, Dong et al. 2022 JAFC; repo cites 'Ohsu et al. (2025)'. The stored EC50 constants are at |
| `soladoye_2020_low_temp_euc_window_v1` | `METADATA-MISMATCH` — Topic matches (sous-vide beef EUC / umami) but authors are Hwang, Ismail & Joo 2020; repo cites 'Soladoye et al. (2020)' in both the registry and computational_priors. |
| `cui_2022_mushroom_gmp_euc_window_v1` | `METADATA-MISMATCH` — Topic broadly matches (shiitake drying: volatiles + taste properties, i.e. nucleotide/EUC context) but authors are Liu, Bau, Jin et al. 2022; repo cites 'Cui et al. (2022)' in both |
| `martins_2001_maillard_kinetics_modelling_v1` | `METADATA-MISMATCH` — Correct paper (Martins, Jongen & van Boekel, TIFS Maillard kinetic-modelling review) but CrossRef dates it 2000; computational_priors.json and benchmark_intake_registry.json both c |
| `blank_1997_rhamnose_proline_hdmf_uplift_v1` | `DEAD` — DOI does not resolve in CrossRef (HTTP 404). |

---

## 6. Full ledger — every anchor

| DOI | Class | Files | Repo ids | Finding |
|---|---|---|---|---|
| `10.1002/recl.19871060202` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `de_bruijn_1987_monosaccharide_alkaline_degradation` | DOI is Chittenden & Regeling 1987, 'An improved procedure for the synthesis of sugar chloroacetates' (Recl Trav Chim Pays-Bas) - a preparative carbohydrate-synthesis paper. Repo uses it as 'de Bruijn 1987 monosaccharide alkaline degradation / alkaline retro-al |
| `10.1006/fstl.1996.0092` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `ilo_1996_maize_sme_lysine_damage` | DOI is Ilo et al. 1996 on apparent viscosity and extrudate properties in twin-screw extrusion of maize grits. Repo id claims 'maize_sme_lysine_damage' and 'shear damage priors'. Extrusion/SME context matches; there is no lysine-damage measurement in this paper |
| `10.1016/S0260-8774(00)00047-9` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `quintas_2000_sucrose_caramelisation` | DOI is Ramaswamy & Zareifard 2000, 'Evaluation of factors influencing tube-flow fluid-to-particle heat transfer coefficient' (J Food Eng). Repo uses it as 'quintas_2000_sucrose_caramelisation' / 'sucrose caramelisation kinetics priors'. No caramelisation conte |
| `10.1016/j.carbpol.2022.120468` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `tran_2023_reducing_sugar_hme` | DOI is Liu et al. 2023, 'Self-healing, antibacterial, and conductive double network hydrogel for strain sensors' (Carbohydr Polym). Repo uses it as 'tran_2023_reducing_sugar_hme' / 'starch reducing-sugar liberation during HME'. Completely unrelated field. |
| `10.1016/j.foodchem.2011.07.037` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `wang_2012_gsh_xylose_sulfur_uplift` | DOI is Koblar et al. 2012, 'Fluoride in teas of different types and forms and the exposure of humans to fluoride with tea and diet'. Repo uses it as 'wang_2012_gsh_xylose_sulfur_uplift' (glutathione/xylose sulfur-uplift prior). Unrelated. |
| `10.1016/j.foodchem.2013.11.082` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `hrncirik_2014_coconut_oil_thermal_profile` | DOI is Mashhadizadeh et al. 2014, solid-phase extraction of trace Ag/Cd/Cu/Hg/Pb using modified Fe3O4 nanoparticles. Repo uses it as 'hrncirik_2014_coconut_oil_thermal_profile'. Unrelated (analytical trace-metal method). |
| `10.1016/j.foodchem.2014.05.097` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `ordoudi_2014_hmf_peak_window` | DOI is Agcam et al. 2014, PEF/heat pasteurisation effects on pectin methylesterase in orange juice. Repo uses it as 'ordoudi_2014_hmf_peak_window' (caramelisation HMF peak-window calibration). Unrelated. |
| `10.1016/j.foodchem.2020.128670` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `comert_gokmen_2021_epicatechin_cysteine_mgo_synergy` | DOI is Azevedo et al. 2021, 'Rhamnolipids-based nanostructured lipid carriers'. Repo uses it as 'comert_gokmen_2021_epicatechin_cysteine_mgo_synergy'. Unrelated (colloid formulation). |
| `10.1016/j.foodchem.2023.136009` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `yeo_mottram_2023_lecithin_crosstalk_anchor` | DOI is Bae et al. 2023, 'Ethanol organosolv lignin as a substitute for commercial antioxidants'. Repo uses it as 'yeo_mottram_2023_lecithin_crosstalk_anchor' (soy lecithin phospholipid crosstalk). Unrelated. |
| `10.1016/j.foodchem.2023.136200` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `poulsen_2023_pbma_cml_cel` | DOI is Weligama Thuppahige et al. 2023 on cassava peel/bagasse for biodegradable food packaging. Repo uses it as 'poulsen_2023_pbma_cml_cel' (commercial PBMA CML/CEL safety anchor). Unrelated. |
| `10.1016/j.foodchem.2025.145815` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `siripitakpong_2026_fft_retention` | DOI is Sharma et al. 2025, a REVIEW on plant-protein isolates and flavour perception. Repo uses it as 'siripitakpong_2026_fft_retention' with notes claiming values 'obtained via molecular docking and homology modeling'. A review contains no such primary dockin |
| `10.1016/j.foodhyd.2024.110509` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `malia_2025_pea_free_sh_crosscheck` | DOI is Gong et al. 2025, 'Effects of esterification and enzymatic modification on the properties of wheat starch and dough'. Repo uses it as 'malia_2025_pea_free_sh_crosscheck' (pea free-SH Ellman crosscheck). Unrelated. |
| `10.1016/j.jfoodeng.2014.06.032` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `luna_aguilera_2014_molten_sugar_color_kinetics` | DOI is Gelaw et al. 2014, ATR-IR microspectroscopy of membrane fouling. Repo uses it as 'luna_aguilera_2014_molten_sugar_color_kinetics' (caramelisation colour kinetics). Unrelated. |
| `10.1016/j.lwt.2006.07.014` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mundt_wedzicha_2007_biscuit_browning` | DOI is Kavitha & Modi 2007 on 5'-IMP degradation vs water activity in a meat model. Repo uses it as 'mundt_wedzicha_2007_biscuit_browning' (biscuit surface browning rate prior). Unrelated. |
| `10.1016/j.lwt.2021.111352` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `huang_2021_sulfur_oav_support` | DOI is Chraibi et al. 2021, antimicrobial mixture-design study of Mentha/Ormenis essential oils vs S. aureus, E. coli, C. tropicalis. Repo uses it as 'huang_2021_sulfur_oav_support' (MFT/FFT sulfur prominence prior). Unrelated. |
| `10.1016/j.lwt.2021.111802` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `yu_2021_corn_hydrolysate_kinetics` | DOI is Patel, Andrade & Rousseau 2021, fat-crystal-stabilised W/O emulsion breakdown during in vitro digestion. Repo uses it as 'yu_2021_corn_hydrolysate_kinetics', kind=quantitative_benchmark. Unrelated. |
| `10.1016/j.lwt.2024.117316` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pruteanu_2023_glucose_phe_arrhenius_browning` | DOI is Bolchini et al. 2025, free-radical scavenging KINETICS of Maillard reaction products in a glucose-glycine system. Repo uses it as 'pruteanu_2023_glucose_phe_arrhenius_browning' (glucose/phenylalanine browning activation energy). Wrong authors, wrong yea |
| `10.1016/j.meatsci.2020.108151` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `sun_2020_solid_matrix_cml_cel_accumulation` | DOI is Jo et al. 2020 on plasma-treated winter mushroom as a nitrite/phosphate replacer in ground ham. Repo uses it as 'sun_2020_solid_matrix_cml_cel_accumulation'. Unrelated. |
| `10.1016/j.mineng.2010.05.003` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `li_2010_phytate_chelation_kinetics` | DOI is Smith, Hadler & Cilliers 2010, 'Flotation bank air addition and distribution for optimal performance' in MINERALS ENGINEERING. Repo uses it as 'li_2010_phytate_chelation_kinetics'. Wrong discipline entirely. |
| `10.1021/acs.jafc.9b06882` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `ref41_ppi_sulfur_binding` | DOI is Zhao et al. 2019, 'Lycopene Prevents DEHP-Induced Leydig Cell Damage with the Nrf2 Antioxidant Signaling Pathway in Mice' - a rodent reproductive-toxicology paper. Repo uses it as 'ref41_ppi_sulfur_binding' (pea-protein sulfur-volatile binding order pri |
| `10.1021/acs.jafc.9b07711` | `TOPIC-MISMATCH` | benchmarks/external_validation/external_validation_bi_2020_raw_pea_hexanal.json, benchmarks/external_validation/external_validation_bi_2020_roasted_pea_hexanal.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `acs_2020_raw_pea_hexanal_baseline`, `jafc_2020_egcg_deoxyosone_trapping`, `external_validation_bi_2020_raw_pea_hexanal`, `external_validation_bi_2020_roasted_pea_hexanal` | DOI is Bi et al. 2020, key aroma compounds in raw and roasted peas - CORRECT for ids acs_2020_raw_pea_hexanal_baseline and both external_validation benchmarks. But the same DOI is ALSO anchored to 'jafc_2020_egcg_deoxyosone_trapping' (EGCG trapping of deoxyoso |
| `10.1021/acsomega.7b00321` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `nashalian_yaylayan_2014_cu_catalyzed_strecker` | DOI is Jefferson et al. 2017 ACS Omega, aqueous metal-enol(ate) complexes of (acetonedicarboxylato)copper. Repo uses it as 'nashalian_yaylayan_2014_cu_catalyzed_strecker'. Cu-enolate chemistry is adjacent but this is not the Nashalian & Yaylayan Strecker study |
| `10.1021/bi00270a010` | `TOPIC-MISMATCH` | lit/arrhenius_params.yml | `beta_elimination_dha` | DOI is Helmerhorst & Stokes 1983 Biochemistry, characterising an acid-stable persulfide-like residue in alkali/sulfhydryl-treated INSULIN. It is mechanistic evidence for beta-elimination of disulfides, not a kinetics paper: the Ea=79.5 kJ/mol (19 kcal/mol) use |
| `10.1021/jf00002a024` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `yeo_shibamoto_1991_wof_hexanal` | DOI is Ramarathnam, Rubin & Diosady 1991, 'Studies on meat flavor. 1. ... uncured and cured pork'. Repo uses it as 'yeo_shibamoto_1991_wof_hexanal' claiming phospholipid hexanal contribution AND lipoxygenase denaturation kinetics. Neither is in this paper. |
| `10.1021/jf00009a012` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `blank_grosch_1991_hdmf_anchor` | DOI is Mattheis et al. 1991, 'Change in apple fruit volatiles after storage in atmospheres inducing anaerobic metabolism'. Repo uses it as 'blank_grosch_1991_hdmf_anchor' (HDMF in boiled beef, kind=quantitative_benchmark). Unrelated. |
| `10.1021/jf00049a015` | `TOPIC-MISMATCH` | benchmarks/quarantined/cys_ribose_150C_Mottram1994.json | `cys_ribose_150C_Mottram1994`, `$.source_metadata.quarantine_reason` | DOI is Zarkadas et al. 1995 beefstock bone protein quality. Repo uses it for the quarantined 'cys_ribose_150C_Mottram1994' benchmark. ALREADY documented in-file by the previous audit; retained here for completeness. |
| `10.1021/jf00058a011` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `glomb_1995_3dg_fragmentation_stoichiometry` | DOI is Cabras et al. 1995, 'Fate of Some Insecticides from Vine to Wine'. Repo uses it as 'glomb_1995_3dg_fragmentation_stoichiometry' (3-deoxyglucosone fragmentation prior). Unrelated. |
| `10.1021/jf00065a009` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `tannenbaum_1985_thiamine_ea` | DOI is Al-Wandawi et al. 1985, 'Tomato processing wastes as essential raw materials source'. Repo uses it as 'tannenbaum_1985_thiamine_ea' (thiamine activation-energy prior). Unrelated. |
| `10.1021/jf00080a032` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `arabshahi_1988_aw_thiamine_kinetics` | DOI is Schaefer & Sandermann 1988, metabolism of PENTACHLOROPHENOL in wheat cell suspension cultures. Repo uses it as 'arabshahi_1988_aw_thiamine_kinetics' (water-activity-dependent thiamine kinetics). Unrelated. |
| `10.1021/jf00083a026` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json | `nakamura_1988_imp_ribose_release` | DOI is Epstein et al. 1988 on sulfamethazine/chloramphenicol/cyromazine drug residues in muscle tissue. Repo uses it as 'nakamura_1988_imp_ribose_release' (IMP-to-ribose delivery anchor). Unrelated. |
| `10.1021/jf00111a008` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `grosch_1982_hexanal_odt` | DOI is Purcell & Walter 1982, 'Stability of amino acids during cooking and processing of sweet potatoes'. Repo uses it as 'grosch_1982_hexanal_odt' - the hexanal ODOUR THRESHOLD anchor used for off-note scoring. Unrelated; the ODT value is unsourced. |
| `10.1021/jf026123f` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `blank_devaud_grosch_2003_g6p_hdmf_uplift_v1`, `blank_devaud_grosch_2003_g6p_hdmf_prior` | DOI is Cerny & Davidek 2003, 'Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction'. Repo uses it as 'blank_devaud_grosch_2003_g6p_hdmf' carrying a g6p_hdmf_multiplier_vs_glucose = 4.4 with 'SIDA with [2H3]-HDMF'. This DOI contain |
| `10.1021/jf047903m` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `martins_van_boekel_2005_ascorbic_amino_browning` | DOI is Adams et al. 2005, 'Thermal Degradation Studies of Food Melanoidins'. Repo uses it as 'martins_van_boekel_2005_ascorbic_amino_browning' / 'ascorbic-basic-amino browning prior'. Melanoidin family is adjacent but there is no ascorbic-acid content and the  |
| `10.1021/jf800268t` | `TOPIC-MISMATCH` | benchmarks/thiamine_cys_xylose_145C_Cerny2008.json | `thiamine_cys_xylose_145C_Cerny2008` | DOI is Limacher, Kerler, Davidek et al. 2008, 'Formation of Furan and Methylfuran by Maillard-Type Reactions'. It is used as source_doi of benchmark thiamine_cys_xylose_145C_Cerny2008.json. The actual Cerny 2008 thiamine/cysteine/xylose paper is 10.1021/jf8017 |
| `10.1080/10408398.2017.1378865` | `TOPIC-MISMATCH` | lit/arrhenius_params.yml | `retro_aldol` | ALjahdali & Carbonero 2017 MRP nutrition/health review. ALREADY FLAGGED in-file. |
| `10.1111/j.1365-2621.1987.tb14251.x` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `buera_1987_maillard_caramelisation_ea` | DOI is Koch, Seib & Hoseney 1987, 'Incorporating L-Ascorbyl 6-Palmitate in Bread and Its Shortening-Sparing and Anti-Firming Effects'. Repo uses it as 'buera_1987_maillard_caramelisation_ea' (browning/caramelisation activation-energy priors). Unrelated. |
| `10.3390/foods10123140` | `TOPIC-MISMATCH` | lit/slr_incorporation_matrix.json | `karolkowski_2021_ppi_ph_release` | DOI is Karolkowski et al. 2021, 'Volatile Compounds in Pulses: A REVIEW' (Foods). Repo uses it as 'karolkowski_2021_ppi_ph_release' - a PPI native-release pH SERIES. A review supplies no such pH-resolved dataset. |
| `10.3390/foods12061349` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc10056349_rubisco_amadori` | DOI is Simkova et al. 2023, 'Variability in Capri Everbearing STRAWBERRY Quality during a Harvest Season'. Repo uses it as 'pmc10056349_rubisco_amadori' (Rubisco amino acid composition and Amadori kinetics). Unrelated. |
| `10.3390/foods13091393` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `xu_2024_soybean_pbma_hexanal` | DOI is Spizzirri & Restuccia 2024, 'Advances in Food Waste Biomass Transformation into High-Value Products' (an editorial/review). Repo uses it as 'xu_2024_soybean_pbma_hexanal' (soy PBMA storage hexanal reference). Unrelated. |
| `10.3390/foods13162590` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc11353891_lentil_deflavoring` | DOI is Streule et al. 2024 on depulping/pod storage/fermentation of Ghanaian COCOA. Repo uses it as 'pmc11353891_lentil_deflavoring' (lentil isolate off-note deflavouring). Unrelated. |
| `10.3390/foods14111881` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc12154226_crosspy_thiol_adduction` | DOI is El Hosry et al. 2025, a general 'Maillard Reaction: Mechanism, Influencing Parameters...' REVIEW. Repo uses it as 'pmc12154226_crosspy_thiol_adduction' - the Family 13 covalent thiol-depletion anchor for CROSSPY-mediated capping in a coffee melanoidin m |
| `10.3390/foods14193453` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `wang_2024_glucosamine_synergy` | DOI is Wu et al. 2025, MRPs from Nile tilapia fish-skin collagen peptides with four REDUCING SUGARS. Repo uses it as 'wang_2024_glucosamine_synergy' claiming 'glucosamine model systems' and glucosamine synergy priors. Wrong authors, wrong year, and glucosamine |
| `10.3390/foods7080126` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc6104182_soybean_fermentation` | DOI is Chambers, Godwin & Terry 2018, 'Recipes for Determining Doneness in POULTRY Do Not Provide Appropriate Information Based on US Government Guidelines'. Repo uses it as 'pmc6104182_soybean_fermentation'. Unrelated. |
| `10.3390/plants13020274` | `TOPIC-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mdpi_plants_2024_hemp_volatiles` | DOI is Chen et al. 2024, 'RNA-Seq-Based WGCNA ... Salt Stress in WHEAT ROOTS' (Plants). Repo uses it as 'mdpi_plants_2024_hemp_volatiles' (hemp protein off-note flavour reference). Wrong field entirely. |
| `10.5194/egusphere-2026-321` | `TOPIC-MISMATCH` | lit/arrhenius_params.yml | `cysteine_thermolysis` | Wu et al. 2026 EGUsphere geoscience preprint (TGA/DSC of geologic materials). ALREADY FLAGGED in-file. |
| `10.1007/s00217-003-0685-6` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `krause_2003_furosine_hydrolysis_yields` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1007/s00217-004-1065-x` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hidalgo_2005_pe_ribose_lysine` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1007/s00726-005-0187-z` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `henle_2005_glycinin_cml` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1007/s11745-003-1082-x` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `kamal_eldin_2003_triolein_scission` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1007/s11746-001-0402-y` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `tan_2001_oil_dsc_oxidation` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/0022-1416(84)90226-9` | `DEAD` | lit/arrhenius_params.yml | `thiamine_degradation` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/0308-8146(91)90013-I` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `farmer_1991_alkyl_thiazoles` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/0891-5849(91)90081-A` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `esterbauer_1991_4hne_kinetics` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/S0260-8774(98)00085-4` | `DEAD` | lit/benchmark_intake_registry.json | `bhandari_1998_beta_cd_aldehyde_binding_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/S0308-8146(98)00174-8` | `DEAD` | benchmarks/quarantined/cys_glucose_150C_Farmer1999.json | `cys_glucose_150C_Farmer1999`, `$.source_metadata.quarantine_reason` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.foodchem.2015.00000` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `fadel_2015_mft_retention` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.foodchem.2022.132009` | `DEAD` | lit/benchmark_intake_registry.json | `zhang_2022_unsaturated_aldehyde_potency_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.foodchem.2025.142222` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc11889959_spi_tvp_volatiles` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.foodres.2021.110345` | `DEAD` | lit/benchmark_intake_registry.json | `comunian_2021_thiamine_encapsulation` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.foodres.2025.001279` | `DEAD` | lit/slr_incorporation_matrix.json | `hao_2025_sph_pyrazines` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1016/j.lwt.2022.113009` | `DEAD` | lit/benchmark_intake_registry.json | `wang_2022_lab_hexanal_cleanup_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/acs.jafc.2c08283` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmid36878579_pe_stoichiometry` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/acs.jafc.5b02009` | `DEAD` | lit/benchmark_intake_registry.json | `pmc9351765_crosspy_trapping_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/acs.jafc.9b05898` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `liu_2020_egcg_arp_kinetics` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/acs.jafc.de_leyn_2019` | `DEAD` | lit/benchmark_intake_registry.json | `de_leyn_2019` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/acs.jafc.liardon_1991` | `DEAD` | lit/benchmark_intake_registry.json | `liardon_1991_r5p_donor_potency` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/bk-1999-0740.ch012` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hauck_tressl_1999_hdmf_non_amino` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf000900i` | `DEAD` | lit/benchmark_intake_registry.json | `blank_2001_epoxydecenal_guardrail` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf001302l` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hofmann_2001_melanoidin_thioether` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf010391+` | `DEAD` | lit/benchmark_intake_registry.json | `mottram_2001_mft_quench_buffering_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf015902m` | `DEAD` | lit/benchmark_intake_registry.json | `brands_2002_mgo_hdmf_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf020082z` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `rawel_2002_cga_cysteine_blocking` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf020582d` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `blank_mottram_2002_ribose_labeling` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf021111m` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `martins_2003_lys_glucose_kinetics` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf025830v` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `martins_boekel_2003_dfg_amadori_degradation` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf035223w` | `DEAD` | lit/benchmark_intake_registry.json | `hidalgo_zamora_2004_4hne_pentylpyrrole` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf0401831` | `DEAD` | lit/benchmark_intake_registry.json | `aliani_2005_donor_potency_nucleotide_context` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf0494452` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `perdiguero_2004_yeast_autolysis_nucleotides` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf053246m` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `davidek_2006_thr_glucose_furan` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf3032779` | `DEAD` | benchmarks/acrylamide_asparagine_glucose_Parker2012.json | `acrylamide_asparagine_glucose_Parker2012` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf502636m` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `marquez_ruiz_2014_oleic_oav_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf960062o` | `DEAD` | benchmarks/quarantined/thiamine_cys_ribose_100C_Hofmann1996.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hofmann_schieberle_grosch_1996`, `$.eligible_references[12].quarantine.quarantine_reason`, `$.entries[167].notes_on_limits`, `$.entries[167].quarantine.quarantine_reason` … | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf960777q` | `DEAD` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `blank_1997_rhamnose_proline_hdmf_uplift_v1`, `blank_1997_rhamnose_proline_hdmf_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1021/jf990111u` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `grosch_1999_c183_propanal_scission` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1039/c9fo00878e` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `finnigan_2019_mycoprotein_rna` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1074/jbc.266.18.11644` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmid_1904866_pentosidine_equivalence_anchor` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1080/17415993.2012.715206` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `tang_2013_thiamine_mft` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1111/j.1745-4530.1993.tb00179.x` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zhang_1993_protein_deamidation_ammonia` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1186/s12934-025-02688-w` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `xu_2025_peptide_hierarchy` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.3989/gya.1051211` | `DEAD` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `bayram_2023_stripped_soybean_oil_hexanal` | DOI does not resolve in CrossRef (HTTP 404). |
| `10.1002/jsfa.10528` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zhu_2021_braised_chicken_ages` | Right paper (Zhu et al., CML/CEL kinetics in chicken braising) but published 2020; repo cites 'Zhu et al. (2021)'. |
| `10.1007/s10068-016-0185-5` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `chen_2016_carbohydrate_meat_cml_cel_sterilization` | Right topic (CML/CEL formation in meat products) but authors are Yu, Gao, Zeng et al. 2016; repo cites 'Chen et al. (2016)'. |
| `10.1007/s10068-024-01633-w` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `choi_2024_garden_pea_acrylamide` | Right topic (acrylamide in air-fryer roasted legumes vs asparagine/free sugars) but authors are Jung, Baek, Ma et al. 2024; repo cites 'Choi et al. (2024)'. |
| `10.1007/s13197-021-05084-7` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `arsa_theerakulkait_2022_pyrazine_ph` | Right topic (pyrazine yield vs pH, rice bran protein hydrolysate) but authors/year are Arsa & Puechkamutr 2021; repo cites 'Arsa & Theerakulkait (2022)'. |
| `10.1016/0163-7827(83)90002-4` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `frankel_1982_c182_hexanal_scission` | Right topic (Frankel, volatile lipid oxidation products) but it is Frankel 1983 sole-author Prog Lipid Res; repo cites 'Frankel et al. (1982)'. |
| `10.1016/j.foodchem.2009.11.049` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `knol_2009_acrylamide_kinetics` | Right paper (Knol et al., fructose/asparagine acrylamide multiresponse kinetics) but published 2010 in Food Chem 120; repo cites 'Knol et al. (2009), Food Chemistry 113:103'. |
| `10.1016/j.foodchem.2022.134998` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `liu_2023_ppi_offnote_baseline`, `liu_cadwallader_2023_pea_aeda`, `liu_2022_ppi_oav_anchors`, `researchgate_2023_pea_aeda` | Correct paper for the PPI off-note baseline (Liu, Cadwallader & Drake 2023, pea protein aroma), but this ONE DOI is reused for FOUR distinct repo ids (liu_2023_ppi_offnote_baseline, liu_cadwallader_2023_pea_aeda, liu_2022_ppi_oav_anchors, researchgate_2023_pea |
| `10.1016/j.foodchem.2024.138795` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `shu_2024_heated_spi` | Right topic (SPI volatile flavour under ultrasonic-thermal treatment) but authors are Kong, Wu, Li et al. 2024; repo cites 'Shu et al. (2024)'. |
| `10.1016/j.foodhyd.2024.110543` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `ince_2024_hexanal_glycinin` | Right topic (11S glycinin - hexanal binding) but the paper is Ince et al. 2025 in Food Hydrocolloids; repo cites 'Ince et al. (2024), Int J Biological Macromolecules' - wrong year and wrong journal. |
| `10.1016/j.foodhyd.2026.112497` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `sun_2026_spi_vsc_retention`, `zhang_2026_spi_vsc_retention` | Correct paper (Sun et al. 2026, VSC-soy protein binding) but slr_incorporation_matrix stores it under paper_id 'zhang_2026_spi_vsc_retention' while its own citation field says Sun et al. Internal id/citation conflict. |
| `10.1016/j.foodres.2016.12.009` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `fratianni_2016_apricot_furosine` | Right paper (Fratianni et al., carotenoid degradation and furosine in dried apricots) but published 2017; repo cites 'Fratianni et al. (2016)'. |
| `10.1016/j.foodres.2022.112187` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `liu_2023_benzoquinone_lysine_kinetics` | Right paper (Liu et al., 4-methylbenzoquinone/lysine kinetics) but published 2023; repo cites 'Liu et al. (2022)'. |
| `10.1016/j.lwt.2010.02.016` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `biondi_2010_oil_microwave_degradation` | Topic matches (microwave heating of vegetable oils, thermal/chemical parameters) but authors are Chiavaro, Rodriguez-Estrada, Vittadini et al. 2010; repo cites 'Biondi et al. (2010)'. |
| `10.1016/j.lwt.2022.113651` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `cui_2022_mushroom_gmp_euc_window_v1`, `cui_2022_mushroom_nucleotide_anchor` | Topic broadly matches (shiitake drying: volatiles + taste properties, i.e. nucleotide/EUC context) but authors are Liu, Bau, Jin et al. 2022; repo cites 'Cui et al. (2022)' in both the registry and computational_priors. |
| `10.1016/j.tifs.2020.01.021` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `yu_2017_cml_cel_meat_review` | Topic matches (CML/CEL in thermally processed meat products review) but authors are Zhu, Huang, Cheng et al. 2020; repo cites 'Yu et al. (2017)'. |
| `10.1016/s0924-2244(01)00022-x` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `martins_2001_maillard_kinetics_modelling_v1`, `martins_2001_maillard_kinetics_modelling` | Correct paper (Martins, Jongen & van Boekel, TIFS Maillard kinetic-modelling review) but CrossRef dates it 2000; computational_priors.json and benchmark_intake_registry.json both cite it as (2001). The repo also stores this DOI in two case variants (uppercase  |
| `10.1021/acs.est.5b02902` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `shirai_2015_bsa_dityrosine_diffusion_limited` | Topic partly matches (dityrosine cross-linking) but the paper is Kampf, Liu, Reinmuth-Selzle et al. 2015 on protein cross-linking upon OZONE exposure (Environ Sci Technol); repo cites 'Shirai et al. (2015)' and claims a BSA diffusion-limited study. |
| `10.1021/acs.jafc.0c04738` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `ding_2020_schiff_base_amadori_emulsion_rates` | Topic matches (Maillard/lipid-oxidation interplay in emulsions) but authors are Troise, Fogliano, Vitaglione et al. 2020; repo cites 'Ding et al. (2020)'. |
| `10.1021/acs.jafc.2c03427` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `bao_2022_xylose_glycylglycine_3_deoxyosone_cleavage` | Topic matches (glycylglycine-catalysed Amadori degradation to deoxyosone) but authors are Cui, Ma, Wang et al. 2022; repo cites 'Bao et al. (2022)'. |
| `10.1021/acs.jafc.2c04919` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `ohsu_2025_kokumi_casr_support_v1`, `ohsu_2025_kokumi_casr_anchor` | Topic matches (kokumi gamma-glutamyl peptides, CaSR/T1R1-T1R3) but the paper is Yang, Liao, Dong et al. 2022 JAFC; repo cites 'Ohsu et al. (2025)'. The stored EC50 constants are attributed to a paper that is not this DOI. |
| `10.1021/acs.jafc.3c08432` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json | `rizzello_2024_fermentation_cleanup`, `acs_jafc_3c08432_crosstalk_cleanup_link` | Correct paper (Flores et al. 2024, fermented texturised pea protein aroma) for id acs_jafc_3c08432_crosstalk_cleanup_link, but the SAME DOI is also anchored to id 'rizzello_2024_fermentation_cleanup' citing 'Rizzello et al. (2024)'. |
| `10.1021/acs.jafc.4c05736` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `huang_2024_dixyl_arp_degradation` | Topic matches (DiXyl-ARP degradation with extra xylose) but authors are Zhang, Cui, Xia et al. 2024; repo cites 'Huang et al. (2024)'. |
| `10.1021/acs.jafc.6b05811` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `li_2017_quercetin_mgo_adduct_browning_mitigation` | Topic matches (quercetin/MGO adducts and alpha-dicarbonyls in lysine/glucose) but authors are Liu, Xia, Lu et al. 2017; repo cites 'Li et al. (2017)'. |
| `10.1021/acs.jafc.7b00264` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `monforte_2018_phenylacetaldehyde_quinones` | Right paper (Monforte et al., gallic acid/glucose/metals in phenylacetaldehyde formation) but published 2017; repo cites 'Monforte et al. (2018)'. |
| `10.1021/acs.jafc.9b07752` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json | `jafc_2019_ref24_polyphenol_thiol_capping` | Topic matches (protein cysteine-phenol adducts in minced beef, i.e. polyphenol-thiol capping) but the paper is Arsad et al. 2020; repo cites it as 'J. Agric. Food Chem. 2019 (Ref. 24)'. |
| `10.1021/bi051902l` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `smagghe_2006_leghemoglobin_oxygen_dissociation` | Right paper (Smagghe et al., ferrous hexacoordinate haemoglobin ligand-binding kinetics incl. leghaemoglobin) but published 2005; repo cites 'Smagghe et al. (2006)'. |
| `10.1021/ef101678s` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `spier_2010_metal_naphthenate_peroxide_decomp` | Topic matches (homogeneous catalysis of hydroperoxide decomposition) but authors/year are West, Adams & Zabarnick 2011 (Energy & Fuels); repo cites 'Spier et al. (2010)'. |
| `10.1021/jf0100797` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mottram_2001_bmfd_retention` | Correct paper (Adams, Mottram, Parker et al. 2001, disulfide interchange between ovalbumin and volatile disulfides) and the citation field says Adams, but the repo paper_id is 'mottram_2001_bmfd_retention'. Also note the matrix is ovalbumin, not soy isolate as |
| `10.1021/jf025561j` | `METADATA-MISMATCH` | lit/arrhenius_params.yml | `amadori` | Davidek, Clety, Aubin et al. 2002 Amadori degradation. ALREADY FLAGGED in-file (comment claims Martins 2003). |
| `10.1021/jf051197n` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `de_vleeschouwer_2006_acrylamide_aqueous`, `claeys_2005_acrylamide_temperature_time` | DOI is Claeys, De Vleeschouwer & Hendrickx 2005 - correct for slr id claeys_2005_acrylamide_temperature_time, but the SAME DOI also anchors 'de_vleeschouwer_2006_acrylamide_aqueous' citing 'De Vleeschouwer et al. (2006)'. Duplicate/misattributed. |
| `10.1021/jf102026c` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zamora_2010_decadienal_asparagine_decarboxylation` | Right paper (asparagine decarboxylation by lipid oxidation products) but first author is Hidalgo (Hidalgo, Delgado, Navarro, Zamora 2010); repo cites 'Zamora et al. (2010)'. Same group - minor. |
| `10.1021/jf3032342` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `jian_2012_ascorbic_ethanolic_degradation` | Topic matches (degradation of ascorbic acid in ethanolic solutions) but authors are Hsu, Tsai, Fu et al. 2012; repo cites 'Jian et al. (2012)'. |
| `10.1021/jf9013394` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `richards_2009_hemoglobin_liposome_oxidation` | Topic matches (haemoglobin-induced lipid oxidation in a liposome model) but authors are Carvajal, Rustad, Mozuraityte et al. 2009; repo cites 'Richards et al. (2009)'. |
| `10.1021/jf9708251` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hendrickx_1998_ascorbic_isobaric_degradation` | Right topic (isobaric-isothermal L-ascorbic acid degradation kinetics) but first authors are Van den Broeck, Ludikhuyze, Weemaes et al. 1998; repo cites 'Hendrickx et al. (1998)' (Hendrickx is a later co-author) and claims a 'squeezed tomato' system not eviden |
| `10.1021/jf970892v` | `METADATA-MISMATCH` | lit/slr_incorporation_matrix.json | `hofmann_schieberle_1997_meat_anchor` | Topic matches exactly (quantification of MFT, FFT, mercaptopentanones in heated meat) but the paper is Kerscher & Grosch 1998; repo cites 'Hofmann & Schieberle (1997/1998)' under paper_id hofmann_schieberle_1997_meat_anchor. |
| `10.1021/jf990120u` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hidalgo_2000_tomato_furosine` | Right paper (Hidalgo & Pompei, HMF and furosine kinetics in tomato products) but published 1999; repo cites 'Hidalgo & Pompei (2000)'. |
| `10.1021/jf9907687` | `METADATA-MISMATCH` | benchmarks/furosine_extrusion_crossover_140C_RamirezJimenez2000.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `ramirez_jimenez_2000_furosine_crossover_benchmark`, `furosine_extrusion_crossover_140C_RamirezJimenez2000` | Correct paper (Ramirez-Jimenez et al. 2000, 'Browning Indicators in Bread'); slr row correctly says toasted_bread, but the intake row records matrix_family='mild_legume_extrudate' for a bread study, and the benchmark file is named furosine_extrusion_crossover. |
| `10.1093/chemse/bjaf043` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mouritsen_2024_umami_thresholds` | Correct paper for the umami-threshold anchor (Tanaka, Itoh & Kondoh 2025) and the citation field says Tanaka, but the repo paper_id is 'mouritsen_2024_umami_thresholds'. |
| `10.1590/1678-457X.08717` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `scielo_brasil_aa_crosslink_hierarchy_anchor` | Right paper (Yu et al., browning kinetics in L-ascorbic acid/basic amino acid systems) but published 2017; repo cites 'Yu et al. (2018)'. Also stored twice in two DOI case variants. |
| `10.1590/1678-457x.08717` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `yu_2018_ascorbic_basic_amino_browning` | Lowercase duplicate of the same Yu et al. 2017 paper; repo cites 2018. |
| `10.3389/fnut.2022.940202` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zhu_2022_acrylamide_lipid_crosstalk` | Topic matches (acrylamide/HMF in glucose-asparagine-linoleic acid kinetic model) but authors are Ma, Long, Li et al. 2022; repo cites 'Zhu et al. (2022)'. |
| `10.3390/antiox9080756` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pereira_2020_metal_pm_haber_weiss_chelation` | Topic adjacent (theoretical study of iron-aminoguanidine complexes, secondary antioxidant activity) but authors are Garcia-Diez & Mora-Diez 2020; repo cites 'Pereira et al. (2020)' for a 'metal PM Haber-Weiss chelation' prior. |
| `10.3390/foods9030251` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json | `soladoye_2020_low_temp_euc_window_v1`, `soladoye_2020_sous_vide_euc_anchor` | Topic matches (sous-vide beef EUC / umami) but authors are Hwang, Ismail & Joo 2020; repo cites 'Soladoye et al. (2020)' in both the registry and computational_priors. |
| `10.3390/molecules25020373` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zamora_2020_pe_dihydropyridine` | Topic matches (carbonyl-amine adducts in PE-fortified soybean oil, radical scavenging) but authors are Goritschnig, Tadus, Koenig et al. 2020; repo cites 'Zamora et al. (2020)'. |
| `10.3390/molecules27196182` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json | `zhao_2022_moromi_precursor_release` | Topic matches (moromi fermentation metabolomics) but authors are Jo, An, Kim et al. 2022; repo cites 'Zhao et al. (2022)'. |
| `10.3390/molecules28073151` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `resconi_2023_pbma_beef_identity_benchmark`, `hernandez_2023_pbma_meat_anchor`, `hernandez_2023_strecker_subanchors` | Correct paper (Hernandez, Woerner, Brooks et al. 2023, sensory/volatiles of PBMAs vs ground beef) - the slr rows cite Hernandez correctly, but the intake row and the benchmark file both call it 'Resconi et al. (2023)'. |
| `10.4014/jmb.2207.07057` | `METADATA-MISMATCH` | lit/benchmark_intake_registry.json | `ahlberg_2021_yeast_extract_grade_anchor` | Topic matches (yeast extract characteristics/production review) but authors are Tao, Yuan, Liu et al. 2022; repo cites 'Ahlberg & Mohammadi (2021)'. |
| `10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9` | `MATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json, lit/slr_incorporation_matrix.json | `van_boekel_2001_maillard_kinetics_review_v1`, `van_boekel_2001_maillard_kinetics_review` | van Boekel 2001 Nahrung Maillard kinetics review - as cited. |
| `10.1002/anie.201300399` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `smuda_glomb_2013_aa_degradation_pathways` | Smuda & Glomb 2013 vitamin C Maillard degradation - as cited. |
| `10.1002/sfp2.1044` | `MATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json, lit/slr_incorporation_matrix.json | `lincoln_2025_polyphenol_crosstalk_v1`, `lincoln_2025`, `lincoln_2025_polyphenol_crosstalk` | Lincoln & Girard 2025 polyphenols in legume protein matrices - as cited. |
| `10.1007/s10068-022-01194-w` | `MATCH` | benchmarks/spi_hvp_xylose_120C_PMC9905368.json, benchmarks/wheat_gluten_hvp_xylose_120C_PMC9905368.json, lit/benchmark_intake_registry.json | `pmc9905368_spi_hvp_xylose_benchmark`, `pmc9905368_wheat_gluten_hvp_xylose_benchmark`, `spi_hvp_xylose_120C_PMC9905368`, `wheat_gluten_hvp_xylose_120C_PMC9905368` | Cho, Park & Kim, HVP + reducing sugar volatiles; CrossRef issue year 2022 vs repo 2023 (print issue) - benign. |
| `10.1016/0146-6380(90)90162-s` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `suzuki_philp_1990_sulfur_melanoidin` | Suzuki & Philp 1990, melanoidins in presence of H2S - as cited. |
| `10.1016/S0924-2244(01)00022-X` | `MATCH` | lit/arrhenius_params.yml | `enolisation` | Martins, Jongen & van Boekel 2000 TIFS Maillard kinetic-modelling review; legitimate estimate anchor (repo says 2001, actual 2000). |
| `10.1016/j.crfs.2025.101173` | `MATCH` | lit/benchmark_intake_registry.json | `li_2025` | Li et al. 2025 heat-induced pea protein modification - as cited. |
| `10.1016/j.fochx.2025.103252` | `MATCH` | lit/benchmark_intake_registry.json | `pmc_12648097_acrylamide_mitigation_anchor` | Zhang et al. 2025 hazards in baked foods review - plausible mitigation anchor as cited. |
| `10.1016/j.foodchem.2004.04.006` | `MATCH` | lit/arrhenius_params.yml | `schiff_condensation` | DOI IDENTITY IS CORRECT: Martins & Van Boekel 2005, 'A kinetic model for the glucose/glycine Maillard reaction pathways', Food Chemistry - the paper the schiff_condensation A_value comment names. NOTE this is a DOI-identity verdict only: a concurrent source-re |
| `10.1016/j.foodchem.2006.11.073` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `serpen_gokmen_2007_ascorbic_redox_kinetics` | Serpen & Gokmen 2007 ascorbic acid redox degradation kinetics - as cited. |
| `10.1016/j.foodchem.2014.09.129` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `sun_2015_ground_beef_ages` | Sun et al. 2015 AGEs in ground beef under pasteurisation - as cited. |
| `10.1016/j.foodchem.2017.12.019` | `MATCH` | lit/slr_incorporation_matrix.json | `troise_2018_soy_amadori` | Troise et al. 2018 free Amadori compounds in soybean products - as cited. |
| `10.1016/j.foodchem.2019.126082` | `MATCH` | benchmarks/pea_isolate_uht_140C_Trikusuma2019.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `trikusuma_2019`, `trikusuma_2019_pea_uht_aroma_panel`, `pea_isolate_uht_140C_Trikusuma2019` | Trikusuma et al. pea protein UHT beverage aroma; CrossRef issue year 2020 vs repo 2019 (online-first) - benign. |
| `10.1016/j.foodchem.2020.126458` | `MATCH` | lit/benchmark_intake_registry.json | `acs_apts_ref24_3dg_arrhenius_anchor` | Yu et al. 2020 furosine/pyrraline kinetics in glucose/lysine - supports the 3-DG Arrhenius anchor as cited. |
| `10.1016/j.foodchem.2020.126500` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zhu_2020_polyphenol_mgo_structure_kinetics` | Zhu et al. 2020 polyphenol structure vs methylglyoxal trapping kinetics - as cited. |
| `10.1016/j.foodchem.2022.132372` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `ishak_2022_phip_hca_kinetics` | Ishak et al. 2022 HCA/PAH kinetics in phenylalanine model - as cited. |
| `10.1016/j.foodchem.2022.134406` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `cga_cys_adduct_sida_2024` | Poojary et al. 2023 caffeic/chlorogenic acid-cysteine adducts - as cited. |
| `10.1016/j.foodchem.2023.137924` | `MATCH` | lit/slr_incorporation_matrix.json | `xu_2023_spi_temporal_release` | Xu et al. off-flavour release from heated SPI; CrossRef issue year 2024 vs repo 2023 (online-first) - benign. |
| `10.1016/j.foodchem.2024.141599` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `nishimura_abe_2024` | Nishimura & Abe, MRPs from soy protein hydrolysates; CrossRef 2025 vs repo 2024 (online-first) - benign. |
| `10.1016/j.foodchem.2025.145652` | `MATCH` | lit/arrhenius_params.yml | `pyrazine_condensation` | Zhan et al. 2025 DFT + experimental pyrazine formation - a defensible estimate anchor for the pyrazine barrier. |
| `10.1016/j.foodchem.2025.146396` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `nguyen_2025_ppi_microencapsulated_oil_stabilization`, `vtechworks_2022_fava_hydrolysis` | Nguyen et al. 2025 mechanistic kinetic model for lipid oxidation in spray-dried microencapsulated oils - as cited for the oil-stabilisation prior. |
| `10.1016/j.foodhyd.2025.111326` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `sun_2025_ppi_vsc_retention` | Sun et al. 2025 pea protein / sulfur flavour compound binding - as cited. |
| `10.1016/j.foodres.2018.06.056` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `voelker_2018_thiamine_salt_degradation` | Voelker et al. 2018 thiamine salt stability kinetics - as cited. |
| `10.1016/j.foodres.2019.03.046` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `comert_gokmen_2019_digestive_mgo_scavenging` | Comert & Gokmen 2019 methylglyoxal scavenging kinetics - as cited. |
| `10.1016/j.freeradbiomed.2019.04.026` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `shu_2019_cysteine_quinone_kinetics` | Shu et al. 2019 quinone-protein adduct kinetics - as cited. |
| `10.1016/j.lwt.2006.09.008` | `MATCH` | lit/slr_incorporation_matrix.json | `cerny_2007_thiamine_partitioning` | Cerny 2007 origin of carbons in sulfur aroma compounds from xylose/cysteine/thiamine - as cited. |
| `10.1016/j.meatsci.2015.04.004` | `MATCH` | lit/slr_incorporation_matrix.json | `watanabe_2015_dmfh_meat` | Watanabe et al. 2015 volatile compounds in aged cooked beef - as cited. |
| `10.1016/s0308-8146(98)00076-4` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mottram_1998_meat_flavour_review` | Mottram 1998 meat flavour review - correctly re-anchored, as its own note states. |
| `10.1021/acs.jafc.0c01761` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zhu_2020_epicatechin_mgo_go_scavenging` | Zhu et al. 2020 epicatechin carbonyl trapping - as cited. |
| `10.1021/acs.jafc.0c01925` | `MATCH` | lit/benchmark_intake_registry.json | `acs_jafc_0c01925_protein_binding_hierarchy` | Anantharamkrishnan et al. 2020 flavour-protein covalent adducts - as cited. |
| `10.1021/acs.jafc.0c04281` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `zha_2020_ppi_glycation_aggregation` | Zha et al. 2020 PPI Maillard glycation/aggregation - as cited. |
| `10.1021/acs.jafc.1c06163` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `jafc_2022_melanoidin_thiol_binding`, `gigl_2021_coffee_thiol_binding` | Gigl, Hofmann & Frank 2021 NMR odorant-melanoidin interactions - as cited (used twice under two ids). |
| `10.1021/acs.jafc.3c02618` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `acs_jafc_3c02618_binding_prior`, `wang_2023_mft_retention` | Wang et al. 2023 aroma-compound binding to myofibrillar/sarcoplasmic protein and collagen - as cited (used under two ids). |
| `10.1021/acs.jafc.3c05991` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `acs_jafc_3c05991_ppi_spi_partitioning`, `jafc_3c05991_hexanal_binding` | Barallat-Perez et al. 2023 flavour binding with commercial food protein isolates - as cited. |
| `10.1021/acs.jafc.5c14296` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hamzalioglu_2026_milk_uht_ages` | Hamzalioglu et al. 2026 multiresponse Maillard kinetics in UHT milk - as cited. |
| `10.1021/acs.jafc.9b04099` | `MATCH` | lit/benchmark_intake_registry.json | `jafc_2019_ref21_pea_gum_arabic_architecture_anchor` | Zha et al. 2019 gum-arabic-mediated glyco-pea protein hydrolysate - as cited (generic 'Ref. 21' citation string). |
| `10.1021/acsfoodscitech.2c00242` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `acs_2022_pba_lysine_loss_benchmark` | Wehrmaker et al. 2022 amino acid modifications during PBMA shearing/sterilisation - as cited. |
| `10.1021/acsfoodscitech.3c00359` | `MATCH` | lit/slr_incorporation_matrix.json | `sen_gokmen_2023_acrylamide_temperature_time` | Sen & Gokmen 2023 multiresponse acrylamide kinetics in low-moisture roasting - as cited. |
| `10.1021/acsfoodscitech.4c00677` | `MATCH` | lit/benchmark_intake_registry.json | `acs_foodscitech_2024_hme_firmness_anchor` | dos Reis Gasparetto Cruz et al. 2025 extruded TVP hamburger analogues - as cited. |
| `10.1021/acsfoodscitech.4c00956` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `takase_2025_lemon_juice_ascorbic_dicarbonyl` | Takase et al. 2025 ascorbic-acid browning/alpha-dicarbonyls in lemon juice - as cited. |
| `10.1021/bm015639p` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `morel_2002_gluten_shear_aggregation` | Morel et al. 2002 heat/shear-mediated wheat gluten aggregation - as cited. |
| `10.1021/bp060389f` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `de_vleeschouwer_2008_acrylamide_aw` | de Vleeschouwer et al. 2008 acrylamide formation/elimination vs water activity - as cited. |
| `10.1021/jf00052a027` | `MATCH` | benchmarks/quarantined/cys_ribose_150C_Mottram1994.json | `$.source_metadata.quarantine_reason` | Mottram & Whitfield 1995 cysteine/ribose/phospholipid volatiles - correct; cited only as the candidate replacement in the quarantine note. |
| `10.1021/jf010789c` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `brands_2002_casein_sugar_melanoidin` | Brands, Wedzicha & van Boekel 2002 melanoidin quantification in sugar-casein - as cited. |
| `10.1021/jf0200826` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mottram_nobrega_2002_furanone_bridge`, `mottram_nobrega_2002` | Mottram & Nobrega 2002 sulfur aroma from cysteine + three ribose forms - as cited. |
| `10.1021/jf034037p` | `MATCH` | lit/benchmark_intake_registry.json | `$.eligible_references[53].citation_aliases[0]` | Yaylayan et al. 2003 thermal decomposition of phosphorylated D-glucoses - as cited. |
| `10.1021/jf0480290` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `rawel_2005_protein_phenolic_binding` | Rawel, Meidtner & Kroll 2005 phenolic-protein binding - as cited. |
| `10.1021/jf050504m` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `knol_2005_acrylamide_kinetics` | Knol et al. 2005 acrylamide kinetic model in glucose-asparagine - as cited. |
| `10.1021/jf060848s` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hidalgo_2006_pe_lysine_antioxidant` | Hidalgo, Leon & Zamora 2006 amino phospholipid antioxidant activity (Rancimat) - as cited. |
| `10.1021/jf062081+` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `munoz_2007_chlorogenic_oxidase_quinone` | Munoz et al. 2007 chlorogenic acid o-quinone kinetics - as cited. |
| `10.1021/jf063024j` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `charissou_2007_cookie_cml_furosine` | Charissou et al. 2007 Maillard indicators in model cookies - as cited. |
| `10.1021/jf070527w` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `hidalgo_2007_decadienal_phenylalanine_styrene` | Hidalgo & Zamora 2007 phenylalanine-to-styrene by 2,4-decadienal - as cited. |
| `10.1021/jf102575r` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `lagrain_2010_cystine_elimination_lanthionine` | Lagrain et al. 2010 cystine beta-elimination and lanthionine in gliadin - as cited. |
| `10.1021/jf3024672` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `rombouts_2012_gluten_crosslinking` | Rombouts et al. 2012 heat-induced gluten cross-linking - as cited. |
| `10.1021/jf801762c` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `cerny_guntz_dubini_2008` | Cerny & Guntz-Dubini 2008 thiamine/cysteine/xylose Maillard - as cited. |
| `10.1021/jf950439o` | `MATCH` | lit/slr_incorporation_matrix.json | `blank_fay_1996_furanones` | Blank & Fay 1996 furanones from pentose sugars - as cited. |
| `10.1021/jf9705983` | `MATCH` | benchmarks/cys_ribose_140C_Hofmann1998.json, lit/slr_incorporation_matrix.json | `hofmann_schieberle_1998_free_mft_fft`, `cys_ribose_140C_Hofmann1998` | Hofmann & Schieberle 1998 quantitative model studies on MFT/FFT precursor systems - as cited, and correct for the cys_ribose_140C_Hofmann1998 benchmark. |
| `10.1021/jp2116033` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `vilanova_2012_pe_schiff_base` | Vilanova et al. 2012 Schiff bases of O-phosphorylethanolamine with PLP - as cited. |
| `10.1021/jp401488j` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `solis_calero_2013_pe_amadori` | Solis-Calero et al. 2013 DFT Amadori rearrangement above a PE surface - as cited. |
| `10.1038/s41598-018-27010-2` | `MATCH` | lit/benchmark_intake_registry.json | `pmc5992167_amadori_pe_burden_anchor` | Kodate et al. 2018 glycated aminophospholipid quantitation in powdered milk - as cited. |
| `10.1039/C4CP05360E` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `solis_calero_2015_pe_glyoxal` | Solis-Calero et al. 2015 DFT carboxymethyl-PE from glyoxal - as cited. |
| `10.1039/D0FO00292E` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `kutzli_2020_pea_maltodextrin_electrospun_glycation` | Kutzli et al. 2020 electrospun pea protein-maltodextrin Maillard glycation - as cited. |
| `10.1046/j.1365-2621.2001.t01-1-00460.x` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `manso_2001_orange_juice_ascorbic_degradation` | Manso et al. 2001 ascorbic acid degradation/browning in orange juice - as cited. |
| `10.1063/5.0288357` | `MATCH` | lit/arrhenius_params.yml | `mutarotation` | Serse et al. 2025 enhanced-sampling kinetics of glucose mutarotation - a defensible anchor for the mutarotation barrier. |
| `10.1073/pnas.0810352106` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `song_2009_benzoquinone_gsh_conjugation` | Song et al. 2009 GSH + quinone nonenzymatic conjugation - as cited (PCB quinones; the entry already flags surrogate transfer limits). |
| `10.1074/jbc.M111.245100` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc3199460_ascorbic_dicarbonyl` | Nemet & Monnier 2011 vitamin C degradation products - as cited. |
| `10.1093/jn/130.4.921S` | `MATCH` | lit/benchmark_intake_registry.json | `yamaguchi_ninomiya_2000_euc_anchor` | Yamaguchi & Ninomiya 2000 'Umami and Food Palatability' - as cited for the EUC constant. |
| `10.1111/1750-3841.17378` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `cao_2024_carp_myoglobin_mrp` | Cao et al. 2024 MRPs inhibiting lipid oxidation via myoglobin stability - as cited. |
| `10.1111/1750-3841.70545` | `MATCH` | lit/benchmark_intake_registry.json | `pmc_2024_pba_cml_cel_ranges_anchor` | Arinzechukwu et al. 2025 plant-based analogue chemical risks - as cited. |
| `10.1111/j.1365-2621.1988.tb13551.x` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `matoba_1988_nucleotide_hydrolysis` | Matoba et al. 1988 thermal degradation of IMP and GMP in aqueous solution - as cited. |
| `10.1155/2015/319505` | `MATCH` | lit/benchmark_intake_registry.json | `pmc_4419266_pe_interfacial_maillard_kinetics` | Solis-Calero et al. 2015 nonenzymatic reactions above phospholipid surfaces - as cited. |
| `10.1186/s13065-021-00773-y` | `MATCH` | lit/benchmark_intake_registry.json | `voelker_2021_thiamine_kinetics` | Voelker et al. 2021 thiamine salt kinetics vs pH/concentration - as cited. |
| `10.1271/bbb.62.893` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `lertsiri_1998_pe_glycation` | Lertsiri et al. 1998 deoxy-D-fructosyl phosphatidylethanolamine - as cited. |
| `10.3389/fnut.2022.1022254` | `MATCH` | lit/benchmark_intake_registry.json, lit/computational_priors.json, lit/slr_incorporation_matrix.json | `frontiers_2022_hcw_aa_arrhenius_v1`, `frontiers_2022_hcw_aa_arrhenius_anchor` | Feng et al. ascorbic acid degradation/browning kinetics in hot-compressed water - as cited. |
| `10.3389/fnut.2022.852225` | `MATCH` | lit/benchmark_intake_registry.json | `asen_2022` | Asen & Aluko 2022 heat-induced pea protein aggregates - as cited. |
| `10.3390/foods10020273` | `MATCH` | lit/slr_incorporation_matrix.json | `wang_2021_peptide_pyrazines` | Wang et al. 2021 pyrazine formation from lysine-containing di/tripeptides - as cited. |
| `10.3390/foods11213505` | `MATCH` | lit/slr_incorporation_matrix.json | `foods_2022_spi_free_sh` | Xie et al. 2022 SPI glycated with galacto-oligosaccharides under HPH; free-SH is plausibly reported. The entry itself already flags 'pending DOI/full-text verification'. |
| `10.3390/foods12061331` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `squeo_2023`, `squeo_2023_pbpi_acrylamide` | Squeo et al. 2023 acrylamide screening in commercial plant-based protein ingredients - as cited. |
| `10.3390/foods12101967` | `MATCH` | benchmarks/cml_cel_commercial_pbma_Foods2023.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `foods_2023_cml_cel_proxy_benchmark`, `acrylamide_heat_trapping_2024`, `foods_2023_pbma_acrylamide_ages`, `cml_cel_commercial_pbma_Foods2023` | Fu et al. 2023 CML/CEL/acrylamide in plant-based meat analogues - as cited (used under four ids and the cml_cel_commercial_pbma benchmark). |
| `10.3390/foods12132543` | `MATCH` | lit/benchmark_intake_registry.json | `wageningen_ref9_hme_rework_hydration_anchor` | Snel et al. 2023 rework potential of soy/pea isolates in HME - as cited. |
| `10.3390/foods12224155` | `MATCH` | lit/slr_incorporation_matrix.json | `laemont_barringer_2023_pyrazine_ph` | Laemont & Barringer 2023 roasted sunflower seed aroma vs pH/sugars/protein - as cited. |
| `10.3390/foods13081257` | `MATCH` | lit/benchmark_intake_registry.json | `pmc11049305_spirulina_offnote_anchor` | Paraskevopoulou et al. 2024 volatile profiling of spirulina supplements - as cited. |
| `10.3390/foods14111940` | `MATCH` | lit/benchmark_intake_registry.json | `pmc12155365_sunflower_roasted_anchor` | Huseynli et al. 2025 sunflower product flavour review - as cited. |
| `10.3390/foods14193295` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `mandelli_2024_plant_milk_hca_quantification` | Mandelli et al. 2025 heterocyclic aromatic amines in plant-based milk - as cited. |
| `10.3390/foods15050912` | `MATCH` | benchmarks/external_validation/external_validation_li_2026_spi_wg_hme_control.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `pmc_2026_hme_hexanal_baseline`, `external_validation_li_2026_spi_wg_hme_control` | Li et al. 2026 enzymatic pretreatment of wheat gluten in HME extrudates - as cited, and correct for the external_validation Li 2026 benchmark. |
| `10.3390/ijms25168668` | `MATCH` | benchmarks/acrylamide_spi_extrusion_130C_ACSRef3.json, lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `acs_ref3_spi_acrylamide_fast_kinetics`, `ma_2024_pbma_extrusion_damage`, `acrylamide_spi_extrusion_130C_ACSRef3` | Ma et al. 2024 extrusion parameters vs CML/CEL/acrylamide in PBMA - as cited (used under two ids and the ACSRef3 benchmark). |
| `10.3390/molecules26134104` | `MATCH` | benchmarks/pea_isolate_40C_PratapSingh2021.json, benchmarks/soy_isolate_40C_PratapSingh2021.json, lit/slr_incorporation_matrix.json | `pratap_singh_2021_native_ppi_spi`, `pea_isolate_40C_PratapSingh2021`, `soy_isolate_40C_PratapSingh2021` | Pratap-Singh et al. 2021 OAV determination in soy/pea/brown-rice proteins - as cited, and correct for both PratapSingh2021 benchmarks. |
| `10.47836/ifrj.28.3.16` | `MATCH` | lit/benchmark_intake_registry.json, lit/slr_incorporation_matrix.json | `yang_2021_ascorbic_glycine_kinetics` | Yang et al. 2021 L-ascorbic acid/glycine browning kinetics - as cited. |

