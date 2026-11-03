# LLM-Digest Echo Census — **NOT a verification result**

> ## ⚠ READ THIS BEFORE READING ANY NUMBER BELOW
> 
> **This report verifies nothing. No paper was opened to produce it.**
> 
> It measures ONE thing: whether a registry entry's numbers appear as text within
> ±30 lines of one of its citation surnames inside the **machine-generated markdown**
> under `data/research_corpus/`. Matching is bare substring matching (`4.5`
> matches `14.52`) and surname matching is unanchored (`Bi` matches `binding`).
> No primary source, DOI resolver, publisher or index is contacted at any point.
> 
> Per this repository's own rule (`data/research_corpus/README.md`):
> 
> > "The deep-research report says so" **is not provenance.**
> 
> **Therefore the DIGEST-ECHO class below is the LAUNDERING CENSUS, not a clean bill
> of health.** A high echo count is a bad sign: it means that for those entries this
> repository can point at nothing but its own LLM research dump. Entries in the
> NO DIGEST ECHO class are not thereby better — they simply have no identified
> upstream at all, which may be worse.
> 
> Nothing in this file may be cited as evidence that a value is verified, sourced,
> anchored or safe to ship. The only verification pass that exists in this repo is
> `results/validation/citation_verification_ledger.md`, and even that checks the DOI's
> identity, not the number.
> 
> *History: until 2026-08-27 this file was titled "Numeric Value Traceability and
> **Verification** Report" and printed its top line as "**Fully Verified** (All values
> matched): 153 (57.5%)". Same computation, opposite claim. Retitled and re-worded by
> Wave T3 of the audit remediation; see `scripts/trace_key_values.py` for the full
> old→new vocabulary map and `docs/history/tasks/audit_remediation.md` § Wave T3.*

---

## Summary Metrics
- **Total Entries Evaluated**: 191
- **DIGEST-ECHO (NOT VERIFICATION) — every number echoes the LLM digest corpus**: 50 (26.2%)
- **PARTIAL DIGEST-ECHO (NOT VERIFICATION) — some numbers echo the corpus**: 37 (19.4%)
- **NO DIGEST ECHO — origin unaccounted for even by the digests**: 104 (54.5%)

---

## Entries with NO / PARTIAL digest echo
The numbers listed for each entry were **not** found near their citation surname in the
LLM digest corpus. This says nothing about whether they are right: it says only that the
digests do not account for them, so their origin is unidentified by this script.

| Entry ID | Citation | Sources | Values NOT echoed anywhere in the digest corpus |
|---|---|---|---|
| `acrylamide_heat_trapping_2024` | Fu et al. (2023) | benchmark_intake_registry.json, safety_reference_payloads.json | 31.81, 186.7 |
| `acs_2020_raw_pea_hexanal_baseline` | Bi et al. (2020), J. Agric. Food Chem. 68:2718 | benchmark_intake_registry.json | 4.5, 20.0, 30.0, 50.0, 160.0, 324.0, 1260.0, 4580.0 |
| `acs_2022_pba_lysine_loss` | ACS Food Science & Technology 2022, 2(2):00242 | safety_reference_payloads.json | 2.78, 3.75, 4.05, 4.43, 15.4, 31.4, 125.5 |
| `acs_2022_pba_lysine_loss_benchmark` | Wehrmaker et al. (2022) | benchmark_intake_registry.json, safety_reference_payloads.json | 15.4, 22.0, 26, 31.4, 125.5 |
| `acs_apts_ref24_3dg_arrhenius_anchor` | Yu et al. (2020) | benchmark_intake_registry.json | 4.02, 14.01, 53.45, 81.7 |
| `acs_foodscitech_2024_hme_firmness_anchor` | Cruz et al. (2025) | benchmark_intake_registry.json | 10.0, 20.0, 25.0, 104.0 |
| `acs_jafc_3c08432_crosstalk_cleanup_link` | Flores et al. (2024), J. Agric. Food Chem. 72:5334 | benchmark_intake_registry.json | 1.42, 3.82, 15, 28.4, 95.0, 382.0 |
| `acs_ref3_spi_acrylamide_fast_kinetics` | Ma et al. (2024) | benchmark_intake_registry.json, safety_reference_payloads.json | 22.36, 25.0, 62.62, 112.0, 130.0 |
| `aliani_2005_donor_potency_nucleotide_context` | Aliani & Farmer (2005), J. Agric. Food Chem. 53:6067 | benchmark_intake_registry.json | 1.4, 3.2, 3.8 |
| `aliani_2005_donor_potency_nucleotide_context_v1` | Aliani & Farmer (2005) | computational_priors.json | 1.4, 3.2, 3.8 |
| `arsa_2022_pyrazine_pH9_optimum` | Arsa & Puechkamutr (2022), Journal of Food Science and Technology 59:890 | flavor_reference_payloads.json | 26.26 |
| `asen_2022` | Asen & Aluko (2022) | benchmark_intake_registry.json | 74.45, 124, 206 |
| `bi_2020_raw_pea_hexanal_point` | Bi et al. (2020), J. Agric. Food Chem. 68:2718-2727 | flavor_reference_payloads.json | 4.5, 280.0, 1260.0 |
| `bi_2020_roasted_pea_hexanal_point` | Bi et al. (2020), J. Agric. Food Chem. 68:2718-2727 | flavor_reference_payloads.json | 4.5, 72.0, 324.0 |
| `blank_2001_epoxydecenal_guardrail_v1` | Blank et al. (2001) | computational_priors.json | 0.07 |
| `blank_mottram_2002_ribose_labeling` | Cerny & Davidek (2004), J. Agric. Food Chem. 52:958 | computational_priors.json, benchmark_intake_registry.json | 90.0 |
| `brands_2002_casein_sugar_melanoidin` | Brands et al. (2002) | computational_priors.json, benchmark_intake_registry.json | 1200.0 |
| `buera_1987_maillard_caramelisation_ea` | Buera, Chirife, Resnik & Wetzler (1987), J. Food Sci. 52:1063 | computational_priors.json, benchmark_intake_registry.json | 25.0, 33.0, 35.0, 36.0 |
| `cerny_guntz_dubini_2008` | Cerny & Guntz-Dubini (2008) | computational_priors.json | 0.0035, 0.018, 145 |
| `cga_cys_adduct_sida_2024` | Poojary et al. (2023) | computational_priors.json, benchmark_intake_registry.json | 90.0, 474.1 |
| `charissou_2007_cookie_cml_furosine` | Charissou et al. (2007) | benchmark_intake_registry.json, safety_reference_payloads.json | 220.0 |
| `davidek_2006_thr_glucose_furan` | Limacher, Kerler, Davidek, Schmalzried & Blank (2008), J. Agric. Food Chem. 56:3639 | computational_priors.json, benchmark_intake_registry.json | 0.045 |
| `davidek_2006_thr_glucose_furan_safety` | Limacher, Kerler, Davidek, Schmalzried & Blank (2008), JAFC 56:3639 | safety_reference_payloads.json | 0.045 |
| `de_vleeschouwer_2008_acrylamide_aw` | De Vleeschouwer et al. (2008) | benchmark_intake_registry.json | 0.82 |
| `ding_2020_schiff_base_amadori_emulsion_rates` | Ding et al. (2020) | computational_priors.json | 0.012 |
| `esterbauer_1991_4hne_kinetics` | Esterbauer, Schaur & Zollner (1991), Free Radic. Biol. Med. 11:81 | computational_priors.json, benchmark_intake_registry.json | 0.0001, 0.001, 1.1, 1.33, 7.4, 8.5, 37 |
| `esterbauer_1991_4hne_kinetics_v1` | Esterbauer et al. (1991) | computational_priors.json | 0.0001, 0.001, 1.33, 7.4, 37.0 |
| `farmer_1991_alkyl_thiazoles` | Farmer & Patterson (1991), Food Chemistry 40:201 | computational_priors.json, benchmark_intake_registry.json | 1.2, 12.0, 50.0, 120 |
| `finnigan_2019_mycoprotein_rna` | Finnigan, Wall, Wilde, Stephens et al. (2019), Curr. Dev. Nutr. 3:nzz021 | benchmark_intake_registry.json | 65.0 |
| `foods_2023_cml_cel_proxy_benchmark` | Fu et al. (2023) | benchmark_intake_registry.json, safety_reference_payloads.json | 16.46, 25.21, 32.0, 38.0, 47.61, 55.0, 86.23 |
| `foods_2023_pba_cml_cel_benchmark` | Fu et al. (2023) | safety_reference_payloads.json | 16.46, 25.21, 47.61, 86.23 |
| `foods_2023_pbma_acrylamide_ages` | Fu et al. (2023) | safety_reference_payloads.json | 31.81, 43.07, 64.0, 68.55, 186.7 |
| `frankel_1982_c182_hexanal_scission` | Frankel (1983), Prog. Lipid Res. 22:1 | computational_priors.json, benchmark_intake_registry.json | 0.88, 14.5, 28.0, 45.0 |
| `frankel_1982_c182_hexanal_scission_v1` | Frankel et al. (1982) | computational_priors.json | 0.88, 14.5, 25.0 |
| `glomb_1995_3dg_fragmentation_stoichiometry` | Glomb & Monnier (1995), J. Biol. Chem. 270:10017 | benchmark_intake_registry.json | 0.06, 0.18, 0.28, 0.41 |
| `glomb_1995_3dg_fragmentation_stoichiometry_v1` | Glomb & Monnier (1995) | computational_priors.json | 0.06, 0.18, 0.28, 0.41 |
| `grosch_1982_hexanal_odt_anchor` | Grosch, W. (1982), in Food Flavours Part A (Morton & MacLeod, eds.), Elsevier, pp. 325-398 | flavor_reference_payloads.json | 4.5 |
| `grosch_1999_c183_propanal_scission` | Grosch (1999) | computational_priors.json | 22.0, 52.0 |
| `hendrickx_1998_ascorbic_isobaric_degradation` | Van den Broeck et al. (1998), J. Agric. Food Chem. 46:2001 | computational_priors.json, benchmark_intake_registry.json | 27.15, 54.8, 246.0 |
| `henle_2005_glycinin_cml` | Henle (2005), Amino Acids 29:313 | benchmark_intake_registry.json | 2.3 |
| `hernandez_2023_2_methylbutanal_panel` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 12.74, 15.56, 18.1, 20.74, 23.51 |
| `hernandez_2023_3_methylbutanal_panel` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 15.69, 16.22, 16.26, 16.61, 16.91 |
| `hernandez_2023_benzaldehyde_panel` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 23.65, 36.67, 62.07, 89.35, 164.42 |
| `hernandez_2023_furfural_ratio_anchor` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 20.0, 987.0, 1093.0 |
| `hernandez_2023_gen_mft_point` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 41.09 |
| `hernandez_2023_methional_panel` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 2.28, 3.77, 4.66, 4.89, 8.38 |
| `hernandez_2023_methylpyrazine_pbma_vs_beef` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 10.0, 13.0, 24.0, 39.0 |
| `hernandez_2023_phenylacetaldehyde_panel` | Hernandez et al. (2023), Molecules 28:3151 | flavor_reference_payloads.json | 6.1, 8.54, 10.1, 11.85, 18.21 |
| `hidalgo_2007_decadienal_phenylalanine_styrene` | Hidalgo et al. (2007) | computational_priors.json, benchmark_intake_registry.json | 8.5 |
| `hofmann_1997_beef_fft_band` | Kerscher & Grosch (1998), Journal of Agricultural and Food Chemistry 46:1954 | flavor_reference_payloads.json | 13.0, 42.0 |
| `hofmann_1997_beef_mft_band` | Kerscher & Grosch (1998), Journal of Agricultural and Food Chemistry 46:1954 | flavor_reference_payloads.json | 7.0, 28.0 |
| `hofmann_schieberle_grosch_1996` | Bolton, Reineccius, Liardon & Huynh-Ba (1994), ACS Symp. Ser. 543:270-278 -- REPLACES the nonexistent "Hofmann et al. (1996)" | benchmark_intake_registry.json | 7.5, 8.0, 11.7, 13.0, 14.7, 60.0, 120.0 |
| `jafc_2019_ref24_polyphenol_thiol_capping` | Arsad et al. (2020), J. Agric. Food Chem. 68:8931 | benchmark_intake_registry.json | 0.62, 1.7, 2.2, 8.1 |
| `jafc_2019_ref24_polyphenol_thiol_capping_v1` | J. Agric. Food Chem. 2019 (Ref. 24) | computational_priors.json | 0.62, 1.7, 2.2, 8.1 |
| `jafc_2020_egcg_deoxyosone_trapping` | Yu, Cui, Tang, Hayat et al. (2020), J. Agric. Food Chem. 68:1714 | benchmark_intake_registry.json | 77.8 |
| `jafc_2020_egcg_deoxyosone_trapping_v1` | J. Agric. Food Chem. 2020 (Ref. 25) | computational_priors.json | 62.8, 77.8, 95.0 |
| `jafc_2022_melanoidin_thiol_binding` | Gigl et al. (2021) | computational_priors.json, benchmark_intake_registry.json | 183.0 |
| `jian_2012_ascorbic_ethanolic_degradation` | Hsu et al. (2012), J. Agric. Food Chem. 60:10696 | computational_priors.json, benchmark_intake_registry.json | 43.3, 96.6 |
| `kamal_eldin_2003_triolein_scission` | Kamal-Eldin, Velasco & Dobarganes (2003), Eur. J. Lipid Sci. Technol. 105:329 | computational_priors.json, benchmark_intake_registry.json | 3.5, 4.2, 38.0 |
| `kamal_eldin_2003_triolein_scission_v1` | Kamal-Eldin et al. (2003) | computational_priors.json | 3.5, 4.2 |
| `krause_2003_furosine_hydrolysis_yields` | Krause, Knoll & Henle (2003), Eur. Food Res. Technol. 216:277 | computational_priors.json, benchmark_intake_registry.json, safety_reference_payloads.json | 0.32, 0.42, 0.46, 0.5, 0.51 |
| `lagrain_2010_cystine_elimination_lanthionine` | Lagrain et al. (2010) | computational_priors.json, benchmark_intake_registry.json | 0.0031, 0.054, 120 |
| `li_2026_spi_wg_hme_1_hexanol_control_point` | Li et al. (2026), Foods 15(5):912 | flavor_reference_payloads.json | 20.04 |
| `li_2026_spi_wg_hme_2_pentylfuran_control_point` | Li et al. (2026), Foods 15(5):912 | flavor_reference_payloads.json | 5625.8 |
| `li_2026_spi_wg_hme_heptanal_control_point` | Li et al. (2026), Foods 15(5):912 | flavor_reference_payloads.json | 89.8 |
| `li_2026_spi_wg_hme_hexanal_control_point` | Li et al. (2026), Foods 15(5):912 | flavor_reference_payloads.json | 605.6 |
| `li_2026_spi_wg_hme_nonanal_control_point` | Li et al. (2026), Foods 15(5):912 | flavor_reference_payloads.json | 72.66 |
| `liu_2020_egcg_arp_kinetics` | Yu, Cui, Tang, Hayat et al. (2020), J. Agric. Food Chem. 68:1714 | computational_priors.json, benchmark_intake_registry.json | 8.5, 77.8 |
| `liu_2022_ppi_oav_anchors` | Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | computational_priors.json, flavor_reference_payloads.json, benchmark_intake_registry.json | 0.002, 0.006, 543.0, 11656.0 |
| `liu_2023_ppi_hexanal_band` | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json | 543.0, 2445.0, 11656.0, 12181.0, 52454.0 |
| `liu_2023_ppi_nonanal_band` | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json | 0.188, 3.42 |
| `liu_2023_ppi_offnote_baseline` | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | benchmark_intake_registry.json | 0.188, 3.42, 6.126, 15.5, 25.5, 1628.0, 2445.0, 52454.0 |
| `liu_cadwallader_2023_pea_aeda` | Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json, benchmark_intake_registry.json | 4096.0 |
| `liu_cadwallader_2023_pea_aeda_anchor` | Liu et al. (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json | 4096.0 |
| `ma_2024_pbma_extrusion_damage` | Ma et al. (2024) | computational_priors.json, benchmark_intake_registry.json, safety_reference_payloads.json | 3.2, 18.5, 39, 101, 150, 160.0 |
| `manso_2001_orange_juice_ascorbic_degradation` | Manso et al. (2001) | computational_priors.json, benchmark_intake_registry.json | 55.6 |
| `martins_2003_lys_glucose_kinetics` | Yu, Zhong, Xie & Guo (2020), Food Chemistry 317:126458 | computational_priors.json, benchmark_intake_registry.json | 0.00014, 0.0024, 28.5, 31.2 |
| `martins_van_boekel_2005_ascorbic_amino_browning` | Yu et al. (2017), Food Sci. Technol. (Campinas) 38:537 | computational_priors.json, benchmark_intake_registry.json | 35.31, 54.94 |
| `mottram_1998_meat_flavour_review` | Mottram (1998) | benchmark_intake_registry.json | 150.0 |
| `mouritsen_2024_umami_thresholds` | Amado, Hanselman, Harmon, Deng et al. (2024), Chemical Senses 49:bjad049 | computational_priors.json, benchmark_intake_registry.json | 29.8, 45.2 |
| `nishimura_abe_2024` | Nishimura & Abe (2024) | benchmark_intake_registry.json | 16.5, 62.5 |
| `perdiguero_2004_yeast_autolysis_nucleotides` | Charpentier, Aussenac, Charpentier & Prome (2005), J. Agric. Food Chem. 53:3000 | benchmark_intake_registry.json | 2.8, 24.0, 55.0 |
| `pereira_2020_metal_pm_haber_weiss_chelation` | Garcia-Diez & Mora-Diez (2020), Antioxidants 9:756 | computational_priors.json, benchmark_intake_registry.json | -58.9, -35.8 |
| `pmc10056349_rubisco_amadori` | Balfany, Gutierrez, Moncada & Komarnytsky (2023), Nutrients 15:1327 | computational_priors.json, benchmark_intake_registry.json | 1.9, 2.2, 9.5 |
| `pmc11049305_spirulina_1_octen_3_ol_oav_point` | Paraskevopoulou et al. (2024) | flavor_reference_payloads.json | 10.0 |
| `pmc11049305_spirulina_beta_ionone_oav_floor` | Paraskevopoulou et al. (2024) | flavor_reference_payloads.json | 100.0 |
| `pmc11049305_spirulina_hexanal_oav_band` | Paraskevopoulou et al. (2024) | flavor_reference_payloads.json | 8.0 |
| `pmc11049305_spirulina_offnote_anchor` | Paraskevopoulou et al. (2024) | benchmark_intake_registry.json | 0.04, 0.8, 8.0, 8.4, 10.0, 10.9, 25.7, 47.2, 100.0 |
| `pmc11353891_lentil_deflavoring` | Vurro, De Angelis, Squeo, Caponio, Summo & Pasqualone (2024), Foods 13:2608 | benchmark_intake_registry.json | 78.0 |
| `pmc12155365_sunflower_2_methylbutanal_fd_point` | Huseynli et al. (2025) | flavor_reference_payloads.json | 2048.0 |
| `pmc12155365_sunflower_3_methylbutanal_fd_point` | Huseynli et al. (2025) | flavor_reference_payloads.json | 2048.0 |
| `pmc12155365_sunflower_4_vinylguaiacol_fd_point` | Huseynli et al. (2025) | flavor_reference_payloads.json | 512.0 |
| `pmc12155365_sunflower_roasted_anchor` | Huseynli et al. (2025) | benchmark_intake_registry.json | 0.25, 20.0, 128.0, 140.0, 256.0, 512.0, 2048.0 |
| `pmc3199460_ascorbic_dicarbonyl` | Nemet & Monnier (2011) | computational_priors.json, benchmark_intake_registry.json | 0.5, 4.7, 5.8, 9.1, 28.0 |
| `pmc5992167_amadori_pe_burden_anchor` | Kodate et al. (2018) | benchmark_intake_registry.json | 10.0, 30.0, 73.0, 87.8, 100.0, 326.1, 488.1 |
| `pmc9905368_spi_hvp_xylose_benchmark` | PMC9905368 (2023) | benchmark_intake_registry.json | 0.18, 0.42, 1.88, 120.0 |
| `pmc9905368_wheat_gluten_hvp_xylose_benchmark` | Cho et al. (2023) | benchmark_intake_registry.json | 0.61, 3.44 |
| `pmc_12648097_acrylamide_mitigation` | PMC PMCID:PMC12648097 (Ref. 5) | safety_reference_payloads.json | 0.978 |
| `pmc_12648097_acrylamide_mitigation_anchor` | Zhang et al. (2025) | benchmark_intake_registry.json | 0.978 |
| `pmc_2024_pba_cml_cel_ranges` | Arinzechukwu et al. (2025) | safety_reference_payloads.json | 16.0, 25.0, 110.0 |
| `pmc_2024_pba_cml_cel_ranges_anchor` | Arinzechukwu et al. (2025) | benchmark_intake_registry.json | 16.0, 25.0, 110.0 |
| `pmc_2026_hme_hexanal_baseline` | Li et al. (2026) | computational_priors.json, benchmark_intake_registry.json | 0.4, 0.6, 1.8, 4.6, 7.0, 20.04, 28.8, 29.42, 30.0, 35.0, 40.0, 57.0, 60.0, 89.8, 90.0, 120.0, 140.0, 150.0, 160.0, 221.51, 250.0, 280.0, 605.6 |
| `pmc_4419266_pe_interfacial_maillard_kinetics` | Solís-Calero et al. (2015) | benchmark_intake_registry.json | 82.9, 92.9, 118.0 |
| `pmid_1904866_pentosidine_equivalence_anchor` | PMID 1904866 pentosidine equivalence | computational_priors.json | 7.4, 13.2, 17.0, 37 |
| `ramirez_jimenez_2000_furosine_crossover` | Ramírez-Jiménez et al. (2000) | safety_reference_payloads.json | 140.0 |
| `ramirez_jimenez_2000_furosine_crossover_benchmark` | Ramírez-Jiménez et al. (2000) | benchmark_intake_registry.json, safety_reference_payloads.json | 140.0, 1250.0 |
| `rawel_2002_cga_cysteine_blocking` | Rawel, Rohn, Kruse & Kroll (2002), Food Chemistry 78:443 | computational_priors.json, flavor_reference_payloads.json, benchmark_intake_registry.json | 2.9, 43.0, 85.0 |
| `resconi_2023_pbma_beef_identity_benchmark` | Hernandez et al. (2023), Molecules 28:3151 | benchmark_intake_registry.json | 5.0, 7.65, 11.5, 12.9, 15.72, 17.5, 19.83, 30.45, 31.65, 38.1, 46.8, 50.0, 51.0, 141.0, 170.0, 560.0, 726.0, 1040.0 |
| `researchgate_2023_pea_aeda` | Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json, benchmark_intake_registry.json | 4096.0 |
| `researchgate_2023_pea_aeda_anchor` | Liu et al. (2023), Food Chemistry 406:134998 | flavor_reference_payloads.json | 4096.0 |
| `richards_2009_hemoglobin_liposome_oxidation` | Carvajal et al. (2009), J. Agric. Food Chem. 57:8134 | computational_priors.json, benchmark_intake_registry.json | 0.67, 1.2, 3.8, 37, 56.6, 66.2 |
| `rizzello_2024_fermentation_cleanup` | Flores et al. (2024), J. Agric. Food Chem. 72:5334 | benchmark_intake_registry.json | 1.42, 3.82, 15, 28.4, 382.0 |
| `rombouts_2012_gluten_crosslinking` | Rombouts et al. (2012) | computational_priors.json, benchmark_intake_registry.json | 12.4, 160 |
| `scielo_brasil_aa_crosslink_hierarchy_anchor` | Yu et al. (2017), Food Sci. Technol. (Campinas) 38:537 | computational_priors.json, benchmark_intake_registry.json | 35.31, 50.08, 54.94 |
| `scielo_brasil_aa_crosslink_hierarchy_v1` | SciELO Brasil (Ref. 4) | computational_priors.json | 35.31, 50.08, 54.94 |
| `shu_2019_cysteine_quinone_kinetics` | Shu et al. (2019) | computational_priors.json, benchmark_intake_registry.json | 100000.0, 1000000.0 |
| `smagghe_2006_leghemoglobin_oxygen_dissociation` | Smagghe et al. (2005), Biochemistry 45:561 | computational_priors.json, benchmark_intake_registry.json | 0.024 |
| `solis_calero_2015_pe_glyoxal` | Solís-Calero et al. (2015) | computational_priors.json, benchmark_intake_registry.json | -5.3 |
| `squeo_2023` | Squeo et al. (2023) | benchmark_intake_registry.json | 7, 24, 185, 451, 748 |
| `squeo_2023_pbpi_acrylamide` | Squeo et al. (2023), Foods 12(6):1331 | safety_reference_payloads.json | 0.999, 7.0, 24.0, 185.0, 230.0, 451.0, 748.0 |
| `sun_2025_ppi_vsc_retention` | Sun et al. (2025) | benchmark_intake_registry.json | 16.5, 42.5 |
| `sun_2026_spi_vsc_retention` | Sun et al. (2026) | benchmark_intake_registry.json | -11.97, -10.4, -9.01, 10.43, 14.53, 38.42 |
| `tannenbaum_1985_thiamine_ea` | Mauri, Alzamora, Chirife & Tomio (1989), Int. J. Food Sci. Technol. 24:1 | computational_priors.json, benchmark_intake_registry.json | 78.5, 105.0 |
| `trikusuma_2019` | Trikusuma et al. (2020) | computational_priors.json, benchmark_intake_registry.json | 2.29, 3.1, 5.0, 7, 8.0, 17.2, 24.0, 32.0, 80.0, 140.0, 163.0, 782.0 |
| `uspto_ptacts_2023_yeast_extract_anchor` | Fraser et al. (2018) | benchmark_intake_registry.json | 2000.0, 8000.0 |
| `uspto_ptacts_2023_yeast_extract_mft_oav_band` | US 9,943,096 B2 (2018), Impossible Foods Inc. | flavor_reference_payloads.json | 2000.0, 8000.0 |
| `van_boekel_2001_maillard_kinetics_review` | van Boekel (2001) | benchmark_intake_registry.json | 8.2, 98.4 |
| `van_boekel_2001_maillard_kinetics_review_v1` | van Boekel (2001), DOI:10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9 | computational_priors.json | 0.075, 0.5, 0.65, 0.87, 1.5, 4.6, 5.9, 7.8, 8.2, 10.4, 18.2, 54.7, 71.3, 98.4, 104.0, 118.0, 132.0, 138.0 |
| `voelker_2018_thiamine_salt_degradation` | Voelker et al. (2018) | computational_priors.json, benchmark_intake_registry.json | 96.5, 112.5 |
| `vtechworks_2022_fava_hydrolysis` | Williams (2025) | flavor_reference_payloads.json, benchmark_intake_registry.json | 8192, 16384 |
| `wageningen_ref9_hme_rework_hydration_anchor` | Snel et al. (2023) | benchmark_intake_registry.json | -90.0, -65.0, -63.0, -17.0 |
| `wang_2021_arg_lys_pyrazine_yield` | Wang et al. (2021), Foods 10:273 | flavor_reference_payloads.json | 73.83 |
| `wang_2022_lab_hexanal_cleanup_oav_target` | Tao et al. (2022), Front. Microbiol. 13:1070773 | flavor_reference_payloads.json | 0.47, 0.957 |
| `xu_2025_peptide_hierarchy` | Lee, Lee & Jo (2025), Food Sci. Anim. Resour. 45:1 | computational_priors.json, benchmark_intake_registry.json | 0.14, 2.25, 2.68 |
| `yu_2018_ascorbic_basic_amino_browning` | Yu et al. (2017), Food Sci. Technol. (Campinas) 38:537 | benchmark_intake_registry.json | 35.31, 50.08, 54.94 |
| `zha_2020_ppi_glycation_aggregation` | Zha et al. (2020) | safety_reference_payloads.json | 0.42, 2.8 |
| `zhang_1993_protein_deamidation_ammonia` | Zhang et al. (1993) | computational_priors.json | 0.015 |
| `zhao_2022_moromi_precursor_release` | Jo et al. (2022), Molecules 27:6182 | benchmark_intake_registry.json | 1.99, 2.38, 2.49, 3.02, 3.73, 6.61, 1175.0, 2221.0 |
| `zhu_2020_polyphenol_mgo_structure_kinetics` | Zhu et al. (2020) | computational_priors.json, benchmark_intake_registry.json | 0.027, 0.063 |
| `zhu_2021_braised_chicken_ages` | Zhu et al. (2020), J. Sci. Food Agric. 100:5064 | benchmark_intake_registry.json | 41.52, 42.45, 68.21, 74.87 |
| `zhu_2022_acrylamide_lipid_crosstalk` | Ma et al. (2022), Front. Nutr. 9:940202 | benchmark_intake_registry.json | 12.87, 14.85 |

---

## DIGEST-ECHO entries — every number found only inside LLM-generated text
**NOT VERIFICATION.** For each entry below, every numeric value was located inside
`data/research_corpus/**` within ±30 lines of a citation surname. No paper was
opened. This is the list of entries whose only demonstrated upstream is a machine-
generated digest — i.e. the remediation worklist, not the safe list.

| Entry ID | Citation | Sources | Values echoed in the digest corpus (unverified) |
|---|---|---|---|
| `blank_2001_epoxydecenal_guardrail` | Lin, Fay, Welti & Blank (2001), Lipids 36:749 | benchmark_intake_registry.json | 0.07 |
| `cantre_2007_corned_beef_furosine` | Cantre et al. (2007) | benchmark_intake_registry.json | 95.7 |
| `cao_2024_carp_myoglobin_mrp` | Cao et al. (2024) | computational_priors.json, benchmark_intake_registry.json | 36.95, 68.19 |
| `comert_gokmen_2019_digestive_mgo_scavenging` | Cömert & Gökmen (2019) | computational_priors.json, benchmark_intake_registry.json | 26.6 |
| `fratianni_2016_apricot_furosine` | Fratianni et al. (2017), Food Res. Int. 99:862 | benchmark_intake_registry.json, safety_reference_payloads.json | 83.3 |
| `frontiers_2022_hcw_aa_arrhenius_anchor` | Feng, Yang, Xie, Liu, Peng, Hu & Yu (2023), Frontiers in Nutrition 9:1022254 | benchmark_intake_registry.json | 15.77, 31.7, 47.53, 150.0, 190.0 |
| `frontiers_2022_hcw_aa_arrhenius_v1` | Feng, Yang, Xie, Liu, Peng, Hu & Yu (2023), Frontiers in Nutrition 9:1022254, DOI:10.3389/fnut.2022.1022254 | computational_priors.json | 15.77, 31.7, 47.53, 150.0, 190.0 |
| `gigl_2021_coffee_thiol_binding` | Gigl et al. (2021) | computational_priors.json, benchmark_intake_registry.json | 2.14, 50.0 |
| `grosch_1982_hexanal_odt` | Belitz, Grosch & Schieberle (2009), Food Chemistry, 4th rev. ed., Springer (aroma-threshold tables, Ch. 5) | computational_priors.json, benchmark_intake_registry.json | 1.2, 4.5, 50.0 |
| `hamzalioglu_2026_milk_uht_ages` | Hamzalioglu et al. (2026) | benchmark_intake_registry.json | 52.1 |
| `hidalgo_2000_tomato_furosine` | Hidalgo & Pompei (1999), J. Agric. Food Chem. 48:78 | benchmark_intake_registry.json | 93.9 |
| `hidalgo_2005_pe_ribose_lysine` | Zamora, Nogales & Hidalgo (2004), Eur. Food Res. Technol. 220:459 | computational_priors.json, benchmark_intake_registry.json | 50.0, 66.5 |
| `hidalgo_2006_pe_lysine_antioxidant` | Hidalgo et al. (2006) | computational_priors.json, benchmark_intake_registry.json | 185.0 |
| `hofmann_2001_melanoidin_thioether` | Hofmann & Schieberle (2001), J. Agric. Food Chem. 50:319 | computational_priors.json, benchmark_intake_registry.json | 90.0 |
| `huang_2024_dixyl_arp_degradation` | Zhang et al. (2024), J. Agric. Food Chem. 72:22461 | computational_priors.json, benchmark_intake_registry.json | 85.0, 115.0 |
| `ishak_2022_phip_hca_kinetics` | Ishak et al. (2022) | benchmark_intake_registry.json | 95.36, 114.12 |
| `knol_2005_acrylamide_kinetics` | Knol et al. (2005) | benchmark_intake_registry.json | 52.1, 72.9 |
| `knol_2009_acrylamide_kinetics` | Knol et al. (2010), Food Chemistry 120:1047 | safety_reference_payloads.json | 8 |
| `kocadagli_2016_saponin_acrylamide_modifier` | Kocadağlı & Gökmen (2016) | safety_reference_payloads.json | 10.0 |
| `lertsiri_1998_pe_glycation` | Lertsiri et al. (1998) | computational_priors.json, benchmark_intake_registry.json | 0.05, 0.12 |
| `liu_2023_benzoquinone_lysine_kinetics` | Liu et al. (2023), Food Res. Int. 163:112187 | computational_priors.json, benchmark_intake_registry.json | 19.0 |
| `mdpi_plants_2024_hemp_volatiles` | Chen, Oliveira, Dias & Ismail (2025), Plants 14:274 | flavor_reference_payloads.json, benchmark_intake_registry.json | 35.0 |
| `monforte_2018_phenylacetaldehyde_quinones` | Monforte, Martins & Silva Ferreira (2017), J. Agric. Food Chem. 66:2459 | computational_priors.json, benchmark_intake_registry.json | 31.5, 32.9 |
| `morel_2002_gluten_shear_aggregation` | Morel et al. (2002) | computational_priors.json, benchmark_intake_registry.json | 10.0, 33.7, 100.0 |
| `mottram_1998_meat_flavour_strecker_v1` | Mottram (1998), Food Chemistry | flavor_reference_payloads.json | 38.0 |
| `mottram_2001_bmfd_retention` | Adams, Mottram, Parker & Brown (2001), J. Agric. Food Chem. 49:4333 | flavor_reference_payloads.json, benchmark_intake_registry.json | 0.18, 8.0, 25.0, 98.5 |
| `mottram_2001_bmfd_retention_prior` | Mottram et al. (2001) | computational_priors.json | 8.0, 25.0, 98.5 |
| `mundt_wedzicha_2007_biscuit_browning` | Mundt & Wedzicha (2007), LWT 40:1078 | computational_priors.json, benchmark_intake_registry.json | 105.0 |
| `munoz_2007_chlorogenic_oxidase_quinone` | Munoz et al. (2007) | computational_priors.json, benchmark_intake_registry.json | 2.73 |
| `nashalian_yaylayan_2014_cu_catalyzed_strecker` | Nashalian & Yaylayan (2014), J. Agric. Food Chem. 62:8518 | computational_priors.json, benchmark_intake_registry.json | 0.41 |
| `pmc6104182_soybean_fermentation` | Zhao, Ding, Yao, Cao, Pan & Kong (2018), Front. Microbiol. 9:1872 | benchmark_intake_registry.json | 2.49, 3.73 |
| `pmid36878579_pe_stoichiometry` | Fujimoto, Hayashi, Hamaguchi, Zhan et al. (2023), J. Oleo Sci. 72:311 | computational_priors.json, benchmark_intake_registry.json | 0.5, 10.0, 160.0, 350.0 |
| `poulsen_2023_pbma_cml_cel` | Fu, Ma, Wang, Sun, Chen, Cheng & Liu (2023), Foods 12:1967 | benchmark_intake_registry.json, safety_reference_payloads.json | 16.0, 25.0, 48.0, 86.0 |
| `quintas_2000_sucrose_caramelisation` | Quintas, Guimaraes, Baylina, Brandao & Silva (2007), Innov. Food Sci. Emerg. Technol. 8:306 | computational_priors.json, benchmark_intake_registry.json | 5.0, 48.0, 60.0, 120.0 |
| `serpen_gokmen_2007_ascorbic_redox_kinetics` | Serpen & Gökmen (2007) | computational_priors.json, benchmark_intake_registry.json | 14.0, 88.0 |
| `smuda_glomb_2013_aa_degradation_pathways` | Smuda & Glomb (2013) | computational_priors.json, benchmark_intake_registry.json | 12.0, 31.0, 32.0, 75.0 |
| `solis_calero_2013_pe_amadori` | Solís-Calero et al. (2013) | computational_priors.json, benchmark_intake_registry.json | 8.76, 16.78, 36.7, 70.2 |
| `song_2009_benzoquinone_gsh_conjugation` | Song et al. (2009) | computational_priors.json, benchmark_intake_registry.json | 11.2, 1547.0 |
| `spier_2010_metal_naphthenate_peroxide_decomp` | West, Adams & Zabarnick (2011), Energy & Fuels 25:897 | computational_priors.json, benchmark_intake_registry.json | 11.0, 14.0, 23.0, 26.0 |
| `sun_2015_ground_beef_ages` | Sun et al. (2015) | benchmark_intake_registry.json | 29.21, 61.01 |
| `sun_2020_solid_matrix_cml_cel_accumulation` | Sun, Tang, Wang, Rasco, Lai & Huang (2016), Meat Science 116:1 | computational_priors.json, benchmark_intake_registry.json, safety_reference_payloads.json | 28.0, 40.971, 44.158 |
| `tang_2013_thiamine_mft` | Tang, Jiang, Yuan & Ho (2013), J. Sulfur Chem. 34:326 (CrossRef issue year 2012) | flavor_reference_payloads.json, benchmark_intake_registry.json | 0.05 |
| `vilanova_2012_pe_schiff_base` | Vilanova et al. (2012) | computational_priors.json, benchmark_intake_registry.json | 13.08, 54.7 |
| `wang_2023_mft_retention` | Wang et al. (2023) | flavor_reference_payloads.json, benchmark_intake_registry.json | 7.0, 15.0, 25.0, 42.0, 85.0 |
| `wang_2023_mft_retention_prior` | Wang et al. (2023) | computational_priors.json | 7.0, 25.0, 85.0 |
| `yang_2021_ascorbic_glycine_kinetics` | Yang et al. (2021) | computational_priors.json, benchmark_intake_registry.json | 60.76, 70.16 |
| `yu_2017_cml_cel_meat_review` | Sun, Tang, Wang, Rasco, Lai & Huang (2015), Food Chemistry 172:802 | benchmark_intake_registry.json, safety_reference_payloads.json | 29.21, 61.01 |
| `zamora_2010_decadienal_asparagine_decarboxylation` | Hidalgo, Delgado, Navarro & Zamora (2010), J. Agric. Food Chem. 58:10512 | computational_priors.json, benchmark_intake_registry.json | 22.0, 81.0, 100 |
| `zamora_2020_pe_dihydropyridine` | Goritschnig et al. (2020), Molecules 25:373 | computational_priors.json, benchmark_intake_registry.json | 85.1 |
| `zhu_2020_epicatechin_mgo_go_scavenging` | Zhu et al. (2020) | computational_priors.json, benchmark_intake_registry.json | 0.059, 1.6 |
