# Literature metadata triage (2026-11-05)

*Hand triage of the flags in `results/validation/literature_metadata_check_2026-11-05.md` (PR #70, branch
`openalex-tool`): `scripts/tools/openalex.py verify` over `data/keys/papers.yml` (sha256 dce659c7df71), 713 DOIs, 35
flagged rows (33 `author_not_in_record`, 7 `year_mismatch`, 1 `doi_not_found`; flags overlap), no retractions.*

**Scope of what was checked.** Each flagged DOI was looked up on OpenAlex (title, year, venue, PMID / PMCID) and
Crossref (authors, venue, volume, pages, online and issue dates) on 2026-11-05, about one request per second, no
account. Each was compared with every record that cites it (the `record_ids` of its `papers.yml` row) and with the
identifier the record says the DOI was resolved from. For the wrong-DOI candidates the PDF in `data/articles/` was
read where held (`wang2012.pdf`, `Kocadagli2016.pdf`, `Xu2022.pdf`). **Not checked:** whether the NUMBERS in the
backlog descriptions and the side-layer payloads are in the papers, except where stated. Most of those numbers come
from the LLM research digests (`raw/NN_*.md`), and none of them was changed here.

**Where the fixes live.** `papers.yml` is generated (`scripts/generators/build_paper_registry.py`). Its `citation`
is the intake-registry citation, else the SLR ledger's, else the first `citation` / `source_citation` met in a cited
record. So each fix is made in the record. A wrong DOI is replaced, with a `doi_repair` {old, new, date, basis}. A
wrong citation label gets the convention `results/literature/deep_research_backlog.json` already uses: a
`verified_citation` field plus a `citation_repair` {old, new, date, basis}. The scraped `citation` /
`source_citation` stays unchanged, because `src/deep_research_runtime_queue.py` joins records on that string.
Consequence: **the registry still prints the scraped labels for the 11 relabelled papers**, and a re-run of
`verify` will flag them again until the registry builder prefers `verified_citation` (owner item 5).

## Counts

| class | rows | changed |
|---|---|---|
| (a) WRONG DOI: the DOI resolves to a different paper | **2** | DOI replaced in the record |
| (e) WRONG CITATION LABEL: the DOI is right (it is the PMID / PMCID / volume-page the record gives, and the title matches the description), the scraped author or year is wrong | **11** | `verified_citation` + `citation_repair` added |
| (c) YEAR: label year wrong, author right | **1** | `verified_citation` + `citation_repair` added |
| (b) FALSE POSITIVE | **19** | none |
| (d) UNRESOLVED by the tool | **2** | none (one resolved here, one left unresolvable) |
| **total** | **35** | |

Year flags (7): 5 are wrong years in a wrong or confabulated label (Ordoudi, Huang, Kocadağlı, Noda; Ince is year
only). One was already repaired (Ruiz-Aceituno / Zheng & Ho). One is a parser artefact (`PMID:20420424` read as
"2042"). None of the 7 is an online-first vs issue-year difference: the tool tolerates ±1 year, which covers those
(Cho 2022/2023, Tao 2022/2023).

## (a) Wrong DOIs

| # | record (file: id) | old DOI -> what it is | new DOI -> what it is | evidence | model rows |
|---|---|---|---|---|---|
| 1 | `results/literature/deep_research_backlog.json`: item "Wang, Z. et al. (2012)" (`registry_id` `wang_2012_gsh_xylose_sulfur_uplift`) | `10.1016/j.foodchem.2011.07.044` -> Tsai, Yen & Chou, *Changes in Muscovy duck breast muscle marinated with ginger extract*, Food Chem. 130:316-320. It was matched from the note's "Food Chem 130:316", and the 2026-08-26 resolution already said the subject did not match | `10.1016/j.foodchem.2011.08.079` -> Wang, Yang & Song, *Key meat flavour compounds formation mechanism in a glutathione-xylose Maillard reaction*, Food Chem. 131:280-285 (2012) | Title, authors, year and venue read from page 1 of `data/articles/wang2012.pdf` and confirmed on Crossref; dossier `wang2012_extraction.md`. The item describes GSH + xylose with CAMOLA, Wang, Food Chem., 2012, and this is the only such paper | **No engine, fit or benchmark row.** The `registry_id` names no record anywhere in the tree (its runtime artifact was deleted earlier; `status: RUNTIME_BOUND` is stale). The real paper's ledger rows (`sulfur.yml`: `wang2012_table2_gsh_cys_cysgly_xylose_thiamine_volatiles`, `wang2012_table3_camola_xylose_isotopomers`, `wang2012_table4_cysgly_diketopiperazine`; `lipid.yml`: `wang2012_table2_lipid_type_volatiles_no_lipid`) were read from the PDF and are unaffected. **Values from the wrong paper? No: the item's numbers are not from the duck paper, and only partly from Wang 2012.** Bis(2-methyl-3-furyl) disulfide is not reported. There is no SIDA (one internal standard, 2-methyl-3-heptanone). "MFT 2.25x" and "4.4x" were not traced to Table 2. Recorded as `doi_repair.content_caveat` |
| 2 | `data/lit/process_state_calibrations.json`: `hrncirik_2014_coconut_oil_thermal_profile` | `10.1016/j.foodchem.2013.11.082` -> Mashhadizadeh et al., solid-phase extraction of trace metals on modified Fe3O4 nanoparticles (Food Chem.) | `10.1007/s11746-013-2401-2` -> Hrncirik & Zeelenberg, *Stability of Essential Fatty Acids and Formation of Nutritionally Undesirable Compounds in Baking and Shallow Frying*, J. Am. Oil Chem. Soc. 91:591-598 (online 2013-12-27, issue 2014) | The 2026-08-26 sweep repaired the twin record of the same id in `data/lit/benchmark_intake_registry.json` to this DOI and left this one. Authors, title and venue from Crossref | **No engine, fit or benchmark row.** This file is read only by side modules (`src/deep_research_runtime_queue.py`, `src/extrusion_benchmark_landing.py`, `src/extrusion_external_closure.py`), not by `src/kinetic_core/` or the front door. **Values from the wrong paper? The record's content is in NEITHER paper.** The JAOCS abstract covers sunflower oil, rapeseed oil, margarines and butter (EFA retention, TFA, polymerised TAG). It names no coconut oil and no delta-decalactone. The intake twin's `key_values.journal_volume_page: "62:9666"` is a third, unidentified locator. Content caveat recorded on the record |

The brief's own example of a wrong DOI, "Kocadağlı & Gökmen (2016)" -> Kaspchak 2022, is **not** a wrong DOI. The
DOI is the paper behind PMID 35749878, which the source record cites, and its title (phytic acid, tannic acid and
saponin in an asparagine-glucose Maillard reaction) is the record's subject. The label is what is wrong:
Kocadağlı & Gökmen 2016 (`Kocadagli2016.pdf`, 10.1016/j.foodchem.2016.05.150) is a glucose / wheat-flour kinetic
model, and its text contains neither "saponin" nor "acrylamide". It is row e4 below.

## (e) Wrong citation labels (DOI right) and (c) wrong year

Every new author string was read from the Crossref record of the DOI. Every label fix is a `verified_citation` +
`citation_repair` on the record named here.

| # | record (file: id or citation key) | scraped label | verified (Crossref) | identifier the DOI came from; match | model rows |
|---|---|---|---|---|---|
| e1 | backlog: "Bacher, Eberhardt & Müller (2023)" | Bacher, Eberhardt & Müller (2023) | Rocchi, Wolkers-Rooijackers, Liao, Tempelaars & Smid (2023), Yeast 40:628-639 | PMID 37930115 -> this DOI; title (thiamine production in *S. cerevisiae*) = description | none (backlog only) |
| e2 | backlog: "Deng et al. (2024)" | Deng et al. (2024) | Zhu, Zhu, Wang, Cao, Ren & Zhao (2025), Foods 14:3065 | PMC12427843; PPI + *E. faecalis* = description | none |
| e3 | backlog: "Huang et al. (2021)" (`huang_2021_sulfur_oav_support`, an id that exists nowhere else) | Huang et al. (2021) | Zheng, Wei, Liu, Thakur, Zhang & Wei (2023), Curr. Res. Food Sci. 6:100445 | PMC9868338; cysteine-xylose-glutamate GCxGC = description. Year also wrong | none |
| e4 | `data/lit/safety_reference_payloads.json`: `kocadagli_2016_saponin_acrylamide_modifier`, and its backlog item | Kocadağlı & Gökmen (2016) | Kaspchak, Toci, Menezes, Igarashi-Mafra & Mafra (2022), Food Chem. 394:133518 | PMID 35749878; saponin / asparagine-glucose = description. Year also wrong | **side layer only**: the payload (LOD 10 ug/kg, "approximately 6 C thermal shift", qualitative acrylamide modifier) is read by `src/doe_generator.py` and `src/deep_research_runtime_queue.py`, not by the engine. **Its numbers are digest numbers, unchecked against Kaspchak 2022** (owner item 3) |
| e5 | backlog: "Noda, Tanaka & Takahashi (2003)" | Noda, Tanaka & Takahashi (2003) | Ito & Matsuyama (2021), J. Fungi 7:658 | PMC8399179 and "MDPI Fungi 7:658" in the record; koji enzymes = description. Year also wrong | none |
| e6 | backlog: "Ordoudi et al. (2014 / PMC12484514)" (`ordoudi_2014_hmf_peak_window`, intake source withdrawn 2026-10-09) | Ordoudi et al. (2014 / PMC12484514) | Chan, Mokhtar & Hong (2025), J. Food Sci. Technol. 62:1999-2012 | PMC12484514; the description itself says "Systematic review caramel models". The DOI is the review, not Ordoudi 2014, so its HMF numbers are secondary at best | none |
| e7 | backlog: "Rietjens et al. (2023)" | Rietjens et al. (2023) | Rodríguez Valerón, Mak, Jahn, Arboleya & Sörensen (2023), Foods 12:4297 | PMC10706403; fermented cereal by-products, gamma-glutamyl peptides = description | none |
| e8 | backlog: "Rocha et al. (2024)" | Rocha et al. (2024) | Berzosa, Marín-Sánchez, Álvarez, Sánchez-Gimeno & Raso (2024), Foods 13:1916 | "Foods 13:1916"; PEF extraction of glutathione = description | none |
| e9 | backlog: "Soladoye et al. (2020)" (`soladoye_2020_sous_vide_euc_anchor`, exists nowhere else) | Soladoye et al. (2020) | Hwang, Ismail & Joo (2020), Foods 9:251 | "Foods 9:251"; sous-vide beef umami / EUC = description | none |
| e10 | backlog: "Ahlberg & Mohammadi (2021)" (`ahlberg_2021_yeast_extract_grade_anchor`, exists nowhere else) | Ahlberg & Mohammadi (2021) | Tao, Yuan, Liu, Liu, Zhang, Liu, Jiang, Huang & Wang (2022), J. Microbiol. Biotechnol. 33:151-166 | PMC9998214; yeast-extract review = description | none |
| e11 | `results/literature/slr_incorporation_matrix.json`: `pratap_singh_2021_native_ppi_spi` | Pratap-Singh et al. (2021), Molecules 26:4104 | Singh, Shi, Magreault, Kitts, Jarzębski, Siejak & Pratap-Singh (2021), Molecules 26:4104 | Pratap-Singh is the LAST author; the DOI and its values were verified from the full text on 2026-08-27 and 2026-09-11 (dossier `singh2021_extraction.md`) | benchmarks `pea_isolate_40C_PratapSingh2021`, `soy_isolate_40C_PratapSingh2021`; `binding_constants.yml` `pratap_singh_2021_molecules`; ledger `singh2021_tabs1_protein_powder_volatiles_oav`, `singh2021_tab_s1_semiquant_oav`. **Unaffected**: the DOI is right and these rows already name Singh et al. |
| c1 | backlog: "Ince et al. (2024)" | Ince et al. (2024), "Int J Biol Macromol" | Ince, Condict, Ashton, Stockmann & Kasapis (2026), Food Chem. Mol. Sci. 12:100379 | PMC12969450; soy 11S glycinin + hexanal = description. **Wrong year** (2026, June issue) and wrong journal, not online-first | none |

## (b) False positives

| # | paper_id / DOI | why the flag is false | model rows (DOI right, so unaffected) |
|---|---|---|---|
| b1 | `kamal_eldin_2003_triolein_scission` / 10.1002/ejlt.200390032 | OpenAlex and Crossref spell Kamal‐Eldin with U+2010 HYPHEN; the tool's fold does not map it to "-" (already noted 2026-08-27 in the record) | `computational_priors.json` `kamal_eldin_2003_triolein_scission_v1` (side layer) |
| b2 | `pmc9905368_spi_hvp_xylose_benchmark` / 10.1007/s10068-022-01194-w | label is a PMC id; PMC9905368 -> this DOI (Cho, Park & Kim, Food Sci. Biotechnol. 32:283, online 2022-11-30, issue 2023); twin row already says "Cho et al. (2023)" | none (both benchmarks deleted 2026-10-27) |
| b3 | `doi_10_1016_j_foodchem_2014_03_113` | PubMed id label; PMID 24837915 -> this DOI (Shrestha & De Meulenaer 2014, PE-glucose Maillard products) | none |
| b4 | `weerawatanakorn2015` / 10.1016/j.jfda.2015.02.001 | PMC id label; PMC9351765 -> this DOI; dossier on file | none |
| b5 | `doi_10_1016_s0006_3495_85_83868_6` | PMC id label; PMC1329438 -> this DOI (Ralph, Lange & Redfield 1985). **Content flag**: the paper is about proton exchange of PE in vesicles; the description's "PE glycation k ~ 1e-3 M-1 s-1" is not its subject | none |
| b6 | `hamzalioglu_2026_milk_uht_ages` / 10.1021/acs.jafc.5c14296 | dotless ı and ǧ in Hamzalıoǧlu do not fold to "i" / "g" | none |
| b7 | `zha2019` / 10.1021/acs.jafc.9b04099 | label "ACS JAFC 9b04099" names no author; Zha, Yang, Rao & Chen 2019 matches the dossier | ledger `zha2019_si_tabs1_95_volatiles_glyco_pph` |
| b8 | `zheng1994` / 10.1021/bk-1994-0564.ch011 | scraped "Ruiz-Aceituno et al. (compiled 2016)" was already repaired 2026-08-26 (`verified_citation`: Zheng & Ho 1994); the registry reads the scraped key. Also the year flag | `src/kinetic_core/parameters_sulfur.py` L186, L237, L2338 (`k_cys_h2s`); `arrhenius_params.yml` `cysteine_thermolysis`. DOI and dossier right |
| b9 | `doi_10_1021_jf0510687` | PubMed id label; PMID 16248560 -> this DOI (Martín-Hernández et al. 2005, proteins in lecithins) | none |
| b10 | `doi_10_1021_jf100649z` | PubMed id label; PMID 20420424 -> this DOI (Mesías-García et al. 2010). The year flag is the parser reading "2042" from the PMID | none |
| b11 | `tai1998` / 10.1021/jf971111t | scraped "Lee, Jo & Kim (1998)" already repaired 2026-08-26 (`verified_citation`: Tai & Ho 1998) | none |
| b12 | `ramirez_jimenez_2000_furosine_crossover_benchmark` / 10.1021/jf9907687 | OpenAlex writes "Ramı́rez" (dotless ı + combining acute) | none in the engine; its claimed content was already flagged 2026-08-26 |
| b13-b15 | `solis_calero_2013_pe_amadori`, `solis_calero_2015_pe_glyoxal`, `pmc_4419266_pe_interfacial_maillard_kinetics` | OpenAlex prints "Solís Calero" without the hyphen; Crossref has "Solís-Calero" | `computational_priors.json` PE priors (side layer) |
| b16 | `doi_10_31274_rtd_180813_12611` | dissertation; OpenAlex uses a curly apostrophe (O’Keefe) | `binding_constants.yml` `okeefe_1988_isu_thesis` |
| b17 | `foods_2022_spi_free_sh` / 10.3390/foods11213505 | label "Foods 11(21):3505 (2022)" names no author; venue, volume and article number match (Xie et al. 2022). Its own SLR note still says "pending DOI/full-text verification" of the free-SH numbers | `process_state_calibrations.json` `foods_2022_spi_free_sh` (side layer) |
| b18 | `doi_10_3390_foods14162787` | PMC id label; PMC12386030 -> this DOI (Li et al. 2025, soy-sauce melanoidins) | none |
| b19 | `doi_10_3390_molecules31020287` | PMC id label; PMC12844103 -> this DOI (Huynh et al. 2026, protein-polyphenol review) | none |

## (d) Unresolved

| # | paper_id / DOI | finding | action |
|---|---|---|---|
| d1 | `xu2022` / 10.7506/spkx1002-6630-20220613-128 | **Resolved: the DOI is real.** `https://doi.org/ra/` names its registration agency as **ISTIC** (China), not Crossref, so OpenAlex and Crossref cannot see it. The DOI is printed in `data/articles/Xu2022.pdf` (Food Science (China) 43(19):28-35) and in the dossier identity row | none. Model rows `parameters_matrix.py` `xu2022_the_covalent_pair` and ledger `xu2022_figs1_2_binding_rate_heat_pea`, `xu2022_tab2_tab3_methylpyrazine_klotz_quench` stand |
| d2 | `doi_10_1021_acs_jafc_3c07476` (backlog "Zhao, Z. et al. (2023)") | **Unresolvable.** The digest itself printed this DOI. It resolves to Zhang, Yao, Zhou et al. 2023, JAFC 71:20251, *Light-Colored Maillard Peptides* (browning, saltiness, kokumi). That is neither "Zhao" nor, on its abstract, the described "MW fractionation; 1-5 kDa fraction best for MFT". A Crossref search of JAFC 2022-2024 found no Zhao paper on peptide molecular-weight fractions and MFT. Label and DOI disagree, and the description matches neither, so the intended paper cannot be identified | none; backlog only, no model row |

## For the owner (nothing here was changed)

1. **Hrncirik record (a2):** `hrncirik_2014_coconut_oil_thermal_profile` exists in both `process_state_calibrations.json` and `benchmark_intake_registry.json`. Its content (coconut oil, delta-decalactone OAV, "62:9666") is not in the paper its DOI now names. Under the 2026-11-02 rule (values not in the cited paper are withdrawn), it is a withdrawal candidate.
2. **Wang 2012 backlog item (a1):** the description's 2.25x / 4.4x and SIDA claims are unsourced. The real paper is already in the ledger from the PDF.
3. **Kaspchak 2022 payload (e4):** check LOD 10 ug/kg, the ~6 C DSC shift and the >60 % pyrazine skew against the paper, or withdraw them.
4. **Content mismatches seen in passing:** Ralph 1985 (b5) does not measure PE glycation. Ordoudi (e6) is a 2025 review. Ramírez-Jiménez 2000 was already flagged. Foods 11:3505 (b17) is still marked "pending full-text verification".
5. **Tooling:** (i) The registry builder could prefer `verified_citation` over a scraped `citation` when both exist, so that the 11 relabelled papers stop printing confabulated authors. (ii) `openalex.py`'s `_fold` could map U+2010, U+2019, ı and ǧ, and drop hyphens and spaces in surnames; that would retire b1, b6, b12-b16. (iii) The year regex could skip digits inside PMIDs and PMCIDs (b10). (iv) DOIs registered outside Crossref (ISTIC, d1) could be checked against `doi.org/ra` before `doi_not_found` is raised.
