# What the reaction rules propose, against what the engine models

*What the literature's reaction rules (data/lit/reaction_rules.yml) propose from each lane's reference charge, placed against the engine's own reactions. Steps and products only: no rate, no concentration, and nothing here is read by the engine. Beyond the first step only products that are engine species or registry compounds react further, so the walk stays on the known map.*

26 rules; 128 proposed steps: 64 mechanism known, 22 modelled, 42 proposed; products: 103 new, 4 registry, 30 species.

| placement | meaning |
|---|---|
| modelled | the engine has a reaction with these reactants and this product, with a rate |
| mechanism_known | the rule's source draws this step; the engine has no reaction for it |
| proposed | analogous to a cited rule, on reactants no source shows; no rate, no source |

## sugar glycine (trunk lane, depth 2)

*the trunk's fitted pot (Martins 2005).* Charge: Glc, Gly.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R01_amadori | Glc + Gly | AMA | modelled | r_schiff > r_amadori |
| R02_schiff_base | Glc + Gly | SB | modelled | r_schiff |
| R03_enolisation_1_2 | AMA | Gly + TDG | modelled | r_ama_mgo, r_ama_odg, r_ama_tdg |
| R04_enolisation_2_3 | AMA | ODG + Gly | modelled | r_ama_mgo, r_ama_odg, r_ama_tdg |

## pentose cysteine (sulfur lane, depth 2)

*the sulfur lane's reference pot (Hofmann 1998).* Charge: PENT, Cys.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + N + H2S | modelled | r_cys_h2s |
| R01_amadori | PENT + Cys | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | PENT + Cys | O=C(O)C(CS)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | Cys + CC=O | CC=NC(CS)C(=O)O | mechanism known |  |
| R03_enolisation_1_2 | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | Cys + TDP | mechanism known |  |
| R04_enolisation_2_3 | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | DPO + Cys | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R15_thiazolidine | PENT + Cys | TTCA | mechanism known |  |
| R15_thiazolidine | Cys + CC=O | CC1NC(C(=O)O)CS1 | mechanism known |  |
| R14_hemithioacetal | PENT + Cys | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | Cys + CC=O | CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + CC=O | CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + CC=O | CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O | proposed |  |

Products that are not engine species: `N`; `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `CC=O` (acetaldehyde); `CC(O)SCC(N)C(=O)O`; `CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC1NC(C(=O)O)CS1`; `CC=NC(CS)C(=O)O`; `NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O`; `NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO`; `O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO`.

## pentose cysteine thiamine (sulfur lane, depth 2)

*the thiamine route (Hofmann 1998 Table 8).* Charge: PENT, Cys, THI.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + N + H2S | modelled | r_cys_h2s |
| R23_thiamine_ring_opening | THI | HMP | modelled | r_thi_hmp |
| R24_hmp_to_mft | HMP | MFT | modelled | r_hmp_mft |
| R01_amadori | PENT + Cys | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | PENT + Cys | O=C(O)C(CS)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | Cys + CC=O | CC=NC(CS)C(=O)O | mechanism known |  |
| R03_enolisation_1_2 | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | Cys + TDP | mechanism known |  |
| R04_enolisation_2_3 | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | DPO + Cys | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + HMP | CC(=O)C(CCO)SSCC(N)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + HMP | CC(=O)C(CCO)SSCC(NCC(=O)C(O)C(O)CO)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + HMP | CC(=O)C(CCO)SSCC(N=CC(O)C(O)C(O)CO)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | HMP + HMP | CC(=O)C(CCO)SSC(CCO)C(C)=O | mechanism known |  |
| R15_thiazolidine | PENT + Cys | TTCA | mechanism known |  |
| R15_thiazolidine | Cys + CC=O | CC1NC(C(=O)O)CS1 | mechanism known |  |
| R14_hemithioacetal | PENT + Cys | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + HMP | CC(=O)C(CCO)SC(O)C(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | Cys + CC=O | CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + CC=O | CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + CC=O | CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | CC=O + HMP | CC(=O)C(CCO)SC(C)O | proposed |  |

Products that are not engine species: `N`; `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `CC=O` (acetaldehyde); `CC(=O)C(CCO)SC(C)O`; `CC(=O)C(CCO)SC(O)C(O)C(O)C(O)CO`; `CC(=O)C(CCO)SSC(CCO)C(C)=O`; `CC(=O)C(CCO)SSCC(N)C(=O)O`; `CC(=O)C(CCO)SSCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(=O)C(CCO)SSCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(N)C(=O)O`; `CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC1NC(C(=O)O)CS1`; `CC=NC(CS)C(=O)O`; `NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O`; `NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO`; `O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO`.

## asparagine glucose (acrylamide lane, depth 2)

*the acrylamide lane's pot (De Vleeschouwer).* Charge: Asn, Glc.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R02_schiff_base | Asn + Glc | SBA | modelled | a_asn_glc_sb |
| R17_asparagine_schiff_to_acrylamide | SBA | ACR + NCC(O)C(O)C(O)C(O)CO | modelled | a_sb_int1 > a_int1_acr |
| R01_amadori | Asn + Glc | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | mechanism known |  |
| R03_enolisation_1_2 | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | Asn + TDG | mechanism known |  |
| R04_enolisation_2_3 | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | ODG + Asn | mechanism known |  |

Products that are not engine species: `NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O`; `NCC(O)C(O)C(O)C(O)CO`.

## linoleate hydroperoxides (lipid lane, depth 2)

*the lipid lane's hydroperoxide pool (Frankel 1989).* Charge: LOOH_13_ct, LOOH_9_ct.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R18a_hydroperoxide_scission_alkane_side | LOOH_13_ct/LOOH_13_tt | PENTANE + ME_13_OXO_TRIDECADIENOATE | modelled | lipid_scission_LOOH_13_tt |
| R18a_hydroperoxide_scission_alkane_side | LOOH_9_ct/LOOH_9_tt | DECADIENAL + ME_OCTANOATE | modelled | lipid_scission_LOOH_9_tt |
| R18b_hydroperoxide_scission_aldehyde_side | LOOH_13_ct/LOOH_13_tt | HEXANAL + COC(=O)CCCCCCCC=CCC=O | modelled | lipid_scission_LOOH_13_tt |
| R18b_hydroperoxide_scission_aldehyde_side | LOOH_9_ct/LOOH_9_tt | CCCCCC=CCC=O + ME_9_OXONONANOATE | modelled | lipid_scission_LOOH_9_tt |

Products that are not engine species: `CCCCCC=CCC=O`; `COC(=O)CCCCCCCC=CCC=O`.

## thiol sink probe (sulfur lane, depth 1)

*the two thiols with every carbonyl and thiol partner the pots hold: what could remove them.* Charge: MFT, FFT, MESH, Cys, H2S, PENT, NF, FUR, HMF, MGO, GO, DA, DECADIENAL, HEXANAL, ACR.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + N + H2S | modelled | r_cys_h2s |
| R09_furfural_h2s_to_fft | H2S + FUR | FFT | modelled | r_fur_fft, r_fur_fft_hs |
| R10_norfuraneol_h2s_to_mft | H2S + NF | MFT | modelled | r_nf_mft |
| R12_thiol_oxidation_to_disulfide | MFT + MFT | MFTD | modelled | ch_dimer_mft |
| R12_thiol_oxidation_to_disulfide | MFT + MESH | MMFT | modelled | ch_mmft |
| R12_thiol_oxidation_to_disulfide | FFT + FFT | FFTD | modelled | ch_dimer_fft |
| R13_thiol_michael_addition | Cys + ACR | ACRCYS | modelled | a_acr_cys |
| R27_dicarbonyl_h2s_to_mercaptoketone | H2S + MGO | MP | modelled | r_mgo_mp |
| R01_amadori | Cys + PENT | O=C(O)C(CS)NCC(=O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | Cys + PENT | O=C(O)C(CS)N=CC(O)C(O)C(O)CO | mechanism known |  |
| R02_schiff_base | Cys + FUR | O=C(O)C(CS)N=Cc1ccco1 | mechanism known |  |
| R02_schiff_base | Cys + HMF | O=C(O)C(CS)N=Cc1ccc(CO)o1 | mechanism known |  |
| R02_schiff_base | Cys + MGO | CC(=O)C=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + GO | O=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + DECADIENAL | CCCCCC=CC=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + HEXANAL | CCCCCC=NC(CS)C(=O)O | mechanism known |  |
| R07_strecker | Cys + MGO | CC(N)C=O + O=CCS | mechanism known |  |
| R07_strecker | Cys + DA | CC(=O)C(C)N + O=CCS | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MFT + FFT | Cc1occc1SSCc1ccco1 | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MFT + Cys | Cc1occc1SSCC(N)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | FFT + MESH | CSSCc1ccco1 | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | FFT + Cys | NC(CSSCc1ccco1)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MESH + MESH | CSSC | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MESH + Cys | CSSCC(N)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R13_thiol_michael_addition | MFT + DECADIENAL | CCCCCC=CC(CC=O)Sc1ccoc1C | mechanism known |  |
| R13_thiol_michael_addition | MFT + ACR | Cc1occc1SCCC(N)=O | mechanism known |  |
| R13_thiol_michael_addition | FFT + DECADIENAL | CCCCCC=CC(CC=O)SCc1ccco1 | mechanism known |  |
| R13_thiol_michael_addition | FFT + ACR | NC(=O)CCSCc1ccco1 | mechanism known |  |
| R13_thiol_michael_addition | MESH + DECADIENAL | CCCCCC=CC(CC=O)SC | mechanism known |  |
| R13_thiol_michael_addition | MESH + ACR | CSCCC(N)=O | mechanism known |  |
| R13_thiol_michael_addition | Cys + DECADIENAL | CCCCCC=CC(CC=O)SCC(N)C(=O)O | mechanism known |  |
| R15_thiazolidine | Cys + PENT | TTCA | mechanism known |  |
| R15_thiazolidine | Cys + FUR | O=C(O)C1CSC(c2ccco2)N1 | mechanism known |  |
| R15_thiazolidine | Cys + HMF | O=C(O)C1CSC(c2ccc(CO)o2)N1 | mechanism known |  |
| R15_thiazolidine | Cys + MGO | CC(=O)C1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + GO | O=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + DECADIENAL | CCCCCC=CC=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + HEXANAL | CCCCCC1NC(C(=O)O)CS1 | mechanism known |  |
| R14_hemithioacetal | MFT + PENT | Cc1occc1SC(O)C(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | MFT + FUR | Cc1occc1SC(O)c1ccco1 | proposed |  |
| R14_hemithioacetal | MFT + HMF | Cc1occc1SC(O)c1ccc(CO)o1 | proposed |  |
| R14_hemithioacetal | MFT + MGO | CC(=O)C(O)Sc1ccoc1C | proposed |  |
| R14_hemithioacetal | MFT + GO | Cc1occc1SC(O)C=O | proposed |  |
| R14_hemithioacetal | MFT + DECADIENAL | CCCCCC=CC=CC(O)Sc1ccoc1C | proposed |  |
| R14_hemithioacetal | MFT + HEXANAL | CCCCCC(O)Sc1ccoc1C | proposed |  |
| R14_hemithioacetal | FFT + PENT | OCC(O)C(O)C(O)C(O)SCc1ccco1 | proposed |  |
| R14_hemithioacetal | FFT + FUR | OC(SCc1ccco1)c1ccco1 | proposed |  |
| R14_hemithioacetal | FFT + HMF | OCc1ccc(C(O)SCc2ccco2)o1 | proposed |  |
| R14_hemithioacetal | FFT + MGO | CC(=O)C(O)SCc1ccco1 | proposed |  |
| R14_hemithioacetal | FFT + GO | O=CC(O)SCc1ccco1 | proposed |  |
| R14_hemithioacetal | FFT + DECADIENAL | CCCCCC=CC=CC(O)SCc1ccco1 | proposed |  |
| R14_hemithioacetal | FFT + HEXANAL | CCCCCC(O)SCc1ccco1 | proposed |  |
| R14_hemithioacetal | MESH + PENT | CSC(O)C(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | MESH + FUR | CSC(O)c1ccco1 | proposed |  |
| R14_hemithioacetal | MESH + HMF | CSC(O)c1ccc(CO)o1 | proposed |  |
| R14_hemithioacetal | MESH + MGO | CSC(O)C(C)=O | proposed |  |
| R14_hemithioacetal | MESH + GO | CSC(O)C=O | proposed |  |
| R14_hemithioacetal | MESH + DECADIENAL | CCCCCC=CC=CC(O)SC | proposed |  |
| R14_hemithioacetal | MESH + HEXANAL | CCCCCC(O)SC | proposed |  |
| R14_hemithioacetal | Cys + PENT | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + FUR | NC(CSC(O)c1ccco1)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + HMF | NC(CSC(O)c1ccc(CO)o1)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + MGO | CC(=O)C(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + GO | NC(CSC(O)C=O)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + DECADIENAL | CCCCCC=CC=CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + HEXANAL | CCCCCC(O)SCC(N)C(=O)O | proposed |  |

Products that are not engine species: `CC(=O)C(C)N`; `CC(=O)C(O)SCC(N)C(=O)O`; `CC(=O)C(O)SCc1ccco1`; `CC(=O)C(O)Sc1ccoc1C`; `CC(=O)C1NC(C(=O)O)CS1`; `CC(=O)C=NC(CS)C(=O)O`; `CC(N)C=O`; `CCCCCC(O)SC`; `CCCCCC(O)SCC(N)C(=O)O`; `CCCCCC(O)SCc1ccco1`; `CCCCCC(O)Sc1ccoc1C`; `CCCCCC1NC(C(=O)O)CS1`; `CCCCCC=CC(CC=O)SC`; `CCCCCC=CC(CC=O)SCC(N)C(=O)O`; `CCCCCC=CC(CC=O)SCc1ccco1`; `CCCCCC=CC(CC=O)Sc1ccoc1C`; `CCCCCC=CC=CC(O)SC`; `CCCCCC=CC=CC(O)SCC(N)C(=O)O`; `CCCCCC=CC=CC(O)SCc1ccco1`; `CCCCCC=CC=CC(O)Sc1ccoc1C`; `CCCCCC=CC=CC1NC(C(=O)O)CS1`; `CCCCCC=CC=CC=NC(CS)C(=O)O`; `CCCCCC=NC(CS)C(=O)O`; `CSC(O)C(C)=O`; `CSC(O)C(O)C(O)C(O)CO`; `CSC(O)C=O`; `CSC(O)c1ccc(CO)o1`; `CSC(O)c1ccco1`; `CSCCC(N)=O`; `CSSCC(N)C(=O)O`; `CSSCc1ccco1`; `Cc1occc1SC(O)C(O)C(O)C(O)CO`; `Cc1occc1SC(O)C=O`; `Cc1occc1SC(O)c1ccc(CO)o1`; `Cc1occc1SC(O)c1ccco1`; `Cc1occc1SCCC(N)=O`; `Cc1occc1SSCC(N)C(=O)O`; `Cc1occc1SSCc1ccco1`; `N`; `NC(=O)CCSCc1ccco1`; `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSC(O)C=O)C(=O)O`; `NC(CSC(O)c1ccc(CO)o1)C(=O)O`; `NC(CSC(O)c1ccco1)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `NC(CSSCc1ccco1)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)N=Cc1ccc(CO)o1`; `O=C(O)C(CS)N=Cc1ccco1`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `O=C(O)C1CSC(c2ccc(CO)o2)N1`; `O=C(O)C1CSC(c2ccco2)N1`; `O=CC(O)SCc1ccco1`; `O=CC1NC(C(=O)O)CS1`; `O=CC=NC(CS)C(=O)O`; `O=CCS`; `OC(SCc1ccco1)c1ccco1`; `OCC(O)C(O)C(O)C(O)SCc1ccco1`; `OCc1ccc(C(O)SCc2ccco2)o1`; `CC=O` (acetaldehyde); `CSSC` (dimethyl_disulfide).

## Rules

| rule | what | status | source |
|---|---|---|---|
| R01_amadori | aldose + primary amine -> Amadori compound | net | martins2005_extraction.md: Table 2, steps 1-2 (glucose + glycine -> DFG) |
| R02_schiff_base | aldehyde + primary amine -> Schiff base | established | devleeschouwer2006_extraction.md: the asparagine + glucose Schiff base of the Zyzak mechanism, step 1 of the fitted scheme |
| R03_enolisation_1_2 | Amadori compound -> 3-deoxyosone + amine (1,2-enolisation) | net | martins2005_extraction.md: Table 2, DFG -> 3-deoxyglucosone (step 3) |
| R04_enolisation_2_3 | Amadori compound -> 1-deoxyosone + amine (2,3-enolisation) | net | martins2005_extraction.md: Table 2, DFG -> 1-deoxyglucosone (step 4) |
| R05a_3dg_to_hmf | 3-deoxyhexosone -> HMF (cyclodehydration) | net | kocadagli2016jafc_extraction.md: 3-DG -> HMF step of the fitted scheme (Table 2) |
| R05b_3dp_to_furfural | 3-deoxypentosone -> furfural (cyclodehydration) | net | hofmann1998_reconciliation.md: Table 5, furan-2-aldehyde formed in situ from ribose |
| R06_dicarbonyl_cleavage | alpha-dicarbonyl with an alpha-hydroxyl -> two carbonyl fragments | net | kocadagli2016jafc_extraction.md: glucosone -> glyoxal and 1-DG -> methylglyoxal steps; and hofmann1998_reconciliation.md Table 10, the C2 + C3 fragments of a pentose |
| R07_strecker | alpha-dicarbonyl + alpha-amino acid -> Strecker aldehyde + alpha-aminoketone (+ CO2) | net | hofmann2000_extraction.md: the Strecker degradation with the dicarbonyl donors of Tables 1-3 |
| R08_cysteine_thermolysis | cysteine -> hydrogen sulfide + acetaldehyde + ammonia (+ CO2) | net | zheng1994_extraction.md: Table I, cysteine thermolysis at four pH values (the lane's k_cys_h2s and its barrier) |
| R09_furfural_h2s_to_fft | furfural + H2S -> 2-furfurylthiol | net | hofmann1998_reconciliation.md: Table 3, furan-2-aldehyde + H2S -> FFT, 0.48 mol% |
| R10_norfuraneol_h2s_to_mft | norfuraneol + H2S -> 2-methyl-3-furanthiol | net | hofmann1998_reconciliation.md: Table 4, norfuraneol + H2S -> MFT, 0.19 mol%; whitfield2001_extraction.md Table 1 |
| R11_c2_c3_recombination_to_mft | hydroxyacetaldehyde + 1-mercapto-2-propanone -> 2-methyl-3-furanthiol | net | hofmann1998_reconciliation.md: Table 10, the C2 + C3 recombination pot, 0.24 mol% MFT |
| R12_thiol_oxidation_to_disulfide | two thiols (+ oxidant) -> disulfide | established | kumazawa2003_extraction.md: Table 3, difurfuryl disulfide the major FFT product at 121 C; zhou2023_extraction.md Table 1; mottram2002_extraction.md Table 1 (mixed disulfides) |
| R13_thiol_michael_addition | thiol + alpha,beta-unsaturated carbonyl -> beta-thioether (Michael adduct) | established | devleeschouwer2006_extraction.md: the acrylamide + cysteine elimination channel (k_acr_cys); starkenmann2008_extraction.md for thiol adducts of unsaturated carbonyls |
| R14_hemithioacetal | thiol + aldehyde -> hemithioacetal (reversible) | proposed | hamzalioglu2018_extraction.md: HMF reacting with the thiol group of cysteine (the HMFAD pool); k6b_adduct_kinetics_synthesis.md section 1c for the aldehyde-thiol class |
| R15_thiazolidine | cysteine + aldehyde -> thiazolidine-4-carboxylic acid | established | zhai2020_extraction.md: the xylose-cysteine thiazolidine (TTCA), about 94 % of the group's Cys-Amadori |
| R17_asparagine_schiff_to_acrylamide | asparagine Schiff base -> acrylamide + aminodeoxysugar (+ CO2) | net | devleeschouwer2006_extraction.md: the Schiff base -> acrylamide step of the fitted scheme (k_int1_acr) |
| R18a_hydroperoxide_scission_alkane_side | allylic hydroperoxide -> alkane + 2,4-dienal (beta-scission, side A) | net | schroen2022_extraction.md: the linoleate hydroperoxide product slate (Frankel 1989) the lipid lane's branch fractions come from |
| R18b_hydroperoxide_scission_aldehyde_side | allylic hydroperoxide -> saturated aldehyde + oxo-alkenoate (beta-scission, side B) | net | schroen2022_extraction.md: as R18a: hexanal from the 13-hydroperoxide, methyl 9-oxononanoate from the 9-hydroperoxide |
| R22_dmhf_h2s_to_thiophenone | furaneol + H2S -> 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone | established | shu1988_extraction.md: Table I, the thiophenone from furaneol + H2S (a structural constant; the paper reports area percent only) |
| R23_thiamine_ring_opening | thiamine -> 5-hydroxy-3-mercapto-2-pentanone (thiazole ring opening) | net | hofmann1998_reconciliation.md: Table 8, thiamin -> MFT via 5-hydroxy-3-mercapto-2-pentanone |
| R24_hmp_to_mft | 5-hydroxy-3-mercapto-2-pentanone -> 2-methyl-3-furanthiol (cyclodehydration) | net | hofmann1998_reconciliation.md: Table 8, the thiamin route |
| R25_dideoxypentosone_h2s_to_mft | 1,4-dideoxypentosone + H2S -> 2-methyl-3-furanthiol (the intact-C5 route) | net | hofmann1998_reconciliation.md: Table 6, ribose + H2S -> MFT with the pentose skeleton intact (isotope labelling) |
| R26a_1dp_to_norfuraneol | 1-deoxypentosone -> norfuraneol (cyclisation) | net | hofmann1998_reconciliation.md: Table 5, norfuraneol formed in situ from ribose (54 530 ug/100 mL) |
| R26b_1dg_to_furaneol | 1-deoxyglucosone -> furaneol (cyclisation with the terminal reduction, net) | net | blank1997_extraction.md: furaneol from 1-deoxyglucosone via acetylformoin (the k5b synthesis, the trunk's r_odg_af / r_af_dmhf) |
| R27_dicarbonyl_h2s_to_mercaptoketone | methylglyoxal + H2S -> 1-mercapto-2-propanone | net | hofmann1998_reconciliation.md: Table 7, 2-oxopropanal + H2S (1:1 and 1:2) |
