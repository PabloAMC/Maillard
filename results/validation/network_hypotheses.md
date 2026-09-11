# What the reaction rules propose, against what the engine models

*What the literature's reaction rules (data/lit/reaction_rules.yml) propose from each lane's reference charge, placed against the engine's own reactions. Steps and products only: no rate, no concentration, and nothing here is read by the engine. Beyond the first step only products that are engine species or registry compounds react further, so the walk stays on the known map.*

38 rules; 270 proposed steps: 172 mechanism known, 26 modelled, 72 proposed; products: 203 new, 8 registry, 69 species.

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
| R03_enolisation_1_2 | AMA | Gly + TDG/DGAL | modelled | r_ama_g, r_ama_mgo, r_ama_odg, r_ama_tdg |
| R04_enolisation_2_3 | AMA | ODG + Gly | modelled | r_ama_g, r_ama_mgo, r_ama_odg, r_ama_tdg |

## pentose cysteine (sulfur lane, depth 2)

*the sulfur lane's reference pot (Hofmann 1998).* Charge: PENT, Cys.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + AMMONIA + H2S | modelled | r_cys_h2s |
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
| R36_hydroxyketone_h2s_to_mercaptoketone | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + H2S | O=C(O)C(CS)NCC(=O)C(S)C(O)CO | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | PENT + AMMONIA | N=CC(O)C(O)C(O)CO | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | CC=O + AMMONIA | CC=N | mechanism known |  |
| R14_hemithioacetal | PENT + Cys | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | Cys + CC=O | CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + CC=O | CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + CC=O | CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O | proposed |  |

Products that are not engine species: `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `CC=O` (acetaldehyde); `CC(O)SCC(N)C(=O)O`; `CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC1NC(C(=O)O)CS1`; `CC=N`; `CC=NC(CS)C(=O)O`; `N=CC(O)C(O)C(O)CO`; `NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O`; `NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O`; `O=C(O)C(CS)NCC(=O)C(S)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO`; `O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO`.

## pentose cysteine thiamine (sulfur lane, depth 2)

*the thiamine route (Hofmann 1998 Table 8).* Charge: PENT, Cys, THI.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + AMMONIA + H2S | modelled | r_cys_h2s |
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
| R36_hydroxyketone_h2s_to_mercaptoketone | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + H2S | O=C(O)C(CS)NCC(=O)C(S)C(O)CO | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | PENT + AMMONIA | N=CC(O)C(O)C(O)CO | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | CC=O + AMMONIA | CC=N | mechanism known |  |
| R14_hemithioacetal | PENT + Cys | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)NCC(=O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + O=C(O)C(CS)N=CC(O)C(O)C(O)CO | O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | PENT + HMP | CC(=O)C(CCO)SC(O)C(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | Cys + CC=O | CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)NCC(=O)C(O)C(O)CO + CC=O | CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | O=C(O)C(CS)N=CC(O)C(O)C(O)CO + CC=O | CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | CC=O + HMP | CC(=O)C(CCO)SC(C)O | proposed |  |

Products that are not engine species: `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `CC=O` (acetaldehyde); `CC(=O)C(CCO)SC(C)O`; `CC(=O)C(CCO)SC(O)C(O)C(O)C(O)CO`; `CC(=O)C(CCO)SSC(CCO)C(C)=O`; `CC(=O)C(CCO)SSCC(N)C(=O)O`; `CC(=O)C(CCO)SSCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(=O)C(CCO)SSCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(N)C(=O)O`; `CC(O)SCC(N=CC(O)C(O)C(O)CO)C(=O)O`; `CC(O)SCC(NCC(=O)C(O)C(O)CO)C(=O)O`; `CC1NC(C(=O)O)CS1`; `CC=N`; `CC=NC(CS)C(=O)O`; `N=CC(O)C(O)C(O)CO`; `NC(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)C(=O)O`; `NC(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)C(=O)O`; `O=C(O)C(CS)NCC(=O)C(S)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSC(O)C(O)C(O)C(O)CO)NCC(=O)C(O)C(O)CO`; `O=C(O)C(CSSCC(N=CC(O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CSSCC(NCC(=O)C(O)C(O)CO)C(=O)O)NCC(=O)C(O)C(O)CO`.

## asparagine glucose (acrylamide lane, depth 2)

*the acrylamide lane's pot (De Vleeschouwer).* Charge: Asn, Glc.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R02_schiff_base | Asn + Glc | SBA | modelled | a_asn_glc_sb |
| R17_asparagine_schiff_to_acrylamide | SBA | ACR + NCC(O)C(O)C(O)C(O)CO | modelled | a_sb_int1 > a_int1_acr |
| R01_amadori | Asn + Glc | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | mechanism known |  |
| R03_enolisation_1_2 | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | Asn + TDG/DGAL | mechanism known |  |
| R04_enolisation_2_3 | NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O | ODG + Asn | mechanism known |  |

Products that are not engine species: `NC(=O)CC(NCC(=O)C(O)C(O)C(O)CO)C(=O)O`; `NCC(O)C(O)C(O)C(O)CO`.

## linoleate hydroperoxides (lipid lane, depth 2)

*the lipid lane's hydroperoxide pool (Frankel 1989) and the 10-hydroperoxide it lacks (Miyazaki 2023).* Charge: LOOH_13_ct, LOOH_9_ct, LOOH_10.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R18a_hydroperoxide_scission_alkane_side | LOOH_13_ct/LOOH_13_tt | PENTANE + ME_13_OXO_TRIDECADIENOATE | modelled | lipid_scission_LOOH_13_tt |
| R18a_hydroperoxide_scission_alkane_side | LOOH_9_ct/LOOH_9_tt | DECADIENAL + ME_OCTANOATE | modelled | lipid_scission_LOOH_9_tt |
| R18b_hydroperoxide_scission_aldehyde_side | LOOH_13_ct/LOOH_13_tt | HEXANAL + COC(=O)CCCCCCCC=CCC=O | modelled | lipid_scission_LOOH_13_tt |
| R18b_hydroperoxide_scission_aldehyde_side | LOOH_9_ct/LOOH_9_tt | CCCCCC=CCC=O + ME_9_OXONONANOATE | modelled | lipid_scission_LOOH_9_tt |
| R18a_hydroperoxide_scission_alkane_side | LOOH_10 | CC=CCCCCC + ME_10_OXO_8_DECENOATE | mechanism known |  |
| R29_oleate_hydroperoxide_scission_alkanal | LOOH_10 | CCCCCC=CCC=O + ME_9_OXONONANOATE | mechanism known |  |
| R30_oleate_hydroperoxide_scission_alkenal | LOOH_10 | CC=CCCCCC + ME_10_OXO_8_DECENOATE | mechanism known |  |
| R31_linoleate_hydroperoxide_furyl_route | LOOH_13_ct/LOOH_13_tt | PENTYLFURAN + ME_9_OXONONANOATE | mechanism known |  |
| R31_linoleate_hydroperoxide_furyl_route | LOOH_9_ct/LOOH_9_tt | HEXANAL + ME_8_FURYL_OCTANOATE | mechanism known |  |
| R32_linoleate_10_hydroperoxide_to_octenol | LOOH_10 | OCTEN3OL + ME_10_OXO_8_DECENOATE | mechanism known |  |

Products that are not engine species: `CC=CCCCCC`; `CCCCCC=CCC=O`; `COC(=O)CCCCCCCC=CCC=O`.

## lipid maillard cross (lipid lane, depth 1)

*the fatty aldehydes an isolate's lipid makes, with the ammonia and hydrogen sulfide the Maillard side supplies: the cross products no lane names (Zamora 2020, Zhou 2000, Du 2023).* Charge: DECADIENAL, HEXANAL, AMMONIA, H2S, Cys.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + AMMONIA + H2S | modelled | r_cys_h2s |
| R02_schiff_base | DECADIENAL + Cys | CCCCCC=CC=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | HEXANAL + Cys | CCCCCC=NC(CS)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R13_thiol_michael_addition | DECADIENAL + Cys | CCCCCC=CC(CC=O)SCC(N)C(=O)O | mechanism known |  |
| R15_thiazolidine | DECADIENAL + Cys | CCCCCC=CC=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | HEXANAL + Cys | CCCCCC1NC(C(=O)O)CS1 | mechanism known |  |
| R33_dienal_ammonia_to_alkylpyridine | DECADIENAL + AMMONIA | PENTYLPYRIDINE | mechanism known |  |
| R34_dienal_h2s_to_alkylthiophene | DECADIENAL + H2S | HEXYLTHIOPHENE | mechanism known |  |
| R35_dienal_h2s_to_alkylthiapyran | DECADIENAL + H2S | PENTYLTHIAPYRAN | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | DECADIENAL + AMMONIA | CCCCCC=CC=CC=N | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | HEXANAL + AMMONIA | HEXANAL_IMINE | mechanism known |  |
| R14_hemithioacetal | DECADIENAL + Cys | CCCCCC=CC=CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | HEXANAL + Cys | CCCCCC(O)SCC(N)C(=O)O | proposed |  |

Products that are not engine species: `CCCCCC(O)SCC(N)C(=O)O`; `CCCCCC1NC(C(=O)O)CS1`; `CCCCCC=CC(CC=O)SCC(N)C(=O)O`; `CCCCCC=CC=CC(O)SCC(N)C(=O)O`; `CCCCCC=CC=CC1NC(C(=O)O)CS1`; `CCCCCC=CC=CC=N`; `CCCCCC=CC=CC=NC(CS)C(=O)O`; `CCCCCC=NC(CS)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `CC=O` (acetaldehyde).

## lipid maillard thiazoles (lipid lane, depth 3)

*Elmore 1997's pot: a lipid alkanal, a Maillard hydroxyketone, ammonia and hydrogen sulfide give the 2-alkyl-3-thiazolines and, oxidised, the registry's 2-alkyl-4-methylthiazoles; the dienal in the same charge shows where the H2S goes instead (Farmer 1990, Mottram 2002).* Charge: HEXANAL, ACETOL, ACETOIN, AMMONIA, H2S, DECADIENAL.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R12_thiol_oxidation_to_disulfide | MP + MP | CC(=O)CSSCC(C)=O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MP + MERCAPTOBUTANONE | CC(=O)CSSC(C)C(C)=O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MERCAPTOBUTANONE + MERCAPTOBUTANONE | CC(=O)C(C)SSC(C)C(C)=O | mechanism known |  |
| R13_thiol_michael_addition | DECADIENAL + MP | CCCCCC=CC(CC=O)SCC(C)=O | mechanism known |  |
| R13_thiol_michael_addition | DECADIENAL + MERCAPTOBUTANONE | CCCCCC=CC(CC=O)SC(C)C(C)=O | mechanism known |  |
| R33_dienal_ammonia_to_alkylpyridine | AMMONIA + DECADIENAL | PENTYLPYRIDINE | mechanism known |  |
| R34_dienal_h2s_to_alkylthiophene | H2S + DECADIENAL | HEXYLTHIOPHENE | mechanism known |  |
| R35_dienal_h2s_to_alkylthiapyran | H2S + DECADIENAL | PENTYLTHIAPYRAN | mechanism known |  |
| R36_hydroxyketone_h2s_to_mercaptoketone | ACETOL + H2S | MP | mechanism known |  |
| R36_hydroxyketone_h2s_to_mercaptoketone | ACETOIN + H2S | MERCAPTOBUTANONE | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | HEXANAL + AMMONIA | HEXANAL_IMINE | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | AMMONIA + DECADIENAL | CCCCCC=CC=CC=N | mechanism known |  |
| R38_mercaptoketone_aldimine_to_thiazoline | MP + HEXANAL_IMINE | THIAZOLINE_2PE_4ME | mechanism known |  |
| R38_mercaptoketone_aldimine_to_thiazoline | MP + CCCCCC=CC=CC=N | CCCCCC=CC=CC1N=C(C)CS1 | mechanism known |  |
| R38_mercaptoketone_aldimine_to_thiazoline | MERCAPTOBUTANONE + HEXANAL_IMINE | THIAZOLINE_2PE_45DM | mechanism known |  |
| R38_mercaptoketone_aldimine_to_thiazoline | MERCAPTOBUTANONE + CCCCCC=CC=CC=N | CCCCCC=CC=CC1N=C(C)C(C)S1 | mechanism known |  |
| R39_thiazoline_to_thiazole | THIAZOLINE_2PE_4ME | CCCCCc1nc(C)cs1 | mechanism known |  |
| R39_thiazoline_to_thiazole | THIAZOLINE_2PE_45DM | CCCCCc1nc(C)c(C)s1 | mechanism known |  |
| R14_hemithioacetal | HEXANAL + MP | CCCCCC(O)SCC(C)=O | proposed |  |
| R14_hemithioacetal | HEXANAL + MERCAPTOBUTANONE | CCCCCC(O)SC(C)C(C)=O | proposed |  |
| R14_hemithioacetal | DECADIENAL + MP | CCCCCC=CC=CC(O)SCC(C)=O | proposed |  |
| R14_hemithioacetal | DECADIENAL + MERCAPTOBUTANONE | CCCCCC=CC=CC(O)SC(C)C(C)=O | proposed |  |

Products that are not engine species: `CCCCCC=CC=CC=N`; `CC(=O)C(C)SSC(C)C(C)=O`; `CC(=O)CSSC(C)C(C)=O`; `CC(=O)CSSCC(C)=O`; `CCCCCC(O)SC(C)C(C)=O`; `CCCCCC(O)SCC(C)=O`; `CCCCCC=CC(CC=O)SC(C)C(C)=O`; `CCCCCC=CC(CC=O)SCC(C)=O`; `CCCCCC=CC=CC(O)SC(C)C(C)=O`; `CCCCCC=CC=CC(O)SCC(C)=O`; `CCCCCC=CC=CC1N=C(C)C(C)S1`; `CCCCCC=CC=CC1N=C(C)CS1`; `CCCCCc1nc(C)c(C)s1`; `CCCCCc1nc(C)cs1` (2_pentyl_4_methylthiazole).

## oleate hydroperoxides (lipid lane, depth 2)

*the four oleate hydroperoxides the lipid lane lumps as LOOH_OL with no edge (nonanal is a declared hold-out): Cao 2020's routes.* Charge: OL_8_OOH, OL_9_OOH, OL_10_OOH, OL_11_OOH.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R18a_hydroperoxide_scission_alkane_side | OL_8_OOH | ME_HEPTANOATE + UNDECENAL_2E | mechanism known |  |
| R18a_hydroperoxide_scission_alkane_side | OL_9_OOH | ME_OCTANOATE + DECENAL_2E | mechanism known |  |
| R18a_hydroperoxide_scission_alkane_side | OL_10_OOH | CCCCCCCC + ME_10_OXO_8_DECENOATE | mechanism known |  |
| R18a_hydroperoxide_scission_alkane_side | OL_11_OOH | CCCCCCC + COC(=O)CCCCCCCC=CC=O | mechanism known |  |
| R29_oleate_hydroperoxide_scission_alkanal | OL_8_OOH | DECANAL + ME_8_OXOOCTANOATE | mechanism known |  |
| R29_oleate_hydroperoxide_scission_alkanal | OL_9_OOH | NONANAL + ME_9_OXONONANOATE | mechanism known |  |
| R29_oleate_hydroperoxide_scission_alkanal | OL_10_OOH | NONANAL + ME_9_OXONONANOATE | mechanism known |  |
| R29_oleate_hydroperoxide_scission_alkanal | OL_11_OOH | OCTANAL + ME_10_OXODECANOATE | mechanism known |  |
| R30_oleate_hydroperoxide_scission_alkenal | OL_8_OOH | ME_HEPTANOATE + UNDECENAL_2E | mechanism known |  |
| R30_oleate_hydroperoxide_scission_alkenal | OL_9_OOH | ME_OCTANOATE + DECENAL_2E | mechanism known |  |
| R30_oleate_hydroperoxide_scission_alkenal | OL_10_OOH | CCCCCCCC + ME_10_OXO_8_DECENOATE | mechanism known |  |
| R30_oleate_hydroperoxide_scission_alkenal | OL_11_OOH | CCCCCCC + COC(=O)CCCCCCCC=CC=O | mechanism known |  |

Products that are not engine species: `CCCCCCC`; `CCCCCCCC`; `COC(=O)CCCCCCCC=CC=O`.

## strecker to pyrazines (trunk lane, depth 2)

*the small dicarbonyls with an amino acid: the Strecker aldehydes, the aminoketones and the pyrazines they condense to (no lane has a pyrazine).* Charge: MGO, GO, DA, Ala, Cys.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + AMMONIA + H2S | modelled | r_cys_h2s |
| R27_dicarbonyl_h2s_to_mercaptoketone | MGO + H2S | MP | modelled | r_mgo_mp |
| R28_aminoketone_condensation_to_pyrazine | AKM + AKM | DMP | modelled | r_akm_dmp |
| R02_schiff_base | MGO + Ala | CC(=O)C=NC(C)C(=O)O | mechanism known |  |
| R02_schiff_base | MGO + Cys | CC(=O)C=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | GO + Ala | CC(N=CC=O)C(=O)O | mechanism known |  |
| R02_schiff_base | GO + Cys | O=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | MGO + AKM | CC(=O)C=NCC(C)=O | mechanism known |  |
| R02_schiff_base | MGO + CC(=O)C(C)N | CC(=O)C=NC(C)C(C)=O | mechanism known |  |
| R02_schiff_base | GO + AKM | CC(=O)CN=CC=O | mechanism known |  |
| R02_schiff_base | GO + CC(=O)C(C)N | CC(=O)C(C)N=CC=O | mechanism known |  |
| R02_schiff_base | Ala + CC(N=CC=O)C(=O)O | CC(N=CC=NC(C)C(=O)O)C(=O)O | mechanism known |  |
| R02_schiff_base | Ala + O=CC=NC(CS)C(=O)O | CC(N=CC=NC(CS)C(=O)O)C(=O)O | mechanism known |  |
| R02_schiff_base | Ala + CC=O | CC=NC(C)C(=O)O | mechanism known |  |
| R02_schiff_base | Ala + O=CCS | CC(N=CCS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + CC(N=CC=O)C(=O)O | CC(N=CC=NC(CS)C(=O)O)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + O=CC=NC(CS)C(=O)O | O=C(O)C(CS)N=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + CC=O | CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | Cys + O=CCS | O=C(O)C(CS)N=CCS | mechanism known |  |
| R02_schiff_base | CC(N=CC=O)C(=O)O + AKM | CC(=O)CN=CC=NC(C)C(=O)O | mechanism known |  |
| R02_schiff_base | CC(N=CC=O)C(=O)O + CC(=O)C(C)N | CC(=O)C(C)N=CC=NC(C)C(=O)O | mechanism known |  |
| R02_schiff_base | O=CC=NC(CS)C(=O)O + AKM | CC(=O)CN=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | O=CC=NC(CS)C(=O)O + CC(=O)C(C)N | CC(=O)C(C)N=CC=NC(CS)C(=O)O | mechanism known |  |
| R02_schiff_base | AKM + CC=O | CC=NCC(C)=O | mechanism known |  |
| R02_schiff_base | AKM + O=CCS | CC(=O)CN=CCS | mechanism known |  |
| R02_schiff_base | CC=O + CC(=O)C(C)N | CC=NC(C)C(C)=O | mechanism known |  |
| R02_schiff_base | O=CCS + CC(=O)C(C)N | CC(=O)C(C)N=CCS | mechanism known |  |
| R07_strecker | MGO + Ala | AKM + CC=O | mechanism known |  |
| R07_strecker | MGO + Cys | AKM + O=CCS | mechanism known |  |
| R07_strecker | DA + Ala | CC(=O)C(C)N + CC=O | mechanism known |  |
| R07_strecker | DA + Cys | CC(=O)C(C)N + O=CCS | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + CC(=O)C=NC(CS)C(=O)O | CC(=O)C=NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=CC=NC(CS)C(=O)O | NC(CSSCC(N=CC=O)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + O=CCS | NC(CSSCC=O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | CC(=O)C=NC(CS)C(=O)O + CC(=O)C=NC(CS)C(=O)O | CC(=O)C=NC(CSSCC(N=CC(C)=O)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | CC(=O)C=NC(CS)C(=O)O + O=CC=NC(CS)C(=O)O | CC(=O)C=NC(CSSCC(N=CC=O)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | CC(=O)C=NC(CS)C(=O)O + O=CCS | CC(=O)C=NC(CSSCC=O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=CC=NC(CS)C(=O)O + O=CC=NC(CS)C(=O)O | O=CC=NC(CSSCC(N=CC=O)C(=O)O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=CC=NC(CS)C(=O)O + O=CCS | O=CC=NC(CSSCC=O)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | O=CCS + O=CCS | O=CCSSCC=O | mechanism known |  |
| R15_thiazolidine | MGO + Cys | CC(=O)C1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | GO + Cys | O=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + CC(N=CC=O)C(=O)O | CC(N=CC1NC(C(=O)O)CS1)C(=O)O | mechanism known |  |
| R15_thiazolidine | Cys + O=CC=NC(CS)C(=O)O | O=C(O)C(CS)N=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + CC=O | CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + O=CCS | O=C(O)C1CSC(CS)N1 | mechanism known |  |
| R28_aminoketone_condensation_to_pyrazine | AKM + CC(=O)C(C)N | Cc1cnc(C)c(C)n1 | mechanism known |  |
| R28_aminoketone_condensation_to_pyrazine | CC(=O)C(C)N + CC(=O)C(C)N | Cc1nc(C)c(C)nc1C | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | MGO + AMMONIA | CC(=O)C=N | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | GO + AMMONIA | N=CC=O | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | CC(N=CC=O)C(=O)O + AMMONIA | CC(N=CC=N)C(=O)O | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | O=CC=NC(CS)C(=O)O + AMMONIA | N=CC=NC(CS)C(=O)O | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | CC=O + AMMONIA | CC=N | mechanism known |  |
| R37_alkanal_ammonia_to_aldimine | O=CCS + AMMONIA | N=CCS | mechanism known |  |
| R14_hemithioacetal | MGO + Cys | CC(=O)C(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | GO + Cys | NC(CSC(O)C=O)C(=O)O | proposed |  |
| R14_hemithioacetal | MGO + CC(=O)C=NC(CS)C(=O)O | CC(=O)C=NC(CSC(O)C(C)=O)C(=O)O | proposed |  |
| R14_hemithioacetal | MGO + O=CC=NC(CS)C(=O)O | CC(=O)C(O)SCC(N=CC=O)C(=O)O | proposed |  |
| R14_hemithioacetal | MGO + O=CCS | CC(=O)C(O)SCC=O | proposed |  |
| R14_hemithioacetal | GO + CC(=O)C=NC(CS)C(=O)O | CC(=O)C=NC(CSC(O)C=O)C(=O)O | proposed |  |
| R14_hemithioacetal | GO + O=CC=NC(CS)C(=O)O | O=CC=NC(CSC(O)C=O)C(=O)O | proposed |  |
| R14_hemithioacetal | GO + O=CCS | O=CCSC(O)C=O | proposed |  |
| R14_hemithioacetal | Cys + CC(N=CC=O)C(=O)O | CC(N=CC(O)SCC(N)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + O=CC=NC(CS)C(=O)O | NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + CC=O | CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + O=CCS | NC(CSC(O)CS)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(=O)C=NC(CS)C(=O)O + CC(N=CC=O)C(=O)O | CC(=O)C=NC(CSC(O)C=NC(C)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(=O)C=NC(CS)C(=O)O + O=CC=NC(CS)C(=O)O | CC(=O)C=NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(=O)C=NC(CS)C(=O)O + CC=O | CC(=O)C=NC(CSC(C)O)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(=O)C=NC(CS)C(=O)O + O=CCS | CC(=O)C=NC(CSC(O)CS)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(N=CC=O)C(=O)O + O=CC=NC(CS)C(=O)O | CC(N=CC(O)SCC(N=CC=O)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | CC(N=CC=O)C(=O)O + O=CCS | CC(N=CC(O)SCC=O)C(=O)O | proposed |  |
| R14_hemithioacetal | O=CC=NC(CS)C(=O)O + O=CC=NC(CS)C(=O)O | O=CC=NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O | proposed |  |
| R14_hemithioacetal | O=CC=NC(CS)C(=O)O + CC=O | CC(O)SCC(N=CC=O)C(=O)O | proposed |  |
| R14_hemithioacetal | O=CC=NC(CS)C(=O)O + O=CCS | O=CC=NC(CSC(O)CS)C(=O)O | proposed |  |
| R14_hemithioacetal | O=CC=NC(CS)C(=O)O + O=CCS | O=CCSC(O)C=NC(CS)C(=O)O | proposed |  |
| R14_hemithioacetal | CC=O + O=CCS | CC(O)SCC=O | proposed |  |
| R14_hemithioacetal | O=CCS + O=CCS | O=CCSC(O)CS | proposed |  |

Products that are not engine species: `CC(=O)C(C)N`; `CC(=O)C(O)SCC(N)C(=O)O`; `CC(=O)C1NC(C(=O)O)CS1`; `CC(=O)C=NC(C)C(=O)O`; `CC(=O)C=NC(CS)C(=O)O`; `CC(N=CC=O)C(=O)O`; `NC(CSC(O)C=O)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `O=CC1NC(C(=O)O)CS1`; `O=CC=NC(CS)C(=O)O`; `O=CCS`; `CC=O` (acetaldehyde); `CC(=O)C(C)N=CC=NC(C)C(=O)O`; `CC(=O)C(C)N=CC=NC(CS)C(=O)O`; `CC(=O)C(C)N=CC=O`; `CC(=O)C(C)N=CCS`; `CC(=O)C(O)SCC(N=CC=O)C(=O)O`; `CC(=O)C(O)SCC=O`; `CC(=O)C=N`; `CC(=O)C=NC(C)C(C)=O`; `CC(=O)C=NC(CSC(C)O)C(=O)O`; `CC(=O)C=NC(CSC(O)C(C)=O)C(=O)O`; `CC(=O)C=NC(CSC(O)C=NC(C)C(=O)O)C(=O)O`; `CC(=O)C=NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O`; `CC(=O)C=NC(CSC(O)C=O)C(=O)O`; `CC(=O)C=NC(CSC(O)CS)C(=O)O`; `CC(=O)C=NC(CSSCC(N)C(=O)O)C(=O)O`; `CC(=O)C=NC(CSSCC(N=CC(C)=O)C(=O)O)C(=O)O`; `CC(=O)C=NC(CSSCC(N=CC=O)C(=O)O)C(=O)O`; `CC(=O)C=NC(CSSCC=O)C(=O)O`; `CC(=O)C=NCC(C)=O`; `CC(=O)CN=CC=NC(C)C(=O)O`; `CC(=O)CN=CC=NC(CS)C(=O)O`; `CC(=O)CN=CC=O`; `CC(=O)CN=CCS`; `CC(N=CC(O)SCC(N)C(=O)O)C(=O)O`; `CC(N=CC(O)SCC(N=CC=O)C(=O)O)C(=O)O`; `CC(N=CC(O)SCC=O)C(=O)O`; `CC(N=CC1NC(C(=O)O)CS1)C(=O)O`; `CC(N=CC=N)C(=O)O`; `CC(N=CC=NC(C)C(=O)O)C(=O)O`; `CC(N=CC=NC(CS)C(=O)O)C(=O)O`; `CC(N=CCS)C(=O)O`; `CC(O)SCC(N)C(=O)O`; `CC(O)SCC(N=CC=O)C(=O)O`; `CC(O)SCC=O`; `CC1NC(C(=O)O)CS1`; `CC=N`; `CC=NC(C)C(=O)O`; `CC=NC(C)C(C)=O`; `CC=NC(CS)C(=O)O`; `CC=NCC(C)=O`; `N=CC=NC(CS)C(=O)O`; `N=CC=O`; `N=CCS`; `NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O`; `NC(CSC(O)CS)C(=O)O`; `NC(CSSCC(N=CC=O)C(=O)O)C(=O)O`; `NC(CSSCC=O)C(=O)O`; `O=C(O)C(CS)N=CC1NC(C(=O)O)CS1`; `O=C(O)C(CS)N=CC=NC(CS)C(=O)O`; `O=C(O)C(CS)N=CCS`; `O=C(O)C1CSC(CS)N1`; `O=CC=NC(CSC(O)C=NC(CS)C(=O)O)C(=O)O`; `O=CC=NC(CSC(O)C=O)C(=O)O`; `O=CC=NC(CSC(O)CS)C(=O)O`; `O=CC=NC(CSSCC(N=CC=O)C(=O)O)C(=O)O`; `O=CC=NC(CSSCC=O)C(=O)O`; `O=CCSC(O)C=NC(CS)C(=O)O`; `O=CCSC(O)C=O`; `O=CCSC(O)CS`; `O=CCSSCC=O`; `Cc1cnc(C)c(C)n1` (trimethylpyrazine); `Cc1nc(C)c(C)nc1C` (tetramethylpyrazine).

## thiol sink probe (sulfur lane, depth 1)

*the two thiols with every carbonyl and thiol partner the pots hold: what could remove them.* Charge: MFT, FFT, MESH, Cys, H2S, PENT, NF, FUR, HMF, MGO, GO, DA, DECADIENAL, HEXANAL, ACR.

| rule | reactants | products | placement | engine reaction |
|---|---|---|---|---|
| R08_cysteine_thermolysis | Cys | CC=O + AMMONIA + H2S | modelled | r_cys_h2s |
| R09_furfural_h2s_to_fft | H2S + FUR | FFT | modelled | r_fur_fft, r_fur_fft_hs |
| R10_norfuraneol_h2s_to_mft | H2S + NF | MFT | modelled | r_nf_mft |
| R12_thiol_oxidation_to_disulfide | MFT + MFT | MFTD | modelled | ch_dimer_mft |
| R12_thiol_oxidation_to_disulfide | MFT + MSH/MESH | MMFT | modelled | ch_mmft |
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
| R07_strecker | Cys + MGO | AKM + O=CCS | mechanism known |  |
| R07_strecker | Cys + DA | CC(=O)C(C)N + O=CCS | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MFT + FFT | Cc1occc1SSCc1ccco1 | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MFT + Cys | Cc1occc1SSCC(N)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | FFT + MSH/MESH | CSSCc1ccco1 | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | FFT + Cys | NC(CSSCc1ccco1)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MSH/MESH + MSH/MESH | DMDS | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | MSH/MESH + Cys | CSSCC(N)C(=O)O | mechanism known |  |
| R12_thiol_oxidation_to_disulfide | Cys + Cys | NC(CSSCC(N)C(=O)O)C(=O)O | mechanism known |  |
| R13_thiol_michael_addition | MFT + DECADIENAL | CCCCCC=CC(CC=O)Sc1ccoc1C | mechanism known |  |
| R13_thiol_michael_addition | MFT + ACR | Cc1occc1SCCC(N)=O | mechanism known |  |
| R13_thiol_michael_addition | FFT + DECADIENAL | CCCCCC=CC(CC=O)SCc1ccco1 | mechanism known |  |
| R13_thiol_michael_addition | FFT + ACR | NC(=O)CCSCc1ccco1 | mechanism known |  |
| R13_thiol_michael_addition | MSH/MESH + DECADIENAL | CCCCCC=CC(CC=O)SC | mechanism known |  |
| R13_thiol_michael_addition | MSH/MESH + ACR | CSCCC(N)=O | mechanism known |  |
| R13_thiol_michael_addition | Cys + DECADIENAL | CCCCCC=CC(CC=O)SCC(N)C(=O)O | mechanism known |  |
| R15_thiazolidine | Cys + PENT | TTCA | mechanism known |  |
| R15_thiazolidine | Cys + FUR | O=C(O)C1CSC(c2ccco2)N1 | mechanism known |  |
| R15_thiazolidine | Cys + HMF | O=C(O)C1CSC(c2ccc(CO)o2)N1 | mechanism known |  |
| R15_thiazolidine | Cys + MGO | CC(=O)C1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + GO | O=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + DECADIENAL | CCCCCC=CC=CC1NC(C(=O)O)CS1 | mechanism known |  |
| R15_thiazolidine | Cys + HEXANAL | CCCCCC1NC(C(=O)O)CS1 | mechanism known |  |
| R34_dienal_h2s_to_alkylthiophene | H2S + DECADIENAL | HEXYLTHIOPHENE | mechanism known |  |
| R35_dienal_h2s_to_alkylthiapyran | H2S + DECADIENAL | PENTYLTHIAPYRAN | mechanism known |  |
| R36_hydroxyketone_h2s_to_mercaptoketone | H2S + NF | CC1=C(S)C(=O)CO1 | mechanism known |  |
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
| R14_hemithioacetal | MSH/MESH + PENT | CSC(O)C(O)C(O)C(O)CO | proposed |  |
| R14_hemithioacetal | MSH/MESH + FUR | CSC(O)c1ccco1 | proposed |  |
| R14_hemithioacetal | MSH/MESH + HMF | CSC(O)c1ccc(CO)o1 | proposed |  |
| R14_hemithioacetal | MSH/MESH + MGO | CSC(O)C(C)=O | proposed |  |
| R14_hemithioacetal | MSH/MESH + GO | CSC(O)C=O | proposed |  |
| R14_hemithioacetal | MSH/MESH + DECADIENAL | CCCCCC=CC=CC(O)SC | proposed |  |
| R14_hemithioacetal | MSH/MESH + HEXANAL | CCCCCC(O)SC | proposed |  |
| R14_hemithioacetal | Cys + PENT | NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + FUR | NC(CSC(O)c1ccco1)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + HMF | NC(CSC(O)c1ccc(CO)o1)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + MGO | CC(=O)C(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + GO | NC(CSC(O)C=O)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + DECADIENAL | CCCCCC=CC=CC(O)SCC(N)C(=O)O | proposed |  |
| R14_hemithioacetal | Cys + HEXANAL | CCCCCC(O)SCC(N)C(=O)O | proposed |  |

Products that are not engine species: `CC(=O)C(C)N`; `CC(=O)C(O)SCC(N)C(=O)O`; `CC(=O)C(O)SCc1ccco1`; `CC(=O)C(O)Sc1ccoc1C`; `CC(=O)C1NC(C(=O)O)CS1`; `CC(=O)C=NC(CS)C(=O)O`; `CC1=C(S)C(=O)CO1`; `CCCCCC(O)SC`; `CCCCCC(O)SCC(N)C(=O)O`; `CCCCCC(O)SCc1ccco1`; `CCCCCC(O)Sc1ccoc1C`; `CCCCCC1NC(C(=O)O)CS1`; `CCCCCC=CC(CC=O)SC`; `CCCCCC=CC(CC=O)SCC(N)C(=O)O`; `CCCCCC=CC(CC=O)SCc1ccco1`; `CCCCCC=CC(CC=O)Sc1ccoc1C`; `CCCCCC=CC=CC(O)SC`; `CCCCCC=CC=CC(O)SCC(N)C(=O)O`; `CCCCCC=CC=CC(O)SCc1ccco1`; `CCCCCC=CC=CC(O)Sc1ccoc1C`; `CCCCCC=CC=CC1NC(C(=O)O)CS1`; `CCCCCC=CC=CC=NC(CS)C(=O)O`; `CCCCCC=NC(CS)C(=O)O`; `CSC(O)C(C)=O`; `CSC(O)C(O)C(O)C(O)CO`; `CSC(O)C=O`; `CSC(O)c1ccc(CO)o1`; `CSC(O)c1ccco1`; `CSCCC(N)=O`; `CSSCC(N)C(=O)O`; `CSSCc1ccco1`; `Cc1occc1SC(O)C(O)C(O)C(O)CO`; `Cc1occc1SC(O)C=O`; `Cc1occc1SC(O)c1ccc(CO)o1`; `Cc1occc1SC(O)c1ccco1`; `Cc1occc1SCCC(N)=O`; `Cc1occc1SSCC(N)C(=O)O`; `Cc1occc1SSCc1ccco1`; `NC(=O)CCSCc1ccco1`; `NC(CSC(O)C(O)C(O)C(O)CO)C(=O)O`; `NC(CSC(O)C=O)C(=O)O`; `NC(CSC(O)c1ccc(CO)o1)C(=O)O`; `NC(CSC(O)c1ccco1)C(=O)O`; `NC(CSSCC(N)C(=O)O)C(=O)O`; `NC(CSSCc1ccco1)C(=O)O`; `O=C(O)C(CS)N=CC(O)C(O)C(O)CO`; `O=C(O)C(CS)N=Cc1ccc(CO)o1`; `O=C(O)C(CS)N=Cc1ccco1`; `O=C(O)C(CS)NCC(=O)C(O)C(O)CO`; `O=C(O)C1CSC(c2ccc(CO)o2)N1`; `O=C(O)C1CSC(c2ccco2)N1`; `O=CC(O)SCc1ccco1`; `O=CC1NC(C(=O)O)CS1`; `O=CC=NC(CS)C(=O)O`; `O=CCS`; `OC(SCc1ccco1)c1ccco1`; `OCC(O)C(O)C(O)C(O)SCc1ccco1`; `OCc1ccc(C(O)SCc2ccco2)o1`; `CC=O` (acetaldehyde).

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
| R07_strecker | alpha-dicarbonyl + alpha-amino acid -> Strecker aldehyde + alpha-aminoketone (+ CO2); the amine lands on the carbonyl carbon that was the aldehyde, or either carbonyl of a diketone | net | hofmann2000_extraction.md: the Strecker degradation with the dicarbonyl donors of Tables 1-3 |
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
| R28_aminoketone_condensation_to_pyrazine | two alpha-aminoketones -> 2,5-disubstituted pyrazine (condensation, dehydration, oxidation; net) | net | zhou2023_extraction.md: the paper's mechanism (Fig. 4): the open-chain aminoketone condenses to the pyrazines, and cysteine's Strecker route supplies more alpha-aminoketones (section on pyrazine promotion) |
| R29_oleate_hydroperoxide_scission_alkanal | mono-ene allylic hydroperoxide -> alkanal + oxo-ester (beta-scission on the hydroperoxide carbon; Cao 2020's "B-scission") | net | cao2020_extraction.md: Figure 6 II-V: 8-OOH -> decanal, 9-OOH and 10-OOH -> nonanal, 11-OOH -> octanal, each with its oxo-glyceride; Table 2 levels at 120-180 C; chen2017_extraction.md corroborates 9-OOH -> nonanal on the free acid |
| R30_oleate_hydroperoxide_scission_alkenal | mono-ene allylic hydroperoxide -> 2-alkenal + alkane-ended ester (scission on the far side; Cao 2020's "A-scission") | net | cao2020_extraction.md: Figure 6 II-III: 8-OOH -> 2-undecenal, 9-OOH -> 2-decenal; Table 2; chen2017_extraction.md: 2-undecenal from 8-OOH and 2-decenal + octanoic acid from 9-OOH, stated |
| R31_linoleate_hydroperoxide_furyl_route | conjugated-diene hydroperoxide -> 2-alkylfuran + oxo-ester (alkoxyl cyclisation onto the diene's far carbon, oxygen, scission; net) | net | miyazaki2023_extraction.md: Figure 2, route LA-13-F: 13-HpODE -> 2-pentylfuran + 9-oxononanoic acid (24 -> 25 -> 26 -> 22 + 23); the mirror route LA-9-F: 9-HpODE -> hexanal + the furan-bearing C12 acid; Table 1: 2-pentylfuran from the 9-, 10- and 13-hydroperoxides and NOT from the 12-hydroperoxide; corroborated by yao2024_extraction.md: the dihydrofuran intermediate (methyl 9-hydroperoxy-9-(5-pentyl-2,5-dihydrofuran-2-yl)nonanoate) detected and quantified at 0.59 mmol/kg in neat methyl linoleate at 180 C, 1 h |
| R32_linoleate_10_hydroperoxide_to_octenol | non-conjugated 10-hydroperoxide -> 10-oxo-8-enoate + 1-octen-3-ol (scission to the 2-octenyl radical, allyl shift, oxygen, reduction; net) | net | miyazaki2023_extraction.md: Figure 3, route LA-10-B: 10-HpODE -> 10-oxo-8-decenoic acid + 2-octenyl radical (40 <-> 41) -> 1-octen-3-ol (the dominant product, area 144.8 M), 1-octen-3-one, 2-octenal, 2-octen-1-ol; yao2024_extraction.md draws the same scission for the 10-hydroperoxide (methyl 10-oxo-8-decenoate + 1-octen-3-ol, both detected at 180 C) |
| R33_dienal_ammonia_to_alkylpyridine | 2,4-alkadienal + ammonia -> 2-alkylpyridine (imine, cyclisation, aromatisation; net) | established | zhou2000_extraction.md: Tables 1-3 and the isotope experiments: 5-13C-2,4-decadienal gives 2-13C-2-pentylpyridine and 15N-ammonia gives the ring nitrogen, in a defatted soy slurry at room temperature, pH 9 > 7 > 4.5; zamora2020_extraction.md Table 2: 2-pentylpyridine 14.68 +/- 0.58 umol per mmol glutamine from 2,4-decadienal at 180 C (authentic standard), the C6 to C10 dienal homologue series; du2023_extraction.md Table 1: 2-butylpyridine 191 ug/L from (E,E)-2,4-nonadienal in a cysteine + glucose pot at 150 C |
| R34_dienal_h2s_to_alkylthiophene | 2,4-alkadienal + H2S -> 2-alkylthiophene (1,4-addition, S onto C1, dehydration; net) | established | farmer1990_extraction.md: Figure 1 middle branch (H2S adds S to C4 and H to C5, S closes on C1, water leaves; 2-hexylthiophene from 2,4-decadienal) and Table 2: 2-hexylthiophene 0 / 0 / 184 / 1220 / 436 relative area for no lipid / triglyceride / lecithin / PC / PE; mottram2002b_extraction.md Figure 3 bottom branch and Table I: 2-pentylthiophene 21, 2-hexylthiophene 9 ng per 0.5 mmol methyl linoleate at 140 C; elmore1997_extraction.md: 2-pentylthiophene 21 % of total area from (E,E)-2,4-nonadienal + acetoin + (NH4)2S; whitfield1988_extraction.md Table: 2-pentylthiophene 95.5 ng only with cysteine + lecithin, and the ratio argument that rejects the furan -> thiophene exchange |
| R35_dienal_h2s_to_alkylthiapyran | 2,4-alkadienal + H2S -> 2-alkyl-2H-thiapyran (1,6-addition, S onto C1, dehydration; net) | established | farmer1990_extraction.md: Figure 1 right branch (H2S adds S to C5 and H to C4, S closes on C1, water leaves) and Table 2: 2-pentyl-2H-thiapyran 0 / 35 / 3150 / 34700 / 12500 relative area, the largest lipid-dependent product in every phospholipid pot; mottram2002b_extraction.md Figure 3 top branch and Table I: 2-ethylthiapyran 399 ng per 0.5 mmol methyl linolenate, and the text that a pure dienal + H2S pot gives up to 100 times more thiapyran than thiophene |
| R36_hydroxyketone_h2s_to_mercaptoketone | alpha-hydroxyketone + H2S -> alpha-mercaptoketone + H2O | established | elmore1997_extraction.md: Figure 1, first step: substitution of the hydroxyketone's -OH by -SH (E97-A); the mercaptoketones themselves are the compounds Farmer 1990 Table 2 rows 44 to 47 measure in cysteine + ribose pots (2-mercapto-3-butanone 0.97 / 0.45 / 0.51 / 0.48 with the four lipids) |
| R37_alkanal_ammonia_to_aldimine | alkanal + NH3 -> aldimine + H2O | established | elmore1997_extraction.md: Figure 1, the imine of the alkanal (E97-B), drawn as the N-donor of the 3-thiazoline's C2; the alkanal series C4 to C10 in Tables 1 and 2 |
| R38_mercaptoketone_aldimine_to_thiazoline | alpha-mercaptoketone + aldimine -> 2-alkyl-3-thiazoline + H2O (S-attack on the imine carbon, condensation; net) | established | elmore1997_extraction.md: Figure 1 (E97-C): the mercaptoketone's S attacks the imine carbon, the amine condenses onto the ketone; Table 1: 4,5-dimethyl-2-pentyl-3-thiazoline 36.7 % of area from acetoin + hexanal + (NH4)2S, 4-methyl-2-pentyl-3-thiazoline 18.3 % from acetol + hexanal (the 4-isomer favoured over the 5-isomer 5.5 to 18.3 versus 1.5 to 4.4 %); 3-thiazolines 15 to 42 % of area across the C4 to C10 alkanals |
| R39_thiazoline_to_thiazole | 3-thiazoline -> thiazole (oxidation; net, terminal) | established | elmore1997_extraction.md: Figure 1 [O] arrow (E97-D); Tables 1 and 2: thiazole / 3-thiazoline area ratio 0.002 to 0.05 in the hydroxyketone pots (hexanal + acetol 0.7 / 18.3) because the ammonium-sulfide pot is reducing, about 1 in the butanedione pots (pentanal 12.8 / 13.4, hexanal 16.8 / 15.7, heptanal 6.9 / 8.0), which go by the dione route (E97-E) without an oxidation step; whitfield1988_extraction.md: 2-acetylthiazole rises 14.8 -> 25.2 ng with lecithin while every other thiazole falls |
