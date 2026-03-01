#!/usr/bin/env python
# encoding: utf-8

name = "Lipid_Schiff_Base"
shortDesc = "Condensation of an aliphatic lipid aldehyde with a free amino acid to form a non-volatile Schiff base."
longDesc = """
Pathway D: Lipid-Maillard Crosstalk & Off-Flavour Trapping
A crucial mitigation strategy in plant-based matrices. Polyunsaturated fat oxidation 
yields foul-tasting volatile aldehydes like hexanal and nonanal. 
These electrophilic aldehydes condense with the nucleophilic amines of added free 
amino acids (originally intended for the Maillard cascade) to form stable, 
non-volatile Schiff bases, physically trapping the off-flavor.
R-CHO + R'-NH2 => R-CH=N-R' + H2O
"""

recipe(name="Lipid_Schiff_Base",
       recipe=[
           ['FORM_BOND', '*1', 1, '*3'], # C-N bond formation between lipid aldehyde and amine
           ['BREAK_BOND', '*1', 2, '*2'], # Break C=O double bond (becomes single)
           ['BREAK_BOND', '*3', 1, '*4'], # Break N-H bond
           ['BREAK_BOND', '*3', 1, '*5'], # Break N-H bond
           ['FORM_BOND', '*2', 1, '*4'], # Form O-H
           ['FORM_BOND', '*2', 1, '*5'], # Form O-H (Water formation)
           ['BREAK_BOND', '*1', 1, '*2'], # Cleave C-O bond entirely
           ['CHANGE_BOND', '*1', 1, '*3'], # C-N single becomes double (C=N)
       ])

entry(
    index = 0,
    label = "AliphaticAldehyde",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2 *2 O u0 {1,D}
3 *3 C u0 {1,S} # Aliphatic chain attachment (C4, C6, C9 etc)
4 *4 H u0 {1,S} 
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "PrimaryAmine",
    group = 
"""
1 *3 N u0 {2,S} {3,S}
2 *4 H u0 {1,S}
3 *5 H u0 {1,S}
""",
    kinetics = None,
)

tree(
    {'AliphaticAldehyde': [],
     'PrimaryAmine': []}
)
