#!/usr/bin/env python
# encoding: utf-8

name = "Schiff_Base_Formation"
shortDesc = "Nucleophilic addition of a primary amine to a carbonyl group (aldehyde/ketone) followed by dehydration to form a Schiff base (imine) and water."
longDesc = """
Maillard Reaction - Phase 1: 
Condensation of reducing sugars (in their open-chain carbonyl form) with 
the free amino groups of amino acids, peptides, or proteins.
Reaction: R(C=O)R' + R''NH2 <=> R(C=NR'')R' + H2O
"""

recipe(name="Schiff_Base_Formation",
       recipe=[
           ['FORM_BOND', '*1', 1, '*3'], # C-N bond formation
           ['BREAK_BOND', '*1', 2, '*2'], # Break C=O double bond (becomes single)
           ['BREAK_BOND', '*3', 1, '*4'], # Break N-H bond
           ['BREAK_BOND', '*3', 1, '*5'], # Break N-H bond
           ['FORM_BOND', '*2', 1, '*4'], # Form O-H
           ['FORM_BOND', '*2', 1, '*5'], # Form O-H (Water formation)
           ['BREAK_BOND', '*1', 1, '*2'], # Cleave C-O bond entirely
           ['CHANGE_BOND', '*1', 1, '*3'], # C-N single becomes double (C=N)
       ])

# Define the fundamental reactant templates
entry(
    index = 0,
    label = "Carbonyl",
    group = 
"""
1 *1 C u0 {2,D}
2 *2 O u0 {1,D}
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
    {'Carbonyl': [],
     'PrimaryAmine': []}
)
