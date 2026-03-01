#!/usr/bin/env python
# encoding: utf-8

name = "Amadori_Rearrangement"
shortDesc = "Isomerisation of an N-substituted aldosylamine to a 1-amino-1-deoxy-2-ketose via a 1,2-enaminol intermediate."
longDesc = """
Maillard Reaction - Phase 1:
The crucial irreversible step turning the Schiff base (glycosylamine) into the Amadori Product.
Requires an alpha-proton on the C2 carbon to shift to the C1 carbon, while the C1=N double bond shifts to a C2=O double bond.
R-CH(OH)-CH=NR' <=> R-C(=O)-CH2-NHR'
"""

recipe(name="Amadori_Rearrangement",
       recipe=[
           ['CHANGE_BOND', '*1', -1, '*2'], # C1=N -> C1-N
           ['CHANGE_BOND', '*3', 1, '*4'], # C2-O -> C2=O
           ['BREAK_BOND', '*4', 1, '*5'], # O-H break (alpha proton)
           ['FORM_BOND', '*1', 1, '*5'], # C1 binds the shifting H
           ['FORM_BOND', '*2', 1, '*6'], # N picks up a proton from solvent (simplified as intramolecular here for RMG matching)
       ])

entry(
    index = 0,
    label = "Aldosylamine",
    group = 
"""
1 *1 C u0 {2,D} {3,S}
2 *2 N u0 {1,D} 
3 *3 C u0 {1,S} {4,S}
4 *4 O u0 {3,S} {5,S}
5 *5 H u0 {4,S}
6 *6 H u0 {2,S} # Assuming secondary amine from amino acid
""",
    kinetics = None,
)

tree(
    {'Aldosylamine': []}
)
