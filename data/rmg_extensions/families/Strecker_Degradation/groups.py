#!/usr/bin/env python
# encoding: utf-8

name = "Strecker_Degradation"
shortDesc = "Reaction of an alpha-dicarbonyl with an amino acid to form an aldehyde and an aminoketone."
longDesc = """
Maillard Reaction - Phase 3 (Advanced):
Crucial flavor-generating step. An alpha-dicarbonyl (from sugar fragmentation/enolisation) 
reacts with a free amino acid. The amino acid undergoes decarboxylation (-CO2) to yield the 
Strecker aldehyde (e.g. Valine -> 2-methylpropanal) and an alpha-aminoketone which goes on
to form pyrazines.
R-C(=O)-C(=O)-R' + R''-CH(NH2)-COOH => R''-CHO + R-C(=O)-CH(NH2)-R' + CO2
"""

recipe(name="Strecker_Degradation",
       recipe=[
           ['BREAK_BOND', '*2', 1, '*4'], # Breakdown of alpha-amino acid C-C(O)OH bond
           ['CHANGE_BOND', '*4', 1, '*5'], # Formation of CO2 (C=O)
           ['BREAK_BOND', '*4', 1, '*6'], # O-H cleavage on carboxyl 
           ['FORM_BOND', '*1', 1, '*3'],  # C-N bond formation between dicarbonyl and amine
           ['BREAK_BOND', '*3', 1, '*7'], # N-H break on amine
           ['FORM_BOND', '*6', 1, '*7'],  # Simplified proton transfer (water/solvent mediated)
       ])

entry(
    index = 0,
    label = "AlphaDicarbonyl",
    group = 
"""
1 *1 C u0 {2,S} {3,D}
2 *8 C u0 {1,S} {4,D}
3 *3 O u0 {1,D}
4 *4 O u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "AlphaAminoAcid",
    group = 
"""
1 *2 C u0 {2,S} {3,S} {4,S}
2 *3 N u0 {1,S} {7,S}
3 *4 C u0 {1,S} {5,D} {6,S} # Carboxyl Carbon
4 *5 O u0 {3,D}
5 *6 O u0 {3,S} {8,S}
6 *7 H u0 {2,S} # Amine proton
7 *8 H u0 {5,S} # Carboxyl proton
""",
    kinetics = None,
)

tree(
    {'AlphaAminoAcid': [],
     'AlphaDicarbonyl': []}
)
