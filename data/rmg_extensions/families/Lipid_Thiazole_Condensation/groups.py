#!/usr/bin/env python
# encoding: utf-8

name = "Lipid_Thiazole_Condensation"
shortDesc = "Condensation of an aliphatic lipid aldehyde with H2S and ammonia to form an alkylthiazole."
longDesc = """
Pathway D: Lipid-Maillard Crosstalk (Synergy)
While lipid oxidation aldehydes (like hexanal) are typically off-flavors, 
when they encounter the abundant H2S and NH3 generated from cysteine degradation, 
they undergo complex condensation to form long-chain alkylthiazoles.
These compounds are critical for providing species-specific "fatty meat" 
background notes (e.g. roasted chicken or beef fat).
"""

recipe(name="Lipid_Thiazole_Condensation",
       recipe=[
           # A highly simplified concerted macro-step for RMG network generation
           ['FORM_BOND', '*1', 1, '*4'], # C(aldehyde)-S
           ['FORM_BOND', '*4', 1, '*2'], # S-C(aldol fragment)
           ['FORM_BOND', '*2', 1, '*5'], # C-N
           ['FORM_BOND', '*5', 1, '*3'], # N-C
           ['FORM_BOND', '*3', 1, '*1'], # C-C(aldehyde) closing the 5-membered ring
           ['CHANGE_BOND', '*1', 1, '*3'], # Forming double bonds in thiazole ring
           ['CHANGE_BOND', '*2', 1, '*5'],
           # Elimination of waters 
           ['BREAK_BOND', '*1', 2, '*6'], # C=O breaks
           ['BREAK_BOND', '*3', 2, '*7'], # C=O breaks
           ['BREAK_BOND', '*4', 1, '*8'], # S-H breaks
           ['BREAK_BOND', '*5', 1, '*9'], # N-H breaks
       ])

entry(
    index = 0,
    label = "AliphaticAldehyde",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {10,S}
2 *6 O u0 {1,D}
3 *3 C u0 {1,S} {7,D} # Simplified alpha-dicarbonyl intermediate from lipid/sugar 
7 *7 O u0 {3,D}
10 *10 C u0 {1,S} # Alkyl chain tail from the lipid
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "Ammonia",
    group = 
"""
1 *5 N u0 {2,S} {3,S} {4,S}
2 *9 H u0 {1,S}
3 *11 H u0 {1,S}
4 *12 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "HydrogenSulfide",
    group = 
"""
1 *4 S u0 {2,S} {3,S}
2 *8 H u0 {1,S}
3 *13 H u0 {1,S}
""",
    kinetics = None,
)

tree(
    {'AliphaticAldehyde': [],
     'Ammonia': [],
     'HydrogenSulfide': []}
)
