#!/usr/bin/env python
# encoding: utf-8

name = "Enolisation"
shortDesc = "Dehydration of Amadori/Heyns products via 1,2- or 2,3-enolisation."
longDesc = """
Maillard Reaction - Phase 2 (Intermediate):
Crucial branching point in the Maillard network:
- 1,2-Enolisation (low pH favoured) leads to 3-deoxyosones and furfurals.
- 2,3-Enolisation (high pH favoured) leads to 1-deoxyosones and pyruvaldehyde/diacetyl.
RMG will currently generate both; our Tier 1 XTBScreener/Ranker will apply the pH modifiers.
"""

# 1,2-Enolisation (loss of water from C3)
recipe(name="1_2_Enolisation",
       recipe=[
           ['BREAK_BOND', '*2', 1, '*4'], # C3-OH break
           ['BREAK_BOND', '*1', 1, '*3'], # C2-H break
           ['FORM_BOND', '*3', 1, '*4'],  # Form Water
           ['CHANGE_BOND', '*1', 1, '*2'], # C2-C3 single becomes double (enol)
       ])

# 2,3-Enolisation (amine elimination from C1)
recipe(name="2_3_Enolisation",
       recipe=[
           ['BREAK_BOND', '*1', 1, '*6'], # C1-N break
           ['BREAK_BOND', '*2', 1, '*3'], # C2-H break 
           ['FORM_BOND', '*3', 1, '*6'],  # Amine leaves carrying the proton
           ['CHANGE_BOND', '*1', 1, '*2'], # C1-C2 single becomes double (enol)
       ])

entry(
    index = 0,
    label = "Amadori_C1_C2_C3",
    group = 
"""
1 *1 C u0 {2,S} {5,D} # C2 (carbonyl)
2 *2 C u0 {1,S} {3,S} # C3
3 *3 H u0 {2,S}
4 *4 O u0 {2,S} {7,S} # OH on C3
5 *5 O u0 {1,D}       # C2=O
6 *6 N u0 {1,S}       # Amine on C1 (for 2,3 enolisation)
7 *7 H u0 {4,S}
""",
    kinetics = None,
)

tree(
    {'Amadori_C1_C2_C3': []}
)
