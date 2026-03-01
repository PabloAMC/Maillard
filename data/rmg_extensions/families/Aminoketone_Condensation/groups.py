#!/usr/bin/env python
# encoding: utf-8

name = "Aminoketone_Condensation"
shortDesc = "Self-condensation of two alpha-aminoketones to form a dihydropyrazine and water."
longDesc = """
Maillard Reaction - Phase 3 (Pyrazines):
Two alpha-aminoketones (products of Strecker degradation) condense together
to form a cyclic dihydropyrazine, which then oxidises to a pyrazine.
This is the primary route to roasted, nutty flavours (e.g. 2,3-dimethylpyrazine).
2 R-C(=O)-CH(NH2)-R' => Dihydropyrazine + 2 H2O
"""

# Note: In RMG, this is typically modelled as a sequence:
# 1. Amine attacks carbonyl (Schiff base)
# 2. Intramolecular amine attacks second carbonyl (cyclization)
# For this targeted library, we provide a concerted macro-step to ensure pyrazines 
# are robustly generated in the graph without deep combinatorics.

recipe(name="Aminoketone_Condensation",
       recipe=[
           ['BREAK_BOND', '*1', 2, '*2'], # C=O double bond breaks on AK1
           ['BREAK_BOND', '*4', 2, '*5'], # C=O double bond breaks on AK2
           ['BREAK_BOND', '*3', 1, '*7'], # N-H breaks on AK1
           ['BREAK_BOND', '*3', 1, '*8'], # N-H breaks on AK1
           ['BREAK_BOND', '*6', 1, '*9'], # N-H breaks on AK2
           ['BREAK_BOND', '*6', 1, '*10'],# N-H breaks on AK2
           ['FORM_BOND', '*1', 1, '*6'],  # C(AK1)-N(AK2) bond
           ['FORM_BOND', '*4', 1, '*3'],  # C(AK2)-N(AK1) bond
           ['CHANGE_BOND', '*1', 1, '*6'], # C-N to C=N double bond
           ['CHANGE_BOND', '*4', 1, '*3'], # C-N to C=N double bond
           ['BREAK_BOND', '*1', 1, '*2'], # C-O bond fully breaks on AK1
           ['BREAK_BOND', '*4', 1, '*5'], # C-O bond fully breaks on AK2
           ['FORM_BOND', '*2', 1, '*7'],  # O(AK1) picks up H -> H2O
           ['FORM_BOND', '*2', 1, '*8'],  # O(AK1) picks up H -> H2O
           ['FORM_BOND', '*5', 1, '*9'],  # O(AK2) picks up H -> H2O
           ['FORM_BOND', '*5', 1, '*10'], # O(AK2) picks up H -> H2O
       ])

entry(
    index = 0,
    label = "AlphaAminoketone_1",
    group = 
"""
1 *1 C u0 {2,D} {11,S}
2 *2 O u0 {1,D}
3 *3 N u0 {11,S} {7,S} {8,S}
7 *7 H u0 {3,S}
8 *8 H u0 {3,S}
11 *11 C u0 {1,S} {3,S}
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "AlphaAminoketone_2",
    group = 
"""
1 *4 C u0 {5,D} {12,S}
2 *5 O u0 {4,D}
3 *6 N u0 {12,S} {9,S} {10,S}
9 *9 H u0 {6,S}
10 *10 H u0 {6,S}
12 *12 C u0 {4,S} {6,S}
""",
    kinetics = None,
)

tree(
    {'AlphaAminoketone_1': [],
     'AlphaAminoketone_2': []}
)
