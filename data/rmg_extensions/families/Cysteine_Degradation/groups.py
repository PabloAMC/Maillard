#!/usr/bin/env python
# encoding: utf-8

name = "Cysteine_Degradation"
shortDesc = "Thermal degradation of cysteine to yield hydrogen sulfide (H2S), ammonia, and pyruvaldehyde."
longDesc = """
Maillard Reaction - Phase 3 (Sulfur Pathway Starter):
Cysteine is the primary source of sulfur in meat flavors. Under thermal stress, 
it undergoes beta-elimination to form H2S, NH3, and a reactive carbonyl (pyruvaldehyde).
The H2S then proceeds to add to furfurals (via Thiol_Addition) to make meaty/roasted thiols.
HS-CH2-CH(NH2)-COOH => H2S + NH3 + CH3-C(=O)-CHO (simplified macro-step for the network generator)
"""

recipe(name="Cysteine_Degradation",
       recipe=[
           ['BREAK_BOND', '*1', 1, '*2'], # C-S bond breaks
           ['BREAK_BOND', '*3', 1, '*4'], # C-N bond breaks
           ['BREAK_BOND', '*5', 1, '*6'], # C-C(carboxyl) bond breaks (decarboxylation)
           ['FORM_BOND', '*2', 1, '*7'],  # S picks up a proton -> H2S (from solvent/intramolecular)
           ['FORM_BOND', '*4', 1, '*8'],  # N picks up a proton -> NH3
           ['CHANGE_BOND', '*5', 1, '*9'], # C-O to C=O (forming pyruvaldehyde backbone)
       ])

entry(
    index = 0,
    label = "Cysteine",
    group = 
"""
1 *1 C u0 {2,S} {3,S} 
2 *2 S u0 {1,S} {7,S} 
3 *3 C u0 {1,S} {4,S} {5,S} # Alpha carbon
4 *4 N u0 {3,S} {8,S} {10,S}
5 *5 C u0 {3,S} {6,S} {9,D} # Carboxyl carbon
6 *6 O u0 {5,S} 
7 *7 H u0 {2,S}
8 *8 H u0 {4,S}
9 *9 O u0 {5,D}
10 *10 H u0 {4,S}
""",
    kinetics = None,
)

tree(
    {'Cysteine': []}
)
