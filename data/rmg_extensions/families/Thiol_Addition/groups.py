#!/usr/bin/env python
# encoding: utf-8

name = "Thiol_Addition"
shortDesc = "Nucleophilic addition of a thiol (R-SH or H2S) to a carbonyl group."
longDesc = """
Maillard Reaction - Phase 3 (Sulfur/Meaty):
H2S (from cysteine/thiamine degradation) or thiols add to aldehydes and ketones.
Particularly important for the addition of H2S to furfural to form 2-furfurylthiol (FFT), 
a defining impact odorant in roasted meat and coffee.
R-CHO + H2S <=> R-CH(OH)(SH) -> R-CH(=S) + H2O
"""

recipe(name="Thiol_Addition",
       recipe=[
           ['FORM_BOND', '*1', 1, '*3'], # S attacks Carbonyl Carbon
           ['BREAK_BOND', '*1', 2, '*2'], # C=O double bond breaks to single
           ['BREAK_BOND', '*3', 1, '*4'], # S-H bond breaks
           ['FORM_BOND', '*2', 1, '*4'],  # O picks up the H to become hydroxyl
       ])

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
    label = "Thiol",
    group = 
"""
1 *3 S u0 {2,S}
2 *4 H u0 {1,S}
""",
    kinetics = None,
)

tree(
    {'Carbonyl': [],
     'Thiol': []}
)
