#!/usr/bin/env python
# encoding: utf-8

name = "Retro_Aldol_Fragmentation"
shortDesc = "Base-catalysed retro-aldol cleavage of sugars and Amadori products."
longDesc = """
Maillard Reaction - Phase 2 (Fragmentation):
C-C bond cleavage generating small reactive carbonyls (e.g. glycolaldehyde, pyruvaldehyde).
Crucial for generating the C2/C3 building blocks needed for meaty and roasted flavours (thiazoles/pyrazines).
R-CH(OH)-C(=O)-CH2(OH) => R-CHO + CH3-C(=O)-CH2(OH)
"""

recipe(name="Retro_Aldol_Fragmentation",
       recipe=[
           ['BREAK_BOND', '*1', 1, '*2'], # C-C bond breaking
           ['CHANGE_BOND', '*2', 1, '*3'], # C-O single becomes double
           ['BREAK_BOND', '*3', 1, '*4'], # O-H bond breaks
           ['FORM_BOND', '*1', 1, '*4'], # C absorbs proton (enol tautomerization simplified)
       ])

entry(
    index = 0,
    label = "Aldol_Adduct",
    group = 
"""
1 *1 C u0 {2,S}
2 *2 C u0 {1,S} {3,S}
3 *3 O u0 {2,S} {4,S}
4 *4 H u0 {3,S}
""",
    kinetics = None,
)

tree(
    {'Aldol_Adduct': []}
)
