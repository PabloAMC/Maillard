#!/usr/bin/env python
# encoding: utf-8

name = "Heyns_Rearrangement"
shortDesc = "Isomerisation of a ketosylamine to a 2-amino-2-deoxyaldose."
longDesc = """
Maillard Reaction - Phase 1:
The ketose equivalent of the Amadori rearrangement (e.g. Fructose + Amine).
Requires a proton shift from C1 to C2, moving the C2=N double bond to a C1=O double bond.
R-C(=NR')-CH2(OH) <=> R-CH(NHR')-CH=O
"""

recipe(name="Heyns_Rearrangement",
       recipe=[
           ['CHANGE_BOND', '*2', -1, '*1'], # C2=N -> C2-N
           ['CHANGE_BOND', '*3', 1, '*4'], # C1-O -> C1=O
           ['BREAK_BOND', '*4', 1, '*5'], # O-H break at C1
           ['FORM_BOND', '*2', 1, '*5'], # C2 binds shifting H
           ['FORM_BOND', '*1', 1, '*6'], # N picks up a proton (intramolecular shortcut for RMG graph logic here)
       ])

entry(
    index = 0,
    label = "Ketosylamine",
    group = 
"""
1 *1 N u0 {2,D} 
2 *2 C u0 {1,D} {3,S}
3 *3 C u0 {2,S} {4,S}
4 *4 O u0 {3,S} {5,S}
5 *5 H u0 {4,S}
6 *6 H u0 {1,S} # Amine secondary proton
""",
    kinetics = None,
)

tree(
    {'Ketosylamine': []}
)
