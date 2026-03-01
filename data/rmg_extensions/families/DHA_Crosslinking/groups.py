#!/usr/bin/env python
# encoding: utf-8

name = "DHA_Crosslinking"
shortDesc = "Nucleophilic addition of a lysine epsilon-amino group to a dehydroalanine (DHA) residue."
longDesc = """
Pathway E: Dehydroalanine (DHA) Competing Pathway
The highly reactive double bond of DHA acts as a powerful electrophile. 
It rapidly undergoes Michael addition with the epsilon-amino group of a free lysine 
residue (or the sulfhydryl of cysteine). 
This crosslinking reaction forms lysinoalanine (LAL) or lanthionine (LAN), building 
textural firmness but permanently consuming the lysine that was needed for the 
Maillard Strecker and pyrazine pathways.
CH2=C(NH2)-COOH + R-NH2 => R-NH-CH2-CH(NH2)-COOH
"""

recipe(name="DHA_Lysinoalanine_Formation",
       recipe=[
           ['FORM_BOND', '*1', 1, '*4'],  # C(beta of DHA)-N(epsilon of Lys) bond formation
           ['CHANGE_BOND', '*1', -1, '*2'], # C(beta)=C(alpha) double bond becomes single
           ['BREAK_BOND', '*4', 1, '*5'],  # N(Lys)-H bond breaks
           ['FORM_BOND', '*2', 1, '*5'],   # C(alpha of DHA) picks up the proton
       ])

entry(
    index = 0,
    label = "Dehydroalanine",
    group = 
"""
1 *1 C u0 {2,D} {6,S} {7,S} # Beta carbon
2 *2 C u0 {1,D} {3,S} {8,S} # Alpha carbon
3 *3 N u0 {2,S} # Alpha amine
6 *6 H u0 {1,S}
7 *7 H u0 {1,S}
8 *8 C u0 {2,S} # Carboxyl attachment
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "LysineEpsilonAmine",
    group = 
"""
1 *4 N u0 {2,S} {3,S} {9,S} # Epsilon amine
2 *5 H u0 {1,S}
3 *10 H u0 {1,S}
9 *9 C u0 {1,S} # Alkyl chain 
10 *11 H u0 {3,S}
""",
    kinetics = None,
)

tree(
    {'Dehydroalanine': [],
     'LysineEpsilonAmine': []}
)
