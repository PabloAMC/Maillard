#!/usr/bin/env python
# encoding: utf-8

name = "Beta_Elimination"
shortDesc = "Thermal beta-elimination of serine or cysteine to form dehydroalanine (DHA)."
longDesc = """
Pathway E: Dehydroalanine (DHA) Competing Pathway
Under severe thermal load and extrusion shear, serine and cysteine residues 
undergo beta-elimination. This expels water or H2S, respectively, and leaves 
behind highly reactive dehydroalanine (DHA) residues. 
This is the critical pre-requisite step for the massive protein cross-linking 
that starves the Maillard reaction of lysine.
X-CH2-CH(NH2)-COOH => CH2=C(NH2)-COOH + HX  (where X = OH or SH)
"""

recipe(name="Beta_Elimination_Cysteine",
       recipe=[
           ['BREAK_BOND', '*1', 1, '*2'], # C(beta)-S bond breaks
           ['BREAK_BOND', '*3', 1, '*4'], # C(alpha)-H bond breaks
           ['FORM_BOND', '*2', 1, '*4'],  # S picks up H -> H2S
           ['CHANGE_BOND', '*1', 1, '*3'], # C(alpha)-C(beta) becomes double bond (DHA)
       ])

recipe(name="Beta_Elimination_Serine",
       recipe=[
           ['BREAK_BOND', '*1', 1, '*5'], # C(beta)-O bond breaks
           ['BREAK_BOND', '*3', 1, '*4'], # C(alpha)-H bond breaks
           ['FORM_BOND', '*5', 1, '*4'],  # O picks up H -> H2O
           ['CHANGE_BOND', '*1', 1, '*3'], # C(alpha)-C(beta) becomes double bond (DHA)
       ])

entry(
    index = 0,
    label = "Cysteine",
    group = 
"""
1 *1 C u0 {2,S} {3,S} # Beta carbon
2 *2 S u0 {1,S} {6,S}
3 *3 C u0 {1,S} {4,S} {7,S} # Alpha carbon
4 *4 H u0 {3,S} # Alpha proton
6 *6 H u0 {2,S} # Sulfhydryl proton
7 *7 N u0 {3,S} # Amine group attached to alpha
""",
    kinetics = None,
)

entry(
    index = 1,
    label = "Serine",
    group = 
"""
1 *1 C u0 {2,S} {3,S} # Beta carbon
2 *5 O u0 {1,S} {6,S}
3 *3 C u0 {1,S} {4,S} {7,S} # Alpha carbon
4 *4 H u0 {3,S} # Alpha proton
6 *6 H u0 {2,S} # Hydroxyl proton
7 *7 N u0 {3,S} # Amine group attached to alpha
""",
    kinetics = None,
)

tree(
    {'Cysteine': [],
     'Serine': []}
)
