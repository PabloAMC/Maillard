#!/usr/bin/env python
# encoding: utf-8

name = "Schiff_Base_Formation"
shortDesc = "Kinetics limits for the condensation of amines with sugars"
longDesc = """
A placeholder for reaction rates of Maillard initial condensation.
Rate is heavily dependent on pH, but RMG-native kinetics don't cleanly 
support external pH scaling out of the box, so we provide baseline aqueous 
favourable kinetics to ensure the pathway is generated and let Tier 1/2 score it.
"""

# Very generic baseline generic rate 
entry(
    index = 0,
    label = "Carbonyl;PrimaryAmine",
    kinetics = ArrheniusEP(
        A = (1e5, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (15.0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Baseline estimated rate for sugar-amine condensation",
)
