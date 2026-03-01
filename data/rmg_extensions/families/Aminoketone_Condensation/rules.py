#!/usr/bin/env python
# encoding: utf-8

name = "Aminoketone_Condensation"
shortDesc = "Kinetics limits for Pyrazine formation"
longDesc = "Baseline aqueous thermodynamics for dihydropyrazine ring closure."

entry(
    index = 0,
    label = "AlphaAminoketone_1;AlphaAminoketone_2",
    kinetics = ArrheniusEP(
        A = (1e4, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (20.0, 'kcal/mol'), 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Pyrazine formation rate",
)
