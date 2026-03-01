#!/usr/bin/env python
# encoding: utf-8

name = "Strecker_Degradation"
shortDesc = "Kinetics limits for Strecker Degradation"
longDesc = "Baseline aqueous thermodynamics for Strecker aldehyde generation."

entry(
    index = 0,
    label = "AlphaDicarbonyl;AlphaAminoAcid",
    kinetics = ArrheniusEP(
        A = (1e5, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (28.0, 'kcal/mol'), 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Strecker condensation and decarboxylation rate",
)
