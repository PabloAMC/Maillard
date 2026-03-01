#!/usr/bin/env python
# encoding: utf-8

name = "Beta_Elimination"
shortDesc = "Kinetics limits for Dehydroalanine generation"
longDesc = "Baseline thermodynamics for high-heat DHA formation from Ser/Cys."

entry(
    index = 0,
    label = "Cysteine",
    kinetics = ArrheniusEP(
        A = (1e5, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (32.0, 'kcal/mol'), # Requires significant thermal/shear stress 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Cys beta-elimination rate",
)

entry(
    index = 1,
    label = "Serine",
    kinetics = ArrheniusEP(
        A = (1e5, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (35.0, 'kcal/mol'), # Typically slightly harder to eliminate OH than SH 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Ser beta-elimination rate",
)
