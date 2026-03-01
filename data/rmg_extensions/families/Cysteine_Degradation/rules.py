#!/usr/bin/env python
# encoding: utf-8

name = "Cysteine_Degradation"
shortDesc = "Kinetics limits for Cysteine beta-elimination"
longDesc = "Baseline aqueous thermodynamics for H2S generation from Cys."

entry(
    index = 0,
    label = "Cysteine",
    kinetics = ArrheniusEP(
        A = (1e4, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (28.0, 'kcal/mol'), # Beta-elimination barrier
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Cysteine -> H2S degradation rate",
)
