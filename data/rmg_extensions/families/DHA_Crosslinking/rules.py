#!/usr/bin/env python
# encoding: utf-8

name = "DHA_Crosslinking"
shortDesc = "Kinetics limits for Lysinoalanine/Lanthionine crosslinking"
longDesc = "Baseline aqueous thermodynamics for DHA Michael additions."

entry(
    index = 0,
    label = "Dehydroalanine;LysineEpsilonAmine",
    kinetics = ArrheniusEP(
        A = (1e6, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (16.0, 'kcal/mol'), # Michael additions to DHA are extremely fast/favorable
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic DHA-Lys condensation rate",
)
