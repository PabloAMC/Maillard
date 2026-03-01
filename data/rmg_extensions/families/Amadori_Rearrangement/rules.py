#!/usr/bin/env python
# encoding: utf-8

name = "Amadori_Rearrangement"
shortDesc = "Kinetics limits for Amadori base-catalysed isomerisation"
longDesc = "Baseline aqueous thermodynamics allowing RMG to traverse into the Amadori product."

entry(
    index = 0,
    label = "Aldosylamine",
    kinetics = ArrheniusEP(
        A = (1e6, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (20.0, 'kcal/mol'),
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Amadori rate",
)
