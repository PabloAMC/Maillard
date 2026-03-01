#!/usr/bin/env python
# encoding: utf-8

name = "Retro_Aldol_Fragmentation"
shortDesc = "Kinetics limits for Aldol cleavage"
longDesc = "Baseline aqueous thermodynamics for Retro-Aldol fragmentation."

entry(
    index = 0,
    label = "Aldol_Adduct",
    kinetics = ArrheniusEP(
        A = (1e5, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (28.0, 'kcal/mol'), # Generally higher barrier than enolisation unless highly basic
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Retro-Aldol rate",
)
