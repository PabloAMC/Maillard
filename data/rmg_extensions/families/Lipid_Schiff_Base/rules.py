#!/usr/bin/env python
# encoding: utf-8

name = "Lipid_Schiff_Base"
shortDesc = "Kinetics limits for Lipid Aldehyde Trapping"
longDesc = "Baseline aqueous thermodynamics allowing RMG to model hexanal suppression."

entry(
    index = 0,
    label = "AliphaticAldehyde;PrimaryAmine",
    kinetics = ArrheniusEP(
        A = (1e5, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (14.0, 'kcal/mol'), # Highly reactive electrophile
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Schiff Base trapping rate",
)
