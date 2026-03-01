#!/usr/bin/env python
# encoding: utf-8

name = "Thiol_Addition"
shortDesc = "Kinetics limits for Thiol Addition"
longDesc = "Baseline aqueous thermodynamics for H2S/thiol addition to carbonyls."

entry(
    index = 0,
    label = "Carbonyl;Thiol",
    kinetics = ArrheniusEP(
        A = (1e6, 'cm^3/(mol*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (12.0, 'kcal/mol'), # Thiols are strong nucleophiles, often faster than amines
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Thiol addition rate",
)
