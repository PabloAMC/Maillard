#!/usr/bin/env python
# encoding: utf-8

name = "Lipid_Thiazole_Condensation"
shortDesc = "Kinetics limits for Alkylthiazole formation"
longDesc = "Baseline aqueous thermodynamics allowing RMG to form lipid-derived thiazoles."

entry(
    index = 0,
    label = "AliphaticAldehyde;Ammonia;HydrogenSulfide",
    kinetics = ArrheniusEP(
        A = (1e3, 'cm^6/(mol^2*s)'),
        n = 1.0,
        alpha = 0,
        E0 = (25.0, 'kcal/mol'), 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Alkylthiazole condensation rate",
)
