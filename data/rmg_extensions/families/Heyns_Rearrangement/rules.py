#!/usr/bin/env python
# encoding: utf-8

name = "Heyns_Rearrangement"
shortDesc = "Kinetics limits for Heyns base-catalysed isomerisation"
longDesc = "Baseline aqueous thermodynamics allowing RMG to traverse into the Heyns product."

entry(
    index = 0,
    label = "Ketosylamine",
    kinetics = ArrheniusEP(
        A = (1e6, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (22.0, 'kcal/mol'), # Heyns is often slightly slower than Amadori 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic Heyns rate",
)
