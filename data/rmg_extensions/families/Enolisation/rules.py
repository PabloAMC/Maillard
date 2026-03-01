#!/usr/bin/env python
# encoding: utf-8

name = "Enolisation"
shortDesc = "Kinetics limits for enolisation branching"
longDesc = "Baseline aqueous thermodynamics allowing RMG to traverse the enolisation forks."

# 1,2-Enolisation baseline
entry(
    index = 0,
    label = "Amadori_C1_C2_C3",
    kinetics = ArrheniusEP(
        A = (1e5, 's^-1'),
        n = 1.0,
        alpha = 0,
        E0 = (25.0, 'kcal/mol'), 
        Tmin = (300, 'K'),
        Tmax = (600, 'K'),
    ),
    rank = 0,
    shortDesc = "Generic 1,2-enolisation rate (acid favoured)",
)
