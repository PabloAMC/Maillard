"""
The hypothesis layer: what the literature's reaction rules PROPOSE over the engine's species, placed
against what the engine MODELS. It produces steps and products, never rates or concentrations.

Firewall: src/kinetic_core never imports this package (tests/unit/test_network_hypotheses.py). The
layer reads the engine's species and reaction lists; nothing flows back.
"""
