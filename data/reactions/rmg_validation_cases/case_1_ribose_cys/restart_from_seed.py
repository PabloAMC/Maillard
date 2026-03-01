restartFromSeed(path='seed')

database(
    thermoLibraries=['primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    # Using families relevant for generic aqueous phase and condensation
    kineticsFamilies=['Schiff_Base_Formation', 'Amadori_Rearrangement', 'Heyns_Rearrangement', 'Enolisation', 'Retro_Aldol_Fragmentation', 'Thiol_Addition', 'Strecker_Degradation', 'Aminoketone_Condensation', 'Cysteine_Degradation', 'Lipid_Schiff_Base', 'Lipid_Thiazole_Condensation', 'Beta_Elimination', 'DHA_Crosslinking', 'Sugar_Ring_Opening'],
    kineticsEstimator='rate rules',
)

species(
    label='ribose',
    reactive=True,
    structure=SMILES('O=CC(O)C(O)C(O)CO')
)

species(
    label='cysteine',
    reactive=True,
    structure=SMILES('N[C@@H](CS)C(=O)O')
)

species(
    label='water',
    reactive=True,
    structure=SMILES('O')
)

# Liquid-phase reactor (150°C, typical Maillard cooking/extrusion temperature)
liquidReactor(
    temperature=(423.15,'K'),
    initialConcentrations={'ribose': (1.0, 'mol/l'), 'cysteine': (1.0, 'mol/l'), 'water': (55.0, 'mol/l')},
    terminationConversion={'ribose': 0.9},
    terminationTime=(10,'hr'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)

options(
    saveEdgeSpecies=False,
    saveSimulationProfiles=True,
    generatePlots=True,
)
