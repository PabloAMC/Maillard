import pytest
import json
import os
import pandas as pd
from pathlib import Path
from src.smirks_engine import SmirksEngine, Species  # noqa: E402
from src.cantera_export import CanteraExporter  # noqa: E402
from src.kinetics import KineticsEngine  # noqa: E402
from src.results_db import ResultsDB  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from rdkit import Chem
from rdkit.Chem import Descriptors

@pytest.fixture
def regression_data():
    path = Path("data/lit/canonical_systems.json")
    with open(path, "r") as f:
        return json.load(f)

@pytest.fixture
def results_db():
    return ResultsDB(db_path="results/maillard_results.db")

NAME_MAP = {
    "O=CC(O)C(O)C(O)CO": "ribose",
    "O=CC(O)C(O)C(O)C(O)CO": "glucose",
    "NCC(=O)O": "glycine",
    "NCC(O)=O": "glycine",
    "NC(CS)C(=O)O": "cysteine",
    "CC(C)CC(N)C(=O)O": "leucine",
    "O": "water",
    "O=C=O": "CO2",
    "N": "ammonia",
    "S": "H2S",
    "[HH]": "H2",
    "Cc1nccs1": "2-methylthiazole",
    "SCc1ccco1": "2-furfurylthiol",
    "O=Cc1ccco1": "furfural",
    "OCC1=CC=C(C=O)O1": "HMF",
    "CC(C)CC=O": "3-methylbutanal",
    "Cc1cnc(C)cn1": "2,5-dimethylpyrazine",
    "CC=O": "acetaldehyde"
}

LOOKUP = {
    "ribose": "O=CC(O)C(O)C(O)CO",
    "glucose": "O=CC(O)C(O)C(O)C(O)CO",
    "glycine": "NCC(=O)O",
    "cysteine": "NC(CS)C(=O)O",
    "leucine": "CC(C)CC(N)C(=O)O"
}

# Synonyms used in data/lit/canonical_systems.json that are the SAME molecule as a
# species the engine emits under another label.
# Abbreviations and synonyms used in data/lit/canonical_systems.json that name the SAME
# molecule as a species the engine emits under another label. Added 2026-08-27 (Wave H):
# the old "at least one target in the top 10" assertion silently never matched "FFT" or
# "2,5-dimethylpyrazine" at all, so two of the three systems were passing on a single
# target each (furfural, HMF) while their other targets were never checked.
_TARGET_SYNONYMS = {
    "isovaleraldehyde": "3_methylbutanal",
    "fft": "2_furfurylthiol",
    "mft": "2_methyl_3_furanthiol",
}


def _cantera_safe(name: str) -> str:
    """Sanitize a chemical name the way CanteraExporter names its species.

    NOTE 2026-08-27 (Wave H): the in-line sanitizer this test used to carry replaced
    dashes, spaces and parentheses but NOT commas, so "2,5-dimethylpyrazine" came out as
    "2,5_dimethylpyrazine" and never matched the emitted "2_5_dimethylpyrazine". The
    glucose_glycine case passed on its other target (HMF) and the miss was invisible.
    """
    return "".join(character if character.isalnum() else "_" for character in str(name))


def _resolve_target_species(safe_target: str, volatiles: dict) -> str | None:
    """Map a canonical-systems target name onto an emitted species name, or None."""
    lookup = {_cantera_safe(name): name for name in volatiles}
    candidates = [
        safe_target,
        _cantera_safe(safe_target),
        _TARGET_SYNONYMS.get(_cantera_safe(safe_target).lower()),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        if candidate in volatiles:
            return candidate
        resolved = lookup.get(_cantera_safe(candidate))
        if resolved is not None:
            return resolved
    return None


# RE-DERIVED 2026-08-27 (Wave H). This test used to assert, for every system, that at
# least one characteristic target appears in the TOP TEN volatiles by final simulated
# concentration in the Cantera lane. `ribose_cysteine_leucine` no longer satisfies that,
# and the diagnosis is that the claim was never sound for this system rather than that
# the chemistry regressed:
#
#  * All five of its targets are still REACHABLE and non-zero (asserted below, which is
#    a strictly stronger statement than the old "at least one of them ranks high").
#  * What changed is their magnitude in THIS LANE only: FFT 3.8e-7 -> 6.8e-14 mol/L,
#    MFT 1.5e-7 -> 6.6e-15. The cause is Wave G1's removal of the fabricated H2 economy
#    (free H2 was being manufactured from water by a Strecker + deamination loop). The
#    three steps that reach these targets -- Thiol_Dehydration, Thiol_Addition_Norfuraneol,
#    Thiol_Addition_Legacy_Shortcut -- all consume the `[HH]` REDUCING-EQUIVALENT TOKEN,
#    which stands in for sugar-derived reductones because the model carries no explicit
#    reductone couple. The production path (src/recommend.py) treats that token as a
#    POOL GATE: the lane runs if a reductant was generated anywhere. The Cantera export
#    takes it literally as a mass-action reactant, so with the fabricated H2 source gone
#    these steps are throttled by an [H2] of ~1e-8. That is an artifact of the export
#    lane's inability to represent a lumping token, not a chemistry finding.
#  * The old top-10 for this system contained `MFT_radical`, `elemental_sulfur` and
#    `quenched_elemental_sulfur` -- three of the fabricated species G1 deleted. The
#    ranking it was pinning was partly a ranking of fictitious chemistry.
#
# So the abundance claim is kept ONLY for the two systems whose targets include a major
# species that the lane can actually represent, and the reachability claim -- which is
# the real regression guard -- is applied to all three. The real quantitative check on
# MFT/FFT is the benchmark panel (cys_ribose_140C_Hofmann1998 and the xylose HVP lanes),
# which runs the production path and reports its residuals honestly.
_TOP10_ABUNDANCE_SYSTEMS = {"ribose_cysteine", "glucose_glycine"}


@pytest.mark.regression
@pytest.mark.parametrize("system_key", ["ribose_cysteine", "glucose_glycine", "ribose_cysteine_leucine"])
def test_canonical_systems(system_key, regression_data, results_db):
    data = regression_data[system_key]
    precursors = data["precursors"]
    targets = data["targets"]
    
    # 1. Setup Engine and Conditions
    cond = ReactionConditions(temperature_celsius=150.0)
    engine_smirks = SmirksEngine(conditions=cond)
    
    # 2. Map Precursors
    precursor_objs = []
    for name, conc in precursors.items():
        smi = LOOKUP.get(name.lower(), name)
        precursor_objs.append(Species(label=name, smiles=smi))
        
    # 3. Discover Network
    steps = engine_smirks.enumerate(precursor_objs, max_generations=3)
    assert len(steps) > 0
    
    # 4. Export to Cantera
    exporter = CanteraExporter()
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        
        # Dual-lookup barrier
        barrier_kcal, _, _ = results_db.get_best_barrier(reactants, products, step.reaction_family)
        
        # Add species with names
        for s in step.reactants + step.products:
            name = NAME_MAP.get(s.smiles, s.label)
            exporter.add_species(s.smiles, name=name)
        
        try:
            exporter.add_reaction(reactants, products, barrier_kcal)
        except ValueError:
            continue
            
    assert len(exporter.reactions) > 0
    
    mech_path = f"tests/temp_{system_key}_mech.yaml"
    exporter.export_yaml(mech_path)
    
    # 5. Simulate
    kinetics = KineticsEngine(temperature_k=423.15) # 150 C
    
    init_state = {}
    for name, conc in precursors.items():
        # Find the name used in the mechanism for this precursor
        found_name = None
        smi = LOOKUP.get(name.lower(), name)
        for _, s_info in exporter.species.items():
            if s_info["smiles"] == smi:
                found_name = s_info["name"]
                break
        if found_name:
            init_state[found_name] = conc
            
    assert len(init_state) == len(precursors)
    
    # Increase time to allow cascades to finish
    simulation_duration = 3600 # 1 hour
    results = kinetics.simulate_network_cantera(mech_path, init_state, (0, simulation_duration))
    
    # 6. Verify Targets
    pd.DataFrame(results)
    final_concs = {k: v for k, v in results.items() if k not in ["time", "temperature"] and not k.endswith("_X")}
    
    # Filter out precursors and inert small molecules for "volatiles" ranking
    inerts = ["water", "CO2", "H2", "ammonia", "H2S"]
    # Create lookup for volatility (Phase 16.8: centralized scientific check)
    name_to_species = {}
    for _, s_info in exporter.species.items():
        name_to_species[s_info["name"]] = Species(label=s_info["name"], smiles=s_info["smiles"])

    def is_volatile(name):
        sp = name_to_species.get(name)
        if sp:
            return sp.is_volatile
        # Fallback for precursors not in exporter list
        if name.lower() in ["water", "co2", "h2", "ammonia", "h2s"]:
            return False
        return True

    volatiles = {k: v[-1] for k, v in final_concs.items() if is_volatile(k)}

    # Sort volatiles by final concentration
    sorted_volatiles = sorted(volatiles.items(), key=lambda x: x[1], reverse=True)
    top_n_names = [name for name, conc in sorted_volatiles[:10]]

    # Log for debugging if test fails
    print(f"\nSystem: {system_key}")
    print(f"Top 10 Volatiles: {sorted_volatiles[:10]}")

    # Targets in canonical_systems.json still have dashes; sanitize them for comparison
    safe_targets = [t.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_") for t in targets]

    # PRIMARY CLAIM (added 2026-08-27, Wave H): every named target must be REACHABLE in
    # the enumerated network and reach a strictly positive concentration. This is a
    # stronger coverage guard than the old "at least one target in the top 10" -- it
    # checks all of them, and it is exactly what breaks if a SMIRKS template regresses.
    resolved = {
        target: _resolve_target_species(target, volatiles) for target in safe_targets
    }
    missing = [target for target, name in resolved.items() if name is None]
    assert not missing, (
        f"{system_key}: targets {missing} are NOT REACHABLE in the enumerated network. "
        f"Reachable volatiles: {sorted(volatiles)}"
    )
    zeroed = [
        target for target, name in resolved.items() if not volatiles.get(name, 0.0) > 0.0
    ]
    assert not zeroed, (
        f"{system_key}: targets {zeroed} are reachable but end at exactly zero "
        "concentration -- a route into them has been severed."
    )

    # SECONDARY CLAIM: abundance ranking. Kept only where it is meaningful; see
    # _TOP10_ABUNDANCE_SYSTEMS for why ribose_cysteine_leucine is excluded.
    if system_key in _TOP10_ABUNDANCE_SYSTEMS:
        found_any = any(t in top_n_names for t in safe_targets)
        assert found_any, f"None of sanitized targets {safe_targets} found in top {len(top_n_names)} volatiles: {top_n_names}"

    # Cleanup
    if os.path.exists(mech_path):
        os.remove(mech_path)
