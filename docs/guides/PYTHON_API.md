# Programmatic Usage (Python API)

While Maillard is primarily interacted with via the CLI, you can easily embed it into custom Python notebooks, internal automated LIMS systems, or design loops.

## The Two Core Objects

The framework is driven by two main objects: `ReactionConditions` (defining the environment) and `MaillardPipeline` (defining the formulation and scoring).

### 1. `ReactionConditions`
Defines the physical boundaries of the reaction.

```python
from src.conditions import ReactionConditions

conditions = ReactionConditions(
    pH=6.5, 
    temperature_celsius=140.0,
    water_activity=0.95,
    protein_type="pea_iso"
)
```

### 2. `MaillardPipeline`
Defines what sensory tag you are optimizing for (e.g., `meaty`, `roasted`) and what you want to penalize (e.g., `beany`, `bitter`).

```python
from src.pipeline import MaillardPipeline

designer = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
```

## Running a Single Evaluation

To predict the chemical outputs of a specific formulation, pass a formulation dictionary to `evaluate_single`. 

```python
formulation = {
    "name": "Pea Isolate Test",
    "sugars": ["ribose", "glucose"],
    "amino_acids": ["cysteine", "leucine"],
    "lipids": ["hexanal"], # Pre-existing off-notes to mask
    "additives": [],
    "molar_ratios": {"ribose": 0.5, "glucose": 0.2, "cysteine": 0.2, "leucine": 0.1},
    "ph": conditions.pH,
    "temp": conditions.temperature_celsius,
    "aw": conditions.water_activity,
    "time_minutes": 45.0,
    "protein_type": conditions.protein_type,
    "denaturation_state": 0.6,
    "catalyst": None
}

result = designer.evaluate_single(formulation, conditions)

print(f"Target Score: {result.target_score}")
print("Predicted Targets:")
for target in result.targets:
    tag = target['target'].label if target['target'] else "Unknown"
    print(f"- {tag}: {target.get('sensory', '')}")
```

## Verifying the Scientific Trust Envelope

Never trust a programmatic output blindly. Always run the input combinations through the `DomainOfValidityChecker` to determine if your prediction is backed by literature benchmarks or is merely an uncalibrated directional extrapolation.

```python
from src.usability_reports import DomainOfValidityChecker, build_confidence_package

checker = DomainOfValidityChecker(target_tag="meaty")

precursors = formulation["sugars"] + formulation["amino_acids"] + formulation["lipids"]
warnings = checker.check(
    precursor_names=precursors,
    protein_type=formulation["protein_type"],
    temp_c=formulation["temp"],
    ph=formulation["ph"]
)

if not warnings:
    print("All inputs are within the rigorously validated envelope.")
else:
    for w in warnings:
        print(f"Warning: {w.description}")
```

> **Note:** See `docs/notebooks/1_Formulation_Screening_Example.ipynb` for a fully runnable Jupyter implementation of this pipeline.
