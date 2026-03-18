# Glossary for Scientists

This document translates the computational framework's internal terminology into standard food science concepts.

| Term | What it means in practice |
|---|---|
| **Benchmark-validated** (or "Strict-ready") | The model prediction is anchored to a specific, quantified, peer-reviewed literature experiment. These predictions can be treated quantitatively (in ppb) with high confidence. |
| **FAST proxy / FAST mode** | The primary, laptop-speed prediction mode. It uses heuristically calibrated energy barriers based on literature rather than running heavy, hours-long Quantum Mechanics calculations for every molecule. |
| **Validated Envelope** (or "Validated Domain") | The physical boundary (pH, temp, time, ingredients) where the model has literature backing. Going outside this envelope (e.g., extremely high pressure or untested protein types) means the output is a hypothesis, not a proven fact. |
| **Observable Projection** | The conversion layer that takes theoretical chemical yields in a beaker and corrects them for what a sensory panel would actually experience (incorporating pH-dependent headspace release, matrix trapping, and degradation). |
| **Target-ranking Validation** | A benchmark where the model doesn't just predict whether a compound forms, but correctly predicts its rank order (e.g., concluding that hexanal will be far more abundant than an oxidized lipid note, matching the lab data). |
| **Intake Check** | A basic test that just verifies the software runs and produces the right compounds, but doesn't guarantee the absolute predicted concentrations are correct yet. |
| **Matrix Retargeting / Matrix Calibration**| Adjusting theoretical free-system predictions (water/buffer) downward to account for physical trapping inside structured plant proteins like pea or soy isolates. |
