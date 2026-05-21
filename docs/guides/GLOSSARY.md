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
| **VoI (Value of Information)** | A priority score ranking which wet-lab experiment is most valuable to run next. It targets compounds where the model is highly uncertain (wide prediction bounds) *and* which have a high sensory impact (far above odour detection threshold). |
| **ODT (Odour Detection Threshold)** | The lowest concentration of a volatile compound that can be detected by the human nose (typically measured in micrograms per kilogram of water/matrix). Used to calculate how decision-relevant a compound's yield is. |
| **Order of Magnitude (or Dex)** | Logarithmic steps representing a ten-fold ($10\times$) difference. For example, a change of 1 order of magnitude (or $1\text{ dex}$) is a $10\times$ difference in concentration (e.g., from 10 ppb to 100 ppb). A change of 2 orders of magnitude is a $100\times$ difference (e.g., from 10 ppb to 1000 ppb). |

| **SIDA (Stable Isotope Dilution Assay)** | A highly precise mass spectrometry method using isotopically labeled standards to measure absolute concentration (ppb) in complex plant matrices. This is the preferred method for generating new benchmark data. |
| **p5 / p95 (Confidence Whiskers)** | The statistical bounds of a prediction. There is a 90% probability that the actual concentration will fall between the lower bound ($p_5$) and the upper bound ($p_{95}$). |
| **Melanoidin Trapping & Retention** | Physical or chemical trapping of volatile aroma molecules (especially sulfur-containing thiols) by complex brown polymers (melanoidins) formed during late-stage Maillard browning. |
