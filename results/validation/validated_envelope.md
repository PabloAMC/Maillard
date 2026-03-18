# Validated Envelope

Target tag: meaty
Supported benchmarks: 8/8
Strict-ready benchmarks: acrylamide_asparagine_glucose_Parker2012, cys_glucose_150C_Farmer1999, cys_ribose_140C_Hofmann1998, cys_ribose_150C_Mottram1994
Matrix-only executable benchmarks: pea_isolate_40C_PratapSingh2021, soy_isolate_40C_PratapSingh2021

## Warnings
- Matrix benchmarks are executable intake/headspace checks, but remain outside the strict gate and target snapshots.
- Benchmark-facing concentrations still use the FAST observable projection; Cantera remains diagnostic-reference-only.
- Peptide-bound and intact-protein accessibility remain outside the validated precursor envelope.
- The validated plant-matrix scope is currently limited to pea/soy matrix-only systems and not yet a broad release gate.

## Next Priorities
- Expose matrix-state and projection explainability in user-facing artifacts.
- Promote domain-of-applicability warnings into routine CLI/report outputs.
- Replace bulk matrix retention with broader compound-aware observability across plant-matrix systems.
