# Validated Envelope

Target tag: meaty
Supported benchmarks: 14/14
Strict-ready benchmarks: none
Matrix-only executable benchmarks: Pea isolate, 40 C (Pratap Singh, 2021), Pea isolate UHT, 140 C (Trikusuma, 2020), Soy isolate, 40 C (Pratap Singh, 2021)

## Warnings
- Matrix benchmarks are executable intake/headspace checks, but remain outside the strict gate and target snapshots.
- Benchmark-facing concentrations still use the FAST observable projection; Cantera remains diagnostic-reference-only.
- Peptide-bound and intact-protein accessibility remain outside the validated precursor envelope.
- The validated plant-matrix scope is currently limited to pea/soy matrix-only systems and not yet a broad release gate.

## Next Priorities
- Expose matrix-state and projection explainability in user-facing artifacts.
- Promote domain-of-applicability warnings into routine CLI/report outputs.
- Replace bulk matrix retention with broader compound-aware observability across plant-matrix systems.
