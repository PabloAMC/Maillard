# Lessons Learned

## 2026-03-11: Architectural Review Corrections
- **Cantera Physics**: Using `ideal-gas` for condensed phases introduces thermodynamic bias in reverse rates (P-V work and entropy errors). Always prefer `IdealSolidSolution` or forward-only kinetics for aquatic/food matrices.
- **Thermo Format**: Cantera NASA7 requires two temperature ranges (14 coeffs). Single-range fitting leads to parse errors or unphysical extrapolation.
- **Complexity Assessment**: Don't assume `max(barriers)` is sufficient for linear cascades; use sequential bottleneck approximation `1/Σ(1/kᵢ)`.
- **Pre-exponential Factors**: Silent fallbacks on `NaN` pre-exponentials mask calibration gaps. Flag sources explicitly (Literature vs. Estimated).
- **State Mutation**: Mutating a shared optimizer grid breaks thread safety and reproducibility. Use stateless evaluation methods.

## 2026-03-11: Cantera Condensed-Phase (R.1) Lessons
- **`ideal-condensed` requires molar volumes**: Every species needs `equation-of-state: {model: constant-volume, molar-volume: X}`. Girolami's atom contribution method provides reasonable estimates for organic intermediates.
- **`IdealSolidSolution` is deprecated**: Use `ideal-condensed` (the current canonical name in Cantera 3.x).
- **Reactor compatibility**: `MoleReactor` (constant-volume) fails with `ideal-condensed` because it tries to set density on an incompressible phase. Use `ConstPressureMoleReactor` instead.
- **Label-based filters are fragile**: The original deamination template used label matching (`"amino-" + "deoxyosone"/"acetone"`) which missed `"amino-pyruvaldehyde"`. SMARTS structural gates are more robust and generic.
- **Test mocking after refactor**: When moving object creation from `__init__` to per-call (thread safety), monkeypatching must target the class method, not an instance attribute.
