# Changelog — 2026-03-17

Resumen de cambios introducidos en la rama `real_tests` (cierre Phase D):

- Headspace:
  - La calibración cuantitativa pea/soy se trasladó al modelo de headspace.
  - `HeadspaceModel.get_matrix_benchmark_headspace_factor()` aplica los factores observables pea/soy y la corrección de pH ya existente.
  - Tests unitarios de headspace actualizados y validados.

- Matrix accessibility:
  - `src/matrix_correction.py` expone ventanas de literatura para `pea_iso` y `soy_iso` en `ACCESSIBILITY_LITERATURE_WINDOWS`.
  - Endpoints de lisina/cisteína ajustados a valores conservadores y propagados a `resolve_matrix_correction()`.
  - La ruta de formulación ahora infiere `denaturation_state` desde temperatura/tiempo/pH cuando el usuario no lo fija explícitamente, y expone `effective_denaturation_state` en los resultados.
  - Tests unitarios/integración para la jerarquía `free > soy_iso > pea_iso` actualizados.

- Benchmarking / Policies:
  - `src/benchmark_validation.py` separa yields base de la observabilidad de headspace; la ruta `matrix_only` usa yields base y deja la observabilidad en `HeadspaceModel`.
  - `thermodynamic_gating_policy` y la auditoría permanecen en su sitio; no se hizo gating benchmark-facing sin materialidad.

- Validación Docker:
  - `scientific-fast` (lane rápida): 91 passed, 3 xfailed.
  - `kinetics-validation` (Cantera): 26 passed.
  - Subsets enfocados para headspace / matrix accessibility: green (25 passed en la pasada focalizada tras la inferencia térmica).

Notas operativas:
- La siguiente prioridad de usabilidad ya no es calibrar más constantes a ciegas, sino exponer mejor el estado efectivo de matriz y el dominio de validez en artefactos/reportes orientados al usuario.
- Se mantiene la política de dejar benchmarks `matrix_only` fuera de la puerta `strict_ready` hasta disponer de calibración absoluta reproducible.

Autores: Cambios aplicados por automatización desde la rama `real_tests` (edición 2026-03-17).
