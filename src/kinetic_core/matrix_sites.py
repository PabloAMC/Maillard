"""
The protein matrix as chemistry: reactive sites charged per gram of protein, and the declared
binding of aldehydes and HMF to them (results/validation/matrix_sites_prereg.md, 2026-09-08).

Two things happen when a spec states a protein loading:

1. The sulfur lane's protein-disulfide pool ``PROT_SS`` is charged (disulfide sites x g/L), so the
   thiol-to-protein exchange channel the engine has carried inert since B2.1 runs with a real pool.
2. Aldehydes and HMF are bound after integration by a pseudo-first-order factor over the thermal
   programme: fraction bound = 1 - exp(-sum_segments k2(T) * [sites] * t), with k2 and its
   temperature dependence DECLARED from the adduct dossiers (brackets, not fits) and the bracket's
   corners priced as an interval width, the way the lipid lane prices its Q10.

No number here is fitted. Every rate is a bracket with its dossier anchor; every site density is a
count from a dossier divided by a molar mass. A matrix with no sites on file charges nothing.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Tuple

from src import data_access, data_paths

R_KJ = 8.314462618e-3
CELSIUS = 273.15

#: The declared binding classes: k2 bracket (M^-1 s^-1) at its reference temperature, the activation
#: energy band (kJ/mol), the site pool it consumes, and the anchor. Brackets and bands are the corpus's,
#: from the adduct-kinetics synthesis (k6b_adduct_kinetics_synthesis.md section 1c and 1a).
BINDING_CLASSES: Dict[str, Dict[str, Any]] = {
    "saturated_aldehyde_amine": {
        "k2_bracket": (6.0e-6, 2.5e-5), "t_ref_c": 20.0, "ea_band_kj_mol": (15.0, 20.0), "pool": "amine",
        "source": "k6b_adduct_kinetics_synthesis.md 1c: hexanal + Na-caseinate/whey lysine <= 2.5e-5 M^-1 s^-1 at 20 C (meynier2004) and "
                  "hexanal + BLG 1e-6 to 1e-5 (anantharamkrishnan2020b); Ea 15-20 kJ/mol from decanal + BLG (shepelev2024, K6b 1a)",
    },
    "unsaturated_aldehyde_amine": {
        "k2_bracket": (5.3e-5, 7.9e-5), "t_ref_c": 20.0, "ea_band_kj_mol": (15.0, 20.0), "pool": "amine",
        "source": "k6b_adduct_kinetics_synthesis.md 1c: trans-2-hexenal + Na-caseinate/whey lysine 5.3-7.9e-5 M^-1 s^-1 at 20 C (meynier2004); "
                  "Ea band as for the saturated aldehydes (no enal Ea is measured)",
    },
    "hmf_thiol": {
        "k2_bracket": (3.95 / 86400.0, 23.3 / 86400.0), "t_ref_c": 25.0, "ea_band_kj_mol": (11.6, 29.6), "pool": "free_thiol",
        "source": "hamzalioglu2018_extraction.md: HMF + cysteine second-order 3.95 / 5.15 / 23.3 M^-1 day^-1 at 5 / 25 / 50 C aqueous pH 3.5, "
                  "Ea 29.6 kJ/mol aqueous and 11.6 in coffee (the bracket spans the three temperatures)",
    },
    "hmf_amine": {
        "k2_bracket": (0.088 / 86400.0 / 0.02, 0.160 / 86400.0 / 0.02), "t_ref_c": 25.0, "ea_band_kj_mol": (10.0, 12.3), "pool": "amine",
        "source": "hamzalioglu2018_extraction.md: HMF + lysine pseudo-first-order 0.088 / 0.090 / 0.160 day^-1 at 5 / 25 / 50 C with 20 mM lysine, "
                  "Ea 10.0 (aqueous) to 12.3 (coffee) kJ/mol",
    },
}

#: Engine species key -> binding classes it undergoes (a compound may bind through more than one).
BINDING_OF_SPECIES: Dict[str, Tuple[str, ...]] = {
    "HEXANAL": ("saturated_aldehyde_amine",),
    "NONANAL": ("saturated_aldehyde_amine",),          # the hexanal bracket, declared: chain-length effect not applied
    "DECADIENAL": ("unsaturated_aldehyde_amine",),
    "ME_13_OXO_TRIDECADIENOATE": ("unsaturated_aldehyde_amine",),
    "ME_9_OXONONANOATE": ("saturated_aldehyde_amine",),
    "FUR": ("saturated_aldehyde_amine",),               # furfural's aldehyde; declared with the hexanal bracket
    "HMF": ("hmf_thiol", "hmf_amine"),
}


@dataclass(frozen=True)
class Sites:
    """Site densities in mmol per gram of protein, with their provenance."""

    free_thiol: float
    disulfide: float
    amine: float
    source: str
    amine_band: Tuple[float, float] = (1.0, 1.0)


@dataclass(frozen=True)
class ChargedSites:
    """The pools in mmol/L for one spec."""

    matrix: str
    protein_g_per_l: float
    free_thiol: float
    disulfide: float
    amine: float
    source: str
    amine_band: Tuple[float, float]

    def as_dict(self) -> Dict[str, Any]:
        return {"matrix": self.matrix, "protein_g_per_l": self.protein_g_per_l,
                "pools_mmol_per_l": {"free_thiol": self.free_thiol, "disulfide": self.disulfide, "amine": self.amine},
                "amine_available_band": list(self.amine_band), "source": self.source}


def matrices() -> Dict[str, Sites]:
    raw = data_access.load_yaml(data_paths.PROTEIN_MATRICES)["matrices"]
    out = {}
    for key, m in raw.items():
        # the densities are the dossier's counts over the molar mass, computed here so the file's
        # rounded convenience values can never drift from the counts (a test asserts they agree)
        counts, mass = m["per_monomer"], float(m["molar_mass_g_per_mol"])
        per_g = lambda n: float(n) / mass * 1000.0  # noqa: E731
        out[key] = Sites(per_g(counts.get("free_cysteine", 0)), per_g(counts.get("disulfide", 0)), per_g(counts.get("lysine", 0)),
                         str(m["source"]), tuple(float(x) for x in m.get("amine_available_band", [1.0, 1.0])))
    return out


class MatrixSpecError(ValueError):
    pass


def resolve(process) -> Tuple[Optional[ChargedSites], Optional[str]]:
    """(the charged pools, or None) and a note for the answer when nothing is charged and why.

    Precedence: a spec's own `protein_sites` (mmol per gram) over the table; either needs
    `protein_g_per_l`. A named matrix without a loading is an error, not a default.
    """
    g_per_l = getattr(process, "protein_g_per_l", None)
    own = getattr(process, "protein_sites", None)
    matrix = str(getattr(process, "matrix", "water") or "water")
    table = matrices()
    if own:
        if g_per_l is None:
            raise MatrixSpecError("protein_sites given without protein_g_per_l: the loading is not defaulted")
        sites = Sites(float(own.get("free_thiol_mmol_per_g", 0.0)), float(own.get("disulfide_mmol_per_g", 0.0)),
                      float(own.get("amine_mmol_per_g", 0.0)), "stated in the spec (protein_sites)")
    elif matrix in table:
        if g_per_l is None:
            raise MatrixSpecError(f"matrix {matrix!r} has site densities on file but the spec states no protein_g_per_l")
        sites = table[matrix]
    else:
        if g_per_l:
            return None, (f"protein_g_per_l stated but matrix {matrix!r} has no site densities on file and the spec "
                          "states no protein_sites: nothing is charged (state free_thiol / disulfide / amine in mmol per gram)")
        return None, None
    g = float(g_per_l)
    if g < 0:
        raise MatrixSpecError("protein_g_per_l must be non-negative")
    return ChargedSites(matrix, g, sites.free_thiol * g, sites.disulfide * g, sites.amine * g, sites.source, sites.amine_band), None


def _k2(cls: Mapping[str, Any], k_ref: float, ea_kj: float, temp_c: float) -> float:
    t_ref = cls["t_ref_c"] + CELSIUS
    return k_ref * math.exp(-ea_kj / R_KJ * (1.0 / (temp_c + CELSIUS) - 1.0 / t_ref))


def bound_fraction(species_key: str, charged: ChargedSites, segments) -> Optional[Dict[str, Any]]:
    """The declared bound fraction of ``species_key`` over the thermal programme (a sequence of
    (minutes, temperature C)), at the bracket's centre and at its corners. ``None`` when the species
    has no binding class or its pool is empty. Pseudo-first-order: the sites are in large excess over
    a trace volatile and are not depleted."""
    classes = BINDING_OF_SPECIES.get(species_key)
    if not classes:
        return None
    exps = {"centre": 0.0, "lo": 0.0, "hi": 0.0}
    used = []
    for name in classes:
        cls = BINDING_CLASSES[name]
        pool_mmol = getattr(charged, cls["pool"])
        if pool_mmol <= 0:
            continue
        pool_m = pool_mmol * 1e-3
        band_lo, band_hi = charged.amine_band if cls["pool"] == "amine" else (1.0, 1.0)
        k_lo, k_hi = cls["k2_bracket"]
        k_c = math.sqrt(k_lo * k_hi)
        ea_lo, ea_hi = cls["ea_band_kj_mol"]
        ea_c = 0.5 * (ea_lo + ea_hi)
        for minutes, temp_c in segments:
            seconds = float(minutes) * 60.0
            exps["centre"] += _k2(cls, k_c, ea_c, temp_c) * pool_m * math.sqrt(band_lo * band_hi) * seconds
            exps["lo"] += min(_k2(cls, k_lo, ea, temp_c) for ea in (ea_lo, ea_hi)) * pool_m * band_lo * seconds
            exps["hi"] += max(_k2(cls, k_hi, ea, temp_c) for ea in (ea_lo, ea_hi)) * pool_m * band_hi * seconds
        used.append(name)
    if not used:
        return None
    frac = {k: 1.0 - math.exp(-v) for k, v in exps.items()}
    remaining = {k: 1.0 - v for k, v in frac.items()}
    # the interval half-width in decades on what REMAINS, from the bracket's corners
    extra = 0.5 * abs(math.log10(max(remaining["lo"], 1e-12) / max(remaining["hi"], 1e-12)))
    return {"classes": used, "bound_fraction": frac["centre"], "bound_fraction_corners": [frac["lo"], frac["hi"]],
            "remaining_fraction": remaining["centre"], "extra_decades": extra,
            "sources": [BINDING_CLASSES[n]["source"] for n in used]}
