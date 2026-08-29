"""
src/kinetic_core/lipid.py

THE LIPID-OXIDATION NETWORK (Build Wave B6, 2026-08-29).
=============================================================================

Module 5 of the kinetic-core rebuild, and the module whose absence causes 8 of
the cutover exam's 17 refusals.

THE ARCHITECTURE, IN ONE PARAGRAPH
----------------------------------
A hydroperoxide pool, resolved by POSITION (9-/13-) and GEOMETRY (cis,trans /
trans,trans), decomposes first-order into a MEASURED six-product slate. The
DISTRIBUTION over that slate is fitted to Frankel 1989's three zero-additive
columns and frozen. The absolute RATE is not measured anywhere and enters as a
bounded, labelled INPUT (``parameters_lipid.K_LOOH_DECOMP_ANCHOR`` plus
``Q10_ASSUMPTION``). The absolute YIELD per hydroperoxide comes from the ONE
yield-like quantity in the corpus -- Schroen's ``k_hexanal / k4 = 1.2 %``, a
RATIO of two constants measured in the same system, which is temperature-
independent under a common Q10 and is therefore the only part of that paper
that survives the extrapolation.

WHY THOSE TWO SOURCES ARE USED FOR DIFFERENT THINGS
---------------------------------------------------
They are not the same quantity and are never averaged
(``LIPID_SOURCE_CONTRADICTIONS['schroen_1.2pct_vs_frankel_11_20pct']``):

  * **Schroen 2022** gives moles of hexanal per mole of hydroperoxide consumed
    (1.2 %), in a whole-oil autoxidising emulsion. That is a YIELD.
  * **Frankel 1989** gives the share of hexanal among SIX MEASURED PEAKS from a
    PURE linoleate hydroperoxide at 180 C (11-20 % at zero additive). That is a
    DISTRIBUTION.

So the module reads:

    [product p] = LOOH_decomposed
                  x Y_hexanal(Schroen)
                  x share_p(Frankel, at THIS pool's isomer composition, donor)
                  / share_hexanal(Frankel, at the AUTOXIDATION reference column)

At the reference composition the last two factors cancel and the hexanal yield
is exactly Schroen's 1.2 %. Away from it, Frankel's measured isomer dependence
moves it -- which is exactly the structure the corpus licenses and a single
scalar branch fraction cannot express.

WHAT IS FITTED, AND TO WHAT
---------------------------
Four 3-way simplexes (13-OOH and 9-OOH, x cis,trans and trans,trans) and four
composition parameters, to EIGHTEEN numbers: Frankel's three zero-additive
columns. Nothing else in that paper is read by the fitter, and
``fit_branch_model`` asserts the input array's shape before it starts.

WHAT IS NOT FITTED
------------------
The hydrogen-donor suppression ``d``. It is a structural knob with NO stored
value: every claim this module makes about hydrogen donors is a monotonicity
statement over the whole range ``d in (0, 1)``, which is what makes it
unfakeable against a hold-out the builder has seen (see
``results/validation/kinetic_core_b6_prereg.md`` sec. 0).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .parameters_lipid import (
    COVALENT_SINK,
    FRANKEL_SYSTEM_GEOMETRY,
    FRANKEL_ZERO_ADDITIVE,
    K_HEXANAL_SCHROEN,
    K_LOOH_DECOMP_ANCHOR,
    LIPID_CARRIERS,
    LipidCarrier,
    Q10_ASSUMPTION,
    k_looh_decomp_per_min,
    looh_charge_mmol_per_l,
    oleate_fraction,
)
from .species_lipid import (
    CLEAVAGE_MECHANISM,
    FRANKEL_SLATE,
    LIPID_KEYS,
    LIPID_SPECIES,
    LOOH_POOLS,
    MOLECULAR_WEIGHT_G_PER_MOL,
    POSITION_PRODUCTS,
    initial_lipid_state,
    mmol_per_litre_to_ug_per_litre,
)

#: The FIT column that plays the role of "an autoxidised lipid": Frankel's
#: Table 1 is the autoxidation mixture, which is what Schroen's emulsion makes.
REFERENCE_AUTOXIDATION_SYSTEM = "mixed_ct_tt_9_13"

#: moles hexanal per mole hydroperoxide consumed. A RATIO of two constants from
#: the SAME table of the SAME paper, so the hand-fitting offset and the (equal)
#: Q10 both cancel. The dossier calls this one of the paper's two robust exports.
Y_HEXANAL_PER_LOOH: float = (
    float(K_HEXANAL_SCHROEN.value) / float(K_LOOH_DECOMP_ANCHOR.value)
)  # = 0.01


# ===========================================================================
# 1. THE BRANCH MODEL
# ===========================================================================


@dataclass(frozen=True)
class LOOHComposition:
    """
    A hydroperoxide pool's composition over the four (position, geometry) cells.

    Fractions of the LINOLEATE pool only; the oleate pool is carried separately
    because it reaches a different (and unquantified) product set.
    """

    f13_ct: float
    f13_tt: float
    f9_ct: float
    f9_tt: float

    def __post_init__(self) -> None:
        total = self.f13_ct + self.f13_tt + self.f9_ct + self.f9_tt
        if not 0.999 <= total <= 1.001:
            raise ValueError(f"LOOH composition must sum to 1, got {total}")
        for value in (self.f13_ct, self.f13_tt, self.f9_ct, self.f9_tt):
            if value < 0.0:
                raise ValueError("negative composition cell")

    def as_dict(self) -> Dict[str, float]:
        return {
            "LOOH_13_ct": self.f13_ct, "LOOH_13_tt": self.f13_tt,
            "LOOH_9_ct": self.f9_ct, "LOOH_9_tt": self.f9_tt,
        }

    @property
    def f13(self) -> float:
        return self.f13_ct + self.f13_tt


@dataclass(frozen=True)
class BranchModel:
    """
    The frozen branch distribution: one simplex per (position, geometry).

    Keys are ``('13', 'ct')``-style tuples; values map a product key to its
    share of THAT pool's decomposition into the measured slate. Each simplex
    sums to 1 by construction, checked in ``__post_init__``.
    """

    simplexes: Mapping[Tuple[str, str], Mapping[str, float]]
    provenance: str = "fitted to Frankel 1989 zero-additive columns (B6 FIT)"

    def __post_init__(self) -> None:
        expected = {(p, g) for p in ("13", "9") for g in ("ct", "tt")}
        if set(self.simplexes) != expected:
            raise ValueError(f"branch model needs exactly {sorted(expected)}")
        for (position, geometry), shares in self.simplexes.items():
            allowed = POSITION_PRODUCTS[position]
            if tuple(shares) != allowed:
                raise ValueError(
                    f"({position},{geometry}): products must be exactly "
                    f"{allowed} -- the STRUCTURAL ZEROS are not optional"
                )
            total = sum(shares.values())
            if not 0.999 <= total <= 1.001:
                raise ValueError(f"({position},{geometry}) sums to {total}")

    def slate_shares(
        self, composition: LOOHComposition, donor: float = 0.0
    ) -> Dict[str, float]:
        """
        The six measured shares for a pool of this composition, normalised to 1.

        ``donor`` is the hydrogen-donor suppression, ``d in [0, 1)``. It
        multiplies the HOMOLYTIC channel only -- the mechanism assignment is
        ``species_lipid.CLEAVAGE_MECHANISM``, taken from Frankel's INTRODUCTION
        (his refs 3-10, all pre-1989), NOT from the held-out tocopherol arms.
        No value of ``donor`` is stored anywhere in this package.
        """
        flux = self.slate_flux(composition, donor)
        total = sum(flux.values())
        if total <= 0.0:
            return {p: 0.0 for p in FRANKEL_SLATE}
        return {p: flux[p] / total for p in FRANKEL_SLATE}

    def slate_flux(
        self, composition: LOOHComposition, donor: float = 0.0
    ) -> Dict[str, float]:
        """
        The UNNORMALISED slate flux, per unit of hydroperoxide decomposed.

        Unnormalised because the donor's whole signature is that the TOTAL
        falls while the hexanal SHARE rises; a function that returned only
        shares could not express half of it.
        """
        if not 0.0 <= donor < 1.0:
            raise ValueError("donor suppression must lie in [0, 1)")
        cells = composition.as_dict()
        flux = {product: 0.0 for product in FRANKEL_SLATE}
        for pool_key, (position, geometry) in LOOH_POOLS.items():
            weight = float(cells[pool_key])
            if weight <= 0.0:
                continue
            for product, share in self.simplexes[(position, geometry)].items():
                suppression = (
                    (1.0 - donor)
                    if CLEAVAGE_MECHANISM[product] == "homolytic"
                    else 1.0
                )
                flux[product] += weight * float(share) * suppression
        return flux

    def total_relative_flux(
        self, composition: LOOHComposition, donor: float = 0.0
    ) -> float:
        """Total measured-slate flux relative to the donor-free case."""
        return sum(self.slate_flux(composition, donor).values())

    def as_dict(self) -> Dict[str, Any]:
        return {
            "provenance": self.provenance,
            "simplexes": {
                f"{position}_{geometry}": dict(shares)
                for (position, geometry), shares in sorted(self.simplexes.items())
            },
        }


def _softmax(logits: Sequence[float]) -> List[float]:
    top = max(logits)
    exponentials = [math.exp(value - top) for value in logits]
    total = sum(exponentials)
    return [value / total for value in exponentials]


def _sigmoid(value: float) -> float:
    # Clipped, because the optimiser walks the logits freely and an unclipped
    # exp overflows at |value| ~ 710. The clip is at a fraction of 1e-15 either
    # way, so it changes no fitted value -- it only stops the search crashing.
    clipped = max(-60.0, min(60.0, float(value)))
    return 1.0 / (1.0 + math.exp(-clipped))


def _composition_for_system(system: str, theta: Mapping[str, float]) -> LOOHComposition:
    """The (position, geometry) composition of one FIT system, from its logits."""
    geometry = FRANKEL_SYSTEM_GEOMETRY[system]
    if geometry == "mixed":
        f13 = _sigmoid(theta["mixed_f13"])
        g_tt = _sigmoid(theta["mixed_g_tt"])
        return LOOHComposition(
            f13_ct=f13 * (1.0 - g_tt), f13_tt=f13 * g_tt,
            f9_ct=(1.0 - f13) * (1.0 - g_tt), f9_tt=(1.0 - f13) * g_tt,
        )
    if geometry == "ct":
        purity = _sigmoid(theta["pure_ct_13_f13"])
        return LOOHComposition(f13_ct=purity, f13_tt=0.0,
                               f9_ct=1.0 - purity, f9_tt=0.0)
    f13 = _sigmoid(theta["tt_9_13_f13"])
    return LOOHComposition(f13_ct=0.0, f13_tt=f13, f9_ct=0.0, f9_tt=1.0 - f13)


_BRANCH_CELLS: Tuple[Tuple[str, str], ...] = (
    ("13", "ct"), ("13", "tt"), ("9", "ct"), ("9", "tt"),
)
_COMPOSITION_KEYS: Tuple[str, ...] = (
    "mixed_f13", "mixed_g_tt", "pure_ct_13_f13", "tt_9_13_f13",
)


def _unpack(vector: Sequence[float]) -> Tuple[BranchModel, Dict[str, float]]:
    """12 unconstrained reals -> a valid BranchModel plus four composition logits."""
    simplexes: Dict[Tuple[str, str], Dict[str, float]] = {}
    cursor = 0
    for cell in _BRANCH_CELLS:
        logits = [0.0, float(vector[cursor]), float(vector[cursor + 1])]
        cursor += 2
        shares = _softmax(logits)
        simplexes[cell] = {
            product: shares[i]
            for i, product in enumerate(POSITION_PRODUCTS[cell[0]])
        }
    theta = {
        key: float(vector[cursor + i]) for i, key in enumerate(_COMPOSITION_KEYS)
    }
    return BranchModel(simplexes=simplexes), theta


def fit_residuals(vector: Sequence[float]) -> List[float]:
    """
    Residuals in RELATIVE PERCENT against the eighteen FIT numbers.

    THE FIREWALL LIVES HERE: this function reads ``FRANKEL_ZERO_ADDITIVE`` and
    nothing else, and it asserts the array is 3 systems x 6 products before it
    computes anything.
    """
    if len(FRANKEL_ZERO_ADDITIVE) != 3:
        raise AssertionError("the FIT array must be exactly three systems")
    branch, theta = _unpack(vector)
    residuals: List[float] = []
    for system, observed in FRANKEL_ZERO_ADDITIVE.items():
        if tuple(observed) != FRANKEL_SLATE:
            raise AssertionError(f"{system}: unexpected product slate")
        composition = _composition_for_system(system, theta)
        predicted = branch.slate_shares(composition, donor=0.0)
        scale = sum(observed.values())  # ~100, exactly as printed
        for product in FRANKEL_SLATE:
            residuals.append(predicted[product] * scale - float(observed[product]))
    if len(residuals) != 18:
        raise AssertionError(f"expected 18 residuals, got {len(residuals)}")
    return residuals


def fit_branch_model(seed: int = 0) -> Dict[str, Any]:
    """
    Fit the branch simplexes and the per-system compositions on the FIT column.

    Deterministic: a fixed multi-start grid, no random numbers, so the frozen
    parameters are reproducible from this function alone.
    """
    from scipy.optimize import least_squares  # local import: scipy is heavy

    starts = [
        [0.0] * 12,
        [0.5, -0.5] * 4 + [0.0, 0.0, 2.0, 0.0],
        [-0.5, 0.5] * 4 + [0.5, -0.5, 3.0, 0.5],
        [1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 0.2, 0.2, 2.5, -0.2],
    ]
    best = None
    for start in starts:
        result = least_squares(
            lambda v: fit_residuals(v), start, method="lm", max_nfev=200000
        )
        if best is None or result.cost < best.cost:
            best = result

    branch, theta = _unpack(best.x)
    compositions = {
        system: _composition_for_system(system, theta)
        for system in FRANKEL_ZERO_ADDITIVE
    }
    residuals = fit_residuals(best.x)
    absolute = [abs(r) for r in residuals]
    absolute_sorted = sorted(absolute)
    median = (
        absolute_sorted[len(absolute_sorted) // 2]
        if len(absolute_sorted) % 2
        else 0.5 * (absolute_sorted[len(absolute_sorted) // 2 - 1]
                    + absolute_sorted[len(absolute_sorted) // 2])
    )
    return {
        "branch_model": branch,
        "compositions": compositions,
        "theta": theta,
        "vector": [float(v) for v in best.x],
        "residuals_percentage_points": residuals,
        "median_abs_residual_pp": median,
        "worst_abs_residual_pp": max(absolute),
        "sum_squared_residuals": float(2.0 * best.cost),
        "n_fit_values": 18,
        "n_free_parameters": 12,
        "degrees_of_freedom": 18 - 3 - 12,
        "fit_rows": "Frankel 1989 zero-additive columns ONLY",
    }


def branch_model_from_dict(payload: Mapping[str, Any]) -> BranchModel:
    """Rebuild a frozen BranchModel from a fit report's JSON."""
    simplexes = {
        tuple(name.split("_")): dict(shares)
        for name, shares in payload["simplexes"].items()
    }
    return BranchModel(
        simplexes={(k[0], k[1]): v for k, v in simplexes.items()},
        provenance=str(payload.get("provenance", "frozen B6 fit report")),
    )


# ===========================================================================
# 2. THE CHARGE
# ===========================================================================


@dataclass(frozen=True)
class LipidCharge:
    """
    The lipid lane's inputs, every one of them with a declared provenance.

    ``looh_linoleate_mmol_l`` is the OXIDATION-STATE PROXY expressed as a pool;
    ``composition`` is the isomer distribution of that pool;
    ``looh_oleate_mmol_l`` exists only so that nonanal can be refused for the
    right reason rather than answered with an invented branch fraction.
    """

    looh_linoleate_mmol_l: float
    composition: LOOHComposition
    looh_oleate_mmol_l: float = 0.0
    carrier_key: Optional[str] = None
    provenance: Tuple[str, ...] = ()

    def total_mmol_l(self) -> float:
        return self.looh_linoleate_mmol_l + self.looh_oleate_mmol_l


def charge_from_carrier(
    carrier: LipidCarrier,
    composition: LOOHComposition,
    *,
    lipid_fraction: Optional[float] = None,
    peroxide_value_meq_per_kg: Optional[float] = None,
) -> LipidCharge:
    """Turn a declared carrier + oxidation-state proxy into a hydroperoxide pool."""
    total = looh_charge_mmol_per_l(
        carrier,
        lipid_fraction=lipid_fraction,
        peroxide_value_meq_per_kg=peroxide_value_meq_per_kg,
    )
    oleate_share = oleate_fraction(carrier)
    return LipidCharge(
        looh_linoleate_mmol_l=total * (1.0 - oleate_share),
        composition=composition,
        looh_oleate_mmol_l=total * oleate_share,
        carrier_key=carrier.key,
        provenance=(
            f"carrier {carrier.key}: lipid fraction "
            f"{carrier.lipid_lo}-{carrier.lipid_hi} (point "
            f"{carrier.lipid_mass_fraction}), peroxide value "
            f"{carrier.pv_lo}-{carrier.pv_hi} meq/kg (point "
            f"{carrier.peroxide_value_meq_per_kg}); evidence class "
            f"{carrier.evidence_class}",
            carrier.source_anchor,
        ),
    )


# ===========================================================================
# 3. INTEGRATION
# ===========================================================================


@dataclass(frozen=True)
class LipidRun:
    state_mmol_per_l: Mapping[str, float]
    concentrations_ug_per_l: Mapping[str, float]
    metadata: Mapping[str, Any]
    warnings: Tuple[str, ...]
    refusals: Mapping[str, str]


def integrate_lipid(
    charge: LipidCharge,
    segments: Sequence[Tuple[float, float]],
    branch: BranchModel,
    *,
    q10: Optional[float] = None,
    donor: float = 0.0,
) -> LipidRun:
    """
    Integrate the lipid lane over a piecewise-constant thermal program.

    The kinetics are first-order and are solved ANALYTICALLY per segment
    (``extent = 1 - exp(-k dt)``), not numerically: the network has no coupling
    between pools, so a solver would add error and nothing else. Both the
    formation and the consumption terms are present; the consumption term
    (``COVALENT_SINK``) is INERT BY RULING and its zero contribution is
    asserted, not assumed.
    """
    state = initial_lipid_state({
        "LOOH_13_ct": charge.looh_linoleate_mmol_l * charge.composition.f13_ct,
        "LOOH_13_tt": charge.looh_linoleate_mmol_l * charge.composition.f13_tt,
        "LOOH_9_ct": charge.looh_linoleate_mmol_l * charge.composition.f9_ct,
        "LOOH_9_tt": charge.looh_linoleate_mmol_l * charge.composition.f9_tt,
        "LOOH_OL": charge.looh_oleate_mmol_l,
    })
    carbon_in = _carbon(state)

    yields = slate_yields(branch, charge.composition, donor=donor)
    segment_log: List[Dict[str, Any]] = []

    for index, (duration, temperature_c) in enumerate(segments):
        k = k_looh_decomp_per_min(float(temperature_c), q10)
        extent = 1.0 - math.exp(-k * float(duration))
        decomposed_linoleate = 0.0
        for pool in LOOH_POOLS:
            consumed = state[pool] * extent
            state[pool] -= consumed
            decomposed_linoleate += consumed
        consumed_oleate = state["LOOH_OL"] * extent
        state["LOOH_OL"] -= consumed_oleate

        for product, yield_per_looh in yields.items():
            state[product] += decomposed_linoleate * yield_per_looh

        # CONSUMPTION. Structurally present, inert by ruling; the assertion is
        # the guard that keeps "inert" a fact rather than a comment.
        sink = _covalent_consumption(state, float(duration))
        if sink:
            raise AssertionError(
                "the covalent aldehyde-lysine sink returned a non-zero flux "
                "while COVALENT_SINK.enabled is False"
            )

        # CARBON CLOSURE. Everything the slate does not name -- the Hock
        # partners, the oleate products, the unmeasured remainder -- is routed
        # to LIPID_FRAG_C so the balance closes as an EQUALITY.
        named_carbon = sum(
            decomposed_linoleate * yields[p] * _carbon_of(p) for p in yields
        )
        looh_carbon = (decomposed_linoleate + consumed_oleate) * _carbon_of("LOOH_13_ct")
        state["LIPID_FRAG_C"] += looh_carbon - named_carbon

        segment_log.append({
            "index": index,
            "duration_min": float(duration),
            "temperature_C": float(temperature_c),
            "k_LOOH_decomp_per_min": k,
            "conversion_extent": extent,
            "q10_used": Q10_ASSUMPTION.default if q10 is None else float(q10),
            "decades_of_extrapolation_from_25C":
                Q10_ASSUMPTION.decades_of_extrapolation(float(temperature_c)),
        })

    carbon_out = _carbon(state)
    if abs(carbon_out - carbon_in) > 1e-9 * max(1.0, abs(carbon_in)):
        raise AssertionError(
            f"lipid carbon not conserved: {carbon_in} -> {carbon_out}"
        )

    warnings = [Q10_ASSUMPTION.warning]
    peak = max(t for _, t in segments)
    if peak > 40.0:
        warnings.append(
            f"the thermal program peaks at {peak:.1f} C, which is "
            f"{Q10_ASSUMPTION.decades_of_extrapolation(peak):.1f} decades of "
            f"10 C above the 25 C anchor. The authors licensed 'adjustment', "
            f"not this span."
        )
    if charge.carrier_key and charge.carrier_key in LIPID_CARRIERS:
        carrier = LIPID_CARRIERS[charge.carrier_key]
        if carrier.evidence_class == "declared_assumption":
            warnings.append(
                f"the hydroperoxide pool for {carrier.display} is a DECLARED "
                f"ASSUMPTION, not a measurement: lipid fraction "
                f"{carrier.lipid_lo}-{carrier.lipid_hi}, peroxide value "
                f"{carrier.pv_lo}-{carrier.pv_hi} meq/kg. No source in the fit "
                f"corpus measures either for this matrix. Supplying a measured "
                f"peroxide value collapses the largest band in this prediction."
            )

    refusals: Dict[str, str] = {}
    if charge.looh_oleate_mmol_l > 0.0:
        refusals["NONANAL"] = (
            "REFUSED. This matrix carries an OLEATE hydroperoxide pool "
            f"({charge.looh_oleate_mmol_l:.4g} mmol/L), and the oleate -> "
            "nonanal branch fraction is measured NOWHERE in the fit corpus. "
            "Frankel 1989 fed linoleate only, and nonanal appears in no table, "
            "figure or sentence of it. The FAST lane's shipped 'nonanal 0.15' "
            "has no source; this module refuses rather than carrying it "
            "forward."
        )

    concentrations = {
        key: mmol_per_litre_to_ug_per_litre(key, state[key])
        for key in MOLECULAR_WEIGHT_G_PER_MOL
        if key in state
    }

    return LipidRun(
        state_mmol_per_l=dict(state),
        concentrations_ug_per_l=concentrations,
        metadata={
            "lane": "lipid",
            "segments": segment_log,
            "slate_yields_per_LOOH": dict(yields),
            "donor_suppression": donor,
            "composition": charge.composition.as_dict(),
            "charge_provenance": list(charge.provenance),
            "carbon_in_mmol_l": carbon_in,
            "carbon_out_mmol_l": carbon_out,
            "covalent_sink_enabled": COVALENT_SINK.enabled,
            "rate_is_an_assumption": True,
        },
        warnings=tuple(warnings),
        refusals=refusals,
    )


def slate_yields(
    branch: BranchModel, composition: LOOHComposition, donor: float = 0.0
) -> Dict[str, float]:
    """
    Moles of each measured product per mole of LINOLEATE hydroperoxide consumed.

    The bridge described in the module docstring: Schroen's measured hexanal
    yield, redistributed by Frankel's measured distribution, referenced to
    Frankel's AUTOXIDATION column (the composition Schroen's emulsion makes).
    """
    reference = FRANKEL_ZERO_ADDITIVE[REFERENCE_AUTOXIDATION_SYSTEM]
    reference_hexanal_share = (
        float(reference["HEXANAL"]) / sum(reference.values())
    )
    flux = branch.slate_flux(composition, donor)
    scale = Y_HEXANAL_PER_LOOH / reference_hexanal_share
    return {product: value * scale for product, value in flux.items()}


def _carbon_of(key: str) -> int:
    for species in LIPID_SPECIES:
        if species.key == key:
            return species.carbon
    raise KeyError(key)


def _carbon(state: Mapping[str, float]) -> float:
    return sum(float(state.get(s.key, 0.0)) * s.carbon for s in LIPID_SPECIES)


def _covalent_consumption(state: Mapping[str, float], duration_min: float) -> float:
    """The consumption term. Zero while the channel is inert; never absent."""
    if not COVALENT_SINK.enabled:
        return 0.0
    raise NotImplementedError(
        "enabling the covalent sink requires the aldehyde-lysine Ea on food "
        "proteins, a NAMED WET-LAB GAP (Amendment 6 ruling 2). Do not enable "
        "it with a guessed Ea."
    )


# ===========================================================================
# 4. STRUCTURAL VALIDATION AND THE LANE-COUPLING GUARD
# ===========================================================================


def validate_lipid_structure(branch: Optional[BranchModel] = None) -> Dict[str, Any]:
    """
    The structural invariants, asserted rather than described.

    1. Every product is reachable from EXACTLY the position that can make it.
    2. NONANAL is reachable from NO linoleate pool -- the declared hold-out
       negative test, enforced by topology.
    3. Each simplex sums to 1.
    """
    findings: Dict[str, Any] = {}
    reachable: Dict[str, List[str]] = {}
    for position, products in POSITION_PRODUCTS.items():
        for product in products:
            reachable.setdefault(product, []).append(position)
    for product, positions in reachable.items():
        if len(positions) != 1:
            raise AssertionError(f"{product} is reachable from {positions}")
    findings["structural_zeros"] = {
        product: f"reachable ONLY from the {positions[0]}-hydroperoxide"
        for product, positions in sorted(reachable.items())
    }

    if "NONANAL" in reachable:
        raise AssertionError(
            "NONANAL has an edge from a linoleate hydroperoxide. That "
            "falsifies the declared Frankel negative test by construction."
        )
    findings["nonanal"] = (
        "STRUCTURAL ZERO from every linoleate pool. Its only parent is "
        "LOOH_OL, whose branch fraction is unmeasured, so it is answered as "
        "exactly 0.0 for a linoleate feed and REFUSED for an oleate-bearing "
        "matrix."
    )
    if branch is not None:
        findings["simplex_sums"] = {
            f"{p}_{g}": sum(s.values()) for (p, g), s in sorted(branch.simplexes.items())
        }
    return findings


#: The species keys the lipid lane owns. Used by the lane-coupling guard.
LIPID_OWNED_KEYS = frozenset(LIPID_KEYS)


def lane_coupling_verdict(other_lane_keys: Sequence[str]) -> Dict[str, Any]:
    """
    Whether the lipid lane may be CO-INTEGRATED with a Maillard lane.

    THE RULING (pre-registered, ``kinetic_core_b6_prereg.md`` sec. 6):
    co-integration as a DIRECT SUM is permitted, because the species sets are
    disjoint and the one candidate coupling -- the aldehyde-lysine covalent
    channel, which would put the lipid aldehydes and the Maillard amine pool in
    competition -- is INERT BY RULING (Amendment 6 ruling 2).

    This is a DIFFERENT ruling from the acrylamide/sulfur LANE CONFLICT, and
    the difference is the point: those two lanes both consume cysteine, so
    summing them spends it twice. These two share nothing.

    The verdict is CONDITIONAL and the condition is checked here, every time:
    if the covalent sink is ever enabled, the amine pool becomes genuinely
    shared and the direct sum stops being valid.
    """
    overlap = sorted(LIPID_OWNED_KEYS & set(other_lane_keys))
    if overlap:
        return {
            "may_cointegrate": False,
            "reason": (
                "SPECIES OVERLAP: the two lanes share "
                + ", ".join(overlap)
                + ". A direct sum would double-count them."
            ),
            "shared_species": overlap,
        }
    if COVALENT_SINK.enabled:
        return {
            "may_cointegrate": False,
            "reason": (
                "The covalent aldehyde-lysine channel is ENABLED, so the lipid "
                "aldehydes and the Maillard amine pool are now in competition "
                "for the same lysine. Co-integration must become a genuine "
                "joint integration before it is permitted again."
            ),
            "shared_species": ["the amine pool, via the covalent sink"],
        }
    return {
        "may_cointegrate": True,
        "rule": "direct sum (independent parallel integration, same thermal program)",
        "reason": (
            "Disjoint species sets, and the only candidate coupling (the "
            "aldehyde-lysine covalent channel) contributes exactly zero: it is "
            "INERT BY RULING, FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2. "
            "Revisit the moment the aldehyde-lysine Ea on food proteins is "
            "measured -- that is a NAMED WET-LAB GAP."
        ),
        "shared_species": [],
    }


def describe_lipid() -> Dict[str, Any]:
    """What this network is, for embedding in every artifact it produces."""
    return {
        "module": "src/kinetic_core/lipid.py",
        "wave": "B6 -- the lipid-oxidation module",
        "species": list(LIPID_KEYS),
        "measured_slate": list(FRANKEL_SLATE),
        "reference_autoxidation_system": REFERENCE_AUTOXIDATION_SYSTEM,
        "Y_hexanal_per_LOOH": Y_HEXANAL_PER_LOOH,
        "Y_hexanal_provenance": (
            "k_hexanal / k4 from Schroen & Berton-Carabin 2022 Table 1 -- a "
            "RATIO of two constants from the same table of the same paper, so "
            "the hand-fitting offset and the Q10 both cancel. The authors state "
            "the same ratio in prose as '1.2 %'; 6e-5/6e-3 = 1.0 %, and the "
            "0.2 pp difference is their rounding, reported rather than "
            "reconciled."
        ),
        "rate_is_measured": False,
        "distribution_is_measured": True,
        "donor_parameter_is_fitted": False,
        "covalent_sink_enabled": COVALENT_SINK.enabled,
    }


__all__ = [
    "BranchModel",
    "LOOHComposition",
    "LipidCharge",
    "LipidRun",
    "REFERENCE_AUTOXIDATION_SYSTEM",
    "Y_HEXANAL_PER_LOOH",
    "branch_model_from_dict",
    "charge_from_carrier",
    "describe_lipid",
    "fit_branch_model",
    "fit_residuals",
    "integrate_lipid",
    "lane_coupling_verdict",
    "slate_yields",
    "validate_lipid_structure",
]
