"""
src/kinetic_core/parameters_lipid_b28.py

WAVE B28's LIPID DATA, DELIBERATELY IN ITS OWN MODULE.
======================================================

Frankel, Neff & Selke 1981 measures the two products the lipid lane used to refuse. It is a
DIFFERENT paper from the one wave B6 fitted the lane on, on a different substrate, at a different
injector temperature, against a different denominator -- and this file exists so that separation is
physical rather than a comment.

THERE IS ALSO A HARD REASON. `tests/unit/test_kinetic_core_b6.py` runs a literal firewall over the
B6 lipid package: twelve values that appear ONLY in Frankel 1989's alpha-tocopherol and
1,4-cyclohexadiene columns, which are that wave's HOLD-OUT, and whose appearance anywhere in the
package would mean a hold-out number had leaked into a parameter. One of them is "9.7", and Frankel
1981's photosensitized methyl octanoate share is also 9.7. That is a coincidence between two
papers, not a leak -- and the wrong ways to resolve it are to write the number as 9.70 so the regex
misses it, or to strike it from the firewall's list. Both would trade a real guard for a cosmetic
pass. Putting this paper's numbers in a file the B6 firewall does not cover keeps that guard exactly
as strict as it was, and states in the open why the collision is not a leak.

Nothing here enters the B6 objective.
"""
from __future__ import annotations

from typing import Mapping, Tuple

# ===========================================================================
# WAVE B28 (2026-09-09) -- THE TWO PRODUCTS THIS LANE USED TO REFUSE
# ===========================================================================
# Frankel, Neff & Selke 1981 (Lipids 16:279-285): pure hydroperoxides of methyl
# oleate, linoleate and linolenate, autoxidised and photosensitized, thermolysed
# neat in a GC injector port at 210 C. Same laboratory and same first author as
# the 1989 slate this lane is FITTED on; a different substrate, a different
# temperature and a different denominator.
#
# THE DENOMINATOR IS THE THING TO WATCH. 1989's shares are fractions of SIX
# measured peaks. 1981's are fractions of a whole ~20-peak chromatogram that
# includes 6-12 % explicitly unidentified. A 1981 share is therefore
# systematically smaller than a 1989 share for the same product, and the two
# must never be pooled as printed. Nothing below enters the B6 objective.

FRANKEL1981_ANCHOR = (
    "Frankel, Neff & Selke (1981) Lipids 16:279-285, Tables II (oleate, p. 281) "
    "and III (linoleate, p. 282); data/articles/frankel1981.pdf, read by Wave B28"
)
FRANKEL1981_DOSSIER = "data/lit/extraction_dossiers/frankel1981_extraction.md"

#: Table II, p. 281. Relative % of the whole volatile chromatogram.
#: THE AUTOXIDATION COLUMN IS NOT AN INDEPENDENT MEASUREMENT. It is footnoted
#: "Data from ref. 21" = Selke, Frankel & Neff 1978, i.e. that single
#: determination republished. Anyone who reads 15 % in both papers and calls it
#: a replicate pair is reading one number twice. The photosensitized column is
#: the only independent oleate determination in the corpus.
FRANKEL1981_OLEATE_SLATE: Mapping[str, Mapping[str, float]] = {
    "oleate_autoxidised": {
        "HEPTANE": 4.4, "OCTANE": 2.7, "HEPTANAL": 0.5, "OCTANOL_1": 0.4,
        "OCTANAL": 11.0, "ME_HEPTANOATE": 1.5, "HEPTANOL_1": 0.4,
        "NONANAL": 15.0, "ME_OCTANOATE": 5.0, "NONENAL_2": 0.5, "DECANAL": 3.9,
        "ME_NONANOATE": 1.5, "DECENAL_2": 5.4, "UNDECENAL_2": 1.7,
        "ME_8_OXOOCTANOATE": 3.5, "ME_9_OXONONANOATE": 15.0,
        "ME_10_OXODECANOATE": 12.0, "ME_10_OXO_8_DECENOATE": 3.4,
        "ME_11_OXO_9_UNDECENOATE": 5.8, "UNIDENTIFIED": 6.4,
    },
    "oleate_photosensitized": {
        "HEPTANE": 4.6, "OCTANE": 10.0, "HEPTANAL": 0.5, "OCTANOL_1": 1.0,
        "OCTANAL": 3.8, "ME_HEPTANOATE": 4.9, "HEPTANOL_1": 0.4,
        "NONANAL": 10.0, "ME_OCTANOATE": 9.7, "NONENAL_2": 0.7, "DECANAL": 2.0,
        "ME_NONANOATE": 0.8, "DECENAL_2": 12.0, "UNDECENAL_2": 7.1,
        "ME_8_OXOOCTANOATE": 3.0, "ME_9_OXONONANOATE": 11.0,
        "ME_10_OXODECANOATE": 1.7, "ME_10_OXO_8_DECENOATE": 5.0,
        "ME_11_OXO_9_UNDECENOATE": 4.6, "UNIDENTIFIED": 6.7,
    },
}
#: Which of the two columns is a fresh measurement, keyed so no reader has to
#: remember the footnote.
FRANKEL1981_OLEATE_PROVENANCE: Mapping[str, str] = {
    "oleate_autoxidised": "[C] Selke 1978 republished -- NOT an independent replicate",
    "oleate_photosensitized": "[M] measured here; the corpus's only independent oleate column",
}

#: Table III, p. 282. The linoleate column, carried for ONE purpose: the
#: cross-laboratory check against the 1989 slate this lane is fitted on, and the
#: alkylfuran ratio below. It is NOT a second fit source.
FRANKEL1981_LINOLEATE_SLATE: Mapping[str, Mapping[str, float]] = {
    "linoleate_autoxidised": {
        "PENTANE": 9.9, "HEXANAL": 15.0, "ME_OCTANOATE": 15.0,
        "DECADIENAL": 14.0, "ME_9_OXONONANOATE": 19.0, "PENTYLFURAN": 2.4,
    },
    "linoleate_photosensitized": {
        "PENTANE": 4.3, "HEXANAL": 17.0, "ME_OCTANOATE": 7.6,
        "DECADIENAL": 4.3, "ME_9_OXONONANOATE": 22.0, "PENTYLFURAN": 0.6,
    },
}
#: The C14 oxo-ester that is 20 % of the 1989 slate is ABSENT from 1981, and the
#: reason is analytical, not chemical: "no authentic references were available".
#: Any comparison of the two slates must drop it from BOTH or it is meaningless.
FRANKEL1981_MISSING_FROM_1981 = (
    "ME_13_OXO_TRIDECADIENOATE: not identified in 1981 for want of an authentic "
    "reference compound, not absent from the chemistry. It is 20 % of the 1989 "
    "mixed-geometry slate, so it must be dropped from BOTH slates before they "
    "are compared, and the comparison renormalised."
)

#: THE FORM 2-PENTYLFURAN SHIPS IN, and the reason it is a ratio.
#: 1981's shares carry 1981's ~20-peak denominator and the lane's carry 1989's
#: six-peak one. Dividing by hexanal, which is in the same column and is already
#: a modelled species, cancels the denominator exactly.
PENTYLFURAN_PER_HEXANAL: Mapping[str, float] = {
    "linoleate_autoxidised": 2.4 / 15.0,       # 0.160
    "linoleate_photosensitized": 0.6 / 17.0,   # 0.0353
}
PENTYLFURAN_ORIGIN_UNASSIGNED = (
    "The paper does NOT assign 2-pentylfuran to a hydroperoxide isomer. Its "
    "Table III Origin column reads '?', its Results say the origin 'is not well "
    "established', and its Discussion only speculates a 10-hydroperoxide "
    "intermediate on a 1966 citation. So this is a share of an autoxidised "
    "linoleate pool's volatile slate and NOT a mechanistic branch from a named "
    "isomer. No position gets a structural zero for it, and none is invented."
)

#: THE DECLARED ASSUMPTION AN ABSOLUTE NONANAL ANSWER RESTS ON, and it is the
#: reason the refusal changes state rather than disappearing.
#: Frankel 1981 prints PEAK AREAS: no internal standard, no response factors, no
#: replicates, no error. The lane turns its LINOLEATE distribution into moles per
#: hydroperoxide by anchoring on Schroen's separately measured hexanal yield.
#: There is no such anchor for oleate anywhere in the corpus.
OLEATE_MOLAR_ANCHOR_ASSUMPTION = (
    "DECLARED ASSUMPTION, not a measurement: that the named-product molar yield "
    "per mole of OLEATE hydroperoxide equals the measured one per mole of "
    "LINOLEATE hydroperoxide. Nothing in the corpus measures an absolute yield "
    "from an oleate hydroperoxide. The band below is propagated into every "
    "absolute nonanal answer and a warning is mandatory on each one. The "
    "precedent is the soy-paste protein loading in the matrix layer, which is "
    "carried the same way for the same reason."
)
#: The band, as a multiplicative factor on the linoleate anchor. Wide on purpose:
#: an oleate hydroperoxide is a MONOENE's, it has no bis-allylic hydrogen, and
#: its decomposition is slower and less productive of volatiles than a diene's.
OLEATE_MOLAR_ANCHOR_BAND: Tuple[float, float] = (0.2, 1.0)
OLEATE_MOLAR_ANCHOR_CENTRE: float = 0.45   # geometric mean of the band

