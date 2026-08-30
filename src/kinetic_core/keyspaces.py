"""
src/kinetic_core/keyspaces.py -- THE THREE KEY-SPACES, IN ONE PLACE (Wave Q1).

WAVE: Q1 (quality pass). No exam of its own: this module contains no science
and no arithmetic. It is a behaviour-preserving consolidation of a mapping that
had been written out by hand in five places.

WHY THIS EXISTS
---------------
One compound wears three different names inside this codebase, and the layer
boundaries are exactly where they change:

  1. DISPLAY NAME  -- what the user typed and what a report prints:
                      ``"2-methyl-3-furanthiol (MFT)"``, ``"hexanal"``,
                      ``"2,4-decadienal"``.
  2. SPECIES KEY   -- what the integrator's state vector is keyed by:
                      ``"MFT"``, ``"HEXANAL"``, ``"DECADIENAL"``.
                      The display -> species map is per-request and lives on
                      ``EnvelopeDeclaration.mapped_targets``, because which
                      display names resolve at all depends on the lane.
  3. B4 KEY        -- what the B4 threshold/OAV tables are keyed by:
                      ``"MFT"`` (unchanged), but ``"hexanal"``,
                      ``"tt_2_4_decadienal"`` for the lipid lane, whose B4
                      records were transcribed under the source papers' own
                      names. The species -> B4 map is ``B4_COMPOUND_KEY``, and
                      a species in ``NO_B4_RECORD`` has NO B4 key at all.

THE BUG CLASS THIS CLOSES
-------------------------
Handing a key from one space to a table keyed by another does not raise -- it
returns ``None`` or ``NoMeasuredThreshold``, which every caller reads as "this
compound has no measured odour threshold". That is a FALSE REFUSAL: it prints
the absence of evidence for a compound the corpus does in fact have a measured
threshold for. It has been found twice by hand already (the B5 cutover found
``.oav()`` being handed display names; Wave Q1 found ``predict_core`` looking
up ``per_species["HEXANAL"]`` in a table keyed ``"hexanal"``, which silently
dropped the OAV of every lipid compound from every report). Both were silent.

So the rule is: NEVER index a B4 table with a key you did not get from
``b4_key()``, and never treat a ``None`` from it as "no threshold" -- ask
``no_b4_reason()``, which returns the RECORDED reason why the compound has no
B4 record, so a report can print why instead of printing nothing.

DECLARED GAPS
-------------
* The display -> species direction is not invertible here: a species key can
  be reached from several display names (aliases), and this module deliberately
  does not choose one. A caller that has only a species key and wants a label
  should print the species key.
* ``NO_B4_RECORD`` is a lipid-lane table. If another lane ever gains a product
  with no B4 structural record, it must be added there, not special-cased in a
  caller.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional

from .species_lipid import B4_COMPOUND_KEY, NO_B4_RECORD

__all__ = [
    "CompoundKeys",
    "b4_key",
    "keys_for",
    "no_b4_reason",
    "species_key_for",
]


@dataclass(frozen=True)
class CompoundKeys:
    """One compound's name in all three key-spaces at once."""

    #: What the user typed / what a report prints.
    display: str
    #: What the integrator's state vector is keyed by.
    species: str
    #: What the B4 threshold/OAV tables are keyed by. ``None`` -- and ONLY
    #: ``None`` -- means the compound has no B4 structural record; the reason
    #: is in :attr:`no_b4_reason`.
    b4: Optional[str]
    #: The RECORDED reason there is no B4 record, or ``None`` when there is one.
    no_b4_reason: Optional[str] = None

    @property
    def has_b4_record(self) -> bool:
        return self.b4 is not None


def species_key_for(display: str, mapped_targets: Mapping[str, str]) -> str:
    """
    DISPLAY NAME -> SPECIES KEY, via one request's declaration.

    ``mapped_targets`` is ``EnvelopeDeclaration.mapped_targets`` (or the same
    dict read back out of ``as_dict()``). A display name the declaration does
    not carry is returned unchanged, which is what every existing caller did.
    """
    return dict(mapped_targets or {}).get(display, display)


def b4_key(species: str) -> Optional[str]:
    """
    SPECIES KEY -> B4 KEY, or ``None`` when the species has no B4 record.

    Never returns the species key for a ``NO_B4_RECORD`` species: that is the
    exact confusion that makes a table lookup answer "no measured threshold"
    for a reason it cannot state.
    """
    if species in NO_B4_RECORD:
        return None
    return B4_COMPOUND_KEY.get(species, species)


def no_b4_reason(species: str) -> Optional[str]:
    """The recorded reason ``species`` has no B4 record, or ``None``."""
    return NO_B4_RECORD.get(species)


def keys_for(display: str, mapped_targets: Mapping[str, str]) -> CompoundKeys:
    """All three keys for one display name, under one request's declaration."""
    species = species_key_for(display, mapped_targets)
    return CompoundKeys(
        display=display,
        species=species,
        b4=b4_key(species),
        no_b4_reason=no_b4_reason(species),
    )
