from src.chem_utils import Species  # noqa: E402
from src.reaction_templates import _thiol_addition  # noqa: E402


def test_thiohemiacetal_pathway():
    """Furfural -> thiohemiacetal -> 2-furfurylthiol.

    UPDATED 2026-08-27 (audit remediation, Wave G1 fix 7). Step 2 used to read

        thiohemiacetal + H2S -> FFT + [S] + H2O

    which balances only by ejecting an atom of ELEMENTAL SULFUR with no
    mechanistic counterpart — and that fictitious `[S]` was the sole reason the
    `Radical_Crosstalk` family existed (it mopped the [S] up by consuming MFT).
    Converting an aryl thiohemiacetal to the thiol is a net REDUCTION, so the
    step now consumes reducing equivalents (H2 as the documented reductone
    token) and the pathway is pool-gated on them.
    """
    furfural = Species(label="furfural", smiles="O=Cc1ccco1")
    h2s = Species(label="H2S", smiles="S")
    h2 = Species(label="H2", smiles="[HH]")
    water = Species(label="water", smiles="O")

    steps = _thiol_addition([furfural, h2s, h2, water])

    step1 = next(s for s in steps if s.reaction_family == "Thiohemiacetal_Formation")
    assert "OC(S)c1ccco1" in [p.smiles for p in step1.products]

    step2 = next(s for s in steps if s.reaction_family == "Thiol_Dehydration")
    assert "SCc1ccco1" in [p.smiles for p in step2.products]
    assert "[S]" not in [p.smiles for p in step2.products]
    assert "[HH]" in [r.smiles for r in step2.reactants]


def test_no_reductant_means_no_thiohemiacetal_lane():
    """Without reducing equivalents in the pool the reduction cannot run."""
    furfural = Species(label="furfural", smiles="O=Cc1ccco1")
    h2s = Species(label="H2S", smiles="S")

    steps = _thiol_addition([furfural, h2s])
    assert [s for s in steps if s.reaction_family == "Thiol_Dehydration"] == []
