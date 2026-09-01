"""Wave I (2026-08-27) chemistry / citation remediation regressions — unit tier.

Everything here is fast: SMARTS matching, single-template calls, file contents.
The network-enumeration and benchmark-prediction regressions for the same wave
live in tests/scientific/test_wave_i_network_chemistry.py.

The findings pinned here, in the numbering of the cold-start red-team audit:

  H4  (fix 8)  the 2[H] reducing-equivalent token that gates the MFT step had
               exactly one producer reachable from a cysteine/sugar system, the
               pyrazine aromatisation.
  H3  (fix 9)  the `Lipid_Schiff_Base` SMARTS was much broader than its own
               comment claimed.
  H5  (fix 10) physically broken Arrhenius pairs, an ordering inversion between
               the two barrier tables, and two docstrings claiming a disabled
               feature was enabled. ANNOTATION ONLY — no runtime value moved,
               and this file pins that nothing moved.
  H1  (fix 4)  the flagship mechanism DOI resolved to the wrong paper.
      (fix 12) the demoted MFT shortcut fired for pentoses despite a
               hexose-only docstring.
      (fix 18) the acrylamide template did not consume the sugar.
"""

import importlib
import math
import re
import sys
from pathlib import Path

import pytest
import yaml
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chem_utils import Species  # noqa: E402
from src.reaction_templates import (  # noqa: E402
    _acrylamide_formation,
    _mft_pathway,
    _thiol_reductant_pool,
)
from src.smirks_engine import _SMIRKS_RULES  # noqa: E402


# ── helpers ────────────────────────────────────────────────────────────────

def _formula(species_list):
    counts: dict[str, int] = {}
    for sp in species_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(sp.smiles))
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return counts


def _lipid_schiff_smarts() -> str:
    for name, family, smirks, _gate in _SMIRKS_RULES:
        if family == "Lipid_Schiff_Base":
            return smirks
    raise AssertionError("the Lipid_Schiff_Base rule disappeared from _SMIRKS_RULES")


# ══ FIX 4 — flagship DOI repair ════════════════════════════════════════════

_BAD_DOI = "10.1021/jf60200a038"
_GOOD_DOI = "10.1021/jf60199a045"

_DOI_REPAIR_SITES = (
    "src/reaction_templates.py",
    "src/smirks_engine.py",
    "src/curated_pathways.py",
    "src/barrier_constants.py",
    "scripts/generators/refit_sulfur_barriers_hofmann.py",
)


@pytest.mark.parametrize("rel", _DOI_REPAIR_SITES)
def test_flagship_doi_is_repaired_everywhere(rel):
    """10.1021/jf60200a038 is live but resolves to a gossypol/rat-feeding paper.

    The real van den Ouweland & Peer 1975 JAFC paper -- norfuraneol + H2S ->
    2-methyl-3-furanthiol, the study the whole sulfur branch rests on -- is
    10.1021/jf60199a045. A live-but-wrong anchor is the hardest citation defect
    to notice, which is exactly why it is pinned rather than left to review.
    """
    text = (ROOT / rel).read_text(encoding="utf-8")
    if rel == "src/reaction_templates.py":
        # The canonical repair record legitimately quotes the old DOI once, on
        # its `old:` line. Nowhere else may.
        occurrences = [
            line for line in text.splitlines()
            if _BAD_DOI in line and not line.strip().startswith("#   old:")
        ]
        assert not occurrences, occurrences
    else:
        assert _BAD_DOI not in text, f"{rel} still cites the gossypol paper"
    assert _GOOD_DOI in text, f"{rel} lost the corrected van den Ouweland anchor"


def test_doi_repair_record_is_complete_and_canonical():
    """One canonical record, carrying old/new/date/basis, as the gate demands."""
    text = (ROOT / "src/reaction_templates.py").read_text(encoding="utf-8")
    block = text.split("doi_repair 2026-08-27", 1)
    assert len(block) == 2, "the canonical doi_repair record is missing"
    record = block[1][:2200]
    for key in ("old:", "new:", "date:", "basis:"):
        assert key in record, f"repair record is missing {key!r}"
    assert _BAD_DOI in record and _GOOD_DOI in record
    assert "gossypol" in record, "the basis must say what the old DOI actually resolves to"


def test_citation_gate_scans_code_and_prose():
    sys.path.insert(0, str(ROOT / "scripts"))
    gate = importlib.import_module("ci.citation_gate")

    for glob in ("src/**/*.py", "docs/**/*.md", "README.md", "AUDIT.md"):
        assert glob in gate.TEXT_SCAN_GLOBS, f"{glob} is not scanned for DOIs"

    # The structured baseline still ratchets down and is still empty.
    assert gate.WAIVERS == ()

    (
        blocking,
        waived,
        stale,
        carried,
        stale_carried,
        dois,
        files_scanned,
    ) = gate.run_offline_detailed()
    assert files_scanned > 0
    assert not blocking, [str(v) for v in blocking]
    assert not waived
    assert not stale
    assert not stale_carried, stale_carried

    # 2026-08-27 (Wave I fix 4): turning the text surface on immediately found a
    # fabricated DOI in docs/slr_benchmark_evaluation.md §7.3 (Hao et al. 2025,
    # an impossible zero-padded Elsevier article number). It was carried in
    # TEXT_SURFACE_WAIVERS for a few hours only, because docs/** sat outside the
    # editing scope of the sub-task that widened the scan, and was then FIXED AT
    # SOURCE rather than waived: the reference is withdrawn in the document.
    # Both lists are therefore empty again, and this assertion is what stops a
    # future widening from quietly parking a defect in the carry-over instead of
    # fixing it.
    assert gate.TEXT_SURFACE_WAIVERS == (), (
        "the text-surface carry-over is a ratchet that only shrinks; a new entry "
        "is a new UNFIXED citation defect and needs a named owner, not a waiver"
    )
    assert carried == [], [str(v) for v in carried]


def test_text_doi_extraction_survives_prose_punctuation():
    sys.path.insert(0, str(ROOT / "scripts"))
    gate = importlib.import_module("ci.citation_gate")

    cases = {
        "see 10.1021/jf60199a045.": "10.1021/jf60199a045",
        "(DOI 10.1021/jf9705983)": "10.1021/jf9705983",
        "`10.1016/j.foodchem.2004.04.006`": "10.1016/j.foodchem.2004.04.006",
        # markdown-escaped, as in data/Gemini_Deep_Research/* (the docs/research/archives/
        # duplicates of those same files were deleted 2026-08-28, Wave S5)
        r"10.1016/0891-5849\\": "10.1016/0891-5849",
        # URL fragment is part of the link, not the DOI
        "10.1080/10942912.2019.1573830#:~:text=x": "10.1080/10942912.2019.1573830",
    }
    for line, expected in cases.items():
        found = [gate._trim_text_doi(m.group(0)) for m in gate.TEXT_DOI_RE.finditer(line)]
        assert expected in found, (line, found)
        for token in found:
            out: list = []
            gate._check_doi_token(token, "x", "#L1", out)
            assert not out, [str(v) for v in out]


def test_text_surface_still_catches_a_confabulated_doi():
    sys.path.insert(0, str(ROOT / "scripts"))
    gate = importlib.import_module("ci.citation_gate")

    for bad in (
        "10.1021/acs.jafc.de_leyn_2019",
        "10.1016/j.foodchem.2015.00000",
        "10.1016/j.foodres.2025.001279",
        "10.xxxx/placeholder",
    ):
        out: list = []
        gate._check_doi_token(bad, "src/fake.py", "#L1", out)
        assert out, f"{bad} passed the text-surface check"
        assert out[0].check == "confabulation-signature"


# ══ FIX 8 — MFT H2 gating ══════════════════════════════════════════════════

def test_thiol_reductant_pool_emits_a_balanced_cysteine_to_cystine_step():
    """`2 cysteine -> cystine + 2[H]`: the decoupled source of the 2[H] token.

    Before Wave I the only producer of `[HH]` reachable from a ribose/cysteine
    system was `Aminoketone_Condensation`, so MFT was a downstream dependent of
    pyrazine chemistry (red-team H4). Cysteine is a PRECURSOR, so this donor
    cannot be switched off by an unrelated lane.
    """
    cysteine = Species("L-Cysteine", "NC(CS)C(=O)O")
    steps = _thiol_reductant_pool([cysteine])
    assert len(steps) == 1
    step = steps[0]
    assert step.reaction_family == "Thiol_Oxidation"
    assert [r.smiles for r in step.reactants] == ["NC(CS)C(=O)O"] * 2
    assert any(p.smiles == "[HH]" for p in step.products)
    assert _formula(step.reactants) == _formula(step.products), (
        _formula(step.reactants), _formula(step.products)
    )


def test_thiol_reductant_pool_is_silent_without_cysteine():
    assert _thiol_reductant_pool([Species("glycine", "NCC(=O)O")]) == []


# ══ FIX 9 — Lipid_Schiff_Base SMARTS ═══════════════════════════════════════

_LIPID_SCHIFF_MATCHES = {
    # genuine lipid / Strecker C3+ aldehydes: must still match
    "hexanal": ("CCCCCC=O", True),
    "nonanal": ("CCCCCCCCC=O", True),
    "propanal": ("CCC=O", True),
    "methional": ("CSCCC=O", True),
    # C2: the comment always said C3+
    "acetaldehyde": ("CC=O", False),
    # alpha-hydroxy carbonyls: the comment always said "excludes sugars"
    "glycolaldehyde": ("OCC=O", False),
    "glyceraldehyde": ("O=CC(O)CO", False),
    "glucose_open_chain": ("O=CC(O)C(O)C(O)C(O)CO", False),
    # alpha-amino aldehyde: already excluded before Wave I, must stay excluded
    "2-aminoethanal": ("NCC=O", False),
}


@pytest.mark.parametrize("name,payload", sorted(_LIPID_SCHIFF_MATCHES.items()))
def test_lipid_schiff_base_smarts_matches_only_what_its_comment_claims(name, payload):
    smiles, should_match = payload
    aldehyde_pattern = _lipid_schiff_smarts().split(">>")[0].split(".")[0]
    patt = Chem.MolFromSmarts(aldehyde_pattern)
    assert patt is not None, "the Lipid_Schiff_Base aldehyde SMARTS no longer parses"
    got = Chem.MolFromSmiles(smiles).HasSubstructMatch(patt)
    assert got is should_match, (
        f"{name} ({smiles}): Lipid_Schiff_Base match={got}, expected {should_match}. "
        "The rule's comment claims C3+ aliphatic aldehydes with no alpha-hydroxyl; "
        "see the Wave I fix 9 note in src/smirks_engine.py."
    )


# ══ FIX 12 — legacy MFT shortcut is hexose-only ════════════════════════════

def _shortcut_steps(deoxyosone_label, deoxyosone_smiles):
    pool = [
        Species(deoxyosone_label, deoxyosone_smiles),
        Species("H2S", "S"),
        Species("H2", "[HH]"),
    ]
    return _mft_pathway(pool)


def test_legacy_mft_shortcut_no_longer_fires_for_pentoses():
    """Its own docstring said hexose-only; until Wave I it accepted pentoses too.

    A pentose reaches MFT through the accepted three-step norfuraneol route, so
    the demoted lump was double-producing pentose MFT and flattering the
    pentose >> hexose ordering with a contribution the Wave G1 review had
    already labelled "NOT a mechanism".
    """
    assert _shortcut_steps("ribose-deoxyosone-3", "O=CC(=O)CC(O)CO") == []


def test_legacy_mft_shortcut_still_gives_hexoses_a_path_to_mft():
    """RESTRICTED, not deleted: the hexose furanone branch needs an extra
    reduction, so deleting the lump would make MFT unreachable from glucose."""
    steps = _shortcut_steps("glucose-deoxyosone-3", "O=CC(=O)CC(O)C(O)CO")
    assert [s.reaction_family for s in steps] == ["Thiol_Addition_Hexose_Legacy_Shortcut"]
    assert _formula(steps[0].reactants) == _formula(steps[0].products)


# ══ FIX 18 — acrylamide sugar-spectator ════════════════════════════════════

_ASPARAGINE = Species("L-Asparagine", "NC(CC(N)=O)C(=O)O")


@pytest.mark.parametrize(
    "label,smiles",
    [
        ("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO"),
        ("D-Ribose", "O=CC(O)C(O)C(O)CO"),
        ("D-Fructose", "OCC(=O)C(O)C(O)C(O)CO"),
    ],
)
def test_acrylamide_template_consumes_the_sugar(label, smiles):
    """The sugar used to be regenerated ("conserved here for simplicity"), i.e.
    formally catalytic, which contradicted src/safety.py::predict_acrylamide --
    that model is explicitly second order, dA/dt = kf*[Asn]*[Sugar]."""
    steps = _acrylamide_formation(Species(label, smiles), _ASPARAGINE)
    assert len(steps) == 1
    step = steps[0]

    product_smiles = [p.smiles for p in step.products]
    assert smiles not in product_smiles, "the sugar is still a spectator"
    assert "C=CC(=O)N" in product_smiles
    assert any("deoxyosone" in p.label for p in step.products), (
        "the consumed sugar must leave as a named species, not vanish"
    )
    assert _formula(step.reactants) == _formula(step.products), (
        _formula(step.reactants), _formula(step.products)
    )


# ══ FIX 10 — Arrhenius honesty. ANNOTATION ONLY: pin that nothing moved. ════

_R_KJ = 8.314


def _half_life(prefactor, ea_kj, temp_k):
    # PyYAML parses "1.0e11" (no exponent sign) as a STRING, which is how
    # src.barrier_constants reads it too (it calls float()). Mirror that.
    k = float(prefactor) * math.exp(-float(ea_kj) * 1000.0 / (_R_KJ * float(temp_k)))
    return math.log(2.0) / k


def _squash(text: str) -> str:
    """Collapse whitespace so a phrase can be asserted across a wrapped comment."""
    return " ".join(text.replace("#", " ").split())


@pytest.fixture(scope="module")
def arrhenius():
    with open(ROOT / "data/lit/arrhenius_params.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)["arrhenius_data"]


def test_wave_i_changed_no_arrhenius_runtime_value(arrhenius):
    """Fix 10 was explicitly annotate-only. If any of these move, the flags in
    the YAML and in src/barrier_constants.py are quoting arithmetic for numbers
    that no longer exist."""
    from src.barrier_constants import FAST_BARRIERS

    assert arrhenius["amadori"]["Ea_kj_mol"] == 59.0
    assert float(arrhenius["amadori"]["A_value"]) == 1.0e11
    assert arrhenius["thiol_addition"]["Ea_kj_mol"] == 29.0
    assert math.isnan(float(arrhenius["thiol_addition"]["A_value"]))
    assert arrhenius["beta_elimination_dha"]["Ea_kj_mol"] == 79.496
    assert math.isnan(float(arrhenius["beta_elimination_dha"]["A_value"]))
    assert arrhenius["schiff_condensation"]["Ea_kj_mol"] == 97.0
    assert float(arrhenius["schiff_condensation"]["A_value"]) == 1.5e11

    assert FAST_BARRIERS["schiff_condensation"][0] == 15.0
    assert FAST_BARRIERS["amadori_rearrangement"][0] == 23.0
    assert FAST_BARRIERS["lipid_homolysis"][0] == 42.0


def test_amadori_pair_really_does_give_a_sub_millisecond_half_life(arrhenius):
    """The flag claims ~0.13 ms at 150 C. Recompute it rather than trust it."""
    t_half = _half_life(
        arrhenius["amadori"]["A_value"], arrhenius["amadori"]["Ea_kj_mol"], 423.15
    )
    assert t_half == pytest.approx(1.33e-4, rel=0.05), t_half
    assert t_half < 1e-3, "the Amadori rearrangement takes minutes to hours in reality"


def test_lipid_homolysis_is_kinetically_dead_on_cooking_timescales():
    """42 kcal/mol through this repo's own Arrhenius helper: ~1.1e3 years at 150 C."""
    from src.barrier_constants import arrhenius_rate_constant

    k = arrhenius_rate_constant(42.0, 423.15, "lipid_homolysis")
    t_half_years = math.log(2.0) / k / 3.15576e7
    assert t_half_years == pytest.approx(1086.0, rel=0.05), t_half_years
    assert t_half_years > 100.0


def test_schiff_amadori_ordering_is_inverted_between_the_two_tables(arrhenius):
    """The inversion is REAL and UNFIXED. This test pins its existence so that
    'fixing' one table without the other cannot pass silently: the day the two
    lanes agree, this test fails and must be deleted with a note saying which
    table won."""
    from src.barrier_constants import arrhenius_rate_constant

    fast_schiff = math.log(2.0) / arrhenius_rate_constant(15.0, 423.15, "schiff_condensation")
    fast_amadori = math.log(2.0) / arrhenius_rate_constant(23.0, 423.15, "amadori_rearrangement")
    assert fast_schiff < fast_amadori, "FAST_BARRIERS: Schiff is the faster step"

    yaml_schiff = _half_life(
        arrhenius["schiff_condensation"]["A_value"],
        arrhenius["schiff_condensation"]["Ea_kj_mol"],
        423.15,
    )
    yaml_amadori = _half_life(
        arrhenius["amadori"]["A_value"], arrhenius["amadori"]["Ea_kj_mol"], 423.15
    )
    assert yaml_amadori < yaml_schiff, "arrhenius_params.yml: Amadori is the faster step"


@pytest.mark.parametrize(
    "key", ["amadori", "thiol_addition", "beta_elimination_dha"]
)
def test_broken_arrhenius_pairs_carry_a_loud_wave_i_flag(arrhenius, key):
    flag = arrhenius[key].get("audit_flag", "")
    assert "2026-08-27 (Wave I fix 10" in flag, f"{key} has no Wave I audit flag"
    assert "VALUE UNCHANGED" in flag or "VALUES UNCHANGED" in flag


def test_ordering_inversion_is_flagged_in_both_tables():
    yml = (ROOT / "data/lit/arrhenius_params.yml").read_text(encoding="utf-8")
    src = (ROOT / "src/barrier_constants.py").read_text(encoding="utf-8")
    for text, where in ((yml, "arrhenius_params.yml"), (src, "barrier_constants.py")):
        squashed = _squash(text)
        assert "ORDERING INVERSION" in squashed, f"{where} does not flag the inversion"
        assert "OPEN OWNER ITEM" in squashed, f"{where} does not mark the inversion unresolved"


# ══ FIX 10b — the lipid-oxidation cap docstrings ═══════════════════════════

def test_lipid_oxidation_cap_actually_ships_disabled():
    import json

    import src.lipid_oxidation as lox

    calib = json.loads(
        (ROOT / "data/lit/lipid_oxidation_calibration.json").read_text(encoding="utf-8")
    )
    assert calib["kinetics"]["max_conversion_fraction"] is None
    assert lox._kinetics()["max_conversion_fraction"] is None
    # Disabled means the multiplier is exactly 1.0, i.e. the unbounded linear model.
    assert lox._conversion_factor(1.0, 1e6, None) == 1.0


def test_lipid_oxidation_docstrings_do_not_claim_the_cap_is_active():
    """Both docstrings used to assert "the hydroperoxide extent is capped at
    physical 100% conversion". It is not, as shipped."""
    import src.lipid_oxidation as lox

    module_doc = lox.__doc__ or ""
    fn_doc = lox.predict_hexanal_generation.__doc__ or ""

    lying = re.compile(r"extent is capped at physical 100% conversion(?! \(S27\)\".)")
    for doc, where in ((module_doc, "module"), (fn_doc, "predict_hexanal_generation")):
        # The phrase may only appear as a QUOTATION of the old false claim.
        for line in doc.splitlines():
            if lying.search(line):
                assert '"' in line or "used to" in doc or "CORRECTED" in doc, (
                    f"{where} docstring still asserts the cap is active"
                )
        assert "SHIPPED DISABLED" in doc or "ships DISABLED" in doc, (
            f"{where} docstring does not say the cap ships disabled"
        )

    # The candour is the point: the headline-preservation clause must survive.
    squashed = _squash(module_doc)
    assert "so the headline stays 37/48" in squashed, (
        "the module docstring must quote the JSON rationale INCLUDING its "
        "headline-preservation sentence, not paraphrase it away"
    )
    assert "regressed the headline 37/48 -> 31/48" in squashed
