"""
src/pathway_extractor.py — Extract elementary steps from RMG-Py output.

Parses RMG's Chemkin-format reaction outputs and adjacency lists (species dictionaries)
to generate a list of elementary reaction steps suitable for Tier 1 xTB screening.

Audit remediation 2026-08-27 (item 3.x, extractor failure modes):
  * species that cannot be resolved no longer vanish from the step — the step is
    marked `source_quality="unresolved_species"` and excluded from energy
    computation by the ranker;
  * the equation splitter understands `<=>` / `=>` and tolerates `=` inside
    species labels (RMG labels are frequently SMILES-like, e.g. `C=CC=O`);
  * `_parse_species_dictionary` parses real RMG adjacency lists instead of
    silently returning an empty dictionary, and raises a clear
    `UnsupportedSpeciesFormatError` when a block is in neither supported form.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import Descriptors
from functools import cached_property


class UnsupportedSpeciesFormatError(ValueError):
    """Raised when a species_dictionary block is in no format we can parse."""


@dataclass
class Species:
    label: str
    smiles: str

    @cached_property
    def is_volatile(self) -> bool:
        """
        Scientific heuristic for aroma volatility.
        Volatiles are generally small (MW < 160) and lack excessive polarity
        (H-bond donors <= 1).
        """
        if not self.smiles:
            return False

        # Inerts and precursors that aren't aroma compounds
        inerts = {"water", "h2o", "h2", "co2", "ammonia", "h2s"}
        if self.label.lower() in inerts or self.smiles.lower() in inerts:
            return False

        try:
            mol = Chem.MolFromSmiles(self.smiles)
            if not mol:
                return False

            mw = Descriptors.MolWt(mol)
            h_donors = Descriptors.NumHDonors(mol)

            # Thresholds aligned with empirical flavor chemistry
            return mw < 160 and h_donors <= 1
        except Exception:
            return False

@dataclass
class ElementaryStep:
    reactants: List[Species]
    products: List[Species]
    reaction_family: Optional[str] = None
    rate_constant_k: Optional[float] = None
    source_quality: str = "heuristic" # "literature", "estimated_tst", "heuristic", "unresolved_species"
    barrier_uncertainty_kcal: float = 5.0 # Default heuristic uncertainty
    barrier_kcal_mol: Optional[float] = None
    reversible: bool = False
    direction: str = "forward"  # "forward" | "reverse"
    unresolved_species: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        reacts = " + ".join([r.label for r in self.reactants])
        prods = " + ".join([p.label for p in self.products])
        arrow = "<=>" if self.reversible else "->"
        suffix = f" [UNRESOLVED: {', '.join(self.unresolved_species)}]" if self.unresolved_species else ""
        return f"{reacts} {arrow} {prods} [{self.reaction_family or 'unknown'}]{suffix}"


# --- RMG adjacency-list parsing ---------------------------------------------

_BOND_ORDERS = {
    "S": Chem.BondType.SINGLE,
    "D": Chem.BondType.DOUBLE,
    "T": Chem.BondType.TRIPLE,
    "Q": Chem.BondType.QUADRUPLE,
    "B": Chem.BondType.AROMATIC,
}

# e.g. "1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}"  (optionally with a "*1" reaction
# centre label between the index and the element symbol)
_ADJACENCY_LINE_RE = re.compile(
    r"^(?P<index>\d+)\s+"
    r"(?:\*\d+\s+)?"
    r"(?P<element>[A-Z][a-z]?)\s+"
    r"u(?P<unpaired>\d+)\s+"
    r"p(?P<lone>\d+)\s+"
    r"c(?P<charge>[+-]?\d+)"
    r"(?P<bonds>(?:\s+\{\d+,[SDTQB]\})*)\s*$"
)
_BOND_RE = re.compile(r"\{(\d+),([SDTQB])\}")
_SMILES_COMMENT_RE = re.compile(r'SMILES="([^"]+)"')


def adjacency_list_to_smiles(lines: List[str]) -> str:
    """Convert an RMG adjacency-list block into a canonical SMILES string.

    RMG adjacency lists carry every hydrogen explicitly, so implicit-H counting
    is switched off. Raises UnsupportedSpeciesFormatError for anything the
    grammar above does not cover, or for a graph RDKit refuses to sanitize.
    """
    parsed: Dict[int, Dict[str, object]] = {}
    for line in lines:
        match = _ADJACENCY_LINE_RE.match(line.strip())
        if match is None:
            raise UnsupportedSpeciesFormatError(f"Unparseable adjacency-list line: {line.strip()!r}")
        parsed[int(match.group("index"))] = {
            "element": match.group("element"),
            "unpaired": int(match.group("unpaired")),
            "charge": int(match.group("charge")),
            "bonds": _BOND_RE.findall(match.group("bonds") or ""),
        }

    if not parsed:
        raise UnsupportedSpeciesFormatError("Empty adjacency list")

    mol = Chem.RWMol()
    index_to_idx: Dict[int, int] = {}
    for index in sorted(parsed):
        spec = parsed[index]
        try:
            atom = Chem.Atom(str(spec["element"]))
        except Exception as exc:  # unknown element symbol
            raise UnsupportedSpeciesFormatError(f"Unknown element {spec['element']!r}") from exc
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(int(spec["unpaired"]))
        atom.SetFormalCharge(int(spec["charge"]))
        index_to_idx[index] = mol.AddAtom(atom)

    for index in sorted(parsed):
        for partner_str, order in parsed[index]["bonds"]:  # type: ignore[union-attr]
            partner = int(partner_str)
            if partner not in index_to_idx:
                raise UnsupportedSpeciesFormatError(
                    f"Adjacency list references undefined atom {partner}"
                )
            if partner <= index:
                continue  # each bond is listed on both atoms
            bond_type = _BOND_ORDERS[order]
            mol.AddBond(index_to_idx[index], index_to_idx[partner], bond_type)
            if bond_type == Chem.BondType.AROMATIC:
                mol.GetAtomWithIdx(index_to_idx[index]).SetIsAromatic(True)
                mol.GetAtomWithIdx(index_to_idx[partner]).SetIsAromatic(True)

    try:
        molecule = mol.GetMol()
        Chem.SanitizeMol(molecule)
        # RMG lists every hydrogen explicitly; collapse them so the SMILES
        # matches the canonical spelling used everywhere else in the repo.
        molecule = Chem.RemoveHs(molecule)
        return Chem.MolToSmiles(molecule)
    except Exception as exc:
        raise UnsupportedSpeciesFormatError(f"RDKit could not sanitize the adjacency list: {exc}") from exc


def _split_on_plus(side: str) -> List[str]:
    """Split a reaction side on '+', ignoring '+' inside square brackets.

    Charged RMG/SMILES labels such as [NH4+] must not be torn in half.
    """
    tokens: List[str] = []
    current: List[str] = []
    depth = 0
    for char in side:
        if char == "[":
            depth += 1
        elif char == "]":
            depth = max(0, depth - 1)
        if char == "+" and depth == 0:
            tokens.append("".join(current))
            current = []
            continue
        current.append(char)
    tokens.append("".join(current))
    return [token.strip() for token in tokens if token.strip()]


class PathwayExtractor:
    """Parses RMG chemkin outputs into explicit reaction objects."""

    def __init__(self, rmg_output_dir: Path, *, emit_reverse_steps: bool = False):
        self.rmg_output_dir = Path(rmg_output_dir)
        self.chemkin_dir = self.rmg_output_dir / "chemkin"
        self.species_dict: dict[str, Species] = {}
        self.elementary_steps: List[ElementaryStep] = []
        # When True, a reversible equation (A<=>B, or bare Chemkin '=') emits a
        # forward AND a reverse step. Off by default so consumers that expect
        # one step per equation keep their current behaviour; every step is
        # flagged with `reversible` either way.
        self.emit_reverse_steps = bool(emit_reverse_steps)

    def _parse_species_dictionary(self):
        """Parse species_dictionary.txt into Species objects with SMILES.

        Supports two block formats:
          1. an `// SMILES="..."` annotation comment (used by the repo fixtures);
          2. a genuine RMG adjacency list, which is converted through RDKit.
        A block in neither format raises UnsupportedSpeciesFormatError rather
        than being dropped — silently returning {} made every downstream
        reaction unresolvable while looking like a clean parse.
        """
        dict_path = self.chemkin_dir / "species_dictionary.txt"
        if not dict_path.exists():
            raise FileNotFoundError(f"Missing {dict_path}")

        raw_blocks = dict_path.read_text().split("\n\n")
        for raw_block in raw_blocks:
            lines = [line.rstrip() for line in raw_block.splitlines() if line.strip()]
            if not lines:
                continue
            self._ingest_species_block(lines)

    def _ingest_species_block(self, lines: List[str]) -> None:
        label: Optional[str] = None
        smiles: Optional[str] = None
        adjacency: List[str] = []

        for line in lines:
            stripped = line.strip()
            comment_match = _SMILES_COMMENT_RE.search(stripped)
            if comment_match:
                smiles = comment_match.group(1)
                continue
            if stripped.startswith("//") or stripped.startswith("!"):
                continue
            if stripped.lower().startswith("multiplicity"):
                continue
            if label is None:
                label = stripped
                continue
            if _ADJACENCY_LINE_RE.match(stripped):
                adjacency.append(stripped)

        if label is None:
            return

        if smiles:
            self.species_dict[label] = Species(label, smiles)
            return

        if adjacency:
            self.species_dict[label] = Species(label, adjacency_list_to_smiles(adjacency))
            return

        raise UnsupportedSpeciesFormatError(
            f"Species block {label!r} in species_dictionary.txt has neither a "
            'SMILES annotation (// SMILES="...") nor a parseable RMG adjacency list.'
        )

    # --- reaction equations -------------------------------------------------

    def _resolved_count(self, side: str) -> Tuple[int, int]:
        tokens = _split_on_plus(side)
        resolved = sum(1 for token in tokens if token in self.species_dict)
        return resolved, len(tokens)

    def _split_reaction_equation(self, equation: str) -> Optional[Tuple[str, str, bool]]:
        """Return (left, right, reversible) for a Chemkin equation.

        `=` cannot simply be split on: RMG species labels routinely contain it
        (`C=CC=O`). The separator is chosen as the `=` position whose two sides
        resolve to the most known species.
        """
        for arrow, reversible in (("<=>", True), ("=>", False)):
            if arrow in equation:
                left, right = equation.split(arrow, 1)
                return left.strip(), right.strip(), reversible

        best: Optional[Tuple[Tuple[int, int], str, str]] = None
        for position, char in enumerate(equation):
            if char != "=":
                continue
            left = equation[:position].strip()
            right = equation[position + 1:].strip()
            if not left or not right:
                continue
            left_resolved, left_total = self._resolved_count(left)
            right_resolved, right_total = self._resolved_count(right)
            if left_total == 0 or right_total == 0:
                continue
            score = (left_resolved + right_resolved, -(left_total + right_total))
            if best is None or score > best[0]:
                best = (score, left, right)

        if best is None:
            return None
        # Bare '=' is reversible in the Chemkin convention.
        return best[1], best[2], True

    def _resolve_side(self, side: str) -> Tuple[List[Species], List[str]]:
        """Resolve one side of an equation, keeping unresolved labels visible."""
        species: List[Species] = []
        unresolved: List[str] = []
        for token in _split_on_plus(side):
            known = self.species_dict.get(token)
            if known is not None:
                species.append(known)
            else:
                # Placeholder with an empty SMILES: the label is preserved so
                # the step can report exactly what it could not resolve.
                species.append(Species(token, ""))
                unresolved.append(token)
        return species, unresolved

    def _parse_chem_inp(self):
        """Parse chem_annotated.inp to extract reaction definitions."""
        inp_path = self.chemkin_dir / "chem_annotated.inp"
        if not inp_path.exists():
            raise FileNotFoundError(f"Missing {inp_path}")

        in_reactions_block = False

        with open(inp_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("REACTIONS"):
                    in_reactions_block = True
                    continue
                if line.startswith("END") and in_reactions_block:
                    break

                if in_reactions_block and line and not line.startswith("!"):
                    # Typical line: A+B=C+D  1.0E3  0.0  10.0
                    rxn_part = line.split()[0] # get just the equation part
                    if "=" not in rxn_part:
                        continue

                    split = self._split_reaction_equation(rxn_part)
                    if split is None:
                        continue
                    left, right, reversible = split

                    reactants, unresolved_r = self._resolve_side(left)
                    products, unresolved_p = self._resolve_side(right)
                    if not reactants or not products:
                        continue

                    unresolved = unresolved_r + unresolved_p
                    # Unresolved species FAIL the step instead of being dropped:
                    # dropping them silently rewrote the stoichiometry (A + B -> C
                    # became A -> C) and the resulting energies were meaningless.
                    source_quality = "unresolved_species" if unresolved else "heuristic"

                    self.elementary_steps.append(
                        ElementaryStep(
                            reactants=reactants,
                            products=products,
                            source_quality=source_quality,
                            reversible=reversible,
                            direction="forward",
                            unresolved_species=unresolved,
                        )
                    )
                    if reversible and self.emit_reverse_steps:
                        self.elementary_steps.append(
                            ElementaryStep(
                                reactants=list(products),
                                products=list(reactants),
                                source_quality=source_quality,
                                reversible=True,
                                direction="reverse",
                                unresolved_species=unresolved,
                            )
                        )

    def run(self) -> List[ElementaryStep]:
        """Execute parsing and return list of steps."""
        self._parse_species_dictionary()
        self._parse_chem_inp()
        return self.elementary_steps
