"""
src/sensory.py

Phase C: Full Sensory Prediction Model.
Maps predicted volatile concentrations to aroma descriptors using a psychophysical mixing model.
"""

import logging
import re
import yaml
from pathlib import Path
from typing import Dict, List, Any, Set, Tuple, Optional
from src.matrix_correction import ProteinType, resolve_compound_matrix_retention, resolve_matrix_correction
from src.headspace import HeadspaceModel  # noqa: E402
from src.text_utils import normalize_compound_name


ROOT = Path(__file__).resolve().parents[1]

logger = logging.getLogger(__name__)


# Odour-activity value (conc / ODT) at which the near-threshold onset ramp reaches 1.0.
# Below this the Stevens power law is scaled by a soft ramp so that intensity rises
# continuously from 0 at OAV = 1 instead of jumping to ~1.0 the instant the threshold is
# crossed (see SensoryPredictor._threshold_onset_ramp). One doubling of the threshold is
# deliberately narrow: it removes the discontinuity the optimizer trips over without
# touching any behaviour above OAV = 2.
THRESHOLD_ONSET_BAND_OAV = 2.0

# Absolute anchor for export_qda_profile: the category intensity that is reported as a
# QDA score of 10. Intensity 10 corresponds to a single compound at OAV = 100 under the
# 0.5 Stevens exponent (i.e. two orders of magnitude above its odour threshold), which is
# the top of the range this model is ever exercised over. It is a *fixed* reference, so
# two different formulations can be compared against each other on the exported scale.
QDA_FULL_SCALE_INTENSITY = 10.0


def _token_match(text: Any, keyword: str) -> bool:
    """
    True when `keyword` occurs as a whole token (or token sequence) inside `text`.

    Descriptor strings in the species YAMLs are free text ("Intense meaty",
    "meaty base note", "Cooked potato"), so exact list membership -- the old test --
    silently missed most of them. Token matching makes 'intense meaty' satisfy 'meaty'
    while still refusing spurious substring hits such as 'creamy' for 'ream'.
    """
    haystack = re.sub(r"[^a-z0-9]+", " ", str(text).lower())
    needle = re.sub(r"[^a-z0-9]+", " ", str(keyword).lower()).strip()
    if not needle:
        return False
    pattern = r"(?<![a-z0-9])" + re.escape(needle) + r"(?![a-z0-9])"
    return re.search(pattern, haystack) is not None


def _alias_candidates(name: str) -> List[str]:
    """
    Explicit aliases for a canonical species name.

    Handles the repo's 'Common name (ABBREV)' convention, including nested parentheses:
      '2-Furfurylthiol (FFT)'            -> ['2-Furfurylthiol', 'FFT']
      'Methional (3-(methylthio)propanal)' -> ['Methional', '3-(methylthio)propanal']
    These are *declared* aliases, not fuzzy prefixes: 'dimethylpyrazine' is deliberately
    not an alias of anything, because it is ambiguous between 2,3- and 2,5-.
    """
    text = str(name).strip()
    if not text.endswith(")"):
        return []
    depth = 0
    open_idx = -1
    for idx in range(len(text) - 1, -1, -1):
        ch = text[idx]
        if ch == ")":
            depth += 1
        elif ch == "(":
            depth -= 1
            if depth == 0:
                open_idx = idx
                break
    if open_idx <= 0:
        return []
    head = text[:open_idx].strip()
    inner = text[open_idx + 1:-1].strip()
    return [a for a in (head, inner) if a]

class SensoryDatabase:
    """
    Unified database for aroma compounds, off-flavours, and toxic markers.
    Loads from YAML files and provides normalized ODT lookups.
    """
    
    def __init__(self, data_dir: str = "data/species"):
        data_path = Path(data_dir)
        if not data_path.is_absolute():
            data_path = ROOT / data_path
        self.data_dir = data_path
        self.compounds = {}  # key: name, value: data
        self.smiles_map = {}
        self.tags = {}
        self.chemical_to_descriptor = {} # key: smiles, value: {odt_ppb, descriptor}
        self.alias_map = {}  # key: normalized alias, value: entry
        self.ambiguous_aliases = set()  # aliases claimed by >1 compound; never resolved
        self.unresolved_identifiers = set()  # identifiers find_entry could not resolve

        self._load_all()
        self._build_alias_index()

    def _load_all(self):
        # 1. Load targets
        files = ["desirable_targets.yml", "off_flavour_targets.yml", "toxic_markers.yml"]
        for fname in files:
            fpath = self.data_dir / fname
            if not fpath.exists():
                continue
            
            with open(fpath, "r") as f:
                data = yaml.safe_load(f)
                for c in data.get("compounds", []):
                    name = c.get("name")
                    # Convert ug/kg (ppb) to ppm
                    threshold_ppb = c.get("odour_threshold_ug_per_kg")
                    threshold_ppm = threshold_ppb / 1000.0 if threshold_ppb is not None else None
                    
                    entry = {
                        "name": name,
                        "threshold_ppm": threshold_ppm,
                        "descriptors": [d.strip().lower() for d in c.get("sensory_desc", "").split(",")],
                        "smiles": c.get("smiles"),
                        "priority": c.get("priority", "medium"),
                        "source": fname
                    }
                    self.compounds[name] = entry
                    if entry["smiles"]:
                        can = entry["smiles"]
                        try:
                            from rdkit import Chem
                            m = Chem.MolFromSmiles(can)
                            if m: can = Chem.MolToSmiles(m)
                        except ImportError:
                            pass

                        self.smiles_map[entry["smiles"]] = entry
                        self.smiles_map[can] = entry
                        
                        desc_dict = {
                            "odt": threshold_ppb,
                            "descriptor": entry["descriptors"][0] if entry["descriptors"] else "unknown"
                        }
                        self.chemical_to_descriptor[entry["smiles"]] = desc_dict
                        self.chemical_to_descriptor[can] = desc_dict

        # 2. Load tags (radar categories)
        tags_path = self.data_dir / "sensory_tags.yml"
        if tags_path.exists():
            with open(tags_path, "r") as f:
                self.tags = yaml.safe_load(f).get("tags", {})

    def _build_alias_index(self) -> None:
        """
        Build the normalized name/alias -> entry index used by find_entry.

        Aliases are declared, not guessed. An alias claimed by more than one compound is
        dropped into `ambiguous_aliases` and never resolves, so a caller asking for
        'dimethylpyrazine' gets None rather than whichever of 2,3- / 2,5- happened to be
        first in dict order.
        """
        claims: Dict[str, Set[str]] = {}

        def claim(alias: str, name: str) -> None:
            key = normalize_compound_name(alias)
            if not key:
                return
            claims.setdefault(key, set()).add(name)

        for name in self.compounds:
            claim(name, name)
        for name in self.compounds:
            for alias in _alias_candidates(name):
                claim(alias, name)

        for key, owners in claims.items():
            if len(owners) == 1:
                self.alias_map[key] = self.compounds[next(iter(owners))]
            else:
                self.ambiguous_aliases.add(key)

    def find_entry(self, identifier: str) -> Optional[Dict]:
        """
        Resolve an identifier to a compound entry by exact name, SMILES, or declared alias.

        Matching is exact after normalization (case/punctuation insensitive) plus the
        parenthetical aliases built in `_build_alias_index`. There is deliberately NO
        fuzzy/substring fallback: the old one resolved any partial name to the first
        compound whose name merely *contained* it and then filed the caller's
        concentration under that arbitrary compound's name. Unknown identifiers now
        return None (warned once each) instead of silently becoming another compound.
        """
        if identifier in self.compounds:
            return self.compounds[identifier]
        if identifier in self.smiles_map:
            return self.smiles_map[identifier]

        target_norm = normalize_compound_name(identifier)
        if not target_norm:
            return None

        entry = self.alias_map.get(target_norm)
        if entry is not None:
            return entry

        for smiles, smiles_entry in self.smiles_map.items():
            if target_norm == normalize_compound_name(smiles):
                return smiles_entry

        if identifier not in self.unresolved_identifiers:
            self.unresolved_identifiers.add(identifier)
            reason = (
                "ambiguous (matches more than one compound)"
                if target_norm in self.ambiguous_aliases
                else "not found in the sensory database"
            )
            logger.warning(
                "sensory: identifier %r is %s; it will not be scored as an aroma compound.",
                identifier,
                reason,
            )
        return None


class SensoryPredictor:
    """
    Predicts sensory profiles using psychophysical mixing.
    """
    
    def __init__(self, database: Optional[SensoryDatabase] = None, headspace: Optional[HeadspaceModel] = None):
        self.db = database or SensoryDatabase()
        self.headspace = headspace
        # Stevens-law exponent for perceived odour intensity, I = OAV^exponent.
        # HONEST PROVENANCE: 0.5 is NOT fitted here. Published Stevens exponents for
        # odorants scatter roughly over 0.2-0.7 depending on compound and method; 0.5 sits
        # inside that range and is the value the unit tests and benchmark calibration were
        # written against, so it was adopted as the shipped default rather than derived.
        # (A stale comment elsewhere in this file used to claim 0.55; there is and was no
        # 0.55 anywhere in the code -- the docstring was simply wrong.)
        self.exponent = 0.5
        # Per-call record of compounds that were NOT scored, and why. Populated by
        # predict_profile; consumers can surface it instead of silently losing compounds.
        self.unscored_compounds: List[Dict[str, Any]] = []
        self._family_cache: Dict[str, Set[str]] = {}
        self.synergy_pow = 1.3  # Group synergy factor (1.0 = additive, 2.0 = vector-sum)
        
        # Expert synergy pairs (Hofmann & Schieberle 2000)
        # (Compound A, Compound B) -> multiplier
        self.synergy_pairs = {
            ("2-furfurylthiol", "methional"): 1.8, # Meaty synergy
            ("2-methyl-3-furanthiol", "methional"): 2.2, # Stronger meaty synergy
            ("2,5-dimethylpyrazine", "pyrazine"): 1.3, # Roasted synergy boost
            ("4-hydroxy-2,5-dimethyl-3(2h)-furanone", "furaneol"): 1.2 # caramel-sweet synergy
        }

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _resolve_odt_ppb(entry: Dict[str, Any]) -> Optional[float]:
        """
        Odour threshold in ppb for an entry, or None when the entry has no threshold data.

        Returning None is meaningful: the compound is NOT scoreable as an aroma (this is
        the case for every toxic marker, which have no `odour_threshold_ug_per_kg` field).
        The previous code read a non-existent `threshold_ppb` key with a 0.1 default, so
        every threshold-less compound was silently handed a *highly potent* fabricated
        threshold; and the `if entry.get("threshold_ppm")` guard that was supposed to
        override it also dropped a legitimate 0.0. Both are fixed by explicit None tests.
        """
        ppm = entry.get("threshold_ppm")
        if ppm is not None:
            return float(ppm) * 1000.0
        ppb = entry.get("threshold_ppb")
        if ppb is not None:
            return float(ppb)
        return None

    @staticmethod
    def _threshold_onset_ramp(oav: float) -> float:
        """
        Soft onset factor in [0, 1] for odour-activity values just above threshold.

        Without it the model is discontinuous at every ODT: OAV = 1 scores 0 (the
        compound is dropped) while OAV = 1 + eps scores ~1.0, which makes the optimizer
        objective jump at every threshold crossing. The ramp is logarithmic in OAV and
        reaches 1.0 at THRESHOLD_ONSET_BAND_OAV, so intensity now rises continuously from
        0 at threshold, and NOTHING above the band changes.

        It also matches the psychophysics better than a step: detection near the ODT is
        probabilistic (the ODT is by definition the 50%-detection point), so a partial
        perceived intensity in the first doubling is the physically honest reading.
        """
        import math

        if oav <= 1.0:
            return 0.0
        if oav >= THRESHOLD_ONSET_BAND_OAV:
            return 1.0
        return math.log(oav) / math.log(THRESHOLD_ONSET_BAND_OAV)

    def _family_members(self, keyword: str) -> Set[str]:
        """
        Canonical names of every database compound belonging to an aroma family.

        Membership is the union of (a) the radar tag list in sensory_tags.yml and
        (b) a token-aware match of the keyword against the compound name and its free-text
        descriptors. Both the family *total* and the family the mask is *applied to* now
        use this one predicate. Previously the total used
        `keyword in name.lower() or keyword in descriptors` (substring on the name,
        exact list membership on descriptors) while the application used descriptor
        membership alone -- so e.g. MFT ("Intense meaty, boiled meat, broth-like") was
        never counted and never masked.
        """
        key = str(keyword).lower().strip()
        if key in self._family_cache:
            return self._family_cache[key]

        members: Set[str] = set()
        for display_name in self.db.tags.get(key, []) or []:
            entry = self.db.find_entry(display_name)
            if entry:
                members.add(entry["name"])
        for name, entry in self.db.compounds.items():
            if _token_match(name, key) or any(_token_match(d, key) for d in entry.get("descriptors", [])):
                members.add(name)

        self._family_cache[key] = members
        return members

    def predict_profile(self, 
                        concentration_dict_ppb: Dict[str, float], 
                        protein_type: str = "free",
                        denaturation_state: float = 0.5,
                        temp_c: Optional[float] = None,
                        fat_fraction: float = 0.0,
                        protein_fraction: float = 0.0,
                        extrusion_process: Optional[Dict[str, object]] = None) -> Dict[str, Tuple[float, float]]:
        """
        Calculate perceived intensity for each compound using Stevens' Power Law.
        Returns {name: (intensity, intensity_uncertainty)}

        Includes matrix ODT corrections, a continuous near-threshold onset, and
        beany -> meaty perceptual masking. Compounds that cannot be scored (unknown
        identifier, or no odour threshold on record) are omitted from the result and
        recorded on `self.unscored_compounds` with a reason.
        """
        p_type = ProteinType(protein_type)
        # Use matrix volatile retention factor to adjust ODT
        bulk_retention = resolve_matrix_correction(p_type, denaturation_state).volatile_retention

        # 1. Headspace Partitioning (optional) - convert ppb to ppm for old model if needed
        # Actually let's stay in ppb for consistency with new safety/benchmarks.
        if self.headspace and temp_c is not None:
            effective_concs = self.headspace.predict_headspace(
                concentration_dict_ppb,
                temp_c,
                fat_fraction,
                protein_fraction,
                protein_type=protein_type,
                denaturation_state=denaturation_state,
                extrusion_process=extrusion_process,
            )
        else:
            effective_concs = concentration_dict_ppb

        # 2. Perceived Intensity calculation
        results = {}
        self.unscored_compounds = []
        for compound, conc_data in effective_concs.items():
            if isinstance(conc_data, tuple):
                conc, unc_ratio = conc_data
            else:
                conc, unc_ratio = conc_data, 0.25

            entry = self.db.find_entry(compound)
            if entry is None:
                self.unscored_compounds.append({
                    "identifier": compound,
                    "reason": "unresolved_identifier",
                    "note": "no exact name/SMILES/alias match in the sensory database",
                })
                continue

            odt_base = self._resolve_odt_ppb(entry)
            if odt_base is None or odt_base <= 0.0:
                # No threshold data => not scoreable as an aroma. This is the honest
                # outcome for the toxic markers, which are reported by the safety layer
                # on a concentration basis; inventing a potent default ODT for them
                # would fabricate a huge aroma intensity out of nothing.
                self.unscored_compounds.append({
                    "identifier": compound,
                    "name": entry["name"],
                    "reason": "no_odour_threshold",
                    "note": "compound has no odour threshold; not scoreable as an aroma",
                })
                continue

            # Matrix Effect: ODT_matrix = ODT_water / Retention
            retention = resolve_compound_matrix_retention(
                compound,
                protein_type=p_type,
                denaturation_state=denaturation_state,
            ) if protein_type != "free" else bulk_retention
            odt_matrix = odt_base / max(0.01, retention)

            # Stevens' Power Law with a continuous near-threshold onset:
            #   I = OAV^exponent * ramp(OAV),  ramp -> 0 at OAV = 1, ramp = 1 above the band.
            oav = conc / odt_matrix
            if oav > 1.0:
                intensity = pow(oav, self.exponent) * self._threshold_onset_ramp(oav)
                if intensity <= 0.0:
                    continue
                i_unc = intensity * self.exponent * unc_ratio
                results[entry["name"]] = (intensity, i_unc)

        # 3. Perceptual Masking
        # Family membership uses one consistent, token-aware predicate for both the
        # aggregate totals and the set the mask is applied to (see _family_members).
        masked_results = {k: v for k, v in results.items()}

        beany_members = self._family_members("beany")
        meaty_members = self._family_members("meaty")

        beany_total = sum(v[0] for k, v in results.items() if k in beany_members)
        meaty_total = sum(v[0] for k, v in results.items() if k in meaty_members)

        if beany_total > 0 and meaty_total > 0:
            # Masking coefficient: meaty is reduced by beany presence
            k_mask = 0.35 # 35% max reduction
            mask_factor = 1.0 - k_mask * (beany_total / (beany_total + meaty_total))
            for k, (v, u) in list(masked_results.items()):
                if k in meaty_members and k not in beany_members:
                    masked_results[k] = (v * mask_factor, u * mask_factor)

        return masked_results

    def export_qda_profile(self, intensity_results: Dict[str, Tuple[float, float]]) -> Dict[str, float]:
        """
        Map category intensities onto an ABSOLUTE 0-10 QDA-style scale.

        The scale is anchored to the fixed reference QDA_FULL_SCALE_INTENSITY: a category
        intensity of 10 (one compound two orders of magnitude above its ODT) reports 10,
        half that reports 5, and scores are clipped at 10.

        This used to divide by the profile's own maximum, which meant the top note was
        ALWAYS exactly 10.0 no matter how weak the formulation, while still being labeled
        a 0-10 QDA scale -- an absolute-sounding label on a purely relative quantity, and
        two formulations were not comparable. Weak profiles now report weak numbers.
        """
        radar = self.get_radar_data_from_intensities(intensity_results)
        return {
            k: min(10.0, (v[0] / QDA_FULL_SCALE_INTENSITY) * 10.0)
            for k, v in radar.items()
            if v[0] > 0
        }

    def get_radar_data_from_intensities(self, compound_intensities: Dict[str, Tuple[float, float]]) -> Dict[str, Tuple[float, float]]:
        """
        Helper to group already calculated intensities into radar categories.
        Returns {category: (score, uncertainty)}
        """
        # 1. Group intensities by category
        radar = {tag: (0.0, 0.0) for tag in self.db.tags.keys()}
        
        def normalize(s):
            import re
            return re.sub(r"[^a-z0-9]+", "", s.lower())

        for tag, search_names in self.db.tags.items():
            intensities = []
            uncertainties = []
            matched_compounds = set()
            
            norm_search_names = [normalize(sn) for sn in search_names]
            
            for c_name, (intensity, unc) in compound_intensities.items():
                c_norm = normalize(c_name)
                for sn_norm in norm_search_names:
                    if sn_norm in c_norm and c_name not in matched_compounds:
                        intensities.append(intensity)
                        uncertainties.append(unc)
                        matched_compounds.add(c_name)
                        break
            
            if not intensities:
                continue

            # Identify which compounds are active for this tag to check for synergy
            active_for_tag = [c.lower() for c in matched_compounds]
            
            # Specialized synergy check.
            # A synergy pair requires TWO DISTINCT compounds. The old test asked only
            # whether some active name contained `a` and some active name contained `b`,
            # which the same compound could satisfy on its own: 2,5-dimethylpyrazine
            # contains both "2,5-dimethylpyrazine" and "pyrazine", so it synergized with
            # itself for a free +30% on the roasted axis. The HDMF/furaneol pair has the
            # same latent shape.
            synergy_boost = 1.0
            for (a, b), boost in self.synergy_pairs.items():
                matches_a = {name for name in active_for_tag if a in name}
                matches_b = {name for name in active_for_tag if b in name}
                if matches_a and matches_b and len(matches_a | matches_b) >= 2:
                    synergy_boost = max(synergy_boost, boost)
            
            if tag.lower() in ["meaty", "sulfury", "sulfurous"]:
                # High additivity for meaty/sulfur (β=1.1)
                group_sum = sum([pow(i, 1.1) for i in intensities])
                score = pow(group_sum, 1.0/1.1) * synergy_boost
            elif tag.lower() == "roasted":
                # Moderate additivity for roasted (β=1.2)
                group_sum = sum([pow(i, 1.2) for i in intensities])
                score = pow(group_sum, 1.0/1.2) * synergy_boost
            else:
                # Standard partial additivity (β=1.3)
                group_sum = sum([pow(i, self.synergy_pow) for i in intensities])
                score = pow(group_sum, 1.0/self.synergy_pow)
            
            total_unc = sum([u**2 for u in uncertainties])**0.5
            radar[tag] = (score, total_unc)
                        
        return radar

    def get_radar_data(self, 
                       concentration_dict_ppb: Dict[str, float],
                       protein_type: str = "free",
                       temp_c: Optional[float] = None,
                       fat_fraction: float = 0.0,
                       protein_fraction: float = 0.0,
                       extrusion_process: Optional[Dict[str, object]] = None) -> Dict[str, Tuple[float, float]]:
        """
        High-level entry for radar categories.
        """
        compound_intensities = self.predict_profile(
            concentration_dict_ppb, 
            protein_type=protein_type,
            temp_c=temp_c, 
            fat_fraction=fat_fraction, 
            protein_fraction=protein_fraction,
            extrusion_process=extrusion_process,
        )
        return self.get_radar_data_from_intensities(compound_intensities)

    def get_dominant_notes(self, radar_profile: Dict[str, Tuple[float, float]], top_n: int = 3) -> List[tuple]:
        """Return the top categories by score."""
        sorted_notes = sorted(radar_profile.items(), key=lambda x: x[1][0], reverse=True)
        return sorted_notes[:top_n]

