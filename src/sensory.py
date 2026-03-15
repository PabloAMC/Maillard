"""
src/sensory.py

Phase C: Full Sensory Prediction Model.
Maps predicted volatile concentrations to aroma descriptors using a psychophysical mixing model.
"""

import yaml
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional
from src.matrix_correction import ProteinType, MATRIX_CORRECTIONS
from src.headspace import HeadspaceModel  # noqa: E402

class SensoryDatabase:
    """
    Unified database for aroma compounds, off-flavours, and toxic markers.
    Loads from YAML files and provides normalized ODT lookups.
    """
    
    def __init__(self, data_dir: str = "data/species"):
        self.data_dir = Path(data_dir)
        self.compounds = {}  # key: name, value: data
        self.smiles_map = {}
        self.tags = {}
        self.chemical_to_descriptor = {} # key: smiles, value: {odt_ppb, descriptor}
        
        self._load_all()

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

    def find_entry(self, identifier: str) -> Optional[Dict]:
        """Lookup by name or SMILES with robustness to underscores/spaces."""
        if identifier in self.compounds:
            return self.compounds[identifier]
        if identifier in self.smiles_map:
            return self.smiles_map[identifier]
        
        def normalize(s):
            import re
            return re.sub(r"[^a-z0-9]+", "", s.lower())

        target_norm = normalize(identifier)
        if not target_norm:
             return None

        # Robust match by normalized name
        for name, entry in self.compounds.items():
            if target_norm == normalize(name):
                return entry
            if target_norm in normalize(name): # Fuzzy fallback
                return entry
        return None


class SensoryPredictor:
    """
    Predicts sensory profiles using psychophysical mixing.
    """
    
    def __init__(self, database: Optional[SensoryDatabase] = None, headspace: Optional[HeadspaceModel] = None):
        self.db = database or SensoryDatabase()
        self.headspace = headspace
        # Perceived intensity scales with exponent ~0.5 (Stevens' Law).
        # Historically this was set to 1.0 for linear scaling in formulation design,
        # but unit tests and benchmark calibration expect 0.5.
        self.exponent = 0.5
        self.synergy_pow = 1.3  # Group synergy factor (1.0 = additive, 2.0 = vector-sum)
        
        # Expert synergy pairs (Hofmann & Schieberle 2000)
        # (Compound A, Compound B) -> multiplier
        self.synergy_pairs = {
            ("2-furfurylthiol", "methional"): 1.8, # Meaty synergy
            ("2-methyl-3-furanthiol", "methional"): 2.2, # Stronger meaty synergy
            ("2,5-dimethylpyrazine", "pyrazine"): 1.3, # Roasted synergy boost
            ("4-hydroxy-2,5-dimethyl-3(2h)-furanone", "furaneol"): 1.2 # caramel-sweet synergy
        }

    def predict_profile(self, 
                        concentration_dict_ppb: Dict[str, float], 
                        protein_type: str = "free",
                        temp_c: Optional[float] = None,
                        fat_fraction: float = 0.0,
                        protein_fraction: float = 0.0) -> Dict[str, Tuple[float, float]]:
        """
        Calculate perceived intensity for each compound using Stevens' Power Law.
        Returns {name: (intensity, intensity_uncertainty)}
        
        Updated with Matrix ODT corrections and Perceptual Masking.
        """
        import math
        
        p_type = ProteinType(protein_type)
        # Use matrix volatile retention factor to adjust ODT
        retention = MATRIX_CORRECTIONS[p_type].volatile_retention

        # 1. Headspace Partitioning (optional) - convert ppb to ppm for old model if needed
        # Actually let's stay in ppb for consistency with new safety/benchmarks.
        if self.headspace and temp_c is not None:
            effective_concs = self.headspace.predict_headspace(
                concentration_dict_ppb, temp_c, fat_fraction, protein_fraction
            )
        else:
            effective_concs = concentration_dict_ppb

        # 2. Perceived Intensity calculation
        results = {}
        for compound, conc_data in effective_concs.items():
            if isinstance(conc_data, tuple):
                conc, unc_ratio = conc_data
            else:
                conc, unc_ratio = conc_data, 0.25
            
            entry = self.db.find_entry(compound)
            if entry:
                # ODT from DB is usually in ppb. 
                # If unavailable, fallback to 0.1 ppb
                odt_base = entry.get("threshold_ppb", 0.1) 
                if entry.get("threshold_ppm"):
                    odt_base = entry["threshold_ppm"] * 1000.0
                
                # Matrix Effect: ODT_matrix = ODT_water / Retention
                odt_matrix = odt_base / max(0.01, retention)
                
                # Stevens' Power Law: I = (C / ODT)^0.55
                if conc > odt_matrix:
                    intensity = pow(conc / odt_matrix, self.exponent)
                    i_unc = intensity * self.exponent * unc_ratio
                    results[entry["name"]] = (intensity, i_unc)

        # 3. Perceptual Masking (Fix 4)
        # We group results by descriptor for masking logic
        masked_results = {k: v for k, v in results.items()}
        
        # Determine aggregate intensity for masking families
        def get_family_intensity(keyword):
            return sum(v[0] for k, v in results.items() if keyword in k.lower() or keyword in self.db.find_entry(k)["descriptors"])

        beany_total = get_family_intensity("beany")
        meaty_total = get_family_intensity("meaty")
        
        if beany_total > 0 and meaty_total > 0:
            # Masking coefficient: meaty is reduced by beany presence
            k_mask = 0.35 # 35% max reduction
            mask_factor = 1.0 - k_mask * (beany_total / (beany_total + meaty_total))
            for k, (v, u) in masked_results.items():
                entry = self.db.find_entry(k)
                if "meaty" in entry["descriptors"]:
                    masked_results[k] = (v * mask_factor, u * mask_factor)

        return masked_results

    def export_qda_profile(self, intensity_results: Dict[str, Tuple[float, float]]) -> Dict[str, float]:
        """
        Normalizes intensities to a 0-10 QDA scale for export.
        """
        import math
        # Aggregated categories
        categories = {}
        # This is similar to get_radar_data but returns 0-10 scores
        radar = self.get_radar_data_from_intensities(intensity_results)
        
        # Normalize radar scores to 10
        max_score = max([v[0] for v in radar.values()]) if radar else 1.0
        qda = {k: min(10.0, (v[0] / max_score) * 10.0) for k, v in radar.items() if v[0] > 0}
        return qda

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
            
            # Specialized synergy check
            synergy_boost = 1.0
            for (a, b), boost in self.synergy_pairs.items():
                has_a = any(a in name for name in active_for_tag)
                has_b = any(b in name for name in active_for_tag)
                if has_a and has_b:
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
                       protein_fraction: float = 0.0) -> Dict[str, Tuple[float, float]]:
        """
        High-level entry for radar categories.
        """
        compound_intensities = self.predict_profile(
            concentration_dict_ppb, 
            protein_type=protein_type,
            temp_c=temp_c, 
            fat_fraction=fat_fraction, 
            protein_fraction=protein_fraction
        )
        return self.get_radar_data_from_intensities(compound_intensities)

    def get_dominant_notes(self, radar_profile: Dict[str, Tuple[float, float]], top_n: int = 3) -> List[tuple]:
        """Return the top categories by score."""
        sorted_notes = sorted(radar_profile.items(), key=lambda x: x[1][0], reverse=True)
        return sorted_notes[:top_n]


if __name__ == "__main__":
    predictor = SensoryPredictor()
    # Mock concentration profile in ppm (1 ppb = 0.001 ppm)
    mock_conc = {
        "2-Furfurylthiol (FFT)": 0.01, # 10 ppb
        "Methional": 0.005,           # 5 ppb
        "Hexanal": 0.05               # 50 ppb
    }
    
    radar = predictor.get_radar_data(mock_conc)
    print("Radar Category Scores (± confidence):")
    for note, (score, unc) in predictor.get_dominant_notes(radar):
        print(f"  - {note}: {score:.4f} ± {unc:.4f}")
