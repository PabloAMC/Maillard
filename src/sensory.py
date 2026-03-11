"""
src/sensory.py

Phase C: Full Sensory Prediction Model.
Maps predicted volatile concentrations to aroma descriptors using a psychophysical mixing model.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional
from src.headspace import HeadspaceModel

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
                    if entry["smiles"]:
                        self.smiles_map[entry["smiles"]] = entry

        # 2. Load tags (radar categories)
        tags_path = self.data_dir / "sensory_tags.yml"
        if tags_path.exists():
            with open(tags_path, "r") as f:
                self.tags = yaml.safe_load(f).get("tags", {})

    def find_entry(self, identifier: str) -> Optional[Dict]:
        """Lookup by name or SMILES."""
        if identifier in self.compounds:
            return self.compounds[identifier]
        if identifier in self.smiles_map:
            return self.smiles_map[identifier]
        
        # Fuzzy match by name
        for name, entry in self.compounds.items():
            if identifier.lower() in name.lower():
                return entry
        return None


class SensoryPredictor:
    """
    Predicts sensory profiles using psychophysical mixing.
    """
    
    def __init__(self, database: Optional[SensoryDatabase] = None, headspace: Optional[HeadspaceModel] = None):
        self.db = database or SensoryDatabase()
        self.headspace = headspace
        self.exponent = 0.5  # Stevens' Law exponent for odorants
        self.synergy_pow = 1.3  # Group synergy factor (1.0 = additive, 2.0 = vector-sum)
        
        # Expert synergy pairs (Hofmann & Schieberle 2000)
        # (Compound A, Compound B) -> multiplier
        self.synergy_pairs = {
            ("2-furfurylthiol", "methional"): 1.8, # Meaty synergy
            ("2-methyl-3-furanthiol", "methional"): 2.0, # Stronger meaty synergy
            ("4-hydroxy-2,5-dimethyl-3(2h)-furanone", "furaneol"): 1.2 # caramel-sweet synergy
        }

    def predict_profile(self, 
                        concentration_dict_ppm: Dict[str, float], 
                        temp_c: Optional[float] = None,
                        fat_fraction: float = 0.0,
                        protein_fraction: float = 0.0) -> Dict[str, Tuple[float, float]]:
        """
        Calculate perceived intensity for each compound.
        Returns {name: (intensity, intensity_uncertainty)}
        """
        # 1. Headspace Partitioning (optional)
        if self.headspace and temp_c is not None:
            effective_concs = self.headspace.predict_headspace(
                concentration_dict_ppm, temp_c, fat_fraction, protein_fraction
            )
        else:
            effective_concs = concentration_dict_ppm

        # 2. Perceived Intensity calculation
        results = {}
        for compound, conc_data in effective_concs.items():
            # Handle both float concs and (conc, unc) tuples
            if isinstance(conc_data, tuple):
                conc, unc_ratio = conc_data # unc_ratio is roughly fractional uncertainty on conc
            else:
                conc, unc_ratio = conc_data, 0.2
            
            entry = self.db.find_entry(compound)
            if entry and entry["threshold_ppm"]:
                oav = conc / entry["threshold_ppm"]
                # Apply psychophysical scaling
                intensity = pow(oav, self.exponent) if oav > 0 else 0
                
                # Propagate uncertainty: dI = I * exp * (dc/c)
                # For Maillard, unc is mostly in the barrier (Ea/RT error), leading to large log-uncertainty
                i_unc = intensity * self.exponent * unc_ratio
                results[entry["name"]] = (intensity, i_unc)
        return results

    def get_radar_data(self, 
                       concentration_dict_ppm: Dict[str, float],
                       temp_c: Optional[float] = None,
                       fat_fraction: float = 0.0,
                       protein_fraction: float = 0.0) -> Dict[str, Tuple[float, float]]:
        """
        Groups perceived intensities into radar categories.
        Returns {category: (score, uncertainty)}
        """
        compound_intensities = self.predict_profile(
            concentration_dict_ppm, temp_c, fat_fraction, protein_fraction
        )
        radar = {tag: (0.0, 0.0) for tag in self.db.tags.keys()}
        
        # 1. Group intensities by category
        category_contributions = {tag: [] for tag in self.db.tags.keys()}
        category_uncertainties = {tag: [] for tag in self.db.tags.keys()}
        for tag, search_names in self.db.tags.items():
            for s_name in search_names:
                for c_name, (intensity, unc) in compound_intensities.items():
                    if s_name.lower() in c_name.lower():
                        category_contributions[tag].append(intensity)
                        category_uncertainties[tag].append(unc)
        
        # 2. Apply Synergy Synergy Model
        # Formula: I_group = (Σ I_i^β)^(1/β) where β relates to synergy/suppression
        for tag, intensities in category_contributions.items():
            if not intensities:
                continue
            
            # Filter duplicates if any
            unique_ints = list(set(intensities))
            uncertainties = category_uncertainties[tag]
            
            if tag.lower() in ["meaty", "sulfurous"]:
                # Specialized synergy check for local pairs
                synergy_boost = 1.0
                active_compounds = [c.lower() for c in compound_intensities.keys()]
                for (a, b), boost in self.synergy_pairs.items():
                    if any(a in c for c in active_compounds) and any(b in c for c in active_compounds):
                        synergy_boost = max(synergy_boost, boost)
                
                group_sum = sum([pow(i, 1.1) for i in unique_ints])
                score = pow(group_sum, 1.0/1.1) * synergy_boost
            else:
                # Standard partial additivity
                group_sum = sum([pow(i, self.synergy_pow) for i in unique_ints])
                score = pow(group_sum, 1.0/self.synergy_pow)
            
            # Simple sum of squares for propagated uncertainty
            total_unc = sum([u**2 for u in uncertainties])**0.5
            radar[tag] = (score, total_unc)
                        
        return radar
                        
        return radar

    def get_dominant_notes(self, radar_profile: Dict[str, float], top_n: int = 3) -> List[tuple]:
        """Return the top categories by score."""
        sorted_notes = sorted(radar_profile.items(), key=lambda x: x[1], reverse=True)
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
