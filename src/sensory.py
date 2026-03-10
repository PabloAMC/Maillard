"""
src/sensory.py

Phase C: Full Sensory Prediction Model.
Maps predicted volatile concentrations to aroma descriptors using a psychophysical mixing model.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional

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
    
    def __init__(self, database: Optional[SensoryDatabase] = None):
        self.db = database or SensoryDatabase()
        self.exponent = 0.5  # Stevens' Law exponent for odorants

    def predict_profile(self, concentration_dict_ppm: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate perceived intensity for each compound using Stevens' Law.
        I = (C / Threshold)^n
        Returns raw OAV_eff scores per compound name.
        """
        results = {}
        for compound, conc in concentration_dict_ppm.items():
            entry = self.db.find_entry(compound)
            if entry and entry["threshold_ppm"]:
                oav = conc / entry["threshold_ppm"]
                # Apply psychophysical scaling
                intensity = pow(oav, self.exponent) if oav > 0 else 0
                results[entry["name"]] = intensity
        return results

    def get_radar_data(self, concentration_dict_ppm: Dict[str, float]) -> Dict[str, float]:
        """
        Groups perceived intensities into radar categories (meaty, roasted, etc.)
        Sums intensities within each category.
        """
        compound_intensities = self.predict_profile(concentration_dict_ppm)
        radar = {tag: 0.0 for tag in self.db.tags.keys()}
        
        # Map compounds to categories using tags
        for tag, search_names in self.db.tags.items():
            for s_name in search_names:
                # Find if any compound in results matches the tag search name
                for c_name, intensity in compound_intensities.items():
                    if s_name.lower() in c_name.lower():
                        radar[tag] += intensity
                        
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
    print("Radar Category Scores:")
    for note, score in predictor.get_dominant_notes(radar):
        print(f"  - {note}: {score:.4f}")
