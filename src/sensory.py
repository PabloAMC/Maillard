"""
src/sensory.py

Phase 12: Sensory profile prediction.
Maps predicted volatile concentrations to aroma descriptors 
based on odor detection thresholds (ODT).
"""

from typing import Dict, Any

# Odor Detection Thresholds (ODT) in water/liquid phase (estimate in ppm)
# Data curated from standard flavor chemistry literature (Fenaroli's, etc.)
SENSORY_DB = {
    "FFT": {
        "name": "2-Furfurylthiol",
        "threshold_ppm": 0.00001, # Extremely potent
        "descriptors": ["meaty", "roasted", "sulfury"],
        "smiles": "CC1=CC=C(S1)CS" # Simplified
    },
    "furfural": {
        "name": "Furfural",
        "threshold_ppm": 3.0,
        "descriptors": ["caramel", "sweet", "bready"],
        "smiles": "C1=COC(=C1)C=O"
    },
    "2,3-dimethylpyrazine": {
        "name": "2,3-Dimethylpyrazine",
        "threshold_ppm": 2.5,
        "descriptors": ["roasted", "nutty", "cocoa"],
        "smiles": "CC1=NC=CN=C1C"
    },
    "3-methylbutanal": {
        "name": "3-Methylbutanal",
        "threshold_ppm": 0.2,
        "descriptors": ["malty", "sweaty", "cocoa"],
        "smiles": "CC(C)CC=O"
    },
    "methional": {
        "name": "Methional",
        "threshold_ppm": 0.0002,
        "descriptors": ["potato", "savory", "sulfury"],
        "smiles": "CSCCC=O"
    }
}

class SensoryPredictor:
    """
    Predicts sensory profiles from kinetic simulation data.
    """
    
    def __init__(self, database: Dict[str, Any] = SENSORY_DB):
        self.db = database

    def predict_profile(self, concentration_dict_ppm: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate the Odor Activity Value (OAV) for each descriptor.
        OAV = Concentration / Threshold
        Returns a normalized radar-style profile.
        """
        # Initialize all known descriptors to 0
        all_descriptors = set()
        for entry in self.db.values():
            all_descriptors.update(entry["descriptors"])
        scores = {desc: 0.0 for desc in all_descriptors}
        
        for compound, conc in concentration_dict_ppm.items():
            # Check if we have sensory data for this species
            # We check both the key and the 'name' field
            entry = None
            if compound in self.db:
                entry = self.db[compound]
            else:
                for k, v in self.db.items():
                    if v["name"] == compound:
                        entry = v
                        break
            
            if entry:
                oav = conc / entry["threshold_ppm"]
                # Distribute OAV across descriptors
                for desc in entry["descriptors"]:
                    scores[desc] = scores.get(desc, 0.0) + oav
                    
        return scores

    def get_dominant_notes(self, profile: Dict[str, float], top_n: int = 3) -> list:
        """Return the top descriptors by score."""
        sorted_notes = sorted(profile.items(), key=lambda x: x[1], reverse=True)
        return sorted_notes[:top_n]

if __name__ == "__main__":
    predictor = SensoryPredictor()
    # Mock concentration profile in ppm
    mock_conc = {"FFT": 0.001, "furfural": 10.0, "2,3-dimethylpyrazine": 0.5}
    profile = predictor.predict_profile(mock_conc)
    print("Predicted Sensory Profile (OAV scores):")
    for note, score in predictor.get_dominant_notes(profile):
        print(f"  - {note}: {score:.2f}")
