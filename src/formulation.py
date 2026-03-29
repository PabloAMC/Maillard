from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

@dataclass
class Formulation:
    """
    Typed schema representing a Maillard precursor formulation and its
    associated process conditions.
    """
    name: str
    sugars: List[str] = field(default_factory=list)
    amino_acids: List[str] = field(default_factory=list)
    lipids: List[str] = field(default_factory=list)
    additives: List[str] = field(default_factory=list)
    interventions: List[Any] = field(default_factory=list)
    molar_ratios: Dict[str, float] = field(default_factory=dict)
    
    # Process conditions (often embedded in formulation grid)
    ph: Optional[float] = None  # if None, pipeline falls back to ReactionConditions.pH
    temperature: float = 120.0  # Celsius
    water_activity: float = 0.8
    time_minutes: float = 60.0
    
    # Structural metadata
    catalyst: Optional[str] = None
    thiamine_availability: Optional[bool] = None
    denaturation_state: Optional[float] = None
    protein_type: str = "free"  # e.g. "soy_isolate", "pea_isolate", "free"
    matrix_type: Optional[str] = None
    notes: str = ""

    def get(self, key: str, default: Any = None) -> Any:
        alias_map = {
            "temp": "temperature",
            "temperature": "temperature",
            "aw": "water_activity",
            "water_activity": "water_activity",
            "ph": "ph",
        }
        return getattr(self, alias_map.get(key, key), default)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Formulation:
        """
        Creates a Formulation instance from a dictionary, handling legacy
        key mappings (e.g. 'temp' -> 'temperature', 'aw' -> 'water_activity').
        """
        # Extract name (required)
        name = data.get("name", "unnamed_formulation")
        
        # Extract components with fallbacks
        sugars = list(data.get("sugars", []))
        amino_acids = list(data.get("amino_acids", []))
        lipids = list(data.get("lipids", []))
        additives = list(data.get("additives", []))
        interventions = list(data.get("interventions", []))
        molar_ratios = dict(data.get("molar_ratios", {}))
        
        # Handle legacy condition keys — use None when key not present so pipeline falls back to conditions.pH
        _ph_raw = data.get("ph", data.get("PH"))
        ph = float(_ph_raw) if _ph_raw is not None else None
        temperature = float(data.get("temperature", data.get("temp", data.get("TEMP", 120.0))))
        water_activity = float(data.get("water_activity", data.get("aw", data.get("AW", 0.8))))
        time_minutes = float(data.get("time_minutes", 60.0))
        
        # Structural fields
        catalyst = data.get("catalyst")
        thiamine_availability = data.get("thiamine_availability")
        denaturation_state = data.get("denaturation_state")
        protein_type = data.get("protein_type", data.get("protein", "free"))
        matrix_type = data.get("matrix_type", data.get("matrix", None))
        notes = data.get("notes", "")
        
        return cls(
            name=name,
            sugars=sugars,
            amino_acids=amino_acids,
            lipids=lipids,
            additives=additives,
            interventions=interventions,
            molar_ratios=molar_ratios,
            ph=ph,
            temperature=temperature,
            water_activity=water_activity,
            time_minutes=time_minutes,
            catalyst=catalyst,
            thiamine_availability=thiamine_availability,
            denaturation_state=denaturation_state,
            protein_type=protein_type,
            matrix_type=matrix_type,
            notes=notes
        )

    def to_dict(self) -> Dict[str, Any]:
        """Returns a dictionary representation for serialization."""
        return {
            "name": self.name,
            "sugars": self.sugars,
            "amino_acids": self.amino_acids,
            "lipids": self.lipids,
            "additives": self.additives,
            "interventions": self.interventions,
            "molar_ratios": self.molar_ratios,
            "ph": self.ph,
            "temperature": self.temperature,
            "water_activity": self.water_activity,
            "time_minutes": self.time_minutes,
            "catalyst": self.catalyst,
            "thiamine_availability": self.thiamine_availability,
            "denaturation_state": self.denaturation_state,
            "protein_type": self.protein_type,
            "matrix_type": self.matrix_type,
            "notes": self.notes
        }
