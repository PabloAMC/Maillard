from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from src.input_normalization import normalize_aliases


_FORMULATION_INPUT_ALIASES = {
    # "pH" is the natural spelling a scientist types; keep it first so the
    # obvious key works without the caller having to know the internal casing.
    "ph": ("pH", "PH"),
    "temperature": ("temp", "TEMP"),
    "water_activity": ("aw", "AW"),
    "protein_type": ("protein",),
    "matrix_type": ("matrix",),
}

# The keys `src.pipeline.MaillardPipeline.evaluate_all` actually reads off a
# formulation mapping. `to_dict()` MUST emit these, otherwise a serialized
# Formulation silently loses its process conditions when fed back to the
# pipeline (it fell back to the global ReactionConditions instead).
CANONICAL_PIPELINE_CONDITION_KEYS = ("ph", "temp", "aw", "time_minutes")

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
            "TEMP": "temperature",
            "temperature": "temperature",
            "aw": "water_activity",
            "AW": "water_activity",
            "water_activity": "water_activity",
            "ph": "ph",
            "pH": "ph",
            "PH": "ph",
            "protein": "protein_type",
            "matrix": "matrix_type",
        }
        return getattr(self, alias_map.get(key, key), default)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Formulation:
        """
        Creates a Formulation instance from a dictionary, handling legacy
        key mappings (e.g. 'temp' -> 'temperature', 'aw' -> 'water_activity').
        """
        normalized = normalize_aliases(data, _FORMULATION_INPUT_ALIASES)

        # Extract name (required)
        name = normalized.get("name", "unnamed_formulation")
        
        # Extract components with fallbacks
        sugars = list(normalized.get("sugars", []))
        amino_acids = list(normalized.get("amino_acids", []))
        lipids = list(normalized.get("lipids", []))
        additives = list(normalized.get("additives", []))
        interventions = list(normalized.get("interventions", []))
        molar_ratios = dict(normalized.get("molar_ratios", {}))
        
        # Handle legacy condition keys — use None when key not present so pipeline falls back to conditions.pH
        _ph_raw = normalized.get("ph")
        ph = float(_ph_raw) if _ph_raw is not None else None
        temperature = float(normalized.get("temperature", 120.0))
        water_activity = float(normalized.get("water_activity", 0.8))
        time_minutes = float(normalized.get("time_minutes", 60.0))
        
        # Structural fields
        catalyst = normalized.get("catalyst")
        thiamine_availability = normalized.get("thiamine_availability")
        denaturation_state = normalized.get("denaturation_state")
        protein_type = normalized.get("protein_type", "free")
        matrix_type = normalized.get("matrix_type")
        notes = normalized.get("notes", "")
        
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
        """Returns a dictionary representation for serialization.

        Schema note (2026-08-27): this used to emit only ``temperature`` /
        ``water_activity`` while the pipeline reads ``temp`` / ``aw``, so a
        round-tripped formulation quietly lost its own temperature and water
        activity and inherited the global ``ReactionConditions`` instead. The
        canonical pipeline-read keys are now emitted alongside the dataclass
        field names; ``from_dict`` accepts either spelling and prefers the
        dataclass names, so the round trip is lossless in both directions.
        """
        return {
            "name": self.name,
            "sugars": self.sugars,
            "amino_acids": self.amino_acids,
            "lipids": self.lipids,
            "additives": self.additives,
            "interventions": self.interventions,
            "molar_ratios": self.molar_ratios,
            "ph": self.ph,
            # canonical pipeline-read keys
            "temp": self.temperature,
            "aw": self.water_activity,
            # dataclass-native names (kept so from_dict round-trips exactly)
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
