import re

def normalize_compound_name(name: str) -> str:
    """
    Strips all non-alphanumeric characters and converts to lowercase.
    Used for strict matching of benchmark and target compound identifiers.
    """
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def normalize_name_underscored(name: str) -> str:
    """Lowercase, collapse whitespace/dashes into single underscores."""
    lowered = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    return "_".join(part for part in lowered.split("_") if part)


def normalize_name_spaced(name: str) -> str:
    """Lowercase, collapse underscores/dashes into single spaces."""
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())
