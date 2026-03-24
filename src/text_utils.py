import re

def normalize_compound_name(name: str) -> str:
    """
    Strips all non-alphanumeric characters and converts to lowercase.
    Used for strict matching of benchmark and target compound identifiers.
    """
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())
