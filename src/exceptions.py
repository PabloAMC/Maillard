"""
src/exceptions.py

Domain-specific exceptions for the Maillard framework.
"""

class MaillardError(Exception):
    """Base exception for all domain-specific errors in the Maillard framework."""
    pass

class FormulationInputError(MaillardError):
    """Raised when pipeline inputs violate schema or logical constraints."""
    pass

class BenchmarkError(MaillardError):
    """Raised when benchmark artifacts are missing, unparseable, or invalid."""
    pass

class KineticsError(MaillardError):
    """Raised when numerical ODE integration fails or network builds fail."""
    pass

class UnbalancedReactionError(KineticsError):
    """Raised when an elementary step violates mass/atom conservation."""
    pass

class ConfigurationError(MaillardError):
    """Raised when environment or path configurations are invalid."""
    pass
