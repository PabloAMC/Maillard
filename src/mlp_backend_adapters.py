from __future__ import annotations

import importlib
import os
from dataclasses import dataclass
from typing import Optional, Protocol

from src.mlp_adoption_contract import MLPModelCandidate


@dataclass(frozen=True)
class BackendAvailability:
    available: bool
    backend_available: bool
    reason: str


class MLPBackendAdapter(Protocol):
    candidate_id: str
    model_family: str
    model_name: str
    proposed_role: str

    def probe_availability(self) -> BackendAvailability:
        ...

    def optimize_geometry(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        ...

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        ...


class _BaseAdapter:
    def __init__(self, candidate: MLPModelCandidate):
        self.candidate_id = candidate.candidate_id
        self.model_family = candidate.model_family
        self.model_name = candidate.model_name
        self.backend_locator = candidate.backend_locator
        self.proposed_role = candidate.proposed_role

    def optimize_geometry(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        raise ImportError(f"{self.model_family} geometry backend is unavailable: {self.probe_availability().reason}")

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        raise ImportError(f"{self.model_family} TS-seed backend is unavailable: {self.probe_availability().reason}")

    def _resolve_backend_locator(self) -> Optional[str]:
        if not self.backend_locator:
            return None
        if self.backend_locator.startswith("env:"):
            env_var = self.backend_locator.split(":", 1)[1].strip()
            resolved = os.environ.get(env_var)
            if not resolved:
                raise ImportError(f"Backend locator environment variable is not set: {env_var}")
            return resolved
        return self.backend_locator

    def _load_python_backend_module(self):
        resolved = self._resolve_backend_locator()
        if not resolved or not resolved.startswith("python:"):
            return None
        module_name = resolved.split(":", 1)[1].strip()
        if not module_name:
            raise ImportError(f"Invalid python backend locator for {self.model_family}: {resolved}")
        return importlib.import_module(module_name)


class MACEBackendAdapter(_BaseAdapter):
    def __init__(self, candidate: MLPModelCandidate, *, model_family: str):
        super().__init__(candidate)
        self.model_family = model_family

    def _build_optimizer(self):
        from src.mlp_optimizer import MLPOptimizer

        return MLPOptimizer(
            model_name=self.model_name,
            model_family=self.model_family,
            backend_locator=self.backend_locator,
            device="cpu",
        )

    def probe_availability(self) -> BackendAvailability:
        try:
            self._build_optimizer()
        except ImportError as exc:
            return BackendAvailability(False, False, str(exc))
        except Exception as exc:
            return BackendAvailability(False, True, f"{self.model_family} backend import succeeded but initialization failed: {exc}")
        return BackendAvailability(True, True, "evaluated")

    def optimize_geometry(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        optimizer = self._build_optimizer()
        return optimizer.optimize_geometry(xyz_string, fmax=fmax, max_steps=max_steps)

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        optimizer = self._build_optimizer()
        return optimizer.optimize_ts(xyz_string, fmax=fmax, max_steps=max_steps)


class AIMNet2BackendAdapter(_BaseAdapter):
    def probe_availability(self) -> BackendAvailability:
        try:
            plugin_module = self._load_python_backend_module()
        except ImportError as exc:
            return BackendAvailability(False, False, str(exc))
        if plugin_module is not None:
            if not hasattr(plugin_module, "probe_backend") or not hasattr(plugin_module, "prepare_ts_seed"):
                return BackendAvailability(False, False, "Configured AIMNet2 python backend does not expose probe_backend and prepare_ts_seed")
            available, reason = plugin_module.probe_backend(self.model_name, self.backend_locator)
            return BackendAvailability(bool(available), bool(available), str(reason))
        try:
            __import__("aimnet2calc")
        except ImportError:
            try:
                __import__("aimnet2")
            except ImportError as exc:
                return BackendAvailability(False, False, f"AIMNet2 backend not installed: {exc}")
        return BackendAvailability(False, True, "AIMNet2 adapter contract exists, but repo-local ASE/TS hook is not wired yet")

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        plugin_module = self._load_python_backend_module()
        if plugin_module is None or not hasattr(plugin_module, "prepare_ts_seed"):
            return super().prepare_ts_seed(xyz_string, fmax=fmax, max_steps=max_steps)
        return plugin_module.prepare_ts_seed(xyz_string, self.model_name, self.backend_locator, fmax, max_steps)


class OrbMolBackendAdapter(_BaseAdapter):
    def probe_availability(self) -> BackendAvailability:
        try:
            plugin_module = self._load_python_backend_module()
        except ImportError as exc:
            return BackendAvailability(False, False, str(exc))
        if plugin_module is not None:
            if not hasattr(plugin_module, "probe_backend") or not hasattr(plugin_module, "prepare_ts_seed"):
                return BackendAvailability(False, False, "Configured OrbMol python backend does not expose probe_backend and prepare_ts_seed")
            available, reason = plugin_module.probe_backend(self.model_name, self.backend_locator)
            return BackendAvailability(bool(available), bool(available), str(reason))
        for module_name in ("orb_models", "orb_models.forcefield", "orb_models.ase"):
            try:
                __import__(module_name)
                return BackendAvailability(False, True, "OrbMol adapter contract exists, but repo-local ASE/TS hook is not wired yet")
            except ImportError:
                continue
        return BackendAvailability(False, False, "OrbMol backend not installed in the active environment")

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        plugin_module = self._load_python_backend_module()
        if plugin_module is None or not hasattr(plugin_module, "prepare_ts_seed"):
            return super().prepare_ts_seed(xyz_string, fmax=fmax, max_steps=max_steps)
        return plugin_module.prepare_ts_seed(xyz_string, self.model_name, self.backend_locator, fmax, max_steps)


class ReactOTBackendAdapter(_BaseAdapter):
    """Adapter for React-OT (deepprinciple/react-ot).

    React-OT is installed in a *secondary* conda env (see
    ``scripts/setup_react_ot_env.sh``) and reached via ``subprocess`` from
    the plugin module ``src.react_ot_backend``. It requires both reactant
    and product geometries, so it cannot satisfy the single-seed
    ``prepare_ts_seed`` contract used by the TS-seed-benchmark validator;
    ``probe_availability`` therefore reports ``available=False`` with an
    explicit contract-gap reason whenever the adapter is consulted from
    that lane. The pilot wrapper
    ``scripts/recover_ts_react_ot_seed.py`` calls
    ``predict_ts_from_reactant_product`` from the plugin module directly.
    """

    def probe_availability(self) -> BackendAvailability:
        try:
            plugin_module = self._load_python_backend_module()
        except ImportError as exc:
            return BackendAvailability(False, False, str(exc))
        if plugin_module is None:
            return BackendAvailability(
                False,
                False,
                "react_ot adapter requires backend_locator='python:src.react_ot_backend'",
            )
        if not hasattr(plugin_module, "probe_backend"):
            return BackendAvailability(
                False,
                False,
                "Configured React-OT python backend does not expose probe_backend",
            )
        try:
            ok, reason = plugin_module.probe_backend(self.model_name, self.backend_locator)
        except Exception as exc:  # noqa: BLE001
            return BackendAvailability(False, False, f"react_ot probe_backend raised: {exc}")
        if not ok:
            return BackendAvailability(False, False, str(reason))
        return BackendAvailability(
            False,
            True,
            "react_ot backend present, but TS-seed benchmark contract requires reactant+product; "
            "use scripts/recover_ts_react_ot_seed.py",
        )

    def prepare_ts_seed(self, xyz_string: str, *, fmax: float = 0.05, max_steps: int = 100) -> str:
        raise NotImplementedError(
            "React-OT requires reactant and product XYZ inputs. Use the pilot wrapper "
            "scripts/recover_ts_react_ot_seed.py instead of the TS-seed benchmark contract."
        )


class UnsupportedBackendAdapter(_BaseAdapter):
    def __init__(self, candidate: MLPModelCandidate, reason: str):
        super().__init__(candidate)
        self._reason = reason

    def probe_availability(self) -> BackendAvailability:
        return BackendAvailability(False, False, self._reason)


def build_candidate_adapter(candidate: MLPModelCandidate) -> MLPBackendAdapter:
    if candidate.model_family in {"mace_mp", "mace_off"}:
        return MACEBackendAdapter(candidate, model_family=candidate.model_family)
    if candidate.model_family == "mace_omol":
        return MACEBackendAdapter(candidate, model_family="mace_omol")
    if candidate.model_family == "aimnet2":
        return AIMNet2BackendAdapter(candidate)
    if candidate.model_family == "orbmol":
        return OrbMolBackendAdapter(candidate)
    if candidate.model_family == "react_ot":
        return ReactOTBackendAdapter(candidate)
    return UnsupportedBackendAdapter(candidate, f"No explicit backend adapter is registered for {candidate.model_family}")


__all__ = [
    "BackendAvailability",
    "MLPBackendAdapter",
    "build_candidate_adapter",
]