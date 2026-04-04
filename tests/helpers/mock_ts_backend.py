from __future__ import annotations


def probe_backend(model_name: str, backend_locator: str | None):
    return True, f"mock backend ready for {model_name}"


def prepare_ts_seed(
    xyz_string: str,
    model_name: str,
    backend_locator: str | None,
    fmax: float,
    max_steps: int,
) -> str:
    return xyz_string