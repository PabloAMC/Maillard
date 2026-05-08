# React-OT external model

This directory holds provenance and (optionally) the local checkpoint for the
React-OT diffusion-based transition-state generator
([deepprinciple/react-ot](https://github.com/deepprinciple/react-ot)).

React-OT is **not** wired into the runtime barrier path. It is treated strictly
as an offline geometric seed generator for the eligible CHON targets listed in
`provenance.json`. Energies produced by the model are never propagated as
runtime barrier authority; every promising seed must clear downstream Sella
DFT and imaginary-mode validation before any barrier change is considered.

## Files

- `provenance.json` — pinned upstream ref, secondary env metadata,
  applicability scope, trust posture.
- `sb-pretrained.ckpt` — *not committed*. Download into this directory after
  running `scripts/setup_react_ot_env.sh`.

## Workflow

1. Build the secondary conda env in the container:
   ```
   ./scripts/docker_maillard.sh react-ot-setup
   ```
2. Download the upstream checkpoint into this folder (the setup script prints
   the exact `wget` command).
3. Run the smoke gate on a small CHON pair:
   ```
   ./scripts/docker_maillard.sh react-ot-smoke
   ```
4. If the smoke gate succeeds, run the pilot wrapper on the eligible target
   set:
   ```
   ./scripts/docker_maillard.sh react-ot-pilot
   ```
5. If the CPU smoke gate fails for `torch_geometric` / CUDA reasons, fall back
   to Google Colab:
   ```
   python scripts/prepare_react_ot_colab_bundle.py
   ```
   Upload the emitted bundle to `notebooks/react_ot_colab_gpu.ipynb`, run the
   notebook on a GPU runtime, then import the downloaded archive locally:
   ```
   ./scripts/docker_maillard.sh react-ot-import-colab PATH/TO/react_ot_colab_artifacts.zip
   ```

## Open the Notebook

- IDE notebook file: `notebooks/react_ot_colab_gpu.ipynb`
- GitHub view: https://github.com/PabloAMC/Maillard/blob/qm-barriers/notebooks/react_ot_colab_gpu.ipynb
- Google Colab: https://colab.research.google.com/github/PabloAMC/Maillard/blob/qm-barriers/notebooks/react_ot_colab_gpu.ipynb
