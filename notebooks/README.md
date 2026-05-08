# React-OT Colab Notebook

- Local notebook in this workspace: [react_ot_colab_gpu.ipynb](react_ot_colab_gpu.ipynb)
- GitHub view: https://github.com/PabloAMC/Maillard/blob/qm-barriers/notebooks/react_ot_colab_gpu.ipynb
- Open directly in Google Colab: https://colab.research.google.com/github/PabloAMC/Maillard/blob/qm-barriers/notebooks/react_ot_colab_gpu.ipynb

Recommended local flow:

1. Generate the upload bundle:
   ```bash
   python scripts/prepare_react_ot_colab_bundle.py
   ```
2. Open the notebook in Colab from the link above.
3. Download `react_ot_colab_artifacts.zip` from Colab.
4. Import the artifacts back into the repo:
   ```bash
   python scripts/import_react_ot_colab_artifacts.py PATH/TO/react_ot_colab_artifacts.zip
   ```