# syntax=docker/dockerfile:1

FROM condaforge/miniforge3:latest

SHELL ["/bin/bash", "-lc"]

WORKDIR /workspace

COPY environment.yml ./environment.yml

RUN python - <<'PY'
from pathlib import Path

source = Path('environment.yml')
filtered = Path('environment.docker.yml')
lines = source.read_text().splitlines()
filtered.write_text('\n'.join(line for line in lines if 'pytorch::pytorch=' not in line) + '\n')
PY

RUN conda create -n maillard python=3.12 -y && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate maillard && \
    conda install -y -c conda-forge jax jaxlib wget xz jupyter ipywidgets && \
    pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu torch && \
    conda env update -n maillard --file environment.docker.yml && \
    wget -q https://github.com/grimme-lab/xtbiff/releases/download/v1.1/xtbiff.tar.xz && \
    tar -xf xtbiff.tar.xz && \
    mv xtbiff "$CONDA_PREFIX/bin/xtbiff" && \
    chmod +x "$CONDA_PREFIX/bin/xtbiff" && \
    rm -f xtbiff.tar.xz environment.docker.yml

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate maillard && \
    python - <<'PY'
from pathlib import Path
import os

patches = [
    (
        Path(os.environ['CONDA_PREFIX']) / 'lib/python3.12/site-packages/e3nn/o3/_wigner.py',
        "torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'))",
        "torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'), weights_only=False)",
    ),
    (
        Path(os.environ['CONDA_PREFIX']) / 'lib/python3.12/site-packages/mace/calculators/mace.py',
        "torch.load(f=model_path, map_location=device)",
        "torch.load(f=model_path, map_location=device, weights_only=False)",
    ),
]

for path, old, new in patches:
    if not path.exists():
        continue
    text = path.read_text()
    if new in text:
        continue
    if old in text:
        path.write_text(text.replace(old, new))
PY

ENV PATH=/opt/conda/envs/maillard/bin:$PATH
ENV CONDA_DEFAULT_ENV=maillard
ENV MKL_INTERFACE_LAYER=LP64

COPY . .

CMD ["bash"]
