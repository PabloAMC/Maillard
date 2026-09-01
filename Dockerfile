# syntax=docker/dockerfile:1
#
# Maillard container. 2026-09-01 (cleaning branch, Phase 1a): the QM/DFT lane and its
# build-time patches (CPU torch, jax, xtbiff, e3nn/MACE torch.load shims, pyGSM AST
# rewrites) were removed together with the lane. The image is now a plain conda env.

FROM condaforge/miniforge3:latest

SHELL ["/bin/bash", "-lc"]

WORKDIR /workspace

COPY environment.yml ./environment.yml

RUN conda env create -n maillard --file environment.yml && \
    conda clean -afy && \
    rm -f environment.yml

ENV PATH=/opt/conda/envs/maillard/bin:$PATH
ENV CONDA_DEFAULT_ENV=maillard
ENV MKL_INTERFACE_LAYER=LP64

COPY . .

CMD ["bash"]
