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

# ----------------------------------------------------------------------------
# pyGSM (Zimmerman group) — used as TS-guess generator (DE-GSM at xTB level).
# Pinned commit + post-install patches for Python 3.12 / NumPy 2.x / setuptools 82
# compatibility. The xTB barriers from pyGSM are diagnostic only; final barriers
# always come from DFT. Keep these patches in sync with src/gsm_backend.py.
# ----------------------------------------------------------------------------
ARG PYGSM_REF=295dc46db57373828e310373496ffac2e815dbf3
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate maillard && \
    pip install --no-cache-dir --no-deps \
        "git+https://github.com/ZimmermanGroup/pyGSM.git@${PYGSM_REF}" && \
    python - <<'PY'
"""Reproducible compatibility patches for the pinned pyGSM commit.

1. setuptools >= 70 dropped pkg_resources.parse_version → use packaging.version.parse.
2. NumPy 2.0 removed ndarray.tostring() → use ndarray.tobytes() (identical bytes).
3. NumPy 2.0 forbids np.sum(<generator>) → wrap as np.sum([...], axis=0).
4. NumPy 2.0 forbids implicit ndarray→scalar coercion in `"%f" % arr`.
   Rewrite every ``BinOp(str_const, %, expr)`` in pyGSM so the RHS is wrapped
   with ``_pygsm_scalarize`` (defined in pyGSM/__init__.py), which converts
   length-1 ndarrays to Python scalars before the % operator runs.

All four are pure-syntax shims with no semantic change. They are applied
in-place at image-build time so every container starts from the same state.
"""
import ast
import re
import sys
from pathlib import Path

site = Path(sys.prefix) / "lib/python3.12/site-packages/pyGSM"
if not site.exists():
    raise SystemExit(f"pyGSM not found at {site}")

# (1) pkg_resources → packaging.version
topology = site / "coordinate_systems/topology.py"
src = topology.read_text()
old = "from pkg_resources import parse_version"
new = "from packaging.version import parse as parse_version"
if old in src:
    topology.write_text(src.replace(old, new))

# (2) ndarray.tostring → ndarray.tobytes
for path in site.rglob("*.py"):
    text = path.read_text()
    if ".tostring(" in text:
        path.write_text(text.replace(".tostring(", ".tobytes("))

# (3) np.sum(<generator>) → np.sum([<generator>], axis=0)
def _wrap_np_sum_generators(source: str) -> str:
    out = []
    i = 0
    while True:
        idx = source.find("np.sum(", i)
        if idx < 0:
            out.append(source[i:])
            return "".join(out)
        out.append(source[i:idx])
        # Find matching closing paren of np.sum(
        depth = 0
        j = idx + len("np.sum(")
        start_arg = j
        while j < len(source):
            c = source[j]
            if c == "(":
                depth += 1
            elif c == ")":
                if depth == 0:
                    break
                depth -= 1
            j += 1
        if j >= len(source):
            # unbalanced; bail out
            out.append(source[idx:])
            return "".join(out)
        arg = source[start_arg:j]
        stripped = arg.strip()
        # Only rewrite if it's a bare generator expression (contains " for "
        # at depth 0, NOT already wrapped in [] or [ ... ]).
        is_generator = (
            " for " in arg
            and not stripped.startswith("[")
            and not stripped.startswith("(")
        )
        if is_generator:
            out.append(f"np.sum([{arg}], axis=0)")
        else:
            out.append(source[idx : j + 1])
        i = j + 1

for path in site.rglob("*.py"):
    text = path.read_text()
    if "np.sum(" not in text:
        continue
    new_text = _wrap_np_sum_generators(text)
    if new_text != text:
        path.write_text(new_text)

# (4a) Inject the runtime scalarize helper into pyGSM/__init__.py.
SCALARIZE_HELPER = '''
# >>> NumPy-2.x compatibility shim injected by Maillard Dockerfile <<<
def _pygsm_scalarize(_x):
    """Coerce length-1 ndarrays / 0-d arrays to Python scalars.

    NumPy 2.0 forbids implicit array→scalar conversion when used as a
    ``%`` format operand.  pyGSM has 199 such call sites; rather than
    patching each, every ``str %`` expression in the package is rewritten
    at build time to ``str % _pygsm_scalarize(rhs)``.

    Tuples are recursed element-wise so multi-arg formats keep working.
    """
    try:
        import numpy as _np
    except Exception:
        _np = None
    if isinstance(_x, tuple):
        return tuple(_pygsm_scalarize(_e) for _e in _x)
    if _np is not None and isinstance(_x, _np.ndarray) and _x.size == 1:
        return _x.item()
    return _x
# <<< end NumPy-2.x compatibility shim >>>
'''
init_py = site / "__init__.py"
init_text = init_py.read_text() if init_py.exists() else ""
if "_pygsm_scalarize" not in init_text:
    init_py.write_text(init_text + SCALARIZE_HELPER)

# (4b) AST rewrite: wrap RHS of every ``str_const % X`` BinOp with
# ``_pygsm_scalarize(X)``, and the sole arg of ``float(X)`` / ``int(X)``
# calls.  Both rewrites are needed because NumPy 2.x forbids implicit
# ndarray→scalar coercion in either path.  We use ast.unparse to round-trip.
class _ScalarizeRewriter(ast.NodeTransformer):
    def __init__(self) -> None:
        self.changed = False

    def _wrap(self, node: ast.AST) -> ast.Call:
        new_call = ast.Call(
            func=ast.Name(id="_pygsm_scalarize", ctx=ast.Load()),
            args=[node],
            keywords=[],
        )
        ast.copy_location(new_call, node)
        return new_call

    @staticmethod
    def _is_already_wrapped(node: ast.AST) -> bool:
        return (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "_pygsm_scalarize")

    def visit_BinOp(self, node: ast.BinOp) -> ast.AST:
        self.generic_visit(node)
        if (isinstance(node.op, ast.Mod)
                and isinstance(node.left, ast.Constant)
                and isinstance(node.left.value, str)
                and not self._is_already_wrapped(node.right)):
            node.right = self._wrap(node.right)
            self.changed = True
        return node

    def visit_Call(self, node: ast.Call) -> ast.AST:
        self.generic_visit(node)
        if (isinstance(node.func, ast.Name)
                and node.func.id in {"float", "int"}
                and len(node.args) == 1
                and not node.keywords
                and not self._is_already_wrapped(node.args[0])):
            node.args[0] = self._wrap(node.args[0])
            self.changed = True
        return node


_IMPORT_LINE = "from pyGSM import _pygsm_scalarize  # NumPy-2 compat shim"

for path in site.rglob("*.py"):
    if path.name == "__init__.py" and path.parent == site:
        continue  # don't rewrite the file that defines the helper
    text = path.read_text()
    # Skip files with neither `%` formatting nor float()/int() conversions.
    if ("%" not in text) and ("float(" not in text) and ("int(" not in text):
        continue
    try:
        tree = ast.parse(text)
    except SyntaxError:
        continue
    rewriter = _ScalarizeRewriter()
    new_tree = rewriter.visit(tree)
    if not rewriter.changed:
        continue
    # Insert the import after any module docstring and __future__ imports.
    insert_idx = 0
    body = new_tree.body
    if (body and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)):
        insert_idx = 1
    while (insert_idx < len(body)
           and isinstance(body[insert_idx], ast.ImportFrom)
           and body[insert_idx].module == "__future__"):
        insert_idx += 1
    import_node = ast.ImportFrom(
        module="pyGSM",
        names=[ast.alias(name="_pygsm_scalarize", asname=None)],
        level=0,
    )
    body.insert(insert_idx, import_node)
    ast.fix_missing_locations(new_tree)
    new_src = ast.unparse(new_tree)
    path.write_text(new_src + "\n")

# Smoke-import to fail the build early if any patch is incomplete.
import pyGSM  # noqa: F401
from pyGSM import _pygsm_scalarize
from pyGSM.molecule import Molecule  # noqa: F401
from pyGSM.growing_string_methods import DE_GSM  # noqa: F401
from pyGSM.utilities import math_utils
import numpy as np
math_utils.orthogonalize(np.eye(3))
# Validate the scalarize shim end-to-end.
assert "%1.4f" % _pygsm_scalarize(np.array([1.2345])) == "1.2345"
assert ("%d %1.2f" % _pygsm_scalarize((np.int64(3), np.array([0.5]))) == "3 0.50")
assert float(_pygsm_scalarize(np.array([3.14]))) == 3.14
assert int(_pygsm_scalarize(np.array([7]))) == 7
print("pyGSM patches applied OK")
PY

ENV PATH=/opt/conda/envs/maillard/bin:$PATH
ENV CONDA_DEFAULT_ENV=maillard
ENV MKL_INTERFACE_LAYER=LP64

COPY . .

CMD ["bash"]
