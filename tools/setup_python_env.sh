#!/usr/bin/env bash
# tools/setup_python_env.sh — ambiente Python di progetto (uv), idempotente.
# R1: lettura immagini BigTIFF (.btf), parquet, h5, figure. R2 aggiungerà Cellpose/StarDist (CUDA).
set -euo pipefail
cd "$(dirname "$0")/.."
UV=${UV:-$HOME/.local/bin/uv}
[ -d .venv ] || "$UV" venv .venv --python 3.12
if [ -f tools/requirements.lock ]; then
  "$UV" pip install --python .venv/bin/python -r tools/requirements.lock
else
  "$UV" pip install --python .venv/bin/python -r tools/requirements-r1.txt
  "$UV" pip freeze --python .venv/bin/python > tools/requirements.lock
fi
.venv/bin/python -c "import tifffile, imagecodecs, pyarrow, h5py, skimage, numpy, pandas, matplotlib; print('python env OK', tifffile.__version__, pyarrow.__version__)"
