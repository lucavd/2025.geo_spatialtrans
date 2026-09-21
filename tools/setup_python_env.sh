#!/usr/bin/env bash
# tools/setup_python_env.sh — ambienti Python di progetto (uv), idempotente.
# R1: .venv  -> tifffile, imagecodecs, pyarrow, h5py, scikit-image (requirements.lock)
# R2: .venv  += torch + Cellpose (CUDA, requirements-r2-cellpose.txt -> requirements-cellpose.lock)
#     .venv-stardist -> tensorflow[and-cuda] + stardist (requirements-r2-stardist.txt -> requirements-stardist.lock)
# Due venv perche' torch e tensorflow pinnano versioni diverse delle librerie nvidia-* (cuDNN).
set -euo pipefail
cd "$(dirname "$0")/.."
UV=${UV:-$HOME/.local/bin/uv}

install_env () {  # $1 venv dir, $2 requirements, $3 lockfile
  [ -d "$1" ] || "$UV" venv "$1" --python 3.12
  if [ -f "$3" ]; then
    "$UV" pip install --python "$1/bin/python" -r "$3"
  else
    "$UV" pip install --python "$1/bin/python" -r "$2"
    "$UV" pip freeze --python "$1/bin/python" > "$3"
  fi
}

install_env .venv tools/requirements-r1.txt tools/requirements.lock
.venv/bin/python -c "import tifffile, imagecodecs, pyarrow, h5py, skimage, numpy, pandas, matplotlib; print('python env OK', tifffile.__version__, pyarrow.__version__)"

if [ "${R2:-1}" = "1" ]; then
  install_env .venv tools/requirements-r2-cellpose.txt tools/requirements-cellpose.lock
  .venv/bin/python -c "import torch, cellpose; print('cellpose env OK', cellpose.version, 'torch', torch.__version__, 'cuda', torch.cuda.is_available())"
  install_env .venv-stardist tools/requirements-r2-stardist.txt tools/requirements-stardist.lock
  tools/py_stardist.sh -c "import tensorflow as tf, stardist; print('stardist env OK', stardist.__version__, 'tf', tf.__version__, 'gpus', len(tf.config.list_physical_devices('GPU')))"
fi
