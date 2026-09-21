#!/usr/bin/env bash
# tools/py_stardist.sh — esegue python del venv StarDist con le librerie CUDA del venv nel loader path.
# tensorflow[and-cuda] 2.21 installa nvidia-*-cu12 (12.9) ma non li trova a runtime ("Cannot dlopen some GPU
# libraries", 0 GPU): con LD_LIBRARY_PATH sulle cartelle nvidia/*/lib la GPU compare (verificato R2, 2026-09-20).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SP="$ROOT/.venv-stardist/lib/python3.12/site-packages"
export LD_LIBRARY_PATH="$(ls -d "$SP"/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export TF_CPP_MIN_LOG_LEVEL=${TF_CPP_MIN_LOG_LEVEL:-2}
exec "$ROOT/.venv-stardist/bin/python" "$@"
