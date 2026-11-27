#!/usr/bin/env bash
# tools/sonda/S9_C_run.sh - S9, catena C: tools/sonda/S9_celltypist.py unit <arm> <run> sulle 50 unita' di results/sonda/S3b/area_values.csv
# (V 10, T 10, B 20, P 10). Rilanciabile (salta le unita' con il JSON). Uso: bash tools/sonda/S9_C_run.sh [W] ; W al piu' 4.
set -uo pipefail
cd ~/2025.geo_spatialtrans
W=${1:-4}; [ "$W" -gt 4 ] && W=4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export S9C_OUT=${S9C_OUT:-/mnt/micron/geo_spatialtrans/S9_sonda/C}
mkdir -p $S9C_OUT/logs
one() {
  local arm=$1 run=$2
  .venv-enact/bin/python -B tools/sonda/S9_celltypist.py unit $arm $run > $S9C_OUT/logs/${arm}_${run}.log 2>&1 && echo "done $arm $run" || echo "FAIL $arm $run"
}
export -f one
tail -n +2 results/sonda/S3b/area_values.csv | cut -d, -f1,5 | tr ',' ' ' | xargs -P $W -n 2 bash -c 'one "$@"' _
