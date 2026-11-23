#!/usr/bin/env bash
# tools/sonda/S3_run.sh - S3: corsa completa (docs/sonda/REGOLA_S3.md). Liste (tools/sonda/S3_jobs.py); catena su V, T, P (in quest'ordine);
# taratura di B su V (tools/sonda/S3_Bparams.py); simulazioni B, 10 aree x 2 semi (tools/sonda/S3_sim_B.py, 10 processi); catena su B;
# tabella, h, D. Processi della catena dal file $S3_ROOT/WORKERS (al piu' 20). Rilanciabile: ogni passo salta cio' che esiste.
# Uso: bash tools/sonda/S3_run.sh
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S3_ROOT=/mnt/micron/geo_spatialtrans/S3_sonda OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
J=$S3_ROOT/jobs; python3 tools/sonda/S3_jobs.py
cat $J/V.tsv $J/T.tsv $J/P.tsv > $J/VTP.tsv
echo "chain VTP start $(date)"; bash tools/sonda/S3_chain.sh $J/VTP.tsv; echo "chain VTP end $(date)"
.venv/bin/python -B tools/sonda/S3_Bparams.py || { echo "BPARAMS FAIL"; exit 1; }
.venv/bin/python -B tools/sonda/S3_sim_B.py 10 all 0,1 || { echo "SIM B FAIL"; exit 1; }
echo "chain B start $(date)"; bash tools/sonda/S3_chain.sh $J/B.tsv; echo "chain B end $(date)"
python3 tools/sonda/S3_table.py && .venv/bin/python -B tools/sonda/S3_h.py
echo "run end $(date)"
