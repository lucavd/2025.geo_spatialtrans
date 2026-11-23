#!/usr/bin/env bash
# tools/sonda/S4_run.sh - S4: simulazione dei bracci e catena (docs/sonda/REGOLA_S4.md). Rilanciabile: ogni passo salta cio' che esiste.
# Uso: bash tools/sonda/S4_run.sh <bracci|all> <i,..|all> <nome lista>   (processi della catena dal file $S4_ROOT/WORKERS, al piu' 20)
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S4_ROOT=/mnt/micron/geo_spatialtrans/S4_sonda; export S3_ROOT=$S4_ROOT S3_REF=/mnt/micron/geo_spatialtrans/S3_sonda/ref
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $S4_ROOT; [ -f $S4_ROOT/WORKERS ] || echo 10 > $S4_ROOT/WORKERS
A=$1; I=$2; N=$3; W=$(cat $S4_ROOT/WORKERS)
AA=$([ "$A" = all ] && echo Y,G,U,A,AE || echo $A); II=$([ "$I" = all ] && echo 0,1,2,3,4,5,6,7,8,9 || echo $I)
echo "sim start $(date)"; .venv/bin/python -B tools/sonda/S4_sim_arms.py $W $AA $II 0 || { echo "SIM FAIL"; exit 1; }
python3 tools/sonda/S4_jobs.py $AA $II $N
echo "chain start $(date)"; bash tools/sonda/S3_chain.sh $S4_ROOT/jobs/$N.tsv; echo "chain end $(date)"
