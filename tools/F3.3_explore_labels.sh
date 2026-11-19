#!/usr/bin/env bash
# tools/F3.3_explore_labels.sh - F3.3 E-2: etichette trasferite su 30 aree vere di taratura e 30 simulate (k = 0), 4 processi
cd ~/2025.geo_spatialtrans
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans; O=$M/F3.3/explore/E2; mkdir -p $O
{ for d in $M/F3.1/main/real/*; do t=$(basename $d); echo "$d ${t:0:2} real $O/real_$t.json"; done
  for d in $M/F3.1/main/sim/*_k0; do t=$(basename $d); echo "$d ${t:0:2} sim $O/sim_$t.json"; done; } > $O/jobs.txt
xargs -P 4 -L 1 bash -c '[ -s "$4" ] || .venv-b2c/bin/python tools/F3.3_explore_labels.py "$1" "$2" "$3" "$4" > "$4.log" 2>&1' _ < $O/jobs.txt
ls $O/*.json | wc -l
