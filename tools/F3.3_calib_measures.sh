#!/usr/bin/env bash
# tools/F3.3_calib_measures.sh - F3.3: misure con le due etichettature sulle 30 aree vere e 60 simulate di F3.1 (4 processi)
cd ~/2025.geo_spatialtrans
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans; O=$M/F3.3/calib; mkdir -p $O
{ for d in $M/F3.1/main/real/*; do echo "$d real $O/real_$(basename $d).json"; done
  for d in $M/F3.1/main/sim/*_k[01]; do echo "$d sim $O/sim_$(basename $d).json"; done; } > $O/jobs.txt
xargs -P 4 -L 1 bash -c '[ -s "$3" ] || .venv-b2c/bin/python tools/F3.3_calib_measures.py "$1" "$2" "$3" > "$3.log" 2>&1' _ < $O/jobs.txt
ls $O/*.json | wc -l; grep -h ' dev ' $O/*.log | awk '{print $4}' | sort -g | tail -1
