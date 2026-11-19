#!/usr/bin/env bash
# tools/F3.3_sim_cg.sh - F3.3: bin2cell, RCTD e misure sulle 20 simulazioni con i geni comuni (ramo condizionale di DEC-097, punto 4)
set -uo pipefail; cd ~/2025.geo_spatialtrans; W=${1:-4}
M=/mnt/micron/geo_spatialtrans; REF=$M/F3.1/main/ref; mkdir -p $M/F3.3/sim_cg_measures
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
.venv/bin/python tools/F3.3_sim_cg.py prepare > $M/F3.3/log_sim_cg.log 2>&1
one() { local d=$1; local b=$(basename $d); local A=${b:0:2}
  [ -f $d/b2c_label.json ] || .venv-b2c/bin/python tools/F3.1_b2c.py $d $REF $A sim label > $d/b2c_label.log 2>&1 || { echo "B2C FAIL $b"; return; }
  [ -f $d/rctd.json ] || Rscript --vanilla tools/F3.1_rctd.R $d $REF $A sim > $d/rctd.log 2>&1 || { echo "RCTD FAIL $b"; return; }
  local o=$M/F3.3/sim_cg_measures/$b.json; [ -s $o ] || .venv-b2c/bin/python tools/F3.3_calib_measures.py $d sim $o > $o.log 2>&1 || echo "MEASURE FAIL $b"
  echo "done $b"; }
export -f one; export M REF
ls -d $M/F3.3/sim_cg/*/ | sed 's#/$##' | xargs -P $W -L 1 bash -c 'one $0' > $M/F3.3/log_sim_cg_runs.log 2>&1
grep -c '^done' $M/F3.3/log_sim_cg_runs.log
