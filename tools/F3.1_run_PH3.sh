#!/usr/bin/env bash
# tools/F3.1_run_PH3.sh - F3.1, analisi post hoc PH-3 (dichiarata dopo gli esiti; non cambia gli esiti): le 30 aree con la semina a compartimenti
# (livelli Q6 di S1.6), 1 seme per area (stessi semi di k = 0 del pannello principale), profili di S2.2; stessi metodi e misure. Rilanciabile.
set -uo pipefail
W=${1:-6}; R=/mnt/micron/geo_spatialtrans/F3.1/main
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
TAGS=$(ls /mnt/micron/geo_spatialtrans/S1.6/layers | grep _U6.rds | sed 's/_U6.rds//' | paste -sd,)
F31_LAYER=Q6 Rscript --vanilla tools/F3.1_sim.R $R/simq6_rds 20261114 $TAGS $W 0 > $R/logs_sim_q6.log 2>&1
one() {  # $1 rds  $2 out_root
  local f=$1 o=$2; local b=$(basename $f .rds); local T=${b%_k*}; local A=${T%%_*}; local d=$o/$b
  if [ ! -f $d/export.json ]; then
    local SIDE=$(.venv/bin/python -c "import pandas as pd; r=pd.read_csv('results/R2/R2_rois_checked.csv'); r=r[(r.archetype+'_'+r.roi_id)=='$T'].iloc[0]; print(r.side_px*r.um_per_px)")
    Rscript --vanilla tools/F3.1_export_sim.R $f $d $SIDE > /dev/null 2>&1 || { echo "EXPORT FAIL $b"; return; }
    .venv/bin/python tools/F3.1_reg_sim.py $d > /dev/null 2>&1 || { echo "REG FAIL $b"; return; }
  fi
  for L in label label_shift; do [ -f $d/b2c_$L.json ] || .venv-b2c/bin/python tools/F3.1_b2c.py $d $R/ref $A sim $L > $d/b2c_$L.log 2>&1 || echo "B2C FAIL $b $L"; done
  [ -f $d/rctd.json ] || Rscript --vanilla tools/F3.1_rctd.R $d $R/ref $A sim > $d/rctd.log 2>&1 || echo "RCTD FAIL $b"
  echo "done $o $b"
}
export -f one; export R
ls $R/simq6_rds/*.rds | sed "s#\$# $R/simq6#" > $R/jobs_PH3.txt
xargs -a $R/jobs_PH3.txt -P $W -L 1 bash -c 'one $0 $1' > $R/logs_PH3.log 2>&1
grep -c "^done" $R/logs_PH3.log; grep -v "^done" $R/logs_PH3.log | head -20 || true; grep -h "ERROR" $R/logs_sim_q6.log | head || true
