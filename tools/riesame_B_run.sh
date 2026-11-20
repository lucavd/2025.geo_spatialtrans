#!/usr/bin/env bash
# tools/riesame_B_run.sh - riesame, Parte B: simulazioni del braccio B (tools/riesame_B_sim.R, 10 aree x 2 semi, seme base 20261124) ed
# esportazione nel formato di S (tools/F3.1_export_sim.R e tools/F3.1_reg_sim.py invariati; lato dell'area uguale a quello
# dell'esportazione di S della stessa area). Poi la catena: bash tools/riesame_panel_run.sh <B>/jobs/panel_B.tsv 2 (e common_B.tsv).
# Rilanciabile: salta cio' che esiste. Uso: bash tools/riesame_B_run.sh [W=2] [tags=tutte] [k=0,1]
set -uo pipefail
cd ~/2025.geo_spatialtrans
W=${1:-2}; TAGS=${2:-A2_r1,A2_r2,A2_r3,A2_r4,A2_r5,A3_r1,A3_r2,A3_r3,A3_r4,A3_r5}; KS=${3:-0,1}
M=/mnt/micron/geo_spatialtrans; R=$M/riesame/B; mkdir -p $R/sim_rds $R/sim
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla tools/riesame_B_sim.R $R/sim_rds 20261124 $TAGS $W $KS >> $R/log_B_sim.log 2>&1
for f in $R/sim_rds/*_k[01].rds; do
  b=$(basename $f .rds); T=${b%_k*}; d=$R/sim/$b
  [ -f $d/export.json ] && python3 -c "import json,sys; sys.exit(0 if 'n_reg' in json.load(open('$d/export.json')) else 1)" && continue
  SIDE=$(python3 -c "import json; print(json.load(open('$M/F3.1/main/sim/${T}_k0/export.json'))['side_um'])")
  Rscript --vanilla tools/F3.1_export_sim.R $f $d $SIDE >> $R/log_B_export.log 2>&1 || { echo "EXPORT FAIL $b"; continue; }
  .venv/bin/python tools/F3.1_reg_sim.py $d >> $R/log_B_export.log 2>&1 || echo "REG FAIL $b"
  echo "exported $b"
done
