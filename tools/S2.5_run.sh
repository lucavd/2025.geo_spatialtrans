#!/usr/bin/env bash
# tools/S2.5_run.sh - S2.5: calcolo pre-registrato (docs/steps/S2.5.md). Fasi rilanciabili (un file per unita').
#   tools/S2.5_run.sh main   : md5 prima, G-1 (test + mutanti), lato vero, pannello main e depth_off, G-2, misure, esiti, regola R-2
#   tools/S2.5_run.sh r2     : pannello main_r2, G-2, misure, esiti
set -euo pipefail
PH=${1:?main o r2}
M=/mnt/micron/geo_spatialtrans/S2.5; RD=$M/real; PD=$M/panel; RS=results/S2.5
mkdir -p $RD $PD $RS; echo 6 > $PD/WORKERS; echo 3 > $RD/WORKERS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
log() { echo "[$(date +%H:%M:%S)] $*"; }
if [ "$PH" = main ]; then
  [ -s $RS/S2.5_md5_before.txt ] || tools/S2.5_md5.sh before $RS
  log "G-1 test"; R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.5.R > $RS/S2.5_G1_tests.txt 2>&1 || true
  ( .venv/bin/python -B tools/S2.5_mutants.py $RS/S2.5_mutants.txt > $M/mutants.log 2>&1 || true ) &
  MUTPID=$!
  log "lato vero"; .venv/bin/python -B tools/S2.5_real.py areas $RD > $RD/areas.log 2>&1
  .venv/bin/python -B tools/S2.5_real.py summary $RD > $RD/summary.log 2>&1; .venv/bin/python -B tools/S2.5_real.py bl219 $RD > $RD/bl219.log 2>&1
  cp $RD/S2.5_real_areas.csv $RD/S2.5_real_refs.csv $RD/S2.5_real_R4_reproduction.csv $RD/S2.5_BL219_share.csv $RD/S2.5_BL219_summary.json $RS/
  for V in main depth_off; do log "pannello $V"; Rscript --vanilla tools/S2.5_panel.R $V $PD 6 > $PD/run_$V.log 2>&1; done
  log "G-2"; Rscript --vanilla tools/S2.5_panel.R check $PD > $RS/S2.5_G2.txt 2>&1 || true; cp $PD/G2.csv $RS/S2.5_G2.csv
  for V in main depth_off; do log "misure $V"; .venv/bin/python -B tools/S2.5_stats.py $PD $V > $PD/stats_$V.log 2>&1; cp $PD/S2.5_sim_$V.csv $RS/; done
  for V in main depth_off; do .venv/bin/python -B tools/S2.5_analysis.py verdicts $RD $PD $V $RS > $RS/S2.5_analysis_$V.txt 2>&1; done
  .venv/bin/python -B tools/S2.5_analysis.py r2 $RD $PD $RS > $RS/S2.5_R2.txt 2>&1
  wait $MUTPID || true
  log "fine fase main"
elif [ "$PH" = r2 ]; then
  log "pannello main_r2"; Rscript --vanilla tools/S2.5_panel.R main_r2 $PD 6 > $PD/run_main_r2.log 2>&1
  log "G-2"; Rscript --vanilla tools/S2.5_panel.R check $PD > $RS/S2.5_G2.txt 2>&1 || true; cp $PD/G2.csv $RS/S2.5_G2.csv
  .venv/bin/python -B tools/S2.5_stats.py $PD main_r2 > $PD/stats_main_r2.log 2>&1; cp $PD/S2.5_sim_main_r2.csv $RS/
  .venv/bin/python -B tools/S2.5_analysis.py verdicts $RD $PD main_r2 $RS > $RS/S2.5_analysis_main_r2.txt 2>&1
  log "fine fase r2"
fi
