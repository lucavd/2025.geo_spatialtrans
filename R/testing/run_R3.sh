#!/usr/bin/env bash
# R/testing/run_R3.sh — R3 completo (comando unico). Uso: bash R/testing/run_R3.sh [stage...]
# stage: synthetic mutants real null calib analysis figures repro (default: tutti in quest'ordine)
set -u
cd "$(dirname "$0")/../.." || exit 1
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
R="Rscript --vanilla"; PY=.venv/bin/python; mkdir -p results/R3 logs
STAGES=${*:-"synthetic mutants real null calib analysis figures repro"}
for s in $STAGES; do
  echo "== $s $(date -Is)"
  case $s in
    synthetic) $R R/testing/test_R3.R synthetic ;;
    mutants)   bash tools/R3_mutants.sh ;;
    real)      $R tools/R3_run.R real && $PY tools/R3_containment.py all > logs/R3_containment.log && $R R/testing/test_R3.R real ;;
    null)      $R tools/R3_run.R null ;;
    calib)     $R tools/R3_run.R calib ;;
    analysis)  $R tools/R3_analysis.R ;;
    figures)   $R tools/R3_figures.R ;;
    repro)     D=/mnt/micron/geo_spatialtrans/R3_rerun; rm -rf "$D"; mkdir -p "$D"
               $R tools/R3_run.R real "$D" > logs/R3_rerun.log 2>&1
               sed "s#^R3 = .*#R3 = \"$D\"#" tools/R3_containment.py > /tmp/R3_containment_rerun.py && $PY /tmp/R3_containment_rerun.py all >> logs/R3_rerun.log
               $R tools/R3_run.R null "$D" >> logs/R3_rerun.log 2>&1
               $R tools/R3_run.R calib "$D" >> logs/R3_rerun.log 2>&1
               $R R/testing/test_R3.R repro /mnt/micron/geo_spatialtrans/R3 "$D" ;;
  esac
done
