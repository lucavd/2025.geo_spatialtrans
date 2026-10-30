#!/usr/bin/env bash
# tools/R4_run.sh [outdir] — R4 completo: stadio per ROI e per sezione (5 + 5 processi), analisi, post hoc (5 processi), riassunto. Rilanciabile.
set -euo pipefail
OUT=${1:-/mnt/micron/geo_spatialtrans/R4/run1}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT/logs"
for k in A1 A23 A4 A5 A6; do
  .venv/bin/python tools/R4_roi.py $k "$OUT" > "$OUT/logs/roi_$k.log" 2>&1 &
  .venv/bin/python tools/R4_section.py $k "$OUT" > "$OUT/logs/section_$k.log" 2>&1 &
done
wait
for k in A1 A23 A4 A5 A6; do .venv/bin/python tools/R4_posthoc.py roi $k "$OUT" > "$OUT/logs/ph_roi_$k.log" 2>&1 & done
OMP_NUM_THREADS=4 .venv/bin/python tools/R4_analysis.py "$OUT" > "$OUT/logs/analysis.log" 2>&1
wait
.venv/bin/python tools/R4_posthoc.py summary "$OUT" > "$OUT/logs/ph_summary.log" 2>&1
echo "R4_run done: $OUT"
