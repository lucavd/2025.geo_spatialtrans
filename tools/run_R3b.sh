#!/usr/bin/env bash
# tools/run_R3b.sh — R3b: catena completa (rilanciabile; le unita' gia' scritte si saltano). Uso (root del repo):
#   bash tools/run_R3b.sh            tutto: maschera CO, G-4, catena R3/S1.4/S1.4d, CP-1, raster, finestre, P-1/P-2, analisi, post hoc
#   bash tools/run_R3b.sh check      solo il cancello G-4 sulla maschera CO (ricalcolo dai file)
#   bash tools/run_R3b.sh regress    solo la regressione G-3 (parte R; le parti Python sono scritte da R3b_p1p2.py e R3b_windows.py)
# Log in /mnt/micron/geo_spatialtrans/R3b/*.log; worker nel file WORKERS (default 6; server condiviso con MMM, limite termico).
set -u
O=/mnt/micron/geo_spatialtrans/R3b; mkdir -p "$O"; [ -f "$O/WORKERS" ] || echo 6 > "$O/WORKERS"; W=$(cat "$O/WORKERS")
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=".venv/bin/python -B"; RS="Rscript --vanilla"
case "${1:-all}" in
  check)   $PY tools/R3b_consensus.py check "$W" > "$O/check.log" 2>&1; tail -1 "$O/check.log"; exit 0 ;;
  regress) $RS tools/R3b_bench.R regress > "$O/regress.log" 2>&1; tail -25 "$O/regress.log"; exit 0 ;;
esac
$PY tools/R3b_consensus.py build "$W" > "$O/build.log" 2>&1; tail -1 "$O/build.log"
$PY tools/R3b_consensus.py check "$W" > "$O/check.log" 2>&1; tail -1 "$O/check.log"
$RS tools/R3b_bench.R r3 > "$O/r3.log" 2>&1; tail -1 "$O/r3.log"
$PY tools/R3b_nuc.py > "$O/nuc.log" 2>&1; tail -1 "$O/nuc.log"
for st in real loo phi gen regress cp1; do $RS tools/R3b_bench.R $st > "$O/$st.log" 2>&1; tail -1 "$O/$st.log"; done
for s in consensus_own stardist_he_own stardist_he_sr cp1; do $PY tools/S1.4c_raster.py sim "$O/gen_in_$s" "$O/sim_out_$s" "$W" > "$O/sim_$s.log" 2>&1; tail -1 "$O/sim_$s.log"; done
$PY tools/R3b_windows.py > "$O/windows.log" 2>&1; tail -2 "$O/windows.log"
$PY tools/R3b_p1p2.py > "$O/p1p2.log" 2>&1; tail -2 "$O/p1p2.log"
$RS tools/R3b_export.R > "$O/export.log" 2>&1; tail -1 "$O/export.log"
$PY tools/R3b_analysis.py > "$O/analysis.log" 2>&1; tail -6 "$O/analysis.log"
# post hoc (dopo i verdetti)
$PY tools/R3b_posthoc.py > "$O/posthoc.log" 2>&1
$RS tools/R3b_posthoc_area.R > "$O/posthoc_area.log" 2>&1
$RS tools/R3b_posthoc_sdarea.R posthoc > "$O/ph_sdarea.log" 2>&1
$PY tools/S1.4c_raster.py sim "$O/gen_in_consensus_sdarea" "$O/sim_out_consensus_sdarea" "$W" > "$O/sim_sdarea.log" 2>&1
$PY tools/R3b_posthoc_sdarea_windows.py > "$O/ph_sdarea_windows.log" 2>&1
$RS tools/R3b_refs.R > "$O/refs.log" 2>&1
echo RUN_DONE
