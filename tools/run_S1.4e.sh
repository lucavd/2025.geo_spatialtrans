#!/usr/bin/env bash
# tools/run_S1.4e.sh - S1.4e: catena del calcolo (docs/steps/S1.4e.md). Uso (root del repo): bash tools/run_S1.4e.sh
# Stadi rilanciabili (ogni unita' gia' scritta e' saltata). Worker: /mnt/micron/geo_spatialtrans/S1.4e/WORKERS (R) e NPY (Python).
set -eu
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=/mnt/micron/geo_spatialtrans/S1.4e; mkdir -p "$O/logs"; NPY=${NPY:-20}
RS="Rscript --vanilla tools/S1.4e_bench.R"; PY=.venv/bin/python
[ -f "$O/md5_before.csv" ] || $RS md5before > "$O/logs/md5before.log" 2>&1
for st in genB genA; do $RS $st > "$O/logs/$st.log" 2>&1; tail -2 "$O/logs/$st.log"; done
for a in K2 ER ET; do $PY tools/S1.4c_raster.py sim "$O/B/nuc/$a" "$O/B/raster/$a" $NPY > "$O/logs/raster_$a.log" 2>&1; tail -1 "$O/logs/raster_$a.log"; done
$PY tools/S1.4e_measure.py pairs B $NPY > "$O/logs/pairsB.log" 2>&1; tail -1 "$O/logs/pairsB.log"
$PY tools/S1.4e_measure.py pairs A $NPY > "$O/logs/pairsA.log" 2>&1; tail -1 "$O/logs/pairsA.log"
$PY tools/S1.4e_measure.py windows > "$O/logs/windows.log" 2>&1; tail -8 "$O/logs/windows.log"
$PY tools/S1.4e_measure.py shape $NPY > "$O/logs/shape.log" 2>&1; tail -8 "$O/logs/shape.log"
$PY tools/S1.4e_measure.py cp1 > "$O/logs/cp1.log" 2>&1; tail -8 "$O/logs/cp1.log"
$PY tools/S1.4e_measure.py g4 > "$O/logs/g4.log" 2>&1; tail -6 "$O/logs/g4.log"
$RS check > "$O/logs/check.log" 2>&1; tail -10 "$O/logs/check.log"
$RS md5after > "$O/logs/md5after.log" 2>&1; tail -2 "$O/logs/md5after.log"
echo "catena S1.4e finita"
