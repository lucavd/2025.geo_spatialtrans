#!/usr/bin/env bash
# tools/run_S1.4f.sh - S1.4f: catena del calcolo (docs/steps/S1.4f.md). Uso (root del repo): bash tools/run_S1.4f.sh
# Stadi rilanciabili (ogni unita' gia' scritta e' saltata). Worker: /mnt/micron/geo_spatialtrans/S1.4f/WORKERS (R) e NPY (Python).
# Ordine pre-registrato: taratura -> lambda -> banco B (tre leggi) -> regola O -> solo allora raster, contatti e coppie di EP, banco A.
set -eu
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=/mnt/micron/geo_spatialtrans/S1.4f; mkdir -p "$O/logs"; NPY=${NPY:-20}
RS="Rscript --vanilla tools/S1.4f_bench.R"; PY=.venv/bin/python
st() { echo "== $1 $(date +%H:%M:%S)"; }
[ -f "$O/md5_before.csv" ] || { st md5before; $RS md5before > "$O/logs/md5before.log" 2>&1; tail -1 "$O/logs/md5before.log"; }
st calB;      $RS calB > "$O/logs/calB.log" 2>&1; tail -2 "$O/logs/calB.log"
st shapecal;  [ -f results/S1.4f/S1.4f_shape_cal.csv ] || $PY tools/S1.4f_measure.py shape cal $NPY > "$O/logs/shape_cal.log" 2>&1; tail -1 "$O/logs/shape_cal.log"
st calibrate; [ -f results/S1.4f/S1.4f_lambda.csv ] || $PY tools/S1.4f_calib.py calibrate > "$O/logs/calibrate.log" 2>&1; tail -40 "$O/logs/calibrate.log"
st genB;      $RS genB > "$O/logs/genB.log" 2>&1; tail -2 "$O/logs/genB.log"
st shapeB;    [ -f results/S1.4f/S1.4f_shape_B.csv ] || $PY tools/S1.4f_measure.py shape B $NPY > "$O/logs/shape_B.log" 2>&1; tail -1 "$O/logs/shape_B.log"
st ruleO;     [ -f results/S1.4f/S1.4f_ruleO.csv ] || $PY tools/S1.4f_calib.py ruleO > "$O/logs/ruleO.log" 2>&1; tail -6 "$O/logs/ruleO.log"
W=$($PY -c "import pandas as pd; print(pd.read_csv('results/S1.4f/S1.4f_ruleO.csv').winner.iloc[0])"); echo "EP = $W"
ln -sfn "$W" "$O/B/nuc/EP"
st rasterEP;  $PY tools/S1.4c_raster.py sim "$O/B/nuc/EP" "$O/B/raster/EP" $NPY > "$O/logs/raster_EP.log" 2>&1; tail -1 "$O/logs/raster_EP.log"
st pairsB;    $PY tools/S1.4f_measure.py pairs B $NPY > "$O/logs/pairsB.log" 2>&1; tail -1 "$O/logs/pairsB.log"
st windows;   $PY tools/S1.4f_measure.py windows > "$O/logs/windows.log" 2>&1; tail -8 "$O/logs/windows.log"
st genA;      $RS genA > "$O/logs/genA.log" 2>&1; tail -2 "$O/logs/genA.log"
for b in A A_dd; do st pairs$b; $PY tools/S1.4f_measure.py pairs $b $NPY > "$O/logs/pairs$b.log" 2>&1; tail -1 "$O/logs/pairs$b.log"; done
st cp1;       $PY tools/S1.4f_measure.py cp1 > "$O/logs/cp1.log" 2>&1; tail -14 "$O/logs/cp1.log"
st g4;        $PY tools/S1.4f_measure.py g4 > "$O/logs/g4.log" 2>&1; tail -25 "$O/logs/g4.log"
st check;     $RS check > "$O/logs/check.log" 2>&1; tail -16 "$O/logs/check.log"
st md5after;  $RS md5after > "$O/logs/md5after.log" 2>&1; tail -2 "$O/logs/md5after.log"
st analysis;  $PY tools/S1.4f_analysis.py > "$O/logs/analysis.log" 2>&1 || true; tail -80 "$O/logs/analysis.log"
echo "catena S1.4f finita $(date +%H:%M:%S)"
