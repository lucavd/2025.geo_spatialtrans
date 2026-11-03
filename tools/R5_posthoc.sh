#!/usr/bin/env bash
# tools/R5_posthoc.sh <outdir> — R5, analisi post hoc (dichiarate dopo i risultati) sulla stessa esecuzione, per G-4 anche sul post hoc.
# PH-1 sintetico CP-3; PH-2 sfocatura interna (griglia + 0.25 µm; scarto 0 in intestino e cuore; coda della sezione); PH-2s controprova sintetica;
# PH-3 intervallo raggiungibile; PH-4 regola R-4 (maschere spostate dello scarto); PH-5 bordi con nucleo stretto (0.7 µm). Processi: <outdir>/WORKERS.
set -euo pipefail
cd "$(dirname "$0")/.."
O=$1; NW=$(cat "$O/WORKERS"); PY=.venv/bin/python; L="$O/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
A16="A1_r1,A1_r2,A1_r3,A1_r4,A1_r5,A6_r1,A6_r2,A6_r3,A6_r4,A6_r5"
A2345="A2_r1,A2_r2,A2_r3,A2_r4,A2_r5,A3_r1,A3_r2,A3_r3,A3_r4,A3_r5,A4_f1,A4_f2,A4_p1,A4_p2,A4_p3,A5_r1,A5_r2,A5_r3,A5_r4,A5_r5"
$PY tools/R5_ph_cp3_synth.py "$NW" "$O/tables/R5_PH1_cp3_synth.csv" > "$L/ph1.log" 2>&1
$PY tools/R5_ph2_inner.py "$O" "$NW" > "$L/ph2.log" 2>&1
R5_PH2_GRID=0.25 $PY tools/R5_ph2_inner.py "$O" "$NW" >> "$L/ph2.log" 2>&1
R5_PH2_DELTA0=1 R5_PH2_TAGS=$A16 $PY tools/R5_ph2_inner.py "$O" "$NW" > "$L/ph2_d0.log" 2>&1
R5_PH2_TAIL=1 R5_PH2_TAGS=$A2345 $PY tools/R5_ph2_inner.py "$O" "$NW" > "$L/ph2_tail.log" 2>&1
$PY tools/R5_ph2_synth.py "$O" "$NW" > "$L/ph2_synth.log" 2>&1
$PY tools/R5_ph3_reach.py "$O" "$NW" > "$L/ph3.log" 2>&1
$PY tools/R5_ph4_r4shift.py "$O" > "$L/ph4.log" 2>&1
R5_PH5_SIG=0.7 $PY tools/R5_ph5_narrowcore.py "$O" "$NW" > "$L/ph5.log" 2>&1
for sub in ph2 ph2_d0 ph2_tail; do $PY tools/R5_ph2_summary.py "$O" $sub > "$L/summary_$sub.log" 2>&1; done
$PY tools/R5_ph2_synth_summary.py "$O" > "$L/summary_ph2_synth.log" 2>&1
$PY tools/R5_ph4_summary.py "$O" > "$L/summary_ph4.log" 2>&1
$PY tools/R5_ph5_summary.py "$O" 0.7 > "$L/summary_ph5.log" 2>&1
echo "R5_posthoc done $(date -Is)" >> "$L/status.log"
