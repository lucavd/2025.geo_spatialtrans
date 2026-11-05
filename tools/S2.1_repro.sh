#!/usr/bin/env bash
# tools/S2.1_repro.sh - S2.1, cancello G-3: seconda esecuzione completa in un'altra cartella con un numero di processi diverso (3 invece di 5),
# poi md5 delle verita' geometriche e delle tabelle di risultato contro la prima. Uso (root del repo): tools/S2.1_repro.sh
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans; R1=$M/S2.1/run1; R2=$M/S2.1/run2; mkdir -p "$R2/results"; echo 3 > "$R2/WORKERS"
export S21_RUN=S2.1/run2 S21_RES=$R2/results
Rscript --vanilla tools/S2.1_panel.R run > "$R2/panel_run.log" 2>&1
Rscript --vanilla tools/S2.1_panel.R check > "$R2/panel_check.log" 2>&1
Rscript --vanilla tools/S2.1_panel.R desc > "$R2/panel_desc.log" 2>&1
for K in A1 A23 A4 A5 A6; do Rscript --vanilla tools/S2.1_cp.R run $K > "$R2/cp_$K.log" 2>&1; done
Rscript --vanilla tools/S2.1_cp.R collect > "$R2/cp_collect.log" 2>&1
.venv/bin/python tools/S2.1_desc_real.py "$R2/results/S2.1_desc_real_roi.csv" > "$R2/desc_real.log" 2>&1
n=0; bad=0
for f in "$R1"/grid/*.rds; do b=$(basename "$f"); n=$((n+1)); [ "$(md5sum < "$f")" = "$(md5sum < "$R2/grid/$b")" ] || { bad=$((bad+1)); echo "diverso: $b"; }; done
for t in S2.1_G2.csv S2.1_desc_sim_roi.csv S2.1_desc_real_roi.csv S2.1_CP1_roi.csv S2.1_CP1_grid.csv; do n=$((n+1))
  [ "$(md5sum < results/S2.1/$t)" = "$(md5sum < "$R2/results/$t")" ] || { bad=$((bad+1)); echo "diverso: $t"; }; done
if [ "$bad" -eq 0 ]; then echo "PASS G-3 riproducibilita': $n file identici (30 verita' geometriche e 5 tabelle)"; else echo "FAIL G-3: $bad file diversi su $n"; fi
