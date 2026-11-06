#!/usr/bin/env bash
# S2.3, cancello G-3: banco (preparazione e verifica di R-1) e pannello (variante main) rifatti in un'altra cartella con un numero di
# processi diverso (3 invece di 6); confronto dei contenuti (i tempi di calcolo esclusi). Uso: tools/S2.3_repro.sh -> results/S2.3/S2.3_G3.txt
set -uo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans/S2.3; R=$M/repro; mkdir -p $R/bench $R/panel; echo 3 > $R/bench/WORKERS; echo 3 > $R/panel/WORKERS
.venv/bin/python tools/S2.3_bench.py prep $R/bench /mnt/micron/geo_spatialtrans/S2.0/run1_wfix/params_flat.json 1 > $R/prep.log 2>&1
cp $M/bench/combos_verify.json $R/bench/
.venv/bin/python tools/S2.3_bench.py grid $R/bench $R/bench/combos_verify.json 10 5 > $R/grid.log 2>&1
ln -sfn $M/panel/hd $R/panel/hd
Rscript --vanilla tools/S2.3_panel.R main $R/panel 3 > $R/panel.log 2>&1
.venv/bin/python tools/S2.3_panel_stats.py $R/panel main > $R/stats.log 2>&1
.venv/bin/python tools/S2.3_repro_compare.py
