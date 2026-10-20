#!/usr/bin/env bash
# tools/run_S1.4d.sh — S1.4d: md5 degli ingressi, test e mutanti, banco, raster, finestre, analisi (rilanciabile).
# Uso (root del repo): bash tools/run_S1.4d.sh   -> log in /mnt/micron/geo_spatialtrans/S1.4d/*.log
# Worker: /mnt/micron/geo_spatialtrans/S1.4d/WORKERS (default 6; server condiviso con MMM, limite termico).
set -u
OUT=/mnt/micron/geo_spatialtrans/S1.4d; mkdir -p "$OUT"; [ -f "$OUT/WORKERS" ] || echo 6 > "$OUT/WORKERS"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans
md5_inputs() { ( cd $M && find S1.4/real S1.4/loo S1.4/null S1.4/perf_S.rds S1.4b/loo S1.4c/gen S1.4c/sim_out R3/real -type f \( -name '*.rds' -o -name '*.parquet' \) | sort
                 ls S1.4c/gen_in/*_K2.parquet ) | (cd $M && xargs md5sum) > "$1"
               ( md5sum results/R2/R2_rois_checked.csv results/R2/tissue_masks/*_valid_ds4.png results/R2b/R2b_windows.csv \
                 results/R2b/annotations/R2b_annotations_Luca.json results/R2c/R2c_dmax.csv results/R2c/R2c_posthoc_ci.csv \
                 results/R2c/R2c_k2_windows.csv docs/BIO_REFERENCES.md ) >> "$1"; }
[ -f "$OUT/md5_before.txt" ] || { md5_inputs "$OUT/md5_before.txt"; echo "md5 prima: $(wc -l < $OUT/md5_before.txt) file"; }
git show db4afea:R/04b4_derive_cell_geometry.R > "$OUT/04b4_db4afea.R"
S14D_LOG=$OUT/tests bash R/testing/run_S1.4d.sh test > "$OUT/gate_G1.log" 2>&1; tail -1 "$OUT/gate_G1.log"
S14D_LOG=$OUT/tests bash R/testing/run_S1.4d.sh mutants > "$OUT/gate_G2.log" 2>&1; tail -1 "$OUT/gate_G2.log"
W=$(cat "$OUT/WORKERS")
for st in gen legacy null; do Rscript --vanilla tools/S1.4d_bench.R $st "$OUT" > "$OUT/$st.log" 2>&1; tail -1 "$OUT/$st.log"; done
.venv/bin/python tools/S1.4c_raster.py sim "$OUT/gen_in_tissue" "$OUT/sim_out_tissue" "$W" > "$OUT/sim_tissue.log" 2>&1; tail -1 "$OUT/sim_tissue.log"
.venv/bin/python tools/S1.4c_raster.py sim "$OUT/gen_in_sd" "$OUT/sim_out_sd" "$W" > "$OUT/sim_sd.log" 2>&1; tail -1 "$OUT/sim_sd.log"
.venv/bin/python tools/S1.4d_windows.py > "$OUT/windows.log" 2>&1; tail -3 "$OUT/windows.log"
Rscript --vanilla tools/S1.4d_bench.R repro "$OUT" > "$OUT/repro.log" 2>&1; tail -1 "$OUT/repro.log"
[ -f "$OUT/perf/perf_S.rds" ] || { /usr/bin/time -v Rscript --vanilla tools/S1.4d_bench.R perf "$OUT" > "$OUT/perf.log" 2>&1; grep -E 'PERF|Maximum resident' "$OUT/perf.log"; }
md5_inputs "$OUT/md5_after.txt"
if cmp -s "$OUT/md5_before.txt" "$OUT/md5_after.txt"; then echo "md5 ingressi invariati: $(wc -l < $OUT/md5_after.txt) file"; else echo "md5 ingressi CAMBIATI"; diff "$OUT/md5_before.txt" "$OUT/md5_after.txt" | head; fi
Rscript --vanilla tools/S1.4d_analysis.R > "$OUT/analysis.log" 2>&1; tail -5 "$OUT/analysis.log"
echo RUN_DONE
