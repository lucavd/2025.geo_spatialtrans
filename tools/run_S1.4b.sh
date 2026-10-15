#!/usr/bin/env bash
# tools/run_S1.4b.sh — S1.4b: tutti gli stadi ufficiali (rilanciabile: le unita' gia' scritte si saltano).
# Uso (root del repo): bash tools/run_S1.4b.sh   -> log in /mnt/micron/geo_spatialtrans/S1.4b/*.log
set -u
OUT=/mnt/micron/geo_spatialtrans/S1.4b; mkdir -p "$OUT"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
for m in none MU1 MU2 MU3 MU4; do S14B_MUTANT=$m Rscript --vanilla R/testing/test_S1.4b.R > "$OUT/test_$m.log" 2>&1; tail -1 "$OUT/test_$m.log"; done
for st in prep loo mut; do Rscript --vanilla tools/S1.4b_bench.R $st "$OUT" > "$OUT/$st.log" 2>&1; tail -2 "$OUT/$st.log"; done
.venv/bin/python tools/S1.4b_arbiter.py "$OUT/arb_in" results/S1.4b/S1.4b_arbiter.csv > "$OUT/arbiter.log" 2>&1; tail -1 "$OUT/arbiter.log"
.venv/bin/python tools/S1.4b_arbiter.py "$OUT/arb_in" results/S1.4b/S1.4b_arbiter_MU5.csv --mutant MU5 > "$OUT/arbiter_MU5.log" 2>&1; tail -1 "$OUT/arbiter_MU5.log"
.venv/bin/python tools/S1.4b_arbiter.py "$OUT/arb_in_mut" results/S1.4b/S1.4b_arbiter_mut.csv > "$OUT/arbiter_mut.log" 2>&1; tail -1 "$OUT/arbiter_mut.log"
.venv/bin/python tools/S1.4b_cp1.py all "$OUT/cp1" > "$OUT/cp1.log" 2>&1; tail -1 "$OUT/cp1.log"
/usr/bin/time -v Rscript --vanilla tools/S1.4b_bench.R perf "$OUT" > "$OUT/perf.log" 2>&1; grep -E 'PERF|Maximum resident|Elapsed' "$OUT/perf.log"
echo RUN_DONE
