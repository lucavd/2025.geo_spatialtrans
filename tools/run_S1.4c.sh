#!/usr/bin/env bash
# tools/run_S1.4c.sh — S1.4c: test, stadi ufficiali e analisi (rilanciabile: le unita' gia' scritte si saltano).
# Uso (root del repo): bash tools/run_S1.4c.sh   -> log in /mnt/micron/geo_spatialtrans/S1.4c/*.log
set -u
OUT=/mnt/micron/geo_spatialtrans/S1.4c; mkdir -p "$OUT"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
for m in none MU1 MU2 MU5; do S14C_MUTANT=$m Rscript --vanilla R/testing/test_S1.4c.R > "$OUT/test_R_$m.log" 2>&1; tail -1 "$OUT/test_R_$m.log"; done
Rscript --vanilla R/testing/test_S1.4c_decision.R > "$OUT/test_decision.log" 2>&1; tail -1 "$OUT/test_decision.log"
for m in none MU3 MU4; do S14C_MUTANT=$m .venv/bin/python tools/S1.4c_raster.py test > "$OUT/test_py_$m.log" 2>&1; tail -1 "$OUT/test_py_$m.log"; done
Rscript --vanilla tools/S1.4c_bench.R prep "$OUT" > "$OUT/prep.log" 2>&1; tail -1 "$OUT/prep.log"
.venv/bin/python tools/S1.4c_raster.py real "$OUT/real_shape" > "$OUT/real.log" 2>&1; tail -1 "$OUT/real.log"
.venv/bin/python tools/S1.4c_raster.py cp2 "$OUT/cp2" 20 > "$OUT/cp2.log" 2>&1; tail -1 "$OUT/cp2.log"
Rscript --vanilla tools/S1.4c_bench.R gen "$OUT" > "$OUT/gen.log" 2>&1; tail -1 "$OUT/gen.log"
.venv/bin/python tools/S1.4c_raster.py sim "$OUT/gen_in" "$OUT/sim_out" 30 > "$OUT/sim.log" 2>&1; tail -1 "$OUT/sim.log"
Rscript --vanilla tools/S1.4c_bench.R obs "$OUT" > "$OUT/obs.log" 2>&1; tail -1 "$OUT/obs.log"
Rscript --vanilla tools/S1.4c_bench.R mut "$OUT" > "$OUT/mut.log" 2>&1; tail -2 "$OUT/mut.log"
/usr/bin/time -v Rscript --vanilla tools/S1.4c_bench.R perf "$OUT" > "$OUT/perf.log" 2>&1; grep -E 'PERF|Maximum resident|Elapsed' "$OUT/perf.log"
echo RUN_DONE
