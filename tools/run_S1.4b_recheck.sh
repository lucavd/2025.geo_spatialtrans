#!/usr/bin/env bash
# tools/run_S1.4b_recheck.sh — S1.4b: dopo le correzioni D-4 (ricerca dei vicini per posizione) e D-5 (doppi vettorizzati):
# test, ricalcolo di prep+loo e confronto con gli ufficiali, arbitro ufficiale rieseguito, arbitro dei mutanti, post hoc PH-A/PH-B, C-perf.
set -u
OUT=/mnt/micron/geo_spatialtrans/S1.4b; RC=/mnt/micron/geo_spatialtrans/S1.4b_recheck; mkdir -p "$RC"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
for m in none MU1 MU2 MU3 MU4; do S14B_MUTANT=$m Rscript --vanilla R/testing/test_S1.4b.R > "$OUT/test2_$m.log" 2>&1; tail -1 "$OUT/test2_$m.log"; done
for st in prep loo; do Rscript --vanilla tools/S1.4b_bench.R $st "$RC" > "$RC/$st.log" 2>&1; tail -1 "$RC/$st.log"; done
Rscript --vanilla tools/S1.4b_recheck.R
.venv/bin/python tools/S1.4b_arbiter.py "$OUT/arb_in" /tmp/s14b_arb_rerun.csv > "$OUT/arbiter_rerun.log" 2>&1; tail -1 "$OUT/arbiter_rerun.log"
cmp results/S1.4b/S1.4b_arbiter.csv /tmp/s14b_arb_rerun.csv && echo ARB_RERUN_IDENTICAL
.venv/bin/python tools/S1.4b_arbiter.py "$OUT/arb_in_mut" results/S1.4b/S1.4b_arbiter_mut.csv > "$OUT/arbiter_mut.log" 2>&1; tail -1 "$OUT/arbiter_mut.log"
.venv/bin/python tools/S1.4b_posthoc_arbiter.py "$OUT/arb_in" results/S1.4b/S1.4b_arbiter.csv results/S1.4b/S1.4b_posthoc_arbiter > "$OUT/posthoc_arbiter.log" 2>&1; tail -3 "$OUT/posthoc_arbiter.log"
/usr/bin/time -v Rscript --vanilla tools/S1.4b_bench.R perf "$OUT" > "$OUT/perf.log" 2>&1; grep -E 'PERF|Maximum resident|Elapsed' "$OUT/perf.log"
echo RECHECK_DONE
