#!/usr/bin/env bash
# R/testing/run_S1.5.sh - S1.5: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5.sh test      -> test_S1.5.R; regressione degli step passati (esiti identici a quelli committati); test OKF
#   bash R/testing/run_S1.5.sh mutants   -> M1..M9: ognuno deve far fallire almeno un'asserzione di test_S1.5.R che passa senza mutante
# Stampa PASS/FAIL per asserzione; ultima riga GATE ... PASS|FAIL; exit 1 se il cancello fallisce.
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15_LOG:-/mnt/micron/geo_spatialtrans/S1.5/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  if [ -n "$(git status --porcelain -- results/S1.1 results/S1.2 results/S1.3 results/S1.4)" ]; then
    echo "FAIL modifiche locali in results/S1.1-S1.4: il test le sovrascriverebbe; fermo (revisione FND-10)"; exit 2; fi
  S15_MUTANT=none Rscript --vanilla R/testing/test_S1.5.R > "$LOG/test_S1.5.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.5.log"; tail -1 "$LOG/test_S1.5.log"
  Rscript --vanilla R/testing/test_S1.1.R > "$LOG/test_S1.1.log" 2>&1; echo "test_S1.1 exit $?"
  Rscript --vanilla R/testing/test_S1.2.R > "$LOG/test_S1.2.log" 2>&1; echo "test_S1.2 exit $?"
  S13_BACKEND=geos S13_MUTANT=none Rscript --vanilla R/testing/test_S1.3.R > "$LOG/test_S1.3.log" 2>&1; echo "test_S1.3 exit $?"
  S14_MUTANT=none Rscript --vanilla R/testing/test_S1.4.R > "$LOG/test_S1.4.log" 2>&1; echo "test_S1.4 exit $?"
  S14D_MUTANT=none Rscript --vanilla R/testing/test_S1.4d.R > "$LOG/test_S1.4d.log" 2>&1 || fail=1
  tail -1 "$LOG/test_S1.4d.log"
  Rscript --vanilla tools/S1.5_regress_cmp.R results/S1.1/S1.1_test_results.csv results/S1.2/S1.2_test_results.csv \
    results/S1.3/S1.3_test_geos_none.csv results/S1.4/S1.4_test_none.csv || fail=1
  python3 -B tools/test_okf.py > "$LOG/test_okf.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_okf.log" | tail -3; tail -1 "$LOG/test_okf.log"
  git status --porcelain -- results/S1.1 results/S1.2 results/S1.3 results/S1.4 > "$LOG/regress_git_status.txt"
  git checkout -- results/S1.1 results/S1.2 results/S1.3 results/S1.4   # gli esiti committati restano il riferimento
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | cut -d' ' -f2 | sort -u > "$LOG/base_fail.txt"
  n=0
  for m in M1 M2 M3 M4 M5 M6 M7 M8 M9; do
    S15_MUTANT=$m Rscript --vanilla R/testing/test_S1.5.R > "$LOG/test_S1.5_$m.log" 2>&1
    newf=$(grep '^FAIL' "$LOG/test_S1.5_$m.log" | cut -d' ' -f2 | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
  done
  echo "mutanti rilevati $n/9"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5.sh test|mutants"; exit 2; fi
exit $fail
