#!/usr/bin/env bash
# R/testing/run_S1.5.sh - S1.5: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5.sh test      -> test_S1.5.R; regressione degli step passati (esiti identici a quelli committati); test OKF
#   bash R/testing/run_S1.5.sh mutants   -> M1..M9 (copie modificate del codice, R/testing/mutants_S1.5.py, da S1.6): ognuno deve far fallire
#                                          almeno un'asserzione di test_S1.5.R che passa senza mutante
# Stampa PASS/FAIL per asserzione; ultima riga GATE ... PASS|FAIL; exit 1 se il cancello fallisce.
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15_LOG:-/mnt/micron/geo_spatialtrans/S1.5/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  if [ -n "$(git status --porcelain -- results/S1.1 results/S1.2 results/S1.3 results/S1.4)" ]; then
    echo "FAIL modifiche locali in results/S1.1-S1.4: fermo (revisione FND-10)"; exit 2; fi
  Rscript --vanilla R/testing/test_S1.5.R > "$LOG/test_S1.5.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.5.log"; tail -1 "$LOG/test_S1.5.log"
  # S1.6 (BL-165): le tabelle dei test di S1.1-S1.4 si scrivono nella cartella dei log, non nei file versionati
  GEO_TEST_OUT="$LOG/reg_S1.1" Rscript --vanilla R/testing/test_S1.1.R > "$LOG/test_S1.1.log" 2>&1; echo "test_S1.1 exit $?"
  GEO_TEST_OUT="$LOG/reg_S1.2" Rscript --vanilla R/testing/test_S1.2.R > "$LOG/test_S1.2.log" 2>&1; echo "test_S1.2 exit $?"
  GEO_TEST_OUT="$LOG/reg_S1.3" S13_BACKEND=geos S13_MUTANT=none Rscript --vanilla R/testing/test_S1.3.R > "$LOG/test_S1.3.log" 2>&1; echo "test_S1.3 exit $?"
  GEO_TEST_OUT="$LOG/reg_S1.4" S14_MUTANT=none Rscript --vanilla R/testing/test_S1.4.R > "$LOG/test_S1.4.log" 2>&1; echo "test_S1.4 exit $?"
  S14D_MUTANT=none Rscript --vanilla R/testing/test_S1.4d.R > "$LOG/test_S1.4d.log" 2>&1 || fail=1
  tail -1 "$LOG/test_S1.4d.log"
  Rscript --vanilla tools/S1.5_regress_cmp.R results/S1.1/S1.1_test_results.csv="$LOG/reg_S1.1/S1.1_test_results.csv" \
    results/S1.2/S1.2_test_results.csv="$LOG/reg_S1.2/S1.2_test_results.csv" results/S1.3/S1.3_test_geos_none.csv="$LOG/reg_S1.3/S1.3_test_geos_none.csv" \
    results/S1.4/S1.4_test_none.csv="$LOG/reg_S1.4/S1.4_test_none.csv" || fail=1
  python3 -B tools/test_okf.py > "$LOG/test_okf.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_okf.log" | tail -3; tail -1 "$LOG/test_okf.log"
  git status --porcelain -- results/S1.1 results/S1.2 results/S1.3 results/S1.4 > "$LOG/regress_git_status.txt"
  if [ -s "$LOG/regress_git_status.txt" ]; then echo "FAIL i test hanno modificato file versionati in results/S1.1-S1.4 (BL-165)"; fail=1; fi
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | cut -d' ' -f2 | sort -u > "$LOG/base_fail.txt"
  # S1.6 (BL-155): mutanti come copie modificate del codice (R/testing/mutants_S1.5.py), niente interruttori nel codice di produzione
  n=0; MS=$(python3 -B tools/mutate.py R/testing/mutants_S1.5.py --list); nm=$(echo $MS | wc -w)
  for m in $MS; do
    T=$(mktemp -d)/$m
    if ! python3 -B tools/mutate.py R/testing/mutants_S1.5.py $m "$T" > /dev/null; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; continue; fi
    S15_RSRC="$T" Rscript --vanilla R/testing/test_S1.5.R > "$LOG/test_S1.5_$m.log" 2>&1
    newf=$(grep '^FAIL' "$LOG/test_S1.5_$m.log" | cut -d' ' -f2 | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/test_S1.5_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
    rm -rf "$(dirname "$T")"
  done
  echo "mutanti rilevati $n/$nm"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5.sh test|mutants"; exit 2; fi
exit $fail
