#!/usr/bin/env bash
# R/testing/run_S1.5b.sh - S1.5b: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5b.sh test      -> test_S1.5b.R, poi la batteria e la regressione di S1.5 (run_S1.5.sh test)
#   bash R/testing/run_S1.5b.sh mutants   -> M1..M7: ognuno deve far fallire almeno un'asserzione di test_S1.5b.R che passa senza mutante
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15B_LOG:-/mnt/micron/geo_spatialtrans/S1.5b/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  S15B_MUTANT=none Rscript --vanilla R/testing/test_S1.5b.R > "$LOG/test_S1.5b.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.5b.log"; tail -1 "$LOG/test_S1.5b.log"
  grep -q '^FAIL' "$LOG/test_S1.5b.log" && fail=1
  S15_LOG="$LOG/S1.5" bash R/testing/run_S1.5.sh test > "$LOG/run_S1.5.log" 2>&1
  grep -E '^FAIL|test_S1|GATE|identic|RISULTATO|PASS.*OKF' "$LOG/run_S1.5.log" | tail -12
  # S1.5 chiuso con 81/83: falliscono solo E10 ed E10b (DEC-011)
  nf=$(grep -c '^FAIL' "$LOG/S1.5/test_S1.5.log" || true); ne=$(grep -E '^FAIL E10b? ' "$LOG/S1.5/test_S1.5.log" | wc -l)
  echo "test_S1.5: FAIL $nf (attesi 2: E10, E10b; trovati fra questi $ne)"
  [ "$nf" = 2 ] && [ "$ne" = 2 ] || fail=1
  # le altre parti della batteria di S1.5 devono passare davvero (revisione RC-02): test_S1.4d, regressioni S1.1-S1.4, test OKF
  grep -q '^FAIL' "$LOG/S1.5/test_S1.4d.log" && { echo "FAIL test_S1.4d"; fail=1; }
  grep -qE '^RESULT [0-9]+/[0-9]+ PASS' "$LOG/S1.5/test_S1.4d.log" || { echo "FAIL test_S1.4d senza riga RESULT"; fail=1; }
  nreg=$(grep -c '^PASS regressione' "$LOG/run_S1.5.log" || true); [ "$nreg" = 4 ] || { echo "FAIL regressioni S1.1-S1.4: $nreg/4"; fail=1; }
  grep -q '^FAIL' "$LOG/S1.5/test_okf.log" && { echo "FAIL test OKF"; fail=1; }
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/S1.5/test_okf.log" || { echo "FAIL test OKF senza RESULT completo"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5b.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | cut -d' ' -f2 | sort -u > "$LOG/base_fail.txt"
  n=0
  for m in M1 M2 M3 M4 M5 M6 M7 M8; do
    S15B_MUTANT=$m Rscript --vanilla R/testing/test_S1.5b.R > "$LOG/test_S1.5b_$m.log" 2>&1
    newf=$(grep '^FAIL' "$LOG/test_S1.5b_$m.log" | cut -d' ' -f2 | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
  done
  echo "mutanti rilevati $n/8"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5b.sh test|mutants"; exit 2; fi
exit $fail
