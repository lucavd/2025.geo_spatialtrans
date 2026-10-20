#!/usr/bin/env bash
# R/testing/run_S1.4d.sh — S1.4d: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.4d.sh test      -> test_S1.4d.R, test_S1.4.R (modo "scale", esiti = quelli committati), test_okf.py
#   bash R/testing/run_S1.4d.sh mutants   -> M1..M7: ognuno deve far fallire test_S1.4d.R
# Stampa PASS/FAIL per asserzione; ultima riga GATE ... PASS|FAIL; exit 1 se il cancello fallisce.
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S14D_LOG:-/mnt/micron/geo_spatialtrans/S1.4d/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  S14D_MUTANT=none Rscript --vanilla R/testing/test_S1.4d.R > "$LOG/test_S1.4d.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.4d.log"; tail -1 "$LOG/test_S1.4d.log"
  S14_MUTANT=none Rscript --vanilla R/testing/test_S1.4.R > "$LOG/test_S1.4.log" 2>&1
  t=$(tail -1 "$LOG/test_S1.4.log"); echo "$t"
  case "$t" in *" 0 FAIL "*) echo "PASS S1.4 (modo scale) senza FAIL";; *) echo "FAIL S1.4 (modo scale)"; fail=1;; esac
  if git diff --quiet -- results/S1.4; then echo "PASS S1.4 esiti riscritti identici a quelli committati"; else echo "FAIL S1.4 esiti diversi da quelli committati"; git diff --stat -- results/S1.4; fail=1; fi
  python3 -B tools/test_okf.py > "$LOG/test_okf.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_okf.log"; tail -1 "$LOG/test_okf.log"
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  n=0
  for m in M1 M2 M3 M4 M5 M6 M7; do
    S14D_MUTANT=$m Rscript --vanilla R/testing/test_S1.4d.R > "$LOG/test_S1.4d_$m.log" 2>&1; rc=$?
    nf=$(grep -c '^FAIL' "$LOG/test_S1.4d_$m.log" || true)
    if [ $rc -ne 0 ] && [ "${nf:-0}" -gt 0 ]; then echo "PASS mutante $m rilevato ($nf FAIL)"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato (exit $rc)"; fail=1; fi
  done
  echo "mutanti rilevati $n/7"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.4d.sh test|mutants"; exit 2; fi
exit $fail
