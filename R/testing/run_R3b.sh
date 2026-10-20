#!/usr/bin/env bash
# R/testing/run_R3b.sh — R3b: cancelli G-1 (test) e G-2 (mutanti) del costruttore della maschera CO.
# Uso (root del repo): bash R/testing/run_R3b.sh test | mutants     (log per esecuzione in ${R3B_LOG:-/tmp/R3b_tests})
set -u
LOG=${R3B_LOG:-/tmp/R3b_tests}; mkdir -p "$LOG"
PY=.venv/bin/python
case "${1:-test}" in
  test)
    R3B_MUTANT=none $PY -B tools/R3b_consensus_test.py > "$LOG/test_none.log" 2>&1; rc=$?
    grep -E '^(PASS|FAIL)' "$LOG/test_none.log"; tail -1 "$LOG/test_none.log"
    [ $rc -eq 0 ] && echo "G-1: PASS" || echo "G-1: FAIL" ;;
  mutants)
    det=0; tot=0
    for M in M1 M2 M3 M4 M5; do
      tot=$((tot+1)); R3B_MUTANT=$M $PY -B tools/R3b_consensus_test.py > "$LOG/test_$M.log" 2>&1
      nf=$(grep -c '^FAIL' "$LOG/test_$M.log" || true)
      if [ "${nf:-0}" -gt 0 ]; then det=$((det+1)); echo "$M rilevato ($nf asserzioni FAIL: $(grep '^FAIL' "$LOG/test_$M.log" | cut -d' ' -f2-3 | tr '\n' ';'))"; else echo "$M NON rilevato"; fi
    done
    echo "G-2: $det/$tot mutanti rilevati" ;;
esac
# mutanti del revisore (dopo la revisione avversariale): bash R/testing/run_R3b.sh reviewer
if [ "${1:-}" = "reviewer" ]; then
  for f in /mnt/micron/geo_spatialtrans/R3b/review_code/mutants/X*.py; do
    R3B_MUTANT=none R3B_CONSENSUS_FILE=$f $PY -B tools/R3b_consensus_test.py > "$LOG/rv_$(basename $f .py).log" 2>&1
    nf=$(grep -c '^FAIL' "$LOG/rv_$(basename $f .py).log" || true); echo "$(basename $f .py): ${nf:-0} FAIL"
  done
fi
