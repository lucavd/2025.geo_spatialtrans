#!/usr/bin/env bash
# R/testing/run_S1.6.sh - S1.6: cancelli G-1 (test), G-2 (mutanti), G-4 (ambiente). Uso (root del repo):
#   bash R/testing/run_S1.6.sh test      -> test_S1.6.R e tools/test_calibration_fit.py (tutte PASS), poi la batteria di S1.4f (che include S1.4e,
#                                           S1.5e, S1.5d, S1.5c, S1.5b, S1.5, S1.4d, S1.1-S1.4, OKF): FAIL ammessi solo T13 di S1.5d, T12 di S1.5c,
#                                           E10 ed E10b di S1.5 (DEC-011); alla fine nessun file versionato modificato (BL-165)
#   bash R/testing/run_S1.6.sh mutants   -> mutanti di S1.6 (N1-N6), S1.5 (M1-M9) e S1.5b (M1-M8) come copie modificate del codice (tools/mutate.py)
#   bash R/testing/run_S1.6.sh env       -> G-4: Imports di DESCRIPTION = pacchetti usati in R/*.R; restore da zero del lockfile in una libreria vuota;
#                                           test_S1.6.R con quella libreria; controllo OKF
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S16_LOG:-/mnt/micron/geo_spatialtrans/S1.6/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
fails_of() { grep '^FAIL' "$1" 2>/dev/null | awk '{print $2}' | sort -u | tr '\n' ' ' | sed 's/ $//'; }
if [ "$mode" = test ]; then
  st0=$(git status --porcelain --untracked-files=no)
  Rscript --vanilla R/testing/test_S1.6.R > "$LOG/test_S1.6.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL|RESULT)' "$LOG/test_S1.6.log"; grep -q '^FAIL' "$LOG/test_S1.6.log" && fail=1
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/test_S1.6.log" || fail=1
  .venv/bin/python -B tools/test_calibration_fit.py > "$LOG/test_calibration_fit.log" 2>&1 || fail=1
  tail -1 "$LOG/test_calibration_fit.log"; grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/test_calibration_fit.log" || fail=1
  S14F_LOG="$LOG/S1.4f" bash R/testing/run_S1.4f.sh test > "$LOG/run_S1.4f.log" 2>&1
  grep -E '^GATE|PASS regressioni|FAIL' "$LOG/run_S1.4f.log" | sed 's/^/[S1.4f] /' | tail -8
  grep -q '^GATE G-1 PASS' "$LOG/run_S1.4f.log" || { echo "FAIL batteria di S1.4f (regressioni)"; fail=1; }
  E="$LOG/S1.4f/S1.4e/S1.5e"
  grep -E '^(PASS|FAIL) regressione' "$E/../run_S1.5e.log" 2>/dev/null | sed 's/^/[S1.5e] /'
  st1=$(git status --porcelain --untracked-files=no)
  if [ "$st0" = "$st1" ]; then echo "PASS nessun file versionato modificato dalla batteria (BL-165)"; else echo "FAIL la batteria ha modificato file versionati:"; diff <(echo "$st0") <(echo "$st1"); fail=1; fi
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.6.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | awk '{print $2}' | sort -u > "$LOG/base_fail.txt"
  n=0; MS=$(python3 -B tools/mutate.py R/testing/mutants_S1.6.py --list); nm=$(echo $MS | wc -w); declare -A TD
  for m in $MS; do
    T=$(mktemp -d); TD[$m]=$T
    if ! python3 -B tools/mutate.py R/testing/mutants_S1.6.py $m "$T/src" > /dev/null; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; touch "$LOG/skip_$m"; continue; fi
    rm -f "$LOG/skip_$m"; ( S16_RSRC="$T/src" Rscript --vanilla R/testing/test_S1.6.R > "$LOG/test_S1.6_$m.log" 2>&1 ) &
  done
  wait
  for m in $MS; do
    [ -f "$LOG/skip_$m" ] && continue
    newf=$(grep '^FAIL' "$LOG/test_S1.6_$m.log" | awk '{print $2}' | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/test_S1.6_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante S1.6 $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante S1.6 $m NON rilevato"; fail=1; fi
    rm -rf "${TD[$m]}"
  done
  echo "mutanti di S1.6 rilevati $n/$nm"
  for s in S1.5 S1.5b; do
    L2="$LOG/S1.4f/S1.4e/S1.5e/S1.5d/S1.5c/S1.5b"; [ "$s" = S1.5 ] && L2="$L2/S1.5"
    if [ "$s" = S1.5 ]; then S15_LOG="$L2" bash R/testing/run_S1.5.sh mutants > "$LOG/mutants_$s.log" 2>&1 || fail=1
    else S15B_LOG="$L2" bash R/testing/run_S1.5b.sh mutants > "$LOG/mutants_$s.log" 2>&1 || fail=1; fi
    grep -E 'mutante|rilevati|GATE' "$LOG/mutants_$s.log" | sed "s/^/[$s] /"
  done
  [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
elif [ "$mode" = env ]; then
  Rscript --vanilla tools/S1.6_env_check.R desc > "$LOG/env_desc.log" 2>&1 || fail=1; cat "$LOG/env_desc.log" | tail -6
  rm -rf "$LOG/restore_lib"; mkdir -p "$LOG/restore_lib"
  Rscript --vanilla tools/S1.6_env_check.R restore "$LOG/restore_lib" > "$LOG/env_restore.log" 2>&1 || fail=1; tail -6 "$LOG/env_restore.log"
  S16_LIB="$LOG/restore_lib" Rscript --vanilla tools/S1.6_env_check.R test "$LOG/restore_lib" > "$LOG/env_test.log" 2>&1 || fail=1; grep -E '^(FAIL|RESULT)' "$LOG/env_test.log"
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/env_test.log" || fail=1
  python3 -B tools/okf.py check --bundle docs > "$LOG/env_okf.log" 2>&1 || fail=1; tail -1 "$LOG/env_okf.log"
  [ $fail = 0 ] && echo "GATE G-4 PASS" || echo "GATE G-4 FAIL"
else echo "uso: run_S1.6.sh test|mutants|env"; exit 2; fi
exit $fail
