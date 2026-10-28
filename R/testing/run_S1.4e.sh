#!/usr/bin/env bash
# R/testing/run_S1.4e.sh - S1.4e: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.4e.sh test      -> test_S1.4e.R e tools/S1.4e_pairs_test.py (tutte PASS), poi la batteria di S1.5e
#                                            (che include S1.5d, S1.5c, S1.5b, S1.5, S1.4d, S1.1-S1.4, OKF con i soli FAIL noti): deve dire GATE G-1 PASS
#   bash R/testing/run_S1.4e.sh mutants   -> M1..M6 come copie modificate del codice: ognuno deve far fallire almeno un'asserzione che passa senza mutante
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S14E_LOG:-/mnt/micron/geo_spatialtrans/S1.4e/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  Rscript --vanilla R/testing/test_S1.4e.R > "$LOG/test_S1.4e.log" 2>&1 || fail=1
  .venv/bin/python tools/S1.4e_pairs_test.py > "$LOG/test_S1.4e_py.log" 2>&1 || fail=1
  for f in "$LOG/test_S1.4e.log" "$LOG/test_S1.4e_py.log"; do
    grep -E '^(PASS|FAIL|RESULT)' "$f"; grep -q '^FAIL' "$f" && fail=1; grep -qE '^RESULT ([0-9]+)/\1 PASS' "$f" || fail=1
  done
  S15E_LOG="$LOG/S1.5e" bash R/testing/run_S1.5e.sh test > "$LOG/run_S1.5e.log" 2>&1
  grep -E '^(PASS|FAIL) regressione|^GATE' "$LOG/run_S1.5e.log" | sed 's/^/[S1.5e] /'
  grep -q '^GATE G-1 PASS' "$LOG/run_S1.5e.log" && echo "PASS regressioni (batteria di S1.5e con i soli FAIL noti)" || { echo "FAIL regressioni"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  for f in test_S1.4e.log test_S1.4e_py.log; do [ -f "$LOG/$f" ] || { echo "manca $LOG/$f: eseguire prima 'test'"; exit 2; }; done
  declare -A FILE PAT KIND
  FILE[M1]=R/04b4_derive_cell_geometry.R; KIND[M1]=R;  PAT[M1]='s/if \(any\(tt\)\) th\[tt\] <- \.dcg_territory_axis\(terr\[tt\]\)/if (any(tt)) th[tt] <- .dcg_territory_axis(terr[tt]) + pi \/ 2/'
  FILE[M2]=tools/S1.4e_bench_lib.R;       KIND[M2]=R;  PAT[M2]='s/ct\$nucleus_aspect <- s14e_aspect\(A\)/ct\$nucleus_aspect <- 1/'
  FILE[M3]=tools/S1.4e_bench_lib.R;       KIND[M3]=R;  PAT[M3]='s/random_seed = L\$metadata\$stage_seeds\[\["size"\]\]/random_seed = L\$metadata\$stage_seeds[["size"]] + 1L/'
  FILE[M4]=tools/S1.4e_measure.py;        KIND[M4]=py; PAT[M4]='s/np\.cos\(theta\) \* ux \+ np\.sin\(theta\) \* uy/np.sin(theta) * ux + np.cos(theta) * uy/'
  FILE[M5]=tools/S1.4e_measure.py;        KIND[M5]=py; PAT[M5]='s/^CHEB = 2 /CHEB = 3 /m'
  FILE[M6]=tools/S1.4e_bench_lib.R;       KIND[M6]=R;  PAT[M6]='s/keep_co <- !inwin\(co\$x, co\$y\) \| inex\(co\$x, co\$y\)/keep_co <- rep(TRUE, length(co\$x))/'
  MS="M1 M2 M3 M4 M5 M6"; declare -A TD; n=0
  for m in $MS; do
    T=$(mktemp -d); TD[$m]=$T; mkdir -p "$T/R" "$T/tools"; cp R/04b*.R "$T/R/"; cp tools/S1.4e_bench_lib.R tools/S1.4e_measure.py tools/S1.4c_raster.py tools/R2c_candidates.py "$T/tools/"
    perl -0pi -e "${PAT[$m]}" "$T/${FILE[$m]}"
    if cmp -s "$T/${FILE[$m]}" "${FILE[$m]}"; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; touch "$LOG/skip_$m"; continue; fi
    rm -f "$LOG/skip_$m"
    if [ "${KIND[$m]}" = R ]; then ( S14E_RSRC="$T" Rscript --vanilla R/testing/test_S1.4e.R > "$LOG/mut_$m.log" 2>&1 ) &
    else ( S14E_PYSRC="$T" .venv/bin/python tools/S1.4e_pairs_test.py > "$LOG/mut_$m.log" 2>&1 ) & fi
  done
  wait
  for m in $MS; do
    [ -f "$LOG/skip_$m" ] && continue
    base="$LOG/test_S1.4e.log"; [ "${KIND[$m]}" = py ] && base="$LOG/test_S1.4e_py.log"
    newf=$(grep '^FAIL' "$LOG/mut_$m.log" | cut -d' ' -f2- | sort -u | comm -23 - <(grep '^FAIL' "$base" | cut -d' ' -f2- | sort -u) | cut -c1-60 | tr '\n' ';')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/mut_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
    rm -rf "${TD[$m]}"
  done
  echo "mutanti rilevati $n/6"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.4e.sh test|mutants"; exit 2; fi
exit $fail
