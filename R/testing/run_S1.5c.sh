#!/usr/bin/env bash
# R/testing/run_S1.5c.sh - S1.5c: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5c.sh test      -> test_S1.5c.R, poi la batteria di S1.5b (che include S1.5, S1.4d, S1.1-S1.4, OKF)
#   bash R/testing/run_S1.5c.sh mutants   -> M1..M8 come copie modificate del codice in una cartella temporanea: ognuno deve far fallire
#                                            almeno un'asserzione di test_S1.5c.R che passa senza mutante
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15C_LOG:-/mnt/micron/geo_spatialtrans/S1.5c/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  Rscript --vanilla R/testing/test_S1.5c.R > "$LOG/test_S1.5c.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.5c.log"; tail -1 "$LOG/test_S1.5c.log"
  grep -q '^FAIL' "$LOG/test_S1.5c.log" && fail=1
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/test_S1.5c.log" || fail=1
  S15B_LOG="$LOG/S1.5b" bash R/testing/run_S1.5b.sh test > "$LOG/run_S1.5b.log" 2>&1 || fail=1
  grep -E 'RESULT|GATE|test_S1.5:' "$LOG/run_S1.5b.log" | tail -6
  grep -q 'GATE G-1 PASS' "$LOG/run_S1.5b.log" || { echo "FAIL batteria di S1.5b"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5c.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | cut -d' ' -f2 | sort -u > "$LOG/base_fail.txt"
  declare -A FILE PAT
  FILE[M1]=R/04b2_seed_centroids.R; PAT[M1]='s/eta\[, t\] <- eta\[, t\] \+ s\[t\] \* /eta[, t] <- eta[, t] + 0 * /'
  FILE[M2]=R/04b2_seed_centroids.R; PAT[M2]='s/wmax <- 1\.1 \* max\(wg\)/wmax <- mean(wg)/'
  FILE[M3]=R/04b2_seed_centroids.R; PAT[M3]='s/a <- \.sc_solve_a\(eta_g, f, rho\)/a <- numeric(K)/'
  FILE[M4]=R/04b2_seed_centroids.R; PAT[M4]='s/\(di \+ dj\) \* 0\.5/pmax(di, dj)/'
  FILE[M5]=R/04b2_seed_centroids.R; PAT[M5]='s/stats::rnorm\(2L \* M, 0, 1 \/ ell\)/stats::rnorm(2L * M, 0, ell)/'
  FILE[M6]=R/04b2_seed_centroids.R; PAT[M6]='s/keep <- u_thin < pr\$w \/ wmax/keep <- u_thin <= 1/'
  FILE[M7]=R/04b2_seed_centroids.R; PAT[M7]='s/tl <- \.sc_draw_type\(pr\$p\[keep, , drop = FALSE\], u_type\[keep\]\)/tl <- .sc_draw_type(matrix(f, sum(keep), K, byrow = TRUE), u_type[keep])/'
  FILE[M8]=tools/S1.5c_bench_lib.R; PAT[M8]='s/fold = c\("loo", "pooled"\)/fold = c("pooled", "loo")/'
  n=0
  for m in M1 M2 M3 M4 M5 M6 M7 M8; do
    T=$(mktemp -d); mkdir -p "$T/R" "$T/tools"; cp R/04b*.R "$T/R/"; cp tools/S1.5c_bench_lib.R "$T/tools/"
    perl -pi -e "${PAT[$m]}" "$T/${FILE[$m]}"
    if cmp -s "$T/${FILE[$m]}" "${FILE[$m]}"; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; rm -rf "$T"; continue; fi
    S15C_RSRC="$T" Rscript --vanilla R/testing/test_S1.5c.R > "$LOG/test_S1.5c_$m.log" 2>&1
    newf=$(grep '^FAIL' "$LOG/test_S1.5c_$m.log" | cut -d' ' -f2 | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/test_S1.5c_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
    rm -rf "$T"
  done
  echo "mutanti rilevati $n/8"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5c.sh test|mutants"; exit 2; fi
exit $fail
