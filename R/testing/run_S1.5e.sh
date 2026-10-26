#!/usr/bin/env bash
# R/testing/run_S1.5e.sh - S1.5e: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5e.sh test      -> test_S1.5e.R (tutte PASS), poi la batteria di S1.5d (che include S1.5c, S1.5b, S1.5, S1.4d, S1.1-S1.4, OKF):
#                                            FAIL ammessi solo T13 di S1.5d (D-6 di S1.5d), T12 di S1.5c (D-10 di S1.5d), E10 ed E10b di S1.5 (DEC-011)
#   bash R/testing/run_S1.5e.sh mutants   -> M1..M9 come copie modificate del codice in cartelle temporanee, in parallelo: ognuno deve far fallire
#                                            almeno un'asserzione di test_S1.5e.R che passa senza mutante
# Dopo le regressioni ripristina i due file che test_S1.1.R e test_S1.2.R riscrivono con tempi nuovi (BL-165, D-13 di S1.5d).
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15E_LOG:-/mnt/micron/geo_spatialtrans/S1.5e/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
fails_of() { grep '^FAIL' "$1" 2>/dev/null | awk '{print $2}' | sort -u | tr '\n' ' ' | sed 's/ $//'; }
if [ "$mode" = test ]; then
  Rscript --vanilla R/testing/test_S1.5e.R > "$LOG/test_S1.5e.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL|RESULT)' "$LOG/test_S1.5e.log"
  grep -q '^FAIL' "$LOG/test_S1.5e.log" && fail=1
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/test_S1.5e.log" || fail=1
  S15D_LOG="$LOG/S1.5d" bash R/testing/run_S1.5d.sh test > "$LOG/run_S1.5d.log" 2>&1
  git checkout -- results/S1.1/S1.1_per_input.csv results/S1.2/S1.2_test_results.csv 2>/dev/null
  D="$LOG/S1.5d"; B="$D/S1.5c/S1.5b"; S="$B/S1.5"
  while read -r f allowed; do
    [ -f "$f" ] || { echo "FAIL regressione: manca $f"; fail=1; continue; }
    grep -qE '^(RESULT|RISULTATO)' "$f" || { echo "FAIL regressione: $f senza RESULT (batteria interrotta)"; fail=1; continue; }
    got=$(fails_of "$f"); [ "$allowed" = "-" ] && allowed=""
    if [ "$got" = "$allowed" ]; then echo "PASS regressione $(basename "$f"): $(grep -E '^(RESULT|RISULTATO)' "$f" | tail -1) FAIL ammessi [${allowed}]"
    else echo "FAIL regressione $(basename "$f"): FAIL trovati [${got}] ammessi [${allowed}]"; fail=1; fi
  done <<LIST
$D/test_S1.5d.log T13
$D/S1.5c/test_S1.5c.log T12
$B/test_S1.5b.log -
$S/test_S1.5.log E10 E10b
$S/test_S1.4d.log -
$S/test_okf.log -
LIST
  for k in 1 2 3 4; do grep -q "test_S1.$k exit 0" "$B/run_S1.5.log" || { echo "FAIL regressione: test_S1.$k non terminato"; fail=1; }; done
  nreg=$(grep -c '^PASS regressione results' "$B/run_S1.5.log"); nbad=$(grep -c '^FAIL regressione' "$B/run_S1.5.log")
  if [ "$nreg" = 4 ] && [ "$nbad" = 0 ]; then echo "PASS regressioni S1.1-S1.4 identiche a HEAD (4/4)"; else echo "FAIL regressioni S1.1-S1.4: $nreg PASS, $nbad FAIL"; fail=1; fi
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5e.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | awk '{print $2}' | sort -u > "$LOG/base_fail.txt"
  declare -A FILE PAT
  FILE[M1]=R/04b2_seed_centroids.R;  PAT[M1]='s/if \(het\$void_mode == "remove" && !is\.null\(fv\)\) \{/if (FALSE && het\$void_mode == "remove" && !is.null(fv)) {/'
  FILE[M2]=R/04b2_seed_centroids.R;  PAT[M2]='s/v <- mean\(!k_in\); void_v\[i\] <- v/v <- het\$void_frac; void_v[i] <- v/'
  FILE[M3]=R/04b2_seed_centroids.R;  PAT[M3]='s/fv <- if \(het\$void_aspect > 1\) \.sc_rff_draw_aniso\(/fv <- if (FALSE) .sc_rff_draw_aniso(/'
  FILE[M4]=R/04b2_seed_centroids.R;  PAT[M4]='s/th <- stats::runif\(1L, 0, pi\)/th <- 0 * stats::runif(1L, 0, pi)/'
  FILE[M5]=R/04b2_seed_centroids.R;  PAT[M5]='s/z\[, 2\] \* sqrt\(k\) \/ ell\)/z[, 2] \/ (ell * sqrt(k)))/'
  FILE[M6]=R/04b2_seed_centroids.R;  PAT[M6]='s/if \(!is\.null\(fw\) && !is\.null\(cap\)\) h <-/if (FALSE) h <-/'
  FILE[M7]=tools/S1.5e_bench_lib.R;  PAT[M7]='s/het_for_e <- function\(A, roi, family, fold = c\("loo", "pooled"\)\)/het_for_e <- function(A, roi, family, fold = c("pooled", "loo"))/'
  FILE[M8]=R/04b_cell_placement.R;   PAT[M8]='s/labels = labels, heterogeneity = htab, quotas = quotas/labels = labels, heterogeneity = if (is.null(htab)) NULL else htab[, intersect(names(htab), c("void_frac", "void_ell_um", "dens_sigma", "dens_ell_um")), drop = FALSE], quotas = quotas/'
  FILE[M9]=R/04b2_seed_centroids.R;  PAT[M9]='s/        quota <- q2\n/        quota <- quota\n/'
  # M10-M13: mutanti del revisore del codice (RC-5), aggiunti dopo la revisione
  FILE[M10]=R/04b2_seed_centroids.R; PAT[M10]='s/cap\$s <- \.sc_cap_scale\(h0, cap\$H\)/cap\$s <- 1/'
  FILE[M11]=R/04b2_seed_centroids.R; PAT[M11]='s/lam_bar <- n \/ \(rdf\$area_um2\[i\] \* mean\(k_in\)\)/lam_bar <- n \/ rdf\$area_um2[i]/'
  FILE[M12]=R/04b2_seed_centroids.R; PAT[M12]='s/lam_cap <- 0\.8 \* 0\.547/lam_cap <- 1.0 * 0.547/'
  FILE[M13]=R/04b2_seed_centroids.R; PAT[M13]='s/out_rdf\$n_target <- plan\$n_region - pl\$n_removed/out_rdf\$n_target <- plan\$n_region/'
  MS="M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13"; declare -A TD
  for m in $MS; do
    T=$(mktemp -d); TD[$m]=$T; mkdir -p "$T/R" "$T/tools"; cp R/04b*.R "$T/R/"; cp tools/S1.5e_bench_lib.R "$T/tools/"
    perl -0pi -e "${PAT[$m]}" "$T/${FILE[$m]}"
    if cmp -s "$T/${FILE[$m]}" "${FILE[$m]}"; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; touch "$LOG/skip_$m"; continue; fi
    rm -f "$LOG/skip_$m"; ( S15E_RSRC="$T" Rscript --vanilla R/testing/test_S1.5e.R > "$LOG/test_S1.5e_$m.log" 2>&1 ) &
  done
  wait
  n=0
  for m in $MS; do
    [ -f "$LOG/skip_$m" ] && continue
    newf=$(grep '^FAIL' "$LOG/test_S1.5e_$m.log" | awk '{print $2}' | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/test_S1.5e_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
    rm -rf "${TD[$m]}"
  done
  echo "mutanti rilevati $n/13 (pre-registrati M1-M9; M10-M13 aggiunti dopo la revisione)"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5e.sh test|mutants"; exit 2; fi
exit $fail
