#!/usr/bin/env bash
# R/testing/run_S1.5d.sh - S1.5d: cancello G-1 (test) e G-2 (mutanti). Uso (root del repo):
#   bash R/testing/run_S1.5d.sh test      -> test_S1.5d.R, poi la batteria di S1.5c (che include S1.5b, S1.5, S1.4d, S1.1-S1.4, OKF)
#   bash R/testing/run_S1.5d.sh mutants   -> M1..M11 come copie modificate del codice in una cartella temporanea: ognuno deve far fallire
#                                            almeno un'asserzione di test_S1.5d.R che passa senza mutante
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LOG=${S15D_LOG:-/mnt/micron/geo_spatialtrans/S1.5d/tests}; mkdir -p "$LOG"
mode=${1:-test}; fail=0
if [ "$mode" = test ]; then
  Rscript --vanilla R/testing/test_S1.5d.R > "$LOG/test_S1.5d.log" 2>&1 || fail=1
  grep -E '^(PASS|FAIL)' "$LOG/test_S1.5d.log"; tail -1 "$LOG/test_S1.5d.log"
  grep -q '^FAIL' "$LOG/test_S1.5d.log" && fail=1
  grep -qE '^RESULT ([0-9]+)/\1 PASS' "$LOG/test_S1.5d.log" || fail=1
  S15C_LOG="$LOG/S1.5c" bash R/testing/run_S1.5c.sh test > "$LOG/run_S1.5c.log" 2>&1 || fail=1
  grep -E 'RESULT|GATE' "$LOG/run_S1.5c.log" | tail -6
  grep -q 'GATE G-1 PASS' "$LOG/run_S1.5c.log" || { echo "FAIL batteria di S1.5c"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-1 PASS" || echo "GATE G-1 FAIL"
elif [ "$mode" = mutants ]; then
  base="$LOG/test_S1.5d.log"; [ -f "$base" ] || { echo "manca $base: eseguire prima 'test'"; exit 2; }
  grep '^FAIL' "$base" | cut -d' ' -f2 | sort -u > "$LOG/base_fail.txt"
  declare -A FILE PAT
  FILE[M1]=R/04b2_seed_centroids.R;  PAT[M1]='s/      if \(!is\.null\(het\)\) pr\$w <- pr\$w \* \.sc_het_weight\(het, fv, fw, cx, cy\)\n//'
  FILE[M2]=R/04b2_seed_centroids.R;  PAT[M2]='s/stats::qnorm\(1 - h\$void_frac\)/stats::qnorm(h\$void_frac)/'
  FILE[M3]=R/04b2_seed_centroids.R;  PAT[M3]='s/a <- \.sc_solve_a\(eta_g\[k_in, , drop = FALSE\], f, rho, h = hg\[k_in\]\)/a <- .sc_solve_a(eta_g[k_in, , drop = FALSE], f, rho)/'
  FILE[M4]=R/04b2_seed_centroids.R;  PAT[M4]='s/fw <- \.sc_rff_draw\(het\$dens_ell_um, M\)/fw <- .sc_rff_draw(1 \/ het\$dens_ell_um, M)/'
  FILE[M5]=R/04b2_seed_centroids.R;  PAT[M5]='s/if \(!is\.null\(fw\)\) ell_g <- min\(ell_g, het\$dens_ell_um\)/if (FALSE) ell_g <- min(ell_g, het\$dens_ell_um)/'
  FILE[M6]=tools/S1.5d_bench_lib.R;  PAT[M6]='s/fold = c\("loo", "pooled"\)/fold = c("pooled", "loo")/g'
  FILE[M7]=R/04b2_seed_centroids.R;  PAT[M7]='s/if \(quotas && quota\[tj\] == 0L\)/if (quota[tj] == 0L)/'
  FILE[M8]=R/04b_cell_placement.R;   PAT[M8]='s/labels = labels, heterogeneity = htab, quotas = quotas, verbose = FALSE\)/labels = labels, heterogeneity = NULL, quotas = quotas, verbose = FALSE)/'
  FILE[M9]=R/04b4_derive_cell_geometry.R; PAT[M9]='s/key\[use\], log\(a_t\[use\]\), stats::pnorm/key[use], log(a_t[use]) - log(a_ref[use]), stats::pnorm/'
  FILE[M10]=R/04b0_tissue_presets.R; PAT[M10]='s/paste0\(preset, ":", out\$cell_types\$size_group\)/paste0(preset, ":", rev(out\$cell_types\$size_group))/'
  FILE[M11]=R/04b4_derive_cell_geometry.R; PAT[M11]='s/use <- !is\.na\(key\) & !abs_mode/use <- rep(FALSE, length(key))/'
  n=0
  for m in M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11; do
    T=$(mktemp -d); mkdir -p "$T/R" "$T/tools"; cp R/04b*.R "$T/R/"; cp tools/S1.5d_bench_lib.R "$T/tools/"
    perl -0pi -e "${PAT[$m]}" "$T/${FILE[$m]}"
    if cmp -s "$T/${FILE[$m]}" "${FILE[$m]}"; then echo "FAIL mutante $m: sostituzione non applicata"; fail=1; rm -rf "$T"; continue; fi
    S15D_RSRC="$T" Rscript --vanilla R/testing/test_S1.5d.R > "$LOG/test_S1.5d_$m.log" 2>&1
    newf=$(grep '^FAIL' "$LOG/test_S1.5d_$m.log" | cut -d' ' -f2 | sort -u | comm -23 - "$LOG/base_fail.txt" | tr '\n' ' ')
    if [ -z "$newf" ] && ! grep -q '^RESULT' "$LOG/test_S1.5d_$m.log"; then newf="(arresto dei test)"; fi
    if [ -n "$newf" ]; then echo "PASS mutante $m rilevato da: $newf"; n=$((n+1)); else echo "FAIL mutante $m NON rilevato"; fail=1; fi
    rm -rf "$T"
  done
  echo "mutanti rilevati $n/11 (M1-M8 dell'eterogeneita', M9-M11 dei residui)"; [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
else echo "uso: run_S1.5d.sh test|mutants"; exit 2; fi
exit $fail
