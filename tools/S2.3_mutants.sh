#!/usr/bin/env bash
# S2.3, cancello G-1: mutanti di R/05b3_place_and_move_transcripts.R; ognuno deve far fallire almeno un test di R/testing/test_S2.3.R.
# Uso (root del repo): tools/S2.3_mutants.sh  -> results/S2.3/S2.3_mutants.txt
set -uo pipefail
cd "$(dirname "$0")/.."
F=R/05b3_place_and_move_transcripts.R; OUT=results/S2.3/S2.3_mutants.txt; mkdir -p results/S2.3; : > $OUT
declare -A MUT
MUT[M01_nucleus_density_inverted]='s/pn <- ifelse(A_n > 0, r \* A_n \/ (r \* A_n + A_c), 0)/pn <- ifelse(A_n > 0, A_n \/ r \/ (A_n \/ r + A_c), 0)/'
MUT[M02_distance_wrong_point]='s/    d <- Dm\[cbind(pa, qa)\]/    d <- Dm[cbind(pmin(nr, pa + 8L), qa)]/'
MUT[M03_tail_exponential_radius]='s/stats::rgamma(nT, shape = 2, scale = params\$ell)/stats::rgamma(nT, shape = 1, scale = params$ell)/'
MUT[M04_field_variance_wrong]='s/shape = 1 \/ s2, scale = s2/shape = 2 \/ s2, scale = s2 \/ 2/g'
MUT[M05_field_ignored_in_acceptance]='s/ua < Fm\[cbind(pa, qa)\] \/ Fmax\[cell\] }/ua < 1 }/'
MUT[M06_lost_dropped]='s/    lost <- tabulate(match(cell\[!inside\], idx), length(idx))/    lost <- integer(length(idx))/'
MUT[M07_sigma_on_radius]='s/sg <- ifelse(wide\[!tail\], params\$sigma_e, params\$sigma_in)/sg <- ifelse(wide[!tail], params$sigma_e, params$sigma_in) \/ sqrt(2)/'
MUT[M08_seed_ignored_in_chunk]='s/^    assign(".Random.seed", streams\[\[k\]\], envir = globalenv())$/    invisible(NULL)/'
MUT[M09_rows_cols_swapped]='s/bin <- ifelse(inside, ii \* gn + jj + 1, NA_real_)/bin <- ifelse(inside \& jj < gnr \& ii < gn, jj * gn + ii + 1, NA_real_)/'
MUT[M10_rng_not_restored]='s/^  on.exit({ do.call(RNGkind/  if (FALSE) on.exit({ do.call(RNGkind/'
MUT[M11_frame_open_ignored]='s/^  if (frame_open) {$/  if (FALSE) {/'
for m in $(printf "%s\n" "${!MUT[@]}" | sort); do
  T=$(mktemp -d); mkdir -p $T/R; cp R/04b*.R R/05b1*.R R/05b2*.R $T/R/; sed -e "${MUT[$m]}" $F > $T/R/05b3_place_and_move_transcripts.R
  if cmp -s $F $T/R/05b3_place_and_move_transcripts.R; then echo "$m NON APPLICATO" | tee -a $OUT; rm -rf $T; continue; fi
  S23_RSRC=$T Rscript --vanilla R/testing/test_S2.3.R > $T/log.txt 2>&1
  nf=$(grep -c '^FAIL' $T/log.txt || true); er=$(grep -c '^Error\|Errore' $T/log.txt || true)
  if [ "${nf:-0}" -gt 0 ] || [ "${er:-0}" -gt 0 ]; then echo "$m RILEVATO ($nf FAIL, $er errori: $(grep '^FAIL' $T/log.txt | cut -d' ' -f2 | tr '\n' ' '))" | tee -a $OUT
  else echo "$m NON RILEVATO" | tee -a $OUT; fi
  rm -rf $T
done
