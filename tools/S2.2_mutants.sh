#!/usr/bin/env bash
# S2.2, cancello G-1: mutanti di R/05b2_simulate_cell_profiles.R; ognuno deve far fallire almeno un test di R/testing/test_S2.2.R.
# Uso (root del repo): tools/S2.2_mutants.sh  -> results/S2.2/S2.2_mutants.txt
set -uo pipefail
cd "$(dirname "$0")/.."
F=R/05b2_simulate_cell_profiles.R; OUT=results/S2.2/S2.2_mutants.txt; : > $OUT
declare -A MUT
MUT[M1_depth_ignored]='s/round(depth_factor \* stats::approx/round(1 * stats::approx/'
MUT[M2_subtype_shifted]='s/prob = sub\$share\[r\])\]/prob = rev(sub$share[r]))]/'
MUT[M3_seed_ignored_in_chunk]='s/^    assign(".Random.seed", streams\[\[k\]\], envir = globalenv())$/    invisible(NULL)/'
MUT[M4_gamma_shape_inverted]='s/shape = shape\[z\], rate = shape\[z\]/shape = 1 \/ shape[z], rate = 1 \/ shape[z]/'
MUT[M5_genes_permuted]='s/j = sums\$gene_index, x/j = rev(seq_len(nrow(genes)))[sums$gene_index], x/'
MUT[M6_totals_misassigned]='s/cell <- rep(seq_along(w), Nc\[w\])/cell <- rep(seq_along(w), rev(Nc[w]))/'
MUT[M7_subtypes_always_on]='s/^  if (!use_subtypes) {   /  if (FALSE) {   /'
MUT[M9_subtype_totals_from_type]='s/profiles\$totals_subtype\[\[sub\$subtype\[s_\]\]\]/profiles$totals[[sub$cell_type[s_]]]/'
MUT[M8_rng_not_restored]='s/^  on.exit({ do.call(RNGkind/  if (FALSE) on.exit({ do.call(RNGkind/'
for m in $(printf "%s\n" "${!MUT[@]}" | sort); do
  T=$(mktemp -d); mkdir -p $T/R; cp R/04b*.R R/05b1*.R $T/R/; sed -e "${MUT[$m]}" $F > $T/R/05b2_simulate_cell_profiles.R
  if cmp -s $F $T/R/05b2_simulate_cell_profiles.R; then echo "$m NON APPLICATO" | tee -a $OUT; rm -rf $T; continue; fi
  S22_RSRC=$T Rscript --vanilla R/testing/test_S2.2.R > $T/log.txt 2>&1
  nf=$(grep -c '^FAIL' $T/log.txt || true); er=$(grep -c 'Error\|Errore' $T/log.txt || true)
  if [ "${nf:-0}" -gt 0 ] || [ "${er:-0}" -gt 0 ]; then echo "$m RILEVATO ($nf FAIL, $er errori: $(grep '^FAIL' $T/log.txt | cut -d' ' -f2 | tr '\n' ' '))" | tee -a $OUT
  else echo "$m NON RILEVATO" | tee -a $OUT; fi
  rm -rf $T
done
