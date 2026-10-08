#!/usr/bin/env bash
# tools/R3_mutants.sh — C-R3.7: ogni mutante deve far fallire almeno un check C sintetico.
cd "$(dirname "$0")/.." || exit 1
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
echo "mutant,n_fail,detected" > results/R3/R3_mutants.csv
for M in M1 M2 M3 M4 M5; do
  R3_MUTANT=$M Rscript --vanilla R/testing/test_R3.R synthetic > results/R3/R3_mutant_$M.log 2>&1
  nf=$(grep -c '^FAIL' results/R3/R3_mutant_$M.log || true)
  det=$([ "${nf:-0}" -gt 0 ] && echo TRUE || echo FALSE)
  echo "$M,${nf:-0},$det" >> results/R3/R3_mutants.csv
  echo "$M FAIL=${nf:-0} rilevato=$det"
done
