#!/usr/bin/env bash
# tools/S1.5_review_mutants.sh - S1.5, controprova dopo la revisione: base + mutanti del revisore R01-R19 + mutanti dell'autore M1-M9
# (copia di review_code/run_mutants.sh con R19 e la cartella postreview). Uso dalla root del repo: bash tools/S1.5_review_mutants.sh
RC=/mnt/micron/geo_spatialtrans/S1.5/postreview; mkdir -p $RC/logs
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
.venv/bin/python tools/S1.5_review_mutants.py || exit 2
run1() { t=$1; m=$2; tree=$3; cd $RC/tree/$tree && S15_MUTANT=$m /usr/bin/time -f "%e s %M kB" Rscript --vanilla R/testing/test_S1.5.R > $RC/logs/$t.log 2>&1; echo "exit $?" >> $RC/logs/$t.log; }
export -f run1; export RC
{ echo "base none base"; for k in $(seq -w 1 19); do echo "R$k none R$k"; done; for k in 1 2 3 4 5 6 7 8 9; do echo "M$k M$k base"; done; } | xargs -P 8 -L 1 bash -c 'run1 $0 $1 $2'
.venv/bin/python tools/S1.5_review_mutants_analyse.py
