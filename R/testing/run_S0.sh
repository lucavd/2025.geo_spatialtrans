#!/usr/bin/env bash
# R/testing/run_S0.sh — S0: esegue in sequenza tutti i check della sessione S0 e scrive i log in logs/
# Uso (dalla root del repo): bash R/testing/run_S0.sh
set -u
cd "$(dirname "$0")/../.."
mkdir -p logs results
# Isolamento della libreria di progetto garantito dal codice, non dal fatto che le site-library di sistema siano vuote
# (rilievo del revisore S0): tutti i passi tranne la controprova 2 girano con R_LIBS_SITE che punta a una directory inesistente.
export R_LIBS_SITE=/nonexistent
export R_LIBS_USER=/nonexistent   # anche la libreria utente (~/R/x86_64-pc-linux-gnu-library/4.6 contiene 65 pacchetti)
echo "== $(date -Is) test_S0_library"
Rscript --vanilla R/testing/test_S0_library.R > logs/S0_test_library.log 2>&1; echo "exit=$?" | tee -a logs/S0_test_library.log
echo "== $(date -Is) full_test run1 (seed 42)"
Rscript --vanilla R/testing/full_test.R --tag=S0_run1 > logs/S0_full_test_run1.log 2>&1; echo "exit=$?" | tee -a logs/S0_full_test_run1.log
echo "== $(date -Is) full_test run2 (seed 42)"
Rscript --vanilla R/testing/full_test.R --tag=S0_run2 > logs/S0_full_test_run2.log 2>&1; echo "exit=$?" | tee -a logs/S0_full_test_run2.log
echo "== $(date -Is) full_test seed43 (controprova 1)"
Rscript --vanilla R/testing/full_test.R --seed=43 --tag=S0_seed43 > logs/S0_full_test_seed43.log 2>&1; echo "exit=$?" | tee -a logs/S0_full_test_seed43.log
echo "== $(date -Is) controprova 2: script originale (dev) senza .libPaths() di progetto"
git show dev:R/testing/full_test.R > /tmp/full_test_dev_orig.R
R_LIBS_SITE= R_LIBS_USER= Rscript --vanilla /tmp/full_test_dev_orig.R > logs/S0_controprova2_nolib.log 2>&1; echo "exit=$?" | tee -a logs/S0_controprova2_nolib.log
echo "== $(date -Is) test_S0_baseline"
Rscript --vanilla R/testing/test_S0_baseline.R > logs/S0_test_baseline.log 2>&1; echo "exit=$?" | tee -a logs/S0_test_baseline.log
echo "== $(date -Is) checksums (file tracciato: results/S0_checksums.md5)"
md5sum results/full_simulation_data_S0_run1.rds results/full_simulation_data_S0_run2.rds results/full_simulation_data_S0_seed43.rds > results/S0_checksums.md5
[ -f results/full_simulation_data_S0.rds ] && md5sum results/full_simulation_data_S0.rds >> results/S0_checksums.md5
cat results/S0_checksums.md5
echo "== $(date -Is) done"
