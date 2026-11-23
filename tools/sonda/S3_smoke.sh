#!/usr/bin/env bash
# tools/sonda/S3_smoke.sh - check smoke (S3, docs/sonda/REGOLA_S3.md): circa 10 minuti, configurazione ridotta, dati di sviluppo.
# 1. test rimasti: tools/sonda/S1_test_cc.py (metro della registrazione di S1);
# 2. catena da capo a fondo su un ingresso ridotto: quadrato di 250 um dentro l'area H2 tumour 1 (ritaglio H&E, StarDist, esportazione,
#    bin2cell, RCTD, misure) e un B sullo stesso quadrato (parametri tarati della corsa completa); uscite in una cartella da buttare;
# 3. C-3: misure identiche (1e-9) a results/sonda/S3/smoke_ref.json (con --init lo scrive);
# 4. metro sulla tabella completa (results/sonda/S3/area_values.csv, h.csv): D, controlli C-1, C-2, C-4, C-5, identita' con
#    results/sonda/S3/D.json; scrive results/misura.json (A1 = D(B)). Esce con 0 solo se tutto passa ("SMOKE PASS").
# Uso: bash tools/sonda/S3_smoke.sh [--init]
# S3b: il metro legge S3_TABLE, S3_H e confronta con S3_DREF (predefiniti: tabelle di S3); vedi tools/sonda/S3b_smoke.sh.
set -uo pipefail
cd ~/2025.geo_spatialtrans
S=/mnt/micron/geo_spatialtrans/S3_sonda/smoke; rm -rf $S; mkdir -p $S/jobs
export S3_ROOT=$S S3_BPAR=/mnt/micron/geo_spatialtrans/S3_sonda/B OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
fail=0; t0=$(date +%s)
OMP_NUM_THREADS=1 .venv/bin/python -B tools/sonda/S1_test_cc.py > $S/S1test.log 2>&1
grep -q 'RESULT PASS' $S/S1test.log && echo "ok   S1_test_cc" || { echo "FAIL S1_test_cc"; fail=1; }
.venv/bin/python -B tools/sonda/S3_export.py crop H2 smoke > $S/chain.log 2>&1 \
 && TF_NUM_INTEROP_THREADS=2 TF_NUM_INTRAOP_THREADS=2 bash tools/py_stardist.sh tools/sonda/S3_stardist.py SMOKE_H2_tumour_1 >> $S/chain.log 2>&1 \
 && .venv/bin/python -B tools/sonda/S3_export.py export H2 smoke >> $S/chain.log 2>&1 \
 && .venv/bin/python -B tools/sonda/S3_sim_B.py 1 0 0 smoke >> $S/chain.log 2>&1 || { echo "FAIL export/sim"; fail=1; }
printf 'SMOKE_V\tSMOKE_H2_tumour_1\t%s\ttumour\treal\tH2\t-\nSMOKE_B\tSMOKE_H2_tumour_1_k0\t%s\ttumour\tsimr\tB\t0\n' $S/runs/SMOKE_V/SMOKE_H2_tumour_1 $S/runs/SMOKE_B/SMOKE_H2_tumour_1_k0 > $S/jobs/smoke.tsv
bash tools/sonda/S3_chain.sh $S/jobs/smoke.tsv 2 >> $S/chain.log 2>&1
.venv/bin/python -B tools/sonda/S3_smoke_check.py $S ${1:-} || fail=1
.venv/bin/python -B tools/sonda/S3_D.py --table ${S3_TABLE:-results/sonda/S3/area_values.csv} --h ${S3_H:-results/sonda/S3/h.csv} --out $S/D.json --per $S/per_measure.csv --misura results/misura.json --provv ${S3_PROVV:-results/sonda/S3/misura_provvisoria.json} > $S/D.log 2>&1 || { echo "FAIL S3_D"; fail=1; }
.venv/bin/python -B tools/sonda/S3_smoke_check.py $S meter || fail=1
echo "elapsed $(( $(date +%s) - t0 )) s"
[ $fail -eq 0 ] && echo "SMOKE PASS" || echo "SMOKE FAIL"
exit $fail
