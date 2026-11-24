#!/usr/bin/env bash
# tools/sonda/S6_run.sh - S6 (docs/sonda/REGOLA_S6.md): md5 prima, simulazione dei bracci con i tipi non copiati, cancelli della simulazione,
# catena (Yalt, YH1, poi Yshuf), tabella, md5 dopo, D. Rilanciabile: ogni passo salta cio' che esiste. Processi della catena da $S6_ROOT/WORKERS.
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S6_ROOT=/mnt/micron/geo_spatialtrans/S6_sonda; export S3_ROOT=$S6_ROOT S3_REF=/mnt/micron/geo_spatialtrans/S3_sonda/ref
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $S6_ROOT/jobs results/sonda/S6; [ -f $S6_ROOT/WORKERS ] || echo 10 > $S6_ROOT/WORKERS
[ -f results/sonda/S6/md5_before.txt ] || bash tools/sonda/S6_md5.sh results/sonda/S6/md5_before.txt
echo "sim start $(date)"; .venv/bin/python -B tools/sonda/S6_sim_types.py 10 Yself,Yalt,YH1,Yshuf all > $S6_ROOT/sim.log 2>&1 || { echo "SIM FAIL"; exit 1; }
.venv/bin/python -B tools/sonda/S6_gates_sim.py || { echo "GATES SIM FAIL"; exit 1; }
: > $S6_ROOT/jobs/chain.tsv
for a in Yalt YH1 Yshuf; do for t in tumour stroma; do for k in 1 2 3 4 5; do u=H2_${t}_${k}
  printf '%s\t%s\t%s\t%s\tsimr\t%s\t0\n' $a $u $S6_ROOT/runs/$a/${u}_k0 $t $a >> $S6_ROOT/jobs/chain.tsv; done; done; done
echo "chain start $(date)"; bash tools/sonda/S3_chain.sh $S6_ROOT/jobs/chain.tsv; echo "chain end $(date)"
.venv/bin/python -B tools/sonda/S3_table.py $S6_ROOT/measures results/sonda/S6/area_values.csv
bash tools/sonda/S6_md5.sh results/sonda/S6/md5_after.txt
cmp -s results/sonda/S6/md5_before.txt results/sonda/S6/md5_after.txt && echo "G-MD5 PASS" || echo "G-MD5 FAIL"
.venv/bin/python -B tools/sonda/S6_D.py
