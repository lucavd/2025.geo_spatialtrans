#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/mnt/micron/geo_spatialtrans/G1_review/rev_registration
PY=/home/user/2025.geo_spatialtrans/.venv/bin/python
cd $OUT
: > tasks2.txt
for p in pos:3.1:0 pos:0:4.3 pos:-3.6:3.45; do for r in R1 R2 R3; do echo "$r $p 200" >> tasks2.txt; done; done
( $PY g1_he.py > logs/he.log 2>&1 ) &
cat tasks2.txt | xargs -P 3 -L 1 bash -c 'f=$(echo "$0_$1_$2" | tr ":" "_"); '$PY' g1_tiles.py $0 $1 $2 > logs/$f.log 2>&1'
wait
bash md5.sh md5_after.txt
diff md5_before.txt md5_after.txt && echo MD5_SAME > md5_compare.txt || echo MD5_DIFF > md5_compare.txt
sensors | grep -i tctl > tctl_end2.txt || true
