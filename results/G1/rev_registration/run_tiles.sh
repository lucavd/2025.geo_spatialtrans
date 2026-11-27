#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/mnt/micron/geo_spatialtrans/G1_review/rev_registration
PY=/home/user/2025.geo_spatialtrans/.venv/bin/python
cd $OUT && mkdir -p tiles logs
: > tasks.txt
for p in main pos:3:0 pos:0:4 pos:-3.5:3.5 warp neg_swap neg_shuffle halfA halfB; do
  for r in R1 R2 R3; do echo "$r $p 200" >> tasks.txt; done; done
for p in halfA halfB; do for r in R1 R2 R3; do echo "$r $p 100" >> tasks.txt; done; done
cat tasks.txt | xargs -P 4 -L 1 bash -c 'f=$(echo "$0_$1_$2" | tr ":" "_"); '$PY' g1_tiles.py $0 $1 $2 > logs/$f.log 2>&1'
cat logs/*.log > tiles_all.log
