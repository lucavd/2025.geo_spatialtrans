#!/usr/bin/env bash
# tools/S2.8_md5.sh before|after - S2.8, cancello G-4: impronte md5 dei file letti dallo step (sola lettura).
set -euo pipefail
cd "$(dirname "$0")/.."
MODE=${1:?before|after}; OUT=results/S2.8/S2.8_md5_${MODE}.txt; mkdir -p results/S2.8
M=/mnt/micron/geo_spatialtrans
TAGS="A1_r1 A2_r1 A3_r1 A4_f1 A5_r1 A6_r1"
{
  for t in $TAGS; do echo "$M/S1.6/layers/${t}_U6.rds"; done
  ls $M/S2.2/profiles/A[1-6]_* 2>/dev/null
  echo "$M/S2.7/tasks_P1.json"
  for t in $TAGS; do for k in 0 1 2; do
    ls $M/S2.7/p1/new/${t}_k${k}.md5 $M/S2.7/p1/new/${t}_k${k}.json $M/S2.7/p1/new/${t}_k${k}*.bin 2>/dev/null
  done
    ls $M/S2.7/p1/old/${t}_k2.json $M/S2.7/p1/old/${t}_k2*.bin 2>/dev/null
    for k in 0 1; do ls $M/S2.5/panel/main_r2/${t}_k${k}.json $M/S2.5/panel/main_r2/${t}_k${k}*.bin 2>/dev/null; done
  done
  echo results/R2/R2_rois.csv
} | sort -u | while read -r f; do md5sum "$f"; done > "$OUT"
echo "$(wc -l < "$OUT") file in $OUT"
