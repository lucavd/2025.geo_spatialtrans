#!/usr/bin/env bash
# Revisione avversariale S1.1 — 07_perf.sh: rilancio di C11 sotto /usr/bin/time -v (sola lettura del repo)
set -u
REV=/mnt/micron/geo_spatialtrans/S1.1/review
cd /home/user/2025.geo_spatialtrans
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
echo "run,id,perf_line,elapsed_wall,max_rss_kb" > $REV/07_perf.csv
for run in 1 2 3; do
 for id in P1_4000_c2_labels P2_4000_c4_quant4; do
  /usr/bin/time -v Rscript --vanilla tools/S1.1_perf.R "$id" > $REV/07_perf_${id}_r${run}.out 2> $REV/07_perf_${id}_r${run}.time
  line=$(grep '^PERF,' $REV/07_perf_${id}_r${run}.out | sed 's/^PERF,//' | tr ',' ';')
  wall=$(grep 'Elapsed (wall clock)' $REV/07_perf_${id}_r${run}.time | awk '{print $NF}')
  kb=$(grep 'Maximum resident set size' $REV/07_perf_${id}_r${run}.time | awk '{print $NF}')
  echo "$run,$id,$line,$wall,$kb" >> $REV/07_perf.csv
 done
done
cat $REV/07_perf.csv
git status --porcelain | head

