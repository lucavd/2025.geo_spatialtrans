#!/usr/bin/env bash
# S1.2 — rifa' tutti i numeri del report di seed_centroids() (dalla root del repo)
#   bash R/testing/run_S1.2.sh            (~1 h su 40 core; le simulazioni di check B saltano le unita' gia' fatte)
set -euo pipefail
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
LOG=/mnt/micron/geo_spatialtrans/S1.2/logs; mkdir -p "$LOG"
Rscript --vanilla R/testing/test_S1.2.R | tee "$LOG/test.log"
/usr/bin/time -v Rscript --vanilla tools/S1.2_perf.R > "$LOG/perf.out" 2> "$LOG/perf.time"
line=$(grep '^PERF,' "$LOG/perf.out" | sed 's/^PERF,//')
kb=$(grep 'Maximum resident set size' "$LOG/perf.time" | awk '{print $NF}')
echo "n_target,n_cells,n_attempts,elapsed_s,max_rss_mb" > results/S1.2/S1.2_perf.csv
echo "$line,$((kb / 1024))" >> results/S1.2/S1.2_perf.csv
cat results/S1.2/S1.2_perf.csv
Rscript --vanilla tools/S1.2_checkB.R all | tee "$LOG/checkB.log"
Rscript --vanilla tools/S1.2_analysis.R | tee "$LOG/analysis.log"
