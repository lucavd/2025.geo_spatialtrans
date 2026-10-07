#!/usr/bin/env bash
# S1.1 — rifà tutti i numeri del report di extract_regions() (dalla root del repo)
#   bash R/testing/run_S1.1.sh
set -euo pipefail
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
LOG=/mnt/micron/geo_spatialtrans/S1.1/logs; mkdir -p "$LOG"
Rscript --vanilla tools/S1.1_inputs.R                        # salta gli input gia' presenti
Rscript --vanilla R/testing/test_S1.1.R | tee "$LOG/test.log"
echo "id,n_px,n_components,n_regions,elapsed_extract_s,rss_after_load_mb,max_rss_process_mb" > results/S1.1/S1.1_perf.csv
for id in P1_4000_c2_labels P2_4000_c4_quant4; do
  /usr/bin/time -v Rscript --vanilla tools/S1.1_perf.R "$id" > "$LOG/perf_$id.out" 2> "$LOG/perf_$id.time"
  line=$(grep '^PERF,' "$LOG/perf_$id.out" | sed 's/^PERF,//')
  kb=$(grep 'Maximum resident set size' "$LOG/perf_$id.time" | awk '{print $NF}')
  echo "$line,$((kb / 1024))" >> results/S1.1/S1.1_perf.csv
done
cat results/S1.1/S1.1_perf.csv
Rscript --vanilla tools/S1.1_figures.R
