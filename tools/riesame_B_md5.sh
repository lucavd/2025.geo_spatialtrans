#!/usr/bin/env bash
# tools/riesame_B_md5.sh - riesame, Parte B (braccio Catena): md5 di ogni file letto, prima e dopo. Uso: bash tools/riesame_B_md5.sh before|after
set -euo pipefail
cd ~/2025.geo_spatialtrans
M=/mnt/micron/geo_spatialtrans; OUT=$M/riesame/B/md5; mkdir -p $OUT
{
  for A in A2 A3; do for k in 1 2 3 4 5; do
    for f in bins.parquet X.parquet cells.parquet genes.tsv b2c_label.json rctd.json rctd_blocks.csv; do echo $M/F3.1/main/real/${A}_r$k/$f; done
    for s in 0 1; do for f in bins.parquet X.parquet cells.parquet genes.tsv origin.parquet profiles.parquet background.parquet export.json b2c_label.json rctd.json rctd_blocks.csv; do echo $M/F3.1/main/sim/${A}_r${k}_k$s/$f; done; done
    echo $M/S1.6/layers/${A}_r${k}_U6.rds; echo $M/S2.3/layers/${A}_r${k}_regions.tsv
  done; done
  cut -f3 $M/riesame/B/jobs/panel_T.tsv | while read d; do for f in bins.parquet X.parquet cells.parquet genes.tsv b2c_label.json rctd.json rctd_blocks.csv; do echo $d/$f; done; done
  for A in A2 A3; do echo $M/F3.1/main/real/../../../F3.3/runs/${A}_r1_cg/rctd_blocks.csv; done | xargs -n1 readlink -f
  find $M/F3.1/main/ref -type f; find $M/S2.2/profiles -type f -name 'A[23]_*'
  echo $M/riesame/xenium/panel_ids_in_common.txt
  ls $M/F3.3/measures/*.json
  for f in results/F3.3/F3.3_area_values.csv results/F3.3/F3.3_crc_common_genes.txt results/F3.3/F3.3_prereg_numbers.csv results/F3.1/F3.1_intervals.csv \
           results/R7/coverage/R7_kappa_final_v4.csv results/R7/coverage/R7_cv_factor.csv results/R2/R2_rois_checked.csv; do echo $PWD/$f; done
} | sort -u > $OUT/files_$1.txt
xargs -a $OUT/files_$1.txt -d '\n' -P 2 -n 50 md5sum | sort -k2 > $OUT/md5_$1.txt
wc -l < $OUT/md5_$1.txt
