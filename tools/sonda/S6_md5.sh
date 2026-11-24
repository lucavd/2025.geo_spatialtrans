#!/usr/bin/env bash
# tools/sonda/S6_md5.sh - S6: md5 dei file letti (V e T omologhe di S3b, Y di S4, parametri di B, tabelle e D di S3b/S3c/S4, mappa, riferimento).
# Uso: bash tools/sonda/S6_md5.sh <out>
cd ~/2025.geo_spatialtrans; S=/mnt/micron/geo_spatialtrans/S3_sonda; T=/mnt/micron/geo_spatialtrans/S3b_sonda/runs/T; Y=/mnt/micron/geo_spatialtrans/S4_sonda/runs/Y
{ for d in $S/runs/V/H2_*; do md5sum $d/bins.parquet $d/X.parquet $d/genes.tsv $d/labels_transferred.parquet $d/export.json; done
  for d in $T/H1_*; do md5sum $d/bins.parquet $d/labels_transferred.parquet $d/export.json; done
  for d in $Y/H2_*; do md5sum $d/bins.parquet $d/X.parquet $d/cells.parquet; done
  md5sum $S/B/params_* results/sonda/S3b/area_values.csv results/sonda/S3b/h.csv results/sonda/S3c/D.json results/sonda/S4/area_values.csv results/sonda/S4/D.json
  md5sum results/sonda/S3b/overlap.json results/sonda/S3b/areas_H1_hom.csv
  find $S/ref -type f | sort | xargs md5sum; } | sed -e "s#$S#S3#" -e "s#$T#S3bT#" -e "s#$Y#S4Y#" > $1
wc -l < $1
