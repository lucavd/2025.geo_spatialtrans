#!/usr/bin/env bash
# tools/sonda/S4_md5.sh - S4: md5 dei file letti (V, parametri di B, tabelle e D di S3b/S3c, riferimento). Uso: bash tools/sonda/S4_md5.sh <out>
cd ~/2025.geo_spatialtrans; S=/mnt/micron/geo_spatialtrans/S3_sonda
{ for d in $S/runs/V/H2_*; do md5sum $d/bins.parquet $d/X.parquet $d/genes.tsv $d/labels_transferred.parquet $d/exp_labels_expanded.npy $d/cells.parquet; done
  md5sum $S/B/params_* results/sonda/S3b/area_values.csv results/sonda/S3b/h.csv results/sonda/S3c/D.json results/sonda/S3c/per_measure.csv
  find $S/ref -type f | sort | xargs md5sum; } | sed "s#$S#S3#" > $1
wc -l < $1
