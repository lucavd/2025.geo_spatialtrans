#!/usr/bin/env bash
# tools/R7_md5.sh before|after — R7, cancello G-4: md5 di tutto ciò che lo step legge.
set -euo pipefail
cd "$(dirname "$0")/.."
tag=${1:?before|after}
out=results/R7/R7_md5_${tag}.txt
{
  ls results/R2/R2_rois_checked.csv results/R2/R2_nuclei_all.parquet results/R2/R2_roi_summary.csv
  ls results/R2/tissue_masks/*
  ls results/R2b/*.csv results/R2c/*.csv
  ls results/R4/R4_roi.csv results/R4/R4_refs.csv results/R4/R4_predictions.csv
  ls /mnt/micron/geo_spatialtrans/R4/run1/roi/*
  ls docs/BIO_REFERENCES.md
  ls results/*/*verdicts*.csv results/*/*predictions*.csv 2>/dev/null
  ls results/R7/inventory/*.csv
} | sort -u | xargs md5sum > "$out"
echo "$out: $(wc -l < "$out") file"
if [ "$tag" = after ] && [ -f results/R7/R7_md5_before.txt ]; then
  if diff -q results/R7/R7_md5_before.txt "$out" >/dev/null; then echo "G-4 PASS: identici"; else echo "G-4 FAIL"; diff results/R7/R7_md5_before.txt "$out" | head; fi
fi
