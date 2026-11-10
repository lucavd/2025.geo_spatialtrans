#!/usr/bin/env bash
# tools/S2.6_md5.sh - S2.6, G-4: md5 di ogni file letto dallo step, prima e dopo; git diff di R/. Uso: tools/S2.6_md5.sh before|after [outdir]
set -euo pipefail
W=${1:?before o after}; OUT=${2:-results/S2.6}; mkdir -p "$OUT"
M=/mnt/micron/geo_spatialtrans; D=$M/data_real
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 spatial/tissue_positions.parquet; do
    p="$D/$ds/binned_outputs_x/binned_outputs/square_008um/$f"; [ -e "$p" ] && echo "$p"
  done
done
ls $M/S2.2/profiles/*
ls R/04b*.R R/05b*.R
} | sort | xargs md5sum > "$OUT/S2.6_md5_${W}.txt"
git diff --stat -- R/ ':!R/testing/test_S2.6.R' > "$OUT/S2.6_gitdiff_R_${W}.txt"
wc -l < "$OUT/S2.6_md5_${W}.txt"
if [ "$W" = after ]; then
  diff "$OUT/S2.6_md5_before.txt" "$OUT/S2.6_md5_after.txt" > "$OUT/S2.6_md5_diff.txt" && echo "G-4 md5 identici" || { echo "G-4 md5 DIFFERENZE"; exit 1; }
  [ ! -s "$OUT/S2.6_gitdiff_R_after.txt" ] && echo "G-4 R/ invariato" || { echo "G-4 R/ MODIFICATO"; exit 1; }
fi
