#!/usr/bin/env bash
# tools/S2.5_md5.sh - S2.5, G-4 (a): md5 di ogni file letto dallo step, prima e dopo. Uso: tools/S2.5_md5.sh before|after [outdir]
set -euo pipefail
W=${1:?before o after}; OUT=${2:-results/S2.5}; mkdir -p "$OUT"
M=/mnt/micron/geo_spatialtrans; D=$M/data_real
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for s in square_002um square_008um; do
    for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 spatial/tissue_positions.parquet; do [ -e "$D/$ds/binned_outputs_x/binned_outputs/$s/$f" ] && echo "$D/$ds/binned_outputs_x/binned_outputs/$s/$f"; done
  done
done
ls results/R2/R2_rois_checked.csv results/R4/R4_roi.csv results/R6/R6_dexp.csv results/R7/coverage/R7_kappa_final_v4.csv results/R7/coverage/R7_cv_factor.csv results/R7/Mc/R7_refs_adopted.csv
ls $M/S1.6/layers/*_U6.rds $M/S2.3/panel/hd/*.rds $M/S2.3/layers/*.tsv $M/S2.3/e2/real/*.json $M/S2.2/profiles/*
ls $M/S2.0/run1_wfix/params_flat.json
ls $M/S2.4/panel/main/A2_r1_k0.bin $M/S2.4/panel/main/A6_r1_k0.bin
} | sort | xargs md5sum > "$OUT/S2.5_md5_${W}.txt"
wc -l < "$OUT/S2.5_md5_${W}.txt"
if [ "$W" = after ]; then diff "$OUT/S2.5_md5_before.txt" "$OUT/S2.5_md5_after.txt" > "$OUT/S2.5_md5_diff.txt" && echo "G-4 (a) identici" || { echo "G-4 (a) DIFFERENZE"; exit 1; }; fi
