#!/usr/bin/env bash
# tools/R5_md5.sh before|after — G-2: md5 di ogni file letto da R5 (dati 10x, immagini, tabelle di R2, R3b, R4, R6). Pre-registrazione: docs/steps/R5.md
set -euo pipefail
cd "$(dirname "$0")/.."
M=/mnt/micron/geo_spatialtrans; D=$M/data_real
OUT=results/R5/R5_md5_$1.txt
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for s in square_002um square_008um; do
    for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 spatial/tissue_positions.parquet spatial/scalefactors_json.json; do
      echo "$D/$ds/binned_outputs_x/binned_outputs/$s/$f"; done; done
  ls $D/$ds/${ds}_tissue_image.* ; done
ls $M/R2/masks/*_cellpose_rgb_native.npz $M/R3b/masks/*_consensus_native.npz $M/R3b/nuclei/*_consensus_native.parquet results/R2/nuclei/*_cellpose_rgb_native.parquet
ls results/R2/R2_rois_checked.csv results/R6/R6_dexp.csv $M/R4/run1/posthoc/*_cells.parquet results/R4/R4_PH1_roi.csv results/R4/R4_PH4_parity.csv
} | sort -u | xargs md5sum > "$OUT"
echo "$(wc -l < "$OUT") files -> $OUT"
