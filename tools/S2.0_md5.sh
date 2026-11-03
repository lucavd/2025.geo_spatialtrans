#!/usr/bin/env bash
# tools/S2.0_md5.sh before|after — G-2 di S2.0: md5 di ogni file letto (dati 10x dei 5 dataset, immagini H&E, maschere e nuclei di R2 e R3b,
# tabelle di R2, R4, R5, R6, derivati di R5 run3: sezioni, tessere, PH-5, parametri, lumi) e dei risultati committati di R2, R4, R5, R6.
# Pre-registrazione: docs/steps/S2.0.md (afef06b). Come tools/R5_md5.sh piu' i derivati di R5.
set -euo pipefail
cd "$(dirname "$0")/.."
M=/mnt/micron/geo_spatialtrans; D=$M/data_real; R5=$M/R5/run3
OUT=results/S2.0/S2.0_md5_$1.txt; mkdir -p results/S2.0
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for s in square_002um square_008um; do
    for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 spatial/tissue_positions.parquet spatial/scalefactors_json.json; do
      echo "$D/$ds/binned_outputs_x/binned_outputs/$s/$f"; done; done
  ls $D/$ds/${ds}_tissue_image.* ; done
ls $M/R2/masks/*_cellpose_rgb_native.npz $M/R3b/masks/*_consensus_native.npz $M/R3b/nuclei/*_consensus_native.parquet results/R2/nuclei/*_cellpose_rgb_native.parquet
ls results/R2/R2_rois_checked.csv $M/R4/run1/posthoc/*_cells.parquet
ls $R5/section/*.npz $R5/tiles/*/*.npz $R5/tiles/*/*.json $R5/ph5_s0.7/*.json $R5/params.json $R5/kernels.json
find results/R2 results/R4 results/R5 results/R6 -type f
} | sort -u | xargs md5sum > "$OUT"
echo "$(wc -l < "$OUT") files -> $OUT"
