#!/usr/bin/env bash
# R6 gate G-2: md5 of every real-data file read by the step. Usage: tools/R6_md5.sh before|after
set -euo pipefail
cd ~/2025.geo_spatialtrans
TAG=${1:?before|after}
D=/mnt/micron/geo_spatialtrans/data_real
OUT=results/R6/R6_md5_${TAG}.txt
: > "$OUT"
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  md5sum $D/$ds/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5 >> "$OUT"
  md5sum $D/$ds/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet >> "$OUT"
  [ -f $D/$ds/${ds}_barcode_mappings.parquet ] && md5sum $D/$ds/${ds}_barcode_mappings.parquet >> "$OUT"
  S=$D/$ds/segmented_outputs_x/segmented_outputs
  [ -f $S/cell_types/10x_Cloud/cell_types.csv ] && md5sum $S/cell_types/10x_Cloud/cell_types.csv >> "$OUT"
  [ -f $S/nucleus_segmentations.geojson ] && md5sum $S/nucleus_segmentations.geojson >> "$OUT"
done
md5sum results/R2/R2_rois_checked.csv results/S1.5/S1.5_real_refs_roi.csv >> "$OUT"
for a in A1 A2 A3 A5; do for r in r1 r2 r3 r4 r5; do
  md5sum results/R2/nuclei/${a}_${r}_cellpose_rgb_native.parquet /mnt/micron/geo_spatialtrans/R2/masks/${a}_${r}_cellpose_rgb_native.npz >> "$OUT"; done; done
for f in /mnt/micron/geo_spatialtrans/R3b/nuclei/A4_*_consensus_native.parquet /mnt/micron/geo_spatialtrans/R3b/nuclei/A6_*_consensus_native.parquet \
         /mnt/micron/geo_spatialtrans/R3b/masks/A4_*_consensus_native.npz /mnt/micron/geo_spatialtrans/R3b/masks/A6_*_consensus_native.npz; do md5sum $f >> "$OUT"; done
for f in /mnt/micron/geo_spatialtrans/R6/PanglaoDB*; do [ -f "$f" ] && md5sum "$f" >> "$OUT"; done
wc -l "$OUT"
