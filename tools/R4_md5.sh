#!/usr/bin/env bash
# tools/R4_md5.sh before|after [outfile] — cancello G-2 di R4: md5 di ogni file letto. Uso dalla root del repo.
set -euo pipefail
WHEN=${1:?before|after}; OUT=${2:-results/R4/R4_md5_${WHEN}.txt}
M=/mnt/micron/geo_spatialtrans; D=$M/data_real
mkdir -p "$(dirname "$OUT")"
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for s in square_002um square_008um; do
    b=$D/$ds/binned_outputs_x/binned_outputs/$s
    echo "$b/filtered_feature_bc_matrix.h5"; echo "$b/raw_feature_bc_matrix.h5"; echo "$b/spatial/tissue_positions.parquet"
  done
  echo "$D/$ds/${ds}_metrics_summary.csv"
done
ls $M/R2/masks/A1_r?_cellpose_rgb_native.npz $M/R2/masks/A2_r?_cellpose_rgb_native.npz $M/R2/masks/A3_r?_cellpose_rgb_native.npz $M/R2/masks/A5_r?_cellpose_rgb_native.npz
ls $M/R3b/masks/A4_*_consensus_native.npz $M/R3b/masks/A6_*_consensus_native.npz $M/R3b/masks/A4_*_consensus_map.parquet $M/R3b/masks/A6_*_consensus_map.parquet
ls $M/R3b/nuclei/A4_*_consensus_native.parquet $M/R3b/nuclei/A6_*_consensus_native.parquet
ls results/R2/nuclei/A[1235]_*_cellpose_rgb_native.parquet
ls $M/R6/counts/*_exp.npz $M/R6/counts/*_nuc.npz $M/R6/counts/*_shift.npz $M/R6/counts/*_nuclei.parquet
ls $M/R6/labels/A?_exp_umi20_seed20261022_nuclei_ph1.parquet
echo results/R2/R2_rois_checked.csv; echo results/R6/R6_dexp.csv
ls $M/R6/counts/*_genes.tsv
} | sort -u | xargs md5sum > "$OUT"
echo "$(wc -l < "$OUT") files -> $OUT"
