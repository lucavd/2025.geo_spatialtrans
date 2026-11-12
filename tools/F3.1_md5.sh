#!/usr/bin/env bash
# tools/F3.1_md5.sh - F3.1 (G-4): md5 dei file letti dallo step, prima e dopo. Uso: tools/F3.1_md5.sh before|after
set -euo pipefail
cd ~/2025.geo_spatialtrans
M=/mnt/micron/geo_spatialtrans; OUT=results/F3.1; mkdir -p $OUT
{
  for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
    S=$M/data_real/$ds/binned_outputs_x/binned_outputs/square_002um
    echo $S/raw_feature_bc_matrix.h5; echo $S/filtered_feature_bc_matrix.h5; echo $S/spatial/tissue_positions.parquet
  done
  ls $M/R2/masks/*_cellpose_rgb_native.npz $M/R3b/masks/*_consensus_native.npz $M/R3b/nuclei/*_consensus_native.parquet
  ls results/R2/nuclei/*_cellpose_rgb_native.parquet results/R2/R2_rois_checked.csv
  ls $M/R6/counts/*_exp.npz $M/R6/counts/*_nuclei.parquet $M/R6/counts/*_genes.tsv $M/R6/labels/*_exp_umi20_seed20261022_nuclei_ph1.parquet
  ls $M/S1.6/layers/*_U6.rds $M/S2.2/profiles/* $M/S2.5/real/S2.5_real_areas.csv
  ls $M/S2.0/run1_wfix/params_flat.json $M/S2.3/layers/*_regions.tsv
  ls $M/S2.0/run1_wfix/section/*.npz 2>/dev/null || true
} | sort -u > $OUT/F3.1_md5_files.txt
xargs -a $OUT/F3.1_md5_files.txt -P 6 -n 20 md5sum | sort -k2 > $OUT/F3.1_md5_$1.txt
wc -l < $OUT/F3.1_md5_$1.txt
if [ "$1" = after ]; then diff <(cat $OUT/F3.1_md5_before.txt) <(cat $OUT/F3.1_md5_after.txt) && echo "G-4 md5 identical"; fi
