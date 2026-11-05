#!/usr/bin/env bash
# S2.2, cancello G-4: md5 dei file letti dallo step (dati veri, maschere, etichette di R6, tabelle dei passi precedenti).
# Uso: tools/S2.2_md5.sh before|after   -> results/S2.2/S2.2_md5_<fase>.txt
set -euo pipefail
cd "$(dirname "$0")/.."
PH=${1:?before|after}; M=/mnt/micron/geo_spatialtrans; OUT=results/S2.2; mkdir -p $OUT
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  md5sum $M/data_real/$ds/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5
  md5sum $M/data_real/$ds/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet
done
md5sum $M/R2/masks/*_cellpose_rgb_native.npz $M/R3b/masks/*_consensus_native.npz
md5sum $M/R6/labels/*_exp_umi20_seed20261022_nuclei_ph1.parquet $M/R6/counts/*_exp.npz $M/R6/counts/*_nuc.npz $M/R6/counts/*_nuclei.parquet $M/R6/counts/*_genes.tsv
md5sum results/R6/R6_dexp.csv results/R2/R2_rois_checked.csv results/R4/*.csv results/R6/*.csv
} | sed "s#$M#<micron>#" > $OUT/S2.2_md5_$PH.txt
wc -l $OUT/S2.2_md5_$PH.txt
