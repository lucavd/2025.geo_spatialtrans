#!/bin/bash
# md5 of the extracted files used here vs the same members inside the downloaded archive
OUT=/mnt/micron/geo_spatialtrans/G1_review/rev_registration
T=/mnt/micron/geo_spatialtrans/data_lung/X1/visium/Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1_binned_outputs.tar.gz
mkdir -p $OUT/tarcheck && cd $OUT/tarcheck
tar -xzf $T --wildcards '*square_002um/filtered_feature_bc_matrix.h5' '*square_002um/spatial/tissue_positions.parquet' '*square_002um/spatial/scalefactors_json.json' '*square_008um/spatial/tissue_positions.parquet'
find . -type f -exec md5sum {} \; | sort -k2 > $OUT/md5_tar_members.txt
cd $OUT && rm -rf tarcheck
