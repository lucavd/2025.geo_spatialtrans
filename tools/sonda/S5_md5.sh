#!/usr/bin/env bash
# tools/sonda/S5_md5.sh - S5 (docs/sonda/REGOLA_S5.md, cancello G-MD5): md5 dei file letti in sola lettura: dati di X1 (Visium HD dopo Xenium e
# Xenium), feature di H1 e H2, composita di S1 per X1, esecuzioni esistenti V, B, P (S3) e T (S3b), riferimento, aree omologhe, tabelle di S3b,
# S3c, S4 e orient.json di S5. Uso: bash tools/sonda/S5_md5.sh <out>
cd ~/2025.geo_spatialtrans; G=/mnt/micron/geo_spatialtrans; X=$G/data_lung/X1
{ md5sum $X/visium/binned/binned_outputs/square_002um/filtered_feature_bc_matrix.h5 $X/visium/binned/binned_outputs/square_002um/raw_feature_bc_matrix.h5 \
         $X/visium/binned/binned_outputs/square_002um/spatial/tissue_positions.parquet $X/visium/binned/binned_outputs/square_002um/spatial/scalefactors_json.json \
         $X/visium/binned/binned_outputs/square_008um/spatial/tissue_positions.parquet $X/visium/binned/binned_outputs/square_008um/spatial/scalefactors_json.json \
         $X/visium/*_tissue_image.btf $X/xenium/*_transcripts.parquet $X/xenium/*_nucleus_boundaries.parquet $X/xenium/*_cells.parquet \
         $X/xenium/*_gene_panel.json $X/xenium/*_he_imagealignment.csv \
         $G/data_lung/H2/visium/binned/binned_outputs/square_002um/filtered_feature_bc_matrix.h5 $G/data_lung/H1/visium/binned/binned_outputs/square_002um/filtered_feature_bc_matrix.h5 \
         $G/data_lung/H2/visium/binned/binned_outputs/square_008um/spatial/tissue_positions.parquet \
         $G/S1_sonda/X1/composite_parts.npy $G/S1_sonda/X1/setup.json $G/S5_sonda/explore/areas_hom_X1.csv
  find $G/S3_sonda/runs/V $G/S3_sonda/runs/B $G/S3_sonda/runs/P $G/S3b_sonda/runs/T $G/S3_sonda/ref -type f | sort | xargs md5sum
  md5sum results/sonda/S3b/area_values.csv results/sonda/S3b/h.csv results/sonda/S3c/D.json results/sonda/S3c/per_measure.csv results/sonda/S4/D.json \
         results/sonda/S4/area_values.csv results/sonda/S5/orient.json docs/sonda/REGOLA_S5.md; } | sed "s#$G/##" > $1
wc -l < $1
