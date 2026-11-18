#!/usr/bin/env bash
# tools/F3.2_md5.sh - F3.2 (G-2 b): md5 dei dati di taratura letti dallo step, prima e dopo. Uso: tools/F3.2_md5.sh before|after
# Controllo storico in piu': i file gia' registrati da F3.1 (results/F3.1/F3.1_md5_after.txt) e da R1 (data/real/checksums.md5) devono avere lo stesso md5.
set -euo pipefail
cd ~/2025.geo_spatialtrans
M=/mnt/micron/geo_spatialtrans; OUT=results/F3.2; mkdir -p $OUT
{
  for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
    B=$M/data_real/$ds/binned_outputs_x/binned_outputs
    echo $B/square_002um/filtered_feature_bc_matrix.h5; echo $B/square_002um/spatial/tissue_positions.parquet
    echo $B/square_008um/filtered_feature_bc_matrix.h5; echo $B/square_008um/spatial/tissue_positions.parquet
    echo $B/square_008um/spatial/scalefactors_json.json; echo $B/square_008um/analysis/clustering/gene_expression_graphclust/clusters.csv
    ls $M/data_real/$ds/${ds}_tissue_image.* | grep -E '\.(btf|tif)$'
  done
  ls $M/S2.2/profiles/*
  ls results/R2/R2_rois_checked.csv results/R2b/R2b_windows.csv results/R5/R5_tiles.csv results/R7/Mc/R7_refs_adopted.csv results/R4/R4_roi.csv inst/extdata/size_resid.csv
} | sort -u > $OUT/F3.2_md5_files.txt
xargs -a $OUT/F3.2_md5_files.txt -P 6 -n 10 md5sum | sort -k2 > $OUT/F3.2_md5_$1.txt
wc -l < $OUT/F3.2_md5_$1.txt
python3 - $OUT/F3.2_md5_$1.txt <<'EOF'
import sys, os
cur = dict(l.split()[::-1] for l in open(sys.argv[1]) if l.strip())
hist = {}
for l in open("results/F3.1/F3.1_md5_after.txt"):
    if l.strip(): m, p = l.split(); hist[p] = (m, "F3.1")
for l in open("data/real/checksums.md5"):
    if l.strip():
        m, p = l.split(maxsplit=1); p = p.strip().lstrip("*")
        p = p[2:] if p.startswith("./") else p
        p = p[len("datasets/"):] if p.startswith("datasets/") else p
        full = os.path.join("/mnt/micron/geo_spatialtrans/data_real", p)
        hist[full] = (m, "R1")
n = bad = 0
for p, m in cur.items():
    if p in hist:
        n += 1
        if hist[p][0] != m: bad += 1; print("DIVERSO", hist[p][1], p)
print(f"confronto storico: {n} file gia' registrati, {bad} diversi")
EOF
if [ "$1" = after ]; then diff $OUT/F3.2_md5_before.txt $OUT/F3.2_md5_after.txt && echo "G-2 b md5 identici prima e dopo"; fi
