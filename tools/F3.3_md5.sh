#!/usr/bin/env bash
# tools/F3.3_md5.sh - F3.3, cancello G-2: md5 di ogni file letto, prima e dopo. Uso: bash tools/F3.3_md5.sh before|after
set -euo pipefail
cd ~/2025.geo_spatialtrans
M=/mnt/micron/geo_spatialtrans; OUT=$M/F3.3/md5; mkdir -p $OUT
{
  for s in $(tail -n +2 results/F3.3/F3.3_areas.csv | cut -d, -f6 | sort -u); do :; done
  python3 - <<'EOF'
import pandas as pd
A = pd.read_csv("results/F3.3/F3.3_areas.csv")
for (d, img, sq), _ in A.groupby(["ds_dir", "image", "sq"]):
    print(img)
    for b in ("square_002um", "square_008um"):
        for f in ("raw_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5", "spatial/tissue_positions.parquet", "spatial/scalefactors_json.json"):
            print(f"{sq}/{b}/{f}")
    g = f"{d}/segmented_outputs_x/segmented_outputs/nucleus_segmentations.geojson"
    import os
    if os.path.exists(g): print(g)
EOF
  find $M/R2/rois $M/R2/masks $M/R3b/masks $M/R3b/nuclei $M/S2.3/layers $M/F3.1/main/ref -type f
  find $M/F3.1/main/real $M/F3.1/main/sim -type f \( -name 'bins.parquet' -o -name 'X.parquet' -o -name 'cells.parquet' -o -name 'rctd_blocks.csv' -o -name 'rctd.json' -o -name 'b2c_label.json' -o -name 'export.json' -o -name 'genes.tsv' \)
  find results/R2/tissue_masks results/R2/nuclei results/F3.1 results/F3.2 results/R2/R2_rois_checked.csv data/holdout/F3.2_download_manifest.csv -type f
} | sort -u > $OUT/files_$1.txt
xargs -a $OUT/files_$1.txt -d '\n' -P 4 -n 50 md5sum | sort -k2 > $OUT/md5_$1.txt
wc -l $OUT/md5_$1.txt
