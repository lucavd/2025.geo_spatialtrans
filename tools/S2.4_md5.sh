#!/usr/bin/env bash
# S2.4, cancello G-4 (a): md5 dei file letti dallo step (dati veri, uscite dei passi precedenti). Uso: tools/S2.4_md5.sh before|after
set -euo pipefail
cd "$(dirname "$0")/.."
PH=${1:?before|after}; M=/mnt/micron/geo_spatialtrans; OUT=results/S2.4; mkdir -p $OUT
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  md5sum $M/data_real/$ds/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet
done
.venv/bin/python - <<'EOF' | xargs -r md5sum
import sys; sys.path.insert(0, "tools")
import R5_tiles as T5
for k in ("A1", "A23", "A4", "A5", "A6"): print(T5.image_info(T5.DS[k])[0])
EOF
md5sum $M/R5/run3/section/*.npz $M/R5/run3/kernels.json $M/R5/run3/params.json
md5sum $M/R2/masks/*_cellpose_rgb_native.npz $M/R3b/masks/*_consensus_native.npz $M/R3b/nuclei/*_consensus_native.parquet
md5sum $M/R4/run1/posthoc/*_cells.parquet
md5sum $M/S2.0/run1/params_flat.json $M/S2.0/run1_wfix/params_flat.json $M/S2.0/run1/stripes/*_factors.npz $M/S2.0/run1_wfix/stripes/*_factors.npz
md5sum $M/S1.6/layers/*_U6.rds
find $M/S2.2/profiles -type f | sort | xargs md5sum
find $M/S2.3/panel/hd $M/S2.3/layers -type f | sort | xargs md5sum
find $M/S2.3/panel/main -name '*.bin' | sort | xargs md5sum
md5sum results/R2/R2_rois_checked.csv results/R5/*.csv results/S2.0/S2.0_PH5_sections.csv results/S2.0/S2.0_PH5_calib.csv results/S2.3/S2.3_R2_fits_PH1_nocorr.csv results/S2.3/S2.3_R1_choice.csv
} | sed "s#$M#<micron>#" > $OUT/S2.4_md5_$PH.txt
wc -l $OUT/S2.4_md5_$PH.txt
