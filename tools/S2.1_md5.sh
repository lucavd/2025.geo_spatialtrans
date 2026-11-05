#!/usr/bin/env bash
# tools/S2.1_md5.sh - S2.1, cancello G-4: md5 dei file letti, prima e dopo. Uso (root del repo): tools/S2.1_md5.sh before|after
# before_extra / after_extra: nuclei del segmentatore letti dal descrittivo (aggiunti prima di leggerli). after: confronta con before e controlla che results/R2, results/R3, results/S1.6 siano invariati rispetto al commit.
set -euo pipefail
STAGE=${1:?before o after}
M=/mnt/micron/geo_spatialtrans; D=$M/data_real; OUT=$M/S2.1; mkdir -p "$OUT"
list() {
  ls "$M"/S1.6/layers/*_U6.rds
  for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
    for s in square_002um square_008um; do
      echo "$D/$ds/binned_outputs_x/binned_outputs/$s/spatial/tissue_positions.parquet"
      echo "$D/$ds/binned_outputs_x/binned_outputs/$s/spatial/scalefactors_json.json"
    done
    if [ "$ds" != Visium_HD_Mouse_Small_Intestine ]; then
      echo "$D/$ds/${ds}_barcode_mappings.parquet"
      echo "$D/$ds/segmented_outputs_x/segmented_outputs/nucleus_segmentations.geojson"
      echo "$D/$ds/segmented_outputs_x/segmented_outputs/cell_segmentations.geojson"
    fi
  done
  ls "$M"/R2/masks/*_cellpose_rgb_native.npz
  ls "$M"/R3b/masks/A4_*_consensus_native.npz "$M"/R3b/masks/A6_*_consensus_native.npz
  echo results/R2/R2_rois_checked.csv
  ls results/R2/tissue_masks/*_valid_ds4.png
}
list_extra() { ls "$M"/R3b/nuclei/A4_*_consensus_native.parquet "$M"/R3b/nuclei/A6_*_consensus_native.parquet; ls results/R2/nuclei/*_cellpose_rgb_native.parquet; }
case "$STAGE" in before_extra|after_extra) list_extra | sort | xargs md5sum > "$OUT/md5_$STAGE.txt";; *) list | sort | xargs md5sum > "$OUT/md5_$STAGE.txt";; esac
echo "$(wc -l < "$OUT/md5_$STAGE.txt") file in $OUT/md5_$STAGE.txt"
if [ "$STAGE" = after_extra ]; then
  if diff -q "$OUT/md5_before_extra.txt" "$OUT/md5_after_extra.txt" > /dev/null; then echo "PASS G-4 md5 dei nuclei del descrittivo identici"; else echo "FAIL G-4 md5 dei nuclei del descrittivo diversi"; fi
fi
if [ "$STAGE" = after ]; then
  if diff -q "$OUT/md5_before.txt" "$OUT/md5_after.txt" > /dev/null; then echo "PASS G-4 md5 identici prima e dopo"; else echo "FAIL G-4 md5 diversi"; diff "$OUT/md5_before.txt" "$OUT/md5_after.txt" | head; fi
  ch=$(git status --porcelain results/R2 results/R3 results/S1.6 | wc -l)
  if [ "$ch" -eq 0 ]; then echo "PASS G-4 results/R2, R3, S1.6 invariati rispetto al commit"; else echo "FAIL G-4 $ch file cambiati in results/R2, R3, S1.6"; fi
fi
