#!/usr/bin/env bash
# tools/S2.7_md5.sh - S2.7, G-4: md5 di ogni file letto dallo step, prima e dopo; file di R/ cambiati rispetto all'inizio dello step.
# Uso (root del repo): tools/S2.7_md5.sh before|after [outdir]
set -euo pipefail
W=${1:?before o after}; OUT=${2:-results/S2.7}; mkdir -p "$OUT"
M=/mnt/micron/geo_spatialtrans; D=$M/data_real; BASE=54b52f2
{
for ds in Visium_HD_Mouse_Small_Intestine Visium_HD_Human_Colon_Cancer Visium_HD_Human_Lymph_Node_FFPE Visium_HD_6p5mm_Mouse_Brain Visium_HD_6p5mm_Human_Heart; do
  for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 spatial/tissue_positions.parquet; do
    p="$D/$ds/binned_outputs_x/binned_outputs/square_008um/$f"; [ -e "$p" ] && echo "$p"
  done
done
ls $M/S2.2/profiles/*
for t in A1_r1 A2_r1 A3_r1 A4_f1 A5_r1 A6_r1; do echo "$M/S1.6/layers/${t}_U6.rds"; echo "$M/S2.3/panel/hd/${t}.rds"; ls $M/S2.5/panel/main/${t}_k[01]*; done
ls $M/S2.6/layers/pc_*_2000um_*.rds $M/S2.6/sections/sec_*.rds
for a in A1 A4 A6; do ls $M/S2.6/q3/one_$a/*; done
echo results/S2.5/S2.5_sim_main.csv; echo results/R7/Mc/R7_Mc_refs.csv; ls results/S2.6/S2.6_section_*.csv
} | sort -u | xargs md5sum > "$OUT/S2.7_md5_${W}.txt"
git diff --name-only $BASE -- R/ | sort > "$OUT/S2.7_R_changed_${W}.txt"
wc -l < "$OUT/S2.7_md5_${W}.txt"
if [ "$W" = after ]; then
  diff "$OUT/S2.7_md5_before.txt" "$OUT/S2.7_md5_after.txt" > "$OUT/S2.7_md5_diff.txt" && echo "G-4 md5 identici" || { echo "G-4 md5 DIFFERENZE"; exit 1; }
  ALLOWED="R/05b2_simulate_cell_profiles.R R/05b3_place_and_move_transcripts.R R/05b6_simulate_visium_hd_pieces.R R/testing/fixtures/S2.7_cell624286.rds R/testing/test_S2.7.R"
  BAD=$(comm -23 "$OUT/S2.7_R_changed_after.txt" <(echo $ALLOWED | tr ' ' '\n' | sort) || true)
  [ -z "$BAD" ] && echo "G-4 R/ cambiato solo nei file dichiarati" || { echo "G-4 R/ file non dichiarati: $BAD"; exit 1; }
fi
