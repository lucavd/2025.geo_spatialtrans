#!/usr/bin/env bash
# S2.2: scarica gli atlanti a sonde (Chromium Flex) della gara sulla fonte dei profili e ne registra md5 e dimensione.
# Uso: tools/S2.2_atlas_fetch.sh crc   (gli altri atlanti si aggiungono con i loro link)
set -euo pipefail
cd "$(dirname "$0")/.."
A=/mnt/micron/geo_spatialtrans/S2.2/atlas; mkdir -p $A results/S2.2
case "${1:?atlante}" in
  crc)
    D=$A/crc; mkdir -p $D
    B=https://cf.10xgenomics.com/samples/cell-exp/8.0.0/HumanColonCancer_Flex_Multiplex/HumanColonCancer_Flex_Multiplex
    for s in _count_filtered_feature_bc_matrix.h5 _aggregation.csv; do [ -s $D/crc$s ] || curl -sSfL --retry 3 -o $D/crc$s "$B$s"; done
    G=https://raw.githubusercontent.com/10XGenomics/HumanColonCancer_VisiumHD/main/MetaData/SingleCell_MetaData.csv.gz
    [ -s $D/SingleCell_MetaData.csv.gz ] || curl -sSfL --retry 3 -o $D/SingleCell_MetaData.csv.gz "$G"
    ;;
  *) echo "atlante sconosciuto"; exit 1;;
esac
md5sum $A/*/* | sed "s#$A#<atlas>#" > results/S2.2/S2.2_atlas_md5.txt
ls -l $A/*/
