#!/usr/bin/env bash
# tools/sonda/S3_download.sh - S3: scarica le due sezioni solo Visium HD del blocco FFPE di polmone 10x (CG000709): H1 =
# Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1, H2 = ..._Experiment2 (nomi fissati in SC2, docs/STATO.md di 9613c16).
# File: metrics_summary.csv, spatial.tar.gz, tissue_image.btf (H&E a piena risoluzione), binned_outputs.tar.gz (estratti solo
# square_002um e square_008um, come X1 e X2). Dimensione di ogni file controllata contro content-length del CDN (S3_sizes.tsv);
# md5 in data_lung/MD5_download_H.txt. Nessuna lettura di espressione. Rilanciabile: salta i file completi.
# Uso: bash tools/sonda/S3_download.sh
set -uo pipefail
D=/mnt/micron/geo_spatialtrans/data_lung; B=https://cf.10xgenomics.com/samples/spatial-exp/3.0.0
: > $D/S3_sizes.tsv
for e in 1 2; do
  n=Visium_HD_Human_Lung_Cancer_HD_Only_Experiment$e; o=$D/H$e/visium; mkdir -p $o/binned
  for s in metrics_summary.csv spatial.tar.gz tissue_image.btf binned_outputs.tar.gz; do
    f=$o/${n}_$s; u=$B/$n/${n}_$s
    exp=$(curl -sIL --max-time 30 $u | grep -i '^content-length' | tail -1 | tr -dc 0-9)
    for try in 1 2 3; do
      [ -f $f ] && [ "$(stat -c %s $f)" = "$exp" ] && break
      curl -sL --retry 5 -C - -o $f $u || echo "curl rc=$? try $try $s"
    done
    got=$(stat -c %s $f 2>/dev/null || echo 0)
    printf 'H%s\t%s\t%s\t%s\t%s\n' $e $s $exp $got $([ "$exp" = "$got" ] && echo ok || echo MISMATCH) >> $D/S3_sizes.tsv
  done
  if [ ! -f $o/binned/.extracted ]; then
    tar -xzf $o/${n}_binned_outputs.tar.gz -C $o/binned --wildcards 'binned_outputs/square_002um/*' 'binned_outputs/square_008um/*' && touch $o/binned/.extracted
  fi
done
cd $D && md5sum H1/visium/*_metrics_summary.csv H1/visium/*_spatial.tar.gz H1/visium/*_binned_outputs.tar.gz H1/visium/*_tissue_image.btf \
  H2/visium/*_metrics_summary.csv H2/visium/*_spatial.tar.gz H2/visium/*_binned_outputs.tar.gz H2/visium/*_tissue_image.btf > MD5_download_H.txt
cat S3_sizes.tsv MD5_download_H.txt
