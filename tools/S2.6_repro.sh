#!/usr/bin/env bash
# tools/S2.6_repro.sh - S2.6, G-3: un'esecuzione a pezzi (cuore, pezzi di 504 um, k = 0) rifatta con 2 processi invece di 4 in un'altra
# cartella, misure rifatte; md5 di UMI, matrice a 8 um e tabelle delle misure confrontati. Uso: tools/S2.6_repro.sh [outfile]
set -euo pipefail
M=/mnt/micron/geo_spatialtrans/S2.6; R=$M/repro; OUT=${1:-results/S2.6/S2.6_G3_repro.txt}; mkdir -p $R
J=$M/q2/P504_A6_k0/run.json
L=$M/layers/$(python3 -c "import json;print(json.load(open('$J'))['layer'])"); S=$(python3 -c "import json;print(json.load(open('$J'))['seed'])")
rm -rf $R/P504_A6_k0
Rscript --vanilla tools/S2.6_pieces.R pieces $L 504 280 $S $R P504_A6_k0 0 2 > $R/log.txt 2>&1
Rscript --vanilla tools/S2.6_measures.R $R/P504_A6_k0 504 >> $R/log.txt 2>&1
{ for f in u2.bin own8.bin dup8.bin measures_T504.csv pairs_T504.csv; do
    a=$(md5sum < $M/q2/P504_A6_k0/$f | cut -c1-32); b=$(md5sum < $R/P504_A6_k0/$f | cut -c1-32); echo "$f $a $b $([ $a = $b ] && echo IDENTICO || echo DIVERSO)"; done
  Rscript --vanilla -e '.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); a <- readRDS(commandArgs(TRUE)[1]); b <- readRDS(commandArgs(TRUE)[2]); cat("c8.rds", if (identical(a, b)) "IDENTICO" else "DIVERSO", "\n")' $M/q2/P504_A6_k0/c8.rds $R/P504_A6_k0/c8.rds
} > $OUT
cat $OUT; ! grep -q DIVERSO $OUT
