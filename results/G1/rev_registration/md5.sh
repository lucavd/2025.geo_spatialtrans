#!/bin/bash
# usage: md5.sh <outfile>
D=/mnt/micron/geo_spatialtrans/data_lung
OUT=/mnt/micron/geo_spatialtrans/G1_review/rev_registration
cd $D
while read f; do md5sum "$f"; done < $OUT/files.txt > $OUT/$1
