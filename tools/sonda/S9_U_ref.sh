#!/usr/bin/env bash
# tools/sonda/S9_U_ref.sh - S9 catena U, cancello G-REF: S9_reference.py HLCA in una cartella di prova (deve riprodurre S3_sonda/ref, md5)
# e riferimento LuCA in S9_sonda/U/ref. md5 dei due H5AD prima e dopo (G-MD5). Uso: bash tools/sonda/S9_U_ref.sh
set -uo pipefail
cd ~/2025.geo_spatialtrans
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans; U=$M/S9_sonda/U; OUT=results/sonda/S9/U; mkdir -p $U/md5 $OUT
H=$M/ref_lung/HLCA_core_066943a2.h5ad; L=$M/ref_lung/LuCA_core_46e0287b.h5ad
[ -f $U/md5/h5ad_before.md5 ] || md5sum $H $L > $U/md5/h5ad_before.md5
echo "ref start $(date)"
( [ -f $U/ref_hlca_check/ref_L.json ] || .venv/bin/python -B tools/sonda/S9_reference.py HLCA $U/ref_hlca_check > $U/ref_hlca_check.log 2>&1 ) &
( [ -f $U/ref/ref_L.json ] || .venv/bin/python -B tools/sonda/S9_reference.py LuCA $U/ref > $U/ref_luca.log 2>&1 ) &
wait
echo "ref end $(date)"
md5sum $H $L > $U/md5/h5ad_after_ref.md5
( cd $M/S3_sonda/ref && md5sum * ) > $U/md5/ref_S3.md5; ( cd $U/ref_hlca_check && md5sum * ) > $U/md5/ref_hlca_check.md5
( cd $U/ref && md5sum * ) > $U/md5/ref_luca.md5
cp $U/md5/*.md5 $OUT/; cp $U/ref/ref_L.json $OUT/ref_luca.json 2>/dev/null; cp $U/ref_hlca_check/ref_L.json $OUT/ref_hlca_check.json 2>/dev/null
tail -3 $U/ref_hlca_check.log $U/ref_luca.log
diff $U/md5/ref_S3.md5 $U/md5/ref_hlca_check.md5; diff $U/md5/h5ad_before.md5 $U/md5/h5ad_after_ref.md5 && echo "h5ad md5 unchanged"
true
