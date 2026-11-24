#!/usr/bin/env bash
# tools/sonda/S5_run.sh - S5 (docs/sonda/REGOLA_S5.md): corsa completa, rilanciabile (ogni passo salta cio' che esiste; un file per unita').
#   prep   md5 prima (G-MD5); copie ristrette Vp, Tp, Bp, Pp (tools/sonda/S5_subset.py, G-SUB); XB: prep, correzione affine locale e metro
#          (tools/sonda/S5_xenium.py reg, G-XR), esportazione, diradamento; XV: ritaglio, StarDist, esportazione, restrizione (tools/sonda/S5_xv.py)
#   pilot  catena di REGOLA_S3 (tools/sonda/S3_chain.sh, invariato) su un'unita' per braccio, con i tempi
#   all    catena sulle 70 esecuzioni (processi dal file $S5_ROOT/WORKERS, al piu' 12)
#   post   tabella (S5_table), h_p (S5_h), D_p e A3 (S5_D), figura (S5_figure), md5 dopo, cancelli (S5_gates)
# Uso: bash tools/sonda/S5_run.sh <prep|pilot|all|post>
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S5_ROOT=/mnt/micron/geo_spatialtrans/S5_sonda; export S3_ROOT=$S5_ROOT S3_REF=/mnt/micron/geo_spatialtrans/S3_sonda/ref
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $S5_ROOT/logs results/sonda/S5; [ -f $S5_ROOT/WORKERS ] || echo 12 > $S5_ROOT/WORKERS
W=$(cat $S5_ROOT/WORKERS); [ "$W" -gt 12 ] && W=12; L=$S5_ROOT/logs; PY=.venv/bin/python
case $1 in
prep)
  [ -s results/sonda/S5/md5_before.txt ] || bash tools/sonda/S5_md5.sh results/sonda/S5/md5_before.txt
  [ -f results/sonda/S5/gate_sub.csv ] || $PY -B tools/sonda/S5_subset.py > $L/subset.log 2>&1 || { echo "SUBSET FAIL"; exit 1; }
  [ -f $S5_ROOT/xen/X1_stroma_5/prep.json ] || $PY -B tools/sonda/S5_xenium.py prep > $L/xen_prep.log 2>&1 || { echo "PREP FAIL"; exit 1; }
  $PY -B tools/sonda/S5_xenium.py reg $W >> $L/xen_reg.log 2>&1 || { echo "REG FAIL"; exit 1; }
  $PY -B tools/sonda/S5_xenium.py export $W >> $L/xen_export.log 2>&1 || { echo "EXPORT XB FAIL"; exit 1; }
  $PY -B tools/sonda/S5_xenium.py thin >> $L/xen_thin.log 2>&1 || { echo "THIN FAIL"; exit 1; }
  [ -f $S5_ROOT/seg/X1_stroma_5/crop.json ] || $PY -B tools/sonda/S5_xv.py crop > $L/xv_crop.log 2>&1 || { echo "CROP FAIL"; exit 1; }
  TF_NUM_INTEROP_THREADS=4 TF_NUM_INTRAOP_THREADS=4 bash tools/py_stardist.sh tools/sonda/S3_stardist.py all >> $L/xv_stardist.log 2>&1 || { echo "STARDIST FAIL"; exit 1; }
  $PY -B tools/sonda/S5_xv.py export >> $L/xv_export.log 2>&1 || { echo "EXPORT XV FAIL"; exit 1; }
  $PY -B tools/sonda/S5_xv.py restrict >> $L/xv_restrict.log 2>&1 || { echo "RESTRICT XV FAIL"; exit 1; }
  python3 tools/sonda/S5_jobs.py ;;
pilot|all)
  python3 tools/sonda/S5_jobs.py
  echo "chain $1 start $(date)"; bash tools/sonda/S3_chain.sh $S5_ROOT/jobs/$1.tsv $W; echo "chain $1 end $(date)" ;;
post)
  python3 tools/sonda/S5_table.py && $PY -B tools/sonda/S5_h.py && $PY -B tools/sonda/S5_D.py && $PY -B tools/sonda/S5_figure.py || { echo "POST FAIL"; exit 1; }
  bash tools/sonda/S5_md5.sh results/sonda/S5/md5_after.txt && $PY -B tools/sonda/S5_gates.py || { echo "GATES FAIL"; exit 1; } ;;
esac
