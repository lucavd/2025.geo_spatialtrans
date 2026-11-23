#!/usr/bin/env bash
# tools/sonda/S3_chain.sh - S3: catena di misura sul polmone (docs/sonda/REGOLA_S3.md, «Catena»): per ogni esecuzione
# bin2cell (tools/F3.1_b2c.py, invariato), RCTD doppietto (tools/F3.1_rctd.R, invariato, seme 20261118) con il riferimento HLCA (sigla L),
# misure M-1...M-6 con etichette trasferite (tools/sonda/S3_measures.py). Rilanciabile: salta cio' che esiste; un file per unita'.
# Lista (TSV, senza intestazione): arm unit run_dir tissue side(real|simr) set k
# Uso: bash tools/sonda/S3_chain.sh <jobs.tsv> [W] ; W predefinito = contenuto di $S3_ROOT/WORKERS (massimo 20).
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S3_ROOT=${S3_ROOT:-/mnt/micron/geo_spatialtrans/S3_sonda}; export S3_REF=${S3_REF:-/mnt/micron/geo_spatialtrans/S3_sonda/ref}
JOBS=$1; W=${2:-$(cat $S3_ROOT/WORKERS 2>/dev/null || echo 4)}; [ "$W" -gt 20 ] && W=20
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
one() {
  IFS=$'\t' read -r arm unit d tissue side set k <<< "$1"
  local o=$S3_ROOT/measures/$arm/$unit.json; mkdir -p $S3_ROOT/measures/$arm
  local t0=$(date +%s.%N)
  [ -f $d/b2c_label.json ] || .venv-b2c/bin/python -B tools/F3.1_b2c.py $d $S3_REF L $side label > $d/b2c_label.log 2>&1 || { echo "B2C FAIL $arm $unit"; return; }
  local t1=$(date +%s.%N)
  [ -f $d/rctd.json ] || Rscript --vanilla tools/F3.1_rctd.R $d $S3_REF L $side > $d/rctd.log 2>&1 || { echo "RCTD FAIL $arm $unit"; return; }
  local t2=$(date +%s.%N)
  [ -s $o ] || .venv-b2c/bin/python -B tools/sonda/S3_measures.py $d $tissue $o arm=$arm unit=$unit run=$unit set=$set k=$k > $o.log 2>&1 || { echo "MEASURE FAIL $arm $unit"; return; }
  local t3=$(date +%s.%N)
  [ -f $d/chain_times.tsv ] || printf 'b2c\trctd\tmeasures\n%.1f\t%.1f\t%.1f\n' $(echo "$t1-$t0" | bc) $(echo "$t2-$t1" | bc) $(echo "$t3-$t2" | bc) > $d/chain_times.tsv
  echo "done $arm $unit"
}
export -f one
tr -d '\r' < $JOBS | grep -v '^#' | grep -v '^$' | xargs -d '\n' -P $W -I{} bash -c 'one "$@"' _ {}
