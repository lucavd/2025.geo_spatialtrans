#!/usr/bin/env bash
# tools/riesame_panel_run.sh - riesame, Parte B: catena di misura condivisa (bracci V, T, S, B e, dopo, Xenium).
#
# README. Data un'esportazione di un'area nel formato di F3.1/F3.3 (cartella con bins.parquet, X.parquet, genes.tsv, cells.parquet;
# per il simulato anche origin/profiles/background), la catena:
#   (a) restringe i conteggi a un insieme di geni (tools/riesame_panel_restrict.py: panel = 408 geni del pannello Xenium dentro i
#       18 024 comuni; common = 18 024 comuni; all = nessuna restrizione);
#   (b) esegue bin2cell come F3.1 (tools/F3.1_b2c.py, .venv-b2c, bin2cell 0.3.5) e RCTD come F3.1 (tools/F3.1_rctd.R, modalita'
#       doppietto, seme 20261118) con il riferimento di F3.1 (/mnt/micron/geo_spatialtrans/F3.1/main/ref): F3.1_rctd.R usa le sole
#       righe del riferimento dei geni dell'area, cioe' il riferimento ristretto agli stessi geni (scelta C-2 di docs/riesame/DECISIONI_B.md);
#   (c) calcola M-1...M-7 con le etichette trasferite di F3.3 (tools/riesame_panel_measures.py, funzioni di tools/F3_3_common.py).
# Uscite: <B>/<arm>/runs/<run>_<genes>/ (una cartella per esecuzione) e <B>/<arm>/measures/<run>_<genes>.json; la tabella per area si
# costruisce con tools/riesame_panel_table.py. Rilanciabile: salta cio' che esiste.
#
# Comando per una sola area (dalla radice del repository, sul server):
#   printf 'V\tpanel\t/mnt/micron/geo_spatialtrans/F3.1/main/real/A2_r1\tA2\treal\tr1\tA2_r1\tA2_r1\tIN\tA2_IN\tIN\tA2_r1\n' > /tmp/one.tsv
#   bash tools/riesame_panel_run.sh /tmp/one.tsv 1
# Colonne della lista (TSV, senza intestazione): arm genes src_dir A side(real|sim) roi(rN per le aree vere di taratura, - altrove)
#   unit run src set level u2. Per un'esportazione nuova (per esempio Xenium) basta una riga con arm, src_dir e A del tessuto.
# Uso: bash tools/riesame_panel_run.sh <jobs.tsv> [W=2] [B_root=/mnt/micron/geo_spatialtrans/riesame/B]
set -uo pipefail
cd ~/2025.geo_spatialtrans
JOBS=$1; W=${2:-2}; BR=${3:-/mnt/micron/geo_spatialtrans/riesame/B}; REF=/mnt/micron/geo_spatialtrans/F3.1/main/ref
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BR REF
one() {
  IFS=$'\t' read -r arm genes src A side roi unit run srcl set level u2 <<< "$1"
  local d=$BR/$arm/runs/${run}_${genes}; local o=$BR/$arm/measures/${run}_${genes}.json; mkdir -p $d $BR/$arm/measures
  local bs=$side; [ "$side" = "sim" ] && [ "$genes" != "all" ] && bs=simr
  local t0=$(date +%s.%N)
  [ -f $d/chain_restrict.json ] || .venv/bin/python -B tools/riesame_panel_restrict.py $src $d $genes $side > $d/restrict.log 2>&1 || { echo "RESTRICT FAIL $run $genes"; return; }
  local t1=$(date +%s.%N)
  [ -f $d/b2c_label.json ] || .venv-b2c/bin/python -B tools/F3.1_b2c.py $d $REF $A $bs label > $d/b2c_label.log 2>&1 || { echo "B2C FAIL $run $genes"; return; }
  local t2=$(date +%s.%N)
  [ -f $d/rctd.json ] || Rscript --vanilla tools/F3.1_rctd.R $d $REF $A $bs > $d/rctd.log 2>&1 || { echo "RCTD FAIL $run $genes"; return; }
  local t3=$(date +%s.%N)
  [ -s $o ] || .venv-b2c/bin/python -B tools/riesame_panel_measures.py $d $A $roi $o arm=$arm genes=$genes unit=$unit run=$run src=$srcl set=$set level=$level u2=$u2 > $o.log 2>&1 || { echo "MEASURE FAIL $run $genes"; return; }
  local t4=$(date +%s.%N)
  [ -f $d/chain_times.tsv ] || printf 'restrict\tb2c\trctd\tmeasures\n%.1f\t%.1f\t%.1f\t%.1f\n' $(echo "$t1-$t0" | bc) $(echo "$t2-$t1" | bc) $(echo "$t3-$t2" | bc) $(echo "$t4-$t3" | bc) > $d/chain_times.tsv
  echo "done $arm $run $genes"
}
export -f one
tr -d '\r' < $JOBS | grep -v '^#' | grep -v '^$' | xargs -d '\n' -P $W -I{} bash -c 'one "$@"' _ {}
