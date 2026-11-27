#!/usr/bin/env bash
# tools/sonda/S9_U_chain.sh - S9 catena U: tools/sonda/S3_chain.sh (invariato) su S9_sonda/U/jobs/all.tsv con S3_ROOT = S9_sonda/U,
# S3_REF = S9_sonda/U/ref, W processi (al piu' 20). Guardia termica: Tctl ogni 60 s in $U/thermal.log; due letture consecutive > 90 C ->
# arresto del gruppo di processi e rilancio con W - 4 (la catena e' rilanciabile: salta cio' che esiste). Uso: bash tools/sonda/S9_U_chain.sh [W]
set -uo pipefail
cd ~/2025.geo_spatialtrans
U=/mnt/micron/geo_spatialtrans/S9_sonda/U; W=${1:-20}; [ "$W" -gt 20 ] && W=20
export S3_ROOT=$U S3_REF=$U/ref OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
tctl() { sensors 2>/dev/null | awk '/Tctl/ {gsub(/[+°C]/, "", $2); print $2; exit}'; }
echo "start $(date +%s) W=$W" >> $U/thermal.log
while true; do
  setsid bash tools/sonda/S3_chain.sh $U/jobs/all.tsv $W >> $U/chain_run.log 2>&1 &
  CP=$!; hot=0; stopped=0
  while kill -0 $CP 2>/dev/null; do
    sleep 60; t=$(tctl); echo "$(date +%s) W=$W Tctl=$t load=$(cut -d' ' -f1 /proc/loadavg)" >> $U/thermal.log
    if awk -v t="$t" 'BEGIN{exit !(t > 90)}'; then hot=$((hot + 1)); else hot=0; fi
    if [ $hot -ge 2 ] && [ $W -gt 4 ]; then
      echo "$(date +%s) HOT: stop group $CP, W $W -> $((W - 4))" >> $U/thermal.log
      kill -TERM -- -$CP 2>/dev/null; sleep 5; kill -KILL -- -$CP 2>/dev/null; wait $CP 2>/dev/null; W=$((W - 4)); stopped=1; break
    fi
  done
  wait $CP 2>/dev/null
  [ $stopped -eq 1 ] && continue
  break
done
echo "end $(date +%s) W=$W" >> $U/thermal.log
