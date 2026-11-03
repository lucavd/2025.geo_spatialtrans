#!/usr/bin/env bash
# tools/R5_chain2.sh — R5 dopo la revisione numerica (ottimizzatore corretto): due esecuzioni complete in parallelo (run3, run4), ciascuna con
# calcolo e post hoc; poi G-4 (tools/R5_repro2.sh) e md5 dopo (G-2). Processi per esecuzione: <outdir>/WORKERS.
cd "$(dirname "$0")/.."
M=/mnt/micron/geo_spatialtrans/R5; S=$M/chain2_status.log
for O in $M/run3 $M/run4; do
  ( mkdir -p "$O/logs"; tools/R5_run.sh "$O" > "$O/logs/run.log" 2>&1 && tools/R5_posthoc.sh "$O" > "$O/logs/posthoc.log" 2>&1 \
      && echo "$O done $(date -Is)" >> "$S" || echo "$O FAILED $(date -Is)" >> "$S" ) &
done
wait
tools/R5_repro2.sh $M/run3 $M/run4 >> "$S" 2>&1
tools/R5_md5.sh after >> "$S" 2>&1
echo "chain2 done $(date -Is)" >> "$S"
