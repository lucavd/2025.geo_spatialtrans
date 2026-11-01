#!/usr/bin/env bash
# tools/R5_chain.sh — catena notturna di R5 (run1 gia' avviato): attende la fine degli adattamenti, completa run1, poi G-4 (run2), poi md5 dopo.
cd "$(dirname "$0")/.."
A=/mnt/micron/geo_spatialtrans/R5/run1
while ! grep -q '^done' "$A/logs/fit.log" 2>/dev/null; do sleep 120; done
echo "fit done $(date -Is)" >> "$A/logs/status.log"
tools/R5_run.sh "$A" >> "$A/logs/status.log" 2>&1 || echo "run1 FAILED $(date -Is)" >> "$A/logs/status.log"
cp "$A"/tables/*.csv "$A"/kernels.json "$A"/params.json results/R5/ 2>/dev/null || true
tools/R5_repro.sh >> "$A/logs/status.log" 2>&1 || echo "repro FAILED $(date -Is)" >> "$A/logs/status.log"
tools/R5_md5.sh after >> "$A/logs/status.log" 2>&1
echo "chain done $(date -Is)" >> "$A/logs/status.log"
