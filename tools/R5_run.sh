#!/usr/bin/env bash
# tools/R5_run.sh <outdir> — R5, esecuzione completa e rilanciabile (ogni tappa salta cio' che esiste gia').
# Tappe: sezioni, tessere di bordo, adattamenti, stime per sezione, lumi (CP-2), rapporto nucleo/citoplasma (P-5, CP-1), tabelle.
# Processi: <outdir>/WORKERS (24). Pre-registrazione: docs/steps/R5.md (c5c9fce, ceca77b).
set -euo pipefail
cd "$(dirname "$0")/.."
O=$1; mkdir -p "$O/logs"; [ -f "$O/WORKERS" ] || echo 24 > "$O/WORKERS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=.venv/bin/python
for k in A1 A23 A4 A5 A6; do [ -f "$O/section/${k}_G3.json" ] || $PY tools/R5_section.py $k "$O" >> "$O/logs/section.log" 2>&1; done
pids=""; for k in A1 A23 A4 A5 A6; do if [ ! -f "$O/tiles/${k}_candidates.csv" ]; then $PY tools/R5_tiles.py $k "$O" > "$O/logs/tiles_$k.log" 2>&1 & pids="$pids $!"; fi; done
for p in $pids; do wait $p; done
$PY tools/R5_fit.py "$O" >> "$O/logs/fit.log" 2>&1
$PY tools/R5_fit_s0.py "$O" >> "$O/logs/fit.log" 2>&1
$PY tools/R5_analysis.py edges "$O" > "$O/logs/analysis_edges.log" 2>&1
pids=""; for k in A1 A23 A4 A5 A6; do if [ ! -f "$O/lumens/${k}_lumens.csv" ]; then $PY tools/R5_lumens.py $k "$O" "$O/params.json" > "$O/logs/lumens_$k.log" 2>&1 & pids="$pids $!"; fi; done
for p in $pids; do wait $p; done
$PY tools/R5_nc.py "$O" "$O/kernels.json" >> "$O/logs/nc.log" 2>&1
$PY tools/R5_analysis.py nc "$O" > "$O/logs/analysis_nc.log" 2>&1
$PY tools/R5_analysis.py lumens "$O" > "$O/logs/analysis_lumens.log" 2>&1
echo "R5_run done $(date -Is)" >> "$O/logs/status.log"
