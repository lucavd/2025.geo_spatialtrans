#!/usr/bin/env bash
# tools/R4_repro.sh — cancello G-4: seconda esecuzione completa (stadi, analisi e post hoc) in un'altra cartella, md5 di tutte le tabelle.
set -euo pipefail
A=${1:-/mnt/micron/geo_spatialtrans/R4/run1}; B=${2:-/mnt/micron/geo_spatialtrans/R4/run2}
tools/R4_run.sh "$B"
(cd "$A/tables" && md5sum *.csv | sort -k2) > results/R4/R4_md5_run1.txt
(cd "$B/tables" && md5sum *.csv | sort -k2) > results/R4/R4_md5_run2.txt
if diff results/R4/R4_md5_run1.txt results/R4/R4_md5_run2.txt > /dev/null; then echo "G-4 PASS: $(wc -l < results/R4/R4_md5_run1.txt) tabelle identiche"; else echo "G-4 FAIL"; diff results/R4/R4_md5_run1.txt results/R4/R4_md5_run2.txt || true; fi
