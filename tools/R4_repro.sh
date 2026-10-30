#!/usr/bin/env bash
# tools/R4_repro.sh — cancello G-4: seconda esecuzione completa in un'altra cartella, md5 delle tabelle confrontati.
set -euo pipefail
A=/mnt/micron/geo_spatialtrans/R4/run1; B=/mnt/micron/geo_spatialtrans/R4/run2
tools/R4_run.sh "$B"
(cd "$A/tables" && md5sum *.csv | sort -k2) > /tmp/R4_md5_run1.txt
(cd "$B/tables" && md5sum *.csv | sort -k2) > /tmp/R4_md5_run2.txt
if diff /tmp/R4_md5_run1.txt /tmp/R4_md5_run2.txt > /dev/null; then echo "G-4 PASS: $(wc -l < /tmp/R4_md5_run1.txt) tabelle identiche"; else echo "G-4 FAIL"; diff /tmp/R4_md5_run1.txt /tmp/R4_md5_run2.txt || true; fi
