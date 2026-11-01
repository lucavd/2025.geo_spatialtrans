#!/usr/bin/env bash
# tools/R5_repro.sh — G-4: seconda esecuzione completa in un'altra cartella con lo stesso seme; tabelle con md5 identici.
set -euo pipefail
cd "$(dirname "$0")/.."
A=/mnt/micron/geo_spatialtrans/R5/run1; B=/mnt/micron/geo_spatialtrans/R5/run2
tools/R5_run.sh "$B"
{ for f in $(cd "$A/tables" && ls *.csv) kernels.json params.json; do
    pa=$A/tables/$f; pb=$B/tables/$f; [ -f "$pa" ] || { pa=$A/$f; pb=$B/$f; }
    ma=$(md5sum < "$pa" | cut -c1-32); mb=$(md5sum < "$pb" | cut -c1-32)
    echo "$f $ma $mb $([ "$ma" = "$mb" ] && echo identical || echo DIFFERENT)"; done; } > results/R5/R5_G4.txt
d=$(grep -c DIFFERENT results/R5/R5_G4.txt || true); t=$(wc -l < results/R5/R5_G4.txt)
echo "G-4: $((t - d))/$t tables identical" >> results/R5/R5_G4.txt
