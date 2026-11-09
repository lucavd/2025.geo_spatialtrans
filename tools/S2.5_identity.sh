#!/usr/bin/env bash
# tools/S2.5_identity.sh - S2.5, G-4 (c): A2_r1 e A6_r1 con i semi dei passi e i blocchi di 2 000 cellule del pannello di S2.4 danno totali per
# quadratino identici a quelli di S2.4 (stessa catena, ora dentro simulate_visium_hd()). Uso: tools/S2.5_identity.sh <panel_dir> [out]
set -euo pipefail
PD=${1:?cartella del pannello}; OUT=${2:-results/S2.5/S2.5_G4_identity.txt}; mkdir -p "$(dirname "$OUT")"
Rscript --vanilla tools/S2.5_panel.R s24_identity "$PD" 1 > "$PD/run_s24_identity.log" 2>&1
: > "$OUT"; ok=1
for t in A2_r1 A6_r1; do
  a=$(md5sum < "$PD/s24_identity/${t}_k0.bin" | cut -d' ' -f1); b=$(md5sum < /mnt/micron/geo_spatialtrans/S2.4/panel/main/${t}_k0.bin | cut -d' ' -f1)
  if [ "$a" = "$b" ]; then echo "IDENTICO $t $a" >> "$OUT"; else echo "DIVERSO $t $a $b" >> "$OUT"; ok=0; fi
done
cat "$OUT"; [ $ok = 1 ]
