#!/usr/bin/env bash
# S2.4, cancello G-4 (c): due aree del pannello di S2.3 rifatte con il codice di S2.4 (seconda coda assente) devono dare quadratini identici.
# Le altre aree sono saltate copiando i loro .json di S2.3 nella cartella temporanea (lo script di S2.3 salta le aree gia' fatte).
set -euo pipefail
cd "$(dirname "$0")/.."
S=/mnt/micron/geo_spatialtrans/S2.3/panel; T=/mnt/micron/geo_spatialtrans/S2.4/s23_identity; rm -rf $T; mkdir -p $T/main
ln -s $S/hd $T/hd
KEEP="A2_r1 A5_r3"
for f in $S/main/*_k?.json; do b=$(basename $f .json); t=${b%_k?}; case " $KEEP " in *" $t "*) ;; *) cp $f $T/main/;; esac; done
Rscript --vanilla tools/S2.3_panel.R main $T 2 > $T/run.log 2>&1 || true
ok=0; n=0
for t in $KEEP; do for k in 0 1; do n=$((n+1)); a=$(md5sum < $S/main/${t}_k$k.bin); b=$(md5sum < $T/main/${t}_k$k.bin); [ "$a" = "$b" ] && ok=$((ok+1)) && echo "IDENTICO ${t}_k$k" || echo "DIVERSO ${t}_k$k"; done; done
echo "G-4 (c): $ok/$n file identici" | tee results/S2.4/S2.4_G4c.txt
