#!/usr/bin/env bash
# tools/S2.5_repro.sh - S2.5, G-3: lato vero, pannello main con k = 0, misure del simulato rifatti in un'altra cartella con altri processi
# (vero 2 invece di 3, pannello 3 invece di 6); confronto md5 dei file per quadratino e dei JSON senza i tempi. Uso: tools/S2.5_repro.sh [out]
set -euo pipefail
M=/mnt/micron/geo_spatialtrans/S2.5; R2=$M/repro; OUT=${1:-results/S2.5/S2.5_G3.txt}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $R2/real $R2/panel; echo 2 > $R2/real/WORKERS; echo 3 > $R2/panel/WORKERS
.venv/bin/python -B tools/S2.5_real.py areas $R2/real > $R2/real/areas.log 2>&1
.venv/bin/python -B tools/S2.5_real.py summary $R2/real > $R2/real/summary.log 2>&1
S25_KS=0 Rscript --vanilla tools/S2.5_panel.R main $R2/panel 3 > $R2/panel/run_main.log 2>&1
.venv/bin/python -B tools/S2.5_stats.py $R2/panel main > $R2/panel/stats_main.log 2>&1
cmpjson() { python3 - "$1" "$2" <<'EOF'
import json, sys
drop = {"sec", "elapsed", "elapsed_total_s", "elapsed_sim_s", "peak_mb"}
a, b = (json.load(open(p)) for p in sys.argv[1:3])
for d in (a, b):
    for k in list(d):
        if k in drop: d.pop(k)
print("ok" if a == b else "diverso")
EOF
}
n=0; bad=0; : > $OUT
for f in $M/real/real/*.json; do t=$(basename $f); n=$((n+1)); r=$(cmpjson $f $R2/real/real/$t); [ "$r" = ok ] || { bad=$((bad+1)); echo "DIVERSO real $t" >> $OUT; }; done
for f in $M/real/real/*_blocks8.parquet; do t=$(basename $f); n=$((n+1)); cmp -s $f $R2/real/real/$t || { bad=$((bad+1)); echo "DIVERSO real $t" >> $OUT; }; done
n=$((n+1)); cmp -s $M/real/S2.5_real_refs.csv $R2/real/S2.5_real_refs.csv || { bad=$((bad+1)); echo "DIVERSO S2.5_real_refs.csv" >> $OUT; }
for f in $M/panel/main/*_k0.bin $M/panel/main/*_k0.g2.bin $M/panel/main/*_k0.u8.bin $M/panel/main/*_k0.g8.bin; do t=$(basename $f); n=$((n+1)); cmp -s $f $R2/panel/main/$t || { bad=$((bad+1)); echo "DIVERSO panel $t" >> $OUT; }; done
for f in $M/panel/main/*_k0.json $M/panel/stats_main/*_k0.json; do d=$(basename $(dirname $f)); t=$(basename $f); n=$((n+1)); r=$(cmpjson $f $R2/panel/$d/$t); [ "$r" = ok ] || { bad=$((bad+1)); echo "DIVERSO $d $t" >> $OUT; }; done
echo "G-3: $n file confrontati, $bad diversi" >> $OUT; cat $OUT; [ $bad = 0 ]
