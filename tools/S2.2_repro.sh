#!/usr/bin/env bash
# S2.2, cancello G-3: gara e pannello rifatti in un'altra cartella con un numero di processi diverso; confronto con la prima esecuzione.
# Uso (root del repo): tools/S2.2_repro.sh  -> results/S2.2/S2.2_G3_repro.txt
set -uo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans/S2.2; R2=$M/repro; rm -rf $R2; mkdir -p $R2/race $R2/panel
for K in A1 A2 A3 A4 A5 A6; do echo $K; done | xargs -P 2 -I{} sh -c "S22_RACE_OUT=$R2/race .venv/bin/python tools/S2.2_race.py {} > $R2/race_{}.log 2>&1"
Rscript --vanilla tools/S2.2_panel.R check $R2/panel 2 > $R2/panel.log 2>&1
OUT=results/S2.2/S2.2_G3_repro.txt
{
echo "gara: file CSV e JSON per area (prima esecuzione: 4 processi; seconda: 2 processi, cartella $R2/race)"
n=0; d=0; for f in $M/race/*.csv $M/race/*_real.json; do b=$(basename $f); n=$((n+1)); cmp -s $f $R2/race/$b || { d=$((d+1)); echo "DIVERSO $b"; }; done; echo "gara: $n file, $d diversi"
.venv/bin/python - <<EOF
import glob, os, numpy as np
n = d = 0
for f in sorted(glob.glob("$M/race/*_p5.npz")):
    g = "$R2/race/" + os.path.basename(f); a, b = np.load(f), np.load(g); n += 1
    if sorted(a.files) != sorted(b.files) or any(not np.array_equal(a[k], b[k]) for k in a.files): d += 1; print("DIVERSO", os.path.basename(f))
print(f"gara, indice degli zeri (npz, confronto dei contenuti): {n} file, {d} diversi")
EOF
echo "pannello: file RDS (prima esecuzione: 4 processi; seconda: 2 processi, cartella $R2/panel)"
n=0; d=0; for f in $M/panel/*.rds; do b=$(basename $f); n=$((n+1)); cmp -s $f $R2/panel/$b || { d=$((d+1)); echo "DIVERSO $b"; }; done; echo "pannello: $n file, $d diversi"
cmp -s results/S2.2/S2.2_G2_panel.csv $R2/panel/S2.2_G2_panel.csv && echo "tabella G-2: identica a meno dei tempi?" || echo "tabella G-2: diversa (contiene i secondi di calcolo, atteso)"
cmp -s results/S2.2/S2.2_G2_recovery.csv $R2/panel/S2.2_G2_recovery.csv && echo "tabella dei profili ricostruiti: identica" || echo "tabella dei profili ricostruiti: DIVERSA"
} > $OUT
cat $OUT
