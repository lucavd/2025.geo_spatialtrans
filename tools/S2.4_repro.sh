#!/usr/bin/env bash
# S2.4, cancello G-3: R-1, R-3 e pannello (variante main, seme k = 0) rifatti in un'altra cartella con un numero di processi diverso.
# Confronto: md5 dei .bin del pannello e delle tabelle (R-1 calibrazione, curva e generatore; R-3 adattamenti e fasce, senza la colonna dei
# secondi), e dei .json senza i campi di tempo e di memoria. Uso: tools/S2.4_repro.sh
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD; M=/mnt/micron/geo_spatialtrans/S2.4; T=$M/repro; rm -rf $T; mkdir -p $T
rsync -a --exclude .git --exclude renv/library --exclude results/S2.4 ./ $T/repo/; ln -s $ROOT/renv/library $T/repo/renv/library; mkdir -p $T/repo/results/S2.4
cp results/S2.4/S2.4_E1_real.csv results/S2.4/S2.4_E12_summary.csv results/S2.4/S2.4_E2_sim.csv results/S2.4/S2.4_R1_generator.csv results/S2.4/S2.4_R3_fits.csv $T/repo/results/S2.4/
cd $T/repo
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python -B tools/S2.4_R1_calib.py 3 > $T/r1.log 2>&1
mkdir -p $T/R3 && echo 2 > $T/R3/WORKERS && .venv/bin/python -B tools/S2.4_R3_glass.py fit $T/R3 > $T/r3.log 2>&1
S24_KS=0 Rscript --vanilla tools/S2.4_panel.R main $T/panel 3 > $T/panel.log 2>&1
cd $ROOT
{
for f in S2.4_R1_calib.csv S2.4_R1_curve.csv S2.4_R1_generator.csv; do a=$(md5sum < results/S2.4/$f); b=$(md5sum < $T/repo/results/S2.4/$f); [ "$a" = "$b" ] && echo "IDENTICO $f" || echo "DIVERSO $f"; done
for f in S2.4_R3_fits.csv S2.4_P5_glass.csv; do a=$(python3 -c "import pandas as pd,sys; d=pd.read_csv(sys.argv[1]); print(d.drop(columns=[c for c in ('sec',) if c in d]).to_csv(index=False))" results/S2.4/$f | md5sum); b=$(python3 -c "import pandas as pd,sys; d=pd.read_csv(sys.argv[1]); print(d.drop(columns=[c for c in ('sec',) if c in d]).to_csv(index=False))" $T/repo/results/S2.4/$f | md5sum); [ "$a" = "$b" ] && echo "IDENTICO $f (senza secondi)" || echo "DIVERSO $f"; done
n=0; ok=0; okj=0
for f in $M/panel/main/*_k0.bin; do b=$(basename $f); n=$((n+1)); [ "$(md5sum < $f)" = "$(md5sum < $T/panel/main/$b)" ] && ok=$((ok+1))
  j=${b%.bin}.json; s1=$(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); [d.pop(k) for k in list(d) if k.startswith('elapsed') or k == 'peak_mb']; print(json.dumps(d, sort_keys=True))" $M/panel/main/$j | md5sum)
  s2=$(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); [d.pop(k) for k in list(d) if k.startswith('elapsed') or k == 'peak_mb']; print(json.dumps(d, sort_keys=True))" $T/panel/main/$j | md5sum); [ "$s1" = "$s2" ] && okj=$((okj+1)); done
echo "pannello: .bin identici $ok/$n; .json identici senza tempi e memoria $okj/$n"
} | tee results/S2.4/S2.4_G3.txt
