#!/usr/bin/env bash
# S2.3, cancello G-3: banco (preparazione e verifica di R-1) e pannello (variante main) rifatti in un'altra cartella con un numero di
# processi diverso (3 invece di 6); confronto dei contenuti (i tempi di calcolo esclusi). Uso: tools/S2.3_repro.sh -> results/S2.3/S2.3_G3.txt
set -uo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/mnt/micron/geo_spatialtrans/S2.3; R=$M/repro; mkdir -p $R/bench $R/panel; echo 3 > $R/bench/WORKERS; echo 3 > $R/panel/WORKERS
.venv/bin/python tools/S2.3_bench.py prep $R/bench /mnt/micron/geo_spatialtrans/S2.0/run1_wfix/params_flat.json 1 > $R/prep.log 2>&1
cp $M/bench/combos_verify.json $R/bench/
.venv/bin/python tools/S2.3_bench.py grid $R/bench $R/bench/combos_verify.json 10 5 > $R/grid.log 2>&1
ln -sfn $M/panel/hd $R/panel/hd
Rscript --vanilla tools/S2.3_panel.R main $R/panel 3 > $R/panel.log 2>&1
.venv/bin/python tools/S2.3_panel_stats.py $R/panel main > $R/stats.log 2>&1
.venv/bin/python - <<'EOF' > results/S2.3/S2.3_G3.txt
import json, glob, os, pickle, hashlib, numpy as np
M = "/mnt/micron/geo_spatialtrans/S2.3"; R = f"{M}/repro"
def strip(d):
    d = dict(d); d.pop("sec", None); d.pop("elapsed_function_s", None); d.pop("elapsed_total_s", None); return d
n = ok = 0
for f in sorted(glob.glob(f"{R}/bench/grid/*.json")):
    a = strip(json.load(open(f))); b = strip(json.load(open(f.replace(R + "/bench", M + "/bench")))); n += 1; ok += a == b
print(f"banco, verifica: {ok}/{n} file identici (tempi esclusi)")
n = ok = 0
for f in sorted(glob.glob(f"{R}/bench/prep/*.pkl")):
    a = pickle.load(open(f, "rb")); b = pickle.load(open(f.replace(R + "/bench", M + "/bench"), "rb")); n += 1
    same = a["r"] == b["r"] and all(np.array_equal(a["Sr"][k], b["Sr"][k]) if isinstance(a["Sr"][k], np.ndarray) else a["Sr"][k] == b["Sr"][k] for k in a["Sr"]) and a["obs"] == b["obs"]
    ok += same
print(f"banco, preparazione: {ok}/{n} aree con r, densita' e statistiche osservate identici")
n = ok = 0
for f in sorted(glob.glob(f"{R}/panel/main/*.bin")):
    n += 1; ok += hashlib.md5(open(f, "rb").read()).hexdigest() == hashlib.md5(open(f.replace(R + "/panel", M + "/panel"), "rb").read()).hexdigest()
print(f"pannello: {ok}/{n} file dei conteggi con md5 identico")
n = ok = 0
for f in sorted(glob.glob(f"{R}/panel/main/stats/*.json")):
    a = strip(json.load(open(f))); b = strip(json.load(open(f.replace(R + "/panel", M + "/panel")))); n += 1; ok += a == b
print(f"pannello, statistiche: {ok}/{n} file identici (tempi esclusi)")
EOF
cat results/S2.3/S2.3_G3.txt
