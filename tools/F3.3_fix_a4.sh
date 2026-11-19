#!/usr/bin/env bash
# tools/F3.3_fix_a4.sh - F3.3, correzione dopo G-3: StarDist deterministico (F33_SD_DET=1) sulle 2 aree di verifica di 500 um del linfonodo,
# rifatte da StarDist in poi nella cartella principale e in quella di riproducibilita'; poi confronto.
set -uo pipefail; cd ~/2025.geo_spatialtrans; M=/mnt/micron/geo_spatialtrans; REF=$M/F3.1/main/ref
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 F33_SD_DET=1
for ROOT in $M/F3.3 $M/F3.3/repro; do
  export F33_ROOT=$ROOT; L=A4_follicolo_1,A4_follicolo_2; [ $ROOT = $M/F3.3/repro ] && L=A4_follicolo_1
  for a in ${L//,/ }; do rm -f $ROOT/seg/$a/stardist_he_native.npz $ROOT/seg/$a/consensus_* ; rm -rf $ROOT/runs/$a; rm -f $ROOT/measures/$a.json $ROOT/cache/runs_${a}_exp.npy $M/F3.3/cache/runs_${a}_exp.npy; done
  tools/py_stardist.sh tools/F3.3_stardist.py $L > $ROOT/log_fix_sd.log 2>&1
  .venv/bin/python tools/F3.3_consensus.py $L > $ROOT/log_fix_co.log 2>&1
  .venv/bin/python tools/F3.3_export.py panel $L > $ROOT/log_fix_exp.log 2>&1
  for a in ${L//,/ }; do d=$ROOT/runs/$a
    .venv-b2c/bin/python tools/F3.1_b2c.py $d $REF A4 real label > $d/b2c_label.log 2>&1; Rscript --vanilla tools/F3.1_rctd.R $d $REF A4 real > $d/rctd.log 2>&1
    F33_CACHE=$ROOT/cache .venv-b2c/bin/python tools/F3.3_holdout_measures.py $d $ROOT/measures/$a.json > /dev/null 2>&1; echo "done $d"; done
done
.venv/bin/python - <<'EOF'
import json, hashlib, numpy as np
A, B = "/mnt/micron/geo_spatialtrans/F3.3", "/mnt/micron/geo_spatialtrans/F3.3/repro"; a = "A4_follicolo_1"
js = lambda p: {k: v for k, v in json.load(open(p)).items() if k not in ("sec", "rctd_sec")}
ok = {f: js(f"{A}/runs/{a}/{f}") == js(f"{B}/runs/{a}/{f}") for f in ("export.json", "b2c_label.json", "rctd.json")}
for f in ("stardist_he_native.npz", "consensus_native.npz"): ok[f] = np.array_equal(np.load(f"{A}/seg/{a}/{f}")["labels"], np.load(f"{B}/seg/{a}/{f}")["labels"])
ok["measures"] = json.load(open(f"{A}/measures/{a}.json")) == json.load(open(f"{B}/measures/{a}.json")); print(ok, all(ok.values()))
d = json.load(open("results/F3.3/F3.3_G3.json"))
for x in d:
    if x["run"] == a: x.update(identical=all(ok.values()), diff=[k for k, v in ok.items() if not v], note="rifatta con StarDist deterministico (F33_SD_DET=1)")
json.dump(d, open("results/F3.3/F3.3_G3.json", "w")); print(sum(x["identical"] for x in d), "su", len(d))
EOF
