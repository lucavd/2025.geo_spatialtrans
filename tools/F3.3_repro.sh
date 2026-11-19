#!/usr/bin/env bash
# tools/F3.3_repro.sh - F3.3, cancello G-3: un'area per sezione di verifica rifatta da zero in un'altra cartella (F33_ROOT) con 2 processi
# invece di 4; confronto dei JSON senza i tempi, delle maschere e delle regioni. Uso: bash tools/F3.3_repro.sh
set -uo pipefail; cd ~/2025.geo_spatialtrans
M=/mnt/micron/geo_spatialtrans; export F33_ROOT=$M/F3.3/repro; REF=$M/F3.1/main/ref; mkdir -p $F33_ROOT/runs $F33_ROOT/measures
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
NON=A1_mucosa_1,A4_follicolo_1,A5_corteccia_1,A6_miocardio_1,MouseBrain_FFPE_SR401_corteccia_1
CRC=A23_tumore_1,CRC_P1_SR300_tumore_2,CRC_P5_SR300_tumore_1,CRC_6p5mm_SR410_tumore_2
.venv/bin/python tools/F3.3_segment.py $NON,$CRC > $F33_ROOT/log_seg.log 2>&1
tools/py_stardist.sh tools/F3.3_stardist.py A4_follicolo_1,A6_miocardio_1 > $F33_ROOT/log_sd.log 2>&1
.venv/bin/python tools/F3.3_consensus.py A4_follicolo_1,A6_miocardio_1 > $F33_ROOT/log_co.log 2>&1
Rscript --vanilla tools/F3.3_regions.R $NON,$CRC > $F33_ROOT/log_reg.log 2>&1
.venv/bin/python tools/F3.3_export.py panel $NON > $F33_ROOT/log_exp.log 2>&1
.venv/bin/python tools/F3.3_export.py common $CRC >> $F33_ROOT/log_exp.log 2>&1
one() { local d=$1; local A=$(python3 -c "import json;print(json.load(open('$d/export.json'))['archetype'])")
  .venv-b2c/bin/python tools/F3.1_b2c.py $d $REF $A real label > $d/b2c_label.log 2>&1; Rscript --vanilla tools/F3.1_rctd.R $d $REF $A real > $d/rctd.log 2>&1
  F33_CACHE=$F33_ROOT/cache .venv-b2c/bin/python tools/F3.3_holdout_measures.py $d $F33_ROOT/measures/$(basename $d).json > /dev/null 2>&1; echo "done $d"; }
export -f one; export REF F33_ROOT
ls -d $F33_ROOT/runs/*/ | sed 's#/$##' | xargs -P 2 -L 1 bash -c 'one $0' > $F33_ROOT/log_runs.log 2>&1
python3 - <<'EOF'
import json, glob, os, hashlib, numpy as np
A, B = "/mnt/micron/geo_spatialtrans/F3.3", "/mnt/micron/geo_spatialtrans/F3.3/repro"; rows = []
def js(p): j = json.load(open(p)); [j.pop(k, None) for k in ("sec", "rctd_sec")]; return j
def md5(p): return hashlib.md5(np.load(p)["labels"].tobytes()).hexdigest()
for d in sorted(glob.glob(f"{B}/runs/*")):
    run = os.path.basename(d); area = run.replace("_cg", ""); ok = {}
    for f in ("export.json", "b2c_label.json", "rctd.json"): ok[f] = js(f"{A}/runs/{run}/{f}") == js(f"{d}/{f}")
    s1, s2 = js(f"{A}/seg/{area}/seg.json"), js(f"{B}/seg/{area}/seg.json"); ok["seg.json"] = s1 == s2
    for npz in glob.glob(f"{B}/seg/{area}/*_native.npz"): ok[os.path.basename(npz)] = md5(npz) == md5(f"{A}/seg/{area}/{os.path.basename(npz)}")
    ok["regions"] = open(f"{A}/seg/{area}/regions.tsv").read() == open(f"{B}/seg/{area}/regions.tsv").read()
    m1, m2 = json.load(open(f"{A}/measures/{run}.json")), json.load(open(f"{B}/measures/{run}.json")); ok["measures"] = m1 == m2
    rows.append((run, all(ok.values()), {k: v for k, v in ok.items() if not v}))
for r in rows: print("PASS" if r[1] else "FAIL", r[0], r[2])
print(sum(r[1] for r in rows), "su", len(rows), "identiche")
json.dump([{"run": r[0], "identical": r[1], "diff": list(r[2])} for r in rows], open("results/F3.3/F3.3_G3.json", "w"))
EOF
