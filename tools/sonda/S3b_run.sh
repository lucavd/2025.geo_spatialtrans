#!/usr/bin/env bash
# tools/sonda/S3b_run.sh - S3b (docs/sonda/REGOLA_S3b.md): aree omologhe di H1 e cancelli G-R, G-O, G-C; se passano: ritaglio H&E, StarDist,
# esportazione e catena di REGOLA_S3 sulle 10 aree di T; tabella (V, B, P di S3 + T nuovo), h con kappa = 1, D. Radice
# /mnt/micron/geo_spatialtrans/S3b_sonda, processi dal file WORKERS (8). Rilanciabile: ogni passo salta cio' che esiste.
set -uo pipefail
cd ~/2025.geo_spatialtrans
export S3_ROOT=/mnt/micron/geo_spatialtrans/S3b_sonda S3_REF=/mnt/micron/geo_spatialtrans/S3_sonda/ref S3_AREAS_H1=$HOME/2025.geo_spatialtrans/results/sonda/S3b/areas_H1_hom.csv S3_NOGUARD=H1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $S3_ROOT/jobs; [ -f $S3_ROOT/WORKERS ] || echo 8 > $S3_ROOT/WORKERS
echo "areas $(date)"; [ -f results/sonda/S3b/gates_areas.json ] || .venv/bin/python -B tools/sonda/S3b_areas_T.py || { echo "AREAS FAIL"; exit 1; }
python3 -c "import json,sys; g=json.load(open('results/sonda/S3b/gates_areas.json')); ok=all(g[k]['pass_'] for k in ('G-R','G-O','G-C')); print('gates', {k: g[k]['pass_'] for k in ('G-R','G-O','G-C','CP-2')}); sys.exit(0 if ok else 3)" || { echo "GATE FAIL: stop (REGOLA_S3b)"; exit 3; }
echo "crop $(date)"; [ -f $S3_ROOT/crop.done ] || { .venv/bin/python -B tools/sonda/S3_export.py crop H1 && touch $S3_ROOT/crop.done; } || { echo "CROP FAIL"; exit 1; }
echo "stardist $(date)"; [ -f $S3_ROOT/stardist.done ] || { TF_NUM_INTEROP_THREADS=4 TF_NUM_INTRAOP_THREADS=4 bash tools/py_stardist.sh tools/sonda/S3_stardist.py all && touch $S3_ROOT/stardist.done; } || { echo "STARDIST FAIL"; exit 1; }
echo "export $(date)"; [ -f $S3_ROOT/export.done ] || { .venv/bin/python -B tools/sonda/S3_export.py export H1 && touch $S3_ROOT/export.done; } || { echo "EXPORT FAIL"; exit 1; }
python3 - <<'EOF'
import os, csv
R = os.environ["S3_ROOT"]; A = list(csv.DictReader(open(os.environ["S3_AREAS_H1"])))
with open(f"{R}/jobs/T.tsv", "w") as f:
    for t in ("tumour", "stroma"):
        for a in sorted([x for x in A if x["tissue"] == t], key=lambda x: int(x["k"])):
            u = f"H1_{t}_{a['k']}"; f.write("\t".join(("T", u, f"{R}/runs/T/{u}", t, "real", "H1", "-")) + "\n")
EOF
echo "chain T start $(date)"; bash tools/sonda/S3_chain.sh $S3_ROOT/jobs/T.tsv $(cat $S3_ROOT/WORKERS); echo "chain T end $(date)"
python3 tools/sonda/S3b_table.py && .venv/bin/python -B tools/sonda/S3b_h.py || { echo "TABLE/H FAIL"; exit 1; }
.venv/bin/python -B tools/sonda/S3_D.py --table results/sonda/S3b/area_values.csv --h results/sonda/S3b/h.csv --out results/sonda/S3b/D.json --per results/sonda/S3b/per_measure.csv --misura results/misura.json --provv results/sonda/S3b/misura_provvisoria.json || { echo "D FAIL"; exit 1; }
echo "run end $(date)"
