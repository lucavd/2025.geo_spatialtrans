#!/usr/bin/env bash
# tools/F3.1_repro.sh - F3.1, cancello G-3: un'area per tessuto, simulata (k = 0) e vera, rifatta da zero in un'altra cartella con 3 processi
# invece di 6; confronto dei JSON (b2c_label, b2c_label_shift, rctd; vero anche b2c_label_dimg) senza i tempi. Uso: tools/F3.1_repro.sh
set -uo pipefail
R=/mnt/micron/geo_spatialtrans/F3.1/main; P=/mnt/micron/geo_spatialtrans/F3.1/repro; mkdir -p $P
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
TAGS=A1_r1,A2_r1,A3_r1,A4_f1,A5_r1,A6_r1
.venv/bin/python tools/F3.1_reference.py $P/ref 20261113 A1 A2 A3 A4 A5 A6 > $P/ref.log 2>&1
Rscript --vanilla tools/F3.1_sim.R $P/sim_rds 20261114 $TAGS 3 0 > $P/sim.log 2>&1
( .venv/bin/python tools/F3.1_export_real.py $P/real A1_r1 > $P/x1.log 2>&1 ) & ( .venv/bin/python tools/F3.1_export_real.py $P/real A2_r1 A3_r1 > $P/x2.log 2>&1 ) &
( .venv/bin/python tools/F3.1_export_real.py $P/real A4_f1 A5_r1 A6_r1 > $P/x3.log 2>&1 ) & wait
for T in ${TAGS//,/ }; do
  A=${T%%_*}; SIDE=$(.venv/bin/python -c "import pandas as pd; r=pd.read_csv('results/R2/R2_rois_checked.csv'); r=r[(r.archetype+'_'+r.roi_id)=='$T'].iloc[0]; print(r.side_px*r.um_per_px)")
  Rscript --vanilla tools/F3.1_export_sim.R $P/sim_rds/${T}_k0.rds $P/sim/${T}_k0 $SIDE > /dev/null 2>&1; .venv/bin/python tools/F3.1_reg_sim.py $P/sim/${T}_k0 > /dev/null
done
job() { d=$1; s=$2; t=$(basename $d); A=${t%%_*}
  for L in label label_shift $( [ $s = real ] && echo label_dimg ); do .venv-b2c/bin/python tools/F3.1_b2c.py $d $P/ref $A $s $L > $d/b2c_$L.log 2>&1; done
  Rscript --vanilla tools/F3.1_rctd.R $d $P/ref $A $s > $d/rctd.log 2>&1; }
export -f job; export P
{ ls -d $P/sim/* | sed 's/$/ sim/'; ls -d $P/real/* | sed 's/$/ real/'; } | xargs -P 3 -L 1 bash -c 'job $0 $1'
python3 - <<'EOF'
import json, glob, os, filecmp
R = "/mnt/micron/geo_spatialtrans/F3.1/main"; P = "/mnt/micron/geo_spatialtrans/F3.1/repro"; bad = 0; n = 0
for side in ["sim", "real"]:
    for d in sorted(glob.glob(f"{P}/{side}/*")):
        for f in ["b2c_label.json", "b2c_label_shift.json", "rctd.json"] + (["b2c_label_dimg.json"] if side == "real" else []):
            a = json.load(open(f"{d}/{f}")); b = json.load(open(f"{R}/{side}/{os.path.basename(d)}/{f}"))
            for x in (a, b): x.pop("sec", None); x.pop("UMI_min_sigma", None)   # tempo; soglia registrata solo nelle esecuzioni non "reuse"
            n += 1; ok = a == b
            if not ok:
                bad += 1; diff = {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)}; print("DIFF", side, os.path.basename(d), f, list(diff.items())[:4])
for f in sorted(glob.glob(f"{P}/ref/*")):
    n += 1; same = filecmp.cmp(f, f"{R}/ref/{os.path.basename(f)}", shallow=False); bad += (not same)
    if not same: print("DIFF ref", os.path.basename(f))
print(f"G-3: {n - bad} identici su {n}")
EOF
