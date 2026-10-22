#!/usr/bin/env bash
# R6 gate G-3: rerun clustering + labelling (exp, UMI>=20, seed 20261022) into a scratch dir and compare labels with the primary run.
set -euo pipefail
cd ~/2025.geo_spatialtrans
L=/mnt/micron/geo_spatialtrans/R6; export OMP_NUM_THREADS=2
mkdir -p $L/repro/labels; rm -f $L/repro/labels/*
sed "s#OUT = f\"{M}/labels\"#OUT = f\"{M}/repro/labels\"#" tools/R6_cluster.py > $L/repro/R6_cluster_repro.py
for a in A1 A2 A3 A4 A5 A6; do
  .venv/bin/python $L/repro/R6_cluster_repro.py $a exp 20 20261022 > $L/repro/cluster_$a.log 2>&1 &
done; wait
.venv/bin/python - <<'EOF'
import pandas as pd, glob, hashlib
L="/mnt/micron/geo_spatialtrans/R6"; rows=[]
for a in ["A1","A2","A3","A4","A5","A6"]:
    p=pd.read_parquet(f"{L}/labels/{a}_exp_umi20_seed20261022_nuclei.parquet"); q=pd.read_parquet(f"{L}/repro/labels/{a}_exp_umi20_seed20261022_nuclei.parquet")
    same_cl=(p.cluster.fillna(-1).astype(str).values==q.cluster.fillna(-1).astype(str).values).all(); same_ty=(p.cell_type.values==q.cell_type.values).all()
    h1=hashlib.md5(pd.util.hash_pandas_object(p[["roi_id","label","cluster","cell_type"]].astype(str)).values).hexdigest(); h2=hashlib.md5(pd.util.hash_pandas_object(q[["roi_id","label","cluster","cell_type"]].astype(str)).values).hexdigest()
    rows.append(dict(archetype=a,n=len(p),same_clusters=same_cl,same_types=same_ty,md5_primary=h1,md5_rerun=h2,identical=h1==h2)); print(rows[-1])
pd.DataFrame(rows).to_csv("results/R6/R6_G3_repro.csv",index=False)
EOF
