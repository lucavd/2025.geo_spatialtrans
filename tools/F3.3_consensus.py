"""tools/F3.3_consensus.py - F3.3: maschera di consenso di R3b (nuclei di Space Ranger con la forma di StarDist,
tools/R3b_consensus.build_consensus senza mutanti) per le aree 'consensus', con la tabella dei nuclei e il flag keep di R2
(R2_common.nuclei_table con le maschere di tessuto e bolle dell'area). Uso: .venv/bin/python tools/F3.3_consensus.py <aree|IN|holdout>"""
import sys, os, json, numpy as np, pandas as pd, tifffile
from pathlib import Path
os.environ["R3B_MUTANT"] = "none"; sys.path.insert(0, "tools")
import R2_common as C, R3b_consensus as R3
OUT = Path("/mnt/micron/geo_spatialtrans/F3.3/seg")
A = pd.read_csv("results/F3.3/F3.3_areas.csv"); A = A[A.seg == "consensus"]; sel = sys.argv[1]
A = A[A.level == "IN"] if sel == "IN" else A[A.level != "IN"] if sel == "holdout" else A[A.area.isin(sel.split(","))]
for _, r in A.iterrows():
    o = OUT / r.area
    sr = np.load(o / "spaceranger_native.npz")["labels"]; sd = np.load(o / "stardist_he_native.npz")["labels"]
    lab, tab = R3.build_consensus(sr, sd)
    rgb = tifffile.imread(o / "roi.tif"); tm, bm = C.masks_for(rgb, float(r.upp), 1)
    nuc = C.nuclei_table(lab.astype(np.uint32), float(r.upp), tm, bm, r.A, r.area, "consensus", 1)
    np.savez_compressed(o / "consensus_native.npz", labels=lab.astype(np.uint32)); nuc.to_parquet(o / "consensus_native.parquet", index=False)
    tab.to_parquet(o / "consensus_map.parquet", index=False)
    j = json.load(open(o / "seg.json")); j.update(n_consensus=int(lab.max()), n_consensus_keep=int(nuc.keep.sum())); json.dump(j, open(o / "seg.json", "w"))
    print(r.area, j, flush=True)
