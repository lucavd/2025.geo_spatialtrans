#!/usr/bin/env python
"""review/code/07_kclip_scan.py — proprieta' geometriche di K-1/K-2 su tutti i 300 insiemi (6 archetipi x 5 fold x 10 semi):
nuclei multiparte, con buchi, vuoti; K-1 contenuto in K-2 (stesso centro, raggio K-2 >= K-1, stesso territorio);
area(K-1) <= area(K-2); tetto area(K-2) <= phi * a_t (1 + 1e-6); K-2 contro il proprio territorio di S1.4 (dentro, distanza dal bordo).
Uso (root del repo): .venv/bin/python results/S1.4c/review/code/07_kclip_scan.py [ncores]
"""
import sys, glob, os
import numpy as np, pandas as pd, shapely
from multiprocessing import Pool
S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; OUTD = "results/S1.4c/review/code"
gen = pd.read_csv(SC + "/exp/gen.csv")
def one(f1):
    b = os.path.basename(f1)[:-11]; A, roi, fold, s = b.split("_"); k = int(fold[4:])
    K1 = shapely.from_wkb(pd.read_parquet(f1).wkb.to_numpy()); K2 = shapely.from_wkb(pd.read_parquet(f1.replace("_K1.", "_K2.")).wkb.to_numpy())
    cc = pd.read_parquet(f"{SC}/exp/cells_{A}_{roi}.parquet"); phi = float(gen[(gen.archetype == A) & (gen.fold == k)].phi.iloc[0])
    a1, a2 = shapely.area(K1), shapely.area(K2); at = cc.a_t.to_numpy()
    nh = lambda G: int(sum(sum(len(p.interiors) for p in shapely.get_parts(g)) > 0 for g in G))
    cont = shapely.covers(shapely.buffer(K2, 1e-6), K1)
    unc = a2 < phi * at - 1e-3   # tetto non attivo: il tetto usa phi * area del territorio rientrato (a_t - perimetro * 1e-5)
    return dict(set=b, archetype=A, n=len(K1), multipart_K1=int(np.sum(shapely.get_num_geometries(K1) > 1)), multipart_K2=int(np.sum(shapely.get_num_geometries(K2) > 1)),
                holes_K1=nh(K1), holes_K2=nh(K2), empty=int(np.sum(shapely.is_empty(K1)) + np.sum(shapely.is_empty(K2))),
                K1_not_in_K2=int(np.sum(~cont)), area_K1_gt_K2=int(np.sum(a1 > a2 * (1 + 1e-9))), uncapped=int(unc.sum()), K1_not_in_K2_uncapped=int(np.sum(~cont & unc)), area_K1_gt_K2_uncapped=int(np.sum((a1 > a2 * (1 + 1e-6)) & unc)),
                cap_violations=int(np.sum(a2 > phi * at * (1 + 1e-6))), frac_K2_eq_K1=float(np.mean(np.abs(a2 - a1) <= 1e-9 * a2)),
                max_ratio_K2_at=float(np.max(a2 / at)), phi=phi)
if __name__ == "__main__":
    fs = sorted(f for f in glob.glob(f"{S}/gen_in/*_K1.parquet") if "_s00_" not in f)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 20) as p: rows = p.map(one, fs, chunksize=1)
    d = pd.DataFrame(rows); d.to_csv(f"{OUTD}/07_kclip_scan.csv", index=False)
    g = d.groupby("archetype").agg(sets=("set", "size"), multipart_K1=("multipart_K1", "sum"), multipart_K2=("multipart_K2", "sum"), holes_K1=("holes_K1", "sum"),
        holes_K2=("holes_K2", "sum"), empty=("empty", "sum"), K1_not_in_K2=("K1_not_in_K2", "sum"), area_K1_gt_K2=("area_K1_gt_K2", "sum"), uncapped=("uncapped", "sum"), K1_not_in_K2_uncapped=("K1_not_in_K2_uncapped", "sum"), area_K1_gt_K2_uncapped=("area_K1_gt_K2_uncapped", "sum"), n=("n", "sum"),
        cap_violations=("cap_violations", "sum"), max_ratio_K2_at=("max_ratio_K2_at", "max"), phi_max=("phi", "max"))
    g.to_csv(f"{OUTD}/07_kclip_scan_summary.csv"); print(len(d), "insiemi"); print(g.to_string())
