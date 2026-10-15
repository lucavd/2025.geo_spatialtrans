
import sys; sys.argv = ["x"]
exec(open("results/S1.4b/review/code/adv_arbiter.py").read().split("if __name__")[0])
import pandas as pd, numpy as np, shapely
f = "/tmp/code_adv/C2_near_exact.parquet"; d = pd.read_parquet(f); reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read())
cells = shapely.from_wkb(d.wkb.to_numpy()); x, y, w = d.x.to_numpy(), d.y.to_numpy(), d.w.to_numpy()
bb = reg.bounds; box = (bb[0]-100, bb[1]-100, bb[2]+100, bb[3]+100)
bf = shapely.intersection(np.array(bf_cells(x, y, w, box), dtype=object), reg)
aR, aB = shapely.area(cells), shapely.area(bf)
offs = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 0]; t = 150
for k, o in enumerate(offs, 1):
    for eqw in (True, False):
        i = 10 * k + int(eqw) - 1  # R 1-based index 10k+eqw -> 0-based
        print(f"off={o:g} eqw={eqw}: base i={i} aR={aR[i]:.6g} aB={aB[i]:.6g} | twin j={t} aR={aR[t]:.6g} aB={aB[t]:.6g}  hiddenR={d.hidden[i]},{d.hidden[t]}")
        t += 1
