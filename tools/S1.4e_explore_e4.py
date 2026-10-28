#!/usr/bin/env python
"""tools/S1.4e_explore_e4.py - S1.4e, esplorazione dichiarata (docs/steps/S1.4e.md), E4: contatto coppia per coppia dei nuclei
simulati nati dai punti di Luca (controprova CP-1 di R3b: linfonodo e cuore, nucleo ritagliato e riportato all'area giusta, 10 semi)
sulle coppie giudicate da Luca in R2c. Contatto simulato = distanza euclidea fra i due poligoni <= 2 px (0.547 um): approssima
il contatto di Chebyshev <= 2 px del raster (tools/S1.4c_raster.py), dichiarato. Nucleo del punto = poligono che contiene il punto.
Uscite: results/S1.4e/explore/S1.4e_E4_pairs.csv, S1.4e_E4_summary.csv."""
import json, glob, re
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
R2B = ROOT / "results/R2b"; R2C = ROOT / "results/R2c"; OUT = ROOT / "results/S1.4e/explore"
GEN = Path("/mnt/micron/geo_spatialtrans/R3b/gen_in_cp1")
BINS = [0, 5, 6, 7, 8, 10, np.inf]; LAB = ["<5", "5-6", "6-7", "7-8", "8-10", ">=10"]

def main():
    wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"]
    it = pd.read_csv(R2C / "R2c_items_Luca.csv"); a = pd.DataFrame(json.load(open(R2C / "annotations/R2c_answers_Luca.json"))["answers"])
    J = it.merge(a, on="item", validate="1:1"); J = J[(J.kind == "pair") & J.win_id.str[:2].isin(["A4", "A6"])]
    rows = []
    for wid, g in J.groupby("win_id"):
        w = wins.loc[wid]; upp = float(w.um_per_px); P = pts[pts.win_id == wid].set_index("idx")
        files = sorted(glob.glob(str(GEN / f"{w.archetype}_{w.roi_id}_fold*_s*_K2.parquet")))
        assert len(files) == 10, (wid, len(files))
        for f in files:
            sd = re.search(r"_s(\d+)_K2", f).group(1)
            geo = shapely.from_wkb(pd.read_parquet(f).wkb.to_numpy()); tree = shapely.STRtree(geo)
            q = tree.query(shapely.points(P[["x_um", "y_um"]].to_numpy()), predicate="within")
            own = {}
            for pi_, gi in zip(*q):
                own.setdefault(int(P.index[pi_]), []).append(int(gi))
            for r in g.itertuples():
                ki, kj = own.get(int(r.i), []), own.get(int(r.j), [])
                if len(ki) != 1 or len(kj) != 1:
                    rows.append(dict(win_id=wid, seed=sd, i=r.i, j=r.j, d_um=r.d_um, answer=r.answer, gap_um=np.nan, sim_contact=np.nan)); continue
                gap = float(shapely.distance(geo[ki[0]], geo[kj[0]]))
                rows.append(dict(win_id=wid, seed=sd, i=r.i, j=r.j, d_um=r.d_um, answer=r.answer, gap_um=gap, sim_contact=float(gap <= 2 * upp)))
    E = pd.DataFrame(rows); E["archetype"] = E.win_id.str[:2]; E["bin"] = pd.cut(E.d_um, BINS, right=False, labels=LAB)
    E.to_csv(OUT / "S1.4e_E4_pairs.csv", index=False)
    print("coppie senza nucleo unico:", int(E.sim_contact.isna().sum()), "su", len(E))
    S = []
    for (A, b), x in E.groupby(["archetype", "bin"], observed=True):
        jd = x[x.answer != "unsure"].drop_duplicates(["win_id", "i", "j"])
        S.append(dict(archetype=A, bin=b, n_pairs=len(jd), p_touch_manual=float((jd.answer == "touch").mean()),
                      p_contact_sim=float(x.sim_contact.mean()), gap_med_sim=float(x.gap_um.median())))
    S = pd.DataFrame(S); S.to_csv(OUT / "S1.4e_E4_summary.csv", index=False); print(S.round(3).to_string())

if __name__ == "__main__":
    main()
