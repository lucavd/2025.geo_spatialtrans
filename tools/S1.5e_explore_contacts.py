#!/usr/bin/env python
"""tools/S1.5e_explore_contacts.py - S1.5e, esplorazione dichiarata PRIMA della pre-registrazione (docs/steps/S1.5e.md).
(E3b) Da dove viene l'eccesso di contatti dei bracci con i vuoti? Per nucleo focale: distanza dal vicino piu' vicino (centro-centro) e
contatto. Vero: focali di R2c (60 per tessuto), contatto dalle risposte di Luca (variante primaria di tools/R2c_analysis.focal_s1),
distanza dal vicino fra i punti manuali di R2b della stessa finestra. Simulato: focali come tools/S1.5d_windows.py (distanza >= d_max
dai lati e dalle esclusioni), contatto = n_c2 > 0 del raster, distanza dal vicino fra i centri del raster; bracci U, N, RV, 10 semi.
Per classe di distanza: frazione dei focali e frazione a contatto; scomposizione s1_cf = somma_b P_vero(b) * s1_sim(b) (contatti che il
modello del nucleo darebbe con le distanze vere). Uscite: results/S1.5e/S1.5e_explore_contacts_*.csv.
Uso: .venv/bin/python tools/S1.5e_explore_contacts.py"""
import importlib.util, json, sys
import numpy as np, pandas as pd, shapely
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
K = mod("r2c_cand", ROOT / "tools/R2c_candidates.py"); AN = mod("r2c_an", ROOT / "tools/R2c_analysis.py")
R2B = ROOT / "results/R2b"; R2C = ROOT / "results/R2c"; RES = ROOT / "results/S1.5e"; RAS = Path("/mnt/micron/geo_spatialtrans/S1.5d/raster")
ARMS = ["U", "N", "RV"]; SEEDS = [f"s{k:02d}" for k in range(1, 11)]
BINS = [0, 5, 6, 7, 8, 10, np.inf]

def nn_dist(P, Q):
    d, _ = cKDTree(P).query(Q, k=2); return d[:, 1]

def binned(nn, y):
    b = np.digitize(nn, BINS[1:-1]); out = {}
    for k in range(len(BINS) - 1):
        m = b == k; out[f"frac_b{k}"] = float(m.mean()); out[f"s1_b{k}"] = float(np.nanmean(y[m])) if m.any() else np.nan; out[f"n_b{k}"] = int(m.sum())
    return out

def main():
    RES.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(R2C / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]; WJ = {w["win_id"]: w for w in wjson}
    ans = AN.load_answers("Luca"); ans = ans[ans.kind == "pair"][["win_id", "i", "j", "answer"]]   # come R2c_analysis.s1_table
    foc = pd.read_csv(R2C / "R2c_focal.csv"); foc = foc[foc.focal].copy()
    foc["s1"] = AN.focal_s1(ans, foc, "primary")
    per = []
    for wr in wins.itertuples():
        w = WJ[wr.win_id]; upp = w["um_per_px"]; mg = w["margin_px"]; p = np.array(w["points"], float).reshape(-1, 2)
        Pm = np.c_[w["x0_um"] + (p[:, 0] - mg) * upp, w["y0_um"] + (p[:, 1] - mg) * upp]
        f = foc[foc.win_id == wr.win_id]
        if len(f):
            nn = nn_dist(Pm, Pm[f.idx.values.astype(int)])
            for i, (d, y) in enumerate(zip(nn, f.s1.values)):
                per.append(dict(source="manual", archetype=wr.archetype, win_id=wr.win_id, seed="", nn_um=d, contact=y))
        dmax = float(dm.d_max_um[wr.archetype]); x0, y0 = float(wr.x0_um), float(wr.y0_um); s = int(wr.side_px) * float(wr.um_per_px)
        ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
        for arm in ARMS:
            for sd in SEEDS:
                t = pd.read_parquet(RAS / arm / f"{wr.archetype}_{wr.roi_id}_{sd}.parquet", columns=["cx", "cy", "n_c2"])
                P = np.c_[t.cx.values, t.cy.values]
                de = np.minimum.reduce([P[:, 0] - x0, x0 + s - P[:, 0], P[:, 1] - y0, y0 + s - P[:, 1]]); keep = de >= dmax
                if u is not None: keep &= shapely.distance(shapely.points(P), u) >= dmax
                if not keep.any(): continue
                nn = nn_dist(P, P[keep]); y = (t.n_c2.values[keep] > 0).astype(float)
                per += [dict(source=arm, archetype=wr.archetype, win_id=wr.win_id, seed=sd, nn_um=a, contact=b) for a, b in zip(nn, y)]
    d = pd.DataFrame(per); d.to_csv(RES / "S1.5e_explore_contacts_focal.csv", index=False)
    rows = []
    for (A, src), g in d.groupby(["archetype", "source"]):
        rows.append(dict(archetype=A, source=src, n=len(g), s1=float(np.nanmean(g.contact)), nn_median=float(g.nn_um.median()),
                         **binned(g.nn_um.values, g.contact.values)))
    S = pd.DataFrame(rows)
    nb = len(BINS) - 1
    for i, r in S.iterrows():
        m = S[(S.archetype == r.archetype) & (S.source == "manual")]
        if len(m):
            m = m.iloc[0]; s1b = np.array([r[f"s1_b{k}"] for k in range(nb)]); pm = np.array([m[f"frac_b{k}"] for k in range(nb)])
            ok = np.isfinite(s1b); S.loc[i, "s1_cf_manual_nn"] = float(np.sum(pm[ok] * s1b[ok]) / pm[ok].sum())
    S.to_csv(RES / "S1.5e_explore_contacts_summary.csv", index=False)
    chk = pd.read_csv(R2C / "R2c_s1.csv"); chk = chk[chk.set == "Luca"].set_index("archetype").s1_primary
    for A in S.archetype.unique():
        v = S[(S.archetype == A) & (S.source == "manual")].s1.iloc[0]
        print(A, "s1 manuale ricostruito", round(v, 6), "R2c", round(float(chk[A]), 6), "uguale" if abs(v - chk[A]) < 1e-9 else "DIVERSO")

if __name__ == "__main__":
    main()
