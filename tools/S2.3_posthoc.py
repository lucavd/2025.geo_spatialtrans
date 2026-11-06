"""tools/S2.3_posthoc.py — S2.3, analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti pre-registrati).
PH-3: i «lumi» delle regioni dello Step 1 sono vuoti? frazione dei quadratini delle fasce 2-8 um dentro i lumi (stessa selezione di P-9) il cui
      centro allineato cade nel tessuto secondo la maschera H&E di R5 (densita' ottica, include il tessuto pallido; BL-085).
PH-4: il calo radiale (P-7) sulla geometria vera: rapporto R sul banco (nuclei veri, campo di R-1, 3 repliche nuove) contro il vero, sui
      quadratini in tessuto di Space Ranger, distanza dal nucleo vero piu' vicino (maschera allineata).
Uso: .venv/bin/python tools/S2.3_posthoc.py <benchdir> -> results/S2.3/S2.3_PH3_lumens.csv, S2.3_PH4_radial.csv"""
import sys, os, json, pickle
import numpy as np, pandas as pd, shapely
from multiprocessing import Pool
sys.path.insert(0, "tools")
import importlib.util
_s = importlib.util.spec_from_file_location("b", "tools/S2.3_bench.py"); BN = importlib.util.module_from_spec(_s); _s.loader.exec_module(BN)
RL, NC, S = BN.RL, BN.NC, BN.S
BD = sys.argv[1]
CH = pd.read_csv("results/S2.3/S2.3_R1_choice.csv").set_index("archetype")

def one(Rw):
    tag = f"{Rw.archetype}_{Rw.roi_id}"; A = Rw.archetype; upp = Rw.um_per_px
    P = pickle.load(open(f"{BD}/prep/{tag}.pkl", "rb"))
    Dd = NC.roi_data(Rw); Mp = S.nuclear_maps(Dd["mask"], upp)
    # PH-4: R sul banco
    rr, cc = S.aligned_px(P["pos"], P["di"], upp)
    dn = S.lookup(Mp["d_out"], rr, cc, np.inf)
    obs = S.profile_ratio(dn, P["u"])
    s2, L = float(CH.loc[A, "s2"]), float(CH.loc[A, "L"])
    sims = []
    for k in range(3):
        rng_f = np.random.default_rng([S.SEED, P["idx"], 900, k, 1]); U = np.random.default_rng([S.SEED, P["idx"], 900, k, 2]).uniform(size=(1, len(P["u"])))
        y = BN.simulate(P, s2, L, rng_f, U)[0]; sims.append(S.profile_ratio(dn, y))
    r4 = {"tag": tag, "archetype": A, "R_obs": obs, "R_sim": float(np.mean(sims)), "ratio": float(np.mean(sims)) / obs}
    # PH-3: lumi dello Step 1 contro la maschera H&E di R5
    rowsA, colsA, posA = RL.box_bins(Rw, False)
    posA_al = posA - P["di"][None, :]
    side = float(Rw.side_px * upp); R = RL.region_geom(tag)
    tis = RL.tissue_lookup(P["fm"], posA_al)
    holes = []
    for poly in getattr(R, "geoms", [R]):
        for ring in poly.interiors:
            h = shapely.Polygon(ring)
            if h.area < RL.LUMEN_MIN_UM2: continue
            x0, y0, x1, y1 = h.bounds
            if x0 <= 0 or y0 <= 0 or x1 >= side or y1 >= side: continue
            holes.append(h)
    x, y = posA_al[:, 1], posA_al[:, 0]; nb = 0; nt = 0
    for h in holes:
        x0, y0, x1, y1 = h.bounds
        near = np.flatnonzero((x >= x0) & (x <= x1) & (y >= y0) & (y <= y1))
        ins = shapely.contains_xy(h, x[near], y[near]); d = shapely.distance(shapely.points(x[near], y[near]), h.exterior)
        m = ins & (d >= 2) & (d < 8); nb += int(m.sum()); nt += int(tis[near][m].sum())
    r3 = {"tag": tag, "archetype": A, "n_lumens": len(holes), "n_band_bins": nb, "n_band_bins_HE_tissue": nt, "frac_HE_tissue": nt / nb if nb else np.nan}
    return r3, r4

if __name__ == "__main__":
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    with Pool(6) as p: out = p.map(one, [Rw for _, Rw in rois.iterrows()])
    T3 = pd.DataFrame([o[0] for o in out]); T4 = pd.DataFrame([o[1] for o in out])
    T3.to_csv("results/S2.3/S2.3_PH3_lumens.csv", index=False); T4.to_csv("results/S2.3/S2.3_PH4_radial.csv", index=False)
    g3 = T3.groupby("archetype")[["n_lumens", "n_band_bins", "n_band_bins_HE_tissue"]].sum(); g3["frac"] = g3.n_band_bins_HE_tissue / g3.n_band_bins
    print(g3)
    g4 = T4.groupby("archetype").ratio.agg(["mean", "std", "count"]); g4["lo2se"] = g4["mean"] - 2 * g4["std"] / np.sqrt(g4["count"]); g4["hi2se"] = g4["mean"] + 2 * g4["std"] / np.sqrt(g4["count"])
    print(g4.round(4)); print(T4.groupby("archetype")[["R_obs", "R_sim"]].mean().round(4))
