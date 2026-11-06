"""tools/S2.3_bench.py — S2.3: banco della taratura del campo di variabilita' fine (geometria vera, modello completo) e esplorazione E-1.
Definizioni: docs/steps/S2.3.md («Il banco della taratura», regola R-1, cancello G-4 c). Prototipo di S2.0 (S2.0_disp.py, PH-2) con i
fattori di cattura della sezione, la regola C2 con i parametri corretti (run1_wfix) e il campo F sulla sottogriglia di 0.5 um.
Uso:
  .venv/bin/python tools/S2.3_bench.py prep <outdir> [params.json] [aligned 0|1]   r invertito e densita' rifinite per area (una volta)
  .venv/bin/python tools/S2.3_bench.py grid <outdir> <combos.json> <rep0> <nrep>   simulazioni per area x combinazione -> <outdir>/grid/
  .venv/bin/python tools/S2.3_bench.py check_s20 <outdir>                          G-4 (c): ritrova S2.0 PH-2 (parametri di run1, s2 = 0)
  .venv/bin/python tools/S2.3_bench.py e1 <outdir>                                 E-1: scala del campo su tessuto sintetico omogeneo
Processi dal file <outdir>/WORKERS."""
import sys, os, json, time, traceback, pickle
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import importlib.util
_spec = importlib.util.spec_from_file_location("s23real", "tools/S2.3_real.py"); RL = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(RL)
NC, C, S0, S = RL.NC, RL.C, RL.S0, RL.S
from scipy import stats as _stats
NREP_S20 = 5
COARSE_L = [0.5, 0.75, 1.0, 1.5, 2.0]
COARSE_S2 = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

def rois():
    return pd.read_csv("results/R2/R2_rois_checked.csv")

def prep_one(args):
    Rw, OUT, PJF, aligned = args
    tag = f"{Rw.archetype}_{Rw.roi_id}"; fo = f"{OUT}/prep/{tag}.pkl"
    if os.path.exists(fo): return fo
    t0 = time.time()
    PJ = json.load(open(PJF)); prm = dict(PJ[NC.KEY[Rw.archetype]]["C2"])
    Dd, IX, IXf, Dobs, am, acell, S_, U5, base = NC.prep(tag, Rw)
    rows, cols, pos = RL.box_bins(Rw, True); assert len(rows) == len(Dd["u"])
    c2 = RL.stripes_c(Rw.archetype, rows, cols); Dd = dict(Dd); Dd["c"] = c2           # come S2.0 PH-2
    Ui = np.random.default_rng([S0.SEED, int(Rw.name)]).uniform(size=U5.shape)
    di = NC.d_image(Dd, prm)
    fm = RL.roi_tissue(Rw, Dd)
    blur = S0.make_blur("C2", prm, Dd["pos"], di, fm)
    r, Sr = S0.invert_r_refined_b(base["obs_star"], S_, blur, Dd["c"], Ui, IX, Dd["L"], Dd["keep"], Dobs, acell, am)
    r_src = "inverted"
    if not np.isfinite(r):
        r = base["obs_star"]; r_src = "observed_DEC034"
        Sr = S0.refine_b(S_, r, blur, Dd["c"], Dd["L"], Dd["keep"], Dobs, acell, am)
    upp = Dd["upp"]
    if aligned:
        rr, cc = S.aligned_px(Dd["pos"], di, upp); tis = RL.tissue_lookup(fm, Dd["pos"] - di[None, :])
    else:
        rr, cc = Dd["rr"], Dd["cc"]; tis = RL.tissue_lookup(fm, Dd["pos"])
    Mp = S.nuclear_maps(Dd["mask"], upp)
    G = S.groups_from_maps(Mp, rr, cc, Dd["keep"], tis, rows, cols)
    lab0 = S0.deep_labels(Dd["mask"], Dd["rr"], Dd["cc"], upp, 1.0)                    # gruppi di S2.0 (per G-4 c)
    obs = {k: S.group_stats(Dd["u"], g, c2, rows, cols) for k, g in G.items()}
    out = {"tag": tag, "archetype": Rw.archetype, "idx": int(Rw.name), "r": float(r), "r_source": r_src, "Sr": Sr, "prm": prm, "di": di,
           "pos": Dd["pos"], "c": c2, "rows": rows, "cols": cols, "u": Dd["u"], "groups": G, "lab_s20": lab0, "keep": Dd["keep"],
           "fm": fm, "obs": obs, "sec": round(time.time() - t0, 1), "aligned": bool(aligned), "params_file": PJF}
    pickle.dump(out, open(fo + ".tmp", "wb")); os.replace(fo + ".tmp", fo)
    return fo

_BL = {}
def blur_of(P):
    k = P["tag"]
    if k not in _BL:
        _BL.clear(); _BL[k] = S0.make_blur("C2", P["prm"], P["pos"], P["di"], P["fm"])
    return _BL[k]

def simulate(P, s2, L, rng_f, U):
    """lambda = c * blur(F * sorgente); conteggi per inversione della Poisson con uniformi U (n_rep x n_bin), come S2.0."""
    src = C.nc_field(P["Sr"], P["r"])
    Y = []
    for k in range(U.shape[0]):
        F = S.gamma_field(src.shape, s2, L, rng_f) if s2 > 0 else 1.0
        lam = P["c"] * np.maximum(blur_of(P)(src * F), 1e-12)
        Y.append(_stats.poisson.ppf(U[k], lam))
    return Y

def grid_one(args):
    tag, OUT, ci, s2, L, rep0, nrep = args
    fo = f"{OUT}/grid/{tag}__c{ci:03d}_r{rep0:02d}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        P = pickle.load(open(f"{OUT}/prep/{tag}.pkl", "rb"))
        reps = []
        for k in range(rep0, rep0 + nrep):
            rng_f = np.random.default_rng([S.SEED, P["idx"], ci, k, 1])
            U = np.random.default_rng([S.SEED, P["idx"], ci, k, 2]).uniform(size=(1, len(P["u"])))
            y = simulate(P, s2, L, rng_f, U)[0]
            reps.append({g: S.group_stats(y, P["groups"][g], P["c"], P["rows"], P["cols"]) for g in P["groups"]})
        res = {"tag": tag, "archetype": P["archetype"], "combo": ci, "s2": s2, "L": L, "rep0": rep0, "nrep": nrep, "reps": reps,
               "obs": P["obs"], "sec": round(time.time() - t0, 1)}
    except Exception:
        res = {"tag": tag, "combo": ci, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

def check_s20(OUT):
    """G-4 (c): con i parametri di run1 (quelli di PH-2), gruppi di S2.0 non allineati, s2 = 0 e le uniformi di S2.0, il banco deve ridare
    psi osservato e simulato salvati da S2.0 in run1/disp_stripes."""
    rows_ = []
    for _, Rw in rois().iterrows():
        tag = f"{Rw.archetype}_{Rw.roi_id}"
        P = pickle.load(open(f"{OUT}/prep/{tag}.pkl", "rb"))
        ref = json.load(open(f"/mnt/micron/geo_spatialtrans/S2.0/run1/disp_stripes/{tag}.json"))
        po = S0.psi_index(P["u"], P["lab_s20"], P["c"], P["keep"])[0]
        Ud = np.random.default_rng([S0.SEED, 1000 + P["idx"]]).uniform(size=(NREP_S20, len(P["u"])))
        Y = simulate(P, 0.0, 1.0, None, Ud)
        ps = float(np.mean([S0.psi_index(y, P["lab_s20"], P["c"], P["keep"])[0] for y in Y]))
        rows_.append({"tag": tag, "r_bench": P["r"], "r_s20": ref["r"], "psi_obs_bench": po, "psi_obs_s20": ref["psi_obs"],
                      "psi_sim_bench": ps, "psi_sim_s20": ref["psi_sim"],
                      "rel_obs": abs(po / ref["psi_obs"] - 1), "rel_sim": abs(ps / ref["psi_sim"] - 1), "rel_r": abs(P["r"] / ref["r"] - 1)})
    T = pd.DataFrame(rows_); T.to_csv(f"{OUT}/check_s20.csv", index=False)
    print(T[["tag", "rel_r", "rel_obs", "rel_sim"]].to_string())
    print(f"G-4 c: max rel psi_obs {T.rel_obs.max():.2e}, psi_sim {T.rel_sim.max():.2e}, r {T.rel_r.max():.2e}")

def e1(OUT):
    """E-1: tessuto sintetico omogeneo 400 x 400 um (tutto tessuto, lontano dal bordo), sorgente costante con 40 UMI per quadratino,
    sfocatura stretta sola (0.7 um), Poisson; per ogni (s2, L) della griglia grossolana: varianza in piu' per quadratino
    v = (psi - 1) / m e rho1, rho2, gruppi di 8 x 8 quadratini; 3 repliche."""
    n_sub = 800; nb = n_sub // 4
    rows, cols = np.divmod(np.arange(nb * nb), nb)
    pos = np.c_[(rows + 0.5) * 2.0, (cols + 0.5) * 2.0]
    fm = np.ones((n_sub, n_sub))
    prm = {"sigma_e": 5.0, "D": 1e-3, "w": 0.0, "ell": 30.0}
    blur = S0.make_blur("C2", prm, pos, np.zeros(2), fm, pad_um=20.0)
    grp = (rows // 8) * 1000 + cols // 8
    src = np.full((n_sub, n_sub), 40.0 / 16.0)
    out = []
    for s2 in [0.0] + COARSE_S2:
        for L in COARSE_L:
            for k in range(3):
                rng = np.random.default_rng([S.SEED, 999, int(s2 * 1e4), int(L * 100), k])
                F = S.gamma_field(src.shape, s2, L, rng) if s2 > 0 else 1.0
                lam = np.maximum(blur(src * F), 1e-12)
                y = rng.poisson(lam).astype(float)
                st = S.group_stats(y, grp, np.ones(len(y)), rows, cols)
                out.append({"s2": s2, "L": L, "rep": k, "v": (st["psi"] - 1) / st["mbar"], "rho1": st["rho1"], "rho2": st["rho2"], "mbar": st["mbar"]})
            if s2 == 0.0: break
    T = pd.DataFrame(out).groupby(["s2", "L"]).mean(numeric_only=True).reset_index()
    T.to_csv(f"{OUT}/E1_field_scale.csv", index=False); print(T.to_string())

if __name__ == "__main__":
    mode, OUT = sys.argv[1], sys.argv[2]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 4
    t0 = time.time()
    if mode == "prep":
        PJF = sys.argv[3] if len(sys.argv) > 3 else "/mnt/micron/geo_spatialtrans/S2.0/run1_wfix/params_flat.json"
        aligned = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        os.makedirs(f"{OUT}/prep", exist_ok=True)
        jobs = [(Rw, OUT, PJF, aligned) for _, Rw in rois().iterrows()]
        with Pool(min(nw, len(jobs))) as p:
            for i, f in enumerate(p.imap_unordered(prep_one, jobs, chunksize=1), 1): print(f"{i}/{len(jobs)} {time.time()-t0:.0f}s {os.path.basename(f)}", flush=True)
    elif mode == "grid":
        combos = json.load(open(sys.argv[3])); rep0, nrep = int(sys.argv[4]), int(sys.argv[5])
        os.makedirs(f"{OUT}/grid", exist_ok=True)
        tags = [f"{Rw.archetype}_{Rw.roi_id}" for _, Rw in rois().iterrows()]
        only = os.environ.get("S23_ONLY", "").split(",") if os.environ.get("S23_ONLY") else None
        jobs = [(t, OUT, int(cb["combo"]), float(cb["s2"]), float(cb["L"]), rep0, nrep) for t in tags for cb in combos
                if (only is None or t.split("_")[0] in only) and (cb.get("archetype") in (None, t.split("_")[0]))]
        jobs = [j for j in jobs if not os.path.exists(f"{OUT}/grid/{j[0]}__c{j[2]:03d}_r{rep0:02d}.json")]
        print(f"grid jobs {len(jobs)} workers {nw}", flush=True)
        with Pool(min(nw, max(1, len(jobs)))) as p:
            for i, f in enumerate(p.imap(grid_one, jobs, chunksize=4), 1):
                if i % 20 == 0 or i == len(jobs): print(f"{i}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    elif mode == "check_s20":
        check_s20(OUT)
    elif mode == "e1":
        e1(OUT)
    print(f"done {time.time()-t0:.0f}s", flush=True)
