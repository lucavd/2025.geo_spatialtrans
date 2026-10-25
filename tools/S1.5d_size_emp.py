"""tools/S1.5d_size_emp.py - S1.5d, esplorazione dichiarata: residui del modello di taglia per gruppo (asimmetria, code) e candidato M3 a residui
empirici: log A_N = a + b log A_T + r, r estratto dai residui veri del gruppo nella stessa classe di log A_T (decili delle ROI di stima).
Lascia-fuori-un-ROI sui territori veri, come tools/S1.5d_size_refs.py. Uscite: results/S1.5d/S1.5d_size_resid.csv, S1.5d_size_check_m3.csv"""
import importlib.util, numpy as np, pandas as pd
from scipy.stats import skew
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); TC = pd.read_csv("results/S1.5b/S1.5b_type_calibration.csv")
R, C = [], []
for A in cal.ARCH:
    _, d = cal.load(A); e = d[d.interior.astype(bool)].copy(); tc = TC[TC.archetype == A]
    e["g"] = e.type.map(dict(zip(tc.cell_type, tc.size_group))); e = e[e.g.notna()]
    for gi, gg in e.groupby("g"):
        x = np.log(gg.a_t.values); y = np.log(gg.a_n.values); b = np.polyfit(x, y, 1); r = y - np.polyval(b, x); s = r.std(ddof=2)
        R.append(dict(archetype=A, size_group=int(gi), n=len(gg), sigma=s, skew=skew(r), q99_z=np.quantile(r / s, 0.99), q999_z=np.quantile(r / s, 0.999),
                      q01_z=np.quantile(r / s, 0.01), gauss_q99=2.326))
        for roi in sorted(gg.roi_id.unique()):
            tr, te = gg[gg.roi_id != roi], gg[gg.roi_id == roi]
            x = np.log(tr.a_t.values); y = np.log(tr.a_n.values); b = np.polyfit(x, y, 1); res = y - np.polyval(b, x)
            edges = np.quantile(x, np.linspace(0, 1, 11)); edges[0], edges[-1] = -np.inf, np.inf
            bi = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, 9); pools = [res[bi == k] for k in range(10)]
            xt = np.log(te.a_t.values); bt = np.clip(np.searchsorted(edges, xt, side="right") - 1, 0, 9); at = te.a_t.values
            rng = np.random.default_rng(20261025); an3 = np.empty((40, len(xt)))
            for k in range(10):
                w = bt == k
                if w.any(): an3[:, w] = np.exp(np.polyval(b, xt[w]) + rng.choice(pools[k], size=(40, w.sum())))
            for t in sorted(te.type.unique()):
                k = (te.type == t).values
                if k.sum() < 30: continue
                C.append(dict(archetype=A, size_group=int(gi), roi_id=roi, cell_type=t, real_nc095=np.mean(te.a_n.values[k] / at[k] >= 0.95),
                              M3_nc095=np.mean(an3[:, k] / at[k] >= 0.95), real_q90=np.quantile(te.a_n.values[k], 0.9), M3_q90=np.quantile(an3[:, k], 0.9),
                              real_med=np.median(te.a_n.values[k]), M3_med=np.median(an3[:, k])))
pd.DataFrame(R).to_csv("results/S1.5d/S1.5d_size_resid.csv", index=False); pd.DataFrame(C).to_csv("results/S1.5d/S1.5d_size_check_m3.csv", index=False); print("ok")
