#!/usr/bin/env python
"""review/code/nonoverlap_opt.py — revisore "code": quanto e' conservativa nonoverlap_scale() (una passata)?
Per A2 e A4, 5 fold, seme s = 1 (x, y, a_tgt, d_reg esportati da vc_poi.R) confronta il troncamento (frazione s < 1 su
tutte le celle, come B-4) e la perdita d'area di:
 onepass  = s_i = min(1, min_j g d_ij/(R_i+R_j), g d_reg/R_i)  (ricalcolato qui; deve coincidere con s_bench)
 greedy   = iterativo: per ogni coppia violata si riduce solo il nucleo piu' grande quanto basta (fino a convergenza)
 lp_r     = LP max sum R_i s_i      s.t. s_i R_i + s_j R_j <= g d_ij (coppie in conflitto), s_i R_i <= g d_reg_i, 0<=s<=1
 lp_a     = LP max sum R_i^2 s_i    (stessi vincoli)
 mvc      = MILP min #troncati: t_i + t_j >= 1 per coppia in conflitto; t forzato se R_i > g d_reg_i o se R_i > g d_ij
            (limite inferiore esatto della frazione troncata con nuclei circolari, centri fissi e sola riduzione)
Uso: python nonoverlap_opt.py <dir_parquet> <out.csv>"""
import sys, glob, os, time
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
K32 = np.sqrt(2 * np.pi / (32 * np.sin(2 * np.pi / 32))); G = 0.99

def ks(a, b):
    a = np.sort(a); b = np.sort(b); v = np.unique(np.concatenate([a, b]))
    return np.max(np.abs(np.searchsorted(a, v, "right") / len(a) - np.searchsorted(b, v, "right") / len(b)))

def run(f):
    d = pd.read_parquet(f); x, y = d.x.to_numpy(), d.y.to_numpy(); n = len(x)
    R = K32 * np.sqrt(d.a_tgt.to_numpy() / np.pi); dreg = d.d_reg.to_numpy()
    tr = cKDTree(np.column_stack([x, y])); pr = tr.query_pairs(2 * R.max() / G, output_type="ndarray")
    i, j = pr[:, 0], pr[:, 1]; dij = np.hypot(x[i] - x[j], y[i] - y[j]); bad = R[i] + R[j] > G * dij
    i, j, dij = i[bad], j[bad], dij[bad]; m = len(i)
    out = {}
    # onepass
    f_ = G * dij / (R[i] + R[j]); sp = np.ones(n); np.minimum.at(sp, i, f_); np.minimum.at(sp, j, f_)
    s1 = np.minimum(sp, np.minimum(1, G * dreg / R)); out["onepass"] = s1
    # greedy: riduci il piu' grande (in raggio corrente) quanto basta
    s = np.minimum(1, G * dreg / R)
    for it in range(1000):
        r = s * R; v = r[i] + r[j] - G * dij; vi = np.where(v > 1e-12)[0]
        if not len(vi): break
        for k in vi[np.argsort(-v[vi])]:
            a, b = i[k], j[k]; ra, rb = s[a] * R[a], s[b] * R[b]
            if ra + rb <= G * dij[k] + 1e-12: continue
            big, sml = (a, b) if ra >= rb else (b, a); rs = s[sml] * R[sml]
            if rs < G * dij[k]: s[big] = (G * dij[k] - rs) / R[big]
            else: f2 = G * dij[k] / (ra + rb); s[a] *= f2; s[b] *= f2
    out["greedy"] = s.copy()
    # LP
    rows = np.concatenate([np.arange(m), np.arange(m)]); cols = np.concatenate([i, j]); vals = np.concatenate([R[i], R[j]])
    Aub = coo_matrix((vals, (rows, cols)), shape=(m, n)).tocsr(); bub = G * dij
    ub = np.minimum(1, G * dreg / R)
    for nm, cvec in (("lp_r", -R), ("lp_a", -R**2)):
        t0 = time.time(); res = linprog(cvec, A_ub=Aub, b_ub=bub, bounds=np.column_stack([np.zeros(n), ub]), method="highs")
        out[nm] = np.clip(res.x, 0, 1) if res.status == 0 else np.full(n, np.nan)
    # MVC
    forced = (R > G * dreg)
    np.logical_or.at(forced, i, R[i] > G * dij); np.logical_or.at(forced, j, R[j] > G * dij)
    Acov = coo_matrix((np.ones(2 * m), (rows, cols)), shape=(m, n)).tocsr()
    lb = forced.astype(float)
    t0 = time.time()
    res = milp(np.ones(n), constraints=[LinearConstraint(Acov, lb=np.ones(m), ub=np.inf)], integrality=np.ones(n), bounds=Bounds(lb, np.ones(n)),
               options=dict(time_limit=600, mip_rel_gap=0))
    tmv = time.time() - t0
    rec = dict(file=os.path.basename(f), n=n, n_conflict_pairs=m, frac_in_conflict=float(np.mean(np.isin(np.arange(n), np.concatenate([i, j])) | (R > G * dreg))),
               onepass_equals_bench=bool(np.allclose(s1, d.s_bench.to_numpy(), rtol=0, atol=1e-12)), vc_trunc=float(d.vc_clamped.mean()))
    for nm, sv in out.items():
        rr = sv * R; ok = np.all(rr[i] + rr[j] <= G * dij + 1e-9) and np.all(rr <= G * dreg + 1e-9)
        rec[f"{nm}_feasible"] = bool(ok); rec[f"{nm}_trunc"] = float(np.mean(sv < 1 - 1e-9)); rec[f"{nm}_loss"] = float(np.mean(1 - sv**2))
        ev = d.ev.to_numpy(); a_fin = sv**2 * d.a_tgt.to_numpy()
        rec[f"{nm}_ks_an"] = float(ks(a_fin[ev], d.a_n_real.to_numpy()[ev])); rec[f"{nm}_ks_nc"] = float(ks(a_fin[ev] / d.a_t.to_numpy()[ev], d.a_n_real.to_numpy()[ev] / d.a_t.to_numpy()[ev]))
    rec["mvc_status"] = int(res.status); rec["mvc_trunc"] = float(res.x.round().mean()) if res.x is not None else np.nan
    rec["mvc_lp_bound"] = np.nan; rec["mvc_gap"] = float(getattr(res, "mip_gap", np.nan)) if res.x is not None else np.nan; rec["mvc_time_s"] = tmv
    rec["mvc_forced_frac"] = float(forced.mean())
    return rec

if __name__ == "__main__":
    rows = [run(f) for f in sorted(glob.glob(os.path.join(sys.argv[1], "*_s1.parquet")))]
    df = pd.DataFrame(rows); df.to_csv(sys.argv[2], index=False); print(df.T.to_string())
