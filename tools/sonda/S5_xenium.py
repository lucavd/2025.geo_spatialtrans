"""tools/sonda/S5_xenium.py - S5 (docs/sonda/REGOLA_S5.md): braccio XB, Xenium in quadratini sulla sezione X1, aree omologhe di questa regola.
Modi:
  prep           per area (quadrato + 240 um di margine, solo meta' di sviluppo x < m - 50 um): UMI del Visium HD di X1 (filtrata) sui geni di G_p,
                 trascritti Xenium con qv >= 20 sui geni di G_p portati sulla griglia da 2 um con la composita di S1 (S1_registration.load_T, variante
                 scelta), poligoni dei nuclei Xenium (nucleus_boundaries.parquet, stessa trasformazione). <S5_ROOT>/xen/<unit>/{vis,tx,nuc}.parquet
  reg [W]        metro di S1 invariato (S1_registration: coarse_shift, tiles_run, L = 200 um, finestra +-12 um, quadrica, validita' >= 500 conteggi):
                 riquadri prima della correzione; residuo sulla scacchiera come S1 stage_main (stima sui pari, valutazione sui dispari e viceversa);
                 affine a 6 parametri centri -> centri + s su tutti i riquadri validi e non di bordo; correzione p = (I + A)^-1 (q - t) dei trascritti;
                 riquadri rifatti dopo la correzione (e, per il segno, con il segno opposto q + s(q)); controllo positivo di S1 sui trascritti corretti.
                 Cancello G-XR per area. <S5_ROOT>/xen/<unit>/reg.json, tiles_*.csv, positive.csv
  export         quadratini da 2 um di X1 con centro nel quadrato dell'area (in_tissue, in_reg come S3_export); trascritti corretti al quadratino del
                 centro piu' vicino; label = nucleo Xenium (poligono corretto) che contiene il centro del quadratino, tutti i nuclei (XD-9), numero =
                 posizione (da 1) della cellula in cells.parquet di Xenium; cells.parquet con cell_type vuoto. <S5_ROOT>/xen/<unit>/raw/
  thin           diradamento (XD-11): f per tessuto = media sulle 5 aree di Vp degli UMI per quadratino in_reg / la stessa media di XB grezzo;
                 ogni conteggio (quadratino, gene) -> Binomiale(n, f), default_rng([20261130, i]) nell'ordine di X.parquet (quadratino, gene).
                 Se f > 1 XB resta alla sua profondita'. <S5_ROOT>/runs/XB/<unit>/, results/sonda/S5/thinning.json
Limite dichiarato (REGOLA_S5): la correzione usa l'espressione del Visium HD di X1; nessuna misura M-1...M-6 entra nella stima.
Uso: .venv/bin/python -B tools/sonda/S5_xenium.py <prep|reg|export|thin> [W]"""
import os, sys, json, time, shutil, numpy as np, pandas as pd, pyarrow.parquet as pq, scipy.sparse as sp
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S5_common as S, S1_registration as S1
C = S.C; mode = sys.argv[1]; W = int(sys.argv[2]) if len(sys.argv) > 2 else 1
MB = int(round(S.MARGIN_UM / 2)); INJ = [(1.0, 0.0), (0.0, 2.5), (1.5, -1.5)]

def xdir(u): return f"{S.ROOT}/xen/{u['unit']}"
def box(u):
    c0 = int(round((u["x0"] - S.MARGIN_UM) / 2)); c1 = int(round((u["x0"] + u["side"] + S.MARGIN_UM) / 2))
    r0 = int(round((u["y0"] - S.MARGIN_UM) / 2)); r1 = int(round((u["y0"] + u["side"] + S.MARGIN_UM) / 2)); return c0, c1, r0, r1
def dev_ok_col(col_int, xmax): return (np.asarray(col_int) + 0.5) * 2.0 < xmax

def prep():
    T, s = S1.load_T("X1"); Tinv = np.linalg.inv(T); G = S.gp(); gi = {e: i for i, e in enumerate(G)}; s2e = S.sym2ens(); xmax = S.xmax_dev()
    A = S.areas_x1(); tp = C.positions("X1", "002").astype({"array_col": np.int64, "array_row": np.int64}); h5 = f"{C.bdir('X1', '002')}/filtered_feature_bc_matrix.h5"
    ids, _ = C.features(h5); gsel = pd.Series(np.arange(len(ids)), index=ids).reindex(G).to_numpy().astype(np.int64)
    bpos = pd.Series(np.arange(len(C.barcodes(h5))), index=C.barcodes(h5)); info = {}
    for u in A:
        os.makedirs(xdir(u), exist_ok=True); c0, c1, r0, r1 = box(u)
        b = tp[(tp.array_col >= c0) & (tp.array_col < c1) & (tp.array_row >= r0) & (tp.array_row < r1)]
        b = b[dev_ok_col(b.array_col.to_numpy(), xmax)]                              # meta' di sviluppo
        col = bpos.reindex(b.barcode.to_numpy()).to_numpy(); have = ~np.isnan(col); b = b[have]
        X, _ = C.read_cols(h5, col[have].astype(np.int64)); X = X[:, gsel].tocoo()
        V = pd.DataFrame({"g": X.col.astype(np.int64), "c": b.array_col.to_numpy()[X.row] - c0, "r": b.array_row.to_numpy()[X.row] - r0, "n": X.data.astype(np.int64)})
        V.to_parquet(f"{xdir(u)}/vis.parquet", index=False)
        corners = np.array([[c0, r0], [c1, r0], [c0, r1], [c1, r1]], float); xb = S1.apply_h(Tinv, corners)
        info[u["unit"]] = dict(c0=c0, c1=c1, r0=r0, r1=r1, xbox=[(xb.min(0) - 30).tolist(), (xb.max(0) + 30).tolist()], vis_umi=int(V.n.sum()),
                               vis_bins_xmax_um=float(((b.array_col.to_numpy() + 0.5) * 2).max()), vis_n_bins=int(len(b)))
        print(u["unit"], "visium umi", int(V.n.sum()), flush=True)
    pf = pq.ParquetFile(f"{S.XPRE}_transcripts.parquet"); acc = {k: [] for k in info}; n_tot = n_q = 0
    for bt in pf.iter_batches(columns=["feature_name", "x_location", "y_location", "qv"], batch_size=5_000_000):
        df = bt.to_pandas(); n_tot += len(df)
        if len(df) and isinstance(df.feature_name.iloc[0], bytes): df["feature_name"] = df.feature_name.str.decode("utf-8")
        df = df[df.qv.to_numpy() >= 20]; n_q += len(df)
        g = df.feature_name.astype(str).map(s2e).map(gi); df = df[g.notna().to_numpy()]; g = g[g.notna()].to_numpy().astype(np.int64)
        for name, I in info.items():
            lo, hi = I["xbox"]; x = df.x_location.to_numpy(); y = df.y_location.to_numpy()
            k = (x >= lo[0]) & (x <= hi[0]) & (y >= lo[1]) & (y <= hi[1])
            if k.any(): acc[name].append(pd.DataFrame({"g": g[k], "x": x[k], "y": y[k]}))
    for u in A:
        I = info[u["unit"]]; d = pd.concat(acc[u["unit"]]); a = S1.apply_h(T, d[["x", "y"]].to_numpy(float))
        k = dev_ok_col(np.floor(a[:, 0] + 0.5), xmax)                                  # trascritti portati nella meta' tenuta da parte: esclusi
        tx = pd.DataFrame({"g": d.g.to_numpy()[k], "c": a[k, 0] - I["c0"], "r": a[k, 1] - I["r0"]}); tx.to_parquet(f"{xdir(u)}/tx.parquet", index=False)
        I.update(tx_n=int(len(tx)), tx_dropped_holdout=int((~k).sum())); print(u["unit"], "xenium tx", len(tx), flush=True)
    cells = pq.read_table(f"{S.XPRE}_cells.parquet", columns=["cell_id"]).to_pandas().cell_id.astype(str).to_numpy()
    lab = pd.Series(np.arange(1, len(cells) + 1), index=cells)
    nb = pq.read_table(f"{S.XPRE}_nucleus_boundaries.parquet", columns=["cell_id", "vertex_x", "vertex_y"]).to_pandas()
    a = S1.apply_h(T, nb[["vertex_x", "vertex_y"]].to_numpy(float)); nb["c"] = a[:, 0]; nb["r"] = a[:, 1]; nb["label"] = lab.reindex(nb.cell_id.astype(str).to_numpy()).to_numpy()
    assert not np.isnan(nb.label.to_numpy(float)).any()
    for u in A:
        I = info[u["unit"]]; cc0, cc1 = u["x0"] / 2 - 30, (u["x0"] + u["side"]) / 2 + 30; rr0, rr1 = u["y0"] / 2 - 30, (u["y0"] + u["side"]) / 2 + 30
        ins = (nb.c.to_numpy() >= cc0) & (nb.c.to_numpy() <= cc1) & (nb.r.to_numpy() >= rr0) & (nb.r.to_numpy() <= rr1)
        L = np.unique(nb.label.to_numpy()[ins]); q = nb[nb.label.isin(L)].copy(); q["c"] -= I["c0"]; q["r"] -= I["r0"]
        q[["label", "cell_id", "c", "r"]].assign(label=q.label.astype(np.int64)).to_parquet(f"{xdir(u)}/nuc.parquet", index=False)
        I.update(nuc_n=int(len(L))); json.dump(dict(u, **I, transcripts_total=n_tot, transcripts_qv20=n_q), open(f"{xdir(u)}/prep.json", "w"), indent=1)
        print(u["unit"], "nuclei", len(L), flush=True)

def p95(v): return float(np.percentile(v, 95)) if len(v) else float("nan")
def s_norm(t):
    v = t[t.valid.astype(bool)]; e = v.edge.astype(bool).to_numpy()
    return np.where(e, np.inf, np.hypot(v.sx_um.to_numpy(float), v.sy_um.to_numpy(float))), v
def checker(t):
    """residuo locale di S1 stage_main: affine stimata sui riquadri pari della scacchiera, valutata sui dispari, e viceversa."""
    v = t[t.valid.astype(bool)].copy(); e = v.edge.astype(bool).values; res_loc = np.full(len(v), 1e6); par = ((v.tx + v.ty) % 2).values; good = ~e
    Aa = np.c_[v.cx_um.values, v.cy_um.values, np.ones(len(v))]
    for f in (0, 1):
        fit = (par == f) & good; ev = par != f
        if fit.sum() >= 6:
            cxy = np.linalg.lstsq(Aa[fit], v[["sx_um", "sy_um"]].values[fit], rcond=None)[0]; pr = Aa @ cxy
            rl = np.hypot(v.sx_um.values - pr[:, 0], v.sy_um.values - pr[:, 1]); res_loc[ev] = np.where(good[ev], rl[ev], 1e6)
    return res_loc
def to_um(c): return 2.0 * (np.asarray(c, float) - MB) + 1.0        # centro del riquadro tx = (tx + 0.5) L, come cx_um di S1.tiles_run
def to_bin(u): return (np.asarray(u, float) - 1.0) / 2.0 + MB
def correct(c, r, coef, sign=+1):
    """p = (I + A)^-1 (q - t) (sign +1, modello q = p + s(p)); sign -1: q + s(q) (segno opposto, solo per la verifica)."""
    q = np.c_[to_um(c), to_um(r)]; Lin = coef[:2, :].T; t = coef[2, :]
    p = np.linalg.solve(np.eye(2) + Lin, (q - t).T).T if sign > 0 else q + q @ Lin.T + t
    return to_bin(p[:, 0]), to_bin(p[:, 1])

def reg_one(u):
    t0 = time.time(); I = json.load(open(f"{xdir(u)}/prep.json")); V = pd.read_parquet(f"{xdir(u)}/vis.parquet"); X = pd.read_parquet(f"{xdir(u)}/tx.parquet")
    nc, nr = I["c1"] - I["c0"], I["r1"] - I["r0"]; side = u["side"]; sb = int(round(side / 2)); L = S.L_UM
    Vp = S1.Pts(V.g.values, V.c.values.astype(float), V.r.values.astype(float), V.n.values); Xp = S1.Pts(X.g.values, X.c.values, X.r.values, np.ones(len(X), np.float32))
    cx, cy, ccv = S1.coarse_shift(Xp, Vp, nc, nr, sb, MB); t = S1.tiles_run(Xp, Vp, side, L, MB, (cx, cy))
    rl = checker(t); r_before, v = s_norm(t); good = ~v.edge.astype(bool).to_numpy()
    Aa = np.c_[v.cx_um.values, v.cy_um.values, np.ones(len(v))]; coef = np.linalg.lstsq(Aa[good], v[["sx_um", "sy_um"]].values[good], rcond=None)[0]
    out = dict(unit=u["unit"], coarse_bins=[int(cx), int(cy)], coarse_cc=float(ccv), n_tiles=int(len(t)), n_valid=int(len(v)), n_edge=int((~good).sum()), n_fit=int(good.sum()),
               s_before_median_um=float(np.median(r_before)), s_before_p95_um=p95(r_before), checker_p95_um=p95(rl), checker_median_um=float(np.median(rl)),
               affine_coef=coef.tolist())
    res = {}
    for nm, sg in (("after", +1), ("wrongsign", -1)):
        c2, r2 = correct(X.c.values, X.r.values, coef, sg); Xc = S1.Pts(X.g.values, c2, r2, np.ones(len(X), np.float32))
        cx2, cy2, cc2 = S1.coarse_shift(Xc, Vp, nc, nr, sb, MB); t2 = S1.tiles_run(Xc, Vp, side, L, MB, (cx2, cy2)); t2.to_csv(f"{xdir(u)}/tiles_{nm}.csv", index=False)
        rn, v2 = s_norm(t2); res[nm] = (Xc, (cx2, cy2), t2)
        out[f"{nm}_coarse_bins"] = [int(cx2), int(cy2)]; out[f"{nm}_n_valid"] = int(len(v2)); out[f"{nm}_n_edge"] = int(v2.edge.astype(bool).sum())
        out[f"{nm}_s_median_um"] = float(np.median(rn)) if len(rn) else float("nan"); out[f"{nm}_s_p95_um"] = p95(rn)
    Xc, (cx2, cy2), t2 = res["after"]; pos = []
    for inj in INJ:
        ti = S1.tiles_run(Xc.shifted(inj[0] / 2, inj[1] / 2), Vp, side, L, MB, (cx2, cy2)); k = S1.okmask(t2, ti)
        dx = float(np.median(ti.sx_um[k] - t2.sx_um[k])); dy = float(np.median(ti.sy_um[k] - t2.sy_um[k]))
        pos.append(dict(unit=u["unit"], inj_x=inj[0], inj_y=inj[1], rec_x=dx, rec_y=dy, n=int(k.sum()), ok=bool(abs(dx - inj[0]) <= 0.3 and abs(dy - inj[1]) <= 0.3)))
    pd.DataFrame(pos).to_csv(f"{xdir(u)}/positive.csv", index=False); t.to_csv(f"{xdir(u)}/tiles_before.csv", index=False)
    out.update(positive=pos, pos_ok=bool(all(p["ok"] for p in pos)), sign_ok=bool(out["after_s_median_um"] < out["wrongsign_s_median_um"]),
               gxr_checker=bool(out["checker_p95_um"] <= 2.0), gxr_after=bool(out["after_s_median_um"] <= 0.5), sec=round(time.time() - t0, 1))
    out["gxr_pass"] = bool(out["gxr_checker"] and out["gxr_after"] and out["pos_ok"])
    json.dump(out, open(f"{xdir(u)}/reg.json", "w"), indent=1); print(json.dumps({k: v for k, v in out.items() if k not in ("positive", "affine_coef")}), flush=True)
    return out

def export_one(u):
    t0 = time.time(); I = json.load(open(f"{xdir(u)}/prep.json")); rg = json.load(open(f"{xdir(u)}/reg.json")); coef = np.array(rg["affine_coef"])
    m = S.m_x1(); tp = C.positions("X1", "002").astype({"array_col": np.int64, "array_row": np.int64}); G = S.gp(); d = f"{xdir(u)}/raw"; os.makedirs(d, exist_ok=True)
    x = (tp.array_col.to_numpy() + 0.5) * 2.0; y = (tp.array_row.to_numpy() + 0.5) * 2.0
    s = (x >= u["x0"]) & (x < u["x0"] + u["side"]) & (y >= u["y0"]) & (y < u["y0"] + u["side"])
    b = tp[s].copy(); xm = C.guard((b.array_col.to_numpy() + 0.5) * 2.0, m)
    b["x_um"] = (b.array_col.to_numpy() + 0.5) * 2.0 - u["x0"]; b["y_um"] = (b.array_row.to_numpy() + 0.5) * 2.0 - u["y0"]
    b = b.sort_values(["array_row", "array_col"]).reset_index(drop=True); side = u["side"]
    bx, by = b.x_um.to_numpy(), b.y_um.to_numpy(); it = b.in_tissue.to_numpy() == 1
    in_reg = it & (bx >= C.FRAME_UM) & (bx <= side - C.FRAME_UM) & (by >= C.FRAME_UM) & (by <= side - C.FRAME_UM)
    rmin, cmin = int(b.array_row.min()), int(b.array_col.min()); grid = -np.ones((int(b.array_row.max()) - rmin + 1, int(b.array_col.max()) - cmin + 1), np.int64)
    grid[b.array_row.to_numpy() - rmin, b.array_col.to_numpy() - cmin] = np.arange(len(b))
    def lookup(cabs, rabs):
        ci = np.floor(cabs + 0.5).astype(np.int64) - cmin; ri = np.floor(rabs + 0.5).astype(np.int64) - rmin
        ok = (ci >= 0) & (ci < grid.shape[1]) & (ri >= 0) & (ri < grid.shape[0]); out = np.full(len(cabs), -1, np.int64); out[ok] = grid[ri[ok], ci[ok]]; return out
    X = pd.read_parquet(f"{xdir(u)}/tx.parquet"); c2, r2 = correct(X.c.values, X.r.values, coef, +1); bi = lookup(c2 + I["c0"], r2 + I["r0"]); k = bi >= 0
    Xm = sp.csr_matrix((np.ones(k.sum(), np.int64), (bi[k], X.g.values[k])), shape=(len(b), len(G))); Xm.sum_duplicates()
    # nuclei: poligono corretto; centro del quadratino (col, row interi = centri nella griglia continua) dentro il poligono
    from matplotlib.path import Path
    N = pd.read_parquet(f"{xdir(u)}/nuc.parquet"); nc2, nr2 = correct(N.c.values, N.r.values, coef, +1); N["ca"] = nc2 + I["c0"]; N["ra"] = nr2 + I["r0"]
    label = np.zeros(len(b), np.int64); conflicts = 0; cid = {}
    for lab_, q in N.groupby("label", sort=True):
        pc, pr = q.ca.to_numpy(), q.ra.to_numpy(); c_lo, c_hi = int(np.ceil(pc.min())), int(np.floor(pc.max())); r_lo, r_hi = int(np.ceil(pr.min())), int(np.floor(pr.max()))
        if c_hi < c_lo or r_hi < r_lo: continue
        cc, rr = np.meshgrid(np.arange(c_lo, c_hi + 1), np.arange(r_lo, r_hi + 1)); cc = cc.ravel(); rr = rr.ravel()
        inside = Path(np.c_[pc, pr]).contains_points(np.c_[cc, rr]); bb = lookup(cc[inside].astype(float), rr[inside].astype(float)); bb = bb[bb >= 0]
        if len(bb): conflicts += int((label[bb] > 0).sum()); label[bb] = int(lab_); cid[int(lab_)] = str(q.cell_id.iloc[0])
    bins = pd.DataFrame({"bin_id": np.arange(len(b), dtype=np.int64), "row": b.array_row.to_numpy().astype(np.int64), "col": b.array_col.to_numpy().astype(np.int64),
                         "x_um": bx, "y_um": by, "in_tissue": it, "label": label, "in_reg": in_reg})
    bins.to_parquet(f"{d}/bins.parquet", index=False); S.write_X(f"{d}/X.parquet", Xm)
    with open(f"{d}/genes.tsv", "w") as f: f.write("\n".join(G) + "\n")
    labs = np.unique(label[label > 0])
    pd.DataFrame({"label": labs, "cell_type": pd.Series([None] * len(labs), dtype=object), "xenium_cell_id": [cid[int(l)] for l in labs]}).to_parquet(f"{d}/cells.parquet", index=False)
    info = dict(u, section="X1", arm="XB", side_um=side, n_bins=int(len(bins)), n_tissue=int(it.sum()), n_reg=int(in_reg.sum()), n_nuclei=int(len(labs)),
                umi_total=int(Xm.sum()), umi_reg=int(Xm[in_reg].sum()), tx_in_region=int(len(X)), tx_in_area=int(k.sum()), label_conflict_bins=conflicts,
                frac_bins_labelled=float((label > 0).mean()), m_um=m, xmax_um=xm, n_genes=len(G), sec_export=round(time.time() - t0, 1))
    json.dump(info, open(f"{d}/export.json", "w")); print(json.dumps(info), flush=True); return info

def umi_per_reg_bin(d):
    bins = pd.read_parquet(f"{d}/bins.parquet"); T = pd.read_parquet(f"{d}/X.parquet"); reg = bins.in_reg.to_numpy()
    return float(T.n.to_numpy()[reg[T.bin.to_numpy()]].sum() / reg.sum())

def thin():
    A = S.areas_x1(); res = {"rule": "docs/sonda/REGOLA_S5.md (XD-11)", "seed": S.SEED_THIN, "by_tissue": {}}
    for tis in ("tumour", "stroma"):
        vp = {k: umi_per_reg_bin(f"{S.ROOT}/runs/Vp/H2_{tis}_{k}") for k in range(1, 6)}
        xb = {k: umi_per_reg_bin(f"{S.ROOT}/xen/X1_{tis}_{k}/raw") for k in range(1, 6)}
        f = float(np.mean(list(vp.values())) / np.mean(list(xb.values())))
        res["by_tissue"][tis] = dict(f=f, Vp_umi_per_reg_bin=vp, XB_raw_umi_per_reg_bin=xb, Vp_mean=float(np.mean(list(vp.values()))), XB_raw_mean=float(np.mean(list(xb.values()))),
                                     thinned=bool(f <= 1.0))
    for u in A:
        s = f"{xdir(u)}/raw"; d = f"{S.ROOT}/runs/XB/{u['unit']}"; f = res["by_tissue"][u["tissue"]]["f"]
        if os.path.exists(f"{d}/export.json"): print("present", u["unit"]); continue
        os.makedirs(d, exist_ok=True)
        for fn in ("bins.parquet", "cells.parquet", "genes.tsv"): shutil.copy(f"{s}/{fn}", f"{d}/{fn}")
        T = pd.read_parquet(f"{s}/X.parquet"); rng = np.random.default_rng([S.SEED_THIN, u["i"]])
        n2 = rng.binomial(T.n.to_numpy(), f) if f <= 1.0 else T.n.to_numpy()
        k = n2 > 0; Y = pd.DataFrame({"bin": T.bin.to_numpy()[k], "gene": T.gene.to_numpy()[k], "n": n2[k].astype(np.int64)}); Y.to_parquet(f"{d}/X.parquet", index=False)
        info = dict(json.load(open(f"{s}/export.json")), f=f, thin_seed=[S.SEED_THIN, u["i"]], umi_total_raw=int(T.n.sum()), umi_total=int(Y.n.sum()))
        json.dump(info, open(f"{d}/export.json", "w")); print("thin", u["unit"], f, int(T.n.sum()), int(Y.n.sum()), flush=True)
    for u in A: res.setdefault("XB_umi_per_reg_bin_after", {})[u["unit"]] = umi_per_reg_bin(f"{S.ROOT}/runs/XB/{u['unit']}")
    json.dump(res, open(f"{S.RES}/thinning.json", "w"), indent=1); print(json.dumps({t: v["f"] for t, v in res["by_tissue"].items()}))

if __name__ == "__main__":
    os.makedirs(f"{S.ROOT}/xen", exist_ok=True); A = S.areas_x1()
    sel = os.environ.get("S5_UNITS"); A = [u for u in A if sel is None or u["unit"] in sel.split(",")]
    if mode == "prep": prep()
    elif mode in ("reg", "export"):
        fn = reg_one if mode == "reg" else export_one; todo = [u for u in A if not os.path.exists(f"{xdir(u)}/reg.json" if mode == "reg" else f"{xdir(u)}/raw/export.json")]
        if W > 1 and len(todo) > 1:
            from multiprocessing import Pool
            with Pool(min(W, len(todo))) as p: p.map(fn, todo)
        else:
            for u in todo: fn(u)
    elif mode == "thin": thin()
