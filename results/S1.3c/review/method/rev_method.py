"""results/S1.3c/review/method/rev_method.py — verifiche del revisore «Metodo» (S1.3c). Sola lettura sugli input.
Uso (root del repo): .venv/bin/python results/S1.3c/review/method/rev_method.py <edge|mask|reg|indep|claims|all>
Uscite: results/S1.3c/review/method/out/*.csv"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage, stats
from scipy.spatial import cKDTree
import importlib.util
ST = sys.argv[1]; sys.argv = sys.argv[:1]
ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location("e", ROOT / "tools/S1.3c_ecm.py"); E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
OUTD = Path(__file__).resolve().parent / "out"; OUTD.mkdir(exist_ok=True)
S13 = Path("/mnt/micron/geo_spatialtrans/S1.3c")

def ctx(r):
    b = pd.read_parquet(S13 / "bins" / f"{r.dataset}_roi_bins.parquet"); b["umi"] = b.umi_com
    return b

def ratio(d, u):
    a = u[(d >= 4) & (d < 8)]; b = u[(d >= 0) & (d < 2)]
    return float(a.mean() / b.mean()) if len(a) > 50 and len(b) > 50 and b.mean() > 0 else np.nan

def setup(r):
    A = r.archetype; terr_m, mask_m = E.PRIMARY[A]
    s, valid4, px4, vb, nm, dmap, dnuc = E.load_roi_context(r, ctx(r), mask_m)
    gx, gy, _ = E.generators(A, r.roi_id, terr_m, valid4, px4)
    terr = cKDTree(np.c_[gx, gy]).query(np.c_[s.x.values, s.y.values])[1]
    return s, valid4, px4, vb, nm, dmap, dnuc, terr, len(gx)

# ---------------------------------------------------------------- (a) gradiente vs bordo del tessuto
def edge_one(args):
    i, r = args
    s, valid4, px4, vb, nm, dmap, dnuc, terr, ng = setup(r)
    u = s.umi.values.astype(float)
    yy = np.clip((s.y / px4).astype(int), 0, valid4.shape[0] - 1).values; xx = np.clip((s.x / px4).astype(int), 0, valid4.shape[1] - 1).values
    de = (ndimage.distance_transform_edt(valid4) * px4)[yy, xx]           # distanza dal bordo della maschera valida (tessuto, lumi, bolle)
    dr = np.minimum.reduce([s.x.values, s.y.values, r.side_um - s.x.values, r.side_um - s.y.values])   # distanza dal bordo del ROI
    tsum = np.bincount(terr[vb], u[vb], minlength=ng); tn = np.bincount(terr[vb], minlength=ng)
    tm = np.where(tn > 0, tsum / np.maximum(tn, 1), np.nan)[terr]
    rec = {"archetype": r.archetype, "roi_id": r.roi_id, "n_valid": int(vb.sum())}
    sets = {"all": vb, "edge20": vb & (de >= 20), "edge20_roi20": vb & (de >= 20) & (dr >= 20), "edge40_roi20": vb & (de >= 40) & (dr >= 20)}
    for k, m in sets.items():
        rec[f"R_{k}"] = ratio(dnuc[m], u[m]); rec[f"n_{k}"] = int(m.sum())
        ok = m & (tm > 0)
        rec[f"Rcomp_{k}"] = ratio(dnuc[ok], tm[ok])                          # R atteso sotto il nullo (densita' uniforme nel territorio)
        rec[f"Rwithin_{k}"] = ratio(dnuc[ok], u[ok] / tm[ok])                # R a parita' di territorio
    # quota di bin vicini al bordo per classe di d_nuc
    for lab, (lo, hi) in {"d0_2": (0, 2), "d4_8": (4, 8)}.items():
        m = vb & (dnuc >= lo) & (dnuc < hi); rec[f"frac_edge<20_{lab}"] = float((de[m] < 20).mean()) if m.sum() else np.nan
        rec[f"median_dedge_{lab}"] = float(np.median(de[m])) if m.sum() else np.nan
    rec["spearman_dnuc_dedge"] = float(stats.spearmanr(dnuc[vb], de[vb])[0])
    rec["spearman_umi_dedge"] = float(stats.spearmanr(u[vb], de[vb])[0])
    for lab, (lo, hi) in {"de0_10": (0, 10), "de10_20": (10, 20), "de20_40": (20, 40), "de40+": (40, 1e9)}.items():
        m = vb & (de >= lo) & (de < hi); rec[f"R_{lab}"] = ratio(dnuc[m], u[m]); rec[f"umi_{lab}"] = float(u[m].mean()) if m.sum() else np.nan
    # placebo CP-1 identico (stesso seme) ristretto ai bin lontani dal bordo
    upp = r.um_per_px; rng = np.random.default_rng(20261013 + int(i))
    ang = rng.uniform(0, 2 * np.pi); mod = rng.uniform(30, 60)
    nm_p = np.roll(nm, (int(round(mod * np.sin(ang) / upp)), int(round(mod * np.cos(ang) / upp))), axis=(0, 1))
    d_p = (ndimage.distance_transform_edt(~nm_p) * upp)[s.ri.values, s.ci.values]
    m = sets["edge20_roi20"]; rec["Rplacebo_all"] = ratio(d_p[vb], u[vb]); rec["Rplacebo_edge20_roi20"] = ratio(d_p[m], u[m])
    # spostamento ottimo di C-2 (pre-registrazione: con FAIL si usa lo spostamento ottimo)
    j = json.load(open(S13 / "roi" / f"{r.archetype}_{r.roi_id}_primary.json"))
    dx, dy = j["C2_best_dx_um"], j["C2_best_dy_um"]
    rr = np.clip(np.round(s.ri.values + dy / upp).astype(int), 0, nm.shape[0] - 1); cc = np.clip(np.round(s.ci.values + dx / upp).astype(int), 0, nm.shape[1] - 1)
    dsh = dmap[rr, cc]; rec.update({"C2_dx": dx, "C2_dy": dy, "C2_shift_bins": j["C2_max_abs_shift_bins"], "R_C2shift_all": ratio(dsh[vb], u[vb]),
                                    "R_C2shift_edge20_roi20": ratio(dsh[m], u[m]), "R_profile_json": j["by_k"]["3"]["R_profile"]})
    return rec

def edge():
    from multiprocessing import Pool
    rois = E.rois_table()
    with Pool(10) as p: res = p.map(edge_one, list(rois.iterrows()))
    pd.DataFrame(res).to_csv(OUTD / "edge_R.csv", index=False)

# ---------------------------------------------------------------- (b)(c) maschera R2 e lumi
def mask_one(args):
    i, r = args
    import tifffile
    from skimage import color
    A = r.archetype
    s = E.roi_bins(ctx(r), r)
    z = np.load(E.R2M / f"{A}_{r.roi_id}_masks.npz"); t4, b4 = z["tissue_ds4"], z["bubble_ds4"]; px4 = r.side_um / t4.shape[1]
    yy = np.clip((s.y / px4).astype(int), 0, t4.shape[0] - 1).values; xx = np.clip((s.x / px4).astype(int), 0, t4.shape[1] - 1).values
    dt = (ndimage.distance_transform_edt(~t4) * px4)[yy, xx]             # distanza dal tessuto R2 (0 dentro)
    it = s.in_tissue.values == 1; T, B = t4[yy, xx], b4[yy, xx]
    rgb = tifffile.imread(f"/mnt/micron/geo_spatialtrans/R2/rois/{A}_{r.roi_id}.tif")[..., :3]
    f = 4; h, w = (rgb.shape[0] // f) * f, (rgb.shape[1] // f) * f
    rgb4 = rgb[:h, :w].reshape(h // f, f, w // f, f, 3).mean(axis=(1, 3)) / 255.0
    hsv = color.rgb2hsv(rgb4); hed = color.rgb2hed(np.clip(rgb4, 1e-3, 1))
    y4 = np.clip(yy, 0, rgb4.shape[0] - 1); x4 = np.clip(xx, 0, rgb4.shape[1] - 1)
    sat, val, hem, eos = hsv[y4, x4, 1], hsv[y4, x4, 2], hed[y4, x4, 0], hed[y4, x4, 1]
    nucs = {}
    for meth in ("cellpose_rgb", "stardist_he"):
        try:
            lab = np.load(E.R2M / f"{A}_{r.roi_id}_{meth}_native.npz")["labels"]; nucs[meth] = lab[s.ri.values, s.ci.values] > 0
        except FileNotFoundError: nucs[meth] = np.zeros(len(s), bool)
    u = s.umi.values.astype(float)
    cats = {"tissue_valid": T & ~B, "lume_CP2": ~T & ~B & it, "lume_near<=8um": ~T & ~B & it & (dt <= 8), "lume_far>8um": ~T & ~B & it & (dt > 8),
            "lume_far>16um": ~T & ~B & it & (dt > 16), "bubble_in_tissue1": B & it, "off_tissue0": ~it, "off_tissue0_out_t4": ~it & ~T}
    rows = []
    for c, m in cats.items():
        n = int(m.sum())
        rec = {"archetype": A, "roi_id": r.roi_id, "cat": c, "n": n, "frac_roi": n / len(s)}
        if n:
            rec.update({"umi_mean": u[m].mean(), "umi_median": float(np.median(u[m])), "frac_umi0": float((u[m] == 0).mean()),
                        "sat_median": float(np.median(sat[m])), "val_median": float(np.median(val[m])), "hem_median": float(np.median(hem[m])),
                        "eos_median": float(np.median(eos[m])), "nuc_cp": float(nucs["cellpose_rgb"][m].mean()), "nuc_sd": float(nucs["stardist_he"][m].mean())})
        rows.append(rec)
    # profilo per fascia di distanza dal tessuto: lumi (in_tissue=1) vs esterno (in_tissue=0)
    prof = []
    for lab, (lo, hi) in {"0-2": (0, 2), "2-4": (2, 4), "4-8": (4, 8), "8-16": (8, 16), "16-32": (16, 32), ">32": (32, 1e9)}.items():
        for c, base in (("lume", ~T & ~B & it), ("esterno", ~T & ~B & ~it)):
            m = base & (dt > lo) & (dt <= hi)
            prof.append({"archetype": A, "roi_id": r.roi_id, "band": lab, "cat": c, "n": int(m.sum()), "umi_mean": u[m].mean() if m.sum() else np.nan,
                         "sat_median": float(np.median(sat[m])) if m.sum() else np.nan, "eos_median": float(np.median(eos[m])) if m.sum() else np.nan})
    m = cats["lume_far>8um"]
    rho = float(stats.spearmanr(u[m], eos[m])[0]) if m.sum() > 50 else np.nan
    rho_s = float(stats.spearmanr(u[m], sat[m])[0]) if m.sum() > 50 else np.nan
    rows.append({"archetype": A, "roi_id": r.roi_id, "cat": "spearman_lumefar_umi_vs_eos", "n": int(m.sum()), "umi_mean": rho, "sat_median": rho_s})
    return rows, prof

def mask():
    from multiprocessing import Pool
    rois = E.rois_table()
    with Pool(10) as p: res = p.map(mask_one, list(rois.iterrows()))
    pd.DataFrame([x for a, b in res for x in a]).to_csv(OUTD / "mask_categories.csv", index=False)
    pd.DataFrame([x for a, b in res for x in b]).to_csv(OUTD / "lumen_vs_exterior_profile.csv", index=False)

# ---------------------------------------------------------------- C-2 su finestra larga (+-30 um)
def reg_one(args):
    i, r = args
    A = r.archetype; _, mask_m = E.PRIMARY[A]
    s, valid4, px4, vb, nm, dmap, dnuc = E.load_roi_context(r, ctx(r), mask_m)
    upp = r.um_per_px; hh = int(round(1.0 / upp)); u = s.umi.values.astype(float)
    S = np.pad(np.cumsum(np.cumsum(nm.astype(np.int32), 0), 1), ((1, 0), (1, 0)))
    out = []
    for dx in range(-30, 31, 2):
        for dy in range(-30, 31, 2):
            rr = np.round(s.ri.values + dy / upp).astype(int); cc = np.round(s.ci.values + dx / upp).astype(int)
            ok = (rr - hh >= 0) & (cc - hh >= 0) & (rr + hh <= nm.shape[0]) & (cc + hh <= nm.shape[1]) & vb
            nf = (S[rr[ok] + hh, cc[ok] + hh] - S[rr[ok] - hh, cc[ok] + hh] - S[rr[ok] + hh, cc[ok] - hh] + S[rr[ok] - hh, cc[ok] - hh]) / (2 * hh) ** 2
            out.append((dx, dy, float(np.corrcoef(nf, u[ok])[0, 1]) if nf.std() > 0 else np.nan))
    g = pd.DataFrame(out, columns=["dx", "dy", "corr"]); g["a"] = g["corr"].abs()
    b = g.loc[g.a.idxmax()]; n10 = g[(g.dx.abs() <= 10) & (g.dy.abs() <= 10)]; b10 = n10.loc[n10.a.idxmax()]
    return {"archetype": A, "roi_id": r.roi_id, "best30_dx": int(b.dx), "best30_dy": int(b.dy), "abs30": b.a, "best10_dx": int(b10.dx), "best10_dy": int(b10.dy),
            "abs10": b10.a, "abs0": float(g[(g.dx == 0) & (g.dy == 0)].a.iloc[0]), "corr0_signed": float(g[(g.dx == 0) & (g.dy == 0)]["corr"].iloc[0])}

def reg():
    from multiprocessing import Pool
    rois = E.rois_table()
    with Pool(10) as p: res = p.map(reg_one, list(rois.iterrows()))
    pd.DataFrame(res).to_csv(OUTD / "C2_wide.csv", index=False)

# ---------------------------------------------------------------- ricalcolo indipendente di phi/E (senza window_table/permute_within/analyse)
def indep():
    rois = E.rois_table(); rows = []
    for _, r in rois[rois.roi_id.isin(["r1"]) & rois.archetype.isin(["A3", "A6", "A2"])].iterrows():
        s, valid4, px4, vb, nm, dmap, dnuc, terr, ng = setup(r)
        thr = json.load(open(S13 / "bins" / f"{r.dataset}_bg.json"))["bg"]["com"]["3"]["thr"]
        u = s.umi.values.astype(np.int64); ar, ac = s.array_row.values, s.array_col.values
        key = (ar // 3) * 100000 + (ac // 3)
        uk, inv = np.unique(key, return_inverse=True)
        nall = np.bincount(inv); nval = np.bincount(inv, vb.astype(int)); comp = (nall == 9) & (nval == 9)
        def phi(vals): return float((np.bincount(inv, vals)[comp] <= thr).mean())
        po = phi(u); ph = []
        for seed in range(1, 21):
            g = np.random.default_rng(1000 + seed); u2 = u.copy(); idx = np.where(vb)[0]
            order = np.argsort(terr[idx], kind="stable"); tt = terr[idx][order]
            cuts = np.flatnonzero(np.diff(tt)) + 1
            for blk in np.split(order, cuts):
                u2[idx[blk]] = u[idx[g.permutation(blk)]]
            ph.append(phi(u2))
        j = json.load(open(S13 / "roi" / f"{r.archetype}_{r.roi_id}_primary.json"))["by_k"]["3"]
        rows.append({"roi": f"{r.archetype}_{r.roi_id}", "n_win": int(comp.sum()), "n_win_json": j["n_tissue_windows"], "phi_obs": po, "phi_obs_json": j["phi_obs"],
                     "phi_null(semi diversi)": float(np.mean(ph)), "phi_null_json": j["phi_null"], "E": po - np.mean(ph), "E_json": j["E"]})
    pd.DataFrame(rows).to_csv(OUTD / "indep_E.csv", index=False)

# ---------------------------------------------------------------- numeri del testo dai CSV
def claims():
    RS = ROOT / "results/S1.3c"; out = {}
    P = pd.read_csv(RS / "S1.3c_PH3_calibration.csv"); P["rec"] = P.E / P.f_zeroed
    M = P.groupby(["archetype", "type", "f_target"]).rec.median().reset_index()
    out["PH3_rec_min"] = M.loc[M.rec.idxmin()].to_dict(); out["PH3_rec_max"] = M.loc[M.rec.idxmax()].to_dict()
    out["PH3_rec_by_arch_minmax"] = M.groupby("archetype").rec.agg(["min", "max"]).round(3).to_dict()
    pr = pd.read_csv(RS / "S1.3c_PH1_profile.csv"); lv = pd.read_csv(RS / "S1.3c_PH1_levels.csv")
    pr = pr.merge(lv[["archetype", "roi_id", "umi_bg_far_per_bin"]], on=["archetype", "roi_id"])
    pr["x_bg"] = pr.umi_mean / pr.umi_bg_far_per_bin
    q = pr[pr.band == "16-32"].groupby("archetype").x_bg.agg(["median", "min", "max"]).round(2)
    out["PH1_16_32_times_bg"] = q.to_dict()
    out["PH1_profile_by_arch_median"] = pr.groupby(["archetype", "band"]).umi_mean.median().unstack().round(2).to_dict()
    cp = pd.read_csv(RS / "S1.3c_cp_roi.csv"); c1 = cp[(cp.cp == "CP-1") & (cp.k == 3)]
    out["CP1_R_placebo_in_0.9_1.1_by_arch"] = c1.assign(ok=c1.R_placebo.between(0.9, 1.1)).groupby("archetype").ok.sum().to_dict()
    out["CP1_BA_nan_or_0.5"] = c1.groupby("archetype").apply(lambda g: int((g.BA_real.isna() | (g.BA_real.round(3) == 0.5)).sum())).to_dict()
    out["CP1_frac_changed_min"] = float(c1.frac_changed.min())
    out["CP1_R_real_range_A1A5A6"] = c1[c1.archetype.isin(["A1", "A5", "A6"])].groupby("archetype").R_real.agg(["min", "median", "max"]).round(3).to_dict()
    out["CP1_R_placebo_range"] = c1.groupby("archetype").R_placebo.agg(["min", "median", "max"]).round(3).to_dict()
    ro = pd.read_csv(RS / "S1.3c_roi.csv"); pk = ro[(ro.variant == "primary") & (ro.k == 3)]
    out["Emax_primary_k3"] = pk.loc[pk.E.idxmax(), ["archetype", "roi_id", "E"]].to_dict(); out["Emax_all_variants"] = ro.loc[ro.E.idxmax(), ["variant", "k", "archetype", "roi_id", "E"]].to_dict()
    out["n_roi_E_ge_0.05"] = int((ro.E >= 0.05).sum()); out["n_rows_roi_csv"] = len(ro)
    out["BA_A4_primary_k3"] = pk[pk.archetype == "A4"].BA.round(3).tolist()
    Path(OUTD / "claims.json").write_text(json.dumps(out, indent=1, default=str))

if __name__ == "__main__":
    for st, fn in (("claims", claims), ("indep", indep), ("edge", edge), ("mask", mask), ("reg", reg)):
        if ST in (st, "all"): fn(); print("done", st, flush=True)
