"""tools/S1.3c_posthoc.py — analisi post hoc PH-1..PH-4 di S1.3c (dichiarate in results/S1.3c/S1.3c_dev_notes.md).
Uso: .venv/bin/python tools/S1.3c_posthoc.py <ph1|ph2|ph3|ph4|all>   (uscite in /mnt/micron/geo_spatialtrans/S1.3c/posthoc)"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
import importlib.util
ST = sys.argv[1]; sys.argv = sys.argv[:1]
spec = importlib.util.spec_from_file_location("e", Path(__file__).resolve().parent / "S1.3c_ecm.py"); E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
PO = E.OUT / "posthoc"; PO.mkdir(exist_ok=True)
BANDS = {"0-2": (0, 2), "2-4": (2, 4), "4-8": (4, 8), "8-16": (8, 16), "16-32": (16, 32), ">32": (32, 1e9)}

def ctx(r):
    b = pd.read_parquet(E.OUT / "bins" / f"{r.dataset}_roi_bins.parquet"); b["umi"] = b.umi_com
    return b

def lumen_bins(r, s):
    z = np.load(E.R2M / f"{r.archetype}_{r.roi_id}_masks.npz"); t4, b4 = z["tissue_ds4"], z["bubble_ds4"]
    px4 = r.side_um / t4.shape[1]
    yy = np.clip((s.y / px4).astype(int), 0, t4.shape[0] - 1); xx = np.clip((s.x / px4).astype(int), 0, t4.shape[1] - 1)
    dt = ndimage.distance_transform_edt(~t4) * px4
    lum = (~t4[yy, xx]) & (~b4[yy, xx]) & (s.in_tissue.values == 1); tis = t4[yy, xx] & (~b4[yy, xx])
    return lum, tis, dt[yy, xx], t4, px4

def ph1():
    rows, prof = [], []
    for _, r in E.rois_table().iterrows():
        s = E.roi_bins(ctx(r), r); lum, tis, d, _, _ = lumen_bins(r, s); u = s.umi.values; off = s.in_tissue.values == 0
        for lab, (lo, hi) in BANDS.items():
            m = lum & (d > lo) & (d <= hi)
            prof.append({"archetype": r.archetype, "roi_id": r.roi_id, "band": lab, "n": int(m.sum()), "umi_mean": float(u[m].mean()) if m.sum() else np.nan})
        bg = json.load(open(E.OUT / "bins" / f"{r.dataset}_bg.json"))["bg"]["com"]["3"]["mean"] / 9
        rows.append({"archetype": r.archetype, "roi_id": r.roi_id, "n_tissue_bins": int(tis.sum()), "umi_tissue_mean": float(u[tis].mean()),
                     "n_lumen_bins": int(lum.sum()), "umi_lumen_mean": float(u[lum].mean()) if lum.sum() else np.nan,
                     "n_off_bins": int(off.sum()), "umi_off_mean": float(u[off].mean()) if off.sum() else np.nan,
                     "lumen_frac_roi": float(lum.mean()), "umi_bg_far_per_bin": bg})
    pd.DataFrame(rows).to_csv(PO / "PH1_levels.csv", index=False); pd.DataFrame(prof).to_csv(PO / "PH1_profile.csv", index=False)

def ph2_one(r):
    A = r.archetype; _, mask_m = E.PRIMARY[A]
    s, valid4, px4, vb, nm, dmap, dnuc = E.load_roi_context(r, ctx(r), mask_m)
    upp = r.um_per_px; h = int(round(1.0 / upp)); u = s.umi.values.astype(float)
    half = ((np.floor(s.x / 100) + np.floor(s.y / 100)) % 2).values.astype(int)
    rows = []
    for dx in range(-10, 11, 2):
        for dy in range(-10, 11, 2):
            nf = E.nuc_fraction_shift(nm, s.ri.values, s.ci.values, dx / upp, dy / upp, h)
            rec = {"archetype": A, "roi_id": r.roi_id, "dx": dx, "dy": dy}
            for hh, sel in (("all", np.ones(len(s), bool)), ("h0", half == 0), ("h1", half == 1)):
                ok = vb & np.isfinite(nf) & sel
                rec[f"corr_{hh}"] = float(np.corrcoef(nf[ok], u[ok])[0, 1]) if ok.sum() > 10 and nf[ok].std() > 0 else np.nan
            rows.append(rec)
    return rows

def ph2():
    from multiprocessing import Pool
    rois = E.rois_table()
    with Pool(15) as p: res = p.map(ph2_one, [r for _, r in rois.iterrows()])
    S = pd.DataFrame([x for rr in res for x in rr]); S.to_csv(PO / "PH2_surfaces.csv", index=False)
    out = []
    for (A, roi), g in S.groupby(["archetype", "roi_id"], sort=False):
        rec = {"archetype": A, "roi_id": roi}
        for hh in ("all", "h0", "h1"):
            i = g[f"corr_{hh}"].abs().idxmax(); rec[f"best_{hh}"] = (int(g.dx[i]), int(g.dy[i])); rec[f"abs_{hh}"] = float(abs(g[f"corr_{hh}"][i]))
            rec[f"abs0_{hh}"] = float(abs(g[(g.dx == 0) & (g.dy == 0)][f"corr_{hh}"].iloc[0]))
        b0, b1 = np.array(rec["best_h0"]), np.array(rec["best_h1"])
        rec["halves_dist_um"] = float(np.hypot(*(b0 - b1))); rec["gain_all"] = rec["abs_all"] / max(rec["abs0_all"], 1e-9)
        out.append(rec)
    pd.DataFrame(out).to_csv(PO / "PH2_summary.csv", index=False)

def ph3_one(r):
    A = r.archetype; terr_m, mask_m = E.PRIMARY[A]; k = 3
    thr = json.load(open(E.OUT / "bins" / f"{r.dataset}_bg.json"))["bg"]["com"]["3"]["thr"]
    s, valid4, px4, vb, nm, dmap, dnuc = E.load_roi_context(r, ctx(r), mask_m)
    gx, gy, _ = E.generators(A, r.roi_id, terr_m, valid4, px4)
    terr = cKDTree(np.c_[gx, gy]).query(np.c_[s.x.values, s.y.values])[1]
    tree = cKDTree(np.c_[s.x.values, s.y.values]); R = np.sqrt(3 * E.T_MED[A] / np.pi)
    rows = []
    for f in (0.05, 0.10, 0.20, 0.30):
        rq = float(np.percentile(dnuc[vb], 100 * (1 - f))); zA = vb & (dnuc > rq)
        rng = np.random.default_rng(20261012); inside = np.zeros(len(s), bool)
        while (inside & vb).sum() < f * vb.sum():
            cx, cy = rng.uniform(0, r.side_um, 2)
            if not valid4[min(int(cy / px4), valid4.shape[0] - 1), min(int(cx / px4), valid4.shape[1] - 1)]: continue
            inside[tree.query_ball_point([cx, cy], R)] = True
        zB = vb & inside
        for typ, z in (("A", zA), ("B", zB)):
            s2 = s.copy(); s2.loc[z, "umi"] = 0
            o, _ = E.analyse(s2, vb, dnuc, terr, thr, k, A)
            rows.append({"archetype": A, "roi_id": r.roi_id, "type": typ, "f_target": f, "f_zeroed": float(z.sum() / vb.sum()),
                         **{kk: o[kk] for kk in ("phi_obs", "phi_null", "E", "BA", "large_frac", "verdict_roi")}})
    return rows

def ph3():
    from multiprocessing import Pool
    with Pool(15) as p: res = p.map(ph3_one, [r for _, r in E.rois_table().iterrows()])
    C = pd.DataFrame([x for rr in res for x in rr]); C.to_csv(PO / "PH3_calibration.csv", index=False)
    M = C.groupby(["archetype", "type", "f_target"]).E.median().reset_index()
    mins = M[M.E >= 0.05].groupby(["archetype", "type"]).f_target.min().rename("f_min_detectable").reset_index()
    mins.to_csv(PO / "PH3_min_detectable.csv", index=False)

def ph4():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import tifffile
    rois = E.rois_table(); sel = pd.concat([rois[rois.archetype == "A5"], rois[rois.archetype != "A5"].groupby("archetype").head(1)])
    for _, r in sel.iterrows():
        s = E.roi_bins(ctx(r), r); lum, tis, d, t4, px4 = lumen_bins(r, s)
        W = 250.0; gx = np.floor(s.x[lum] / W); gy = np.floor(s.y[lum] / W)
        if lum.sum() == 0: continue
        cnt = pd.Series(1, index=pd.MultiIndex.from_arrays([gx, gy])).groupby(level=[0, 1]).size().idxmax()
        x0, y0 = cnt[0] * W, cnt[1] * W
        rgb = tifffile.imread(f"/mnt/micron/geo_spatialtrans/R2/rois/{r.archetype}_{r.roi_id}.tif")
        c0, r0 = int(x0 / r.um_per_px), int(y0 / r.um_per_px); n = int(W / r.um_per_px)
        crop = rgb[r0:r0 + n, c0:c0 + n]
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), facecolor="white")
        for ax in axes: ax.imshow(crop, extent=(x0, x0 + W, y0 + W, y0)); ax.set_xlim(x0, x0 + W); ax.set_ylim(y0 + W, y0)
        tm = t4[int(y0 / px4):int((y0 + W) / px4), int(x0 / px4):int((x0 + W) / px4)].astype(float)
        axes[1].contour(np.linspace(x0, x0 + W, tm.shape[1]), np.linspace(y0, y0 + W, tm.shape[0]), tm, levels=[0.5], colors="cyan", linewidths=0.8)
        m = lum & (s.x >= x0).values & (s.x < x0 + W).values & (s.y >= y0).values & (s.y < y0 + W).values
        sc = axes[1].scatter(s.x[m], s.y[m], c=np.log1p(s.umi[m]), s=2, cmap="magma", vmin=0, vmax=np.log1p(60), marker="s")
        plt.colorbar(sc, ax=axes[1], fraction=0.04, label="log(1+UMI) bin «lume»")
        axes[0].set_title(f"{r.archetype} {r.roi_id} — H&E (x {x0:.0f}–{x0+W:.0f} µm, y {y0:.0f}–{y0+W:.0f} µm)", fontsize=9)
        axes[1].set_title(f"contorno maschera R2 (ciano); bin in_tissue=1 fuori maschera: {int(m.sum())}", fontsize=9)
        fig.tight_layout(); fig.savefig(PO / f"PH4_{r.archetype}_{r.roi_id}.png", dpi=110); plt.close(fig)

if __name__ == "__main__":
    for st, fn in (("ph1", ph1), ("ph2", ph2), ("ph3", ph3), ("ph4", ph4)):
        if ST in (st, "all"): fn(); print("done", st, flush=True)
