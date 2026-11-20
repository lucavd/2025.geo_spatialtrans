"""tools/riesame_xen_export.py - Riesame Parte B, bracci X e X1 (docs/riesame/REGOLA_B.md; scelte XD-1...XD-12 di docs/riesame/DECISIONI_B_XENIUM.md).
Esportazioni nel formato di tools/F3.3_export.py (aree di verifica): X.parquet (bin, gene, n; 0-based), bins.parquet (quadratini di V di F3.1,
label = nucleo Xenium sotto il centro del quadratino), cells.parquet (label, tipi vuoti, xenium_cell_id), genes.tsv (408 Ensembl), export.json.
Fasi: cache  -> trascritti qv >= 20 dei 408 geni, portati in um del riquadro di ogni area (finestra allargata di 100 um), un file per area;
      X      -> braccio X;   X1 -> diradamento binomiale unico per tessuto verso la media di V (results/riesame/xen_targets.csv).
Uso: .venv/bin/python -B tools/riesame_xen_export.py <cache|X|X1> [aree separate da virgola]"""
import sys, os, json, time, numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq, shapely
from shapely import STRtree
M = "/mnt/micron/geo_spatialtrans"; XP = f"{M}/data_xenium/P2/outs"; OUTB = f"{M}/riesame/B"; F = f"{M}/F3.1/main/real"; CACHE = f"{OUTB}/_cache"
AREAS = [f"A2_r{i}" for i in range(1, 6)] + [f"A3_r{i}" for i in range(1, 6)]
QV_MIN = 20.0; MARGIN_UM = 100.0; SEED_X1 = 20261124
PANEL = pd.read_csv(f"{M}/riesame/xenium/panel_ids_in_common.txt", header=None)[0].astype(str).tolist()
AR = pd.read_csv("results/F3.3/F3.3_areas.csv").set_index("area")
_REG = {}
def REGJ():
    if not _REG: _REG.update(json.load(open(f"{OUTB}/registration.json")))
    return _REG
def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)
def panel_map():
    """simbolo Xenium -> indice nei 408 (via identificativo Ensembl del pannello)."""
    gp = json.load(open(f"{XP}/gene_panel.json")); m = {}
    for t in gp["payload"]["targets"]:
        d = t["type"]["data"]
        if t["type"].get("descriptor") == "gene" and d.get("id") in PANEL: m[d["name"]] = PANEL.index(d["id"])
    return m
def to_um(area, x, y):
    """coordinate Xenium (um) -> um del riquadro dell'area (x = colonna, y = riga), come x_um, y_um dei quadratini di V."""
    B = np.asarray(REGJ()["areas"][area]["B"]); r = AR.loc[area]
    col, row = (np.c_[x, y, np.ones(len(x))] @ B).T
    return (col - r.c0) * r.upp, (row - r.r0) * r.upp
def side(area): r = AR.loc[area]; return float(r.side_px * r.upp)

def stage_cache(areas):
    os.makedirs(CACHE, exist_ok=True); pm = panel_map(); assert len(set(pm.values())) == 408, len(set(pm.values()))
    t = pq.read_table(f"{XP}/transcripts.parquet", columns=["feature_name", "x_location", "y_location", "qv", "transcript_id"])
    fn = t.column("feature_name").to_pandas().astype(str).to_numpy(); qv = t.column("qv").to_numpy(); x = t.column("x_location").to_numpy(); y = t.column("y_location").to_numpy()
    tid = t.column("transcript_id").to_numpy(); del t
    gi = pd.Series(fn).map(pm).to_numpy(); inp = ~pd.isna(gi)
    stats = {"n_transcripts_all": int(len(fn)), "n_panel408": int(inp.sum()), "n_panel408_qv20": int((inp & (qv >= QV_MIN)).sum()),
             "n_all_qv20": int((qv >= QV_MIN).sum()), "n_feature_names": int(len(np.unique(fn))), "n_symbols_mapped": len(pm)}
    keep = inp & (qv >= QV_MIN); x, y, gi, qv_k, tid = x[keep], y[keep], gi[keep].astype(np.int16), qv[keep], tid[keep]
    for a in areas:
        xu, yu = to_um(a, x, y); s = side(a)
        w = (xu >= -MARGIN_UM) & (xu < s + MARGIN_UM) & (yu >= -MARGIN_UM) & (yu < s + MARGIN_UM)
        pq.write_table(pa.table({"x_um": xu[w], "y_um": yu[w], "gene": gi[w], "transcript_id": tid[w]}), f"{CACHE}/tx_{a}.parquet")
        stats[f"n_window_{a}"] = int(w.sum()); log("cache", a, int(w.sum()))
    # quota qv >= 20 sui 408 geni per area (riquadro), per il controllo
    json.dump(stats, open(f"{CACHE}/tx_stats.json", "w"), indent=1); log(stats)

_nuc = {}
def nuclei():
    if not _nuc:
        ce = pd.read_parquet(f"{XP}/cells.parquet", columns=["cell_id", "x_centroid", "y_centroid"])
        lab = pd.Series(np.arange(1, len(ce) + 1), index=ce.cell_id.astype(str))
        nb = pd.read_parquet(f"{XP}/nucleus_boundaries.parquet"); cb = pd.read_parquet(f"{XP}/cell_boundaries.parquet")
        _nuc.update(ce=ce, lab=lab, nb=nb, cb=cb)
    return _nuc
def polys(df, area):
    """poligoni (uno per cell_id e label_id) in um del riquadro; ritorna (geometrie, cell_id)."""
    xu, yu = to_um(area, df.vertex_x.to_numpy(float), df.vertex_y.to_numpy(float)); s = side(area)
    key = df.cell_id.astype(str).to_numpy() + "|" + df.label_id.astype(str).to_numpy() if "label_id" in df else df.cell_id.astype(str).to_numpy()
    k, first, inv = np.unique(key, return_index=True, return_inverse=True)
    cx = np.bincount(inv, weights=xu) / np.bincount(inv); cy = np.bincount(inv, weights=yu) / np.bincount(inv)
    sel = (cx > -60) & (cx < s + 60) & (cy > -60) & (cy < s + 60)
    m = sel[inv]; xu, yu, inv2 = xu[m], yu[m], inv[m]
    o = np.argsort(inv2, kind="stable"); xu, yu, inv2 = xu[o], yu[o], inv2[o]
    _, ring = np.unique(inv2, return_inverse=True)
    lr = shapely.linearrings(np.c_[xu, yu], indices=ring); g = shapely.make_valid(shapely.polygons(lr))
    cid = df.cell_id.astype(str).to_numpy()[first][np.flatnonzero(sel)]
    return g, cid

def v_area(a):
    d = f"{F}/{a}"; b = pd.read_parquet(f"{d}/bins.parquet"); return b
def grid_lut(b):
    """affine (row, col) -> (y_um, x_um) dai quadratini di V; inversa per assegnare un punto al quadratino."""
    D = np.c_[b.row.to_numpy(float), b.col.to_numpy(float), np.ones(len(b))]
    Bm, *_ = np.linalg.lstsq(D, np.c_[b.x_um.to_numpy(), b.y_um.to_numpy()], rcond=None)
    res = np.hypot(*(np.c_[b.x_um, b.y_um] - D @ Bm).T)
    L = Bm[:2].T; t = Bm[2]; Li = np.linalg.inv(L)
    r0, c0 = int(b.row.min()), int(b.col.min()); lut = -np.ones((int(b.row.max()) - r0 + 1, int(b.col.max()) - c0 + 1), np.int64)
    lut[b.row.to_numpy() - r0, b.col.to_numpy() - c0] = np.arange(len(b))
    def assign(x, y):
        rc = (np.c_[x, y] - t) @ Li.T; r = np.rint(rc[:, 0]).astype(np.int64) - r0; c = np.rint(rc[:, 1]).astype(np.int64) - c0
        ok = (r >= 0) & (r < lut.shape[0]) & (c >= 0) & (c < lut.shape[1]); out = -np.ones(len(x), np.int64); out[ok] = lut[r[ok], c[ok]]
        return out
    return assign, {"grid_affine": Bm.tolist(), "grid_fit_resid_max_um": float(res.max()), "bin_pitch_um": float(np.sqrt(abs(np.linalg.det(L))))}

def write_export(out, b, X, cells, info):
    os.makedirs(out, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(X[["bin", "gene", "n"]].astype({"bin": np.int64, "gene": np.int64, "n": np.int32}), preserve_index=False), f"{out}/X.parquet")
    b.to_parquet(f"{out}/bins.parquet", index=False); cells.to_parquet(f"{out}/cells.parquet", index=False)
    pd.Series(PANEL).to_csv(f"{out}/genes.tsv", index=False, header=False)
    u = np.bincount(X.bin.to_numpy(), weights=X.n.to_numpy(), minlength=len(b)); reg = b.in_reg.to_numpy()
    info.update({"n_bins": int(len(b)), "n_in_filtered": int(b.in_filtered.sum()), "n_reg": int(reg.sum()), "n_reg_not_filtered": int((reg & ~b.in_filtered.to_numpy()).sum()),
                 "umi": int(X.n.sum()), "umi_reg": int(u[reg].sum()), "mean_umi_reg": float(u[reg].mean()), "n_genes": len(PANEL),
                 "n_lab_bins": int((b.label > 0).sum()), "n_lab_bins_dimg": int((b.label_dimg > 0).sum()), "n_nuclei_keep": int(len(cells))})
    json.dump(info, open(f"{out}/export.json", "w")); log(json.dumps({k: v for k, v in info.items() if k != "grid_affine"}))

def stage_X(areas):
    N = nuclei()
    for a in areas:
        t0 = time.time(); r = AR.loc[a]; b = v_area(a); assign, gi = grid_lut(b)
        tx = pd.read_parquet(f"{CACHE}/tx_{a}.parquet"); bi = assign(tx.x_um.to_numpy(), tx.y_um.to_numpy()); ok = bi >= 0
        X = pd.DataFrame({"bin": bi[ok], "gene": tx.gene.to_numpy()[ok].astype(np.int64)}).groupby(["bin", "gene"]).size().reset_index(name="n")
        # nuclei sotto il centro del quadratino
        g, cid = polys(N["nb"], a); tree = STRtree(g)
        pts = shapely.points(b.x_um.to_numpy(), b.y_um.to_numpy()); pi, gj = tree.query(pts, predicate="within")
        lab = np.zeros(len(b), np.int64); labs = N["lab"].reindex(cid[gj]).to_numpy().astype(np.int64)
        o = np.lexsort((labs, pi)); pi, labs = pi[o], labs[o]; first = np.r_[True, pi[1:] != pi[:-1]]; lab[pi[first]] = labs[first]
        n_multi = int((~first).sum())
        # cellule Xenium (per il controllo della quota di tessuto, XD-6)
        gc, _ = polys(N["cb"], a); pc, _ = STRtree(gc).query(pts, predicate="within"); in_cell = np.zeros(len(b), bool); in_cell[pc] = True
        bb = b.copy(); bb["label"] = lab; bb["label_dimg"] = lab
        ul = np.unique(lab[lab > 0]); inv = pd.Series(N["ce"].cell_id.astype(str).to_numpy(), index=np.arange(1, len(N["ce"]) + 1))
        cells = pd.DataFrame({"label": ul}); cells["cell_type"] = None; cells["classifiable"] = None; cells["cell_type_r6"] = None
        cells["xenium_cell_id"] = inv.reindex(ul).to_numpy()
        # centroidi dei nuclei nel riquadro (densita', XD-6)
        cxy = N["ce"]; cxu, cyu = to_um(a, cxy.x_centroid.to_numpy(float), cxy.y_centroid.to_numpy(float)); s = side(a)
        has_nuc = N["lab"].index.isin(np.unique(N["nb"].cell_id.astype(str)))
        n_c_box = int(((cxu >= 0) & (cxu < s) & (cyu >= 0) & (cyu < s) & has_nuc).sum())
        tis = b.in_tissue.to_numpy(); reg = b.in_reg.to_numpy()
        info = {"tag": a, "area": a, "archetype": a[:2], "level": "XEN", "section": "P2_xenium", "genes": "panel408", "arm": "X", "side_um": s, "d_img_um": [0.0, 0.0],
                "source": "GSE280314 P2 Xenium, qv>=20", "registration_mode": REGJ()["mode"], **gi,
                "n_tx_window": int(len(tx)), "n_tx_in_box_bins": int(ok.sum()), "n_bins_multi_nucleus": n_multi,
                "n_xen_nuclei_centroid_in_box": n_c_box, "xen_nuclei_per_mm2": n_c_box / (s / 1000) ** 2,
                "share_bins_in_xen_cell_all": float(in_cell.mean()), "share_bins_in_xen_cell_in_tissue": float(in_cell[tis].mean()),
                "share_bins_in_xen_cell_in_reg": float(in_cell[reg].mean()), "share_in_tissue_v": float(tis.mean()), "share_in_reg_v": float(reg.mean()),
                "sec": None}
        np.save(f"{CACHE}/in_cell_{a}.npy", in_cell)
        info["sec"] = round(time.time() - t0, 1); write_export(f"{OUTB}/X/{a}", bb, X, cells, info)

def stage_X1(areas):
    tg = pd.read_csv("results/riesame/xen_targets.csv").set_index("area")
    ex = {a: json.load(open(f"{OUTB}/X/{a}/export.json")) for a in AREAS}
    f = {}
    for t in ("A2", "A3"):
        aa = [a for a in AREAS if a.startswith(t)]
        f[t] = float(tg.loc[aa, "v_mean408_reg"].mean() / np.mean([ex[a]["mean_umi_reg"] for a in aa]))
    json.dump({"factor": f, "seed_base": SEED_X1, "target": {t: float(tg.loc[[a for a in AREAS if a.startswith(t)], "v_mean408_reg"].mean()) for t in f}},
              open(f"{OUTB}/X1_factors.json", "w"), indent=1); log("X1 factors", f)
    for a in areas:
        t0 = time.time(); ft = f[a[:2]]
        if ft > 1: log("X1 not possible (factor > 1)", a, ft); continue
        d = f"{OUTB}/X/{a}"; X = pd.read_parquet(f"{d}/X.parquet"); b = pd.read_parquet(f"{d}/bins.parquet"); cells = pd.read_parquet(f"{d}/cells.parquet")
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED_X1, AREAS.index(a)])))
        X = X.assign(n=rng.binomial(X.n.to_numpy(), ft)); X = X[X.n > 0].reset_index(drop=True)
        info = {k: v for k, v in ex[a].items() if k not in ("umi", "umi_reg", "mean_umi_reg", "sec")}
        info.update({"tag": a, "arm": "X1", "thin_factor": ft, "seed": [SEED_X1, AREAS.index(a)], "seed_kind": "numpy SeedSequence -> PCG64"})
        info["sec"] = round(time.time() - t0, 1); write_export(f"{OUTB}/X1/{a}", b, X, cells, info)

if __name__ == "__main__":
    st = sys.argv[1]; areas = sys.argv[2].split(",") if len(sys.argv) > 2 else AREAS
    {"cache": stage_cache, "X": stage_X, "X1": stage_X1}[st](areas)
