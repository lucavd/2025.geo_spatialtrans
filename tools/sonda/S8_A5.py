"""tools/sonda/S8_A5.py - S8, A5 (docs/sonda/REGOLA_S8.md): cinque regole di assegnazione dei quadratini (ENACT 0.2.5: naive, weighted_by_area,
weighted_by_gene, weighted_by_cluster, chiamate senza modifiche tramite ENACT.assign_bins_to_cells; bin2cell da S8_b2c.py) sul simulato banale B
(10 aree di H2 x 2 semi) e sul Visium HD di X1 (XVfull, 10 aree); verita' per cellula (B: territori ricostruiti dal seme; X1: trascritti Xenium
qv >= 20 della cellula abbinata), coseno sui 284 geni di G_p, cancelli G-ENACT, G-TRUTH, G-MATCH; riepilogo T5, S5c, A5, CP-shuf.
Uso: .venv-enact/bin/python -B tools/sonda/S8_A5.py unit <arm> <unit> [out_root] | summary | list"""
import sys, os, json, time, logging, glob, numpy as np, pandas as pd, scipy.sparse as sp
import geopandas as gpd, shapely, anndata as ad, scanpy as sc, pyarrow.dataset as pds, pyarrow as pa
from scipy.spatial import cKDTree
from scipy.stats import kendalltau
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C, S5_common as S5
from enact.pipeline import ENACT
M = C.M; OUT = os.environ.get("S8_OUT", f"{M}/S8_sonda"); RES = f"{C.REPO}/results/sonda/S8"; SEED = 20261201; BIN = 2.0
EM = {"E-naive": "naive", "E-area": "weighted_by_area", "E-gene": "weighted_by_gene", "E-cluster": "weighted_by_cluster"}
METHODS = list(EM) + ["b2c"]; UNITS = S5.UNITS
logging.basicConfig(level=logging.WARNING)

def units():
    out = []
    for i, (t, k) in enumerate(UNITS):
        for s in (0, 1): out.append(dict(arm="B", unit=f"H2_{t}_{k}_k{s}", d=f"{M}/S3_sonda/runs/B/H2_{t}_{k}_k{s}", i=i, s=s, tissue=t))
    for i, (t, k) in enumerate(UNITS):
        out.append(dict(arm="X1", unit=f"X1_{t}_{k}", d=f"{M}/S5_sonda/runs/XVfull/X1_{t}_{k}", xb=f"{M}/S5_sonda/runs/XB/X1_{t}_{k}", i=i, s=0, tissue=t))
    return out

def load(d):
    bins = pd.read_parquet(f"{d}/bins.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
    T = pd.read_parquet(f"{d}/X.parquet"); X = sp.csr_matrix((T.n.to_numpy().astype(np.float64), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
    return bins, genes, X

def group_sum(ids, Xrows, labels):
    """somma delle righe di Xrows per etichetta (ids: etichetta di ogni riga), nell'ordine di labels"""
    pos = pd.Series(np.arange(len(labels)), index=labels).reindex(ids).to_numpy()
    k = ~np.isnan(pos); P = sp.csr_matrix((np.ones(k.sum()), (pos[k].astype(np.int64), np.flatnonzero(k))), shape=(len(labels), len(ids)))
    return np.asarray((P @ Xrows).todense())

def nuclei_gdf(bins, tb):
    lab = bins.label.to_numpy(); b = tb[lab[tb] > 0]; L = lab[b]
    boxes = shapely.box(bins.x_um.to_numpy()[b] - 1, bins.y_um.to_numpy()[b] - 1, bins.x_um.to_numpy()[b] + 1, bins.y_um.to_numpy()[b] + 1)
    df = pd.DataFrame({"id": L, "g": boxes}); geoms = df.groupby("id").g.apply(lambda s: shapely.union_all(s.to_numpy()))
    g = gpd.GeoDataFrame({"id": geoms.index.to_numpy()}, geometry=list(geoms.to_numpy())).reset_index(drop=True)
    cen = g.geometry.centroid; g["cell_x"] = cen.x.to_numpy(); g["cell_y"] = cen.y.to_numpy()
    g["id"] = ["c%d" % v for v in g.id.to_numpy()]; return g          # ENACT usa identificativi testuali

def enact_obj(method, wd):
    so = object.__new__(ENACT); so.logger = logging.getLogger("enact"); so.patch_size = 1e12; so.bin_representation = "polygon"
    so.configs = {"params": {"chunks_to_run": [], "n_clusters": 4, "n_pcs": 250}}; so.bin_to_cell_method = method
    so.cell_chunks_dir = f"{wd}/chunks/cells"; so.bin_chunks_dir = f"{wd}/chunks/bins"
    so.bin_assign_dir = f"{wd}/{method}/assign"; so.cell_ix_lookup_dir = f"{wd}/{method}/lookup"; return so

def run_enact(bins, genes, X, tb, gp_ids, wd, seed):
    a = ad.AnnData(X[tb]); a.obs_names = ["b" + str(v) for v in bins.bin_id.to_numpy()[tb]]; a.var_names = genes
    a.obs["index"] = a.obs_names.to_numpy(); a.obs["in_tissue"] = 1
    a.obs["pxl_col_in_fullres"] = bins.x_um.to_numpy()[tb]; a.obs["pxl_row_in_fullres"] = bins.y_um.to_numpy()[tb]
    # geni come load_visiumhd_dataset di ENACT (use_hvg, n_hvg 1000) con G_p come marcatori
    an = a.copy(); sc.pp.normalize_total(an); sc.pp.log1p(an); sc.pp.highly_variable_genes(an, n_top_genes=1000)
    mask = an.var["highly_variable"].copy(); n_hvg = int(mask.sum()); mask.loc[list(set(gp_ids) & set(genes))] = True; a = a[:, mask.to_numpy()].copy()
    nuc = nuclei_gdf(bins, tb); so0 = enact_obj("naive", wd)
    cells = so0.expand_nuclei_with_voronoi(nuc.copy(), expansion_size=2 * BIN)
    so0.split_df_to_chunks(df=cells.copy(), x_col="cell_x", y_col="cell_y", output_dir=so0.cell_chunks_dir)
    bg = so0.generate_bins_gdf(a, BIN); so0.split_df_to_chunks(df=bg, x_col="pxl_col_in_fullres", y_col="pxl_row_in_fullres", output_dir=so0.bin_chunks_dir)
    res = {}
    for m, em in EM.items():
        np.random.seed(seed); so = enact_obj(em, wd); so.assign_bins_to_cells(a, None)
        ch = os.listdir(so.bin_assign_dir); assert ch == ["patch_0_0.csv"], ch
        cnt = pd.read_csv(f"{so.bin_assign_dir}/patch_0_0.csv", index_col=0); lk = pd.read_csv(f"{so.cell_ix_lookup_dir}/patch_0_0.csv", index_col=0)
        assert len(cnt) == len(lk); cnt.index = lk["id"].astype(str).str[1:].astype(np.int64).to_numpy(); res[m] = cnt
    return a, cells, res, n_hvg

def enact_gate(a, cells, res):
    """G-ENACT: totali attesi calcolati a parte (STRtree di shapely, intersezioni vettoriali)"""
    xs, ys = a.obs.pxl_col_in_fullres.to_numpy(), a.obs.pxl_row_in_fullres.to_numpy(); boxes = shapely.box(xs - 1, ys - 1, xs + 1, ys + 1)
    cg = cells.geometry.to_numpy(); tree = shapely.STRtree(cg); bi, ci = tree.query(boxes, predicate="intersects")
    nper = np.bincount(bi, minlength=len(boxes)); rs = np.asarray(a.X.sum(axis=1)).ravel(); inter = nper > 0; uniq = nper == 1
    exp = {"E-naive": rs[uniq].sum(), "E-cluster": rs[inter].sum()}
    sh = nper[bi] > 1; ar = shapely.area(shapely.intersection(boxes[bi[sh]], cg[ci[sh]])) / (BIN * BIN)
    exp["E-area"] = rs[uniq].sum() + (ar * rs[bi[sh]]).sum()
    ids = cells.id.to_numpy(); ub = np.flatnonzero(uniq); cu = np.zeros(len(boxes), np.int64) - 1
    cu[bi[~sh]] = ci[~sh]; U = group_sum(cu[ub], a.X[ub], np.arange(len(cg)))        # celle x geni dai soli quadratini propri
    Xc = a.X.tocsr(); lost = 0.0
    sb = np.unique(bi[sh])
    for b in sb:
        cc = ci[bi == b]; z = (U[cc].sum(axis=0) == 0); row = Xc[b].toarray().ravel(); lost += row[z].sum()
    exp["E-gene"] = rs[inter].sum() - lost
    got = {m: float(res[m].to_numpy().sum()) for m in EM}
    rel = {m: abs(got[m] - exp[m]) / max(exp[m], 1) for m in EM}
    # E-cluster (errore di pre-registrazione, REGOLA_S8 «Correzione del cancello G-ENACT»): una cellula senza quadratini propri non ha profilo;
    # weight_by_cluster le da' peso 1/n e le altre si dividono comunque il quadratino: eccesso <= somma_b conteggi_b * k_b / n_b
    has_u = np.zeros(len(cg), bool); has_u[ci[~sh]] = True; exc_max = 0.0; n_bk = 0
    for b in sb:
        cc = ci[bi == b]; kb = int((~has_u[cc]).sum())
        if 0 < kb < len(cc): exc_max += rs[b] * kb / len(cc); n_bk += 1
    exc = got["E-cluster"] - exp["E-cluster"]; tol = 1e-9 * max(exp["E-cluster"], 1)
    ok_cl = (abs(exc) <= tol) if n_bk == 0 else (-tol <= exc <= exc_max + tol)
    ok = all(rel[m] < 1e-9 for m in ("E-naive", "E-area", "E-gene")) and ok_cl
    return dict(expected={m: float(v) for m, v in exp.items()}, got=got, rel=rel, cluster_excess=float(exc), cluster_excess_max=float(exc_max),
                n_bins_cell_without_own=n_bk, ok=bool(ok), n_shared_bins=int(len(sb)), n_inter_bins=int(inter.sum()))

def truth_B(u, bins, X, tb, gpi):
    ck = json.load(open(f"{u['d']}/B_check.json")); par = json.load(open(f"{C.ROOT}/B/params_{u['tissue']}.json"))
    x, y = bins.x_um.to_numpy(), bins.y_um.to_numpy(); rng = np.random.default_rng(ck["seed"])
    N = int(rng.poisson(par["lambda_um2"] * 4.0 * len(tb))); pick = rng.integers(0, len(tb), N)
    cx = x[tb[pick]] + rng.uniform(-1, 1, N); cy = y[tb[pick]] + rng.uniform(-1, 1, N)
    dist, own = cKDTree(np.c_[cx, cy]).query(np.c_[x[tb], y[tb]]); keep = np.unique(own)
    newid = np.zeros(N, np.int64); newid[keep] = np.arange(1, len(keep) + 1); cell = newid[own]
    lab = np.zeros(len(bins), np.int64); lab[tb[dist <= par["r_um"]]] = cell[dist <= par["r_um"]]
    labels = np.arange(1, len(keep) + 1); Tt = group_sum(cell, X[tb][:, gpi], labels)
    tot = float(X[tb].sum())
    g = dict(label_equal=bool(np.array_equal(lab, bins.label.to_numpy())), N=N, n_cells=int(len(keep)), umi_total=tot, umi_total_check=ck["umi_total"])
    g["ok"] = g["label_equal"] and abs(tot - ck["umi_total"]) < 0.5 and N == ck["n_centres"]
    return labels, Tt, g

def truth_X1(u, bins, gp_ids):
    xb = pd.read_parquet(f"{u['xb']}/bins.parquet"); assert np.array_equal(xb.bin_id.to_numpy(), bins.bin_id.to_numpy())
    s, x = bins.label.to_numpy(), xb.label.to_numpy(); k = (s > 0) & (x > 0)
    pr = pd.DataFrame({"s": s[k], "x": x[k]}).value_counts().rename("n").reset_index()
    ss = pd.Series(s[s > 0]).value_counts(); xs = pd.Series(x[x > 0]).value_counts()
    bs = pr.sort_values(["s", "n", "x"], ascending=[True, False, True]).drop_duplicates("s"); bx = pr.sort_values(["x", "n", "s"], ascending=[True, False, True]).drop_duplicates("x")
    mt = bs.merge(bx, on=["s", "x", "n"]); mt = mt[mt.n >= 0.5 * np.minimum(ss.reindex(mt.s).to_numpy(), xs.reindex(mt.x).to_numpy())]
    cid = pd.read_parquet(f"{u['xb']}/cells.parquet").set_index("label").xenium_cell_id; mt["cell_id"] = cid.reindex(mt.x).to_numpy()
    g = dict(n_stardist=int(len(ss)), n_xenium=int(len(xs)), n_matched=int(len(mt)), one_to_one=bool(mt.x.is_unique and mt.s.is_unique and mt.cell_id.is_unique))
    s2e = S5.sym2ens(); gix = pd.Series(np.arange(len(gp_ids)), index=gp_ids)
    dset = pds.dataset(f"{S5.XPRE}_transcripts.parquet"); typ = dset.schema.field("cell_id").type
    want = [c.encode() for c in mt.cell_id] if pa.types.is_binary(typ) else list(mt.cell_id)
    t = dset.to_table(columns=["cell_id", "feature_name", "qv"], filter=(pds.field("cell_id").isin(want)) & (pds.field("qv") >= 20)).to_pandas()
    for col in ("cell_id", "feature_name"):
        if len(t) and isinstance(t[col].iloc[0], bytes): t[col] = t[col].str.decode("utf-8")
    gi = t.feature_name.map(s2e).map(gix); t = t[gi.notna().to_numpy()]; gi = gi[gi.notna()].to_numpy().astype(np.int64)
    rowpos = pd.Series(np.arange(len(mt)), index=mt.cell_id.to_numpy()).reindex(t.cell_id.to_numpy()).to_numpy().astype(np.int64)
    Tt = np.zeros((len(mt), len(gp_ids))); np.add.at(Tt, (rowpos, gi), 1.0)
    g["tx_panel_qv20"] = int(len(t)); g["ok"] = g["one_to_one"]
    return mt.s.to_numpy().astype(np.int64), Tt, g

def cosine(A, T):
    num = (A * T).sum(axis=1); den = np.sqrt((A * A).sum(axis=1) * (T * T).sum(axis=1))
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)

def run_unit(arm, unit, out_root=OUT):
    u = [v for v in units() if v["arm"] == arm and v["unit"] == unit][0]; od = f"{out_root}/{arm}/{unit}"; os.makedirs(od, exist_ok=True)
    if os.path.exists(f"{od}/unit.json"): return "present"
    t0 = time.time(); bins, genes, X = load(u["d"]); tb = np.flatnonzero(bins.in_tissue.to_numpy()); gp_ids = S5.gp()
    gix = pd.Series(np.arange(len(genes)), index=genes); gpi = gix.reindex(gp_ids).to_numpy(); assert not np.isnan(gpi.astype(float)).any(); gpi = gpi.astype(np.int64)
    seed = SEED + 100 * u["i"] + u["s"] + (0 if arm == "B" else 50)
    a, cells, res, n_hvg = run_enact(bins, genes, X, tb, gp_ids, f"{od}/enact", seed); gate_e = enact_gate(a, cells, res)
    if arm == "B": labels, Tt, gate_t = truth_B(u, bins, X, tb, gpi)
    else: labels, Tt, gate_t = truth_X1(u, bins, gp_ids)
    # cellule della catena: nucleo con almeno meta' dei quadratini in in_reg
    lab = bins.label.to_numpy(); reg = bins.in_reg.to_numpy(); L, inv = np.unique(lab[lab > 0], return_inverse=True)
    nb = np.bincount(inv, minlength=len(L)); nr = np.bincount(inv, weights=reg[lab > 0].astype(float), minlength=len(L)); inreg = set(L[nr >= 0.5 * nb])
    A = {}
    for m in EM: A[m] = res[m].reindex(index=labels, columns=list(gp_ids)).fillna(0.0).to_numpy()
    exp = np.load(f"{od}/exp_b2c.npy"); A["b2c"] = group_sum(exp, X[:, gpi], labels)
    keep = np.array([l in inreg for l in labels]) & (Tt.sum(axis=1) > 0)
    for m in METHODS: keep &= A[m].sum(axis=1) > 0
    rng = np.random.default_rng([SEED, u["i"]]); perm = rng.permutation(int(keep.sum()))
    df = pd.DataFrame({"label": labels[keep], "truth_panel": Tt[keep].sum(axis=1)})
    for m in METHODS:
        df[f"panel_{m}"] = A[m][keep].sum(axis=1); df[f"cos_{m}"] = cosine(A[m][keep], Tt[keep]); df[f"cosshuf_{m}"] = cosine(A[m][keep], Tt[keep][perm])
    df.to_parquet(f"{od}/cells.parquet", index=False)
    info = dict(u, n_hvg=n_hvg, n_genes_enact=int(a.shape[1]), n_cells_truth=int(len(labels)), n_cells_inreg=int(sum(l in inreg for l in labels)), n_kept=int(keep.sum()),
                score={m: float(df[f"cos_{m}"].mean()) for m in METHODS}, score_shuf={m: float(df[f"cosshuf_{m}"].mean()) for m in METHODS},
                g_enact=gate_e, g_truth=gate_t, seed_np=seed, sec=round(time.time() - t0, 1))
    json.dump(info, open(f"{od}/unit.json", "w"), indent=1, default=str); return f"{arm} {unit} kept {keep.sum()} {info['sec']} s"

def taub(a, b):
    t = kendalltau(a, b).statistic; return float(t) if np.isfinite(t) else 0.0

def summary():
    os.makedirs(RES, exist_ok=True); rows = []
    for u in units():
        j = json.load(open(f"{OUT}/{u['arm']}/{u['unit']}/unit.json"))
        for m in METHODS: rows.append(dict(arm=u["arm"], unit=u["unit"], i=u["i"], s=u["s"], tissue=u["tissue"], method=m, score=j["score"][m], score_shuf=j["score_shuf"][m], n_kept=j["n_kept"]))
    pa_ = pd.DataFrame(rows); pa_.to_csv(f"{RES}/per_unit.csv", index=False)
    S = {arm: pa_[pa_.arm == arm].groupby(["i", "method"]).score.mean().unstack()[METHODS].sort_index().to_numpy() for arm in ("B", "X1")}
    Sh = pa_[pa_.arm == "X1"].groupby(["i", "method"]).score_shuf.mean().unstack()[METHODS].sort_index().to_numpy()
    import itertools
    TU, ST = list(range(5)), list(range(5, 10)); splits = [(list(a) + list(b), [x for x in TU if x not in a] + [x for x in ST if x not in b]) for a in itertools.combinations(TU, 3) for b in itertools.combinations(ST, 2)]
    def stats(h1, h2):
        m = lambda Z, h: Z[h].mean(axis=0)
        t5 = taub(m(S["X1"], h1), m(S["X1"], h2)); s5 = 0.5 * (taub(m(S["B"], h1), m(S["X1"], h2)) + taub(m(S["B"], h2), m(S["X1"], h1)))
        return t5, s5, s5 - t5, taub(m(Sh, h1), m(Sh, h2))
    P = np.array([stats(h1, h2) for h1, h2 in splits]); est = P.mean(axis=0)
    rng = np.random.default_rng(SEED); BS = []
    for _ in range(2000):
        h1, h2 = splits[rng.integers(len(splits))]
        r = lambda h: [int(v) for t in (TU, ST) for v in rng.choice([x for x in h if x in t], size=len([x for x in h if x in t]), replace=True)]
        BS.append(stats(r(h1), r(h2)))
    BS = np.array(BS); se = BS.std(axis=0, ddof=1)
    names = ["T5", "S5c", "A5", "T5_shuf"]; out = {n: float(est[k]) for k, n in enumerate(names)}; out.update({f"{n}_se": float(se[k]) for k, n in enumerate(names)})
    out["tau_full_B_X1"] = taub(S["B"].mean(axis=0), S["X1"].mean(axis=0)); out["methods"] = METHODS
    out["score_B"] = dict(zip(METHODS, S["B"].mean(axis=0).tolist())); out["score_X1"] = dict(zip(METHODS, S["X1"].mean(axis=0).tolist())); out["score_X1_shuf"] = dict(zip(METHODS, Sh.mean(axis=0).tolist()))
    T5, T5s, A5, A5s = out["T5"], out["T5_se"], out["A5"], out["A5_se"]
    if T5 - 2 * T5s <= 0 <= T5 + 2 * T5s: out["esito"] = "non misurabile a questa profondita"
    elif A5 - 2 * A5s >= -0.2: out["esito"] = "regge"
    elif A5 + 2 * A5s < -0.2: out["esito"] = "non regge"
    else: out["esito"] = "incerta"
    out["cp_shuf_contains_0"] = bool(out["T5_shuf"] - 2 * out["T5_shuf_se"] <= 0 <= out["T5_shuf"] + 2 * out["T5_shuf_se"])
    out["n_splits"] = len(splits); out["B_boot"] = 2000; out["seed"] = SEED
    json.dump(out, open(f"{RES}/A5.json", "w"), indent=1); pd.DataFrame(BS, columns=names).to_csv(f"{RES}/boot.csv", index=False)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, dict)}, indent=0))

if __name__ == "__main__":
    if sys.argv[1] == "unit": print(run_unit(sys.argv[2], sys.argv[3], *(sys.argv[4:5] or [OUT])))
    elif sys.argv[1] == "summary": summary()
    elif sys.argv[1] == "list":
        for u in units(): print(u["arm"], u["unit"], u["d"], sep="\t")
