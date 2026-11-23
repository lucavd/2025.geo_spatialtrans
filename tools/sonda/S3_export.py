"""tools/sonda/S3_export.py - S3: aree vere nel formato di F3.1 (bins.parquet, X.parquet, genes.tsv, cells.parquet) e braccio P.
Modi:
  crop   <sec> [smoke]  ritaglio dell'H&E (area + 50 pixel) per StarDist: <ROOT>/seg/<unit>/roi.tif, crop.json
  export <sec> [smoke]  quadratini da 2 um con centro nel quadrato dell'area; conteggi dalla matrice grezza (sole colonne dell'area; quadratini assenti = zero);
                        geni = identificativi Ensembl delle feature della matrice filtrata (pannello, come F3.1); label = nucleo StarDist sotto il pixel del centro del quadratino;
                        in_reg = in_tissue e centro a >= 25 um dal bordo; uscita <ROOT>/runs/<V|T>/<unit>/
  perm   H2             braccio P (controllo C-4): vettori di conteggio dei quadratini in tissue permutati fra i quadratini in tissue,
                        default_rng([20261126, i]), i = 0..9 (tumore 1..5, stroma 1..5); nuclei fermi. <ROOT>/runs/P/<unit>/
smoke: un'unica area ridotta, il quadrato di 250 um con angolo (x0 + 250, y0 + 250) dell'area H2 tumour 1 (check smoke, C-3), sotto
<root>/runs/SMOKE_V/. La radice si cambia con S3_ROOT. Uso: .venv/bin/python -B tools/sonda/S3_export.py <modo> <sec> [smoke]"""
import sys, os, json, time, shutil, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow as pa, pyarrow.parquet as pq, tifffile
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
ROOT = os.environ.get("S3_ROOT", C.ROOT); MARGIN = 50
mode, sec = sys.argv[1], sys.argv[2]; smoke = len(sys.argv) > 3 and sys.argv[3] == "smoke"
ARM = {"H2": "V", "H1": "T"}

def units():
    a = C.areas(sec)
    if smoke:
        r = a[(a.tissue == "tumour") & (a.k == 1)].iloc[0]
        return [dict(unit=f"SMOKE_{sec}_tumour_1", tissue="tumour", k=1, i=0, x0=r.x0 + 250.0, y0=r.y0 + 250.0, side=250.0, arm="SMOKE_" + ARM[sec])]
    out = []
    for i, (_, r) in enumerate(a.sort_values(["tissue", "k"], key=lambda s: s.map({"tumour": 0, "stroma": 1}) if s.name == "tissue" else s).iterrows()):
        out.append(dict(unit=C.unit_name(sec, r.tissue, int(r.k)), tissue=r.tissue, k=int(r.k), i=i, x0=float(r.x0), y0=float(r.y0), side=float(r.side), arm=ARM[sec]))
    return out

def area_bins(tp, u, m):
    x = (tp.array_col.to_numpy() + 0.5) * 2.0; y = (tp.array_row.to_numpy() + 0.5) * 2.0
    s = (x >= u["x0"]) & (x < u["x0"] + u["side"]) & (y >= u["y0"]) & (y < u["y0"] + u["side"])
    b = tp[s].copy(); xa = (b.array_col.to_numpy() + 0.5) * 2.0; xm = C.guard(xa, m) if C.guard_on(sec) else float(xa.max())
    b["x_um"] = (b.array_col.to_numpy() + 0.5) * 2.0 - u["x0"]; b["y_um"] = (b.array_row.to_numpy() + 0.5) * 2.0 - u["y0"]
    return b.sort_values(["array_row", "array_col"]).reset_index(drop=True), xm

if mode in ("crop", "export"):
    m = C.split_m(sec); tp = C.positions(sec, "002"); upp = C.upp(sec)
    if mode == "export":
        h5 = f"{C.bdir(sec, '002')}/raw_feature_bc_matrix.h5"; ids_raw, _ = C.features(h5)
        ids, _ = C.features(f"{C.bdir(sec, '002')}/filtered_feature_bc_matrix.h5")      # geni del pannello: feature della filtrata, come F3.1
        gsel = pd.Series(np.arange(len(ids_raw)), index=ids_raw).reindex(ids).to_numpy(); assert not np.isnan(gsel).any(); gsel = gsel.astype(np.int64)
        bpos = pd.Series(np.arange(len(C.barcodes(h5))), index=C.barcodes(h5))
    for u in units():
        t0 = time.time(); sd = f"{ROOT}/seg/{u['unit']}"; os.makedirs(sd, exist_ok=True)
        b, xm = area_bins(tp, u, m)
        if mode == "crop":
            if os.path.exists(f"{sd}/crop.json"): print("present", u["unit"]); continue
            pr, pc = b.pxl_row_in_fullres.to_numpy(), b.pxl_col_in_fullres.to_numpy()
            r0, r1 = int(np.floor(pr.min())) - MARGIN, int(np.ceil(pr.max())) + MARGIN + 1; c0, c1 = int(np.floor(pc.min())) - MARGIN, int(np.ceil(pc.max())) + MARGIN + 1
            rgb = C.read_region(C.image(sec), max(r0, 0), r1, max(c0, 0), c1)
            tifffile.imwrite(f"{sd}/roi.tif", rgb, photometric="rgb", compression="zlib")
            json.dump(dict(u, section=sec, r0=max(r0, 0), r1=r1, c0=max(c0, 0), c1=c1, shape=list(rgb.shape), upp=upp, m_um=m, xmax_um=xm), open(f"{sd}/crop.json", "w"))
            print("crop", u["unit"], rgb.shape, round(time.time() - t0, 1), flush=True); continue
        d = f"{ROOT}/runs/{u['arm']}/{u['unit']}"
        if os.path.exists(f"{d}/export.json"): print("present", u["unit"]); continue
        os.makedirs(d, exist_ok=True); cj = json.load(open(f"{sd}/crop.json")); lab = np.load(f"{sd}/stardist.npz")["labels"]
        rr = np.rint(b.pxl_row_in_fullres.to_numpy()).astype(np.int64) - cj["r0"]; cc = np.rint(b.pxl_col_in_fullres.to_numpy()).astype(np.int64) - cj["c0"]
        assert rr.min() >= 0 and cc.min() >= 0 and rr.max() < lab.shape[0] and cc.max() < lab.shape[1]
        label = lab[rr, cc].astype(np.int64)
        side = u["side"]; x, y = b.x_um.to_numpy(), b.y_um.to_numpy(); it = b.in_tissue.to_numpy() == 1
        in_reg = it & (x >= C.FRAME_UM) & (x <= side - C.FRAME_UM) & (y >= C.FRAME_UM) & (y <= side - C.FRAME_UM)
        bins = pd.DataFrame({"bin_id": np.arange(len(b), dtype=np.int64), "row": b.array_row.to_numpy().astype(np.int64), "col": b.array_col.to_numpy().astype(np.int64),
                             "x_um": x, "y_um": y, "in_tissue": it, "label": label, "in_reg": in_reg})
        col = bpos.reindex(b.barcode.to_numpy()).to_numpy(); have = ~np.isnan(col)   # quadratini assenti dalla grezza: zero conteggi
        Xp, nr = C.read_cols(h5, col[have].astype(np.int64)); Xp = Xp[:, gsel]
        Sel = sp.csr_matrix((np.ones(have.sum()), (np.flatnonzero(have), np.arange(have.sum()))), shape=(len(b), int(have.sum())))
        X = (Sel @ Xp).tocsr().astype(np.int64); T = X.tocoo()
        pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": T.row.astype(np.int64), "gene": T.col.astype(np.int64), "n": T.data.astype(np.int64)}), preserve_index=False), f"{d}/X.parquet")
        bins.to_parquet(f"{d}/bins.parquet", index=False)
        with open(f"{d}/genes.tsv", "w") as f: f.write("\n".join(ids) + "\n")
        labs = np.unique(label[label > 0]); pd.DataFrame({"label": labs, "cell_type": pd.Series([None] * len(labs), dtype=object)}).to_parquet(f"{d}/cells.parquet", index=False)
        info = dict(u, section=sec, side_um=side, n_bins=int(len(bins)), n_tissue=int(it.sum()), n_reg=int(in_reg.sum()), n_nuclei=int(len(labs)),
                    nuclei_stardist_crop=int(lab.max()), umi_total=int(T.data.sum()), umi_reg=int(np.asarray(X[in_reg].sum())), m_um=m, xmax_um=xm,
                    n_ranges=int(nr), n_genes=int(len(ids)), n_missing_barcodes=int((~have).sum()), n_missing_in_tissue=int((~have & it).sum()),
                    sec_export=round(time.time() - t0, 1))
        json.dump(info, open(f"{d}/export.json", "w")); print(json.dumps(info), flush=True)

if mode == "perm":
    for u in units():
        s = f"{ROOT}/runs/V/{u['unit']}"; d = f"{ROOT}/runs/P/{u['unit']}"
        if os.path.exists(f"{d}/export.json"): print("present", u["unit"]); continue
        os.makedirs(d, exist_ok=True)
        for f in ("bins.parquet", "genes.tsv", "cells.parquet"): shutil.copy(f"{s}/{f}", f"{d}/{f}")
        bins = pd.read_parquet(f"{s}/bins.parquet"); T = pd.read_parquet(f"{s}/X.parquet"); ng = len(open(f"{s}/genes.tsv").read().split())
        X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), ng))
        it = np.flatnonzero(bins.in_tissue.to_numpy()); rng = np.random.default_rng([C.SEED_P, u["i"]]); src = np.arange(len(bins)); src[it] = it[rng.permutation(len(it))]
        Y = X[src].tocoo()
        pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": Y.row.astype(np.int64), "gene": Y.col.astype(np.int64), "n": Y.data.astype(np.int64)}), preserve_index=False), f"{d}/X.parquet")
        info = dict(json.load(open(f"{s}/export.json")), arm="P", perm_seed=[C.SEED_P, u["i"]], umi_total_check=int(Y.data.sum()))
        assert info["umi_total_check"] == info["umi_total"]; json.dump(info, open(f"{d}/export.json", "w")); print("perm", u["unit"], flush=True)
