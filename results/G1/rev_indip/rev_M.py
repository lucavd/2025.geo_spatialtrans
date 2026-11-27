"""Independent bottom-up recomputation of M-1..M-6 per area from chain outputs (bins, counts, bin2cell expansion, RCTD blocks).
Transferred labels recomputed from the HLCA reference (multinomial ML, pseudocount 0.5, >= 20 UMI). Own code."""
import os, sys, glob, json, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
from multiprocessing import Pool
G = "/mnt/micron/geo_spatialtrans/"; OUT = G + "G1_review/rev_indip/"; REF = G + "S3_sonda/ref/"
def runs():
    L = []
    for t in ("tumour", "stroma"):
        for k in range(1, 6):
            L.append(("V", t, k, "-", f"{G}S3_sonda/runs/V/H2_{t}_{k}"))
            L.append(("T", t, k, "-", f"{G}S3b_sonda/runs/T/H1_{t}_{k}"))
            for s in (0, 1): L.append(("B", t, k, str(s), f"{G}S3_sonda/runs/B/H2_{t}_{k}_k{s}"))
            L.append(("Y", t, k, "0", f"{G}S4_sonda/runs/Y/H2_{t}_{k}_k0"))
            for a in ("Yalt", "YH1", "Yshuf"): L.append((a, t, k, "0", f"{G}S6_sonda/runs/{a}/H2_{t}_{k}_k0"))
    return L
FILES = ["bins.parquet", "X.parquet", "genes.tsv", "exp_labels_expanded.npy", "cells.parquet", "rctd_blocks.csv", "labels_transferred.parquet"]
def prep():
    rc = pd.read_csv(REF + "ref_L_cells.csv"); rg = pd.read_csv(REF + "ref_L_genes.tsv", header=None)[0].astype(str).to_numpy()
    M = sio.mmread(REF + "ref_L_counts.mtx").tocsc()          # genes x cells
    assert M.shape == (len(rg), len(rc)), M.shape
    types = sorted(pd.read_csv(REF + "markers_L.csv").cell_type.unique())
    S = np.vstack([np.asarray(M[:, np.flatnonzero(rc.cell_type.to_numpy() == ty)].sum(axis=1)).ravel() for ty in types])
    np.savez(OUT + "ref_type_sums.npz", S=S, genes=rg, types=np.array(types))
def one(r):
    arm, t, k, seed, d = r
    ref = np.load(OUT + "ref_type_sums.npz", allow_pickle=True); S, rg, types = ref["S"], ref["genes"].astype(str), list(ref["types"].astype(str))
    bins = pd.read_parquet(d + "/bins.parquet"); genes = pd.read_csv(d + "/genes.tsv", header=None)[0].astype(str).to_numpy()
    T = pd.read_parquet(d + "/X.parquet"); nb, ng = len(bins), len(genes)
    X = sp.coo_matrix((T.n.to_numpy().astype(float), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(nb, ng)).tocsr()
    exp = np.load(d + "/exp_labels_expanded.npy"); assert len(exp) == nb
    cells = pd.read_parquet(d + "/cells.parquet"); umi = np.asarray(X.sum(1)).ravel()
    gpos = {g: i for i, g in enumerate(genes)}
    # --- transferred labels (own implementation)
    keep_ref = np.array([g in gpos for g in rg]); col = np.array([gpos[g] for g in rg[keep_ref]])
    P = S[:, keep_ref] + 0.5; LP = np.log(P / P.sum(1, keepdims=True))
    def agg(labels):  # counts summed over the expanded territory of each label -> csr (len(labels) x genes)
        lab = np.asarray(labels, np.int64); m = pd.Series(np.arange(len(lab)), index=lab)
        rows = m.reindex(exp).to_numpy(); ok = ~np.isnan(rows)
        A = sp.csr_matrix((np.ones(ok.sum()), (rows[ok].astype(int), np.flatnonzero(ok))), shape=(len(lab), nb))
        return (A @ X).tocsr()
    labs = cells.label.to_numpy(np.int64); Cl = agg(labs); ul = np.asarray(Cl.sum(1)).ravel()
    ll = np.asarray(Cl[:, col] @ LP.T); my = np.array(types, dtype=object)[ll.argmax(1)]; my[ul < 20] = None
    srt = np.sort(ll, 1); gap = srt[:, -1] - srt[:, -2]
    tl = pd.read_parquet(d + "/labels_transferred.parquet").set_index("label").reindex(labs)
    same = (pd.Series(my).fillna("NA").to_numpy() == tl.cell_type.fillna("NA").to_numpy())
    ctype = pd.Series(my, index=labs)
    # --- measures
    reg = bins.in_reg.to_numpy(bool); nuc = bins.label.to_numpy(np.int64)
    out = {"arm": arm, "A": t, "area": k, "seed": seed, "dir": d, "n_labels": int(len(labs)), "label_agree": int(same.sum()),
           "min_ll_gap_disagree": float(gap[~same].min()) if (~same).any() else None}
    out["M1"] = umi[reg & (exp > 0)].sum() / umi[reg].sum()
    nn = pd.DataFrame({"l": nuc[nuc > 0], "r": reg[nuc > 0]}).groupby("l").r.agg(["size", "sum"])
    cs = nn.index.to_numpy()[nn["sum"].to_numpy() >= 0.5 * nn["size"].to_numpy()]; out["n_cells"] = int(len(cs))
    C = agg(cs); u = np.asarray(C.sum(1)).ravel(); out["M2"] = float(u.mean())
    mk = pd.read_csv(REF + "markers_L.csv"); mk = mk[mk.gene.isin(list(gpos))]
    MU = np.column_stack([np.asarray(C[:, [gpos[g] for g in mk.gene[mk.cell_type == ty]]].sum(1)).ravel() for ty in types])
    tot = MU.sum(1); ct = ctype.reindex(cs).to_numpy()
    own = np.array([MU[i, types.index(c)] if c in types else np.nan for i, c in enumerate(ct)])
    typed = np.array([c in types for c in ct]) & (tot >= 5); out["n_typed"] = int(typed.sum())
    out["M3"] = float(np.mean(1 - own[typed] / tot[typed]))
    b = pd.read_csv(d + "/rctd_blocks.csv"); ok = b.spot_class.notna() & (b.spot_class != "reject")
    out["M4"] = float((b.spot_class[ok] == "singlet").mean())
    nt = ctype.reindex(b.nuc_label.to_numpy()).to_numpy(); nt = np.where(b.nuc_bins.to_numpy() >= 4, nt, None)
    hn = ok.to_numpy() & np.array([x in types for x in nt]); out["n_nuc_blocks"] = int(hn.sum())
    out["M5"] = float((b["first"].to_numpy()[hn] == nt[hn]).mean())
    bo = b[ok]; first = dict(zip(zip(bo.r8, bo.c8), bo["first"])); eq = n = 0
    for (r8, c8), f in first.items():
        for nbk in ((r8, c8 + 1), (r8 + 1, c8)):
            if nbk in first: n += 1; eq += (first[nbk] == f)
    out["M6"] = eq / n; out["n_pairs"] = n
    return out
if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    if sys.argv[1] == "md5":
        import hashlib
        fl = [REF + f for f in ("ref_L_cells.csv", "ref_L_genes.tsv", "ref_L_counts.mtx", "markers_L.csv")] + [r[4] + "/" + f for r in runs() for f in FILES]
        with open(OUT + sys.argv[2], "w") as fh:
            for f in fl: fh.write(hashlib.md5(open(f, "rb").read()).hexdigest() + "  " + f + "\n")
        print("md5", len(fl))
    elif sys.argv[1] == "prep": prep(); print("prep done")
    else:
        with Pool(4) as p: res = p.map(one, runs(), chunksize=1)
        pd.DataFrame(res).to_csv(OUT + "rev_M_cells.csv", index=False); print("done", len(res))
