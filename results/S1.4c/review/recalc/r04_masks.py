#!/usr/bin/env python
"""results/S1.4c/review/recalc/r04_masks.py — revisione avversariale S1.4c (recalc).
S-1 reale, n_bnd/n_c2 per nucleo (S-2) e coppie a contatto (CP-0) direttamente dalle maschere native, con un metodo
diverso dall'operatore del banco: ciclo per etichetta su bounding box (ndimage.find_objects) e dilatazione binaria
5x5 (Chebyshev 2) invece dei filtri max/min sull'immagine intera. Appaiamento A4/A6 (BL-060) rifatto con shapely
(centroide StarDist dentro un solo territorio Space Ranger; territorio con un solo centroide). CP-0: Delaunay di scipy
sui centroidi ricalcolati dalla maschera; AUC Mann-Whitney in double; score come implementato ((r_i+r_j)/d) e come
scritto nella pre-registrazione (r_i+r_j-d).
Uso: .venv/bin/python r04_masks.py [ncores]
"""
import sys, os, json, numpy as np, pandas as pd, shapely
from scipy import ndimage as ndi
from scipy.spatial import Delaunay
from scipy.stats import rankdata
from multiprocessing import Pool
B = "/mnt/micron/geo_spatialtrans/S1.4c/"; D = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump/"
OUTD = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/masks/"; os.makedirs(OUTD, exist_ok=True)
R3 = "/mnt/micron/geo_spatialtrans/R3/real/"; MASKS = "/mnt/micron/geo_spatialtrans/R2/masks/"
SEG = dict(A1="cellpose_rgb", A2="cellpose_rgb", A3="cellpose_rgb", A4="stardist_he", A5="cellpose_rgb", A6="stardist_he")
K5 = np.ones((5, 5), bool); CROSS = ndi.generate_binary_structure(2, 1)

def per_label(lab, d=2):
    """Per ogni etichetta > 0: n_px, somme r/c, n_bnd (4-vicino diverso o fuori immagine), n_cd (pixel di bordo con
    un'altra etichetta entro Chebyshev d), insieme dei partner a contatto. Ciclo esplicito su bounding box."""
    H, W = lab.shape; K = np.ones((2 * d + 1, 2 * d + 1), bool); sl = ndi.find_objects(lab)
    rows = []; pairs = []
    for i, s in enumerate(sl):
        if s is None: continue
        L = i + 1
        r0 = max(s[0].start - d - 1, 0); r1 = min(s[0].stop + d + 1, H); c0 = max(s[1].start - d - 1, 0); c1 = min(s[1].stop + d + 1, W)
        cr = lab[r0:r1, c0:c1]; me = cr == L
        # bordo: pixel di L con un 4-vicino non L; i lati dell'immagine contano come diverso
        pad = np.pad(me, 1, constant_values=False)
        if r0 > 0: pad[0, :] = np.pad(lab[r0 - 1, c0:c1] == L, 1, constant_values=False)
        if r1 < H: pad[-1, :] = np.pad(lab[r1, c0:c1] == L, 1, constant_values=False)
        if c0 > 0: pad[1:-1, 0] = lab[r0:r1, c0 - 1] == L
        if c1 < W: pad[1:-1, -1] = lab[r0:r1, c1] == L
        inner = pad[1:-1, 1:-1] & pad[:-2, 1:-1] & pad[2:, 1:-1] & pad[1:-1, :-2] & pad[1:-1, 2:]
        bnd = me & ~inner
        oth = (cr > 0) & ~me
        near = ndi.binary_dilation(oth, structure=K) if oth.any() else np.zeros_like(me)
        dil = ndi.binary_dilation(me, structure=K)
        part = np.unique(cr[dil & oth])
        rr, cc = np.nonzero(me)
        rows.append((L, me.sum(), (rr + r0).sum(), (cc + c0).sum(), bnd.sum(), (bnd & near).sum(), len(part)))
        pairs += [(min(L, int(p)), max(L, int(p))) for p in part]
    df = pd.DataFrame(rows, columns=["label", "n_px", "sr", "sc", "n_bnd", "n_c", "n_part"])
    pr = pd.DataFrame(sorted(set(pairs)), columns=["label_i", "label_j"])
    return df, pr

def auc(s, y):
    y = np.asarray(y, bool); r = rankdata(s); n1 = float(y.sum()); n0 = float((~y).sum())
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def pair_A4A6(A, roi, cc):
    """BL-060 rifatto: nucleo StarDist il cui centroide cade in esattamente un territorio SR; territorio appaiato se
    contiene esattamente un centroide StarDist."""
    sg = pd.read_parquet(R3 + f"{A}_{roi}_stardist_he_gen.parquet")
    T = shapely.from_wkb(pd.read_parquet(D + f"terr_{A}_{roi}.parquet").wkb.to_numpy()); n = len(T)
    tree = shapely.STRtree(T); pts = shapely.points(sg.x.to_numpy(), sg.y.to_numpy())
    pi, ti = tree.query(pts, predicate="intersects")
    cnt_p = np.bincount(pi, minlength=len(sg)); one = cnt_p[pi] == 1
    h1 = np.full(len(sg), -1); h1[pi[one]] = ti[one]
    cnt_t = np.bincount(h1[h1 >= 0], minlength=n); matched = cnt_t == 1
    lab = np.full(n, -1, np.int64)
    for p_ in np.nonzero(h1 >= 0)[0]:
        t = h1[p_]
        if matched[t] and lab[t] < 0: lab[t] = int(sg.label.iloc[p_])
    return matched, lab

def one(args):
    A, roi = args; seg = SEG[A]
    o = OUTD + f"{A}_{roi}.json"
    if os.path.exists(o): return json.load(open(o))
    gen = pd.read_parquet(R3 + f"{A}_{roi}_{seg}_gen.parquet"); upp = float(gen.um_per_px.iloc[0])
    lab = np.load(MASKS + f"{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    keep = np.zeros(int(lab.max()) + 1, bool); keep[gen.label.to_numpy().astype(np.int64)] = True; lab[~keep[lab]] = 0
    df, pr = per_label(lab, 2); df = df.set_index("label")
    df["cx"] = df.sc / df.n_px * upp; df["cy"] = df.sr / df.n_px * upp; df["area"] = df.n_px * upp ** 2
    cc = pd.read_parquet(D + f"cells_{A}_{roi}.parquet"); ev = cc.ev_prep.to_numpy().astype(bool)
    res = dict(archetype=A, roi=roi, seg=seg, n_lab=len(df), n_gen=len(gen))
    if A in ("A4", "A6"):
        matched, sl = pair_A4A6(A, roi, cc)
        res["pair_matched_identical"] = bool((matched == cc.matched.to_numpy().astype(bool)).all())
        pl = cc.shape_label.to_numpy()
        res["pair_label_identical"] = bool((sl[matched] == pl[matched]).all())
        ev_lab = sl[ev]
    else:
        ev_lab = gen.label.to_numpy()[ev]
    # confronto con l'operatore del banco (real_shape) e con R3
    rs = pd.read_parquet(B + f"real_shape/{A}_{roi}_{seg}_shape.parquet").set_index("label").loc[df.index]
    res["max_dcx"] = float(np.abs(df.cx - rs.cx).max()); res["max_darea"] = float(np.abs(df.area - rs.area).max())
    res["n_bnd_diff"] = int((df.n_bnd != rs.n_bnd).sum()); res["n_c2_diff"] = int((df.n_c != rs.n_c2).sum())
    e = df.loc[ev_lab]
    res["n_ev"] = int(len(e)); res["s1_mine"] = float((e.n_part > 0).mean()); res["s1_from_nc2_mine"] = float((e.n_c > 0).mean())
    res["s2_median_mine"] = float(np.median(e.n_c / e.n_bnd))
    # coppie
    rp = pd.read_parquet(B + f"real_shape/{A}_{roi}_{seg}_pairs.parquet")
    a = set(map(tuple, pr.to_numpy())); b = set(map(tuple, rp.to_numpy()))
    res["n_pairs_mine"] = len(a); res["n_pairs_rep"] = len(b); res["pairs_symdiff"] = len(a ^ b)
    # CP-0: Delaunay scipy su tutti i centroidi della maschera, lati con almeno un estremo ev
    labs = df.index.to_numpy(); P = df[["cx", "cy"]].to_numpy()
    tri = Delaunay(P).simplices
    E = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]]), axis=1), axis=0)
    evs = np.isin(labs, ev_lab); E = E[evs[E[:, 0]] | evs[E[:, 1]]]
    dd = np.hypot(*(P[E[:, 0]] - P[E[:, 1]]).T); rr = np.sqrt(df.area.to_numpy() / np.pi)
    lo = np.minimum(labs[E[:, 0]], labs[E[:, 1]]); hi = np.maximum(labs[E[:, 0]], labs[E[:, 1]])
    M = int(labs.max()) + 1; key = set((pr.label_i.to_numpy() * M + pr.label_j.to_numpy()).tolist())
    y = np.fromiter(((l * M + h) in key for l, h in zip(lo, hi)), bool, len(lo))
    s_ratio = (rr[E[:, 0]] + rr[E[:, 1]]) / dd; s_diff = rr[E[:, 0]] + rr[E[:, 1]] - dd
    res.update(n_edges=int(len(E)), frac_contact=float(y.mean()), auc_ratio=float(auc(s_ratio, y)), auc_diff=float(auc(s_diff, y)))
    np.savez_compressed(OUTD + f"{A}_{roi}_cp0.npz", s_ratio=s_ratio, s_diff=s_diff, y=y)
    json.dump(res, open(o, "w")); return res

if __name__ == "__main__":
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    with Pool(nc) as p:
        res = pd.DataFrame(p.map(one, list(zip(rois.archetype, rois.roi_id)), chunksize=1))
    refs = pd.read_csv("results/S1.4c/S1.4c_real_refs.csv")[["archetype", "roi", "n_ev", "s1", "s2_median"]]
    c0 = pd.read_csv("results/S1.4c/S1.4c_cp0_roi.csv").rename(columns={"n_pairs": "cp0_n_rep", "frac_contact": "cp0_fc_rep", "auc": "cp0_auc_rep"})
    res = res.merge(refs, on=["archetype", "roi"], suffixes=("", "_rep")).merge(c0, on=["archetype", "roi"])
    pooled = []
    for A in SEG:
        z = [np.load(OUTD + f"{A}_{r}_cp0.npz") for r in rois.roi_id[rois.archetype == A]]
        y = np.concatenate([q["y"] for q in z])
        pooled.append(dict(archetype=A, n_edges=len(y), frac_contact=float(y.mean()), auc_ratio=auc(np.concatenate([q["s_ratio"] for q in z]), y),
                           auc_diff=auc(np.concatenate([q["s_diff"] for q in z]), y)))
    pooled = pd.DataFrame(pooled).merge(pd.read_csv("results/S1.4c/S1.4c_cp0.csv"), on="archetype", suffixes=("", "_rep"))
    res.to_csv(OUTD + "../r04_masks_roi.csv", index=False); pooled.to_csv(OUTD + "../r04_cp0_pooled.csv", index=False)
    pd.set_option("display.width", 250)
    print(res[["archetype", "roi", "n_ev", "n_ev_rep", "s1_mine", "s1", "s1_from_nc2_mine", "s2_median_mine", "s2_median", "n_bnd_diff", "n_c2_diff", "pairs_symdiff", "max_dcx", "max_darea"]].to_string())
    print(res[[c for c in res.columns if c.startswith("pair_")] + ["archetype", "roi"]].dropna().to_string())
    print(res[["archetype", "roi", "n_edges", "cp0_n_rep", "frac_contact", "cp0_fc_rep", "auc_ratio", "cp0_auc_rep", "auc_diff"]].to_string())
    print(pooled.to_string())
    print("max|s1 mio - riportato|", float(np.abs(res.s1_mine - res.s1).max()), "max|auc_ratio - rep| per ROI", float(np.abs(res.auc_ratio - res.cp0_auc_rep).max()))
