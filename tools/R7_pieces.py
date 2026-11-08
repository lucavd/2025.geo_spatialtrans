"""tools/R7_pieces.py — R7: geometria (M-1), pezzi disgiunti (M-2), cancelli G-2 e G-3, valori per area e per pezzo
(P-3, M-e), quota di nuclei condivisi nelle validazioni lasciando fuori un'area (M-d), due stadi contro uno (M-f).
Uso: .venv/bin/python tools/R7_pieces.py  -> results/R7/pieces/*.csv"""
import os, sys, json, itertools
import numpy as np, pandas as pd, pyarrow.parquet as pq, h5py
from PIL import Image
sys.path.insert(0, "tools")
import R7_common as C
import R4_common as R4
OUT = "results/R7/pieces"; os.makedirs(OUT, exist_ok=True)
R4RUN = "/mnt/micron/geo_spatialtrans/R4/run1"
DR = "data/real/datasets"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; AIDX = {a: i for i, a in enumerate(ARCH)}
E0 = {("A1", "r3", "r5"): .345, ("A1", "r2", "r4"): .122, ("A1", "r3", "r4"): .028, ("A2", "r4", "r5"): .077,
      ("A3", "r1", "r2"): .278, ("A3", "r2", "r5"): .278, ("A3", "r3", "r5"): .240, ("A3", "r3", "r4"): .120,
      ("A4", "f2", "p2"): .384, ("A4", "p2", "p3"): .278, ("A4", "f2", "p1"): .131, ("A4", "f2", "p3"): .082, ("A4", "f1", "f2"): .002,
      ("A5", "r2", "r5"): .301}
R4COLS = ["mean_umi2", "median_umi2", "q90_umi2", "q99_umi2", "cv_umi2", "frac_empty2", "mean_umi8", "median_umi8", "q90_umi8",
          "q99_umi8", "cv_umi8", "frac_empty8", "mean_genes2", "mean_genes8", "sparsity2", "sparsity8",
          "share_nucleus", "share_ring", "share_acellular", "share_other_nucleus", "mean_umi2_nucleus", "mean_umi2_ring",
          "mean_umi2_acellular", "frac_empty2_nucleus", "frac_empty2_ring", "frac_empty2_acellular"]
P3_R4 = ["mean_umi2", "mean_umi8", "frac_empty2", "mean_genes8", "rho_obs"]

rois = pd.read_csv("results/R2/R2_rois_checked.csv")
R4roi = pd.read_csv("results/R4/R4_roi.csv"); R4refs = pd.read_csv("results/R4/R4_refs.csv")
nuc_all = pq.read_table("results/R2/R2_nuclei_all.parquet").to_pandas()
nuc_all = nuc_all[(nuc_all.method == "cellpose_rgb") & (nuc_all.scale == 1) & nuc_all.keep]
summ = pd.read_csv("results/R2/R2_roi_summary.csv"); summ = summ[(summ.method == "cellpose_rgb") & (summ.scale == 1)]

def vals_2um(u, s, ng):
    return {"mean_umi2": float(u.mean()), "frac_empty2": float((u == 0).mean()), "rho_obs": R4.density_ratio(u, s), "mean_genes2": float(ng.mean())}
def vals_8um(u8, g8):
    return {"mean_umi8": float(u8.mean()), "mean_genes8": float(g8.mean())}
def vals_nuc(area, n, tissue_um2):
    return {"nuc_area_median": float(np.median(area)), "nuc_density": float(n / (tissue_um2 * 1e-6))}

Frows, G2, G3, V, LEAK = [], [], [], [], []
tp8_cache = {}
for A in ARCH:
    g = rois[rois.archetype == A].reset_index(drop=True); ids = list(g.roi_id); n = len(g)
    boxes = g[["c0", "r0", "c1", "r1"]].to_numpy(float); upp = g.um_per_px.to_numpy()
    Fb = C.box_overlap_matrix(boxes)
    # ---- 2 um: bin per area dai file di R4, posizione globale in pixel interi
    bins = []
    for k, rid in enumerate(ids):
        b = np.load(f"{R4RUN}/roi/{A}_{rid}_bins.npz")
        cc = np.rint(b["x_um"] / upp[k]).astype(np.int64); rr = np.rint(b["y_um"] / upp[k]).astype(np.int64)
        bins.append(dict(b=b, gx=boxes[k, 0] + cc, gy=boxes[k, 1] + rr))
    # F_bin: bin unici dell'unione, appartenenza per contenimento
    keys = np.concatenate([bb["gx"].astype(np.int64) * 10**7 + bb["gy"].astype(np.int64) for bb in bins])
    uk, first = np.unique(keys, return_index=True)
    ux = (uk // 10**7).astype(float); uy = (uk % 10**7).astype(float)
    mem = (ux[:, None] >= boxes[None, :, 0]) & (ux[:, None] < boxes[None, :, 2]) & (uy[:, None] >= boxes[None, :, 1]) & (uy[:, None] < boxes[None, :, 3])
    # bin che stanno in un file di area ma che per contenimento risultano altrove (devono essere 0)
    origin_ok = all(mem[np.searchsorted(uk, bb["gx"].astype(np.int64) * 10**7 + bb["gy"].astype(np.int64)), k].all() for k, bb in enumerate(bins))
    # i bin in tessuto di un'area devono comparire anche nel file dell'altra area quando cadono nella sua casella (stessa lista Space Ranger)
    count_in_files = np.zeros((len(uk), n), bool)
    for k, bb in enumerate(bins):
        count_in_files[np.searchsorted(uk, bb["gx"].astype(np.int64) * 10**7 + bb["gy"].astype(np.int64)), k] = True
    files_vs_boxes = int((count_in_files != mem).sum())
    Fbin = C.unit_overlap_matrix(mem)
    piece_u = C.assign_pieces(ux, uy, boxes)
    # ---- nuclei: centroide globale
    nuc = []
    for k, rid in enumerate(ids):
        d = nuc_all[(nuc_all.archetype == A) & (nuc_all.roi_id == rid)]
        nuc.append(dict(area=d.area_um2.to_numpy(), gx=boxes[k, 0] + d.x_um.to_numpy() / upp[k], gy=boxes[k, 1] + d.y_um.to_numpy() / upp[k]))
    nin = np.zeros((n, n))
    for i in range(n):
        x, y = nuc[i]["gx"], nuc[i]["gy"]
        for j in range(n):
            nin[i, j] = np.sum((x >= boxes[j, 0]) & (x < boxes[j, 2]) & (y >= boxes[j, 1]) & (y < boxes[j, 3]))
    Fnuc = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        nij = (nin[i, j] + nin[j, i]) / 2; Fnuc[i, j] = Fnuc[j, i] = nij / np.sqrt(nin[i, i] * nin[j, j])
    for i, j in itertools.combinations(range(n), 2):
        Frows.append(dict(archetype=A, i=ids[i], j=ids[j], F_box=Fb[i, j], F_bin=Fbin[i, j], F_nuc=Fnuc[i, j]))
        e0 = E0.get((A, ids[i], ids[j]), 0.0)
        G2.append(dict(check="G-2a F_box vs E-0/BL-215", archetype=A, pair=f"{ids[i]}-{ids[j]}", value=Fb[i, j], expected=e0, ok=abs(Fb[i, j] - e0) <= 0.005))
    G2.append(dict(check="G-2b bin dei file di R4 dentro la propria area", archetype=A, pair="", value=float(origin_ok), expected=1.0, ok=bool(origin_ok)))
    G2.append(dict(check="G-2c bin in tessuto: file per area = contenimento nelle caselle", archetype=A, pair="", value=files_vs_boxes, expected=0, ok=files_vs_boxes == 0))
    G2.append(dict(check="G-2d ogni bin unico in un solo pezzo, unione completa", archetype=A, pair="", value=int((piece_u < 0).sum()), expected=0,
                   ok=bool((piece_u >= 0).all() and mem[np.arange(len(uk)), piece_u].all())))
    # ---- 8 um: coordinate dei bin nello stesso ordine di R4_roi.py
    ds = g.dataset[0]
    if ds not in tp8_cache:
        sq = f"{DR}/{ds}/binned_outputs_x/binned_outputs/square_008um"
        with h5py.File(f"{sq}/filtered_feature_bc_matrix.h5", "r") as h: bc8 = h["matrix"]["barcodes"][:].astype(str)
        tp8 = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
        tp8_cache[ds] = tp8[tp8.barcode.isin(set(bc8))]
    tp8 = tp8_cache[ds]
    # ---- valori per area e per pezzo
    mask_area = {}
    for k, rid in enumerate(ids):
        b = bins[k]["b"]; u, s, ng = b["umi"], b["s"], b["ngen"]
        pk = C.assign_pieces(bins[k]["gx"], bins[k]["gy"], boxes) == k
        R = g.iloc[k]
        b8 = tp8[(tp8.pxl_row_in_fullres >= R.r0) & (tp8.pxl_row_in_fullres < R.r1) & (tp8.pxl_col_in_fullres >= R.c0) & (tp8.pxl_col_in_fullres < R.c1)]
        if len(b8) != len(b["umi8"]): raise SystemExit(f"{A}_{rid}: 8 um {len(b8)} contro {len(b['umi8'])}")
        p8 = C.assign_pieces(np.floor(b8.pxl_col_in_fullres.to_numpy()), np.floor(b8.pxl_row_in_fullres.to_numpy()), boxes) == k
        # tessuto valido (maschera di R2 a 1/4), centro del pixel in coordinate globali
        m = np.array(Image.open(f"results/R2/tissue_masks/{A}_{rid}_valid_ds4.png")) > 0
        jj, ii = np.meshgrid(np.arange(m.shape[1]), np.arange(m.shape[0]))
        mx = boxes[k, 0] + 4 * jj[m] + 2; my = boxes[k, 1] + 4 * ii[m] + 2
        pm = C.assign_pieces(mx, my, boxes) == k
        px_um2 = (4 * upp[k]) ** 2
        nk = nuc[k]; pn = C.assign_pieces(nk["gx"], nk["gy"], boxes) == k
        va = {**vals_2um(u, s, ng), **vals_8um(b["umi8"], b["ngen8"]), **vals_nuc(nk["area"], len(nk["area"]), m.sum() * px_um2)}
        vp = {**vals_2um(u[pk], s[pk], ng[pk]), **vals_8um(b["umi8"][p8], b["ngen8"][p8]), **vals_nuc(nk["area"][pn], pn.sum(), pm.sum() * px_um2)}
        for q in va:
            V.append(dict(archetype=A, roi_id=rid, quantity=q, area_value=va[q], piece_value=vp[q]))
        fr = dict(bins2=pk.mean(), bins8=p8.mean(), nuclei=pn.mean(), tissue=pm.mean())
        V.append(dict(archetype=A, roi_id=rid, quantity="piece_share_nuclei", area_value=1.0, piece_value=fr["nuclei"]))
        V.append(dict(archetype=A, roi_id=rid, quantity="piece_share_bins2", area_value=1.0, piece_value=fr["bins2"]))
        # G-3: riproduzione di R4_roi.csv
        r4 = R4roi[(R4roi.archetype == A) & (R4roi.roi_id == rid)].iloc[0]
        for q in ["mean_umi2", "frac_empty2", "rho_obs", "mean_genes2", "mean_umi8", "mean_genes8"]:
            G3.append(dict(check="G-3a valore per area = R4_roi.csv", archetype=A, roi_id=rid, quantity=q, ours=va[q], ref=float(r4[q]),
                           rel=abs(va[q] - r4[q]) / abs(r4[q]), ok=abs(va[q] - r4[q]) <= 1e-9 * abs(r4[q])))
        s2 = summ[(summ.archetype == A) & (summ.roi_id == rid)].iloc[0]
        for q, col in (("nuc_area_median", "area_median"), ("nuc_density", "density_per_mm2")):
            G3.append(dict(check="info: nuclei per area contro R2_roi_summary (maschera a 1/4 per la densità)", archetype=A, roi_id=rid, quantity=q,
                           ours=va[q], ref=float(s2[col]), rel=abs(va[q] - s2[col]) / abs(s2[col]), ok=True))
        # M-d: quota dei nuclei dell'area che stanno anche in un'altra area dello stesso tessuto
        other = np.zeros(len(nk["gx"]), bool)
        for j in range(n):
            if j != k: other |= (nk["gx"] >= boxes[j, 0]) & (nk["gx"] < boxes[j, 2]) & (nk["gy"] >= boxes[j, 1]) & (nk["gy"] < boxes[j, 3])
        LEAK.append(dict(archetype=A, held_out=rid, n_nuclei=len(other), share_in_training_areas=float(other.mean())))
    if A == "A6":
        Va = pd.DataFrame([v for v in V if v["archetype"] == "A6"])
        G2.append(dict(check="G-2e cuore: valori sui pezzi identici a quelli sulle aree", archetype=A, pair="", value=int((Va.area_value != Va.piece_value).sum()),
                       expected=0, ok=bool((Va.area_value == Va.piece_value).all())))
    print(A, "ok", flush=True)

# G-3b: intervalli registrati di R4 (uno stadio) dai valori per area con i semi di R4
rng_for = lambda a, q: np.random.default_rng([R4.SEED, AIDX[a], q])
for a in ARCH:
    Ra = R4roi[R4roi.archetype == a]
    jobs = [("rho_obs nucleus/ring UMI per bin (uncorrected)", "rho_obs", 3)] + [(col, col, 100 + qi) for qi, col in enumerate(R4COLS)]
    for name, col, qq in jobs:
        e, lo, hi, _ = R4.boot_one_stage(Ra[col].values, rng=rng_for(a, qq))
        rr = R4refs[(R4refs.archetype == a) & (R4refs.quantity == name)].iloc[0]
        ok = all(abs(x - y) <= 1e-9 * max(abs(y), 1e-12) for x, y in ((e, rr.estimate), (lo, rr.lo), (hi, rr.hi)))
        G3.append(dict(check="G-3b intervallo di R4_refs.csv ricalcolato", archetype=a, roi_id="", quantity=name, ours=lo, ref=rr.lo,
                       rel=max(abs(lo - rr.lo), abs(hi - rr.hi)) / max(abs(rr.hi - rr.lo), 1e-12), ok=ok))

# M-f: bootstrap a due stadi (aree, cellule) contro uno stadio sulle mediane per area, prima riga a due stadi di R4
MF = []
for a in ARCH:
    cells = []
    for rid in rois[rois.archetype == a].roi_id:
        ce = pd.read_parquet(f"{R4RUN}/roi/{a}_{rid}_cells.parquet"); ce = ce[ce.interior] if "interior" in ce else ce
        cells.append(ce.umi_exp.to_numpy(float))
    e2, lo2, hi2, _ = R4.boot_two_stage(cells, np.median, rng=np.random.default_rng([C.SEED, AIDX[a], 1]))
    med = np.array([np.median(x) for x in cells])
    e1, lo1, hi1, _ = R4.boot_one_stage(med, stat=np.mean, rng=np.random.default_rng([C.SEED, AIDX[a], 2]))
    MF.append(dict(archetype=a, quantity="median UMI per interior cell (nucleus + ring)", two_stage_lo=lo2, two_stage_hi=hi2, two_stage_est=e2,
                   one_stage_lo=lo1, one_stage_hi=hi1, one_stage_est=e1, width_ratio=(hi2 - lo2) / (hi1 - lo1)))

pd.DataFrame(Frows).to_csv(f"{OUT}/R7_F_matrices.csv", index=False)
pd.DataFrame(G2).to_csv(f"{OUT}/R7_G2.csv", index=False)
pd.DataFrame(G3).to_csv(f"{OUT}/R7_G3.csv", index=False)
pd.DataFrame(V).to_csv(f"{OUT}/R7_values_area_piece.csv", index=False)
pd.DataFrame(LEAK).to_csv(f"{OUT}/R7_Md_leakage.csv", index=False)
pd.DataFrame(MF).to_csv(f"{OUT}/R7_Mf_two_stage.csv", index=False)
g2 = pd.DataFrame(G2); g3 = pd.DataFrame(G3); g3c = g3[~g3.check.str.startswith("info")]
print(f"G-2 {'PASS' if g2.ok.all() else 'FAIL'} ({int(g2.ok.sum())}/{len(g2)}); G-3 {'PASS' if g3c.ok.all() else 'FAIL'} ({int(g3c.ok.sum())}/{len(g3c)})")
