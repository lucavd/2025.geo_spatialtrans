"""tools/S2.5_real.py — S2.5: valori veri delle 30 aree per il confronto del simulatore completo (docs/steps/S2.5.md).
Due definizioni del tessuto, stesse misure:
  reg : quadratini con il centro (allineato con lo scarto della sezione, come S2.3 e S2.4) nelle regioni dello Step 1 e a >= 25 um dal
        bordo dell'area (come S2.4 E-1 e il lato simulato); blocchi di 8 um con almeno 8 dei 16 quadratini in `reg`;
  sr  : quadratini e blocchi marcati in tessuto da Space Ranger nel riquadro dell'area (definizione di R4, B-270...B-331).
Misure: UMI medio, frazione di vuoti, geni rilevati medi a 2 e 8 um (matrici raw di Space Ranger, tutte le feature, come R4).
Rapporto nucleo/corona (statistica di R4 PH-1 con la maschera spostata): dal lato vero di S2.3 (panel_star, definizione reg).
Anisotropia del campo dei blocchi di 8 um (tools/S2_5_common.py aniso, definizione reg).
Uso: .venv/bin/python tools/S2.5_real.py areas <outdir> | summary <outdir> | bl219 <outdir>"""
import sys, os, json, time, glob, traceback
import numpy as np, pandas as pd, pyarrow.parquet as pq, h5py
sys.path.insert(0, "tools")
import S2_5_common as CM
import importlib.util
def _load(name, path):
    sp = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m

M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
S23REAL = f"{M}/S2.3/e2/real"
FRAME_UM = 25.0
MIN_CHILD8 = 8
B = 2000; SEED = 20261121
METRICS = ["mean_umi2", "frac_empty2", "mean_genes2", "mean_umi8", "frac_empty8", "mean_genes8"]

def colsum_h5(p, keep_ids):
    """UMI e geni rilevati per colonna della matrice raw, sui soli geni `keep_ids` (le feature della matrice filtrata = pannello di sonde,
    come R4; la raw ha anche i geni fuori pannello, S2.5 E-0)."""
    with h5py.File(p, "r") as h:
        m = h["matrix"]; ip = m["indptr"][:]; bc = m["barcodes"][:].astype(str); dat = m["data"][:]; ind = m["indices"][:]
        ids = m["features"]["id"][:].astype(str)
    keep = np.isin(ids, keep_ids); kk = keep[ind]
    cum = np.concatenate([[0], np.cumsum(np.where(kk, dat, 0), dtype=np.int64)]); cnz = np.concatenate([[0], np.cumsum(kk, dtype=np.int64)])
    cs = cum[ip[1:]] - cum[ip[:-1]]; nz = cnz[ip[1:]] - cnz[ip[:-1]]   # somme per colonna (anche con colonne vuote in coda)
    return pd.DataFrame({"barcode": bc, "umi": cs, "ngen": nz})

def filtered_ids(p):
    with h5py.File(p, "r") as h: return h["matrix"]["features"]["id"][:].astype(str)

stats = CM.stats

def job_section(key, OUT):
    import shapely
    RL = _load("s23real", "tools/S2.3_real.py")
    t0 = time.time()
    sq = f"{D}/{DS[key]}/binned_outputs_x/binned_outputs"
    fid = filtered_ids(f"{sq}/square_002um/filtered_feature_bc_matrix.h5")
    c2 = colsum_h5(f"{sq}/square_002um/raw_feature_bc_matrix.h5", fid); c8 = colsum_h5(f"{sq}/square_008um/raw_feature_bc_matrix.h5", fid)
    tp2 = pq.read_table(f"{sq}/square_002um/spatial/tissue_positions.parquet").to_pandas().merge(c2, on="barcode", how="left")
    tp8 = pq.read_table(f"{sq}/square_008um/spatial/tissue_positions.parquet").to_pandas().merge(c8, on="barcode", how="left")
    for t in (tp2, tp8): t["umi"] = t["umi"].fillna(0).astype(np.int64); t["ngen"] = t["ngen"].fillna(0).astype(np.int64)
    t8 = tp8.set_index(["array_row", "array_col"])
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype.map(KEY) == key]
    PJ = json.load(open(f"{RL.S20}/params_flat.json")); r4 = pd.read_csv("results/R4/R4_roi.csv")
    for _, Rw in rois.iterrows():
        A = Rw.archetype; tag = f"{A}_{Rw.roi_id}"; fo = f"{OUT}/real/{tag}.json"
        if os.path.exists(fo): continue
        try:
            prm = PJ[RL.NC.KEY[A]]["C2"]; Dd = RL.NC.roi_data(Rw); di = RL.NC.d_image(Dd, prm)
            upp = Rw.um_per_px; side = float(Rw.side_px * upp)
            m = (tp2.pxl_row_in_fullres >= Rw.r0) & (tp2.pxl_row_in_fullres < Rw.r1) & (tp2.pxl_col_in_fullres >= Rw.c0) & (tp2.pxl_col_in_fullres < Rw.c1)
            b = tp2[m]
            pos = np.c_[(b.pxl_row_in_fullres.to_numpy() - Rw.r0) * upp, (b.pxl_col_in_fullres.to_numpy() - Rw.c0) * upp] - np.asarray(di)[None, :]
            R = RL.region_geom(tag)
            inreg = shapely.contains_xy(R, pos[:, 1], pos[:, 0])
            fr = (pos[:, 0] >= FRAME_UM) & (pos[:, 0] <= side - FRAME_UM) & (pos[:, 1] >= FRAME_UM) & (pos[:, 1] <= side - FRAME_UM)
            P = inreg & fr
            u2, g2 = b.umi.to_numpy(), b.ngen.to_numpy()
            res = {"tag": tag, "archetype": A, "roi_id": Rw.roi_id, "side_um": side, "d_img_um": [float(v) for v in di]}
            for k, v in stats(u2[P], g2[P]).items(): res[f"reg_{k}2"] = v
            sr = b.in_tissue.to_numpy() == 1
            for k, v in stats(u2[sr], g2[sr]).items(): res[f"sr_{k}2"] = v
            # blocchi di 8 um: (array_row // 4, array_col // 4) = indici del blocco di Space Ranger (identita' verificata in S2.5 E-0)
            blk = pd.DataFrame({"r8": b.array_row.to_numpy() // 4, "c8": b.array_col.to_numpy() // 4, "P": P, "u": u2})
            gb = blk.groupby(["r8", "c8"]).agg(nP=("P", "sum"), nch=("P", "size"), u=("u", "sum"))
            sel = gb[(gb.nP >= MIN_CHILD8)]
            full = gb[gb.nch == 16]
            look = t8.reindex(full.index)
            res["check8_n_full_blocks"] = int(len(full)); res["check8_identical"] = int((look.umi.to_numpy() == full.u.to_numpy()).sum())
            s8 = t8.reindex(sel.index)
            for k, v in stats(s8.umi.to_numpy(), s8.ngen.to_numpy()).items(): res[f"reg_{k}8"] = v
            sel.reset_index()[["r8", "c8"]].to_parquet(f"{OUT}/real/{tag}_blocks8.parquet", index=False)
            res["reg_n8_partial_children"] = int((sel.nch < 16).sum())
            # anisotropia del campo dei blocchi di 8 um in tessuto (definizione reg); immagine sugli indici (riga, colonna) del blocco
            r0, c0 = gb.index.get_level_values(0).min(), gb.index.get_level_values(1).min()
            nr8, nc8 = gb.index.get_level_values(0).max() - r0 + 1, gb.index.get_level_values(1).max() - c0 + 1
            U8 = np.zeros((nr8, nc8)); T8 = np.zeros((nr8, nc8), bool)
            ri, ci = sel.index.get_level_values(0) - r0, sel.index.get_level_values(1) - c0
            U8[ri, ci] = s8.umi.to_numpy(); T8[ri, ci] = True
            res["reg_aniso8"], res["reg_angle8"] = CM.aniso(U8, T8)
            m8 = (tp8.pxl_row_in_fullres >= Rw.r0) & (tp8.pxl_row_in_fullres < Rw.r1) & (tp8.pxl_col_in_fullres >= Rw.c0) & (tp8.pxl_col_in_fullres < Rw.c1) & (tp8.in_tissue == 1)
            for k, v in stats(tp8.umi[m8].to_numpy(), tp8.ngen[m8].to_numpy()).items(): res[f"sr_{k}8"] = v
            j = json.load(open(f"{S23REAL}/{tag}.json")); res["reg_star"] = float(j["panel_star"]); res["reg_star_n_bins"] = int(j["panel_n_bins"])
            q = r4[(r4.archetype == A) & (r4.roi_id == Rw.roi_id)].iloc[0]
            res["r4_check"] = {k: float(q[k]) for k in METRICS}
            res["sec"] = round(time.time() - t0, 1)
        except Exception:
            res = {"tag": tag, "error": traceback.format_exc()}
        json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo)
        print(f"[{key}] {tag} {res.get('sec')}", flush=True)
    return key

corrected_interval = CM.corrected_interval
SEED = 20261121

def job_bl219(OUT):
    """BL-219: quota degli UMI veri dello stroma (A3, definizione reg, blocchi di 8 um) sui geni del pannello senza profilo nell'atlante."""
    import shapely
    RL = _load("s23real", "tools/S2.3_real.py")
    sq = f"{D}/{DS['A23']}/binned_outputs_x/binned_outputs"
    prof = "/mnt/micron/geo_spatialtrans/S2.2/profiles"
    genes = pd.read_csv(f"{prof}/A3_genes.tsv", sep="\t"); sums = pd.read_csv(f"{prof}/A3_sums.tsv.gz", sep="\t")
    nz = np.zeros(len(genes), bool); nz[sums.gene_index.to_numpy() - 1] = True          # gene_index 1-based (R)
    missing = genes.gene_id[~nz].to_numpy()
    with h5py.File(f"{sq}/square_008um/raw_feature_bc_matrix.h5", "r") as h:
        m = h["matrix"]; ip = m["indptr"][:]; bc = m["barcodes"][:].astype(str); dat = m["data"][:]; ind = m["indices"][:]; ids = m["features"]["id"][:].astype(str)
    fid = filtered_ids(f"{sq}/square_002um/filtered_feature_bc_matrix.h5")
    inpanel = np.isin(ids, fid); ismiss = np.isin(ids, missing)
    tp8 = pq.read_table(f"{sq}/square_008um/spatial/tissue_positions.parquet").to_pandas(); col = pd.Series(np.arange(len(bc)), index=bc)
    out = []
    for f in sorted(glob.glob(f"{OUT}/real/A3_*.json")):
        j = json.load(open(f)); tag = j["tag"]
        blk = pd.read_parquet(f"{OUT}/real/{tag}_blocks8.parquet")
        b = tp8.set_index(["array_row", "array_col"]).reindex(pd.MultiIndex.from_arrays([blk.r8, blk.c8])).barcode
        cols = col.reindex(b.dropna()).dropna().astype(int).to_numpy()
        tot = 0; mis = 0
        for cc in cols:
            sl = slice(ip[cc], ip[cc + 1]); g = ind[sl]; d = dat[sl]
            tot += int(d[inpanel[g]].sum()); mis += int(d[ismiss[g]].sum())
        out.append({"tag": tag, "n_blocks": int(len(cols)), "umi_panel": tot, "umi_missing_genes": mis, "share": mis / tot if tot else np.nan})
    R = pd.DataFrame(out); R.to_csv(f"{OUT}/S2.5_BL219_share.csv", index=False)
    print("geni senza profilo nell'atlante (stroma):", len(missing), "di cui nel pannello:", int(np.isin(missing, fid).sum()))
    ci = corrected_interval(R.share.to_numpy(float), "A3", np.random.default_rng([SEED, 219]))
    print(R.to_string(index=False)); print(ci)
    json.dump({"n_missing": int(len(missing)), **ci}, open(f"{OUT}/S2.5_BL219_summary.json", "w"))

if __name__ == "__main__":
    mode, OUT = sys.argv[1], sys.argv[2]
    os.makedirs(f"{OUT}/real", exist_ok=True)
    if mode == "areas":
        from multiprocessing import Pool
        nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 3
        with Pool(min(nw, 5)) as p: print(p.starmap(job_section, [(k, OUT) for k in DS], chunksize=1))
    elif mode == "bl219":
        job_bl219(OUT)
    elif mode == "summary":
        J = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/real/*.json"))])
        assert "error" not in J.columns or J["error"].isna().all(), J.get("error")
        J.drop(columns=["r4_check", "d_img_um"]).to_csv(f"{OUT}/S2.5_real_areas.csv", index=False)
        # riproduzione dei valori per area di R4 con la definizione sr
        chk = []
        for _, r in J.iterrows():
            for mtr in METRICS:
                chk.append({"tag": r.tag, "metric": mtr, "r4": r.r4_check[mtr], "sr": r[f"sr_{mtr}"], "absdiff": abs(r.r4_check[mtr] - r[f"sr_{mtr}"])})
        pd.DataFrame(chk).to_csv(f"{OUT}/S2.5_real_R4_reproduction.csv", index=False)
        rows = []
        for A, g in J.groupby("archetype"):
            ai = ["A1", "A2", "A3", "A4", "A5", "A6"].index(A)
            for qi, mtr in enumerate(METRICS + ["star", "aniso8"]):
                for di_, dfn in enumerate(("reg", "sr")):
                    col = f"{dfn}_{mtr}"
                    if col not in g: continue
                    ci = corrected_interval(g[col].to_numpy(float), A, np.random.default_rng([SEED, ai, qi, di_]))
                    rows.append({"archetype": A, "metric": mtr, "definition": dfn, **ci})
        S = pd.DataFrame(rows); S.to_csv(f"{OUT}/S2.5_real_refs.csv", index=False)
        print(pd.DataFrame(chk).groupby("metric").absdiff.max())
        print(S.round(5).to_string(index=False))
