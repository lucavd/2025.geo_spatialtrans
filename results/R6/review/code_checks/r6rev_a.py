
import pandas as pd, numpy as np, scipy.sparse as sp
M="/mnt/micron/geo_spatialtrans/R6"
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
bad=0; rows=[]
for A in ["A1","A2","A3","A4","A5","A6"]:
    rr = rois[rois.archetype==A].roi_id.tolist()
    lab = pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei.parquet")
    lab1 = pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    order_ok = rr == list(lab.roi_id.unique())
    cat = pd.concat([pd.read_parquet(f"{M}/counts/{A}_{r}_nuclei.parquet") for r in rr], ignore_index=True)
    same_rows = (len(cat)==len(lab)) and (cat.label.values==lab.label.values).all() and (cat.roi_id.values==lab.roi_id.values).all() and np.allclose(cat.umi_exp.values, lab.umi.values)
    same_ph1 = (len(lab1)==len(lab)) and (lab1.label.values==lab.label.values).all() and (lab1.roi_id.values==lab.roi_id.values).all() and (lab1.cell_type_prereg.values==lab.cell_type.values).all()
    for r in rr:
        t = pd.read_parquet(f"{M}/counts/{A}_{r}_nuclei.parquet")
        for fp in ["exp","nuc","shift"]:
            C = sp.load_npz(f"{M}/counts/{A}_{r}_{fp}.npz")
            rs = np.asarray(C.sum(axis=1)).ravel()
            ok = C.shape[0]==len(t) and np.array_equal(rs, t[f"umi_{fp}"].values)
            if not ok: bad+=1
            rows.append((A,r,fp,C.shape[0],len(t),ok, float(np.abs(rs-t[f"umi_{fp}"].values).max()) if C.shape[0]==len(t) else np.nan))
    print(A, "rois csv order == lab.roi_id.unique():", order_ok, "| concat counts == labels parquet (label, roi_id, umi):", same_rows, "| ph1 aligned:", same_ph1, "| n", len(lab), "| label unique within roi:", not cat.duplicated(["roi_id","label"]).any())
df = pd.DataFrame(rows, columns=["A","roi","fp","nrow_npz","nrow_parquet","ok","maxabsdiff"])
print("ROI x footprint checked:", len(df), "failures:", bad)
print(df[df.fp=="exp"].to_string())
