
import pandas as pd, numpy as np, os
M="/mnt/micron/geo_spatialtrans/R6"
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
NON = {"non classificabile","non risolto"}
print("==== (d) PH-1 cluster counts and labels")
for A in ["A1","A2","A3","A4","A5","A6"]:
    base=f"{M}/labels/{A}_exp_umi20_seed20261022"
    lab = pd.read_parquet(base+"_nuclei_ph1.parquet"); cl = pd.read_csv(base+"_clusters_ph1.csv"); cl0 = pd.read_csv(base+"_clusters.csv")
    cnt = lab[lab.cluster.notna()].groupby("cluster").size(); cnt.index = cnt.index.astype(int)
    cl = cl.set_index("cluster"); cl0 = cl0.set_index("cluster")
    ok_n = cnt.reindex(cl.index).equals(cl.n_nuclei.astype(int)) ; ok_n0 = cnt.reindex(cl0.index).equals(cl0.n_nuclei.astype(int))
    # type per cluster unique and equal to type_ph1 / type
    tpc = lab[lab.cluster.notna()].groupby("cluster").cell_type.agg(lambda s: s.unique().tolist()); tpc.index = tpc.index.astype(int)
    uni = all(len(v)==1 for v in tpc.values)
    same = all(tpc.loc[k][0] == cl.type_ph1.loc[k] for k in cl.index)
    tpc0 = lab[lab.cluster.notna()].groupby("cluster").cell_type_prereg.agg(lambda s: s.unique().tolist()); tpc0.index = tpc0.index.astype(int)
    same0 = all(tpc0.loc[k][0] == cl0.type.loc[k] for k in cl0.index)
    nan_cl = lab.cluster.isna(); nc_ok = (lab.cell_type[nan_cl] == "non classificabile").all() and (lab.cell_type_prereg[nan_cl] == "non classificabile").all() and (~lab.classifiable[nan_cl]).all() and (lab.classifiable[~nan_cl]).all()
    nc_only_nan = (lab.cell_type == "non classificabile").sum() == nan_cl.sum()
    print(f"{A}: clusters csv {len(cl)} / parquet {cnt.size}; n_nuclei per cluster csv==parquet: ph1 {ok_n}, prereg {ok_n0}; type unique per cluster {uni}; type_ph1 matches {same}; prereg type matches {same0}; NaN cluster -> non classificabile & !classifiable {nc_ok}; 'non classificabile' only for NaN cluster {nc_only_nan}; n NaN {int(nan_cl.sum())} of {len(lab)}; changed ph1 vs prereg {(lab.cell_type!=lab.cell_type_prereg).mean():.3f}")
print("\n==== (e) territory join coverage and duplicates")
for A in ["A1","A2","A3","A4","A5","A6"]:
    lab = pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei.parquet")
    parts=[]
    for roi in lab.roi_id.unique():
        t = pd.read_parquet(f"{M}/terr/{A}_{roi}.parquet"); t["roi_id"]=roi; parts.append(t)
    t = pd.concat(parts, ignore_index=True)
    t["kx"]=t.x.round(4); t["ky"]=t.y.round(4); lab["kx"]=lab.x_um.round(4); lab["ky"]=lab.y_um.round(4)
    dup_t = t.duplicated(["roi_id","kx","ky"]).sum(); dup_l = lab.duplicated(["roi_id","kx","ky"]).sum()
    mg = lab.merge(t[["roi_id","kx","ky","a_t","idx"]], on=["roi_id","kx","ky"], how="left")
    cov = mg.a_t.notna().mean(); 
    # exact join on label == idx?
    mg2 = lab.merge(t[["roi_id","idx","a_t"]], left_on=["roi_id","label"], right_on=["roi_id","idx"], how="left")
    unm = mg[mg.a_t.isna()]
    print(f"{A}: nuclei {len(lab)} territories {len(t)}; join (x,y 4 dp) coverage {cov:.4f} unmatched {int(mg.a_t.isna().sum())}; rows after merge {len(mg)} (dup keys terr {dup_t}, lab {dup_l}); join on label==idx coverage {mg2.a_t.notna().mean():.4f}; label==idx identical ordering within roi: {(mg.idx.fillna(-1).astype(int)==lab.label).mean():.4f}")
    if len(unm): print("   unmatched per roi:", unm.roi_id.value_counts().to_dict(), "; sample", unm[["roi_id","label","x_um","y_um"]].head(3).to_dict("records"))
