
import pandas as pd, numpy as np, os, sys
M="/mnt/micron/geo_spatialtrans/R6"
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
src = open("tools/R6_analysis.py").read().split("# ================= main")[0]
# execute only the function definitions (no side effects except reading refs csv)
exec(src, globals())
NON = {"non classificabile","non risolto"}
for SUF in ["", "_ph1"]:
    print(f"\n==== variant '{SUF or 'prereg'}'")
    refs_t = pd.read_csv(f"results/R6/R6_refs_type{SUF}.csv"); pred = pd.read_csv(f"results/R6/R6_predictions{SUF}.csv")
    labs = {}
    for A in ["A1","A2","A3","A4","A5","A6"]:
        lab = pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei{SUF}.parquet"); labs[A]=lab
        # refs_type: fraction (denominator all nuclei), median area, median ecc, n
        r = refs_t[refs_t.archetype==A].set_index("cell_type")
        mine = lab.groupby("cell_type").agg(n=("label","size"), median_area_nuc=("area_um2","median"), median_ecc_nuc=("eccentricity","median"))
        mine["fraction"] = mine.n/len(lab)
        j = r.join(mine, rsuffix="_mine")
        dn = (j.n != j.n_mine).sum(); dfr = np.abs(j.fraction-j.fraction_mine).max(); da = np.abs(j.median_area_nuc-j.median_area_nuc_mine).max(); de = np.abs(j.median_ecc_nuc-j.median_ecc_nuc_mine).max()
        print(f"{A}: refs_type rows {len(r)} types {len(mine)}; n mismatches {dn}; max|dfraction| {dfr:.2e}; max|darea| {da:.2e}; max|decc| {de:.2e}; fraction sums to {r.fraction.sum():.4f} (incl. nontype {r.fraction.reindex(list(NON)).fillna(0).sum():.3f})")
    # predictions
    p = pred.set_index(["id","archetype"]).estimate
    out = []
    l5 = labs["A5"]; est = l5.area_um2[l5.cell_type=="neuroni"].median()/l5.area_um2[l5.cell_type.isin(GLIA)].median(); out.append(("P-1","A5",p.get(("P-1","A5")), est))
    for A,lab in labs.items(): out.append(("P-2",A,p.get(("P-2",A)), ratio_max_min(lab)))
    l6 = labs["A6"]; d = l6[~l6.cell_type.isin(NON)]; out.append(("P-3","A6",p.get(("P-3","A6")), d.eccentricity[d.cell_type=="cardiomiociti"].median()-d.eccentricity[d.cell_type!="cardiomiociti"].median()))
    l4 = labs["A4"]; out.append(("P-4","A4",p.get(("P-4","A4")), (l4.cell_type[l4.roi_id.isin(["f1","f2"])]=="linfociti B").mean()-(l4.cell_type[l4.roi_id.isin(["p1","p2","p3"])]=="linfociti B").mean()))
    DOM = {"A1": ("epiteliale", None), "A2": ("epiteliale tumorale", None), "A3": (["fibroblasti","muscolo liscio"], None), "A4": ("linfociti B", ["f1","f2"]), "A5": ("neuroni", None), "A6": ("cardiomiociti", None)}
    for A,(ty,rr) in DOM.items():
        lab = labs[A] if rr is None else labs[A][labs[A].roi_id.isin(rr)]
        out.append(("P-5",A,p.get(("P-5",A)), lab.cell_type.isin(ty if isinstance(ty,list) else [ty]).mean()))
    for A in ["A5","A6"]:
        lab = add_sr(labs[A], A); print(f"   {A}: rows before/after add_sr {len(labs[A])}/{len(lab)}; sr_cell_id notna {int(labs[A].sr_cell_id.notna().sum())}; merged sr_type notna {int(lab.sr_type.notna().sum())}; lowumi {int(lab.sr_lowumi.fillna(False).sum())}")
        out.append(("P-6",A,p.get(("P-6",A)), agreement_informative(lab))); out.append(("P-6 (all matched)",A,p.get(("P-6 (all matched)",A)), agreement_all(lab)))
        if SUF=="": 
            ct = pd.read_csv(SRFILE[A]); ct["sr_type"]=[sr_map(A,c,f) for c,f in zip(ct.coarse_cell_type.astype(str), ct.fine_cell_type.astype(str))]
            print(f"   {A} sr_map coarse->type:"); print(pd.crosstab(ct.coarse_cell_type, ct.sr_type).to_string())
    o = pd.DataFrame(out, columns=["id","A","csv_estimate","recomputed"]); o["absdiff"]=(o.csv_estimate-o.recomputed).abs()
    print(o.to_string()); print("max abs diff:", o.absdiff.max())
