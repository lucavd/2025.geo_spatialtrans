import pandas as pd, numpy as np, glob, os
from scipy.stats import spearmanr
pd.set_option("display.width",250)
fs=[f for f in glob.glob("results/R3/review/code/cells/*.parquet") if "_jit" not in f]
L=[]
for f in fs:
    A,roi,*m=os.path.basename(f)[:-8].split("_"); d=pd.read_parquet(f); d["archetype"]=A; d["roi_id"]=roi; d["method"]="_".join(m); L.append(d)
X=pd.concat(L); I=X[X.interior].copy()
I["eq_r"]=np.sqrt(I.area/np.pi); I["nc"]=I.area_nuc/I.area; I["ratio"]=np.sqrt(I.nc)
def axd(a,b):
    d=np.abs(a-b)%np.pi; return np.minimum(d,np.pi-d)
rows=[]
for (A,me),d in I.groupby(["archetype","method"]):
    rows.append(dict(archetype=A,method=me,n_interior=len(d),area_median=d.area.median(),eq_r_median=d.eq_r.median(),nc_median=d.nc.median(),
      ratio_median=d.ratio.median(),ecc_T_median=d.ecc_T.median(),spearman_nuc_terr=spearmanr(d.area_nuc,d.area)[0],frac_cut=np.mean(d.frac_out>0.05),
      max_abs_eccN_minus_regionprops=np.nanmax(np.abs(d.ecc_N-d.ecc_pq)),
      p95_thetaN_vs_regionprops_deg=np.nanquantile(np.degrees(axd(d.theta_N,d.theta_pq))[d.ecc_N>0.3],0.95),
      n_elig_pixel=int(((d.ecc_N>=0.8)&(d.ecc_T>=0.5)).sum()), n_elig_regionprops=int(((d.ecc_pq>=0.8)&(d.ecc_T>=0.5)).sum())))
mine=pd.DataFrame(rows)
th=pd.read_csv("results/R3/R3_cells_summary.csv")
M=mine.merge(th,on=["archetype","method"],suffixes=("","_R3"))
out=[]
for c in ["n_interior","area_median","eq_r_median","nc_median","ratio_median","ecc_T_median","spearman_nuc_terr","frac_cut"]:
    M[c+"_reldiff"]=(M[c]-M[c+"_R3"]).abs()/M[c+"_R3"].abs()
print(M[["archetype","method"]+[c for c in M.columns if c.endswith("_reldiff")]].to_string())
print(mine[["archetype","method","frac_cut","max_abs_eccN_minus_regionprops","p95_thetaN_vs_regionprops_deg","n_elig_pixel","n_elig_regionprops"]].round(4).to_string())
M.to_csv("results/R3/review/code/rc_compare_archetype.csv",index=False)
