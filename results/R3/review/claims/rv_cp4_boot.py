
import numpy as np, pandas as pd, glob
RV="results/R3/review/claims"; OUT="/mnt/micron/geo_spatialtrans/R3"
cc=pd.concat([pd.read_parquet(f) for f in glob.glob(f"{OUT}/calib/*.parquet")]); cc=cc[cc.interior]
cc["A"]=cc.win_id.str[:2]; cc["k"]=np.where(cc.source=="manual_Luca","m","s")
rng=np.random.default_rng(1); cv=lambda a: a.std(ddof=1)/a.mean(); rows=[]
for A in sorted(cc.A.unique()):
    m=cc.area[(cc.A==A)&(cc.k=="m")].values; s=cc.area[(cc.A==A)&(cc.k=="s")].values
    b=[cv(rng.choice(s,len(s)))-cv(rng.choice(m,len(m))) for _ in range(2000)]
    rows.append(dict(archetype=A, diff=cv(s)-cv(m), ci_lo=np.quantile(b,.025), ci_hi=np.quantile(b,.975)))
o=pd.DataFrame(rows); print(o.round(3).to_string()); o.to_csv(f"{RV}/rv_cp4_bootstrap.csv",index=False)
V=pd.read_csv("results/R3/R3_verdicts.csv"); pc=V.prediction_confirmed
print("rows",len(V),"PASS",(V.assertion=="PASS").sum(),"FAIL",(V.assertion=="FAIL").sum(),"registrato",(V.assertion=="registrato").sum(),
      "pred TRUE",(pc==True).sum(),"pred FALSE",(pc==False).sum(),"pred NA",pc.isna().sum())
print(V[pc==False][["id","archetype"]].to_string())

