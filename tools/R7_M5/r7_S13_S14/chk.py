
import pandas as pd
R="/home/user/2025.geo_spatialtrans/results/"
v=pd.read_csv(R+"S1.4e/S1.4e_verdicts.csv"); s=pd.read_csv(R+"R2c/R2c_s1.csv"); s=s[s.set=="Luca"].set_index("archetype")
kb=1.979389361922006
def vd(lo,hi,a,b): return "confermata" if lo>=a and hi<=b else ("smentita" if hi<a or lo>b else "indeterminata")
sc={}
for arm in ["K2","ER","ET"]:
    n=0; det=[]
    for _,r in v[(v.arm==arm)&(v.id.isin(["P-1","P-2"]))].iterrows():
        lo,hi=r.est-2*r.se, r.est+2*r.se
        if r.id=="P-1":
            c=s.loc[r.archetype,"s1_primary"]; a=c-kb*(c-s.loc[r.archetype,"lo_primary"]); b=c+kb*(s.loc[r.archetype,"hi_primary"]-c)
        else: a,b=r.lo,r.hi
        x=vd(lo,hi,a,b); n+= x=="confermata"; det.append("%s%s:%s"%(r.id,r.archetype,x[:4]))
    sc[arm]=n; print(arm,n," ".join(det))
print(sc)
