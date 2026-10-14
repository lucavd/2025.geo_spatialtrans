
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="results/S1.3c/review/method/out/"
d=pd.read_csv(D+"edge_R.csv")
cols=[("R_all","R osservato (= R_profile)"),("R_edge20_roi20","R, bin a ≥20 µm da bordo tessuto e ROI"),("Rcomp_all","R atteso dal nullo (densità uniforme nel territorio)"),("Rwithin_all","R a parità di territorio"),("Rwithin_edge20_roi20","R a parità di territorio, ≥20 µm dai bordi")]
fig,axes=plt.subplots(1,6,figsize=(16,4.2),sharey=True,facecolor="white")
for ax,A in zip(axes,[f"A{i}" for i in range(1,7)]):
    g=d[d.archetype==A]
    for j,(c,l) in enumerate(cols):
        ax.scatter(np.full(len(g),j)+np.linspace(-0.15,0.15,len(g)),g[c],s=18,label=l if A=="A1" else None,color=f"C{j}")
    ax.axhline(1,color="grey",lw=0.8,ls="--"); ax.set_title(A); ax.set_xticks(range(len(cols))); ax.set_xticklabels(["oss","≥20µm","nullo","terr","terr≥20"],rotation=60,fontsize=8)
axes[0].set_ylabel("R = UMI medio a d_nuc 4–8 µm / 0–2 µm")
fig.legend(loc="lower center",ncol=3,fontsize=8,frameon=False); fig.tight_layout(rect=(0,0.12,1,1))
fig.savefig(D+"F_rev_R_decomposition.png",dpi=120)
print("ok")
