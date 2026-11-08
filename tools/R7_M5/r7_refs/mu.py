
import pandas as pd
d = pd.read_csv("/tmp/r7_refs/R7_Mc_refs.csv"); o = pd.read_csv("/tmp/r7_refs/R7_Mc_refs_v2.csv"); u = pd.read_csv("/tmp/r7_refs/use.csv")
ids = u[u.margin_use.astype(str).str.upper() == "TRUE"].id.tolist(); print("margin TRUE ids", ids)
x = d[d.id.isin(ids)].merge(o[["id","metric","lo_new","hi_new"]], on=["id","metric"], how="left", suffixes=("","_v2"))
pd.set_option("display.width", 250)
print(x[["id","archetype","interval_construction","metric","lo","hi","lo_new_v2","hi_new_v2","kappa","lo_new","hi_new","lo_stress","hi_stress","lo_geomA5","hi_geomA5"]].round(4).to_string(max_colwidth=38))
print(u[u.id.isin(ids)][["id","margin_evidence"]].to_string(max_colwidth=200))
