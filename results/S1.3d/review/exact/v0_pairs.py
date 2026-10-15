"""review/exact/v0_pairs.py — punto (c): caratterizza ogni coppia violata da V0 (167 righe o1_V0): area esatta della coppia (rx), aree esatte
delle due celle; il valore di V0 coincide con l'area di una delle due celle? Insiemi e copertura (coverage json v0_valid)."""
import sys, os, json, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np, rx
D = "/mnt/micron/geo_spatialtrans/S1.3d"
V = pd.read_csv(f"{D}/review_exact/viol_recheck_all.csv"); V = V[(V.vname == "o1_V0") & V.what.isin(["pair_ij", "pair_ji"])]
rows = []
for key, d in V.groupby("key"):
    ids = set(d.i.astype(int)) | set(d.j.astype(int)); cells = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want=ids)
    cv = json.load(open(f"{D}/coverage/{key}.json")) if os.path.exists(f"{D}/coverage/{key}.json") else {}
    for _, r in d.iterrows():
        i, j = int(r.i), int(r.j); ai, aj = rx.area(cells[i]), rx.area(cells[j])
        rows.append(dict(key=key, what=r.what, i=i, j=j, exact=r.exact, cand=r.cand, a_i=float(ai), a_j=float(aj), tau=r.tau,
                         cand_minus_ai=abs(r.cand - float(ai)), cand_minus_aj=abs(r.cand - float(aj)),
                         cand_is_whole_cell=min(abs(r.cand - float(ai)), abs(r.cand - float(aj))) <= r.tau,
                         whole_cell_is=("i" if abs(r.cand - float(ai)) <= abs(r.cand - float(aj)) else "j"),
                         v0_valid=(cv.get("V0") or {}).get("valid")))
R = pd.DataFrame(rows)
def setof(k):
    for p in ("U2_", "CP1P_", "CP1T_"):
        if k.startswith(p): return p[:-1]
    return "H" if int(k.split("_")[-1]) >= 61 else ("N" if k[0] == "A" else "?")
R["set"] = R.key.map(setof); R.to_csv(f"{D}/review_exact/v0_pairs.csv", index=False)
print(R.groupby(["set", "v0_valid"], dropna=False).agg(n=("i", "size"), keys=("key", "nunique"), whole=("cand_is_whole_cell", "sum"), max_exact=("exact", "max"),
      min_cand=("cand", "min"), max_dev=("cand_minus_ai", lambda x: np.nan)).to_string())
R["dev"] = np.minimum(R.cand_minus_ai, R.cand_minus_aj)
print("max |cand - area cella intera| fra le 'whole':", R[R.cand_is_whole_cell].dev.max(), " non-whole:", (~R.cand_is_whole_cell).sum())
print(R[~R.cand_is_whole_cell].to_string())
# V0 N copertura valida: le 20 chiavi
O1 = pd.read_csv("/home/user/2025.geo_spatialtrans/results/S1.3d/S1.3d_o1_long.csv"); v = O1[(O1.cand == "V0") & (O1.set == "N") & (O1.v0_valid == True) & (O1.silent == True)]
k20 = sorted(v.key); sub = R[R.key.isin(k20)]
print("20 chiavi:", len(k20), "; coppie violate in esse:", len(sub), "; chiavi con coppie:", sub.key.nunique(), "; tutte 'cella intera':", bool(sub.cand_is_whole_cell.all()), "; esatto max:", sub.exact.max())
print("chiavi delle 20 senza coppie violate:", sorted(set(k20) - set(sub.key)))
