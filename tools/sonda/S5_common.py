"""tools/sonda/S5_common.py - S5 (docs/sonda/REGOLA_S5.md): percorsi, pannello G_p, aree omologhe di X1 e letture comuni ai bracci XB e XV.
G_p = identificativi Ensembl del pannello Xenium di X1 (gene_panel.json) presenti fra le feature di square_002um/filtered_feature_bc_matrix.h5
di H2, nell'ordine di quelle feature (= ordine di genes.tsv delle esecuzioni di V). Aree: /mnt/micron/geo_spatialtrans/S5_sonda/explore/areas_hom_X1.csv
(quadrati da 750 um in um di array di X1; si usano x0, y0, side, tissue, k), unita' X1_<tessuto>_<k>, i = 0..9 (tumore 1..5, stroma 1..5).
Meta' tenuta da parte di X1 come REGOLA_S1: m = mediana di x sui quadratini da 8 um in tissue (3300 um); si usa solo x < m - 50 um."""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
C.SEC["X1"] = f"{C.M}/data_lung/X1/visium"
M = C.M; ROOT = f"{M}/S5_sonda"; REPO = C.REPO; RES = f"{REPO}/results/sonda/S5"; REF = f"{M}/S3_sonda/ref"
XEN = f"{M}/data_lung/X1/xenium"; XPRE = f"{XEN}/Xenium_V1_Human_Lung_Cancer_FFPE"
AREAS_X1 = f"{ROOT}/explore/areas_hom_X1.csv"
SEED_THIN = 20261130; MARGIN_UM = 240.0; L_UM = 200.0
SRC = {"V": f"{M}/S3_sonda/runs/V", "T": f"{M}/S3b_sonda/runs/T", "B": f"{M}/S3_sonda/runs/B", "P": f"{M}/S3_sonda/runs/P"}
DST = {"V": "Vp", "T": "Tp", "B": "Bp", "P": "Pp"}
UNITS = [("tumour", k) for k in range(1, 6)] + [("stroma", k) for k in range(1, 6)]

def sym2ens():
    panel = json.load(open(f"{XPRE}_gene_panel.json"))["payload"]["targets"]; out = {}
    for t in panel:
        dd = t.get("type", {}).get("data", {})
        if str(dd.get("id", "")).startswith("ENSG"): out[dd["name"]] = dd["id"]
    return out

def gp():
    """G_p nell'ordine delle feature della filtrata di H2."""
    ids, _ = C.features(f"{C.bdir('H2', '002')}/filtered_feature_bc_matrix.h5"); ens = set(sym2ens().values())
    g = [x for x in ids if x in ens]; assert len(g) == 284, len(g)
    return np.array(g)

def areas_x1():
    a = pd.read_csv(AREAS_X1); out = []
    for i, (t, k) in enumerate(UNITS):
        r = a[(a.tissue == t) & (a.k == k)].iloc[0]
        out.append(dict(unit=f"X1_{t}_{k}", tissue=t, k=k, i=i, x0=float(r.x0), y0=float(r.y0), side=float(r.side)))
    return out

def m_x1(): return C.split_m("X1")
def xmax_dev(): return m_x1() - 50.0

def write_X(path, X):
    """X csr quadratini x geni -> X.parquet (bin, gene, n) ordinato per quadratino e gene, come S3_export."""
    import pyarrow as pa, pyarrow.parquet as pq
    T = X.tocoo(); o = np.lexsort((T.col, T.row))
    pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": T.row[o].astype(np.int64), "gene": T.col[o].astype(np.int64), "n": T.data[o].astype(np.int64)}), preserve_index=False), path)
