"""tools/sonda/S3_common.py - S3: percorsi e letture comuni (docs/sonda/REGOLA_S3.md).
Metà tenuta da parte come REGOLA_S1: m = mediana di x = (array_col + 0.5) * 8 um sui quadratini da 8 um in tissue; sviluppo x < m - 50 um.
read_cols() legge da una matrice 10x (h5, CSC geni x quadratini) le sole colonne richieste, per intervalli di indptr: l'espressione degli
altri quadratini non si carica. guard() verifica che ogni quadratino usato stia nella metà di sviluppo e restituisce il x massimo."""
import os, glob, json, hashlib, numpy as np, pandas as pd, pyarrow.parquet as pq, h5py, scipy.sparse as sp
M = "/mnt/micron/geo_spatialtrans"; ROOT = f"{M}/S3_sonda"; REF = f"{ROOT}/ref"; REPO = os.path.expanduser("~/2025.geo_spatialtrans")
SEC = {s: f"{M}/data_lung/{s}/visium" for s in ("H1", "H2")}
SEED_B, SEED_P, SEED_REF = 20261125, 20261126, 20261123
FRAME_UM = 25.0

def bdir(sec, size): return f"{SEC[sec]}/binned/binned_outputs/square_{size}um"
def image(sec): return glob.glob(f"{SEC[sec]}/*_tissue_image.btf")[0]
def positions(sec, size="002"): return pq.read_table(f"{bdir(sec, size)}/spatial/tissue_positions.parquet").to_pandas()
def upp(sec): return float(json.load(open(f"{bdir(sec, '008')}/spatial/scalefactors_json.json"))["microns_per_pixel"])

def split_m(sec):
    tp8 = positions(sec, "008"); x8 = (tp8.array_col.to_numpy() + 0.5) * 8.0
    return float(np.median(x8[tp8.in_tissue.to_numpy() == 1]))

def guard(x_um_abs, m):
    xm = float(np.max(x_um_abs)); assert xm < m - 50.0, f"holdout guard: x {xm} >= m - 50 = {m - 50.0}"
    return xm

def features(h5path):
    with h5py.File(h5path, "r") as h:
        f = h["matrix/features"]; return f["id"][:].astype(str), f["name"][:].astype(str)

def barcodes(h5path):
    with h5py.File(h5path, "r") as h: return h["matrix/barcodes"][:].astype(str)

def read_cols(h5path, cols):
    """colonne (quadratini) richieste, nell'ordine dato; restituisce csr quadratini x geni (int64) e il numero di intervalli letti."""
    cols = np.asarray(cols, np.int64); order = np.argsort(cols, kind="stable"); sc = cols[order]
    with h5py.File(h5path, "r") as h:
        g = h["matrix"]; ng = int(g["shape"][0]); ip = g["indptr"]; D = g["data"]; I = g["indices"]
        brk = np.flatnonzero(np.diff(sc) != 1) + 1; starts = np.r_[0, brk]; ends = np.r_[brk, len(sc)]
        data, ind, ptr = [], [], [0]
        for a, b in zip(starts, ends):
            c0, c1 = int(sc[a]), int(sc[b - 1]) + 1; p = ip[c0:c1 + 1]
            data.append(D[p[0]:p[-1]]); ind.append(I[p[0]:p[-1]]); ptr.append(np.diff(p))
        nnz = np.concatenate(ptr[1:]) if len(ptr) > 1 else np.zeros(0, np.int64)
        X = sp.csr_matrix((np.concatenate(data).astype(np.int64), np.concatenate(ind), np.r_[0, np.cumsum(nnz)]), shape=(len(sc), ng))
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    return X[inv], len(starts)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()

def read_region(path, r0, r1, c0, c1):
    """regione [r0:r1, c0:c1] del livello 0 dell'H&E senza caricare l'immagine (come tools/R1_verify.py)."""
    import tifffile, zarr
    store = tifffile.imread(path, aszarr=True)
    try:
        z = zarr.open(store, mode="r"); z = z["0"] if isinstance(z, zarr.Group) else z
        with tifffile.TiffFile(path) as tf: axes = tf.series[0].axes
        if axes.startswith("YX"): arr = np.asarray(z[r0:r1, c0:c1])
        else: arr = np.moveaxis(np.asarray(z[:, r0:r1, c0:c1]), 0, -1)
    finally:
        store.close()
    return arr[..., :3]

def areas(sec): return pd.read_csv(f"{REPO}/results/sonda/S3/areas_{sec}.csv")
def unit_name(sec, tissue, k): return f"{sec}_{tissue}_{k}"
