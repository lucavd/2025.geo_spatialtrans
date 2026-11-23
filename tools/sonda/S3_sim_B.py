"""tools/sonda/S3_sim_B.py - S3: braccio banale B (docs/sonda/REGOLA_S3.md, «Braccio banale B»), sulla griglia da 2 um delle aree di V.
Per area i (0..9: tumore 1..5, stroma 1..5) e seme k: default_rng([20261125, i, k]). N ~ Poisson(lambda * 4 * quadratini in tissue); centri =
quadratini in tissue estratti con ripetizione + scarto uniforme di +-1 um; territorio = Voronoi discreto (quadratino in tissue al centro piu'
vicino, cKDTree); cellule senza quadratini tolte; nucleo = quadratini del territorio a distanza <= r dal centro; tipo dalla composizione;
n ~ Poisson(rho * quadratini del territorio), geni ~ Multinomiale(n, profilo del tipo), ogni trascritto in un quadratino del territorio a caso.
Fuori tessuto zero; nessuno spostamento, striscia o fondo. Uscita: $S3_ROOT/runs/B/<unit>_k<k>/ nel formato di F3.1 e B_check.json.
Uso: .venv/bin/python -B tools/sonda/S3_sim_B.py <W> [i,i,..] [k,k] [smoke]"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow as pa, pyarrow.parquet as pq
from scipy.spatial import cKDTree
from multiprocessing import Pool
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
ROOT = os.environ.get("S3_ROOT", C.ROOT); PAR = os.environ.get("S3_BPAR", f"{C.ROOT}/B")
UNITS = [("tumour", k) for k in range(1, 6)] + [("stroma", k) for k in range(1, 6)]

def sim(i, k, src=None, dst=None):
    tissue, kk = UNITS[i]; src = src or f"{ROOT}/runs/V/H2_{tissue}_{kk}"; dst = dst or f"{ROOT}/runs/B/H2_{tissue}_{kk}_k{k}"
    if os.path.exists(f"{dst}/export.json"): return f"present {dst}"
    t0 = time.time(); os.makedirs(dst, exist_ok=True)
    par = json.load(open(f"{PAR}/params_{tissue}.json")); z = np.load(f"{PAR}/params_{tissue}.npz"); P = z["profiles"]; types = list(z["types"])
    genes = open(f"{src}/genes.tsv").read().split(); assert genes == list(z["genes"])
    bins = pd.read_parquet(f"{src}/bins.parquet"); x, y = bins.x_um.to_numpy(), bins.y_um.to_numpy()
    tb = np.flatnonzero(bins.in_tissue.to_numpy()); rng = np.random.default_rng([C.SEED_B, i, k])
    N = int(rng.poisson(par["lambda_um2"] * 4.0 * len(tb)))
    pick = rng.integers(0, len(tb), N); cx = x[tb[pick]] + rng.uniform(-1, 1, N); cy = y[tb[pick]] + rng.uniform(-1, 1, N)
    dist, own = cKDTree(np.c_[cx, cy]).query(np.c_[x[tb], y[tb]])
    keep = np.unique(own); newid = np.zeros(N, np.int64); newid[keep] = np.arange(1, len(keep) + 1); n = len(keep)
    cell = newid[own]                                        # cellula (1..n) di ogni quadratino in tissue
    label = np.zeros(len(bins), np.int64); label[tb[dist <= par["r_um"]]] = cell[dist <= par["r_um"]]
    ctype = rng.choice(len(types), size=n, p=np.array([par["composition"][t] for t in types]))
    size = np.bincount(cell, minlength=n + 1)[1:]
    ncnt = rng.poisson(par["rho_per_bin"] * size)
    order = np.argsort(cell, kind="stable"); start = np.r_[0, np.cumsum(size)[:-1]]; tb_sorted = tb[order]
    rows, cols = [], []
    for t in range(len(types)):
        w = np.flatnonzero(ctype == t)
        for c0 in range(0, len(w), 200):
            ww = w[c0:c0 + 200]; G = rng.multinomial(ncnt[ww], P[t])          # celle x geni
            ci, gi = np.nonzero(G); m = G[ci, gi]
            cc = np.repeat(ww[ci], m); gg = np.repeat(gi, m)
            off = (rng.random(len(cc)) * size[cc]).astype(np.int64)
            rows.append(tb_sorted[start[cc] + off]); cols.append(gg)
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    X = sp.csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(len(bins), len(genes))); X.sum_duplicates(); T = X.tocoo()
    out = bins.copy(); out["label"] = label
    out.to_parquet(f"{dst}/bins.parquet", index=False)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": T.row.astype(np.int64), "gene": T.col.astype(np.int64), "n": T.data.astype(np.int64)}), preserve_index=False), f"{dst}/X.parquet")
    with open(f"{dst}/genes.tsv", "w") as f: f.write("\n".join(genes) + "\n")
    pd.DataFrame({"label": np.arange(1, n + 1), "cell_type": np.array(types, dtype=object)[ctype]}).to_parquet(f"{dst}/cells.parquet", index=False)
    reg = bins.in_reg.to_numpy(); umi = np.asarray(X.sum(axis=1)).ravel()
    nb = np.bincount(label[label > 0], minlength=n + 1)[1:]
    chk = {"unit": os.path.basename(dst), "i": i, "k": k, "seed": [C.SEED_B, i, k], "src": src, "n_centres": N, "n_cells": int(n), "n_dropped": int(N - n),
           "lambda_target_mm2": par["lambda_mm2"], "lambda_realised_mm2": n / (len(tb) * 4.0) * 1e6, "n_nuc_target": par["n_nuc_mean"],
           "n_nuc_realised": float(nb[nb > 0].mean()), "n_cells_without_nucleus": int((nb == 0).sum()), "rho_target": par["rho_per_bin"],
           "rho_realised_reg": float(umi[reg].sum() / reg.sum()), "umi_total": int(umi.sum()), "umi_out_of_tissue": int(umi[~bins.in_tissue.to_numpy()].sum()),
           "composition_realised": {types[t]: float((ctype == t).mean()) for t in range(len(types))}, "sec": round(time.time() - t0, 1)}
    json.dump(chk, open(f"{dst}/B_check.json", "w"), indent=1)
    json.dump(dict(chk, arm="B", tissue=tissue, side_um=float(json.load(open(f"{src}/export.json"))["side_um"]), n_bins=int(len(bins)), n_reg=int(reg.sum())), open(f"{dst}/export.json", "w"))
    return f"{chk['unit']} cells {n} lambda {chk['lambda_realised_mm2']:.1f}/{chk['lambda_target_mm2']:.1f} rho {chk['rho_realised_reg']:.3f}/{chk['rho_target']:.3f} nnuc {chk['n_nuc_realised']:.2f}/{chk['n_nuc_target']:.2f} {chk['sec']} s"

if __name__ == "__main__":
    W = int(sys.argv[1]); I = [int(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 and sys.argv[2] != "all" else list(range(10))
    K = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0, 1]
    if len(sys.argv) > 4 and sys.argv[4] == "smoke":
        print(sim(0, 0, src=f"{ROOT}/runs/SMOKE_V/SMOKE_H2_tumour_1", dst=f"{ROOT}/runs/SMOKE_B/SMOKE_H2_tumour_1_k0")); sys.exit(0)
    jobs = [(i, k) for i in I for k in K]
    with Pool(W) as p:
        for r in p.starmap(sim, jobs): print(r, flush=True)
