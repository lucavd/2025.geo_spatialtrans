"""tools/sonda/S3_Bparams.py - S3: parametri del braccio banale B tarati su H2 (docs/sonda/REGOLA_S3.md, «Braccio banale B»).
Per tessuto, sommando numeratori e denominatori sulle 5 aree di V dopo la catena (bins, X, nuclei allargati di bin2cell, etichette trasferite):
cellule = nuclei con almeno meta' dei quadratini in in_reg (F3.3); lambda = cellule / (quadratini in_reg * 4 um^2); n_nuc = quadratini medi
per nucleo; composizione = quote dei tipi trasferiti non nulli; profilo del tipo = somma dei conteggi dei nuclei allargati delle sue cellule,
normalizzata; rho = UMI nei quadratini in_reg / quadratini in_reg. Uscite: $S3_ROOT/B/params_<tissue>.npz e .json, results/sonda/S3/B_params.json.
Uso: .venv/bin/python -B tools/sonda/S3_Bparams.py"""
import sys, os, json, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
ROOT = os.environ.get("S3_ROOT", C.ROOT); os.makedirs(f"{ROOT}/B", exist_ok=True); allp = {}
for tissue in ("tumour", "stroma"):
    ncell = nreg = nnuc = umi = 0.0; comp = {}; prof = {}; genes = None; per_area = []
    for k in range(1, 6):
        d = f"{ROOT}/runs/V/H2_{tissue}_{k}"
        bins = pd.read_parquet(f"{d}/bins.parquet"); g = open(f"{d}/genes.tsv").read().split(); genes = g if genes is None else genes
        assert g == genes; T = pd.read_parquet(f"{d}/X.parquet")
        X = sp.csr_matrix((T.n.to_numpy().astype(float), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(g)))
        exp = np.load(f"{d}/exp_labels_expanded.npy"); tl = pd.read_parquet(f"{d}/labels_transferred.parquet").set_index("label").cell_type
        reg = bins.in_reg.to_numpy(); nuc = bins.label.to_numpy().astype(np.int64)
        ids, inv = np.unique(nuc[nuc > 0], return_inverse=True)
        nb_all = np.bincount(inv, minlength=len(ids)); nb_reg = np.bincount(inv, weights=reg[nuc > 0].astype(float), minlength=len(ids))
        cs = ids[nb_reg >= 0.5 * nb_all]; nbc = nb_all[nb_reg >= 0.5 * nb_all]
        ty = tl.reindex(cs).to_numpy(); u = float(np.asarray(X[reg].sum()))
        ncell += len(cs); nreg += reg.sum(); nnuc += nbc.sum(); umi += u
        per_area.append({"unit": f"H2_{tissue}_{k}", "cells": int(len(cs)), "n_reg": int(reg.sum()), "lambda_mm2": len(cs) / (reg.sum() * 4) * 1e6,
                         "n_nuc_mean": float(nbc.mean()), "rho": u / reg.sum(), "typed": int(pd.notna(ty).sum())})
        pos = pd.Series(np.arange(len(cs)), index=cs); sel = np.isin(exp, cs)
        I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(cs), len(bins)))
        Cc = (I @ X).tocsr()
        for t in pd.unique(ty[pd.notna(ty)]):
            w = np.flatnonzero(ty == t); comp[t] = comp.get(t, 0) + len(w)
            prof[t] = prof.get(t, 0) + np.asarray(Cc[w].sum(axis=0)).ravel()
    types = sorted(comp); n_t = np.array([comp[t] for t in types], float); P = np.vstack([prof[t] for t in types]); P = P / P.sum(axis=1, keepdims=True)
    par = {"tissue": tissue, "lambda_um2": ncell / (nreg * 4.0), "lambda_mm2": ncell / (nreg * 4.0) * 1e6, "n_nuc_mean": nnuc / ncell,
           "r_um": float(np.sqrt(4.0 * (nnuc / ncell) / np.pi)), "rho_per_bin": umi / nreg, "types": types,
           "composition": dict(zip(types, (n_t / n_t.sum()).tolist())), "n_cells": int(ncell), "n_typed": int(n_t.sum()), "n_reg": int(nreg), "per_area": per_area}
    np.savez_compressed(f"{ROOT}/B/params_{tissue}.npz", profiles=P, types=np.array(types), genes=np.array(genes))
    json.dump(par, open(f"{ROOT}/B/params_{tissue}.json", "w"), indent=1); allp[tissue] = par
    print(tissue, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in par.items() if k not in ("per_area",)}, flush=True)
os.makedirs(f"{C.REPO}/results/sonda/S3", exist_ok=True); json.dump(allp, open(f"{C.REPO}/results/sonda/S3/B_params.json", "w"), indent=1)
