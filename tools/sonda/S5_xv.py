"""tools/sonda/S5_xv.py - S5 (docs/sonda/REGOLA_S5.md): braccio XV, Visium HD dopo Xenium su X1 (descrittivo), stesse aree di XB.
Modi:
  crop, export  tools/sonda/S3_export.py invariato, eseguito sulla sezione X1 (ritaglio dell'H&E di X1 per StarDist; quadratini da 2 um con centro
                nel quadrato dell'area, conteggi della matrice grezza, label = nucleo StarDist, in_reg come S3), con la sola aggiunta della voce
                "X1": "XVfull" al dizionario ARM (sostituzione testuale controllata, una occorrenza). Aree: S3_AREAS_X1 = areas_hom_X1.csv; guardia
                della meta' di sviluppo di S3_common attiva (m di X1). Radice S3_ROOT = /mnt/micron/geo_spatialtrans/S5_sonda.
  restrict      runs/XVfull/<unit> -> runs/XV/<unit>, X.parquet e genes.tsv ristretti a G_p come tools/sonda/S5_subset.py, con lo stesso controllo
                (results/sonda/S5/gate_sub_xv.csv, fuori dal conteggio 50/50 di G-SUB).
StarDist fra crop ed export: S3_ROOT=<S5_ROOT> tools/py_stardist.sh tools/sonda/S3_stardist.py all. Uso: .venv/bin/python -B tools/sonda/S5_xv.py <modo>"""
import os, sys, json, shutil, numpy as np, pandas as pd
os.environ.setdefault("S3_ROOT", "/mnt/micron/geo_spatialtrans/S5_sonda"); os.environ.setdefault("S3_AREAS_X1", "/mnt/micron/geo_spatialtrans/S5_sonda/explore/areas_hom_X1.csv")
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S5_common as S
mode = sys.argv[1]
if mode in ("crop", "export"):
    p = f"{S.REPO}/tools/sonda/S3_export.py"; src = open(p).read(); old = 'ARM = {"H2": "V", "H1": "T"}'
    assert src.count(old) == 1; src = src.replace(old, 'ARM = {"H2": "V", "H1": "T", "X1": "XVfull"}')
    sys.argv = [p, mode, "X1"]; exec(compile(src, p, "exec"), {"__name__": "__main__", "__file__": p})
elif mode == "restrict":
    G = S.gp(); Gs = set(G); rows = []
    for u in S.areas_x1():
        s = f"{S.ROOT}/runs/XVfull/{u['unit']}"; d = f"{S.ROOT}/runs/XV/{u['unit']}"
        genes = open(f"{s}/genes.tsv").read().split(); keep = np.array([g in Gs for g in genes]); gnew = [g for g in genes if g in Gs]; assert gnew == list(G)
        remap = -np.ones(len(genes), np.int64); remap[np.flatnonzero(keep)] = np.arange(keep.sum())
        T = pd.read_parquet(f"{s}/X.parquet"); k = keep[T.gene.to_numpy()]
        exp = pd.DataFrame({"bin": T.bin.to_numpy()[k].astype(np.int64), "gene": remap[T.gene.to_numpy()[k]], "n": T.n.to_numpy()[k].astype(np.int64)})
        if not os.path.exists(f"{d}/export.json"):
            os.makedirs(d, exist_ok=True)
            for f in ("bins.parquet", "cells.parquet"): shutil.copy(f"{s}/{f}", f"{d}/{f}")
            exp.to_parquet(f"{d}/X.parquet", index=False)
            with open(f"{d}/genes.tsv", "w") as f: f.write("\n".join(gnew) + "\n")
            json.dump(dict(json.load(open(f"{s}/export.json")), arm="XV", panel="G_p", n_genes=len(gnew), umi_total_panel=int(exp.n.sum())), open(f"{d}/export.json", "w"))
        W = pd.read_parquet(f"{d}/X.parquet"); gw = np.array(open(f"{d}/genes.tsv").read().split())
        a = pd.DataFrame({"bin": T.bin.to_numpy(), "gid": np.array(genes)[T.gene.to_numpy()], "n": T.n.to_numpy()}); a = a[a.gid.isin(Gs)].sort_values(["bin", "gid"]).reset_index(drop=True)
        b = pd.DataFrame({"bin": W.bin.to_numpy(), "gid": gw[W.gene.to_numpy()], "n": W.n.to_numpy()}).sort_values(["bin", "gid"]).reset_index(drop=True)
        same = bool(len(a) == len(b) and (a.to_numpy() == b.to_numpy()).all())
        rows.append(dict(arm="XV", unit=u["unit"], n_entries=len(b), umi_panel=int(b.n.sum()), umi_orig_all=int(T.n.sum()), frac_panel=float(b.n.sum() / max(T.n.sum(), 1)), pass_=same))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{S.RES}/gate_sub_xv.csv", index=False)
