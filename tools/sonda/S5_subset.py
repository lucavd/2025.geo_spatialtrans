"""tools/sonda/S5_subset.py - S5 (docs/sonda/REGOLA_S5.md, «Bracci»): copie ristrette al pannello G_p delle esecuzioni esistenti di V (H2,
S3_sonda/runs/V), T (H1 omologo, S3b_sonda/runs/T), B (2 semi, S3_sonda/runs/B) e P (S3_sonda/runs/P). X.parquet e genes.tsv ristretti a G_p
(ordine dei geni di genes.tsv, indici rinumerati 0..283); bins.parquet (quadratini e nuclei) e cells.parquet copiati; export.json copiato con
le voci del pannello. Nessuna nuova simulazione. Cancello G-SUB: per ogni esecuzione i conteggi scritti sono l'originale ristretto a G_p
(stesse terne quadratino, gene, n e stessa somma). Uscite: <S5_ROOT>/runs/{Vp,Tp,Bp,Pp}/<unit>/, results/sonda/S5/gate_sub.csv.
Uso: .venv/bin/python -B tools/sonda/S5_subset.py"""
import os, sys, json, glob, shutil, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S5_common as S
G = S.gp(); Gs = set(G); rows = []
for arm in ("V", "T", "B", "P"):
    for s in sorted(glob.glob(f"{S.SRC[arm]}/*")):
        unit = os.path.basename(s); d = f"{S.ROOT}/runs/{S.DST[arm]}/{unit}"
        genes = open(f"{s}/genes.tsv").read().split(); keep = np.array([g in Gs for g in genes]); gnew = [g for g in genes if g in Gs]
        assert gnew == list(G), (arm, unit, "ordine di G_p diverso da genes.tsv")
        remap = -np.ones(len(genes), np.int64); remap[np.flatnonzero(keep)] = np.arange(keep.sum())
        T = pd.read_parquet(f"{s}/X.parquet"); k = keep[T.gene.to_numpy()]
        exp = pd.DataFrame({"bin": T.bin.to_numpy()[k].astype(np.int64), "gene": remap[T.gene.to_numpy()[k]], "n": T.n.to_numpy()[k].astype(np.int64)})
        if not os.path.exists(f"{d}/export.json"):
            os.makedirs(d, exist_ok=True)
            for f in ("bins.parquet", "cells.parquet"): shutil.copy(f"{s}/{f}", f"{d}/{f}")
            exp.to_parquet(f"{d}/X.parquet", index=False)
            with open(f"{d}/genes.tsv", "w") as f: f.write("\n".join(gnew) + "\n")
            info = json.load(open(f"{s}/export.json")) if os.path.exists(f"{s}/export.json") else {}
            info.update(panel="G_p", n_genes=len(gnew), source=s, umi_total_panel=int(exp.n.sum()))
            json.dump(info, open(f"{d}/export.json", "w"))
        # G-SUB: rilettura indipendente dal file scritto, confronto per identificativo Ensembl
        W = pd.read_parquet(f"{d}/X.parquet"); gw = np.array(open(f"{d}/genes.tsv").read().split())
        a = pd.DataFrame({"bin": T.bin.to_numpy(), "gid": np.array(genes)[T.gene.to_numpy()], "n": T.n.to_numpy()}); a = a[a.gid.isin(Gs)]
        b = pd.DataFrame({"bin": W.bin.to_numpy(), "gid": gw[W.gene.to_numpy()], "n": W.n.to_numpy()})
        a = a.sort_values(["bin", "gid"]).reset_index(drop=True); b = b.sort_values(["bin", "gid"]).reset_index(drop=True)
        same = bool(len(a) == len(b) and (a.bin.to_numpy() == b.bin.to_numpy()).all() and (a.gid.to_numpy() == b.gid.to_numpy()).all() and (a.n.to_numpy() == b.n.to_numpy()).all())
        same_files = all(S.C.md5(f"{s}/{f}") == S.C.md5(f"{d}/{f}") for f in ("bins.parquet", "cells.parquet"))
        rows.append(dict(arm=S.DST[arm], unit=unit, n_entries=len(b), umi_panel=int(b.n.sum()), umi_orig_restricted=int(a.n.sum()), umi_orig_all=int(T.n.sum()),
                         frac_panel=float(b.n.sum() / max(T.n.sum(), 1)), entries_identical=same, bins_cells_copied=same_files, pass_=bool(same and same_files and a.n.sum() == b.n.sum())))
        print(rows[-1], flush=True)
R = pd.DataFrame(rows); os.makedirs(S.RES, exist_ok=True); R.to_csv(f"{S.RES}/gate_sub.csv", index=False)
print("G-SUB", int(R.pass_.sum()), "/", len(R))
