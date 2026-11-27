"""tools/sonda/S9_celltypist.py - S9, catena C (docs/sonda/REGOLA_S9.md, «Catena C»): CellTypist 1.6.5 (Human_Lung_Atlas.pkl, 61 etichette)
al posto del classificatore trasferito e di RCTD. Uscite di bin2cell della catena di oggi lette in sola lettura dalle cartelle originali.
- cellule: conteggi del nucleo allargato come F3_3_common.transfer_labels, >= 20 UMI, tipo = livello 2 (HLCA ann_level_2) dell'etichetta
  con probabilita' massima; nessun majority_voting.
- blocchi da 8 um (blocks.csv, blocks_counts.mtx): >= 20 UMI; tipi = livelli 2 delle etichette con p >= 0.5; assegnato se >= 1 tipo;
  primo = livello 2 dell'etichetta con probabilita' massima.
- M-1, M-2, M-3: F3_3_common.measures (invariata, chiamata con i tipi CellTypist; riferimento e marcatori di S3_sonda/ref).
  M-4c, M-5c, M-6c: stesso codice di F3_3_common.measures per M4/M5/M6 con le chiamate CellTypist dei blocchi.
Ingresso a CellTypist: normalize_total(1e4) + log1p su tutti i geni, nomi di Space Ranger (filtered_feature_bc_matrix.h5).
Uso (.venv-enact/bin/python -B):  mapping | refctl | unit <arm> <run> | gcode <arm> <run> | table | md5 <tag> | md5cmp <tagA> <tagB>
Uscite per unita' in $S9C_OUT (predefinito /mnt/micron/geo_spatialtrans/S9_sonda/C), risultati in results/sonda/S9/C."""
import sys, os, json, time, glob, hashlib, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
REPO = os.path.expanduser("~/2025.geo_spatialtrans"); sys.path.insert(0, f"{REPO}/tools"); import F3_3_common as L
M = "/mnt/micron/geo_spatialtrans"; L.REF = f"{M}/S3_sonda/ref"; A_REF = "L"
OUT = os.environ.get("S9C_OUT", f"{M}/S9_sonda/C"); RES = f"{REPO}/results/sonda/S9/C"
MODEL = os.path.expanduser("~/.celltypist/data/models/Human_Lung_Atlas.pkl"); MODEL_MD5 = "bb198782a612415e56f6565fa9401d21"
HLCA = f"{M}/ref_lung/HLCA_core_066943a2.h5ad"
SR_H5 = {s: f"{M}/data_lung/{s}/visium/binned/binned_outputs/square_002um/filtered_feature_bc_matrix.h5" for s in ("H1", "H2")}
UMI_MIN = 20; P_THR = 0.5
RUNDIR = {"V": f"{M}/S3_sonda/runs/V", "B": f"{M}/S3_sonda/runs/B", "P": f"{M}/S3_sonda/runs/P", "T": f"{M}/S3b_sonda/runs/T"}
TABLE1 = f"{REPO}/results/sonda/S3b/area_values.csv"
META = ["arm", "genes", "A", "unit", "run", "set", "k"]
COLS = META + ["M1", "M2", "M3", "M4", "M5", "M6", "n_cells", "n_typed", "n_nuc_blocks", "n_pairs", "umi_reg", "rctd_n_blocks", "rctd_n_run",
               "rctd_n_genes", "rctd_sec"]
UNIT_FILES = ["bins.parquet", "X.parquet", "genes.tsv", "cells.parquet", "blocks.csv", "blocks_counts.mtx", "exp_labels_expanded.npy",
              "rctd_blocks.csv", "b2c_label.json", "labels_transferred.parquet"]

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()

def table1(): return pd.read_csv(TABLE1, dtype={c: str for c in META})

def sr_names(genes):
    import h5py
    with h5py.File(SR_H5["H2"], "r") as h: f = h["matrix/features"]; nm = pd.Series(f["name"][:].astype(str), index=f["id"][:].astype(str))
    out = nm.reindex(genes).to_numpy(); assert not pd.isna(out).any(); return out.astype(str)

def load_model():
    from celltypist import models
    assert md5(MODEL) == MODEL_MD5; return models.Model.load(MODEL)

def level2_map():
    t = pd.read_csv(f"{RES}/label_map_level2.csv"); return pd.Series(t.level2.to_numpy(), index=t.label.to_numpy())

def ct_predict(Xc, names, model):
    """CellTypist sui conteggi Xc (csr, n x geni dell'area): normalize_total 1e4 + log1p su tutti i geni, best match senza majority voting.
    Ritorna (etichetta di probabilita' massima, matrice di probabilita' n x 61 nell'ordine di model.cell_types, decisione, n geni comuni)."""
    import anndata as ad, scanpy as sc, celltypist
    a = ad.AnnData(sp.csr_matrix(Xc, dtype=np.float64)); a.var_names = pd.Index(names); a.var_names_make_unique()
    feat = set(model.features); ov0 = set(names) & feat; ov1 = set(a.var_names) & feat; assert ov0 == ov1, "make_unique changed overlap"
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    res = celltypist.annotate(a, model=model, mode="best match")  # 1.6.5: nessun majority voting (argomento assente)
    lab = np.asarray(model.cell_types); P = res.probability_matrix[lab].to_numpy(); Dm = res.decision_matrix[lab].to_numpy()
    best = lab[Dm.argmax(axis=1)]; assert (best == res.predicted_labels["predicted_labels"].astype(str).to_numpy()).all()
    return best, P, Dm, len(ov0)

def cell_counts(X, exp, labels):
    """stessa matrice di F3_3_common.transfer_labels: conteggi del nucleo allargato per ogni etichetta di <labels>."""
    labels = np.asarray(labels, np.int64); pos = pd.Series(np.arange(len(labels)), index=labels); sel = np.isin(exp, labels)
    I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(labels), X.shape[0]))
    C = (I @ X).tocsr(); u = np.asarray(C.sum(axis=1)).ravel(); return labels, C, u

def measures_c(m, cell_types, types):
    """M-4c, M-5c, M-6c: codice di F3_3_common.measures (M4, M5, M6) con m = chiamate dei blocchi (colonne r8, c8, nuc_label, nuc_bins,
    assigned, n_types, first). assigned sostituisce 'run & spot_class != reject', n_types == 1 sostituisce 'singlet'."""
    ok = m.assigned.to_numpy(bool); r = {}
    r["M4"] = float((m.n_types[ok] == 1).mean())
    nt = pd.Series(np.where(m.nuc_bins.to_numpy() >= 4, cell_types.reindex(m.nuc_label).to_numpy(), None), index=m.index)
    hn = ok & nt.notna().to_numpy() & nt.isin(types).to_numpy()
    r["M5"] = float((m["first"].to_numpy()[hn] == nt.to_numpy()[hn]).mean()); r["n_nuc_blocks"] = int(hn.sum())
    key = (m.r8.astype(str) + " " + m.c8.astype(str)).to_numpy(); ft = pd.Series(np.where(ok, m["first"].to_numpy(), None), index=key)
    nb = np.concatenate([ft.reindex((m.r8.astype(str) + " " + (m.c8 + 1).astype(str)).to_numpy()).to_numpy(),
                         ft.reindex(((m.r8 + 1).astype(str) + " " + m.c8.astype(str)).to_numpy()).to_numpy()])
    sf = np.concatenate([ft.to_numpy(), ft.to_numpy()]); pp = pd.notna(nb) & pd.notna(sf)
    r["M6"] = float((nb[pp] == sf[pp]).mean()); r["n_pairs"] = int(pp.sum())
    return r

def rctd_calls(d):
    """chiamate di RCTD della catena di oggi nel formato di measures_c (per G-CODE)."""
    m = pd.read_csv(f"{d}/rctd_blocks.csv"); run = m.spot_class.notna().to_numpy(); ok = run & (m.spot_class != "reject").to_numpy()
    return m.assign(assigned=ok, n_types=np.where(m.spot_class == "singlet", 1, 2))

def load_unit(arm, run):
    d = f"{RUNDIR[arm]}/{run}"; bins, genes, X = L.load_area(d); assert os.path.exists(f"{d}/exp_labels_expanded.npy")
    exp = np.load(f"{d}/exp_labels_expanded.npy"); cells = pd.read_parquet(f"{d}/cells.parquet"); return d, bins, genes, X, exp, cells

def row_of(arm, run):
    t = table1(); r = t[(t.arm == arm) & (t.run == run)]; assert len(r) == 1, (arm, run, len(r)); return r.iloc[0]

def cmd_unit(arm, run):
    o = f"{OUT}/{arm}/{run}.json"
    if os.path.exists(o) and os.path.getsize(o) > 0: print("skip", o); return
    os.makedirs(f"{OUT}/{arm}", exist_ok=True); t0 = time.time(); meta = row_of(arm, run)
    d, bins, genes, X, exp, cells = load_unit(arm, run); names = sr_names(genes); model = load_model(); l2 = level2_map()
    types = sorted(pd.read_csv(f"{L.REF}/markers_{A_REF}.csv").cell_type.unique())
    # cellule
    labels, C, u = cell_counts(X, exp, cells.label.to_numpy()); okc = u >= UMI_MIN; t1 = time.time()
    bc, Pc, Dc, ng = ct_predict(C[okc], names, model)
    ct = np.full(len(labels), None, dtype=object); ct[okc] = l2.reindex(bc).to_numpy()
    lab_best = np.full(len(labels), None, dtype=object); lab_best[okc] = bc
    cdf = pd.DataFrame({"label": labels, "umi_exp_all": u, "run": okc, "best_label": lab_best, "cell_type": ct})
    cdf.loc[okc, "best_prob"] = Pc.max(axis=1); cell_types = pd.Series(ct, index=labels)
    # blocchi
    blk = pd.read_csv(f"{d}/blocks.csv"); Xb = sio.mmread(f"{d}/blocks_counts.mtx").T.tocsr(); assert Xb.shape == (len(blk), len(genes))
    ub = np.asarray(Xb.sum(axis=1)).ravel(); assert np.array_equal(ub, blk.umi.to_numpy(float)); okb = ub >= UMI_MIN
    bb, Pb, Db, _ = ct_predict(Xb[okb], names, model); t2 = time.time()
    lab = np.asarray(model.cell_types); l2lab = l2.reindex(lab).to_numpy(); hit = Pb >= P_THR
    tset = ["|".join(sorted(set(l2lab[h]))) for h in hit]
    ntyp = np.zeros(len(blk), int); ntyp[okb] = [len(s.split("|")) if s else 0 for s in tset]
    first = np.full(len(blk), None, dtype=object); first[okb] = l2.reindex(bb).to_numpy()
    tys = np.full(len(blk), None, dtype=object); tys[okb] = tset
    bdf = blk[["block", "r8", "c8", "x", "y", "umi", "nuc_label", "nuc_bins"]].assign(run=okb, n_types=ntyp, assigned=ntyp >= 1, types=tys,
                                                                                        first=first)
    bdf.loc[okb, "best_label"] = bb; bdf.loc[okb, "best_prob"] = Pb.max(axis=1); bdf.loc[okb, "n_labels_p05"] = hit.sum(axis=1)
    # misure
    r1 = L.measures(d, A_REF, cell_types, exp, bins, genes, X); rc = measures_c(bdf, cell_types, types)
    bj = json.load(open(f"{d}/b2c_label.json"))
    row = {**{k: meta[k] for k in META}, "M1": r1["M1"], "M2": r1["M2"], "M3": r1["M3"], **{k: rc[k] for k in ("M4", "M5", "M6")},
           "n_cells": r1["n_cells"], "n_typed": r1["n_typed"], "n_nuc_blocks": rc["n_nuc_blocks"], "n_pairs": rc["n_pairs"],
           "umi_reg": float(bj["umi_reg"]), "rctd_n_blocks": int(len(blk)), "rctd_n_run": int(okb.sum()), "rctd_n_genes": int(ng),
           "rctd_sec": round(t2 - t1, 3),
           "diag": {"n_cells_ct": int(okc.sum()), "n_blocks_assigned": int((ntyp >= 1).sum()), "n_types_dist": {str(k): int(v) for k, v in
                    zip(*np.unique(ntyp[okb], return_counts=True))}, "n_ties_p05_cells": int((Pc == P_THR).sum()),
                    "n_ties_p05_blocks": int((Pb == P_THR).sum()), "n_blocks_ge_p05_strict_differs": int(((Pb >= P_THR) != (Pb > P_THR)).any(1).sum()),
                    "rctd_M4_same_code": r1["M4"], "rctd_M6_same_code": r1["M6"], "cell_type_share": pd.Series(ct[okc]).value_counts(normalize=True).round(5).to_dict(),
                    "first_share": pd.Series(first[okb]).value_counts(normalize=True).round(5).to_dict(), "sec_total": round(time.time() - t0, 1)}}
    for nm, obj in (("cells", cdf), ("blocks", bdf)):
        obj.to_parquet(f"{OUT}/{arm}/{run}_{nm}.parquet.tmp", index=False); os.replace(f"{OUT}/{arm}/{run}_{nm}.parquet.tmp", f"{OUT}/{arm}/{run}_{nm}.parquet")
    for nm, obj in (("cells_prob", Pc), ("blocks_prob", Pb)):
        np.save(f"{OUT}/{arm}/{run}_{nm}.tmp.npy", obj); os.replace(f"{OUT}/{arm}/{run}_{nm}.tmp.npy", f"{OUT}/{arm}/{run}_{nm}.npy")
    json.dump(row, open(o + ".tmp", "w"), indent=1, default=float); os.replace(o + ".tmp", o)
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "diag"}), flush=True)

def cmd_gcode(arm, run):
    """G-CODE: tipi trasferiti (F3_3_common.transfer_labels, ricalcolati senza scrivere) e chiamate RCTD al posto di CellTypist."""
    d, bins, genes, X, exp, cells = load_unit(arm, run); types = sorted(pd.read_csv(f"{L.REF}/markers_{A_REF}.csv").cell_type.unique())
    tl, _ = L.transfer_labels(X, exp, cells.label.to_numpy(), A_REF, genes); old = pd.read_parquet(f"{d}/labels_transferred.parquet")
    same_tl = bool(old.label.equals(tl.label) and (old.cell_type.fillna("NA").to_numpy() == tl.cell_type.fillna("NA").to_numpy()).all())
    cts = pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy())
    r1 = L.measures(d, A_REF, cts, exp, bins, genes, X); rc = measures_c(rctd_calls(d), cts, types); t = row_of(arm, run)
    got = {"M1": r1["M1"], "M2": r1["M2"], "M3": r1["M3"], **{k: rc[k] for k in ("M4", "M5", "M6", "n_nuc_blocks", "n_pairs")},
           "n_cells": r1["n_cells"], "n_typed": r1["n_typed"]}
    diff = {k: abs(float(got[k]) - float(t[k])) for k in got}
    out = {"arm": arm, "run": run, "same_transferred_labels": same_tl, "got": got, "table": {k: float(t[k]) for k in got}, "absdiff": diff,
           "max_absdiff_M": max(diff[k] for k in ("M1", "M2", "M3", "M4", "M5", "M6")), "pass": bool(same_tl and max(diff.values()) < 1e-9)}
    os.makedirs(f"{RES}/gcode", exist_ok=True); json.dump(out, open(f"{RES}/gcode/{arm}_{run}.json", "w"), indent=1); print(json.dumps(out))

def cmd_mapping():
    import h5py
    model = load_model()
    def cat(h, k): g = h["obs"][k]; return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))
    with h5py.File(HLCA, "r") as h: fl = cat(h, "ann_finest_level"); l2 = cat(h, "ann_level_2")
    ct = pd.crosstab(pd.Series(np.asarray(fl), name="fin"), pd.Series(np.asarray(l2), name="l2")); rows = []
    for lb in model.cell_types:
        if lb not in ct.index: rows.append({"label": lb, "level2": None, "n_cells": 0}); continue
        v = ct.loc[lb].sort_values(ascending=False); n = int(v.sum())
        rows.append({"label": lb, "level2": v.index[0], "n_cells": n, "n_majority": int(v.iloc[0]), "purity": float(v.iloc[0] / n),
                     "second": v.index[1] if v.iloc[1] > 0 else None, "n_second": int(v.iloc[1])})
    T = pd.DataFrame(rows); os.makedirs(RES, exist_ok=True); T.to_csv(f"{RES}/label_map_level2.csv", index=False)
    mk = sorted(pd.read_csv(f"{L.REF}/markers_{A_REF}.csv").cell_type.unique()); l2s = sorted(T.level2.dropna().unique())
    chk = {"n_labels": int(len(T)), "n_mapped": int(T.level2.notna().sum()), "n_hlca_cells": int(ct.to_numpy().sum()),
           "finest_not_in_model": sorted(set(ct.index) - set(model.cell_types)), "level2_types": l2s, "marker_types": mk,
           "level2_not_in_markers": sorted(set(l2s) - set(mk)), "markers_not_reached": sorted(set(mk) - set(l2s)),
           "min_purity": float(T.purity.min()), "n_purity_lt_0.9": int((T.purity < 0.9).sum()), "model_md5": md5(MODEL)}
    json.dump(chk, open(f"{RES}/label_map_check.json", "w"), indent=1); print(json.dumps(chk))

def cmd_refctl():
    """G-CT (c): CellTypist sulle cellule del riferimento di S3 (S3_sonda/ref, HLCA), accuratezza di livello 2."""
    model = load_model(); l2 = level2_map()
    Xr = sio.mmread(f"{L.REF}/ref_{A_REF}_counts.mtx").tocsr(); rg = pd.read_csv(f"{L.REF}/ref_{A_REF}_genes.tsv", header=None)[0].astype(str).to_numpy()
    rc = pd.read_csv(f"{L.REF}/ref_{A_REF}_cells.csv"); assert Xr.shape == (len(rg), len(rc)), Xr.shape
    best, P, Dm, ng = ct_predict(Xr.T.tocsr(), sr_names(rg), model); pred = l2.reindex(best).to_numpy(); ok = pred == rc.cell_type.to_numpy()
    per = rc.assign(pred=pred, ok=ok).groupby("cell_type").agg(n=("ok", "size"), accuracy=("ok", "mean")).reset_index()
    per.to_csv(f"{RES}/refctl_per_type.csv", index=False)
    pd.crosstab(rc.cell_type.rename("true"), pd.Series(pred, name="pred")).to_csv(f"{RES}/refctl_confusion.csv")
    out = {"n_cells": int(len(rc)), "n_genes_ref": int(len(rg)), "n_genes_common_model": int(ng), "accuracy": float(ok.mean()),
           "balanced_accuracy": float(per.accuracy.mean()), "min_type_accuracy": float(per.accuracy.min()), "threshold": 0.85,
           "pass": bool(ok.mean() >= 0.85)}
    json.dump(out, open(f"{RES}/refctl.json", "w"), indent=1); print(json.dumps(out)); print(per.to_string(index=False))

def cmd_table():
    t = table1(); rows = []
    for f in sorted(glob.glob(f"{OUT}/*/*.json")): r = json.load(open(f)); r.pop("diag"); rows.append(r)
    U = pd.DataFrame(rows)[COLS]
    for c_ in META: U[c_] = U[c_].astype(str)
    U = U.sort_values(["arm", "A", "run"]).reset_index(drop=True)
    assert U[META].reset_index(drop=True).equals(t[META].sort_values(["arm", "A", "run"]).reset_index(drop=True)), "rows differ from S3b table"
    U.to_csv(f"{RES}/area_values_C.csv", index=False)
    print(U.groupby(["arm", "A"]).size().to_dict(), "missing M:", int(U[["M1", "M2", "M3", "M4", "M5", "M6"]].isna().sum().sum()))

def inputs():
    t = table1(); fs = [f"{RUNDIR[a]}/{r}/{x}" for a, r in zip(t.arm, t.run) for x in UNIT_FILES]
    fs += [f"{L.REF}/{x}" for x in ("markers_L.csv", "ref_L_counts.mtx", "ref_L_genes.tsv", "ref_L_cells.csv")]
    return fs + [HLCA, MODEL, SR_H5["H1"], SR_H5["H2"], TABLE1, f"{REPO}/tools/F3_3_common.py"]

def cmd_md5(tag):
    from concurrent.futures import ThreadPoolExecutor
    fs = inputs()
    with ThreadPoolExecutor(4) as ex: h = list(ex.map(md5, fs))
    os.makedirs(RES, exist_ok=True); pd.DataFrame({"file": fs, "md5": h}).to_csv(f"{RES}/md5_{tag}.csv", index=False); print(tag, len(fs))

def cmd_md5cmp(a, b):
    A = pd.read_csv(f"{RES}/md5_{a}.csv"); B = pd.read_csv(f"{RES}/md5_{b}.csv"); same = A.equals(B)
    print(json.dumps({"n_files": int(len(A)), "identical": bool(same), "n_diff": int((A.md5 != B.md5).sum()) if len(A) == len(B) else None}))

if __name__ == "__main__":
    a = sys.argv[1:]
    {"unit": lambda: cmd_unit(a[1], a[2]), "gcode": lambda: cmd_gcode(a[1], a[2]), "mapping": cmd_mapping, "refctl": cmd_refctl,
     "table": cmd_table, "md5": lambda: cmd_md5(a[1]), "md5cmp": lambda: cmd_md5cmp(a[1], a[2])}[a[0]]()
