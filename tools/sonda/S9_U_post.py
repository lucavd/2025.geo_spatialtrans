"""tools/sonda/S9_U_post.py - S9 catena U (docs/sonda/REGOLA_S9.md) dopo tools/sonda/S3_chain.sh:
1. G-B2C: per ogni unita', uscite di bin2cell di U contro la catena di oggi (V, B, P: S3_sonda/runs; T: S3b_sonda/runs/T).
   md5 identici: exp_labels_expanded.npy, blocks.csv, blocks_counts.mtx (non dipendono dal riferimento).
   b2c_label.json e cells_label.parquet contengono anche quantita' dei marcatori del riferimento (types, foreign_*, n_typed_*) e il tempo (sec):
   confronto esatto dei soli campi/colonne che non dipendono dal riferimento; md5 del file intero riportato (atteso diverso).
2. Controllo della costruzione: S3_table.py sulle misure di oggi (V, B, P di S3 + T di S3b) deve dare results/sonda/S3b/area_values.csv.
3. area_values_U.csv = S3_table.py sulle misure di U (V, T omologo, B, P); h_U.csv = tools/sonda/S3b_h.py (invariato).
4. G-MD5: md5 dei file di esportazione letti (sorgenti e copie) dopo la catena, contro exports_before.tsv.
Uso: .venv/bin/python -B tools/sonda/S9_U_post.py"""
import os, sys, json, csv, glob, shutil, hashlib, subprocess, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans"; U = f"{M}/S9_sonda/U"; REPO = os.path.expanduser("~/2025.geo_spatialtrans"); OUT = f"{REPO}/results/sonda/S9/U"
os.makedirs(OUT, exist_ok=True); os.chdir(REPO)
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
def orig(arm, unit): return f"{M}/S3b_sonda/runs/T/{unit}" if arm == "T" else f"{M}/S3_sonda/runs/{arm}/{unit}"
REFDEP = lambda k: k in ("types", "sec") or k.startswith("foreign_") or k.startswith("n_typed_")
jobs = [r for r in csv.reader(open(f"{U}/jobs/all.tsv"), delimiter="\t") if r]
rows = []
for arm, unit, d, *_ in jobs:
    o = orig(arm, unit); r = {"arm": arm, "unit": unit}
    for f in ("exp_labels_expanded.npy", "blocks.csv", "blocks_counts.mtx"):
        a, b = f"{d}/{f}", f"{o}/{f}"; r[f"md5eq_{f}"] = bool(os.path.exists(a) and md5(a) == md5(b))
    a, b = f"{d}/b2c_label.json", f"{o}/b2c_label.json"
    if os.path.exists(a):
        ja, jb = json.load(open(a)), json.load(open(b)); ka = {k for k in ja if not REFDEP(k)}; kb = {k for k in jb if not REFDEP(k)}
        r["b2c_json_refindep_eq"] = bool(ka == kb and all(ja[k] == jb[k] for k in ka)); r["b2c_json_n_fields_compared"] = len(ka)
        r["b2c_json_md5eq_whole"] = md5(a) == md5(b)
        r["b2c_json_refdep_fields_differ"] = ",".join(sorted(k for k in set(ja) | set(jb) if REFDEP(k) and ja.get(k) != jb.get(k)))
    else: r["b2c_json_refindep_eq"] = False
    a, b = f"{d}/cells_label.parquet", f"{o}/cells_label.parquet"
    if os.path.exists(a):
        ta, tb = pd.read_parquet(a), pd.read_parquet(b); ca = [c for c in ta.columns if not c.startswith("foreign_")]; cb = [c for c in tb.columns if not c.startswith("foreign_")]
        r["cells_label_refindep_eq"] = bool(ca == cb and ta[ca].equals(tb[cb])); r["cells_label_cols_compared"] = ";".join(ca); r["cells_label_md5eq_whole"] = md5(a) == md5(b)
    else: r["cells_label_refindep_eq"] = False
    r["b2c_identical"] = all(r[k] for k in ("md5eq_exp_labels_expanded.npy", "md5eq_blocks.csv", "md5eq_blocks_counts.mtx", "b2c_json_refindep_eq", "cells_label_refindep_eq"))
    rows.append(r)
G = pd.DataFrame(rows); G.to_csv(f"{OUT}/gate_b2c.csv", index=False)
# 2. controllo della costruzione della tabella
tc = f"{U}/table_check/measures"
for arm in ("V", "B", "P"):
    os.makedirs(f"{tc}/{arm}", exist_ok=True)
    for f in glob.glob(f"{M}/S3_sonda/measures/{arm}/*.json"): shutil.copy2(f, f"{tc}/{arm}/")
os.makedirs(f"{tc}/T", exist_ok=True)
for f in glob.glob(f"{M}/S3b_sonda/measures/T/*.json"): shutil.copy2(f, f"{tc}/T/")
subprocess.run(["python3", "tools/sonda/S3_table.py", tc, f"{U}/table_check/area_values_check.csv"], check=True)
ref = pd.read_csv("results/sonda/S3b/area_values.csv"); chk = pd.read_csv(f"{U}/table_check/area_values_check.csv")
table_check = {"equal_frame": bool(ref.equals(chk)), "md5_equal": md5("results/sonda/S3b/area_values.csv") == md5(f"{U}/table_check/area_values_check.csv"),
               "n_rows": [len(ref), len(chk)]}
# 3. tabella di U e margini
subprocess.run(["python3", "tools/sonda/S3_table.py", f"{U}/measures", f"{OUT}/area_values_U.csv"], check=True)
T = pd.read_csv(f"{OUT}/area_values_U.csv")
fmt = {"same_columns": list(T.columns) == list(ref.columns), "same_arm_unit_k": bool(T[["arm", "A", "unit", "run", "set", "k"]].astype(str).reset_index(drop=True).equals(ref[["arm", "A", "unit", "run", "set", "k"]].astype(str).reset_index(drop=True))),
       "counts": T.groupby(["arm", "A"]).size().astype(int).to_dict(), "missing_M": int(T[["M1", "M2", "M3", "M4", "M5", "M6"]].isna().sum().sum())}
fmt["counts"] = {f"{a}|{b}": v for (a, b), v in fmt["counts"].items()}
subprocess.run([".venv/bin/python", "-B", "tools/sonda/S3b_h.py", f"{OUT}/area_values_U.csv", f"{OUT}/h_U.csv"], check=True)
# 4. md5 dopo
B = pd.read_csv(f"{U}/md5/exports_before.tsv", sep="\t"); B["md5_src_after"] = [md5(p) for p in B.src]; B["md5_copy_after"] = [md5(p) for p in B["copy"]]
B["unchanged"] = (B.md5_src == B.md5_src_after) & (B["md5_copy"] == B["md5_copy_after"]) & (B.md5_src == B.md5_copy)
B.to_csv(f"{OUT}/md5_exports.tsv", sep="\t", index=False)
summ = {"G-B2C": {"n_units": int(len(G)), "n_identical": int(G.b2c_identical.sum()), "failing": G.loc[~G.b2c_identical, "unit"].tolist(),
                  "n_md5eq_exp": int(G["md5eq_exp_labels_expanded.npy"].sum()), "n_md5eq_blocks_csv": int(G["md5eq_blocks.csv"].sum()), "n_md5eq_blocks_mtx": int(G["md5eq_blocks_counts.mtx"].sum()),
                  "n_b2c_json_refindep_eq": int(G.b2c_json_refindep_eq.sum()), "n_cells_label_refindep_eq": int(G.cells_label_refindep_eq.sum()),
                  "n_b2c_json_md5eq_whole": int(G.get("b2c_json_md5eq_whole", pd.Series(dtype=bool)).sum()), "n_cells_label_md5eq_whole": int(G.get("cells_label_md5eq_whole", pd.Series(dtype=bool)).sum())},
        "table_check_vs_S3b": table_check, "format_vs_S3b": fmt,
        "G-MD5_exports": {"n_files": int(len(B)), "n_unchanged": int(B.unchanged.sum())}}
json.dump(summ, open(f"{OUT}/post_summary.json", "w"), indent=1, default=str); print(json.dumps(summ, default=str))
