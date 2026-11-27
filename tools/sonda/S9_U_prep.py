"""tools/sonda/S9_U_prep.py - S9 catena U (docs/sonda/REGOLA_S9.md): cartelle per unita' in S9_sonda/U/runs/<arm>/<unit>/ con COPIA delle sole
esportazioni (file scritti da tools/sonda/S3_export.py e S3_sim_B.py: bins.parquet, X.parquet, genes.tsv, cells.parquet, export.json; B anche
B_check.json), mai le uscite della catena. Bracci: V, B, P = S3_sonda/runs, T = S3b_sonda/runs/T (aree omologhe di H1). md5 sorgente = copia.
Liste per tools/sonda/S3_chain.sh (formato di S3_sonda/jobs): S9_sonda/U/jobs/{pilot,all}.tsv. Uso: python3 tools/sonda/S9_U_prep.py"""
import os, csv, json, shutil, hashlib
M = "/mnt/micron/geo_spatialtrans"; U = f"{M}/S9_sonda/U"
EXPORT = ["bins.parquet", "X.parquet", "genes.tsv", "cells.parquet", "export.json"]
CHAIN = {"b2c_label.json", "b2c_label.log", "blocks.csv", "blocks_counts.mtx", "cells_label.parquet", "chain_times.tsv", "exp_labels_expanded.npy",
         "labels_transferred.parquet", "rctd_blocks.csv", "rctd.json", "rctd.log"}
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
src_jobs = [f"{M}/S3_sonda/jobs/V.tsv", f"{M}/S3b_sonda/jobs/T.tsv", f"{M}/S3_sonda/jobs/B.tsv", f"{M}/S3_sonda/jobs/P.tsv"]
rows = []
for j in src_jobs:
    for r in csv.reader(open(j), delimiter="\t"):
        if r and not r[0].startswith("#"): rows.append(r)
assert [r[0] for r in rows].count("V") == 10 and [r[0] for r in rows].count("T") == 10 and [r[0] for r in rows].count("B") == 20 and [r[0] for r in rows].count("P") == 10, len(rows)
os.makedirs(f"{U}/jobs", exist_ok=True); os.makedirs(f"{U}/md5", exist_ok=True); man = []; unexpected = []; out = []
for arm, unit, d, tissue, side, st, k in rows:
    files = EXPORT + (["B_check.json"] if arm == "B" else [])
    have = set(os.listdir(d)); extra = sorted(have - set(files) - CHAIN)
    if extra: unexpected.append((arm, unit, extra))
    dst = f"{U}/runs/{arm}/{unit}"; os.makedirs(dst, exist_ok=True)
    for f in files:
        s = f"{d}/{f}"; ms = md5(s); t = f"{dst}/{f}"
        if not os.path.exists(t): shutil.copy2(s, t)
        mt = md5(t); assert ms == mt, (s, ms, mt)
        man.append((arm, unit, s, t, ms, mt))
    out.append((arm, unit, dst, tissue, side, st, k))
with open(f"{U}/jobs/all.tsv", "w") as f:
    for r in out: f.write("\t".join(r) + "\n")
with open(f"{U}/jobs/pilot.tsv", "w") as f:
    f.write("\t".join([r for r in out if r[0] == "V" and r[1] == "H2_tumour_1"][0]) + "\n")
with open(f"{U}/md5/exports_before.tsv", "w") as f:
    f.write("arm\tunit\tsrc\tcopy\tmd5_src\tmd5_copy\n")
    for r in man: f.write("\t".join(r) + "\n")
print(json.dumps({"n_units": len(out), "n_files": len(man), "unexpected_files": unexpected,
                  "by_arm": {a: sum(1 for r in out if r[0] == a) for a in "VTBP"}}))
