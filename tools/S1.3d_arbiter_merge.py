"""tools/S1.3d_arbiter_merge.py — S1.3d: riunisce i blocchi dell'arbitro (S) nei file standard <key>_o1*.parquet/json e <key>_o2_s*.parquet.
Somme in virgola mobile dei valori esatti per blocco (errore ~1e-9 µm² << tau_U ~ 10 µm²; dichiarato in dev_notes).
Uso: .venv/bin/python tools/S1.3d_arbiter_merge.py <truthdir> <key> <K>"""
import sys, os, json, glob
import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
truth, key, K = sys.argv[1], sys.argv[2], int(sys.argv[3])
def cat(suffix):
    fs = sorted(glob.glob(os.path.join(truth, f"{key}_chunk*of{K:03d}_{suffix}.parquet"))); assert len(fs) == K, (suffix, len(fs))
    return pa.concat_tables([pq.read_table(f) for f in fs])
out = {"key": key}
if os.path.exists(os.path.join(truth, key + "_o1_areg.json")):
    P = cat("o1_pairs").sort_by("q"); C = cat("o1_cells").sort_by("cell")
    pq.write_table(P.select(["i", "j", "a", "zero"]), os.path.join(truth, key + "_o1_pairs.parquet")); pq.write_table(C, os.path.join(truth, key + "_o1_cells.parquet"))
    A = json.load(open(os.path.join(truth, key + "_o1_areg.json")))["A_reg"]
    S = float(pc.sum(C.column("a")).as_py()); Sov = float(pc.sum(P.column("a")).as_py()); Iin = float(pc.sum(C.column("inside")).as_py())
    U = S - Sov; I = Iin - Sov
    out["o1"] = dict(sum_a=S, A_reg=A, sum_ov=Sov, union=U, inside=I, symd=(U - I) + (A - I), approx_bound=Sov, n_local=None, t_s=None, merged_from=K)
    json.dump(out, open(os.path.join(truth, key + "_o1.json"), "w"))
for s in (2, 8):
    fs = glob.glob(os.path.join(truth, f"{key}_chunk*of{K:03d}_o2_s{s}.parquet"))
    if fs: T = cat(f"o2_s{s}").sort_by([("i", "ascending"), ("ix", "ascending"), ("iy", "ascending")]); pq.write_table(T, os.path.join(truth, f"{key}_o2_s{s}.parquet")); out[f"n{s}"] = T.num_rows
if "n2" in out: json.dump({"key": key, "o2": {"n2": out["n2"], "n8": out.get("n8"), "agg_mismatch": None, "merged_from": K}}, open(os.path.join(truth, key + "_o2.json"), "w"))
print(json.dumps(out))
