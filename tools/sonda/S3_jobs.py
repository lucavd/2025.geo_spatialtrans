"""tools/sonda/S3_jobs.py - S3: liste delle esecuzioni per tools/sonda/S3_chain.sh: V (H2), T (H1), P (H2 permutato), B (H2, 2 semi).
Uso: python3 tools/sonda/S3_jobs.py  ->  $S3_ROOT/jobs/{V,T,P,B}.tsv"""
import os, csv
R = os.environ.get("S3_ROOT", "/mnt/micron/geo_spatialtrans/S3_sonda"); os.makedirs(f"{R}/jobs", exist_ok=True)
A = {s: list(csv.DictReader(open(os.path.expanduser(f"~/2025.geo_spatialtrans/results/sonda/S3/areas_{s}.csv")))) for s in ("H2", "H1")}
def rows(arm, sec, side="real", ks=(None,)):
    out = []
    for t in ("tumour", "stroma"):
        for a in sorted([x for x in A[sec] if x["tissue"] == t], key=lambda x: int(x["k"])):
            u = f"{sec}_{t}_{a['k']}"
            for k in ks:
                uu = u if k is None else f"{u}_k{k}"
                out.append((arm, uu, f"{R}/runs/{arm}/{uu}", t, side, {"V": "H2", "T": "H1", "P": "P", "B": "B"}[arm], "-" if k is None else str(k)))
    return out
for arm, sec, side, ks in (("V", "H2", "real", (None,)), ("T", "H1", "real", (None,)), ("P", "H2", "real", (None,)), ("B", "H2", "simr", (0, 1))):
    with open(f"{R}/jobs/{arm}.tsv", "w") as f:
        for r in rows(arm, sec, side, ks): f.write("\t".join(r) + "\n")
    print(arm, len(rows(arm, sec, side, ks)))
