"""tools/sonda/S4_jobs.py - S4: lista delle esecuzioni della catena (tools/sonda/S3_chain.sh) per i bracci Y, G, U, A, AE, seme 0.
Unita' senza il seme (H2_<tessuto>_<k>), cartella <unit>_k0. Uso: python3 tools/sonda/S4_jobs.py [bracci,..] [i,i,..] -> $S4_ROOT/jobs/<nome>.tsv"""
import os, sys
R = os.environ.get("S4_ROOT", "/mnt/micron/geo_spatialtrans/S4_sonda"); os.makedirs(f"{R}/jobs", exist_ok=True)
arms = (sys.argv[1] if len(sys.argv) > 1 else "Y,G,U,A,AE").split(","); ii = [int(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else list(range(10))
name = sys.argv[3] if len(sys.argv) > 3 else "all"
UNITS = [("tumour", k) for k in range(1, 6)] + [("stroma", k) for k in range(1, 6)]
with open(f"{R}/jobs/{name}.tsv", "w") as f:
    for a in arms:
        for i in ii:
            t, k = UNITS[i]; u = f"H2_{t}_{k}"; f.write("\t".join((a, u, f"{R}/runs/{a}/{u}_k0", t, "simr", a, "0")) + "\n")
print(name, len(arms) * len(ii))
