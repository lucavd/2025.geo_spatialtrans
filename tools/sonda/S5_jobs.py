"""tools/sonda/S5_jobs.py - S5 (docs/sonda/REGOLA_S5.md, «Calcolo»): liste per tools/sonda/S3_chain.sh, 70 esecuzioni: Vp 10, Tp 10, Bp 20 (2 semi),
Pp 10, XB 10, XV 10. Colonne: arm unit run_dir tissue side set k; arm = nome del braccio nella tabella (V, T, B, P, XB, XV; genes = panel), cartelle
<S5_ROOT>/runs/{Vp,Tp,Bp,Pp,XB,XV}. pilot = un'unita' per braccio (tumore 1, B seme 0). Uso: python3 tools/sonda/S5_jobs.py -> <S5_ROOT>/jobs/{all,pilot}.tsv"""
import os
R = "/mnt/micron/geo_spatialtrans/S5_sonda"; os.makedirs(f"{R}/jobs", exist_ok=True)
UNITS = [("tumour", k) for k in range(1, 6)] + [("stroma", k) for k in range(1, 6)]
rows = []
for t, k in UNITS:   # B per primo: RCTD piu' lento in S3
    for s in (0, 1): rows.append(("B", f"H2_{t}_{k}_k{s}", f"{R}/runs/Bp/H2_{t}_{k}_k{s}", t, "simr", "B", str(s)))
for arm, d, sec, st in (("V", "Vp", "H2", "H2"), ("T", "Tp", "H1", "H1"), ("P", "Pp", "H2", "P"), ("XB", "XB", "X1", "X1"), ("XV", "XV", "X1", "X1")):
    for t, k in UNITS: rows.append((arm, f"{sec}_{t}_{k}", f"{R}/runs/{d}/{sec}_{t}_{k}", t, "real", st, "-"))
pilot = [r for r in rows if r[1].split("_")[1:3] == ["tumour", "1"] and r[6] in ("-", "0")]
for name, rr in (("all", rows), ("pilot", pilot)):
    with open(f"{R}/jobs/{name}.tsv", "w") as f:
        for r in rr: f.write("\t".join(r) + "\n")
    print(name, len(rr))
