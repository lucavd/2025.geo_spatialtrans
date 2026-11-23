"""tools/sonda/S3_smoke_check.py - verifiche del check smoke di S3 (tools/sonda/S3_smoke.sh).
  <smoke_dir> [--init]  C-3: misure della catena ridotta (V e B) uguali entro 1e-9 a results/sonda/S3/smoke_ref.json (--init: lo scrive)
  <smoke_dir> meter     metro: controlli C-1, C-2, C-4, C-5 passati e D identico (1e-9) a results/sonda/S3/D.json; results/misura.json presente"""
import os, sys, json, numpy as np
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); S = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else ""
REF = "results/sonda/S3/smoke_ref.json"; KEYS = ["M1", "M2", "M3", "M4", "M5", "M6", "n_cells", "n_typed", "n_nuc_blocks", "n_pairs", "umi_reg", "rctd_n_blocks", "rctd_n_run"]
ok = True
if mode != "meter":
    cur = {}
    for arm, u in (("SMOKE_V", "SMOKE_H2_tumour_1"), ("SMOKE_B", "SMOKE_H2_tumour_1_k0")):
        f = f"{S}/measures/{arm}/{u}.json"
        if not os.path.exists(f): print(f"FAIL chain: {f} missing"); sys.exit(1)
        r = [x for x in json.load(open(f)) if x["labels"] == "tl"][0]; cur[arm] = {k: r[k] for k in KEYS}
    print("ok   chain end to end", {a: {k: round(v, 4) for k, v in d.items() if k.startswith("M")} for a, d in cur.items()})
    if mode == "--init" or not os.path.exists(REF):
        json.dump(cur, open(REF, "w"), indent=1); print("init smoke_ref.json")
    else:
        ref = json.load(open(REF)); dev = max(abs(cur[a][k] - ref[a][k]) / max(abs(ref[a][k]), 1e-12) for a in ref for k in KEYS)
        if dev <= 1e-9: print(f"ok   C-3 chain identity, max rel dev {dev:.2e}")
        else: print(f"FAIL C-3 chain identity, max rel dev {dev:.3e}"); ok = False
else:
    d = json.load(open(f"{S}/D.json")); r = json.load(open("results/sonda/S3/D.json"))
    for c in ("C1_synthetic", "C2_colon_identity", "C4_positive", "C5_negative"):
        p = d.get(c, {}).get("pass"); print(("ok  " if p else "FAIL"), c); ok &= bool(p)
    dev = max(abs(d["by"][g]["D"][k] - r["by"][g]["D"][k]) / max(abs(r["by"][g]["D"][k]), 1e-12) for g in r["by"] for k in r["by"][g]["D"])
    print(("ok  " if dev <= 1e-9 else "FAIL"), f"meter identity vs results/sonda/S3/D.json, max rel dev {dev:.2e}"); ok &= dev <= 1e-9
    f = "results/misura.json" if os.path.exists("results/misura.json") else "results/sonda/S3/misura_provvisoria.json"
    m = json.load(open(f)); print("ok  " if f == "results/misura.json" else "FAIL", f, "A1 = D(B) =", m["valore"], "se", m["se"]); ok &= f == "results/misura.json"
sys.exit(0 if ok else 1)
