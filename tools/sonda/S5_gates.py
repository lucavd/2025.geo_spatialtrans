"""tools/sonda/S5_gates.py - S5 (docs/sonda/REGOLA_S5.md, «Cancelli»): raccoglie G-XR (per area: residuo sulla scacchiera p95 <= 2 um, mediana di |s|
dei riquadri rifatti dopo la correzione <= 0,5 um, controllo positivo entro 0,3 um per asse; piu' la verifica del segno), G-SUB (50/50,
results/sonda/S5/gate_sub.csv), G-HALF (x massimo di ogni quadratino di X1 usato da correzione, XB, XV < m - 50 um), G-MD5 (md5_before = md5_after)
e G-MET / CP-XB da D.json; scrive results/sonda/S5/gates.json e la sezione "gates" di results/sonda/S5/D.json (nessun numero di D cambia).
Uso: .venv/bin/python -B tools/sonda/S5_gates.py"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S5_common as S
os.chdir(S.REPO); A = S.areas_x1(); m = S.m_x1(); xmax = m - 50.0; g = {}
xr = []
for u in A:
    r = json.load(open(f"{S.ROOT}/xen/{u['unit']}/reg.json"))
    xr.append({k: r[k] for k in ("unit", "checker_p95_um", "after_s_median_um", "wrongsign_s_median_um", "s_before_median_um", "n_valid", "n_edge", "after_n_valid", "after_n_edge",
                                 "pos_ok", "sign_ok", "gxr_checker", "gxr_after", "gxr_pass")} | {"pos_max_err_um": max(max(abs(p["rec_x"] - p["inj_x"]), abs(p["rec_y"] - p["inj_y"])) for p in r["positive"])})
X = pd.DataFrame(xr); X.to_csv("results/sonda/S5/gate_xr.csv", index=False)
g["G-XR"] = {"pass": bool(X.gxr_pass.all()), "n_checker_ok": int(X.gxr_checker.sum()), "n_after_ok": int(X.gxr_after.sum()), "n_pos_ok": int(X.pos_ok.sum()),
             "n_sign_ok": int(X.sign_ok.sum()), "checker_p95_max_um": float(X.checker_p95_um.max()), "after_s_median_max_um": float(X.after_s_median_um.max()),
             "pos_max_err_um": float(X.pos_max_err_um.max()), "wrongsign_s_median_min_um": float(X.wrongsign_s_median_um.min()), "per_area": xr}
Sb = pd.read_csv("results/sonda/S5/gate_sub.csv"); g["G-SUB"] = {"pass": bool(Sb.pass_.all() and len(Sb) == 50), "n_pass": int(Sb.pass_.sum()), "n": int(len(Sb))}
half = {}
for u in A:
    half[u["unit"]] = {"correction_vis": json.load(open(f"{S.ROOT}/xen/{u['unit']}/prep.json"))["vis_bins_xmax_um"],
                       "XB": json.load(open(f"{S.ROOT}/runs/XB/{u['unit']}/export.json"))["xmax_um"],
                       "XV": json.load(open(f"{S.ROOT}/runs/XVfull/{u['unit']}/export.json"))["xmax_um"],
                       "XV_crop": json.load(open(f"{S.ROOT}/seg/{u['unit']}/crop.json"))["xmax_um"]}
mx = max(max(v.values()) for v in half.values())
g["G-HALF"] = {"pass": bool(mx < xmax), "m_um": m, "limit_um": xmax, "max_x_um": mx, "n_exceptions": int(sum(x >= xmax for v in half.values() for x in v.values())), "per_area": half}
b, a = "results/sonda/S5/md5_before.txt", "results/sonda/S5/md5_after.txt"
if os.path.exists(a):
    B = open(b).read().splitlines(); Af = open(a).read().splitlines(); diff = sorted(set(B) ^ set(Af))
    g["G-MD5"] = {"pass": bool(B == Af and len(B) > 0), "n_files": len(B), "n_diff_lines": len(diff), "diff": diff[:20]}
else: g["G-MD5"] = {"pass": None, "note": "md5_after mancante"}
D = json.load(open("results/sonda/S5/D.json")) if os.path.exists("results/sonda/S5/D.json") else None
if D is not None:
    g["G-MET"] = {"pass": D["G_MET_pass"], "C1p": D["C1p"]["pass"], "C4": D["C4"]["pass"], "C5": D["C5"]["pass"]}; g["CP-XB"] = {"pass": D["CP_XB"]["pass"]}
else:
    g["G-MET"] = {"pass": None, "note": "D.json mancante: catena non completata"}; g["CP-XB"] = {"pass": None, "note": "D.json mancante: catena non completata"}
g["G-REV"] = {"pass": None, "note": "ricalcolo indipendente da un sotto-agente: non eseguito in questo calcolo"}
g["all_computed_pass"] = bool(all(v["pass"] for k, v in g.items() if isinstance(v, dict) and k != "G-REV" and v.get("pass") is not None))
g["not_computed"] = [k for k, v in g.items() if isinstance(v, dict) and v.get("pass") is None]
json.dump(g, open("results/sonda/S5/gates.json", "w"), indent=1, default=float)
if D is not None:
    D["gates"] = {k: ({kk: vv for kk, vv in v.items() if kk != "per_area"} if isinstance(v, dict) else v) for k, v in g.items()}
    json.dump(D, open("results/sonda/S5/D.json", "w"), indent=1, default=float)
print({k: (v["pass"] if isinstance(v, dict) else v) for k, v in g.items()})
