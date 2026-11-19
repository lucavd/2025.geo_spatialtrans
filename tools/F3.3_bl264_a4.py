"""tools/F3.3_bl264_a4.py - F3.3, chiusura della discrepanza della revisione (marcatori presenti in due liste nel linfonodo, BL-264):
M-3 con i marcatori doppi esclusi (M3u) sui quarti di follicolo dentro campione e simulati, da confrontare con M3u delle aree di verifica di
500 um, per vedere se la definizione cambia lo scarto fuori meno dentro e simulato meno vero. Uso: .venv-b2c/bin/python tools/F3.3_bl264_a4.py"""
import sys, glob, json, os, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "tools"); import F3_3_common as L
M = "/mnt/micron/geo_spatialtrans"; CACHE = f"{M}/F3.3/cache"; REF = f"{M}/F3.1/main/ref"
mk = pd.read_csv(f"{REF}/markers_A4.csv"); types = sorted(mk.cell_type.unique()); dup = set(mk.gene[mk.gene.duplicated(keep=False)])
def m3(d, cache, labf, box, drop):
    bins, genes, X = L.load_area(d); exp = np.load(cache); tl = pd.read_parquet(labf); ct = pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy())
    reg = bins.in_reg.to_numpy(); nuc = bins["label"].to_numpy().astype(np.int64); xb, yb = bins.x_um.to_numpy(), bins.y_um.to_numpy()
    ids, inv = np.unique(nuc[nuc > 0], return_inverse=True); na = np.bincount(inv, minlength=len(ids)); nr = np.bincount(inv, weights=reg[nuc > 0].astype(float), minlength=len(ids))
    cs = nr >= 0.5 * na
    if box is not None:
        cx = np.bincount(inv, weights=xb[nuc > 0], minlength=len(ids)) / na; cy = np.bincount(inv, weights=yb[nuc > 0], minlength=len(ids)) / na
        cs &= (cx >= box[0]) & (cx < box[2]) & (cy >= box[1]) & (cy < box[3])
    cell = ids[cs]; cty = ct.reindex(cell).to_numpy(); gidx = pd.Series(np.arange(len(genes)), index=genes); m_ = mk[mk.gene.isin(gidx.index)]
    Mm = np.zeros((len(genes), len(types)))
    for k, ty in enumerate(types):
        g = m_.gene[(m_.cell_type == ty) & (~m_.gene.isin(dup) if drop else True)]; Mm[gidx.reindex(g).dropna().astype(int).to_numpy(), k] = 1
    sel = np.isin(exp, cell); pos = pd.Series(np.arange(len(cell)), index=cell)
    I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(cell), len(bins)))
    MU = np.asarray((I @ X) @ Mm); tot = MU.sum(axis=1); own = np.full(len(cell), np.nan)
    for k, ty in enumerate(types): own[cty == ty] = MU[cty == ty, k]
    ok = np.isin(cty, types) & (tot >= 5); return float(np.nanmean((1 - own / np.where(tot > 0, tot, np.nan))[ok]))
rows = []
Q = [(f"q{i}{j}", (i * 500.0, j * 500.0, (i + 1) * 500.0, (j + 1) * 500.0)) for j in (0, 1) for i in (0, 1)]
for side, pat in (("real", f"{M}/F3.1/main/real/A4_f[12]"), ("sim", f"{M}/F3.1/main/sim/A4_f[12]_k[01]")):
    for d in sorted(glob.glob(pat)):
        run = os.path.basename(d)
        for q, box in Q:
            rows.append({"src": "IN" if side == "real" else "SIM", "run": run, "part": q,
                         "M3": m3(d, f"{CACHE}/{side}_{run}_exp.npy", f"{CACHE}/labels/{side}_{run}_tl.parquet", box, False),
                         "M3u": m3(d, f"{CACHE}/{side}_{run}_exp.npy", f"{CACHE}/labels/{side}_{run}_tl.parquet", box, True)})
for d in sorted(glob.glob(f"{M}/F3.3/runs/A4_follicolo_*")):
    run = os.path.basename(d); j = json.load(open(f"{M}/F3.3/measures/{run}.json"))[0]
    rows.append({"src": "OUT", "run": run, "part": "whole", "M3": j["M3"], "M3u": j["M3u"]})
D = pd.DataFrame(rows); D.to_csv("results/F3.3/F3.3_bl264_a4.csv", index=False)
S = D.copy(); S["unit"] = np.where(S.src == "SIM", S.run.str[:5] + "_" + S.part, S.run + "_" + S.part); S = S.groupby(["src", "unit"])[["M3", "M3u"]].mean().reset_index()
m = S.groupby("src")[["M3", "M3u"]].mean(); print(m.round(4))
print("fuori meno dentro:", (m.loc["OUT"] - m.loc["IN"]).round(4).to_dict(), "simulato meno dentro:", (m.loc["SIM"] - m.loc["IN"]).round(4).to_dict())
