"""tools/F3.1_PH5_aniso.py - F3.1, indagine PH-5 (dichiarata il 2026-11-18): anisotropia del campo dei tipi di RCTD e del campo degli UMI a 8 um.
Per 8 direzioni (0, 22.5, ..., 157.5 gradi) e distanze 1-3 blocchi (spostamento arrotondato sulla griglia): quota delle coppie di blocchi non scartati con lo
stesso primo tipo, e correlazione di log(1 + UMI) fra blocchi; anisotropia = massimo meno minimo su theta della media sulle 3 distanze. Per area; vero,
simulato (U6, media dei semi), compartimenti (Q6). Uso: .venv/bin/python tools/F3.1_PH5_aniso.py <main_dir> <out_dir>"""
import sys, os, glob, numpy as np, pandas as pd
main, out = sys.argv[1:3]; TH = np.arange(8) * np.pi / 8; rows = []
def offsets(th):
    return sorted({(int(round(d * np.sin(th))), int(round(d * np.cos(th)))) for d in (1, 2, 3)} - {(0, 0)})
def run(d, side):
    m = pd.read_csv(f"{d}/rctd_blocks.csv", usecols=["r8", "c8", "umi", "spot_class", "first"])
    key = pd.Series(np.arange(len(m)), index=m.r8.astype(np.int64) * 100000 + m.c8)
    ok = m.spot_class.notna().to_numpy() & (m.spot_class != "reject").to_numpy(); ft = m["first"].astype(str).to_numpy(); lu = np.log1p(m.umi.to_numpy())
    same, corr = [], []
    for th in TH:
        s_, c_ = [], []
        for dr, dc in offsets(th):
            j = key.reindex((m.r8.astype(np.int64) + dr) * 100000 + m.c8 + dc).to_numpy(); has = ~np.isnan(j); i = np.flatnonzero(has); j = j[has].astype(int)
            pr = ok[i] & ok[j]; s_.append((ft[i][pr] == ft[j][pr]).mean()); c_.append(np.corrcoef(lu[i], lu[j])[0, 1])
        same.append(np.mean(s_)); corr.append(np.mean(c_))
    same, corr = np.array(same), np.array(corr)
    b = os.path.basename(d)
    rows.append({"side": side, "run": b, "tag": b.split("_k")[0], "A": b[:2], "aniso_type": same.max() - same.min(), "theta_type": float(np.degrees(TH[same.argmax()])),
                 "aniso_umi": corr.max() - corr.min(), "theta_umi": float(np.degrees(TH[corr.argmax()])), "same_mean": same.mean()})
for side, pat in [("real", f"{main}/real/*"), ("sim", f"{main}/sim/*_k*"), ("simq6", f"{main}/simq6/*_k0")]:
    for d in sorted(glob.glob(pat)):
        if os.path.exists(f"{d}/rctd_blocks.csv"): run(d, side)
T = pd.DataFrame(rows); T.to_csv(f"{out}/F3.1_PH5_aniso.csv", index=False)
A_ = T.groupby(["side", "A", "tag"])[["aniso_type", "aniso_umi"]].mean().reset_index()
S = A_.groupby(["A", "side"]).agg(type_mean=("aniso_type", "mean"), type_se=("aniso_type", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
                                  umi_mean=("aniso_umi", "mean"), umi_se=("aniso_umi", lambda x: x.std(ddof=1) / np.sqrt(len(x))), n=("tag", "size"))
S.to_csv(f"{out}/F3.1_PH5_aniso_summary.csv"); pd.set_option("display.width", 200); print(S.round(4).to_string())
# differenza vero - simulato (U6) per tessuto, con 2 errori standard (aree indipendenti sui due lati)
D = []
for A in sorted(A_.A.unique()):
    r = A_[(A_.A == A) & (A_.side == "real")]; s = A_[(A_.A == A) & (A_.side == "sim")]
    for col in ["aniso_type", "aniso_umi"]:
        dlt = r[col].mean() - s[col].mean(); se = np.sqrt(r[col].var(ddof=1) / len(r) + s[col].var(ddof=1) / len(s))
        D.append({"A": A, "measure": col, "real_minus_sim": dlt, "se": se, "beyond_2se": bool(dlt - 2 * se > 0)})
D = pd.DataFrame(D); D.to_csv(f"{out}/F3.1_PH5_aniso_diff.csv", index=False); print(D.round(4).to_string(index=False))
