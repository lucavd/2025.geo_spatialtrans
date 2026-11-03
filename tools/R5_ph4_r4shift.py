"""tools/R5_ph4_r4shift.py — R5, regola R-4 (P-4 smentita nell'intestino): sensibilita' di R4 PH-1 (rapporto nucleo/corona pesato per cellula,
diviso lo stesso con la maschera spostata di 20 µm) con le maschere dei nuclei spostate dello scarto stimato della sezione (K2, mediana delle tessere).
Convenzione del modello in avanti: lambda(x) = K*T(x - delta) -> la maschera allineata ai conteggi e' mask(x - delta) (sign = +1);
sign = -1 e' il controllo (spostamento opposto). Per R6 (impronte): quota dei bin del nucleo che restano allo stesso nucleo dopo lo spostamento
e correlazione di Spearman degli UMI per nucleo (le tipizzazioni di R6 non sono ricalcolate).
Uso: .venv/bin/python tools/R5_ph4_r4shift.py <outdir> [A_roi ...] -> <outdir>/tables/R5_PH4_r4shift.csv"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT = sys.argv[1]; ONLY = sys.argv[2:]
sys.argv = ["R5_nc.py", OUT, f"{OUT}/kernels.json"]
sys.path.insert(0, "tools")
import importlib.util
spec = importlib.util.spec_from_file_location("nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(spec); spec.loader.exec_module(NC)
import R5_common as C
from R6_common import footprint_labels

K = json.load(open(f"{OUT}/kernels.json")); rois = pd.read_csv("results/R2/R2_rois_checked.csv")
rows = []
for _, Rw in rois.iterrows():
    tag = f"{Rw.archetype}_{Rw.roi_id}"
    if ONLY and tag not in ONLY: continue
    Dd = NC.roi_data(Rw); kern = K[NC.KEY[Rw.archetype]]["K2"]; upp = Dd["upp"]
    di = NC.d_image(Dd, kern); spx = np.round(di / upp).astype(int)
    s20 = int(round(20.0 / upp)); rr, cc, u, keep = Dd["rr"], Dd["cc"], Dd["u"], Dd["keep"]
    base_nuc = Dd["L"]["nuc"]
    for sign in (0, 1, -1):
        m2 = np.roll(np.roll(Dd["mask"], sign * spx[0], 0), sign * spx[1], 1)
        ms = np.roll(np.roll(m2, s20, 0), s20, 1)
        L = {"nuc": footprint_labels(m2, rr, cc, upp, NC.D_EXP, "nuc"), "exp": footprint_labels(m2, rr, cc, upp, NC.D_EXP, "exp"),
             "nuc_shift": footprint_labels(ms, rr, cc, upp, NC.D_EXP, "nuc"), "exp_shift": footprint_labels(ms, rr, cc, upp, NC.D_EXP, "exp")}
        st = C.star_fast(u, C.stat_index(L, keep))
        inn = np.isin(base_nuc, keep)
        same = float(np.mean(L["nuc"][inn] == base_nuc[inn]))
        t0 = pd.Series(u[inn]).groupby(base_nuc[inn]).sum(); sel = np.isin(L["nuc"], keep)
        t1 = pd.Series(u[sel]).groupby(L["nuc"][sel]).sum().reindex(t0.index).fillna(0)
        rows.append({"tag": tag, "archetype": Rw.archetype, "sign": sign, "d_img_row_um": float(di[0]), "d_img_col_um": float(di[1]),
                     "shift_px_row": int(sign * spx[0]), "shift_px_col": int(sign * spx[1]), "star": float(st),
                     "frac_nuc_bins_same_label": same, "spearman_nuc_umi": float(spearmanr(t0.values, t1.values)[0])})
    print(tag, [round(r["star"], 3) for r in rows[-3:]], flush=True)
os.makedirs(f"{OUT}/tables", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{OUT}/tables/R5_PH4_r4shift.csv", index=False)
print("done")
