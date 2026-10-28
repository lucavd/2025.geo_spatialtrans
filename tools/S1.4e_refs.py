#!/usr/bin/env python
"""tools/S1.4e_refs.py - S1.4e: intervalli di rilevanza calcolati sui dati veri PRIMA di ogni simulazione nuova (docs/steps/S1.4e.md).
R-1 (per P-2) contatto manuale per coppia, per classe di distanza fra i punti di Luca (coppie giudicate 'pair' di R2c, 'unsure' escluse):
    stima e intervallo bootstrap a due stadi al 95 % (finestre con reinserimento, poi coppie con reinserimento dentro la finestra),
    B = 2000, seme 20261028. Classi: < 6 um e 6-10 um.
R-2 (per P-3) allineamento: quota dei nuclei allungati (allungamento >= 1.2) con il vicino piu' vicino a meno di 5 um che hanno
    l'asse maggiore di traverso al vicino (psi > 45 gradi); per ROI (segmentatore di forma, E5); intervallo = [minimo, massimo] dei 5 ROI.
R-3 (per P-4) eccentricita' mediana del nucleo per ROI (stesso segmentatore, keep, scala nativa); intervallo = [minimo, massimo] dei 5 ROI.
Uscita: results/S1.4e/S1.4e_refs.csv"""
import json
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/S1.4e"
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
B = 2000; SEED = 20261028
CLS = {"lt6": (0, 6), "6to10": (6, 10)}

def boot2(g, rng):
    w = g.win_id.unique(); out = np.empty(B)
    by = {k: v.touch.to_numpy() for k, v in g.groupby("win_id")}
    for b in range(B):
        ws = rng.choice(w, len(w), replace=True); y = np.concatenate([rng.choice(by[k], len(by[k]), replace=True) for k in ws])
        out[b] = y.mean()
    return out

def main():
    rng = np.random.default_rng(SEED)
    it = pd.read_csv(ROOT / "results/R2c/R2c_items_Luca.csv"); a = pd.DataFrame(json.load(open(ROOT / "results/R2c/annotations/R2c_answers_Luca.json"))["answers"])
    J = it.merge(a, on="item", validate="1:1"); J = J[(J.kind == "pair") & (J.answer != "unsure")].copy()
    J["archetype"] = J.win_id.str[:2]; J["touch"] = (J.answer == "touch").astype(float)
    rows = []
    for A in sorted(J.archetype.unique()):
        for c, (lo, hi) in CLS.items():
            g = J[(J.archetype == A) & (J.d_um >= lo) & (J.d_um < hi)]
            if len(g) < 10:
                rows.append(dict(ref="R-1", archetype=A, cls=c, n=len(g), n_win=g.win_id.nunique(), est=np.nan, lo=np.nan, hi=np.nan)); continue
            bs = boot2(g, rng)
            rows.append(dict(ref="R-1", archetype=A, cls=c, n=len(g), n_win=g.win_id.nunique(), est=float(g.touch.mean()),
                             lo=float(np.quantile(bs, .025)), hi=float(np.quantile(bs, .975))))
    E = pd.read_csv(OUT / "explore/S1.4e_E5_nuclei.csv.gz")
    for A, g in E.groupby("archetype"):
        el = g[(g.aspect >= 1.2) & (g.nn_um < 5)]
        per = el.groupby("roi_id").psi.apply(lambda s: float(np.mean(s > 45)))
        rows.append(dict(ref="R-2", archetype=A, cls="nn<5", n=len(el), n_win=len(per), est=float(per.mean()), lo=float(per.min()), hi=float(per.max())))
    rois = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")
    for A in sorted(rois.archetype.unique()):
        med = []
        for roi in rois.roi_id[rois.archetype == A]:
            d = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{roi}_{SHAPE_SEG[A]}_native.parquet", columns=["eccentricity", "scale", "keep"])
            med.append(float(d.eccentricity[(d.scale == 1) & d.keep].median()))
        rows.append(dict(ref="R-3", archetype=A, cls="all", n=len(med), n_win=len(med), est=float(np.mean(med)), lo=float(min(med)), hi=float(max(med))))
    R = pd.DataFrame(rows); R.to_csv(OUT / "S1.4e_refs.csv", index=False); print(R.round(3).to_string())

if __name__ == "__main__":
    main()
