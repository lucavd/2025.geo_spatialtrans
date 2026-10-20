#!/usr/bin/env python
"""tools/R3b_posthoc.py — R3b, analisi POST HOC (scritta dopo aver visto P-2: CO nel linfonodo 0.75 contro l'atteso 0.45-0.65).
Non cambia previsioni, intervalli, verdetti ne' la regola. Domande: (1) i contatti della maschera CO sui focali di R2c vengono dai
poligoni Space Ranger rimasti (contigui) o dai nuclei StarDist? (2) sugli stessi focali, quanto concorda ogni maschera con il
giudizio di Luca (tabella 2x2)?
Uscite: results/R3b/R3b_posthoc_partners.csv, R3b_posthoc_agreement.csv.  Uso: .venv/bin/python tools/R3b_posthoc.py
"""
import importlib.util, sys
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "tools"))
def _load(name, f):
    spec = importlib.util.spec_from_file_location(name, ROOT / f); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
K = _load("r2c_cand", "tools/R2c_candidates.py"); B = _load("r3b", "tools/R3b_consensus.py"); A_ = _load("r2c_ana", "tools/R2c_analysis.py")
R3B = Path("/mnt/micron/geo_spatialtrans/R3b"); RES = ROOT / "results/R3b"; TIS = ["A4", "A6"]


def main():
    foc = pd.read_csv(ROOT / "results/R2c/R2c_focal.csv"); foc = foc[foc.focal & foc.archetype.isin(TIS)][["win_id", "archetype", "idx"]].reset_index(drop=True)
    L = A_.load_answers("Luca"); P = L[L.kind == "pair"][["win_id", "i", "j", "answer"]]
    foc["manual"] = A_.focal_s1(P, foc, "primary")
    sf = pd.read_csv(RES / "R3b_p2_focal.csv").merge(foc, on=["win_id", "archetype", "idx"], validate="m:1")
    wins = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id"); rois = K.S.rois_table()
    ag = []
    for (A, s), g in sf[sf.matched].groupby(["archetype", "seg"]):
        g = g[np.isfinite(g.manual)]; y, m = g.manual.astype(bool), g.seg_contact.astype(bool)
        ag.append(dict(archetype=A, mask=s, n=len(g), manual_touch=int(y.sum()), mask_touch=int(m.sum()), both=int((y & m).sum()),
                       mask_only=int((~y & m).sum()), manual_only=int((y & ~m).sum()), neither=int((~y & ~m).sum()),
                       manual_s1_same_nuclei=float(y.mean()), mask_s1=float(m.mean())))
    pd.DataFrame(ag).to_csv(RES / "R3b_posthoc_agreement.csv", index=False)
    rows = []
    for wid, g in sf[(sf.seg == "consensus") & sf.matched].groupby("win_id"):
        A, roi = wins.loc[wid].archetype, wins.loc[wid].roi_id
        rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
        lab = B.cp2_filter(np.load(R3B / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"], upp, K.S.read_valid(A, roi), side_um)
        pr = K.S.contact_pairs(lab, 2); mp = pd.read_parquet(R3B / "masks" / f"{A}_{roi}_consensus_map.parquet").set_index("label").src
        part = {}
        for i, j in zip(pr.label_i.astype(int), pr.label_j.astype(int)):
            part.setdefault(i, []).append(j); part.setdefault(j, []).append(i)
        for r in g.itertuples():
            p = part.get(r.label, []); ps = [mp[q] for q in p]
            rows.append(dict(win_id=wid, archetype=A, idx=r.idx, label=r.label, focal_src=mp[r.label], manual=r.manual,
                             contact=bool(p), n_partners=len(p), any_sd_partner="sd" in ps, any_sr_partner="sr" in ps,
                             contact_sd_sd=bool(mp[r.label] == "sd" and "sd" in ps)))
    pt = pd.DataFrame(rows); pt.to_csv(RES / "R3b_posthoc_partners.csv", index=False)
    pd.set_option("display.width", 200)
    print(pd.DataFrame(ag).round(3).to_string(index=False))
    print(pt.groupby("archetype").agg(n=("idx", "size"), focal_sr=("focal_src", lambda x: (x == "sr").sum()), contact=("contact", "mean"),
                                      via_sr=("any_sr_partner", "mean"), via_sd=("any_sd_partner", "mean"), sd_sd=("contact_sd_sd", "mean")).round(3).to_string())


if __name__ == "__main__":
    main()
