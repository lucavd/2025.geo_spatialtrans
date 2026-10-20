#!/usr/bin/env python
"""tools/R3b_p1p2.py — R3b, previsioni P-1 (nuclei della maschera / punti di Luca nelle finestre) e P-2 (contatti nella maschera
sui nuclei focali di R2c), per Space Ranger (riferimento), StarDist e CO, in linfonodo (A4) e cuore (A6).
Codice di R2b (tools/R2b_match.py::analyse_window: appaiamento punto-oggetto, oggetti keep con centroide nella finestra e fuori
dalle esclusioni) e di R2c (tools/R2c_refs.py: contatto = etichetta in una coppia di contact_pairs(lab, 2) sulle etichette
filtrate come cp2). Stima di P-2 come R2c (media sui focali appaiati); intervallo bootstrap a due stadi (finestre, poi focali),
B = 2000, seme 20261021 (funzione copiata da tools/R2c_posthoc.py::two_stage). Regressione G-3: con SR e SD riproduce n_obj, tp
di R2b_matching.csv, appaiamenti di R2b_points.csv e seg_contact di R2c_seg_focal.csv.
Uscite: results/R3b/R3b_p1_windows.csv, R3b_p2_focal.csv, R3b_p1p2_summary.csv, R3b_G3_p1p2.csv
Uso (root del repo): .venv/bin/python tools/R3b_p1p2.py
"""
import importlib.util, json, sys
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "tools"))
def _load(name, f):
    spec = importlib.util.spec_from_file_location(name, ROOT / f); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
MB = _load("r2b_match", "tools/R2b_match.py"); K = _load("r2c_cand", "tools/R2c_candidates.py"); B = _load("r3b", "tools/R3b_consensus.py")
R3B = Path("/mnt/micron/geo_spatialtrans/R3b"); RES = ROOT / "results/R3b"; RES.mkdir(parents=True, exist_ok=True)
SEGS = ["spaceranger", "stardist_he", "consensus"]; TIS = ["A4", "A6"]; SEED = 20261021; NB = 2000


def two_stage(win, y, B_=NB, seed=SEED):
    rng = np.random.default_rng(seed); W = list(dict.fromkeys(win)); g = [y[win == w] for w in W]; bs = []
    for _ in range(B_):
        pick = rng.choice(len(W), len(W), replace=True)
        bs.append(np.nanmean(np.concatenate([rng.choice(g[k], len(g[k]), replace=True) for k in pick])))
    return tuple(np.percentile(bs, [2.5, 97.5]))


def label_sets(A, roi):
    sets, keeps = {}, {}
    for s in SEGS:
        if s == "consensus":
            lab = np.load(R3B / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"]
            d = pd.read_parquet(R3B / "nuclei" / f"{A}_{roi}_consensus_native.parquet", columns=["label", "keep"])
        else:
            lab = MB.load_labels(A, roi, s, 1)
            d = pd.read_parquet(MB.NUC_DIR / f"{A}_{roi}_{s}_native.parquet", columns=["label", "keep"])
        sets[s] = (lab, MB.objects_table(lab)); keeps[s] = set(d.label[d.keep].astype(int))
    return sets, keeps


def main():
    wins = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id")
    wj_all = json.load(open(ROOT / "results/R2b/annotations/R2b_annotations_Luca.json"))["windows"]
    rois = K.S.rois_table()
    mat, pts, con = [], [], {}
    for wj in wj_all:
        if wj["archetype"] not in TIS:
            continue
        wrow = wins.loc[wj["win_id"]]; A, roi = wj["archetype"], wj["roi_id"]
        sets, keeps = label_sets(A, roi)
        row, mrows, prow, _ = MB.analyse_window(wj, wrow, "Luca", sets, keeps)
        mat += mrows; pts.append(prow)
        rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
        valid = K.S.read_valid(A, roi)
        for s in SEGS:
            if (A, roi, s) not in con:
                lab = K.cp2_lab(A, roi, s, rois)[0] if s != "consensus" else B.cp2_filter(sets[s][0], upp, valid, side_um)
                pr = K.S.contact_pairs(lab, 2); con[(A, roi, s)] = set(pr.label_i.astype(int)) | set(pr.label_j.astype(int))
    mat = pd.DataFrame(mat); pts = pd.concat(pts, ignore_index=True)
    mat["ratio"] = mat.n_obj / mat.n_manual
    mat.to_csv(RES / "R3b_p1_windows.csv", index=False)
    # P-2 sui focali di R2c
    foc = pd.read_csv(ROOT / "results/R2c/R2c_focal.csv"); foc = foc[foc.focal & foc.archetype.isin(TIS)]
    F = foc.merge(pts, on=["win_id", "idx"], how="left", validate="1:1")
    rows = []
    for r in F.itertuples():
        A = r.archetype; roi = wins.loc[r.win_id].roi_id
        for s in SEGS:
            m = bool(getattr(r, f"match_{s}")); l = int(getattr(r, f"label_{s}")) if m else -1
            rows.append(dict(win_id=r.win_id, archetype=A, idx=r.idx, seg=s, matched=m, label=l,
                             seg_contact=float(l in con[(A, roi, s)]) if m else np.nan))
    sf = pd.DataFrame(rows); sf.to_csv(RES / "R3b_p2_focal.csv", index=False)
    # sommario
    summ = []
    for (A, s), g in mat.groupby(["archetype", "method"]):
        q = sf[(sf.archetype == A) & (sf.seg == s)].reset_index(drop=True)
        lo, hi = two_stage(q.win_id.to_numpy(), q.seg_contact.to_numpy())
        summ.append(dict(archetype=A, mask=s, n_windows=len(g), n_manual=int(g.n_manual.sum()), n_obj=int(g.n_obj.sum()),
                         p1_ratio_mean=g.ratio.mean(), p1_ratio_sd=g.ratio.std(ddof=1), p1_se=g.ratio.std(ddof=1) / np.sqrt(len(g)),
                         p1_ratio_pooled=g.n_obj.sum() / g.n_manual.sum(), recall=g.tp_points.sum() / g.n_manual.sum(),
                         precision=g.tp.sum() / g.n_obj.sum(), p2_n_focal=len(q), p2_n_matched=int(q.matched.sum()),
                         p2_est=float(np.nanmean(q.seg_contact)), p2_lo=lo, p2_hi=hi))
    summ = pd.DataFrame(summ); summ.to_csv(RES / "R3b_p1p2_summary.csv", index=False)
    # regressione G-3
    ref_m = pd.read_csv(ROOT / "results/R2b/R2b_matching.csv"); ref_m = ref_m[(ref_m.rater == "Luca") & ref_m.archetype.isin(TIS)]
    ref_p = pd.read_csv(ROOT / "results/R2b/R2b_points.csv"); ref_p = ref_p[ref_p.rater == "Luca"]
    ref_s = pd.read_csv(ROOT / "results/R2c/R2c_seg_focal.csv"); ref_s = ref_s[ref_s.archetype.isin(TIS)]
    g3 = []
    for s in ["spaceranger", "stardist_he"]:
        a = mat[mat.method == s].set_index("win_id").sort_index(); b = ref_m[ref_m.method == s].set_index("win_id").sort_index()
        for c in ["n_manual", "n_obj", "tp", "tp_points", "n_obj_unkept", "n_obj_in_excl"]:
            g3.append(dict(check=f"R2b_matching {s} {c}", equal=bool(np.array_equal(a[c].to_numpy(), b[c].to_numpy())), n=len(a)))
        pa = pts.set_index(["win_id", "idx"]); pb = ref_p[ref_p.win_id.isin(pts.win_id)].set_index(["win_id", "idx"]).loc[pa.index]
        g3.append(dict(check=f"R2b_points {s} match/label", n=len(pa),
                       equal=bool(np.array_equal(pa[f"match_{s}"].astype(bool).to_numpy(), pb[f"match_{s}"].astype(bool).to_numpy())
                                  and np.array_equal(pa[f"label_{s}"].astype(int).to_numpy(), pb[f"label_{s}"].fillna(0).astype(int).to_numpy()))))
        x = sf[sf.seg == s].set_index(["win_id", "idx"]).sort_index(); y = ref_s[ref_s.seg == s].set_index(["win_id", "idx"]).sort_index()
        g3.append(dict(check=f"R2c_seg_focal {s} seg_contact", n=len(x),
                       equal=bool(len(x) == len(y) and np.array_equal(x.matched.to_numpy(), y.matched.astype(bool).to_numpy())
                                  and np.allclose(x.seg_contact.to_numpy(), y.seg_contact.astype(float).to_numpy(), equal_nan=True, rtol=0, atol=0))))
    g3 = pd.DataFrame(g3); g3.to_csv(RES / "R3b_G3_p1p2.csv", index=False)
    print(summ.round(4).to_string(index=False)); print(g3.to_string(index=False))


if __name__ == "__main__":
    main()
