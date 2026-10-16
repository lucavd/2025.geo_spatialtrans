#!/usr/bin/env python
"""tools/R2c_refs.py — R2c: valori di confronto indipendenti dalle risposte manuali.
(1) S-1 dei segmentatori sui nuclei focali appaiati: contatto = etichetta con n_c2 > 0 nello stadio cp2 di S1.4c (cellpose_rgb,
    stardist_he; identico all'operatore di R2c per C-3) e, in A4/A6, Space Ranger 4 dalle coppie a contatto di R2c (descrittivo).
(2) S-1 di K-2 (S1.4c, sim_out, semi s01..s10, ROI di test del proprio fold): nuclei simulati con centroide nella finestra a
    distanza >= d_max dal bordo e dalle esclusioni (stessa regola dei focali); S-1 = frazione con n_c2 > 0.
Output: results/R2c/R2c_seg_focal.csv, results/R2c/R2c_k2_windows.csv
Uso: .venv/bin/python tools/R2c_refs.py
"""
import importlib.util, json
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
RES = ROOT / "results/R2c"; R2B = ROOT / "results/R2b"
SIM = Path("/mnt/micron/geo_spatialtrans/S1.4c/sim_out"); CP2 = Path("/mnt/micron/geo_spatialtrans/S1.4c/cp2")
SEEDS = [f"s{k:02d}" for k in range(1, 11)]


def main():
    wins = pd.read_csv(R2B / "R2b_windows.csv"); rois = K.S.rois_table()
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"]
    foc = pd.read_csv(RES / "R2c_focal.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    F = foc[foc.focal].merge(pts, on=["win_id", "idx"], how="left", validate="1:1")
    rows = []; sr_cache = {}
    for wr in wins.itertuples():
        f = F[F.win_id == wr.win_id]
        for seg in ["cellpose_rgb", "stardist_he"] + (["spaceranger"] if wr.archetype in K.SR_ARCH else []):
            if seg == "spaceranger":
                key = (wr.archetype, wr.roi_id)
                if key not in sr_cache:
                    lab, _ = K.cp2_lab(wr.archetype, wr.roi_id, seg, rois); pr = K.S.contact_pairs(lab, 2)
                    sr_cache[key] = set(pr.label_i.astype(int)) | set(pr.label_j.astype(int))
                con = sr_cache[key]; ev = None
            else:
                cp = pd.read_parquet(CP2 / f"{wr.archetype}_{wr.roi_id}_{seg}_cp2.parquet").set_index("label")
                con = set(cp.index[cp.n_c2 > 0].astype(int)); ev = cp.ev
            for r in f.itertuples():
                m = str(getattr(r, f"match_{seg}")) == "True"
                l = int(getattr(r, f"label_{seg}")) if m else -1
                rows.append(dict(win_id=r.win_id, archetype=wr.archetype, idx=r.idx, seg=seg, matched=m, label=l,
                                 seg_contact=(l in con) if m else np.nan,
                                 ev=(bool(ev.get(l, False)) if (m and ev is not None) else np.nan)))
    sf = pd.DataFrame(rows); sf.to_csv(RES / "R2c_seg_focal.csv", index=False)
    # ---- K-2
    kr = []
    for wr in wins.itertuples():
        d = float(dm.d_max_um[wr.archetype]); upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um)
        s = int(wr.side_px) * upp
        ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
        fs = sorted(SIM.glob(f"{wr.archetype}_{wr.roi_id}_fold*_s??_K2.parquet"))
        assert len(fs) == 10 and [p.name.split("_")[3] for p in fs] == SEEDS, (wr.win_id, [p.name for p in fs])
        for p in fs:
            t = pd.read_parquet(p)
            de = np.asarray(np.minimum.reduce([t.cx - x0, x0 + s - t.cx, t.cy - y0, y0 + s - t.cy]), float)
            keep = de >= d
            if u is not None:
                keep &= shapely.distance(shapely.points(t[["cx", "cy"]].to_numpy()), u) >= d
            kr.append(dict(win_id=wr.win_id, archetype=wr.archetype, seed=p.name.split("_")[3], n=int(keep.sum()),
                           n_contact=int((t.n_c2[keep] > 0).sum())))
    kr = pd.DataFrame(kr); kr["s1"] = kr.n_contact / kr.n
    kr.to_csv(RES / "R2c_k2_windows.csv", index=False)
    a = kr.groupby(["archetype", "seed"])[["n", "n_contact"]].sum(); a["s1"] = a.n_contact / a.n
    print(a.groupby("archetype").s1.agg(["median", "min", "max"]).round(3).to_string())
    g = sf[sf.matched].groupby(["archetype", "seg"]).agg(n_focal_matched=("idx", "size"), s1=("seg_contact", "mean"))
    print(g.round(3).to_string())


if __name__ == "__main__":
    main()
