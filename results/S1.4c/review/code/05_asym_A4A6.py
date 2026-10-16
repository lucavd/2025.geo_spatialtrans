#!/usr/bin/env python
"""review/code/05_asym_A4A6.py — effetto quantitativo delle asimmetrie reale/simulato in A4, A6 (Space Ranger + StarDist).

Fatti: i simulati stanno su tutti i generatori Space Ranger (n_SR); la maschera reale di forma e' StarDist (n_SD < n_SR);
B-8 reale usa a_n (StarDist appaiato) o a_own (area SR) per i non appaiati; B-7 reale usa il generatore SR e il Voronoi SR.
Alternative (ciascuna applicata in modo coerente a reali, metro inter-ROI e, dove serve, simulati):
  B8-R1  raggi dei non appaiati da a_own * mediana(a_n / a_own) degli appaiati (area SR ricalibrata a StarDist)
  B8-R2  insieme dei vicini = sole celle appaiate, nei reali e nei simulati (simulati: nuclei delle celle appaiate)
  B8-R3  reali = operatore sulla sola maschera StarDist (centroidi e aree StarDist, tutti i nuclei StarDist); simulati invariati
  B7-R2  Voronoi delle sole celle appaiate, nei reali (generatori SR) e nei simulati (centroidi raster)
  B7-R3  reali = centroidi StarDist, Voronoi di tutti i centroidi StarDist, area StarDist; simulati invariati
  S-RM   S-1/S-2/S-3 con le sole celle appaiate nei reali (etichette StarDist appaiate) e nei simulati (ri-rasterizzati)
Uso (root del repo): .venv/bin/python results/S1.4c/review/code/05_asym_A4A6.py [ncores]
"""
import os, sys, glob, json, importlib.util
import numpy as np, pandas as pd, shapely
from multiprocessing import Pool
def load(name, path):
    s = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
O2 = load("o2", "results/S1.4c/review/code/02_obs_recompute.py"); rb = load("rb", "tools/S1.4c_raster.py")
S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; EXP = SC + "/exp"; OUTD = "results/S1.4c/review/code"
R3 = "/mnt/micron/geo_spatialtrans/R3/real"; MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); AS = ["A4", "A6"]
RO = {A: list(rois.roi_id[rois.archetype == A]) for A in AS}
ks = O2.ks

def real_alt(A, roi):
    cc = pd.read_parquet(f"{EXP}/cells_{A}_{roi}.parquet"); ev = cc.ev_prep.to_numpy(); m = (cc.matched & np.isfinite(cc.a_n) & (cc.a_n > 0)).to_numpy()
    x, y = cc.x.to_numpy(), cc.y.to_numpy(); a_n = cc.a_n.to_numpy(); a_own = cc.a_own.to_numpy()
    r0 = np.sqrt(np.where(m, a_n, a_own) / np.pi)
    c = np.median(a_n[m] / a_own[m]); r1 = np.sqrt(np.where(m, a_n, a_own * c) / np.pi)
    out = dict(n_sr=len(cc), n_matched=int(m.sum()), n_ev=int(ev.sum()), calib=float(c))
    out["b8_R0"] = float(np.mean(O2.circle_hits(x, y, r0)[ev]))
    out["b8_R1"] = float(np.mean(O2.circle_hits(x, y, r1)[ev]))
    h = np.zeros(len(x), bool); idm = np.nonzero(m)[0]; h[idm] = O2.circle_hits(x[idm], y[idm], np.sqrt(a_n[idm] / np.pi)); out["b8_R2"] = float(np.mean(h[ev]))
    sg = pd.read_parquet(f"{R3}/{A}_{roi}_stardist_he_gen.parquet"); out["n_sd"] = len(sg)
    pos = pd.Series(np.arange(len(sg)), index=sg.label.to_numpy()); evl = cc.shape_label.to_numpy()[ev]; ie = pos.loc[evl].to_numpy()
    hs = O2.circle_hits(sg.x.to_numpy(), sg.y.to_numpy(), np.sqrt(sg.area_nuc.to_numpy() / np.pi)); out["b8_R3"] = float(np.mean(hs[ie]))
    reg = O2.regions()[(A, roi)]
    out["b7_R0"] = float(np.mean(np.sqrt(a_n[ev] / np.pi) > cc.d_gen.to_numpy()[ev]))
    xm = np.where(m, x, np.nan); ym = np.where(m, y, np.nan); _, dm, _ = O2.voronoi_obs(xm, ym, reg)
    out["b7_R2"] = float(np.mean(np.sqrt(a_n[ev] / np.pi) > dm[ev]))
    _, ds, _ = O2.voronoi_obs(sg.x.to_numpy(), sg.y.to_numpy(), reg)
    out["b7_R3"] = float(np.mean(np.sqrt(sg.area_nuc.to_numpy()[ie] / np.pi) > ds[ie]))
    # S-RM: maschera StarDist con le sole etichette appaiate (a una cella SR)
    lab = np.load(f"{MASKS}/{A}_{roi}_stardist_he_native.npz")["labels"].astype(np.int64)
    keep = np.zeros(int(lab.max()) + 1, bool); keep[cc.shape_label.to_numpy()[m].astype(np.int64)] = True; lab[~keep[lab]] = 0
    upp = float(sg.um_per_px.iloc[0]); mm = rb.measure(lab, upp, frac_out=False).set_index("label").loc[evl]
    out.update(S1_RM=float(np.mean(mm.n_c2 > 0)), _s2=(mm.n_c2 / mm.n_bnd).to_numpy().tolist(), _s3=mm.sol.to_numpy().tolist())
    rs = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_stardist_he_shape.parquet").set_index("label").loc[evl]
    out.update(S1_R0=float(np.mean(rs.n_c2 > 0)), _s2_0=(rs.n_c2 / rs.n_bnd).to_numpy().tolist(), _s3_0=rs.sol.to_numpy().tolist())
    return (A, roi), out

def sim_alt(f):
    b = os.path.basename(f)[:-8]; A, roi, fold, s, V = b.split("_")
    cc = pd.read_parquet(f"{EXP}/cells_{A}_{roi}.parquet"); ev = cc.ev_prep.to_numpy(); m = (cc.matched & np.isfinite(cc.a_n) & (cc.a_n > 0)).to_numpy(); n = len(cc)
    so = pd.read_parquet(f"{S}/sim_out/{b}.parquet")
    def v(col):
        z = np.full(n, np.nan); z[so.label.to_numpy() - 1] = so[col].to_numpy(); return z
    cx, cy, area = v("cx"), v("cy"), v("area")
    ok = np.isfinite(cx) & (area > 0); idm = np.nonzero(ok & m)[0]
    h = np.zeros(n, bool); h[idm] = O2.circle_hits(cx[idm], cy[idm], np.sqrt(area[idm] / np.pi))
    reg = O2.regions()[(A, roi)]
    _, dm, _ = O2.voronoi_obs(np.where(m, cx, np.nan), np.where(m, cy, np.nan), reg)
    e = ev & np.isfinite(cx) & np.isfinite(dm)
    # ri-rasterizzazione dei soli nuclei delle celle appaiate (operatore del banco)
    t = pd.read_parquet(f"{S}/gen_in/{b}.parquet"); keep = m[t.idx.to_numpy() - 1]
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side = int(rw.side_px)
    lab, _ = rb.rasterize(shapely.from_wkb(t.wkb.to_numpy()[keep]), t.idx.to_numpy()[keep], (side, side), upp)
    mm = rb.measure(lab, upp, frac_out=False).set_index("label"); ie = np.nonzero(ev)[0] + 1; ie = ie[np.isin(ie, mm.index)]; z = mm.loc[ie]
    return dict(archetype=A, roi=roi, fold=int(fold[4:]), seed_idx=int(s[1:]), variant={"Kc": "K-c", "K1": "K-1", "K2": "K-2"}[V],
                b8_M=float(np.mean(h[ev & ok])), b7_M=float(np.mean(np.sqrt(area[e] / np.pi) > dm[e])),
                s1_M=float(np.mean(z.n_c2 > 0)), L_s2=(z.n_c2 / z.n_bnd).to_numpy().tolist(), L_s3=z.sol.to_numpy().tolist())

if __name__ == "__main__":
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    with Pool(nc) as p:
        RA = dict(p.starmap(real_alt, [(A, r) for A in AS for r in RO[A]]))
        fs = sorted(f for A in AS for f in glob.glob(f"{S}/sim_out/{A}_*_s[01][0-9]_K[c12].parquet") if "_s00_" not in f)
        SA = p.map(sim_alt, fs, chunksize=1)
    pd.DataFrame([dict(archetype=k[0], roi=k[1], **{a: b for a, b in v.items() if not a.startswith("_")}) for k, v in RA.items()]).to_csv(f"{OUTD}/05_real_alternatives.csv", index=False)
    sa = pd.DataFrame(SA); obs = pd.read_csv("results/S1.4c/S1.4c_obs_raw.csv")
    key = ["archetype", "roi", "fold", "seed_idx", "variant"]; J = obs.merge(sa, on=key)
    # distanze per file nelle alternative
    J["B8_R1"] = [abs(r.b8 - RA[(r.archetype, r.roi)]["b8_R1"]) for r in J.itertuples()]
    J["B8_R2"] = [abs(r.b8_M - RA[(r.archetype, r.roi)]["b8_R2"]) for r in J.itertuples()]
    J["B8_R3"] = [abs(r.b8 - RA[(r.archetype, r.roi)]["b8_R3"]) for r in J.itertuples()]
    J["B7_R2"] = [abs(r.b7_M - RA[(r.archetype, r.roi)]["b7_R2"]) for r in J.itertuples()]
    J["B7_R3"] = [abs(r.b7 - RA[(r.archetype, r.roi)]["b7_R3"]) for r in J.itertuples()]
    J["S1_RM"] = [abs(r.s1_M - RA[(r.archetype, r.roi)]["S1_RM"]) for r in J.itertuples()]
    J["S2_RM"] = [ks(r.L_s2, RA[(r.archetype, r.roi)]["_s2"]) for r in J.itertuples()]
    J["S3_RM"] = [ks(r.L_s3, RA[(r.archetype, r.roi)]["_s3"]) for r in J.itertuples()]
    J.drop(columns=["L_s2", "L_s3"]).to_csv(f"{OUTD}/05_sim_alternatives.csv", index=False)
    loo = pd.read_csv(EXP + "/loo.csv")
    def yard(A, kr, fn): return np.array([fn(A, kr, j) for j in RO[A] if j != kr])
    alts = {"B8_R1": lambda A, a, b: abs(RA[(A, a)]["b8_R1"] - RA[(A, b)]["b8_R1"]), "B8_R2": lambda A, a, b: abs(RA[(A, a)]["b8_R2"] - RA[(A, b)]["b8_R2"]),
            "B8_R3": lambda A, a, b: abs(RA[(A, a)]["b8_R3"] - RA[(A, b)]["b8_R3"]), "B7_R2": lambda A, a, b: abs(RA[(A, a)]["b7_R2"] - RA[(A, b)]["b7_R2"]),
            "B7_R3": lambda A, a, b: abs(RA[(A, a)]["b7_R3"] - RA[(A, b)]["b7_R3"]), "S1_RM": lambda A, a, b: abs(RA[(A, a)]["S1_RM"] - RA[(A, b)]["S1_RM"]),
            "S2_RM": lambda A, a, b: ks(RA[(A, a)]["_s2"], RA[(A, b)]["_s2"]), "S3_RM": lambda A, a, b: ks(RA[(A, a)]["_s3"], RA[(A, b)]["_s3"]),
            "B8": lambda A, a, b: abs(RA[(A, a)]["b8_R0"] - RA[(A, b)]["b8_R0"]), "B7": lambda A, a, b: abs(RA[(A, a)]["b7_R0"] - RA[(A, b)]["b7_R0"])}
    rows = []
    for A in AS:
        for V in ["K-c", "K-1", "K-2"]:
            for met, fn in alts.items():
                fv = []
                for k in range(1, 6):
                    kr = loo.roi_test[(loo.archetype == A) & (loo.fold == k)].iloc[0]; z = J[(J.archetype == A) & (J.variant == V) & (J.fold == k)]
                    d = np.median(z[met]); yd = yard(A, kr, fn)
                    fv.append("PASS" if d <= np.median(yd) else ("WARN" if d <= yd.max() else "FAIL"))
                av = "PASS" if fv.count("PASS") >= 4 else ("FAIL" if fv.count("FAIL") >= 2 else "WARN")
                rows.append(dict(archetype=A, variant=V, metric=met, verdict=av, folds=",".join(fv)))
    vt = pd.DataFrame(rows); vt.to_csv(f"{OUTD}/05_alt_verdicts.csv", index=False)
    # punteggi sul pannello {B1,B2,B7,B8} sostituendo in A4/A6 B-7 e/o B-8 con le alternative
    rep = pd.read_csv("results/S1.4c/S1.4c_verdicts.csv"); rep = rep[rep.metric.isin(["B1", "B2", "B7", "B8"])]
    sc = []
    for b7, b8 in [("B7", "B8"), ("B7", "B8_R1"), ("B7", "B8_R2"), ("B7", "B8_R3"), ("B7_R2", "B8"), ("B7_R3", "B8"), ("B7_R2", "B8_R2"), ("B7_R3", "B8_R3")]:
        s = {}
        for V in ["K-c", "K-1", "K-2"]:
            base = rep[(rep.variant == V) & ~((rep.archetype.isin(AS)) & rep.metric.isin(["B7", "B8"]))]
            alt = vt[(vt.variant == V) & vt.metric.isin([b7, b8])]
            s[V] = int((base.verdict == "PASS").sum() + (alt.verdict == "PASS").sum())
        sc.append(dict(B7_A4A6=b7, B8_A4A6=b8, **s))
    sc = pd.DataFrame(sc); sc.to_csv(f"{OUTD}/05_alt_scores.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.read_csv(f"{OUTD}/05_real_alternatives.csv").to_string()); print(vt.to_string()); print(sc.to_string())
