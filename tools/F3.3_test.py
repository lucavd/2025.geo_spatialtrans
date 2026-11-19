"""tools/F3.3_test.py - F3.3, cancello G-1: batteria di prove, PASS/FAIL per asserzione.
(a) ritagli delle 30 aree di taratura uguali a R2; (b) maschere dei nuclei (Cellpose, Space Ranger, StarDist, consenso) e flag keep uguali
a R2 e R3b; (c) maschere valide a 1/4 e regioni uguali a R2 e ai livelli di S1.6; (d) esportazione uguale a F3.1 e misure di F3_3_common
uguali a F3.1 sulle 90 esecuzioni; (e) etichettatore su dati sintetici e mutanti; (f) statistiche su dati sintetici.
Uso: .venv-b2c/bin/python tools/F3.3_test.py [a,b,c,d,e,f]"""
import sys, os, json, glob, numpy as np, pandas as pd, tifffile
sys.path.insert(0, "tools")
M = "/mnt/micron/geo_spatialtrans"; SEG = f"{M}/F3.3/seg"; RUNS = f"{M}/F3.3/runs"
parts = sys.argv[1].split(",") if len(sys.argv) > 1 else list("abcdef")
res = []
def check(part, name, ok, detail=""):
    res.append((part, name, bool(ok), detail)); print(f"{'PASS' if ok else 'FAIL'} [{part}] {name} {detail}", flush=True)
A = pd.read_csv("results/F3.3/F3.3_areas.csv"); IN = A[A.level == "IN"]
if "a" in parts:
    for _, r in IN.iterrows():
        a = tifffile.imread(f"{SEG}/{r.area}/roi.tif"); b = tifffile.imread(f"{M}/R2/rois/{r.area}.tif")
        check("a", f"ritaglio {r.area}", a.shape == b.shape and np.array_equal(a, b))
if "b" in parts:
    KEYCOLS = ["label", "area_um2", "eq_diam_um", "x_um", "y_um", "in_tissue", "in_bubble", "flag_small", "flag_large", "keep"]
    for _, r in IN.iterrows():
        o = f"{SEG}/{r.area}"
        if r.seg == "cellpose_rgb":
            a = np.load(f"{o}/cellpose_rgb_native.npz")["labels"]; b = np.load(f"{M}/R2/masks/{r.area}_cellpose_rgb_native.npz")["labels"]
            check("b", f"Cellpose {r.area}", np.array_equal(a, b))
            ta = pd.read_parquet(f"{o}/cellpose_rgb_native.parquet"); tb = pd.read_parquet(f"results/R2/nuclei/{r.area}_cellpose_rgb_native.parquet")
            check("b", f"nuclei e keep {r.area}", len(ta) == len(tb) and all(np.array_equal(ta[c].to_numpy(), tb[c].to_numpy()) for c in KEYCOLS))
        else:
            for meth, ref in [("spaceranger", f"{M}/R2/masks/{r.area}_spaceranger_native.npz"), ("stardist_he", f"{M}/R2/masks/{r.area}_stardist_he_native.npz"),
                              ("consensus", f"{M}/R3b/masks/{r.area}_consensus_native.npz")]:
                a = np.load(f"{o}/{meth}_native.npz")["labels"]; b = np.load(ref)["labels"]
                check("b", f"{meth} {r.area}", np.array_equal(a.astype(np.int64), b.astype(np.int64)))
            ta = pd.read_parquet(f"{o}/consensus_native.parquet"); tb = pd.read_parquet(f"{M}/R3b/nuclei/{r.area}_consensus_native.parquet")
            check("b", f"nuclei e keep consenso {r.area}", len(ta) == len(tb) and all(np.array_equal(ta[c].to_numpy(), tb[c].to_numpy()) for c in KEYCOLS))
if "c" in parts:
    from skimage.io import imread
    for _, r in IN.iterrows():
        a = imread(f"{SEG}/{r.area}/valid_ds4.png"); b = imread(f"results/R2/tissue_masks/{r.area}_valid_ds4.png")
        check("c", f"maschera valida {r.area}", np.array_equal(a, b))
        ra = pd.read_csv(f"{SEG}/{r.area}/regions.tsv", sep="\t"); rb = pd.read_csv(f"{M}/S2.3/layers/{r.area}_regions.tsv", sep="\t")
        check("c", f"regioni {r.area}", len(ra) == len(rb) and (ra.wkb.to_numpy() == rb.wkb.to_numpy()).all(), f"{len(ra)} contro {len(rb)}")
if "d" in parts:
    for _, r in IN.iterrows():
        a, b = f"{RUNS}/{r.area}", f"{M}/F3.1/main/real/{r.area}"
        if not os.path.exists(f"{a}/bins.parquet"): check("d", f"esportazione {r.area}", False, "manca"); continue
        ba, bb = pd.read_parquet(f"{a}/bins.parquet"), pd.read_parquet(f"{b}/bins.parquet")
        xa, xb = pd.read_parquet(f"{a}/X.parquet"), pd.read_parquet(f"{b}/X.parquet")
        ca, cb = pd.read_parquet(f"{a}/cells.parquet"), pd.read_parquet(f"{b}/cells.parquet")
        ga, gb = open(f"{a}/genes.tsv").read(), open(f"{b}/genes.tsv").read()
        okb = ba.shape == bb.shape and all(np.array_equal(ba[c].to_numpy(), bb[c].to_numpy()) for c in bb.columns)
        okx = xa.shape == xb.shape and all(np.array_equal(xa[c].to_numpy(), xb[c].to_numpy()) for c in xb.columns)
        okc = ca.shape == cb.shape and np.array_equal(ca.label.to_numpy(), cb.label.to_numpy()) and (ca.cell_type.fillna("NA").to_numpy() == cb.cell_type.fillna("NA").to_numpy()).all()
        check("d", f"esportazione {r.area}", okb and okx and okc and ga == gb, f"quadratini {okb} conteggi {okx} nuclei {okc} geni {ga == gb}")
    devs = []
    for f in glob.glob(f"{M}/F3.3/calib/*.json"): devs.append(json.load(open(f))[0]["g1d_max_rel_dev"])
    check("d", "misure di F3_3_common uguali a F3.1 sulle 90 esecuzioni", len(devs) == 90 and max(devs) <= 1e-9, f"n = {len(devs)}, scarto massimo {max(devs):.1e}")
if "e" in parts:
    import scipy.sparse as sp, F3_3_common as L
    rng = np.random.default_rng(1); G = 400; K = 4; types = [f"t{k}" for k in range(K)]
    P = rng.dirichlet(np.full(G, 0.2), size=K)
    def synth(n, depth, prof):
        y = rng.integers(0, K, n); C = np.vstack([rng.multinomial(depth, prof[k]) for k in y]); return sp.csr_matrix(C), y
    C, y = synth(400, 500, P); Lg = np.log(P)
    pred = np.asarray(C @ Lg.T).argmax(axis=1); check("e", "riconoscimento a 500 UMI", (pred == y).mean() >= 0.99, f"{(pred == y).mean():.3f}")
    Psh = P[rng.permutation(K)]; predm = np.asarray(C @ np.log(Psh).T).argmax(axis=1)
    check("e", "mutante: profili rimescolati", (predm == y).mean() <= 1 / K + 0.05 or (predm == y).mean() < 0.5, f"{(predm == y).mean():.3f}")
    # fuori area: conteggio delle cellule del riferimento usate sulle aree di taratura
    rc = pd.read_csv(f"{M}/F3.1/main/ref/ref_A1_cells.csv"); _, _, _, n_out = L.ref_profiles("A1", pd.read_csv(f"{M}/F3.1/main/real/A1_r1/genes.tsv", header=None)[0].astype(str).to_numpy(), "r1")
    check("e", "profili fuori area (A1_r1)", n_out == int((rc.roi_id.astype(str) != "r1").sum()) and n_out < len(rc), f"{n_out} su {len(rc)}")
    _, _, _, n_all = L.ref_profiles("A1", pd.read_csv(f"{M}/F3.1/main/real/A1_r1/genes.tsv", header=None)[0].astype(str).to_numpy(), None)
    check("e", "mutante: profili non fuori area rilevati", n_all != n_out, f"{n_all} contro {n_out}")
if "f" in parts:
    import F3_3_stats as F
    rng = np.random.default_rng(7)
    def mk_cells(shift):
        cells = []
        for j in range(12):
            sd = 1.0; h = 1.6; n_in, n_out = 5, int(rng.integers(1, 6)); mu = rng.normal(0, 1)
            cells.append({"in": rng.normal(mu, sd, n_in), "out": {"s": rng.normal(mu + shift * h * rng.choice([-1, 1]), sd, n_out)}, "h": h, "sd": sd, "df": n_in - 1 + n_out - 1,
                          "sim": mu + rng.normal(0, 0.5), "good": True})
        return cells
    e0 = [F.corrected(mk_cells(0.0), F.stat_S, "perm", 200, 11 + i)["est"] for i in range(100)]
    check("f", "S* senza spostamento, permutazione", abs(np.mean(e0)) <= 0.15, f"media {np.mean(e0):.3f}")
    e0p = [F.corrected(mk_cells(0.0), F.stat_S, "param", 200, 11 + i)["est"] for i in range(100)]
    check("f", "S* senza spostamento, modello t", abs(np.mean(e0p)) <= 0.15, f"media {np.mean(e0p):.3f}")
    g0 = [F.corrected(mk_cells(0.0), F.stat_G, "perm", 200, 31 + i)["est"] for i in range(100)]
    check("f", "G* senza spostamento", abs(np.mean(g0)) <= 0.15, f"media {np.mean(g0):.3f}")
    e1 = [F.corrected(mk_cells(1.0), F.stat_S, "param", 200, 51 + i)["est"] for i in range(100)]
    check("f", "S* con spostamento di 1 h (statistica scritta; fallita il 2026-11-19, DEC-098)", abs(np.mean(e1) - 1) <= 0.25, f"media {np.mean(e1):.3f}")
    for k, fun in (("QS", F.stat_QS), ("QG", F.stat_QG)):
        q0 = [fun(mk_cells(0.0)) for _ in range(300)]
        check("f", f"{k} senza spostamento (quadrato medio)", abs(np.mean(q0)) <= 0.15, f"media {np.mean(q0):.3f} h^2")
        q1 = [F.q_estimate(mk_cells(1.0), fun, 200, 71 + i)["est"] for i in range(60)]
        check("f", f"{k} con spostamento di 1 h", abs(np.median(q1) - 1) <= 0.25, f"mediana {np.median(q1):.3f}")
    check("f", "regola V", F.vrule(0.1, 0.2, 0.5) == "confermata" and F.vrule(1.0, 0.2, 0.5) == "smentita" and F.vrule(0.4, 0.2, 0.5) == "indeterminata"
          and F.vrule(-1.0, 0.2, 0.5) == "smentita")
    vc = F.var_components([rng.normal(0, 1, 6) + 3 * i for i in range(3)], 500, 3)
    check("f", "componenti di varianza: rapporto alto quando i pazienti differiscono", vc["rho"] > 1.5, f"rho {vc['rho']:.2f}")
os.makedirs("results/F3.3", exist_ok=True)
pd.DataFrame(res, columns=["part", "assertion", "pass", "detail"]).to_csv(f"results/F3.3/F3.3_G1_{'_'.join(parts)}.csv", index=False)
print(f"{sum(r[2] for r in res)} PASS, {sum(not r[2] for r in res)} FAIL")
