"""tools/sonda/S3_areas.py - S3: 5 aree di tumore e 5 di stroma per sezione (docs/sonda/REGOLA_S3.md, «Aree»).
Quadrati di lato 1000 um (poi 750) con angolo su un reticolo da 100 um, interamente nella metà di sviluppo; copertura = quota dei
quadratini da 8 um in tissue sul numero atteso (lato/8)^2; punteggio r = log2((e + 1e-4) / (s + 1e-4)) sui quadratini in tissue
(matrice filtrata da 8 um, sole colonne della metà di sviluppo). Tumore: r massimo, 5 volte, togliendo i candidati che intersecano;
stroma: r minimo fra i disgiunti da tutte le aree scelte. Parità: x minore, poi y minore. Rilassamenti: (1000, 0.9), (1000, 0.8),
(750, 0.9), (750, 0.8). Uso: .venv/bin/python -B tools/sonda/S3_areas.py H2 H1  ->  results/sonda/S3/areas_<sec>.csv"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
E = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"]; S = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A3"]
RELAX = [(1000.0, 0.9), (1000.0, 0.8), (750.0, 0.9), (750.0, 0.8)]; LAT = 100.0
OUT = f"{C.REPO}/results/sonda/S3"; os.makedirs(OUT, exist_ok=True)

def integral(G):
    I = np.zeros((G.shape[0] + 1, G.shape[1] + 1)); I[1:, 1:] = G.cumsum(0).cumsum(1); return I

def boxsum(I, r0, c0, n):
    return I[r0 + n, c0 + n] - I[r0, c0 + n] - I[r0 + n, c0] + I[r0, c0]

def run(sec):
    m = C.split_m(sec); tp = C.positions(sec, "008"); h5 = f"{C.bdir(sec, '008')}/filtered_feature_bc_matrix.h5"
    ids, names = C.features(h5); bc = C.barcodes(h5)
    R, Cc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1
    x8 = (tp.array_col.to_numpy() + 0.5) * 8.0
    tis = np.zeros((R, Cc)); tis[tp.array_row.to_numpy(), tp.array_col.to_numpy()] = tp.in_tissue.to_numpy()
    # espressione: sole colonne in tissue della metà di sviluppo
    dev = tp[(tp.in_tissue.to_numpy() == 1) & (x8 < m - 50.0)]
    pos = pd.Series(np.arange(len(bc)), index=bc); col = pos.reindex(dev.barcode.to_numpy()).to_numpy()
    assert not np.isnan(col).any(); xmax = C.guard((dev.array_col.to_numpy() + 0.5) * 8.0, m)
    X, nruns = C.read_cols(h5, col.astype(np.int64))
    gE = [g for g in E if g in set(names)]; gS = [g for g in S if g in set(names)]
    nm = pd.Series(np.arange(len(names)), index=names)
    tot = np.asarray(X.sum(1)).ravel(); ue = np.asarray(X[:, nm[gE].to_numpy()].sum(1)).ravel(); us = np.asarray(X[:, nm[gS].to_numpy()].sum(1)).ravel()
    G = {k: np.zeros((R, Cc)) for k in ("tot", "e", "s")}
    rr, cc = dev.array_row.to_numpy(), dev.array_col.to_numpy()
    G["tot"][rr, cc] = tot; G["e"][rr, cc] = ue; G["s"][rr, cc] = us
    I = {k: integral(v) for k, v in G.items()}; It = integral(tis)
    out = None
    for lev, (side, cov) in enumerate(RELAX):
        n = int(round(side / 8.0)); cand = []
        for x0 in np.arange(0.0, m - 50.0 - side + 1e-9, LAT):
            c0 = int(np.ceil(x0 / 8.0 - 0.5))
            for y0 in np.arange(0.0, R * 8.0 - side + 1e-9, LAT):
                r0 = int(np.ceil(y0 / 8.0 - 0.5))
                if r0 + n > R or c0 + n > Cc: continue
                cv = boxsum(It, r0, c0, n) / n ** 2
                if cv < cov: continue
                t, e, s = boxsum(I["tot"], r0, c0, n), boxsum(I["e"], r0, c0, n), boxsum(I["s"], r0, c0, n)
                e, s = e / t, s / t
                cand.append((x0, y0, cv, e, s, float(np.log2((e + 1e-4) / (s + 1e-4)))))
        cd = pd.DataFrame(cand, columns=["x0", "y0", "coverage", "e", "s", "r"])
        chosen = []
        def free(row): return all(abs(row.x0 - q["x0"]) >= side or abs(row.y0 - q["y0"]) >= side for q in chosen)
        for tissue, asc in (("tumour", False), ("stroma", True)):
            srt = cd.sort_values(["r", "x0", "y0"], ascending=[asc, True, True])
            k = 0
            for _, row in srt.iterrows():
                if k == 5: break
                if free(row):
                    k += 1; flag = (row.e <= row.s) if tissue == "tumour" else (row.s <= row.e)
                    chosen.append(dict(section=sec, tissue=tissue, k=k, x0=row.x0, y0=row.y0, side=side, coverage=row.coverage, e=row.e, s=row.s, r=row.r,
                                       flag_wrong_side=bool(flag), relax_level=lev, coverage_min=cov, n_candidates=len(cd), m_um=m))
        if len(chosen) == 10: out = pd.DataFrame(chosen); break
        print(sec, "relax level", lev, "chosen", len(chosen), "candidates", len(cd), flush=True)
    assert out is not None, f"{sec}: fewer than 10 disjoint areas at every relaxation"
    out["xmax_area_um"] = out.x0 + out.side; assert (out.xmax_area_um <= m - 50.0).all()
    out.to_csv(f"{OUT}/areas_{sec}.csv", index=False)
    info = dict(section=sec, m_um=m, xmax_expression_read_um=xmax, n_bins_read=int(len(dev)), n_ranges=int(nruns), genes_E=gE, genes_S=gS,
                missing_E=sorted(set(E) - set(gE)), missing_S=sorted(set(S) - set(gS)), n_genes=int(len(names)), names_unique=bool(len(set(names)) == len(names)))
    json.dump(info, open(f"{OUT}/areas_{sec}.json", "w"), indent=1)
    print(json.dumps(info)); print(out[["tissue", "k", "x0", "y0", "side", "coverage", "e", "s", "r", "flag_wrong_side"]].round(4).to_string(index=False), flush=True)

for s in sys.argv[1:]: run(s)
