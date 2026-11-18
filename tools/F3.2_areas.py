#!/usr/bin/env python3
"""tools/F3.2_areas.py - F3.2: mappa delle aree usate (G-4 a), classi dei gruppi di Space Ranger e filtro H&E (regola automatica del
file dello step), prova della regola sulle aree di R2 (G-4 d), estrazione delle aree di verifica (L1 sulle sezioni di taratura, L2/L3 sulle altre).
Coordinate: micrometri nel riferimento dell'immagine H&E a piena risoluzione (x = colonna * um/px, y = riga * um/px).
Uso (radice del repository, .venv/bin/python):
  tools/F3.2_areas.py used                   -> results/F3.2/F3.2_used_areas.csv
  tools/F3.2_areas.py classes <label>        -> results/F3.2/classes/<label>_clusters.csv (+ classe per quadratino di 8 um su micron)
  tools/F3.2_areas.py validate               -> results/F3.2/F3.2_rule_validation.csv (aree di R2, ritaglio di R1)
  tools/F3.2_areas.py select <label>         -> results/F3.2/F3.2_areas_<label>.csv   (prima regola, fallita: solo per memoria)
  tools/F3.2_areas.py validate2              -> results/F3.2/F3.2_rule2_validation.csv (seconda regola, G-4 d')
  tools/F3.2_areas.py select2 <label>        -> results/F3.2/F3.2_areas2_<label>.csv (fino a 10 aree per classe, in ordine di estrazione)
  tools/F3.2_areas.py sheet [classe]         -> results/F3.2/sheets/sheet_<classe>.png (foglio della seconda regola, non usato: regola fallita)
  tools/F3.2_areas.py pool <label>           -> results/F3.2/F3.2_pool_<label>.csv (candidate senza classe per il confronto a vista)
  tools/F3.2_areas.py pool_sheet <label>     -> results/F3.2/sheets/pool_<label>.png
Nessuna grandezza confrontata con il simulatore si calcola qui (vedi il file dello step, regola automatica, punto 4).
"""
import sys, os, json, math, csv, time
import numpy as np, pandas as pd, h5py, tifffile
import scipy.sparse as sp
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R1_verify import read_region
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/F3.2"; OUT.mkdir(parents=True, exist_ok=True)
MIC = Path("/mnt/micron/geo_spatialtrans")
SEAL = MIC / "F3.2_sealed"                       # classi per quadratino e maschere H&E dei dati di verifica (fuori git)
SEED = 20261118; MARGIN_UM = 280.0; SIDE_UM = 1000.0; TILE_UM = 50.0; TEST_TILE_UM = 250.0
CAL = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
       "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
TISSUE_OF = {"A1": "intestine", "A23": "crc", "A4": "lymph", "A5": "brain", "A6": "heart"}
HOLD = {"CRC_P1_SR300": ("Visium_HD_Human_Colon_Cancer_P1", "crc"), "CRC_P5_SR300": ("Visium_HD_Human_Colon_Cancer_P5", "crc"),
        "CRC_6p5mm_SR410": ("Visium_HD_6p5mm_Human_Colon_Cancer", "crc"), "MouseBrain_FFPE_SR401": ("Visium_HD_Mouse_Brain", "brain"),
        "HumanKidney_FFPE_SR401": ("Visium_HD_Human_Kidney_FFPE", "none"), "RatLiver_6p5mm_SR410": ("Visium_HD_6p5mm_Rat_Liver", "none")}
# pannello del file dello step: classe -> (insieme positivo, insiemi di esclusione)
PANEL = {
  "intestine": {"mucosa": (["Epcam", "Vil1", "Lgals4"], [["Acta2", "Des", "Myh11"]])},
  "crc": {"tumore": (["EPCAM", "KRT20", "CEACAM5"], [["CA1", "SLC26A3", "AQP8"], ["COL1A1", "COL3A1", "DCN"]]),
          "stroma": (["COL1A1", "COL3A1", "DCN"], [["EPCAM", "KRT20"]])},
  "lymph": {"follicolo": (["MS4A1", "CR2", "CD19"], [["CD3E", "CD3D", "CCR7"]]),
            "paracorticale": (["CD3E", "CD3D", "CCR7"], [["MS4A1", "CR2", "CD19"]])},
  "brain": {"corteccia": (["Slc17a7", "Satb2"], [["Mbp", "Plp1"], ["Ppp1r1b", "Penk"], ["Prox1"]])},
  "heart": {"miocardio": (["TNNT2", "MYH7", "ACTC1"], [["ADIPOQ", "PLIN1"]])},
  "none": {"tessuto": (None, [])},
}
R2_CLASS = {"A1": "mucosa", "A2": "tumore", "A3": "stroma", "A5": "corteccia", "A6": "miocardio"}

def paths(label):
    if label in CAL:
        d = MIC / "data_real" / CAL[label]; ds = CAL[label]; tissue = TISSUE_OF[label]
        bo = d / "binned_outputs_x/binned_outputs"
    else:
        ds, tissue = HOLD[label]; d = MIC / ("data_bridge" if label.endswith("P2_SR300") else "data_holdout") / label
        bo = d / "binned_outputs_x/binned_outputs"
    img = [p for p in d.glob(f"{ds}_tissue_image.*") if p.suffix in (".btf", ".tif", ".tiff")]
    return dict(dir=d, ds=ds, tissue=tissue, sq2=bo / "square_002um", sq8=bo / "square_008um", img=img[0] if img else None)

def upp_of(P):
    return json.load(open(P["sq8"] / "spatial/scalefactors_json.json"))["microns_per_pixel"]

def gap(a, b):
    """distanza fra due rettangoli (x0, y0, x1, y1) in um (0 se si toccano o si sovrappongono)"""
    dx = max(b[0] - a[2], a[0] - b[2], 0.0); dy = max(b[1] - a[3], a[1] - b[3], 0.0)
    return math.hypot(dx, dy)

# ---------------------------------------------------------------- mappa delle aree usate
def used():
    rows = []
    r2 = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")
    key_of_ds = {v: k for k, v in CAL.items()}
    for _, r in r2.iterrows():
        u = r.um_per_px
        rows.append(dict(section=key_of_ds[r.dataset], set="R2_ROI", id=f"{r.archetype}_{r.roi_id}", x0=r.c0 * u, y0=r.r0 * u, x1=r.c1 * u, y1=r.r1 * u,
                         source="results/R2/R2_rois_checked.csv", note=r.label))
    w = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv")
    for _, r in w.iterrows():
        u = r.um_per_px
        rows.append(dict(section=key_of_ds[r.dataset], set="R2b_WINDOW", id=r.win_id, x0=r.c0_img * u, y0=r.r0_img * u,
                         x1=(r.c0_img + r.side_px) * u, y1=(r.r0_img + r.side_px) * u, source="results/R2b/R2b_windows.csv",
                         note=f"dentro {r.archetype}_{r.roi_id}: {bool(r.inside_roi)}; usate anche da R2c, R3b, S1.4e-S1.5e"))
    t5 = pd.read_csv(ROOT / "results/R5/R5_tiles.csv")
    for key, g in t5.groupby("key"):
        P = paths(key); u = upp_of(P)
        pos = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet", columns=["array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        A = pos.set_index(["array_row", "array_col"])
        nr, nc = int(pos.array_row.max()) + 1, int(pos.array_col.max()) + 1
        for t in g.tile:
            R0, C0 = (int(v) for v in t[1:].split("_"))
            rr = [max(R0 - 75, 0), min(R0 + 150 + 75 - 1, nr - 1)]; cc = [max(C0 - 75, 0), min(C0 + 150 + 75 - 1, nc - 1)]
            pts = A.loc[[(a, b) for a in rr for b in cc]]
            x = pts.pxl_col_in_fullres.values * u; y = pts.pxl_row_in_fullres.values * u
            rows.append(dict(section=key, set="R5_EDGE_TILE", id=f"{key}_{t}", x0=x.min() - 1, y0=y.min() - 1, x1=x.max() + 1, y1=y.max() + 1,
                             source="results/R5/R5_tiles.csv", note="tessera 300 um piu' margine di 150 um letto dall'adattamento (R5, S2.0)"))
    for key in CAL:
        for fam, src in [("stripes", "tools/S2.0_stripes.py, S2.0_ph5_stripes_hp.py (sezione intera)"),
                         ("glass_background_second_tail", "tools/S2.3_glass.py, S2.4 R3 (vetro della sezione intera)"),
                         ("section_measures", "tools/R4_section.py, R4_posthoc.py (alternanza), S2.6/S2.7_section_measures.py")]:
            rows.append(dict(section=key, set="SECTION", id=f"{key}_{fam}", x0=np.nan, y0=np.nan, x1=np.nan, y1=np.nan, source=src,
                             note="sezione intera: per questi parametri le aree nuove della stessa sezione non sono fuori campione"))
    d = pd.DataFrame(rows); d.to_csv(OUT / "F3.2_used_areas.csv", index=False)
    print(d.groupby(["set"]).size().to_string()); print("->", OUT / "F3.2_used_areas.csv")

# ---------------------------------------------------------------- classi dei gruppi di Space Ranger
def read_h5(sq8):
    f = h5py.File(sq8 / "filtered_feature_bc_matrix.h5", "r")["matrix"]
    X = sp.csc_matrix((f["data"][:], f["indices"][:], f["indptr"][:]), shape=tuple(f["shape"][:]))
    bc = [b.decode() for b in f["barcodes"][:]]; names = [n.decode() for n in f["features/name"][:]]
    return X, bc, names

def classes(label, quiet=False):
    P = paths(label); panel = PANEL[P["tissue"]]
    cl = pd.read_csv(P["sq8"] / "analysis/clustering/gene_expression_graphclust/clusters.csv")
    cl.columns = ["barcode", "cluster"]
    pos = pd.read_parquet(P["sq8"] / "spatial/tissue_positions.parquet")
    pos = pos[pos.in_tissue == 1][["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"]]
    if P["tissue"] == "none":
        out = pos.assign(klass="tessuto", cluster=-1); tab = pd.DataFrame([dict(cluster=-1, klass="tessuto")])
    else:
        X, bc, names = read_h5(P["sq8"])
        lib = np.asarray(X.sum(axis=0)).ravel()
        genes = sorted({g for c in panel.values() for s in [c[0]] + c[1] for g in s})
        missing = [g for g in genes if g not in names]
        gi = {g: names.index(g) for g in genes if g in names}
        bidx = {b: i for i, b in enumerate(bc)}
        cl = cl[cl.barcode.isin(bidx)]
        cols = np.array([bidx[b] for b in cl.barcode])
        L = np.log1p(np.asarray(X[[gi[g] for g in gi]][:, cols].todense()) / np.maximum(lib[cols], 1) * 1e4)   # geni x quadratini
        M = pd.DataFrame(L.T, columns=list(gi)); M["cluster"] = cl.cluster.values
        mean_g = M.groupby("cluster").mean()                                   # media per gruppo, per gene
        def sc(s): return mean_g[[g for g in s if g in mean_g]].mean(axis=1)
        sets = {}
        for k, (posg, excl) in panel.items():
            sets[k + ":pos"] = sc(posg)
            for j, e in enumerate(excl): sets[f"{k}:excl{j}"] = sc(e)
        Z = pd.DataFrame({k: (v - v.mean()) / v.std(ddof=0) for k, v in sets.items()})
        lab = []
        for c in Z.index:
            best, bz = "altro", -np.inf
            for k in panel:
                z = Z.loc[c, k + ":pos"]
                ex = [Z.loc[c, x] for x in Z.columns if x.startswith(k + ":excl")]
                if z >= 1 and all(e < 1 for e in ex) and z > bz: best, bz = k, z
            lab.append(best)
        tab = Z.copy(); tab["klass"] = lab; tab["n_bins"] = M.groupby("cluster").size(); tab = tab.reset_index()
        tab["missing_genes"] = ";".join(missing)
        out = cl.merge(tab[["cluster", "klass"]], on="cluster").merge(pos, on="barcode")
    (OUT / "classes").mkdir(exist_ok=True); tab.to_csv(OUT / "classes" / f"{label}_clusters.csv", index=False)
    SEAL.mkdir(exist_ok=True); out.to_parquet(SEAL / f"{label}_bin8_class.parquet")
    if not quiet: print(label, tab[["cluster", "klass"]].groupby("klass").size().to_dict())
    return out, upp_of(P)

# ---------------------------------------------------------------- filtro H&E
def he_tiles(img, y0px, x0px, y1px, x1px, tile_px):
    """per le tessere di lato tile_px dentro [y0, y1) x [x0, x1): True = da scartare (vetro > 1/2 o grigio senza colore > 1/2)"""
    reg = read_region(img, y0px, y1px, x0px, x1px)
    if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
    reg = reg[..., :3].astype(np.float32)
    od = -np.log10((reg + 1.0) / 256.0).sum(axis=-1)
    mx, mn = reg.max(axis=-1), reg.min(axis=-1); sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-6), 0.0)
    glass = od < 0.15; grey = (sat < 0.08) & ~glass
    H, W = (reg.shape[0] // tile_px) * tile_px, (reg.shape[1] // tile_px) * tile_px
    def frac(m): return m[:H, :W].reshape(H // tile_px, tile_px, W // tile_px, tile_px).mean(axis=(1, 3))
    return (frac(glass) > 0.5) | (frac(grey) > 0.5)

def he_reject_fraction(img, upp, x0, y0, side_um):
    t = int(round(TILE_UM / upp)); n = int(round(side_um / TILE_UM))
    c0, r0 = int(round(x0 / upp)), int(round(y0 / upp))
    bad = he_tiles(img, r0, c0, r0 + n * t, c0 + n * t, t)
    return float(bad.mean())

# ---------------------------------------------------------------- prova sulle aree di R2 e sul ritaglio di R1
def validate():
    r2 = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv"); key_of_ds = {v: k for k, v in CAL.items()}
    res = []
    cache = {}
    for _, r in r2.iterrows():
        key = key_of_ds[r.dataset]
        if key not in cache: cache[key] = classes(key)
        b8, u = cache[key]
        x = b8.pxl_col_in_fullres * u; y = b8.pxl_row_in_fullres * u
        x0, y0, x1, y1 = r.c0 * u, r.r0 * u, r.c1 * u, r.r1 * u
        inside = b8[(x >= x0) & (x < x1) & (y >= y0) & (y < y1)]
        share = inside.klass.value_counts(normalize=True)
        want = R2_CLASS.get(r.archetype) or ("follicolo" if r.roi_id.startswith("f") else "paracorticale")
        maj = share.index[0] if len(share) else "nessuno"
        P = paths(key)
        rej = he_reject_fraction(P["img"], u, x0, y0, SIDE_UM)
        res.append(dict(kind="R2_ROI", id=f"{r.archetype}_{r.roi_id}", section=key, label=r.label, class_expected=want, class_majority=maj,
                        share_expected=float(share.get(want, 0.0)), n_bins8=len(inside), he_reject_frac=rej,
                        ok=(maj == want), meets_80=float(share.get(want, 0.0)) >= 0.8, he_ok=rej <= 0.05))
        print(res[-1], flush=True)
    # ritaglio 1 del linfonodo di R1 (Federica Pezzuto: «sicuramente non linfociti»), ricostruito con la stessa estrazione di R1_verify.py
    P = paths("A4"); u = upp_of(P)
    pos = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet"); it = pos[pos.in_tissue == 1]
    rng = np.random.default_rng(42); s = it.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
    half = int(round(400.0 / u / 2)); r0c, c0c = max(int(s.pxl_row_in_fullres) - half, 0), max(int(s.pxl_col_in_fullres) - half, 0)
    rej = he_reject_fraction(P["img"], u, c0c * u, r0c * u, 400.0)
    res.append(dict(kind="R1_crop", id="A4_crop1", section="A4", label="non rappresentativo (F. Pezzuto, R1)", class_expected="scartato",
                    class_majority=f"centro riga {int(s.pxl_row_in_fullres)} colonna {int(s.pxl_col_in_fullres)}", share_expected=np.nan,
                    n_bins8=0, he_reject_frac=rej, ok=rej > 0.05, meets_80=np.nan, he_ok=rej <= 0.05))
    d = pd.DataFrame(res); d.to_csv(OUT / "F3.2_rule_validation.csv", index=False)
    n_ok = int(d[d.kind == "R2_ROI"].ok.sum()); crop = bool(d[d.kind == "R1_crop"].ok.iloc[0])
    print(f"G-4 d: classe giusta in {n_ok}/30 aree di R2 (soglia 28); ritaglio di R1 scartato: {crop}; aree di R2 con H&E ok: {int(d[d.kind=='R2_ROI'].he_ok.sum())}/30")

# ---------------------------------------------------------------- estrazione
def select(label):
    b8, u = classes(label)
    P = paths(label); img = P["img"]; panel = PANEL[P["tissue"]]
    level = "L1" if label in CAL else "L2/L3"
    used_df = pd.read_csv(OUT / "F3.2_used_areas.csv")
    used_r = [tuple(v) for v in used_df[(used_df.section == label) & used_df.x0.notna()][["x0", "y0", "x1", "y1"]].values] if label in CAL else []
    # tessera di prova di StarDist (G-3): 250 um centrata sulla mediana dei quadratini in tessuto; conta come area usata
    cx, cy = float(np.median(b8.pxl_col_in_fullres)) * u, float(np.median(b8.pxl_row_in_fullres)) * u
    test = (cx - TEST_TILE_UM / 2, cy - TEST_TILE_UM / 2, cx + TEST_TILE_UM / 2, cy + TEST_TILE_UM / 2)
    blocked = used_r + [test]
    # griglia di tessere di 50 um sull'involucro del tessuto
    x = b8.pxl_col_in_fullres.values * u; y = b8.pxl_row_in_fullres.values * u
    gx0, gy0 = math.floor(x.min() / TILE_UM) * TILE_UM, math.floor(y.min() / TILE_UM) * TILE_UM
    nx, ny = int((x.max() - gx0) // TILE_UM) + 1, int((y.max() - gy0) // TILE_UM) + 1
    ix = ((x - gx0) // TILE_UM).astype(int); iy = ((y - gy0) // TILE_UM).astype(int)
    tiss = np.zeros((ny, nx)); np.add.at(tiss, (iy, ix), 1)
    # filtro H&E per tessere, a bande
    t = int(round(TILE_UM / u)); bad = np.ones((ny, nx), bool)
    for j in range(ny):
        r0 = int(round((gy0 + j * TILE_UM) / u)); c0 = int(round(gx0 / u))
        try:
            bt = he_tiles(img, r0, c0, r0 + t, c0 + nx * t, t)[0]; bad[j, :len(bt)] = bt
        except Exception as e:
            print("banda", j, "non letta:", e, flush=True)
    np.save(SEAL / f"{label}_he_bad_tiles.npy", bad)
    n = int(round(SIDE_UM / TILE_UM))
    def boxsum(A):
        S = np.zeros((A.shape[0] + 1, A.shape[1] + 1)); S[1:, 1:] = A.cumsum(0).cumsum(1)
        return S[n:, n:] - S[:-n, n:] - S[n:, :-n] + S[:-n, :-n]          # somma su ogni quadrato n x n (angolo in alto a sinistra)
    B = boxsum(bad.astype(float)) / (n * n); T = boxsum(tiss)
    rows = []
    for ci, klass in enumerate(sorted(panel)):
        K = np.zeros((ny, nx)); kk = (b8.klass.values == klass); np.add.at(K, (iy[kk], ix[kk]), 1)
        F = np.where(T > 0, boxsum(K) / np.maximum(T, 1), 0)
        ok = (F >= 0.8) & (B <= 0.05) & (T > 0)
        cand = np.argwhere(ok)                                            # (j, i) angolo in alto a sinistra, passo 50 um
        rng = np.random.default_rng([SEED, sorted(list(CAL) + list(HOLD)).index(label), ci])
        cand = cand[rng.permutation(len(cand))]
        got = []
        for j, i in cand:
            box = (gx0 + i * TILE_UM, gy0 + j * TILE_UM, gx0 + i * TILE_UM + SIDE_UM, gy0 + j * TILE_UM + SIDE_UM)
            if any(gap(box, b) < MARGIN_UM for b in blocked + got): continue
            got.append(box)
            rows.append(dict(label=label, level=level, klass=klass, area_id=f"{label}_{klass}_{len(got)}", x0=box[0], y0=box[1], x1=box[2], y1=box[3],
                             class_share=float(F[j, i]), he_reject_frac=float(B[j, i]), um_per_px=u, n_admissible=int(len(cand))))
            if len(got) == 5: break
        print(label, klass, "ammessi", len(cand), "estratte", len(got), flush=True)
    rows.append(dict(label=label, level="test_tile", klass="stardist_test", area_id=f"{label}_stardist_test", x0=test[0], y0=test[1], x1=test[2], y1=test[3],
                     class_share=np.nan, he_reject_frac=np.nan, um_per_px=u, n_admissible=0))
    d = pd.DataFrame(rows); d.to_csv(OUT / f"F3.2_areas_{label}.csv", index=False); print("->", OUT / f"F3.2_areas_{label}.csv")

# ================================================================ seconda regola (file dello step, «Regola automatica, seconda versione»; DEC-092)
LN2 = math.log(2.0); GREY_MAX = 0.02; COV_DEFAULT = 0.90
R2_TISSUE = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
REF_SECTION = {"intestine": "A1", "crc": "A23", "lymph": "A4", "brain": "A5", "heart": "A6"}

def set_list(tissue):
    out = []
    for k, (posg, excl) in PANEL[tissue].items():
        if posg is None: continue
        for s in [posg] + excl:
            if tuple(s) not in [tuple(x) for x in out]: out.append(s)
    return out

def bin_scores(label):
    """quadratini di 8 um in tessuto con x, y (um) e il punteggio di ogni insieme del pannello; cache sigillata su micron"""
    P = paths(label); u = upp_of(P); f = SEAL / f"{label}_bin8_setscores.parquet"
    if f.exists(): return pd.read_parquet(f), u
    pos = pd.read_parquet(P["sq8"] / "spatial/tissue_positions.parquet")
    pos = pos[pos.in_tissue == 1][["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"]].reset_index(drop=True)
    out = pd.DataFrame({"x": pos.pxl_col_in_fullres.values * u, "y": pos.pxl_row_in_fullres.values * u})
    sets = set_list(P["tissue"])
    if sets:
        X, bc, names = read_h5(P["sq8"]); bidx = pd.Series(np.arange(len(bc)), index=bc)
        cols = bidx.reindex(pos.barcode).to_numpy(); lib = np.asarray(X.sum(axis=0)).ravel()[cols]
        for s in sets:
            gi = [names.index(g) for g in s if g in names]
            L = np.log1p(np.asarray(X[gi][:, cols].todense()) / np.maximum(lib, 1) * 1e4)
            out["|".join(s)] = L.mean(axis=0)
    SEAL.mkdir(exist_ok=True); out.to_parquet(f); return out, u

def area_stats(B, box, sets, upp):
    """rapporti rho degli insiemi, copertura del tessuto (quadratini di 8 um in tessuto / attesi) nel riquadro box (um)"""
    m = (B.x >= box[0]) & (B.x < box[2]) & (B.y >= box[1]) & (B.y < box[3])
    exp_bins = ((box[2] - box[0]) / 8.0) * ((box[3] - box[1]) / 8.0)
    rho = {k: float(B.loc[m, k].mean() / B[k].mean()) for k in sets}
    return rho, float(m.sum() / exp_bins)

def grey_fraction(img, upp, box):
    t = int(round(TILE_UM / upp)); c0, r0 = int(round(box[0] / upp)), int(round(box[1] / upp))
    n = int(round((box[2] - box[0]) / TILE_UM)); reg = read_region(img, r0, r0 + n * t, c0, c0 + n * t)
    if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
    reg = reg[..., :3].astype(np.float32); od = -np.log10((reg + 1.0) / 256.0).sum(axis=-1)
    mx, mn = reg.max(axis=-1), reg.min(axis=-1); sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-6), 0.0)
    grey = (sat < 0.08) & (od >= 0.15); H, W = (grey.shape[0] // t) * t, (grey.shape[1] // t) * t
    return float((grey[:H, :W].reshape(H // t, t, W // t, t).mean(axis=(1, 3)) > 0.5).mean())

def classify(rho, refs):
    """classe k se max |ln(rho / ref_k)| <= ln 2 e minore della distanza dalle altre classi; altrimenti 'altro'"""
    d = {k: max(abs(math.log(max(rho[s], 1e-9) / max(r[s], 1e-9))) for s in r) for k, r in refs.items()}
    k = min(d, key=d.get)
    ok = d[k] <= LN2 and all(d[k] < d[j] for j in d if j != k)
    return (k if ok else "altro"), d

def r2_table():
    r2 = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv"); rows = []
    for _, r in r2.iterrows():
        u = r.um_per_px; want = R2_CLASS.get(r.archetype) or ("follicolo" if r.roi_id.startswith("f") else "paracorticale")
        rows.append(dict(id=f"{r.archetype}_{r.roi_id}", archetype=r.archetype, section=R2_TISSUE[r.archetype], klass=want,
                         box=(r.c0 * u, r.r0 * u, r.c1 * u, r.r1 * u), label=r.label))
    return rows

def r2_refs(tissue):
    """rapporti rho delle aree di R2 della sezione di taratura del tessuto, per classe e per area"""
    sec = REF_SECTION[tissue]; B, u = bin_scores(sec); sets = ["|".join(s) for s in set_list(tissue)]
    per = []
    for a in r2_table():
        if a["section"] != sec: continue
        rho, cov = area_stats(B, a["box"], sets, u); per.append(dict(a, rho=rho, cov=cov))
    return per, sets

def ref_of(per, sets, exclude=None):
    refs = {}
    for k in sorted({p["klass"] for p in per}):
        ps = [p for p in per if p["klass"] == k and p["id"] != exclude]
        refs[k] = {s: float(np.mean([p["rho"][s] for p in ps])) for s in sets}
    return refs

def cov_min(tissue):
    if tissue == "none": return COV_DEFAULT
    per, _ = r2_refs(tissue); return min(p["cov"] for p in per) - 0.05

def validate2():
    res = []
    for tissue in ["intestine", "crc", "lymph", "brain", "heart"]:
        per, sets = r2_refs(tissue); cmin = cov_min(tissue); P = paths(REF_SECTION[tissue]); u = upp_of(P)
        for p in per:
            refs = ref_of(per, sets, exclude=p["id"])                     # riferimento senza l'area in esame
            cls, d = classify(p["rho"], refs)
            grey = grey_fraction(P["img"], u, p["box"])
            rec = dict(kind="R2_ROI", id=p["id"], tissue=tissue, class_expected=p["klass"], class_assigned=cls, own_ok=(cls == p["klass"]),
                       coverage=p["cov"], cov_min=cmin, grey_frac=grey, admitted=(p["cov"] >= cmin and grey <= GREY_MAX))
            for k, v in d.items(): rec[f"d_{k}"] = v
            others = [k for k in refs if k != p["klass"]]
            rec["negatives"] = len(others); rec["negatives_rejected"] = sum(cls != k for k in others)
            res.append(rec); print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    # ritaglio 1 del linfonodo di R1
    P = paths("A4"); u = upp_of(P); B, _ = bin_scores("A4")
    pos = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet"); it = pos[pos.in_tissue == 1]
    rng = np.random.default_rng(42); s = it.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
    half = int(round(400.0 / u / 2)); r0c, c0c = max(int(s.pxl_row_in_fullres) - half, 0), max(int(s.pxl_col_in_fullres) - half, 0)
    box = (c0c * u, r0c * u, c0c * u + 400.0, r0c * u + 400.0)
    _, cov = area_stats(B, box, [], u); grey = grey_fraction(P["img"], u, box); cmin = cov_min("lymph")
    res.append(dict(kind="R1_crop", id="A4_crop1", tissue="lymph", class_expected="scartato", class_assigned="", own_ok=np.nan,
                    coverage=cov, cov_min=cmin, grey_frac=grey, admitted=(cov >= cmin and grey <= GREY_MAX), negatives=0, negatives_rejected=0))
    d = pd.DataFrame(res); d.to_csv(OUT / "F3.2_rule2_validation.csv", index=False)
    r = d[d.kind == "R2_ROI"]
    print(f"G-4 d': (i) classe propria {int(r.own_ok.sum())}/30 (soglia 28); (ii) negativi scartati {int(r.negatives_rejected.sum())}/{int(r.negatives.sum())} "
          f"(soglia 14/15); (iii) ammesse da copertura e artefatti {int(r.admitted.sum())}/30 (soglia 28); (iv) ritaglio di R1 scartato: {not bool(d[d.kind=='R1_crop'].admitted.iloc[0])}")

def select2(label, n_max=10):
    P = paths(label); tissue = P["tissue"]; B, u = bin_scores(label); img = P["img"]
    sets = ["|".join(s) for s in set_list(tissue)]
    if tissue != "none":
        per, _ = r2_refs(tissue); refs = ref_of(per, sets)
    else:
        refs = {"tessuto": {}}
    cmin = cov_min(tissue); level = "L1" if label in CAL else ("L3" if tissue == "none" else "L2")
    used_df = pd.read_csv(OUT / "F3.2_used_areas.csv")
    blocked = [tuple(v) for v in used_df[(used_df.section == label) & used_df.x0.notna()][["x0", "y0", "x1", "y1"]].values] if label in CAL else []
    cx, cy = float(np.median(B.x)), float(np.median(B.y))
    test = (cx - TEST_TILE_UM / 2, cy - TEST_TILE_UM / 2, cx + TEST_TILE_UM / 2, cy + TEST_TILE_UM / 2); blocked.append(test)
    gx0, gy0 = math.floor(B.x.min() / TILE_UM) * TILE_UM, math.floor(B.y.min() / TILE_UM) * TILE_UM
    nx, ny = int((B.x.max() - gx0) // TILE_UM) + 1, int((B.y.max() - gy0) // TILE_UM) + 1
    ix = ((B.x.values - gx0) // TILE_UM).astype(int); iy = ((B.y.values - gy0) // TILE_UM).astype(int)
    n = int(round(SIDE_UM / TILE_UM)); exp_bins = (SIDE_UM / 8.0) ** 2
    def boxsum(A):
        S = np.zeros((A.shape[0] + 1, A.shape[1] + 1)); S[1:, 1:] = A.cumsum(0).cumsum(1)
        return S[n:, n:] - S[:-n, n:] - S[n:, :-n] + S[:-n, :-n]
    cnt = np.zeros((ny, nx)); np.add.at(cnt, (iy, ix), 1); C = boxsum(cnt)
    R = {}
    for s in sets:
        a = np.zeros((ny, nx)); np.add.at(a, (iy, ix), B[s].values); R[s] = boxsum(a) / np.maximum(C, 1) / B[s].mean()
    cov = C / exp_bins
    rows = []; taken = []
    for ci, klass in enumerate(sorted(refs)):
        if tissue != "none":
            D = {k: np.max(np.stack([np.abs(np.log(np.maximum(R[s], 1e-9) / max(r[s], 1e-9))) for s in sets]), axis=0) for k, r in refs.items()}
            others = [D[j] for j in D if j != klass]
            okc = (D[klass] <= LN2) & np.all(np.stack([D[klass] < o for o in others]), axis=0) if others else (D[klass] <= LN2)
        else:
            okc = np.ones_like(cov, bool)
        cand = np.argwhere(okc & (cov >= cmin))
        rng = np.random.default_rng([SEED, sorted(list(CAL) + list(HOLD)).index(label), ci, 2]); cand = cand[rng.permutation(len(cand))]
        got = 0
        for j, i in cand:
            box = (gx0 + i * TILE_UM, gy0 + j * TILE_UM, gx0 + i * TILE_UM + SIDE_UM, gy0 + j * TILE_UM + SIDE_UM)
            if any(gap(box, b) < MARGIN_UM for b in blocked + taken): continue
            try: grey = grey_fraction(img, u, box)
            except Exception as e: print("H&E non letto", box, e, flush=True); continue
            if grey > GREY_MAX: continue
            taken.append(box); got += 1
            rows.append(dict(label=label, level=level, klass=klass, order=got, area_id=f"{label}_{klass}_{got}", x0=box[0], y0=box[1], x1=box[2], y1=box[3],
                             coverage=float(cov[j, i]), cov_min=cmin, grey_frac=grey, d_class=float(D[klass][j, i]) if tissue != "none" else np.nan,
                             um_per_px=u, n_admissible=int(len(cand))))
            if got == n_max: break
        print(label, klass, "ammessi", len(cand), "estratte", got, flush=True)
    rows.append(dict(label=label, level="test_tile", klass="stardist_test", order=0, area_id=f"{label}_stardist_test", x0=test[0], y0=test[1], x1=test[2], y1=test[3],
                     coverage=np.nan, cov_min=np.nan, grey_frac=np.nan, d_class=np.nan, um_per_px=u, n_admissible=0))
    pd.DataFrame(rows).to_csv(OUT / f"F3.2_areas2_{label}.csv", index=False); print("->", OUT / f"F3.2_areas2_{label}.csv")

def sheet(klass_filter=None):
    """fogli per il confronto a vista: per classe, le aree estratte (fino a 10 per sezione) accanto alle 5 aree di R2 della classe"""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    files = sorted(OUT.glob("F3.2_areas2_*.csv")); A = pd.concat([pd.read_csv(f) for f in files]); A = A[A.level != "test_tile"]
    th = ROOT / "results/R2/rois_thumbs"; out = OUT / "sheets"; out.mkdir(exist_ok=True)
    for klass, g in A.groupby("klass"):
        if klass_filter and klass != klass_filter: continue
        refs = [a for a in r2_table() if a["klass"] == klass]
        n = len(g); ncol = 5; nrow = 1 + math.ceil(n / ncol)
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.4 * nrow), facecolor="white")
        for ax in axs.ravel(): ax.set_axis_off()
        for k, a in enumerate(refs[:5]):
            im = plt.imread(th / f"{a['id']}.png"); axs[0, k].imshow(im[int(im.shape[0]*0.07):, :]); axs[0, k].set_title(f"RIFERIMENTO {a['id']}", fontsize=8, color="darkgreen")
        for k, (_, r) in enumerate(g.sort_values(["label", "order"]).iterrows()):
            P = paths(r.label); u = r.um_per_px; t = int(round(SIDE_UM / u)); c0, r0 = int(round(r.x0 / u)), int(round(r.y0 / u))
            reg = read_region(P["img"], r0, r0 + t, c0, c0 + t)
            if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
            s = t // 4 * 4; th4 = reg[:s, :s, :3].reshape(s // 4, 4, s // 4, 4, 3).mean(axis=(1, 3)).astype(np.uint8)
            ax = axs[1 + k // ncol, k % ncol]; ax.imshow(th4); ax.set_title(f"#{k+1}  {r.area_id}", fontsize=7)
        fig.suptitle(f"Classe: {klass}. Prima riga: aree di R2 approvate. Sotto: aree estratte. Somiglia / non somiglia?", fontsize=10)
        fig.tight_layout(); fig.savefig(out / f"sheet_{klass}.png", dpi=90); plt.close(fig); print(out / f"sheet_{klass}.png", n)

# ================================================================ selezione a vista (dopo il fallimento delle due regole; file dello step, DEC-092)
def pool(label):
    """aree candidate senza classe: copertura e artefatti della seconda regola, 280 um dalle aree usate e fra loro; 2 x 5 x (numero di classi)"""
    P = paths(label); tissue = P["tissue"]; B, u = bin_scores(label); img = P["img"]
    n_cls = len(PANEL[tissue]); n_want = 2 * 5 * n_cls; cmin = cov_min(tissue)
    level = "L1" if label in CAL else ("L3" if tissue == "none" else "L2")
    used_df = pd.read_csv(OUT / "F3.2_used_areas.csv")
    blocked = [tuple(v) for v in used_df[(used_df.section == label) & used_df.x0.notna()][["x0", "y0", "x1", "y1"]].values] if label in CAL else []
    cx, cy = float(np.median(B.x)), float(np.median(B.y))
    test = (cx - TEST_TILE_UM / 2, cy - TEST_TILE_UM / 2, cx + TEST_TILE_UM / 2, cy + TEST_TILE_UM / 2); blocked.append(test)
    gx0, gy0 = math.floor(B.x.min() / TILE_UM) * TILE_UM, math.floor(B.y.min() / TILE_UM) * TILE_UM
    nx, ny = int((B.x.max() - gx0) // TILE_UM) + 1, int((B.y.max() - gy0) // TILE_UM) + 1
    ix = ((B.x.values - gx0) // TILE_UM).astype(int); iy = ((B.y.values - gy0) // TILE_UM).astype(int)
    n = int(round(SIDE_UM / TILE_UM)); exp_bins = (SIDE_UM / 8.0) ** 2
    cnt = np.zeros((ny, nx)); np.add.at(cnt, (iy, ix), 1)
    S = np.zeros((ny + 1, nx + 1)); S[1:, 1:] = cnt.cumsum(0).cumsum(1); C = S[n:, n:] - S[:-n, n:] - S[n:, :-n] + S[:-n, :-n]
    cov = C / exp_bins; cand = np.argwhere(cov >= cmin)
    rng = np.random.default_rng([SEED, sorted(list(CAL) + list(HOLD)).index(label), 3]); cand = cand[rng.permutation(len(cand))]
    rows, taken = [], []
    for j, i in cand:
        box = (gx0 + i * TILE_UM, gy0 + j * TILE_UM, gx0 + i * TILE_UM + SIDE_UM, gy0 + j * TILE_UM + SIDE_UM)
        if any(gap(box, b) < MARGIN_UM for b in blocked + taken): continue
        try: grey = grey_fraction(img, u, box)
        except Exception as e: print("H&E non letto", box, e, flush=True); continue
        if grey > GREY_MAX: continue
        taken.append(box)
        rows.append(dict(label=label, level=level, order=len(taken), area_id=f"{label}_c{len(taken):02d}", x0=box[0], y0=box[1], x1=box[2], y1=box[3],
                         coverage=float(cov[j, i]), cov_min=cmin, grey_frac=grey, um_per_px=u, n_admissible=int(len(cand))))
        if len(taken) == n_want: break
    rows.append(dict(label=label, level="test_tile", order=0, area_id=f"{label}_stardist_test", x0=test[0], y0=test[1], x1=test[2], y1=test[3],
                     coverage=np.nan, cov_min=np.nan, grey_frac=np.nan, um_per_px=u, n_admissible=0))
    pd.DataFrame(rows).to_csv(OUT / f"F3.2_pool_{label}.csv", index=False)
    print(label, "ammessi", len(cand), "candidate", len(taken), "di", n_want, flush=True)

def pool_sheet(label):
    """foglio per il confronto a vista: in alto le aree di R2 approvate di ogni classe del tessuto, sotto le candidate numerate"""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    P = paths(label); tissue = P["tissue"]; A = pd.read_csv(OUT / f"F3.2_pool_{label}.csv"); A = A[A.level != "test_tile"]
    refs = [a for a in r2_table() if tissue != "none" and a["section"] == REF_SECTION[tissue]]
    klasses = sorted({a["klass"] for a in refs}); ncol = 5
    nref_rows = len(klasses); nrow = nref_rows + math.ceil(len(A) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.25 * nrow), facecolor="white")
    axs = np.atleast_2d(axs)
    for ax in axs.ravel(): ax.set_axis_off()
    th = ROOT / "results/R2/rois_thumbs"
    for ri, k in enumerate(klasses):
        for ci, a in enumerate([a for a in refs if a["klass"] == k][:5]):
            im = plt.imread(th / f"{a['id']}.png"); axs[ri, ci].imshow(im[int(im.shape[0] * 0.07):, :])
            axs[ri, ci].set_title(f"RIFERIMENTO {k.upper()}\n{a['id']}", fontsize=8, color="darkgreen")
    for k, (_, r) in enumerate(A.sort_values("order").iterrows()):
        u = r.um_per_px; t = int(round(SIDE_UM / u)); c0, r0 = int(round(r.x0 / u)), int(round(r.y0 / u))
        reg = read_region(P["img"], r0, r0 + t, c0, c0 + t)
        if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
        s = t // 4 * 4; th4 = reg[:s, :s, :3].reshape(s // 4, 4, s // 4, 4, 3).mean(axis=(1, 3)).astype(np.uint8)
        ax = axs[nref_rows + k // ncol, k % ncol]; ax.imshow(th4); ax.set_title(f"#{k + 1}", fontsize=12, fontweight="bold")
    q = (" / ".join(klasses) + " / nessuna") if klasses else "tessuto adatto? sì / no"
    fig.suptitle(f"{label}: per ogni candidata numerata, a quale riga di riferimento somiglia? ({q})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98)); (OUT / "sheets").mkdir(exist_ok=True)
    fig.savefig(OUT / "sheets" / f"pool_{label}.png", dpi=80); plt.close(fig); print(OUT / "sheets" / f"pool_{label}.png", len(A))

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "used": used()
    elif cmd == "classes": classes(sys.argv[2])
    elif cmd == "validate": validate()
    elif cmd == "select": select(sys.argv[2])
    elif cmd == "validate2": validate2()
    elif cmd == "select2": select2(sys.argv[2])
    elif cmd == "sheet": sheet(sys.argv[2] if len(sys.argv) > 2 else None)
    elif cmd == "pool": pool(sys.argv[2])
    elif cmd == "pool_sheet": pool_sheet(sys.argv[2])
