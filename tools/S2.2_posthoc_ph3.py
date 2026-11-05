"""S2.2, ANALISI POST HOC PH-3 (dichiarata dopo PH-1b/c): il sottotipo si estrae dalle quote dentro il decimo del totale della cellula
(tabella stimata sulle 4 aree di stima), cioè un totale per sottotipo; equivale a dare a ogni sottotipo la sua distribuzione dei totali.
Come PH-1: braccio «oracolo» in cui ogni cellula vera riceve come sottotipo il suo gruppo di R6 (la disposizione dei sottotipi nello spazio resta quella vera),
fonte R6 impronta intera, modelli M e G. Funzioni copiate senza modifiche da tools/S2.2_race.py. Testo originale della gara:
 (regole R-1, R-2), controprova CP-1 e CP-2, previsione P-5, su un tessuto.
Banco: geometria e totali veri (quadratini di valutazione di tools/S2.2_real.py), composizione simulata da ogni candidata;
ogni area vera è tenuta fuori a turno dalla stima dei profili di R6. Uso: .venv/bin/python tools/S2.2_race.py <A>
Scrive /mnt/micron/geo_spatialtrans/S2.2/race/<A>_<roi>.csv (una riga per braccio e replica; rilanciabile: salta le aree
già scritte), <A>_<roi>_real.json (misure vere sull'universo comune), <A>_<roi>_p5.npz. Pre-registrazione: docs/steps/S2.2.md."""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "tools")
import S22_common as S
import R4_common as R4C

M = "/mnt/micron/geo_spatialtrans"; REAL = f"{M}/S2.2/real"; OUT = f"{M}/S2.2/posthoc"; os.makedirs(OUT, exist_ok=True)
A = sys.argv[1]
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
ATLAS = {"A2": "crc", "A3": "crc", "A4": "ln", "A5": "brain"}.get(A)
NULL_OTHER = {"A1": "A5", "A5": "A1", "A2": "A6", "A3": "A6", "A4": "A2", "A6": "A4"}
MIN_SUB = 30; SEED0 = 20261106; REPS = 3; SEED_M4 = 20261106
ARMS = ["exp_M", "exp_G", "nuc_M", "nuc_G", "atlas_M", "atlas_G", "null_pooled", "null_other",
        "exp_M_nosub", "exp_G_nosub", "nuc_M_nosub", "nuc_G_nosub", "atlas_M_nosub", "atlas_G_nosub"]
TYPES = S.PRESET_TYPES[A]
rois_tab = pd.read_csv("results/R2/R2_rois_checked.csv"); rois_tab = rois_tab[rois_tab.archetype == A].reset_index(drop=True)
ROIS = list(rois_tab.roi_id)
t0 = time.time()

# ------------------------------------------------------------------ universo dei geni
vis = pd.read_csv(f"{M}/R6/counts/{KEY[A]}_genes.tsv", sep="\t")
if ATLAS and os.path.exists(f"{M}/S2.2/atlas/{ATLAS}_X.npz"):
    atl_genes = pd.read_csv(f"{M}/S2.2/atlas/{ATLAS}_genes.tsv", sep="\t")
    pos_atl = pd.Series(np.arange(len(atl_genes)), index=atl_genes.gene_id)
    in_both = vis.gene_id.isin(atl_genes.gene_id).to_numpy()
    GV = np.flatnonzero(in_both)                                   # colonne di Visium nell'universo
    GA = pos_atl.reindex(vis.gene_id[in_both]).to_numpy()          # colonne dell'atlante, stesso ordine
    HAVE_ATLAS = True
else:
    GV = np.arange(len(vis)); GA = None; HAVE_ATLAS = False
G = len(GV)
print(f"[{A}] universo {G} geni (Visium {len(vis)}); atlante {'si' if HAVE_ATLAS else 'no'}", flush=True)


# ------------------------------------------------------------------ cellule di stima delle candidate
def r6_cells(arch, roi, fp):
    """Nuclei di R6 di un'area con tipo del preset (classificabili), conteggi sull'universo."""
    tab = pd.read_parquet(f"{M}/R6/counts/{arch}_{roi}_nuclei.parquet")
    C = sp.load_npz(f"{M}/R6/counts/{arch}_{roi}_{fp}.npz").tocsr()
    lab = pd.read_parquet(f"{M}/R6/labels/{arch}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    lab = lab[lab.roi_id == roi].set_index("label")
    ty = lab.cell_type.reindex(tab.label).replace(S.MERGE).to_numpy()
    cl = lab.cluster.reindex(tab.label).to_numpy()
    ok = np.isin(ty, S.PRESET_TYPES[arch]) & pd.notna(cl)
    return C[ok], ty[ok], np.array([f"{t}|{c}" for t, c in zip(ty[ok], cl[ok])])


def merge_small(ty, sub):
    s = pd.Series(sub); n = s.map(s.value_counts())
    return np.where(n.to_numpy() < MIN_SUB, np.char.add(np.asarray(ty, dtype=str), "|rest"), sub)


def build_source(X, ty, sub, use_sub=True):
    """Profili per sottotipo (somme / totale), quote dei sottotipi nel tipo, phi per gene (momenti, dentro i sottotipi)."""
    X = sp.csr_matrix(X, dtype=np.float64)
    if not use_sub:
        sub = np.asarray(ty, dtype=str)
    sub = np.asarray(sub, dtype=str)
    us, inv = np.unique(sub, return_inverse=True)
    I = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(us), len(inv)))
    sums = np.asarray((I @ X).todense()); tot = sums.sum(axis=1)
    P = sums / np.where(tot > 0, tot, 1)[:, None]
    Nc = np.asarray(X.sum(axis=1)).ravel()
    sx2 = np.asarray(X.multiply(X).sum(axis=0)).ravel()
    sxN = np.asarray((I @ sp.diags(Nc) @ X).todense())            # sottotipi x geni: sum_c N_c x_cg
    sN2 = np.asarray(I @ (Nc ** 2)).ravel(); sN = np.asarray(I @ Nc).ravel()
    num = sx2 - 2 * (P * sxN).sum(axis=0) + (P ** 2 * sN2[:, None]).sum(axis=0) - (P * sN[:, None]).sum(axis=0)
    den = (P ** 2 * sN2[:, None]).sum(axis=0)
    phi = np.where(den > 0, np.maximum(num, 0) / np.where(den > 0, den, 1), 0.0)
    sub_type = np.array([u.split("|")[0] for u in us])
    shares = {}
    for t in np.unique(ty):
        m = sub_type == t
        cnt = np.bincount(inv, minlength=len(us))[m].astype(float)
        shares[t] = (np.flatnonzero(m), cnt / cnt.sum())
    return dict(P=P, phi=phi, shares=shares, names=us)


def pooled_source(X):
    X = sp.csr_matrix(X, dtype=np.float64); s = np.asarray(X.sum(axis=0)).ravel()
    return dict(P=(s / s.sum())[None, :], phi=np.zeros(X.shape[1]), shares=None, names=np.array(["pooled"]))


# atlante (fisso, senza aree)
if HAVE_ATLAS:
    ac = pd.read_parquet(f"{M}/S2.2/atlas/{ATLAS}_cells.parquet")
    AX = sp.load_npz(f"{M}/S2.2/atlas/{ATLAS}_X.npz").tocsr()[:, GA]
    ok = ac.cell_type.isin(TYPES).to_numpy()
    ac = ac[ok].reset_index(drop=True); AX = AX[ok]
    a_sub = merge_small(ac.cell_type.to_numpy(), (ac.cell_type + "|" + ac.subtype_raw.astype(str)).to_numpy())
    SRC_ATLAS = {True: build_source(AX, ac.cell_type.to_numpy(), a_sub, True), False: build_source(AX, ac.cell_type.to_numpy(), a_sub, False)}
    missing = sorted(set(TYPES) - set(ac.cell_type))
    print(f"[{A}] atlante {ATLAS}: {AX.shape[0]} cellule, tipi mancanti {missing}", flush=True)
# altro tessuto (nullo b): profilo medio di R6 impronta intera, 5 aree, sui geni dell'universo presenti nell'altro tessuto
o = NULL_OTHER[A]; vis_o = pd.read_csv(f"{M}/R6/counts/{KEY[o]}_genes.tsv", sep="\t"); pos_o = pd.Series(np.arange(len(vis_o)), index=vis_o.gene_id)
rois_o = pd.read_csv("results/R2/R2_rois_checked.csv").query("archetype == @o").roi_id
Xo = sp.vstack([r6_cells(o, r, "exp")[0] for r in rois_o]).tocsc()
so = np.asarray(Xo.sum(axis=0)).ravel()
po = np.zeros(G); idx = pos_o.reindex(vis.gene_id.iloc[GV]).to_numpy(); okg = ~np.isnan(idx)
po[okg] = so[idx[okg].astype(int)]
SRC_OTHER = dict(P=(po / po.sum())[None, :], phi=np.zeros(G), shares=None, names=np.array(["other"]))


# ------------------------------------------------------------------ simulazione sul banco
def simulate(src, noise, cells_ty, Nc, rng):
    """Conteggi per cellula (lista di geni per trascritto, in ordine di cellula). cells_ty: tipo; Nc: totale (sull'universo)."""
    n = len(Nc); P = src["P"]
    if src["shares"] is None:
        subidx = np.zeros(n, np.int64)
    else:
        subidx = np.empty(n, np.int64)
        for t in np.unique(cells_ty):
            m = np.flatnonzero(cells_ty == t)
            if t not in src["shares"]:
                raise RuntimeError(f"tipo {t} assente dalla candidata")
            ids, pr = src["shares"][t]
            subidx[m] = ids[rng.choice(len(ids), size=len(m), p=pr)]
    off = np.r_[0, np.cumsum(Nc)]
    genes = np.empty(off[-1], np.int64)
    if noise == "M":
        for s in np.unique(subidx):
            cs = np.flatnonzero(subidx == s)
            T = int(Nc[cs].sum())
            if T == 0:
                continue
            draws = rng.choice(G, size=T, p=P[s])
            slots = np.repeat(off[cs], Nc[cs]) + (np.arange(T) - np.repeat(np.r_[0, np.cumsum(Nc[cs])[:-1]], Nc[cs]))
            genes[slots] = draws
    else:
        phi = src["phi"]; pos = phi > 0
        shape = np.where(pos, 1.0 / np.where(pos, phi, 1.0), 1.0)
        for c in range(n):
            if Nc[c] == 0:
                continue
            lam = P[subidx[c]].copy()
            nzp = np.flatnonzero((lam > 0) & pos)
            lam[nzp] *= rng.gamma(shape[nzp], 1.0 / shape[nzp])
            lam /= lam.sum()
            x = rng.multinomial(int(Nc[c]), lam)
            g = np.flatnonzero(x)
            genes[off[c]:off[c + 1]] = np.repeat(g, x[g])
    return genes, subidx


def to_bins(genes, Nc, cell_of_bin, nb, n_bins, rng):
    """Distribuisce i trascritti di ogni cellula fra i suoi quadratini con i totali veri (ipergeometrica multivariata)."""
    cell_t = np.repeat(np.arange(len(Nc)), Nc)
    order = np.lexsort((rng.random(len(genes)), cell_t))
    g = genes[order]
    border = np.argsort(cell_of_bin, kind="stable")
    bin_t = np.repeat(border, nb[border])
    return sp.csr_matrix((np.ones(len(g), np.int32), (bin_t, g)), shape=(n_bins, G))


def cell_matrix(genes, Nc):
    cell_t = np.repeat(np.arange(len(Nc)), Nc)
    return sp.csr_matrix((np.ones(len(genes), np.int32), (cell_t, genes)), shape=(len(Nc), G))



for k, roi in enumerate(ROIS):
    fout = f"{OUT}/{A}_{roi}_ph3.csv"
    if os.path.exists(fout):
        continue
    Xr = sp.load_npz(f"{REAL}/{A}_{roi}_X.npz").tocsc()[:, GV].tocsr()
    bz = np.load(f"{REAL}/{A}_{roi}_bins.npz"); rows, cols, lab = bz["rows"], bz["cols"], bz["label"]
    cells = pd.read_parquet(f"{REAL}/{A}_{roi}_cells.parquet").reset_index(drop=True)
    pos_c = pd.Series(np.arange(len(cells)), index=cells.label); cob = pos_c.reindex(lab).to_numpy()
    nb = np.asarray(Xr.sum(axis=1)).ravel().astype(np.int64)
    Nc = np.bincount(cob, weights=nb, minlength=len(cells)).astype(np.int64)
    cty = cells.ptype.to_numpy(); ctype = dict(zip(cells.label, cells.ptype)); cxy = {l: (x, y) for l, x, y in zip(cells.label, cells.x_um, cells.y_um)}
    side = rois_tab.side_px[k] * rois_tab.um_per_px[k]
    interior = ((cells.x_um >= R4C.INTERIOR_UM) & (cells.y_um >= R4C.INTERIOR_UM) & (cells.x_um <= side - R4C.INTERIOR_UM) & (cells.y_um <= side - R4C.INTERIOR_UM)).to_numpy()
    Xc = Xr.tocoo(); Creal = sp.csr_matrix((Xc.data.astype(float), (cob[Xc.row], Xc.col)), shape=(len(cells), G))
    p5_sets = {}
    for t in np.unique(cty):
        m = np.flatnonzero(interior & (cty == t))
        if len(m) >= R4C.MIN_CELLS:
            gsel = R4C.eligible_genes(Creal[m])
            if len(gsel): p5_sets[t] = (m, gsel)
    def p5_stat(Cm):
        out = {}
        for t, (m, gsel) in p5_sets.items():
            Y = Cm[m]; lib = np.asarray(Y.sum(axis=1)).ravel().astype(float); sf = lib / lib.mean(); ok = sf > 0
            out[t] = R4C.nb_sf_excess(Y[ok][:, gsel], sf[ok])
        return out
    train = [r for r in ROIS if r != roi]
    parts = [r6_cells(A, r, "exp") for r in train]
    X_ = sp.vstack([p[0] for p in parts]).tocsc()[:, GV]; ty_ = np.concatenate([p[1] for p in parts]); sb_ = merge_small(ty_, np.concatenate([p[2] for p in parts]))
    src = build_source(X_, ty_, sb_, True)
    # tabella P(sottotipo | tipo, decimo del totale) sulle cellule di stima
    Ntr = np.asarray(X_.sum(axis=1)).ravel()
    name_pos = {n: i for i, n in enumerate(src["names"])}
    sub_tr = np.array([name_pos[s] for s in sb_])
    cuts = {t: np.quantile(Ntr[ty_ == t], np.linspace(0, 1, 11)[1:-1]) for t in np.unique(ty_)}
    dec_tr = np.array([np.searchsorted(cuts[t], n_, side="right") for t, n_ in zip(ty_, Ntr)])
    dec_ho = np.array([np.searchsorted(cuts[t], n_, side="right") for t, n_ in zip(cty, Nc)])
    rows_out, p5 = [], {}
    for noise in ("M", "G"):
        for r in range(REPS):
            rng = np.random.default_rng([SEED0 + 1000 * r + 970 + (noise == "G"), k])
            orc = np.empty(len(Nc), np.int64)
            for t in np.unique(cty):
                for d in range(10):
                    w = np.flatnonzero((cty == t) & (dec_ho == d))
                    if not len(w): continue
                    tab = np.bincount(sub_tr[(ty_ == t) & (dec_tr == d)], minlength=len(src["names"])).astype(float)
                    if tab.sum() == 0:
                        ids, pr = src["shares"][t]; tab = np.zeros(len(src["names"])); tab[ids] = pr
                    orc[w] = rng.choice(len(tab), size=len(w), p=tab / tab.sum())
            P_ = src["P"]
            off = np.r_[0, np.cumsum(Nc)]; genes = np.empty(off[-1], np.int64)
            if noise == "M":
                for s in np.unique(orc):
                    cs = np.flatnonzero(orc == s); T = int(Nc[cs].sum())
                    if T == 0: continue
                    draws = rng.choice(G, size=T, p=P_[s])
                    slots = np.repeat(off[cs], Nc[cs]) + (np.arange(T) - np.repeat(np.r_[0, np.cumsum(Nc[cs])[:-1]], Nc[cs]))
                    genes[slots] = draws
            else:
                phi = src["phi"]; pos = phi > 0; shape = np.where(pos, 1.0 / np.where(pos, phi, 1.0), 1.0)
                for c_ in range(len(Nc)):
                    if Nc[c_] == 0: continue
                    lam = P_[orc[c_]].copy(); nzp = np.flatnonzero((lam > 0) & pos)
                    lam[nzp] *= rng.gamma(shape[nzp], 1.0 / shape[nzp]); lam /= lam.sum()
                    x = rng.multinomial(int(Nc[c_]), lam); g_ = np.flatnonzero(x); genes[off[c_]:off[c_ + 1]] = np.repeat(g_, x[g_])
            Xs = to_bins(genes, Nc, cob, nb, len(nb), rng)
            m_ = S.bin_metrics(Xs, rows, cols, lab, ctype, cxy, SEED_M4); m_.update(archetype=A, roi_id=roi, arm=f"exp_{noise}_subtot", rep=r)
            rows_out.append(m_)
            for t, v in p5_stat(cell_matrix(genes, Nc)).items(): p5[f"exp_{noise}_subtot|{r}|{t}"] = v
    np.savez_compressed(f"{OUT}/{A}_{roi}_ph3_p5.npz", **p5)
    pd.DataFrame(rows_out).to_csv(fout, index=False)
    print(f"[{A}_{roi}] PH-3 ({time.time()-t0:.0f}s)", flush=True)
