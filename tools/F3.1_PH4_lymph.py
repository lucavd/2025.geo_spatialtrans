"""tools/F3.1_PH4_lymph.py - F3.1, indagine PH-4 (dichiarata il 2026-11-18): linfociti B del linfonodo chiamati plasmacellule.
Classificatore: tipo di massima verosimiglianza multinomiale contro i profili medi per tipo del riferimento di F3.1 (stessi dati di RCTD; profilo
p_t = (somma dei conteggi + 0.5) / (totale + 0.5 G)). Applicato a: (a) profili veri generati delle cellule simulate (profiles.parquet), ridotti a caso
(seme 20261120) a totali estratti dagli UMI dei nuclei B veri, cosi' da confrontare quantita' di informazione simili; (b) conteggi dei nuclei veri con la
corona di 8 um (R6, impronta exp); (c) conteggi del solo nucleo dopo spostamento e cattura (quadratini con l'etichetta del nucleo, X.parquet), sul vero,
sul simulato con la semina uniforme (k = 0) e con i compartimenti (Q6). Per tipo vero (o di R6) e, per i B simulati, per sottotipo.
Uso: .venv/bin/python tools/F3.1_PH4_lymph.py <main_dir> <out_dir>"""
import sys, os, glob, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
sys.path.insert(0, "tools"); import S22_common as S
main, out = sys.argv[1:3]; M = "/mnt/micron/geo_spatialtrans"; A = "A4"; rng = np.random.default_rng(20261120)
Xr = sio.mmread(f"{main}/ref/ref_{A}_counts.mtx").tocsc(); rc = pd.read_csv(f"{main}/ref/ref_{A}_cells.csv")
rg = pd.read_csv(f"{main}/ref/ref_{A}_genes.tsv", header=None)[0].astype(str).to_numpy()
types = sorted(rc.cell_type.unique()); G = len(rg)
P = np.vstack([np.asarray(Xr[:, np.flatnonzero(rc.cell_type.to_numpy() == t)].sum(axis=1)).ravel() for t in types])
logP = np.log((P + 0.5) / (P.sum(axis=1, keepdims=True) + 0.5 * G))              # tipi x geni
gpos = pd.Series(np.arange(G), index=rg)
def classify(C):                                                                  # C: cellule x geni (ordine del riferimento)
    return np.array(types)[np.asarray(C @ logP.T).argmax(axis=1)]
def to_ref(C, genes):                                                             # riporta le colonne all'ordine del riferimento
    idx = gpos.reindex(genes).to_numpy(); ok = ~np.isnan(idx); C = C.tocsc()[:, np.flatnonzero(ok)]
    return sp.csr_matrix((C.tocoo().data, (C.tocoo().row, idx[ok].astype(int)[C.tocoo().col])), shape=(C.shape[0], G))
rows = []
def tally(src, area, truth, called, extra=None):
    d = pd.DataFrame({"truth": truth, "called": called})
    if extra is not None: d["sub"] = extra
    for (t, cl), n in d.groupby(["truth", "called"]).size().items(): rows.append({"source": src, "area": area, "truth": t, "called": cl, "n": int(n)})
    if extra is not None:
        for (sb, cl), n in d[d.truth == "linfociti B"].groupby(["sub", "called"]).size().items():
            rows.append({"source": src + "_subtype", "area": area, "truth": sb, "called": cl, "n": int(n)})
# (b) nuclei veri con la corona e totali dei B veri
lab = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
r6g = pd.read_csv(f"{M}/R6/counts/A4_genes.tsv", sep="\t").gene_id.astype(str).to_numpy()
Btot = []
for roi in ["f1", "f2", "p1", "p2", "p3"]:
    nuc = pd.read_parquet(f"{M}/R6/counts/{A}_{roi}_nuclei.parquet"); E = sp.load_npz(f"{M}/R6/counts/{A}_{roi}_exp.npz").tocsr()
    L = lab[lab.roi_id == roi].set_index("label"); ty = L.cell_type.reindex(nuc.label).replace(S.MERGE).to_numpy()
    ok = np.isin(ty, types) & (np.asarray(E.sum(axis=1)).ravel() >= 20)
    C = to_ref(E[np.flatnonzero(ok)], r6g); tally("real_exp", roi, ty[ok], classify(C))
    Btot.append(np.asarray(E[np.flatnonzero(ok & (ty == "linfociti B"))].sum(axis=1)).ravel())
Btot = np.concatenate(Btot)
# (a) profili veri generati del simulato, ridotti
for d in sorted(glob.glob(f"{main}/sim/A4_*_k0")):
    cells = pd.read_parquet(f"{d}/cells.parquet").drop_duplicates("label").set_index("label")
    pf = pd.read_parquet(f"{d}/profiles.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
    ids = np.sort(pf.label.unique()); pos = pd.Series(np.arange(len(ids)), index=ids)
    C = sp.csr_matrix((pf.n.to_numpy().astype(float), (pos.reindex(pf.label).to_numpy(), pf.gene.to_numpy())), shape=(len(ids), len(genes)))
    tgt = rng.choice(Btot, size=len(ids)); tot = np.asarray(C.sum(axis=1)).ravel(); C = C.tocoo()
    keep = rng.binomial(C.data.astype(np.int64), np.clip(tgt[C.row] / np.maximum(tot[C.row], 1), 0, 1))
    Ct = sp.csr_matrix((keep, (C.row, C.col)), shape=C.shape)
    ty = cells.cell_type.reindex(ids).to_numpy(); sb = cells.subtype.reindex(ids).to_numpy()
    tally("sim_profile_thinned", os.path.basename(d), ty, classify(to_ref(Ct, genes)), sb)
    tally("sim_profile_full", os.path.basename(d), ty, classify(to_ref(C.tocsr(), genes)), sb)
# (c) solo nucleo dopo spostamento e cattura
for side, pat in [("real_nucleus_bins", f"{main}/real/A4_*"), ("sim_nucleus_bins", f"{main}/sim/A4_*_k0"), ("simq6_nucleus_bins", f"{main}/simq6/A4_*_k0")]:
    for d in sorted(glob.glob(pat)):
        b = pd.read_parquet(f"{d}/bins.parquet", columns=["label", "in_reg"]); T = pd.read_parquet(f"{d}/X.parquet")
        genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
        X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(b), len(genes)))
        lb = b.label.to_numpy(); sel = lb > 0; ids, inv = np.unique(lb[sel], return_inverse=True)
        I = sp.csr_matrix((np.ones(sel.sum()), (inv, np.flatnonzero(sel))), shape=(len(ids), len(b))); C = I @ X
        cells = pd.read_parquet(f"{d}/cells.parquet").drop_duplicates("label").set_index("label")
        ty = cells.cell_type.reindex(ids).to_numpy(); ok = np.isin(ty, types) & (np.asarray(C.sum(axis=1)).ravel() >= 20)
        sb = cells.subtype.reindex(ids).to_numpy()[ok] if "subtype" in cells.columns else None
        tally(side, os.path.basename(d), ty[ok], classify(to_ref(C[np.flatnonzero(ok)], genes)), sb)
T = pd.DataFrame(rows); T.to_csv(f"{out}/F3.1_PH4_lymph_calls.csv", index=False)
S_ = T[~T.source.str.endswith("_subtype")].groupby(["source", "truth", "called"]).n.sum().reset_index()
S_["share"] = S_.n / S_.groupby(["source", "truth"]).n.transform("sum")
W = S_.pivot_table(index=["source", "truth"], columns="called", values="share", fill_value=0).round(3)
pd.set_option("display.width", 250); print(W.to_string())
sb = T[T.source == "sim_profile_thinned_subtype"].groupby(["truth", "called"]).n.sum().unstack(fill_value=0)
sb["share_plasma"] = (sb.get("plasmacellule", 0) / sb.sum(axis=1)).round(3); print(sb.sort_values("share_plasma", ascending=False).head(15).to_string())
# composizione per area (Q3)
comp = []
for side, pat in [("real", f"{main}/real/A4_*"), ("sim", f"{main}/sim/A4_*_k0")]:
    for d in sorted(glob.glob(pat)):
        ct = pd.read_parquet(f"{d}/cells.parquet").cell_type; ct = ct[ct.isin(types)]
        for t, v in ct.value_counts(normalize=True).items(): comp.append({"side": side, "area": os.path.basename(d)[:5], "type": t, "share": v})
C3 = pd.DataFrame(comp).pivot_table(index=["type"], columns=["side", "area"], values="share").round(3); C3.to_csv(f"{out}/F3.1_PH4_composition.csv"); print(C3.to_string())
