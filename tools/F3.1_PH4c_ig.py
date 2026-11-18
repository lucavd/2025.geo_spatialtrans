"""tools/F3.1_PH4c_ig.py - F3.1, PH-4 seguito (2026-11-18): da dove viene il segnale B/plasmacellula nei blocchi non B del linfonodo simulato.
(1) quota di UMI dei geni delle immunoglobuline (IGKC, JCHAIN, IGHA1, IGHA2, IGHE, IGHG1, IGHG3, IGHD, IGHM, IGLC1, IGLC7) per tipo: riferimento di F3.1,
profili veri generati (simulato k = 0, Q6), quadratini del nucleo (vero, simulato k = 0, Q6); (2) quota sul totale dell'area; (3) solo simulato, da origin.parquet:
quota degli UMI dei quadratini del nucleo che viene dalla cellula stessa e da altre cellule B o plasmacellule.
Uso: .venv/bin/python tools/F3.1_PH4c_ig.py <main_dir> <out_dir>"""
import sys, os, glob, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
main, out = sys.argv[1:3]; A = "A4"
IG = ["IGKC", "JCHAIN", "IGHA1", "IGHA2", "IGHE", "IGHG1", "IGHG3", "IGHD", "IGHM", "IGLC1", "IGLC7"]
gt = pd.read_csv("/mnt/micron/geo_spatialtrans/R6/counts/A4_genes.tsv", sep="\t"); igid = set(gt.gene_id[gt.gene_name.isin(IG)].astype(str))
print("geni Ig trovati:", len(igid))
rows = []
def add(src, area, ty, C, genes):
    isig = np.isin(genes, list(igid)); tot = np.asarray(C.sum(axis=1)).ravel(); ig = np.asarray(C[:, np.flatnonzero(isig)].sum(axis=1)).ravel()
    d = pd.DataFrame({"t": ty, "tot": tot, "ig": ig}).groupby("t")[["tot", "ig"]].sum()
    for t, r in d.iterrows(): rows.append({"source": src, "area": area, "type": t, "umi": r.tot, "ig": r.ig})
    if len(d) > 1: rows.append({"source": src, "area": area, "type": "TUTTE", "umi": tot.sum(), "ig": ig.sum()})
Xr = sio.mmread(f"{main}/ref/ref_{A}_counts.mtx").tocsr(); rc = pd.read_csv(f"{main}/ref/ref_{A}_cells.csv")
rg = pd.read_csv(f"{main}/ref/ref_{A}_genes.tsv", header=None)[0].astype(str).to_numpy(); types = sorted(rc.cell_type.unique())
add("reference", "all", rc.cell_type.to_numpy(), Xr.T.tocsr() if Xr.shape[0] == len(rg) else Xr, rg)
orow = []
for side, pat in [("real", f"{main}/real/A4_*"), ("sim", f"{main}/sim/A4_*_k0"), ("simq6", f"{main}/simq6/A4_*_k0")]:
    for d in sorted(glob.glob(pat)):
        area = os.path.basename(d); b = pd.read_parquet(f"{d}/bins.parquet", columns=["bin_id", "label"]); T = pd.read_parquet(f"{d}/X.parquet")
        genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
        X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(b), len(genes)))
        cells = pd.read_parquet(f"{d}/cells.parquet").drop_duplicates("label").set_index("label")
        add(side + "_area_all_bins", area, np.repeat("TUTTE_AREA", len(b)), X, genes)
        lb = b.label.to_numpy(); sel = lb > 0; ids, inv = np.unique(lb[sel], return_inverse=True)
        I = sp.csr_matrix((np.ones(sel.sum()), (inv, np.flatnonzero(sel))), shape=(len(ids), len(b))); C = (I @ X).tocsr()
        ty = cells.cell_type.reindex(ids).to_numpy(); ok = np.isin(ty, types)
        add(side + "_nucleus_bins", area, ty[ok], C[np.flatnonzero(ok)], genes)
        if side != "real":
            pf = pd.read_parquet(f"{d}/profiles.parquet"); pids = np.sort(pf.label.unique()); pos = pd.Series(np.arange(len(pids)), index=pids)
            P = sp.csr_matrix((pf.n.to_numpy().astype(float), (pos.reindex(pf.label).to_numpy(), pf.gene.to_numpy())), shape=(len(pids), len(genes)))
            add(side + "_profiles", area, cells.cell_type.reindex(pids).to_numpy(), P, genes)
            og = pd.read_parquet(f"{d}/origin.parquet"); bl = pd.Series(b.label.to_numpy(), index=b.bin_id.to_numpy())
            og["host"] = bl.reindex(og.bin_id).to_numpy(); og = og[og.host > 0]
            og["host_type"] = cells.cell_type.reindex(og.host).to_numpy(); og["src_type"] = cells.cell_type.reindex(og.label).to_numpy()
            og["own"] = og.label == og.host; og["from_Bpl"] = (~og.own) & og.src_type.isin(["linfociti B", "plasmacellule"])
            for t, g in og.groupby("host_type"):
                orow.append({"side": side, "area": area, "type": t, "umi": g.n.sum(), "own": g.n[g.own].sum(), "from_B_or_plasma": g.n[g.from_Bpl].sum()})
R = pd.DataFrame(rows); R.to_csv(f"{out}/F3.1_PH4c_ig_rows.csv", index=False)
S = R.groupby(["source", "type"])[["umi", "ig"]].sum(); S["ig_share"] = (S.ig / S.umi).round(4)
pd.set_option("display.width", 220); W = S.ig_share.unstack("source"); print(W.to_string()); W.to_csv(f"{out}/F3.1_PH4c_ig_share.csv")
O = pd.DataFrame(orow); O.to_csv(f"{out}/F3.1_PH4c_origin_rows.csv", index=False)
Og = O.groupby(["side", "type"])[["umi", "own", "from_B_or_plasma"]].sum(); Og["own_share"] = (Og.own / Og.umi).round(3); Og["B_pl_share"] = (Og.from_B_or_plasma / Og.umi).round(3)
print(Og[["own_share", "B_pl_share"]].to_string()); Og.to_csv(f"{out}/F3.1_PH4c_origin.csv")
