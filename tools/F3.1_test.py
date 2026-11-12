"""tools/F3.1_test.py - F3.1, cancello G-1: batteria di prove; stampa PASS/FAIL per asserzione, esce con 1 se una fallisce.
Uso: .venv-b2c/bin/python tools/F3.1_test.py <main_dir>   (dalla radice del repository)"""
import sys, os, json, glob, subprocess, tempfile, numpy as np, pandas as pd, scipy.sparse as sp, h5py
main = sys.argv[1]; FAIL = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok: FAIL.append(name)
def X_of(d, nb, ng):
    T = pd.read_parquet(f"{d}/X.parquet"); return sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(nb, ng))
# (a) esportazione simulata: UMI per quadratino = origine + fondo (tutte le aree simulate presenti)
for d in sorted(glob.glob(f"{main}/sim/*_k*"))[:60]:
    if not os.path.isdir(d) or not os.path.exists(f"{d}/X.parquet"): continue
    b = pd.read_parquet(f"{d}/bins.parquet"); g = pd.read_csv(f"{d}/genes.tsv", header=None)[0]
    u = np.asarray(X_of(d, len(b), len(g)).sum(axis=1)).ravel()
    og = pd.read_parquet(f"{d}/origin.parquet").groupby("bin_id").n.sum(); bg = pd.read_parquet(f"{d}/background.parquet")
    bgs = bg.set_index("bin_id")[bg.columns[-1]] if "bin_id" in bg.columns else pd.Series(dtype=float)
    exp = og.reindex(b.bin_id).fillna(0).to_numpy() + bgs.reindex(b.bin_id).fillna(0).to_numpy()
    check(f"a-sim {os.path.basename(d)}: UMI = origine + fondo per quadratino", np.array_equal(u, exp), f"max diff {np.abs(u-exp).max():.3g}")
# (a) esportazione vera: 300 quadratini a caso confrontati con la matrice grezza, letti in ordine di colonna (un'area per sezione)
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
rng = np.random.default_rng(20261116)
for A, ds in DS.items():
    d = sorted(glob.glob(f"{main}/real/{A}_*"))[0]; b = pd.read_parquet(f"{d}/bins.parquet"); g = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
    X = X_of(d, len(b), len(g)); sq = f"/mnt/micron/geo_spatialtrans/data_real/{ds}/binned_outputs_x/binned_outputs/square_002um"
    tp = pd.read_parquet(f"{sq}/spatial/tissue_positions.parquet")
    key = tp.set_index(["array_row", "array_col"]).barcode
    pick = rng.choice(len(b), 300, replace=False); bcs = key.reindex(list(zip(b.row.to_numpy()[pick], b.col.to_numpy()[pick]))).to_numpy().astype(str)
    with h5py.File(f"{sq}/raw_feature_bc_matrix.h5", "r") as h:
        m = h["matrix"]; allbc = m["barcodes"][:].astype(str); ip = m["indptr"]; ids = m["features"]["id"][:].astype(str)
        col = pd.Series(np.arange(len(allbc)), index=allbc); gpos = pd.Series(np.arange(len(g)), index=g)
        ok = True; nchk = 0
        cidx = np.array([col.get(bc, -1) for bc in bcs]); o = np.argsort(cidx); pick, bcs = pick[o], bcs[o]
        for k, bc in zip(pick, bcs):
            if bc not in col.index: continue
            c = col[bc]; s0, s1 = ip[c], ip[c + 1]; ind = m["indices"][s0:s1]; dat = m["data"][s0:s1]
            gid = ids[ind]; keep = np.isin(gid, g); v = np.zeros(len(g)); v[gpos.reindex(gid[keep]).to_numpy()] = dat[keep]
            ok &= np.array_equal(X[k].toarray().ravel(), v); nchk += 1
    check(f"a-real {os.path.basename(d)}: conteggi = matrice grezza sul pannello", ok and nchk > 250, f"{nchk} quadratini")
# (b) continuita' con S2.5 sulle 30 aree vere
s25 = pd.read_csv("/mnt/micron/geo_spatialtrans/S2.5/real/S2.5_real_areas.csv").set_index("tag")
for d in sorted(glob.glob(f"{main}/real/*")):
    e = json.load(open(f"{d}/export.json")); t = e["tag"]
    ok = e["n_reg"] == int(s25.loc[t, "reg_n2"]) and abs(e["umi_reg"] / e["n_reg"] - s25.loc[t, "reg_mean_umi2"]) <= 1e-9 * s25.loc[t, "reg_mean_umi2"]
    check(f"b {t}: quadratini in tessuto e UMI medio = S2.5", ok, f"{e['n_reg']} / {e['umi_reg']/e['n_reg']:.6f} vs {s25.loc[t,'reg_n2']} / {s25.loc[t,'reg_mean_umi2']:.6f}")
# (c) somma per cellula = bin_to_cell() di bin2cell (un'area vera per tessuto, nuclei senza allargamento)
import anndata as ad, bin2cell as b2c
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    d = sorted(glob.glob(f"{main}/real/{A}_*"))[0]; b = pd.read_parquet(f"{d}/bins.parquet"); g = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str)
    X = X_of(d, len(b), len(g)).astype(np.float32); a = ad.AnnData(X); a.obs_names = b.bin_id.astype(str).to_numpy(); a.var_names = g.to_numpy()
    a.obs["labels"] = b.label.to_numpy(); a.obsm["spatial"] = np.c_[b.x_um.to_numpy(), b.y_um.to_numpy()]
    a.uns["spatial"] = {"lib": {"scalefactors": {"spot_diameter_fullres": 1.0}}}
    a.obs["array_row"] = b.row.to_numpy(); a.obs["array_col"] = b.col.to_numpy()
    cd = b2c.bin_to_cell(a, labels_key="labels", spatial_keys=["spatial"])
    lab = b.label.to_numpy(); sel = lab > 0
    own = pd.DataFrame({"l": lab[sel], "u": np.asarray(X.sum(axis=1)).ravel()[sel]}).groupby("l").u.sum()
    b2 = pd.Series(np.asarray(cd.X.sum(axis=1)).ravel(), index=cd.obs_names.astype(int))
    check(f"c {os.path.basename(d)}: UMI per cellula = bin_to_cell()", np.allclose(own.reindex(b2.index).to_numpy(), b2.to_numpy()), f"{len(b2)} cellule")
# (d) griglia costruita: 2 cellule, 2 geni marcatori, risposte note
with tempfile.TemporaryDirectory() as td:
    nr = nc = 40; rr, cc = np.divmod(np.arange(nr * nc), nc)
    lab = np.zeros(nr * nc, int); lab[(rr >= 10) & (rr < 14) & (cc >= 10) & (cc < 14)] = 1; lab[(rr >= 10) & (rr < 14) & (cc >= 26) & (cc < 30)] = 2
    near1 = (rr >= 8) & (rr < 16) & (cc >= 8) & (cc < 16); near2 = (rr >= 8) & (rr < 16) & (cc >= 24) & (cc < 32)
    gA = np.where(near1, 9, 0) + np.where(near2, 1, 0); gB = np.where(near2, 9, 0) + np.where(near1, 1, 0)   # 10 % estranei per costruzione
    T = pd.DataFrame({"bin": np.r_[np.arange(nr*nc), np.arange(nr*nc)], "gene": np.r_[np.zeros(nr*nc, int), np.ones(nr*nc, int)], "n": np.r_[gA, gB]})
    T = T[T.n > 0]; T.to_parquet(f"{td}/X.parquet")
    pd.DataFrame({"bin_id": np.arange(nr*nc), "row": rr, "col": cc, "x_um": cc*2+1., "y_um": rr*2+1., "in_reg": np.ones(nr*nc, bool), "label": lab}).to_parquet(f"{td}/bins.parquet")
    pd.Series(["gA", "gB"]).to_csv(f"{td}/genes.tsv", index=False, header=False)
    pd.DataFrame({"label": [1, 2], "cell_type": ["tA", "tB"]}).to_parquet(f"{td}/cells.parquet")
    os.makedirs(f"{td}/ref", exist_ok=True); pd.DataFrame({"cell_type": ["tA", "tB"], "gene": ["gA", "gB"], "rank": [1, 1], "lfc": [3., 3.]}).to_csv(f"{td}/ref/markers_A9.csv", index=False)
    og = pd.DataFrame({"bin_id": np.r_[np.flatnonzero(near1), np.flatnonzero(near2)], "label": np.r_[np.ones(near1.sum(), int), 2*np.ones(near2.sum(), int)], "n": 10})
    og.to_parquet(f"{td}/origin.parquet")
    r = subprocess.run([sys.executable, "tools/F3.1_b2c.py", td, f"{td}/ref", "A9", "sim", "label"], capture_output=True, text=True)
    ok = r.returncode == 0
    if ok:
        j = json.load(open(f"{td}/b2c_label.json"))
        # allargamento fino a distanza euclidea 2 (in quadratini) dal nucleo 4x4: il quadrato 8x8 meno 3 quadratini per angolo = 52 quadratini,
        # 520 UMI per cellula; gli UMI sono nei due quadrati 8x8 (1 280): quota assegnata 1 040 / 1 280 = 0.8125
        check("d griglia: UMI assegnati = 0.8125", abs(j["share_assigned_exp"] - 0.8125) < 1e-12, str(j["share_assigned_exp"]))
        check("d griglia: UMI per cellula = 520", abs(j["umi_per_cell_exp"] - 520) < 1e-9, str(j["umi_per_cell_exp"]))
        check("d griglia: marcatori estranei = 0.1", abs(j["foreign_exp"] - 0.1) < 1e-12, str(j["foreign_exp"]))
        check("d griglia: purezza vera = 1, recupero vero = 0.8125", abs(j["truth_purity_exp"] - 1) < 1e-12 and abs(j["truth_recall_exp"] - 0.8125) < 1e-12,
              f'{j["truth_purity_exp"]} {j["truth_recall_exp"]}')
        check("d griglia: nucleo solo, UMI assegnati = 16*10*2/(64*10*2)", abs(j["share_assigned_nuc"] - 0.25) < 1e-12, str(j["share_assigned_nuc"]))
    else: check("d griglia: lo script gira", False, r.stderr[-300:])
# (e) mutanti e ricalcolo delle misure di RCTD dalle tabelle per blocco
for d in sorted(glob.glob(f"{main}/real/*")) + sorted(glob.glob(f"{main}/sim/*_k0"))[:6]:
    if not os.path.exists(f"{d}/rctd_blocks.csv"): continue
    m = pd.read_csv(f"{d}/rctd_blocks.csv"); j = json.load(open(f"{d}/rctd.json")); ref = pd.read_csv(sorted(glob.glob(f"{main}/ref/ref_{d.split('/')[-1][:2]}_cells.csv"))[0])
    ok_ = m.spot_class.notna() & (m.spot_class != "reject"); hn = ok_ & m.nuc_type.notna() & m.nuc_type.isin(ref.cell_type)
    ag = float((m["first"][hn] == m.nuc_type[hn]).mean()); sg = float((m.spot_class[ok_] == "singlet").mean())
    check(f"e {os.path.basename(d)}: concordanza e singoletti ricalcolati", abs(ag - j["agree_nucleus"]) < 1e-9 and abs(sg - j["frac_singlet"]) < 1e-9)
    rs = np.random.default_rng(1); nt = m.nuc_type[hn].to_numpy().copy(); rs.shuffle(nt)
    p = m["first"][hn].value_counts(normalize=True); q = m.nuc_type[hn].value_counts(normalize=True)   # caso: tipi di RCTD e del nucleo indipendenti sugli stessi blocchi; chance = float((p.reindex(q.index).fillna(0) * q).sum())
    agm = float((m["first"][hn].to_numpy() == nt).mean())
    check(f"e {os.path.basename(d)}: mutante tipi rimescolati al livello del caso", abs(agm - chance) < 0.03, f"{agm:.3f} vs {chance:.3f}")
for f in sorted(glob.glob(f"{main}/sim/*_k0/b2c_label_shift.json"))[:30]:
    j = json.load(open(f)); check(f"e {f.split('/')[-2]}: nuclei spostati, purezza vera < 0.1", j["truth_purity_exp"] < 0.1, f"{j['truth_purity_exp']:.3f}")
print(f"\n{'FAIL' if FAIL else 'PASS'}: {len(FAIL)} asserzioni fallite"); sys.exit(1 if FAIL else 0)
