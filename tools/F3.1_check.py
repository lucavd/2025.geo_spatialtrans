"""tools/F3.1_check.py - F3.1, cancello G-2: invarianti per area (vero, simulato, profili dal nucleo) e del riferimento.
(1) UMI in tessuto e quota di UMI nei nuclei ricalcolati dai file contro b2c_label.json; (2) blocchi di 8 um = somma dei 16 figli;
(3) blocchi passati a RCTD = blocchi con almeno 8 figli in tessuto; (4) riferimento: 25 <= cellule per tipo <= 3000.
Uso: .venv-b2c/bin/python tools/F3.1_check.py <main_dir>; stampa PASS/FAIL, esce con 1 se una fallisce."""
import sys, os, json, glob, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
main = sys.argv[1]; FAIL = []; n_ok = 0
def check(name, ok, detail=""):
    global n_ok
    if ok: n_ok += 1
    else: print(f"FAIL {name} [{detail}]", flush=True); FAIL.append(name)
dirs = sorted(glob.glob(f"{main}/real/*")) + sorted(glob.glob(f"{main}/sim/*_k*")) + sorted(glob.glob(f"{main}/simnuc/*_k*"))
for d in dirs:
    nm = "/".join(d.split("/")[-2:]); b = pd.read_parquet(f"{d}/bins.parquet"); g = pd.read_csv(f"{d}/genes.tsv", header=None)[0]
    T = pd.read_parquet(f"{d}/X.parquet"); X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(b), len(g)))
    u = np.asarray(X.sum(axis=1)).ravel(); reg = b.in_reg.to_numpy(); j = json.load(open(f"{d}/b2c_label.json"))
    check(f"{nm} (1) UMI in tessuto", float(u[reg].sum()) == j["umi_reg"], f"{u[reg].sum()} vs {j['umi_reg']}")
    lab = b.label.to_numpy(); sh = float(u[reg & (lab > 0)].sum() / u[reg].sum())
    check(f"{nm} (1) quota di UMI nei nuclei ricalcolata dalle etichette", abs(sh - j["share_assigned_nuc"]) < 1e-12, f"{sh} vs {j['share_assigned_nuc']}")
    key = (b.row.to_numpy() // 4).astype(np.int64) * 100000 + b.col.to_numpy() // 4
    bk, inv = np.unique(key, return_inverse=True); nreg = np.bincount(inv, weights=reg.astype(float), minlength=len(bk))
    blk = pd.read_csv(f"{d}/blocks.csv"); XB = sio.mmread(f"{d}/blocks_counts.mtx").tocsc()
    check(f"{nm} (3) blocchi = blocchi con >= 8 figli in tessuto", len(blk) == int((nreg >= 8).sum()), f"{len(blk)} vs {(nreg >= 8).sum()}")
    ub = np.bincount(inv, weights=u, minlength=len(bk))[nreg >= 8]
    check(f"{nm} (2) UMI del blocco = somma dei figli", np.array_equal(np.asarray(XB.sum(axis=0)).ravel(), ub) and np.array_equal(blk.umi.to_numpy(), ub))
    r = json.load(open(f"{d}/rctd.json")); check(f"{nm} (3) blocchi in RCTD", r["n_blocks"] == len(blk), f"{r['n_blocks']} vs {len(blk)}")
for f in sorted(glob.glob(f"{main}/ref/ref_*_cells.csv")):
    vc = pd.read_csv(f).cell_type.value_counts(); check(f"{os.path.basename(f)} (4) 25..3000 per tipo", bool(((vc >= 25) & (vc <= 3000)).all()), vc.to_dict())
print(f"{'FAIL' if FAIL else 'PASS'}: {n_ok} controlli passati, {len(FAIL)} falliti su {len(dirs)} aree"); sys.exit(1 if FAIL else 0)
