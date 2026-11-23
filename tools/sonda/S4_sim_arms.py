"""tools/sonda/S4_sim_arms.py - S4: bracci a componenti veri (docs/sonda/REGOLA_S4.md). Il banale B di tools/sonda/S3_sim_B.py con un
interruttore per componente; ogni componente e' B (come S3_sim_B) oppure V (preso dall'area omologa di H2, la sezione di taratura):
  geo  B: centri Poisson con la densita' di B, territori di Voronoi discreti, nucleo = disco di raggio r
       V: nuclei veri (etichette StarDist di V), centri = baricentri, territorio = quadratino in tissue al baricentro piu' vicino,
          poi i quadratini del nucleo vero tornano alla sua cellula
  typ  B: tipo estratto dalla composizione, indipendente
       V: tipo della cellula vera tipizzata piu' vicina (baricentri; con geo V e' la cellula stessa se tipizzata)
  dens B: n ~ Poisson(rho x quadratini del territorio), ogni trascritto in un quadratino del territorio a caso
       V: ogni quadratino in tissue riceve esattamente gli UMI di V in quel quadratino, attribuiti alla cellula del territorio
  prof B: geni ~ Multinomiale(n, profilo medio del tipo, parametri di B)
       V: geni ~ Multinomiale(n, profilo della cellula vera sorgente = conteggi del suo nucleo allargato da bin2cell in V, normalizzati);
          se quella cellula ha < 20 UMI si usa il profilo medio del tipo
Bracci: Y (typ V), G (geo V), U (dens V), A (geo, typ, dens V), AE (A + prof V). Seme default_rng([20261129, id del braccio, i, k]).
Fuori tessuto zero conteggi; nessuno spostamento, striscia o fondo. Uscite nel formato di F3.1 in $S4_ROOT/runs/<braccio>/<unit>_k<k>/
(bins, X, genes, cells) piu' S4_check.json e cells_src.parquet (cellula vera sorgente, solo per controllo).
Uso: .venv/bin/python -B tools/sonda/S4_sim_arms.py <W> <bracci,..|all> [i,i,..|all] [k]"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow as pa, pyarrow.parquet as pq
from scipy.spatial import cKDTree
from multiprocessing import Pool
SEED_S4 = 20261129; MIN_UMI_PROF = 20
ARMS = {"Y": dict(geo="B", typ="V", dens="B", prof="B"), "G": dict(geo="V", typ="B", dens="B", prof="B"),
        "U": dict(geo="B", typ="B", dens="V", prof="B"), "A": dict(geo="V", typ="V", dens="V", prof="B"),
        "AE": dict(geo="V", typ="V", dens="V", prof="V")}
ARM_ID = {a: j for j, a in enumerate(ARMS)}
UNITS = [("tumour", k) for k in range(1, 6)] + [("stroma", k) for k in range(1, 6)]
S3 = "/mnt/micron/geo_spatialtrans/S3_sonda"; ROOT = os.environ.get("S4_ROOT", "/mnt/micron/geo_spatialtrans/S4_sonda")
PAR = os.environ.get("S3_BPAR", f"{S3}/B")


def real_cells(x, y, lab):
    """etichette vere > 0 -> (ids ordinati, baricentri (n, 2), indice della cellula di ogni quadratino o -1)."""
    m = lab > 0; ids, inv = np.unique(lab[m], return_inverse=True)
    cnt = np.bincount(inv, minlength=len(ids)).astype(float)
    cxy = np.c_[np.bincount(inv, weights=x[m], minlength=len(ids)) / cnt, np.bincount(inv, weights=y[m], minlength=len(ids)) / cnt]
    pos = np.full(len(lab), -1, np.int64); pos[np.flatnonzero(m)] = inv
    return ids, cxy, pos


def simulate(x, y, in_tissue, lab_v, umi_v, Xexp, real_type, par, P, types, spec, rng):
    """Nucleo del simulatore (puro, testabile).
    x, y: centri dei quadratini (um); in_tissue: bool; lab_v: etichette dei nuclei veri per quadratino (0 = nessuno);
    umi_v: UMI veri per quadratino; Xexp: csr (cellule vere x geni) conteggi dei nuclei allargati, righe nell'ordine di ids di real_cells;
    real_type: array (cellule vere) di indici di tipo in `types` o -1 se non tipizzata; par: parametri di B; P: profili (tipi x geni).
    Restituisce label (nucleo per quadratino), righe/colonne dei trascritti, tipo per cellula, sorgente per cellula, diagnostica."""
    nb = len(x); tb = np.flatnonzero(in_tissue); ids, cxy, pos = real_cells(x, y, lab_v)
    if spec["geo"] == "B":
        N = int(rng.poisson(par["lambda_um2"] * 4.0 * len(tb)))
        pick = rng.integers(0, len(tb), N); cx = x[tb[pick]] + rng.uniform(-1, 1, N); cy = y[tb[pick]] + rng.uniform(-1, 1, N)
        dist, own = cKDTree(np.c_[cx, cy]).query(np.c_[x[tb], y[tb]])
        keep = np.unique(own); newid = np.zeros(N, np.int64); newid[keep] = np.arange(1, len(keep) + 1); n = len(keep)
        cell = newid[own]; centres = np.c_[cx, cy][keep]
        label = np.zeros(nb, np.int64); label[tb[dist <= par["r_um"]]] = cell[dist <= par["r_um"]]
        src_self = np.full(n, -1, np.int64)
    else:
        n = len(ids); centres = cxy
        _, own = cKDTree(cxy).query(np.c_[x[tb], y[tb]]); cell = own + 1
        nuc_t = pos[tb] >= 0; cell[nuc_t] = pos[tb][nuc_t] + 1           # il nucleo vero resta nella sua cellula
        label = np.zeros(nb, np.int64); label[tb[nuc_t]] = pos[tb][nuc_t] + 1
        src_self = np.arange(n)
    # sorgente vera: cellula vera tipizzata piu' vicina (con geo V la cellula stessa, se tipizzata)
    typed = np.flatnonzero(real_type >= 0)
    _, nn = cKDTree(cxy[typed]).query(centres); src = typed[nn]
    if spec["geo"] == "V":
        ok = real_type[src_self] >= 0; src[ok] = src_self[ok]
    if spec["typ"] == "B":
        ctype = rng.choice(len(types), size=n, p=np.array([par["composition"][t] for t in types]))
    else:
        ctype = real_type[src]
    size = np.bincount(cell, minlength=n + 1)[1:]
    if spec["dens"] == "B":
        ncnt = rng.poisson(par["rho_per_bin"] * size)
    else:
        m_b = np.zeros(nb, np.int64); m_b[tb] = umi_v[tb]
        ncnt = np.bincount(cell, weights=m_b[tb], minlength=n + 1)[1:].astype(np.int64)
    # geni per cellula
    gi_all, cc_all, rows_all = [], [], []
    order = np.argsort(cell, kind="stable"); start = np.r_[0, np.cumsum(size)[:-1]]; tb_sorted = tb[order]
    if spec["prof"] == "B":
        for t in range(len(types)):
            w = np.flatnonzero(ctype == t)
            for c0 in range(0, len(w), 200):
                ww = w[c0:c0 + 200]; Gm = rng.multinomial(ncnt[ww], P[t]); ci, gi = np.nonzero(Gm); mm = Gm[ci, gi]
                cc = np.repeat(ww[ci], mm); cc_all.append(cc); gi_all.append(np.repeat(gi, mm))
                if spec["dens"] == "B":                                   # stesso ordine di estrazione di S3_sim_B
                    off = (rng.random(len(cc)) * size[cc]).astype(np.int64); rows_all.append(tb_sorted[start[cc] + off])
        n_fallback = 0
    else:
        tot = np.asarray(Xexp.sum(axis=1)).ravel(); n_fallback = 0
        for c in range(n):
            if ncnt[c] == 0: continue
            s = src[c]
            if tot[s] >= MIN_UMI_PROF:
                row = Xexp.getrow(s); p = row.data / tot[s]; g = rng.multinomial(ncnt[c], p); gi = row.indices[g > 0]; mm = g[g > 0]
            else:
                n_fallback += 1; g = rng.multinomial(ncnt[c], P[ctype[c]]); gi = np.flatnonzero(g); mm = g[gi]
            cc_all.append(np.full(int(mm.sum()), c, np.int64)); gi_all.append(np.repeat(gi, mm))
    cc = np.concatenate(cc_all) if cc_all else np.zeros(0, np.int64); gg = np.concatenate(gi_all) if gi_all else np.zeros(0, np.int64)
    # quadratino di ogni trascritto
    if spec["dens"] == "B" and spec["prof"] == "B":
        rows = np.concatenate(rows_all) if rows_all else np.zeros(0, np.int64)
    elif spec["dens"] == "B":
        off = (rng.random(len(cc)) * size[cc]).astype(np.int64); rows = tb_sorted[start[cc] + off]
    else:
        o = np.lexsort((rng.random(len(cc)), cc)); cc = cc[o]; gg = gg[o]          # geni in ordine casuale dentro la cellula
        ob = np.argsort(cell, kind="stable"); bins_by_cell = tb[ob]
        rows = np.repeat(bins_by_cell, m_b[bins_by_cell])                            # quadratini ripetuti per i loro UMI, per cellula
        assert len(rows) == len(cc), (len(rows), len(cc))
    diag = {"n_cells": int(n), "n_real_cells": int(len(ids)), "n_real_typed": int(len(typed)), "n_fallback_profile": int(n_fallback),
            "share_type_equals_source": float(np.mean(ctype == real_type[src])) if n else None}
    return label, rows, gg, ctype, src, ids, diag


def run(arm, i, k=0):
    spec = ARMS[arm]; tissue, kk = UNITS[i]; src_dir = f"{S3}/runs/V/H2_{tissue}_{kk}"; dst = f"{ROOT}/runs/{arm}/H2_{tissue}_{kk}_k{k}"
    if os.path.exists(f"{dst}/export.json"): return f"present {dst}"
    t0 = time.time(); os.makedirs(dst, exist_ok=True)
    par = json.load(open(f"{PAR}/params_{tissue}.json")); z = np.load(f"{PAR}/params_{tissue}.npz"); P = z["profiles"]; types = list(z["types"])
    genes = open(f"{src_dir}/genes.tsv").read().split(); assert genes == list(z["genes"])
    bins = pd.read_parquet(f"{src_dir}/bins.parquet"); x, y = bins.x_um.to_numpy(), bins.y_um.to_numpy(); it = bins.in_tissue.to_numpy().astype(bool)
    T = pd.read_parquet(f"{src_dir}/X.parquet")
    Xv = sp.csr_matrix((T.n.to_numpy().astype(np.int64), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
    umi_v = np.asarray(Xv.sum(axis=1)).ravel().astype(np.int64)
    lab_v = bins.label.to_numpy().astype(np.int64); ids, _, _ = real_cells(x, y, lab_v)
    exp = np.load(f"{src_dir}/exp_labels_expanded.npy"); sel = np.isin(exp, ids); posi = pd.Series(np.arange(len(ids)), index=ids)
    Ie = sp.csr_matrix((np.ones(sel.sum()), (posi.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(ids), len(bins)))
    Xexp = (Ie @ Xv).tocsr()
    tl = pd.read_parquet(f"{src_dir}/labels_transferred.parquet").set_index("label").cell_type.reindex(ids)
    tix = {t: j for j, t in enumerate(types)}; real_type = np.array([tix.get(t, -1) if isinstance(t, str) else -1 for t in tl.to_numpy()], np.int64)
    rng = np.random.default_rng([SEED_S4, ARM_ID[arm], i, k])
    label, rows, gg, ctype, src, ids, diag = simulate(x, y, it, lab_v, umi_v, Xexp, real_type, par, P, types, spec, rng)
    X = sp.csr_matrix((np.ones(len(rows), np.int64), (rows, gg)), shape=(len(bins), len(genes))); X.sum_duplicates(); C = X.tocoo()
    umi = np.asarray(X.sum(axis=1)).ravel()
    # identita' dei componenti veri
    chk_id = {}
    if spec["dens"] == "V": chk_id["dens_identity_bins_equal"] = bool(np.array_equal(umi[it], umi_v[it]) and umi[~it].sum() == 0)
    if spec["geo"] == "V":
        chk_id["geo_identity_nuclei_equal"] = bool(np.array_equal(label[it] > 0, lab_v[it] > 0) and len(np.unique(label[label > 0])) == len(np.unique(lab_v[it & (lab_v > 0)])))
    if spec["typ"] == "V": chk_id["typ_identity_share"] = diag["share_type_equals_source"]
    assert all(v is True or v == 1.0 for v in chk_id.values()), chk_id
    out = bins.copy(); out["label"] = label; out.to_parquet(f"{dst}/bins.parquet", index=False)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": C.row.astype(np.int64), "gene": C.col.astype(np.int64), "n": C.data.astype(np.int64)}), preserve_index=False), f"{dst}/X.parquet")
    with open(f"{dst}/genes.tsv", "w") as f: f.write("\n".join(genes) + "\n")
    n = diag["n_cells"]
    pd.DataFrame({"label": np.arange(1, n + 1), "cell_type": np.array(types, dtype=object)[ctype]}).to_parquet(f"{dst}/cells.parquet", index=False)
    pd.DataFrame({"label": np.arange(1, n + 1), "src_real_label": ids[src]}).to_parquet(f"{dst}/cells_src.parquet", index=False)
    reg = bins.in_reg.to_numpy(); nbn = np.bincount(label[label > 0], minlength=n + 1)[1:]
    chk = {"unit": os.path.basename(dst), "arm": arm, "spec": spec, "i": i, "k": k, "seed": [SEED_S4, ARM_ID[arm], i, k], "src": src_dir, **diag, **chk_id,
           "lambda_target_mm2": par["lambda_mm2"], "lambda_realised_mm2": n / (it.sum() * 4.0) * 1e6, "n_nuc_target": par["n_nuc_mean"],
           "n_nuc_realised": float(nbn[nbn > 0].mean()), "n_cells_without_nucleus": int((nbn == 0).sum()), "rho_target": par["rho_per_bin"],
           "rho_realised_reg": float(umi[reg].sum() / reg.sum()), "rho_V_reg": float(umi_v[reg].sum() / reg.sum()), "umi_total": int(umi.sum()),
           "umi_out_of_tissue": int(umi[~it].sum()), "composition_realised": {types[t]: float((ctype == t).mean()) for t in range(len(types))},
           "sec": round(time.time() - t0, 1)}
    json.dump(chk, open(f"{dst}/S4_check.json", "w"), indent=1, default=float)
    json.dump(dict(chk, tissue=tissue, side_um=float(json.load(open(f"{src_dir}/export.json"))["side_um"]), n_bins=int(len(bins)), n_reg=int(reg.sum())),
              open(f"{dst}/export.json", "w"), default=float)
    return f"{arm} {chk['unit']} cells {n} lam {chk['lambda_realised_mm2']:.0f}/{chk['lambda_target_mm2']:.0f} rho {chk['rho_realised_reg']:.3f}/{chk['rho_target']:.3f} {chk_id} {chk['sec']} s"


def _one(a): return run(*a)

if __name__ == "__main__":
    W = int(sys.argv[1]); arms = list(ARMS) if sys.argv[2] == "all" else sys.argv[2].split(",")
    ii = list(range(10)) if len(sys.argv) < 4 or sys.argv[3] == "all" else [int(v) for v in sys.argv[3].split(",")]
    k = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    jobs = [(a, i, k) for a in arms for i in ii]
    with Pool(W) as p:
        for r in p.imap_unordered(_one, jobs): print(r, flush=True)
