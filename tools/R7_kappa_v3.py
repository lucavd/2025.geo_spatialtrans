"""tools/R7_kappa_v3.py — R7, correzione dopo la revisione avversariale (G-5):
(1) kappa ricalibrati con 10^6 repliche (prima 5 x 20 000); (2) fattore per gli intervalli «stima ± 2 errori standard»
(kappa_se), che mancava: lo stesso allargamento attorno alla stima che porta la copertura a 0.95; (3) fattori per il lato
simulato (aree simulate con semi indipendenti: nessuna zona comune, solo poche unità). Matrici: F_box, F_bin, F_nuc per
tessuto; correlazione implicita (M-4) per campo; identità per 5 e 4 unità.
-> results/R7/coverage/R7_kappa_v3_all.csv, R7_kappa_final_v3.csv"""
import sys, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R7_common as C
from scipy import stats
NBIG = 1_000_000; CH = 20_000
F = pd.read_csv("results/R7/pieces/R7_F_matrices.csv"); PR = pd.read_csv("results/R7/cp1/R7_cp1_pairs.csv")
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
def mat_from(df, A, col):
    ids = list(rois[rois.archetype == A].roi_id); M = np.eye(5)
    for _, r in df.iterrows():
        M[ids.index(r.i), ids.index(r.j)] = M[ids.index(r.j), ids.index(r.i)] = r[col]
    return M
cases = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    for col in ["F_box", "F_bin", "F_nuc"]: cases.append((A, col, mat_from(F[F.archetype == A], A, col)))
    for fld in ["umi", "nuclei"]: cases.append((A, f"implied_{fld}", mat_from(PR[(PR.archetype == A) & (PR.field == fld)], A, "R_implied")))
cases += [("sim", "identity_5", np.eye(5)), ("window", "identity_4", np.eye(4))]
def run(ix):
    A, nm, R = cases[ix]; n = len(R); L = np.linalg.cholesky(R)
    rng = np.random.default_rng([C.SEED, 3000 + ix]); nb, ns, rng_mm = [], [], []
    for _ in range(NBIG // CH):
        X = rng.standard_normal((CH, n)) @ L.T; lo, hi = C.exact_boot_ci(X)
        nb.append(C.needed_kappa(X, 0.0, lo, hi))
        se = X.std(1, ddof=1) / np.sqrt(n); ns.append(np.abs(X.mean(1)) / (2 * se))
        rng_mm.append(X.max(1) - X.min(1))
    nb, ns = np.concatenate(nb), np.concatenate(ns)
    kb, ks = float(np.quantile(nb, 0.95)), float(np.quantile(ns, 0.95))
    # errore Monte Carlo del quantile (binomiale sugli ordini statistici)
    def qse(v):
        v = np.sort(v); m = len(v); d = int(1.96 * np.sqrt(m * .95 * .05)); i = int(.95 * m)
        return float((v[i + d] - v[i - d]) / (2 * 1.96))
    Xi = rng.standard_normal((200_000, n)); kr = float(np.mean(Xi.max(1) - Xi.min(1)) / np.mean(np.concatenate(rng_mm)))
    return dict(archetype=A, matrix=nm, n=n, rho_bar=C.rho_bar(R), kappa_boot=kb, kappa_boot_se=qse(nb), kappa_se=ks, kappa_se_se=qse(ns), kappa_range=kr)
if __name__ == "__main__":
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as p: rows = p.map(run, range(len(cases)))
    D = pd.DataFrame(rows); D.to_csv("results/R7/coverage/R7_kappa_v3_all.csv", index=False)
    print(D.round(4).to_string(index=False))
    t4 = stats.t.ppf(.975, 4) / 2; print("controllo: kappa_se identità 5 =", D[D.matrix == "identity_5"].kappa_se.iloc[0], "atteso t4/2 =", t4)
    # tabella finale: R-2 per tipo; stroma e linfonodo dalla struttura misurata (R-1 dopo P-2 smentita), box = massimo dei due campi
    out = []
    for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
        g = D[D.archetype == A].set_index("matrix")
        if A in ("A3", "A4"): src = {"bin": "implied_umi", "nuc": "implied_nuclei"}
        else: src = {"bin": "F_bin", "nuc": "F_nuc", "box": "F_box"}
        r = dict(archetype=A, unit="area")
        for t in ("bin", "nuc"):
            for k in ("kappa_boot", "kappa_se", "kappa_range"): r[f"{k}_{t}"] = g.loc[src[t], k]
        for k in ("kappa_boot", "kappa_se", "kappa_range"):
            r[f"{k}_box"] = g.loc["F_box", k] if A not in ("A3", "A4") else max(g.loc["implied_umi", k], g.loc["implied_nuclei", k])
        r["source"] = "struttura misurata (R-1)" if A in ("A3", "A4") else "geometria (R-2)"
        out.append(r)
    for A, nm, lab in (("sim", "identity_5", "sim_5"), ("window", "identity_4", "window_4")):
        g = D[D.matrix == nm].iloc[0]; r = dict(archetype="all", unit=lab, source="unità indipendenti")
        for t in ("box", "bin", "nuc"): r.update({f"kappa_boot_{t}": g.kappa_boot, f"kappa_se_{t}": g.kappa_se, f"kappa_range_{t}": 1.0})
        out.append(r)
    pd.DataFrame(out).to_csv("results/R7/coverage/R7_kappa_final_v3.csv", index=False); print(pd.DataFrame(out).drop(columns="source").round(4).to_string(index=False))
