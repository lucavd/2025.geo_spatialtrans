"""tools/R7_kappa_final_v4.py — R7: tabella definitiva dei fattori dopo la revisione avversariale.
kappa_boot e kappa_se da R7_kappa_v3_all.csv (10^6 repliche); kappa_range ricalcolato con le stesse estrazioni per
numeratore e denominatore (correzione C-3); regola R-1 applicata per tessuto secondo gli esiti di P-2 rifatta
(results/R7/cp1_v2/R7_P2.csv): tessuto con almeno un campo smentito -> correlazione implicita (bin = campo UMI,
nuc = campo nuclei, box = il maggiore dei due); altrimenti matrici geometriche (R-2). -> results/R7/coverage/R7_kappa_final_v4.csv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C
D = pd.read_csv("results/R7/coverage/R7_kappa_v3_all.csv"); P2 = pd.read_csv("results/R7/cp1_v2/R7_P2.csv")
F = pd.read_csv("results/R7/pieces/R7_F_matrices.csv"); PR = pd.read_csv("results/R7/cp1_v2/R7_cp1_pairs.csv")
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
def mat_from(df, A, col):
    ids = list(rois[rois.archetype == A].roi_id); M = np.eye(5)
    for _, r in df.iterrows(): M[ids.index(r.i), ids.index(r.j)] = M[ids.index(r.j), ids.index(r.i)] = r[col]
    return M
mats = {}
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    for col in ["F_box", "F_bin", "F_nuc"]: mats[(A, col)] = mat_from(F[F.archetype == A], A, col)
    for fld in ["umi", "nuclei"]: mats[(A, f"implied_{fld}")] = mat_from(PR[(PR.archetype == A) & (PR.field == fld)], A, "R_implied")
kr = {k: C.kappa_range(R, np.random.default_rng([C.SEED, 5000 + i]), nrep=1_000_000) for i, (k, R) in enumerate(mats.items())}
D["kappa_range_crn"] = [kr.get((a, m), 1.0) for a, m in zip(D.archetype, D.matrix)]
D.to_csv("results/R7/coverage/R7_kappa_v4_all.csv", index=False)
out = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    g = D[D.archetype == A].set_index("matrix"); p = P2[P2.archetype == A]
    fired = bool((p.verdict == "smentita").any()) if len(p) else False
    src = {"bin": "implied_umi", "nuc": "implied_nuclei"} if fired else {"bin": "F_bin", "nuc": "F_nuc", "box": "F_box"}
    r = dict(archetype=A, unit="area", R1_fired=fired, P2=";".join(f"{f}:{v}" for f, v in zip(p.field, p.verdict)))
    for k, kk in (("kappa_boot", "kappa_boot"), ("kappa_se", "kappa_se"), ("kappa_range", "kappa_range_crn")):
        for t in ("bin", "nuc"): r[f"{k}_{t}"] = g.loc[src[t], kk]
        r[f"{k}_box"] = max(g.loc["implied_umi", kk], g.loc["implied_nuclei", kk]) if fired else g.loc["F_box", kk]
    r["source"] = "struttura misurata dentro le aree (R-1: P-2 smentita)" if fired else "geometria delle zone comuni (R-2)"
    out.append(r)
for nm, lab in (("identity_5", "sim_5"), ("identity_4", "window_4")):
    g = D[D.matrix == nm].iloc[0]; r = dict(archetype="all", unit=lab, R1_fired=False, P2="", source="unità indipendenti (solo poche unità)")
    for t in ("box", "bin", "nuc"): r.update({f"kappa_boot_{t}": g.kappa_boot, f"kappa_se_{t}": g.kappa_se, f"kappa_range_{t}": 1.0})
    out.append(r)
T = pd.DataFrame(out); T.to_csv("results/R7/coverage/R7_kappa_final_v4.csv", index=False)
print(T.drop(columns=["source"]).round(4).to_string(index=False))
