
# results/R3/review/claims/rv_verdicts.py — ricalcolo indipendente dei verdetti R3 dal testo pre-registrato
# (non usa R3_analysis.R). Uso: cd ~/2025.geo_spatialtrans && .venv/bin/python results/R3/review/claims/rv_verdicts.py
import numpy as np, pandas as pd, glob, os
from scipy.stats import spearmanr
R = "results/R3"; OUT = "/mnt/micron/geo_spatialtrans/R3"; RV = "results/R3/review/claims"
PRIM = dict(A1="cellpose_rgb", A2="cellpose_rgb", A3="cellpose_rgb", A4="spaceranger", A5="cellpose_rgb", A6="spaceranger")
cells = pd.read_parquet(f"{OUT}/cells_all.parquet")
cells = cells[cells.method == cells.archetype.map(PRIM)]
ci = cells[cells.interior].copy()
ci["eqr"] = np.sqrt(ci.area / np.pi); ci["ncr"] = ci.area_nuc / ci.area; ci["rr"] = np.sqrt(ci.ncr)
rows = []
def add(id_, A, stat, mine, claimed_val, my_assert, claimed_assert, pred, note=""):
    rows.append(dict(id=id_, archetype=A, statistic=stat, value_recomputed=mine, value_claimed=claimed_val,
                     assertion_recomputed=my_assert, assertion_claimed=claimed_assert, predicted=pred, note=note))
V = pd.read_csv(f"{R}/R3_verdicts.csv")
def cl(id_, A):
    z = V[(V.id == id_) & (V.archetype == A)]
    return (z.value.iloc[0], z.assertion.iloc[0]) if len(z) == 1 else (np.nan, f"n={len(z)}")
pf = lambda b: "PASS" if b else "FAIL"
# B-R3.1
for A in PRIM:
    m = ci.eqr[ci.archetype == A].median(); c = cl("B-R3.1", A)
    add("B-R3.1", A, "eq_r mediano pool interne", m, *[c[0]], pf(5 <= m <= 25), c[1], "FAIL" if A == "A4" else "PASS")
# B-R3.2
cat = dict(A1=.45, A2=.45, A3=.45, A4=.55)
for A in PRIM:
    m = ci.rr[ci.archetype == A].median(); c = cl("B-R3.2", A)
    a = pf(abs(m / cat[A] - 1) <= .2) if A in cat else "registrato"
    add("B-R3.2", A, "rapporto raggi mediano", m, c[0], a, c[1], dict(A1="FAIL", A2="PASS", A3="FAIL", A4="PASS").get(A, "registrare"))
nc = {A: ci.ncr[ci.archetype == A].median() for A in PRIM}
c = cl("B-R3.3a", "A4"); add("B-R3.3a", "A4", "N/C mediano", nc["A4"], c[0], pf(.8 <= nc["A4"] <= .9), c[1], "FAIL")
c = cl("B-R3.3b", "A3"); add("B-R3.3b", "A3", "N/C A3 < A1 e < A2", nc["A3"], c[0], pf(nc["A3"] < nc["A1"] and nc["A3"] < nc["A2"]), c[1], "PASS",
                           f"A1={nc['A1']:.4f} A2={nc['A2']:.4f}")
# B-R3.4 (pool) + per-ROI robustness
nsup = 0; roi_sp = []
for A in PRIM:
    d = ci[ci.archetype == A]; s = spearmanr(d.area_nuc, d.area).statistic; nsup += s >= .3
    per = [spearmanr(g.area_nuc, g.area).statistic for _, g in d.groupby("roi_id")]
    roi_sp.append(dict(archetype=A, pool=s, roi_min=min(per), roi_max=max(per), roi_n_ge_03=sum(np.array(per) >= .3)))
    c = cl("B-R3.4", A); add("B-R3.4", A, "Spearman pool", s, c[0], pf(s >= .3), c[1], "registrare",
                             f"per ROI {min(per):.3f}..{max(per):.3f}, ROI>=0.3: {sum(np.array(per) >= .3)}/5")
c = cl("B-R3.4", "tutti"); add("B-R3.4", "tutti", "archetipi sostenuti", nsup, c[0], pf(nsup <= 2), c[1], "PASS")
pd.DataFrame(roi_sp).to_csv(f"{RV}/rv_B34_spearman_roi.csv", index=False)
# CP-1 from null summary + roi summary (not from cp1 csv)
ns = pd.read_csv(f"{R}/R3_null_summary.csv"); rr = pd.read_csv(f"{R}/R3_roi_summary.csv")
rr = rr[rr.method == rr.archetype.map(PRIM)]
cp1 = []
for _, r in rr.iterrows():
    cs = ns[(ns.archetype == r.archetype) & (ns.roi_id == r.roi_id) & (ns.model == "CSR")]
    rs = ns[(ns.archetype == r.archetype) & (ns.roi_id == r.roi_id) & (ns.model == "RSA")]
    pos = "sotto" if r.cv_loc < cs.cv_loc.min() else ("sopra" if r.cv_loc > cs.cv_loc.max() else "dentro")
    posg = "sotto" if r.cv < cs.cv.min() else ("sopra" if r.cv > cs.cv.max() else "dentro")
    cp1.append(dict(archetype=r.archetype, roi_id=r.roi_id, n_real=r.n, n_csr_min=cs.n.min(), n_csr_max=cs.n.max(),
                    n_rsa_min=rs.n.min(), n_rsa_max=rs.n.max(), nrep_csr=len(cs), nrep_rsa=len(rs),
                    nfail_csr=cs.n_failed.sum(), nfail_rsa=rs.n_failed.sum(), conserv_max=max(cs.conserv.max(), rs.conserv.max()),
                    cv_loc=r.cv_loc, pos_loc=pos, pos_glob=posg, relerr=abs(rs.cv_loc.median() - r.cv_loc) / r.cv_loc,
                    rsa_cvloc_min=rs.cv_loc.min(), rsa_cvloc_max=rs.cv_loc.max(),
                    pos_loc_rsa="sotto" if r.cv_loc < rs.cv_loc.min() else ("sopra" if r.cv_loc > rs.cv_loc.max() else "dentro"),
                    ecc_real=r.median_ecc_T, ecc_rsa=rs.median_ecc_T.median(),
                    frac_int_real=r.frac_interior, frac_int_csr=cs.frac_interior.median(), frac_int_rsa=rs.frac_interior.median()))
cp1 = pd.DataFrame(cp1); cp1.to_csv(f"{RV}/rv_cp1_roi.csv", index=False)
for A in PRIM:
    d = cp1[cp1.archetype == A]
    if A in ("A4", "A5"):
        k = (d.pos_loc == "sotto").sum(); c = cl("CP-1a", A); add("CP-1a", A, "ROI CV_loc sotto CSR", k, c[0], pf(k >= 4), c[1], "PASS",
                                                                   ",".join(d.pos_loc) + " | CV globale: " + ",".join(d.pos_glob))
    elif A == "A6":
        k = (d.pos_loc == "sopra").sum(); c = cl("CP-1a", A); add("CP-1a", A, "ROI CV_loc sopra CSR", k, c[0], pf(k >= 4), c[1], "PASS",
                                                                   ",".join(d.pos_loc) + " | CV globale: " + ",".join(d.pos_glob))
    else:
        c = cl("CP-1a", A); add("CP-1a", A, "posizione CV_loc", (d.pos_loc == "sopra").sum(), c[0], "registrato", c[1], "registrare", ",".join(d.pos_loc))
    k = (d.relerr <= .10).sum(); c = cl("CP-1b", A)
    if A in ("A4", "A5", "A6"): add("CP-1b", A, "ROI RSA entro 10%", k, c[0], pf(k >= 4), c[1], "PASS", " ".join(f"{v:.3f}" for v in d.relerr))
    else: add("CP-1b", A, "ROI RSA fuori 10%", 5 - k, c[0], pf(5 - k >= 3), c[1], "PASS", " ".join(f"{v:.3f}" for v in d.relerr))
# CP-2 from per-cell table (independent permutation RNG)
def adiff(a, b): d = np.abs(a - b) % np.pi; return np.minimum(d, np.pi - d)
rng = np.random.default_rng(12345); cp2 = []
for (A, roi), g in ci.groupby(["archetype", "roi_id"]):
    e = g[(g.ecc_N >= .8) & (g.ecc_T >= .5) & np.isfinite(g.theta_N)]
    dth = np.degrees(adiff(e.theta_T.values, e.theta_N.values)); obs = np.median(dth)
    perm = np.array([np.median(np.degrees(adiff(e.theta_T.values, rng.permutation(e.theta_N.values)))) for _ in range(2000)])
    cp2.append(dict(archetype=A, roi_id=roi, n_elig=len(e), med_dth=obs, p=(perm <= obs).mean(), p_plus1=((perm <= obs).sum() + 1) / 2001,
                    ecc_real_parquet=g.ecc_T.median()))
cp2 = pd.DataFrame(cp2).merge(cp1[["archetype", "roi_id", "ecc_real", "ecc_rsa"]]); cp2["d_ecc"] = cp2.ecc_real_parquet - cp2.ecc_rsa
cp2.to_csv(f"{RV}/rv_cp2_roi.csv", index=False)
d6 = cp2[cp2.archetype == "A6"]; k = (d6.d_ecc >= .05).sum(); c = cl("CP-2a", "A6"); add("CP-2a", "A6", "ROI d_ecc>=0.05", k, c[0], pf(k >= 4), c[1], "PASS")
for A in PRIM:
    d = cp2[cp2.archetype == A]; k = ((d.med_dth <= 35) & (d.p < .01)).sum(); c = cl("CP-2b", A)
    add("CP-2b", A, "ROI dtheta<=35 & p<0.01", k, c[0], pf(k >= 4) if A in ("A1", "A6") else "registrato", c[1],
        "PASS" if A in ("A1", "A6") else "registrare", " ".join(f"{m:.1f}/{p:.3f}/n{n}" for m, p, n in zip(d.med_dth, d.p, d.n_elig)))
d4 = cp2[cp2.archetype == "A4"]; k = (d4.d_ecc.abs() < .05).sum(); c = cl("CP-2c", "A4"); add("CP-2c", "A4", "ROI |d_ecc|<0.05", k, c[0], pf(k >= 4), c[1], "PASS")
# CP-3
cp3 = []
for A in PRIM:
    d = ci[ci.archetype == A]; v = (d.frac_out > .05).mean(); c = cl("CP-3", A)
    thr = {"A1": pf(v >= .10), "A5": pf(v <= .05)}.get(A, "registrato")
    add("CP-3", A, "frazione tagliati >5%", v, c[0], thr, c[1], {"A1": "PASS", "A5": "PASS"}.get(A, "registrare"),
        f"NA frac_out={d.frac_out.isna().sum()}; >10%: {(d.frac_out > .10).mean():.3f}; >20%: {(d.frac_out > .20).mean():.3f}")
    cp3.append(dict(archetype=A, cut5=v, cut10=(d.frac_out > .10).mean(), cut20=(d.frac_out > .20).mean(), roi_min=d.groupby('roi_id').apply(lambda z: (z.frac_out > .05).mean()).min(),
                    roi_max=d.groupby('roi_id').apply(lambda z: (z.frac_out > .05).mean()).max()))
pd.DataFrame(cp3).to_csv(f"{RV}/rv_cp3.csv", index=False)
# CP-4 from calib parquets: pooled CV (text) + per-window CV (sensitivity)
cc = pd.concat([pd.read_parquet(f) for f in glob.glob(f"{OUT}/calib/*.parquet")]); cc = cc[cc.interior]
cc["archetype"] = cc.win_id.str[:2]; cc["kind"] = np.where(cc.source == "manual_Luca", "manual", "seg")
cv = lambda v: np.std(v, ddof=1) / np.mean(v)
cp4 = []; npool = 0
for A in PRIM:
    m = cc[(cc.archetype == A) & (cc.kind == "manual")]; s = cc[(cc.archetype == A) & (cc.kind == "seg")]
    dif = cv(s.area) - cv(m.area); npool += dif > 0
    pw = [(cv(s.area[s.win_id == w]) - cv(m.area[m.win_id == w])) for w in sorted(set(m.win_id) & set(s.win_id))]
    # normalised by window mean (removes between-window density differences)
    norm = lambda z: z.groupby("win_id").area.transform(lambda a: a / a.mean())
    dn = cv(norm(s)) - cv(norm(m))
    cp4.append(dict(archetype=A, n_m=len(m), n_s=len(s), cv_m=cv(m.area), cv_s=cv(s.area), diff_pool=dif, diff_win_median=np.median(pw),
                    n_win_pos=int(np.sum(np.array(pw) > 0)), n_win=len(pw), diff_pool_winnorm=dn))
    c = cl("CP-4", A); add("CP-4", A, "CV_seg-CV_man pool", dif, c[0], pf(dif > 0), c[1], "registrare",
                           f"per finestra mediana {np.median(pw):.3f} ({int(np.sum(np.array(pw) > 0))}/{len(pw)} >0); pool normalizzato per finestra {dn:.3f}")
cp4 = pd.DataFrame(cp4); cp4.to_csv(f"{RV}/rv_cp4.csv", index=False)
c = cl("CP-4", "tutti"); add("CP-4", "tutti", "archetipi CV_seg>CV_man", npool, c[0], pf(npool >= 5), c[1], "PASS",
                             f"per-finestra: {(cp4.diff_win_median > 0).sum()}/6; normalizzato: {(cp4.diff_pool_winnorm > 0).sum()}/6")
out = pd.DataFrame(rows)
out["value_match"] = [ (abs(a - b) <= 1e-4 * max(1, abs(b))) if pd.notna(b) else False for a, b in zip(out.value_recomputed.astype(float), out.value_claimed.astype(float))]
out["assert_match"] = out.assertion_recomputed == out.assertion_claimed
out.to_csv(f"{RV}/rv_verdicts_recomputed.csv", index=False)
print(out[["id", "archetype", "value_recomputed", "value_claimed", "assertion_recomputed", "assertion_claimed", "value_match", "assert_match"]].to_string())
print(cp1[["archetype","roi_id","n_real","n_csr_min","n_csr_max","n_rsa_min","n_rsa_max","nrep_csr","nrep_rsa","nfail_csr","nfail_rsa","conserv_max"]].to_string())
print(cp4.to_string())

