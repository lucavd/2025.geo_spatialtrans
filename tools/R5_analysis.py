"""tools/R5_analysis.py — R5, tabelle di sintesi. Uso:
  edges  <outdir> : stime per sezione dalle tessere (P-1...P-4, R-1...R-4, CP-3, S-0...S-4, D-2, D-3), kernels.json per il rapporto, params.json per i lumi
  nc     <outdir> : P-5, CP-1, S-5 dalle ROI
  lumens <outdir> : CP-2
Tabelle in <outdir>/tables/. Bootstrap per blocchi sulle unita' (tessere, ROI, tessere con lumi), B = 2000, seme 20261031.
Pre-registrazione: docs/steps/R5.md (c5c9fce, ceca77b)."""
import sys, os, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R5_common as C

MODE, OUT = sys.argv[1], sys.argv[2]; TB = f"{OUT}/tables"; os.makedirs(TB, exist_ok=True)
KEYS = os.environ.get("R5_KEYS", "A1 A23 A4 A5 A6").split()
ARCH_OF = {"A1": ["A1"], "A23": ["A2", "A3"], "A4": ["A4"], "A5": ["A5"], "A6": ["A6"]}
TISSUE = {"A1": "intestino", "A23": "tumore del colon e stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
rng_for = lambda key, q: np.random.default_rng([C.SEED, KEYS.index(key) if key in KEYS else 9, q])

def boot_stat(n, stat, B=C.B, rng=None):
    idx = rng.integers(0, n, size=(B, n))
    reps = np.array([stat(ix) for ix in idx])
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps

def load_tiles(key):
    rows, errs = [], 0
    for f in sorted(glob.glob(f"{OUT}/tiles/{key}/t*.json")):
        d = json.load(open(f))
        if "error" in d: errs += 1; continue
        r = {"key": key, "tile": d["tile"], "b_far": d["b_far"]}
        for m in ("K1", "K2", "K3", "S0_bfixed", "S1_altmask", "S2_parity_glass", "S3_glass", "S4_200um", "CP3_shift", "SW_w05"):
            if d.get(m) is None: continue
            for k in ("sigma", "sigma_r", "sigma_c", "w", "ell", "dr", "dc", "a", "b", "nll", "converged", "n_bins", "nit"):
                r[f"{m}_{k}"] = d[m][k]
        I = np.array(d["fisher_delta_K2"]); r.update({"I_rr": I[0, 0], "I_rc": I[0, 1], "I_cc": I[1, 1]})
        z = np.load(f.replace(".json", ".npz"))               # tessuto H&E nei bin fuori dall'array (rilievo F14)
        occ = C.bin_occupancy(z["f_sub"].astype(float)); r["tissue_outside_array"] = bool(((~z["inarr_bin"]) & (occ > 0)).any())
        rows.append(r)          # senza tempi di calcolo: le tabelle devono essere identiche fra esecuzioni (G-4)
    return pd.DataFrame(rows), errs

def tail_beyond(r, sigma, w, ell):
    return (1 - w) * np.exp(-r ** 2 / (2 * sigma ** 2)) + w * (1 + r / ell) * np.exp(-r / ell)

if MODE == "edges":
    allT, out, sens, cp3, desc, kern, params = [], [], [], [], [], {}, {}
    per_key_sigma_reps = {}
    for key in KEYS:
        T, errs = load_tiles(key); allT.append(T); n = len(T)
        bg = json.load(open(f"{OUT}/tiles/{key}_background.json"))
        med = lambda col: (lambda ix: np.median(T[col].to_numpy()[ix]))
        # P-1 sigma
        e = float(T.K2_sigma.median()); lo, hi, reps = boot_stat(n, med("K2_sigma"), rng=rng_for(key, 1)); per_key_sigma_reps[key] = reps
        out.append({"key": key, "prediction": "P-1", "quantity": "sigma K2 (um)", "estimate": e, "lo": lo, "hi": hi, "interval": "[1, 4]",
                    "verdict": C.verdict(lo, hi, 1.0, 4.0), "n_tiles": n, "n_errors": errs})
        # P-2 w (previsione: sopra 0.01; [0, 0.01] = trascurabile)
        e = float(T.K2_w.median()); lo, hi, _ = boot_stat(n, med("K2_w"), rng=rng_for(key, 2))
        v = "confermata" if lo > 0.01 else ("smentita" if hi <= 0.01 else "indeterminata")
        out.append({"key": key, "prediction": "P-2", "quantity": "tail mass w K2", "estimate": e, "lo": lo, "hi": hi, "interval": "[0, 0.01] trascurabile",
                    "verdict": v, "n_tiles": n, "n_errors": errs})
        # P-3 anisotropia
        ratio = (T.K3_sigma_r / T.K3_sigma_c).to_numpy()
        e = float(np.median(ratio)); lo, hi, _ = boot_stat(n, lambda ix: np.median(ratio[ix]), rng=rng_for(key, 3))
        out.append({"key": key, "prediction": "P-3", "quantity": "sigma_row/sigma_col K3", "estimate": e, "lo": lo, "hi": hi, "interval": "[0.8, 1.25]",
                    "verdict": C.verdict(lo, hi, 0.8, 1.25), "n_tiles": n, "n_errors": errs})
        # P-4 scarto: modulo del vettore delle mediane
        dr, dc = T.K2_dr.to_numpy(), T.K2_dc.to_numpy()
        nrm = lambda ix: float(np.hypot(np.median(dr[ix]), np.median(dc[ix])))
        e = nrm(np.arange(n)); lo, hi, _ = boot_stat(n, nrm, rng=rng_for(key, 4))
        v = C.verdict(lo, hi, 0.0, 1.0) if key != "A5" else "senza attesa (" + C.verdict(lo, hi, 0.0, 1.0) + ")"
        out.append({"key": key, "prediction": "P-4", "quantity": "|delta| (um), median vector", "estimate": e, "lo": lo, "hi": hi, "interval": "[0, 1]",
                    "verdict": v, "n_tiles": n, "n_errors": errs, "dr_median": float(np.median(dr)), "dc_median": float(np.median(dc))})
        # scarto pesato con l'informazione di Fisher (sensibilita')
        I = np.zeros((2, 2)); v_ = np.zeros(2)
        for r in T.itertuples():
            Ii = np.array([[r.I_rr, r.I_rc], [r.I_rc, r.I_cc]])
            if np.all(np.linalg.eigvalsh(Ii) > 0): I += Ii; v_ += Ii @ np.array([r.K2_dr, r.K2_dc])
        dF = np.linalg.solve(I, v_)
        # sensibilita'
        for m, lab in (("K2", "primary (as pre-registered: b free, rho from all tissue bins)"), ("S0_bfixed", "S-0 b fixed from far glass"), ("S1_altmask", "S-1 alternative H&E mask"),
                       ("S2_parity_glass", "S-2 parity from glass"), ("S3_glass", "S-3 glass-only likelihood"), ("S4_200um", "S-4 central 200 um"), ("SW_w05", "S-W w <= 0.5 (bound of the first version)"), ("K1", "K1 Gaussian only")):
            if f"{m}_sigma" not in T: continue
            sub = T.dropna(subset=[f"{m}_sigma"])
            sens.append({"key": key, "variant": lab, "n_tiles": len(sub), "sigma": sub[f"{m}_sigma"].median(), "w": sub[f"{m}_w"].median() if f"{m}_w" in sub else 0.0,
                         "ell": sub[f"{m}_ell"].median(), "dr": sub[f"{m}_dr"].median(), "dc": sub[f"{m}_dc"].median(), "a": sub[f"{m}_a"].median(), "b": sub[f"{m}_b"].median()})
        sens.append({"key": key, "variant": "Fisher-weighted delta (K2)", "n_tiles": n, "dr": dF[0], "dc": dF[1]})
        # robustezza (dopo la revisione del codice, rilievi F06, F14, F16, F21)
        good = T.K2_converged.astype(bool) & (T.K2_w < C.BOUNDS["w"][1] - 0.001) & (T.K2_dr.abs() < 5.99) & (T.K2_dc.abs() < 5.99) & (T.K2_ell > 5.01) & (T.K2_sigma > 0.21) & (T.K2_sigma < 9.99)
        sub = T[good]
        sens.append({"key": key, "variant": "K2 converged and not at bounds", "n_tiles": int(good.sum()), "sigma": sub.K2_sigma.median(), "w": sub.K2_w.median(),
                     "ell": sub.K2_ell.median(), "dr": sub.K2_dr.median(), "dc": sub.K2_dc.median(), "a": sub.K2_a.median(), "b": sub.K2_b.median()})
        sub = T[~T.tissue_outside_array.astype(bool)]
        sens.append({"key": key, "variant": "tiles without tissue outside the array", "n_tiles": len(sub), "sigma": sub.K2_sigma.median(), "w": sub.K2_w.median(),
                     "ell": sub.K2_ell.median(), "dr": sub.K2_dr.median(), "dc": sub.K2_dc.median(), "a": sub.K2_a.median(), "b": sub.K2_b.median()})
        k3ok = (T.K3_nit > 0) & (T.K3_sigma_r > 0.21) & (T.K3_sigma_c > 0.21) & (T.K3_sigma_r < 9.99) & (T.K3_sigma_c < 9.99)
        rr3 = (T.K3_sigma_r / T.K3_sigma_c)[k3ok].to_numpy()
        lo3, hi3, _ = boot_stat(len(rr3), lambda ix: np.median(rr3[ix]), rng=rng_for(key, 8)) if len(rr3) > 1 else (np.nan, np.nan, None)
        sens.append({"key": key, "variant": "P-3 ratio, K3 moved from start and not at bounds", "n_tiles": int(k3ok.sum()), "ratio_K3": float(np.median(rr3)) if len(rr3) else np.nan,
                     "ratio_lo": lo3, "ratio_hi": hi3})
        nrm_t = np.hypot(dr, dc); lo4, hi4, _ = boot_stat(n, lambda ix: np.median(nrm_t[ix]), rng=rng_for(key, 9))
        sens.append({"key": key, "variant": "P-4 alternative: median of per-tile |delta|", "n_tiles": n, "abs_delta": float(np.median(nrm_t)), "abs_delta_lo": lo4, "abs_delta_hi": hi4})
        sens.append({"key": key, "variant": "fits at bounds or not converged (counts)", "n_tiles": n, "n_K2_w_at_bound": int((T.K2_w >= C.BOUNDS["w"][1] - 0.001).sum()), "n_K2_w_above_0.5": int((T.K2_w > 0.5).sum()), "n_K2_delta_at_bound": int(((T.K2_dr.abs() >= 5.99) | (T.K2_dc.abs() >= 5.99)).sum()), "n_K2_ell_at_lower_bound": int((T.K2_ell <= 5.01).sum()),
                     "n_K2_sigma_at_bound": int(((T.K2_sigma <= 0.21) | (T.K2_sigma >= 9.99)).sum()), "n_K2_not_converged": int((~T.K2_converged.astype(bool)).sum()),
                     "n_K3_nit0": int((T.K3_nit == 0).sum()), "n_tissue_outside_array": int(T.tissue_outside_array.astype(bool).sum())})
        # CP-3
        ddr = (T.CP3_shift_dr - T.K2_dr).to_numpy(); ddc = (T.CP3_shift_dc - T.K2_dc).to_numpy()
        lr, hr, _ = boot_stat(n, lambda ix: np.median(ddr[ix]), rng=rng_for(key, 5)); lc, hc, _ = boot_stat(n, lambda ix: np.median(ddc[ix]), rng=rng_for(key, 6))
        s0, s1 = T.K2_sigma.median(), T.CP3_shift_sigma.median(); w0, w1 = T.K2_w.median(), T.CP3_shift_w.median()
        ok = abs(np.median(ddr) + 2) <= 0.25 and abs(np.median(ddc) - 2) <= 0.25 and abs(s1 - s0) <= 0.15 + 0.05 * s0 and abs(w1 - w0) <= 0.01
        cp3.append({"key": key, "d_dr_median": float(np.median(ddr)), "d_dr_lo": lr, "d_dr_hi": hr, "d_dc_median": float(np.median(ddc)), "d_dc_lo": lc, "d_dc_hi": hc,
                    "sigma_K2": s0, "sigma_CP3": s1, "w_K2": w0, "w_CP3": w1, "pass": bool(ok)})
        # D-2, D-3 e nuclei per le tappe successive
        s, w, l = float(T.K2_sigma.median()), float(T.K2_w.median()), float(T.K2_ell.median())
        d = {"key": key, "sigma": s, "w": w, "ell": l, "b_tiles_median": float(T.K2_b.median()), "b_far": bg["b_far"], "b_far_rule": bg["rule"], "a_median": float(T.K2_a.median()),
             "mean_displacement_um": (1 - w) * s * np.sqrt(np.pi / 2) + w * 2 * l}
        for r in (2, 4, 8, 16, 32):
            d[f"frac_beyond_{r}um"] = float(tail_beyond(r, s, w, l))
        d["delta_spread_iqr_dr"] = float(np.subtract(*np.percentile(dr, [75, 25]))); d["delta_spread_iqr_dc"] = float(np.subtract(*np.percentile(dc, [75, 25])))
        desc.append(d)
        k2 = {"sigma_r": s, "sigma_c": s, "w": w, "ell": l, "dr": float(np.median(dr)), "dc": float(np.median(dc))}
        k1 = {"sigma_r": float(T.K1_sigma.median()), "sigma_c": float(T.K1_sigma.median()), "w": 0.0, "ell": 30.0, "dr": float(T.K1_dr.median()), "dc": float(T.K1_dc.median())}
        rg = rng_for(key, 7); draws = []
        for j in range(20):
            ix = rg.integers(0, n, n)
            draws.append({"sigma_r": float(np.median(T.K2_sigma.to_numpy()[ix])), "sigma_c": float(np.median(T.K2_sigma.to_numpy()[ix])),
                          "w": float(np.median(T.K2_w.to_numpy()[ix])), "ell": float(np.median(T.K2_ell.to_numpy()[ix])),
                          "dr": float(np.median(dr[ix])), "dc": float(np.median(dc[ix]))})
        kern[key] = {"K2": k2, "K1": k1, "draws": draws}
        params[key] = {"sigma": s, "w": w, "ell": l, "dr": k2["dr"], "dc": k2["dc"], "a": float(T.K2_a.median()), "b": float(T.K2_b.median())}
    # R-2: rapporto fra il sigma piu' grande e il piu' piccolo
    R = np.vstack([per_key_sigma_reps[k] for k in KEYS]); ests = [float(pd.DataFrame(allT[i]).K2_sigma.median()) for i in range(len(KEYS))]
    ratio_reps = R.max(axis=0) / R.min(axis=0); lo, hi = np.percentile(ratio_reps, [2.5, 97.5])
    rule2 = "un solo sigma" if hi < 1.25 else ("sigma per tessuto" if lo > 1.25 else "un solo sigma come predefinito, per tessuto come opzione")
    out.append({"key": "all", "prediction": "R-2", "quantity": "max/min sigma", "estimate": max(ests) / min(ests), "lo": lo, "hi": hi, "interval": "1.25", "verdict": rule2})
    pd.concat(allT).to_csv(f"{TB}/R5_tiles.csv", index=False)
    pd.DataFrame(out).to_csv(f"{TB}/R5_predictions_edges.csv", index=False)
    pd.DataFrame(sens).to_csv(f"{TB}/R5_sensitivity_edges.csv", index=False)
    pd.DataFrame(cp3).to_csv(f"{TB}/R5_CP3.csv", index=False)
    pd.DataFrame(desc).to_csv(f"{TB}/R5_kernel_descriptive.csv", index=False)
    json.dump(kern, open(f"{OUT}/kernels.json", "w"), indent=1); json.dump(params, open(f"{OUT}/params.json", "w"), indent=1)
    print(pd.DataFrame(out).round(4).to_string())

elif MODE == "nc":
    rows, fits, nerr = {}, [], 0
    for f in sorted(glob.glob(f"{OUT}/nc/A?_*__*.json")):
        d = json.load(open(f))
        if "error" in d: nerr += 1; print("ERROR", d["tag"], d["kernel"]); continue
        rows[d["tag"]] = {k: d[k] for k in ("tag", "archetype", "roi_id", "n_bins", "n_cells", "obs_star", "obs_star_far", "G3c_identical_to_R4")}
        fits.append({k: d[k] for k in ("tag", "archetype", "obs_star", "obs_star_far", "kernel", "r_true", "pred_star_far", "sigma_r", "sigma_c", "w", "ell", "dr", "dc")})
    Rr, F = pd.DataFrame(list(rows.values())), pd.DataFrame(fits)
    print("nc errors", nerr)
    r4 = pd.read_csv("results/R4/R4_PH1_roi.csv"); r4["tag"] = r4.archetype + "_" + r4.roi_id
    Rr = Rr.merge(r4[["tag", "cellweighted_star"]], on="tag", how="left"); Rr["obs_equals_R4"] = np.isclose(Rr.obs_star, Rr.cellweighted_star, rtol=1e-12, atol=0)
    Rr.to_csv(f"{TB}/R5_nc_roi.csv", index=False); F.to_csv(f"{TB}/R5_nc_fits.csv", index=False)
    expected = {"A1": "uniforme", "A2": "due densita (sotto)", "A3": "due densita (sotto)", "A4": "due densita (sotto)", "A5": "due densita (sopra)", "A6": "due densita (sopra)"}
    out = []
    for a in ["A1", "A2", "A3", "A4", "A5", "A6"]:
        Fa = F[F.archetype == a]; tags = sorted(Fa.tag.unique())
        P = Fa[Fa.kernel == "K2"].set_index("tag").reindex(tags)
        Dr = Fa[Fa.kernel.str.startswith("draw")].pivot(index="tag", columns="kernel", values="r_true").reindex(tags)
        Dp = Fa[Fa.kernel.str.startswith("draw")].pivot(index="tag", columns="kernel", values="pred_star_far").reindex(tags)
        obs = P.obs_star.to_numpy(); obsf = P.obs_star_far.to_numpy()
        rg = np.random.default_rng([C.SEED, 600 + ["A1", "A2", "A3", "A4", "A5", "A6"].index(a)])
        n, J = len(tags), Dr.shape[1]
        ix = rg.integers(0, n, size=(C.B, n)); jj = rg.integers(0, J, size=C.B)
        Drn, Dpn = Dr.to_numpy(), Dp.to_numpy()
        with np.errstate(all="ignore"):
            rr = np.array([np.nanmean(Drn[ix[b], jj[b]]) for b in range(C.B)])
            cpr = np.array([np.nanmean(Dpn[ix[b], jj[b]]) / np.mean(obsf[ix[b]]) for b in range(C.B)])
        n_fail = int(P.r_true.isna().sum()); frac_fail_draws = float(np.mean(~np.isfinite(Drn)))
        r_obs = float(np.mean(obs)); rK1 = Fa[Fa.kernel == "K1"].set_index("tag").reindex(tags).r_true
        cp1_name = "CP-1" if a in ("A2", "A4", "A5", "A6") else "CP-1 (descrittiva)"
        # revisione numerica F-8: (1) come scritto = media di TUTTE le ROI (e delle repliche): se una ROI o una replica non si inverte la
        # stima o l'intervallo non esistono ("non calcolabile come scritto"); (2) sensibilita' dichiarata = sole ROI riuscite (nanmean)
        with np.errstate(all="ignore"):
            rr_all = np.array([np.mean(Drn[ix[b], jj[b]]) for b in range(C.B)])
            cp_all = np.array([np.mean(Dpn[ix[b], jj[b]]) / np.mean(obsf[ix[b]]) for b in range(C.B)])
        for variant, est, reps, est_cp, reps_cp in (
                ("as written (all ROIs)", float(np.mean(P.r_true.to_numpy())), rr_all, float(np.mean(P.pred_star_far.to_numpy()) / np.mean(obsf)), cp_all),
                ("sensitivity: successful ROIs only (nanmean)", float(np.nanmean(P.r_true)) if n_fail < n else np.nan, rr,
                 float(np.nanmean(P.pred_star_far) / np.mean(obsf[P.pred_star_far.notna().to_numpy()])) if n_fail < n else np.nan, cpr)):
            fin = np.isfinite(reps); fin_cp = np.isfinite(reps_cp)
            ok_int = fin.all() if variant.startswith("as written") else fin.any()
            lo, hi = (np.nanpercentile(reps, [2.5, 97.5]) if ok_int else (np.nan, np.nan))
            if not np.isfinite(est) or not ok_int:
                got, v, farther = "non calcolabile", "non calcolabile come scritto" if variant.startswith("as written") else "non calcolabile", np.nan
            else:
                cls = C.verdict(lo, hi, 0.909, 1.1); side = "sotto" if est < 1 else "sopra"
                got = {"confermata": "uniforme", "smentita": f"due densita ({side})", "indeterminata": "a cavallo"}[cls]
                farther = bool(abs(np.log(est)) > abs(np.log(r_obs)))
                v = "confermata" if got == expected[a] else ("indeterminata" if got == "a cavallo" else "smentita")
                if v == "confermata" and expected[a] != "uniforme" and not farther:
                    v = "smentita"     # pre-registrato: due densita' con r_vero piu' lontano da 1 dell'osservato (rilievo F09)
            out.append({"archetype": a, "prediction": "P-5", "variant": variant, "r_obs_mean": r_obs, "estimate": est, "lo": lo, "hi": hi,
                        "median_reps": float(np.nanmedian(reps)) if fin.any() else np.nan, "frac_nan_reps": float(np.mean(~fin)),
                        "class": got, "expected": expected[a], "verdict": v, "farther_from_1": farther, "n_roi": n, "n_roi_failed": n_fail,
                        "frac_failed_draws": frac_fail_draws, "n_draws": J, "r_true_K1_mean": float(np.mean(rK1)) if rK1.notna().all() else np.nan})
            ok_cp = fin_cp.all() if variant.startswith("as written") else fin_cp.any()
            lo2, hi2 = (np.nanpercentile(reps_cp, [2.5, 97.5]) if ok_cp else (np.nan, np.nan))
            v2 = C.verdict(lo2, hi2, 0.909, 1.1) if (np.isfinite(est_cp) and ok_cp) else ("non calcolabile come scritto" if variant.startswith("as written") else "non calcolabile")
            out.append({"archetype": a, "prediction": cp1_name, "variant": variant, "r_obs_mean": float(np.mean(obsf)), "estimate": est_cp, "lo": lo2, "hi": hi2,
                        "median_reps": float(np.nanmedian(reps_cp)) if fin_cp.any() else np.nan, "frac_nan_reps": float(np.mean(~fin_cp)),
                        "class": "pred/obs", "expected": "[0.909, 1.1]", "verdict": v2, "n_roi": n, "n_roi_failed": n_fail, "n_draws": J})
    pd.DataFrame(out).to_csv(f"{TB}/R5_predictions_nc.csv", index=False)
    print(pd.DataFrame(out).round(4).to_string())

elif MODE == "lumens":
    out = []
    for key in KEYS:
        f = f"{OUT}/lumens/{key}_lumens.csv"
        if not os.path.exists(f) or os.path.getsize(f) < 5: continue
        L = pd.read_csv(f)
        if L.empty: continue
        g = L.groupby("tile").sum(numeric_only=True); n = len(g)
        for lo_, hi_ in ((2, 4), (4, 8), (8, 16), (16, 32)):
            o, p = g[f"obs_{lo_}_{hi_}"].to_numpy(), g[f"pred_{lo_}_{hi_}"].to_numpy()
            if p.sum() <= 0: continue
            e = o.sum() / p.sum(); lo, hi, _ = boot_stat(n, lambda ix: o[ix].sum() / max(p[ix].sum(), 1e-12), rng=rng_for(key, 40 + lo_))
            out.append({"key": key, "band_um": f"{lo_}-{hi_}", "n_lumens": len(L), "n_tiles": n, "n_bins": int(g[f"n_{lo_}_{hi_}"].sum()),
                        "obs_per_bin": o.sum() / max(g[f"n_{lo_}_{hi_}"].sum(), 1), "pred_per_bin": p.sum() / max(g[f"n_{lo_}_{hi_}"].sum(), 1),
                        "ratio": e, "lo": lo, "hi": hi, "verdict": C.verdict(lo, hi, 0.8, 1.25)})
    pd.DataFrame(out).to_csv(f"{TB}/R5_CP2_lumens.csv", index=False)
    print(pd.DataFrame(out).round(3).to_string())
