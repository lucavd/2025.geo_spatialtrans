"""tools/S2.4_analysis.py — S2.4: esiti di P-1...P-4, P-6, P-7, CP-1, CP-2 e misure (docs/steps/S2.4.md).
Statistiche dei simulati con le stesse funzioni dei veri (tools/S2.4_explore.py job_sim) in <PD>/stats_<variante>/sim/.
P-1...P-3: differenza simulato (media dei due semi) - vero per area, media sulle 5 aree +- 2 errori standard (sd / sqrt(n));
P-4, P-6: rapporto simulato / vero dei valori aggregati sulle aree del tessuto, errore standard con bootstrap delle aree appaiate
(B = 2000, seme 20261115); P-7: tempo e memoria per area (media dei due semi), media +- 2 errori standard.
Uso: .venv/bin/python tools/S2.4_analysis.py stats <PD> | analyse <PD>"""
import sys, os, json, glob, subprocess
import numpy as np, pandas as pd
B = 2000; SEED = 20261115
TIS = {"A1": "intestino", "A2": "tumore", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
def verdict(e, se, lo, hi):
    if not np.isfinite(e) or not np.isfinite(se): return "non valutabile"
    a, b = e - 2 * se, e + 2 * se
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
def load_stats(PD, VAR):
    fs = sorted(glob.glob(f"{PD}/stats_{VAR}/sim/*.json"))
    S = pd.DataFrame([json.load(open(f)) for f in fs]); S["tag"] = S["tag"].astype(str); return S
def pooled_ratio(sim, real, num, den, rng):
    """rapporto (sum num / sum den)_sim / (sum num / sum den)_real con bootstrap delle aree appaiate."""
    tags = sorted(set(sim.tag) & set(real.tag)); s = sim.set_index("tag").loc[tags]; r = real.set_index("tag").loc[tags]
    def val(ix):
        ss, rr = s.iloc[ix], r.iloc[ix]
        a = ss[num].sum() / ss[den].sum() if ss[den].sum() > 0 else np.nan; b = rr[num].sum() / rr[den].sum() if rr[den].sum() > 0 and rr[num].sum() > 0 else np.nan
        return a / b if np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan
    e = val(np.arange(len(tags))); bs = np.array([val(rng.integers(0, len(tags), len(tags))) for _ in range(B)])
    return e, float(np.nanstd(bs, ddof=1)) if np.isfinite(bs).sum() > 10 else np.nan, int(len(tags))
if __name__ == "__main__":
    mode, PD = sys.argv[1], sys.argv[2]
    if mode == "stats":
        for V in ("main", "stripes_off", "tail2_off"):
            subprocess.run([".venv/bin/python", "-B", "tools/S2.4_explore.py", "sim", f"{PD}/stats_{V}", PD, V], check=True)
        sys.exit(0)
    E12 = pd.read_csv("results/S2.4/S2.4_E12_summary.csv").set_index("A")
    R = pd.read_csv("results/S2.4/S2.4_E1_real.csv"); R["tag"] = R["tag"].astype(str)
    Sm = load_stats(PD, "main"); So = load_stats(PD, "stripes_off"); St = load_stats(PD, "tail2_off")
    num = Sm.select_dtypes("number").columns
    Sa = Sm.groupby(["tag", "archetype"])[list(num)].mean().reset_index()
    rng = np.random.default_rng(SEED); rows = []
    for A in sorted(TIS):
        r = R[R.archetype == A].set_index("tag"); s = Sa[Sa.archetype == A].set_index("tag"); tags = sorted(set(r.index) & set(s.index))
        for pid, q in (("P-1", "odd_over_even"), ("P-2", "sd_repro_hp"), ("P-3", "lag1_repro_hp")):
            for ax in ("rows", "cols"):
                col = f"stripes_{ax}_{q}"; d = (s.loc[tags, col] - r.loc[tags, col]).to_numpy(float)
                hw = float(E12.loc[A, f"real_{col}_maxdiff"]); e = float(np.mean(d)); se = float(np.std(d, ddof=1) / np.sqrt(len(d)))
                rows.append({"id": pid, "A": A, "tessuto": TIS[A], "quantita": f"{q} {ax}", "sim": float(s.loc[tags, col].mean()), "real": float(r.loc[tags, col].mean()),
                             "stima": e, "se": se, "lo": -hw, "hi": hw, "n": len(d), "esito": verdict(e, se, -hw, hw)})
        sA = Sm[Sm.archetype == A]
        for b in ("8_16", "16_32", "32_64"):
            sim2 = sA.groupby("tag")[[f"lumen_band_{b}_sum", f"lumen_band_{b}_n", "lumen_ref_sum", "lumen_ref_n"]].sum().reset_index()
            sim2["bm"] = sim2[f"lumen_band_{b}_sum"]; sim2["rm"] = sim2["lumen_ref_sum"]
            sim2["bn"] = sim2[f"lumen_band_{b}_n"]; sim2["rn"] = sim2["lumen_ref_n"]
            rr = R[R.archetype == A].copy(); rr["bm"] = rr[f"lumen_band_{b}_sum"]; rr["rm"] = rr["lumen_ref_sum"]; rr["bn"] = rr[f"lumen_band_{b}_n"]; rr["rn"] = rr["lumen_ref_n"]
            for d_ in (sim2, rr): d_["band_mean_x_refn"] = d_["bm"] * d_["rn"]; d_["ref_x_bandn"] = d_["rm"] * d_["bn"]
            # rapporto di fascia = (sum band / sum band_n) / (sum ref / sum ref_n); per il bootstrap si usa la forma aggregata
            def val(ix, ss=sim2.reset_index(drop=True), rr_=rr.reset_index(drop=True)):
                a = ss.iloc[ix]; c_ = rr_.iloc[ix]
                fa = (a.bm.sum() / a.bn.sum()) / (a.rm.sum() / a.rn.sum()) if a.bn.sum() > 0 else np.nan
                fc = (c_.bm.sum() / c_.bn.sum()) / (c_.rm.sum() / c_.rn.sum()) if c_.bn.sum() > 0 and c_.bm.sum() > 0 else np.nan
                return fa / fc if np.isfinite(fa) and np.isfinite(fc) else np.nan
            sim2 = sim2.set_index("tag").loc[sorted(rr.tag)].reset_index(); rr = rr.sort_values("tag")
            n_ = len(rr); e = val(np.arange(n_), sim2, rr.reset_index(drop=True))
            bs = np.array([val(rng.integers(0, n_, n_), sim2, rr.reset_index(drop=True)) for _ in range(B)]); se = float(np.nanstd(bs, ddof=1)) if np.isfinite(bs).sum() > 10 else np.nan
            judged = A in ("A1", "A2", "A6")
            rows.append({"id": "P-4", "A": A, "tessuto": TIS[A], "quantita": f"fascia {b.replace('_', '-')} um", "sim": np.nan, "real": np.nan, "stima": e, "se": se,
                         "lo": 0.8, "hi": 1.25, "n": int((rr.bn > 0).sum()), "esito": verdict(e, se, 0.8, 1.25) if judged else "descrittivo"})
        e, se, n_ = pooled_ratio(Sa[Sa.archetype == A].assign(z=lambda d: d.frac_zero_tissue * d.n_tissue_bins), R[R.archetype == A].assign(z=lambda d: d.frac_zero_tissue * d.n_tissue_bins), "z", "n_tissue_bins", rng)
        rows.append({"id": "P-6", "A": A, "tessuto": TIS[A], "quantita": "quadratini vuoti in tessuto", "sim": float(Sa[Sa.archetype == A].frac_zero_tissue.mean()),
                     "real": float(R[R.archetype == A].frac_zero_tissue.mean()), "stima": e, "se": se, "lo": 0.8, "hi": 1.25, "n": n_, "esito": verdict(e, se, 0.8, 1.25)})
    # P-7 tempo e memoria
    J = pd.DataFrame([json.load(open(f)) for f in glob.glob(f"{PD}/main/*_k?.json")])
    Ja = J.groupby(["tag", "archetype"])[["elapsed_total_s", "peak_mb", "M", "elapsed_place_s", "elapsed_profiles_s", "elapsed_capture_s", "elapsed_stripes_s", "n_generated", "n_captured", "n_background", "n_lost"]].mean().reset_index()
    for A in sorted(TIS):
        j = Ja[Ja.archetype == A]
        for q, hi in (("elapsed_total_s", 300.0), ("peak_mb", 8192.0)):
            v = j[q].to_numpy(float); e = float(v.mean()); se = float(v.std(ddof=1) / np.sqrt(len(v)))
            rows.append({"id": "P-7", "A": A, "tessuto": TIS[A], "quantita": q, "sim": e, "real": np.nan, "stima": e, "se": se, "lo": 0.0, "hi": hi, "n": len(v), "esito": verdict(e, se, 0.0, hi)})
    P = pd.DataFrame(rows); P.to_csv("results/S2.4/S2.4_predictions.csv", index=False)
    # CP-1: strisce spente
    cp = []
    for A in sorted(TIS):
        o = So[So.archetype == A]
        sd = max(o.stripes_rows_sd_repro_hp.mean(), o.stripes_cols_sd_repro_hp.mean()); pr = max(abs(o.stripes_rows_odd_over_even.mean() - 1), abs(o.stripes_cols_odd_over_even.mean() - 1))
        cp.append({"id": "CP-1", "A": A, "sd_fine_max": sd, "parity_dev_max": pr, "ok": bool(sd < 0.05 and pr <= 0.01)})
    # CP-2: seconda coda spenta, tumore, stessi semi di main k = 0
    m0 = Sm[(Sm.k == 0) & (Sm.archetype == "A2")].set_index("tag"); t0 = St[St.archetype == "A2"].set_index("tag"); tags = sorted(set(m0.index) & set(t0.index))
    for b in ("16_32", "32_64"):
        def fv(d, ix): d = d.iloc[ix]; return (d[f"lumen_band_{b}_sum"].sum() / d[f"lumen_band_{b}_n"].sum()) / (d.lumen_ref_sum.sum() / d.lumen_ref_n.sum()) if d[f"lumen_band_{b}_n"].sum() > 0 else np.nan
        mm, tt = m0.loc[tags], t0.loc[tags]; e = fv(mm, np.arange(len(tags))) - fv(tt, np.arange(len(tags)))
        bs = []
        for _ in range(B): ix = rng.integers(0, len(tags), len(tags)); bs.append(fv(mm, ix) - fv(tt, ix))
        se = float(np.nanstd(bs, ddof=1)); cp.append({"id": "CP-2", "A": "A2", "band": b, "main_minus_off": e, "se": se, "main": fv(mm, np.arange(len(tags))), "off": fv(tt, np.arange(len(tags))), "ok": bool(e > 2 * se)})
    CP = pd.DataFrame(cp); CP.to_csv("results/S2.4/S2.4_CP12.csv", index=False)
    # misure
    S23 = pd.read_csv("results/S2.4/S2.4_E2_sim.csv"); S23a = S23.groupby("archetype").mean_umi_tissue.mean()
    Jo = pd.DataFrame([json.load(open(f)) for f in glob.glob(f"{PD}/stripes_off/*_k0.json")]).set_index("tag")
    J0 = J[J.k == 0].set_index("tag")
    Ms = []
    for A in sorted(TIS):
        tg = sorted(set(J0[J0.archetype == A].index) & set(Jo.index))
        Ms.append({"A": A, "M_mean": float(Ja[Ja.archetype == A].M.mean()), "M_min": float(J[J.archetype == A].M.min()), "M_max": float(J[J.archetype == A].M.max()),
                   "umi_tissue_S24_over_S23": float(Sa[Sa.archetype == A].mean_umi_tissue.mean() / S23a[A]),
                   "umi_tissue_S24_over_real": float(Sa[Sa.archetype == A].mean_umi_tissue.mean() / R[R.archetype == A].mean_umi_tissue.mean()),
                   "place_time_ratio_M_over_M1": float((J0.loc[tg].elapsed_place_s / Jo.loc[tg].elapsed_place_s).mean()),
                   "total_time_ratio_M_over_M1": float((J0.loc[tg].elapsed_total_s / Jo.loc[tg].elapsed_total_s).mean()),
                   "captured_over_target_R4": float((J[J.archetype == A].n_captured / J[J.archetype == A].n_target_R4).mean()),
                   "lost_frac": float((J[J.archetype == A].n_lost / J[J.archetype == A].n_generated).mean())})
    MS = pd.DataFrame(Ms); MS.to_csv("results/S2.4/S2.4_measures.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    print(P.round(4).to_string(index=False)); print(CP.round(4).to_string(index=False)); print(MS.round(3).to_string(index=False))
