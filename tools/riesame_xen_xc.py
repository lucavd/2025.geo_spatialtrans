"""tools/riesame_xen_xc.py - Riesame Parte B, bracci XC e XC2 (docs/riesame/REGOLA_B.md; scelte XD-13...XD-22 di docs/riesame/DECISIONI_B_XENIUM.md).
Trascritti Xenium di P2 (qv >= 20, 408 geni; cache di tools/riesame_xen_export.py) -> spostamento dello Step 2 con i parametri del preset del
tessuto (regola della zona di bordo e seconda coda, come place_and_move_transcripts() e simulate_visium_hd_pieces(); campo e r spenti perche'
le nascite sono quelle vere) -> strisce (simulate_stripes() del pacchetto, via tools/riesame_xen_params.R) con tetto M*, cattura binomiale
min(1, k c) con k unico per tessuto e braccio, fondo Poisson(b_p c). XC2: lunghezze dello spostamento (sigma_in, sigma_e, ell, ell2) raddoppiate.
Fasi: prep (maschere e semi) | move (arrivi per area, seme, braccio) | depth (k per tessuto) | capture (esportazioni).
Uso: .venv/bin/python -B tools/riesame_xen_xc.py <prep|move|depth|capture> [aree]"""
import sys, os, json, time, importlib.util, numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq, shapely
from scipy import ndimage, optimize
sys.path.insert(0, "tools")
_s = importlib.util.spec_from_file_location("xexp", "tools/riesame_xen_export.py"); XE = importlib.util.module_from_spec(_s); _s.loader.exec_module(XE)
M = XE.M; OUTB = XE.OUTB; CACHE = XE.CACHE; W = f"{OUTB}/_xc"; AREAS = XE.AREAS; KS = (0, 1); ARMS = ("XC", "XC2")
SEED_XC = 20261125; SEED_STRIPES = 20261126; STEP = 0.5
PRM = json.load(open("results/riesame/xen_preset_params.json"))
def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)
def regions(a):
    t = pd.read_csv(f"{M}/S2.3/layers/{a}_regions.tsv", sep="\t"); return shapely.union_all(shapely.from_wkb([bytes.fromhex(h) for h in t.wkb]))
def params(a, arm):
    tp = dict(PRM[a[:2]]["transcript"]); cp = PRM[a[:2]]["capture"]; p = {k: tp[k] for k in ("sigma_in", "sigma_e", "D", "w", "ell")}
    p["w2"] = cp["w2"]; p["ell2"] = cp["ell2"]
    if arm == "XC2":
        for k in ("sigma_in", "sigma_e", "ell", "ell2"): p[k] *= 2.0
    return p
def rng_for(a, k):
    ss = np.random.SeedSequence([SEED_XC, AREAS.index(a), k]); mv, cap, bg = ss.spawn(3)
    return [np.random.Generator(np.random.PCG64(s)) for s in (mv, cap, bg)]
def stripe_seed(a, k): return int(np.random.SeedSequence([SEED_STRIPES, AREAS.index(a), k]).generate_state(1)[0] % (2**31 - 1))
def fold(x, lo, hi): Wd = hi - lo; t = np.mod(x - lo, 2 * Wd); return lo + Wd - np.abs(t - Wd)

def stage_prep(areas):
    os.makedirs(W, exist_ok=True); rows = []
    for a in areas:
        b = pd.read_parquet(f"{XE.F}/{a}/bins.parquet"); R = regions(a); shapely.prepare(R)
        inr = shapely.contains_xy(R, b.x_um.to_numpy(), b.y_um.to_numpy())
        pq.write_table(pa.table({"bin_id": b.bin_id.to_numpy(), "row": b.row.to_numpy(), "col": b.col.to_numpy(), "in_regions": inr}), f"{W}/stripes_in_{a}.parquet")
        # maschera del tessuto sulla sottogriglia di 0.5 um nel riquadro [0, lato]; fuori dal riquadro tessuto (frame_open), distanza dal bordo
        s = XE.side(a); o = -XE.MARGIN_UM; n = int(np.ceil((s + 2 * XE.MARGIN_UM) / STEP)); cc = o + (np.arange(n) + 0.5) * STEP
        Xg, Yg = np.meshgrid(cc, cc); fr = (Xg >= 0) & (Xg <= s) & (Yg >= 0) & (Yg <= s)
        Mt = np.ones((n, n), bool); Mt[fr] = shapely.contains_xy(R, Xg[fr], Yg[fr])          # fuori dal riquadro: tessuto (frame_open)
        Dm = ndimage.distance_transform_edt(Mt) * STEP if not Mt.all() else np.full((n, n), np.inf)
        np.save(f"{W}/dist_{a}.npy", Dm.astype(np.float32))
        for k in KS: rows.append({"area": a, "A": a[:2], "k": k, "seed": stripe_seed(a, k)})
        log("prep", a, "in_regions", float(inr.mean()), "tissue subpx", float(Mt.mean()), "dist max", float(np.nanmax(Dm[np.isfinite(Dm)])) if np.isfinite(Dm).any() else None)
    sd = pd.DataFrame(rows); f = f"{W}/stripes_seeds.csv"
    if os.path.exists(f): sd = pd.concat([pd.read_csv(f), sd]).drop_duplicates(["area", "k"], keep="last")
    sd.to_csv(f, index=False)

def stage_move(areas):
    for a in areas:
        b = pd.read_parquet(f"{XE.F}/{a}/bins.parquet"); assign, _ = XE.grid_lut(b); s = XE.side(a)
        tx = pd.read_parquet(f"{CACHE}/tx_{a}.parquet"); x0 = tx.x_um.to_numpy(); y0 = tx.y_um.to_numpy(); gene = tx.gene.to_numpy().astype(np.int64)
        Dm = np.load(f"{W}/dist_{a}.npy"); n = Dm.shape[0]; o = -XE.MARGIN_UM
        ins = (x0 >= 0) & (x0 < s) & (y0 >= 0) & (y0 < s)
        pi = np.clip(np.floor((y0 - o) / STEP).astype(np.int64), 0, n - 1); qi = np.clip(np.floor((x0 - o) / STEP).astype(np.int64), 0, n - 1)
        d = Dm[pi, qi].astype(float)                               # distanza dal vetro (sottogriglia di 0.5 um sulla finestra, riquadro aperto)
        for k in KS:
            for arm in ARMS:
                fo = f"{W}/arr_{arm}_{a}_k{k}.parquet"
                if os.path.exists(fo): continue
                t0 = time.time(); p = params(a, arm); rmv = rng_for(a, k)[0]; nt = len(x0)
                u1 = rmv.random(nt); u2 = rmv.random(nt); tail = u1 < p["w"]; wide = ~tail & (u2 < np.exp(-d / p["D"]))
                dx = np.zeros(nt); dy = np.zeros(nt); nT = int(tail.sum())
                rho = rmv.gamma(2.0, p["ell"], nT); th = rmv.uniform(0, 2 * np.pi, nT); dx[tail] = rho * np.cos(th); dy[tail] = rho * np.sin(th)
                g = ~tail; sg = np.where(wide[g], p["sigma_e"], p["sigma_in"]); dx[g] = rmv.standard_normal(int(g.sum())) * sg; dy[g] = rmv.standard_normal(int(g.sum())) * sg
                x1 = x0 + dx; y1 = y0 + dy
                sec = rmv.random(nt) < p["w2"]; nS = int(sec.sum())
                rho2 = rmv.gamma(2.0, p["ell2"], nS); th2 = rmv.uniform(0, 2 * np.pi, nS)
                x1[sec] = fold(x0[sec] + rho2 * np.cos(th2), 0.0, s); y1[sec] = fold(y0[sec] + rho2 * np.sin(th2), 0.0, s)
                bi = assign(x1, y1); ok = bi >= 0
                A = pd.DataFrame({"bin": bi[ok], "gene": gene[ok]}).groupby(["bin", "gene"]).size().reset_index(name="n")
                pq.write_table(pa.Table.from_pandas(A.astype({"bin": np.int64, "gene": np.int64, "n": np.int32}), preserve_index=False), fo + ".tmp"); os.replace(fo + ".tmp", fo)
                st = {"area": a, "k": k, "arm": arm, "params": p, "n_births_window": nt, "n_births_in_box": int(ins.sum()), "n_arrived_box": int(ok.sum()),
                      "share_tail": float(tail.mean()), "share_wide": float(wide.mean()), "share_second_tail": float(sec.mean()),
                      "mean_disp_um_nonsecond": float(np.hypot(dx, dy)[~sec].mean()), "sec": round(time.time() - t0, 1)}
                json.dump(st, open(fo.replace(".parquet", ".json"), "w")); log(json.dumps({k_: v for k_, v in st.items() if k_ != "params"}))

def _load_ct(a, k):
    c = pd.read_parquet(f"{W}/stripes_{a}_k{k}.parquet"); b = pd.read_parquet(f"{XE.F}/{a}/bins.parquet")
    assert (c.bin_id.to_numpy() == b.bin_id.to_numpy()).all()
    return np.minimum(c.c.to_numpy(), PRM[a[:2]]["Mstar"]), b.in_reg.to_numpy()
def stage_depth(areas):
    tg = pd.read_csv("results/riesame/xen_targets.csv").set_index("area")
    s_share = {t: float(tg.loc[[a for a in AREAS if a.startswith(t)], ["s_panel_share_k0", "s_panel_share_k1"]].to_numpy().mean()) for t in ("A2", "A3")}
    out = {"s_panel_share": s_share, "b_panel": {t: PRM[t]["capture"]["b"] * s_share[t] for t in s_share}, "Mstar": {t: PRM[t]["Mstar"] for t in s_share}, "k": {}}
    for arm in ARMS:
        for t in ("A2", "A3"):
            aa = [a for a in AREAS if a.startswith(t)]; target = float(tg.loc[aa, "v_mean408_reg"].mean()); bp = out["b_panel"][t]; parts = []
            for a in aa:
                for k in KS:
                    ct, reg = _load_ct(a, k); A = pd.read_parquet(f"{W}/arr_{arm}_{a}_k{k}.parquet")
                    arr = np.bincount(A.bin.to_numpy(), weights=A.n.to_numpy(), minlength=len(ct)); parts.append((arr[reg], ct[reg]))
            def mean_exp(kk): return np.mean([(np.minimum(1.0, kk * c) * ar).sum() / len(ar) + bp * c.mean() for ar, c in parts])
            kmax = 1.0 / min(c.min() for _, c in parts)
            if mean_exp(kmax * 10) < target: out["k"][f"{arm}_{t}"] = {"k": None, "target": target, "max_reachable": float(mean_exp(kmax * 10))}; continue
            kk = optimize.brentq(lambda v: mean_exp(v) - target, 1e-9, kmax * 10, xtol=1e-14)
            trunc_bins = float(np.mean(np.concatenate([kk * c > 1 for _, c in parts])))
            lost = sum(((kk * c - np.minimum(1, kk * c)) * ar).sum() for ar, c in parts) / sum((kk * c * ar).sum() for ar, c in parts)
            out["k"][f"{arm}_{t}"] = {"k": kk, "k_times_Mstar": kk * PRM[t]["Mstar"], "target": target, "check_mean": float(mean_exp(kk)),
                                      "share_reg_bins_truncated": trunc_bins, "share_expected_transcripts_truncated": float(lost),
                                      "mean_arrived_per_reg_bin": float(np.mean([ar.mean() for ar, _ in parts]))}
            log(arm, t, out["k"][f"{arm}_{t}"])
    json.dump(out, open(f"{OUTB}/XC_depth.json", "w"), indent=1)

def stage_capture(areas):
    dp = json.load(open(f"{OUTB}/XC_depth.json"))
    for a in areas:
        t = a[:2]; bp = dp["b_panel"][t]
        tx = pd.read_parquet(f"{CACHE}/tx_{a}.parquet", columns=["gene"]); comp = np.bincount(tx.gene.to_numpy().astype(np.int64), minlength=408).astype(float); comp /= comp.sum()
        xd = f"{OUTB}/X/{a}"; b = pd.read_parquet(f"{xd}/bins.parquet"); cells = pd.read_parquet(f"{xd}/cells.parquet"); ex = json.load(open(f"{xd}/export.json"))
        for arm in ARMS:
            kk = dp["k"][f"{arm}_{t}"]["k"]
            if kk is None: log("skip", arm, a, "target not reachable"); continue
            for k in KS:
                out = f"{OUTB}/{arm}/{a}_k{k}"
                if os.path.exists(f"{out}/export.json"): continue
                t0 = time.time(); ct, reg = _load_ct(a, k); _, rcap, rbg = rng_for(a, k)
                A = pd.read_parquet(f"{W}/arr_{arm}_{a}_k{k}.parquet"); p = np.minimum(1.0, kk * ct[A.bin.to_numpy()])
                kept = rcap.binomial(A.n.to_numpy(), p); A = A.assign(n=kept)[kept > 0]
                K = rbg.poisson(bp * ct); nb = int(K.sum()); gb = rbg.choice(408, size=nb, p=comp); bb = np.repeat(np.arange(len(ct)), K)
                X = pd.concat([A[["bin", "gene", "n"]], pd.DataFrame({"bin": bb, "gene": gb, "n": 1})]).groupby(["bin", "gene"], as_index=False).n.sum()
                mv = json.load(open(f"{W}/arr_{arm}_{a}_k{k}.json"))
                info = {kx: v for kx, v in ex.items() if kx not in ("umi", "umi_reg", "mean_umi_reg", "sec", "arm", "tag")}
                info.update({"tag": f"{a}_k{k}", "arm": arm, "k_seed": k, "seeds": {"numpy": [SEED_XC, AREAS.index(a), k], "stripes_R": stripe_seed(a, k)},
                             "depth_k": kk, "b_panel": bp, "Mstar": dp["Mstar"][t], "displacement": mv["params"], "n_births_window": mv["n_births_window"],
                             "n_arrived_box": mv["n_arrived_box"], "n_captured": int(kept.sum()), "n_background": nb,
                             "share_bins_ct_over_1_over_k": float((kk * ct > 1).mean())})
                info["sec"] = round(time.time() - t0, 1); XE.write_export(out, b, X, cells, info)

if __name__ == "__main__":
    st = sys.argv[1]; areas = sys.argv[2].split(",") if len(sys.argv) > 2 else AREAS
    {"prep": stage_prep, "move": stage_move, "depth": stage_depth, "capture": stage_capture}[st](areas)
