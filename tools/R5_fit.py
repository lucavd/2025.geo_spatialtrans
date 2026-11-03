"""tools/R5_fit.py — R5, adattamenti per tessera (resumable: un json per tessera). Uso: .venv/bin/python tools/R5_fit.py <outdir> [KEY ...]
Numero di processi dal file <outdir>/WORKERS (letto all'avvio). Per ogni tessera: K2 (primaria, fondo libero, come pre-registrato), K1, K3, CP-3
(maschera spostata di +2 µm righe, -2 µm colonne), S-0 (fondo fissato dal vetro lontano), S-1 (maschera alternativa), S-2 (parita' dal vetro),
S-3 (solo vetro), S-4 (200 µm centrali); informazione di Fisher dello scarto; profili per distanza dal bordo (D-1).
Pre-registrazione: docs/steps/R5.md (c5c9fce, ceca77b); deviazioni D-1, D-2 dopo G-1."""
import sys, os, json, glob, time, traceback
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C

OUT = sys.argv[1]; KEYS = sys.argv[2:] or ["A1", "A23", "A4", "A5", "A6"]
TB, MB = 150, 75
TILE = (slice(MB, MB + TB), slice(MB, MB + TB))
TILE200 = (slice(MB + 25, MB + 125), slice(MB + 25, MB + 125))
DBINS = np.arange(-20, 152, 2.0)

def region_P(y, f, rows, cols, valid, pf_from_glass=False):
    P = C.prepare_region(y, f, rows, cols)
    if pf_from_glass:
        cls = C.parity_class(rows, cols)
        pf = C.parity_factors(y.astype(float), cls, (P["occ"] == 0) & (P["dist"] >= 10) & valid)
        P = C.prepare_region(y, f, rows, cols)          # densita' ricalcolata con i fattori del vetro
        c = pf[cls]; rho = C.rho_field(y / c, P["occ"] >= 1)     # come la primaria: tutti i bin interamente di tessuto (rilievo F13 corretto)
        P.update({"pf": pf, "c": c, "rho": rho, "Tft": C.source_ft(C.source_sub(f, rho), P["pshape"])})
    return P

def fisher_delta(P, fit, y, sel, b_fixed, h=0.1):
    """Matrice di informazione (2 x 2) dello scarto: derivate seconde centrali della NLL profilata (a profilato, b fissato)."""
    def nll(dr, dc):
        g = C.expected_bins(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], fit["sigma_r"], fit["sigma_c"], fit["w"], fit["ell"], dr, dc)[sel]
        return C.fit_a(y[sel].astype(float), g, P["c"][sel], b_fixed)[2]
    d0r, d0c = fit["dr"], fit["dc"]; f0 = nll(d0r, d0c)
    frr = (nll(d0r + h, d0c) - 2 * f0 + nll(d0r - h, d0c)) / h ** 2
    fcc = (nll(d0r, d0c + h) - 2 * f0 + nll(d0r, d0c - h)) / h ** 2
    frc = (nll(d0r + h, d0c + h) - nll(d0r + h, d0c - h) - nll(d0r - h, d0c + h) + nll(d0r - h, d0c - h)) / (4 * h * h)
    return [[float(frr), float(frc)], [float(frc), float(fcc)]]

def profile(P, y, sel, fits):
    d = P["dist"][sel]; k = np.digitize(d, DBINS)
    out = {"edges": DBINS.tolist(), "n": np.bincount(k, minlength=len(DBINS) + 1).tolist(),
           "y": np.bincount(k, weights=y[sel], minlength=len(DBINS) + 1).tolist()}
    for nm, ft in fits.items():
        lam = C.predict(ft, P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], P["c"])[sel]
        out["pred_" + nm] = np.bincount(k, weights=lam, minlength=len(DBINS) + 1).tolist()
    out["rho_tile_mean"] = float(P["rho"][TILE][P["occ"][TILE] >= 1].mean()) if (P["occ"][TILE] >= 1).any() else float("nan")
    return out

def run_tile(args):
    key, path, bfar = args
    fout = path.replace(".npz", ".json")
    if os.path.exists(fout):
        return fout
    t0 = time.time()
    try:
        z = np.load(path)
        y = z["U"].astype(float); f = z["f_sub"].astype(float); fa = z["f_alt"].astype(float)
        R0, C0 = int(z["R0"]), int(z["C0"])
        rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
        valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
        P = region_P(y, f, rows, cols, valid)
        sel = C.likelihood_mask(P["occ"], P["dist"], TILE) & valid
        res = {"key": key, "tile": os.path.basename(path)[:-4], "b_far": bfar, "pf": P["pf"].tolist()}
        # revisione numerica F-1: K2 con griglia di partenza (C.fit_k2_screened); K3 e sensibilita' partono da K2 e dalla seconda soluzione
        k2, k1 = C.fit_k2_screened(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, None)
        S2 = [k2] + ([dict(k2, **k2["alt"])] if k2.get("alt") else [])
        st3 = [dict(k2, sigma_r=k2["sigma"], sigma_c=k2["sigma"]), dict(k2, sigma_r=k2["sigma"] * 1.3, sigma_c=k2["sigma"] / 1.3)]
        k3 = C.fit_model("K3", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, st3, None)
        res.update({"K2": k2, "K1": k1, "K3": k3, "fisher_delta_K2": fisher_delta(P, k2, y, sel, k2["b"]), "profile": profile(P, y, sel, {"K1": k1, "K2": k2, "K3": k3})})
        # S-0: fondo fissato dal vetro lontano della sezione (variante scartata in G-1, qui come sensibilita')
        res["S0_bfixed"] = C.fit_model("K2", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, S2, bfar, polish=False)
        res["S0_fixed_ok"] = True
        # S-3: solo vetro; S-4: 200 µm centrali
        sel3 = C.likelihood_mask(P["occ"], P["dist"], TILE, glass_only=True) & valid
        res["S3_glass"] = C.fit_model("K2", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel3, S2, None, polish=False)
        sel4 = C.likelihood_mask(P["occ"], P["dist"], TILE200) & valid
        res["S4_200um"] = C.fit_model("K2", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel4, S2, None, polish=False) if sel4.sum() > 500 else None
        # CP-3: maschera spostata di +2 µm lungo le righe e -2 µm lungo le colonne (4 sottopixel)
        fs = np.zeros_like(f); fs[4:, :-4] = f[:-4, 4:]
        Ps = region_P(y, fs, rows, cols, valid)
        sels = C.likelihood_mask(Ps["occ"], Ps["dist"], TILE) & valid
        st = dict(k2); st["dr"] = k2["dr"] - 2.0; st["dc"] = k2["dc"] + 2.0
        res["CP3_shift"] = C.fit_model("K2", Ps["Tft"], Ps["kr"], Ps["kc"], Ps["pshape"], Ps["region_bins"], y, Ps["c"], sels, [st, k2], None)
        # S-1: maschera alternativa
        Pa = region_P(y, fa, rows, cols, valid)
        sela = C.likelihood_mask(Pa["occ"], Pa["dist"], TILE) & valid
        res["S1_altmask"] = C.fit_model("K2", Pa["Tft"], Pa["kr"], Pa["kc"], Pa["pshape"], Pa["region_bins"], y, Pa["c"], sela, S2, None, polish=False)
        # S-2: parita' dal vetro
        Pg = region_P(y, f, rows, cols, valid, pf_from_glass=True)
        res["S2_parity_glass"] = C.fit_model("K2", Pg["Tft"], Pg["kr"], Pg["kc"], Pg["pshape"], Pg["region_bins"], y, Pg["c"], sel, S2, None, polish=False)
        res["pf_glass"] = Pg["pf"].tolist(); res["S2_fixed_ok"] = True
        # S-W (revisione numerica F-2): limite w <= 0.5 della prima versione, come sensibilita'
        # se K2 ha w <= 0.5 la soluzione e' gia' ammissibile e coincide con quella vincolata: si copia
        wb = C.BOUNDS["w"]; C.BOUNDS["w"] = (0.0, min(0.5, wb[1]))
        if k2["w"] <= 0.5:
            res["SW_w05"] = dict(k2, copied_from_K2=True); C.BOUNDS["w"] = wb
        else:
          try:
            res["SW_w05"] = C.fit_model("K2", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel,
                                        [dict(x, w=min(x["w"], 0.5)) for x in S2] + C.STARTS["K2"], None, polish=True)
          finally:
            C.BOUNDS["w"] = wb
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"key": key, "tile": os.path.basename(path)[:-4], "error": traceback.format_exc()}
    tmp = fout + ".tmp"; json.dump(res, open(tmp, "w")); os.replace(tmp, fout)
    return fout

if __name__ == "__main__":
    jobs = []
    for key in KEYS:
        bfar = json.load(open(f"{OUT}/tiles/{key}_background.json"))["b_far"]
        jobs += [(key, p, bfar) for p in sorted(glob.glob(f"{OUT}/tiles/{key}/t*.npz"))]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 20
    todo = [j for j in jobs if not os.path.exists(j[1].replace(".npz", ".json"))]
    print(f"tiles {len(jobs)} todo {len(todo)} workers {nw}", flush=True)
    t0 = time.time()
    with Pool(nw) as pool:
        for i, f in enumerate(pool.imap_unordered(run_tile, todo, chunksize=1), 1):
            if i % 10 == 0 or i == len(todo):
                print(f"{i}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("done", flush=True)
