"""tools/S2.0_stripes.py — S2.0, M1 (BL-192): strisce della griglia per sezione. Pre-registrazione docs/steps/S2.0.md (afef06b, a8690b3, f2a23d7).
Uso: .venv/bin/python tools/S2.0_stripes.py <outdir> [KEY ...]. Legge /mnt/micron/geo_spatialtrans/R5/run3/section/<KEY>.npz (UMI per bin a 2 µm,
geni della filtrata, DEC-036; in_tissue di Space Ranger). Scrive <outdir>/stripes/<KEY>.json e <KEY>_factors.npz.
P-1: deviazione standard riproducibile dei fattori non periodici di riga e colonna (meta' disgiunte; righe/colonne con >= 200 bin di tessuto per
meta'), bootstrap per blocchi di 25 righe o colonne. P-2: deviazione standard riproducibile fra tessere di 300 µm del logaritmo dei rapporti
dispari/pari (meta' sinistra/destra per le colonne, sopra/sotto per le righe), bootstrap delle tessere. B = 2000, seme 20261103."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, "tools")
import S2_0_common as S

R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
MIN_HALF = 200

def lag1_repro(r1, r2):
    """Autocorrelazione a passo 1 della parte riproducibile: cov(r1[i], r2[i+1]) / cov(r1[i], r2[i])."""
    ok0 = np.isfinite(r1) & np.isfinite(r2)
    ok1 = np.isfinite(r1[:-1]) & np.isfinite(r2[1:])
    if ok0.sum() < 10 or ok1.sum() < 10: return float("nan")
    c0 = np.cov(r1[ok0], r2[ok0])[0, 1]; c1 = np.cov(r1[:-1][ok1], r2[1:][ok1])[0, 1]
    return float(c1 / c0) if c0 > 0 else float("nan")

def run(KEY, out):
    t0 = time.time()
    z = np.load(f"{R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); T = z["in_tissue"].astype(bool)
    al, be, mr, mc = S.stripe_factors(U, T)
    res = {"key": KEY, "n_tissue_bins": int(T.sum()), "shape": list(U.shape), "seed": S.SEED, "B": S.B}
    fac = {"alpha": al, "beta": be}
    for ax, f, other, m, nm in ((0, al, be, mr, "rows"), (1, be, al, mc, "cols")):
        idx = np.arange(len(f))
        mu, npf = S.periodic_split(f, idx)
        h1, h2, n1, n2 = S.half_split_factors(U, T, ax)      # due adattamenti separati sulle meta' disgiunte (come pre-registrato)
        ok = (n1 >= MIN_HALF) & (n2 >= MIN_HALF)
        _, r1 = S.periodic_split(np.where(ok, h1, np.nan), idx); _, r2 = S.periodic_split(np.where(ok, h2, np.nan), idx)
        sd = S.repro_sd(r1, r2); lo, hi, _ = S.block_boot_sd(r1, r2, block=25)
        okf = np.isfinite(npf)
        res[nm] = {"odd_over_even": float(np.exp(mu[1] - mu[0])), "n_eligible": int(ok.sum()), "n_with_factor": int(okf.sum()),
                   "sd_repro": sd, "sd_lo": lo, "sd_hi": hi, "sd_total_np": float(np.nanstd(npf)),
                   "corr_halves": float(np.corrcoef(r1[ok & np.isfinite(r1) & np.isfinite(r2)], r2[ok & np.isfinite(r1) & np.isfinite(r2)])[0, 1]),
                   "lag1_repro": lag1_repro(r1, r2), "lag1_total": float(np.corrcoef(npf[:-1][okf[:-1] & okf[1:]], npf[1:][okf[:-1] & okf[1:]])[0, 1])}
        fac[f"{nm}_half1"], fac[f"{nm}_half2"] = r1, r2
    tp = S.tile_parity(U, T)
    res["tiles"] = {"n_tiles": int(len(tp))}
    for nm, a, b in (("cols", 2, 3), ("rows", 4, 5)):
        sd = S.repro_sd(tp[:, a], tp[:, b])
        lo, hi, _ = S.boot_rows(lambda ii: S.repro_sd(tp[ii, a], tp[ii, b]), len(tp))
        allv = 0.5 * (tp[:, a] + tp[:, b])
        res["tiles"][nm] = {"sd_between_repro": sd, "sd_lo": lo, "sd_hi": hi, "median_ratio": float(np.exp(np.median(allv))),
                            "p10_ratio": float(np.exp(np.percentile(allv, 10))), "p90_ratio": float(np.exp(np.percentile(allv, 90)))}
    fac["tiles"] = tp
    np.savez_compressed(f"{out}/stripes/{KEY}_factors.npz", **fac)
    res["sec"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{out}/stripes/{KEY}.json", "w"), indent=1)
    print(KEY, json.dumps({k: res[k] for k in ("rows", "cols")})[:400], res["tiles"], f"{res['sec']} s", flush=True)


# ------------------------------------------------------------------ E3: precisione dello stimatore di P-1 su tessuto con struttura forte
def e3(out, seeds=4):
    sys.path.insert(0, "tools")
    import importlib.util
    sp = importlib.util.spec_from_file_location("t", "tools/S2.0_test.py"); T = importlib.util.module_from_spec(sp); sp.loader.exec_module(T)
    rows = []
    for amp in (0.2, 0.4, 0.8):
        for sd in (0.0, 0.05):
            for k in range(seeds):
                U, Tm = T.synth_section(sd, 50 + k, amp=amp)
                o = T.stripe_summary(U, Tm)
                rows.append({"amp": amp, "sd_true": sd, "seed": k, "sd_rows": o["rows"]["sd"], "sd_cols": o["cols"]["sd"]})
                print(rows[-1], flush=True)
    import pandas as pd
    pd.DataFrame(rows).to_csv(f"{out}/stripes/E3_precision.csv", index=False)

if __name__ == "__main__":
    out = sys.argv[1]; keys = sys.argv[2:] or ["A1", "A23", "A4", "A5", "A6"]
    os.makedirs(f"{out}/stripes", exist_ok=True)
    if keys == ["e3"]:
        e3(out); sys.exit(0)
    for k in keys:
        run(k, out)
