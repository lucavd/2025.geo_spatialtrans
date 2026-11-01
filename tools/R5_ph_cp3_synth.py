"""tools/R5_ph_cp3_synth.py — R5, post hoc PH-1 (dichiarato dopo aver visto CP-3 nell'intestino): risposta dello stimatore a uno spostamento
noto della maschera (+2 µm righe, -2 µm colonne) su tessere sintetiche con bordo curvo, sigma in {2, 5} µm, w = 0.1, scarto vero 0 e (0, -2.5).
Uso: .venv/bin/python tools/R5_ph_cp3_synth.py <workers>. Scrive results/R5/R5_PH1_cp3_synth.csv."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools"); import R5_common as C
TILE = (slice(75, 225), slice(75, 225))
def one(a):
    s, dc, rep = a
    rng = np.random.default_rng([C.SEED, 77, int(s * 10), int((dc + 5) * 10), rep])
    S = C.synth_tile(rng, s, s, 0.10, 25.0, 0.0, dc, geometry="curve")
    out = {"sigma_true": s, "dc_true": dc, "rep": rep}
    for tag, f in (("orig", S["f_sub"]), ("shift", None)):
        if f is None:
            f = np.zeros_like(S["f_sub"]); f[4:, :-4] = S["f_sub"][:-4, 4:]
        fit, P, sel = C.fit_region(S["y"], f, S["rows"], S["cols"], TILE)
        for k in ("sigma", "w", "dr", "dc"): out[f"{tag}_{k}"] = fit[k]
    return out
if __name__ == "__main__":
    jobs = [(s, dc, r) for s in (2.0, 5.0) for dc in (0.0, -2.5) for r in range(10)]
    with Pool(int(sys.argv[1])) as p: D = pd.DataFrame(p.map(one, jobs, chunksize=1))
    D["resp_dr"] = D.shift_dr - D.orig_dr; D["resp_dc"] = D.shift_dc - D.orig_dc
    D.to_csv("results/R5/R5_PH1_cp3_synth.csv", index=False)
    print(D.groupby(["sigma_true", "dc_true"])[["orig_sigma", "orig_dc", "resp_dr", "resp_dc"]].median().round(3).to_string())
