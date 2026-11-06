"""tools/S2.3_select.py — S2.3, regola R-1 (taratura del campo sul nucleo) e giudizio di P-1...P-3 sul banco.
Uso: .venv/bin/python tools/S2.3_select.py combos_coarse <out.json>
     .venv/bin/python tools/S2.3_select.py combos_fine <benchdir> <out.json>     (griglia fine attorno al minimo grossolano, per tessuto)
     .venv/bin/python tools/S2.3_select.py choose <benchdir> <out_verify.json>    (scelta finale -> results/S2.3/S2.3_R1_choice.csv, S2.3_R1_loss.csv)
     .venv/bin/python tools/S2.3_select.py judge <benchdir>                       (verifica, repliche 10-14 -> results/S2.3/S2.3_P123_bench.csv)"""
import sys, os, json, glob
import numpy as np, pandas as pd
COARSE_L = [0.5, 0.75, 1.0, 1.5, 2.0]
COARSE_S2 = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
GROUPS = ["nucleus", "corona", "acellular"]
MIN_GROUPS = 20

def table(bd, rep0):
    rows = []
    for f in glob.glob(f"{bd}/grid/*_r{rep0:02d}.json"):
        d = json.load(open(f))
        if "error" in d: raise RuntimeError(f"{f}: {d['error'][-300:]}")
        for g in GROUPS:
            o = d["obs"][g]
            sims = [r[g] for r in d["reps"]]
            ok = o["n_groups"] >= MIN_GROUPS and np.isfinite(o["psi"])
            rows.append({"tag": d["tag"], "archetype": d["archetype"], "combo": d["combo"], "s2": d["s2"], "L": d["L"], "group": g, "ok": bool(ok),
                         "psi_obs": o["psi"], "psi_sim": float(np.mean([s["psi"] for s in sims])),
                         "rho1_obs": o["rho1"], "rho1_sim": float(np.mean([s["rho1"] for s in sims])),
                         "rho2_obs": o["rho2"], "rho2_sim": float(np.mean([s["rho2"] for s in sims])), "n_groups_obs": o["n_groups"], "nrep": len(sims)})
    T = pd.DataFrame(rows)
    T["q"] = T.psi_sim / T.psi_obs; T["e1"] = T.rho1_sim - T.rho1_obs; T["e2"] = T.rho2_sim - T.rho2_obs
    return T

def loss(T):
    N = T[(T.group == "nucleus") & T.ok]
    L = N.groupby(["archetype", "combo", "s2", "L"]).agg(q=("q", "mean"), e1=("e1", "mean"), n_roi=("tag", "nunique")).reset_index()
    L["loss"] = (L.q - 1) ** 2 + L.e1 ** 2
    return L

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "combos_coarse":
        cb = [{"combo": i * len(COARSE_L) + j, "s2": s2, "L": L} for i, s2 in enumerate(COARSE_S2) for j, L in enumerate(COARSE_L)]
        json.dump(cb, open(sys.argv[2], "w")); print(len(cb))
    elif mode == "combos_fine":
        bd = sys.argv[2]; L = loss(table(bd, 0)); cb = []
        for ai, A in enumerate(ARCH):
            b = L[L.archetype == A].sort_values("loss").iloc[0]
            k = 0
            for fs in (-0.5, -0.25, 0.0, 0.25, 0.5):
                for fl in (0.85, 1.0, 1.18):
                    if fs == 0.0 and fl == 1.0: continue
                    cb.append({"combo": 100 + 20 * ai + k, "archetype": A, "s2": float(b.s2 * 2 ** fs), "L": float(b.L * fl)}); k += 1
        json.dump(cb, open(sys.argv[3], "w")); print(len(cb))
    elif mode == "choose":
        bd = sys.argv[2]; L = loss(table(bd, 0)); os.makedirs("results/S2.3", exist_ok=True)
        L.sort_values(["archetype", "loss"]).to_csv("results/S2.3/S2.3_R1_loss.csv", index=False)
        ch = L.sort_values("loss").groupby("archetype").head(1).sort_values("archetype")
        ch[["archetype", "combo", "s2", "L", "q", "e1", "loss", "n_roi"]].to_csv("results/S2.3/S2.3_R1_choice.csv", index=False)
        print(ch.to_string())
        cb = [{"combo": 500 + ARCH.index(r.archetype), "archetype": r.archetype, "s2": float(r.s2), "L": float(r.L)} for r in ch.itertuples()]
        json.dump(cb, open(sys.argv[3], "w"))
    elif mode == "judge":
        bd = sys.argv[2]; T = table(bd, 10); T = T[T.combo >= 500]
        iv = pd.read_csv("results/S2.3/S2.3_E2_relevance.csv").set_index(["quantity", "archetype"])
        out = []
        for (A, g), d in T[T.ok].groupby(["archetype", "group"]):
            n = len(d)
            for m, lo, hi in (("q", 0.909, 1.1), ("e1", None, None), ("e2", None, None)):
                x = d[m].to_numpy(); est = x.mean(); se = x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
                if m != "q":
                    h = iv.loc[(f"bench_{g}_rho{m[1]}", A), "interval_halfwidth"]; lo, hi = -h, h
                a, b = est - 2 * se, est + 2 * se
                verdict = "confermata" if (a >= lo and b <= hi) else ("smentita" if (b < lo or a > hi) else "indeterminata")
                out.append({"archetype": A, "group": g, "measure": m, "n_roi": n, "estimate": est, "lo2se": a, "hi2se": b, "int_lo": lo, "int_hi": hi, "verdict": verdict})
        R = pd.DataFrame(out); R.to_csv("results/S2.3/S2.3_P123_bench.csv", index=False)
        T.to_csv("results/S2.3/S2.3_bench_verify_roi.csv", index=False)
        pd.set_option("display.width", 200); print(R.round(4).to_string())
