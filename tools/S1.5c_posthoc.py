"""tools/S1.5c_posthoc.py - S1.5c, analisi post hoc (dopo i verdetti, senza verdetto). PH-1: da dove viene la variabilita' dei territori che manca (P-2)?
Per ROI, celle interne: quantili dell'area del territorio divisa per la mediana (q10, q90, q99), frazione di territori sopra 3 volte la mediana, e CV
dopo aver tolto i territori sopra 3 volte la mediana (CV troncato), reale contro bracci U e vincente (10 semi, media). Stima = media dei 5 ROI +- 2 ES.
Uso: .venv/bin/python tools/S1.5c_posthoc.py <braccio vincente>   Uscita: results/S1.5c/S1.5c_posthoc_terr.csv"""
import sys, glob, importlib.util, numpy as np, pandas as pd
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
WIN = sys.argv[1] if len(sys.argv) > 1 else "N-a"; OUT = "/mnt/micron/geo_spatialtrans/S1.5c"
def stats(a):
    m = np.median(a); r = a / m; k = r <= 3
    return dict(q10=np.quantile(r, 0.1), q90=np.quantile(r, 0.9), q99=np.quantile(r, 0.99), frac_gt3=np.mean(r > 3),
                cv=a.std(ddof=1) / a.mean(), cv_trim=a[k].std(ddof=1) / a[k].mean())
rows = []
for A in cal.ARCH:
    _, d = cal.load(A)
    for r, e in d.groupby("roi_id"):
        rows.append(dict(archetype=A, roi_id=r, source="reale", **stats(e.a_t.values[e.interior.values])))
        for arm in ["U", WIN]:
            st = [stats(c.a_t.values[c.interior.values]) for c in (pd.read_parquet(f) for f in sorted(glob.glob(f"{OUT}/cells/{arm}/{A}_{r}_s*.parquet")))]
            rows.append(dict(archetype=A, roi_id=r, source=arm, **pd.DataFrame(st).mean().to_dict()))
D = pd.DataFrame(rows); D.to_csv("results/S1.5c/S1.5c_posthoc_terr.csv", index=False)
pd.set_option("display.width", 250)
print(D.groupby(["archetype", "source"])[["q10", "q90", "q99", "frac_gt3", "cv", "cv_trim"]].mean().round(3).to_string())
