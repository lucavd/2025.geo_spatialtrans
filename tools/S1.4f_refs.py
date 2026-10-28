#!/usr/bin/env python
"""tools/S1.4f_refs.py - S1.4f: intervalli di rilevanza nuovi calcolati sui dati veri PRIMA di ogni simulazione nuova (docs/steps/S1.4f.md).
R-2 (ripreso da S1.4e, ricalcolato come controllo): quota dei nuclei allungati (>= 1.2) con il vicino a meno di 5 um e asse maggiore di
    traverso al vicino (psi > 45 gradi); per ROI; intervallo = [minimo, massimo] dei 5 ROI. Deve coincidere con results/S1.4e/S1.4e_refs.csv.
R-4 (per la controprova CP-2): la stessa quota con il vicino a 5-6 um e a 6-8 um (classi non usate nella taratura); per ROI; [minimo, massimo].
Dati: results/S1.4e/explore/S1.4e_E5_nuclei.csv.gz (nuclei del segmentatore di forma, 5 ROI per tessuto). Uscita: results/S1.4f/S1.4f_refs.csv"""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/S1.4f"; OUT.mkdir(parents=True, exist_ok=True)
CLS = {"nn<5": (0, 5), "nn5-6": (5, 6), "nn6-8": (6, 8)}

def main():
    E = pd.read_csv(ROOT / "results/S1.4e/explore/S1.4e_E5_nuclei.csv.gz"); rows = []
    for A, g in E.groupby("archetype"):
        for c, (lo, hi) in CLS.items():
            el = g[(g.aspect >= 1.2) & (g.nn_um >= lo) & (g.nn_um < hi)]
            per = el.groupby("roi_id").psi.apply(lambda s: float(np.mean(s > 45)))
            rows.append(dict(ref="R-2" if c == "nn<5" else "R-4", archetype=A, cls=c, n=len(el), n_roi=len(per), est=float(per.mean()),
                             lo=float(per.min()), hi=float(per.max()), roi_values=";".join(f"{k}:{v:.4f}" for k, v in per.items())))
    R = pd.DataFrame(rows)
    old = pd.read_csv(ROOT / "results/S1.4e/S1.4e_refs.csv"); old = old[old.ref == "R-2"].set_index("archetype")
    r2 = R[R.ref == "R-2"].set_index("archetype"); dd = float((r2[["lo", "hi"]] - old.loc[r2.index, ["lo", "hi"]]).abs().max().max())
    print("R-2 contro S1.4e, scarto massimo:", dd); assert dd < 1e-12
    R.to_csv(OUT / "S1.4f_refs.csv", index=False); print(R.drop(columns="roi_values").round(3).to_string())

if __name__ == "__main__":
    main()
