"""tools/S1.5c_type_nuc_area.py - S1.5c: area nucleare mediana reale per tipo (tutti i nuclei con un tipo, etichette R6 PH-1, come la regola dei
gruppi di S1.5b), per la distanza minima per nucleo (scelta 1). Uscita: results/S1.5c/S1.5c_type_nuc_area.csv"""
import importlib.util, pandas as pd
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
rows = []
for A in cal.ARCH:
    lab, _ = cal.load(A); lab = lab[~lab.type.isin(list(cal.UNTYPED))]
    for t, g in lab.groupby("type"): rows.append(dict(archetype=A, cell_type=t, median_area_um2=g.area_um2.median(), n=len(g)))
pd.DataFrame(rows).to_csv("results/S1.5c/S1.5c_type_nuc_area.csv", index=False); print(pd.DataFrame(rows).to_string())
