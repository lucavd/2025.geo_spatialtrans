"""tools/test_calibration_fit.py - S1.6 (BL-169): test di tools/calibration_fit.py. Stampa PASS/FAIL per asserzione e RESULT n/m PASS."""
import importlib.util, numpy as np, pandas as pd
s = importlib.util.spec_from_file_location("cf", "tools/calibration_fit.py"); cf = importlib.util.module_from_spec(s); s.loader.exec_module(cf)
res = []
def ck(name, ok, val=""):
    res.append(bool(ok)); print(("PASS" if ok else "FAIL"), name, val)
# T1: varianza delle sole ROI di addestramento
real = pd.DataFrame(dict(roi_id=["a", "b", "c", "d", "e"], t=[1.0, 2.0, 3.0, 4.0, 100.0]))
sim = pd.DataFrame(dict(combo=["x"] * 5 + ["y"] * 5, roi_id=list("abcde") * 2, t=[2.0, 3, 4, 5, 101] + [1.0, 2, 3, 4, 100]))
L = cf.fold_loss(sim, real, ["t"], ["a", "b", "c", "d"])
var4 = np.var([1, 2, 3, 4], ddof=1)
ck("T1a perdita della piega con la varianza delle 4 ROI di addestramento (la quinta, anomala, non entra)", abs(L["x"] - 4 / var4) < 1e-12 and L["y"] == 0, L.to_dict())
L5 = cf.fold_loss(sim, real, ["t"], list("abcde"))
ck("T1b piega con tutte le ROI: varianza delle 5", abs(L5["x"] - 5 / real.t.var(ddof=1)) < 1e-12)
ck("T1c una combinazione senza tutte le ROI della piega e' esclusa", "z" not in cf.fold_loss(pd.concat([sim, pd.DataFrame(dict(combo=["z"], roi_id=["a"], t=[1.0]))]), real, ["t"], list("abcd")).index)
# T2: bordi alti e bassi
G = pd.DataFrame(dict(void_frac=[0, 0.025, 0.05, 0.1, 0.25], void_ell_um=[np.nan, 8, 16, 32, 64]))
ck("T2a bordo alto della frazione e basso della lunghezza", cf.grid_edges(dict(void_frac=0.25, void_ell_um=8), G, ["void_frac", "void_ell_um"], ("void_frac",)) == ["void_frac_hi", "void_ell_um_lo"])
ck("T2b bordo basso della frazione = il piu' piccolo valore positivo (lo zero e' «spento», non un bordo)",
   cf.grid_edges(dict(void_frac=0.025, void_ell_um=16), G, ["void_frac", "void_ell_um"], ("void_frac",)) == ["void_frac_lo"])
ck("T2c frazione 0 (spento): nessun bordo per quel parametro", cf.grid_edges(dict(void_frac=0.0, void_ell_um=np.nan), G, ["void_frac", "void_ell_um"], ("void_frac",)) == [])
ck("T2d interno: nessun bordo", cf.grid_edges(dict(void_frac=0.05, void_ell_um=16), G, ["void_frac", "void_ell_um"], ("void_frac",)) == [])
# T3: media pesata per area
ck("T3 media pesata per area, non per punti", abs(cf.area_weighted([0.5, 0.0], [100.0, 300.0]) - 0.125) < 1e-15)
print(f"RESULT {sum(res)}/{len(res)} PASS")
