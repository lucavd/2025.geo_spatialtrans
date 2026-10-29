"""tools/calibration_fit.py - S1.6 (BL-169): funzioni comuni per la scelta dei parametri calibrati della semina (vuoti, densita'), corrette
rispetto agli strumenti di S1.5d e S1.5e (che restano com'erano per riprodurre quei numeri). Usato dalle calibrazioni future (BL-168).
- fold_loss: perdita della piega = somma sulle ROI di addestramento e sui bersagli di (simulato - vero)^2 / varianza del bersaglio fra le
  sole ROI di addestramento (S1.5d ed S1.5e usavano la varianza fra tutte le 5 ROI anche nelle pieghe lascia-fuori: revisione RC-8 di S1.5e).
- grid_edges: bordi della griglia toccati dall'ottimo, alti e bassi, per ogni parametro; per un parametro con lo zero come «spento» (frazione
  dei vuoti, ampiezza della densita') il bordo basso e' il piu' piccolo valore positivo (S1.5e segnalava solo i bordi alti: RC-9).
- area_weighted: media di una quantita' per regione pesata per l'area della regione (S1.5e pesava la frazione dei vuoti per i punti della
  griglia: RC-6)."""
import numpy as np, pandas as pd

def fold_loss(sim, real, targets, train_rois):
    """sim: righe combo, roi_id, <bersagli> (media delle repliche); real: righe roi_id, <bersagli>. Restituisce una Series combo -> perdita."""
    r = real[real.roi_id.isin(train_rois)]
    var = r[targets].var(ddof=1)
    s = sim[sim.roi_id.isin(train_rois)].merge(r[["roi_id"] + targets], on="roi_id", suffixes=("", "_real"))
    loss = sum((s[t] - s[t + "_real"]) ** 2 / var[t] for t in targets)
    n = s.assign(L=loss).groupby("combo").agg(L=("L", "sum"), n_roi=("roi_id", "nunique"))
    return n.L[n.n_roi == len(set(train_rois))]

def grid_edges(best, grid, params, zero_off=()):
    """best: riga (dict o Series) dell'ottimo; grid: DataFrame dei candidati; params: nomi dei parametri; zero_off: parametri per cui 0 = spento.
    Restituisce la lista dei bordi toccati, per esempio ['void_frac_hi', 'void_ell_lo']."""
    out = []
    for p in params:
        v = best[p]
        if v is None or (isinstance(v, float) and np.isnan(v)): continue
        vals = np.asarray(sorted(set(pd.to_numeric(grid[p], errors="coerce").dropna())), float)
        if p in zero_off:
            if v == 0: continue
            vals = vals[vals > 0]
        if not len(vals): continue
        if v >= vals.max(): out.append(f"{p}_hi")
        if v <= vals.min(): out.append(f"{p}_lo")
    return out

def area_weighted(values, areas):
    v = np.asarray(values, float); a = np.asarray(areas, float); k = np.isfinite(v) & np.isfinite(a) & (a > 0)
    return float(np.sum(v[k] * a[k]) / np.sum(a[k])) if k.any() else float("nan")
