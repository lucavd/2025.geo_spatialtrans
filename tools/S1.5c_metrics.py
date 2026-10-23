"""tools/S1.5c_metrics.py - S1.5c: metriche spaziali per tipo, stesse funzioni per i nuclei veri e per le simulazioni (definizioni: docs/steps/S1.5c.md).
Ingresso per ROI: coordinate dei centri (um), tipo (stringa; None/«non risolto»/«non classificabile» = senza tipo), indicatore di cella interna.
- nn_rel_t: mediana della distanza dal vicino piu' vicino (di qualsiasi tipo, tutti i centri) delle celle interne del tipo t,
  divisa per la stessa mediana su tutte le celle interne (adimensionale, non dipende dalla densita').
- markconn_t(b): P(vicino con tipo = t | focale di tipo t, distanza in b) / f_t, fra i centri con un tipo; focali = celle interne del tipo t
  (al piu' 1500, campione fisso con seme 20261023); f_t = frazione di t fra i centri con un tipo del ROI; bin 0-10-20-40-80-160 um.
- cv_local4: coefficiente di variazione del numero di vicini entro R = 4 x (mediana della distanza dal vicino piu' vicino delle celle interne),
  focali = celle interne (tutti i centri contati)."""
import numpy as np
from scipy.spatial import cKDTree
UNTYPED = {"non risolto", "non classificabile"}
BINS = np.array([0, 10, 20, 40, 80, 160.])
def roi_metrics(x, y, typ, interior, main_types, seed=20261023):
    P = np.c_[x, y]; tr = cKDTree(P); I = np.asarray(interior, bool)
    tt = np.array([None if (t is None or (isinstance(t, float) and np.isnan(t)) or t in UNTYPED) else t for t in typ], dtype=object)
    typed = np.array([t is not None for t in tt])
    dd, ii = tr.query(P, k=2); nn = dd[:, 1]; m_all = np.median(nn[I])
    out = {"nn_median_all": m_all}
    cnt = np.array([len(v) - 1 for v in tr.query_ball_point(P[I], 4 * m_all)]); out["cv_local4"] = cnt.std(ddof=1) / cnt.mean()
    rng = np.random.default_rng(seed)
    for t in main_types:
        foc = np.where(I & (tt == t))[0]
        if len(foc) < 20: continue
        out[("nn_rel", t)] = np.median(nn[foc]) / m_all
        f_t = np.mean(tt[typed] == t)
        sub = foc if len(foc) <= 1500 else np.sort(rng.choice(foc, 1500, replace=False))
        num = np.zeros(len(BINS) - 1); den = np.zeros(len(BINS) - 1)
        for k, lst in zip(sub, tr.query_ball_point(P[sub], BINS[-1])):
            lst = np.asarray(lst, dtype=np.int64); lst = lst[(lst != k) & typed[lst]]
            if not len(lst): continue
            b = np.minimum(np.digitize(np.hypot(*(P[lst] - P[k]).T), BINS) - 1, len(BINS) - 2)   # d = 160 esatto nell'ultimo bin (revisione RC-08)
            np.add.at(num, b, tt[lst] == t); np.add.at(den, b, 1)
        for kb in range(len(BINS) - 1):
            out[("markconn", t, f"{int(BINS[kb])}-{int(BINS[kb+1])}")] = num[kb] / den[kb] / f_t if den[kb] else np.nan
    return out
