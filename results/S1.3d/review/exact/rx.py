"""review/exact/rx.py — revisore esatto S1.3d: geometria esatta INDIPENDENTE dall'arbitro (non importa exactgeom né S1.3d_arbiter).
Metodo: integrazione per strisce verticali (decomposizione trapezoidale). Fra due ascisse di rottura consecutive (vertici, incroci fra lati,
incroci dei lati con le rette y = y0, y = y1 del bin) la lunghezza della sezione verticale è lineare in x, quindi
∫ = (xb - xa) * L((xa + xb)/2) esattamente (regola del punto medio, esatta per funzioni lineari). Tutto in fractions.Fraction a partire dai
double bit-esatti letti dai parquet. Sezione di un poligono valido = coppie consecutive delle ordinate ordinate (regola pari-dispari).
Diverso dall'arbitro, che usa Green sulla 1-catena (O1) e Sutherland–Hodgman (O2)."""
from fractions import Fraction as F
import numpy as np, pyarrow.parquet as pq, math

def read_polys(path, want=None):
    """parquet (cell, part, ring, k, x, y) -> {cell: [ring, ...]} (tutti gli anelli di tutte le parti; pari-dispari non ne ha bisogno)"""
    t = pq.read_table(path); c = t.column("cell").to_numpy(); p = t.column("part").to_numpy(); r = t.column("ring").to_numpy()
    k = t.column("k").to_numpy(); x = t.column("x").to_numpy(); y = t.column("y").to_numpy()
    sel = np.ones(len(c), bool) if want is None else np.isin(c, np.array(sorted(want)))
    idx = np.nonzero(sel)[0]; idx = idx[np.lexsort((k[idx], r[idx], p[idx], c[idx]))]
    out = {}
    for q in idx:
        out.setdefault(int(c[q]), {}).setdefault((int(p[q]), int(r[q])), []).append((F(float(x[q])), F(float(y[q]))))
    res = {}
    for ci, rings in out.items():
        L = []
        for key in sorted(rings):
            rg = rings[key]
            if len(rg) > 1 and rg[0] == rg[-1]: rg = rg[:-1]          # anello chiuso nel parquet: tolgo il punto di chiusura
            L.append(rg)
        res[ci] = L
    return res

def edges(rings):
    E = []
    for rg in rings:
        n = len(rg)
        for q in range(n):
            a, b = rg[q], rg[(q + 1) % n]
            if a != b: E.append((a, b))
    return E

def shoelace(rings):
    """area con segno della 1-catena (esterni CCW +, buchi CW -) — solo come controllo incrociato, l'area 'vera' è slab"""
    s = F(0)
    for rg in rings:
        n = len(rg)
        for q in range(n):
            (x1, y1), (x2, y2) = rg[q], rg[(q + 1) % n]; s += x1 * y2 - x2 * y1
    return s / 2

def _ys_at(E, x):
    ys = []
    for (a, b) in E:
        ax, bx = a[0], b[0]
        if (ax < x < bx) or (bx < x < ax):
            ys.append(a[1] + (x - ax) * (b[1] - a[1]) / (bx - ax))
    ys.sort(); assert len(ys) % 2 == 0, "sezione dispari: poligono non chiuso/non valido"
    return [(ys[2 * q], ys[2 * q + 1]) for q in range(len(ys) // 2)]

def _len_clip(iv, y0, y1):
    s = F(0)
    for a, b in iv:
        lo = a if (y0 is None or a > y0) else y0; hi = b if (y1 is None or b < y1) else y1
        if hi > lo: s += hi - lo
    return s

def _cross_x_hline(E, yv, xa, xb):
    out = []
    for (a, b) in E:
        if (a[1] - yv) * (b[1] - yv) < 0:
            out.append(a[0] + (yv - a[1]) * (b[0] - a[0]) / (b[1] - a[1]))
    return [x for x in out if xa < x < xb]

def area_box(rings, x0, x1, y0, y1):
    """area esatta di poligono ∩ [x0,x1]x[y0,y1] (bordi del bin come Fraction)"""
    x0, x1, y0, y1 = F(x0), F(x1), F(y0), F(y1)
    E = edges(rings)
    xs = {x0, x1}
    for (a, b) in E:
        if x0 < a[0] < x1: xs.add(a[0])
    xs.update(_cross_x_hline(E, y0, x0, x1)); xs.update(_cross_x_hline(E, y1, x0, x1))
    xs = sorted(xs); A = F(0)
    for q in range(len(xs) - 1):
        xa, xb = xs[q], xs[q + 1]
        A += (xb - xa) * _len_clip(_ys_at(E, (xa + xb) / 2), y0, y1)
    return A

def area(rings):
    E = edges(rings); xs = sorted({a[0] for a, b in E}); A = F(0)
    for q in range(len(xs) - 1):
        xa, xb = xs[q], xs[q + 1]; A += (xb - xa) * _len_clip(_ys_at(E, (xa + xb) / 2), None, None)
    return A

def _seg_x(p, q, r, s):
    """ascissa dell'incrocio proprio/improprio fra i segmenti pq e rs (None se paralleli o disgiunti)"""
    d = (q[0] - p[0]) * (s[1] - r[1]) - (q[1] - p[1]) * (s[0] - r[0])
    if d == 0: return None
    t = ((r[0] - p[0]) * (s[1] - r[1]) - (r[1] - p[1]) * (s[0] - r[0])) / d
    u = ((r[0] - p[0]) * (q[1] - p[1]) - (r[1] - p[1]) * (q[0] - p[0])) / d
    if 0 <= t <= 1 and 0 <= u <= 1: return p[0] + t * (q[0] - p[0])
    return None

def _inter_len(I, J):
    s = F(0); a = b = 0
    while a < len(I) and b < len(J):
        lo = max(I[a][0], J[b][0]); hi = min(I[a][1], J[b][1])
        if hi > lo: s += hi - lo
        if I[a][1] < J[b][1]: a += 1
        else: b += 1
    return s

def area_inter(RA, RB, EB=None):
    """area esatta di A ∩ B (A, B poligoni validi, anche multipli/bucati)"""
    EA = edges(RA); EB = edges(RB) if EB is None else EB
    ax0 = min(min(a[0], b[0]) for a, b in EA); ax1 = max(max(a[0], b[0]) for a, b in EA)
    ay0 = min(min(a[1], b[1]) for a, b in EA); ay1 = max(max(a[1], b[1]) for a, b in EA)
    EBl = [e for e in EB if max(e[0][0], e[1][0]) >= ax0 and min(e[0][0], e[1][0]) <= ax1]
    bx0 = min(min(a[0], b[0]) for a, b in EBl) if EBl else ax1; bx1 = max(max(a[0], b[0]) for a, b in EBl) if EBl else ax0
    lo, hi = max(ax0, bx0), min(ax1, bx1)
    if hi <= lo: return F(0)
    xs = {lo, hi}
    for (a, b) in EA + EBl:
        if lo < a[0] < hi: xs.add(a[0])
    # incroci fra lati di A e di B (solo lati con riquadri sovrapposti)
    EBy = [e for e in EBl if max(e[0][1], e[1][1]) >= ay0 and min(e[0][1], e[1][1]) <= ay1]
    for (p, q) in EA:
        px0, px1 = min(p[0], q[0]), max(p[0], q[0]); py0, py1 = min(p[1], q[1]), max(p[1], q[1])
        for (r, s) in EBy:
            if max(r[0], s[0]) < px0 or min(r[0], s[0]) > px1 or max(r[1], s[1]) < py0 or min(r[1], s[1]) > py1: continue
            x = _seg_x(p, q, r, s)
            if x is not None and lo < x < hi: xs.add(x)
    xs = sorted(xs); A = F(0)
    for k in range(len(xs) - 1):
        xa, xb = xs[k], xs[k + 1]; xm = (xa + xb) / 2
        A += (xb - xa) * _inter_len(_ys_at(EA, xm), _ys_at(EBl, xm))
    return A

def fbox(rings):
    xs = [float(p[0]) for r in rings for p in r]; ys = [float(p[1]) for r in rings for p in r]
    return min(xs), max(xs), min(ys), max(ys)

def pieces(rings, X0, Y0, s, nx, ny):
    """tutti i pezzi non nulli cella x bin (indici ix = floor((x - X0)/s)), esatti; riquadro allargato di 1 bin come margine"""
    x0, x1, y0, y1 = fbox(rings)
    ia, ib = max(0, math.floor((x0 - X0) / s) - 1), min(nx - 1, math.floor((x1 - X0) / s) + 1)
    ja, jb = max(0, math.floor((y0 - Y0) / s) - 1), min(ny - 1, math.floor((y1 - Y0) / s) + 1)
    out = {}
    for ix in range(ia, ib + 1):
        for iy in range(ja, jb + 1):
            a = area_box(rings, F(X0) + s * ix, F(X0) + s * (ix + 1), F(Y0) + s * iy, F(Y0) + s * (iy + 1))
            if a != 0: out[(ix, iy)] = a
    return out

def cmpf(exact, reported):
    """(uguale al double più vicino?, |diff| relativa) — float(Fraction) è arrotondato correttamente"""
    fe = float(exact); eq = (fe == reported)
    den = max(abs(fe), 1e-300); return eq, (0.0 if eq else abs(float(F(reported) - exact)) / den)
