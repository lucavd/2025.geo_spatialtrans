
"""exactgeom.py — geometria esatta su razionali (fractions.Fraction) per poligoni a coordinate double.
Revisore 'exact', S1.3b. Poligono = lista di parti; parte = lista di anelli; anello = lista di (Fraction, Fraction)
senza ripetizione del punto di chiusura. Area dell'intersezione: teorema di Green sulla 1-catena del bordo di P∩Q
(sottosegmenti di dP interni a Q, di dQ interni a P, lati collineari concordi una sola volta)."""
from fractions import Fraction as F
import struct

def fx(h): return F(float.fromhex(h)) if isinstance(h, str) else F(h)

def ulp_key(x):
    i = struct.unpack('<q', struct.pack('<d', float(x)))[0]
    return i if i >= 0 else -(i & 0x7fffffffffffffff)
def ulp_diff(a, b): return abs(ulp_key(a) - ulp_key(b))

def cross(ax, ay, bx, by): return ax * by - ay * bx
def orient(a, b, c): return cross(b[0] - a[0], b[1] - a[1], c[0] - a[0], c[1] - a[1])

def ring_area2(r):
    n = len(r); return sum(r[k][0] * r[(k + 1) % n][1] - r[(k + 1) % n][0] * r[k][1] for k in range(n))
def strip_close(r):
    r = list(r)
    while len(r) > 1 and r[0] == r[-1]: r = r[:-1]
    return r
def normalize(poly):
    """esterni CCW, buchi CW; toglie il punto di chiusura"""
    out = []
    for part in poly:
        np_ = []
        for k, r in enumerate(part):
            r = strip_close(r); a = ring_area2(r)
            if (k == 0 and a < 0) or (k > 0 and a > 0): r = r[::-1]
            np_.append(r)
        out.append(np_)
    return out
def area(poly):
    return sum(sum(ring_area2(r) for r in part) for part in normalize(poly)) / 2   # buchi CW: negativi
def edges(poly):
    E = []
    for part in normalize(poly):
        for r in part:
            n = len(r)
            for k in range(n): E.append((r[k], r[(k + 1) % n]))
    return E
def on_seg(p, a, b):
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
def pip(p, E):
    """'in' | 'out' | 'on' (pari-dispari su tutti gli anelli, esatto)"""
    inside = False; px, py = p
    for a, b in E:
        ay, by_ = a[1], b[1]
        if (ay > py and by_ > py) or (ay < py and by_ < py): continue        # lato tutto sopra/sotto: né 'on' né crossing
        if (a[0] < px and b[0] < px) and ay != by_ and not (ay == py or by_ == py): continue  # tutto a sinistra, non orizzontale
        if on_seg(p, a, b): return 'on'
        if (a[1] > p[1]) != (b[1] > p[1]):
            xint = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if p[0] < xint: inside = not inside
    return 'in' if inside else 'out'

def _splits(EA, EB):
    """parametri t di taglio dei lati di A dovuti ai lati di B (esatti, inclusi sovrapposti collineari)"""
    S = [ {F(0), F(1)} for _ in EA ]
    for ia, (a, b) in enumerate(EA):
        dx, dy = b[0] - a[0], b[1] - a[1]
        minx, maxx, miny, maxy = min(a[0], b[0]), max(a[0], b[0]), min(a[1], b[1]), max(a[1], b[1])
        for c, d in EB:
            if max(c[0], d[0]) < minx or min(c[0], d[0]) > maxx or max(c[1], d[1]) < miny or min(c[1], d[1]) > maxy: continue
            ex, ey = d[0] - c[0], d[1] - c[1]
            den = cross(dx, dy, ex, ey)
            if den != 0:
                t = cross(c[0] - a[0], c[1] - a[1], ex, ey) / den
                u = cross(c[0] - a[0], c[1] - a[1], dx, dy) / den
                if 0 <= t <= 1 and 0 <= u <= 1: S[ia].add(t)
            elif cross(c[0] - a[0], c[1] - a[1], dx, dy) == 0:   # collineari
                L2 = dx * dx + dy * dy
                for q in (c, d):
                    t = ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / L2
                    if 0 <= t <= 1: S[ia].add(t)
    return S

def _subsegs(EA, S):
    for (a, b), ts in zip(EA, S):
        ts = sorted(ts); dx, dy = b[0] - a[0], b[1] - a[1]
        for t0, t1 in zip(ts[:-1], ts[1:]):
            p0 = (a[0] + t0 * dx, a[1] + t0 * dy); p1 = (a[0] + t1 * dx, a[1] + t1 * dy)
            yield p0, p1, (a, b)

def inter_area(P, Q, detail=False):
    """area esatta di P∩Q (Fraction) via Green sulla 1-catena del bordo"""
    EP, EQ = edges(P), edges(Q)
    tot = F(0); pieces = []
    for p0, p1, (a, b) in _subsegs(EP, _splits(EP, EQ)):
        m = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2); s = pip(m, EQ)
        if s == 'in': inc = True
        elif s == 'on':   # lato collineare: incluso se le due orientazioni sono concordi (interni dallo stesso lato)
            inc = None
            for c, d in EQ:
                if on_seg(p0, c, d) and on_seg(p1, c, d):
                    inc = ((b[0] - a[0]) * (d[0] - c[0]) + (b[1] - a[1]) * (d[1] - c[1])) > 0; break
            assert inc is not None
        else: inc = False
        if inc: tot += cross(p0[0], p0[1], p1[0], p1[1]); pieces.append(('P', s, p0, p1))
    for q0, q1, _ in _subsegs(EQ, _splits(EQ, EP)):
        m = ((q0[0] + q1[0]) / 2, (q0[1] + q1[1]) / 2)
        if pip(m, EP) == 'in': tot += cross(q0[0], q0[1], q1[0], q1[1]); pieces.append(('Q', 'in', q0, q1))
    A = tot / 2
    return (A, pieces) if detail else A

def is_convex_ring(poly):
    P = normalize(poly)
    if len(P) != 1 or len(P[0]) != 1: return False
    r = P[0][0]; n = len(r)
    return all(orient(r[k], r[(k + 1) % n], r[(k + 2) % n]) >= 0 for k in range(n))

def sh_clip_area(P, Qconvex):
    """controprova indipendente: Sutherland–Hodgman esatto di P (un anello) contro Q convesso"""
    out = normalize(P)[0][0]; q = normalize(Qconvex)[0][0]; n = len(q)
    for k in range(n):
        c, d = q[k], q[(k + 1) % n]; inp = out; out = []
        if not inp: break
        for m in range(len(inp)):
            cur, prv = inp[m], inp[m - 1]
            ic, ip = orient(c, d, cur) >= 0, orient(c, d, prv) >= 0
            if ic != ip:
                ex, ey = cur[0] - prv[0], cur[1] - prv[1]
                t = orient(c, d, prv) / (orient(c, d, prv) - orient(c, d, cur))
                out.append((prv[0] + t * ex, prv[1] + t * ey))
            if ic: out.append(cur)
    return ring_area2(out) / 2 if len(out) >= 3 else F(0)

def seg_dist2(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]; L2 = dx * dx + dy * dy
    t = F(0) if L2 == 0 else ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2
    t = min(F(1), max(F(0), t)); qx, qy = a[0] + t * dx - p[0], a[1] + t * dy - p[1]
    return qx * qx + qy * qy
def bnd_dist2(p, poly): return min(seg_dist2(p, a, b) for a, b in edges(poly))
