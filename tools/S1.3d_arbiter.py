"""tools/S1.3d_arbiter.py — S1.3d (BL-084): verità esatta (razionali) per O1 e O2.
Pre-registrazione results/S1.3d/S1.3d_preregistration.md (ee6208b); definizioni operative in results/S1.3d/S1.3d_dev_notes.md.
Riusa SENZA modifiche results/S1.3b/review/exact/exactgeom.py (validato in S1.3b; md5 registrato in dev_notes).
Nuovo qui: (1) intersezione esatta cella ∩ regione con indice dei lati della regione (Green sulla 1-catena, come exactgeom.inter_area);
(2) Sutherland–Hodgman esatto per anello contro bin quadrati allineati agli assi (O2). Validati in C-2 contro exactgeom.inter_area.
Uso: .venv/bin/python tools/S1.3d_arbiter.py <geomdir> <truthdir> <key> [--o2]   (key = A_roi_model_rep)"""
import sys, os, json, time, math
from fractions import Fraction as F
import numpy as np, pyarrow.parquet as pq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "results", "S1.3b", "review", "exact"))
import exactgeom as eg

# ---- lettura -------------------------------------------------------------------------------------------------
def polys_from_table(t):
    """parquet (cell, part, ring, k, x, y) -> {cell: [[ring(F,F)...]...]} con coordinate esatte dai double"""
    c, p, r, k = (t.column(n).to_numpy() for n in ("cell", "part", "ring", "k"))
    x, y = t.column("x").to_numpy(), t.column("y").to_numpy()
    o = np.lexsort((k, r, p, c)); out = {}
    for idx in o:
        ci, pi, ri = int(c[idx]), int(p[idx]), int(r[idx])
        parts = out.setdefault(ci, {}); rings = parts.setdefault(pi, {}); rings.setdefault(ri, []).append((F(float(x[idx])), F(float(y[idx]))))
    return {ci: [[parts[pi][ri] for ri in sorted(parts[pi])] for pi in sorted(parts)] for ci, parts in out.items()}

def fbbox(poly):
    xs = [float(q[0]) for part in poly for r in part for q in r]; ys = [float(q[1]) for part in poly for r in part for q in r]
    return min(xs), max(xs), min(ys), max(ys)

# ---- regione indicizzata -----------------------------------------------------------------------------------------
class RegionIndex:
    def __init__(self, reg):
        self.reg = reg; self.E = eg.edges(reg)
        a = np.array([[float(e[0][0]), float(e[1][0]), float(e[0][1]), float(e[1][1])] for e in self.E])
        self.xmin, self.xmax = np.minimum(a[:, 0], a[:, 1]), np.maximum(a[:, 0], a[:, 1])
        self.ymin, self.ymax = np.minimum(a[:, 2], a[:, 3]), np.maximum(a[:, 2], a[:, 3])
    def local(self, bb, m=1e-9):
        x0, x1, y0, y1 = bb
        q = np.nonzero((self.xmax >= x0 - m) & (self.xmin <= x1 + m) & (self.ymax >= y0 - m) & (self.ymin <= y1 + m))[0]
        return [self.E[i] for i in q]
    def pip(self, p):
        """identico a eg.pip(p, E) sui soli lati della striscia orizzontale (gli altri sono saltati da eg.pip comunque)"""
        py = float(p[1]); m = 1e-9 * max(1.0, abs(py))
        q = np.nonzero((self.ymax >= py - m) & (self.ymin <= py + m))[0]
        return eg.pip(p, [self.E[i] for i in q])

def inside_area(cell, RI):
    """area esatta di cell ∩ regione: stessa 1-catena di eg.inter_area(cell, regione), con i lati della regione localizzati"""
    EP = eg.edges(cell); bb = fbbox(cell); EQ = RI.local(bb)
    tot = F(0)
    for p0, p1, (a, b) in eg._subsegs(EP, eg._splits(EP, EQ)):
        mm = ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2); s = RI.pip(mm)
        if s == 'in': inc = True
        elif s == 'on':
            inc = None
            for c, d in EQ:
                if eg.on_seg(p0, c, d) and eg.on_seg(p1, c, d):
                    inc = ((b[0] - a[0]) * (d[0] - c[0]) + (b[1] - a[1]) * (d[1] - c[1])) > 0; break
            assert inc is not None
        else: inc = False
        if inc: tot += eg.cross(p0[0], p0[1], p1[0], p1[1])
    if EQ:
        for q0, q1, _ in eg._subsegs(EQ, eg._splits(EQ, EP)):
            mm = ((q0[0] + q1[0]) / 2, (q0[1] + q1[1]) / 2)
            if eg.pip(mm, EP) == 'in': tot += eg.cross(q0[0], q0[1], q1[0], q1[1])
    return tot / 2, len(EQ)

# ---- Sutherland–Hodgman esatto contro bin assiali -----------------------------------------------------------------
def _clip_half(pts, axis, val, keep_ge):
    val = F(val)                                   # mai float nell'aritmetica (Fraction - float = float)
    out = []; n = len(pts)
    if not n: return out
    for k in range(n):
        cur, prv = pts[k], pts[k - 1]
        ic = (cur[axis] >= val) if keep_ge else (cur[axis] <= val)
        ip = (prv[axis] >= val) if keep_ge else (prv[axis] <= val)
        if ic != ip:
            t = (val - prv[axis]) / (cur[axis] - prv[axis])
            q = [None, None]; q[axis] = F(val); o = 1 - axis; q[o] = prv[o] + t * (cur[o] - prv[o]); out.append((q[0], q[1]))
        if ic: out.append(cur)
    return out
def sh_ring_box(ring, x0, x1, y0, y1):
    p = _clip_half(ring, 0, x0, True); p = _clip_half(p, 0, x1, False); p = _clip_half(p, 1, y0, True); p = _clip_half(p, 1, y1, False)
    return p
def piece_box(cell_n, x0, x1, y0, y1):
    """area esatta (Fraction) e perimetro (float, anelli S-H) di cell ∩ [x0,x1]x[y0,y1]; cell_n già normalizzata (esterni CCW, buchi CW)"""
    A = F(0); per = 0.0
    for part in cell_n:
        for r in part:
            p = sh_ring_box(r, x0, x1, y0, y1)
            if len(p) >= 3:
                A += eg.ring_area2(p) / 2
                per += sum(math.hypot(float(p[(k + 1) % len(p)][0] - p[k][0]), float(p[(k + 1) % len(p)][1] - p[k][1])) for k in range(len(p)))
    return A, per

def o2_truth(cells, gs, s):
    X0, Y0 = gs["X0"], gs["Y0"]; nx = int(gs["nx8"] * 8 // s); ny = int(gs["ny8"] * 8 // s)
    rows = []
    for ci, poly in cells.items():
        pn = eg.normalize(poly); x0, x1, y0, y1 = fbbox(poly)
        ia, ib = max(0, int(math.floor((x0 - X0) / s)) - 1), min(nx - 1, int(math.floor((x1 - X0) / s)) + 1)
        ja, jb = max(0, int(math.floor((y0 - Y0) / s)) - 1), min(ny - 1, int(math.floor((y1 - Y0) / s)) + 1)
        for ix in range(ia, ib + 1):
            for iy in range(ja, jb + 1):
                a, per = piece_box(pn, X0 + s * ix, X0 + s * (ix + 1), Y0 + s * iy, Y0 + s * (iy + 1))
                if a != 0: rows.append((ci, ix, iy, a, per))
    return rows

# ---- O1 -----------------------------------------------------------------------------------------------------------
def o1_truth(cells, reg, pairs, n):
    t0 = time.time(); RI = RegionIndex(reg)
    a = {i: eg.area(cells[i]) if i in cells else F(0) for i in range(1, n + 1)}
    A_reg = eg.area(reg)
    pa = [eg.inter_area(cells[int(i)], cells[int(j)]) for i, j in pairs]
    ins = {}; n_loc = 0
    for i in range(1, n + 1):
        if i not in cells: ins[i] = F(0); continue
        bb = fbbox(cells[i]); loc = RI.local(bb)
        if loc: ins[i], _ = inside_area(cells[i], RI); n_loc += 1
        else:                                      # nessun lato della regione nel riquadro: tutta dentro o tutta fuori
            v = cells[i][0][0][0]; ins[i] = a[i] if RI.pip(v) == 'in' else F(0)
    S = sum(a.values()); Sov = sum(pa); Iin = sum(ins.values())
    U = S - Sov; I = Iin - Sov                     # approssimazione: sovrapposizioni triple trascurate (limite = Sov)
    SD = (U - I) + (A_reg - I)
    return dict(sum_a=float(S), A_reg=float(A_reg), sum_ov=float(Sov), union=float(U), inside=float(I), symd=float(SD),
                approx_bound=float(Sov), n_local=n_loc, t_s=time.time() - t0), a, pa, ins

def write_parquet(path, cols):
    import pyarrow as pa_
    pq.write_table(pa_.table(cols), path)

def main_chunk(geomdir, truthdir, key, k, K, do_o1, do_o2):
    """blocco k di K (celle con cell % K == k, coppie con indice % K == k); riunione con tools/S1.3d_arbiter_merge.py"""
    meta = json.load(open(os.path.join(geomdir, key + "_meta.json")))
    cells = polys_from_table(pq.read_table(os.path.join(geomdir, key + "_cells.parquet")))
    reg = polys_from_table(pq.read_table(os.path.join(geomdir, key + "_region.parquet")))[0]
    mine = {c: p for c, p in cells.items() if c % K == k}; tag = f"{key}_chunk{k:03d}of{K:03d}"
    if do_o1:
        pt = pq.read_table(os.path.join(geomdir, key + "_pairs.parquet")); I, J = pt.column("i").to_pylist(), pt.column("j").to_pylist()
        sel = [q for q in range(len(I)) if q % K == k]; RI = RegionIndex(reg)
        pa = [eg.inter_area(cells[I[q]], cells[J[q]]) for q in sel]
        write_parquet(os.path.join(truthdir, tag + "_o1_pairs.parquet"), {"q": sel, "i": [I[q] for q in sel], "j": [J[q] for q in sel], "a": [float(v) for v in pa], "zero": [v == 0 for v in pa]})
        ids = sorted(mine); a = [eg.area(mine[c]) for c in ids]; ins = []
        for c in ids:
            loc = RI.local(fbbox(mine[c]))
            if loc: ins.append(inside_area(mine[c], RI)[0])
            else: ins.append(eg.area(mine[c]) if RI.pip(mine[c][0][0][0]) == 'in' else F(0))
        write_parquet(os.path.join(truthdir, tag + "_o1_cells.parquet"), {"cell": ids, "a": [float(v) for v in a], "inside": [float(v) for v in ins]})
        if k == 0: json.dump({"A_reg": float(eg.area(reg))}, open(os.path.join(truthdir, key + "_o1_areg.json"), "w"))
    if do_o2:
        for s_ in (2, 8):
            rows = o2_truth(mine, meta["gs"], s_)
            write_parquet(os.path.join(truthdir, f"{tag}_o2_s{s_}.parquet"), {"i": [r[0] for r in rows], "ix": [r[1] for r in rows], "iy": [r[2] for r in rows],
                          "area": [float(r[3]) for r in rows], "perim": [r[4] for r in rows]})

def main():
    if "--chunk" in sys.argv:
        k, K = map(int, sys.argv[sys.argv.index("--chunk") + 1].split("/"))
        return main_chunk(sys.argv[1], sys.argv[2], sys.argv[3], k, K, "--no-o1" not in sys.argv, "--o2" in sys.argv)
    geomdir, truthdir, key = sys.argv[1:4]; do_o2 = "--o2" in sys.argv; do_o1 = "--no-o1" not in sys.argv
    os.makedirs(truthdir, exist_ok=True)
    meta = json.load(open(os.path.join(geomdir, key + "_meta.json")))
    cells = polys_from_table(pq.read_table(os.path.join(geomdir, key + "_cells.parquet")))
    reg = polys_from_table(pq.read_table(os.path.join(geomdir, key + "_region.parquet")))[0]
    out = {"key": key, "n": meta["n"]}
    if do_o1:
        pt = pq.read_table(os.path.join(geomdir, key + "_pairs.parquet")); pairs = list(zip(pt.column("i").to_pylist(), pt.column("j").to_pylist()))
        summ, a, pa, ins = o1_truth(cells, reg, pairs, meta["n"]); out["o1"] = summ
        write_parquet(os.path.join(truthdir, key + "_o1_pairs.parquet"), {"i": [p[0] for p in pairs], "j": [p[1] for p in pairs],
                      "a": [float(v) for v in pa], "zero": [v == 0 for v in pa]})
        ids = list(range(1, meta["n"] + 1))
        write_parquet(os.path.join(truthdir, key + "_o1_cells.parquet"), {"cell": ids, "a": [float(a[i]) for i in ids], "inside": [float(ins[i]) for i in ids]})
    if do_o2:
        t0 = time.time(); r2 = o2_truth(cells, meta["gs"], 2); r8 = o2_truth(cells, meta["gs"], 8)
        # identità interna: somma esatta dei 16 sotto-pezzi = pezzo esatto da 8 µm
        agg = {}
        for ci, ix, iy, a_, _ in r2: agg[(ci, ix // 4, iy // 4)] = agg.get((ci, ix // 4, iy // 4), F(0)) + a_
        d8 = {(ci, ix, iy): a_ for ci, ix, iy, a_, _ in r8}
        keys = set(agg) | set(d8); mism = sum(1 for k_ in keys if agg.get(k_, F(0)) != d8.get(k_, F(0)))
        for s_, rows in ((2, r2), (8, r8)):
            write_parquet(os.path.join(truthdir, f"{key}_o2_s{s_}.parquet"), {"i": [r[0] for r in rows], "ix": [r[1] for r in rows], "iy": [r[2] for r in rows],
                          "area": [float(r[3]) for r in rows], "perim": [r[4] for r in rows]})
        out["o2"] = dict(n2=len(r2), n8=len(r8), agg_mismatch=mism, n_keys8=len(keys), t_s=time.time() - t0)
    json.dump(out, open(os.path.join(truthdir, key + ("_o1" if do_o1 else "") + ("_o2" if do_o2 else "") + ".json"), "w"))
    print(json.dumps(out))

if __name__ == "__main__":
    main()
