"""tools/S1.3d_arbiter_check.py — S1.3d: C-2 (due vie esatte indipendenti coincidono come Fraction), CP-4 (placebo e mutante dell'arbitro).
C-2a cella ∩ regione: inside_area (indice dei lati) vs exactgeom.inter_area(cella, regione) con tutti i lati.
C-2b coppie con una cella convessa: exactgeom.inter_area (Green) vs exactgeom.sh_clip_area (Sutherland–Hodgman).
C-2c pezzi O2: piece_box (S-H assiale) vs exactgeom.inter_area(cella, quadrato del bin) (Green).
CP-4 placebo: celle in ordine permutato (input cambiato, asserito) -> aree e coppie esatte identiche dopo la rimappatura.
CP-4 mutante (M-arb): un vertice tolto da una cella -> l'area esatta di quella cella cambia (rilevato).
Uso: .venv/bin/python tools/S1.3d_arbiter_check.py <outdir> <file_chiavi> <m_per_chiave> <out.json> [--cores 40] [--seed 20261014]"""
import sys, os, json, random, argparse, math
from multiprocessing import Pool
from fractions import Fraction as F
import importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("arb", os.path.join(HERE, "S1.3d_arbiter.py")); arb = importlib.util.module_from_spec(spec); spec.loader.exec_module(arb)
eg = arb.eg; pq = arb.pq

def sq(x0, x1, y0, y1): return [[[(F(x0), F(y0)), (F(x1), F(y0)), (F(x1), F(y1)), (F(x0), F(y1))]]]

def one(args):
    outdir, key, m, seed = args; rng = random.Random(f"{seed}-{key}"); g = os.path.join(outdir, "geom")
    meta = json.load(open(os.path.join(g, key + "_meta.json")))
    cells = arb.polys_from_table(pq.read_table(os.path.join(g, key + "_cells.parquet")))
    reg = arb.polys_from_table(pq.read_table(os.path.join(g, key + "_region.parquet")))[0]
    pt = pq.read_table(os.path.join(g, key + "_pairs.parquet")); pairs = list(zip(pt.column("i").to_pylist(), pt.column("j").to_pylist()))
    RI = arb.RegionIndex(reg); res = {"key": key}
    # C-2a
    loc = [c for c in cells if RI.local(arb.fbbox(cells[c]))]; pick = rng.sample(loc, min(m, len(loc)))
    eq = sum(1 for c in pick if arb.inside_area(cells[c], RI)[0] == eg.inter_area(cells[c], reg)); res["c2a"] = [eq, len(pick)]
    # C-2b
    cand = [(i, j) for i, j in pairs if eg.is_convex_ring(cells[j]) and len(eg.normalize(cells[i])) == 1 and len(eg.normalize(cells[i])[0]) == 1]
    pick = rng.sample(cand, min(m, len(cand)))
    eq = sum(1 for i, j in pick if eg.inter_area(cells[i], cells[j]) == eg.sh_clip_area(cells[i], cells[j])); res["c2b"] = [eq, len(pick)]
    # C-2c
    gs = meta["gs"]; s = 2; X0, Y0 = gs["X0"], gs["Y0"]; tot = 0; eq = 0
    for c in rng.sample(sorted(cells), min(m, len(cells))):
        x0, x1, y0, y1 = arb.fbbox(cells[c]); ix = int(math.floor(((x0 + x1) / 2 - X0) / s)); iy = int(math.floor(((y0 + y1) / 2 - Y0) / s))
        for dx, dy in ((0, 0), (1, 0), (0, 1)):
            bx0, by0 = X0 + s * (ix + dx), Y0 + s * (iy + dy)
            a1, _ = arb.piece_box(eg.normalize(cells[c]), bx0, bx0 + s, by0, by0 + s); a2 = eg.inter_area(cells[c], sq(bx0, bx0 + s, by0, by0 + s))
            tot += 1; eq += (a1 == a2)
    res["c2c"] = [eq, tot]
    # CP-4 placebo: permutazione dell'ordine delle celle
    ids = sorted(cells); perm = ids[:]; rng.shuffle(perm)
    if perm == ids: perm = ids[::-1]                     # permutazione identica (casi da 2-3 celle): ordine inverso
    changed = perm != ids
    newc = {k + 1: cells[perm[k]] for k in range(len(perm))}; back = {perm[k]: k + 1 for k in range(len(perm))}
    sp = rng.sample(pairs, min(m, len(pairs)))
    same_a = all(eg.area(newc[back[c]]) == eg.area(cells[c]) for c in rng.sample(ids, min(m, len(ids))))
    same_p = all(eg.inter_area(newc[back[i]], newc[back[j]]) == eg.inter_area(cells[i], cells[j]) for i, j in sp)
    res["cp4_placebo"] = {"input_changed": changed, "identical": bool(same_a and same_p)}
    # CP-4 mutante: un vertice tolto dalla cella con più vertici (anello esterno)
    c = max(ids, key=lambda q: len(cells[q][0][0])); mut = [[list(r) for r in part] for part in cells[c]]; del mut[0][0][len(mut[0][0]) // 2]
    res["cp4_mutant"] = {"cell": c, "detected": eg.area(mut) != eg.area(cells[c])}
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("keys"); ap.add_argument("m", type=int); ap.add_argument("out")
    ap.add_argument("--cores", type=int, default=40); ap.add_argument("--seed", default="20261014"); a = ap.parse_args()
    keys = [k.strip() for k in open(a.keys) if k.strip()]
    with Pool(a.cores) as p: R = p.map(one, [(a.outdir, k, a.m, a.seed) for k in keys])
    tot = {t: [sum(r[t][0] for r in R), sum(r[t][1] for r in R)] for t in ("c2a", "c2b", "c2c")}
    tot["cp4_placebo"] = [sum(r["cp4_placebo"]["identical"] and r["cp4_placebo"]["input_changed"] for r in R), len(R)]
    tot["cp4_mutant"] = [sum(r["cp4_mutant"]["detected"] for r in R), len(R)]
    json.dump({"totals": tot, "per_key": R}, open(a.out, "w")); print(json.dumps(tot))
