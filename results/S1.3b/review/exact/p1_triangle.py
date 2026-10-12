
import sys; sys.path.insert(0,'.')
exec(open('p1_exact.py').read().split("geo, T, nbad")[0])
geo, T, nb = load(f"{E}/P1_pairs_coords_hex.csv")
import json
rows = json.load(open(f"{E}/P1_exact_results.json"))
for x in rows:
    Pi, Pj = tofrac(geo[(x['case'],'i')]), tofrac(geo[(x['case'],'j')])
    Vi = {v for p in Pi for r in p for v in strip_close(r)}; Vj = {v for p in Pj for r in p for v in strip_close(r)}
    (a,) = tuple(Vi & Vj); d = x['diff'][0]
    vi = (F(float.fromhex(d['vi'][0])), F(float.fromhex(d['vi'][1]))); vj = (F(float.fromhex(d['vj'][0])), F(float.fromhex(d['vj'][1])))
    tri = abs(orient(a, vi, vj)) / 2; A = F(x['inter_exact'])
    print(x['case'], "triangolo(a,vi,vj) == intersezione esatta:", tri == A, "%.4g" % float(A), "lato condiviso %.3f" % float(((a[0]-vi[0])**2+(a[1]-vi[1])**2))**0.5)
