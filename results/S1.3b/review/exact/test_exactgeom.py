
import sys; sys.path.insert(0, '.')
from exactgeom import *
from fractions import Fraction as F
import random
def sq(x0, y0, x1, y1): return [[[(F(x0), F(y0)), (F(x1), F(y0)), (F(x1), F(y1)), (F(x0), F(y1))]]]
e = F(2) ** -52
cases = [
 ("quadrati adiacenti (lato comune)", sq(0,0,1,1), sq(1,0,2,1), F(0)),
 ("quadrati adiacenti per vertice", sq(0,0,1,1), sq(1,1,2,2), F(0)),
 ("sovrapposti [0,2]^2 ∩ [1,3]^2", sq(0,0,2,2), sq(1,1,3,3), F(1)),
 ("identici", sq(0,0,2,2), sq(0,0,2,2), F(4)),
 ("identici, orientazione opposta", sq(0,0,2,2), [[sq(0,0,2,2)[0][0][::-1]]], F(4)),
 ("annidati", sq(0,0,4,4), sq(1,1,2,3), F(2)),
 ("lato parziale comune, interni opposti", sq(0,0,1,1), sq(1,F(1,2),3,3), F(0)),
 ("lato parziale comune, interni concordi", sq(0,0,2,2), sq(1,0,3,1), F(1)),
 ("scheggia 1 ulp: vertice (1-e,1)", sq(0,0,1,1), [[[(F(1),F(0)),(F(2),F(0)),(F(2),F(1)),(1-e,F(1))]]], e/2),
 ("gap 1 ulp: vertice (1+e,1)", sq(0,0,1,1), [[[(F(1),F(0)),(F(2),F(0)),(F(2),F(1)),(1+e,F(1))]]], F(0)),
 ("L concava ∩ quadrato", [[[(F(0),F(0)),(F(3),F(0)),(F(3),F(1)),(F(1),F(1)),(F(1),F(3)),(F(0),F(3))]]], sq(F(1,2),F(1,2),2,2), F(1,2)*F(1,2)+F(1,2)*1+1*F(1,2)),
 ("quadrato con buco ∩ quadrato", [[sq(0,0,4,4)[0][0], sq(1,1,3,3)[0][0]]], sq(0,0,2,2), F(3)),
 ("multipoligono ∩ striscia", [sq(0,0,1,1)[0], sq(2,0,3,1)[0]], sq(F(1,2),0,F(5,2),1), F(1)),
]
ok = True
for name, P, Q, exp in cases:
    a1 = inter_area(P, Q); a2 = inter_area(Q, P)
    sh = sh_clip_area(P, Q) if (len(P) == 1 and len(P[0]) == 1 and is_convex_ring(Q)) else None
    good = a1 == exp and a2 == exp and (sh is None or sh == exp); ok &= good
    print(("PASS" if good else "FAIL"), name, "atteso", exp, "PQ", a1, "QP", a2, "SH", sh)
# casuali vs shapely (non degeneri): poligoni stellati con coordinate double
try:
    import shapely
    from shapely.geometry import Polygon
    import math
    random.seed(1); worst = 0.0
    for it in range(300):
        def star(cx, cy):
            k = random.randint(3, 12); angs = sorted(random.uniform(0, 2*math.pi) for _ in range(k))
            return [(cx + random.uniform(.2, 1)*math.cos(t), cy + random.uniform(.2, 1)*math.sin(t)) for t in angs]
        a = star(0, 0); b = star(random.uniform(-1, 1), random.uniform(-1, 1))
        A, B = Polygon(a), Polygon(b)
        if not (A.is_valid and B.is_valid): continue
        ex = float(inter_area([[[(F(x), F(y)) for x, y in a]]], [[[(F(x), F(y)) for x, y in b]]]))
        worst = max(worst, abs(ex - A.intersection(B).area))
    print(("PASS" if worst < 1e-12 else "FAIL"), "300 coppie stellate casuali vs shapely: max |diff| =", worst); ok &= worst < 1e-12
except ImportError: print("shapely assente")
print("TUTTI PASS" if ok else "QUALCHE FAIL")
