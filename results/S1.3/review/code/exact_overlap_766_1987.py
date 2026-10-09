#!/usr/bin/env python3
"""results/S1.3/review/code/exact_overlap_766_1987.py — compito 5: area ESATTA (aritmetica razionale) dell'intersezione
fra i territori 766 e 1987 di A3 r3 CSR replica 3 (coordinate double esportate da tessellate_voronoi(), WKB -> repr Python).
Sutherland-Hodgman razionale: 766 ritagliato dai semipiani del pentagono convesso 1987."""
from fractions import Fraction as F
a = [[973.3384350022055, 18.85435339788798], [1000.1043756993286, 2.841546762068755], [1000.1043756993286, 0.0], [970.7050669390757, 0.0], [969.4090590704214, 9.895122257082992]]
b = [[969.6669003931302, 32.00545838482873], [970.2013793599356, 32.16445884749549], [1000.1043756993286, 27.32101704924243], [1000.1043756993286, 2.841546762068754], [973.3384350022055, 18.85435339788798]]
A = [(F(x), F(y)) for x, y in a]; B = [(F(x), F(y)) for x, y in b]
def area(P): return sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1] for i in range(len(P))) / 2
def convex(P):
    s = [ (P[(i+1)%len(P)][0]-P[i][0])*(P[(i+2)%len(P)][1]-P[(i+1)%len(P)][1]) - (P[(i+1)%len(P)][1]-P[i][1])*(P[(i+2)%len(P)][0]-P[(i+1)%len(P)][0]) for i in range(len(P))]
    return all(v >= 0 for v in s) or all(v <= 0 for v in s)
sB = 1 if area(B) > 0 else -1
def clip(P, a0, a1):
    out = []
    side = lambda p: sB * ((a1[0] - a0[0]) * (p[1] - a0[1]) - (a1[1] - a0[1]) * (p[0] - a0[0]))   # >= 0 dentro
    for i in range(len(P)):
        p, q = P[i], P[(i + 1) % len(P)]; fp, fq = side(p), side(q)
        if fp >= 0: out.append(p)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq); out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out
P = A
for i in range(len(B)):
    P = clip(P, B[i], B[(i + 1) % len(B)])
    if len(P) < 3: break
ov = abs(area(P)) if len(P) >= 3 else F(0)
print("1987 convesso:", convex(B), " area 766 esatta:", float(abs(area(A))), " area 1987 esatta:", float(abs(area(B))))
print("area ESATTA dell'intersezione 766 ∩ 1987 (µm²):", float(ov), " vertici:", len(P))
