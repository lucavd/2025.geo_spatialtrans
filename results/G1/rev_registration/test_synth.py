import numpy as np, time
from scipy.ndimage import gaussian_filter
import g1_tiles as T
rng = np.random.default_rng(1)
def synth(true_s, seed):
    r = np.random.default_rng(seed)
    n = 250; G = 30                     # bins (2um), region with margin
    # cells: random centres, each with a type -> gene profile
    nc = 2500; cx = r.uniform(0, 2*n, nc); cy = r.uniform(0, 2*n, nc); ty = r.integers(0, 5, nc)
    prof = r.gamma(0.5, 1, (5, G)); prof /= prof.sum(1, keepdims=True)
    # transcripts per cell ~ 60, spread sigma 3 um
    k = r.poisson(60, nc); ci = np.repeat(np.arange(nc), k)
    g = np.array([r.choice(G, p=prof[ty[i]]) for i in ci])
    ux = cx[ci] + r.normal(0, 3, len(ci)); uy = cy[ci] + r.normal(0, 3, len(ci))
    # Visium: thinned copy of same molecules (other draw), binned at 2um
    k2 = r.poisson(25, nc); ci2 = np.repeat(np.arange(nc), k2)
    g2 = np.array([r.choice(G, p=prof[ty[i]]) for i in ci2])
    vx = cx[ci2] + r.normal(0, 3, len(ci2)); vy = cy[ci2] + r.normal(0, 3, len(ci2))
    vc = np.floor(vx/2 + 0.5).astype(int); vr = np.floor(vy/2 + 0.5).astype(int)
    # Xenium shifted by true_s (Xenium = Visium position + s)
    return g, ux + true_s[0], uy + true_s[1], g2, vc, vr, np.ones(len(g2)), G
for s in [(0, 0), (3.3, -1.7), (7.6, 9.1)]:
    g, ux, uy, vg, vc, vr, vn, G = synth(s, 5)
    sc = (int(2*round(s[0]/2)), int(2*round(s[1]/2)))
    t0 = time.time()
    t = T.Tile(g, ux, uy, vg, vc, vr, vn, 75, 75, 100, sc, G)
    out = t.search()
    print("true", s, "coarse", sc, "est", (sc[0]+out["ox"], sc[1]+out["oy"]), "edge", out["edge"], "genes", t.G, "valid", t.valid, round(time.time()-t0, 1), "s")
# wrong tissue
g, ux, uy, *_ = synth((0, 0), 5); _, _, _, vg, vc, vr, vn, G = synth((0, 0), 99)
t = T.Tile(g, ux, uy, vg, vc, vr, vn, 75, 75, 100, (0, 0), G); o = t.search()
print("wrong tissue est", o["ox"], o["oy"], "edge", o["edge"], "peak", round(o["peak"], 3))
