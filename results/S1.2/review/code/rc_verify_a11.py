# rc_verify_a11.py — profilo di densita' vicino ai bordi (200 seed), verifica indipendente
import os, json, math, numpy as np, pandas as pd
OUT = "results/S1.2/review/code/out"; VER = "results/S1.2/review/code/verify"; os.makedirs(VER, exist_ok=True)
rows = []
S = pd.read_csv(os.path.join(OUT, "A11b_square", "centroids_200seeds_near_edges.csv")); ns = S.seed.nunique()
d = np.minimum.reduce([S.x.values, 1000 - S.x.values, S.y.values, 1000 - S.y.values])
ref_k = ((d >= 10) & (d < 30)).sum(); ref_rho = ref_k / (((1000 - 20) ** 2 - (1000 - 60) ** 2) * ns)
for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 7.5), (7.5, 10), (10, 20), (20, 30)]:
    k = ((d >= lo) & (d < hi)).sum(); ar = ((1000 - 2 * lo) ** 2 - (1000 - 2 * hi) ** 2) * ns
    rows.append(dict(config="square_outer_edge", band_um=f"[{lo},{hi})", side="esterno", n=int(k), rel_density=k / ar / ref_rho, se=math.sqrt(k) / ar / ref_rho))
for nm in ["A11b_adjLR", "A11b_adjRL"]:
    D = pd.read_csv(os.path.join(OUT, nm, "centroids_200seeds_near_edges.csv")); D = D[(D.y > 50) & (D.y < 950)]; ns = D.seed.nunique()
    dx = D.x.values - 500; ref = ((np.abs(dx) >= 10) & (np.abs(dx) < 30)).sum() / (2 * 20 * 900 * ns)
    for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 7.5), (7.5, 10), (10, 20)]:
        for side, m in [("sinistra", (dx < -lo) & (dx >= -hi)), ("destra", (dx >= lo) & (dx < hi))]:
            k = m.sum(); ar = (hi - lo) * 900 * ns; rr = sorted(D.region_row.values[m].tolist() and set(D.region_row.values[m].tolist()))
            rows.append(dict(config=nm, band_um=f"[{lo},{hi})", side=f"{side} (riga {rr})", n=int(k), rel_density=k / ar / ref, se=math.sqrt(k) / ar / ref))
T = pd.DataFrame(rows); T.to_csv(os.path.join(VER, "V_A11_boundary_200seeds.csv"), index=False)
print(T.round(4).to_string())
