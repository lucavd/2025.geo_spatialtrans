"""tools/S2.6_PH3_edge.py - S2.6, analisi post hoc PH-3 (dichiarata dopo i risultati): vetro a pezzi / pezzo solo separando i blocchi di vetro
vicini al bordo dell'area di cattura (centro a 50-300 um dal bordo) da quelli lontani (oltre 300 um), stesse fasce di distanza dal tessuto.
Ipotesi: nei pezzi la seconda coda si riflette anche sul bordo dell'area di cattura (frame_open = TRUE), nel pezzo solo si perde.
Uso: .venv/bin/python tools/S2.6_PH3_edge.py [out.csv]"""
import sys, json, numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy import ndimage
SEC = {"A1": "Visium_HD_Mouse_Small_Intestine", "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A6": "Visium_HD_6p5mm_Human_Heart"}
M = "/mnt/micron/geo_spatialtrans/S2.6/q3"; rows = []
def blocks(tag, NR, NC):
    j = json.load(open(f"{M}/{tag}/run.json")); g = j["grid"]
    u2 = np.fromfile(f"{M}/{tag}/u2.bin", dtype="<i4", count=g["n_row"] * g["n_col"]).reshape(g["n_row"], g["n_col"])
    b8 = u2.reshape(g["n_row"] // 4, 4, g["n_col"] // 4, 4).sum(axis=(1, 3)); o = int(round(-g["x0"] / 8))
    return b8[o:o + NR, o:o + NC].astype(float)
for A, s in SEC.items():
    tp = pq.read_table(f"/mnt/micron/geo_spatialtrans/data_real/{s}/binned_outputs_x/binned_outputs/square_008um/spatial/tissue_positions.parquet", columns=["in_tissue", "array_row", "array_col"]).to_pandas()
    NR, NC = tp.array_row.max() + 1, tp.array_col.max() + 1
    tis = np.zeros((NR, NC), bool); tis[tp.array_row[tp.in_tissue == 1], tp.array_col[tp.in_tissue == 1]] = True
    d = ndimage.distance_transform_edt(~tis) * 8
    yc = (np.arange(NR) + 0.5) * 8; xc = (np.arange(NC) + 0.5) * 8
    de = np.minimum(np.minimum(yc, NR * 8 - yc)[:, None], np.minimum(xc, NC * 8 - xc)[None, :])
    one = blocks(f"one_{A}", NR, NC); pcs = blocks(f"pieces_{A}", NR, NC)
    for zone, zm in [("bordo 50-300 um", (de >= 50) & (de < 300)), ("interno > 300 um", de >= 300)]:
        for lo, hi in [(16, 256), (256, 1e9)]:
            w = (~tis) & zm & (d >= lo) & (d < hi)
            rows.append(dict(archetype=A, zone=zone, glass_band=f"{lo}-{int(hi) if hi < 1e8 else 'inf'}", n_blocks=int(w.sum()),
                             umi_one=one[w].sum(), umi_pieces=pcs[w].sum(), ratio=pcs[w].sum() / max(one[w].sum(), 1e-12)))
out = sys.argv[1] if len(sys.argv) > 1 else "results/S2.6/S2.6_PH3_edge.csv"
pd.DataFrame(rows).to_csv(out, index=False); print(pd.DataFrame(rows).to_string())
