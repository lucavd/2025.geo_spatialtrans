"""tools/sonda/S3b_fig_halves.py - S3b: figura per la decisione. Spazio array di H2 (um): tessuto di cui S3 ha letto l'espressione in H2
(metà di sviluppo di H2) e in H1 (metà di sviluppo di H1 riportata su H2 con la registrazione globale), aree di S3 delle due sezioni; residuo locale."""
import os, sys, json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap; from matplotlib.patches import Patch, Polygon
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
E = "/mnt/micron/geo_spatialtrans/S3b_sonda/explore"; I = json.load(open(f"{E}/overlap.json")); Wh = np.array(I["H2hires_to_H1hires_colrow"]); Wi = np.linalg.inv(Wh)
co = {s: np.array(I[s]["array_to_fullres"]).T for s in ("H1", "H2")}; sf = {s: I[s]["hires_scalef"] for s in ("H1", "H2")}; m = {s: I[s]["m_um"] for s in ("H1", "H2")}
def a2h(s, r, c): P = np.c_[r, c, np.ones(len(r))] @ co[s]; return np.c_[P[:, 1], P[:, 0]] * sf[s]
def h2a(s, xy): A = co[s][:2].T; b = co[s][2]; return np.linalg.solve(A, (xy[:, ::-1] / sf[s] - b).T).T
M = {}
for s in ("H1", "H2"):
    tp = C.positions(s, "008"); R, Cc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1; Z = np.zeros((R, Cc), np.uint8); Z[tp.array_row, tp.array_col] = tp.in_tissue; M[s] = Z
R, Cc = M["H2"].shape; rr, cc = [a.ravel() for a in np.meshgrid(np.arange(R), np.arange(Cc), indexing="ij")]
xy1 = (Wh @ np.c_[a2h("H2", rr, cc), np.ones(len(rr))].T).T[:, :2]; rc1 = h2a("H1", xy1); ri, ci = np.rint(rc1[:, 0]).astype(int), np.rint(rc1[:, 1]).astype(int)
ok = (ri >= 0) & (ri < M["H1"].shape[0]) & (ci >= 0) & (ci < M["H1"].shape[1]); t1 = np.zeros(len(rr), bool); t1[ok] = M["H1"][ri[ok], ci[ok]] == 1
t2 = M["H2"][rr, cc] == 1; d2 = (cc + 0.5) * 8 < m["H2"] - 50; d1 = (rc1[:, 1] + 0.5) * 8 < m["H1"] - 50
Z = np.zeros(len(rr)); Z[t2] = 1; Z[t2 & d2] = 2; Z[t2 & t1 & d1] = 3; Z[t2 & t1 & d1 & d2] = 4
fig, ax = plt.subplots(1, 2, figsize=(13, 7.4), gridspec_kw=dict(width_ratios=[1.35, 1])); fig.patch.set_facecolor("white")
cm = ListedColormap(["white", "0.8", "#9ecae1", "#fdae6b", "k"])
ax[0].imshow(Z.reshape(R, Cc), cmap=cm, vmin=0, vmax=4, extent=(0, Cc * 8, R * 8, 0), interpolation="nearest")
for s, col in (("H2", "#08519c"), ("H1", "#a63603")):
    for a in C.areas(s).itertuples():
        sq = np.array([[a.x0, a.y0], [a.x0 + a.side, a.y0], [a.x0 + a.side, a.y0 + a.side], [a.x0, a.y0 + a.side]])
        if s == "H1":
            xy = a2h("H1", sq[:, 1] / 8 - 0.5, sq[:, 0] / 8 - 0.5); xy2 = (Wi @ np.c_[xy, np.ones(4)].T).T[:, :2]; rc = h2a("H2", xy2); sq = np.c_[(rc[:, 1] + 0.5) * 8, (rc[:, 0] + 0.5) * 8]
        ax[0].add_patch(Polygon(sq, closed=True, fill=False, ec=col, lw=1.6, ls="-" if a.tissue == "tumour" else "--"))
ax[0].set_xlabel("x in H2 array frame (um)"); ax[0].set_ylabel("y (um)"); ax[0].set_title("Development halves of H1 and H2 in a common frame")
ax[0].legend(handles=[Patch(fc="#9ecae1", label="expression read in H2 (H2 dev half)"), Patch(fc="#fdae6b", label="expression read in H1 (H1 dev half)"),
    Patch(fc="0.8", label="tissue in neither dev half (thin strip at the split)"), Patch(fc="none", ec="#08519c", label="S3 areas, H2 (solid tumour, dashed stroma)"),
    Patch(fc="none", ec="#a63603", label="S3 areas, H1, mapped onto H2")], loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=7.5, frameon=False)
a_ = 64e-6; ax[0].text(0.02, 0.02, f"tissue in both dev halves: {(t2 & t1 & d1 & d2).sum() * a_:.2f} mm2 of {(t2 & t1).sum() * a_:.1f} mm2", transform=ax[0].transAxes, fontsize=8.5, bbox=dict(fc="white", ec="0.6"))
L = pd.read_csv(f"{E}/local_ecc.csv"); ax[1].hist(L.norm_2, bins=20, color="0.45", ec="white"); ax[1].axvline(L.norm_2.quantile(0.95), color="#a63603", ls="--", label=f"95th pct {L.norm_2.quantile(0.95):.0f} um")
ax[1].set_xlabel("local residual after global H1->H2 registration (um, 1 mm tiles)"); ax[1].set_ylabel("tiles"); ax[1].set_title(f"Registration H1->H2, {len(L)} tiles"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{E}/S3b_halves.png", dpi=130, facecolor="white"); print("ok", Z.shape)
