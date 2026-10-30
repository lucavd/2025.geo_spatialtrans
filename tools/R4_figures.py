"""tools/R4_figures.py — figure del report R4. Uso: .venv/bin/python tools/R4_figures.py [outdir_run]  -> results/R4/figures/*.png
Legge results/R4/*.csv e, per le distribuzioni per bin, <outdir_run>/roi/*_bins.npz."""
import sys, os, glob
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RUN = sys.argv[1] if len(sys.argv) > 1 else "/mnt/micron/geo_spatialtrans/R4/run1"
RES = "results/R4"; FIG = f"{RES}/figures"; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8, "axes.titlesize": 8,
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlelocation": "left", "legend.frameon": False})
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
NAME = {"A1": "mouse intestine", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "mouse brain", "A6": "heart"}
OI = {"blue": "#0072B2", "orange": "#E69F00", "verm": "#D55E00", "green": "#009E73", "grey": "#8C8C8C"}
x = np.arange(6); lab = [NAME[a].replace(" ", "\n", 1) for a in ARCH]
TYPE_EN = {"epiteliale": "epithelial", "immunitario": "immune", "stromale e muscolo liscio": "stromal + smooth muscle", "epiteliale normale": "normal epithelium",
           "epiteliale tumorale": "tumour epithelium", "fibroblasti": "fibroblasts", "endoteliale": "endothelial", "linfociti B": "B cells", "linfociti T": "T cells",
           "plasmacellule": "plasma cells", "mieloidi": "myeloid", "stromali": "stromal", "endoteliali": "endothelial", "neuroni": "neurons", "astrociti": "astrocytes",
           "oligodendrociti": "oligodendrocytes", "microglia": "microglia", "vascolari": "vascular", "cardiomiociti": "cardiomyocytes", "immunitari": "immune"}

def pts(ax, df, xoff, color, label, filled=True, marker="o"):
    df = df.set_index("archetype").reindex(ARCH)
    ax.errorbar(x + xoff, df.estimate, yerr=[df.estimate - df.lo, df.hi - df.estimate], fmt=marker, color=color, ms=5, capsize=2, lw=1,
                mfc=color if filled else "white", label=label)

P = pd.read_csv(f"{RES}/R4_predictions.csv"); PH1 = pd.read_csv(f"{RES}/R4_PH1.csv"); PH2 = pd.read_csv(f"{RES}/R4_PH2.csv")
# F1 nucleo / corona
fig, ax = plt.subplots(figsize=(6.8, 3.1))
ax.axhspan(0.909, 1.1, color="#EEEEEE", zorder=0); ax.axhline(1, color=OI["grey"], lw=0.6, zorder=0)
ax.text(5.45, 1.0, "uniform\n0.909-1.1", fontsize=6.5, color="#555555", va="center", ha="left")
pts(ax, P[P.id == "P-1"], -0.21, OI["blue"], "pre-registered: pooled bins, divided by shifted mask")
pts(ax, P[P.id == "CP-1"], -0.07, OI["blue"], "shifted mask, pooled bins (counter-check, expected 1)", filled=False)
pts(ax, PH1[PH1.quantity == "cellweighted_star"], 0.07, OI["orange"], "post hoc: cell-weighted, divided by shifted mask")
pts(ax, PH1[PH1.quantity == "cellweighted_shift"], 0.21, OI["orange"], "shifted mask, cell-weighted (expected 1)", filled=False, marker="s")
ax.set_xticks(x, lab); ax.set_xlim(-0.5, 6.2); ax.set_ylabel("UMI per 2 µm bin:\nnucleus / its 8 µm ring")
ax.set_title("Nucleus denser than its ring in brain and heart, sparser in colon tumour and lymph node")
ax.legend(loc="upper left", fontsize=6.5); ax.margins(y=0.06); fig.tight_layout(); fig.savefig(f"{FIG}/R4_F1_nucleus_ring.png", dpi=200); plt.close(fig)

# F2 zeri in eccesso
fig, ax = plt.subplots(figsize=(6.8, 3.1))
ax.axhspan(-0.02, 0.02, color="#EEEEEE", zorder=0); ax.axhline(0, color=OI["grey"], lw=0.6)
pts(ax, P[P.id == "P-2"].dropna(subset=["estimate"]), -0.2, OI["blue"], "pre-registered: negative binomial, one mean per gene")
s = PH2[(PH2.kind == "obs") & (PH2.model == "NB with cell size factor, mean >= 2")].dropna(subset=["estimate"])
pts(ax, s, 0.0, OI["orange"], "post hoc: negative binomial with cell size factor")
pts(ax, P[P.id == "CP-2"].dropna(subset=["estimate"]), 0.2, OI["verm"], "positive control: 30 % of non-zero values set to 0", filled=False, marker="v")
ax.text(3, 0.0, "not\nassessable", ha="center", va="center", fontsize=6.5, color="#555555")
ax.set_xticks(x, lab); ax.set_xlim(-0.5, 5.5); ax.set_ylabel("cells with 0 UMI: observed - expected\n(median over type x gene)")
ax.set_title("No excess zeros within cell types once cell size is accounted for, except intestinal subtypes")
ax.legend(loc="upper left", fontsize=6.5); ax.margins(y=0.06); fig.tight_layout(); fig.savefig(f"{FIG}/R4_F2_excess_zeros.png", dpi=200); plt.close(fig)

# F3 UMI per cellula per tipo
T = pd.read_csv(f"{RES}/R4_types.csv")
fig, axs = plt.subplots(1, 6, figsize=(7.4, 3.3), sharey=True)
for ax, a in zip(axs, ARCH):
    t = T[T.archetype == a].sort_values("median_umi_cell", ascending=False).reset_index(drop=True)
    for i, r in t.iterrows():
        c = OI["blue"] if r.in_P4 else OI["grey"]
        ax.errorbar(i, r.median_umi_cell, yerr=[[r.median_umi_cell - r.lo], [r.hi - r.median_umi_cell]], fmt="o", color=c, ms=4, capsize=2, lw=1, mfc=c if r.in_P4 else "white")
    ax.set_yscale("log"); ax.set_xticks(range(len(t)), [TYPE_EN.get(c_, c_) for c_ in t.cell_type], rotation=60, ha="right", fontsize=6)
    ax.set_xlim(-0.6, len(t) - 0.4); ax.set_title(NAME[a], fontsize=7)
    pr = P[(P.id == "P-4") & (P.archetype == a)].iloc[0]
    ax.text(0.97, 0.97, f"max/min\n{pr.estimate:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
axs[0].set_ylabel("median UMI per cell\n(nucleus + 8 µm ring)")
axs[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
fig.suptitle("UMI per cell differ between cell types (filled: types with at least 5 % of interior nuclei)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/R4_F3_umi_per_cell_type.png", dpi=200); plt.close(fig)

# F4 parita' di riga e colonna + ritaglio
W = pd.read_csv(f"{RES}/R4_PH4_parity.csv")
order = ["Visium_HD_Mouse_Small_Intestine", "Visium_HD_Human_Colon_Cancer", "Visium_HD_Human_Lymph_Node_FFPE", "Visium_HD_6p5mm_Mouse_Brain", "Visium_HD_6p5mm_Human_Heart"]
sname = ["mouse\nintestine", "colon cancer\nsection", "lymph\nnode", "mouse\nbrain", "heart"]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axs[0]
for j, ds in enumerate(order):
    sec = W[(W.unit == "section") & (W.dataset == ds)].iloc[0]; rois = W[(W.unit != "section") & (W.dataset == ds)]
    ax.scatter(np.full(len(rois), j - 0.12), rois.col_odd_over_even, s=10, color=OI["blue"], alpha=0.6, lw=0)
    ax.scatter(np.full(len(rois), j + 0.12), rois.row_odd_over_even, s=10, color=OI["orange"], alpha=0.6, lw=0)
    ax.plot([j - 0.22, j - 0.02], [sec.col_odd_over_even] * 2, color=OI["blue"], lw=2); ax.plot([j + 0.02, j + 0.22], [sec.row_odd_over_even] * 2, color=OI["orange"], lw=2)
ax.axhline(1, color=OI["grey"], lw=0.6); ax.set_xticks(range(5), sname, fontsize=6.5); ax.set_ylabel("mean UMI per 2 µm bin, odd / even")
ax.plot([], [], color=OI["blue"], lw=2, label="columns (bar: section; dots: ROIs)"); ax.plot([], [], color=OI["orange"], lw=2, label="rows")
ax.legend(loc="upper right", fontsize=6.5); ax.set_title("Odd rows and columns capture more UMI in 4 of 5 sections"); ax.margins(y=0.08)
ds = "Visium_HD_Human_Colon_Cancer"; b = f"/mnt/micron/geo_spatialtrans/data_real/{ds}/binned_outputs_x/binned_outputs/square_002um"
with h5py.File(f"{b}/filtered_feature_bc_matrix.h5") as f:
    m = f["matrix"]; ip = m["indptr"][:]; dat = m["data"][:]; bc = m["barcodes"][:].astype(str)
cs = np.concatenate([[0], np.cumsum(dat, dtype=np.int64)]); u = cs[ip[1:]] - cs[ip[:-1]]
tp = pq.read_table(f"{b}/spatial/tissue_positions.parquet").to_pandas().merge(pd.DataFrame({"barcode": bc, "umi": u}), on="barcode")
roi = pd.read_csv("results/R2/R2_rois_checked.csv").query("archetype == 'A2' and roi_id == 'r1'").iloc[0]
s_ = tp[(tp.pxl_row_in_fullres >= roi.r0) & (tp.pxl_row_in_fullres < roi.r1) & (tp.pxl_col_in_fullres >= roi.c0) & (tp.pxl_col_in_fullres < roi.c1)]
r0, c0 = int(s_.array_row.median()), int(s_.array_col.median())
cr = s_[(s_.array_row >= r0) & (s_.array_row < r0 + 60) & (s_.array_col >= c0) & (s_.array_col < c0 + 60)]
im = np.full((60, 60), np.nan); im[cr.array_row - r0, cr.array_col - c0] = cr.umi
h = axs[1].imshow(im, cmap="magma", vmin=0, vmax=np.nanpercentile(im, 99), interpolation="nearest")
axs[1].set_xticks([]); axs[1].set_yticks([]); axs[1].set_title("Colon tumour, 120 x 120 µm at 2 µm: stripes")
cb = fig.colorbar(h, ax=axs[1], fraction=0.046, pad=0.02); cb.set_label("UMI per 2 µm bin", fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/R4_F4_parity.png", dpi=200); plt.close(fig)

# F5 distribuzioni per bin a 2 e 8 um (5 ROI raggruppate): curve cumulative
fig, axs = plt.subplots(1, 6, figsize=(7.4, 2.9), sharey=True)
for j, a in enumerate(ARCH):
    fs = sorted(glob.glob(f"{RUN}/roi/{a}_*_bins.npz")); ax = axs[j]
    for key, col, nm in (("umi", OI["blue"], "2 µm"), ("umi8", OI["green"], "8 µm")):
        uu = np.sort(np.concatenate([np.load(f)[key] for f in fs])); yv = np.arange(1, len(uu) + 1) / len(uu)
        ax.step(uu + 1, yv, where="post", color=col, lw=1.2, label=nm)
    ax.set_xscale("log"); ax.set_title(NAME[a], fontsize=7); ax.tick_params(labelsize=6); ax.set_xlabel("UMI per bin + 1", fontsize=7)
axs[0].set_ylabel("fraction of bins\nat or below"); axs[-1].legend(loc="lower right", fontsize=6.5)
fig.suptitle("8 µm bins are almost never empty; 2 µm bins are rarely empty except in the heart\n(UMI per bin in tissue, 5 ROIs per tissue pooled; medians and empty fractions in the reference table)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/R4_F5_umi_per_bin.png", dpi=200); plt.close(fig)
print("figures:", sorted(os.listdir(FIG)))
