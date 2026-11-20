"""tools/riesame_xen_checks.py - Riesame Parte B, braccio Xenium: controlli indipendenti della registrazione e plausibilita' rispetto a V (XD-6
di docs/riesame/DECISIONI_B_XENIUM.md). Solo riportati, nessuna scelta. Per area: densita' dei nuclei Xenium contro i nuclei di V; quota dei
quadratini in una cellula Xenium; UMI per quadratino in_reg (X grezzo contro V, 408 geni); correlazione incrociata (spostamento del picco)
fra le mappe di V e di X a 2 um lisciate (gaussiana 4 um), per i conteggi dei 408 geni e per i quadratini di nucleo; correlazione di Spearman
dei totali per gene (408) fra V e X. Uscita: results/riesame/xen_plausibility.csv. Uso: .venv/bin/python -B tools/riesame_xen_checks.py"""
import json, numpy as np, pandas as pd
from scipy import ndimage, stats
M = "/mnt/micron/geo_spatialtrans"; F = f"{M}/F3.1/main/real"; XB = f"{M}/riesame/B/X"
AREAS = [f"A2_r{i}" for i in range(1, 6)] + [f"A3_r{i}" for i in range(1, 6)]
PANEL = pd.read_csv(f"{M}/riesame/xenium/panel_ids_in_common.txt", header=None)[0].astype(str).tolist()
MAXS = 20
def to_map(b, v):
    r = b.row.to_numpy() - b.row.min(); c = b.col.to_numpy() - b.col.min(); m = np.zeros((r.max() + 1, c.max() + 1)); m[r, c] = v; return m
def xcorr(a, b):
    a = ndimage.gaussian_filter(a, 2.0); b = ndimage.gaussian_filter(b, 2.0); a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
    out = np.zeros((2 * MAXS + 1, 2 * MAXS + 1)); H, W = a.shape
    for i, dr in enumerate(range(-MAXS, MAXS + 1)):
        for j, dc in enumerate(range(-MAXS, MAXS + 1)):
            A = a[max(0, dr):H + min(0, dr), max(0, dc):W + min(0, dc)]; B = b[max(0, -dr):H + min(0, -dr), max(0, -dc):W + min(0, -dc)]
            out[i, j] = np.mean(A * B)
    k = np.unravel_index(out.argmax(), out.shape); return out[MAXS, MAXS], out.max(), (k[0] - MAXS, k[1] - MAXS), out[MAXS + 10, MAXS + 10]
rows = []
for a in AREAS:
    vb = pd.read_parquet(f"{F}/{a}/bins.parquet"); vg = pd.read_csv(f"{F}/{a}/genes.tsv", header=None)[0].astype(str).to_numpy()
    pos = pd.Series(np.arange(len(vg)), index=vg); gi = pos.reindex(PANEL).to_numpy().astype(int)
    VX = pd.read_parquet(f"{F}/{a}/X.parquet"); VX = VX[VX.gene.isin(gi)]; VX["g"] = pd.Series(np.arange(408), index=gi).reindex(VX.gene).to_numpy()
    xb = pd.read_parquet(f"{XB}/{a}/bins.parquet"); XX = pd.read_parquet(f"{XB}/{a}/X.parquet"); ex = json.load(open(f"{XB}/{a}/export.json")); ve = json.load(open(f"{F}/{a}/export.json"))
    assert (xb.bin_id.to_numpy() == vb.bin_id.to_numpy()).all() and (xb.row.to_numpy() == vb.row.to_numpy()).all()
    uv = np.bincount(VX.bin, weights=VX.n, minlength=len(vb)); ux = np.bincount(XX.bin, weights=XX.n, minlength=len(xb)); reg = vb.in_reg.to_numpy()
    c0, cmax, sh, c10 = xcorr(to_map(vb, uv), to_map(xb, ux)); n0, nmax, nsh, n10 = xcorr(to_map(vb, (vb.label > 0).astype(float)), to_map(xb, (xb.label > 0).astype(float)))
    gv = np.bincount(VX.g, weights=VX.n, minlength=408); gx = np.bincount(XX.gene, weights=XX.n, minlength=408)
    s = ex["side_um"] / 1000
    rows.append({"area": a, "v_nuclei_keep_per_mm2": ve["n_nuclei_keep"] / s ** 2, "x_nuclei_per_mm2": ex["xen_nuclei_per_mm2"],
                 "ratio_nuclei_x_over_v": ex["xen_nuclei_per_mm2"] / (ve["n_nuclei_keep"] / s ** 2),
                 "v_share_nucleus_bins": float((vb.label > 0).mean()), "x_share_nucleus_bins": float((xb.label > 0).mean()),
                 "v_share_in_tissue": ex["share_in_tissue_v"], "v_share_in_reg": ex["share_in_reg_v"],
                 "x_share_bins_in_cell_all": ex["share_bins_in_xen_cell_all"], "x_share_bins_in_cell_in_tissue": ex["share_bins_in_xen_cell_in_tissue"],
                 "x_share_bins_in_cell_in_reg": ex["share_bins_in_xen_cell_in_reg"],
                 "v_mean408_reg": float(uv[reg].mean()), "x_mean408_reg": float(ux[reg].mean()), "ratio_umi_x_over_v": float(ux[reg].mean() / uv[reg].mean()),
                 "xcorr_umi_at0": float(c0), "xcorr_umi_peak": float(cmax), "xcorr_umi_peak_shift_bins_row": sh[0], "xcorr_umi_peak_shift_bins_col": sh[1],
                 "xcorr_umi_peak_shift_um": float(2 * np.hypot(*sh)), "xcorr_umi_at_20um_diag": float(c10),
                 "xcorr_nuc_at0": float(n0), "xcorr_nuc_peak": float(nmax), "xcorr_nuc_peak_shift_um": float(2 * np.hypot(*nsh)), "xcorr_nuc_at_20um_diag": float(n10),
                 "spearman_gene_totals": float(stats.spearmanr(gv, gx).correlation), "pearson_log_gene_totals": float(np.corrcoef(np.log1p(gv), np.log1p(gx))[0, 1])})
    print(rows[-1], flush=True)
df = pd.DataFrame(rows); df.to_csv("results/riesame/xen_plausibility.csv", index=False); print(df.round(3).T.to_string())
