#!/usr/bin/env python3
"""tools/F3.2_bridge.py - F3.2, previsioni P-1 (controprova CP-2): Space Ranger 3.0.0 contro 4.0.1 sulla stessa sezione P2 del tumore del colon.
Per ogni area di R2 del tumore (A2) e dello stroma (A3): quadratini di 2 um dentro il riquadro dell'area e in tessuto per la 4.0.1 (stessa scelta di
tools/R4_roi.py), gli stessi codici a barre nella 3.0.0 (matrice grezza). Scarto dichiarato: la grezza della 3.0.0 non contiene
tutti i quadratini (9 468 111 contro 11 222 500); il confronto si fa sui quadratini presenti in entrambe, la quota mancante si riporta a parte. Grandezze: UMI medi per quadratino, quota di
quadratini vuoti, geni rilevati per quadratino. Controllo interno: il lato 4.0.1 deve riprodurre results/R4/R4_roi.csv (tolleranza relativa 1e-9).
Stima per tessuto: media sulle 5 aree della differenza relativa (3.0.0 - 4.0.1) / 4.0.1; errore standard fra aree * kappa di R7; esito con +- 2 ES.
Uso: .venv/bin/python tools/F3.2_bridge.py -> results/F3.2/F3.2_bridge_roi.csv, F3.2_bridge_predictions.csv, F3.2_bridge_image_md5.csv
"""
import sys, json, csv
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/F3.2"
SR4 = Path("/mnt/micron/geo_spatialtrans/data_real/Visium_HD_Human_Colon_Cancer/binned_outputs_x/binned_outputs/square_002um")
SR3 = Path("/mnt/micron/geo_spatialtrans/data_bridge/CRC_P2_SR300/binned_outputs_x/binned_outputs/square_002um")
# previsioni pre-registrate (file dello step, commit 1a7c698): margine = cv_between_areas / 4, kappa di R7
PRED = [("P-1a", "A2", "mean_umi2", "B-280", "mean UMI per 2 um bin in tissue"), ("P-1b", "A2", "frac_empty2", "B-280", "fraction of empty 2 um bins in tissue"),
        ("P-1c", "A2", "mean_genes2", "B-282", "genes detected per 2 um bin"), ("P-1d", "A3", "mean_umi2", "B-290", "mean UMI per 2 um bin in tissue"),
        ("P-1e", "A3", "frac_empty2", "B-290", "fraction of empty 2 um bins in tissue"), ("P-1f", "A3", "mean_genes2", "B-292", "genes detected per 2 um bin")]

class H5:
    def __init__(self, p):
        f = h5py.File(p, "r")["matrix"]; self.f = f
        self.bc = pd.Series(np.arange(f["shape"][1]), index=[b.decode() for b in f["barcodes"][:]])
        self.indptr = f["indptr"][:]; self.n_genes = int(f["shape"][0])
    def has(self, barcodes):
        return np.asarray(pd.Index(barcodes).isin(self.bc.index))
    def totals(self, barcodes):
        """UMI e geni rilevati per codice a barre (letture per blocchi contigui di colonne); i codici a barre devono esistere"""
        cols = self.bc.reindex(barcodes).to_numpy(); miss = 0
        assert not pd.isna(cols).any(), "codici a barre assenti: filtrare prima con has()"
        cols = cols.astype(np.int64)
        order = np.argsort(cols); umi = np.zeros(len(cols)); ng = np.zeros(len(cols))
        data = self.f["data"]
        lo, hi = self.indptr[cols[order].min()], self.indptr[cols[order].max() + 1]
        d = data[lo:hi].astype(np.float64)
        for k in order:
            a, b = self.indptr[cols[k]] - lo, self.indptr[cols[k] + 1] - lo
            umi[k] = d[a:b].sum(); ng[k] = np.count_nonzero(d[a:b])
        return umi, ng, miss

def main():
    rois = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype.isin(["A2", "A3"])]
    r4 = pd.read_csv(ROOT / "results/R4/R4_roi.csv")
    tp = pq.read_table(SR4 / "spatial/tissue_positions.parquet").to_pandas()
    m4 = H5(SR4 / "filtered_feature_bc_matrix.h5"); m3 = H5(SR3 / "raw_feature_bc_matrix.h5")
    rows = []
    for _, R in rois.iterrows():
        box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)]
        bcs = box.barcode[box.barcode.isin(m4.bc.index)].to_numpy()
        u4full, _, _ = m4.totals(bcs)                                     # insieme completo: controllo contro R4
        both = m3.has(bcs); miss3 = int((~both).sum())                   # scarto dichiarato: confronto sui quadratini presenti in entrambe
        u4, g4, _ = m4.totals(bcs[both]); u3, g3, _ = m3.totals(bcs[both])
        ref = r4[(r4.archetype == R.archetype) & (r4.roi_id == R.roi_id)].iloc[0]
        rec = dict(archetype=R.archetype, roi_id=R.roi_id, n_bins=len(bcs), n_bins_both=int(both.sum()), n_missing_in_sr3=miss3,
                   frac_missing_in_sr3=miss3 / len(bcs), sr4full_mean_umi2=u4full.mean(),
                   sr4_mean_umi2=u4.mean(), sr3_mean_umi2=u3.mean(), sr4_frac_empty2=(u4 == 0).mean(), sr3_frac_empty2=(u3 == 0).mean(),
                   sr4_mean_genes2=g4.mean(), sr3_mean_genes2=g3.mean(), r4_mean_umi2=ref.mean_umi2, r4_frac_empty2=ref.frac_empty2,
                   r4_n_bins2=ref.n_bins2, frac_bins_identical_umi=float((u3 == u4).mean()))
        for q in ("mean_umi2", "frac_empty2", "mean_genes2"):
            rec[f"rel_{q}"] = (rec[f"sr3_{q}"] - rec[f"sr4_{q}"]) / rec[f"sr4_{q}"]
        rec["r4_reproduced"] = bool(len(bcs) == ref.n_bins2 and abs(u4full.mean() / ref.mean_umi2 - 1) < 1e-9
                                    and abs((u4full == 0).mean() - ref.frac_empty2) <= 1e-9 * max(ref.frac_empty2, 1e-12) + 1e-15)
        rows.append(rec); print({k: rec[k] for k in ("archetype", "roi_id", "n_bins", "rel_mean_umi2", "rel_frac_empty2", "rel_mean_genes2", "frac_missing_in_sr3", "r4_reproduced")}, flush=True)
    d = pd.DataFrame(rows); d.to_csv(OUT / "F3.2_bridge_roi.csv", index=False)
    r7 = pd.read_csv(ROOT / "results/R7/Mc/R7_refs_adopted.csv")
    pr = []
    for pid, A, q, bid, metric in PRED:
        x = d[d.archetype == A][f"rel_{q}"].to_numpy()
        rr = r7[(r7["id"] == bid) & (r7.metric == metric)].iloc[0]
        margin = float(rr.cv_between_areas) / 4; kappa = float(rr.kappa)
        est = x.mean(); se = x.std(ddof=1) / np.sqrt(len(x)); lo, hi = est - 2 * se * kappa, est + 2 * se * kappa
        verdict = "confermata" if (lo >= -margin and hi <= margin) else ("smentita" if (lo > margin or hi < -margin) else "indeterminata")
        pr.append(dict(prediction=pid, archetype=A, quantity=q, n_areas=len(x), estimate=est, se=se, kappa=kappa, lo=lo, hi=hi,
                       margin=margin, margin_source=f"{bid} cv_between_areas/4 = {rr.cv_between_areas:.4f}/4", verdict=verdict))
        print(pr[-1], flush=True)
    pd.DataFrame(pr).to_csv(OUT / "F3.2_bridge_predictions.csv", index=False)
    inv = pd.read_csv(ROOT / "data/real/inventory/files_10x_pages_2026-09-18.csv")
    im = inv[(inv.dataset_id == "Visium_HD_Human_Colon_Cancer") & inv.filename.str.contains("tissue_image")]
    im = im.assign(version=im.url_download.str.extract(r"/spatial-exp/([0-9.]+)/")[0])[["version", "filename", "md5", "size_as_shown"]]
    im["same_md5_as_4.0.1"] = im.md5 == im[im.version == "4.0.1"].md5.iloc[0]
    im.to_csv(OUT / "F3.2_bridge_image_md5.csv", index=False); print(im.to_string())

if __name__ == "__main__":
    main()
