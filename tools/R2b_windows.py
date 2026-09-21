"""tools/R2b_windows.py — R2b: campionamento con seed delle finestre per la verita' manuale (BL-024, BL-025).

Per archetipo: 4 finestre quadrate a lato adattato alla densita' R2 (~70-110 nuclei attesi), una per ROI su 4 ROI distinti
(scelti a caso fra i 5 di R2), posizione uniforme dentro il ROI con margine, accettata se tessuto valido >= 90 % e 0 % bolle.
Check C-R2b.1: finestra dentro il ROI, ritaglio dal ROI identico al ritaglio letto dall'immagine intera alle coordinate
globali ricostruite. Check C-R2b.4: md5 dell'array della finestra (riproducibilita' con lo stesso seed).
Output: results/R2b/R2b_windows.csv, results/R2b/R2b_windows_check.csv, /mnt/micron/geo_spatialtrans/R2b/windows/<win>.png (fuori git)
Uso: .venv/bin/python tools/R2b_windows.py
"""
import sys, hashlib, numpy as np, pandas as pd
from pathlib import Path
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, rois_table, load_roi, masks_for, NUC_DIR
from R2_rois import img_path, read_region

OUT = ROOT / "results/R2b"; OUT.mkdir(exist_ok=True)
WIN_DIR = Path("/mnt/micron/geo_spatialtrans/R2b/windows"); WIN_DIR.mkdir(parents=True, exist_ok=True)
SEED = 20260921
SIDE_UM = {"A1": 100, "A2": 120, "A3": 150, "A4": 60, "A5": 250, "A6": 250}   # pre-registrati (chat 2026-09-21)
N_WIN = 4; TISSUE_MIN = 0.90; BUBBLE_MAX = 0.0; MARGIN_UM = 20.0; MAX_TRIES = 5000
METHODS = ["cellpose_rgb", "stardist_he", "spaceranger"]

def md5_arr(a):
    return hashlib.md5(np.ascontiguousarray(a).tobytes()).hexdigest()

def count_in_window(archetype, roi_id, method, x0, y0, side_um):
    f = NUC_DIR / f"{archetype}_{roi_id}_{method}_native.parquet"
    if not f.exists(): return np.nan
    d = pd.read_parquet(f, columns=["x_um", "y_um", "keep"])
    d = d[d.keep]
    return int(((d.x_um >= x0) & (d.x_um < x0 + side_um) & (d.y_um >= y0) & (d.y_um < y0 + side_um)).sum())

def main():
    rois = rois_table()
    cons = pd.read_csv(ROOT / "results/R2/R2_consensus_density.csv").set_index(["archetype", "roi_id"])
    rows = []
    for ai, (arch, g) in enumerate(rois.groupby("archetype", sort=True)):
        rng = np.random.default_rng(SEED + ai)
        chosen = sorted(rng.choice(g.roi_id.values, N_WIN, replace=False).tolist())
        rater2 = set(rng.choice(chosen, 2, replace=False).tolist())
        for wi, roi_id in enumerate(chosen, 1):
            r = g[g.roi_id == roi_id].iloc[0]
            upp = float(r.um_per_px); side_um = SIDE_UM[arch]; side = int(round(side_um / upp)); margin = int(round(MARGIN_UM / upp))
            rgb = load_roi(arch, roi_id); H, W = rgb.shape[:2]
            t, b = masks_for(rgb, upp, 1, arch, roi_id)
            accepted = None
            for k in range(1, MAX_TRIES + 1):
                c0 = int(rng.integers(margin, W - side - margin)); r0 = int(rng.integers(margin, H - side - margin))
                tf = float(t[r0:r0 + side, c0:c0 + side].mean()); bf = float(b[r0:r0 + side, c0:c0 + side].mean())
                if tf >= TISSUE_MIN and bf <= BUBBLE_MAX:
                    accepted = (c0, r0, tf, bf, k); break
            assert accepted is not None, f"{arch} {roi_id}: nessuna finestra valida in {MAX_TRIES} tentativi"
            c0, r0, tf, bf, k = accepted
            win = rgb[r0:r0 + side, c0:c0 + side]
            # C-R2b.1: round trip dall'immagine intera alle coordinate globali ricostruite
            C0, R0 = int(r.c0) + c0, int(r.r0) + r0
            full = read_region(img_path(r.dataset), R0, R0 + side, C0, C0 + side)
            identical = bool(full.shape == win.shape and np.array_equal(full, win))
            maxdiff = int(np.abs(full.astype(int) - win.astype(int)).max()) if full.shape == win.shape else -1
            win_id = f"{arch}_w{wi}"
            png = WIN_DIR / f"{win_id}.png"; Image.fromarray(win).save(png, compress_level=6)
            x0_um, y0_um = c0 * upp, r0 * upp
            row = dict(win_id=win_id, archetype=arch, dataset=r.dataset, roi_id=roi_id, roi_label=r.label, side_um=side_um, side_px=side,
                       um_per_px=upp, area_mm2=(side * upp) ** 2 / 1e6, c0_roi=c0, r0_roi=r0, x0_um=x0_um, y0_um=y0_um, c0_img=C0, r0_img=R0,
                       inside_roi=bool(c0 >= 0 and r0 >= 0 and c0 + side <= W and r0 + side <= H),
                       tissue_frac=tf, bubble_frac=bf, n_tries=k, roundtrip_identical=identical, roundtrip_maxdiff=maxdiff,
                       md5_array=md5_arr(win), png=str(png), rater2=roi_id in rater2,
                       density_consensus_roi=float(cons.loc[(arch, roi_id), "density_consensus"]))
            row["expected_consensus"] = row["density_consensus_roi"] * row["area_mm2"]
            for m in METHODS:
                row[f"n_{m}"] = count_in_window(arch, roi_id, m, x0_um, y0_um, side * upp)
            rows.append(row)
            print(f"{win_id} {roi_id} side={side}px tissue={tf:.3f} bubble={bf:.4f} tries={k} roundtrip={identical} "
                  f"nCP={row['n_cellpose_rgb']} nSD={row['n_stardist_he']} nSR={row['n_spaceranger']} exp={row['expected_consensus']:.0f}")
    df = pd.DataFrame(rows); df.to_csv(OUT / "R2b_windows.csv", index=False)
    chk = pd.DataFrame({"win_id": df.win_id,
                        "inside_roi": df.inside_roi, "tissue_ge_0.90": df.tissue_frac >= TISSUE_MIN, "bubble_eq_0": df.bubble_frac <= BUBBLE_MAX,
                        "roundtrip_identical": df.roundtrip_identical})
    chk["C-R2b.1"] = np.where(chk.iloc[:, 1:].all(axis=1), "PASS", "FAIL")
    chk.to_csv(OUT / "R2b_windows_check.csv", index=False)
    print(chk["C-R2b.1"].value_counts().to_dict())

if __name__ == "__main__":
    main()

