#!/usr/bin/env python3
"""R1_verify.py — verifiche della sessione R1 (geo_spatialtrans).

Per ogni dataset del manifest:
  C-R1.1  md5 ricalcolato di ogni file scaricato = md5 pubblicato da 10x (manifest)
  C-R1.2  immagine microscopica H&E apribile (tifffile); scala µm/px letta da
          scalefactors_json.json (microns_per_pixel) e dai tag TIFF di risoluzione
  B-R1.3  µm/px <= 0.5 PASS, 0.5-1 WARN, > 1 FAIL
  CP-1    controprova: passo della griglia 2 µm in pixel full-res (da
          tissue_positions.parquet) × microns_per_pixel = 2.00 µm ± 5 %
  CP-2    ritagli H&E a risoluzione nativa (400 µm) + overview con bin in_tissue
Scrive results/R1/*.csv e results/R1/crops/*.png (README e checksums.md5 sono generati a parte).

Uso: .venv/bin/python tools/R1_verify.py --data data/real/datasets \
        --manifest data/real/inventory/R1_download_manifest.csv --out results/R1 [--skip-md5]
"""
import argparse, csv, hashlib, json, os, re, subprocess, sys, tarfile, time
from pathlib import Path

import numpy as np
import pandas as pd

def md5sum(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

def extract_tar(tar_path, dest, marker):
    """Estrae una volta sola (marker file). Ritorna True se ha estratto ora."""
    if marker.exists():
        return False
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as t:
        t.extractall(dest, filter="data")
    marker.write_text(time.strftime("%FT%TZ", time.gmtime()))
    return True

def find_one(root, pattern):
    hits = sorted(root.rglob(pattern))
    return hits[0] if hits else None

def tiff_metadata(path):
    import tifffile
    out = {"tiff_opens": False}
    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        page = s.pages[0] if hasattr(s, "pages") else tf.pages[0]
        out.update(
            tiff_opens=True, is_bigtiff=tf.is_bigtiff, is_ome=tf.is_ome,
            n_series=len(tf.series), n_levels=len(getattr(s, "levels", [s])),
            shape=str(s.shape), dtype=str(s.dtype), axes=s.axes,
            compression=str(page.compression).split(".")[-1],
            tile=str(getattr(page, "tile", None)),
        )
        # risoluzione dai tag TIFF
        xr = page.tags.get("XResolution"); ru = page.tags.get("ResolutionUnit")
        um_per_px_tag = None
        if xr is not None:
            num, den = xr.value
            if den and num:
                px_per_unit = num / den
                unit = ru.value if ru is not None else None
                unit = getattr(unit, "value", unit)
                if unit == 2:      # inch
                    um_per_px_tag = 25400.0 / px_per_unit
                elif unit == 3:    # centimetre
                    um_per_px_tag = 10000.0 / px_per_unit
                out["res_tag_raw"] = f"{num}/{den} unit={unit}"
        out["um_per_px_tifftag"] = um_per_px_tag
        # OME physical size
        if tf.is_ome and tf.ome_metadata:
            m = re.search(r'PhysicalSizeX="([\d.eE+-]+)"', tf.ome_metadata)
            u = re.search(r'PhysicalSizeXUnit="([^"]+)"', tf.ome_metadata)
            if m:
                out["ome_physical_size_x"] = f"{m.group(1)} {u.group(1) if u else ''}"
        # dimensioni YX del livello 0
        ax = s.axes; shp = s.shape
        out["height_px"] = shp[ax.index("Y")]; out["width_px"] = shp[ax.index("X")]
    return out

def read_region(path, r0, r1, c0, c1):
    """Legge una regione [r0:r1, c0:c1] del livello 0 senza caricare tutta l'immagine."""
    import tifffile, zarr
    store = tifffile.imread(path, aszarr=True)
    try:
        z = zarr.open(store, mode="r")
        if isinstance(z, zarr.Group):
            z = z["0"]
        with tifffile.TiffFile(path) as tf:
            axes = tf.series[0].axes
        if axes.startswith("YX"):
            arr = z[r0:r1, c0:c1]
        elif axes.startswith("SYX") or axes.startswith("CYX"):
            arr = np.moveaxis(z[:, r0:r1, c0:c1], 0, -1)
        else:
            raise ValueError(f"axes {axes} non gestiti")
    finally:
        store.close()
    if arr.ndim == 3 and arr.shape[-1] > 3:
        arr = arr[..., :3]
    return np.asarray(arr)

def grid_pitch_px(pos):
    """Passo della griglia in pixel full-res: least squares (pxl_col,pxl_row) ~ A·(array_col,array_row)+b."""
    p = pos.dropna(subset=["pxl_col_in_fullres", "pxl_row_in_fullres"])   # tutti i bin (rilievo revisore R1)
    X = np.c_[p.array_col.values, p.array_row.values, np.ones(len(p))]
    Y = np.c_[p.pxl_col_in_fullres.values, p.pxl_row_in_fullres.values]
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    A = coef[:2, :].T                      # 2x2
    sv = np.linalg.svd(A, compute_uv=False)
    resid = Y - X @ coef
    ax1 = float(np.hypot(*A[:, 0])); ax2 = float(np.hypot(*A[:, 1]))     # passo lungo array_col e array_row
    rot = float(np.degrees(np.arctan2(A[1, 0], A[0, 0])))                  # orientamento dell'asse array_col nell'immagine
    return dict(pitch_px_sv1=sv[0], pitch_px_sv2=sv[1], pitch_px=float(np.sqrt(abs(np.linalg.det(A)))),
                pitch_px_axis1=ax1, pitch_px_axis2=ax2, grid_rotation_deg=rot,
                rms_resid_px=float(np.sqrt((resid ** 2).mean())), n_fit=len(p))

def status(v, ok, warn=None):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "FAIL"
    if ok(v):
        return "PASS"
    if warn is not None and warn(v):
        return "WARN"
    return "FAIL"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True); ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--skip-md5", action="store_true")
    ap.add_argument("--crop-um", type=float, default=400.0)
    a = ap.parse_args()
    data = Path(a.data); out = Path(a.out); crops = out / "crops"; crops.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(a.manifest).drop_duplicates("url")
    checks, meta_rows, scale_rows, md5_rows = [], [], [], []

    # ---------- C-R1.1 md5 ----------
    for _, r in man.iterrows():
        f = data / r.dataset_id / r.filename
        if not f.exists():
            md5_rows.append(dict(dataset_id=r.dataset_id, filename=r.filename, expected=r.md5, observed="MISSING", bytes=0, status="FAIL")); continue
        obs = "skipped" if a.skip_md5 else md5sum(f)
        st = "PASS" if (a.skip_md5 or obs == r.md5) else "FAIL"
        md5_rows.append(dict(dataset_id=r.dataset_id, filename=r.filename, expected=r.md5, observed=obs, bytes=f.stat().st_size, status=st))
        print(f"[md5] {st} {r.filename}", flush=True)
    md5df = pd.DataFrame(md5_rows); md5df.to_csv(out / "R1_md5.csv", index=False)
    for ds_id, g in md5df.groupby("dataset_id"):
        checks.append(dict(check="C-R1.1", dataset_id=ds_id, value=f"{(g.status=='PASS').sum()}/{len(g)} md5 ok",
                           threshold="tutti", status="PASS" if (g.status == "PASS").all() else "FAIL"))

    # ---------- per dataset ----------
    arch_of = man.groupby("dataset_id").archetype.agg(lambda s: "/".join(sorted(set(s)))).to_dict()
    for ds_id in man.dataset_id.unique():
        d = data / ds_id; arch = arch_of[ds_id]
        print(f"\n=== {ds_id} ({arch})", flush=True)
        # estrazione tar (una volta)
        for tarname, sub in [("spatial.tar.gz", "spatial_toplevel"), ("binned_outputs.tar.gz", "binned_outputs_x"), ("segmented_outputs.tar.gz", "segmented_outputs_x")]:
            tp = find_one(d, f"*{tarname}")
            if tp:
                did = extract_tar(tp, d / sub, d / sub / ".extracted")
                print(f"[tar] {tarname}: {'estratto' if did else 'già presente'}", flush=True)
        # immagine microscopica
        img = find_one(d, "*_tissue_image.btf") or find_one(d, "*_tissue_image.tif")
        cyt = find_one(d, "*_image.tif")
        meta = dict(dataset_id=ds_id, archetype=arch, microscope_image=img.name if img else None)
        if img:
            try:
                meta.update(tiff_metadata(img))
            except Exception as e:
                meta.update(tiff_opens=False, error=str(e)[:200])
        if cyt:
            try:
                cm = tiff_metadata(cyt); meta.update(cytassist_shape=cm["shape"], cytassist_um_per_px_tag=cm["um_per_px_tifftag"])
            except Exception as e:
                meta.update(cytassist_error=str(e)[:200])
        checks.append(dict(check="C-R1.2a", dataset_id=ds_id, value=f"opens={meta.get('tiff_opens')} shape={meta.get('shape')} comp={meta.get('compression')} levels={meta.get('n_levels')}",
                           threshold="immagine apribile", status="PASS" if meta.get("tiff_opens") else "FAIL"))
        # scalefactors (2 µm) e posizioni
        sf2 = find_one(d, "square_002um/spatial/scalefactors_json.json")
        pos2 = find_one(d, "square_002um/spatial/tissue_positions.parquet")
        sf_top = find_one(d / "spatial_toplevel", "scalefactors_json.json") if (d / "spatial_toplevel").exists() else None
        sc = dict(dataset_id=ds_id, archetype=arch)
        upp = None
        for tag, p in [("sf2", sf2), ("sftop", sf_top)]:
            if p:
                j = json.load(open(p)); sc[f"{tag}_keys"] = ",".join(sorted(j.keys()))
                for k in ("microns_per_pixel", "bin_size_um", "spot_diameter_fullres", "tissue_hires_scalef", "tissue_lowres_scalef", "regist_target_img_scalef"):
                    if k in j: sc[f"{tag}_{k}"] = j[k]
                if upp is None and "microns_per_pixel" in j: upp = float(j["microns_per_pixel"])
        sc["um_per_px_scalefactors"] = upp; sc["um_per_px_tifftag"] = meta.get("um_per_px_tifftag")
        checks.append(dict(check="C-R1.2b", dataset_id=ds_id, value=f"microns_per_pixel={upp}", threshold="presente", status="PASS" if upp else "FAIL"))
        checks.append(dict(check="B-R1.3", dataset_id=ds_id, value=f"{upp} µm/px" if upp else None, threshold="<=0.5 PASS; <=1 WARN",
                           status=status(upp, lambda v: v <= 0.5, lambda v: v <= 1.0)))
        # ---------- CP-1: passo della griglia ----------
        if pos2 and upp:
            pos = pd.read_parquet(pos2)
            sc["n_bins_2um"] = len(pos); sc["n_in_tissue_2um"] = int(pos.in_tissue.sum()) if "in_tissue" in pos else None
            gp = grid_pitch_px(pos); sc.update(gp)
            pitch_um = gp["pitch_px"] * upp; sc["pitch_um_from_grid"] = pitch_um
            sc["pitch_um_axis1"] = gp["pitch_px_axis1"] * upp; sc["pitch_um_axis2"] = gp["pitch_px_axis2"] * upp
            aniso = gp["pitch_px_sv1"] / gp["pitch_px_sv2"]; sc["grid_anisotropy_sv1_sv2"] = aniso
            checks.append(dict(check="CP-1", dataset_id=ds_id, value=f"pitch={gp['pitch_px']:.4f} px × {upp} = {pitch_um:.4f} µm (sv ratio {aniso:.4f}, rms {gp['rms_resid_px']:.3f} px)",
                               threshold="2.00 µm ± 5 %", status=status(pitch_um, lambda v: abs(v - 2.0) <= 0.10)))
            # posizioni dentro l'immagine?
            if meta.get("tiff_opens"):
                inside = ((pos.pxl_col_in_fullres >= 0) & (pos.pxl_col_in_fullres < meta["width_px"]) & (pos.pxl_row_in_fullres >= 0) & (pos.pxl_row_in_fullres < meta["height_px"]))
                frac = float(inside.mean()); sc["frac_all_bins_inside_image"] = frac
                it_mask = (pos.in_tissue == 1) if "in_tissue" in pos else np.ones(len(pos), bool)
                frac_it = float(inside[it_mask].mean()); sc["frac_in_tissue_bins_inside_image"] = frac_it
                sc["capture_area_px_est"] = f"{(pos.array_col.max()-pos.array_col.min()+1)*gp['pitch_px']:.0f}x{(pos.array_row.max()-pos.array_row.min()+1)*gp['pitch_px']:.0f}"
                checks.append(dict(check="C-R1.2c", dataset_id=ds_id, value=f"in_tissue: {frac_it:.4f} dentro l'immagine {meta['width_px']}x{meta['height_px']} (tutti i bin: {frac:.4f}; area di cattura stimata {sc['capture_area_px_est']} px)",
                                   threshold=">= 0.99 dei bin in_tissue", status=status(frac_it, lambda v: v >= 0.99)))
            # ---------- CP-2: ritagli nativi ----------
            if meta.get("tiff_opens"):
                try:
                    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
                    it = pos[pos.in_tissue == 1] if "in_tissue" in pos else pos
                    half = int(round(a.crop_um / upp / 2))
                    rng = np.random.default_rng(42)
                    # 3 ritagli: centroide del tessuto + 2 punti casuali in tessuto
                    centers = [(int(it.pxl_row_in_fullres.median()), int(it.pxl_col_in_fullres.median()))]
                    for _ in range(2):
                        s = it.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
                        centers.append((int(s.pxl_row_in_fullres), int(s.pxl_col_in_fullres)))
                    for k, (cr, cc) in enumerate(centers):
                        r0, c0 = max(cr - half, 0), max(cc - half, 0)
                        reg = read_region(img, r0, r0 + 2 * half, c0, c0 + 2 * half)
                        fig, ax = plt.subplots(figsize=(6, 6), dpi=150, facecolor="white")
                        ax.imshow(reg, interpolation="nearest"); ax.set_axis_off()
                        bar_px = 50 / upp
                        ax.plot([20, 20 + bar_px], [reg.shape[0] - 25] * 2, "k-", lw=3); ax.text(20, reg.shape[0] - 35, "50 µm", fontsize=8)
                        ax.set_title(f"{arch} · {ds_id}\ncrop {k} @ (row {cr}, col {cc}), {a.crop_um:.0f} µm, {upp:.4f} µm/px", fontsize=8)
                        fig.savefig(crops / f"{arch.replace('/','-')}_{ds_id}_crop{k}.png", bbox_inches="tight"); plt.close(fig)
                    # overview: hires image + bin in_tissue (sottocampionati)
                    hires = find_one(d / "spatial_toplevel", "tissue_hires_image.png") or find_one(d, "square_008um/spatial/tissue_hires_image.png")
                    sfh = json.load(open(sf_top)) if sf_top else (json.load(open(sf2)) if sf2 else {})
                    if hires and "tissue_hires_scalef" in sfh:
                        from PIL import Image
                        him = np.asarray(Image.open(hires)); s = sfh["tissue_hires_scalef"]
                        sub = it.sample(min(len(it), 100_000), random_state=42)
                        fig, ax = plt.subplots(figsize=(7, 7), dpi=130, facecolor="white")
                        ax.imshow(him); ax.scatter(sub.pxl_col_in_fullres * s, sub.pxl_row_in_fullres * s, s=0.05, c="lime", alpha=0.3)
                        for k, (cr, cc) in enumerate(centers):
                            ax.add_patch(plt.Rectangle(((cc - half) * s, (cr - half) * s), 2 * half * s, 2 * half * s, fill=False, ec="red", lw=1)); ax.text((cc + half) * s, (cr - half) * s, str(k), color="red", fontsize=7)
                        ax.set_axis_off(); ax.set_title(f"{arch} · {ds_id}: bin 2 µm in_tissue (verde) su tissue_hires_image; riquadri = ritagli", fontsize=8)
                        fig.savefig(crops / f"{arch.replace('/','-')}_{ds_id}_overview.png", bbox_inches="tight"); plt.close(fig)
                    checks.append(dict(check="CP-2", dataset_id=ds_id, value=f"{len(centers)} ritagli + overview", threshold="verifica visiva di Luca", status="PENDING"))
                except Exception as e:
                    checks.append(dict(check="CP-2", dataset_id=ds_id, value=f"errore: {str(e)[:150]}", threshold="—", status="FAIL"))
        else:
            checks.append(dict(check="CP-1", dataset_id=ds_id, value="tissue_positions.parquet o microns_per_pixel mancanti", threshold="—", status="FAIL"))
        # Space Ranger version dal web_summary
        ws = find_one(d, "*web_summary.html")
        if ws:
            t = ws.read_text(errors="ignore")
            m = re.search(r'(?:spaceranger|pipeline_version|Pipeline Version)[^0-9]{0,40}([0-9]+\.[0-9]+\.[0-9]+)', t, flags=re.I)
            if m:
                meta["spaceranger_version_websummary"] = m.group(1); meta["spaceranger_version_source"] = "etichetta"
            else:
                from collections import Counter
                cands = Counter(re.findall(r'\b([34]\.[0-9]\.[0-9])\b', t))
                meta["spaceranger_version_websummary"] = cands.most_common(1)[0][0] if cands else None
                meta["spaceranger_version_source"] = f"stringa piu frequente {dict(cands)}" if cands else None
        meta_rows.append(meta); scale_rows.append(sc)
        pd.DataFrame(checks).to_csv(out / "R1_checks.csv", index=False)  # salvataggio incrementale

    pd.DataFrame(meta_rows).to_csv(out / "R1_image_metadata.csv", index=False)
    pd.DataFrame(scale_rows).to_csv(out / "R1_scale.csv", index=False)
    ck = pd.DataFrame(checks); ck.to_csv(out / "R1_checks.csv", index=False)
    print("\n" + ck.to_string(index=False, max_colwidth=90))
    print("\nRIEPILOGO:", ck.status.value_counts().to_dict())
    return 0 if not (ck.status == "FAIL").any() else 1

if __name__ == "__main__":
    sys.exit(main())
