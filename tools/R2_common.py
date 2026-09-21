"""tools/R2_common.py — utilità comuni della sessione R2 (importabile da .venv e .venv-stardist).
Maschera tessuto, maschera bolle di montaggio (BL-017), canale ematossilina, sottocampionamento, tabella nuclei.
"""
import json, numpy as np, pandas as pd, tifffile
from pathlib import Path
from scipy import ndimage as ndi
from skimage import color, morphology, measure, filters

ROOT = Path(__file__).resolve().parents[1]
ROI_DIR = Path("/mnt/micron/geo_spatialtrans/R2/rois")
MASK_DIR = Path("/mnt/micron/geo_spatialtrans/R2/masks"); MASK_DIR.mkdir(parents=True, exist_ok=True)
NUC_DIR = ROOT / "results/R2/nuclei"; NUC_DIR.mkdir(parents=True, exist_ok=True)
SCALES = {1: "native", 2: "ds2", 4: "ds4"}
AREA_MIN_UM2, AREA_MAX_UM2 = 4.0, 400.0     # C-R2.7 (convenzione pre-registrata)

def rois_table():
    return pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")

def load_roi(archetype, roi_id):
    return tifffile.imread(ROI_DIR / f"{archetype}_{roi_id}.tif")

def downsample(rgb, f):
    if f == 1: return rgb
    h, w = (rgb.shape[0] // f) * f, (rgb.shape[1] // f) * f
    return rgb[:h, :w].reshape(h // f, f, w // f, f, rgb.shape[2]).mean(axis=(1, 3)).astype(np.uint8)

def tissue_mask(rgb, upp):
    """Tessuto = pixel saturi (non bianco/grigio di fondo), lisciato ~2 um, chiusura 5 um, oggetti >= 500 um2, buchi < 100 um2 riempiti."""
    hsv = color.rgb2hsv(rgb)
    sat = filters.gaussian(hsv[..., 1], sigma=max(1, 2.0 / upp))
    val = filters.gaussian(hsv[..., 2], sigma=max(1, 2.0 / upp))
    m = (sat > 0.06) & (val < 0.97)
    r = max(1, int(round(5.0 / upp)))
    m = morphology.closing(m, morphology.disk(r))
    px_um2 = upp ** 2
    m = morphology.remove_small_objects(m, max_size=int(500 / px_um2))
    m = morphology.remove_small_holes(m, max_size=int(100 / px_um2))
    return m

def bubble_mask(rgb, upp):
    """Bolle di montaggio (BL-017): anelli grigi (bassa cromaticita', piu' scuri del fondo) -> dilatazione 5 um -> riempimento
    dell'interno (buchi fino a ~250 um di diametro) -> margine di sicurezza 3 um."""
    rgb = rgb.astype(np.int16)
    chroma = rgb.max(axis=2) - rgb.min(axis=2)
    mean = rgb.mean(axis=2)
    gray = (chroma < 18) & (mean < 225) & (mean > 60)
    gray = morphology.remove_small_objects(gray, max_size=int(20 / upp ** 2))
    m = morphology.dilation(gray, morphology.disk(max(1, int(round(5.0 / upp)))))
    m = morphology.remove_small_holes(m, max_size=int(np.pi * 125 ** 2 / upp ** 2))
    m = morphology.dilation(m, morphology.disk(max(1, int(round(3.0 / upp)))))
    return m

def hematoxylin(rgb):
    """Canale ematossilina (deconvoluzione colore Ruifrok), normalizzato ai percentili 1-99.8, uint8 (nuclei chiari)."""
    hed = color.rgb2hed(rgb)
    h = hed[..., 0]
    lo, hi = np.percentile(h, [1, 99.8])
    h = np.clip((h - lo) / max(hi - lo, 1e-6), 0, 1)
    return (h * 255).astype(np.uint8)

def nuclei_table(labels, upp, tmask, bmask, archetype, roi_id, method, scale):
    """Una riga per nucleo: geometria in um, flag di qualita'. Centroide in coordinate native del ROI (um dal bordo)."""
    if labels.max() == 0:
        return pd.DataFrame()
    props = measure.regionprops_table(labels, properties=("label", "area", "centroid", "eccentricity",
                                      "major_axis_length", "minor_axis_length", "solidity", "perimeter", "orientation"))
    d = pd.DataFrame(props)
    d["area_um2"] = d.area * upp ** 2
    d["eq_diam_um"] = 2 * np.sqrt(d.area_um2 / np.pi)
    d["major_um"] = d.major_axis_length * upp; d["minor_um"] = d.minor_axis_length * upp
    d["perimeter_um"] = d.perimeter * upp
    d["x_um"] = d["centroid-1"] * upp; d["y_um"] = d["centroid-0"] * upp
    ri = np.clip(d["centroid-0"].round().astype(int), 0, tmask.shape[0] - 1)
    ci = np.clip(d["centroid-1"].round().astype(int), 0, tmask.shape[1] - 1)
    d["in_tissue"] = tmask[ri, ci]; d["in_bubble"] = bmask[ri, ci]
    d["flag_small"] = d.area_um2 < AREA_MIN_UM2; d["flag_large"] = d.area_um2 > AREA_MAX_UM2
    d["keep"] = d.in_tissue & ~d.in_bubble & ~d.flag_small & ~d.flag_large
    d["archetype"] = archetype; d["roi_id"] = roi_id; d["method"] = method; d["scale"] = scale; d["um_per_px"] = upp
    return d.drop(columns=["area", "major_axis_length", "minor_axis_length", "perimeter", "centroid-0", "centroid-1"])

def save_labels(labels, archetype, roi_id, method, scale):
    np.savez_compressed(MASK_DIR / f"{archetype}_{roi_id}_{method}_{SCALES[scale]}.npz", labels=labels.astype(np.uint32))

def load_labels(archetype, roi_id, method, scale):
    return np.load(MASK_DIR / f"{archetype}_{roi_id}_{method}_{SCALES[scale]}.npz")["labels"]

def masks_for(rgb_native, upp_native, scale, archetype=None, roi_id=None):
    """Maschere tessuto/bolle: calcolate a 1/4 della risoluzione nativa (strutture >> 1 um; 16x piu' veloce), riportate alla
    risoluzione nativa per ripetizione, cache su disco per ROI (MASK_DIR/<A>_<roi>_masks.npz), poi ridotte alla scala richiesta."""
    cache = MASK_DIR / f"{archetype}_{roi_id}_masks.npz" if archetype else None
    if cache is not None and cache.exists():
        z = np.load(cache); t4, b4 = z["tissue_ds4"], z["bubble_ds4"]
    else:
        rgb4 = downsample(rgb_native, 4); upp4 = upp_native * 4
        t4 = tissue_mask(rgb4, upp4); b4 = bubble_mask(rgb4, upp4)
        if cache is not None: np.savez_compressed(cache, tissue_ds4=t4, bubble_ds4=b4)
    H, W = rgb_native.shape[0], rgb_native.shape[1]
    def up(m, f):
        mm = np.repeat(np.repeat(m, f, axis=0), f, axis=1)
        out = np.zeros((H, W), bool); h, w = min(H, mm.shape[0]), min(W, mm.shape[1]); out[:h, :w] = mm[:h, :w]
        if h < H: out[h:, :w] = out[h - 1:h, :w]
        if w < W: out[:, w:] = out[:, w - 1:w]
        return out
    def ds(m, f):
        h, w = (m.shape[0] // f) * f, (m.shape[1] // f) * f
        return m[:h, :w].reshape(h // f, f, w // f, f).mean(axis=(1, 3)) >= 0.5
    if scale == 4:
        return t4[:H // 4, :W // 4], b4[:H // 4, :W // 4]
    t, b = up(t4, 4), up(b4, 4)
    if scale == 1: return t, b
    return ds(t, scale), ds(b, scale)
