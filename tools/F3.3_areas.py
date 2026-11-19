"""tools/F3.3_areas.py - F3.3: tabella delle aree (30 di taratura di R2 e 38 di verifica di F3.2, tessuti nuovi esclusi, DEC: Luca
2026-11-19 «Lasciarli sigillati»). Riquadro in pixel dell'immagine a piena risoluzione: c0 = round(x0 / upp), lato = round(lato_um / upp)
(per le aree di R2 si ritrovano i riquadri di results/R2/R2_rois_checked.csv). Scarto di registrazione: quello della sezione usato in
F3.1 per le sezioni di taratura (export.json), 0 per le altre (DEC-097, punto 3). Uso: python3 tools/F3.3_areas.py"""
import json, glob, numpy as np, pandas as pd
from pathlib import Path
M = Path("/mnt/micron/geo_spatialtrans"); DR, DH = M / "data_real", M / "data_holdout"
SEC = {"A1": ("Visium_HD_Mouse_Small_Intestine", DR), "A23": ("Visium_HD_Human_Colon_Cancer", DR), "A4": ("Visium_HD_Human_Lymph_Node_FFPE", DR),
       "A5": ("Visium_HD_6p5mm_Mouse_Brain", DR), "A6": ("Visium_HD_6p5mm_Human_Heart", DR),
       "CRC_P1_SR300": ("CRC_P1_SR300", DH), "CRC_P5_SR300": ("CRC_P5_SR300", DH), "CRC_6p5mm_SR410": ("CRC_6p5mm_SR410", DH),
       "MouseBrain_FFPE_SR401": ("MouseBrain_FFPE_SR401", DH)}
ARCH_OF = {"mucosa": "A1", "tumore": "A2", "stroma": "A3", "follicolo": "A4", "corteccia": "A5", "miocardio": "A6"}
SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
CRC = {"A23", "CRC_P1_SR300", "CRC_P5_SR300", "CRC_6p5mm_SR410"}
def sec_paths(label):
    ds, root = SEC[label]; d = root / ds
    img = [p for p in d.glob("*_tissue_image.*") if p.suffix in (".btf", ".tif", ".tiff")]; assert len(img) == 1, (label, img)
    sq = d / "binned_outputs_x/binned_outputs"; upp = json.load(open(sq / "square_008um/spatial/scalefactors_json.json"))["microns_per_pixel"]
    return str(d), str(img[0]), str(sq), float(upp)
dimg = {}
for f in glob.glob(str(M / "F3.1/main/real/*/export.json")):
    j = json.load(open(f)); lab = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}[j["archetype"]]
    v = tuple(j["d_img_um"]); assert dimg.setdefault(lab, v) == v, "scarto non costante nella sezione"
rows = []
R2 = pd.read_csv("results/R2/R2_rois_checked.csv")
for _, r in R2.iterrows():
    lab = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}[r.archetype]; d, img, sq, upp = sec_paths(lab)
    assert abs(upp - r.um_per_px) < 1e-12
    rows.append(dict(area=f"{r.archetype}_{r.roi_id}", A=r.archetype, klass=None, level="IN", section=lab, set=f"{lab}_IN", ds_dir=d, image=img, sq=sq, upp=upp,
                     c0=int(r.c0), r0=int(r.r0), side_px=int(r.side_px), side_um=float(r.side_px * upp), seg=SEG[r.archetype], crc=lab in CRC,
                     d_img_r=dimg[lab][0], d_img_c=dimg[lab][1], roi_id=r.roi_id))
S = pd.read_csv("results/F3.2/F3.2_selected_all.csv"); S = S[S.level.isin(["L1", "L2"])]
for _, r in S.iterrows():
    A = ARCH_OF[r.klass]; d, img, sq, upp = sec_paths(r.label); side_px = int(round(r.side_um / upp))
    c0, r0 = int(round(r.x0 / upp)), int(round(r.y0 / upp)); di = dimg.get(r.label, (0.0, 0.0)) if r.level == "L1" else (0.0, 0.0)
    rows.append(dict(area=r.area_id, A=A, klass=r.klass, level=r.level, section=r.label, set=f"{r.label}_{r.level}", ds_dir=d, image=img, sq=sq, upp=upp,
                     c0=c0, r0=r0, side_px=side_px, side_um=float(side_px * upp), seg=SEG[A], crc=r.label in CRC, d_img_r=di[0], d_img_c=di[1], roi_id=None))
T = pd.DataFrame(rows); T["r1"] = T.r0 + T.side_px; T["c1"] = T.c0 + T.side_px
T.to_csv("results/F3.3/F3.3_areas.csv", index=False)
print(T.groupby(["level", "section", "A"]).size().to_string()); print(len(T))
