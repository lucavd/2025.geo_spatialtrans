#!/usr/bin/env python
"""tools/R2c_check_annotator.py — R2c, check C-2 sull'annotatore generato (stampa PASS/FAIL, scrive results/R2c/R2c_c2.csv).
C-2a ritagli nello zip pixel-identici al ROI R2 (con bianco fuori ROI); C-2b sulla parte che cade nella finestra R2b,
identici al PNG della finestra; C-2c marcatori: (mk + c0 + 0.5) * upp = coordinate dei punti di Luca R2b (|err| <= 1e-9 µm);
C-2d marcatori dentro il ritaglio; C-2e cecita': l'HTML non contiene parole che rivelano tipo di item, finestra o segmentatori;
C-2f ogni item del CSV ha un'immagine nello zip e le ripetute puntano alla stessa immagine dell'originale.
Uso: .venv/bin/python tools/R2c_check_annotator.py Luca|Federica
"""
import sys, io, zipfile, re
import numpy as np, pandas as pd
from pathlib import Path
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_roi
RES = ROOT / "results/R2c"; R2B = ROOT / "results/R2b"; OUT = Path("/mnt/micron/geo_spatialtrans/R2c")
BANNED = ["cellpose", "stardist", "spaceranger", "space ranger", "segment", "catch", "trappol", "repeat", "ripetut", "band",
          "fascia", "focal", "pair_id", "kind", "archetyp", "win_id", "d_um", "dmax", "d_max"]


def main(rater):
    res = []
    def ck(name, ok, val=""):
        res.append(dict(rater=rater, check=name, result="PASS" if ok else "FAIL", value=str(val)))
        print(("PASS " if ok else "FAIL ") + name + " " + str(val))
    cr = pd.read_csv(RES / f"R2c_crops_{rater}.csv"); it = pd.read_csv(RES / f"R2c_items_{rater}.csv")
    wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"].set_index(["win_id", "idx"])
    z = zipfile.ZipFile(OUT / f"R2c_annotator_{rater}.zip")
    names = set(z.namelist())
    ck("C-2f item->immagine", all(f"R2c_{rater}/{n}" in names for n in cr.img) and len(cr) == len(it), f"{len(cr)} item")
    rep = cr[cr.kind == "repeat"]
    orig = cr[cr.kind == "pair"].set_index(["win_id", "i", "j"]).img
    ck("C-2f ripetute = stessa immagine", all(orig.loc[(r.win_id, r.i, r.j)] == r.img for r in rep.itertuples()), f"{len(rep)} ripetute")
    rgb = {}; bad_a = 0; bad_b = 0; nb = 0; err = 0.0; seen = set()
    for r in cr.itertuples():
        w = wins.loc[r.win_id]; upp = float(w.um_per_px)
        pa = pts.loc[(r.win_id, r.i)]; pb = pts.loc[(r.win_id, r.j)]
        err = max(err, abs((r.mk_ax + r.c0 + 0.5) * upp - pa.x_um), abs((r.mk_ay + r.r0 + 0.5) * upp - pa.y_um),
                  abs((r.mk_bx + r.c0 + 0.5) * upp - pb.x_um), abs((r.mk_by + r.r0 + 0.5) * upp - pb.y_um))
        if r.img in seen:
            continue
        seen.add(r.img)
        k = (w.archetype, w.roi_id)
        if k not in rgb:
            rgb[k] = load_roi(*k)
        R = rgb[k]; H, W = R.shape[:2]; s = int(r.side_px)
        a = np.asarray(Image.open(io.BytesIO(z.read(f"R2c_{rater}/{r.img}"))).convert("RGB"))
        ref = np.full((s, s, 3), 255, np.uint8)
        a0, b0, a1, b1 = max(r.r0, 0), max(r.c0, 0), min(r.r0 + s, H), min(r.c0 + s, W)
        ref[a0 - r.r0:a1 - r.r0, b0 - r.c0:b1 - r.c0] = R[a0:a1, b0:b1, :3]
        bad_a += int(not np.array_equal(a, ref))
        wp = np.asarray(Image.open(w.png).convert("RGB")); ws = int(w.side_px)
        wr0, wc0 = int(w.r0_roi), int(w.c0_roi)
        y0, x0 = max(r.r0, wr0), max(r.c0, wc0); y1, x1 = min(r.r0 + s, wr0 + ws), min(r.c0 + s, wc0 + ws)
        if y1 > y0 and x1 > x0:
            nb += 1
            bad_b += int(not np.array_equal(a[y0 - r.r0:y1 - r.r0, x0 - r.c0:x1 - r.c0], wp[y0 - wr0:y1 - wr0, x0 - wc0:x1 - wc0]))
    ck("C-2a ritagli = ROI", bad_a == 0, f"{len(seen)} immagini, diverse {bad_a}")
    ck("C-2b ritagli = finestra R2b (parte interna)", bad_b == 0 and nb == len(seen), f"{nb} confrontate, diverse {bad_b}")
    ck("C-2c marcatori -> punti R2b", err <= 1e-9, f"max err {err:.2e} µm")
    ck("C-2d marcatori nel ritaglio", bool(cr.inside.all()), f"{int(cr.inside.sum())}/{len(cr)}")
    html = (OUT / f"R2c_annotator_{rater}.html").read_text().lower()
    n_kw = len(re.findall(r"\bcatch\s*\(", html))
    html = re.sub(r"\bcatch\s*\(", "(", html)      # parola chiave JavaScript try/catch, non un dato (dev notes)
    hits = {b: len(re.findall(re.escape(b), html)) for b in BANNED}
    hits = {k: v for k, v in hits.items() if v}
    wid = re.findall(r"a[1-6]_w[1-4]", html)
    ck("C-2e cecita' (grep)", not hits and not wid, f"parole trovate {hits}, finestre {len(wid)}, catch( JS esclusi {n_kw}")
    out = RES / "R2c_c2.csv"
    old = pd.read_csv(out) if out.exists() else pd.DataFrame()
    if len(old):
        old = old[old.rater != rater]
    pd.concat([old, pd.DataFrame(res)]).to_csv(out, index=False)
    return all(r["result"] == "PASS" for r in res)


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
