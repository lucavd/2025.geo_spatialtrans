#!/usr/bin/env python
"""tools/R2c_annotator.py — R2c: annotatore autonomo si'/no sui contatti fra nuclei (zip: index.html + img/*.png).

Per ogni item di results/R2c/R2c_items_<rater>.csv ritaglia dal ROI R2 un quadrato di lato fisso per archetipo
(R2c_dmax.csv: crop_um) centrato sul punto medio della coppia, PNG senza perdita, nome opaco (hash della coppia: le ripetute
riusano la stessa immagine). Due anelli sui punti di Luca R2b. Convenzione dei punti R2b: pixel k = [k, k+1) (R2b_annotator:
drawImage a 0,0), quindi x_crop = x_um / upp - c0 e il marcatore e' passato a toDisplay nel formato centro-pixel (x_crop - 0.5).
L'HTML non contiene il tipo di item, l'archetipo, la finestra ne' alcuna informazione dei segmentatori (C-2, grep).
Il nucleo logico e' tools/R2c_annotator_core.js (test: node tools/R2c_annotator_core.test.js).
Output: /mnt/micron/geo_spatialtrans/R2c/R2c_annotator_<rater>.zip (fuori git) + results/R2c/R2c_crops_<rater>.csv
Uso: .venv/bin/python tools/R2c_annotator.py Luca|Federica
"""
import sys, io, json, hashlib, zipfile
import numpy as np, pandas as pd
from pathlib import Path
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_roi

RES = ROOT / "results/R2c"; R2B = ROOT / "results/R2b"
OUT = Path("/mnt/micron/geo_spatialtrans/R2c"); OUT.mkdir(parents=True, exist_ok=True)
TOOL_VERSION = "R2c-annotator 1.0"
DISP = 520
SALT = "R2c-20261016"


def crop_rgb(rgb, c0, r0, side):
    """Ritaglio side x side con angolo (r0, c0); fuori dal ROI bianco."""
    H, W = rgb.shape[:2]
    out = np.full((side, side, 3), 255, np.uint8)
    a0, b0 = max(r0, 0), max(c0, 0); a1, b1 = min(r0 + side, H), min(c0 + side, W)
    if a1 > a0 and b1 > b0:
        out[a0 - r0:a1 - r0, b0 - c0:b1 - c0] = rgb[a0:a1, b0:b1, :3]
    return out


def crop_geometry(xa, ya, xb, yb, side):
    mx, my = (xa + xb) / 2, (ya + yb) / 2
    return int(np.floor(mx - side / 2 + 0.5)), int(np.floor(my - side / 2 + 0.5))


HTML = r"""<!doctype html><html lang="it"><head><meta charset="utf-8"><title>R2c — contatti fra nuclei</title>
<style>
body{margin:0;font-family:system-ui,sans-serif;background:#1e1e1e;color:#eee}
#wrap{display:flex;gap:18px;padding:14px}
#left{width:__DISP__px}
canvas{background:#fff;display:block;border:1px solid #555}
#right{max-width:520px;font-size:14px;line-height:1.45}
button{font-size:15px;margin:4px 4px 4px 0;padding:7px 12px;cursor:pointer}
.k{display:inline-block;min-width:16px;padding:0 4px;border:1px solid #888;border-radius:3px;text-align:center;font-family:monospace}
#status{margin:8px 0;font-size:15px}
#cur{color:#ffd54f}
</style></head><body><div id="wrap">
<div id="left"><canvas id="cv" width="__DISP__" height="__DISP__"></canvas>
<div id="status"></div>
<div><button id="b1">1 · Si toccano</button><button id="b2">2 · Separati</button><button id="b3">3 · Non giudicabile</button></div>
<div><button id="bb">&larr; Indietro</button><button id="bm">M · anelli</button><button id="be">Esporta JSON</button>
<label style="font-size:13px;margin-left:6px">Riprendi da file: <input type="file" id="fi" accept=".json"></label></div></div>
<div id="right"><h3 style="margin-top:0">I due nuclei indicati dagli anelli azzurri si toccano?</h3>
<p><span class="k">1</span> <b>Si toccano</b>: fra i due nuclei non c'&egrave; una fessura chiara (citoplasma o fondo), anche se si sovrappongono.</p>
<p><span class="k">2</span> <b>Separati</b>: una fessura chiara li separa lungo tutto il bordo affacciato, oppure c'&egrave; un altro nucleo in mezzo.</p>
<p><span class="k">3</span> <b>Non giudicabile</b>: l'immagine non permette di dirlo, oppure un anello non sta su un nucleo.</p>
<p><span class="k">&larr;</span> (Backspace) torna all'immagine precedente per correggere &middot; <span class="k">M</span> mostra/nasconde gli anelli.</p>
<p>Gli anelli stanno dove avevi cliccato in R2b. La barra bianca in basso misura 5 &micro;m.
Il lavoro si salva da solo in questo browser: puoi chiudere e riprendere quando vuoi. Alla fine (o a ogni pausa) premi
<b>Esporta JSON</b> e invia il file.</p></div></div>
<script>__CORE__</script>
<script>
"use strict";
const DOC = __DOC__;
const C = window.R2cCore;
const KEY = "R2c_" + DOC.doc_id;
const cv = document.getElementById("cv"), ctx = cv.getContext("2d"), D = cv.width;
let st = C.newState(DOC.items.length), showMk = true;
try { const s = JSON.parse(localStorage.getItem(KEY) || "null"); if (s) C.restore(st, s); } catch (e) {}
const cache = {};
function img(src){ if (!cache[src]) { const im = new Image(); im.onload = draw; im.src = src; cache[src] = im; } return cache[src]; }
function save(){ localStorage.setItem(KEY, JSON.stringify({n: st.n, ans: st.ans, t_ms: st.t_ms})); }
const LAB = {touch: "si toccano", separate: "separati", unsure: "non giudicabile"};
function draw(){
  const it = DOC.items[st.cur], im = img(it.src);
  ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, D, D);
  if (im.complete && im.naturalWidth) { ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high"; ctx.drawImage(im, 0, 0, D, D); }
  if (showMk) { ctx.strokeStyle = "#00e5ff"; ctx.lineWidth = 1.6;
    it.mk.forEach(p => { const q = C.toDisplay(p[0], p[1], it.side, D); ctx.beginPath(); ctx.arc(q[0], q[1], 9, 0, 2*Math.PI); ctx.stroke(); }); }
  const L = 5 / it.upp * D / it.side; ctx.fillStyle = "#000"; ctx.fillRect(10, D - 20, L + 4, 10); ctx.fillStyle = "#fff"; ctx.fillRect(12, D - 18, L, 6);
  const a = st.ans[st.cur];
  document.getElementById("status").innerHTML = "Immagine <b>" + (st.cur + 1) + "</b> / " + st.n + " &middot; fatte " + C.progress(st) +
    (a ? " &middot; <span id='cur'>risposta attuale: " + LAB[a] + "</span>" : "");
  if (st.cur + 1 < st.n) img(DOC.items[st.cur + 1].src);
}
function act(k){ C.answer(st, k, Date.now()); save(); draw(); }
document.addEventListener("keydown", e => {
  if (e.target.tagName === "INPUT") return;
  if (e.key === "1" || e.key === "2" || e.key === "3") { act(e.key); e.preventDefault(); }
  else if (e.key === "Backspace" || e.key === "ArrowLeft") { C.back(st, Date.now()); draw(); e.preventDefault(); }
  else if (e.key === "m" || e.key === "M") { showMk = !showMk; draw(); }
});
document.getElementById("b1").onclick = () => act("1");
document.getElementById("b2").onclick = () => act("2");
document.getElementById("b3").onclick = () => act("3");
document.getElementById("bb").onclick = () => { C.back(st, Date.now()); draw(); };
document.getElementById("bm").onclick = () => { showMk = !showMk; draw(); };
document.getElementById("be").onclick = () => {
  const doc = C.exportDoc(st, {tool_version: DOC.tool_version, doc_id: DOC.doc_id, rater: DOC.rater, now_iso: new Date().toISOString()});
  const b = new Blob([JSON.stringify(doc)], {type: "application/json"}), u = URL.createObjectURL(b), a = document.createElement("a");
  a.href = u; a.download = "R2c_answers_" + DOC.rater + ".json"; a.click(); URL.revokeObjectURL(u);
};
document.getElementById("fi").onchange = e => { const f = e.target.files[0]; if (!f) return; const r = new FileReader();
  r.onload = () => { try { const o = JSON.parse(r.result); if (o.doc_id !== DOC.doc_id) { alert("File di un altro annotatore/documento"); return; }
    C.restore(st, {n: o.n_items, ans: o.answers.map(x => x.answer), t_ms: o.answers.map(x => x.t_ms)}); save(); draw(); } catch (err) { alert("File non valido"); } };
  r.readAsText(f); };
C.show(st, Date.now()); draw();
</script></body></html>
"""


def build(rater):
    items = pd.read_csv(RES / f"R2c_items_{rater}.csv")
    wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
    dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"].set_index(["win_id", "idx"])
    rgb = {}; imgs = {}; doc_items = []; rows = []
    for it in items.itertuples():
        w = wins.loc[it.win_id]; A = w.archetype; upp = float(w.um_per_px)
        side = int(round(float(dm.crop_um[A]) / upp))
        pa = pts.loc[(it.win_id, it.i)]; pb = pts.loc[(it.win_id, it.j)]
        xa, ya, xb, yb = pa.x_um / upp, pa.y_um / upp, pb.x_um / upp, pb.y_um / upp
        c0, r0 = crop_geometry(xa, ya, xb, yb, side)
        key = f"{it.win_id}:{it.i}-{it.j}"
        name = "img/i_" + hashlib.sha1((SALT + key).encode()).hexdigest()[:12] + ".png"
        if name not in imgs:
            k = (A, w.roi_id)
            if k not in rgb:
                rgb[k] = load_roi(A, w.roi_id)
            arr = crop_rgb(rgb[k], c0, r0, side)
            buf = io.BytesIO(); Image.fromarray(arr).save(buf, format="PNG", compress_level=9); imgs[name] = buf.getvalue()
        mk = [[xa - c0 - 0.5, ya - r0 - 0.5], [xb - c0 - 0.5, yb - r0 - 0.5]]
        doc_items.append(dict(src=name, side=side, upp=upp, mk=mk))
        rows.append(dict(item=it.item, win_id=it.win_id, i=it.i, j=it.j, kind=it.kind, img=name, c0=c0, r0=r0, side_px=side,
                         mk_ax=mk[0][0], mk_ay=mk[0][1], mk_bx=mk[1][0], mk_by=mk[1][1],
                         inside=bool(0 <= mk[0][0] < side and 0 <= mk[0][1] < side and 0 <= mk[1][0] < side and 0 <= mk[1][1] < side)))
    doc = dict(tool_version=TOOL_VERSION, doc_id=f"R2c_{rater}_v10", rater=rater, items=doc_items)
    core = (ROOT / "tools/R2c_annotator_core.js").read_text()
    html = HTML.replace("__DISP__", str(DISP)).replace("__CORE__", core).replace("__DOC__", json.dumps(doc))
    zp = OUT / f"R2c_annotator_{rater}.zip"
    with zipfile.ZipFile(str(zp) + ".tmp", "w", zipfile.ZIP_STORED) as z:
        z.writestr(f"R2c_{rater}/index.html", html)
        for n, b in sorted(imgs.items()):
            z.writestr(f"R2c_{rater}/{n}", b)
    Path(str(zp) + ".tmp").replace(zp)
    (OUT / f"R2c_annotator_{rater}.html").write_text(html)            # copia per i controlli di cecita' (grep)
    pd.DataFrame(rows).to_csv(RES / f"R2c_crops_{rater}.csv", index=False)
    print(rater, "items", len(doc_items), "immagini", len(imgs), "zip MB", round(zp.stat().st_size / 1e6, 1),
          "marcatori nel ritaglio", int(pd.DataFrame(rows).inside.sum()), "/", len(rows))


if __name__ == "__main__":
    build(sys.argv[1])
