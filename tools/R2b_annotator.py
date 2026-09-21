"""tools/R2b_annotator.py — R2b: genera l'annotatore HTML autonomo per la verita' manuale (punti, contorni, esclusioni).

Per ogni finestra di results/R2b/R2b_windows.csv ritaglia dal ROI R2 la finestra + margine di contesto (MARGIN_UM), incorpora il PNG
in base64 e 10 bersagli (croci) seedati per la scelta non selettiva dei nuclei da contornare (BL-025). L'annotatore non mostra alcuna
segmentazione (conteggio alla cieca). Il nucleo geometrico e' in tools/R2b_annotator_core.js (testato con node, C-R2b.2).
Output: /mnt/micron/geo_spatialtrans/R2b/R2b_annotator_<rater>.html (fuori git; artifact in chat) + results/R2b/R2b_targets.csv
Uso: .venv/bin/python tools/R2b_annotator.py            # tutte le 24 finestre (annotatore principale)
     .venv/bin/python tools/R2b_annotator.py --rater2   # solo le 12 finestre del secondo annotatore (CP-R2b.1)
"""
import sys, io, json, base64, hashlib, numpy as np, pandas as pd
from pathlib import Path
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_roi

OUT_DIR = Path("/mnt/micron/geo_spatialtrans/R2b"); OUT_DIR.mkdir(parents=True, exist_ok=True)
RES = ROOT / "results/R2b"
MARGIN_UM = 10.0          # contesto attorno alla finestra (sbiancato nell'annotatore)
N_TARGETS = 10            # croci per finestra -> 40 contorni per archetipo
TARGET_INSET_UM = 5.0
SEED = 20260921

def targets_for(win_id, side_px, upp, rng):
    """10 bersagli uniformi dentro la finestra (inset 5 um), separati di almeno side/5 (rigetto)."""
    inset = TARGET_INSET_UM / upp; dmin = side_px / 5
    pts = []
    for _ in range(20000):
        p = rng.uniform(inset, side_px - inset, 2)
        if all(np.hypot(*(p - q)) >= dmin for q in pts):
            pts.append(p)
            if len(pts) == N_TARGETS: break
    assert len(pts) == N_TARGETS, win_id
    return np.array(pts)

def build_doc(rater2=False):
    win = pd.read_csv(RES / "R2b_windows.csv")
    if rater2: win = win[win.rater2].copy()
    doc = {"doc_id": "R2b_rater2" if rater2 else "R2b_main", "tool_version": "R2b-annotator 1.0", "windows": []}
    trows = []
    for _, r in win.iterrows():
        upp = float(r.um_per_px); m = int(round(MARGIN_UM / upp)); side = int(r.side_px)
        rgb = load_roi(r.archetype, r.roi_id)
        c0, r0 = int(r.c0_roi) - m, int(r.r0_roi) - m
        assert c0 >= 0 and r0 >= 0 and c0 + side + 2 * m <= rgb.shape[1] and r0 + side + 2 * m <= rgb.shape[0], r.win_id
        img = rgb[r0:r0 + side + 2 * m, c0:c0 + side + 2 * m]
        # controllo: il centro del ritaglio con margine coincide con la finestra salvata da R2b_windows.py
        core = np.asarray(Image.open(r.png)); assert np.array_equal(core, img[m:m + side, m:m + side]), r.win_id
        buf = io.BytesIO(); Image.fromarray(img).save(buf, format="PNG", compress_level=6)
        rng = np.random.default_rng(SEED + int(hashlib.md5(r.win_id.encode()).hexdigest()[:8], 16) % 100000)
        tg = targets_for(r.win_id, side, upp, rng) + m       # in px immagine (con margine)
        for k, (x, y) in enumerate(tg, 1):
            trows.append(dict(win_id=r.win_id, target=k, x_img_px=x, y_img_px=y, x_win_px=x - m, y_win_px=y - m,
                              x_um=r.x0_um + (x - m) * upp, y_um=r.y0_um + (y - m) * upp))
        doc["windows"].append(dict(win_id=r.win_id, archetype=r.archetype, roi_id=r.roi_id, side_px=side, side_um=float(r.side_um),
                                   margin_px=m, um_per_px=upp, x0_um=float(r.x0_um), y0_um=float(r.y0_um),
                                   img_w=int(img.shape[1]), img_h=int(img.shape[0]),
                                   targets=[[float(x), float(y)] for x, y in tg],
                                   png="data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()))
    if not rater2:
        pd.DataFrame(trows).to_csv(RES / "R2b_targets.csv", index=False)
    return doc

HTML = r"""<!DOCTYPE html>
<html lang="it"><head><meta charset="utf-8"><title>R2b — annotatore nuclei (__DOCID__)</title>
<style>
 html,body{margin:0;height:100%;font-family:system-ui,sans-serif;font-size:13px;background:#222;color:#eee}
 #side{position:absolute;left:0;top:0;bottom:0;width:300px;overflow-y:auto;background:#2b2b2b;padding:10px;box-sizing:border-box;border-right:1px solid #444}
 #main{position:absolute;left:300px;top:0;right:0;bottom:0}
 canvas{display:block;cursor:crosshair}
 button{margin:2px;padding:4px 8px;background:#444;color:#eee;border:1px solid #666;border-radius:3px;cursor:pointer}
 button.active{background:#c94;color:#000;font-weight:bold}
 .win{display:block;width:100%;text-align:left;margin:1px 0;padding:3px 6px}
 .win.cur{outline:2px solid #fc3}
 .win.done{color:#8f8}
 #status{position:absolute;left:310px;bottom:6px;background:rgba(0,0,0,.6);padding:3px 8px;border-radius:3px;pointer-events:none}
 h3{margin:8px 0 4px 0;font-size:14px;color:#fc3} details{margin:6px 0} summary{cursor:pointer;color:#fc3}
 input[type=text],textarea{width:100%;box-sizing:border-box;background:#333;color:#eee;border:1px solid #666}
 .small{color:#aaa;font-size:11px} kbd{background:#444;padding:0 4px;border-radius:3px}
</style></head><body>
<div id="side">
 <h3>R2b — conteggio nuclei</h3>
 <div class="small">Documento: __DOCID__ · __NWIN__ finestre · v1.0</div>
 <label>Annotatore: <input type="text" id="rater" placeholder="nome"></label>
 <h3>Modalità</h3>
 <button id="m_points">1 Nuclei</button><button id="m_poly">2 Contorni</button><button id="m_excl">3 Esclusione</button>
 <div id="excl_opts" style="display:none">tipo: <select id="excl_type"><option value="unreadable">illeggibile (bolla/artefatto)</option><option value="off_archetype">fuori archetipo</option></select></div>
 <div><button id="b_undo">Annulla ultimo (Z)</button><button id="b_close">Chiudi poligono (Invio)</button><button id="b_fit">Adatta (F)</button></div>
 <div><label><input type="checkbox" id="done"> finestra completata</label></div>
 <textarea id="note" rows="2" placeholder="note sulla finestra (opzionale)"></textarea>
 <h3>Finestre</h3><div id="winlist"></div>
 <h3>File</h3>
 <button id="b_export">Esporta JSON</button> <label class="small">Importa: <input type="file" id="f_import" accept=".json"></label>
 <div class="small">Il lavoro si salva automaticamente nel browser (localStorage). Esporta comunque il JSON a ogni pausa e inviamelo in chat.</div>
 <details open><summary>Istruzioni</summary>
 <ol class="small">
  <li><b>Riquadro giallo</b> = finestra di conteggio. Fuori dal riquadro c'è solo contesto (sbiancato): non si annota.</li>
  <li><b>Nuclei (1)</b>: un click sul centro di ogni nucleo il cui centro cade dentro il riquadro. Conta ogni profilo nucleare colorato dall'ematossilina, anche pallido, piccolo o allungato; non contare citoplasma, globuli rossi, pigmento, frammenti senza cromatina. Click su un punto già segnato lo rimuove. I nuclei dentro le zone escluse si possono saltare (in analisi vengono ignorati).</li>
  <li><b>Contorni (2)</b>: per ciascuna delle 10 <b>croci azzurre</b> traccia il contorno del nucleo più vicino alla croce (click sui vertici, doppio click o <kbd>Invio</kbd> per chiudere, <kbd>Esc</kbd> annulla). Se entro ~10 µm dalla croce non c'è un nucleo, salta la croce. Contorna il bordo esterno della cromatina.</li>
  <li><b>Esclusione (3)</b>: contorna le zone non giudicabili (bolle, sfocato) come <i>illeggibile</i> e quelle che non appartengono all'archetipo (es. ghiandola tumorale in una finestra di stroma) come <i>fuori archetipo</i>.</li>
  <li>Navigazione: rotella = zoom sul cursore; trascina con tasto destro/centrale o <kbd>spazio</kbd>+trascina = sposta; <kbd>F</kbd> adatta; <kbd>N</kbd>/<kbd>P</kbd> finestra successiva/precedente; tasto destro dentro un poligono chiuso lo elimina.</li>
  <li>Segna <i>finestra completata</i> quando hai finito nuclei e contorni. Alla fine: <b>Esporta JSON</b>.</li>
 </ol></details>
</div>
<div id="main"><canvas id="cv"></canvas></div>
<div id="status"></div>
<script id="core">__CORE__</script>
<script>
const DOC = __DOC__;
const KEY = "r2b_" + DOC.doc_id;
let cur = 0, mode = "points", view = {z:1,px:0,py:0}, draft = [], hover = null, spaceDown = false;
let state = DOC.windows.map(() => ({points:[], polygons:[], exclusions:[], done:false, note:""}));
const imgs = DOC.windows.map(w => { const im = new Image(); im.src = w.png; im.onload = draw; return im; });
const cv = document.getElementById("cv"), ctx = cv.getContext("2d");
const $ = id => document.getElementById(id);
// ---- persistenza
function save(){ localStorage.setItem(KEY, JSON.stringify({rater: $("rater").value, state: state})); refreshList(); }
function load(){ const s = localStorage.getItem(KEY); if(!s) return; try{ const o = JSON.parse(s); if(o.state && o.state.length===state.length) state = o.state; $("rater").value = o.rater||""; }catch(e){} }
// ---- vista
function resize(){ cv.width = $("main").clientWidth; cv.height = $("main").clientHeight; draw(); }
function fit(){ const w = DOC.windows[cur]; view = fitView(w.img_w, w.img_h, cv.width, cv.height); draw(); }
function setWin(i){ cur = Math.max(0, Math.min(DOC.windows.length-1, i)); draft = []; fit(); syncSide(); }
function setMode(m){ mode = m; draft = []; ["points","poly","excl"].forEach(k => $("m_"+k).classList.toggle("active", k===m)); $("excl_opts").style.display = m==="excl" ? "block":"none"; draw(); }
function syncSide(){ const a = state[cur]; $("done").checked = a.done; $("note").value = a.note; refreshList(); }
function refreshList(){
  const L = $("winlist"); L.innerHTML = "";
  DOC.windows.forEach((w,i) => { const a = state[i]; const b = document.createElement("button"); b.className = "win"+(i===cur?" cur":"")+(a.done?" done":"");
    b.textContent = `${w.win_id} · ${w.side_um} µm · ${countInside(a, w.margin_px, w.side_px)} nuclei · ${a.polygons.length}/10 cont.` + (a.exclusions.length? ` · ${a.exclusions.length} escl.`:"") + (a.done?" ✓":"");
    b.onclick = () => setWin(i); L.appendChild(b); });
}
// ---- disegno
function draw(){
  const w = DOC.windows[cur], a = state[cur], im = imgs[cur];
  ctx.setTransform(1,0,0,1,0,0); ctx.fillStyle = "#222"; ctx.fillRect(0,0,cv.width,cv.height);
  ctx.setTransform(view.z,0,0,view.z,view.px,view.py); ctx.imageSmoothingEnabled = true;
  if (im.complete) ctx.drawImage(im, 0, 0);
  const m = w.margin_px, s = w.side_px, W = w.img_w, H = w.img_h;
  ctx.fillStyle = "rgba(255,255,255,0.55)";
  ctx.fillRect(0,0,W,m); ctx.fillRect(0,m+s,W,H-m-s); ctx.fillRect(0,m,m,s); ctx.fillRect(m+s,m,W-m-s,s);
  const lw = 1.5/view.z;
  ctx.lineWidth = lw; ctx.strokeStyle = "#ff0"; ctx.strokeRect(m, m, s, s);
  // esclusioni
  a.exclusions.forEach(e => { const u = e.type==="unreadable"; drawPoly(e.pts, u ? "rgba(255,140,0,0.9)" : "rgba(200,0,255,0.9)", u ? "rgba(255,140,0,0.15)" : "rgba(200,0,255,0.12)"); });
  // contorni
  a.polygons.forEach(pg => drawPoly(pg, "#0f0", "rgba(0,255,0,0.12)"));
  // bersagli
  const r = 7/view.z;
  ctx.strokeStyle = "#0ff"; ctx.lineWidth = lw;
  w.targets.forEach((t,k) => { ctx.beginPath(); ctx.moveTo(t[0]-r,t[1]); ctx.lineTo(t[0]+r,t[1]); ctx.moveTo(t[0],t[1]-r); ctx.lineTo(t[0],t[1]+r); ctx.stroke();
    ctx.fillStyle = "#0ff"; ctx.font = `${11/view.z}px sans-serif`; ctx.fillText(String(k+1), t[0]+r*0.8, t[1]-r*0.8); });
  // punti
  a.points.forEach(p => { ctx.beginPath(); ctx.arc(p.x, p.y, 4/view.z, 0, 2*Math.PI); ctx.fillStyle = insideWindow(p,m,s) ? "rgba(255,0,0,0.85)" : "rgba(255,0,0,0.3)"; ctx.fill(); ctx.strokeStyle = "#fff"; ctx.lineWidth = 1/view.z; ctx.stroke(); });
  // poligono in corso
  if (draft.length){ ctx.strokeStyle = mode==="poly" ? "#0f0" : "#f80"; ctx.lineWidth = lw; ctx.setLineDash([4/view.z, 3/view.z]); ctx.beginPath(); draft.forEach((p,i) => i? ctx.lineTo(p.x,p.y) : ctx.moveTo(p.x,p.y)); if (hover) ctx.lineTo(hover.x, hover.y); ctx.stroke(); ctx.setLineDash([]);
    draft.forEach(p => { ctx.beginPath(); ctx.arc(p.x,p.y,3/view.z,0,2*Math.PI); ctx.fillStyle = "#fff"; ctx.fill(); }); }
  $("status").textContent = `${w.win_id} (${w.archetype}, ROI ${w.roi_id}) · lato ${w.side_um} µm · zoom ${view.z.toFixed(2)}× · ` + (hover ? `x ${((hover.x-m)*w.um_per_px).toFixed(1)} µm, y ${((hover.y-m)*w.um_per_px).toFixed(1)} µm · ` : "") + `modalità: ${mode==="points"?"nuclei":mode==="poly"?"contorni":"esclusione"}`;
}
function drawPoly(pts, stroke, fill){ if (pts.length<2) return; ctx.beginPath(); pts.forEach((p,i) => i? ctx.lineTo(p.x,p.y) : ctx.moveTo(p.x,p.y)); ctx.closePath(); ctx.strokeStyle = stroke; ctx.lineWidth = 1.5/view.z; ctx.stroke(); if (fill){ ctx.fillStyle = fill; ctx.fill(); } }
// ---- interazione
let dragging = false, dragStart = null, moved = false, panBtn = false;
cv.addEventListener("mousedown", e => { dragStart = {x:e.offsetX, y:e.offsetY, px:view.px, py:view.py}; moved = false; panBtn = (e.button!==0) || spaceDown; dragging = true; });
cv.addEventListener("mousemove", e => { hover = toImage(e.offsetX, e.offsetY, view);
  if (dragging && panBtn){ view.px = dragStart.px + (e.offsetX - dragStart.x); view.py = dragStart.py + (e.offsetY - dragStart.y); }
  if (dragging && (Math.abs(e.offsetX-dragStart.x)>3 || Math.abs(e.offsetY-dragStart.y)>3)) moved = true; draw(); });
cv.addEventListener("mouseup", e => { if (dragging && !moved && !panBtn && e.button===0) click(toImage(e.offsetX, e.offsetY, view)); if (dragging && !moved && e.button===2) rightClick(toImage(e.offsetX, e.offsetY, view)); dragging = false; });
cv.addEventListener("mouseleave", () => { dragging = false; hover = null; draw(); });
cv.addEventListener("contextmenu", e => e.preventDefault());
cv.addEventListener("wheel", e => { e.preventDefault(); view = zoomAt(view, e.offsetX, e.offsetY, Math.pow(1.1, -e.deltaY/100), 0.1, 40); draw(); }, {passive:false});
cv.addEventListener("dblclick", e => { if (mode!=="points"){ if (draft.length>=2){ const l = draft[draft.length-1], q = draft[draft.length-2]; if (Math.hypot(l.x-q.x,l.y-q.y) < 2/view.z) draft.pop(); } closePoly(); } });
function click(p){ const a = state[cur];
  if (mode==="points"){ const i = nearestIndex(p, a.points, 6/view.z); if (i>=0) a.points.splice(i,1); else a.points.push(p); save(); }
  else draft.push(p); draw(); }
function rightClick(p){ const a = state[cur]; if (draft.length){ draft = []; draw(); return; }
  if (mode==="poly"){ const i = a.polygons.findIndex(pg => pointInPoly(p, pg)); if (i>=0){ a.polygons.splice(i,1); save(); } }
  if (mode==="excl"){ const i = a.exclusions.findIndex(e => pointInPoly(p, e.pts)); if (i>=0){ a.exclusions.splice(i,1); save(); } } draw(); }
function closePoly(){ if (draft.length<3){ return; } const a = state[cur];
  if (mode==="poly") a.polygons.push(draft); else if (mode==="excl") a.exclusions.push({type: $("excl_type").value, pts: draft});
  draft = []; save(); draw(); }
function undo(){ const a = state[cur]; if (draft.length){ draft.pop(); } else if (mode==="points") a.points.pop(); else if (mode==="poly") a.polygons.pop(); else a.exclusions.pop(); save(); draw(); }
document.addEventListener("keydown", e => { if (e.target.tagName==="INPUT" || e.target.tagName==="TEXTAREA") return;
  if (e.code==="Space"){ spaceDown = true; e.preventDefault(); }
  else if (e.key==="1") setMode("points"); else if (e.key==="2") setMode("poly"); else if (e.key==="3") setMode("excl");
  else if (e.key==="z"||e.key==="Z") undo(); else if (e.key==="Enter") closePoly(); else if (e.key==="Escape"){ draft=[]; draw(); }
  else if (e.key==="f"||e.key==="F") fit(); else if (e.key==="n"||e.key==="N") setWin(cur+1); else if (e.key==="p"||e.key==="P") setWin(cur-1); });
document.addEventListener("keyup", e => { if (e.code==="Space") spaceDown = false; });
$("m_points").onclick = () => setMode("points"); $("m_poly").onclick = () => setMode("poly"); $("m_excl").onclick = () => setMode("excl");
$("b_undo").onclick = undo; $("b_close").onclick = closePoly; $("b_fit").onclick = fit;
$("done").onchange = () => { state[cur].done = $("done").checked; save(); };
$("note").onchange = () => { state[cur].note = $("note").value; save(); };
$("rater").onchange = save;
$("b_export").onclick = () => { const rater = $("rater").value.trim() || "anonimo"; const obj = buildExport(DOC, rater, state, new Date().toISOString());
  const blob = new Blob([JSON.stringify(obj)], {type:"application/json"}); const a = document.createElement("a"); a.href = URL.createObjectURL(blob);
  a.download = `R2b_annotations_${DOC.doc_id}_${rater.replace(/\s+/g,"_")}_${new Date().toISOString().slice(0,10)}.json`; a.click(); };
$("f_import").onchange = e => { const f = e.target.files[0]; if(!f) return; const rd = new FileReader(); rd.onload = () => { try { const r = parseImport(DOC, JSON.parse(rd.result)); state = r.state; if (r.rater) $("rater").value = r.rater; save(); syncSide(); draw(); alert(`Importate ${r.n_windows} finestre`); } catch(err){ alert("File non valido: "+err); } }; rd.readAsText(f); };
window.addEventListener("resize", resize);
load(); resize(); setMode("points"); setWin(0);
</script></body></html>
"""

def main():
    rater2 = "--rater2" in sys.argv
    doc = build_doc(rater2)
    core = (Path(__file__).resolve().parent / "R2b_annotator_core.js").read_text()
    html = (HTML.replace("__DOCID__", doc["doc_id"]).replace("__NWIN__", str(len(doc["windows"])))
                .replace("__CORE__", core).replace("__DOC__", json.dumps(doc)))
    out = OUT_DIR / f"R2b_annotator_{'rater2' if rater2 else 'main'}.html"
    out.write_text(html)
    print(out, f"{out.stat().st_size/1e6:.1f} MB", len(doc["windows"]), "finestre")

if __name__ == "__main__":
    main()
