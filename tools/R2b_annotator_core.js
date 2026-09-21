// tools/R2b_annotator_core.js — funzioni pure dell'annotatore R2b (testate con node in R/testing/test_R2b_core.js).
// Coordinate: "img" = pixel dell'immagine mostrata (finestra + margine di contesto), origine in alto a sinistra, asse y verso il basso.
// "win" = pixel della finestra (img - margin). "um" = micrometri nel sistema del ROI R2 (x0_um + win_px * um_per_px).
const R2B_TOOL_VERSION = "R2b-annotator 1.0";

// Trasformazione schermo <-> immagine con vista {z: zoom, px, py: traslazione in pixel schermo}.
function toImage(sx, sy, view) { return { x: (sx - view.px) / view.z, y: (sy - view.py) / view.z }; }
function toScreen(ix, iy, view) { return { x: ix * view.z + view.px, y: iy * view.z + view.py }; }
// Zoom attorno a un punto schermo (ancora fissa): il punto immagine sotto il cursore non si muove.
function zoomAt(view, sx, sy, factor, zmin, zmax) {
  const z = Math.min(zmax, Math.max(zmin, view.z * factor));
  const f = z / view.z;
  return { z: z, px: sx - (sx - view.px) * f, py: sy - (sy - view.py) * f };
}
function fitView(imgW, imgH, cw, ch) {
  const z = Math.min(cw / imgW, ch / imgH);
  return { z: z, px: (cw - imgW * z) / 2, py: (ch - imgH * z) / 2 };
}
// Punto (img px) dentro la finestra [margin, margin+side) ?
function insideWindow(p, margin, side) {
  return p.x >= margin && p.x < margin + side && p.y >= margin && p.y < margin + side;
}
// Ray casting; vertici [{x,y}]; bordo incluso in modo approssimato.
function pointInPoly(p, poly) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const xi = poly[i].x, yi = poly[i].y, xj = poly[j].x, yj = poly[j].y;
    const cross = ((yi > p.y) !== (yj > p.y)) && (p.x < (xj - xi) * (p.y - yi) / (yj - yi) + xi);
    if (cross) inside = !inside;
  }
  return inside;
}
// Area (shoelace) in px^2, positiva.
function polyArea(poly) {
  let a = 0;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) a += (poly[j].x + poly[i].x) * (poly[j].y - poly[i].y);
  return Math.abs(a) / 2;
}
function polyCentroid(poly) {
  let cx = 0, cy = 0, a = 0;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const f = poly[j].x * poly[i].y - poly[i].x * poly[j].y; a += f; cx += (poly[j].x + poly[i].x) * f; cy += (poly[j].y + poly[i].y) * f;
  }
  if (Math.abs(a) < 1e-9) { const n = poly.length; return { x: poly.reduce((s, p) => s + p.x, 0) / n, y: poly.reduce((s, p) => s + p.y, 0) / n }; }
  return { x: cx / (3 * a), y: cy / (3 * a) };
}
// Indice del punto entro tol (px immagine) dal punto p, altrimenti -1.
function nearestIndex(p, pts, tol) {
  let best = -1, bd = tol * tol;
  for (let i = 0; i < pts.length; i++) { const d = (pts[i].x - p.x) ** 2 + (pts[i].y - p.y) ** 2; if (d < bd) { bd = d; best = i; } }
  return best;
}
function countInside(ann, margin, side) {
  return ann.points.filter(p => insideWindow(p, margin, side)).length;
}
// Esportazione: un oggetto serializzabile con tutte le coordinate in px immagine (float) e i metadati per convertirle.
function buildExport(doc, rater, state, now) {
  return {
    tool_version: R2B_TOOL_VERSION, doc_id: doc.doc_id, rater: rater, exported_at: now,
    windows: doc.windows.map((w, i) => {
      const a = state[i] || { points: [], polygons: [], exclusions: [], done: false, note: "" };
      return {
        win_id: w.win_id, archetype: w.archetype, roi_id: w.roi_id, side_px: w.side_px, margin_px: w.margin_px, um_per_px: w.um_per_px,
        x0_um: w.x0_um, y0_um: w.y0_um, img_w: w.img_w, img_h: w.img_h,
        points: a.points.map(p => [p.x, p.y]),
        polygons: a.polygons.map(pg => pg.map(p => [p.x, p.y])),
        exclusions: a.exclusions.map(e => ({ type: e.type, pts: e.pts.map(p => [p.x, p.y]) })),
        done: !!a.done, note: a.note || "",
        n_points_in_window: countInside(a, w.margin_px, w.side_px)
      };
    })
  };
}
// Importazione: dal JSON esportato allo stato interno; ignora finestre non presenti nel documento.
function parseImport(doc, obj) {
  const idx = {}; doc.windows.forEach((w, i) => idx[w.win_id] = i);
  const state = doc.windows.map(() => ({ points: [], polygons: [], exclusions: [], done: false, note: "" }));
  let n = 0;
  (obj.windows || []).forEach(w => {
    const i = idx[w.win_id]; if (i === undefined) return; n++;
    state[i] = {
      points: (w.points || []).map(q => ({ x: q[0], y: q[1] })),
      polygons: (w.polygons || []).map(pg => pg.map(q => ({ x: q[0], y: q[1] }))),
      exclusions: (w.exclusions || []).map(e => ({ type: e.type, pts: e.pts.map(q => ({ x: q[0], y: q[1] })) })),
      done: !!w.done, note: w.note || ""
    };
  });
  return { state: state, n_windows: n, rater: obj.rater || "" };
}
if (typeof module !== "undefined") module.exports = { R2B_TOOL_VERSION, toImage, toScreen, zoomAt, fitView, insideWindow, pointInPoly, polyArea, polyCentroid, nearestIndex, countInside, buildExport, parseImport };
