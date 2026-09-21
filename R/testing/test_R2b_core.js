// R/testing/test_R2b_core.js — check C-R2b.2: round-trip geometrico ed esportazione/importazione dell'annotatore R2b.
// Uso: node R/testing/test_R2b_core.js   (stampa PASS/FAIL per asserzione, exit 1 se una fallisce)
const core = require("../../tools/R2b_annotator_core.js");
let fails = 0;
function check(name, ok, detail) { console.log((ok ? "PASS" : "FAIL") + "  " + name + (detail ? "  [" + detail + "]" : "")); if (!ok) fails++; }
// PRNG deterministico (LCG) per riproducibilita'
let seed = 20260921; function rnd() { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; }

// 1. schermo <-> immagine: 10 000 punti casuali, viste casuali
let maxErr = 0;
for (let k = 0; k < 10000; k++) {
  const view = { z: 0.2 + rnd() * 12, px: (rnd() - 0.5) * 4000, py: (rnd() - 0.5) * 4000 };
  const ix = rnd() * 1000, iy = rnd() * 1000;
  const s = core.toScreen(ix, iy, view), b = core.toImage(s.x, s.y, view);
  maxErr = Math.max(maxErr, Math.abs(b.x - ix), Math.abs(b.y - iy));
}
check("C-R2b.2a round-trip schermo->immagine->schermo < 1e-6 px", maxErr < 1e-6, "max " + maxErr.toExponential(2));

// 2. zoom ancorato: il punto immagine sotto il cursore non si muove
let maxShift = 0;
for (let k = 0; k < 2000; k++) {
  const view = { z: 0.5 + rnd() * 5, px: rnd() * 500, py: rnd() * 500 };
  const sx = rnd() * 1200, sy = rnd() * 800;
  const before = core.toImage(sx, sy, view);
  const v2 = core.zoomAt(view, sx, sy, rnd() < 0.5 ? 1.1 : 1 / 1.1, 0.1, 40);
  const after = core.toImage(sx, sy, v2);
  maxShift = Math.max(maxShift, Math.abs(after.x - before.x), Math.abs(after.y - before.y));
}
check("C-R2b.2b zoom ancorato al cursore, spostamento < 1e-6 px", maxShift < 1e-6, "max " + maxShift.toExponential(2));

// 3. fitView centra e contiene
const fv = core.fitView(987, 987, 1400, 900);
check("C-R2b.2c fitView: immagine contenuta e centrata", Math.abs(fv.z - 900 / 987) < 1e-12 && Math.abs(fv.py) < 1e-9 && Math.abs(fv.px - (1400 - 900) / 2) < 1e-9);

// 4. geometria poligoni
const sq = [{ x: 10, y: 10 }, { x: 20, y: 10 }, { x: 20, y: 20 }, { x: 10, y: 20 }];
check("C-R2b.2d pointInPoly dentro/fuori quadrato", core.pointInPoly({ x: 15, y: 15 }, sq) && !core.pointInPoly({ x: 25, y: 15 }, sq) && !core.pointInPoly({ x: 15, y: 5 }, sq));
check("C-R2b.2e polyArea quadrato 10x10 = 100", Math.abs(core.polyArea(sq) - 100) < 1e-9);
const tri = [{ x: 0, y: 0 }, { x: 4, y: 0 }, { x: 0, y: 3 }];
check("C-R2b.2f polyArea triangolo 3-4 = 6 (orientamento inverso incluso)", Math.abs(core.polyArea(tri) - 6) < 1e-9 && Math.abs(core.polyArea(tri.slice().reverse()) - 6) < 1e-9);
const cen = core.polyCentroid(sq);
check("C-R2b.2g polyCentroid quadrato = (15,15)", Math.abs(cen.x - 15) < 1e-9 && Math.abs(cen.y - 15) < 1e-9);
check("C-R2b.2h insideWindow: bordi [margin, margin+side)", core.insideWindow({ x: 37, y: 37 }, 37, 365) && !core.insideWindow({ x: 36.99, y: 100 }, 37, 365) && core.insideWindow({ x: 401.99, y: 100 }, 37, 365) && !core.insideWindow({ x: 402, y: 100 }, 37, 365));
check("C-R2b.2i nearestIndex entro tolleranza", core.nearestIndex({ x: 5, y: 5 }, [{ x: 0, y: 0 }, { x: 5.5, y: 5.2 }], 1) === 1 && core.nearestIndex({ x: 50, y: 50 }, [{ x: 0, y: 0 }], 1) === -1);

// 5. esportazione -> JSON -> importazione: stato casuale su un documento sintetico di 3 finestre
const doc = { doc_id: "test", windows: [1, 2, 3].map(i => ({ win_id: "T_w" + i, archetype: "T", roi_id: "r" + i, side_px: 365, margin_px: 37, um_per_px: 0.27376, x0_um: 100 * i, y0_um: 50 * i, img_w: 439, img_h: 439 })) };
const state = doc.windows.map(() => ({
  points: Array.from({ length: 50 }, () => ({ x: rnd() * 439, y: rnd() * 439 })),
  polygons: Array.from({ length: 5 }, () => Array.from({ length: 8 }, () => ({ x: rnd() * 439, y: rnd() * 439 }))),
  exclusions: [{ type: "bubble", pts: Array.from({ length: 6 }, () => ({ x: rnd() * 439, y: rnd() * 439 })) }],
  done: true, note: "n"
}));
const exp = core.buildExport(doc, "tester", state, "2026-09-21T00:00:00Z");
const back = core.parseImport(doc, JSON.parse(JSON.stringify(exp)));
let maxD = 0, nPts = 0, nPoly = 0, nExc = 0;
back.state.forEach((a, i) => {
  a.points.forEach((p, k) => { maxD = Math.max(maxD, Math.abs(p.x - state[i].points[k].x), Math.abs(p.y - state[i].points[k].y)); nPts++; });
  a.polygons.forEach((pg, k) => { pg.forEach((p, m) => maxD = Math.max(maxD, Math.abs(p.x - state[i].polygons[k][m].x), Math.abs(p.y - state[i].polygons[k][m].y))); nPoly++; });
  a.exclusions.forEach((e, k) => { e.pts.forEach((p, m) => maxD = Math.max(maxD, Math.abs(p.x - state[i].exclusions[k].pts[m].x))); nExc += (e.type === "bubble") ? 1 : 0; });
});
check("C-R2b.2j export->JSON->import: 150 punti, 15 poligoni, 3 esclusioni, errore max < 0.5 px", nPts === 150 && nPoly === 15 && nExc === 3 && maxD < 0.5, "max " + maxD.toExponential(2));
check("C-R2b.2k export riporta n_points_in_window coerente con countInside", exp.windows.every((w, i) => w.n_points_in_window === core.countInside(state[i], 37, 365)));
check("C-R2b.2l export contiene i metadati di conversione (um_per_px, x0_um, y0_um, margin_px)", exp.windows.every(w => typeof w.um_per_px === "number" && typeof w.x0_um === "number" && typeof w.y0_um === "number" && w.margin_px === 37));

console.log(fails === 0 ? "ALL PASS" : fails + " FAIL");
process.exit(fails === 0 ? 0 : 1);
