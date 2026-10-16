/* tools/R2c_annotator_core.js — R2c: nucleo dell'annotatore (stato, risposte, tempi, esportazione, ripristino).
   Nessun accesso al DOM: testato con node da tools/R2c_annotator_core.test.js (C-2). */
(function (root) {
  "use strict";
  var KEYS = { "1": "touch", "2": "separate", "3": "unsure" };
  function newState(n) {
    var a = [], t = [];
    for (var i = 0; i < n; i++) { a.push(null); t.push(null); }
    return { n: n, cur: 0, ans: a, t_ms: t, shown_at: null };
  }
  function show(st, now) { st.shown_at = now; return st; }
  function answer(st, key, now) {
    var v = KEYS[key];
    if (!v || st.n === 0) return st;
    st.ans[st.cur] = v;
    st.t_ms[st.cur] = (st.shown_at === null) ? null : Math.max(0, now - st.shown_at);
    if (st.cur < st.n - 1) st.cur += 1;
    st.shown_at = now;
    return st;
  }
  function back(st, now) { if (st.cur > 0) st.cur -= 1; st.shown_at = now; return st; }
  function firstEmpty(st) { var k = st.ans.indexOf(null); return k < 0 ? st.n - 1 : k; }
  function progress(st) { var c = 0; for (var i = 0; i < st.n; i++) if (st.ans[i] !== null) c++; return c; }
  function exportDoc(st, meta) {
    var out = [];
    for (var i = 0; i < st.n; i++) out.push({ item: i + 1, answer: st.ans[i], t_ms: st.t_ms[i] });
    return { tool_version: meta.tool_version, doc_id: meta.doc_id, rater: meta.rater, exported_at: meta.now_iso,
             n_items: st.n, n_done: progress(st), answers: out };
  }
  function restore(st, saved) {
    if (!saved || saved.n !== st.n || !saved.ans || saved.ans.length !== st.n) return st;
    st.ans = saved.ans.slice(); st.t_ms = saved.t_ms.slice(); st.cur = firstEmpty(st);
    return st;
  }
  /* coordinate del marcatore: punto in pixel del ritaglio (centro del pixel = intero) -> pixel del canvas */
  function toDisplay(x_px, y_px, side_px, disp) { var s = disp / side_px; return [(x_px + 0.5) * s, (y_px + 0.5) * s]; }
  var api = { KEYS: KEYS, newState: newState, show: show, answer: answer, back: back, firstEmpty: firstEmpty,
              progress: progress, exportDoc: exportDoc, restore: restore, toDisplay: toDisplay };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.R2cCore = api;
})(this);
