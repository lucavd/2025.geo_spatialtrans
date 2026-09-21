"""tools/R2b_import_model.py — R2b: converte gli output dei modelli visione-linguaggio nel formato dell'annotatore HTML.

Input : uno o piu' file/cartelle JSON prodotti con il prompt results/R2b/R2b_model_prompt.md. Ogni oggetto: {"image": "A4_w2.png",
        "n_nuclei": int|null, "points": [[x_frac, y_frac], ...], "unreadable_fraction": float, "confidence": str, "notes": str}.
        Sono accettati: un array di oggetti, un oggetto singolo, o un oggetto {"windows": [...]}; anche JSON annidati in testo (```json ... ```).
Output: results/R2b/annotations/R2b_annotations_model_<nome>.json nel formato dell'annotatore (points in px immagine, margin_px = 0),
        rater = "model:<nome>"; results/R2b/R2b_model_import_log.csv con, per immagine, n_nuclei dichiarato vs len(points) e anomalie.
Uso: .venv/bin/python tools/R2b_import_model.py --model <nome> <file_o_cartella> [...]
Regole: se len(points) != n_nuclei si tiene n_nuclei come conteggio (colonna n_declared) e i punti come posizioni (appaiamento solo
        se presenti); coordinate fuori [0, 1] vengono ritagliate al bordo e segnalate; "n_nuclei": null -> finestra non valutata.
"""
import sys, json, re, glob, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT

RES = ROOT / "results/R2b"; ANN = RES / "annotations"; ANN.mkdir(exist_ok=True)

def load_objects(path):
    p = Path(path)
    files = sorted(p.rglob("*.json")) if p.is_dir() else [p]
    objs = []
    for f in files:
        txt = f.read_text()
        try: data = json.loads(txt)
        except json.JSONDecodeError:
            m = re.search(r"\{.*\}|\[.*\]", txt, re.S); data = json.loads(m.group(0)) if m else None
        if data is None: print("non parsabile:", f); continue
        if isinstance(data, dict) and "windows" in data: data = data["windows"]
        if isinstance(data, dict): data = [data]
        for o in data: o["_file"] = f.name
        objs += data
    return objs

def main(model, paths):
    wins = pd.read_csv(RES / "R2b_windows.csv").set_index("win_id")
    objs = []
    for p in paths: objs += load_objects(p)
    out, log = [], []
    for wid, w in wins.iterrows():
        cands = [o for o in objs if str(o.get("image", "")).split("/")[-1].replace(".png", "") == wid]
        if not cands:
            log.append(dict(win_id=wid, status="missing")); continue
        o = cands[-1]                                              # l'ultimo file vince
        n_decl = o.get("n_nuclei"); pts = o.get("points") or []
        side = int(w.side_px); arr = np.array(pts, float).reshape(-1, 2) if len(pts) else np.zeros((0, 2))
        oob = int(((arr < 0) | (arr > 1)).any(axis=1).sum()) if len(arr) else 0
        arr = np.clip(arr, 0, 0.999999) * side                    # frazioni -> px immagine (margin 0)
        rec = dict(win_id=wid, archetype=w.archetype, roi_id=w.roi_id, side_px=side, margin_px=0, um_per_px=float(w.um_per_px),
                   x0_um=float(w.x0_um), y0_um=float(w.y0_um), img_w=side, img_h=side,
                   points=arr.round(2).tolist(), polygons=[], exclusions=[], done=n_decl is not None,
                   note=str(o.get("notes", ""))[:200], n_points_in_window=len(arr),
                   n_declared=n_decl, unreadable_fraction=o.get("unreadable_fraction"), confidence=o.get("confidence"), source=o.get("_file"))
        out.append(rec)
        log.append(dict(win_id=wid, status="ok" if n_decl is not None else "null", n_declared=n_decl, n_points=len(arr),
                        mismatch=(n_decl is not None and len(arr) and n_decl != len(arr)), out_of_bounds=oob,
                        unreadable_fraction=o.get("unreadable_fraction"), confidence=o.get("confidence"), n_candidates=len(cands), file=o.get("_file")))
    doc = dict(tool_version="R2b_import_model 1.0", doc_id="R2b_model", rater=f"model:{model}", exported_at=pd.Timestamp.utcnow().isoformat(), windows=out)
    f = ANN / f"R2b_annotations_model_{model}.json"; json.dump(doc, open(f, "w"))
    lg = pd.DataFrame(log); lg["model"] = model
    lf = RES / "R2b_model_import_log.csv"
    if lf.exists(): lg = pd.concat([pd.read_csv(lf).query("model != @model"), lg], ignore_index=True)
    lg.to_csv(lf, index=False)
    print(f, "\n", lg[["win_id", "status", "n_declared", "n_points", "mismatch", "out_of_bounds"]].to_string(index=False))

if __name__ == "__main__":
    a = sys.argv[1:]; model = a[a.index("--model") + 1]; paths = [x for i, x in enumerate(a) if x != "--model" and a[i - 1] != "--model"]
    main(model, paths)
