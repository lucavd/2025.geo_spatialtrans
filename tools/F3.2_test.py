#!/usr/bin/env python3
"""tools/F3.2_test.py - F3.2, batteria dei cancelli G-1...G-4: stampa PASS/FAIL per asserzione; esce con 1 se una fallisce.
Uso (radice del repository): python3 tools/F3.2_test.py [G-1|G-2|G-3|G-4|all]
G-4 (d) della prima regola e' registrato come fallito come scritto (12/30): la batteria lo stampa come tale, non lo nasconde.
"""
import sys, os, csv, json, math, subprocess, tempfile, shutil
from pathlib import Path
import pandas as pd, numpy as np
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/F3.2"; FAIL = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok: FAIL.append(name)

def compare_values(a, b):
    """True se le due tabelle del censimento hanno le stesse righe e gli stessi valori (relativo 1e-9; md5 per le tabelle)"""
    k = ["family", "fun", "archetype", "parameter"]
    if len(a) != len(b): return False
    m = a.merge(b, on=k, how="outer", suffixes=("_a", "_b"), indicator=True)
    if (m["_merge"] != "both").any(): return False
    for _, r in m.iterrows():
        if r.kind_a == "scalar":
            if r.value_a in ("NA", "") or r.value_b in ("NA", ""):
                if r.value_a != r.value_b: return False
                continue
            x, y = float(r.value_a), float(r.value_b)
            if not (x == y or abs(x - y) <= 1e-9 * max(abs(x), abs(y))): return False
        elif r.kind_a == "table":
            if r.values_md5_a != r.values_md5_b or str(r.n_values_a) != str(r.n_values_b): return False
        elif str(r.value_a) != str(r.value_b): return False
    return True

def g1():
    ref = pd.read_csv(OUT / "F3.2_census_values.csv", dtype=str, keep_default_na=False)
    tmp = Path(tempfile.mkdtemp()) / "v.csv"
    r = subprocess.run(["Rscript", "--vanilla", "tools/F3.2_census.R", str(tmp)], cwd=ROOT, capture_output=True, text=True)
    check("G-1 censimento riletto dal pacchetto senza errori", r.returncode == 0 and tmp.exists(), r.stderr[-200:])
    new = pd.read_csv(tmp, dtype=str, keep_default_na=False)
    check("G-1 valori del censimento uguali al pacchetto (1e-9 relativo; md5 delle tabelle)", compare_values(ref, new), f"{len(ref)} righe")
    cen = pd.read_csv(OUT / "F3.2_census.csv", dtype=str, keep_default_na=False)
    check("G-1 nessun valore senza insieme di aree", (cen.area_set != "UNKNOWN").all() and len(cen) == len(ref), f"{(cen.area_set=='UNKNOWN').sum()} senza")
    check("G-1 nessun valore senza sezione", (cen.section != "").all())
    check("G-1 nessun valore senza fonte", (cen.source_id.str.len() > 0).all(), f"{(cen.source_id=='').sum()} senza fonte")
    m1 = ref.drop(index=ref.index[5]); m2 = ref.copy()
    i = m2.index[(m2.kind == "scalar") & ~m2.value.isin(["NA", ""])][7]; m2.loc[i, "value"] = str(float(m2.loc[i, "value"]) * 1.01)
    check("G-1 mutante 1 (riga tolta) rilevato", not compare_values(m1, new))
    check("G-1 mutante 2 (valore cambiato dell'1 %) rilevato", not compare_values(m2, new))

def g2():
    man = list(csv.DictReader(open(ROOT / "data/holdout/F3.2_download_manifest.csv")))
    st = {r["dest"]: r for r in csv.DictReader(open(ROOT / "data/holdout/F3.2_download_status.csv"))}
    check("G-2 (a) ogni file del manifesto scaricato e verificato", all(m["dest"] in st and st[m["dest"]]["status"] == "OK" for m in man),
          f"{sum(m['dest'] in st and st[m['dest']]['status']=='OK' for m in man)}/{len(man)}")
    n_md5 = n_et = 0; bad = []
    for m in man:
        s = st.get(m["dest"])
        if not s: continue
        if m["md5_published"]:
            n_md5 += 1
            if s["md5_observed"] != m["md5_published"]: bad.append(m["filename"])
        else:
            n_et += 1
            if not (s["etag_recomputed"] == m["etag"] == s["etag_at_end"]): bad.append(m["filename"])
        if not Path(m["dest"]).exists() or Path(m["dest"]).stat().st_size != int(s["bytes_observed"]): bad.append(m["filename"] + " (dimensione)")
    check("G-2 (a) md5 uguale al pubblicato o ETag multiparte ricalcolato uguale", not bad, f"md5 {n_md5}, ETag {n_et}, eccezioni {bad[:3]}")
    b = OUT / "F3.2_md5_before.txt"; a = OUT / "F3.2_md5_after.txt"
    if b.exists() and a.exists():
        B = dict(l.split()[::-1] for l in open(b) if l.strip()); A_ = dict(l.split()[::-1] for l in open(a) if l.strip())
        check("G-2 (b) md5 dei dati di taratura letti uguali prima e dopo", B == A_ and len(B) > 0, f"{len(B)} file")
    else:
        check("G-2 (b) md5 dei dati di taratura letti uguali prima e dopo", False, "file md5 mancanti")

def g3():
    cal = ["A1", "A23", "A4", "A5", "A6"]; C = OUT / "compat"
    res = {l: json.load(open(C / f"{l}.json")) for l in cal if (C / f"{l}.json").exists()}
    check("G-3 PASS su tutte le sezioni di taratura (5; la pre-registrazione dice 6, errore dichiarato)", len(res) == 5 and all(r["PASS"] for r in res.values()),
          ", ".join(f"{k}:{'PASS' if v['PASS'] else 'FAIL'}" for k, v in res.items()))
    for f in ["scale10", "rot90", "trunc", "trunc_mid"]:
        det = [l for l in cal if (C / f"{l}_{f}.json").exists() and not json.load(open(C / f"{l}_{f}.json"))["PASS"]]
        check(f"G-3 guasto {f} rilevato", len(det) >= 1, f"rilevato in {len(det)}/5 sezioni: {det}")

def gap(a, b):
    dx = max(b[0] - a[2], a[0] - b[2], 0.0); dy = max(b[1] - a[3], a[1] - b[3], 0.0); return math.hypot(dx, dy)

def g4():
    u = pd.read_csv(OUT / "F3.2_used_areas.csv")
    n_r2 = len(pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")); w = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv"); n_r5 = len(pd.read_csv(ROOT / "results/R5/R5_tiles.csv"))
    cnt = u.groupby("set").size().to_dict()
    check("G-4 (a) mappa: stesse righe dei file sorgente", cnt.get("R2_ROI") == n_r2 and cnt.get("R2b_WINDOW") == len(w) and cnt.get("R5_EDGE_TILE") == n_r5,
          f"R2 {cnt.get('R2_ROI')}/{n_r2}, R2b {cnt.get('R2b_WINDOW')}/{len(w)}, R5 {cnt.get('R5_EDGE_TILE')}/{n_r5}")
    check("G-4 (a) finestre di R2b dentro le aree di R2", bool(w.inside_roi.all()))
    files = sorted(OUT.glob("F3.2_selected_*.csv"))
    if files:
        mins = []; pair = []
        for f in files:
            d = pd.read_csv(f); lab = d.label.iloc[0]
            if lab not in ("A1", "A23", "A4", "A5", "A6"): continue
            used = [tuple(v) for v in u[(u.section == lab) & u.set.isin(["R2_ROI", "R2b_WINDOW"])][["x0", "y0", "x1", "y1"]].values]   # DEC-094
            ar = d[d.level == "L1"]; boxes = [tuple(v) for v in ar[["x0", "y0", "x1", "y1"]].values]
            test = [tuple(v) for v in d[d.level == "test_tile"][["x0", "y0", "x1", "y1"]].values]
            for i, b in enumerate(boxes):
                mins.append(min(gap(b, x) for x in used + test))
                pair += [gap(b, c) for c in boxes[i + 1:]]
        check("G-4 (b) aree di L1 a >= 280 um dalle aree di R2, dalle finestre e dalla tessera di prova (DEC-094)", bool(mins) and min(mins) >= 280 - 1e-6, f"minimo {min(mins):.1f} um su {len(mins)} aree" if mins else "nessuna area")
        check("G-4 (b) aree di L1 a >= 280 um fra loro", (not pair) or min(pair) >= 280 - 1e-6, f"minimo {min(pair):.1f} um" if pair else "")
    else:
        check("G-4 (b) aree di L1 estratte", False, "nessun file F3.2_areas2_*.csv")
    r = subprocess.run([sys.executable, "tools/F3.2_guard.py"], cwd=ROOT, capture_output=True, text=True)
    check("G-4 (c) nessuna occorrenza dei dati di verifica fuori dai lettori ammessi", r.returncode == 0, r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr[-200:])
    tmp = Path(tempfile.mkdtemp()); shutil.copy(ROOT / "tools/S2.0_stripes.py", tmp / "S2.0_stripes.py")
    with open(tmp / "S2.0_stripes.py", "a") as f: f.write('\nX = "/mnt/micron/geo_spatialtrans/data_holdout/CRC_P1_SR300"\n')
    r = subprocess.run([sys.executable, "tools/F3.2_guard.py", str(tmp)], cwd=ROOT, capture_output=True, text=True)
    check("G-4 (c) mutante (percorso piantato in uno script di taratura) rilevato", r.returncode == 1, r.stdout.strip().splitlines()[-1] if r.stdout else "")
    v1 = pd.read_csv(OUT / "F3.2_rule_validation.csv"); r2 = v1[v1.kind == "R2_ROI"]
    check("G-4 (d) prima regola: classe giusta in >= 28 aree di R2 su 30 [fallito come scritto, DEC-092]", int(r2.ok.sum()) >= 28, f"{int(r2.ok.sum())}/30")
    f2 = OUT / "F3.2_rule2_validation.csv"
    if f2.exists():
        v = pd.read_csv(f2); r = v[v.kind == "R2_ROI"]; cr = v[v.kind == "R1_crop"].iloc[0]
        check("G-4 (d') (i) seconda regola: classe propria >= 28/30 (riferimento senza l'area)", int(r.own_ok.sum()) >= 28, f"{int(r.own_ok.sum())}/30")
        check("G-4 (d') (ii) negativi scartati >= 14/15", int(r.negatives_rejected.sum()) >= 14 and int(r.negatives.sum()) == 15, f"{int(r.negatives_rejected.sum())}/{int(r.negatives.sum())}")
        check("G-4 (d') (iii) copertura e artefatti ammettono >= 28/30 aree di R2", int(r.admitted.sum()) >= 28, f"{int(r.admitted.sum())}/30")
        check("G-4 (d') (iv) ritaglio del linfonodo di R1 scartato", not bool(cr.admitted), f"copertura {cr.coverage:.3f}, grigio {cr.grey_frac:.3f}")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name, fn in [("G-1", g1), ("G-2", g2), ("G-3", g3), ("G-4", g4)]:
        if which in ("all", name): fn()
    print(f"\n{'TUTTO PASS' if not FAIL else 'FALLITI: ' + ', '.join(FAIL)}"); sys.exit(1 if FAIL else 0)
