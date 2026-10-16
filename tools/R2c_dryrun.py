#!/usr/bin/env python
"""tools/R2c_dryrun.py — R2c, prova a secco dell'analisi (BL-125a) su risposte SINTETICHE in /tmp/r2c_dry (mai in results/).
Verita' piantata: 'touch' = contatto StarDist della coppia (se non appaiata: d < 0.5 d_max); ripetute identiche; trappole e
fascia 'separate'; Federica = Luca con il 10 % delle coppie invertite (seme fisso). Attese: CP-2 kappa = 1, CP-3 = 1,
D-3 PASS, CP-1 kappa ~0.7-0.85, B-dec che preferisce StarDist dove i segmentatori si separano.
Uso: .venv/bin/python tools/R2c_dryrun.py"""
import os, json, shutil, subprocess, sys
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"; D = Path("/tmp/r2c_dry")
shutil.rmtree(D, ignore_errors=True); (D / "annotations").mkdir(parents=True)
for f in ["R2c_items_Luca.csv", "R2c_items_Federica.csv", "R2c_pairs.csv", "R2c_focal.csv", "R2c_seg_focal.csv", "R2c_k2_windows.csv", "R2c_dmax.csv"]:
    shutil.copy(RES / f, D / f)
pr = pd.read_csv(RES / "R2c_pairs.csv").set_index("pair_id"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype").d_max_um
def truth(pid, kind):
    if kind in ("catch", "band"):
        return "separate"
    r = pr.loc[pid]
    if r.both_stardist_he == True:
        return "touch" if r.contact_stardist_he == True else "separate"
    return "touch" if r.d_um < 0.5 * dm[r.archetype] else "separate"
rng = np.random.default_rng(7)
flip = {"touch": "separate", "separate": "touch"}
for rater in ["Luca", "Federica"]:
    it = pd.read_csv(RES / f"R2c_items_{rater}.csv")
    ans = [truth(p, k) for p, k in zip(it.pair_id, it.kind)]
    if rater == "Federica":
        ans = [flip[a] if (k == "pair" and rng.random() < 0.10) else a for a, k in zip(ans, it.kind)]
        orig = dict((p, a) for p, a, k in zip(it.pair_id, ans, it.kind) if k == "pair")
        ans = [orig[p] if k == "repeat" else a for p, a, k in zip(it.pair_id, ans, it.kind)]
    doc = dict(tool_version="dry", doc_id=f"R2c_{rater}_v10", rater=rater, exported_at="dry", n_items=len(it), n_done=len(it),
               answers=[dict(item=int(i), answer=a, t_ms=1500) for i, a in zip(it["item"], ans)])
    json.dump(doc, open(D / f"annotations/R2c_answers_{rater}.json", "w"))
env = dict(os.environ, R2C_RES=str(D))
p = subprocess.run([sys.executable, str(ROOT / "tools/R2c_analysis.py")], env=env, capture_output=True, text=True)
print(p.stdout[-9000:]); print(p.stderr[-3000:]); sys.exit(p.returncode)
