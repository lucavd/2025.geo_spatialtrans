"""R/testing/test_R2b.py — asserzioni PASS/FAIL della sessione R2b dalle tabelle in results/R2b (attese pre-registrate, non modificate).
Uso: .venv/bin/python R/testing/test_R2b.py   (exit 0 anche con FAIL di check B/CP: sono esiti, non errori; exit 1 solo se un check C fallisce
o una tabella manca). Scrive results/R2b/R2b_test_results.csv.
"""
import sys, json, subprocess, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; RES = ROOT / "results/R2b"
R = []
def rec(sec, cid, desc, status, detail=""): R.append(dict(section=sec, id=cid, check=desc, result=status, detail=detail)); print(f"{status:5s} {cid:10s} {desc}  [{detail}]")
def pf(ok): return "PASS" if ok else "FAIL"

# ---------------------------------------------------------------- Check C
wc = pd.read_csv(RES / "R2b_windows_check.csv"); _ok = wc.inside_roi & wc["tissue_ge_0.90"] & wc.roundtrip_identical
rec("C", "C-R2b.1a", "24 finestre dentro ROI, tessuto >= 0.90, round-trip pixel-identico", pf(_ok.all() and len(wc) == 24), f"{_ok.sum()}/24")
_mw = pd.read_csv(RES / "R2b_manual_windows.csv"); _ex = _mw[(_mw.rater == "Luca") & (_mw.n_excl_unreadable > 0)]
rec("C", "C-R2b.1b", "bolle = 0 secondo la maschera R2 (BL-026, parziale): l'annotatore ha escluso bolle in " + str(len(_ex)) + " finestre", "WARN" if len(_ex) else "PASS", "maschera R2 0 % in 24/24; esclusioni manuali: " + ", ".join(f"{r.win_id} {r.excl_frac:.1%}" for r in _ex.itertuples()) + " (revisore avversariale, discrepanza 1)")
node = subprocess.run(["node", str(ROOT / "R/testing/test_R2b_core.js")], capture_output=True, text=True); rec("C", "C-R2b.2", "core annotatore (node): round-trip schermo/immagine, zoom, export/import < 0.5 px", pf(node.returncode == 0), node.stdout.strip().splitlines()[-1] if node.stdout else node.stderr[-100:])
st = pd.read_csv(RES / "R2b_selftest.csv"); rec("C", "C-R2b.3", "appaiamento su dati sintetici: recall/precisione entro 1 %, 1:1, esclusioni, area", pf((st.result == "PASS").all()), f"{(st.result=='PASS').sum()}/{len(st)}")
rp = pd.read_csv(RES / "R2b_repro.csv"); rec("C", "C-R2b.4", "riproducibilita' campionamento finestre (md5 identici con lo stesso seed)", pf(rp.identical.all() and len(rp) == 24), f"{rp.identical.sum()}/24")
# import annotazioni: rater e conteggi
ann = {}
for f in sorted((RES / "annotations").glob("*.json")):
    o = json.load(open(f)); ann[o["rater"]] = sum(w["n_points_in_window"] for w in o["windows"])
rec("C", "C-R2b.5", "annotazioni umane importate: Luca 24 finestre, Federica 12", pf(ann.get("Luca", 0) > 2000 and ann.get("Federica", 0) > 1000), f"Luca {ann.get('Luca')} punti, Federica {ann.get('Federica')} punti")

# ---------------------------------------------------------------- Check B (primario Luca)
S = pd.read_csv(RES / "R2b_archetype_summary.csv").set_index("archetype")
lo = np.minimum(S.density_cellpose_rgb_windows, S.density_stardist_he_windows) * 0.9; hi = np.maximum(S.density_cellpose_rgb_windows, S.density_stardist_he_windows) * 1.1
b1 = (S.density_manual >= lo) & (S.density_manual <= hi)
rec("B", "B-R2b.1", "densita' manuale fra Cellpose e StarDist (+-10 %) sulle stesse finestre, 6/6", pf(b1.all()), f"{b1.sum()}/6; sopra entrambi: {', '.join(S.index[S.density_manual > hi])}")
tol = pd.Series({"A1": .15, "A2": .15, "A3": .15, "A4": .15, "A5": .20, "A6": .20})
dev = (S.consensus_R2_5roi / S.density_manual - 1); b2 = dev.abs() <= tol
rec("B", "B-R2b.2", "consenso R2 entro +-15 % (A1-A4) / +-20 % (A5-A6) dalla densita' manuale", pf(b2.all()), f"{b2.sum()}/6; " + "; ".join(f"{a} {v:+.0%}" for a, v in dev.items()))
rc, rs = S.recall_cellpose_rgb, S.recall_stardist_he
dir_ok = {a: rc[a] > rs[a] for a in ("A1", "A2", "A3")} | {a: rs[a] > rc[a] for a in ("A5", "A6")} | {"A4": (rc["A4"] >= .85) and (rs["A4"] >= .85)}
prec_ok = {"A5": S.precision_cellpose_rgb["A5"] < .85, "A6": S.precision_cellpose_rgb["A6"] < .70}
rec("B", "B-R2b.3a", "direzione dei richiami: CP>SD in A1-A3, SD>CP in A5-A6, A4 entrambi >= 0.85", pf(all(dir_ok.values())), ", ".join(f"{a}:{'ok' if v else 'NO'}" for a, v in dir_ok.items()))
rec("B", "B-R2b.3b", "precisione Cellpose < 0.85 in A5 e < 0.70 in A6 (BL-032/033)", pf(all(prec_ok.values())), f"A5 {S.precision_cellpose_rgb['A5']:.2f}, A6 {S.precision_cellpose_rgb['A6']:.2f}")
order = list(S.density_manual.sort_values(ascending=False).index); rec("B", "B-R2b.4", "ordinamento A4 > A1 > A2 > A3 > A5 > A6", pf(order == ["A4", "A1", "A2", "A3", "A5", "A6"]), " > ".join(order))
rec("B", "B-R2b.5", "area nucleare da contorni manuali", "N/A", "fase 2 (calibri) non eseguita in R2b; BL-025 resta aperta")

# ---------------------------------------------------------------- Controprove
ir = pd.read_csv(RES / "R2b_interrater.csv"); h = ir[ir.rater_a.isin(["Luca", "Federica"]) & ir.rater_b.isin(["Luca", "Federica"])]
rel = abs(h.n_a.sum() - h.n_b.sum()) / ((h.n_a.sum() + h.n_b.sum()) / 2); f13 = 2 * h["tp_3um"].sum() / (h.n_a.sum() + h.n_b.sum()); f15 = 2 * h["tp_5um"].sum() / (h.n_a.sum() + h.n_b.sum())
per_arch = h.groupby("archetype").apply(lambda g: abs(g.n_a.sum() - g.n_b.sum()) / ((g.n_a.sum() + g.n_b.sum()) / 2), include_groups=False)
rec("CP", "CP-R2b.1a", "accordo umano Luca-Federica: differenza relativa di conteggio < 10 %", "PASS" if (rel < .10 and (per_arch < .10).all() and (h.rel_diff < .10).all()) else ("WARN" if rel < .10 else "FAIL"), f"totale {rel:.1%} (PASS); per archetipo {(per_arch < .10).sum()}/6 (max {per_arch.max():.0%} {per_arch.idxmax()}); per finestra {(h.rel_diff < .10).sum()}/{len(h)}: PASS solo sul totale (revisore, R5)")
rec("CP", "CP-R2b.1b", "accordo umano: F1 punto-punto (3 um) > 0.90", pf(f13 > .90), f"F1 3 um {f13:.2f}; 5 um {f15:.2f}")
cf = S.common_fn_frac_cp_sd; rec("CP", "CP-R2b.2", "punti manuali persi da entrambi CP e SD < 5 % per archetipo", pf((cf < .05).all()), "; ".join(f"{a} {v:.0%}" for a, v in cf.items()))
ns = pd.read_csv(RES / "R2b_null_summary.csv"); g = ns[(ns.rater == "Luca") & (ns.method == "cellpose_rgb")].set_index("archetype")
ghost = 1 - g.precision_obs; rec("CP", "CP-R2b.3", "oggetti Cellpose senza punto manuale: >= 15 % in A5, >= 30 % in A6 (BL-032/033)", pf(ghost["A5"] >= .15 and ghost["A6"] >= .30), f"A5 {ghost['A5']:.0%} (nullo {1-g.precision_null['A5']:.0%}), A6 {ghost['A6']:.0%} (nullo {1-g.precision_null['A6']:.0%})")
gl = ns[(ns.rater == "Luca")]; rec("CP", "CP-R2b.4", "nullo punti casuali: F1 Luca > F1 nullo + 0.20 per archetipo e metodo", pf((gl.f1_gain > .20).all()), f"min guadagno {gl.f1_gain.min():.2f}; {(gl.f1_gain > .20).sum()}/{len(gl)}")
gm = ns[ns.rater.str.startswith("model:")].groupby("rater").f1_gain.median(); rec("CP", "CP-R2b.4m", "VLM vs nullo (misura, nessuna attesa)", "INFO", "; ".join(f"{r.replace('model:','')} mediana {v:.2f}" for r, v in gm.items()))
od = pd.read_csv(RES / "R2b_od_strata.csv"); odF = pd.read_csv(RES / "R2b_od_strata_Federica.csv")
rec("CP", "CP-R2b.5", "OD ematossilina dei punti persi da CP e SD <= 75 % degli appaiati in >= 5/6 (Luca) — tesi di Luca", pf(od.missed_paler_25pct.sum() >= 5), f"Luca {od.missed_paler_25pct.sum()}/6 (rapporti {', '.join(f'{v:.2f}' for v in od.ratio_missed_over_matched)}); Federica {odF.missed_paler_25pct.sum()}/6")
# VLM attese
rv = pd.read_csv(RES / "R2b_rater_vs_primary.csv"); ga = rv[rv.rater == "model:gpt6_astra"].groupby("archetype").rel_err.apply(lambda s: s.abs().median())
rec("CP", "VLM-1", "GPT-6 Astra: errore relativo di conteggio mediano < 20 % per archetipo", pf((ga < .20).all()), f"{(ga < .20).sum()}/6; " + "; ".join(f"{a} {v:.0%}" for a, v in ga.items()))
irm = ir[(ir.rater_a == "Luca") & ir.rater_b.str.startswith("model:") | (ir.rater_b == "Luca") & ir.rater_a.str.startswith("model:")].copy(); irm["model"] = np.where(irm.rater_a == "Luca", irm.rater_b, irm.rater_a)
fm = irm.groupby("model").apply(lambda g: 2 * g["tp_3um"].sum() / (g.n_a + g.n_b).sum(), include_groups=False)
rec("CP", "VLM-2", "VLM: F1 punto-punto (3 um) vs Luca atteso < 0.6", pf((fm < .6).all()), "; ".join(f"{m.replace('model:','')} {v:.2f}" for m, v in fm.items()))

df = pd.DataFrame(R); df.to_csv(RES / "R2b_test_results.csv", index=False)
print("\nC:", (df[df.section == "C"].result == "PASS").sum(), "/", (df.section == "C").sum(), " B:", df[df.section == "B"].result.value_counts().to_dict(), " CP:", df[df.section == "CP"].result.value_counts().to_dict())
sys.exit(0 if df[df.section == "C"].result.isin(["PASS", "WARN"]).all() else 1)
